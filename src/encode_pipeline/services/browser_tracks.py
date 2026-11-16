"""Run-scoped WashU track catalog, derivatives, and range reads."""

from __future__ import annotations

from collections.abc import Callable, Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path

from encode_pipeline.adapters.browser.catalog import (
    BrowserArtifactDescription,
    describe_browser_artifact,
)
from encode_pipeline.adapters.browser.scaffolds import viewer_scaffold
from encode_pipeline.platform.result_generations import (
    validate_artifact_generation,
    validate_artifact_revision,
)
from encode_pipeline.platform.runs import RunArtifactRef
from encode_pipeline.services.artifact_downloads import ArtifactDownloadService
from encode_pipeline.services.browser_bgzf import BrowserIndexError, query_indexed_bed
from encode_pipeline.services.browser_bytes import (
    BrowserByteError,
    ByteRange,
    iter_plan_bytes,
    parse_byte_range,
)
from encode_pipeline.services.browser_derivatives import (
    MAX_SOURCE_BYTES,
    DisplayDerivativeError,
    build_display_derivative,
)
from encode_pipeline.services.browser_store import (
    BrowserStoreError,
    StoredDerivative,
    iter_descriptor,
    open_derivative,
    read_derivative_summary,
    write_derivative,
)
from encode_pipeline.services.run_repositories import ResultGenerationChangedError
from encode_pipeline.services.runs import RunService


_ARTIFACT_SCAN_LIMIT = 500
_REGION_COUNT_LIMIT = 100_001
_DRAW_LIMIT = 100_000


class BrowserTrackError(Exception):
    """Stable browser-track failure without a host path or exception text."""

    def __init__(
        self,
        status: int,
        code: str,
        message: str,
        *,
        path: str,
        reason: str | None = None,
    ) -> None:
        super().__init__(code)
        self.status = status
        self.code = code
        self.message = message
        self.path = path
        self.reason = reason


@dataclass(frozen=True)
class BrowserByteBody:
    """Headers and a one-shot byte iterator. Closing is inside the iterator."""

    status: int
    media_type: str
    headers: dict[str, str]
    chunks: Iterator[bytes] | None
    close: Callable[[], None]


class BrowserTrackService:
    """Describe and read tracks for one authenticated run view."""

    def __init__(
        self,
        *,
        run_service: RunService,
        download_service: ArtifactDownloadService,
        workspace_root: Path,
        reference_profile_service: object,
    ) -> None:
        self._runs = run_service
        self._downloads = download_service
        self._root = workspace_root
        self._references = reference_profile_service

    def session(self, run_id: str) -> dict[str, object]:
        run = self._run(run_id)
        generation = self._generation(run_id)
        reference = self._reference(run_id, run.workflow_id)
        tracks: list[dict[str, object]] = []
        unavailable: list[dict[str, object]] = []
        if generation is not None:
            for artifact in self._artifacts(run_id):
                described = self._describe(run.workflow_id, artifact, reference)
                if described is None:
                    continue
                if described["available"] is not True:
                    unavailable.append(described["unavailable"])
                    continue
                tracks.append(
                    self._track_payload(
                        run_id,
                        artifact,
                        generation,
                        described["description"],
                        reference,
                    )
                )
        return {
            "ok": True,
            "run_id": run_id,
            "workflow_id": run.workflow_id,
            "artifact_generation": generation,
            "reference": reference,
            "annotation_status": (
                "available"
                if any(track["role"] == "annotation" for track in tracks)
                else "absent"
            ),
            "interaction_status": (
                "available"
                if any(track["role"] == "interaction" for track in tracks)
                else "absent"
            ),
            "tracks": tracks,
            "unavailable": unavailable,
            "draw_record_limit": _DRAW_LIMIT,
        }

    def prepare(
        self,
        run_id: str,
        artifact_id: str,
        *,
        generation: str,
        revision: str,
    ) -> dict[str, object]:
        self._expect_tokens(generation, revision)
        run = self._run(run_id)
        artifact = self._artifact(run_id, artifact_id, generation, revision)
        description = describe_browser_artifact(
            workflow_id=run.workflow_id,
            output_type=_text(artifact.metadata, "output_type"),
            name=artifact.name,
        )
        if description is None or not description.available:
            raise BrowserTrackError(
                404,
                "BROWSER_TRACK_NOT_FOUND",
                "This artifact cannot be prepared as a browser track.",
                path="artifact_id",
            )
        if description.source_form != "display_derivative":
            raise BrowserTrackError(
                409,
                "BROWSER_TRACK_PREPARE_CONFLICT",
                "This track is already an original indexed file.",
                path="artifact_id",
                reason="original_file",
            )
        size = _size(artifact)
        if size is not None and size > MAX_SOURCE_BYTES:
            raise BrowserTrackError(
                409,
                "BROWSER_TRACK_PREPARE_CONFLICT",
                "The source is too large to index in this view.",
                path="artifact_id",
                reason="source_too_large",
            )
        source = self._read_source(run_id, artifact_id, generation, revision)
        try:
            derivative = build_display_derivative(
                kind=description.kind,
                source=source,
                source_artifact_id=artifact.artifact_id,
                source_revision=artifact.revision,
                source_generation=generation,
                source_output_type=_text(artifact.metadata, "output_type"),
            )
        except DisplayDerivativeError as exc:
            raise BrowserTrackError(
                409,
                "BROWSER_TRACK_PREPARE_CONFLICT",
                _prepare_message(exc.reason),
                path="artifact_id",
                reason=exc.reason,
            ) from None
        try:
            public = write_derivative(
                self._root,
                run_id,
                artifact_id,
                revision,
                data=derivative.data,
                index=derivative.index,
                provenance=derivative.provenance,
            )
        except BrowserStoreError as exc:
            raise BrowserTrackError(
                409,
                "BROWSER_TRACK_PREPARE_CONFLICT",
                "The display file could not be stored.",
                path="artifact_id",
                reason=exc.reason,
            ) from None
        self._artifact(run_id, artifact_id, generation, revision)
        return self._ready_track(run_id, artifact, generation, description, public)

    def region_summary(
        self,
        run_id: str,
        artifact_id: str,
        *,
        generation: str,
        revision: str,
        chrom: str,
        start: int,
        end: int,
    ) -> dict[str, object]:
        self._expect_tokens(generation, revision)
        if start < 0 or end <= start or end - start > 10_000_000:
            raise BrowserTrackError(
                400,
                "BROWSER_REGION_INVALID",
                "The requested region is invalid.",
                path="region",
            )
        run = self._run(run_id)
        artifact = self._artifact(run_id, artifact_id, generation, revision)
        description = describe_browser_artifact(
            workflow_id=run.workflow_id,
            output_type=_text(artifact.metadata, "output_type"),
            name=artifact.name,
        )
        if description is None or not description.available:
            raise BrowserTrackError(
                404,
                "BROWSER_TRACK_NOT_FOUND",
                "This artifact is not a browser track.",
                path="artifact_id",
            )
        if description.source_form != "display_derivative":
            return {
                "ok": True,
                "kind": description.kind,
                "records": None,
                "truncated": False,
                "complete": False,
                "detail": (
                    "This track is not a record count. A missing file is not "
                    "reported as zero signal."
                    if description.kind == "signal"
                    else "This track does not publish a PET or interval count."
                ),
            }
        indexed = self._load_indexed(run_id, artifact_id, revision)
        try:
            rows, truncated = query_indexed_bed(
                indexed,
                chrom,
                start,
                end,
                limit=_REGION_COUNT_LIMIT,
            )
        except BrowserIndexError as exc:
            raise BrowserTrackError(
                409,
                "BROWSER_TRACK_READ_CONFLICT",
                "The display index could not be queried.",
                path="artifact_id",
                reason="derivative_changed",
            ) from exc
        return {
            "ok": True,
            "kind": description.kind,
            "records": len(rows),
            "truncated": truncated,
            "complete": not truncated,
            "detail": (
                "No records in this window. This is not a zero-signal measurement."
                if not rows and not truncated
                else "The window count stopped at the display limit."
                if truncated
                else "Record count for this window."
            ),
        }

    def open_bytes(
        self,
        run_id: str,
        artifact_id: str,
        *,
        generation: str,
        revision: str,
        index: bool,
        range_header: str | None,
        head: bool,
    ) -> BrowserByteBody:
        self._expect_tokens(generation, revision)
        run = self._run(run_id)
        artifact = self._artifact(run_id, artifact_id, generation, revision)
        description = describe_browser_artifact(
            workflow_id=run.workflow_id,
            output_type=_text(artifact.metadata, "output_type"),
            name=artifact.name,
        )
        if description is None or not description.available:
            raise BrowserTrackError(
                404,
                "BROWSER_TRACK_NOT_FOUND",
                "This artifact is not a browser track.",
                path="artifact_id",
            )
        if description.source_form == "display_derivative":
            return self._open_derivative_bytes(
                run_id,
                artifact_id,
                revision,
                index=index,
                range_header=range_header,
                head=head,
            )
        if index:
            raise BrowserTrackError(
                404,
                "BROWSER_TRACK_INDEX_NOT_FOUND",
                "This track has no separate index file.",
                path="revision",
                reason="index_not_applicable",
            )
        result = self._downloads.prepare(
            run_id,
            artifact_id,
            expected_generation=generation,
            expected_revision=revision,
        )
        if result.is_failure or result.value is None:
            issue = result.issues[0]
            reason = issue.context.get("reason_code") if issue.context else None
            raise BrowserTrackError(
                409
                if issue.code == "RUN_ARTIFACT_DOWNLOAD_CONFLICT"
                else 404
                if issue.code in {"RUN_NOT_FOUND", "RUN_ARTIFACT_NOT_FOUND"}
                else 500,
                "BROWSER_TRACK_READ_CONFLICT"
                if issue.code == "RUN_ARTIFACT_DOWNLOAD_CONFLICT"
                else issue.code,
                "Artifact content is no longer available as indexed."
                if issue.code == "RUN_ARTIFACT_DOWNLOAD_CONFLICT"
                else issue.message,
                path=issue.path or "artifact_id",
                reason=reason if isinstance(reason, str) else None,
            )
        plan = result.value
        try:
            selected = parse_byte_range(range_header, plan.size_bytes)
        except BrowserByteError:
            plan.close()
            return _range_error(plan.size_bytes)
        return _body_from_plan(plan, selected, head)

    def _open_derivative_bytes(
        self,
        run_id: str,
        artifact_id: str,
        revision: str,
        *,
        index: bool,
        range_header: str | None,
        head: bool,
    ) -> BrowserByteBody:
        try:
            stored = open_derivative(
                self._root,
                run_id,
                artifact_id,
                revision,
                index=index,
            )
        except BrowserStoreError as exc:
            raise BrowserTrackError(
                409,
                "BROWSER_TRACK_READ_CONFLICT",
                "The display file is no longer available as indexed.",
                path="revision",
                reason=exc.reason,
            ) from None
        try:
            selected = parse_byte_range(range_header, stored.size)
        except BrowserByteError:
            stored.close()
            return _range_error(stored.size)
        return _body_from_store(stored, selected, head)

    def _load_indexed(self, run_id: str, artifact_id: str, revision: str):
        from encode_pipeline.services.browser_bgzf import IndexedText

        data = self._read_stored(run_id, artifact_id, revision, index=False)
        index = self._read_stored(run_id, artifact_id, revision, index=True)
        return IndexedText(bgzf=data, tabix=index, record_count=0)

    def _read_stored(
        self,
        run_id: str,
        artifact_id: str,
        revision: str,
        *,
        index: bool,
    ) -> bytes:
        stored = open_derivative(
            self._root,
            run_id,
            artifact_id,
            revision,
            index=index,
        )
        try:
            return b"".join(iter_descriptor(stored.descriptor, 0, stored.size))
        finally:
            stored.close()

    def _read_source(
        self,
        run_id: str,
        artifact_id: str,
        generation: str,
        revision: str,
    ) -> bytes:
        result = self._downloads.prepare(
            run_id,
            artifact_id,
            expected_generation=generation,
            expected_revision=revision,
        )
        if result.is_failure or result.value is None:
            raise BrowserTrackError(
                409,
                "BROWSER_TRACK_PREPARE_CONFLICT",
                "Artifact content is no longer available as indexed.",
                path="revision",
                reason="source_changed",
            )
        return b"".join(iter_plan_bytes(result.value, None))

    def _track_payload(
        self,
        run_id: str,
        artifact: RunArtifactRef,
        generation: str,
        description: BrowserArtifactDescription,
        reference: dict[str, object],
    ) -> dict[str, object]:
        mismatch = _reference_mismatch(artifact, reference)
        if mismatch is not None:
            return {
                "artifact_id": artifact.artifact_id,
                "name": artifact.name,
                "output_type": _text(artifact.metadata, "output_type"),
                "revision": artifact.revision,
                "role": description.role,
                "kind": description.kind,
                "wash_u_type": description.wash_u_type,
                "label": description.label,
                "detail": description.detail,
                "score_meaning": description.score_meaning,
                "source_form": description.source_form,
                "preparation": "reference_mismatch",
                "content_path": None,
                "index_path": None,
                "provenance": None,
            }
        if description.source_form == "original":
            return self._ready_track(run_id, artifact, generation, description, None)
        size = _size(artifact)
        if size is not None and size > MAX_SOURCE_BYTES:
            payload = self._ready_track(run_id, artifact, generation, description, None)
            payload["preparation"] = "source_too_large"
            payload["content_path"] = None
            payload["index_path"] = None
            return payload
        existing = read_derivative_summary(
            self._root,
            run_id,
            artifact.artifact_id,
            artifact.revision,
        )
        if (
            existing is not None
            and existing.get("source_revision") == artifact.revision
        ):
            return self._ready_track(
                run_id,
                artifact,
                generation,
                description,
                existing,
            )
        payload = self._ready_track(run_id, artifact, generation, description, None)
        payload["preparation"] = "needed"
        payload["content_path"] = None
        payload["index_path"] = None
        payload["provenance"] = None
        return payload

    def _ready_track(
        self,
        run_id: str,
        artifact: RunArtifactRef,
        generation: str,
        description: BrowserArtifactDescription,
        provenance: dict[str, object] | None,
    ) -> dict[str, object]:
        content = (
            f"/api/v1/runs/{run_id}/browser-tracks/{artifact.artifact_id}"
            f"/bytes/{generation}/{artifact.revision}"
        )
        indexed = description.source_form == "display_derivative"
        ready = description.source_form == "original" or provenance is not None
        return {
            "artifact_id": artifact.artifact_id,
            "name": artifact.name,
            "output_type": _text(artifact.metadata, "output_type"),
            "revision": artifact.revision,
            "role": description.role,
            "kind": description.kind,
            "wash_u_type": description.wash_u_type,
            "label": f"{description.label}: {artifact.name}",
            "detail": description.detail,
            "score_meaning": description.score_meaning,
            "source_form": description.source_form,
            "preparation": "ready" if ready else "needed",
            "content_path": content if ready else None,
            "index_path": content + ".tbi"
            if indexed and provenance is not None
            else None,
            "provenance": provenance,
        }

    def _describe(
        self,
        workflow_id: str,
        artifact: RunArtifactRef,
        reference: dict[str, object],
    ) -> dict[str, object] | None:
        del reference
        description = describe_browser_artifact(
            workflow_id=workflow_id,
            output_type=_text(artifact.metadata, "output_type"),
            name=artifact.name,
        )
        if description is None:
            return None
        if description.available:
            return {"available": True, "description": description}
        return {
            "available": False,
            "unavailable": {
                "artifact_id": artifact.artifact_id,
                "name": artifact.name,
                "output_type": _text(artifact.metadata, "output_type"),
                "reason": description.unavailable_reason,
                "detail": description.detail,
            },
        }

    def _reference(self, run_id: str, workflow_id: str) -> dict[str, object]:
        del workflow_id
        try:
            binding = self._runs.get_run_reference_binding(run_id)
        except KeyError as exc:
            raise BrowserTrackError(
                404,
                "RUN_NOT_FOUND",
                "Run was not found.",
                path="run_id",
            ) from exc
        except Exception as exc:
            raise BrowserTrackError(
                500,
                "BROWSER_TRACK_UNAVAILABLE",
                "The browser view is temporarily unavailable.",
                path="run_id",
            ) from exc
        if binding is None:
            return {"status": "unbound", "viewer_scaffold": None}
        summary = self._summary(binding.revision_id)
        if (
            summary is None
            or summary.profile_id != binding.profile_id
            or summary.public_identity_sha256 != binding.revision_public_identity_sha256
        ):
            return {
                "status": "unavailable",
                "profile_id": binding.profile_id,
                "revision_id": binding.revision_id,
                "identity_sha256": binding.revision_public_identity_sha256,
                "binding_digest": binding.binding_digest,
                "assembly_label": None,
                "viewer_scaffold": None,
            }
        scaffold = viewer_scaffold(summary.assembly)
        return {
            "status": "bound",
            "profile_id": binding.profile_id,
            "revision_id": binding.revision_id,
            "identity_sha256": binding.revision_public_identity_sha256,
            "binding_digest": binding.binding_digest,
            "assembly_label": summary.assembly,
            "viewer_scaffold": scaffold,
        }

    def _summary(self, revision_id: str):
        try:
            return self._references.get_revision_summary(revision_id)
        except Exception:
            return None

    def _artifacts(self, run_id: str) -> tuple[RunArtifactRef, ...]:
        try:
            return self._runs.list_artifacts(run_id, limit=_ARTIFACT_SCAN_LIMIT)
        except KeyError as exc:
            raise BrowserTrackError(
                404,
                "RUN_NOT_FOUND",
                "Run was not found.",
                path="run_id",
            ) from exc
        except Exception as exc:
            raise BrowserTrackError(
                500,
                "BROWSER_TRACK_UNAVAILABLE",
                "The browser view is temporarily unavailable.",
                path="run_id",
            ) from exc

    def _artifact(
        self,
        run_id: str,
        artifact_id: str,
        generation: str,
        revision: str,
    ) -> RunArtifactRef:
        try:
            _found, artifact = self._runs.get_artifact_at_generation(
                run_id,
                artifact_id,
                expected_generation=generation,
            )
        except ResultGenerationChangedError as exc:
            raise BrowserTrackError(
                409,
                "BROWSER_TRACK_READ_CONFLICT",
                "Artifact generation changed. Reload the browser view.",
                path="generation",
                reason="generation_changed",
            ) from exc
        except KeyError as exc:
            self._run(run_id)
            raise BrowserTrackError(
                404,
                "RUN_ARTIFACT_NOT_FOUND",
                "Artifact was not found for this run.",
                path="artifact_id",
            ) from exc
        except Exception as exc:
            raise BrowserTrackError(
                500,
                "BROWSER_TRACK_UNAVAILABLE",
                "The browser view is temporarily unavailable.",
                path="artifact_id",
            ) from exc
        if artifact.revision != revision:
            raise BrowserTrackError(
                409,
                "BROWSER_TRACK_READ_CONFLICT",
                "Artifact content is no longer available as indexed.",
                path="revision",
                reason="revision_changed",
            )
        return artifact

    def _generation(self, run_id: str) -> str | None:
        try:
            state = self._runs.get_result_state(run_id)
        except KeyError as exc:
            raise BrowserTrackError(
                404,
                "RUN_NOT_FOUND",
                "Run was not found.",
                path="run_id",
            ) from exc
        except Exception as exc:
            raise BrowserTrackError(
                500,
                "BROWSER_TRACK_UNAVAILABLE",
                "The browser view is temporarily unavailable.",
                path="run_id",
            ) from exc
        generation = state.artifact_generation
        if generation is None:
            return None
        try:
            return validate_artifact_generation(generation)
        except (TypeError, ValueError) as exc:
            raise BrowserTrackError(
                500,
                "BROWSER_TRACK_UNAVAILABLE",
                "The browser view is temporarily unavailable.",
                path="generation",
            ) from exc

    def _run(self, run_id: str):
        try:
            return self._runs.get_run(run_id)
        except KeyError as exc:
            raise BrowserTrackError(
                404,
                "RUN_NOT_FOUND",
                "Run was not found.",
                path="run_id",
            ) from exc
        except Exception as exc:
            raise BrowserTrackError(
                500,
                "BROWSER_TRACK_UNAVAILABLE",
                "The browser view is temporarily unavailable.",
                path="run_id",
            ) from exc

    @staticmethod
    def _expect_tokens(generation: str, revision: str) -> None:
        try:
            validate_artifact_generation(generation)
            validate_artifact_revision(revision)
        except (TypeError, ValueError) as exc:
            raise BrowserTrackError(
                409,
                "BROWSER_TRACK_READ_CONFLICT",
                "Artifact content is no longer available as indexed.",
                path="generation",
                reason="expectation_invalid",
            ) from exc


def _body_from_plan(
    plan: object, selected: ByteRange | None, head: bool
) -> BrowserByteBody:
    size = int(getattr(plan, "size_bytes"))
    if head:
        plan.close()
        status, headers = _range_headers(size, selected)
        return BrowserByteBody(
            status, "application/octet-stream", headers, None, plan.close
        )

    def chunks() -> Iterator[bytes]:
        yield from iter_plan_bytes(plan, selected)

    status, headers = _range_headers(size, selected)
    return BrowserByteBody(
        status,
        "application/octet-stream",
        headers,
        chunks(),
        plan.close,
    )


def _body_from_store(
    stored: StoredDerivative,
    selected: ByteRange | None,
    head: bool,
) -> BrowserByteBody:
    if selected is None:
        start, length = 0, stored.size
    else:
        start, length = selected.start, selected.length
    status, headers = _range_headers(stored.size, selected)
    if head:
        stored.close()
        return BrowserByteBody(
            status,
            "application/octet-stream",
            headers,
            None,
            stored.close,
        )

    def chunks() -> Iterator[bytes]:
        try:
            if length:
                yield from iter_descriptor(stored.descriptor, start, length)
        finally:
            stored.close()

    return BrowserByteBody(
        status,
        "application/octet-stream",
        headers,
        chunks(),
        stored.close,
    )


def _range_headers(size: int, selected: ByteRange | None) -> tuple[int, dict[str, str]]:
    if selected is None:
        length = size
        status = 200
        headers = {
            "Accept-Ranges": "bytes",
            "Content-Length": str(length),
            "Cache-Control": "private, no-store",
            "X-Content-Type-Options": "nosniff",
        }
        return status, headers
    headers = {
        "Accept-Ranges": "bytes",
        "Content-Length": str(selected.length),
        "Content-Range": f"bytes {selected.start}-{selected.end}/{size}"
        if selected.length
        else f"bytes */{size}",
        "Cache-Control": "private, no-store",
        "X-Content-Type-Options": "nosniff",
    }
    if selected.length == 0 and selected.end < selected.start:
        headers["Content-Range"] = f"bytes */{size}"
        return 200, headers
    return 206, headers


def _range_error(size: int) -> BrowserByteBody:
    return BrowserByteBody(
        416,
        "application/octet-stream",
        {
            "Accept-Ranges": "bytes",
            "Content-Range": f"bytes */{size}",
            "Content-Length": "0",
            "Cache-Control": "private, no-store",
            "X-Content-Type-Options": "nosniff",
        },
        None,
        lambda: None,
    )


def _reference_mismatch(
    artifact: RunArtifactRef,
    reference: dict[str, object],
) -> str | None:
    if reference.get("status") != "bound":
        return None
    assembly = _text(artifact.metadata, "reference_assembly")
    label = reference.get("assembly_label")
    if assembly and isinstance(label, str) and assembly != label:
        return "reference_mismatch"
    binding_id = _text(artifact.metadata, "reference_binding_id")
    if not binding_id:
        return None
    allowed = {
        reference.get("profile_id"),
        reference.get("revision_id"),
        reference.get("binding_digest"),
    }
    if binding_id not in allowed:
        return "reference_mismatch"
    return None


def _text(metadata: object, key: str) -> str:
    if not isinstance(metadata, Mapping):
        return ""
    value = metadata.get(key)
    return value if isinstance(value, str) else ""


def _size(artifact: RunArtifactRef) -> int | None:
    if not isinstance(artifact.metadata, Mapping):
        return None
    value = artifact.metadata.get("size_bytes")
    if isinstance(value, bool) or not isinstance(value, int):
        return None
    return value


def _prepare_message(reason: str) -> str:
    if reason == "source_too_large":
        return "The source is too large to index in this view."
    if reason == "source_has_no_records":
        return "The source has no records to display."
    return "The source could not be indexed for display."
