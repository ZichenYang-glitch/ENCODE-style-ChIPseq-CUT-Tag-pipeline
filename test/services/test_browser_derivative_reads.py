"""Derivative identity, byte ranges, and browser session states."""

from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import sys

import pytest

from encode_pipeline.api.main import create_app
from encode_pipeline.services.browser_bytes import (
    BrowserByteError,
    iter_plan_bytes,
    parse_byte_range,
)
from encode_pipeline.services.browser_derivatives import MAX_SOURCE_BYTES
from encode_pipeline.services.browser_store import (
    BrowserStoreError,
    derivative_directory,
    iter_descriptor,
    open_derivative,
    read_derivative_summary,
    write_derivative,
)
from encode_pipeline.services.browser_tracks import (
    BrowserTrackError,
    BrowserTrackService,
)


_API_DIR = Path(__file__).resolve().parents[1] / "api"
if str(_API_DIR) not in sys.path:
    sys.path.insert(0, str(_API_DIR))
_HELPER_SPEC = importlib.util.spec_from_file_location(
    "browser_track_route_helpers",
    _API_DIR / "test_routes_browser_tracks.py",
)
assert _HELPER_SPEC is not None and _HELPER_SPEC.loader is not None
_helpers = importlib.util.module_from_spec(_HELPER_SPEC)
_HELPER_SPEC.loader.exec_module(_helpers)


def _fds() -> set[int]:
    return {int(name) for name in os.listdir("/proc/self/fd") if name.isdigit()}


def _derivative_descriptors() -> list[str]:
    found: list[str] = []
    for name in os.listdir("/proc/self/fd"):
        try:
            linked = os.readlink(f"/proc/self/fd/{name}")
        except OSError:
            continue
        if "browser-derivatives" in linked:
            found.append(linked)
    return found


def _service(app) -> BrowserTrackService:
    return BrowserTrackService(
        run_service=app.state.run_service,
        download_service=app.state.artifact_download_service,
        workspace_root=app.state.workspace_root,
        reference_profile_service=app.state.reference_profile_service,
    )


def _app(tmp_path: Path):
    app = create_app(
        database_url=f"sqlite:///{tmp_path / 'platform.db'}",
        workspace_root=tmp_path / "workspaces",
    )
    _helpers.seed_test_authentication(app)
    return app


def _client(app):
    return _helpers.ApiTestClient(app)


def _peak(app):
    client = _client(app)
    run_id = _helpers._create_run(client)
    _helpers._bind_reference(client, run_id, assembly="hg38")
    artifact = _helpers._record(
        client,
        run_id,
        "artifact-peak",
        b"chr1\t1500\t1700\tpeak-a\t80\t+\n",
        name="sample_peaks.narrowPeak",
        output_type="macs3_peak",
        revision_kind="content",
        reference_assembly="hg38",
    )
    return run_id, artifact, _helpers._generation(client, run_id)


def _store(
    tmp_path: Path, *, data: bytes = b"display-bytes", index: bytes = b"index-bytes"
):
    provenance = {
        "kind": "peak",
        "source_revision": "rev-1",
        "data_size": len(data),
        "index_size": len(index),
    }
    public = write_derivative(
        tmp_path,
        "run-1",
        "artifact-peak",
        "rev-1",
        data=data,
        index=index,
        provenance=provenance,
    )
    return public


def test_stored_derivative_reads_back_and_releases_its_descriptors(
    tmp_path: Path,
) -> None:
    data = b"display-bytes"
    index = b"index-bytes"
    before = _fds()
    public = _store(tmp_path, data=data, index=index)
    assert "device" not in public
    assert "index_stamp" not in public
    opened = open_derivative(tmp_path, "run-1", "artifact-peak", "rev-1", index=False)
    assert b"".join(iter_descriptor(opened.descriptor, 0, opened.size)) == data
    opened.close()
    opened.close()
    indexed = open_derivative(tmp_path, "run-1", "artifact-peak", "rev-1", index=True)
    assert b"".join(iter_descriptor(indexed.descriptor, 0, indexed.size)) == index
    indexed.close()
    assert _fds() == before


def test_changed_missing_linked_or_invalid_derivatives_are_refused(
    tmp_path: Path,
) -> None:
    before = _fds()
    _store(tmp_path)
    directory = derivative_directory(tmp_path, "run-1", "artifact-peak", "rev-1")
    data_path = directory / "track.gz"
    replacement = directory / "replacement"
    replacement.write_bytes(b"replaced-display")
    os.replace(replacement, data_path)
    with pytest.raises(BrowserStoreError, match="derivative_changed"):
        open_derivative(tmp_path, "run-1", "artifact-peak", "rev-1", index=False)

    _store(tmp_path)
    directory = derivative_directory(tmp_path, "run-1", "artifact-peak", "rev-1")
    summary = json.loads((directory / "provenance.json").read_text(encoding="ascii"))
    summary["data_size"] = int(summary["data_size"]) + 1
    (directory / "provenance.json").write_text(
        json.dumps(summary),
        encoding="ascii",
    )
    with pytest.raises(BrowserStoreError, match="derivative_changed"):
        open_derivative(tmp_path, "run-1", "artifact-peak", "rev-1", index=False)

    _store(tmp_path)
    directory = derivative_directory(tmp_path, "run-1", "artifact-peak", "rev-1")
    (directory / "track.gz").unlink()
    with pytest.raises(BrowserStoreError, match="derivative_missing"):
        open_derivative(tmp_path, "run-1", "artifact-peak", "rev-1", index=False)

    _store(tmp_path)
    directory = derivative_directory(tmp_path, "run-1", "artifact-peak", "rev-1")
    target = directory / "elsewhere"
    target.write_bytes(b"other")
    data_path = directory / "track.gz"
    data_path.unlink()
    data_path.symlink_to(target.name)
    with pytest.raises(BrowserStoreError, match="derivative_missing"):
        open_derivative(tmp_path, "run-1", "artifact-peak", "rev-1", index=False)

    _store(tmp_path)
    directory = derivative_directory(tmp_path, "run-1", "artifact-peak", "rev-1")
    (directory / "provenance.json").write_bytes(b"\xff")
    with pytest.raises(BrowserStoreError, match="derivative_invalid"):
        open_derivative(tmp_path, "run-1", "artifact-peak", "rev-1", index=False)

    _store(tmp_path)
    directory = derivative_directory(tmp_path, "run-1", "artifact-peak", "rev-1")
    summary = json.loads((directory / "provenance.json").read_text(encoding="ascii"))
    summary.pop("device")
    (directory / "provenance.json").write_text(
        json.dumps(summary),
        encoding="ascii",
    )
    with pytest.raises(BrowserStoreError, match="derivative_changed"):
        open_derivative(tmp_path, "run-1", "artifact-peak", "rev-1", index=False)

    summary["index_stamp"] = "not-a-stamp"
    (directory / "provenance.json").write_text(json.dumps(summary), encoding="ascii")
    with pytest.raises(BrowserStoreError, match="derivative_invalid"):
        open_derivative(tmp_path, "run-1", "artifact-peak", "rev-1", index=True)

    link = tmp_path / "linked-root"
    link.symlink_to(tmp_path, target_is_directory=True)
    with pytest.raises(BrowserStoreError, match="derivative_invalid"):
        write_derivative(
            link,
            "run-2",
            "artifact-peak",
            "rev-2",
            data=b"a",
            index=b"b",
            provenance={"kind": "peak"},
        )
    assert _fds() == before


def test_summary_reader_ignores_missing_linked_and_broken_provenance(
    tmp_path: Path,
) -> None:
    assert read_derivative_summary(tmp_path, "missing", "artifact", "rev") is None
    _store(tmp_path)
    directory = derivative_directory(tmp_path, "run-1", "artifact-peak", "rev-1")
    summary = read_derivative_summary(tmp_path, "run-1", "artifact-peak", "rev-1")
    assert summary is not None
    assert "inode" not in summary
    provenance = directory / "provenance.json"
    provenance.unlink()
    provenance.symlink_to(directory / "track.gz")
    assert read_derivative_summary(tmp_path, "run-1", "artifact-peak", "rev-1") is None
    provenance.unlink()
    provenance.write_text("[]", encoding="ascii")
    assert read_derivative_summary(tmp_path, "run-1", "artifact-peak", "rev-1") is None


def test_truncated_open_derivative_is_not_a_short_success(tmp_path: Path) -> None:
    before = _fds()
    payload = b"display-bytes"
    _store(tmp_path, data=payload)
    opened = open_derivative(tmp_path, "run-1", "artifact-peak", "rev-1", index=False)
    os.truncate(
        derivative_directory(tmp_path, "run-1", "artifact-peak", "rev-1") / "track.gz",
        0,
    )
    with pytest.raises(BrowserStoreError, match="derivative_changed"):
        b"".join(iter_descriptor(opened.descriptor, 0, opened.size))
    opened.close()
    assert _fds() == before


def test_derivative_ranges_head_and_early_close_release_descriptors(
    tmp_path: Path,
) -> None:
    app = _app(tmp_path)
    run_id, artifact, generation = _peak(app)
    service = _service(app)
    prepared = service.prepare(
        run_id,
        artifact.artifact_id,
        generation=generation,
        revision=artifact.revision,
    )
    assert prepared["preparation"] == "ready"
    before = _fds()
    full = service.open_bytes(
        run_id,
        artifact.artifact_id,
        generation=generation,
        revision=artifact.revision,
        index=False,
        range_header=None,
        head=False,
    )
    assert full.status == 200
    body = b"".join(full.chunks or ())
    assert body
    assert full.headers["Content-Length"] == str(len(body))

    ranged = service.open_bytes(
        run_id,
        artifact.artifact_id,
        generation=generation,
        revision=artifact.revision,
        index=False,
        range_header="bytes=0-0",
        head=False,
    )
    assert ranged.status == 206
    assert b"".join(ranged.chunks or ()) == body[:1]
    assert ranged.headers["Content-Range"].startswith("bytes 0-0/")

    suffix = service.open_bytes(
        run_id,
        artifact.artifact_id,
        generation=generation,
        revision=artifact.revision,
        index=False,
        range_header="bytes=-1",
        head=False,
    )
    assert suffix.status == 206
    assert b"".join(suffix.chunks or ()) == body[-1:]

    head = service.open_bytes(
        run_id,
        artifact.artifact_id,
        generation=generation,
        revision=artifact.revision,
        index=True,
        range_header="bytes=0-0",
        head=True,
    )
    assert head.status == 206
    assert head.chunks is None
    assert head.headers["Content-Length"] == "1"

    rejected = service.open_bytes(
        run_id,
        artifact.artifact_id,
        generation=generation,
        revision=artifact.revision,
        index=False,
        range_header="bytes=0-1,2-3",
        head=False,
    )
    assert rejected.status == 416
    beyond = service.open_bytes(
        run_id,
        artifact.artifact_id,
        generation=generation,
        revision=artifact.revision,
        index=False,
        range_header=f"bytes={len(body) + 5}-{len(body) + 6}",
        head=False,
    )
    assert beyond.status == 416

    streaming = service.open_bytes(
        run_id,
        artifact.artifact_id,
        generation=generation,
        revision=artifact.revision,
        index=False,
        range_header=None,
        head=False,
    )
    assert streaming.chunks is not None
    next(streaming.chunks)
    streaming.chunks.close()
    assert _derivative_descriptors() == []
    assert _fds() == before

    truncated = service.open_bytes(
        run_id,
        artifact.artifact_id,
        generation=generation,
        revision=artifact.revision,
        index=False,
        range_header=None,
        head=False,
    )
    stored = derivative_directory(
        app.state.workspace_root,
        run_id,
        artifact.artifact_id,
        artifact.revision,
    )
    os.truncate(stored / "track.gz", 0)
    with pytest.raises(BrowserStoreError, match="derivative_changed"):
        b"".join(truncated.chunks or ())
    assert _derivative_descriptors() == []
    assert _fds() == before


def test_session_keeps_reference_and_preparation_states_distinct(
    tmp_path: Path,
) -> None:
    app = _app(tmp_path)
    client = _client(app)
    run_id = _helpers._create_run(client)
    service = _service(app)
    empty = service.session(run_id)
    assert empty["tracks"] == []
    assert empty["unavailable"] == []
    assert empty["reference"]["status"] == "unbound"
    assert empty["reference"]["viewer_scaffold"] is None

    _helpers._bind_reference(client, run_id, assembly="hg38")
    peak = _helpers._record(
        client,
        run_id,
        "artifact-peak",
        b"chr1\t1500\t1700\tpeak-a\t80\t+\n",
        name="sample_peaks.narrowPeak",
        output_type="macs3_peak",
        revision_kind="content",
        reference_assembly="mm10",
    )
    _helpers._record(
        client,
        run_id,
        "artifact-cool",
        b"not-a-matrix",
        name="contacts.mcool",
        output_type="mcool",
    )
    mismatched = service.session(run_id)
    peak_row = next(
        item for item in mismatched["tracks"] if item["artifact_id"] == "artifact-peak"
    )
    assert peak_row["preparation"] == "reference_mismatch"
    assert peak_row["content_path"] is None
    assert peak_row["index_path"] is None
    assert any(
        item["artifact_id"] == "artifact-cool" for item in mismatched["unavailable"]
    )

    recorded_needed = _helpers._record(
        client,
        run_id,
        "artifact-needed",
        b"chr1\t10\t20\tpeak-b\t5\t+\n",
        name="needed.narrowPeak",
        output_type="macs3_peak",
        revision_kind="content",
        reference_assembly="hg38",
    )
    prepared_source = _helpers._record(
        client,
        run_id,
        "artifact-ready",
        b"chr1\t30\t40\tpeak-c\t9\t+\n",
        name="ready.narrowPeak",
        output_type="macs3_peak",
        revision_kind="content",
        reference_assembly="hg38",
    )
    generation = _helpers._generation(client, run_id)
    ready = service.prepare(
        run_id,
        prepared_source.artifact_id,
        generation=generation,
        revision=prepared_source.revision,
    )
    assert ready["content_path"]
    assert ready["index_path"] == ready["content_path"] + ".tbi"
    viewed = service.session(run_id)
    by_id = {item["artifact_id"]: item for item in viewed["tracks"]}
    assert recorded_needed.artifact_id == "artifact-needed"
    assert by_id["artifact-needed"]["preparation"] == "needed"
    assert by_id["artifact-needed"]["content_path"] is None
    assert by_id["artifact-ready"]["preparation"] == "ready"
    assert by_id["artifact-ready"]["content_path"]
    assert by_id["artifact-peak"]["content_path"] is None
    assert peak.artifact_id == "artifact-peak"


def test_unavailable_reference_has_no_scaffold(tmp_path: Path) -> None:
    app = _app(tmp_path)
    client = _client(app)
    run_id = _helpers._create_run(client)
    _helpers._bind_reference(client, run_id, assembly="hg38")
    real = app.state.reference_profile_service

    class _UnreadableReferences:
        def __getattr__(self, name: str):
            return getattr(real, name)

        def get_revision_summary(self, _revision_id: str):
            raise RuntimeError("sqlite:///secret-reference.db")

    service = BrowserTrackService(
        run_service=app.state.run_service,
        download_service=app.state.artifact_download_service,
        workspace_root=app.state.workspace_root,
        reference_profile_service=_UnreadableReferences(),
    )
    viewed = service.session(run_id)
    assert viewed["reference"]["status"] == "unavailable"
    assert viewed["reference"]["viewer_scaffold"] is None
    assert viewed["reference"]["assembly_label"] is None
    assert "secret" not in json.dumps(viewed)


def test_failed_prepare_does_not_publish_a_ready_derivative(tmp_path: Path) -> None:
    app = _app(tmp_path)
    client = _client(app)
    run_id = _helpers._create_run(client)
    _helpers._bind_reference(client, run_id, assembly="hg38")
    empty = _helpers._record(
        client,
        run_id,
        "artifact-empty",
        b"",
        name="empty.narrowPeak",
        output_type="macs3_peak",
        revision_kind="content",
        reference_assembly="hg38",
    )
    generation = _helpers._generation(client, run_id)
    service = _service(app)
    with pytest.raises(BrowserTrackError) as empty_error:
        service.prepare(
            run_id,
            empty.artifact_id,
            generation=generation,
            revision=empty.revision,
        )
    assert empty_error.value.status == 409
    assert empty_error.value.reason == "source_has_no_records"
    assert "no records" in empty_error.value.message
    assert (
        read_derivative_summary(
            app.state.workspace_root,
            run_id,
            empty.artifact_id,
            empty.revision,
        )
        is None
    )
    huge = _helpers._record(
        client,
        run_id,
        "artifact-huge",
        b"chr1\t1\t2\tpeak\t1\t+\n",
        name="huge.narrowPeak",
        output_type="macs3_peak",
        revision_kind="content",
        reference_assembly="hg38",
    )
    generation = _helpers._generation(client, run_id)
    stored = app.state.run_service.get_artifact_at_generation(
        run_id,
        huge.artifact_id,
        expected_generation=generation,
    )[1]
    metadata = dict(stored.metadata)
    metadata["size_bytes"] = MAX_SOURCE_BYTES + 1
    from dataclasses import replace

    app.state.run_service.replace_artifacts(
        run_id,
        tuple(
            replace(item, metadata=metadata)
            if item.artifact_id == huge.artifact_id
            else item
            for item in app.state.run_service.list_artifacts(run_id)
        ),
    )
    current = _helpers._generation(client, run_id)
    with pytest.raises(BrowserTrackError) as huge_error:
        service.prepare(
            run_id,
            huge.artifact_id,
            generation=current,
            revision=huge.revision,
        )
    assert huge_error.value.reason == "source_too_large"
    assert "too large" in huge_error.value.message
    assert (
        read_derivative_summary(
            app.state.workspace_root,
            run_id,
            huge.artifact_id,
            huge.revision,
        )
        is None
    )


def test_repository_faults_stay_public_and_do_not_echo_the_cause(
    tmp_path: Path,
) -> None:
    app = _app(tmp_path)
    client = _client(app)
    run_id = _helpers._create_run(client)
    real = app.state.run_service

    class _Runs:
        def __getattr__(self, name: str):
            return getattr(real, name)

        def get_run(self, requested: str):
            if requested == "missing":
                raise KeyError(requested)
            raise RuntimeError("sqlite:///secret.db")

    service = BrowserTrackService(
        run_service=_Runs(),
        download_service=app.state.artifact_download_service,
        workspace_root=app.state.workspace_root,
        reference_profile_service=app.state.reference_profile_service,
    )
    with pytest.raises(BrowserTrackError) as not_found:
        service.session("missing")
    assert not_found.value.status == 404
    assert "secret" not in not_found.value.message
    with pytest.raises(BrowserTrackError) as unavailable:
        service.session(run_id)
    assert unavailable.value.status == 500
    assert unavailable.value.code == "BROWSER_TRACK_UNAVAILABLE"
    assert "secret" not in unavailable.value.message
    assert "sqlite" not in unavailable.value.message


def test_closed_prepare_and_reference_identity_do_not_publish_reads(
    tmp_path: Path,
) -> None:
    app = _app(tmp_path)
    client = _client(app)
    run_id = _helpers._create_run(client)
    peak = _helpers._record(
        client,
        run_id,
        "artifact-peak",
        b"chr1\t1500\t1700\tpeak-a\t80\t+\n",
        name="sample_peaks.narrowPeak",
        output_type="macs3_peak",
        revision_kind="content",
        reference_assembly="hg38",
    )
    service = _service(app)
    unbound = service.session(run_id)
    unbound_peak = next(
        item for item in unbound["tracks"] if item["artifact_id"] == peak.artifact_id
    )
    assert unbound["reference"]["status"] == "unbound"
    assert unbound_peak["preparation"] == "needed"
    assert unbound_peak["content_path"] is None
    assert unbound_peak["index_path"] is None

    _helpers._bind_reference(client, run_id, assembly="hg38")
    from dataclasses import replace

    foreign = replace(
        peak,
        metadata={**peak.metadata, "reference_binding_id": "not-a-binding"},
    )
    app.state.run_service.replace_artifacts(run_id, (foreign,))
    mismatched = service.session(run_id)
    mismatched_peak = next(
        item for item in mismatched["tracks"] if item["artifact_id"] == peak.artifact_id
    )
    assert mismatched["reference"]["status"] == "bound"
    assert mismatched_peak["preparation"] == "reference_mismatch"
    assert mismatched_peak["content_path"] is None
    assert mismatched_peak["index_path"] is None

    broken = _helpers._record(
        client,
        run_id,
        "artifact-broken",
        b"chr1\tonly-two\n",
        name="broken.narrowPeak",
        output_type="macs3_peak",
        revision_kind="content",
        reference_assembly="hg38",
    )
    generation = _helpers._generation(client, run_id)
    with pytest.raises(BrowserTrackError) as broken_error:
        service.prepare(
            run_id,
            broken.artifact_id,
            generation=generation,
            revision=broken.revision,
        )
    assert broken_error.value.status == 409
    assert broken_error.value.reason == "source_unreadable"
    assert "could not be indexed" in broken_error.value.message
    assert "only-two" not in broken_error.value.message
    assert (
        read_derivative_summary(
            app.state.workspace_root,
            run_id,
            broken.artifact_id,
            broken.revision,
        )
        is None
    )

    ignored_size = _helpers._record(
        client,
        run_id,
        "artifact-size",
        b"chr1\t10\t20\tpeak-b\t5\t+\n",
        name="size.narrowPeak",
        output_type="macs3_peak",
        revision_kind="content",
        reference_assembly="hg38",
    )
    generation = _helpers._generation(client, run_id)
    sized = app.state.run_service.get_artifact_at_generation(
        run_id,
        ignored_size.artifact_id,
        expected_generation=generation,
    )[1]
    app.state.run_service.replace_artifacts(
        run_id,
        tuple(
            replace(item, metadata={**dict(item.metadata), "size_bytes": "12"})
            if item.artifact_id == ignored_size.artifact_id
            else item
            for item in app.state.run_service.list_artifacts(run_id)
        ),
    )
    generation = _helpers._generation(client, run_id)
    with pytest.raises(BrowserTrackError) as sized_error:
        service.prepare(
            run_id,
            sized.artifact_id,
            generation=generation,
            revision=sized.revision,
        )
    assert sized_error.value.status == 409
    assert sized_error.value.reason == "source_changed"
    assert (
        read_derivative_summary(
            app.state.workspace_root,
            run_id,
            sized.artifact_id,
            sized.revision,
        )
        is None
    )

    signal = _helpers._record(
        client,
        run_id,
        "artifact-signal",
        (
            Path(__file__).resolve().parents[1]
            / "fixtures"
            / "browser"
            / "chr1-signal.bw"
        ).read_bytes(),
        name="chr1-signal.bw",
        output_type="cpm_bigwig",
    )
    generation = _helpers._generation(client, run_id)
    with pytest.raises(BrowserTrackError) as original:
        service.prepare(
            run_id,
            signal.artifact_id,
            generation=generation,
            revision=signal.revision,
        )
    assert original.value.reason == "original_file"
    assert (
        read_derivative_summary(
            app.state.workspace_root,
            run_id,
            signal.artifact_id,
            signal.revision,
        )
        is None
    )
    with pytest.raises(BrowserTrackError) as indexed:
        service.open_bytes(
            run_id,
            signal.artifact_id,
            generation=generation,
            revision=signal.revision,
            index=True,
            range_header=None,
            head=True,
        )
    assert indexed.value.status == 404
    assert indexed.value.reason == "index_not_applicable"


def test_store_refusal_and_repository_reads_stay_redacted(tmp_path: Path) -> None:
    app = _app(tmp_path)
    client = _client(app)
    run_id = _helpers._create_run(client)
    _helpers._bind_reference(client, run_id, assembly="hg38")
    peak = _helpers._record(
        client,
        run_id,
        "artifact-peak",
        b"chr1\t30\t40\tpeak-c\t9\t+\n",
        name="ready.narrowPeak",
        output_type="macs3_peak",
        revision_kind="content",
        reference_assembly="hg38",
    )
    cool = _helpers._record(
        client,
        run_id,
        "artifact-cool",
        b"not-a-matrix",
        name="contacts.mcool",
        output_type="mcool",
    )
    generation = _helpers._generation(client, run_id)
    linked = tmp_path / "linked-workspace"
    linked.symlink_to(app.state.workspace_root, target_is_directory=True)
    service = BrowserTrackService(
        run_service=app.state.run_service,
        download_service=app.state.artifact_download_service,
        workspace_root=linked,
        reference_profile_service=app.state.reference_profile_service,
    )
    with pytest.raises(BrowserTrackError) as stored:
        service.prepare(
            run_id,
            peak.artifact_id,
            generation=generation,
            revision=peak.revision,
        )
    assert stored.value.status == 409
    assert stored.value.reason == "derivative_invalid"
    assert stored.value.message == "The display file could not be stored."
    assert "linked-workspace" not in stored.value.message
    assert (
        read_derivative_summary(
            app.state.workspace_root,
            run_id,
            peak.artifact_id,
            peak.revision,
        )
        is None
    )

    direct = _service(app)
    with pytest.raises(BrowserTrackError) as missing_artifact:
        direct.open_bytes(
            run_id,
            "missing-artifact",
            generation=generation,
            revision=peak.revision,
            index=False,
            range_header=None,
            head=False,
        )
    assert missing_artifact.value.status == 404
    assert missing_artifact.value.code == "RUN_ARTIFACT_NOT_FOUND"
    with pytest.raises(BrowserTrackError) as bad_token:
        direct.open_bytes(
            run_id,
            peak.artifact_id,
            generation="not-a-generation",
            revision=peak.revision,
            index=False,
            range_header=None,
            head=False,
        )
    assert bad_token.value.status == 409
    assert bad_token.value.reason == "expectation_invalid"
    with pytest.raises(BrowserTrackError) as region:
        direct.region_summary(
            run_id,
            peak.artifact_id,
            generation=generation,
            revision=peak.revision,
            chrom="chr1",
            start=20,
            end=20,
        )
    assert region.value.status == 400
    assert region.value.code == "BROWSER_REGION_INVALID"
    with pytest.raises(BrowserTrackError) as unavailable:
        direct.region_summary(
            run_id,
            cool.artifact_id,
            generation=generation,
            revision=cool.revision,
            chrom="chr1",
            start=0,
            end=1000,
        )
    assert unavailable.value.status == 404
    prepared = direct.prepare(
        run_id,
        peak.artifact_id,
        generation=generation,
        revision=peak.revision,
    )
    counted = direct.region_summary(
        run_id,
        peak.artifact_id,
        generation=generation,
        revision=peak.revision,
        chrom="chr1",
        start=0,
        end=100,
    )
    empty = direct.region_summary(
        run_id,
        peak.artifact_id,
        generation=generation,
        revision=peak.revision,
        chrom="chr1",
        start=1000,
        end=2000,
    )
    assert prepared["preparation"] == "ready"
    assert counted["records"] == 1
    assert counted["complete"] is True
    assert empty["records"] == 0
    assert "zero-signal" not in str(counted["detail"])
    assert "not a zero-signal" in str(empty["detail"])

    real = app.state.run_service

    class _HiddenArtifacts:
        def __getattr__(self, name: str):
            return getattr(real, name)

        def list_artifacts(self, *_args, **_kwargs):
            raise RuntimeError("sqlite:///secret.db")

    hidden = BrowserTrackService(
        run_service=_HiddenArtifacts(),
        download_service=app.state.artifact_download_service,
        workspace_root=app.state.workspace_root,
        reference_profile_service=app.state.reference_profile_service,
    )
    with pytest.raises(BrowserTrackError) as listing:
        hidden.session(run_id)
    assert listing.value.status == 500
    assert listing.value.code == "BROWSER_TRACK_UNAVAILABLE"
    assert "secret" not in listing.value.message
    assert "sqlite" not in listing.value.message


def test_byte_parser_rejects_empty_and_inverted_ranges() -> None:
    assert parse_byte_range(None, 4) is None
    with pytest.raises(BrowserByteError, match="invalid_range"):
        parse_byte_range("bytes=-", 0)
    with pytest.raises(BrowserByteError, match="unsatisfiable_range"):
        parse_byte_range("bytes=3-1", 4)
    selected = parse_byte_range("bytes=-2", 4)
    assert selected is not None
    assert (selected.start, selected.end) == (2, 3)

    class _Plan:
        def __init__(self) -> None:
            self._content = b"ab"
            self.size_bytes = 4
            self.closed = False

        def close(self) -> None:
            self.closed = True

    plan = _Plan()
    with pytest.raises(BrowserByteError, match="short_read"):
        list(iter_plan_bytes(plan, None))
    assert plan.closed
