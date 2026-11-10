"""Restricted authoring drafts, proposals, and validation orchestration."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from datetime import datetime
import re
from typing import Any

from encode_pipeline.adapters.encode_authoring_projection import (
    explainable_codes,
    explainable_issues,
    project_encode_fields,
)
from encode_pipeline.persistence.authoring import (
    AuthoringSessionRecord,
    AuthoringTurnRecord,
    SqlAlchemyAuthoringRepository,
    TurnCompletion,
    decode_authoring_cursor,
)
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.authoring import (
    ADMISSION_FAILED,
    AUTHENTICATION_REQUIRED,
    AuthoringActor,
    AuthoringField,
    AuthoringRejected,
    BRIDGE_UNAVAILABLE,
    REQUEST_TOO_LARGE,
    CAPABILITY_CHANGED,
    CAPABILITY_UNAVAILABLE,
    CONNECTION_LOST,
    DRAFT_EXPIRED,
    INVALID_REQUEST,
    MAX_MODEL_ANSWER_BYTES,
    MAX_PROJECTION_BYTES,
    MAX_PROPOSAL_FIELDS,
    MAX_QUESTION_BYTES,
    MODEL_OUTPUT_REJECTED,
    NOT_FOUND,
    PROPOSAL_EXPIRED,
    PROPOSAL_INVALIDATED,
    PROJECTION_VERSION,
    SAMPLE_SELECTION_REQUIRED,
    SCOPE_DENIED,
    STALE_DRAFT,
    WORKFLOW_UNAVAILABLE,
    adopt_request_digest,
    apply_field_changes,
    authoring_message,
    canonical_json,
    check_scalar,
    draft_digest,
    normalize_draft_content,
    parse_strict_json_object,
    projected_field,
    projection_identity_digest,
    proposal_digest,
    question_digest,
    sample_count,
    schema_identity_digest,
    sha256_text,
    tool_envelope,
    tool_issue,
    utc_now,
)
from encode_pipeline.platform.data_registry import (
    LEGACY_PROJECT_ID,
    ProjectKind,
    ProjectSampleSelection,
)
from encode_pipeline.platform.input_registry import InputFileRevisionSelection
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.platform.snapshots import ValidatedInputSnapshot
from encode_pipeline.services.authoring_bridge import (
    AuthoringBridge,
    BridgeDisconnected,
    BridgeSessionRequest,
    BridgeTurnRequest,
    UnavailableAuthoringBridge,
)
from encode_pipeline.services.data_registry import DataRegistryService


_INSTANCE_ID = re.compile(r"^[A-Za-z0-9_-]{8,36}$")
Clock = Callable[[], datetime]


@dataclass(frozen=True)
class SessionView:
    sid: str
    project_id: str
    workflow_id: str
    draft_state: str
    draft_revision: int
    expires_at: datetime
    purge_after: datetime
    schema_digest: str
    projection_digest: str | None
    draft: dict[str, Any] | None
    capability: dict[str, Any]
    connection_state: str
    connection_reason_code: str | None
    requires_new_connection: bool
    changed: bool | None = None


@dataclass(frozen=True)
class TurnView:
    turn_id: str
    client_turn_id: str
    status: str
    failure_reason: str | None
    result_kind: str | None
    proposal_id: str | None
    proposal: dict[str, Any] | None
    proposal_state: str | None
    adoptable: bool
    reason_code: str | None
    proposal_expires_at: datetime | None
    receipt: dict[str, Any] | None


@dataclass(frozen=True)
class AdoptionView:
    replayed: bool
    receipt: dict[str, Any]
    draft_state: str
    draft_revision: int
    expires_at: datetime
    purge_after: datetime


@dataclass(frozen=True)
class ProjectPage:
    projects: tuple[tuple[str, str], ...]
    next_cursor: str | None


@dataclass(frozen=True)
class SamplePage:
    project_id: str
    samples: tuple[tuple[str, str, str], ...]
    next_cursor: str | None


@dataclass(frozen=True)
class PreparedValidation:
    draft_revision: int
    result: Result[ValidatedInputSnapshot | None]


class AuthoringService:
    """Owner-scoped drafts. The bridge is optional and never bypasses this service."""

    def __init__(
        self,
        *,
        repository: SqlAlchemyAuthoringRepository,
        registry: Any,
        data_registry: DataRegistryService,
        validated_inputs: Any,
        snapshots: Any,
        bridge: AuthoringBridge | None = None,
        clock: Clock | None = None,
    ) -> None:
        self._repository = repository
        self._registry = registry
        self._data = data_registry
        self._validated = validated_inputs
        self._snapshots = snapshots
        self._bridge = bridge if bridge is not None else UnavailableAuthoringBridge()
        self._clock = clock
        self._capability_tokens: dict[tuple[str, int], str] = {}

    def list_projects(
        self,
        actor: AuthoringActor,
        *,
        cursor: str | None,
        limit: int,
    ) -> ProjectPage:
        self._require_login(actor)
        projects = [
            project
            for project in self._data.list_projects(include_archived=False)
            if project.kind is ProjectKind.USER
            and project.project_id != LEGACY_PROJECT_ID
        ]
        projects.sort(key=lambda project: (project.created_at, project.project_id))
        if cursor:
            created_at, project_id = _decode_pair_cursor(cursor, "prj_")
            projects = [
                project
                for project in projects
                if (
                    utc_now(lambda moment=project.created_at: moment),
                    project.project_id,
                )
                > (created_at, project_id)
            ]
        page = projects[:limit]
        next_cursor = None
        if len(projects) > limit and page:
            last = page[-1]
            next_cursor = _encode_pair_cursor(last.created_at, last.project_id)
        return ProjectPage(
            tuple((project.project_id, project.display_name) for project in page),
            next_cursor,
        )

    def list_samples(
        self,
        actor: AuthoringActor,
        project_id: str,
        *,
        cursor: str | None,
        limit: int,
    ) -> SamplePage:
        self._require_login(actor)
        self._require_visible_project(project_id)
        samples = list(self._data.list_samples(project_id))
        samples.sort(key=lambda view: view.sample.sample_id)
        if cursor:
            _stamp, sample_id = _decode_pair_cursor(cursor, "smp_")
            samples = [view for view in samples if view.sample.sample_id > sample_id]
        page = samples[:limit]
        next_cursor = None
        if len(samples) > limit and page:
            last = page[-1]
            next_cursor = _encode_pair_cursor(
                last.sample.created_at,
                last.sample.sample_id,
            )
        return SamplePage(
            project_id,
            tuple(
                (
                    view.sample.sample_id,
                    view.display_name,
                    view.revision.sample_revision_id,
                )
                for view in page
            ),
            next_cursor,
        )

    def create_session(
        self,
        actor: AuthoringActor,
        *,
        project_id: str,
        workflow_id: str,
    ) -> SessionView:
        self._require_visible_project(project_id)
        schema_digest, projection_digest, _fields = self._identity(workflow_id)
        record = self._repository.create_session(
            actor,
            project_id=project_id,
            workflow_id=workflow_id,
            schema_digest=schema_digest,
            projection_digest=projection_digest,
        )
        return self._session_view(actor, record, changed=False)

    def get_session(self, actor: AuthoringActor, sid: str) -> SessionView:
        return self._session_view(actor, self._repository.get_session(actor, sid))

    def list_sessions(
        self,
        actor: AuthoringActor,
        *,
        project_id: str | None,
        workflow_id: str | None,
        cursor: str | None,
        limit: int,
    ) -> tuple[tuple[SessionView, ...], str | None]:
        cursor_time = None
        cursor_sid = None
        if cursor:
            cursor_time, cursor_sid = decode_authoring_cursor(cursor)
        page = self._repository.list_sessions(
            actor,
            project_id=project_id,
            workflow_id=workflow_id,
            cursor_time=cursor_time,
            cursor_sid=cursor_sid,
            limit=limit,
        )
        return (
            tuple(self._session_view(actor, record) for record in page.records),
            page.next_cursor,
        )

    def save_draft(
        self,
        actor: AuthoringActor,
        sid: str,
        *,
        expected_revision: int,
        content: Mapping[str, Any],
    ) -> SessionView:
        current = self._repository.get_session(actor, sid)
        schema_digest, projection_digest, _fields = self._identity(current.workflow_id)
        normalized = normalize_draft_content(content)
        record, changed = self._repository.save_draft(
            actor,
            sid,
            expected_revision=expected_revision,
            content=normalized,
            schema_digest=schema_digest,
            projection_digest=projection_digest,
        )
        if (
            (
                record.schema_digest != current.schema_digest
                or record.projection_digest != current.projection_digest
            )
            and current.connection_state != "closed"
            and current.bridge_owner_epoch
        ):
            self._repository.fail_connection(
                sid,
                generation=current.generation,
                owner_epoch=current.bridge_owner_epoch,
                reason_code=CAPABILITY_CHANGED,
            )
            record = self._repository.get_session(actor, sid)
        return self._session_view(actor, record, changed=changed)

    def open_connection(self, actor: AuthoringActor, sid: str) -> SessionView:
        current = self._repository.get_session(actor, sid)
        self._require_active(current)
        if not self._bridge.status().available:
            raise AuthoringRejected(BRIDGE_UNAVAILABLE)
        schema_digest, projection_digest, _fields = self._identity(current.workflow_id)
        if projection_digest is None:
            raise AuthoringRejected(CAPABILITY_UNAVAILABLE)
        if (
            current.schema_digest != schema_digest
            or current.projection_digest != projection_digest
        ):
            if current.connection_state != "closed":
                self._repository.fail_connection(
                    sid,
                    generation=current.generation,
                    owner_epoch=current.bridge_owner_epoch or "",
                    reason_code=CAPABILITY_CHANGED,
                )
            current, _changed = self._repository.save_draft(
                actor,
                sid,
                expected_revision=current.revision,
                content=current.content,
                schema_digest=schema_digest,
                projection_digest=projection_digest,
            )
        epoch = self._bridge.owner_epoch()
        claim = self._repository.open_connection(
            actor,
            sid,
            owner_epoch=epoch,
            projection_digest=projection_digest,
        )
        if claim.reused or claim.capability_token is None:
            return self._session_view(actor, claim.record)
        self._capability_tokens[(sid, claim.record.generation)] = claim.capability_token
        request = BridgeSessionRequest(
            sid=sid,
            generation=claim.record.generation,
            owner_epoch=epoch,
            capability_token=claim.capability_token,
            schema_digest=schema_digest,
            projection_digest=projection_digest,
        )
        try:
            admission = self._bridge.admit(request)
            if (
                admission.state != "ready"
                or _INSTANCE_ID.fullmatch(admission.instance_id) is None
            ):
                raise AuthoringRejected(ADMISSION_FAILED)
            record = self._repository.mark_connection_ready(
                actor,
                sid,
                generation=claim.record.generation,
                instance_id=admission.instance_id,
                owner_epoch=epoch,
            )
        except AuthoringRejected:
            self._abandon(request, ADMISSION_FAILED)
            raise
        except BridgeDisconnected:
            self._abandon(request, ADMISSION_FAILED)
            raise AuthoringRejected(ADMISSION_FAILED) from None
        except Exception:
            self._abandon(request, ADMISSION_FAILED)
            raise AuthoringRejected(ADMISSION_FAILED) from None
        return self._session_view(actor, record)

    def close_connection(self, actor: AuthoringActor, sid: str) -> SessionView:
        before = self._repository.get_session(actor, sid)
        record = self._repository.close_connection(actor, sid)
        self._capability_tokens.pop((sid, before.generation), None)
        if before.generation != record.generation and before.instance_id is not None:
            self._bridge.close(
                BridgeSessionRequest(
                    sid=sid,
                    generation=before.generation,
                    owner_epoch=before.bridge_owner_epoch or "",
                    capability_token="",
                    schema_digest=before.schema_digest,
                    projection_digest=before.projection_digest or "",
                )
            )
        return self._session_view(actor, record)

    def start_turn(
        self,
        actor: AuthoringActor,
        sid: str,
        *,
        client_turn_id: str,
        expected_revision: int,
        question: Mapping[str, Any],
    ) -> TurnView:
        current = self._repository.get_session(actor, sid)
        self._require_active(current)
        if not self._bridge.status().available:
            raise AuthoringRejected(BRIDGE_UNAVAILABLE)
        schema_digest, projection_digest, fields = self._identity(current.workflow_id)
        if projection_digest is None or fields is None:
            raise AuthoringRejected(CAPABILITY_UNAVAILABLE)
        self._validate_question(question, fields)
        encoded = canonical_json(dict(question)).encode("utf-8")
        if len(encoded) > MAX_QUESTION_BYTES:
            raise AuthoringRejected(REQUEST_TOO_LARGE)
        digest = question_digest(
            question,
            expected_revision=expected_revision,
            schema_digest=schema_digest,
            projection_digest=projection_digest,
        )
        try:
            turn, replayed = self._repository.start_turn(
                actor,
                sid,
                client_turn_id=client_turn_id,
                request_digest=digest,
                question=question,
                expected_revision=expected_revision,
                schema_digest=schema_digest,
                projection_digest=projection_digest,
                owner_epoch=self._bridge.owner_epoch(),
            )
        except AuthoringRejected as exc:
            if (
                exc.reason_code == CAPABILITY_CHANGED
                and current.connection_state != "closed"
                and current.bridge_owner_epoch
            ):
                self._repository.fail_connection(
                    sid,
                    generation=current.generation,
                    owner_epoch=current.bridge_owner_epoch,
                    reason_code=CAPABILITY_CHANGED,
                )
            raise
        if replayed or turn.status != "generating":
            record = self._repository.get_session(actor, sid)
            return self._turn_view(actor, record, turn)
        issue_codes = self._issue_codes_for(current, question)
        request = BridgeTurnRequest(
            sid=sid,
            turn_id=turn.turn_id,
            generation=turn.generation,
            instance_id=turn.instance_id,
            owner_epoch=turn.bridge_owner_epoch,
            capability_token=self._capability_tokens.get((sid, turn.generation), ""),
            question=dict(question),
            issue_codes=issue_codes,
            schema_digest=schema_digest,
            projection_digest=projection_digest,
            base_revision=turn.base_revision,
        )
        try:
            raw = self._bridge.execute_turn(request)
            completion = self._completion(raw, turn, current, fields, issue_codes)
        except BridgeDisconnected:
            self._repository.fail_connection(
                sid,
                generation=turn.generation,
                owner_epoch=turn.bridge_owner_epoch,
                reason_code=CONNECTION_LOST,
            )
            record, failed = self._repository.get_turn(actor, sid, turn.turn_id)
            return self._turn_view(actor, record, failed)
        except AuthoringRejected as exc:
            if exc.reason_code == MODEL_OUTPUT_REJECTED:
                self._repository.fail_turn(
                    turn.turn_id,
                    generation=turn.generation,
                    reason_code=MODEL_OUTPUT_REJECTED,
                )
                record, failed = self._repository.get_turn(actor, sid, turn.turn_id)
                return self._turn_view(actor, record, failed)
            raise
        completed = self._repository.complete_turn(actor, turn.turn_id, completion)
        record = self._repository.get_session(actor, sid)
        return self._turn_view(actor, record, completed)

    def get_turn(self, actor: AuthoringActor, sid: str, turn_id: str) -> TurnView:
        record, turn = self._repository.get_turn(actor, sid, turn_id)
        return self._turn_view(actor, record, turn)

    def adopt(
        self,
        actor: AuthoringActor,
        sid: str,
        proposal_id: str,
        *,
        expected_revision: int,
        proposal_digest_value: str,
    ) -> AdoptionView:
        request_digest = adopt_request_digest(expected_revision, proposal_digest_value)
        current, turn = self._repository.get_proposal(actor, sid, proposal_id)
        content = None
        content_changed = None
        digest = None
        if turn.decision == "pending" and turn.proposal is not None:
            fields = self._fields_for(current.workflow_id)
            content = apply_field_changes(
                current.content,
                [
                    {
                        "field_token": change["field_token"],
                        "value": change["value"],
                    }
                    for change in turn.proposal.get("changes", [])
                ],
                fields,
            )
            digest = draft_digest(content)
            content_changed = digest != current.draft_digest
        result = self._repository.adopt(
            actor,
            sid,
            proposal_id,
            expected_revision=expected_revision,
            proposal_digest_value=proposal_digest_value,
            request_digest=request_digest,
            content=content,
            content_digest=digest,
            changed=content_changed,
        )
        return self._adoption_view(result.record, result.turn, result.replayed)

    def reject(
        self,
        actor: AuthoringActor,
        sid: str,
        proposal_id: str,
        *,
        proposal_digest_value: str,
    ) -> TurnView:
        turn = self._repository.reject(
            actor,
            sid,
            proposal_id,
            proposal_digest_value=proposal_digest_value,
        )
        record = self._repository.get_session(actor, sid)
        return self._turn_view(actor, record, turn)

    def validate(
        self,
        actor: AuthoringActor,
        sid: str,
        *,
        expected_revision: int,
    ) -> PreparedValidation:
        record = self._repository.load_for_validation(
            actor,
            sid,
            expected_revision=expected_revision,
        )
        cached = self._cached_validation(record)
        if cached is not None:
            return PreparedValidation(record.revision, cached)
        if not record.content["sample_revision_ids"]:
            return PreparedValidation(
                record.revision,
                Result.failure(
                    [
                        Issue(
                            code=SAMPLE_SELECTION_REQUIRED,
                            message=authoring_message(SAMPLE_SELECTION_REQUIRED),
                            path="sample_revision_ids",
                            source="authoring",
                        )
                    ]
                ),
            )
        result = self._validated.validate(
            record.workflow_id,
            WorkflowInputs(
                config=record.content["config"],
                samples=record.content["samples"],
                options=record.content["options"],
            ),
            project_sample_selection=ProjectSampleSelection(
                project_id=record.project_id,
                sample_revision_ids=tuple(record.content["sample_revision_ids"]),
            ),
            input_file_revision_selections=tuple(
                InputFileRevisionSelection(
                    input_use_key=item["input_use_key"],
                    occurrence=item["occurrence"],
                    input_file_revision_ids=tuple(item["input_file_revision_ids"]),
                )
                for item in record.content["input_selections"]
            ),
            reference_profile_revision_id=record.content[
                "reference_profile_revision_id"
            ],
        )
        schema_digest, _projection, _fields = self._identity(record.workflow_id)
        public = [_public_issue(issue) for issue in result.issues]
        codes = tuple(
            issue.code for issue in result.issues if issue.code in explainable_codes()
        )
        snapshot_id = None if result.value is None else result.value.snapshot_id
        self._repository.store_validation(
            actor,
            sid,
            expected_revision=expected_revision,
            schema_digest=schema_digest,
            snapshot_id=snapshot_id if result.is_success else None,
            issue_codes=codes,
            issues=public,
        )
        return PreparedValidation(expected_revision, result)

    def read_tool(
        self,
        *,
        sid: str,
        turn_id: str,
        generation: int,
        instance_id: str,
        owner_epoch: str,
        capability_token: str,
        name: str,
    ) -> dict[str, Any]:
        """Read-only tool body used by the bridge. Failures do not adopt."""

        record = self._repository.peek_session(sid)
        turn = self._repository.peek_turn(sid, turn_id)
        if record is None or turn is None:
            return tool_envelope(ok=False, data=None, issues=[tool_issue(SCOPE_DENIED)])
        if (
            record.connection_capability_digest != sha256_text(capability_token)
            or record.generation != generation
            or turn.generation != generation
            or record.instance_id != instance_id
            or record.bridge_owner_epoch != owner_epoch
            or record.connection_state != "ready"
        ):
            return tool_envelope(ok=False, data=None, issues=[tool_issue(SCOPE_DENIED)])
        if record.revision != turn.base_revision:
            self._repository.fail_turn(
                turn.turn_id,
                generation=turn.generation,
                reason_code=STALE_DRAFT,
            )
            return tool_envelope(ok=False, data=None, issues=[tool_issue(STALE_DRAFT)])
        schema_digest, projection_digest, fields = self._identity(record.workflow_id)
        if (
            fields is None
            or projection_digest != record.projection_digest
            or schema_digest != record.schema_digest
        ):
            self._repository.fail_connection(
                sid,
                generation=generation,
                owner_epoch=owner_epoch,
                reason_code=CAPABILITY_CHANGED,
            )
            return tool_envelope(
                ok=False,
                data=None,
                issues=[tool_issue(CAPABILITY_CHANGED)],
            )
        if name == "get_authoring_schema":
            return tool_envelope(
                ok=True,
                data={
                    "fields": [field.public_dict() for field in fields],
                    "projection_digest": projection_digest,
                    "schema_digest": schema_digest,
                    "workflow_id": record.workflow_id,
                },
            )
        if name != "get_draft_summary":
            return tool_envelope(ok=False, data=None, issues=[tool_issue(SCOPE_DENIED)])
        return tool_envelope(
            ok=True,
            data={
                "fields": [projected_field(record.content, field) for field in fields],
                "issue_codes": list(record.validation_issue_codes),
                "projection_digest": projection_digest,
                "revision": record.revision,
                "sample_count": sample_count(record.content.get("samples")),
                "schema_digest": schema_digest,
            },
        )

    def _completion(
        self,
        raw: str,
        turn: AuthoringTurnRecord,
        record: AuthoringSessionRecord,
        fields: Mapping[str, AuthoringField],
        issue_codes: tuple[str, ...],
    ) -> TurnCompletion:
        parsed = parse_strict_json_object(raw, max_bytes=MAX_MODEL_ANSWER_BYTES)
        expected = {
            "base_revision",
            "changes",
            "kind",
            "reason_codes",
            "schema_digest",
        }
        if set(parsed) != expected:
            raise AuthoringRejected(MODEL_OUTPUT_REJECTED)
        if parsed["schema_digest"] != turn.schema_digest:
            raise AuthoringRejected(MODEL_OUTPUT_REJECTED)
        if parsed["base_revision"] != turn.base_revision:
            raise AuthoringRejected(MODEL_OUTPUT_REJECTED)
        if parsed["kind"] not in {"proposal", "explanation"}:
            raise AuthoringRejected(MODEL_OUTPUT_REJECTED)
        changes = parsed["changes"]
        reasons = parsed["reason_codes"]
        if not isinstance(changes, list) or not isinstance(reasons, list):
            raise AuthoringRejected(MODEL_OUTPUT_REJECTED)
        if len(changes) > MAX_PROPOSAL_FIELDS or len(reasons) > MAX_PROPOSAL_FIELDS:
            raise AuthoringRejected(MODEL_OUTPUT_REJECTED)
        if any(
            not isinstance(code, str) or code not in issue_codes for code in reasons
        ):
            raise AuthoringRejected(MODEL_OUTPUT_REJECTED)
        normalized_changes = []
        for change in changes:
            if not isinstance(change, dict) or set(change) != {"field_token", "value"}:
                raise AuthoringRejected(MODEL_OUTPUT_REJECTED)
            normalized_changes.append(
                {"field_token": change["field_token"], "value": change["value"]}
            )
        if parsed["kind"] == "explanation":
            if normalized_changes:
                raise AuthoringRejected(MODEL_OUTPUT_REJECTED)
            digest = proposal_digest(
                sid=turn.sid,
                base_revision=turn.base_revision,
                schema_digest=turn.schema_digest,
                projection_digest=turn.projection_digest,
                changes=[],
                reason_codes=list(reasons),
            )
            return TurnCompletion("explanation", digest, None, None, None)
        updated = apply_field_changes(record.content, normalized_changes, fields)
        presented = []
        for change in normalized_changes:
            field = fields[change["field_token"]]
            presented.append(
                {
                    "after": projected_field(updated, field),
                    "before": projected_field(record.content, field),
                    "field_token": change["field_token"],
                    "value": change["value"],
                }
            )
        digest = proposal_digest(
            sid=turn.sid,
            base_revision=turn.base_revision,
            schema_digest=turn.schema_digest,
            projection_digest=turn.projection_digest,
            changes=[
                {"field_token": item["field_token"], "value": item["value"]}
                for item in presented
            ],
            reason_codes=list(reasons),
        )
        from secrets import token_hex

        return TurnCompletion(
            "proposal",
            digest,
            f"prop_{token_hex(16)}",
            {"changes": presented, "reason_codes": list(reasons)},
            digest,
        )

    def _validate_question(
        self,
        question: Mapping[str, Any],
        fields: Mapping[str, AuthoringField],
    ) -> None:
        if set(question) != {"field_tokens", "intent", "preferences"}:
            raise AuthoringRejected(INVALID_REQUEST)
        if question["intent"] not in {"suggest", "explain_issues"}:
            raise AuthoringRejected(INVALID_REQUEST)
        tokens = question["field_tokens"]
        preferences = question["preferences"]
        if (
            not isinstance(tokens, list)
            or not isinstance(preferences, list)
            or len(tokens) > MAX_PROPOSAL_FIELDS
            or len(set(tokens)) != len(tokens)
        ):
            raise AuthoringRejected(INVALID_REQUEST)
        for token in tokens:
            if token not in fields:
                raise AuthoringRejected(INVALID_REQUEST)
        for preference in preferences:
            if not isinstance(preference, dict) or set(preference) != {
                "field_token",
                "value",
            }:
                raise AuthoringRejected(INVALID_REQUEST)
            token = preference["field_token"]
            if token not in tokens:
                raise AuthoringRejected(INVALID_REQUEST)
            try:
                check_scalar(fields[token].constraint, preference["value"])
            except AuthoringRejected:
                raise AuthoringRejected(INVALID_REQUEST) from None

    def _identity(
        self,
        workflow_id: str,
    ) -> tuple[str, str | None, dict[str, AuthoringField] | None]:
        try:
            schema = self._registry.get(workflow_id).schema()
        except Exception:
            raise AuthoringRejected(WORKFLOW_UNAVAILABLE) from None
        payload = schema.to_dict()
        schema_digest = schema_identity_digest(workflow_id, payload)
        projected = project_encode_fields(workflow_id, payload)
        if projected is None:
            return schema_digest, None, None
        issues = explainable_issues()
        digest = projection_identity_digest(workflow_id, projected, issues)
        if (
            len(
                canonical_json(
                    {
                        "fields": [field.public_dict() for field in projected],
                        "issues": [issue.public_dict() for issue in issues],
                    }
                ).encode("utf-8")
            )
            > MAX_PROJECTION_BYTES
        ):
            return schema_digest, None, None
        return schema_digest, digest, {field.token: field for field in projected}

    def _fields_for(self, workflow_id: str) -> dict[str, AuthoringField]:
        _schema, _projection, fields = self._identity(workflow_id)
        if fields is None:
            raise AuthoringRejected(CAPABILITY_UNAVAILABLE)
        return fields

    def _session_view(
        self,
        actor: AuthoringActor,
        record: AuthoringSessionRecord,
        *,
        changed: bool | None = None,
    ) -> SessionView:
        now = utc_now(self._clock)
        active = now < record.draft_expires_at
        schema_digest, projection_digest, fields = self._safe_identity(
            record.workflow_id
        )
        capability = _capability(fields, projection_digest)
        ready = (
            record.connection_state == "ready"
            and record.originating_auth_session_ref == actor.session_digest
            and self._bridge.status().available
            and record.bridge_owner_epoch == _safe_epoch(self._bridge)
            and self._repository.login_is_live(actor.user_id, actor.session_digest)
        )
        return SessionView(
            sid=record.sid,
            project_id=record.project_id,
            workflow_id=record.workflow_id,
            draft_state="active" if active else "expired",
            draft_revision=record.revision,
            expires_at=record.draft_expires_at,
            purge_after=record.purge_after,
            schema_digest=schema_digest,
            projection_digest=projection_digest,
            draft=dict(record.content) if active else None,
            capability=capability,
            connection_state=record.connection_state,
            connection_reason_code=record.connection_reason_code,
            requires_new_connection=not ready,
            changed=changed,
        )

    def _turn_view(
        self,
        actor: AuthoringActor,
        record: AuthoringSessionRecord,
        turn: AuthoringTurnRecord,
    ) -> TurnView:
        now = utc_now(self._clock)
        expired = now >= record.draft_expires_at
        adoptable, reason = self._adoptable(actor, record, turn, now)
        hide = expired and turn.decision != "adopted"
        receipt = None
        if turn.decision == "adopted":
            receipt = _receipt(turn)
        proposal = None if hide or turn.proposal is None else dict(turn.proposal)
        return TurnView(
            turn_id=turn.turn_id,
            client_turn_id=turn.client_turn_id,
            status=turn.status,
            failure_reason=turn.failure_reason,
            result_kind=None if hide else turn.result_kind,
            proposal_id=None if hide else turn.proposal_id,
            proposal=proposal,
            proposal_state=turn.decision,
            adoptable=adoptable,
            reason_code=reason,
            proposal_expires_at=None if hide else turn.proposal_expires_at,
            receipt=receipt,
        )

    def _adoption_view(
        self,
        record: AuthoringSessionRecord,
        turn: AuthoringTurnRecord,
        replayed: bool,
    ) -> AdoptionView:
        now = utc_now(self._clock)
        receipt = _receipt(turn)
        if receipt is None:
            raise AuthoringRejected(PROPOSAL_INVALIDATED)
        return AdoptionView(
            replayed=replayed,
            receipt=receipt,
            draft_state="active" if now < record.draft_expires_at else "expired",
            draft_revision=record.revision,
            expires_at=record.draft_expires_at,
            purge_after=record.purge_after,
        )

    def _adoptable(
        self,
        actor: AuthoringActor,
        record: AuthoringSessionRecord,
        turn: AuthoringTurnRecord,
        now: datetime,
    ) -> tuple[bool, str | None]:
        if turn.status != "completed" or turn.result_kind != "proposal":
            return False, turn.failure_reason
        if turn.decision == "adopted":
            return False, None
        if turn.decision == "rejected":
            return False, None
        if turn.decision == "invalidated":
            return False, turn.invalidation_reason or PROPOSAL_INVALIDATED
        if now >= record.draft_expires_at:
            return False, DRAFT_EXPIRED
        if turn.proposal_expires_at is None or now >= turn.proposal_expires_at:
            return False, PROPOSAL_EXPIRED
        if actor.session_digest != turn.originating_auth_session_ref:
            return False, PROPOSAL_INVALIDATED
        if not self._repository.login_is_live(
            actor.user_id,
            turn.originating_auth_session_ref,
        ):
            return False, PROPOSAL_INVALIDATED
        if record.revision != turn.base_revision:
            return False, STALE_DRAFT
        if (
            record.schema_digest != turn.schema_digest
            or record.projection_digest != turn.projection_digest
            or record.connection_state != "ready"
            or record.generation != turn.generation
            or record.instance_id != turn.instance_id
            or record.bridge_owner_epoch != turn.bridge_owner_epoch
        ):
            return False, PROPOSAL_INVALIDATED
        return True, None

    def _cached_validation(
        self,
        record: AuthoringSessionRecord,
    ) -> Result[ValidatedInputSnapshot | None] | None:
        if (
            record.last_snapshot_id is None
            or record.last_validated_revision != record.revision
        ):
            return None
        schema_digest, _projection, _fields = self._identity(record.workflow_id)
        if schema_digest != record.schema_digest:
            return None
        try:
            snapshot = self._snapshots.get_validated_input_snapshot(
                record.last_snapshot_id
            )
        except KeyError:
            return None
        if snapshot.workflow_id != record.workflow_id:
            return None
        if utc_now(self._clock) >= snapshot.expires_at:
            return None
        issues = tuple(_issue_from_public(issue) for issue in record.validation_issues)
        if any(issue.severity == "error" for issue in issues):
            return None
        return Result.success(snapshot, issues)

    def _issue_codes_for(
        self,
        record: AuthoringSessionRecord,
        question: Mapping[str, Any],
    ) -> tuple[str, ...]:
        if question.get("intent") != "explain_issues":
            return ()
        if record.last_validated_revision != record.revision:
            return ()
        allowed = explainable_codes()
        return tuple(code for code in record.validation_issue_codes if code in allowed)

    def _require_active(self, record: AuthoringSessionRecord) -> None:
        if utc_now(self._clock) >= record.draft_expires_at:
            raise AuthoringRejected(DRAFT_EXPIRED)

    def _require_login(self, actor: AuthoringActor) -> None:
        if not self._repository.login_is_live(actor.user_id, actor.session_digest):
            raise AuthoringRejected(AUTHENTICATION_REQUIRED)

    def _require_visible_project(self, project_id: str) -> None:
        try:
            project = self._data.get_project(project_id)
        except KeyError:
            raise AuthoringRejected(NOT_FOUND) from None
        if (
            project.kind is not ProjectKind.USER
            or project.archived_at is not None
            or project.project_id == LEGACY_PROJECT_ID
        ):
            raise AuthoringRejected(NOT_FOUND)

    def _safe_identity(
        self,
        workflow_id: str,
    ) -> tuple[str, str | None, dict[str, AuthoringField] | None]:
        try:
            return self._identity(workflow_id)
        except AuthoringRejected:
            return "", None, None

    def _abandon(self, request: BridgeSessionRequest, reason_code: str) -> None:
        self._capability_tokens.pop((request.sid, request.generation), None)
        self._repository.fail_connection(
            request.sid,
            generation=request.generation,
            owner_epoch=request.owner_epoch,
            reason_code=reason_code,
        )
        try:
            self._bridge.close(request)
        except Exception:
            return


def _capability(
    fields: Mapping[str, AuthoringField] | None,
    projection_digest: str | None,
) -> dict[str, Any]:
    if fields is None or projection_digest is None:
        return {
            "allowed_intents": [],
            "allowed_issue_codes": [],
            "available": False,
            "fields": [],
            "projection_digest": None,
            "projection_version": PROJECTION_VERSION,
            "reason_code": CAPABILITY_UNAVAILABLE,
        }
    return {
        "allowed_intents": ["explain_issues", "suggest"],
        "allowed_issue_codes": [issue.public_dict() for issue in explainable_issues()],
        "available": True,
        "fields": [field.public_dict() for field in fields.values()],
        "projection_digest": projection_digest,
        "projection_version": PROJECTION_VERSION,
        "reason_code": None,
    }


def _receipt(turn: AuthoringTurnRecord) -> dict[str, Any] | None:
    if turn.decision != "adopted" or turn.proposal_id is None:
        return None
    return {
        "adopted_at": turn.adopted_at,
        "base_revision": turn.base_revision,
        "changed": bool(turn.adoption_changed),
        "draft_expires_at": turn.adoption_draft_expires_at,
        "proposal_id": turn.proposal_id,
        "result_revision": turn.result_revision,
    }


def _public_issue(issue: Issue) -> dict[str, Any]:
    payload = issue.to_dict()
    payload["technical_message"] = None
    return payload


def _issue_from_public(payload: Mapping[str, Any]) -> Issue:
    return Issue(
        code=str(payload["code"]),
        message=str(payload["message"]),
        severity=str(payload.get("severity") or "error"),
        path=payload.get("path") if isinstance(payload.get("path"), str) else None,
        source=payload.get("source")
        if isinstance(payload.get("source"), str)
        else None,
        technical_message=None,
        hint=payload.get("hint") if isinstance(payload.get("hint"), str) else None,
        context=payload.get("context")
        if isinstance(payload.get("context"), dict)
        else {},
    )


def _safe_epoch(bridge: AuthoringBridge) -> str | None:
    try:
        if not bridge.status().available:
            return None
        return bridge.owner_epoch()
    except Exception:
        return None


def _encode_pair_cursor(moment: datetime, identifier: str) -> str:
    import base64

    raw = f"{utc_now(lambda: moment).strftime('%Y-%m-%dT%H:%M:%SZ')}\n{identifier}"
    return base64.urlsafe_b64encode(raw.encode()).decode("ascii").rstrip("=")


def _decode_pair_cursor(value: str, prefix: str) -> tuple[datetime, str]:
    import base64

    if not isinstance(value, str) or len(value) > 256:
        raise AuthoringRejected(INVALID_REQUEST)
    padded = value + ("=" * (-len(value) % 4))
    try:
        raw = base64.urlsafe_b64decode(padded.encode("ascii")).decode("utf-8")
        stamp, identifier = raw.split("\n", 1)
        parsed = datetime.strptime(stamp, "%Y-%m-%dT%H:%M:%SZ").replace(
            tzinfo=__import__("datetime").timezone.utc
        )
    except (ValueError, UnicodeError):
        raise AuthoringRejected(INVALID_REQUEST) from None
    if not identifier.startswith(prefix):
        raise AuthoringRejected(INVALID_REQUEST)
    return parsed, identifier
