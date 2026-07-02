"""Agent service for the workflow platform API.

The agent service exposes a chat endpoint that answers questions about
registered workflows. It operates only through read-only inspection tools and
a provider-neutral LLM client; it never executes workflows or mutates inputs.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, cast
from encode_pipeline.api.models import AgentRequest, AgentResponse, AgentToolCall, IssueResponse
from encode_pipeline.platform.results import Issue, IssueSeverity, Result
from encode_pipeline.services.agent_audit import AgentAuditSink, AuditEvent, InMemoryAuditSink
from encode_pipeline.services.agent_output_filter import OutputFilter
from encode_pipeline.services.agent_redaction import RedactionPolicy
from encode_pipeline.services.agent_tools import (
    ReadOnlyToolRegistry,
    build_explain_issues_tool,
    build_get_workflow_schema_tool,
    build_list_workflows_tool,
    build_validate_workflow_input_tool,
)
from encode_pipeline.services.llm_client import LLMClient
from encode_pipeline.services.validation import ValidationService
from encode_pipeline.services.workflow_info import WorkflowInfoService


class AgentService:
    """Chat service backed by read-only workflow tools and an LLM client."""

    def __init__(
        self,
        workflow_info: WorkflowInfoService,
        validation_service: ValidationService,
        llm_client: LLMClient,
        tool_registry: ReadOnlyToolRegistry | None = None,
        redaction_policy: RedactionPolicy | None = None,
        output_filter: OutputFilter | None = None,
        audit_sink: AgentAuditSink | None = None,
    ) -> None:
        self._workflow_info = workflow_info
        self._validation_service = validation_service
        self._llm_client = llm_client
        self._redaction_policy = redaction_policy if redaction_policy is not None else RedactionPolicy()
        self._output_filter = output_filter if output_filter is not None else OutputFilter()
        self._audit_sink = audit_sink if audit_sink is not None else InMemoryAuditSink()
        if tool_registry is None:
            registry = ReadOnlyToolRegistry()
            registry.register(build_list_workflows_tool(workflow_info))
            registry.register(build_get_workflow_schema_tool(workflow_info))
            registry.register(build_validate_workflow_input_tool(validation_service))
            registry.register(build_explain_issues_tool(llm_client))
            self._tool_registry = registry
        else:
            self._tool_registry = tool_registry

    async def chat(self, workflow_id: str, request: AgentRequest) -> AgentResponse:
        """Answer a user message for the given workflow.

        Returns a failure envelope when the workflow is not registered.
        For valid workflows the response includes the LLM-generated
        explanation and a transparent record of the read-only
        ``explain_issues`` tool call that produced it.

        When safety components are configured, the service redacts sensitive
        values and absolute filesystem paths, filters execution-like output,
        and emits audit events for the request, tool call, and response.
        """
        issue_dicts = self._request_issue_dicts(request)
        redacted_issue_dicts, input_redaction = self._redact_issues(issue_dicts)

        self._record_audit_event(
            event_type="AGENT_REQUEST",
            workflow_id=workflow_id,
            request=request,
            issue_dicts=redacted_issue_dicts,
            redaction_counts=input_redaction,
        )

        lookup = self._workflow_info.get_capabilities(workflow_id)
        if lookup.is_failure:
            response = self._build_not_found_response(workflow_id, request, lookup.issues)
            response = self._redact_response(response)
            self._record_audit_event(
                event_type="LLM_RESPONSE",
                workflow_id=workflow_id,
                request=request,
                issue_dicts=redacted_issue_dicts,
                redaction_counts=self._response_redaction_counts(response),
            )
            return response

        explain_tool = self._tool_registry.get("explain_issues")
        tool_output = await explain_tool.handler(issues=redacted_issue_dicts)

        self._record_audit_event(
            event_type="TOOL_CALL",
            workflow_id=workflow_id,
            request=request,
            issue_dicts=redacted_issue_dicts,
            tool_name=explain_tool.name,
            redaction_counts=input_redaction,
        )

        if isinstance(tool_output, Result):
            raw_explanation = tool_output.value if tool_output.is_success else ""
            explanation = str(raw_explanation) if raw_explanation is not None else ""
        else:
            explanation = ""

        filtered, filtered_text, _matched_patterns = self._apply_output_filter(
            explanation
        )
        if filtered:
            self._record_audit_event(
                event_type="OUTPUT_FILTERED",
                workflow_id=workflow_id,
                request=request,
                issue_dicts=redacted_issue_dicts,
                filtered=True,
            )
            explanation = filtered_text

        explanation_result = self._redaction_policy.redact_text(explanation)
        explanation = cast(str, explanation_result.value)
        input_summary = cast(
            dict[str, Any],
            self._redact_value({"issues": redacted_issue_dicts}),
        )
        output_summary = self._redact_text(str(explanation))

        tool_call = AgentToolCall(
            tool_name=explain_tool.name,
            input_summary=input_summary,
            output_summary=output_summary,
        )

        response = AgentResponse(
            ok=True,
            session_id=request.session_id,
            message=explanation,
            tool_calls=[tool_call],
        )

        self._record_audit_event(
            event_type="LLM_RESPONSE",
            workflow_id=workflow_id,
            request=request,
            issue_dicts=redacted_issue_dicts,
            redaction_counts=self._response_redaction_counts(
                response, prior_message_paths=explanation_result.absolute_paths_redacted
            ),
        )
        return response

    def _build_not_found_response(
        self,
        workflow_id: str,
        request: AgentRequest,
        issues: tuple[Issue, ...],
    ) -> AgentResponse:
        issue_responses = [
            IssueResponse(
                code=issue.code,
                message=issue.message,
                severity=(
                    issue.severity.value
                    if isinstance(issue.severity, IssueSeverity)
                    else issue.severity
                ),
                path=issue.path,
                source=issue.source,
                technical_message=issue.technical_message,
                hint=issue.hint,
                context=dict(issue.context),
            )
            for issue in issues
        ]
        return AgentResponse(
            ok=False,
            session_id=request.session_id,
            message=f"Workflow '{workflow_id}' was not found.",
            issues=issue_responses,
        )

    def _request_issue_dicts(self, request: AgentRequest) -> list[dict[str, object]]:
        if request.context is None:
            return []
        return [issue.model_dump() for issue in request.context.current_issues]

    def _redact_issues(
        self, issue_dicts: list[dict[str, object]]
    ) -> tuple[list[dict[str, object]], dict[str, int]]:
        result = self._redaction_policy.redact_issues(
            cast(list[Mapping[str, object]], issue_dicts)
        )
        value = result.value
        assert isinstance(value, list)
        return cast(list[dict[str, object]], value), {
            "sensitive_keys": result.sensitive_keys_redacted,
            "paths": result.absolute_paths_redacted,
        }

    def _redact_text(self, text: str) -> str:
        result = self._redaction_policy.redact_text(text)
        assert isinstance(result.value, str)
        return result.value

    def _redact_value(self, value: object) -> object:
        return self._redaction_policy.redact_value(value).value

    def _apply_output_filter(self, text: str) -> tuple[bool, str, tuple[str, ...]]:
        result = self._output_filter.filter_text(text)
        return result.filtered, result.text, result.matched_patterns

    def _redact_response(self, response: AgentResponse) -> AgentResponse:
        redacted_issue_dicts = cast(
            list[dict[str, object]],
            self._redact_value([issue.model_dump() for issue in response.issues]),
        )
        redacted_issues = [
            IssueResponse(**cast(dict[str, Any], issue))
            for issue in redacted_issue_dicts
        ]
        return response.model_copy(
            update={
                "message": self._redact_text(response.message),
                "issues": redacted_issues,
            }
        )

    def _response_redaction_counts(
        self,
        response: AgentResponse,
        prior_message_paths: int = 0,
    ) -> dict[str, int]:
        message_result = self._redaction_policy.redact_text(response.message)
        issues_result = self._redaction_policy.redact_issues(
            [issue.model_dump() for issue in response.issues]
        )
        tool_input_paths = 0
        tool_input_sensitive = 0
        tool_output_paths = 0
        tool_output_sensitive = 0
        for call in response.tool_calls:
            input_result = self._redaction_policy.redact_value(call.input_summary)
            tool_input_paths += input_result.absolute_paths_redacted
            tool_input_sensitive += input_result.sensitive_keys_redacted
            output_result = self._redaction_policy.redact_text(call.output_summary)
            tool_output_paths += output_result.absolute_paths_redacted
            tool_output_sensitive += output_result.sensitive_keys_redacted
        return {
            "sensitive_keys": (
                message_result.sensitive_keys_redacted
                + issues_result.sensitive_keys_redacted
                + tool_input_sensitive
                + tool_output_sensitive
            ),
            "paths": (
                prior_message_paths
                + issues_result.absolute_paths_redacted
                + tool_input_paths
                + tool_output_paths
            ),
        }

    def _record_audit_event(
        self,
        *,
        event_type: str,
        workflow_id: str,
        request: AgentRequest,
        issue_dicts: list[dict[str, object]],
        tool_name: str | None = None,
        redaction_counts: dict[str, int] | None = None,
        filtered: bool = False,
    ) -> None:
        if redaction_counts is None:
            redaction_counts = {"sensitive_keys": 0, "paths": 0}
        self._audit_sink.record(
            AuditEvent(
                event_type=event_type,
                workflow_id=workflow_id,
                session_id=request.session_id,
                issue_count=len(issue_dicts),
                tool_name=tool_name,
                filtered=filtered,
                redaction_counts=redaction_counts,
                issue_codes=self._issue_codes(issue_dicts),
                issue_sources=self._issue_sources(issue_dicts),
                issue_severities=self._issue_severities(issue_dicts),
            )
        )

    @staticmethod
    def _issue_codes(issue_dicts: list[dict[str, object]]) -> tuple[str, ...]:
        return tuple(str(code) for code in (issue.get("code") for issue in issue_dicts) if code)

    @staticmethod
    def _issue_sources(issue_dicts: list[dict[str, object]]) -> tuple[str, ...]:
        return tuple(str(source) for source in (issue.get("source") for issue in issue_dicts) if source)

    @staticmethod
    def _issue_severities(issue_dicts: list[dict[str, object]]) -> tuple[str, ...]:
        return tuple(str(severity) for severity in (issue.get("severity") for issue in issue_dicts) if severity)
