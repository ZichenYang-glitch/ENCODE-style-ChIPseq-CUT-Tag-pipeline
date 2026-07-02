"""Agent service for the workflow platform API.

The agent service exposes a chat endpoint that answers questions about
registered workflows. It operates only through read-only inspection tools and
a provider-neutral LLM client; it never executes workflows or mutates inputs.
"""

from __future__ import annotations

from encode_pipeline.api.models import AgentRequest, AgentResponse, AgentToolCall, IssueResponse
from encode_pipeline.platform.results import Issue, Result
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
    ) -> None:
        self._workflow_info = workflow_info
        self._validation_service = validation_service
        self._llm_client = llm_client
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
        """
        lookup = self._workflow_info.get_capabilities(workflow_id)
        if lookup.is_failure:
            return self._build_not_found_response(workflow_id, request, lookup.issues)

        issue_dicts = []
        if request.context is not None:
            issue_dicts = [
                issue.model_dump() for issue in request.context.current_issues
            ]

        explain_tool = self._tool_registry.get("explain_issues")
        tool_output = await explain_tool.handler(issues=issue_dicts)

        if isinstance(tool_output, Result):
            explanation = tool_output.value if tool_output.is_success else ""
        else:
            explanation = ""

        tool_call = AgentToolCall(
            tool_name=explain_tool.name,
            input_summary={"issues": issue_dicts},
            output_summary=str(explanation),
        )

        return AgentResponse(
            ok=True,
            session_id=request.session_id,
            message=explanation,
            tool_calls=[tool_call],
        )

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
                severity=issue.severity.value,
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
