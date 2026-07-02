"""Read-only agent tool registry and tool builders.

This module provides the registry used by the agent service to expose
read-only workflow inspection tools. Mutating or execution-oriented tools
are rejected at registration time via a name-prefix deny list.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

from encode_pipeline.platform.adapters import WorkflowInputs, WorkflowMetadata
from encode_pipeline.platform.results import Result
from encode_pipeline.services.llm_client import LLMClient, LLMMessage
from encode_pipeline.services.validation import ValidationService
from encode_pipeline.services.workflow_info import WorkflowInfoService


@dataclass(frozen=True)
class AgentTool:
    """Definition of a tool exposed to the agent runtime."""

    name: str
    description: str
    parameters_schema: Mapping[str, object]
    handler: Callable[..., Any]
    read_only: Literal[True] = True


class ReadOnlyToolRegistry:
    """Registry that only accepts read-only agent tools.

    Tools with names matching a mutation/execution prefix are rejected even
    when marked ``read_only=True``. The registry instance owns its internal
    state; there is no hidden global registry.
    """

    _DENIED_PREFIXES: tuple[str, ...] = (
        "run_",
        "submit_",
        "execute_",
        "delete_",
        "update_",
        "apply_",
        "write_",
        "modify_",
    )

    def __init__(self) -> None:
        self._tools: dict[str, AgentTool] = {}

    def register(self, tool: AgentTool) -> None:
        """Register a read-only tool.

        Raises:
            ValueError: If the tool is not read-only, its name starts with a
                denied prefix, or its name is already registered.
        """
        if tool.read_only is not True:
            raise ValueError(
                f"Tool {tool.name!r} must have read_only=True to be registered"
            )

        if any(tool.name.startswith(prefix) for prefix in self._DENIED_PREFIXES):
            raise ValueError(
                f"Tool {tool.name!r} starts with a denied prefix and cannot be "
                "registered in a read-only registry"
            )

        if tool.name in self._tools:
            raise ValueError(f"Tool {tool.name!r} is already registered")

        self._tools[tool.name] = tool

    def get(self, name: str) -> AgentTool:
        """Return the tool with the given name.

        Raises:
            KeyError: If no tool with that name is registered.
        """
        try:
            return self._tools[name]
        except KeyError as exc:
            raise KeyError(name) from exc

    def list_tools(self) -> list[AgentTool]:
        """Return all registered tools in registration order."""
        return list(self._tools.values())


def build_list_workflows_tool(workflow_info: WorkflowInfoService) -> AgentTool:
    """Build a read-only tool that lists registered workflows."""

    def handler(**kwargs: Any) -> Result[list[WorkflowMetadata]]:
        return Result.success(workflow_info.list_workflows())

    return AgentTool(
        name="list_workflows",
        description="List all registered workflows with their metadata.",
        parameters_schema={"type": "object", "properties": {}, "required": []},
        handler=handler,
    )


def build_get_workflow_schema_tool(workflow_info: WorkflowInfoService) -> AgentTool:
    """Build a read-only tool that fetches a workflow schema."""

    def handler(workflow_id: str, **kwargs: Any) -> Result[object]:
        return workflow_info.get_schema(workflow_id)

    return AgentTool(
        name="get_workflow_schema",
        description="Return the input schema for a registered workflow.",
        parameters_schema={
            "type": "object",
            "properties": {
                "workflow_id": {
                    "type": "string",
                    "description": "Unique workflow identifier.",
                },
            },
            "required": ["workflow_id"],
        },
        handler=handler,
    )


def build_validate_workflow_input_tool(
    validation_service: ValidationService,
) -> AgentTool:
    """Build a read-only tool that validates workflow inputs."""

    def handler(
        workflow_id: str,
        config: Mapping[str, object],
        samples: str | Path | None = None,
        options: Mapping[str, object] | None = None,
        **kwargs: Any,
    ) -> Result[object]:
        inputs = WorkflowInputs(
            config=config,
            samples=samples,
            options=dict(options or {}),
        )
        return validation_service.validate(workflow_id, inputs)

    return AgentTool(
        name="validate_workflow_input",
        description="Validate inputs for a registered workflow.",
        parameters_schema={
            "type": "object",
            "properties": {
                "workflow_id": {
                    "type": "string",
                    "description": "Unique workflow identifier.",
                },
                "config": {
                    "type": "object",
                    "description": "Workflow configuration to validate.",
                },
                "samples": {
                    "type": "string",
                    "description": "Optional path to a sample sheet.",
                },
                "options": {
                    "type": "object",
                    "description": "Optional workflow options.",
                },
            },
            "required": ["workflow_id", "config"],
        },
        handler=handler,
    )


def build_explain_issues_tool(llm_client: LLMClient) -> AgentTool:
    """Build a read-only tool that explains validation issues via an LLM."""

    async def handler(
        issues: list[Mapping[str, object]],
        **kwargs: Any,
    ) -> Result[str]:
        messages = [
            LLMMessage(
                role="system",
                content="You are a helpful workflow assistant.",
            ),
            LLMMessage(
                role="user",
                content=f"Explain these workflow validation issues: {issues}",
            ),
        ]
        response = await llm_client.complete(messages, tools=[])
        return Result.success(response.content)

    return AgentTool(
        name="explain_issues",
        description="Explain workflow validation issues using the LLM.",
        parameters_schema={
            "type": "object",
            "properties": {
                "issues": {
                    "type": "array",
                    "description": "List of validation issues to explain.",
                    "items": {"type": "object"},
                },
            },
            "required": ["issues"],
        },
        handler=handler,
    )
