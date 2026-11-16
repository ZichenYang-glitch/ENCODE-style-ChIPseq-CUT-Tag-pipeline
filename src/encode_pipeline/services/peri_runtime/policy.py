"""Bootstrap settings and the fixed A1 tool policy."""

from __future__ import annotations

from encode_pipeline.services.peri_runtime.identity import (
    CANONICAL_DRAFT_TOOL,
    CANONICAL_SCHEMA_TOOL,
    DRAFT_TOOL,
    SCHEMA_TOOL,
    SERVER_NAME,
)


PORT = 18765
SYNTHETIC_API_KEY = "synthetic-a1-key"
_EMPTY_SCHEMA = {
    "type": "object",
    "properties": {},
    "additionalProperties": False,
}
TOOLS = (
    {
        "name": SCHEMA_TOOL,
        "description": "Read the approved authoring field catalog for this session.",
        "inputSchema": _EMPTY_SCHEMA,
    },
    {
        "name": DRAFT_TOOL,
        "description": "Read the approved scalar projection of the current draft.",
        "inputSchema": _EMPTY_SCHEMA,
    },
)


def settings() -> dict[str, object]:
    disabled = [
        "DefaultSystemPromptMiddleware",
        "LangMiddleware",
        "AgentsMdMiddleware",
        "PluginMiddleware",
        "SkillsMiddleware",
        "SkillPreloadMiddleware",
        "AtMentionMiddleware",
        "ImageMiddleware",
        "GitAttributionMiddleware",
        "TodoMiddleware",
        "HookMiddleware",
        "HumanInTheLoopMiddleware",
        "SubAgentMiddleware",
        "GoalMiddleware",
        "BuiltInSubagents",
        "WebMiddleware",
        "ArtifactMiddleware",
        "CronMiddleware",
        "MonitorMiddleware",
        "WorkspaceMiddleware",
        "WorkflowMiddleware",
    ]
    return {
        "config": {
            "active_alias": "sonnet",
            "compact": {"auto_compact_enabled": False},
            "providers": [
                {
                    "id": "stub",
                    "type": "openai",
                    "apiKey": SYNTHETIC_API_KEY,
                    "baseUrl": f"http://127.0.0.1:{PORT}/v1",
                }
            ],
            "profiles": {
                name: {
                    "provider": "stub",
                    "model": "a1-stub",
                    "max_tokens": 1024,
                    "contextWindow": 32768,
                }
                for name in ["sonnet", "haiku", "opus", "fable"]
            },
            "meta_harness": {
                **dict.fromkeys(disabled, False),
                "PermissionMiddleware": True,
                "McpMiddleware": True,
                "ToolSearch": True,
            },
        }
    }


def process_policy() -> dict[str, object]:
    entries: list[dict[str, object]] = []
    for name in ["SearchExtraTools", "ExecuteExtraTool"]:
        entries.append(
            {
                "name": name,
                "server": None,
                "wire": None,
                "namespace": "meta",
                "builtin": None,
                "aliases": [],
                "parameters": None,
            }
        )
    for tool in TOOLS:
        wire = str(tool["name"])
        canonical = (
            CANONICAL_SCHEMA_TOOL if wire == SCHEMA_TOOL else CANONICAL_DRAFT_TOOL
        )
        entries.append(
            {
                "name": canonical,
                "server": SERVER_NAME,
                "wire": wire,
                "namespace": None,
                "builtin": None,
                "aliases": [],
                "parameters": tool["inputSchema"],
            }
        )
    return {"version": 1, "tools": entries}
