"""Pinned A0-only configuration. No platform or user configuration is consumed."""

PERI_COMMIT = "c7d154c9985cbe08ba25e964cb7468d15014e6d0"
CONTRACT = "helixweave-peri-a0-closeout/1"
PORT = 18765
TOOLS = [
    {
        "name": "lookup",
        "description": "Return one synthetic object in this instance.",
        "inputSchema": {
            "type": "object",
            "properties": {"object_id": {"type": "string", "enum": ["item-1"]}},
            "required": ["object_id"],
            "additionalProperties": False,
        },
    },
    {
        "name": "increment",
        "description": "Increment an in-memory synthetic counter after approval.",
        "inputSchema": {
            "type": "object",
            "properties": {},
            "additionalProperties": False,
        },
    },
]


def settings():
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
                    "apiKey": "synthetic-a0-key",
                    "baseUrl": f"http://127.0.0.1:{PORT}/v1",
                }
            ],
            "profiles": {
                name: {
                    "provider": "stub",
                    "model": "a0-stub",
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


def process_policy():
    entries = []
    for name in ["SearchExtraTools", "ExecuteExtraTool"]:
        entries.append(
            dict(
                name=name,
                server=None,
                wire=None,
                namespace="meta",
                builtin=None,
                aliases=[],
                parameters=None,
            )
        )
    for name in ["lookup", "increment"]:
        entries.append(
            dict(
                name="mcp__fixture__" + name,
                server="fixture",
                wire=name,
                namespace=None,
                builtin=None,
                aliases=[],
                parameters=next(
                    tool["inputSchema"] for tool in TOOLS if tool["name"] == name
                ),
            )
        )
    return {"version": 1, "tools": entries}
