"""Host broker unit tests; these are not evidence about native Peri execution."""

from collections import Counter
import copy
from configuration import TOOLS
import json
import threading
import unittest

from fixtures import MARKERS, hostile_document, project_fixture
from protocol import Client, Rejected, grant_key


def client_shell():
    client = Client.__new__(Client)
    client.instance, client.generation, client.session = "a", 2, "session-a"
    client.server_id = "server-a"
    client.state, client.approve, client.active_turn = "ready", True, True
    client.gate = threading.RLock()
    client.grants, client.used_requests = Counter(), set()
    client.counter = client.lookups = client.tool_arrivals = 0
    client.approvals, client.denied_callbacks, client.events, client.errors = (
        [],
        [],
        [],
        [],
    )
    client.approval_failure = client.callback_failure = None
    client.output_fixture = None
    client.declared_tools = copy.deepcopy(TOOLS)
    client.process = None
    return client


def tool_request(**updates):
    value = {
        "serverId": "server-a",
        "requestId": "request-1",
        "method": "tools/call",
        "params": {"name": "increment", "arguments": {}},
    }
    value.update(updates)
    return value


class ContractTests(unittest.TestCase):
    def test_tool_requires_active_exact_single_permission(self):
        client = client_shell()
        request = tool_request()
        client.callback("mcp/message", request, "a", 2)
        self.assertEqual(client.counter, 0)
        client.grants[grant_key("mcp__fixture__increment", {})] = 1
        client.callback("mcp/message", request, "a", 2)
        self.assertEqual(client.counter, 1)
        client.callback("mcp/message", request, "a", 2)
        self.assertEqual(client.counter, 1)

    def test_scope_generation_unknown_arguments_and_closed_reject(self):
        client = client_shell()
        client.grants[grant_key("mcp__fixture__increment", {})] = 5
        for instance, generation, request in [
            ("b", 2, tool_request()),
            ("a", 1, tool_request()),
            ("a", 2, tool_request(sessionId="session-b")),
            ("a", 2, tool_request(serverId="server-b")),
            (
                "a",
                2,
                tool_request(params={"name": "increment", "arguments": {"extra": 1}}),
            ),
        ]:
            client.callback("mcp/message", request, instance, generation)
        self.assertEqual(client.counter, 0)
        client.invalidate("SYNTHETIC_FAILURE")
        client.callback("mcp/message", tool_request(), "a", 2)
        self.assertEqual(client.counter, 0)
        self.assertFalse(client.grants)
        self.assertEqual(client.state, "closed")

    def test_missing_session_permission_and_no_active_turn_cannot_grant(self):
        client = client_shell()
        permission = {
            "toolCall": {"title": "mcp__fixture__increment", "rawInput": {}},
            "options": [{"optionId": "allow_once"}, {"optionId": "reject_once"}],
        }
        client.callback("session/request_permission", permission, "a", 2)
        self.assertFalse(client.grants)
        permission["sessionId"] = client.session
        client.active_turn = False
        client.callback("session/request_permission", permission, "a", 2)
        self.assertFalse(client.grants)

    def test_closed_projection_drops_free_text_even_when_named_description_or_default(
        self,
    ):
        projected = project_fixture(hostile_document())
        wire = json.dumps(projected)
        for marker in MARKERS:
            self.assertNotIn(marker, wire)
        self.assertEqual(projected["sample_count"], 2)
        self.assertNotIn("default", wire)
        self.assertNotIn("description", wire)

    def test_active_tools_list_rejects_changed_description_before_publication(self):
        client = client_shell()
        client.declared_tools[0]["description"] = "SYNTHETIC_UNAPPROVED_DESCRIPTION"
        reply = client.callback(
            "mcp/message",
            {
                "serverId": client.server_id,
                "requestId": "active-refresh",
                "method": "tools/list",
                "params": {},
            },
            "a",
            2,
        )
        self.assertEqual(reply["error"]["message"], "catalog changed")
        self.assertNotIn("SYNTHETIC_UNAPPROVED_DESCRIPTION", json.dumps(reply))
        self.assertEqual(client.state, "closed")
        self.assertFalse(client.grants)
        self.assertEqual(client.counter, 0)

    def test_projection_does_not_accept_identifier_in_numeric_field(self):
        value = hostile_document()
        value["summary"]["samples"] = MARKERS[2]
        value["schema"]["properties"]["threads"]["maximum"] = MARKERS[0]
        result = project_fixture(value)
        self.assertTrue(result["requires_local_input"])
        self.assertIsNone(result["sample_count"])
        self.assertNotIn(MARKERS[0], json.dumps(result))

    def test_arbitrary_acp_methods_and_unknown_prompt_fields_are_not_exposed(self):
        client = client_shell()
        client.active_turn = False
        for method in [
            "session/load",
            "session/fork",
            "session/resume",
            "session/set_mode",
            "session/update_config",
            "peri/mcp/control",
            "filesystem/read",
        ]:
            with self.assertRaises(Rejected):
                client.invoke(method, {"sessionId": client.session})
        with self.assertRaises(Rejected):
            client.invoke(
                "session/prompt",
                {"sessionId": client.session, "prompt": [], "mode": "bypass"},
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)
