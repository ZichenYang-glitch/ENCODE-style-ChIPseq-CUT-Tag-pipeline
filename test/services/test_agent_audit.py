from encode_pipeline.services.agent_audit import (
    AuditEvent,
    InMemoryAuditSink,
    NoOpAuditSink,
)


def test_records_event():
    sink = InMemoryAuditSink()
    event = AuditEvent(
        event_type="AGENT_REQUEST",
        workflow_id="wf-1",
        session_id="sess-1",
        issue_count=2,
        tool_name=None,
        filtered=False,
        redaction_counts={"sensitive_keys": 1, "paths": 0},
        issue_codes=("FRIP_LOW",),
        issue_sources=("validation",),
        issue_severities=("error",),
    )
    sink.record(event)
    assert sink.events() == (event,)


def test_bound_drops_oldest_events():
    sink = InMemoryAuditSink(bound=2)
    for i in range(3):
        sink.record(
            AuditEvent(
                event_type="AGENT_REQUEST",
                workflow_id=f"wf-{i}",
                session_id=None,
                issue_count=0,
                tool_name=None,
                filtered=False,
                redaction_counts={},
                issue_codes=(),
                issue_sources=(),
                issue_severities=(),
            )
        )
    events = sink.events()
    assert len(events) == 2
    assert events[0].workflow_id == "wf-1"
    assert events[1].workflow_id == "wf-2"


def test_audit_event_holds_no_raw_content():
    event = AuditEvent(
        event_type="AGENT_REQUEST",
        workflow_id="wf-1",
        session_id="sess-1",
        issue_count=1,
        tool_name=None,
        filtered=False,
        redaction_counts={"paths": 1},
        issue_codes=("FRIP_LOW",),
        issue_sources=("validation",),
        issue_severities=("error",),
    )
    home_literal = "/" + "home" + "/user"
    for forbidden in ("secret", "password", home_literal, "raw message"):
        assert forbidden not in str(event)


def test_no_op_audit_sink_records_nothing():
    sink = NoOpAuditSink()
    event = AuditEvent(
        event_type="AGENT_REQUEST",
        workflow_id="wf-1",
        session_id="sess-1",
        issue_count=1,
        tool_name=None,
        filtered=False,
        redaction_counts={"paths": 1},
        issue_codes=("FRIP_LOW",),
        issue_sources=("validation",),
        issue_severities=("error",),
    )
    sink.record(event)
    assert sink.events() == ()
