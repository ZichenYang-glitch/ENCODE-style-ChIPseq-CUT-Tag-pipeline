import pytest


def test_run_status_values_and_terminal_detection():
    from encode_pipeline.platform.runs import RunStatus

    assert RunStatus.CREATED.value == "created"
    assert not RunStatus.CREATED.is_terminal
    assert RunStatus.SUCCEEDED.is_terminal
    assert RunStatus.FAILED.is_terminal
    assert RunStatus.CANCELLED.is_terminal


def test_can_transition_accepts_pr99_path():
    from encode_pipeline.platform.runs import RunStatus, can_transition

    assert can_transition(RunStatus.CREATED, RunStatus.VALIDATING)
    assert can_transition(RunStatus.VALIDATING, RunStatus.PLANNED)
    assert can_transition(RunStatus.VALIDATING, RunStatus.FAILED)
    assert can_transition(RunStatus.PLANNED, RunStatus.QUEUED)
    assert can_transition(RunStatus.QUEUED, RunStatus.RUNNING)
    assert can_transition(RunStatus.RUNNING, RunStatus.SUCCEEDED)
    assert can_transition(RunStatus.RUNNING, RunStatus.FAILED)
    assert can_transition(RunStatus.RUNNING, RunStatus.CANCELLED)


def test_can_transition_rejects_invalid_transitions():
    from encode_pipeline.platform.runs import RunStatus, can_transition

    assert not can_transition(RunStatus.CREATED, RunStatus.RUNNING)
    assert not can_transition(RunStatus.RUNNING, RunStatus.PLANNED)
    assert not can_transition(RunStatus.SUCCEEDED, RunStatus.FAILED)
    assert not can_transition(RunStatus.FAILED, RunStatus.CANCELLED)


def test_require_transition_raises_on_invalid():
    from encode_pipeline.platform.runs import RunStatus, require_transition

    with pytest.raises(ValueError, match="Illegal transition"):
        require_transition(RunStatus.SUCCEEDED, RunStatus.FAILED)
