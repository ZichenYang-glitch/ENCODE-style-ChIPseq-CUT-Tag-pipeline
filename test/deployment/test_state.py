from __future__ import annotations

import os
from pathlib import Path

import pytest

from encode_pipeline.deployment.errors import DeploymentError
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.models import DeploymentState, PLATFORM
from encode_pipeline.deployment.state import STATE_FILENAME, StateStore


FIRST = "sha256-" + "1" * 64
SECOND = "sha256-" + "2" * 64


def test_slot_transitions_preserve_explicit_previous_identity() -> None:
    state = DeploymentState.initial()
    staged = state.stage(PLATFORM, FIRST)
    active = staged.activate(PLATFORM)
    upgraded = active.stage(PLATFORM, SECOND).activate(PLATFORM)
    rolled_back = upgraded.rollback(PLATFORM)

    assert active.components[PLATFORM].active == FIRST
    assert upgraded.components[PLATFORM].active == SECOND
    assert upgraded.components[PLATFORM].previous == FIRST
    assert rolled_back.components[PLATFORM].active == FIRST
    assert rolled_back.components[PLATFORM].previous == SECOND


def test_state_rejects_duplicate_slot_identity() -> None:
    state = DeploymentState.initial().stage(PLATFORM, FIRST).activate(PLATFORM)
    with pytest.raises(DeploymentError) as captured:
        state.stage(PLATFORM, FIRST)
    assert captured.value.issue.code == "DEPLOYMENT_STAGE_CONFLICT"


def test_state_generation_pointer_switch_is_atomic(tmp_path: Path) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    store = StateStore(layout)
    initial = store.initialize()
    replacement = initial.stage(PLATFORM, FIRST)

    store.commit(
        replacement,
        operation="stage-platform",
        expected_current_identity=initial.identity,
    )

    assert store.read() == replacement
    assert os.readlink(layout.current_state_link) == (
        f"generations/{replacement.identity}"
    )
    assert (
        layout.state_generations / replacement.identity / STATE_FILENAME
    ).stat().st_mode & 0o777 == 0o444


def test_interruption_before_pointer_keeps_old_generation_active(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    store = StateStore(layout)
    initial = store.initialize()
    replacement = initial.stage(PLATFORM, FIRST)

    def crash(point: str) -> None:
        if point == "generation-committed":
            raise RuntimeError("injected")

    with pytest.raises(RuntimeError, match="injected"):
        store.commit(
            replacement,
            operation="stage-platform",
            expected_current_identity=initial.identity,
            fault=crash,
        )

    assert store.read() == initial
    assert store.pending_transactions()


def test_interruption_after_pointer_is_detectable_and_recoverable(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    store = StateStore(layout)
    initial = store.initialize()
    replacement = initial.stage(PLATFORM, FIRST)

    def crash(point: str) -> None:
        if point == "pointer-committed":
            raise RuntimeError("injected")

    with pytest.raises(RuntimeError, match="injected"):
        store.commit(
            replacement,
            operation="stage-platform",
            expected_current_identity=initial.identity,
            fault=crash,
        )

    assert store.read() == replacement
    assert store.pending_transactions()


def test_state_file_symlink_is_rejected(tmp_path: Path) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    store = StateStore(layout)
    state = store.initialize()
    generation = layout.state_generations / state.identity
    state_file = generation / STATE_FILENAME
    generation.chmod(0o755)
    state_file.chmod(0o644)
    state_file.unlink()
    state_file.symlink_to("/etc/passwd")

    with pytest.raises(DeploymentError) as captured:
        store.read()
    assert captured.value.issue.code == "DEPLOYMENT_STATE_INVALID"


def test_state_commit_compare_and_swap_rejects_a_stale_base(tmp_path: Path) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    store = StateStore(layout)
    initial = store.initialize()
    first = initial.stage(PLATFORM, FIRST)
    stale = initial.stage("encode-runtime", SECOND)

    store.commit(
        first,
        operation="stage-platform",
        expected_current_identity=initial.identity,
    )
    with pytest.raises(DeploymentError) as captured:
        store.commit(
            stale,
            operation="stage-encode-runtime",
            expected_current_identity=initial.identity,
        )

    assert captured.value.issue.code == "DEPLOYMENT_STATE_CHANGED"
    assert store.read() == first


def test_state_rejects_a_writable_operator_ancestor(tmp_path: Path) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    store = StateStore(layout)
    store.initialize(
        expected_owner_uid=os.getuid(),
        expected_owner_gid=os.getgid(),
    )
    layout.data_root.chmod(0o770)

    with pytest.raises(DeploymentError) as captured:
        store.read(
            expected_owner_uid=os.getuid(),
            expected_owner_gid=os.getgid(),
        )

    assert captured.value.issue.code == "DEPLOYMENT_STATE_INVALID"


def test_state_lock_fails_closed_under_contention(tmp_path: Path) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    store = StateStore(layout)
    store.initialize()

    with store.transaction(exclusive=True):
        with pytest.raises(DeploymentError) as captured:
            store.read()

    assert captured.value.issue.code == "DEPLOYMENT_BUSY"
