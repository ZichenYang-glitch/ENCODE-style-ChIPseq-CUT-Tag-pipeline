from __future__ import annotations

import os
from pathlib import Path

import pytest

import encode_pipeline.deployment.state as state_module
from encode_pipeline.deployment.errors import DeploymentError
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.models import DeploymentState, ENCODE_RUNTIME, PLATFORM
from encode_pipeline.deployment.state import (
    PLATFORM_ENV_FILENAME,
    STATE_FILENAME,
    StateStore,
    parse_platform_environment,
    render_platform_environment,
)


FIRST = "sha256-" + "1" * 64
SECOND = "sha256-" + "2" * 64
API_DIGEST = "a" * 64


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


@pytest.mark.parametrize(
    ("crash_point", "desired", "expected"),
    (
        ("transaction-prepared", "restore-prior", "prior"),
        ("generation-committed", "complete-candidate", "candidate"),
        ("pointer-committed", "restore-prior", "prior"),
    ),
)
def test_pending_state_transaction_recovery_is_identity_bound_and_atomic(
    tmp_path: Path,
    crash_point: str,
    desired: str,
    expected: str,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / crash_point)
    store = StateStore(layout)
    initial = store.initialize()
    candidate = initial.stage(PLATFORM, FIRST)

    def crash(point: str) -> None:
        if point == crash_point:
            raise RuntimeError("injected")

    with pytest.raises(RuntimeError, match="injected"):
        store.commit(
            candidate,
            operation="stage-platform",
            expected_current_identity=initial.identity,
            fault=crash,
        )

    selected = store.recover_pending_transaction(
        prior_state_identity=initial.identity,
        candidate_state_identity=candidate.identity,
        desired=desired,
    )

    assert selected == (initial if expected == "prior" else candidate)
    assert store.read() == selected
    assert store.pending_transactions() == ()
    suffix = "recovered.json" if desired == "restore-prior" else "complete.json"
    receipts = tuple(layout.state_transactions.glob(f"*.{suffix}"))
    assert receipts
    assert (layout.state_generations / initial.identity).is_dir()
    if crash_point != "transaction-prepared":
        assert (layout.state_generations / candidate.identity).is_dir()


def test_pending_state_transaction_recovery_rejects_journal_identity_mismatch(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    store = StateStore(layout)
    initial = store.initialize()
    candidate = initial.stage(PLATFORM, FIRST)

    def crash(point: str) -> None:
        if point == "generation-committed":
            raise RuntimeError("injected")

    with pytest.raises(RuntimeError, match="injected"):
        store.commit(
            candidate,
            operation="stage-platform",
            expected_current_identity=initial.identity,
            fault=crash,
        )

    with pytest.raises(DeploymentError) as captured:
        store.recover_pending_transaction(
            prior_state_identity=initial.identity,
            candidate_state_identity=SECOND,
            desired="complete-candidate",
        )

    assert captured.value.issue.code == "DEPLOYMENT_RECOVERY_IDENTITY_MISMATCH"
    assert store.read() == initial
    assert store.pending_transactions()


def test_pending_state_transaction_recovery_retries_after_pointer_crash(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    store = StateStore(layout)
    initial = store.initialize()
    candidate = initial.stage(PLATFORM, FIRST)

    def commit_crash(point: str) -> None:
        if point == "generation-committed":
            raise RuntimeError("commit-crash")

    with pytest.raises(RuntimeError, match="commit-crash"):
        store.commit(
            candidate,
            operation="stage-platform",
            expected_current_identity=initial.identity,
            fault=commit_crash,
        )

    def recovery_crash(point: str) -> None:
        if point == "recovery-pointer-committed":
            raise RuntimeError("recovery-crash")

    with pytest.raises(RuntimeError, match="recovery-crash"):
        store.recover_pending_transaction(
            prior_state_identity=initial.identity,
            candidate_state_identity=candidate.identity,
            desired="complete-candidate",
            fault=recovery_crash,
        )

    assert store.read() == candidate
    assert store.pending_transactions()
    assert (
        store.recover_pending_transaction(
            prior_state_identity=initial.identity,
            candidate_state_identity=candidate.identity,
            desired="complete-candidate",
        )
        == candidate
    )
    assert store.pending_transactions() == ()


def test_completed_state_commit_can_restore_prior_across_outer_journal_window(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "completed-handoff")
    store = StateStore(layout)
    initial = store.initialize()
    candidate = initial.stage(PLATFORM, FIRST)

    store.commit(
        candidate,
        operation="stage-platform",
        expected_current_identity=initial.identity,
    )

    assert store.pending_transactions() == ()
    assert len(tuple(layout.state_transactions.glob("*.complete.json"))) == 2

    def crash_after_pointer(point: str) -> None:
        if point == "recovery-pointer-committed":
            raise RuntimeError("outer-handoff-recovery-crash")

    with pytest.raises(RuntimeError, match="outer-handoff-recovery-crash"):
        store.recover_pending_transaction(
            prior_state_identity=initial.identity,
            candidate_state_identity=candidate.identity,
            desired="restore-prior",
            fault=crash_after_pointer,
        )

    assert store.read() == initial
    restored = store.recover_pending_transaction(
        prior_state_identity=initial.identity,
        candidate_state_identity=candidate.identity,
        desired="restore-prior",
    )
    assert restored == initial
    assert len(tuple(layout.state_transactions.glob("*.recovered.json"))) == 1
    restored_again = store.recover_pending_transaction(
        prior_state_identity=initial.identity,
        candidate_state_identity=candidate.identity,
        desired="restore-prior",
    )
    assert restored_again == initial


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


def _activate_platform(
    store: StateStore,
    layout: DeploymentLayout,
) -> DeploymentState:
    initial = store.initialize()
    staged = initial.stage(PLATFORM, FIRST)
    store.commit(
        staged,
        operation="stage-platform",
        expected_current_identity=initial.identity,
    )
    active = staged.activate(PLATFORM)
    store.commit(
        active,
        operation="activate-platform",
        expected_current_identity=staged.identity,
        platform_environment=render_platform_environment(
            layout,
            active,
            api_contract_sha256=API_DIGEST,
        ),
    )
    return active


def test_generation_atomically_binds_state_and_nonsecret_environment(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    store = StateStore(
        layout,
        reader_gid=os.getgid(),
        service_gid=os.getgid(),
    )
    active = _activate_platform(store, layout)
    generation = layout.state_generations / active.identity

    assert {entry.name for entry in generation.iterdir()} == {
        STATE_FILENAME,
        PLATFORM_ENV_FILENAME,
    }
    assert generation.stat().st_mode & 0o777 == 0o555
    assert (generation / STATE_FILENAME).stat().st_mode & 0o777 == 0o444
    environment = generation / PLATFORM_ENV_FILENAME
    assert environment.stat().st_mode & 0o777 == 0o440
    assert environment.stat().st_gid == os.getgid()
    parsed = parse_platform_environment(layout, active, environment.read_bytes())
    assert parsed.api_contract_sha256 == API_DIGEST


def test_metadata_reader_verifies_environment_metadata_without_reading_content(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    writer = StateStore(
        layout,
        reader_gid=os.getgid(),
        service_gid=os.getgid(),
    )
    active = _activate_platform(writer, layout)
    original = state_module.read_regular_file

    def guarded(path, **kwargs):
        if path.name == PLATFORM_ENV_FILENAME:
            raise AssertionError("metadata reader must not read platform.env")
        return original(path, **kwargs)

    monkeypatch.setattr(state_module, "read_regular_file", guarded)
    reader = StateStore(
        layout,
        reader_gid=os.getgid(),
        service_gid=os.getgid(),
        verify_environment_content=False,
    )

    assert (
        reader.read(
            expected_owner_uid=os.getuid(),
            expected_owner_gid=os.getgid(),
        )
        == active
    )
    with pytest.raises(DeploymentError) as captured:
        with reader.transaction(
            exclusive=True,
            expected_owner_uid=os.getuid(),
            expected_owner_gid=os.getgid(),
        ):
            pass
    assert captured.value.issue.code == "DEPLOYMENT_STATE_READ_ONLY"


def test_stage_rebinds_only_the_state_identity_in_the_active_environment(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    store = StateStore(layout)
    active = _activate_platform(store, layout)
    old_environment = (
        layout.state_generations / active.identity / PLATFORM_ENV_FILENAME
    ).read_text()
    staged = active.stage(ENCODE_RUNTIME, SECOND)

    store.commit(
        staged,
        operation="stage-encode-runtime",
        expected_current_identity=active.identity,
    )

    new_environment = (
        layout.state_generations / staged.identity / PLATFORM_ENV_FILENAME
    ).read_text()
    old_lines = dict(line.split("=", 1) for line in old_environment.splitlines())
    new_lines = dict(line.split("=", 1) for line in new_environment.splitlines())
    assert new_lines.pop("HELIXWEAVE_DEPLOYMENT_IDENTITY") == staged.identity
    assert old_lines.pop("HELIXWEAVE_DEPLOYMENT_IDENTITY") == active.identity
    assert new_lines == old_lines


@pytest.mark.parametrize(
    "mutator",
    (
        lambda content: content.replace(
            b"ENCODE_PIPELINE_QUEUE_NAME=encode-runs\n",
            b"UNKNOWN_KEY=value\n",
        ),
        lambda content: content.replace(
            b"ENCODE_PIPELINE_QUEUE_NAME=encode-runs\n",
            b"ENCODE_PIPELINE_DATABASE_URL=duplicate\n",
        ),
        lambda content: content.replace(FIRST.encode("ascii"), SECOND.encode("ascii")),
    ),
)
def test_environment_rejects_unknown_duplicate_or_active_platform_drift(
    tmp_path: Path,
    mutator,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    store = StateStore(layout)
    active = _activate_platform(store, layout)
    content = render_platform_environment(
        layout,
        active,
        api_contract_sha256=API_DIGEST,
    ).content

    with pytest.raises(DeploymentError) as captured:
        parse_platform_environment(layout, active, mutator(content))

    assert captured.value.issue.code == "DEPLOYMENT_CONFIGURATION_INVALID"


def test_interruption_before_pointer_cannot_select_candidate_environment(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    store = StateStore(layout)
    active = _activate_platform(store, layout)
    staged = active.stage(ENCODE_RUNTIME, SECOND)
    store.commit(
        staged,
        operation="stage-encode-runtime",
        expected_current_identity=active.identity,
    )
    candidate = staged.activate(ENCODE_RUNTIME)
    candidate_environment = render_platform_environment(
        layout,
        candidate,
        api_contract_sha256=API_DIGEST,
    )

    def crash(point: str) -> None:
        if point == "generation-committed":
            raise RuntimeError("injected")

    with pytest.raises(RuntimeError, match="injected"):
        store.commit(
            candidate,
            operation="activate-encode-runtime",
            expected_current_identity=staged.identity,
            platform_environment=candidate_environment,
            fault=crash,
        )

    assert store.read() == staged
    selected = os.readlink(layout.current_state_link)
    assert selected == f"generations/{staged.identity}"
    assert (
        candidate_environment.content
        == (
            layout.state_generations / candidate.identity / PLATFORM_ENV_FILENAME
        ).read_bytes()
    )
