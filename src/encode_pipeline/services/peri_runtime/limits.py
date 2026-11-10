"""Cgroup checks for one native Peri process, not the API process."""

from __future__ import annotations

from pathlib import Path
import subprocess
import time

from encode_pipeline.services.peri_runtime.identity import MEMORY_MAX_BYTES, PID_MAX


def cgroup_path_is_qualified(root: Path) -> bool:
    try:
        memory = (root / "memory.max").read_text(encoding="utf-8").strip()
        swap = (root / "memory.swap.max").read_text(encoding="utf-8").strip()
        pids = (root / "pids.max").read_text(encoding="utf-8").strip()
        cpu = (root / "cpu.max").read_text(encoding="utf-8").strip()
    except (OSError, ValueError):
        return False
    if memory == "max" or int(memory) > MEMORY_MAX_BYTES or swap != "0":
        return False
    if pids == "max" or int(pids) > PID_MAX:
        return False
    quota, period = cpu.split()
    return quota != "max" and int(quota) <= int(period)


def cgroup_root_for_pid(pid: int) -> Path | None:
    try:
        relative = next(
            line.split("::", 1)[1]
            for line in Path(f"/proc/{pid}/cgroup")
            .read_text(encoding="utf-8")
            .splitlines()
            if line.startswith("0::")
        )
    except (OSError, StopIteration, ValueError):
        return None
    return Path("/sys/fs/cgroup") / relative.lstrip("/")


def cgroup_is_qualified() -> bool:
    """Return whether this process itself is inside the fixed native scope."""

    root = _self_cgroup()
    return root is not None and cgroup_path_is_qualified(root)


def _self_cgroup() -> Path | None:
    try:
        relative = next(
            line.split("::", 1)[1]
            for line in Path("/proc/self/cgroup")
            .read_text(encoding="utf-8")
            .splitlines()
            if line.startswith("0::")
        )
    except (OSError, StopIteration, ValueError):
        return None
    return Path("/sys/fs/cgroup") / relative.lstrip("/")


def _children(pid: int) -> list[int]:
    path = Path(f"/proc/{pid}/task/{pid}/children")
    try:
        text = path.read_text(encoding="utf-8")
    except OSError:
        return []
    return [int(item) for item in text.split() if item.isdigit()]


def _descendants(pid: int) -> list[int]:
    found: list[int] = []
    stack = [pid]
    while stack:
        current = stack.pop()
        children = _children(current)
        found.extend(children)
        stack.extend(children)
    return found


def qualified_native_pid(unit: str, ancestor_pid: int, *, timeout: float = 2.0) -> int:
    """Return the native PID whose own cgroup is this unit's qualified scope."""

    deadline = time.monotonic() + timeout
    while True:
        candidates = [ancestor_pid, *_descendants(ancestor_pid)]
        for pid in reversed(candidates):
            root = cgroup_root_for_pid(pid)
            if root is None or unit not in str(root):
                continue
            if cgroup_path_is_qualified(root):
                return pid
        if time.monotonic() >= deadline:
            break
        time.sleep(0.02)
    raise OSError("native cgroup scope is not qualified")


def stop_scope(unit: str | None) -> None:
    if not unit:
        return
    try:
        subprocess.run(
            ["systemctl", "--user", "stop", unit],
            check=False,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=3,
        )
    except (OSError, subprocess.TimeoutExpired):
        return
