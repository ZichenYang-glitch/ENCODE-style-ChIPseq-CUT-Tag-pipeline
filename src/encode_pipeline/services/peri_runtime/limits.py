"""Host cgroup admission for one qualified Peri runtime."""

from __future__ import annotations

from pathlib import Path

from encode_pipeline.services.peri_runtime.identity import MEMORY_MAX_BYTES, PID_MAX


def cgroup_is_qualified() -> bool:
    """Return whether this process is inside the fixed A1 runtime scope."""

    try:
        relative = next(
            line.split("::", 1)[1]
            for line in Path("/proc/self/cgroup")
            .read_text(encoding="utf-8")
            .splitlines()
            if line.startswith("0::")
        )
        root = Path("/sys/fs/cgroup") / relative.lstrip("/")
        memory = (root / "memory.max").read_text(encoding="utf-8").strip()
        swap = (root / "memory.swap.max").read_text(encoding="utf-8").strip()
        pids = (root / "pids.max").read_text(encoding="utf-8").strip()
        cpu = (root / "cpu.max").read_text(encoding="utf-8").strip()
    except (OSError, StopIteration, ValueError):
        return False
    if memory == "max" or int(memory) > MEMORY_MAX_BYTES or swap != "0":
        return False
    if pids == "max" or int(pids) > PID_MAX:
        return False
    quota, period = cpu.split()
    return quota != "max" and int(quota) <= int(period)
