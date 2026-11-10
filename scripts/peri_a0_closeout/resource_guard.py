"""Require effective aggregate cgroup v2 limits, not merely RLIMIT or RSS samples."""

from pathlib import Path


def inspect_limits():
    relative = next(
        line.split("::", 1)[1]
        for line in Path("/proc/self/cgroup").read_text().splitlines()
        if line.startswith("0::")
    )
    root = Path("/sys/fs/cgroup") / relative.lstrip("/")
    values = {
        name: (root / name).read_text().strip()
        for name in ["memory.max", "memory.swap.max", "pids.max", "cpu.max"]
    }
    assert values["memory.max"] != "max" and int(values["memory.max"]) <= 536870912, (
        values
    )
    assert values["memory.swap.max"] == "0", values
    assert values["pids.max"] != "max" and int(values["pids.max"]) <= 128, values
    quota, period = values["cpu.max"].split()
    assert quota != "max" and int(quota) <= int(period), values
    return {"path": relative, "effective_limits": values}
