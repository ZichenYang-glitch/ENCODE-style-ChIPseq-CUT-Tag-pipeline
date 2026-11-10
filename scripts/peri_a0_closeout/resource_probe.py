"""Exercise only this transient scope's aggregate controllers, with finite loads."""

import errno
import json
import os
from pathlib import Path
import sys
import time
from resource_guard import inspect_limits


def main():
    destination = Path(sys.argv[1])
    assert not destination.exists()
    report = {"limits": inspect_limits(), "passed": False}
    root = Path("/sys/fs/cgroup") / report["limits"]["path"].lstrip("/")
    children = []
    try:
        report["memory"] = {
            "not_repeated": True,
            "evidence": "resource-oom-systemd.txt: first scope stopped with Result=oom-kill",
        }
        before = (root / "pids.events").read_text()
        rejected = False
        for _ in range(32):
            try:
                pid = os.fork()
            except OSError as error:
                assert error.errno == errno.EAGAIN
                rejected = True
                break
            if pid == 0:
                time.sleep(3)
                os._exit(0)
            children.append(pid)
        report["pids"] = {
            "rejected": rejected,
            "children_started": len(children),
            "before": before,
            "after": (root / "pids.events").read_text(),
        }
        assert rejected
        for pid in children:
            os.kill(pid, 9)
        for pid in children:
            os.waitpid(pid, 0)
        children = []
        before = (root / "cpu.stat").read_text()
        stop = time.monotonic() + 0.7
        while time.monotonic() < stop:
            pass
        after = (root / "cpu.stat").read_text()
        report["cpu"] = {"before": before, "after": after}

        def stat(text):
            return dict(line.split() for line in text.splitlines())

        assert int(stat(after)["nr_throttled"]) > int(stat(before)["nr_throttled"])
        report["memory_peak_bytes"] = (root / "memory.peak").read_text().strip()
        report["passed"] = True
    finally:
        for pid in children:
            try:
                os.kill(pid, 9)
            except ProcessLookupError:
                pass
            os.waitpid(pid, 0)
        report["children_reaped"] = True
        destination.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
