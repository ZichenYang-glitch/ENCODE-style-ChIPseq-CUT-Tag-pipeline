"""Bounded sampling of only the process tree started by this experiment."""

from pathlib import Path
import threading
import time


class Observation:
    def __init__(self, pid):
        self.pid = pid
        self.samples = []
        self.pids = set()
        self.stopping = threading.Event()
        self.thread = threading.Thread(target=self._run, daemon=True)
        self.thread.start()

    def _run(self):
        while not self.stopping.is_set():
            pending = [self.pid]
            rss = 0
            while pending:
                pid = pending.pop()
                try:
                    base = Path("/proc") / str(pid)
                    children = (
                        (base / "task" / str(pid) / "children").read_text().split()
                    )
                    pending.extend(map(int, children))
                    fields = dict(
                        line.split(":", 1)
                        for line in (base / "status").read_text().splitlines()
                    )
                    rss += int(fields.get("VmRSS", "0 kB").split()[0])
                    self.pids.add(pid)
                except (FileNotFoundError, ProcessLookupError):
                    pass
            if len(self.samples) < 2000:
                self.samples.append({"at": time.monotonic(), "tree_rss_kib": rss})
            self.stopping.wait(0.1)

    def close(self):
        self.stopping.set()
        self.thread.join(timeout=5)
        return {
            "interval_seconds": 0.1,
            "samples": len(self.samples),
            "peak_sampled_tree_rss_kib": max(
                (s["tree_rss_kib"] for s in self.samples), default=0
            ),
            "observed_pids": sorted(self.pids),
            "still_present": [
                pid for pid in self.pids if Path(f"/proc/{pid}").exists()
            ],
            "boundary": "sampled bwrap/relay/Peri descendants; excludes host stub; not an enforced memory limit",
        }
