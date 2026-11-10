"""Host-independent adversarial probe; reads only synthetic test coordinates."""

import json
import os
from pathlib import Path
import socket


def cannot_read(path):
    try:
        Path(path).read_bytes()
    except (FileNotFoundError, PermissionError, IsADirectoryError):
        return True
    return False


def cannot_connect(address, family):
    with socket.socket(family, socket.SOCK_STREAM) as client:
        client.settimeout(1)
        try:
            client.connect(address)
        except OSError:
            return True
    return False


coordinates = json.loads(Path("/state/probe.json").read_text())
checks = {
    name: cannot_read(path) for name, path in coordinates["forbidden_files"].items()
}
checks["host_tcp_denied"] = cannot_connect(
    ("127.0.0.1", coordinates["tcp_port"]), socket.AF_INET
)
checks["platform_socket_denied"] = cannot_connect(
    coordinates["socket_path"], socket.AF_UNIX
)
checks["abstract_socket_denied"] = cannot_connect(
    "\0" + coordinates["abstract"], socket.AF_UNIX
)
checks["secret_environment_absent"] = "A0_BRIDGE_SECRET" not in os.environ
checks["host_process_hidden"] = not Path(f"/proc/{coordinates['host_pid']}").exists()
checks["private_state_writable"] = True
Path("/state/probe-write").write_text("synthetic")
fds = {}
for name in os.listdir("/proc/self/fd"):
    try:
        fds[name] = os.readlink("/proc/self/fd/" + name)
    except FileNotFoundError:
        pass
checks["no_unexpected_fds"] = set(fds) <= {"0", "1", "2"}
status = Path("/proc/self/status").read_text()
checks["no_new_privs"] = "NoNewPrivs:\t1" in status
checks["no_effective_caps"] = "CapEff:\t0000000000000000" in status
checks["no_host_mounts"] = not Path("/home").exists() and not Path("/run").exists()
print(
    json.dumps(
        {
            "checks": checks,
            "fds": fds,
            "uid": os.getuid(),
            "gid": os.getgid(),
            "namespaces": {
                n: os.readlink("/proc/self/ns/" + n)
                for n in ["user", "pid", "net", "mnt"]
            },
            "environment_keys": sorted(os.environ),
        }
    )
)
raise SystemExit(0 if all(checks.values()) else 1)
