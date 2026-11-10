"""Fixed bubblewrap policy for A0; never falls back to an ordinary subprocess."""

from pathlib import Path
import subprocess
import json
from configuration import process_policy


ROOT = Path(__file__).resolve().parent


def command(binary, state, model_socket, *, probe=None, policy=None, prediction="off"):
    args = [
        "/usr/bin/bwrap",
        "--unshare-all",
        "--unshare-user",
        "--disable-userns",
        "--die-with-parent",
        "--new-session",
        "--cap-drop",
        "ALL",
        "--clearenv",
        "--proc",
        "/proc",
        "--dev",
        "/dev",
        "--tmpfs",
        "/tmp",
        "--dir",
        "/work",
        "--dir",
        "/runtime",
        "--ro-bind",
        str(binary.resolve()),
        "/runtime/peri",
        "--ro-bind",
        str(ROOT / "relay.py"),
        "/runtime/relay.py",
        "--ro-bind",
        "/usr/bin/python3.12",
        "/usr/bin/python3.12",
        "--ro-bind",
        "/usr/lib/python3.12",
        "/usr/lib/python3.12",
        "--ro-bind",
        "/lib/x86_64-linux-gnu",
        "/lib/x86_64-linux-gnu",
        "--ro-bind",
        "/lib64",
        "/lib64",
        "--ro-bind",
        "/etc/ssl/certs/ca-certificates.crt",
        "/etc/ssl/certs/ca-certificates.crt",
        "--bind",
        str(state.resolve()),
        "/state",
        "--ro-bind",
        str(model_socket.resolve()),
        "/model.sock",
        "--setenv",
        "PATH",
        "/usr/bin",
        "--setenv",
        "HOME",
        "/state/home",
        "--setenv",
        "XDG_CONFIG_HOME",
        "/state/config",
        "--setenv",
        "XDG_CACHE_HOME",
        "/state/cache",
        "--setenv",
        "XDG_DATA_HOME",
        "/state/data",
        "--setenv",
        "TMPDIR",
        "/tmp",
        "--setenv",
        "LANG",
        "C.UTF-8",
        "--setenv",
        "TZ",
        "UTC",
        "--setenv",
        "PERI_MCP_BUILTIN",
        "off",
        "--setenv",
        "SSL_CERT_FILE",
        "/etc/ssl/certs/ca-certificates.crt",
        "--chdir",
        "/work",
    ]
    args.extend(
        [
            "--setenv",
            "PERI_PROCESS_TOOL_POLICY",
            json.dumps(process_policy() if policy is None else policy),
            "--setenv",
            "PERI_PREDICTION_POLICY",
            prediction,
            "--setenv",
            "TOKIO_WORKER_THREADS",
            "2",
        ]
    )
    if probe is not None:
        args.extend(["--ro-bind", str(probe.resolve()), "/runtime/probe.py"])
    args.extend(
        [
            "--remount-ro",
            "/",
            "--",
            "/usr/bin/python3.12",
            "-I",
            "-S",
            "/runtime/probe.py" if probe else "/runtime/relay.py",
        ]
    )
    return args


def spawn(binary, state, model_socket, stderr, *, policy=None, prediction="off"):
    for name in ["home", "config", "cache", "data"]:
        (state / name).mkdir(parents=True, exist_ok=True)
    return subprocess.Popen(
        [
            "/usr/bin/prlimit",
            "--nofile=128:128",
            "--core=0:0",
            "--fsize=16777216:16777216",
            "--",
            *command(binary, state, model_socket, policy=policy, prediction=prediction),
        ],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=stderr,
        close_fds=True,
        env={},
        start_new_session=True,
    )
