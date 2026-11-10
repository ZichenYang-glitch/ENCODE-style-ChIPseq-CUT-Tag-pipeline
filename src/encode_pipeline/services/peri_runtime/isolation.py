"""Bubblewrap launch for one Peri process. There is no plain-subprocess fallback."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess

from encode_pipeline.services.peri_runtime.policy import process_policy


BWRAP = Path("/usr/bin/bwrap")
PRLIMIT = Path("/usr/bin/prlimit")
PYTHON = Path("/usr/bin/python3.12")
PYTHON_LIB = Path("/usr/lib/python3.12")
LIB_DIR = Path("/lib/x86_64-linux-gnu")
LIB64 = Path("/lib64")
CA_BUNDLE = Path("/etc/ssl/certs/ca-certificates.crt")
RELAY_SCRIPT = Path(__file__).with_name("relay.py")


def isolation_paths_exist() -> bool:
    return all(
        path.exists()
        for path in (
            BWRAP,
            PRLIMIT,
            PYTHON,
            PYTHON_LIB,
            LIB_DIR,
            LIB64,
            CA_BUNDLE,
            RELAY_SCRIPT,
        )
    )


def command(
    binary: Path,
    state: Path,
    model_socket: Path,
    *,
    policy: dict[str, object] | None = None,
    prediction: str = "off",
) -> list[str]:
    args = [
        str(BWRAP),
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
        str(RELAY_SCRIPT.resolve()),
        "/runtime/relay.py",
        "--ro-bind",
        str(PYTHON),
        "/usr/bin/python3.12",
        "--ro-bind",
        str(PYTHON_LIB),
        "/usr/lib/python3.12",
        "--ro-bind",
        str(LIB_DIR),
        "/lib/x86_64-linux-gnu",
        "--ro-bind",
        str(LIB64),
        "/lib64",
        "--ro-bind",
        str(CA_BUNDLE),
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
        "--setenv",
        "PERI_PROCESS_TOOL_POLICY",
        json.dumps(process_policy() if policy is None else policy),
        "--setenv",
        "PERI_PREDICTION_POLICY",
        prediction,
        "--setenv",
        "TOKIO_WORKER_THREADS",
        "2",
        "--remount-ro",
        "/",
        "--",
        "/usr/bin/python3.12",
        "-I",
        "-S",
        "/runtime/relay.py",
    ]
    return args


def spawn(
    binary: Path,
    state: Path,
    model_socket: Path,
    stderr,
    *,
    policy: dict[str, object] | None = None,
    prediction: str = "off",
) -> subprocess.Popen[bytes]:
    for name in ("home", "config", "cache", "data"):
        (state / name).mkdir(parents=True, exist_ok=True)
    return subprocess.Popen(
        [
            str(PRLIMIT),
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
