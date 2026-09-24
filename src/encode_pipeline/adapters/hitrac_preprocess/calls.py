"""Transparent, fixed-tracPre2 child invocation receipts (not a public runner)."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shlex
import signal
import subprocess
import sys
import time
import uuid
from dataclasses import dataclass


@dataclass(frozen=True)
class CallFailure(Exception):
    code: str
    sample: str | None = None
    stage: str | None = None


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def strict_json(data: str):
    def unique(pairs):
        obj = {}
        for key, value in pairs:
            if key in obj:
                raise ValueError("duplicate key")
            obj[key] = value
        return obj

    return json.loads(
        data,
        object_pairs_hook=unique,
        parse_constant=lambda _: (_ for _ in ()).throw(ValueError("constant")),
    )


def write_exclusive(path: Path, data: dict) -> None:
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
    with os.fdopen(fd, "w") as stream:
        json.dump(data, stream, sort_keys=True, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())


def plan_calls(
    output: Path, reference: Path, samples: list[str], threads: int, mapq: int
) -> list[dict]:
    plan = []

    def add(tool, stage, argv, sample=None):
        plan.append(
            dict(tool=tool, stage=stage, sample=sample, argv=list(map(str, argv)))
        )

    all_files, unique_files = [], []
    for token in samples:
        prefix = output / token / token
        sam, bam = str(prefix) + ".sam", str(prefix) + ".bam"
        all_bed, unique = str(prefix) + "_all.bedpe", str(prefix) + "_unique.bedpe"
        add(
            "bowtie2",
            "align",
            [
                "-p",
                threads,
                "-q",
                "--end-to-end",
                "--very-sensitive",
                "-x",
                reference,
                "-1",
                str(prefix) + "_R1.fastq.gz",
                "-2",
                str(prefix) + "_R2.fastq.gz",
                "-S",
                sam,
            ],
            token,
        )
        add(
            "samtools",
            "view",
            ["view", "-b", "-F", 4, "-@", 2, "-q", mapq, "-o", bam, sam],
            token,
        )
        add(
            "samtools",
            "sort",
            ["sort", "-n", "-@", 2, bam, "-T", prefix, "-o", bam],
            token,
        )
        add("rm", "remove_sam", [sam], token)
        add("bamToBed", "bedpe", ["-bedpe", "-i", bam], token)
        add("gzip", "compress", [all_bed, unique], token)
        all_files.append(all_bed + ".gz")
        unique_files.append(unique + ".gz")
    for stage, files, name in [
        ("qc_all", all_files, "allBedpeQc"),
        ("qc_noBg", unique_files, "uniNonBgBedpeQc"),
    ]:
        add(
            "cLoops2",
            stage,
            [
                "qc",
                "-f",
                ",".join(files),
                "-o",
                output / name,
                "-p",
                min(len(samples), threads),
            ],
        )
    add(
        "rm",
        "remove_qc",
        [output / "allBedpeQc_bedpeQc.txt", output / "uniNonBgBedpeQc_bedpeQc.txt"],
    )
    return plan


def equivalent(actual: list[str], entry: dict) -> bool:
    expected = entry["argv"]
    if entry["stage"] in ("qc_all", "qc_noBg"):
        # glob/sample discovery order is not scientific ordering; do not rewrite argv.
        return (
            len(actual) == len(expected)
            and actual[:2] == expected[:2]
            and sorted(actual[2].split(",")) == sorted(expected[2].split(","))
            and actual[3:] == expected[3:]
        )
    if entry["stage"] == "remove_qc":
        return sorted(actual) == sorted(expected)
    return actual == expected


def classify(tool: str, argv: list[str], plan: list[dict]) -> tuple[int, dict]:
    matched = [
        (i, entry)
        for i, entry in enumerate(plan)
        if entry["tool"] == tool and equivalent(argv, entry)
    ]
    if len(matched) == 1:
        return matched[0]
    # Original isTool uses find_executable, not subprocess. No tool probes are
    # expected in the scientific process; separately labelled probes cannot
    # substitute any planned invocation.
    if argv in (["--version"], ["--help"]):
        return -1, dict(tool=tool, stage="probe", sample=None, argv=argv)
    return -2, dict(tool=tool, stage="unknown", sample=None, argv=argv)


def prepare_shims(
    attempt: Path, tools: dict[str, Path], plan: list[dict], identity: str
) -> tuple[Path, Path]:
    private = attempt / "private"
    receipts = private / "calls"
    receipts.mkdir(mode=0o700)
    shim_dir = private / "bin"
    shim_dir.mkdir(mode=0o700)
    config = private / "call-plan.json"
    doc = dict(
        attempt=str(attempt),
        cwd=str(attempt),
        identity=identity,
        receipts=str(receipts),
        plan=plan,
        tools={
            name: {"path": str(path), "sha256": digest(path)}
            for name, path in tools.items()
        },
    )
    write_exclusive(config, doc)
    config_sha = digest(config)
    implementation = Path(__file__).absolute()
    for tool in tools:
        # A shell shebang executes Python with real argv and no shebang splitting.
        args = [
            sys.executable,
            "-I",
            "-B",
            str(implementation),
            str(config),
            config_sha,
            tool,
        ]
        target = shim_dir / tool
        target.write_text("#!/bin/sh\nexec " + shlex.join(args) + ' "$@"\n')
        target.chmod(0o700)
    return shim_dir, config


def verify_calls(config: Path) -> list[dict]:
    try:
        doc = strict_json(config.read_text())
        expected_sha = digest(config)
        directory = Path(doc["receipts"])
        records = {}
        for path in directory.iterdir():
            if path.is_symlink() or not path.is_file():
                raise ValueError("unexpected file")
            parts = path.name.split(".")
            if (
                len(parts) != 3
                or parts[1] not in ("start", "end")
                or parts[2] != "json"
            ):
                raise ValueError("unexpected receipt")
            if len(parts[0]) != 32 or any(
                x not in "0123456789abcdef" for x in parts[0]
            ):
                raise ValueError("bad id")
            records.setdefault(parts[0], {})[parts[1]] = strict_json(path.read_text())
        counts = [0] * len(doc["plan"])
        complete = []
        for invocation, pair in records.items():
            if "start" not in pair:
                raise ValueError("missing start")
            begin = pair["start"]
            if (
                not isinstance(begin, dict)
                or set(begin)
                != {
                    "invocation",
                    "attempt",
                    "identity",
                    "config_sha256",
                    "tool",
                    "tool_sha256",
                    "argv",
                    "cwd",
                    "pid",
                    "time_ns",
                    "stage",
                    "sample",
                }
                or type(begin["pid"]) is not int
                or begin["pid"] <= 0
                or type(begin["time_ns"]) is not int
                or begin["time_ns"] <= 0
            ):
                raise ValueError("start schema")
            tool, argv = begin["tool"], begin["argv"]
            if not isinstance(argv, list) or not all(isinstance(a, str) for a in argv):
                raise ValueError("argv")
            idx, entry = classify(tool, argv, doc["plan"])
            if (
                begin["invocation"] != invocation
                or begin["config_sha256"] != expected_sha
                or begin["attempt"] != doc["attempt"]
                or begin["cwd"] != doc["cwd"]
                or begin["identity"] != doc["identity"]
                or tool not in doc["tools"]
                or begin["tool_sha256"] != doc["tools"][tool]["sha256"]
                or begin["stage"] != entry["stage"]
                or begin["sample"] != entry["sample"]
            ):
                raise ValueError("receipt binding")
            if idx == -2:
                raise CallFailure("call_unrecognized")
            if "end" not in pair:
                raise CallFailure("call_unfinished", entry["sample"], entry["stage"])
            end = pair["end"]
            if (
                not isinstance(end, dict)
                or set(end) != {"invocation", "start_sha256", "returncode", "time_ns"}
                or type(end["time_ns"]) is not int
                or end["time_ns"] <= 0
            ):
                raise ValueError("end schema")
            if (
                end["invocation"] != invocation
                or end["start_sha256"] != digest(directory / f"{invocation}.start.json")
                or type(end["returncode"]) is not int
            ):
                raise ValueError("end binding")
            if end["returncode"] != 0:
                raise CallFailure("child_failed", entry["sample"], entry["stage"])
            if idx >= 0:
                counts[idx] += 1
            complete.append(begin | {"returncode": end["returncode"]})
        for count, entry in zip(counts, doc["plan"]):
            if count != 1:
                raise CallFailure(
                    "call_missing" if count == 0 else "call_duplicate",
                    entry["sample"],
                    entry["stage"],
                )
        return complete
    except CallFailure:
        raise
    except (OSError, ValueError, KeyError, TypeError, IndexError):
        raise CallFailure("call_records_invalid") from None


def delegate(config: Path, expected_sha: str, tool: str, argv: list[str]) -> int:
    # No diagnostic text is added to scientific stdout/stderr, even on failure.
    # Missing/unfinished receipts are fail-closed in the parent.
    try:
        if digest(config) != expected_sha:
            return 125
        doc = strict_json(config.read_text())
        info = doc["tools"][tool]
        invocation = uuid.uuid4().hex
        directory = Path(doc["receipts"])
        _, entry = classify(tool, argv, doc["plan"])
        begin = dict(
            invocation=invocation,
            attempt=doc["attempt"],
            identity=doc["identity"],
            config_sha256=expected_sha,
            tool=tool,
            tool_sha256=info["sha256"],
            argv=argv,
            cwd=os.getcwd(),
            pid=os.getpid(),
            time_ns=time.time_ns(),
            stage=entry["stage"],
            sample=entry["sample"],
        )
        start = directory / f"{invocation}.start.json"
        write_exclusive(start, begin)
        if entry["stage"] == "unknown" or os.getcwd() != doc["cwd"]:
            code = 126
        elif digest(Path(info["path"])) != info["sha256"]:
            code = 125
        else:
            # Absolute alias path retained: bamToBed must not turn into bare bedtools.
            code = subprocess.call([info["path"], *argv])
        write_exclusive(
            directory / f"{invocation}.end.json",
            dict(
                invocation=invocation,
                start_sha256=digest(start),
                returncode=code,
                time_ns=time.time_ns(),
            ),
        )
        if code < 0:
            if -code not in (signal.SIGKILL, signal.SIGSTOP):
                signal.signal(-code, signal.SIG_DFL)
            os.kill(os.getpid(), -code)
        return code
    except (OSError, ValueError, KeyError, TypeError):
        return 125


if __name__ == "__main__":
    raise SystemExit(
        delegate(Path(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4:])
    )
