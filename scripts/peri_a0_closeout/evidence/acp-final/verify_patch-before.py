"""Run only catalog/policy library tests against the recorded final candidate."""

import argparse
import hashlib
import json
import re
from pathlib import Path
import subprocess
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--pattern", action="append")
    parser.add_argument(
        "--crate",
        default="peri-agent",
        choices=["peri-agent", "peri-acp", "peri-middlewares"],
    )
    parser.add_argument("--online-locked", action="store_true")
    parser.add_argument("--test-opt-level-zero", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    build = json.loads(args.manifest.read_text())
    patch = subprocess.check_output(["git", "diff", "HEAD"], cwd=args.source)
    assert hashlib.sha256(patch).hexdigest() == build["patch_sha256"]
    environment = dict(build["environment"])
    if args.online_locked:
        environment["CARGO_NET_OFFLINE"] = "false"
    report = {
        "environment": environment,
        "patch_sha256": build["patch_sha256"],
        "started_at_unix": time.time(),
        "commands": [],
        "passed": False,
    }
    try:
        patterns = args.pattern or [
            "process_policy_tests",
            "tool_filter_tests",
            "session::tool_catalog::tests",
        ]
        for pattern in patterns:
            command = [
                "./scripts/cargo-rmcp-patched.sh",
                "test",
                "--locked",
                "--release",
                "-p",
                args.crate,
                "--lib",
                pattern,
            ]
            if args.test_opt_level_zero:
                command[1:1] = [
                    "--config",
                    f"profile.release.package.{args.crate}.opt-level=0",
                ]
            with (args.output / (pattern + ".log")).open("xb") as log:
                result = subprocess.run(
                    command,
                    cwd=args.source,
                    env=environment,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    timeout=600,
                )
            summary = re.search(
                r"test result: ok\. (\d+) passed;",
                (args.output / (pattern + ".log")).read_text(),
            )
            count = int(summary[1]) if summary else 0
            report["commands"].append(
                {"argv": command, "exit_code": result.returncode, "passed_tests": count}
            )
            if result.returncode or count == 0:
                break
        report["passed"] = len(report["commands"]) == len(patterns) and all(
            x["exit_code"] == 0 and x["passed_tests"] > 0 for x in report["commands"]
        )
    finally:
        report["finished_at_unix"] = time.time()
        (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
