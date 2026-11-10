"""Run only catalog/policy library tests against the recorded final candidate."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    build = json.loads(args.manifest.read_text())
    patch = subprocess.check_output(["git", "diff", "HEAD"], cwd=args.source)
    assert hashlib.sha256(patch).hexdigest() == build["patch_sha256"]
    report = {
        "patch_sha256": build["patch_sha256"],
        "started_at_unix": time.time(),
        "commands": [],
        "passed": False,
    }
    try:
        for pattern in [
            "process_policy_tests",
            "tool_filter_tests",
            "tool_catalog_test",
        ]:
            command = [
                "./scripts/cargo-rmcp-patched.sh",
                "test",
                "--locked",
                "--release",
                "-p",
                "peri-agent",
                "--lib",
                pattern,
            ]
            with (args.output / (pattern + ".log")).open("xb") as log:
                result = subprocess.run(
                    command,
                    cwd=args.source,
                    env=build["environment"],
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    timeout=600,
                )
            report["commands"].append({"argv": command, "exit_code": result.returncode})
            if result.returncode:
                break
        report["passed"] = len(report["commands"]) == 3 and all(
            x["exit_code"] == 0 for x in report["commands"]
        )
    finally:
        report["finished_at_unix"] = time.time()
        (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
