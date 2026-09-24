#!/usr/bin/env python3
"""Verify an existing Hi-TrAC runtime and record an externally built reference.

This does not download, build, execute science, or approve the new reference
record. Its printed reference digest must be independently reviewed before use.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import runpy
import sys

if not sys.flags.isolated or not sys.flags.no_site:
    raise SystemExit("invoke with Python -I -S")
sys.dont_write_bytecode = True
root = Path(__file__).resolve().parents[1]
provenance = runpy.run_path(str(root / "scripts/source_provenance.py"))
provenance["bootstrap_checkout"](root)

from encode_pipeline.adapters.hitrac_preprocess.admission import (  # noqa: E402
    INDEX_SUFFIXES,
    AdmissionError,
    _fasta_contigs,
    load_reference_binding,
    load_runtime_binding,
    sha256_file,
)


def write_private(path: Path, value: dict) -> None:
    with path.open("x", encoding="utf-8") as handle:
        path.chmod(0o600)
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime-prefix", type=Path, required=True)
    parser.add_argument("--tracpre2", type=Path, required=True)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--index-prefix", type=Path, required=True)
    parser.add_argument("--index-build-record", type=Path, required=True)
    parser.add_argument("--index-build-record-sha256", required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    try:
        if sha256_file(args.index_build_record) != args.index_build_record_sha256:
            raise AdmissionError("reference_build_record_changed")
        record = json.loads(args.index_build_record.read_text(encoding="utf-8"))
        argv = record.get("argv")
        # The record is a separately trusted execution record, not an assertion
        # that this preparation command has run bowtie2-build itself.
        expected_tool = args.runtime_prefix / "bin/bowtie2-build"
        if (
            record.get("exit") != 0
            or isinstance(record.get("exit"), bool)
            or record.get("status") != "completed"
            or not isinstance(argv, list)
            or len(argv) < 3
            or argv[0] != str(expected_tool)
            or argv[-2:] != [str(args.reference_fasta), str(args.index_prefix)]
        ):
            raise AdmissionError("reference_build_record_invalid")
        output = args.output_directory
        output.mkdir(mode=0o700)
        lock = root / "config/hitrac_preprocess/tools.lock.json"
        runtime_path = output / "runtime-binding.json"
        write_private(
            runtime_path,
            {
                "schema_version": "hitrac-runtime-binding-v1",
                "prefix": str(args.runtime_prefix),
                "script": str(args.tracpre2),
                "lock_sha256": sha256_file(lock),
            },
        )
        runtime = load_runtime_binding(runtime_path)
        files = {"fasta": sha256_file(args.reference_fasta)}
        for suffix in INDEX_SUFFIXES:
            files[suffix] = sha256_file(Path(f"{args.index_prefix}{suffix}"))
        reference_path = output / "reference-binding.json"
        write_private(
            reference_path,
            {
                "schema_version": "hitrac-reference-binding-v1",
                "fasta": str(args.reference_fasta),
                "prefix": str(args.index_prefix),
                "files": files,
                "contigs": _fasta_contigs(args.reference_fasta),
                "provenance": {
                    "runtime_lock_sha256": runtime.lock_sha256,
                    "build_record_sha256": args.index_build_record_sha256,
                    "build_argv": argv,
                    "build_exit": 0,
                    "preparation_executed_index_build": False,
                },
            },
        )
        reference_sha = sha256_file(reference_path)
        load_reference_binding(reference_path, reference_sha)
        print(
            json.dumps(
                {
                    "state": "verified_pending_reference_approval",
                    "runtime_binding_sha256": sha256_file(runtime_path),
                    "reference_binding_sha256": reference_sha,
                    "runtime_lock_sha256": runtime.lock_sha256,
                },
                sort_keys=True,
            )
        )
        return 0
    except AdmissionError as exc:
        print(json.dumps({"state": "rejected", "reason_code": exc.reason_code}))
        return 1
    except FileExistsError:
        print(json.dumps({"state": "rejected", "reason_code": "binding_output_exists"}))
        return 1
    except (OSError, ValueError, TypeError, AttributeError):
        print(
            json.dumps(
                {"state": "rejected", "reason_code": "binding_preparation_invalid"}
            )
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
