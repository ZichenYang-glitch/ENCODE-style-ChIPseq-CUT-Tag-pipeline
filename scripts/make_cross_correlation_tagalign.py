#!/usr/bin/env python3
"""Build the cross-correlation tagAlign and its policy/window metadata.

Input is the canonical analysis BAM.  PE uses only the R1 alignment of each
complete valid template; SE uses valid primary alignments.  No additional
duplicate or blacklist filtering is applied here: the analysis BAM already
carries the run's policy.  The read-length estimate follows the pinned
``run_spp.R``: median span of the first 500 tag records, rounded with R's
round-half-to-even rule.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC = _REPO_ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from encode_pipeline.qc.histone import classify_histone_target  # noqa: E402

from qc_template import (  # noqa: E402
    TemplateReader,
    cigar_ref_length,
    classify_pe_template,
)

SAMTOOLS = os.environ.get("SAMTOOLS", "samtools")

_FLAG_UNMAPPED = 0x4
_FLAG_SECONDARY = 0x100
_FLAG_QC_FAIL = 0x200
_FLAG_SUPPLEMENTARY = 0x800
_FLAG_REVERSE = 0x10
_FLAG_FIRST = 0x40


def _r_round(value: float) -> int:
    """Match R's round() for the finite positive values used here.

    Python's round() uses round-half-to-even, as does R's round() for these
    values.  Keep the helper explicit so the method is visible in tests.
    """
    return int(round(value))


def _resolve_target(assay: str, target: str, requested: str):
    if requested not in ("auto", "tf", "histone"):
        raise ValueError(f"unsupported cross-correlation target type: {requested!r}")
    inferred = classify_histone_target(target)["inferred_histone_class"]
    if assay != "chipseq":
        return "unknown", "diagnostic_non_chip", inferred
    if requested == "auto":
        if inferred != "unknown":
            return "histone", "auto_histone", inferred
        return "unknown", "auto_unknown", inferred
    return requested, "explicit", inferred


def _window(read_length: int | None, target_type: str):
    lower = -500
    base = (read_length or 0) + 10
    if target_type == "tf":
        upper = max(base, 50)
        policy = "tf"
    elif target_type == "histone":
        upper = max(base, 100)
        policy = "histone"
    else:
        upper = max(base, 100)
        policy = "diagnostic"
    return lower, upper, policy


def _pe_tags(reader: TemplateReader):
    valid = 0
    invalid = 0
    for template in reader:
        classification = classify_pe_template(template)
        if not classification.valid:
            invalid += 1
            continue
        valid += 1
        first = next(
            (record for record in template.records if record.flag & _FLAG_FIRST),
            None,
        )
        if first is None:
            invalid += 1
            valid -= 1
            continue
        if (
            any(record.flag & _FLAG_QC_FAIL for record in template.records)
            or first.pos < 1
            or cigar_ref_length(first.cigar) <= 0
        ):
            invalid += 1
            valid -= 1
            continue
        start = first.pos - 1
        end = start + cigar_ref_length(first.cigar)
        strand = "-" if first.flag & _FLAG_REVERSE else "+"
        yield first.rname, start, end, strand
    reader.valid_templates = valid
    reader.invalid_templates = invalid


def _se_tags(bam: str, samtools: str):
    process = subprocess.Popen(
        [samtools, "view", bam],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    try:
        assert process.stdout is not None
        for line in process.stdout:
            if not line or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError("BAM alignment record is malformed")
            try:
                flag = int(fields[1])
                pos = int(fields[3])
            except ValueError as exc:
                raise ValueError("BAM alignment record is malformed") from exc
            if flag & (
                _FLAG_UNMAPPED | _FLAG_SECONDARY | _FLAG_SUPPLEMENTARY | _FLAG_QC_FAIL
            ):
                continue
            start = pos - 1
            span = cigar_ref_length(fields[5])
            if pos < 1 or fields[2] == "*" or span <= 0:
                continue
            end = start + span
            strand = "-" if flag & _FLAG_REVERSE else "+"
            yield fields[2], start, end, strand
    finally:
        if process.stdout is not None:
            process.stdout.close()
        stderr = process.stderr.read() if process.stderr is not None else ""
        returncode = process.wait()
        if returncode != 0:
            raise RuntimeError(f"samtools view failed on {bam}: {stderr.strip()}")


def _read_stats(tagalign: Path):
    spans = []
    distribution = {}
    with tagalign.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            try:
                span = int(fields[2]) - int(fields[1])
            except ValueError:
                continue
            if span <= 0:
                raise ValueError("tagAlign spans must be positive")
            distribution[str(span)] = distribution.get(str(span), 0) + 1
            if len(spans) < 500:
                spans.append(span)
    if not spans:
        return None, {}, False
    read_length = _r_round(statistics.median(spans))
    return read_length, dict(sorted(distribution.items())), len(distribution) > 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-bam", required=True)
    parser.add_argument("--layout", required=True, choices=("PE", "SE"))
    parser.add_argument("--sample", required=True)
    parser.add_argument("--assay", required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument(
        "--target-type", default="auto", choices=("auto", "tf", "histone")
    )
    parser.add_argument("--tagalign", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--samtools", default=SAMTOOLS)
    args = parser.parse_args()

    target_type, target_source, inferred_class = _resolve_target(
        args.assay, args.target, args.target_type
    )
    tagalign = Path(args.tagalign)
    tagalign.parent.mkdir(parents=True, exist_ok=True)

    valid = 0
    invalid = 0
    alignments_total = 0
    if not os.path.isfile(args.analysis_bam):
        raise FileNotFoundError("analysis BAM is unavailable")
    if args.layout == "PE":
        with TemplateReader(args.analysis_bam, samtools=args.samtools) as reader:
            with tagalign.open("w", encoding="utf-8") as handle:
                for rname, start, end, strand in _pe_tags(reader):
                    handle.write(f"{rname}\t{start}\t{end}\tN\t1000\t{strand}\n")
            valid = getattr(reader, "valid_templates", 0)
            invalid = getattr(reader, "invalid_templates", 0)
            alignments_total = reader.alignments_total
        status = "ready" if tagalign.stat().st_size else "no_data"
    else:
        with tagalign.open("w", encoding="utf-8") as handle:
            for rname, start, end, strand in _se_tags(args.analysis_bam, args.samtools):
                handle.write(f"{rname}\t{start}\t{end}\tN\t1000\t{strand}\n")
                valid += 1
        status = "ready" if valid else "no_data"

    records = sum(1 for _ in tagalign.open(encoding="utf-8"))
    if status == "ready" and records > 1:
        env = os.environ.copy()
        env["LC_ALL"] = "C"
        subprocess.run(
            ["sort", "-k1,1", "-k2,2n", str(tagalign), "-o", str(tagalign)],
            check=True,
            env=env,
        )
    read_length, distribution, variable = _read_stats(tagalign)
    lower, upper, policy = _window(read_length, target_type)
    if status != "ready":
        policy = "no_data"
    metadata = {
        "sample": args.sample,
        "assay": args.assay,
        "target": args.target,
        "layout": args.layout,
        "input_stage": "analysis_r1" if args.layout == "PE" else "analysis_se",
        "status": status,
        "tagalign_records": records,
        "valid_templates": valid if args.layout == "PE" else None,
        "invalid_templates": invalid if args.layout == "PE" else None,
        "alignments_total": alignments_total if args.layout == "PE" else None,
        "read_length": read_length,
        "read_length_method": "median_first_500_tag_spans_r_round",
        "read_length_distribution": distribution,
        "variable_read_length": bool(variable),
        "target_type": target_type,
        "target_type_source": target_source,
        "target_inferred_histone_class": inferred_class,
        "window_lower": lower,
        "window_upper": upper,
        "window_policy": policy,
    }
    Path(args.metadata).parent.mkdir(parents=True, exist_ok=True)
    Path(args.metadata).write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
