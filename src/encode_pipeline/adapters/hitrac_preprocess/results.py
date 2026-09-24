"""Verified Hi-TrAC result projection and numeric QC; no replacement science."""

from __future__ import annotations

import csv
from dataclasses import replace
from decimal import Decimal, InvalidOperation, ROUND_HALF_EVEN
import io
import math
from pathlib import Path
import re
import shutil
import stat
import subprocess
import uuid

from encode_pipeline.platform.adapters import (
    ExtractedArtifactCandidate,
    ExtractedQcMetricCandidate,
    WorkflowAvailability,
    WorkflowCapabilities,
)
from encode_pipeline.platform.results import Result

from .adapter import HiTracPreprocessAdapter
from .admission import _regular, sha256_file
from .calls import CallFailure, digest, plan_calls, strict_json, verify_calls
from .execution import ATTEMPT, REQUEST, encoded, failure
from .outputs import iter_bedpe, summary_columns, verify_outputs
from .pairs import verify_pairs
from .qualification import implementation_identity

SUMMARY_TYPE = "hitrac_summary"
ALL_TYPE = "hitrac_bedpe_all"
NO_BG_TYPE = "hitrac_bedpe_no_bg"
METRIC_KEYS = (
    "raw_pairs",
    "trimmed_pairs",
    "mapping_ratio",
    "all_pets",
    "all_redundancy",
    "all_cis_ratio",
    "all_close_ratio",
    "all_middle_ratio",
    "all_distal_ratio",
    "no_bg_pets",
    "yield",
    "no_bg_cis_ratio",
    "no_bg_close_ratio",
    "no_bg_middle_ratio",
    "no_bg_distal_ratio",
)
_COUNTS = frozenset({0, 1, 3, 9})


def _require(value):
    if not value:
        raise ValueError("result_contract_invalid")


def _directory(path):
    _require(path.is_absolute() and path.resolve() == path)
    current = Path(path.anchor)
    for part in path.parts[1:]:
        current /= part
        _require(stat.S_ISDIR(current.lstat().st_mode))
    return path


def _json(path):
    value = strict_json(_regular(path).read_text(encoding="utf-8"))
    _require(isinstance(value, dict))
    return value


def parse_summary(content, samples, mapq):
    """Decode the upstream's original 16 physical columns without changing bytes.

    Full scientific correspondence is checked by verify_outputs before projection.
    Here the QC consumer additionally checks its bounded source and frozen IDs.
    """
    rows = csv.reader(io.StringIO(content.decode("ascii")), delimiter="\t")
    _require(next(rows, None) == [""] + summary_columns(mapq))
    result = {}
    for row in rows:
        _require(len(row) == 16 and row[0] in samples and row[0] not in result)
        _require(
            all(
                re.fullmatch(
                    r"-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?", text
                )
                for text in row[1:]
            )
        )
        values = [Decimal(text) for text in row[1:]]
        _require(all(value.is_finite() for value in values))
        for index, value in enumerate(values):
            if index in _COUNTS:
                _require(value >= 0 and value == value.to_integral_value())
            elif index == 10:
                _require(value >= 0)
            else:
                _require(0 <= value <= 1)
        _require(values[0] == samples[row[0]]["raw_pairs"])
        _require(0 < values[1] <= values[0] and values[3] > 0 and values[9] > 0)
        _require(values[9] <= values[3] and values[5] > 0 and values[11] > 0)
        _require(
            math.isclose(
                float(values[10]),
                float(values[9] / values[0]),
                abs_tol=1e-12,
                rel_tol=0,
            )
        )
        for start in (6, 12):
            _require(abs(sum(values[start : start + 3]) - 1) <= Decimal("1e-12"))
        result[row[0]] = values
    _require(set(result) == set(samples))
    return result


def _verified_outputs(binding, inputs, workspace):
    runtime, reference, normalized, samples = binding.verify()
    _require(encoded(normalized.to_dict()) == encoded(inputs.to_dict()))
    workspace = _directory(Path(workspace))
    attempt = _directory(workspace / ATTEMPT)
    private = _directory(attempt / "private")
    output = _directory(attempt / "output")
    mapq, threads = normalized.options["mapq"], normalized.options["threads"]
    expected_plan = {
        "schema_version": "hitrac-private-plan-v1",
        "workspace": str(workspace),
        "runtime_binding": str(binding.runtime.binding),
        "runtime_sha256": binding.runtime.binding_sha256,
        "reference_binding": str(binding.reference),
        "reference_sha256": binding.reference_sha256,
        "inputs": normalized.to_dict(),
        "samples": samples,
        "identity": binding.identity(),
        "timeout": binding.runtime.timeout,
    }
    _require(_json(workspace / REQUEST) == expected_plan)
    complete = _json(attempt / "complete.json")
    identity = implementation_identity() | {
        "runtime_lock_sha256": runtime.lock_sha256,
        "runtime_binding_sha256": runtime.binding_sha256,
        "reference_binding_sha256": reference.binding_sha256,
    }
    _require(
        set(complete) == {"schema_version", "identity", "results", "mapq", "threads"}
    )
    _require(complete["schema_version"] == "hitrac-qualification-complete-v1")
    _require(complete["identity"] == identity == _json(private / "identity.json"))
    _require(type(complete["mapq"]) is int and type(complete["threads"]) is int)
    _require(complete["mapq"] == mapq and complete["threads"] == threads)
    expected_inputs = {
        "runtime_binding_sha256": runtime.binding_sha256,
        "runtime_lock_sha256": runtime.lock_sha256,
        "reference_binding_sha256": reference.binding_sha256,
        "samples": samples,
    }
    _require(_json(private / "inputs.json") == expected_inputs)
    _require(
        _json(private / "outcome.json")
        == {
            "status": "complete",
            "reason_code": None,
            "qualification_sha256": identity["sha256"],
            "samples": list(samples),
        }
    )
    staged_reference = attempt / "reference/genome"
    for suffix, expected in reference.files.items():
        path = (
            attempt / "reference/reference.fa"
            if suffix == "fasta"
            else Path(str(staged_reference) + suffix)
        )
        _require(sha256_file(_regular(path)) == expected)
    for token, sample in samples.items():
        for mate in ("r1", "r2"):
            path = attempt / "input" / (token + "_" + mate.upper() + ".fastq.gz")
            _require(sha256_file(_regular(path)) == sample["sha256"][mate])
    tools = {
        name: runtime.tools[name]
        for name in ("bowtie2", "samtools", "bamToBed", "gzip", "cLoops2", "rm")
    }
    expected_calls = {
        "attempt": str(attempt),
        "cwd": str(attempt),
        "identity": identity["sha256"],
        "receipts": str(private / "calls"),
        "plan": plan_calls(output, staged_reference, list(samples), threads, mapq),
        "tools": {
            name: {"path": str(path), "sha256": digest(path)}
            for name, path in tools.items()
        },
    }
    _require(_json(private / "call-plan.json") == expected_calls)
    _directory(private / "calls")
    verify_calls(private / "call-plan.json")
    command = _json(private / "command.json")
    _require(command["cwd"] == str(attempt))
    _require(
        command["argv"]
        == [
            str(runtime.python),
            "-I",
            "-B",
            str(runtime.script),
            "-fqd",
            str(attempt / "input"),
            "-o",
            str(output),
            "-ref",
            str(staged_reference),
            "-n",
            "1",
            "-p",
            str(threads),
            "-mapq",
            str(mapq),
        ]
    )
    execution = _json(private / "execution.json")
    _require(execution["returncode"] == 0 and execution["reason"] is None)
    scratch = private / ("result-check-" + uuid.uuid4().hex)
    scratch.mkdir(mode=0o700)
    results = verify_outputs(
        output,
        {token: sample["raw_pairs"] for token, sample in samples.items()},
        reference.contigs,
        mapq,
        scratch,
    )
    _require(encoded(results) == encoded(complete["results"]))
    # Re-read the current BAM with the admitted absolute tool. This is a bounded
    # file conversion for provenance, not a second aligner/bedtools/QC execution.
    # subprocess.run owns only this fixed samtools child; no process-wide prctl.
    for token in samples:
        bam = _regular(output / token / (token + ".bam"))
        before = sha256_file(bam)
        sam = scratch / (token + ".sam")
        with (
            sam.open("xb") as stdout,
            (scratch / (token + ".stderr")).open("xb") as stderr,
        ):
            run = subprocess.run(
                [str(runtime.tools["samtools"]), "view", "-h", str(bam)],
                cwd=attempt,
                stdin=subprocess.DEVNULL,
                stdout=stdout,
                stderr=stderr,
                timeout=binding.runtime.timeout,
                check=False,
                env={"PATH": "/usr/bin:/bin", "LANG": "C.UTF-8"},
            )
        _require(run.returncode == 0 and sha256_file(bam) == before)
        with sam.open() as stream:
            verify_pairs(
                stream,
                iter_bedpe(
                    output / token / (token + "_all.bedpe.gz"), reference.contigs
                ),
                scratch / (token + ".sqlite"),
            )
    # Reject concurrent input, tool, implementation or completion changes.
    binding.verify()
    _require(_json(attempt / "complete.json") == complete)
    return output, samples, results, sha256_file(attempt / "complete.json")


def _project(workspace, output, samples, results, completion_sha):
    """Copy only accepted public bytes into the existing results/ boundary.

    A partial private filesystem projection is never indexed. Existing bytes
    must match exactly; nothing is overwritten, removed or linked.
    """
    directory = Path(workspace)
    for component in ("results", "hitrac", completion_sha):
        directory /= component
        try:
            directory.mkdir(mode=0o700)
        except FileExistsError:
            _directory(directory)
    specs = [(results["summary"], SUMMARY_TYPE, {"scope": "run"})]
    for token, sample in samples.items():
        for artifact, kind, collection in zip(
            results["samples"][token]["artifacts"],
            (ALL_TYPE, NO_BG_TYPE),
            ("all", "noBg"),
        ):
            specs.append(
                (
                    artifact,
                    kind,
                    {
                        "scope": "sample",
                        "sample_id": sample["display_id"],
                        "collection": collection,
                    },
                )
            )
    candidates = []
    for record, kind, metadata in specs:
        relative = Path(record["path"])
        _require(not relative.is_absolute() and ".." not in relative.parts)
        source = _regular(output / relative)
        _require(
            source.stat().st_size == record["size_bytes"]
            and digest(source) == record["sha256"]
        )
        destination = directory / relative
        try:
            destination.parent.mkdir(mode=0o700)
        except FileExistsError:
            _directory(destination.parent)
        if destination.exists() or destination.is_symlink():
            _require(digest(_regular(destination)) == record["sha256"])
        else:
            with source.open("rb") as stream, destination.open("xb") as target:
                shutil.copyfileobj(stream, target, 1024 * 1024)
            destination.chmod(0o400)
        _require(
            destination.stat().st_size == record["size_bytes"]
            and digest(destination) == record["sha256"]
        )
        candidates.append(
            ExtractedArtifactCandidate(
                output_type=kind,
                relative_path=destination.relative_to(workspace).as_posix(),
                mime_type="text/tab-separated-values"
                if kind == SUMMARY_TYPE
                else "application/gzip",
                metadata=metadata
                | {"sha256": record["sha256"], "qualification_sha256": completion_sha},
            )
        )
    return tuple(candidates)


class HiTracPreprocessResultsAdapter(HiTracPreprocessAdapter):
    """Server-admitted execution with verified, atomic artifact/QC publication."""

    metadata = replace(
        HiTracPreprocessAdapter.metadata,
        description="Pinned tracPre2 preprocessing with verified BEDPE and sample QC results.",
    )
    capabilities = WorkflowCapabilities(
        supports=(
            *HiTracPreprocessAdapter.capabilities.supports,
            "artifact_extract",
            "qc_summary_extract",
        )
    )

    def requires_atomic_result_publication(self):
        return True

    def execution_availability(self):
        if self._runtime is None:
            return super().execution_availability()
        try:
            self._runtime.verify()
        except (OSError, ValueError, TypeError, KeyError):
            return WorkflowAvailability(
                execution="unavailable", reason_code="WORKFLOW_EXECUTION_UNAVAILABLE"
            )
        return WorkflowAvailability(
            execution="available", reason_code="WORKFLOW_EXECUTION_READY"
        )

    def extract_artifacts(self, inputs, workspace):
        if self._binding is None:
            return failure("HITRAC_EXECUTION_BINDING_REQUIRED")
        try:
            output, samples, results, completion_sha = _verified_outputs(
                self._binding, inputs, workspace
            )
            return Result.success(
                _project(Path(workspace), output, samples, results, completion_sha)
            )
        except (
            OSError,
            ValueError,
            TypeError,
            KeyError,
            csv.Error,
            subprocess.SubprocessError,
            CallFailure,
        ):
            return failure("HITRAC_RESULTS_INVALID", "results")

    def qc_source_output_types(self):
        return (SUMMARY_TYPE,)

    def extract_qc_metrics(self, inputs, sources):
        if self._binding is None:
            return failure("HITRAC_EXECUTION_BINDING_REQUIRED")
        try:
            _, _, normalized, samples = self._binding.verify()
            _require(encoded(normalized.to_dict()) == encoded(inputs.to_dict()))
            _require(isinstance(sources, tuple) and len(sources) == 1)
            document = sources[0]
            _require(document.source.output_type == SUMMARY_TYPE)
            path = Path(document.source.relative_path)
            _require(
                path.parts[:2] == ("results", "hitrac")
                and len(path.parts) == 4
                and path.name == "tracPre_summary.txt"
            )
            _require(
                len(path.parts[2]) == 64
                and all(c in "0123456789abcdef" for c in path.parts[2])
            )
            _require(document.source.metadata == {"scope": "run"})
            mapq = normalized.options["mapq"]
            table = parse_summary(document.content, samples, mapq)
            columns = summary_columns(mapq)
            return Result.success(
                tuple(
                    ExtractedQcMetricCandidate(
                        metric_key=METRIC_KEYS[index],
                        display_name=columns[index],
                        # Public durable QC has at most twelve fractional
                        # digits. Preserve full original TSV and validate science
                        # before this display-only ROUND_HALF_EVEN projection.
                        value=(
                            value.quantize(Decimal("1e-12"), rounding=ROUND_HALF_EVEN)
                            if index not in _COUNTS and value.as_tuple().exponent < -12
                            else value
                        ),
                        unit="count"
                        if index in _COUNTS
                        else ("ratio" if index == 10 else "fraction"),
                        scope="sample",
                        sample_id=samples[token]["display_id"],
                        assay="Hi-TrAC",
                        source_artifact_id=document.source.artifact_id,
                    )
                    for token in samples
                    for index, value in enumerate(table[token])
                )
            )
        except (
            OSError,
            ValueError,
            TypeError,
            KeyError,
            AttributeError,
            csv.Error,
            InvalidOperation,
        ):
            return failure("HITRAC_QC_INVALID", "qc_sources")
