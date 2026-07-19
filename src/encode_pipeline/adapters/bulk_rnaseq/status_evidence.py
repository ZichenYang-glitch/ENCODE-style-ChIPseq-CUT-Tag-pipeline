"""Shared fail-closed sample-status evidence for nf-core/rnaseq 3.26.0."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation, localcontext
from io import BytesIO
import json
import math
from pathlib import PurePosixPath
import re
import struct
from zipfile import BadZipFile, ZipFile


_FASTP_VERSION = "1.0.1"
_FASTQC_VERSION = "0.12.1"
_PLAIN_NUMBER = re.compile(r"(?:0|[1-9][0-9]*)(?:\.[0-9]+)?")
_TABLE_NUMBER = re.compile(r"(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?")
_MAX_STATUS_VALUE_CHARACTERS = 64
_MAX_ZIP_ENTRIES = 256
_MAX_ZIP_ENTRY_BYTES = 16 * 1024 * 1024
_MAX_ZIP_TOTAL_BYTES = 64 * 1024 * 1024
_MAX_SOURCE_BYTES = 16 * 1024 * 1024
_MAX_TEXT_LINES = 20_000
_MAX_LINE_LENGTH = 8_192
_MAX_JSON_DEPTH = 16
_MAX_JSON_NODES = 50_000
_MAX_JSON_STRING = 8_192
_MAX_STAR_LOG_BYTES = 64 * 1024
_MAX_STAR_LINE_CHARACTERS = 1024
_MAX_STAR_TEMPLATE_COUNT = 10**25
_MAX_STATUS_COUNT = 10**25
_CUTADAPT_PERCENT_TOLERANCE = Decimal("0.000000001")
_STAR_PERCENT_FRACTION_TOLERANCE = Decimal("0.00005")
_STAR_TIMESTAMP = re.compile(
    r"(?P<month>Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec) "
    r"(?P<day>0?[1-9]| [1-9]|[12][0-9]|3[01]) "
    r"(?:[01][0-9]|2[0-3]):[0-5][0-9]:(?:[0-5][0-9]|60)"
)
_STAR_MONTH_DAYS = {
    "Jan": 31,
    "Feb": 29,
    "Mar": 31,
    "Apr": 30,
    "May": 31,
    "Jun": 30,
    "Jul": 31,
    "Aug": 31,
    "Sep": 30,
    "Oct": 31,
    "Nov": 30,
    "Dec": 31,
}
_STAR_LOG_FINAL_LAYOUT = (
    "Started job on",
    "Started mapping on",
    "Finished on",
    "Mapping speed, Million of reads per hour",
    "Number of input reads",
    "Average input read length",
    "UNIQUE READS:",
    "Uniquely mapped reads number",
    "Uniquely mapped reads %",
    "Average mapped length",
    "Number of splices: Total",
    "Number of splices: Annotated (sjdb)",
    "Number of splices: GT/AG",
    "Number of splices: GC/AG",
    "Number of splices: AT/AC",
    "Number of splices: Non-canonical",
    "Mismatch rate per base, %",
    "Deletion rate per base",
    "Deletion average length",
    "Insertion rate per base",
    "Insertion average length",
    "MULTI-MAPPING READS:",
    "Number of reads mapped to multiple loci",
    "% of reads mapped to multiple loci",
    "Number of reads mapped to too many loci",
    "% of reads mapped to too many loci",
    "UNMAPPED READS:",
    "Number of reads unmapped: too many mismatches",
    "% of reads unmapped: too many mismatches",
    "Number of reads unmapped: too short",
    "% of reads unmapped: too short",
    "Number of reads unmapped: other",
    "% of reads unmapped: other",
    "CHIMERIC READS:",
    "Number of chimeric reads",
    "% of chimeric reads",
)
_STAR_SECTION_HEADERS = frozenset(
    {
        "UNIQUE READS:",
        "MULTI-MAPPING READS:",
        "UNMAPPED READS:",
        "CHIMERIC READS:",
    }
)


class StatusEvidenceError(ValueError):
    """A fixed status table or its required native evidence is inconsistent."""


class StatusEvidenceContractError(StatusEvidenceError):
    """A status row could not have been emitted by the fixed upstream route."""


@dataclass(frozen=True)
class StatusTable:
    """One optional fixed MultiQC status table."""

    present: bool
    values: Mapping[str, Decimal]


@dataclass(frozen=True)
class ReconciledSampleStatus:
    """Samples whose fixed upstream route actually removes later outputs."""

    trimmed_failed: frozenset[str]
    mapped_failed: frozenset[str]


@dataclass(frozen=True)
class StarLogFinalEvidence:
    """Validated count and printed-percent evidence from fixed STAR 2.7.11b."""

    input_templates: int
    uniquely_mapped_templates: int
    accepted_multimapped_templates: int
    too_many_loci_templates: int
    unmapped_too_many_mismatches_templates: int
    unmapped_too_short_templates: int
    unmapped_other_templates: int
    uniquely_mapped_percent: Decimal


@dataclass(frozen=True)
class CutadaptRowEvidence:
    """One fully validated fixed MultiQC 1.33 Cutadapt row."""

    row_identity: str
    owner: str
    reads_processed: int
    reads_with_adapters: int
    reads_written: int
    bases_processed: int
    quality_trimmed_bases: int
    bases_written: int
    percent_trimmed: Decimal


@dataclass(frozen=True)
class FastpSummaryEvidence:
    """The fixed fastp 1.0.1 counters used for trimming status and QC."""

    input_reads: int
    retained_reads: int
    input_bases: int
    retained_bases: int


def parse_cutadapt_processed_reads(
    content: bytes,
    *,
    expected_row_owners: Mapping[str, str],
) -> Mapping[str, Decimal]:
    """Parse the fixed MultiQC Cutadapt input-read evidence by canonical sample."""
    rows = parse_cutadapt_rows(
        content,
        expected_row_owners=expected_row_owners,
    )
    by_sample: dict[str, list[Decimal]] = {}
    for row in rows:
        by_sample.setdefault(row.owner, []).append(Decimal(row.reads_processed))
    result: dict[str, Decimal] = {}
    for sample, values in by_sample.items():
        if not values or len(set(values)) != 1:
            raise StatusEvidenceError
        result[sample] = values[0]
    return result


def parse_cutadapt_rows(
    content: bytes,
    *,
    expected_row_owners: Mapping[str, str],
) -> tuple[CutadaptRowEvidence, ...]:
    """Validate every status-relevant field in the fixed Cutadapt TSV."""
    text = _decode_bounded_text(content)
    if not text or "\x00" in text or "\r" in text or text.startswith("\ufeff"):
        raise StatusEvidenceError
    lines = text.split("\n")
    if lines[-1] == "":
        lines.pop()
    expected_header = (
        "Sample\tcutadapt_version\tr_processed\tr_with_adapters\tr_written\t"
        "bp_processed\tquality_trimmed\tbp_written\tpercent_trimmed"
    )
    if (
        not expected_row_owners
        or not lines
        or lines[0] != expected_header
        or len(lines) != len(expected_row_owners) + 1
        or any(not line or len(line) > _MAX_LINE_LENGTH for line in lines)
    ):
        raise StatusEvidenceError
    observed: set[str] = set()
    parsed: list[CutadaptRowEvidence] = []
    for line in lines[1:]:
        fields = line.split("\t")
        if (
            len(fields) != 9
            or fields[0] not in expected_row_owners
            or fields[0] in observed
            or fields[1] != "4.0"
        ):
            raise StatusEvidenceError
        counts = tuple(_bounded_nonnegative_int(value) for value in fields[2:8])
        (
            reads_processed,
            reads_with_adapters,
            reads_written,
            bases_processed,
            quality_trimmed_bases,
            bases_written,
        ) = counts
        if (
            len(fields[8]) > _MAX_STATUS_VALUE_CHARACTERS
            or _TABLE_NUMBER.fullmatch(fields[8]) is None
        ):
            raise StatusEvidenceError
        try:
            percent_trimmed = Decimal(fields[8])
        except InvalidOperation as exc:
            raise StatusEvidenceError from exc
        if (
            reads_processed <= 0
            or bases_processed <= 0
            or reads_with_adapters > reads_processed
            # The fixed product surface does not expose Cutadapt's
            # --discard-untrimmed route, so Trim Galore 2.1.0 reports every
            # processed read as written at this stage.
            or reads_written != reads_processed
            or quality_trimmed_bases > bases_processed
            or bases_written > bases_processed
            or quality_trimmed_bases > bases_processed - bases_written
            or not percent_trimmed.is_finite()
            or percent_trimmed < 0
            or percent_trimmed > 100
        ):
            raise StatusEvidenceError
        with localcontext() as context:
            context.prec = 64
            exact_percent = (
                Decimal(bases_processed - bases_written)
                / Decimal(bases_processed)
                * Decimal(100)
            )
        if abs(percent_trimmed - exact_percent) > _CUTADAPT_PERCENT_TOLERANCE:
            raise StatusEvidenceError
        observed.add(fields[0])
        parsed.append(
            CutadaptRowEvidence(
                row_identity=fields[0],
                owner=expected_row_owners[fields[0]],
                reads_processed=reads_processed,
                reads_with_adapters=reads_with_adapters,
                reads_written=reads_written,
                bases_processed=bases_processed,
                quality_trimmed_bases=quality_trimmed_bases,
                bases_written=bases_written,
                percent_trimmed=percent_trimmed,
            )
        )
    if observed != set(expected_row_owners):
        raise StatusEvidenceError
    return tuple(parsed)


def parse_status_table(
    content: bytes | None,
    *,
    kind: str,
    expected_header: str,
    known_samples: frozenset[str],
) -> StatusTable:
    """Parse one exact bounded table without interpreting threshold semantics."""
    if content is None:
        return StatusTable(present=False, values={})
    try:
        text = content.decode("utf-8", errors="strict")
    except UnicodeError as exc:
        raise StatusEvidenceError from exc
    if not text or "\x00" in text or "\r" in text or text.startswith("\ufeff"):
        raise StatusEvidenceError
    lines = text.split("\n")
    if lines[-1] == "":
        lines.pop()
    if (
        not lines
        or lines[0] != expected_header
        or len(lines) > len(known_samples) + 1
        or any(not line for line in lines)
    ):
        raise StatusEvidenceError
    values: dict[str, Decimal] = {}
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) != 2:
            raise StatusEvidenceError
        sample, raw_value = fields
        if (
            sample not in known_samples
            or sample in values
            or not raw_value
            or len(raw_value) > _MAX_STATUS_VALUE_CHARACTERS
            or _TABLE_NUMBER.fullmatch(raw_value) is None
        ):
            raise StatusEvidenceError
        try:
            value = Decimal(raw_value)
        except InvalidOperation as exc:
            raise StatusEvidenceError from exc
        if not value.is_finite() or value < 0:
            raise StatusEvidenceError
        if kind == "trimmed" and value != value.to_integral_value():
            raise StatusEvidenceError
        if kind == "mapped" and value > 100:
            raise StatusEvidenceError
        if kind not in {"trimmed", "mapped"}:
            raise StatusEvidenceError
        values[sample] = value
    return StatusTable(present=True, values=values)


def reconcile_sample_status(
    *,
    known_samples: frozenset[str],
    trimming_enabled: bool,
    trimming_tool: str,
    trim_threshold: Decimal,
    mapped_threshold: Decimal,
    trim_table: StatusTable,
    mapped_table: StatusTable,
    retained_read_evidence: Mapping[str, Decimal],
    trimgalore_input_read_evidence: Mapping[str, Decimal],
    uniquely_mapped_percent_evidence: Mapping[str, Decimal],
) -> ReconciledSampleStatus:
    """Reconcile status rows with the exact native values that created them."""
    if (
        trim_threshold < 0
        or mapped_threshold < 0
        or mapped_threshold > 100
        or trimming_tool not in {"fastp", "trimgalore"}
    ):
        raise StatusEvidenceError
    if not trimming_enabled and trim_table.present:
        raise StatusEvidenceError
    _validate_evidence_samples(known_samples, retained_read_evidence)
    _validate_evidence_samples(known_samples, trimgalore_input_read_evidence)
    _validate_evidence_samples(known_samples, uniquely_mapped_percent_evidence)

    if trimming_enabled:
        trimmed_failed = _reconcile_trim_rows(
            table=trim_table,
            evidence=retained_read_evidence,
            trimgalore_input_reads=trimgalore_input_read_evidence,
            threshold=trim_threshold,
            tool=trimming_tool,
        )
    elif retained_read_evidence:
        raise StatusEvidenceError
    else:
        trimmed_failed = frozenset()
    if trimmed_failed.intersection(mapped_table.values):
        raise StatusEvidenceError
    if any(sample in trimmed_failed for sample in uniquely_mapped_percent_evidence):
        raise StatusEvidenceError
    _reconcile_mapped_rows(
        table=mapped_table,
        evidence=uniquely_mapped_percent_evidence,
        threshold=mapped_threshold,
    )
    return ReconciledSampleStatus(
        trimmed_failed=trimmed_failed,
        mapped_failed=frozenset(mapped_table.values),
    )


def parse_fastp_retained_reads(content: bytes) -> Decimal:
    """Return the exact fixed fastp count used by the upstream status channel."""
    return Decimal(parse_fastp_summary(content).retained_reads)


def parse_fastp_summary(content: bytes) -> FastpSummaryEvidence:
    """Validate the fixed fastp report counters shared by status and QC."""
    try:
        payload = parse_strict_json_document(content)
    except StatusEvidenceError:
        raise
    if (
        not isinstance(payload, Mapping)
        or payload.get("fastp_version") != _FASTP_VERSION
    ):
        raise StatusEvidenceError
    summary = payload.get("summary")
    if not isinstance(summary, Mapping):
        raise StatusEvidenceError
    before = summary.get("before_filtering")
    after = summary.get("after_filtering")
    filtering = payload.get("filtering_result")
    if (
        not isinstance(before, Mapping)
        or not isinstance(after, Mapping)
        or not isinstance(filtering, Mapping)
    ):
        raise StatusEvidenceError
    input_reads = _bounded_json_nonnegative_int(before.get("total_reads"))
    retained_reads = _bounded_json_nonnegative_int(after.get("total_reads"))
    input_bases = _bounded_json_nonnegative_int(before.get("total_bases"))
    retained_bases = _bounded_json_nonnegative_int(after.get("total_bases"))
    passed_filter_reads = _bounded_json_nonnegative_int(
        filtering.get("passed_filter_reads")
    )
    if (
        input_reads <= 0
        or input_bases <= 0
        or retained_reads > input_reads
        or retained_bases > input_bases
        or passed_filter_reads != retained_reads
    ):
        raise StatusEvidenceError
    return FastpSummaryEvidence(
        input_reads=input_reads,
        retained_reads=retained_reads,
        input_bases=input_bases,
        retained_bases=retained_bases,
    )


def parse_strict_json_document(content: bytes) -> object:
    """Parse one bounded JSON document with duplicate-key rejection."""
    text = _decode_bounded_text(content)
    _preflight_json_text(text)

    def reject_constant(_value: str) -> None:
        raise StatusEvidenceError

    def pairs(values: list[tuple[str, object]]) -> dict[str, object]:
        result: dict[str, object] = {}
        for key, value in values:
            if key in result:
                raise StatusEvidenceError
            result[key] = value
        return result

    try:
        value = json.loads(
            text,
            parse_float=Decimal,
            parse_int=int,
            parse_constant=reject_constant,
            object_pairs_hook=pairs,
        )
    except (json.JSONDecodeError, RecursionError, TypeError, ValueError) as exc:
        if isinstance(exc, StatusEvidenceError):
            raise
        raise StatusEvidenceError from exc
    _validate_json_tree(value)
    return value


def parse_fastqc_total_sequences(
    content: bytes,
    *,
    expected_root: str,
    expected_filename: str,
) -> Decimal:
    """Read the exact Basic Statistics count from one bounded fixed FastQC ZIP."""
    try:
        with ZipFile(BytesIO(content)) as archive:
            infos = archive.infolist()
            if len(infos) > _MAX_ZIP_ENTRIES:
                raise StatusEvidenceError
            names: set[str] = set()
            total_size = 0
            for info in infos:
                path = PurePosixPath(info.filename)
                if (
                    info.filename in names
                    or path.as_posix() != info.filename
                    or path.is_absolute()
                    or any(part in {"", ".", ".."} for part in path.parts)
                    or info.file_size < 0
                    or info.file_size > _MAX_ZIP_ENTRY_BYTES
                ):
                    raise StatusEvidenceError
                names.add(info.filename)
                total_size += info.file_size
                if total_size > _MAX_ZIP_TOTAL_BYTES:
                    raise StatusEvidenceError
            member = f"{expected_root}/fastqc_data.txt"
            if member not in names:
                raise StatusEvidenceError
            data = archive.read(member)
    except (BadZipFile, KeyError, OSError, RuntimeError, ValueError) as exc:
        if isinstance(exc, StatusEvidenceError):
            raise
        raise StatusEvidenceError from exc
    try:
        text = data.decode("utf-8", errors="strict")
    except UnicodeError as exc:
        raise StatusEvidenceError from exc
    if not text.startswith(f"##FastQC\t{_FASTQC_VERSION}\n"):
        raise StatusEvidenceError
    in_basic = False
    filename: str | None = None
    total: Decimal | None = None
    for line in text.splitlines()[1:]:
        if line.startswith(">>Basic Statistics\t"):
            if in_basic:
                raise StatusEvidenceError
            in_basic = True
            continue
        if in_basic and line == ">>END_MODULE":
            in_basic = False
            continue
        if not in_basic or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 2:
            raise StatusEvidenceError
        if fields[0] == "Filename":
            if filename is not None:
                raise StatusEvidenceError
            filename = fields[1]
        elif fields[0] == "Total Sequences":
            if total is not None or _PLAIN_NUMBER.fullmatch(fields[1]) is None:
                raise StatusEvidenceError
            total = Decimal(fields[1])
    if (
        filename != expected_filename
        or total is None
        or total != total.to_integral_value()
    ):
        raise StatusEvidenceError
    return total


def parse_star_log_final(content: bytes) -> StarLogFinalEvidence:
    """Validate the complete fixed STAR 2.7.11b Log.final.out contract."""
    if len(content) > _MAX_STAR_LOG_BYTES or b"\x00" in content or b"\r" in content:
        raise StatusEvidenceError
    try:
        text = content.decode("utf-8", errors="strict")
    except UnicodeError as exc:
        raise StatusEvidenceError from exc
    if not text or text.startswith("\ufeff"):
        raise StatusEvidenceError
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if len(lines) != len(_STAR_LOG_FINAL_LAYOUT) or any(
        len(line) > _MAX_STAR_LINE_CHARACTERS for line in lines
    ):
        raise StatusEvidenceError
    values: dict[str, str] = {}
    for line, expected in zip(lines, _STAR_LOG_FINAL_LAYOUT, strict=True):
        if expected in _STAR_SECTION_HEADERS:
            if line != expected:
                raise StatusEvidenceError
            continue
        fields = line.split("|")
        if len(fields) != 2:
            raise StatusEvidenceError
        key, raw = (field.strip() for field in fields)
        if key != expected or not raw or key in values:
            raise StatusEvidenceError
        values[key] = raw

    for key in ("Started job on", "Started mapping on", "Finished on"):
        _validate_star_timestamp(values[key])
    _star_nonnegative_decimal(values["Mapping speed, Million of reads per hour"])
    average_input_length = _star_nonnegative_decimal(
        values["Average input read length"]
    )
    average_mapped_length = _star_nonnegative_decimal(values["Average mapped length"])

    input_templates = _star_nonnegative_int(values["Number of input reads"])
    splice_total = _star_nonnegative_int(values["Number of splices: Total"])
    splice_annotated = _star_nonnegative_int(
        values["Number of splices: Annotated (sjdb)"]
    )
    splice_motifs = tuple(
        _star_nonnegative_int(values[key])
        for key in (
            "Number of splices: GT/AG",
            "Number of splices: GC/AG",
            "Number of splices: AT/AC",
            "Number of splices: Non-canonical",
        )
    )
    for key in (
        "Mismatch rate per base, %",
        "Deletion rate per base",
        "Insertion rate per base",
    ):
        _star_percent(values[key])
    for key in ("Deletion average length", "Insertion average length"):
        _star_nonnegative_decimal(values[key])

    counts = {
        "unique": _star_nonnegative_int(values["Uniquely mapped reads number"]),
        "multi": _star_nonnegative_int(
            values["Number of reads mapped to multiple loci"]
        ),
        "too_many": _star_nonnegative_int(
            values["Number of reads mapped to too many loci"]
        ),
        "too_many_mismatches": _star_nonnegative_int(
            values["Number of reads unmapped: too many mismatches"]
        ),
        "too_short": _star_nonnegative_int(
            values["Number of reads unmapped: too short"]
        ),
        "other": _star_nonnegative_int(values["Number of reads unmapped: other"]),
    }
    chimeric_count = _star_nonnegative_int(values["Number of chimeric reads"])
    accepted_mapped_count = counts["unique"] + counts["multi"]
    reported_percentages = {
        "unique": _star_percent(values["Uniquely mapped reads %"]),
        "multi": _star_percent(values["% of reads mapped to multiple loci"]),
        "too_many": _star_percent(values["% of reads mapped to too many loci"]),
        "too_many_mismatches": _star_percent(
            values["% of reads unmapped: too many mismatches"]
        ),
        "too_short": _star_percent(values["% of reads unmapped: too short"]),
        "other": _star_percent(values["% of reads unmapped: other"]),
    }
    chimeric_percentage = _star_percent(values["% of chimeric reads"])
    if (
        input_templates <= 0
        or average_input_length <= 0
        or any(count > input_templates for count in counts.values())
        or sum(counts.values()) != input_templates
        or chimeric_count > input_templates
        or splice_annotated > splice_total
        or sum(splice_motifs) != splice_total
        or (accepted_mapped_count > 0 and average_mapped_length <= 0)
        or (
            accepted_mapped_count == 0
            and (average_mapped_length != 0 or splice_total != 0)
        )
    ):
        raise StatusEvidenceError
    with localcontext() as context:
        context.prec = 64
        count_percentage_pairs = (
            *(
                (counts[name], percentage)
                for name, percentage in reported_percentages.items()
            ),
            (chimeric_count, chimeric_percentage),
        )
        for count, reported_percentage in count_percentage_pairs:
            exact_fraction = Decimal(count) / Decimal(input_templates)
            reported_fraction = reported_percentage / Decimal(100)
            if (
                abs(exact_fraction - reported_fraction)
                > _STAR_PERCENT_FRACTION_TOLERANCE
            ):
                raise StatusEvidenceError
    return StarLogFinalEvidence(
        input_templates=input_templates,
        uniquely_mapped_templates=counts["unique"],
        accepted_multimapped_templates=counts["multi"],
        too_many_loci_templates=counts["too_many"],
        unmapped_too_many_mismatches_templates=counts["too_many_mismatches"],
        unmapped_too_short_templates=counts["too_short"],
        unmapped_other_templates=counts["other"],
        uniquely_mapped_percent=reported_percentages["unique"],
    )


def parse_star_uniquely_mapped_percent(content: bytes) -> Decimal:
    """Return the validated literal used by the pinned STAR status route."""
    return parse_star_log_final(content).uniquely_mapped_percent


def _star_nonnegative_int(raw: str) -> int:
    if re.fullmatch(r"0|[1-9][0-9]*", raw) is None:
        raise StatusEvidenceError
    value = int(raw)
    if value > _MAX_STAR_TEMPLATE_COUNT:
        raise StatusEvidenceError
    return value


def _star_percent(raw: str) -> Decimal:
    if not raw.endswith("%") or _PLAIN_NUMBER.fullmatch(raw[:-1]) is None:
        raise StatusEvidenceError
    value = Decimal(raw[:-1])
    if not value.is_finite() or value > 100:
        raise StatusEvidenceError
    return value


def _star_nonnegative_decimal(raw: str) -> Decimal:
    if _PLAIN_NUMBER.fullmatch(raw) is None:
        raise StatusEvidenceError
    try:
        value = Decimal(raw)
    except InvalidOperation as exc:
        raise StatusEvidenceError from exc
    if not value.is_finite() or value < 0 or value > Decimal(_MAX_STAR_TEMPLATE_COUNT):
        raise StatusEvidenceError
    return value


def _validate_star_timestamp(raw: str) -> None:
    match = _STAR_TIMESTAMP.fullmatch(raw)
    if match is None:
        raise StatusEvidenceError
    if int(match.group("day")) > _STAR_MONTH_DAYS[match.group("month")]:
        raise StatusEvidenceError


def _bounded_nonnegative_int(raw: str) -> int:
    if (
        len(raw) > _MAX_STATUS_VALUE_CHARACTERS
        or re.fullmatch(r"0|[1-9][0-9]*", raw) is None
    ):
        raise StatusEvidenceError
    try:
        value = int(raw)
    except ValueError as exc:
        raise StatusEvidenceError from exc
    if value > _MAX_STATUS_COUNT:
        raise StatusEvidenceError
    return value


def _bounded_json_nonnegative_int(raw: object) -> int:
    if isinstance(raw, bool) or not isinstance(raw, int):
        raise StatusEvidenceError
    if raw < 0 or raw > _MAX_STATUS_COUNT:
        raise StatusEvidenceError
    return raw


def _decode_bounded_text(content: bytes) -> str:
    if len(content) > _MAX_SOURCE_BYTES or b"\x00" in content:
        raise StatusEvidenceError
    try:
        text = content.decode("utf-8", errors="strict")
    except UnicodeError as exc:
        raise StatusEvidenceError from exc
    lines = text.splitlines()
    if len(lines) > _MAX_TEXT_LINES or any(
        len(line) > _MAX_LINE_LENGTH for line in lines
    ):
        raise StatusEvidenceError
    return text


def _preflight_json_text(text: str) -> None:
    depth = 0
    structural_tokens = 0
    in_string = False
    escaped = False
    string_length = 0
    for character in text:
        if in_string:
            if escaped:
                escaped = False
                string_length += 1
            elif character == "\\":
                escaped = True
            elif character == '"':
                in_string = False
            else:
                string_length += 1
            if string_length > _MAX_JSON_STRING:
                raise StatusEvidenceError
            continue
        if character == '"':
            in_string = True
            string_length = 0
        elif character in "[{":
            depth += 1
            structural_tokens += 1
            if depth > _MAX_JSON_DEPTH:
                raise StatusEvidenceError
        elif character in "]}":
            depth -= 1
            if depth < 0:
                raise StatusEvidenceError
        elif character in ",:":
            structural_tokens += 1
        if structural_tokens > _MAX_JSON_NODES:
            raise StatusEvidenceError
    if in_string or escaped or depth != 0:
        raise StatusEvidenceError


def _validate_json_tree(root: object) -> None:
    nodes = 0
    stack: list[tuple[object, int]] = [(root, 1)]
    while stack:
        value, depth = stack.pop()
        nodes += 1
        if nodes > _MAX_JSON_NODES or depth > _MAX_JSON_DEPTH:
            raise StatusEvidenceError
        if isinstance(value, str):
            if len(value) > _MAX_JSON_STRING or "\x00" in value:
                raise StatusEvidenceError
        elif isinstance(value, Mapping):
            for key, child in value.items():
                if not isinstance(key, str) or len(key) > 256 or "\x00" in key:
                    raise StatusEvidenceError
                stack.append((child, depth + 1))
        elif isinstance(value, list):
            stack.extend((child, depth + 1) for child in value)
        elif isinstance(value, Decimal):
            if not value.is_finite():
                raise StatusEvidenceError
        elif value is not None and not isinstance(value, (bool, int)):
            raise StatusEvidenceError


def _reconcile_mapped_rows(
    *,
    table: StatusTable,
    evidence: Mapping[str, Decimal],
    threshold: Decimal,
) -> None:
    normalized_threshold = _binary32_decimal(threshold)
    normalized_table = {
        sample: _binary32_decimal(value) for sample, value in table.values.items()
    }
    normalized_evidence = {
        sample: _binary32_decimal(value) for sample, value in evidence.items()
    }
    for sample, value in normalized_table.items():
        if value >= normalized_threshold:
            raise StatusEvidenceContractError
        if normalized_evidence.get(sample) != value:
            raise StatusEvidenceError
    for sample, value in normalized_evidence.items():
        if (value < normalized_threshold) != (sample in normalized_table):
            raise StatusEvidenceError


def _reconcile_trim_rows(
    *,
    table: StatusTable,
    evidence: Mapping[str, Decimal],
    trimgalore_input_reads: Mapping[str, Decimal],
    threshold: Decimal,
    tool: str,
) -> frozenset[str]:
    if threshold != threshold.to_integral_value():
        raise StatusEvidenceError
    table_membership_threshold = _binary32_integral(threshold)
    if tool == "trimgalore":
        normalized_table = {
            sample: _exact_integral(value) for sample, value in table.values.items()
        }
        if not set(evidence).issubset(trimgalore_input_reads):
            raise StatusEvidenceError
        route_values: dict[str, Decimal] = {}
        for sample, retained in evidence.items():
            total = _exact_integral(trimgalore_input_reads[sample])
            retained = _exact_integral(retained)
            if retained > total:
                raise StatusEvidenceError
            # Trim Galore 2.1.0 increments both mate counters whenever either
            # mate is too short. Therefore SE records and PE pairs share this
            # exact retained-count identity for the fixed nf-core route.
            route_values[sample] = _binary32_subtract_integrals(
                total,
                total - retained,
            )
        normalized_evidence = route_values
        downstream_threshold = table_membership_threshold
    elif tool == "fastp":
        if trimgalore_input_reads:
            raise StatusEvidenceError
        normalized_table = {
            sample: _exact_integral(value) for sample, value in table.values.items()
        }
        normalized_evidence = {
            sample: _exact_integral(value) for sample, value in evidence.items()
        }
        downstream_threshold = threshold
    else:
        raise StatusEvidenceError

    for sample, value in normalized_table.items():
        if value > table_membership_threshold:
            raise StatusEvidenceContractError
        if normalized_evidence.get(sample) != value:
            raise StatusEvidenceError
    for sample, value in normalized_evidence.items():
        if (value <= table_membership_threshold) != (sample in normalized_table):
            raise StatusEvidenceError
    return frozenset(
        sample
        for sample, value in normalized_evidence.items()
        if value < downstream_threshold
    )


def _exact_integral(value: Decimal) -> Decimal:
    if not value.is_finite() or value < 0 or value != value.to_integral_value():
        raise StatusEvidenceError
    return value


def _binary32_integral(value: Decimal) -> Decimal:
    value = _exact_integral(value)
    return _binary32_decimal(value)


def _binary32_decimal(value: Decimal) -> Decimal:
    if not value.is_finite() or value < 0:
        raise StatusEvidenceError
    try:
        packed = struct.pack(">f", float(value))
        rounded = struct.unpack(">f", packed)[0]
    except (OverflowError, struct.error, ValueError) as exc:
        raise StatusEvidenceError from exc
    if not math.isfinite(rounded) or rounded < 0:
        raise StatusEvidenceError
    result = Decimal.from_float(rounded)
    return result


def _binary32_subtract_integrals(left: Decimal, right: Decimal) -> Decimal:
    left_float = float(_binary32_integral(left))
    right_float = float(_binary32_integral(right))
    # Pinned Nextflow/Groovy promotes Float subtraction to java.lang.Double.
    # Python's binary64 subtraction of the two exact binary32 operands matches it.
    difference = left_float - right_float
    if not math.isfinite(difference) or difference < 0:
        raise StatusEvidenceError
    result = Decimal.from_float(difference)
    if result != result.to_integral_value():
        raise StatusEvidenceError
    return result


def _validate_evidence_samples(
    known_samples: frozenset[str],
    evidence: Mapping[str, Decimal],
) -> None:
    for sample, value in evidence.items():
        if (
            sample not in known_samples
            or not isinstance(value, Decimal)
            or not value.is_finite()
            or value < 0
        ):
            raise StatusEvidenceError
