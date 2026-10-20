"""Platform reference/control admission for the ENCODE adapter.

The platform path is stricter than the legacy CLI:
- every sample genome must resolve to a configured genome resource;
- one experiment must use one Bowtie2 index;
- referenced controls must match treatment genome/layout;
- external control BAMs must expose a compatible, coordinate-sorted header.

This module is intentionally independent of Snakemake and optional QC tools.
"""

from __future__ import annotations

import gzip
import os
import struct

from encode_pipeline.errors import ValidationError
from encode_pipeline.qc.histone import classify_histone_target

__all__ = [
    "validate_reference_control_admission",
]


def _fragment_cpm_candidate(sample, mode):
    """Match the workflow's PR-E1 fragment eligibility predicate."""
    assay = sample.get("assay", "")
    if assay in ("atac", "mnase"):
        return True
    if assay == "cuttag":
        if mode == "on":
            return True
        info = classify_histone_target(
            sample.get("target", ""), sample.get("peak_mode")
        )
        return info["inferred_histone_class"] != "unknown"
    return False


def _read_bam_header(
    path: str, *, error_cls=ValidationError
) -> tuple[str, dict[str, int]]:
    """Read BAM magic, text header, and binary reference dictionary.

    Returns ``(@HD line, {@SQ name: length})``. The text @SQ dictionary and
    the binary reference dictionary must agree. This is an admission check,
    not a full alignment-integrity scan.
    """
    try:
        with gzip.open(path, "rb") as handle:
            magic = handle.read(4)
            if magic != b"BAM\x01":
                raise error_cls(f"External control BAM has no BAM magic: {path}")
            raw_length = handle.read(4)
            if len(raw_length) != 4:
                raise error_cls(f"External control BAM header is truncated: {path}")
            text_length = struct.unpack("<i", raw_length)[0]
            if text_length < 0 or text_length > 16 * 1024 * 1024:
                raise error_cls(f"External control BAM header is invalid: {path}")
            raw_text = handle.read(text_length)
            if len(raw_text) != text_length:
                raise error_cls(f"External control BAM header is truncated: {path}")

            raw_n_ref = handle.read(4)
            if len(raw_n_ref) != 4:
                raise error_cls(
                    f"External control BAM binary header is truncated: {path}"
                )
            n_ref = struct.unpack("<i", raw_n_ref)[0]
            if n_ref < 0 or n_ref > 10_000_000:
                raise error_cls(f"External control BAM has invalid n_ref: {path}")
            binary_refs: list[tuple[str, int]] = []
            binary_seen: set[str] = set()
            for _ in range(n_ref):
                raw_name_length = handle.read(4)
                if len(raw_name_length) != 4:
                    raise error_cls(
                        f"External control BAM reference name is truncated: {path}"
                    )
                name_length = struct.unpack("<i", raw_name_length)[0]
                if name_length <= 0 or name_length > 16 * 1024 * 1024:
                    raise error_cls(
                        f"External control BAM reference name is invalid: {path}"
                    )
                raw_name = handle.read(name_length)
                if len(raw_name) != name_length or not raw_name.endswith(b"\x00"):
                    raise error_cls(
                        f"External control BAM reference name is truncated: {path}"
                    )
                raw_ref_length = handle.read(4)
                if len(raw_ref_length) != 4:
                    raise error_cls(
                        f"External control BAM reference length is truncated: {path}"
                    )
                ref_length = struct.unpack("<i", raw_ref_length)[0]
                if ref_length < 0:
                    raise error_cls(
                        f"External control BAM reference length is invalid: {path}"
                    )
                name = raw_name[:-1].decode("utf-8", errors="strict")
                if name in binary_seen:
                    raise error_cls(
                        f"External control BAM has duplicate reference name: {path}"
                    )
                binary_seen.add(name)
                binary_refs.append((name, ref_length))
    except (OSError, EOFError, gzip.BadGzipFile) as exc:
        raise error_cls(f"External control BAM is not readable: {path}") from exc

    header_text = raw_text.decode("utf-8", errors="replace")
    hd_line = ""
    text_refs: list[tuple[str, int]] = []
    text_seen: set[str] = set()
    for line in header_text.splitlines():
        if line.startswith("@HD"):
            hd_line = line
        elif line.startswith("@SQ"):
            fields = dict(
                field.split(":", 1) for field in line.split("\t")[1:] if ":" in field
            )
            name = fields.get("SN")
            length = fields.get("LN")
            if not name or not length or not length.isdigit():
                raise error_cls(f"External control BAM has invalid @SQ: {path}")
            if name in text_seen:
                raise error_cls(
                    f"External control BAM has duplicate text @SQ name: {path}"
                )
            text_seen.add(name)
            text_refs.append((name, int(length)))
    if not text_refs:
        raise error_cls(f"External control BAM has no @SQ records: {path}")
    if text_refs != binary_refs:
        raise error_cls(
            f"External control BAM text/binary reference dictionary mismatch: {path}"
        )
    return hd_line, dict(binary_refs)


def _read_chrom_sizes(path: str, *, error_cls=ValidationError) -> dict[str, int]:
    sizes: dict[str, int] = {}
    try:
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 2:
                    continue
                name, length = fields[0], fields[1]
                if name and length.isdigit():
                    sizes[name] = int(length)
    except OSError as exc:
        raise error_cls(f"chrom_sizes is not readable: {path}") from exc
    if not sizes:
        raise error_cls(f"chrom_sizes is empty or invalid: {path}")
    return sizes


def _validate_external_control_bam(
    path: str,
    *,
    genome: str,
    resources: dict,
    error_cls=ValidationError,
) -> None:
    if not os.path.isfile(path):
        raise error_cls(f"control_bam file not found: {path}")
    hd_line, sq = _read_bam_header(path, error_cls=error_cls)
    if "SO:coordinate" not in hd_line.split("\t"):
        raise error_cls(
            f"External control BAM is not declared coordinate-sorted: {path}"
        )

    resource = resources.get(genome, {})
    chrom_sizes = resource.get("chrom_sizes") if isinstance(resource, dict) else ""
    if not chrom_sizes:
        raise error_cls(
            "External control validation requires genome_resources["
            f"{genome!r}].chrom_sizes"
        )
    expected = _read_chrom_sizes(chrom_sizes, error_cls=error_cls)
    for name, length in sq.items():
        if name not in expected:
            raise error_cls(
                f"External control BAM contig {name!r} is absent from chrom_sizes"
            )
        if expected[name] != length:
            raise error_cls(
                f"External control BAM contig {name!r} length {length} "
                f"does not match chrom_sizes {expected[name]}"
            )


def validate_reference_control_admission(
    samples: list[dict],
    config: dict,
    *,
    reference_mode: str,
    external_control_policy: str,
    error_cls=ValidationError,
) -> None:
    """Fail-closed platform admission for reference/control identity.

    ``reference_mode=legacy`` follows the historical lenient behavior.
    ``strict`` enforces genome resource membership, experiment index identity,
    control genome/layout identity, and external control BAM admission.
    """
    strict = reference_mode == "strict"
    resources = config.get("genome_resources", {})
    if strict:
        if not isinstance(resources, dict):
            raise error_cls("genome_resources must be a mapping")
        for sample in samples:
            genome = sample.get("genome", "")
            if genome not in resources:
                raise error_cls(
                    f"Sample {sample.get('id', '?')!r}: genome {genome!r} is not "
                    "present in genome_resources"
                )

    use_control = bool(config.get("use_control", False))
    resources = config.get("genome_resources", {})
    if config.get("blacklist_mode") == "analysis":
        if not isinstance(resources, dict):
            raise error_cls("blacklist_mode=analysis requires genome_resources")
        for sample in samples:
            if sample.get("role") == "treatment" or (
                use_control and sample.get("role") == "control"
            ):
                genome = sample.get("genome", "")
                entry = resources.get(genome, {})
                if not isinstance(entry, dict) or not entry.get("blacklist"):
                    raise error_cls(
                        f"Sample {sample.get('id', '?')!r}: "
                        "blacklist_mode=analysis requires a non-empty "
                        f"genome_resources[{genome!r}].blacklist"
                    )
    if (
        config.get("analysis_profile") == "updated"
        and config.get("fragment_cpm") == "on"
        and not any(
            sample.get("role") == "treatment"
            and sample.get("assay") in ("atac", "cuttag", "mnase")
            for sample in samples
        )
    ):
        raise error_cls(
            "fragment_cpm=on requires at least one ATAC, CUT&Tag, or MNase "
            "treatment sample; ChIP-seq remains read-level in PR-E1."
        )
    if config.get("analysis_profile") == "updated" and config.get(
        "tss_fragment_profile"
    ):
        mode = str(config.get("fragment_cpm", "off"))
        candidates = [
            sample
            for sample in samples
            if sample.get("role") == "treatment"
            and sample.get("assay") in ("atac", "cuttag")
            and _fragment_cpm_candidate(sample, mode)
        ]
        if not candidates:
            raise error_cls(
                "tss_fragment_profile=true requires at least one eligible "
                "ATAC or CUT&Tag fragment sample."
            )
        missing_gtf = [
            sample.get("id", "?")
            for sample in candidates
            if not isinstance(resources, dict)
            or not resources.get(sample.get("genome", ""), {}).get("gtf")
        ]
        if missing_gtf:
            raise error_cls(
                "tss_fragment_profile=true requires genome_resources."
                "<genome>.gtf for every eligible sample; missing for: "
                + ", ".join(missing_gtf)
            )
    if use_control and config.get("analysis_profile") == "updated":
        for sample in samples:
            if sample.get("role") == "treatment" and sample.get("control_bam"):
                raise error_cls(
                    f"Sample {sample.get('id', '?')!r}: analysis_profile=updated "
                    "does not yet support external control_bam. Use a FASTQ-based "
                    "control_sample or analysis_profile=legacy."
                )
    exp_treatments: dict[str, list[dict]] = {}
    for sample in samples:
        if sample.get("role") == "treatment":
            exp_treatments.setdefault(sample.get("experiment", ""), []).append(sample)

    by_id = {sample.get("id"): sample for sample in samples}
    external_controls: dict[str, set[str]] = {}

    for exp, rows in exp_treatments.items():
        if strict:
            treatment_indexes = {row.get("bt2_idx", "") for row in rows}
            if len(treatment_indexes) != 1:
                raise error_cls(
                    f"Experiment {exp!r}: treatment rows use multiple Bowtie2 indexes"
                )
        else:
            treatment_indexes = set()
        if not use_control:
            continue

        controls: list[dict] = []
        for row in rows:
            control_sample = row.get("control_sample", "")
            control_bam = row.get("control_bam", "")
            if control_sample:
                control = by_id.get(control_sample)
                if control is None:
                    raise error_cls(
                        f"Sample {row.get('id')!r}: control_sample "
                        f"{control_sample!r} is missing"
                    )
                if strict:
                    if control.get("experiment") != exp:
                        raise error_cls(
                            f"Experiment {exp!r}: control {control_sample!r} "
                            "belongs to a different experiment"
                        )
                    if control.get("genome") != row.get("genome"):
                        raise error_cls(
                            f"Experiment {exp!r}: control {control_sample!r} genome "
                            "does not match treatment"
                        )
                    if control.get("layout") != row.get("layout"):
                        raise error_cls(
                            f"Experiment {exp!r}: control {control_sample!r} layout "
                            "does not match treatment"
                        )
                controls.append(control)
            elif control_bam:
                external_controls.setdefault(control_bam, set()).add(
                    row.get("genome", "")
                )

        if strict:
            control_indexes = {c.get("bt2_idx", "") for c in controls}
            if control_indexes - treatment_indexes:
                raise error_cls(
                    f"Experiment {exp!r}: control Bowtie2 index does not match treatment"
                )

    if external_control_policy == "validate":
        for path, genomes in sorted(external_controls.items()):
            for genome in sorted(genomes):
                _validate_external_control_bam(
                    path,
                    genome=genome,
                    resources=resources if isinstance(resources, dict) else {},
                    error_cls=error_cls,
                )
