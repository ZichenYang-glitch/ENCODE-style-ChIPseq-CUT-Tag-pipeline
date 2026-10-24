"""A2 platform reference/control admission tests."""

from __future__ import annotations

import gzip
from pathlib import Path
import struct

import pytest
import yaml

from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter
from encode_pipeline.cli._validator import main as cli_validator_main
from encode_pipeline.errors import ValidationError
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.samples.admission import (
    _get_reference_contigs,
    _read_fai,
    _validate_blacklist_bed,
    validate_reference_control_admission,
)


HEADER = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n"


def _write_bam(
    path: Path,
    header: str = HEADER,
    refs: list[tuple[str, int]] | None = None,
) -> Path:
    payload = header.encode("utf-8")
    if refs is None:
        refs = []
        for line in header.splitlines():
            if not line.startswith("@SQ"):
                continue
            fields = dict(
                field.split(":", 1) for field in line.split("\t")[1:] if ":" in field
            )
            refs.append((fields["SN"], int(fields["LN"])))
    with gzip.open(path, "wb") as handle:
        handle.write(b"BAM\x01")
        handle.write(struct.pack("<i", len(payload)))
        handle.write(payload)
        handle.write(struct.pack("<i", len(refs)))
        for name, length in refs:
            encoded = name.encode("utf-8") + b"\x00"
            handle.write(struct.pack("<i", len(encoded)))
            handle.write(encoded)
            handle.write(struct.pack("<i", length))
    return path


def _truncate_after_text(path: Path, header: str = HEADER) -> Path:
    payload = header.encode("utf-8")
    with gzip.open(path, "wb") as handle:
        handle.write(b"BAM\x01")
        handle.write(struct.pack("<i", len(payload)))
        handle.write(payload)
    return path


def _binary_mismatch(path: Path, header: str = HEADER) -> Path:
    payload = header.encode("utf-8")
    with gzip.open(path, "wb") as handle:
        handle.write(b"BAM\x01")
        handle.write(struct.pack("<i", len(payload)))
        handle.write(payload)
        handle.write(struct.pack("<i", 1))
        name = b"chr2\x00"
        handle.write(struct.pack("<i", len(name)))
        handle.write(name)
        handle.write(struct.pack("<i", 1000))
    return path


def _write_sizes(path: Path, rows: str = "chr1\t1000\n") -> Path:
    path.write_text(rows, encoding="utf-8")
    return path


def _sample(**overrides):
    row = {
        "id": "S1",
        "role": "treatment",
        "experiment": "EXP1",
        "genome": "hg38",
        "layout": "PE",
        "bt2_idx": "/idx/hg38",
        "control_sample": "",
        "control_bam": "",
    }
    row.update(overrides)
    return row


def _config(chrom_sizes: Path | None = None, use_control: bool = True):
    sizes = str(chrom_sizes) if chrom_sizes else ""
    resources = {
        "hg38": {"chrom_sizes": sizes},
        "mm10": {"chrom_sizes": sizes},
    }
    return {
        "genome_resources": resources,
        "use_control": use_control,
        "blacklist_mode": "off",
    }


def _admit(samples, config, *, reference_mode="strict", external="validate"):
    validate_reference_control_admission(
        samples,
        config,
        reference_mode=reference_mode,
        external_control_policy=external,
        error_cls=ValidationError,
    )


def test_updated_tss_fragment_requires_eligible_gtf_sample():
    config = _config(use_control=False)
    config.update(
        {
            "analysis_profile": "updated",
            "fragment_cpm": "auto",
            "tss_fragment_profile": True,
        }
    )
    with pytest.raises(ValidationError, match="tss_fragment_profile=true"):
        _admit(
            [_sample(assay="chipseq")], config, reference_mode="legacy", external="warn"
        )

    config["genome_resources"]["hg38"]["gtf"] = "/tmp/annotation.gtf"
    _admit([_sample(assay="atac")], config, reference_mode="legacy", external="warn")


def test_updated_tss_fragment_rejects_unknown_cuttag_target_in_auto():
    config = _config(use_control=False)
    config["genome_resources"]["hg38"]["gtf"] = "/tmp/annotation.gtf"
    config.update(
        {
            "analysis_profile": "updated",
            "fragment_cpm": "auto",
            "tss_fragment_profile": True,
        }
    )
    with pytest.raises(ValidationError, match="at least one eligible"):
        _admit(
            [_sample(assay="cuttag", target="UNKNOWN_TARGET", peak_mode="narrow")],
            config,
            reference_mode="legacy",
            external="warn",
        )


def test_updated_tss_fragment_requires_gtf_for_every_candidate():
    config = _config(use_control=False)
    config["genome_resources"]["hg38"]["gtf"] = "/tmp/annotation.gtf"
    config.update(
        {
            "analysis_profile": "updated",
            "fragment_cpm": "auto",
            "tss_fragment_profile": True,
        }
    )
    with pytest.raises(ValidationError, match="missing for: S2"):
        _admit(
            [
                _sample(id="S1", assay="atac", genome="hg38"),
                _sample(id="S2", assay="atac", genome="mm10"),
            ],
            config,
            reference_mode="legacy",
            external="warn",
        )


def test_updated_tss_fragment_on_allows_unknown_cuttag_with_gtf():
    config = _config(use_control=False)
    config["genome_resources"]["hg38"]["gtf"] = "/tmp/annotation.gtf"
    config.update(
        {
            "analysis_profile": "updated",
            "fragment_cpm": "on",
            "tss_fragment_profile": True,
        }
    )
    _admit(
        [_sample(assay="cuttag", target="UNKNOWN_TARGET", peak_mode="narrow")],
        config,
        reference_mode="legacy",
        external="warn",
    )


def test_updated_fragment_cpm_on_requires_eligible_assay():
    config = _config(use_control=False)
    config.update({"analysis_profile": "updated", "fragment_cpm": "on"})
    with pytest.raises(ValidationError, match="fragment_cpm=on requires"):
        _admit(
            [_sample(assay="chipseq")], config, reference_mode="legacy", external="warn"
        )
    _admit([_sample(assay="atac")], config, reference_mode="legacy", external="warn")


def test_updated_rejects_external_control_with_legacy_reference_mode(tmp_path):
    bam = _write_bam(tmp_path / "control.bam")
    config = _config(use_control=True)
    config["analysis_profile"] = "updated"

    with pytest.raises(ValidationError, match="external control_bam"):
        _admit(
            [_sample(control_bam=str(bam))],
            config,
            reference_mode="legacy",
            external="warn",
        )


def test_updated_ignores_external_control_when_use_control_false(tmp_path):
    bam = _write_bam(tmp_path / "control.bam")
    config = _config(use_control=False)
    config["analysis_profile"] = "updated"

    _admit(
        [_sample(control_bam=str(bam))],
        config,
        reference_mode="legacy",
        external="warn",
    )


def test_strict_rejects_genome_absent_from_resources():
    with pytest.raises(ValidationError, match="not present in genome_resources"):
        _admit([_sample(genome="not_a_genome")], _config())


def test_strict_rejects_multiple_indexes_in_experiment():
    samples = [
        _sample(id="S1", bt2_idx="/idx/a"),
        _sample(id="S2", bt2_idx="/idx/b"),
    ]
    with pytest.raises(ValidationError, match="multiple Bowtie2 indexes"):
        _admit(samples, _config())


def test_strict_rejects_control_genome_or_layout_mismatch():
    treatment = _sample(control_sample="C1")
    control = _sample(
        id="C1",
        role="control",
        genome="mm10",
        layout="SE",
        control_sample="",
    )
    with pytest.raises(ValidationError, match="genome does not match"):
        _admit([treatment, control], _config())


def test_strict_rejects_singleton_control_from_other_experiment():
    treatment = _sample(control_sample="C1")
    control = _sample(
        id="C1",
        role="control",
        experiment="EXP2",
        control_sample="",
    )
    with pytest.raises(ValidationError, match="different experiment"):
        _admit([treatment, control], _config())


def test_strict_external_control_accepts_matching_header(tmp_path: Path):
    sizes = _write_sizes(tmp_path / "sizes")
    bam = _write_bam(tmp_path / "control.bam")
    treatment = _sample(control_bam=str(bam))
    _admit([treatment], _config(sizes))


def test_strict_external_control_accepts_real_bgzf_bam(tmp_path: Path):
    sizes = _write_sizes(tmp_path / "sizes")
    bam = Path(__file__).resolve().parents[1] / "fixtures/a2/tiny_control.bgzf"
    treatment = _sample(control_bam=str(bam))
    _admit([treatment], _config(sizes))


def test_strict_external_control_rejects_text_only_truncated_header(tmp_path: Path):
    sizes = _write_sizes(tmp_path / "sizes")
    bam = _truncate_after_text(tmp_path / "control.bam")
    treatment = _sample(control_bam=str(bam))
    with pytest.raises(ValidationError, match="binary header is truncated"):
        _admit([treatment], _config(sizes))


def test_strict_external_control_rejects_text_binary_mismatch(tmp_path: Path):
    sizes = _write_sizes(tmp_path / "sizes")
    bam = _binary_mismatch(tmp_path / "control.bam")
    treatment = _sample(control_bam=str(bam))
    with pytest.raises(
        ValidationError, match="text/binary reference dictionary mismatch"
    ):
        _admit([treatment], _config(sizes))


def test_strict_external_control_rejects_reference_order_mismatch(tmp_path: Path):
    sizes = _write_sizes(tmp_path / "sizes", "chr1\t1000\nchr2\t2000\n")
    header = (
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n"
    )
    bam = _write_bam(
        tmp_path / "control.bam",
        header=header,
        refs=[("chr2", 2000), ("chr1", 1000)],
    )
    treatment = _sample(control_bam=str(bam))
    with pytest.raises(
        ValidationError, match="text/binary reference dictionary mismatch"
    ):
        _admit([treatment], _config(sizes))


def test_strict_external_control_rejects_duplicate_text_reference(tmp_path: Path):
    sizes = _write_sizes(tmp_path / "sizes")
    header = (
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr1\tLN:1000\n"
    )
    bam = _write_bam(
        tmp_path / "control.bam",
        header=header,
        refs=[("chr1", 1000)],
    )
    treatment = _sample(control_bam=str(bam))
    with pytest.raises(ValidationError, match="duplicate text @SQ"):
        _admit([treatment], _config(sizes))


def test_strict_external_control_rejects_wrong_contig(tmp_path: Path):
    sizes = _write_sizes(tmp_path / "sizes")
    bam = _write_bam(
        tmp_path / "control.bam",
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr2\tLN:1000\n",
    )
    treatment = _sample(control_bam=str(bam))
    with pytest.raises(ValidationError, match="absent from chrom_sizes"):
        _admit([treatment], _config(sizes))


def test_strict_external_control_rejects_missing_coordinate_sort(tmp_path: Path):
    sizes = _write_sizes(tmp_path / "sizes")
    bam = _write_bam(
        tmp_path / "control.bam",
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n",
    )
    treatment = _sample(control_bam=str(bam))
    with pytest.raises(ValidationError, match="coordinate-sorted"):
        _admit([treatment], _config(sizes))


def test_strict_external_control_requires_chrom_sizes(tmp_path: Path):
    bam = _write_bam(tmp_path / "control.bam")
    treatment = _sample(control_bam=str(bam))
    with pytest.raises(ValidationError, match="requires genome_resources"):
        _admit([treatment], _config(None))


def test_legacy_mode_keeps_lenient_behavior():
    _admit([_sample(genome="not_a_genome")], _config(), reference_mode="legacy")


def _sha256(value: bytes) -> str:
    import hashlib

    return hashlib.sha256(value).hexdigest()


def _write(path: Path, value: bytes) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(value)
    return _sha256(value)


def _binding(root: Path) -> dict[str, object]:
    resources: dict[str, object] = {}
    for name, content in (
        ("reference_fasta", b">chr1\nACGT\n"),
        ("gtf", b'chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id "g1";\n'),
        ("chrom_sizes", b"chr1\t1000\n"),
        ("blacklist", b"chr1\t1\t2\n"),
    ):
        suffix = {
            "reference_fasta": "fa",
            "gtf": "gtf",
            "chrom_sizes": "sizes",
            "blacklist": "bed",
        }[name]
        path = root / f"reference/{name}.{suffix}"
        resources[name] = {"path": str(path), "sha256": _write(path, content)}

    prefix = root / "reference/bowtie2/GRCh38"
    files: dict[str, str] = {}
    for suffix in (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"):
        path = Path(f"{prefix}{suffix}")
        files[suffix] = _write(path, suffix.encode("ascii"))
    return {
        "schema_version": "encode-reference-binding-v1",
        "assembly": "GRCh38",
        "effective_genome_size": "hs",
        "genome_resources": resources,
        "bowtie2_index": {"prefix": str(prefix), "files": files},
    }


def _inputs(root: Path) -> WorkflowInputs:
    return WorkflowInputs(
        config={
            "threads": 1,
            "reference_mode": "legacy",
            "external_control_policy": "warn",
            "blacklist_mode": "off",
        },
        samples=[
            {
                "sample": "S1",
                "fastq_1": str(root / "S1_R1.fastq.gz"),
                "fastq_2": str(root / "S1_R2.fastq.gz"),
                "layout": "PE",
                "assay": "chipseq",
                "target": "CTCF",
                "peak_mode": "narrow",
            }
        ],
        options={},
    )


def test_platform_admission_overrides_caller_downgrade(tmp_path: Path):
    adapter = EncodeStyleWorkflowAdapter()
    bound = adapter.bind_reference_profile(_inputs(tmp_path), _binding(tmp_path))
    assert bound.is_success
    result = bound.value.adapter.validate(bound.value.inputs)
    assert result.is_success, result.issues
    config = result.value["config"]
    assert config["reference_mode"] == "strict"
    assert config["external_control_policy"] == "validate"
    assert any(
        issue.code == "ENCODE_PLATFORM_ADMISSION_ENFORCED" for issue in result.issues
    )


def test_platform_admission_accepts_updated_execution(tmp_path: Path):
    adapter = EncodeStyleWorkflowAdapter()
    inputs = _inputs(tmp_path)
    inputs.config["analysis_profile"] = "updated"
    bound = adapter.bind_reference_profile(inputs, _binding(tmp_path))
    assert bound.is_success
    result = bound.value.adapter.validate(bound.value.inputs)
    assert result.is_success, result.issues
    assert result.value["config"]["analysis_profile"] == "updated"


def test_bare_bundle_rejects_updated_external_control(tmp_path: Path):
    bam = _write_bam(tmp_path / "control.bam")
    inputs = _inputs(tmp_path)
    inputs.config.update({"analysis_profile": "updated", "use_control": True})
    inputs.samples[0].update(
        {
            "genome": "hg38",
            "bowtie2_index": "/idx/hg38",
            "control_bam": str(bam),
        }
    )

    result = EncodeStyleWorkflowAdapter().validate(inputs)

    assert result.is_failure
    assert any(
        issue.code == "ENCODE_REFERENCE_ADMISSION_INVALID" for issue in result.issues
    )


def test_bare_bundle_updated_ignores_external_control_when_disabled(tmp_path: Path):
    bam = _write_bam(tmp_path / "control.bam")
    inputs = _inputs(tmp_path)
    inputs.config.update({"analysis_profile": "updated", "use_control": False})
    inputs.samples[0].update(
        {
            "genome": "hg38",
            "bowtie2_index": "/idx/hg38",
            "control_bam": str(bam),
        }
    )

    result = EncodeStyleWorkflowAdapter().validate(inputs)

    assert result.is_success, result.issues


def _cli_fixture(tmp_path: Path, *, config: dict) -> Path:
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
        "bowtie2_index\n"
        "S1\t/abs/S1_R1.fq.gz\t/abs/S1_R2.fq.gz\tPE\tchipseq\tCTCF\tnarrow\thg38\t"
        "/idx/hg38\n",
        encoding="utf-8",
    )
    config["samples"] = str(samples)
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config), encoding="utf-8")
    return path


def test_cli_accepts_updated_execution(tmp_path: Path):
    path = _cli_fixture(
        tmp_path,
        config={
            "threads": 1,
            "use_control": False,
            "analysis_profile": "updated",
            "blacklist_mode": "off",
            "genome_resources": {"hg38": {"effective_genome_size": "hs"}},
        },
    )
    assert cli_validator_main(["--config", str(path)]) is None


def test_cli_executes_strict_admission(tmp_path: Path):
    path = _cli_fixture(
        tmp_path,
        config={
            "threads": 1,
            "use_control": False,
            "reference_mode": "strict",
            "genome_resources": {"hs": {"effective_genome_size": "hs"}},
        },
    )
    with pytest.raises(SystemExit):
        cli_validator_main(["--config", str(path)])


def test_validate_blacklist_bed_rules(tmp_path: Path):
    # Valid BED (tab and space)
    bed = tmp_path / "valid.bed"
    bed.write_text("chr1\t10\t50\nchr1 60 100\n", encoding="utf-8")
    _validate_blacklist_bed(str(bed), ref_contigs={"chr1": 200}, genome="hg38")

    # Missing file
    with pytest.raises(ValidationError, match="blacklist file not found"):
        _validate_blacklist_bed(str(tmp_path / "missing.bed"))

    # Fewer than 3 fields
    bad_fields = tmp_path / "bad_fields.bed"
    bad_fields.write_text("chr1\t10\n", encoding="utf-8")
    with pytest.raises(ValidationError, match="fewer than 3 fields"):
        _validate_blacklist_bed(str(bad_fields))

    # Non-integer
    bad_int = tmp_path / "bad_int.bed"
    bad_int.write_text("chr1\tfoo\t100\n", encoding="utf-8")
    with pytest.raises(ValidationError, match="not integers"):
        _validate_blacklist_bed(str(bad_int))

    # Negative start
    bad_neg = tmp_path / "bad_neg.bed"
    bad_neg.write_text("chr1\t-1\t100\n", encoding="utf-8")
    with pytest.raises(ValidationError, match="negative"):
        _validate_blacklist_bed(str(bad_neg))

    # start >= end
    bad_span = tmp_path / "bad_span.bed"
    bad_span.write_text("chr1\t100\t50\n", encoding="utf-8")
    with pytest.raises(ValidationError, match="start=100 >= end=50"):
        _validate_blacklist_bed(str(bad_span))

    # Contig absent
    with pytest.raises(ValidationError, match="absent from reference"):
        _validate_blacklist_bed(str(bed), ref_contigs={"chr2": 200}, genome="hg38")

    # Exceeds contig length
    with pytest.raises(ValidationError, match="exceeds contig length"):
        _validate_blacklist_bed(str(bed), ref_contigs={"chr1": 80}, genome="hg38")

    # Empty BED
    empty_bed = tmp_path / "empty.bed"
    empty_bed.write_text("# comment only\n", encoding="utf-8")
    with pytest.raises(ValidationError, match="empty or contains no valid"):
        _validate_blacklist_bed(str(empty_bed))


def test_reference_contig_validation_and_fai(tmp_path: Path):
    # Matching sizes and fai
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t1000\nchr2\t500\n", encoding="utf-8")
    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\nACGT\n>chr2\nACGT\n", encoding="utf-8")
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t1000\t6\t50\t51\nchr2\t500\t100\t50\t51\n", encoding="utf-8")

    resource = {"chrom_sizes": str(sizes), "reference_fasta": str(fa)}
    contigs = _get_reference_contigs(resource, genome="hg38")
    assert contigs == {"chr1": 1000, "chr2": 500}

    # Length mismatch between FAI and chrom_sizes
    fai.write_text("chr1\t9999\t6\t50\t51\nchr2\t500\t100\t50\t51\n", encoding="utf-8")
    with pytest.raises(ValidationError, match="length mismatch"):
        _get_reference_contigs(resource, genome="hg38")

    # Empty fai
    empty_fai = tmp_path / "empty.fa.fai"
    empty_fai.write_text("", encoding="utf-8")
    with pytest.raises(ValidationError, match="empty or invalid"):
        _read_fai(str(empty_fai))


def test_admission_blacklist_mode_off_vs_analysis(tmp_path: Path):
    # Sample
    samples = [_sample(id="S1", role="treatment", genome="hg38")]
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t1000\n", encoding="utf-8")

    # Invalid blacklist path and invalid content
    nonexistent_path = str(tmp_path / "does_not_exist.bed")
    invalid_bed = tmp_path / "invalid.bed"
    invalid_bed.write_text("chr1\t-5\t100\n", encoding="utf-8")

    # 1. blacklist_mode: off passes with nonexistent blacklist
    cfg_off_nonexistent = {
        "use_control": False,
        "blacklist_mode": "off",
        "genome_resources": {
            "hg38": {
                "chrom_sizes": str(sizes),
                "blacklist": nonexistent_path,
            }
        },
    }
    _admit(samples, cfg_off_nonexistent)

    # 2. blacklist_mode: off passes with invalid BED content
    cfg_off_invalid = {
        "use_control": False,
        "blacklist_mode": "off",
        "genome_resources": {
            "hg38": {
                "chrom_sizes": str(sizes),
                "blacklist": str(invalid_bed),
            }
        },
    }
    _admit(samples, cfg_off_invalid)

    # 3. blacklist_mode: analysis fails with nonexistent blacklist
    cfg_analysis_nonexistent = {
        "use_control": False,
        "blacklist_mode": "analysis",
        "genome_resources": {
            "hg38": {
                "chrom_sizes": str(sizes),
                "blacklist": nonexistent_path,
            }
        },
    }
    with pytest.raises(ValidationError, match="blacklist file not found"):
        _admit(samples, cfg_analysis_nonexistent)

    # 4. blacklist_mode: analysis fails with invalid BED content
    cfg_analysis_invalid = {
        "use_control": False,
        "blacklist_mode": "analysis",
        "genome_resources": {
            "hg38": {
                "chrom_sizes": str(sizes),
                "blacklist": str(invalid_bed),
            }
        },
    }
    with pytest.raises(ValidationError, match="negative"):
        _admit(samples, cfg_analysis_invalid)


def test_admission_organelle_validation(tmp_path: Path):
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t1000\n", encoding="utf-8")
    samples = [_sample(id="S1", role="treatment", genome="hg38", assay="cuttag")]

    # Explicit organelle contig missing from reference
    cfg_explicit_missing = {
        "use_control": False,
        "blacklist_mode": "off",
        "organelle_mode": "auto",
        "genome_resources": {
            "hg38": {
                "chrom_sizes": str(sizes),
                "organelle_chroms": ["chrM_missing"],
            }
        },
    }
    with pytest.raises(ValidationError, match="explicit organelle contig"):
        _admit(samples, cfg_explicit_missing)

    # Auto mode for cuttag resolves to exclude, and none of chrM/MT/chrMT exist in ref_contigs
    cfg_exclude_missing = {
        "use_control": False,
        "blacklist_mode": "off",
        "organelle_mode": "auto",  # cuttag -> exclude
        "genome_resources": {
            "hg38": {
                "chrom_sizes": str(sizes),
            }
        },
    }
    with pytest.raises(ValidationError, match="none of candidate contigs"):
        _admit(samples, cfg_exclude_missing)


def test_control_organelle_mode_uses_treatment_context(tmp_path: Path):
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t1000\n", encoding="utf-8")
    samples = [
        _sample(
            id="T1",
            role="treatment",
            assay="chipseq",
            control_sample="C1",
        ),
        _sample(id="C1", role="control", assay="atac"),
    ]
    cfg = {
        "use_control": True,
        "blacklist_mode": "off",
        "organelle_mode": "auto",
        "genome_resources": {"hg38": {"chrom_sizes": str(sizes)}},
    }
    _admit(samples, cfg)


def test_shared_control_with_same_resolved_organelle_policy(tmp_path: Path):
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t1000\n", encoding="utf-8")
    samples = [
        _sample(id="T1", role="treatment", assay="chipseq", control_sample="C1"),
        _sample(
            id="T2",
            role="treatment",
            assay="mnase",
            experiment="EXP2",
            control_sample="C1",
        ),
        _sample(id="C1", role="control", assay="chipseq"),
    ]
    cfg = {
        "use_control": True,
        "blacklist_mode": "off",
        "organelle_mode": "auto",
        "genome_resources": {"hg38": {"chrom_sizes": str(sizes)}},
    }
    _admit(samples, cfg, reference_mode="legacy")


def test_copies_mode_requires_blacklist_resource(tmp_path: Path):
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t1000\n", encoding="utf-8")
    samples = [_sample(id="S1", role="treatment", genome="hg38")]
    cfg = {
        "use_control": False,
        "blacklist_mode": "copies",
        "genome_resources": {"hg38": {"chrom_sizes": str(sizes), "blacklist": ""}},
    }
    with pytest.raises(ValidationError, match="blacklist_mode=copies"):
        _admit(samples, cfg)


def test_reference_contigs_union_preserves_fai_only_contigs(tmp_path: Path):
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t1000\n", encoding="utf-8")
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n>chrM\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text(
        "chr1\t1000\t0\t50\t51\nchrM\t200\t0\t50\t51\n", encoding="utf-8"
    )
    contigs = _get_reference_contigs(
        {"chrom_sizes": str(sizes), "reference_fasta": str(fasta)}, genome="hg38"
    )
    assert contigs == {"chr1": 1000, "chrM": 200}


def test_reference_contigs_disjoint_sources_fail(tmp_path: Path):
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t1000\n", encoding="utf-8")
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr2\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr2\t500\t0\t50\t51\n", encoding="utf-8")
    with pytest.raises(ValidationError, match="no contigs in common"):
        _get_reference_contigs(
            {"chrom_sizes": str(sizes), "reference_fasta": str(fasta)},
            genome="hg38",
        )


def test_default_blacklist_mode_is_analysis_and_requires_resource(tmp_path: Path):
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t1000\n", encoding="utf-8")
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("chr1\t10\t20\n", encoding="utf-8")
    samples = [_sample(id="S1", role="treatment", genome="hg38")]

    missing = {
        "use_control": False,
        "genome_resources": {"hg38": {"chrom_sizes": str(sizes), "blacklist": ""}},
    }
    with pytest.raises(ValidationError, match="blacklist_mode=analysis"):
        _admit(samples, missing)

    present = {
        "use_control": False,
        "genome_resources": {
            "hg38": {"chrom_sizes": str(sizes), "blacklist": str(blacklist)}
        },
    }
    _admit(samples, present)
