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
    return {"genome_resources": resources, "use_control": use_control}


def _admit(samples, config, *, reference_mode="strict", external="validate"):
    validate_reference_control_admission(
        samples,
        config,
        reference_mode=reference_mode,
        external_control_policy=external,
        error_cls=ValidationError,
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
