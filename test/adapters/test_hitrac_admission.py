"""Fast admission contracts; no local runtime or scientific tool is required."""

import gzip
import json
from pathlib import Path

import pytest

from encode_pipeline.adapters.hitrac_preprocess.admission import (
    INDEX_SUFFIXES,
    AdmissionError,
    ReferenceBinding,
    RuntimeBinding,
    Sample,
    load_reference_binding,
    load_runtime_binding,
    prepare_inputs,
    sha256_file,
    validate_fastq_pair,
)


def write_fastq(path, records):
    with gzip.open(path, "wt", encoding="ascii") as handle:
        handle.write(records)
    return path


def fq(name="id", sequence="ACGT", quality="IIII", plus="+"):
    return f"@{name}\n{sequence}\n{plus}\n{quality}\n"


@pytest.mark.parametrize(
    ("left", "right", "count"),
    [
        (fq(), fq(), 1),
        (fq("id/1"), fq("id/2"), 1),
        (fq("id 1:N:0:ACGT"), fq("id 2:N:0:ACGT"), 1),
        (fq("1_120_-1"), fq("1_120_-1"), 1),
        (fq("id") * 2, fq("id") * 2, 2),
        (fq(sequence="NRYA"), fq(sequence="nrya"), 1),
        ("", "", 0),
    ],
)
def test_fastq_full_pair_contract(tmp_path, left, right, count):
    paths = [
        write_fastq(tmp_path / f"R{mate}.gz", text)
        for mate, text in ((1, left), (2, right))
    ]
    before = [path.read_bytes() for path in paths]
    assert validate_fastq_pair(*paths) == count
    assert [path.read_bytes() for path in paths] == before


@pytest.mark.parametrize(
    ("left", "right", "reason"),
    [
        (fq(), "", "fastq_pair_count_mismatch"),
        ("", fq(), "fastq_pair_count_mismatch"),
        (fq("x"), fq("y"), "fastq_pair_id_mismatch"),
        (fq("x/2"), fq("x/2"), "fastq_mate_role_invalid"),
        (fq("x/1"), fq("x/1"), "fastq_mate_role_invalid"),
        (fq("x 2:N:0:AT"), fq("x 2:N:0:AT"), "fastq_mate_role_invalid"),
        (fq("x 3:N:0:AT"), fq("x 2:N:0:AT"), "fastq_mate_role_invalid"),
        ("@x\nACTG\n+\n", fq(), "fastq_truncated_record"),
        (fq(sequence="ACT"), fq(), "fastq_sequence_quality_invalid"),
        (fq(sequence="AC?T"), fq(), "fastq_sequence_quality_invalid"),
        (fq(quality="II I"), fq(), "fastq_quality_invalid"),
        (fq(plus="-"), fq(), "fastq_structure_invalid"),
        (fq(plus="+other"), fq(), "fastq_plus_header_mismatch"),
    ],
)
def test_fastq_malformed_or_mismatched_pair_is_rejected(tmp_path, left, right, reason):
    a = write_fastq(tmp_path / "R1.gz", left)
    b = write_fastq(tmp_path / "R2.gz", right)
    with pytest.raises(AdmissionError, match=f"^{reason}$") as caught:
        validate_fastq_pair(a, b)
    assert str(tmp_path) not in str(caught.value)


@pytest.mark.parametrize(
    "content", [b"not-gzip", b"", gzip.compress(fq().encode())[:-5]]
)
def test_bad_gzip_and_trailer_are_rejected(tmp_path, content):
    a = tmp_path / "R1.gz"
    a.write_bytes(content)
    b = write_fastq(tmp_path / "R2.gz", fq())
    with pytest.raises(AdmissionError, match="^fastq_gzip_invalid$"):
        validate_fastq_pair(a, b)


def reference_fixture(tmp_path):
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chrA\nACGT\n>chrB\nTGCAN\n")
    prefix = tmp_path / "index"
    files = {"fasta": sha256_file(fasta)}
    for suffix in INDEX_SUFFIXES:
        path = Path(f"{prefix}{suffix}")
        path.write_bytes(f"synthetic-index-{suffix}".encode())
        files[suffix] = sha256_file(path)
    manifest = tmp_path / "reference.json"
    manifest.write_text(
        json.dumps(
            {
                "schema_version": "hitrac-reference-binding-v1",
                "fasta": str(fasta),
                "prefix": str(prefix),
                "files": files,
                "contigs": {"chrA": 4, "chrB": 5},
            }
        )
    )
    return manifest


def test_reference_manifest_is_explicitly_trusted_and_files_rechecked(tmp_path):
    manifest = reference_fixture(tmp_path)
    approved = sha256_file(manifest)
    reference = load_reference_binding(manifest, approved)
    assert reference.contigs == {"chrA": 4, "chrB": 5}
    altered = Path(f"{reference.prefix}.rev.2.bt2")
    altered.write_bytes(b"other-assembly-index")
    with pytest.raises(AdmissionError, match="^reference_file_identity_mismatch$"):
        load_reference_binding(manifest, approved)
    payload = json.loads(manifest.read_text())
    payload["files"][".rev.2.bt2"] = sha256_file(altered)
    manifest.write_text(json.dumps(payload))
    with pytest.raises(AdmissionError, match="^reference_binding_identity_mismatch$"):
        load_reference_binding(manifest, approved)


@pytest.mark.parametrize(
    "kind", ["missing", "empty", "symlink", "fasta", "contigs", "extra", "bt2l"]
)
def test_reference_rejections(tmp_path, kind):
    manifest = reference_fixture(tmp_path)
    payload = json.loads(manifest.read_text())
    index = Path(f"{payload['prefix']}.2.bt2")
    if kind == "missing":
        index.unlink()
        expected = "input_missing"
    elif kind == "empty":
        index.write_bytes(b"")
        payload["files"][".2.bt2"] = sha256_file(index)
        expected = "reference_file_empty"
    elif kind == "symlink":
        index.unlink()
        index.symlink_to(Path(f"{payload['prefix']}.1.bt2"))
        expected = "input_symlink"
    elif kind == "fasta":
        Path(payload["fasta"]).write_text(">chrA\nACGT\n>chrA\nTGCA\n")
        payload["files"]["fasta"] = sha256_file(Path(payload["fasta"]))
        expected = "reference_fasta_invalid"
    elif kind == "contigs":
        payload["contigs"]["chrA"] = 400
        expected = "reference_contigs_mismatch"
    elif kind == "extra":
        payload["files"]["extra"] = "0" * 64
        expected = "reference_index_set_invalid"
    else:
        payload["files"][".2.bt2l"] = payload["files"].pop(".2.bt2")
        expected = "reference_index_set_invalid"
    manifest.write_text(json.dumps(payload))
    with pytest.raises(AdmissionError, match=f"^{expected}$"):
        load_reference_binding(manifest, sha256_file(manifest))


def runtime_stub(tmp_path):
    return RuntimeBinding(
        tmp_path, tmp_path / "python", tmp_path / "tracPre2.py", {}, "a" * 64, "b" * 64
    )


def test_staging_is_private_tokenized_and_source_read_only(tmp_path):
    manifest = reference_fixture(tmp_path)
    reference = load_reference_binding(manifest, sha256_file(manifest))
    folder = tmp_path / "raw inputs with spaces"
    folder.mkdir()
    left, right = (write_fastq(folder / f"r{i}.gz", fq(f"id/{i}")) for i in (1, 2))
    before = [
        (p.read_bytes(), p.stat().st_mode, p.stat().st_mtime_ns) for p in (left, right)
    ]
    attempt = tmp_path / "attempt"
    attempt.mkdir(mode=0o700)
    result = prepare_inputs(
        [Sample("Display Sample", left, right)],
        runtime_stub(tmp_path),
        reference,
        attempt,
    )
    assert result.samples["s000001"]["raw_pairs"] == 1
    assert result.samples["s000001"]["display_id"] == "Display Sample"
    for i, source in ((1, left), (2, right)):
        staged = result.staged_fastq_dir / f"s000001_R{i}.fastq.gz"
        assert staged.read_bytes() == source.read_bytes()
        assert staged.stat().st_mode & 0o777 == 0o400
        assert not staged.is_symlink()
        assert staged.stat().st_ino != source.stat().st_ino
    assert [
        (p.read_bytes(), p.stat().st_mode, p.stat().st_mtime_ns) for p in (left, right)
    ] == before
    with pytest.raises(AdmissionError, match="^attempt_staging_exists$"):
        prepare_inputs(
            [Sample("Display Sample", left, right)],
            runtime_stub(tmp_path),
            reference,
            attempt,
        )
    assert (result.staged_fastq_dir / "s000001_R1.fastq.gz").exists()


@pytest.mark.parametrize(
    "name", ["../sample", "x;touch", "x\ny", "", "-option", "x/child"]
)
def test_unsafe_sample_names_rejected_before_staging(tmp_path, name):
    attempt = tmp_path / "attempt"
    attempt.mkdir(mode=0o700)
    reference = ReferenceBinding(
        tmp_path / "unused", tmp_path / "unused", {}, {}, "c" * 64
    )
    with pytest.raises(AdmissionError, match="^sample_id_invalid$"):
        prepare_inputs(
            [Sample(name, tmp_path / "absent1", tmp_path / "absent2")],
            runtime_stub(tmp_path),
            reference,
            attempt,
        )
    assert not (attempt / "input").exists()


def test_failed_input_remains_private_and_identifies_only_token(tmp_path):
    manifest = reference_fixture(tmp_path)
    reference = load_reference_binding(manifest, sha256_file(manifest))
    a = write_fastq(tmp_path / "a.gz", fq("a"))
    b = write_fastq(tmp_path / "b.gz", fq("b"))
    attempt = tmp_path / "attempt"
    attempt.mkdir(mode=0o700)
    with pytest.raises(AdmissionError, match="^fastq_pair_id_mismatch$") as caught:
        prepare_inputs(
            [Sample("private sample", a, b)], runtime_stub(tmp_path), reference, attempt
        )
    assert caught.value.sample == "s000001"
    assert (attempt / "input/s000001_R1.fastq.gz").read_bytes() == a.read_bytes()
    assert not (attempt / "completed.json").exists()


def test_reference_drift_between_admission_and_staging_fails(tmp_path):
    manifest = reference_fixture(tmp_path)
    reference = load_reference_binding(manifest, sha256_file(manifest))
    a = write_fastq(tmp_path / "a.gz", fq())
    b = write_fastq(tmp_path / "b.gz", fq())
    reference.fasta.write_text(">chrA\nAAAA\n")
    attempt = tmp_path / "attempt"
    attempt.mkdir(mode=0o700)
    with pytest.raises(AdmissionError, match="^reference_changed_during_staging$"):
        prepare_inputs(
            [Sample("sample", a, b)], runtime_stub(tmp_path), reference, attempt
        )


def test_runtime_binding_cannot_supply_its_own_lock(tmp_path):
    binding = tmp_path / "runtime.json"
    binding.write_text(
        json.dumps(
            {
                "schema_version": "hitrac-runtime-binding-v1",
                "prefix": str(tmp_path),
                "script": str(tmp_path / "untrusted.py"),
                "lock_sha256": "0" * 64,
            }
        )
    )
    with pytest.raises(AdmissionError, match="^runtime_lock_mismatch$"):
        load_runtime_binding(binding)


def test_concatenated_gzip_members_form_one_sample_pair(tmp_path):
    first = gzip.compress(fq("a/1").encode()) + gzip.compress(fq("b/1").encode())
    second = gzip.compress(fq("a/2").encode()) + gzip.compress(fq("b/2").encode())
    a, b = tmp_path / "a.gz", tmp_path / "b.gz"
    a.write_bytes(first)
    b.write_bytes(second)
    assert validate_fastq_pair(a, b) == 2


def test_missing_mate_and_symlink_input_reject(tmp_path):
    a = write_fastq(tmp_path / "a.gz", fq())
    with pytest.raises(AdmissionError, match="^input_missing$"):
        validate_fastq_pair(a, tmp_path / "absent.gz")
    link = tmp_path / "link.gz"
    link.symlink_to(a)
    with pytest.raises(AdmissionError, match="^input_symlink$"):
        validate_fastq_pair(link, a)
