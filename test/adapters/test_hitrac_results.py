"""Fast original-summary contracts; scientific correspondence is covered by H2
and explicit real result consumer tests, not asserted from synthetic TSV alone.
"""

from decimal import Decimal

import pytest

from encode_pipeline.adapters.hitrac_preprocess.results import (
    HiTracPreprocessResultsAdapter,
    METRIC_KEYS,
    parse_summary,
)
from encode_pipeline.adapters.hitrac_preprocess.outputs import summary_columns
from encode_pipeline.platform.adapters import QcSummaryExtractingAdapter

VALUES = [
    12,
    10,
    0.8,
    8,
    0.125,
    6 / 7,
    0.5,
    1 / 3,
    1 / 6,
    5,
    5 / 12,
    0.8,
    0.25,
    0.5,
    0.25,
]


def summary(mapq=10, tokens=("s000001", "s000002")):
    return (
        "\t"
        + "\t".join(summary_columns(mapq))
        + "\n"
        + "".join(token + "\t" + "\t".join(map(str, VALUES)) + "\n" for token in tokens)
    ).encode()


@pytest.mark.parametrize("mapq", [0, 10, 30, 255])
def test_summary_original_sixteen_columns_and_all_metrics(mapq):
    samples = {
        token: {"raw_pairs": 12, "display_id": name}
        for token, name in (("s000001", "first"), ("s000002", "second"))
    }
    result = parse_summary(summary(mapq), samples, mapq)
    assert list(result) == list(samples)
    assert len(METRIC_KEYS) == 15
    assert result["s000001"] == list(map(lambda v: Decimal(str(v)), VALUES))
    assert result["s000002"] == result["s000001"]


@pytest.mark.parametrize(
    "tokens",
    [
        ("s000001",),
        ("s000001", "s000001"),
        ("s000001", "s000002", "s000003"),
        ("wrong", "s000002"),
    ],
)
def test_summary_rejects_missing_duplicate_extra_wrong_samples(tokens):
    with pytest.raises(ValueError):
        parse_summary(
            summary(tokens=tokens),
            {"s000001": {"raw_pairs": 12}, "s000002": {"raw_pairs": 12}},
            10,
        )


@pytest.mark.parametrize(
    "index,value",
    [
        (0, "11"),
        (0, "1.5"),
        (1, "-1"),
        (3, "0"),
        (9, "0"),
        (5, "0"),
        (11, "0"),
        (2, "NaN"),
        (2, "inf"),
        (2, "abc"),
        (2, "1.1"),
        (10, "0.9"),
        (6, "0.9"),
    ],
)
def test_summary_rejects_illegal_metrics(index, value):
    lines = summary().decode().splitlines()
    fields = lines[1].split("\t")
    fields[index + 1] = value
    lines[1] = "\t".join(fields)
    from decimal import InvalidOperation

    with pytest.raises((ValueError, InvalidOperation)):
        parse_summary(
            ("\n".join(lines) + "\n").encode(),
            {"s000001": {"raw_pairs": 12}, "s000002": {"raw_pairs": 12}},
            10,
        )


@pytest.mark.parametrize("change", ["header", "width", "encoding"])
def test_summary_rejects_contract_mismatch(change):
    content = summary()
    if change == "header":
        content = content.replace(b"mapq>=10", b"mapq>=20")
    elif change == "width":
        content = content.replace(b"s000001\t", b"s000001\textra\t")
    else:
        content += b"\xff"
    with pytest.raises(ValueError):
        parse_summary(
            content, {"s000001": {"raw_pairs": 12}, "s000002": {"raw_pairs": 12}}, 10
        )


def test_results_capabilities_require_real_binding_and_atomic_publication():
    adapter = HiTracPreprocessResultsAdapter()
    assert isinstance(adapter, QcSummaryExtractingAdapter)
    assert adapter.requires_atomic_result_publication() is True
    assert adapter.execution_availability().execution == "not_configured"
    assert adapter.qc_source_output_types() == ("hitrac_summary",)
    assert adapter.extract_artifacts(None, None).is_failure
    assert adapter.extract_qc_metrics(None, ()).is_failure


@pytest.fixture
def results_case(bound_case):
    case = bound_case
    adapter = HiTracPreprocessResultsAdapter(
        runtime=case["runtime"], binding=case["bound"]._binding
    )
    return case, adapter


# The same admitted-runtime stub as H3; no scientific result is manufactured.
from test_hitrac_adapter import bound_case as shared_bound_case  # noqa: E402

bound_case = shared_bound_case


def qc_document(content, **changes):
    from encode_pipeline.platform.adapters import QcSourceArtifact, QcSourceDocument

    fields = dict(
        artifact_id="source",
        output_type="hitrac_summary",
        relative_path="results/hitrac/" + "a" * 64 + "/tracPre_summary.txt",
        metadata={"scope": "run"},
    )
    fields.update(changes)
    return QcSourceDocument(source=QcSourceArtifact(**fields), content=content)


def small_summary():
    values = [1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0]
    return (
        "\t"
        + "\t".join(summary_columns(10))
        + "\n"
        + "".join(
            token + "\t" + "\t".join(map(str, values)) + "\n"
            for token in ("s000001", "s000002")
        )
    ).encode()


def test_qc_consumer_maps_frozen_order_and_units_without_flagging_science(results_case):
    case, adapter = results_case
    result = adapter.extract_qc_metrics(
        case["bound"]._binding.verify()[2],
        (qc_document(small_summary()),),
    )
    assert result.is_success, result.issues
    assert len(result.value) == 30
    assert [metric.sample_id for metric in result.value] == ["zeta"] * 15 + [
        "alpha"
    ] * 15
    assert [metric.unit for metric in result.value[:15]] == [
        "count",
        "count",
        "fraction",
        "count",
        *["fraction"] * 5,
        "count",
        "ratio",
        *["fraction"] * 4,
    ]
    assert all(
        metric.qc_flag is None and metric.source_artifact_id == "source"
        for metric in result.value
    )


@pytest.mark.parametrize(
    "change",
    [
        {"output_type": "private_bam"},
        {"relative_path": "hitrac-attempt/private/request.json"},
        {"relative_path": "results/hitrac/old/tracPre_summary.txt"},
        {"metadata": {"scope": "sample", "sample_id": "wrong"}},
    ],
)
def test_qc_source_contract_rejections_do_not_leak_values(results_case, change):
    case, adapter = results_case
    result = adapter.extract_qc_metrics(
        case["bound"]._binding.verify()[2], (qc_document(small_summary(), **change),)
    )
    assert result.is_failure
    assert {issue.code for issue in result.issues} == {"HITRAC_QC_INVALID"}
    assert "private_bam" not in repr(result.issues)
    assert "request.json" not in repr(result.issues)


def test_projection_existing_mismatch_and_symlinks_fail_closed(tmp_path):
    from encode_pipeline.adapters.hitrac_preprocess.results import _project
    from hashlib import sha256

    output = tmp_path / "output"
    output.mkdir()
    original = output / "tracPre_summary.txt"
    original.write_bytes(b"original-bytes\n")
    record = {
        "path": original.name,
        "size_bytes": original.stat().st_size,
        "sha256": sha256(original.read_bytes()).hexdigest(),
    }
    result = _project(tmp_path, output, {}, {"summary": record}, "a" * 64)
    projected = tmp_path / result[0].relative_path
    assert projected.read_bytes() == original.read_bytes()
    # Only this task-created synthetic projection is altered by the test.
    projected.chmod(0o600)
    projected.write_bytes(b"different-bytes\n")
    with pytest.raises(ValueError):
        _project(tmp_path, output, {}, {"summary": record}, "a" * 64)
    projected.unlink()
    projected.symlink_to(original)
    with pytest.raises(ValueError):
        _project(tmp_path, output, {}, {"summary": record}, "a" * 64)
    assert original.read_bytes() == b"original-bytes\n"


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("0.7142857142857143", "0.714285714286"),
        ("0.1234567890125", "0.123456789012"),
        ("0.1234567890135", "0.123456789014"),
    ],
)
def test_public_fraction_uses_existing_twelve_decimal_contract(
    results_case, raw, expected
):
    from encode_pipeline.services.qc_summary_indexing import QcSummaryIndexingService

    case, adapter = results_case
    original = small_summary()
    lines = original.decode().splitlines()
    for index in (1, 2):
        values = lines[index].split("\t")
        values[3] = raw
        lines[index] = "\t".join(values)
    content = ("\n".join(lines) + "\n").encode()
    inputs = case["bound"]._binding.verify()[2]
    samples = case["bound"]._binding.verify()[3]
    assert parse_summary(content, samples, 10)["s000001"][2] == Decimal(raw)
    source = qc_document(content)
    result = adapter.extract_qc_metrics(inputs, (source,))
    assert result.is_success, result.issues
    mapped = [metric for metric in result.value if metric.metric_key == "mapping_ratio"]
    assert len(mapped) == 2
    for metric in mapped:
        assert metric.value == Decimal(expected)
        assert abs(metric.value - Decimal(raw)) <= Decimal("0.0000000000005")
        QcSummaryIndexingService._validate_candidate(metric, {"source"})
    assert source.content == content
    assert all(
        metric.value == 1 for metric in result.value if metric.metric_key == "raw_pairs"
    )
