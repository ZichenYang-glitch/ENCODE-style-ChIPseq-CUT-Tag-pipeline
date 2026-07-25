from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
DECISION_PATH = (
    REPOSITORY_ROOT / "docs" / "architecture" / "data-foundation-decision.md"
)


def test_data_foundation_decision_records_required_compatibility_contracts() -> None:
    decision = DECISION_PATH.read_text(encoding="utf-8")

    required_contracts = (
        "# Data Foundation Architecture and Compatibility Decision",
        "Project",
        "Sample",
        "SampleRevision",
        "InputFile",
        "StoragePool",
        "RunSample",
        "ArtifactSample",
        "SQLite is authoritative",
        "JSON manifests",
        "stable internal IDs",
        "created_at",
        "started_at",
        "ended_at",
        "produced_at",
        "published_at",
        "archived_at",
        "reserved Legacy Project",
        "must not infer sample provenance",
        "snapshot binding envelope",
        "unresolved",
        "administrator-only local CLI",
        "must not add an unauthenticated management write API",
        "manifest-last",
        "cross-filesystem",
    )

    missing = [contract for contract in required_contracts if contract not in decision]
    assert not missing, f"decision is missing required contracts: {missing}"


def test_data_foundation_decision_pins_normative_compatibility_matrix() -> None:
    decision = DECISION_PATH.read_text(encoding="utf-8")

    required_rows = (
        "| new HTTP server path after managed-input cutover | forbidden |",
        "| pre-existing snapshot | Legacy + unresolved |",
        "| snapshot consumption | one SQLite transaction |",
        "| snapshot replay | exact ordered binding equality or fail closed |",
        "| managed-input execution | ephemeral adapter view; resolved path not persisted |",
        "| publication move unit | one regular-file Artifact revision |",
        "| publication visibility | valid manifest AND indexed SQLite state |",
        "| publication content identity | separate from existing artifact revision |",
        "| bound sample-scoped artifact | one or more snapshot contributors required |",
        "| sample query date axis | published_at UTC [from, to) |",
        "| cumulative migration evidence | rev08 to exact stage head |",
    )

    missing = [row for row in required_rows if row not in decision]
    assert not missing, f"decision is missing normative matrix rows: {missing}"


def test_data_foundation_decision_rejects_legacy_path_permission() -> None:
    decision = DECISION_PATH.read_text(encoding="utf-8")

    assert (
        "`legacy_v1` describes missing provenance and never grants path-submission"
        in decision
    )
    assert (
        "Only a historical Legacy artifact\nwithout evidence uses "
        "`unresolved_legacy`" in decision
    )
    assert (
        "Every pre-migration run and every pre-migration validation snapshot"
        in decision
    )


def test_architecture_overview_links_the_data_foundation_decision() -> None:
    overview = (
        REPOSITORY_ROOT / "docs" / "architecture" / "platform-overview.md"
    ).read_text(encoding="utf-8")

    assert "(data-foundation-decision.md)" in overview
