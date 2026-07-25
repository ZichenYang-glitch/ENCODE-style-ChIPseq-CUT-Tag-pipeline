import re
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
DECISION_PATH = (
    REPOSITORY_ROOT / "docs" / "architecture" / "data-foundation-decision.md"
)


def _section(decision: str, heading: str, *, level: int = 3) -> str:
    prefix = "#" * level
    marker = f"{prefix} {heading}\n"
    start = decision.index(marker) + len(marker)
    following = re.search(rf"\n#{{1,{level}}} ", decision[start:])
    return (
        decision[start:]
        if following is None
        else decision[start : start + following.start()]
    )


def _assert_contains(
    section_name: str, section: str, contracts: tuple[str, ...]
) -> None:
    normalized = re.sub(r"\s+", " ", section)
    missing = [
        contract
        for contract in contracts
        if re.sub(r"\s+", " ", contract) not in normalized
    ]
    assert not missing, f"{section_name} is missing contracts: {missing}"


def _unconditional_global_path_rejections(text: str) -> list[str]:
    normalized = re.sub(r"\s+", " ", text)
    claims: list[str] = []

    known_global_claims = (
        r"StoragePool/InputFile stage is the mandatory cutover",
        r"from that stage onward, new HTTP and frontend validation rejects "
        r"physical server paths",
        r"whether or not the caller supplies a binding envelope",
        r"new HTTP server path after managed-input cutover",
    )
    claims.extend(
        pattern
        for pattern in known_global_claims
        if re.search(pattern, normalized, flags=re.IGNORECASE)
    )

    path_policy = re.compile(
        r"\b(?:reject\w*|forbid\w*|prohibit\w*|disallow\w*|den(?:y|ies|ied)|"
        r"unavailable)\b|\bcannot\s+(?:submit|accept|use)\b",
        flags=re.IGNORECASE,
    )
    path_subject = re.compile(
        r"\b(?:input|physical|server|frontend)\b[^.!?]{0,40}\bpaths?\b",
        flags=re.IGNORECASE,
    )
    universal_scope = re.compile(
        r"\b(?:all|every|any|global|stage-wide|milestone-wide)\b|"
        r"\b(?:registry|StoragePool/InputFile|persistence|Project/Sample)\s+stage\b",
        flags=re.IGNORECASE,
    )
    qualified_scope = re.compile(
        r"\bqualified input use\b|"
        r"\binput use\b[^.!?]{0,60}\bcompleted managed cutover\b|"
        r"\bonly that qualified\b",
        flags=re.IGNORECASE,
    )
    negated_policy = re.compile(
        r"\b(?:does not|do not|must not|never|no unconditional)\b[^.!?]{0,80}"
        r"\b(?:reject|forbid|prohibit|cutover)\w*\b",
        flags=re.IGNORECASE,
    )

    for sentence in re.split(r"(?<=[.!?])\s+", normalized):
        if (
            path_policy.search(sentence)
            and path_subject.search(sentence)
            and universal_scope.search(sentence)
            and not qualified_scope.search(sentence)
            and not negated_policy.search(sentence)
        ):
            claims.append(sentence)
    return claims


def _compatibility_matrix(decision: str) -> dict[str, str]:
    section = _section(decision, "Normative compatibility matrix")
    rows: dict[str, str] = {}
    for line in section.splitlines():
        if (
            not line.startswith("| ")
            or line.startswith("| Situation ")
            or "---" in line
        ):
            continue
        cells = [cell.strip() for cell in line.strip("|").split("|")]
        assert len(cells) == 2, f"malformed compatibility row: {line}"
        situation, contract = cells
        assert situation not in rows, f"duplicate compatibility row: {situation}"
        rows[situation] = contract
    return rows


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

    rows = _compatibility_matrix(decision)
    expected_rows = {
        "managed cutover identity": (
            "workflow ID + adapter contract version + input use capability"
        ),
        "qualified managed input use": ("opaque revision only; physical path rejected"),
        "new catalog-aware unqualified input use": (
            "explicit transitional_unmanaged_v1; trusted-local path contract preserved"
        ),
        "mixed managed/unmanaged snapshot": (
            "per-use provenance required and exact adapter contract must allow the mix"
        ),
        "fully managed adapter": "every execution-required input use qualified",
        "pre-migration or historical unbound v1 snapshot": "Legacy + unresolved",
        "snapshot consumption": "one SQLite transaction",
        "snapshot replay": (
            "exact ordered per-use mode and evidence equality or fail closed"
        ),
        "managed-input execution": (
            "ephemeral adapter view; resolved path not persisted"
        ),
        "publication move unit": "one regular-file Artifact revision",
        "publication visibility": "valid manifest AND indexed SQLite state",
        "publication content identity": "separate from existing artifact revision",
        "bound sample-scoped artifact": ("one or more snapshot contributors required"),
        "sample query date axis": "published_at UTC [from, to)",
        "cumulative migration evidence": "rev08 to exact stage head",
    }

    assert rows == expected_rows


def test_managed_input_cutover_is_per_adapter_per_use_and_never_global() -> None:
    decision = DECISION_PATH.read_text(encoding="utf-8")
    binding = _section(
        decision,
        "Snapshot-time project, sample, and input binding",
    )
    legacy = _section(
        decision,
        "Legacy compatibility, capability-scoped transition, and conservative migration",
    )
    storage = _section(decision, "Storage pools and input revisions")
    consequences = _section(
        decision,
        "Consequences and non-goals",
        level=2,
    )

    _assert_contains(
        "snapshot binding",
        binding,
        (
            "workflow ID and adapter contract version",
            "stable opaque key",
            "exact input-use capability and closure-contract versions",
            "provenance mode",
            "ordered zero-or-more opaque resource/revision IDs",
            "logical member coordinates",
            "generic closure digest",
            "`managed_revision_v1`",
            "`transitional_unmanaged_v1`",
            "without copying or interpreting their physical paths",
            "complete ordered per-input-use mode",
            "mixed managed/unmanaged snapshot",
            "every execution-required input use must appear exactly once",
            "Missing, duplicated, undeclared, or contract-forbidden combinations "
            "fail closed",
            "fully managed only when every input use required by that execution",
        ),
    )
    _assert_contains(
        "legacy transition",
        legacy,
        (
            "Managed cutover is per-adapter, per-input-use",
            "workflow ID, adapter contract version, and input-use capability",
            "complete closure can be represented by catalog revisions",
            "adapter can validate the scientific layout",
            "worker can independently resolve and reverify the same closure",
            "Only that qualified input use rejects a supplied physical path",
            "Cutover of one use never changes the path policy of another use",
            "`transitional_unmanaged_v1`",
            "existing trusted-local v1 path validation remains available for that use",
            "`legacy_v1` is a provenance classification; it neither grants nor prohibits",
            "Project/Sample stage does not reject physical paths",
        ),
    )
    _assert_contains(
        "storage resolver",
        storage,
        (
            "A request for `managed_revision_v1` fails admission",
            "Once an input use is qualified and cut over, it never falls back",
            "An unqualified input use instead follows only an explicitly permitted "
            "`transitional_unmanaged_v1` resolver",
            "regular-file InputFile support is additive",
            "transitional adapter-private v1 payload",
            "never copied into the binding envelope",
            "redacted from public projections",
        ),
    )
    _assert_contains(
        "consequences",
        consequences,
        (
            "no unconditional global physical-path cutover",
            "regular-file support alone cannot justify a fully managed adapter claim",
            "Project/Sample delivery in particular preserves current trusted-local "
            "adapter validation",
        ),
    )

    decision_contract = _section(decision, "Decision", level=2)
    normative_prose = decision_contract[
        : decision_contract.index("\n### Normative compatibility matrix")
    ]
    normative_scope = "\n".join(
        (
            normative_prose,
            consequences,
        )
    )
    present = _unconditional_global_path_rejections(normative_scope)
    assert not present, f"decision still contains global cutover claims: {present}"

    assert _unconditional_global_path_rejections(
        "All input paths become forbidden after rev10."
    )
    assert _unconditional_global_path_rejections(
        "The registry stage rejects every frontend server path."
    )
    assert not _unconditional_global_path_rejections(
        "Every physical path for that qualified input use is rejected."
    )


def test_managed_resource_ownership_and_closure_qualification_are_explicit() -> None:
    decision = DECISION_PATH.read_text(encoding="utf-8")
    storage = _section(decision, "Storage pools and input revisions")
    sequence = _section(decision, "Migration and compatibility sequence")

    platform_contracts = (
        "opaque resource and revision IDs",
        "Project and StoragePool authorization",
        "generic closure digest",
    )
    adapter_contracts = (
        "Bowtie2 prefix",
        "FASTA sidecars",
        "STAR",
        "Salmon",
        "SortMeRNA",
        "producer, tool, container, and build parameters",
    )
    non_inference_sources = (
        "filename",
        "directory name",
        "sample sheet",
        "configuration path",
        "adapter-private payload",
    )
    _assert_contains(
        "managed closure ownership",
        storage,
        platform_contracts + adapter_contracts + non_inference_sources,
    )

    qualification_contracts = (
        "PR #154",
        "stale",
        "exact execution closure",
        "does not requalify",
        "runtime remains unavailable and fails closed",
        "input use may remain transitional",
        "transitional catalog provenance does not admit a stale execution runtime",
        "exact-HEAD qualification",
    )
    _assert_contains("qualification boundaries", sequence, qualification_contracts)


def test_legacy_classification_does_not_set_per_use_path_policy() -> None:
    decision = DECISION_PATH.read_text(encoding="utf-8")
    legacy = _section(
        decision,
        "Legacy compatibility, capability-scoped transition, and conservative migration",
    )

    _assert_contains(
        "legacy classification",
        legacy,
        (
            "`legacy_v1` is a provenance classification; it neither grants nor "
            "prohibits physical-path submission",
            "For a new catalog-aware snapshot, path acceptance is decided only by "
            "the frozen adapter/input-use capability and provenance mode",
            "historical v1 snapshot has no synthesized per-use capability record",
            "frozen v1 workflow contract and digest",
            "private compatibility resolver",
            "Every pre-migration run and every pre-migration validation snapshot",
        ),
    )


def test_historical_artifact_provenance_remains_unresolved_without_inference() -> None:
    decision = DECISION_PATH.read_text(encoding="utf-8")
    artifact = _section(decision, "Artifact contributor provenance")

    _assert_contains(
        "historical artifact provenance",
        artifact,
        (
            "It does not infer contributors from filenames, output paths, QC labels, "
            "or adapter-private scientific roles",
            "Only a historical Legacy artifact without evidence uses "
            "`unresolved_legacy`",
            "it is never assigned a convenient sample",
        ),
    )


def test_architecture_overview_links_the_data_foundation_decision() -> None:
    overview = (
        REPOSITORY_ROOT / "docs" / "architecture" / "platform-overview.md"
    ).read_text(encoding="utf-8")

    assert "(data-foundation-decision.md)" in overview
