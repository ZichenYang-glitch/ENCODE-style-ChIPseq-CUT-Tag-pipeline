#!/usr/bin/env python3
"""Check the scientific-workflow evidence ledger against repository sources.

Reads files in this checkout only. It does not start a workflow, download
data, open a run, or treat a software qualification as a scientific conclusion.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "docs/development/scientific-workflow-evidence-ledger.md"
QUALIFICATION = (
    ROOT
    / "src/encode_pipeline/contracts/nfcore_rnaseq"
    / "default-execution-qualification-1.1.0.json"
)
RESULTS = (
    ROOT
    / "src/encode_pipeline/contracts/nfcore_rnaseq"
    / "results-contract-3.26.0.json"
)
MANIFEST = (
    ROOT
    / "src/encode_pipeline/contracts/nfcore_rnaseq"
    / "execution-implementation-manifest-1.0.0.json"
)
LOCK = ROOT / "config/hitrac_preprocess/tools.lock.json"
SUMMARY = ROOT / "docs/development/sci7-public-data/sci7-batch-summary.md"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str, failures: list[str]) -> None:
    if not condition:
        failures.append(message)


def main() -> int:
    failures: list[str] = []
    ledger = LEDGER.read_text(encoding="utf-8")
    qualification = json.loads(QUALIFICATION.read_text(encoding="utf-8"))
    results = json.loads(RESULTS.read_text(encoding="utf-8"))
    lock = json.loads(LOCK.read_text(encoding="utf-8"))
    implementation = qualification["execution_implementation"]

    facts = {
        "qualification file": digest(QUALIFICATION),
        "results contract file": digest(RESULTS),
        "manifest file": digest(MANIFEST),
        "Hi-TrAC lock file": digest(LOCK),
        "qualification record": qualification["record_sha256"],
        "execution aggregate": implementation["aggregate_sha256"],
        "manifest record": implementation["manifest_sha256"],
        "runtime identity": qualification["runtime"]["identity_sha256"],
        "results contract record": qualification["result_contract"]["sha256"],
        "upstream commit": results["upstream"]["commit"],
        "tracPre2 script": lock["upstream"]["script_sha256"],
    }
    require(
        qualification["candidate_status"] == "source-enabled",
        "Bulk candidate_status is no longer source-enabled",
        failures,
    )
    require(
        qualification["adapter_command_contract"]["workflow_id"] == "bulk-rnaseq",
        "Bulk qualification workflow_id changed",
        failures,
    )
    require(
        implementation["file_count"] == 119,
        "Bulk closure file_count changed",
        failures,
    )
    require(
        facts["manifest file"] == facts["manifest record"],
        "manifest file bytes no longer match the qualification record",
        failures,
    )
    require(
        facts["results contract file"] == facts["results contract record"],
        "results contract bytes no longer match the qualification record",
        failures,
    )
    require(
        results["workflow_id"] == "bulk-rnaseq"
        and results["upstream"]["release"] == "3.26.0"
        and results["upstream"]["route"] == "star_salmon",
        "Bulk results contract coordinates changed",
        failures,
    )
    require(
        "No production default is frozen." in SUMMARY.read_text(encoding="utf-8"),
        "SCI-7 batch summary no longer states that no production default is frozen",
        failures,
    )

    policy = (ROOT / "src/encode_pipeline/config/scientific_policy.py").read_text(
        encoding="utf-8"
    )
    defaults = (ROOT / "src/encode_pipeline/config/defaults.py").read_text(
        encoding="utf-8"
    )
    config = (ROOT / "config/config.yaml").read_text(encoding="utf-8")
    metadata = (ROOT / "workflow/rules/metadata.smk").read_text(encoding="utf-8")
    hitrac_results = (
        ROOT / "src/encode_pipeline/adapters/hitrac_preprocess/results.py"
    ).read_text(encoding="utf-8")
    conformance = (ROOT / "docs/development/adapter-conformance.md").read_text(
        encoding="utf-8"
    )
    require(
        "SCIENTIFIC_POLICY_REVISION = 7" in policy,
        "ENCODE scientific policy revision is no longer 7",
        failures,
    )
    require(
        'return "retain" if assay == "cuttag" else "exclude"' in metadata,
        "ENCODE auto duplicate policy source changed",
        failures,
    )
    require(
        "mapq: 30" in config and 'duplicate_analysis: "auto"' in config,
        "ENCODE config defaults changed",
        failures,
    )
    expected_insert = {
        "chipseq": 2000,
        "atac": 2000,
        "cuttag": 700,
        "mnase": 500,
    }
    for assay, value in expected_insert.items():
        require(
            f'"{assay}": {value}' in defaults,
            f"Bowtie2 max_insert for {assay} changed",
            failures,
        )
    authoring = (
        ROOT / "src/encode_pipeline/adapters/hitrac_preprocess/authoring.py"
    ).read_text(encoding="utf-8")
    require(
        "DEFAULT_MAPQ = 10" in authoring,
        "Hi-TrAC authoring MAPQ default changed",
        failures,
    )
    for token in (
        'SUMMARY_TYPE = "hitrac_summary"',
        'ALL_TYPE = "hitrac_bedpe_all"',
        'NO_BG_TYPE = "hitrac_bedpe_no_bg"',
        "def qc_source_output_types",
    ):
        require(
            token in hitrac_results,
            f"Hi-TrAC results contract lost {token}",
            failures,
        )
    require(
        "not_configured` until qualified artifact/QC publication is implemented"
        in conformance,
        "adapter-conformance.md lost the publication sentence",
        failures,
    )
    reference_sources = {
        "docs/development/sci7-public-data/minimal-smoke-results.md": (
            "the GRCh38 FASTA matched the official byte counts and MD5",
            "six Bowtie2 index files",
        ),
        "docs/development/sci7-public-data/tf-jund-confirmation.md": (
            "9c3fc2ca260a767530dddb0f26721a6b",
            "e9a2d5a5cd225293646ae298998ad4ea8c11e8f22729bc2d6f5c3dcfefbf8ef8",
        ),
        "docs/development/sci7-public-data/atac-treg-discovery.md": (
            "3bc591be24b77f710b6ba5d41022fc5a",
            "617b10dc7ef90354c3b6af986e45d6d9621242b64ed3a94c9abeac3e45f18c17",
            "six Bowtie2 `.bt2` files",
        ),
        "docs/development/sci7-public-data/reference-identity.md": (
            "gencode.v44.primary_assembly.annotation.gtf.gz",
            "b182a9f3b134b9cc2da566a2b1692557",
            "gencode.vM25.primary_assembly.annotation.gtf.gz",
            "c5125258a0a2c5250ddb4c192abbf4e8",
            "Neither GTF was downloaded.",
        ),
    }
    for relative, needles in reference_sources.items():
        text = (ROOT / relative).read_text(encoding="utf-8")
        for needle in needles:
            require(needle in text, f"{relative} lost {needle}", failures)

    citations = [
        "5f050c81b28cefe8bacd6b499939f46210d166b4f4fd8630a593017a4637603b",
        *facts.values(),
        "source-enabled",
        "No production default is frozen.",
        "encode-style-chipseq-cuttag-atac-mnase",
        "bulk-rnaseq",
        "hitrac-preprocess",
        "bulk_rnaseq.salmon.quant_gene",
        "bulk_rnaseq.salmon.quant_transcript",
        "bulk_rnaseq.salmon.merged.gene_tpm",
        "salmon.mapping_fraction",
        "hitrac_bedpe_all",
        "hitrac_bedpe_no_bg",
        "hitrac_summary",
        'return "retain" if assay == "cuttag" else "exclude"',
        "SCIENTIFIC_POLICY_REVISION = 7",
        "until qualified artifact/QC publication is implemented",
        "IRENE-35",
        "IRENE-21",
        "IRENE-9",
        "IRENE-34",
        "11_17 的计算与核验跟踪在 IRENE-9",
        "试用验收在 IRENE-34",
        "11_17",
        "WashU",
        "9c3fc2ca260a767530dddb0f26721a6b",
        "e9a2d5a5cd225293646ae298998ad4ea8c11e8f22729bc2d6f5c3dcfefbf8ef8",
        "3bc591be24b77f710b6ba5d41022fc5a",
        "617b10dc7ef90354c3b6af986e45d6d9621242b64ed3a94c9abeac3e45f18c17",
        "gencode.v44.primary_assembly.annotation.gtf.gz",
        "gencode.vM25.primary_assembly.annotation.gtf.gz",
        "b182a9f3b134b9cc2da566a2b1692557",
        "c5125258a0a2c5250ddb4c192abbf4e8",
        "minimal-smoke-results.md",
        "tf-jund-confirmation.md",
        "atac-treg-discovery.md",
        "ENCSR213YKJ",
        "ENCSR639PPN",
    ]
    require(
        "工单列表没有 IRENE-9" not in ledger and "工单列表里没有 IRENE-9" not in ledger,
        "ledger still describes IRENE-9 from a transient issue list",
        failures,
    )
    require(
        "大 FASTA、GTF 和 Bowtie2 索引仍未按" not in ledger,
        "ledger still treats the later FASTA and index records as unverified",
        failures,
    )
    for citation in citations:
        require(citation in ledger, f"ledger is missing {citation}", failures)

    if failures:
        print("scientific workflow evidence check failed:", file=sys.stderr)
        for item in failures:
            print(f"- {item}", file=sys.stderr)
        return 1

    print("scientific workflow evidence check passed")
    print(f"qualification {facts['qualification file']}")
    print(
        "bulk aggregate "
        f"{facts['execution aggregate']} files {implementation['file_count']}"
    )
    print(f"results contract {facts['results contract file']}")
    print(f"hitrac lock {facts['Hi-TrAC lock file']}")
    print("software qualifications remain separate from scientific conclusions")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
