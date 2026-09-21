# Stage naming retirement

Status: Accepted and implemented in PR-4, 2026-09-21.

The ENCODE workflow uses semantic names for configuration, rules, and public
artifact metadata. This is a breaking change; there are no compatibility
aliases, duplicate outputs, or automatic rewrites of existing results.

| Retired name | Required name |
| --- | --- |
| Config `stage4b` | `replicate_analysis` |
| Config `stage5` | `chipseq_idr` |
| Rule and output type `stage3_qc_summary` | `project_qc_summary` |
| `results/multiqc/stage3_qc_summary.tsv` | `results/multiqc/project_qc_summary.tsv` |
| Rule and script `stage5b_summary` / `scripts/stage5b_summary.py` | `chipseq_idr_summary` / `scripts/chipseq_idr_summary.py` |
| `<exp>.stage5b.summary.log` | `<exp>.chipseq_idr.summary.log` |
| Public-validation input column `stage5_idr` | `chipseq_idr` |

Users must rename the two configuration keys while retaining their boolean
values. Platform form inputs already use the semantic names and keep their
`{enabled: boolean}` representation; the adapter renders the same-named
scientific boolean keys. The authoring contract advances to 2.0.0. Revalidate
saved drafts before submitting new runs. Imported Bundles containing retired
keys must also be migrated before fresh validation; the upstream Bundle render
contract is unchanged.

Update downstream QC readers, rule targets, script invocations, and manifest
output-type filters to the new names. Existing results and historical release
evidence are retained as-is. ChIP-seq IDR scientific outputs under `06_idr/`,
the separate reproducibility summaries, and scientific gating behavior remain
unchanged.

Development-stage decision (2026-09-21): retain the validator's existing silent
handling of unknown keys, including retired keys, without adding rejection
checks or rejection tests; the platform is not yet in production and configs
migrate with the code in the same PR. No path-source refactor or
command-ownership change is part of this decision.
