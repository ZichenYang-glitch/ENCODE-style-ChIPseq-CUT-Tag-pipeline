# Platform worker smoke profile

This deterministic profile contains one control-only row while `use_control` is
disabled, so the bundled workflow resolves to its empty `all` target. The
SCI-1 identity marker remains a production prerequisite, so the tiny run now
executes `analysis_identity_marker` and `all` (2 of 2 steps) and writes the
versioned `results/.helixweave/analysis_identity.json` marker. It proves the
API → Redis/RQ → independent worker → Snakemake → SQLite lifecycle without
running bioinformatics tools or using scientific data. The scientific workflow
and real-data behavior remain covered by the existing Snakemake test profiles.
