# Platform worker smoke profile

This deterministic profile contains one control-only row while `use_control` is
disabled, so the bundled workflow resolves to its empty `all` target. It proves
the API → Redis/RQ → independent worker → Snakemake → SQLite lifecycle without
running bioinformatics tools or using scientific data. The scientific workflow
and real-data behavior remain covered by the existing Snakemake test profiles.
