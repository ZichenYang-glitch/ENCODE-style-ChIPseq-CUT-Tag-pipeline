# ENCODE command ownership

Accepted in PR-5b, 2026-09-21.

The ENCODE adapter now constructs its own Snakemake `CommandSpec` through
`encode_execution.py`. Source/runtime coordinates are bound during composition
and retained when the adapter binds a Reference Profile. The command service
only delegates; argv, environment, preflight, and `cwd=None` retain their
pre-migration values. No public `CommandSpec`, `WorkspacePlan`, or adapter
protocol fields/signatures change.

Workspace planning adds `config/encode-execution.json`, canonical UTF-8 JSON
with sorted keys, compact separators, and a final newline. Its exact fields
are `schema_version` (`"1.0.0"`), `cores` (integer 1–1024, excluding booleans),
and `workspace_contract_sha256` (SHA-256). The digest binds cores, ordered
workspace directories, and the paths/content hashes of the other planned
files. This follows the bulk adapter's planned-file integrity pattern.

Command construction reparses the planned bytes, rejects duplicate/unknown
keys and invalid types/ranges, and requires exact canonical bytes and a
matching digest. An existing materialized execution file must be a regular
non-symlink file matching the plan byte-for-byte. Planning remains side-effect
free; building directly from an unmaterialized plan remains supported.
Older workspace plans without this file must be replanned, not silently
assigned a cores value.

The registry still authorizes only explicitly trusted registered instances.
Reference-bound ENCODE delegation additionally requires the same adapter type,
metadata, and capabilities as the registered instance. This trusted branch
retains the existing `cwd=None` exception. The generic/bulk command path still
rejects missing cwd and retains all existing workspace checks.

The committed pre-migration command fixture captures eight real CommandBuilder
outputs (default/1/8/1024 cores, with and without the admitted Conda runtime).
Tests expand only temporary path placeholders and compare argv/environment/
preflight bytes as well as the complete serialized command with the adapter
and service outputs. Workspace contract assertions deliberately include the
new private file and its revised file count.
