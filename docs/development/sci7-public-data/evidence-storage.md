# Evidence storage and integrity

Seventeen source-evidence files from the initial SCI-7 preparation commit are
stored as deterministic, lossless gzip members. Their original bytes contained
trailing whitespace that made the full base-to-head `git diff --check` fail.
The bytes were **not** normalized or edited; only the storage representation
changed. This preserves the original evidence while keeping the complete PR
diff auditable.

[`evidence-integrity.tsv`](evidence-integrity.tsv) records, for every compressed
file:

- original and compressed repository-relative paths;
- original byte count and SHA256 of the **decompressed original bytes**;
- SHA256 of the gzip member;
- deterministic gzip settings (level 9, empty original filename, `mtime=0`);
- the official query URL, or a note that the file is a derived evidence index.

The original `queries.tsv` and other query indexes may contain the original
uncompressed `local_response` names. Use `evidence-integrity.tsv` to map those
names to the current `.gz` paths rather than changing the preserved original
query bytes.

## Reading and verification

Decompress a file to standard output without changing the repository:

```bash
gzip -dc docs/development/sci7-public-data/metadata/geo-sra/queries.tsv.gz
```

Verify every compressed member and its original identity:

```bash
python3 - <<'PY'
import csv, gzip, hashlib
from pathlib import Path
root = Path('docs/development/sci7-public-data')
with (root / 'evidence-integrity.tsv').open(newline='') as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))
assert len(rows) == 17
for row in rows:
    payload = (root / row['compressed_path']).read_bytes()
    original = gzip.decompress(payload)
    assert hashlib.sha256(payload).hexdigest() == row['gzip_sha256']
    assert len(original) == int(row['original_bytes'])
    assert hashlib.sha256(original).hexdigest() == row['original_sha256']
print('17 compressed evidence files: decompression and SHA256 OK')
PY
```

The complete-PR check must be run against the committed range, not only against
a clean working tree:

```bash
git diff dbc66d6a66cbc6d2ef95783e90203d9b353135a0..HEAD --check
```

The two pre-existing ENCODE blacklist `.bed.gz` files remain governed by their
official ENCODE size/MD5 records in `reference-identity.md`; they are not part of
this 17-file gzip set.
