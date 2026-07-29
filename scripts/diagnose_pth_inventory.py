#!/usr/bin/env python3
"""Print a path-redacted inventory of executable startup hooks."""

from __future__ import annotations

import base64
import csv
from email.parser import Parser
import hashlib
import json
from pathlib import Path
import re
import sys
import sysconfig


def _identity(path: Path) -> tuple[str, str]:
    message = Parser().parsestr(path.read_text(encoding="utf-8"))
    return message.get("Name", "<missing>"), message.get("Version", "<missing>")


def main() -> int:
    site_root = Path(sysconfig.get_path("purelib")).resolve(strict=True)
    prefix = Path(sys.prefix).resolve(strict=True)
    conda_root = prefix / "conda-meta"
    for hook in sorted(site_root.glob("*.pth")):
        raw = hook.read_bytes()
        lines = tuple(
            line.strip()
            for line in raw.decode("utf-8").splitlines()
            if line.strip() and not line.strip().startswith("#")
        )
        executable = tuple(
            line for line in lines if re.match(r"^import(?:[ \t]|$)", line)
        )
        if not executable:
            continue
        digest = hashlib.sha256(raw).hexdigest()
        urlsafe_digest = (
            base64.urlsafe_b64encode(hashlib.sha256(raw).digest())
            .decode("ascii")
            .rstrip("=")
        )
        print(
            "PTH"
            f" basename={hook.name}"
            f" sha256={digest}"
            f" lines={len(lines)}"
            f" executable={len(executable)}"
            f" path={len(lines) - len(executable)}"
            f" trailing_newline={raw.endswith(bytes([10]))}"
        )
        record_claims: list[str] = []
        for metadata in sorted(site_root.glob("*.dist-info")):
            record = metadata / "RECORD"
            if not record.is_file():
                continue
            rows = csv.reader(record.read_text(encoding="utf-8").splitlines())
            for row in rows:
                if len(row) != 3 or row[0] != hook.name:
                    continue
                name, version = _identity(metadata / "METADATA")
                record_claims.append(
                    f"{name}=={version}"
                    f":hash={row[1] == 'sha256=' + urlsafe_digest}"
                    f":size={row[2] == str(len(raw))}"
                )
        print("RECORD claims=" + (",".join(record_claims) or "<none>"))
        conda_claims: list[str] = []
        if conda_root.is_dir():
            relative = hook.relative_to(prefix).as_posix()
            for inventory in sorted(conda_root.glob("*.json")):
                metadata = json.loads(inventory.read_text(encoding="utf-8"))
                for entry in metadata.get("paths_data", {}).get("paths", []):
                    if entry.get("_path") != relative:
                        continue
                    conda_claims.append(
                        f"{metadata.get('name')}=={metadata.get('version')}"
                        f"={metadata.get('build')}"
                        f":path_type={entry.get('path_type')}"
                        f":sha256={entry.get('sha256') == digest}"
                        f":sha256_in_prefix="
                        f"{entry.get('sha256_in_prefix') == digest}"
                        f":size={entry.get('size_in_bytes') == len(raw)}"
                    )
        print("CONDA claims=" + (",".join(conda_claims) or "<none>"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
