"""Machine-readable supported-deployment schemas."""

from __future__ import annotations

from importlib import resources
import json


def load_schema(name: str) -> dict[str, object]:
    if name not in {
        "deployment-bundle-v1.schema.json",
        "deployment-state-v1.schema.json",
    }:
        raise ValueError("unsupported deployment schema")
    return json.loads(resources.files(__package__).joinpath(name).read_bytes())


__all__ = ["load_schema"]
