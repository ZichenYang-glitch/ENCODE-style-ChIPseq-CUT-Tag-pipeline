"""Artifact catalog query interface."""

from pathlib import Path
from typing import Optional

from encode_pipeline.artifacts.models import (
    Artifact,
    load_artifacts,
)


_CATALOG_PATH = Path(__file__).resolve().parents[3] / "docs" / "architecture" / "artifact-inventory.yaml"


def _default_catalog_path() -> str:
    return str(_CATALOG_PATH)


def load_catalog(path: Optional[str] = None) -> list[Artifact]:
    """Load the artifact inventory as a list of validated ``Artifact`` objects.

    If ``path`` is not provided, uses the default inventory at
    ``docs/architecture/artifact-inventory.yaml`` relative to the repository
    root.
    """
    if path is None:
        path = _default_catalog_path()
    return load_artifacts(path)


def describe(output_type: str, catalog: Optional[list[Artifact]] = None) -> Optional[Artifact]:
    """Return the ``Artifact`` whose ``manifest_output_type`` equals ``output_type``.

    Returns ``None`` if no match is found. If ``catalog`` is omitted, the
    default inventory is loaded.
    """
    if catalog is None:
        catalog = load_catalog()
    for artifact in catalog:
        if artifact.manifest_output_type == output_type:
            return artifact
    return None


def expected_for_assay(assay: str, catalog: Optional[list[Artifact]] = None) -> list[Artifact]:
    """Return artifacts expected for a given assay.

    Matching rules:

    - ``assay_gate == "all"`` always matches.
    - ``assay_gate == "peak_centric"`` matches when assay is one of
      ``chipseq``, ``cuttag``, or ``atac``.
    - ``assay_gate == "idr"`` matches when assay is one of ``chipseq``,
      ``cuttag``, or ``atac`` (IDR is a peak-centric reproducibility path).
    - Otherwise the assay gate must equal the assay.

    If ``catalog`` is omitted, the default inventory is loaded.
    """
    if catalog is None:
        catalog = load_catalog()

    peak_centric = {"chipseq", "cuttag", "atac"}

    results = []
    for artifact in catalog:
        gate = artifact.assay_gate
        if gate == "all":
            results.append(artifact)
        elif gate == "peak_centric" and assay in peak_centric:
            results.append(artifact)
        elif gate == "idr" and assay in peak_centric:
            results.append(artifact)
        elif gate == assay:
            results.append(artifact)
    return results
