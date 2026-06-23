"""Sample record dataclass."""

from dataclasses import dataclass


@dataclass(frozen=True)
class SampleRecord:
    """A validated sample row from the sample sheet.

    This is a read-only, hashable representation of the sample dicts
    produced by load_and_validate_samples(). It exists for type safety
    inside Python helpers; Snakemake rules continue to receive plain dicts.
    """

    id: str
    fq1: str
    fq2: str | None
    layout: str  # "PE" or "SE"
    assay: str  # "chipseq", "cuttag", "atac", "mnase"
    target: str
    peak_mode: str  # "narrow", "broad", "nucleosome"
    genome: str
    bt2_idx: str
    control_bam: str
    role: str  # "treatment" or "control"
    control_sample: str
    experiment: str
    condition: str
    replicate: int
    biological_replicate: int
    technical_replicate: int

    @classmethod
    def from_dict(cls, sample: dict) -> "SampleRecord":
        """Build a SampleRecord from a legacy sample dict."""
        return cls(
            id=sample["id"],
            fq1=sample["fq1"],
            fq2=sample.get("fq2") or None,
            layout=sample["layout"],
            assay=sample["assay"],
            target=sample["target"],
            peak_mode=sample["peak_mode"],
            genome=sample["genome"],
            bt2_idx=sample["bt2_idx"],
            control_bam=sample.get("control_bam", ""),
            role=sample["role"],
            control_sample=sample.get("control_sample", ""),
            experiment=sample["experiment"],
            condition=sample["condition"],
            replicate=int(sample["replicate"]),
            biological_replicate=int(sample["biological_replicate"]),
            technical_replicate=int(sample["technical_replicate"]),
        )
