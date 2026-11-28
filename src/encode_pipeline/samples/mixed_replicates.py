"""Per-experiment MACS3 decisions for an explicitly requested mixed run."""

from dataclasses import dataclass
from typing import Literal

from encode_pipeline.config.cuttag import validate_cuttag_config
from encode_pipeline.errors import ValidationError


@dataclass(frozen=True)
class MixedReplicateRoute:
    method: Literal["idr", "consensus"]
    n_bioreps: int
    run_macs3_consensus: bool

    def as_record(self):
        record = {"route": self.method, "n_bioreps": self.n_bioreps, "caller": "macs3"}
        if not self.run_macs3_consensus:
            record["auxiliary_consensus_skipped"] = "insufficient_biological_replicates"
        return record


def treatment_groups(samples):
    groups = {}
    for sample in samples:
        if sample["role"] == "treatment":
            groups.setdefault(sample["experiment"], []).append(sample)
    return groups


def resolve_mixed_replicate_routes(
    samples, config, idr_flags, *, error_cls=ValidationError
):
    """Return only groups governed by the new policy. Strict runs have no overlay.

    The caller supplies the effective IDR flags from config validation. Both
    raw run records and normalized configs share the same sample admission.
    """
    repro = (config or {}).get("reproducibility") or {}
    if repro.get("mixed_replicate_policy") != "per_group_fallback":
        return {}

    minimum = int(repro.get("consensus", {}).get("min_replicates", 2))
    seacr = validate_cuttag_config(config.get("cuttag", {}), error_cls)["seacr"][
        "enabled"
    ]
    routes = {}
    for experiment, rows in treatment_groups(samples).items():
        first = rows[0]
        assay, mode = first["assay"], first["peak_mode"]
        if assay not in ("chipseq", "atac", "cuttag"):
            continue
        count = len({row["biological_replicate"] for row in rows})
        eligible = mode == "narrow" and idr_flags[f"{assay}_narrow"]
        if eligible and count not in (2, 3):
            raise error_cls(
                f"Experiment {experiment!r}: per_group_fallback supports exactly 2 or 3 "
                f"biological replicates for {assay} narrow, found {count}."
            )
        if count < 2:
            continue  # No consensus branch is scheduled for unrelated singletons.
        if count < minimum:
            if seacr and assay == "cuttag" and first["layout"] == "PE":
                raise error_cls(
                    f"Experiment {experiment!r}: SEACR consensus requires "
                    f"min_replicates={minimum}, but has {count} biological replicates."
                )
            if not eligible:
                raise error_cls(
                    f"Experiment {experiment!r}: MACS3 consensus requires "
                    f"min_replicates={minimum}, but has {count} biological replicates."
                )
        if eligible:
            routes[experiment] = MixedReplicateRoute(
                method="idr" if count == 2 else "consensus",
                n_bioreps=count,
                run_macs3_consensus=count >= minimum,
            )
    if not any(route.method == "idr" for route in routes.values()):
        raise error_cls(
            "per_group_fallback requires at least one eligible narrow-peak IDR experiment "
            "with exactly 2 biological replicates."
        )
    return routes


def mixed_replicate_notice(routes):
    """Describe the chosen methods without disclosing input paths."""
    groups = []
    for experiment, route in sorted(routes.items()):
        suffix = "" if route.run_macs3_consensus else ", auxiliary consensus skipped"
        groups.append(
            f"{experiment} ({route.n_bioreps} biological replicates, {route.method}{suffix})"
        )
    return (
        "Per-group MACS3 reproducibility: "
        + "; ".join(groups)
        + ". IDR and consensus scores use different scales."
    )
