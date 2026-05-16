# replicates.smk — Stage 4b replicate-aware grouped outputs
# ===========================================================
# Produces biological-replicate BAMs, pooled treatment BAMs, and
# pooled control BAMs. All rules run after per-sample final.bam.
#
# Rules:
#   merge_biorep_bam    — merge technical replicates within a biological replicate
#   pool_treatment_bam  — merge treatment biorep BAMs for the experiment
#   pool_control_bam    — merge/symlink control BAM(s) for the experiment


# ---------------------------------------------------------------------------
# Helper: list of (sample_id,) final.bam paths for a (experiment, role, bio_rep) key
# ---------------------------------------------------------------------------

def _biorep_techrep_bams(experiment, role, bio_rep):
    """Return sorted final.bam paths for every tech-rep in a biorep group."""
    key = (experiment, role, bio_rep)
    sids = sorted(BIO_REP_GROUPS.get(key, []))
    return [
        f"{OUTDIR}/{sid}/02_align/{sid}.final.bam"
        for sid in sids
    ]


# ---------------------------------------------------------------------------
# 1. Biological-replicate BAM — merge technical replicates or symlink
# ---------------------------------------------------------------------------

rule merge_biorep_bam:
    output:
        bam = f"{OUTDIR}/experiments/{{experiment}}/02_align/biorep{{bio_rep}}.final.bam",
        bai = f"{OUTDIR}/experiments/{{experiment}}/02_align/biorep{{bio_rep}}.final.bam.bai",
    input:
        lambda wc: _biorep_techrep_bams(wc.experiment, "treatment", int(wc.bio_rep)),
    params:
        bio_rep = "{bio_rep}",
    wildcard_constraints:
        bio_rep = r"\d+",
    threads: THREADS,
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        set -e -o pipefail
        set -- {input:q}
        mkdir -p $(dirname {output.bam:q})

        if [[ "$#" -eq 1 ]]; then
            # Single technical replicate — symlink, then index for correct .bai
            ln -sf "$(readlink -f "$1")" {output.bam:q}
            samtools index -@ {threads} {output.bam:q}
        else
            # Multiple technical replicates — samtools merge + index
            samtools merge -@ {threads} {output.bam:q} "$@"
            samtools index -@ {threads} {output.bam:q}
        fi
        """


# ---------------------------------------------------------------------------
# 2. Pooled treatment BAM — merge all treatment biorep BAMs for an experiment
# ---------------------------------------------------------------------------

def _treatment_biorep_bams(experiment):
    """Return sorted treatment biorep BAM paths for an experiment."""
    bioreps = sorted(set(
        br for (exp, role, br) in BIO_REP_GROUPS
        if exp == experiment and role == "treatment"
    ))
    return [
        f"{OUTDIR}/experiments/{experiment}/02_align/biorep{br}.final.bam"
        for br in bioreps
    ]


rule pool_treatment_bam:
    output:
        bam = f"{OUTDIR}/experiments/{{experiment}}/02_align/{{experiment}}.pooled.final.bam",
        bai = f"{OUTDIR}/experiments/{{experiment}}/02_align/{{experiment}}.pooled.final.bam.bai",
    input:
        lambda wc: _treatment_biorep_bams(wc.experiment),
    threads: THREADS,
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        set -e -o pipefail
        set -- {input:q}
        mkdir -p $(dirname {output.bam:q})
        samtools merge -@ {threads} {output.bam:q} "$@"
        samtools index -@ {threads} {output.bam:q}
        """


# ---------------------------------------------------------------------------
# 3. Pooled control BAM — merge or symlink control(s) for the experiment
# ---------------------------------------------------------------------------

def _pooled_control_inputs(experiment):
    """Return list of control final.bam paths for pooled control BAM.

    1 unique control → [single_path]  (shell symlinks)
    >1 unique controls → [path1, path2, ...]  (shell merges)
    0 controls → empty list (rule not scheduled)
    """
    resolved = _resolve_experiment_controls(experiment)
    return [path for (_stype, path) in resolved]


rule pool_control_bam:
    output:
        bam = f"{OUTDIR}/experiments/{{experiment}}/02_align/{{experiment}}.pooled.control.final.bam",
        bai = f"{OUTDIR}/experiments/{{experiment}}/02_align/{{experiment}}.pooled.control.final.bam.bai",
    input:
        lambda wc: _pooled_control_inputs(wc.experiment),
    threads: THREADS,
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        set -e -o pipefail
        set -- {input:q}
        mkdir -p $(dirname {output.bam:q})

        if [[ "$#" -eq 1 ]]; then
            # Single control BAM → symlink, then index to produce correct .bai
            ln -sf "$(readlink -f "$1")" {output.bam:q}
            samtools index -@ {threads} {output.bam:q}
        else
            # Multiple control BAMs → samtools merge + index
            samtools merge -@ {threads} {output.bam:q} "$@"
            samtools index -@ {threads} {output.bam:q}
        fi
        """
