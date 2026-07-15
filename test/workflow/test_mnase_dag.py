"""MNase-specific completion-target DAG contracts."""

import os
import shutil
import subprocess

from conftest import prepare_profile_workdir


MNASE_COMPLETION_RULES = {
    "mnase_split_mono",
    "mnase_split_sub",
    "mnase_split_di",
    "mnase_dyad_bigwig",
    "mnase_mono_bigwig",
    "mnase_qc_summary",
}


def test_direct_pipeline_done_target_schedules_complete_mnase_outputs(
    profiles_dir,
    snakefile,
    snakemake_executable,
):
    profile_dir = os.path.join(profiles_dir, "mnase_pe_noctrl")
    workdir, config_path = prepare_profile_workdir(profile_dir)
    try:
        env = os.environ.copy()
        env["PYTHONDONTWRITEBYTECODE"] = "1"
        env["XDG_CACHE_HOME"] = os.path.join(workdir, ".cache")
        result = subprocess.run(
            [
                snakemake_executable,
                "-s",
                snakefile,
                "--configfile",
                config_path,
                "--dry-run",
                "results/M1/logs/M1.pipeline.done",
            ],
            cwd=workdir,
            capture_output=True,
            text=True,
            check=False,
            env=env,
        )
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    for rule in MNASE_COMPLETION_RULES:
        assert rule in output
