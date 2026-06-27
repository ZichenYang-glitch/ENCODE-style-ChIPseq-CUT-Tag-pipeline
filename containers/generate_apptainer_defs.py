#!/usr/bin/env python3
"""Generate Apptainer definition files using workflow/envs/*.lock.

The matching workflow/envs/*.yml is read only to determine the conda
environment name (`name:`). The generated definition copies the explicit
lockfile into the container and installs from it for reproducible builds.

Usage:
    python3 containers/generate_apptainer_defs.py
    python3 containers/generate_apptainer_defs.py --envs-dir workflow/envs --output-dir containers
"""

import argparse
from pathlib import Path

import yaml


TEMPLATE = """\
Bootstrap: docker
From: condaforge/miniforge3:24.11.3-0

%files
    {lock_rel} /opt/pipeline/{lock_rel}

%post
    conda create -n {name} --file /opt/pipeline/{lock_rel} -y
    conda clean -afy

%environment
    export PATH=/opt/conda/envs/{name}/bin:$PATH
{runscript}
"""

RUNNER_RUNSCRIPT = """
%runscript
    exec snakemake "$@"
"""


def render_def(envs_dir: Path, lock_path: Path) -> str:
    yml_path = lock_path.with_suffix(".yml")
    if yml_path.exists():
        with open(yml_path) as fh:
            spec = yaml.safe_load(fh)
        name = spec.get("name", lock_path.stem)
    else:
        name = lock_path.stem
    lock_rel = Path("workflow/envs") / lock_path.name
    runscript = RUNNER_RUNSCRIPT if name == "chipseq-runner" else ""
    return TEMPLATE.format(
        name=name,
        lock_rel=lock_rel.as_posix(),
        runscript=runscript,
    )


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate Apptainer definition files from conda env lockfiles."
    )
    parser.add_argument("--envs-dir", default="workflow/envs")
    parser.add_argument("--output-dir", default="containers")
    args = parser.parse_args(argv)

    envs_dir = Path(args.envs_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    generated = 0
    for lock in sorted(envs_dir.glob("*.lock")):
        def_text = render_def(envs_dir, lock)
        out_path = output_dir / f"Apptainer.{lock.stem}.def"
        with open(out_path, "w") as fh:
            # Ensure exactly one trailing newline for standard text-file
            # formatting while avoiding extra blank lines.
            fh.write(def_text.rstrip("\n") + "\n")
        print(f"Generated {out_path}")
        generated += 1

    print(f"Generated {generated} Apptainer definition files in {output_dir}")


if __name__ == "__main__":
    main()
