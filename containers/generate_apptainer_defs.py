#!/usr/bin/env python3
"""Generate Apptainer definition files from workflow/envs/*.yml.

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
    {env_rel} /opt/pipeline/{env_rel}

%post
    conda env create -f /opt/pipeline/{env_rel} -n {name} -y
    conda clean -afy

%environment
    export PATH=/opt/conda/envs/{name}/bin:$PATH
{runscript}
"""

RUNNER_RUNSCRIPT = """
%runscript
    exec snakemake "$@"
"""


def render_def(env_path: Path) -> str:
    with open(env_path) as fh:
        spec = yaml.safe_load(fh)
    name = spec.get("name", env_path.stem)
    env_rel = Path("workflow/envs") / env_path.name
    runscript = RUNNER_RUNSCRIPT if name == "chipseq-runner" else ""
    return TEMPLATE.format(
        name=name,
        env_rel=env_rel.as_posix(),
        runscript=runscript,
    )


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate Apptainer definition files from conda env YAMLs."
    )
    parser.add_argument("--envs-dir", default="workflow/envs")
    parser.add_argument("--output-dir", default="containers")
    args = parser.parse_args(argv)

    envs_dir = Path(args.envs_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    generated = 0
    for yml in sorted(envs_dir.glob("*.yml")):
        def_text = render_def(yml)
        out_path = output_dir / f"Apptainer.{yml.stem}.def"
        with open(out_path, "w") as fh:
            # Ensure exactly one trailing newline for standard text-file
            # formatting while avoiding extra blank lines.
            fh.write(def_text.rstrip("\n") + "\n")
        print(f"Generated {out_path}")
        generated += 1

    print(f"Generated {generated} Apptainer definition files in {output_dir}")


if __name__ == "__main__":
    main()
