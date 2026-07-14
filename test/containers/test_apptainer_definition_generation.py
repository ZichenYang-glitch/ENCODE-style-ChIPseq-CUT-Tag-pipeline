"""Offline behavior contracts for generated Apptainer definitions."""

from pathlib import Path

from containers import generate_apptainer_defs


REPO_ROOT = Path(__file__).resolve().parents[2]


def test_definition_uses_declared_environment_and_runner_entrypoint(tmp_path):
    envs_dir = tmp_path / "envs"
    envs_dir.mkdir()
    lock = envs_dir / "runner.lock"
    lock.write_text("@EXPLICIT\n", encoding="utf-8")
    lock.with_suffix(".yml").write_text(
        "name: chipseq-runner\ndependencies: []\n",
        encoding="utf-8",
    )

    rendered = generate_apptainer_defs.render_def(envs_dir, lock)

    assert "workflow/envs/runner.lock" in rendered
    assert "conda create -n chipseq-runner" in rendered
    assert "export PATH=/opt/conda/envs/chipseq-runner/bin:$PATH" in rendered
    assert '%runscript\n    exec snakemake "$@"' in rendered
    assert rendered.endswith("\n")


def test_definition_falls_back_to_lock_stem_without_runner_entrypoint(tmp_path):
    lock = tmp_path / "analysis.lock"
    lock.write_text("@EXPLICIT\n", encoding="utf-8")

    rendered = generate_apptainer_defs.render_def(tmp_path, lock)

    assert "conda create -n analysis" in rendered
    assert "export PATH=/opt/conda/envs/analysis/bin:$PATH" in rendered
    assert "%runscript" not in rendered


def test_main_generates_every_lock_in_sorted_order_under_requested_output(
    tmp_path, capsys
):
    envs_dir = tmp_path / "envs"
    output_dir = tmp_path / "definitions"
    envs_dir.mkdir()
    for name in ("zeta", "alpha"):
        (envs_dir / f"{name}.lock").write_text("@EXPLICIT\n", encoding="utf-8")
    (envs_dir / "zeta.yml").write_text("name: custom-zeta\n", encoding="utf-8")

    generate_apptainer_defs.main(
        ["--envs-dir", str(envs_dir), "--output-dir", str(output_dir)]
    )

    assert sorted(path.name for path in output_dir.iterdir()) == [
        "Apptainer.alpha.def",
        "Apptainer.zeta.def",
    ]
    assert "conda create -n alpha" in (output_dir / "Apptainer.alpha.def").read_text(
        encoding="utf-8"
    )
    assert "conda create -n custom-zeta" in (
        output_dir / "Apptainer.zeta.def"
    ).read_text(encoding="utf-8")
    assert capsys.readouterr().out.splitlines() == [
        f"Generated {output_dir / 'Apptainer.alpha.def'}",
        f"Generated {output_dir / 'Apptainer.zeta.def'}",
        f"Generated 2 Apptainer definition files in {output_dir}",
    ]


def test_tracked_definitions_match_the_generator():
    envs_dir = REPO_ROOT / "workflow" / "envs"
    for lock in sorted(envs_dir.glob("*.lock")):
        tracked = REPO_ROOT / "containers" / f"Apptainer.{lock.stem}.def"
        rendered = generate_apptainer_defs.render_def(envs_dir, lock)
        assert tracked.read_text(encoding="utf-8") == rendered.rstrip("\n") + "\n"
