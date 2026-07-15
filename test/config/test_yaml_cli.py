"""Pin _load_yaml, _parse_config_minimal, and validator CLI behavior.

These tests characterize the YAML loading and CLI surfaces in
``encode_pipeline.config.validator`` before PR74 extracts them. They are
intentionally behavior-preserving: any output shape that changes in PR74 must
be matched or updated here.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

VALID_SAMPLES_HEADER = "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n"


def _valid_samples(sid="S1"):
    return (
        VALID_SAMPLES_HEADER + f"{sid}\tR1.fq\tR2.fq\tPE\tchipseq\tT\tnarrow\ths\tidx\n"
    )


def _write_text_config(tmp_path, text, samples_text=None):
    """Write raw YAML config text and optional samples TSV, return paths."""
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(text, encoding="utf-8")

    if samples_text is None:
        samples_text = _valid_samples()
    samples_path = tmp_path / "samples.tsv"
    samples_path.write_text(samples_text, encoding="utf-8")

    return cfg_path, samples_path


def _make_minimal_config(tmp_path, config_text, samples_text=None):
    """Write a config that references the generated samples.tsv.

    If *config_text* is empty, injects a ``samples`` line pointing to the
    generated TSV so the config is valid for CLI success paths.
    """
    cfg_path, samples_path = _write_text_config(
        tmp_path, config_text, samples_text=samples_text
    )
    if config_text.strip() == "":
        cfg_path.write_text(f'samples: "{samples_path}"\n', encoding="utf-8")
    return cfg_path, samples_path


# ---------------------------------------------------------------------------
# _parse_config_minimal direct tests
# ---------------------------------------------------------------------------


def test_parse_config_minimal_top_level_scalars(tmp_path):
    from encode_pipeline.config.validator import _parse_config_minimal

    cfg_path, _ = _make_minimal_config(
        tmp_path,
        'outdir: "results"\nthreads: 8\nuse_control: false\nmapq: 30\npi: 3.14\n',
    )
    result = _parse_config_minimal(str(cfg_path))
    assert result == {
        "outdir": "results",
        "threads": 8,
        "use_control": False,
        "mapq": 30,
        "pi": "3.14",
    }


def test_parse_config_minimal_quoted_strings(tmp_path):
    from encode_pipeline.config.validator import _parse_config_minimal

    cfg_path, _ = _make_minimal_config(
        tmp_path,
        "single: 'quoted'\ndouble: \"quoted\"\nunquoted: value\n",
    )
    result = _parse_config_minimal(str(cfg_path))
    assert result == {
        "single": "quoted",
        "double": "quoted",
        "unquoted": "value",
    }


def test_parse_config_minimal_ints_and_bools(tmp_path):
    from encode_pipeline.config.validator import _parse_config_minimal

    cfg_path, _ = _make_minimal_config(
        tmp_path,
        "a: 1\nb: 0\nflag_true: true\nflag_false: false\nflag_yes: yes\n",
    )
    result = _parse_config_minimal(str(cfg_path))
    # Note: "yes" is not parsed as a boolean by the minimal parser.
    assert result == {
        "a": 1,
        "b": 0,
        "flag_true": True,
        "flag_false": False,
        "flag_yes": "yes",
    }


def test_parse_config_minimal_genome_resources(tmp_path):
    from encode_pipeline.config.validator import _parse_config_minimal

    cfg_path, _ = _make_minimal_config(
        tmp_path,
        "genome_resources:\n"
        "  hs:\n"
        "    effective_genome_size: hs\n"
        '    chrom_sizes: ""\n'
        "  mm:\n"
        "    effective_genome_size: mm\n",
    )
    result = _parse_config_minimal(str(cfg_path))
    assert result == {
        "genome_resources": {
            "hs": {
                "effective_genome_size": "hs",
                "chrom_sizes": "",
            },
            "mm": {
                "effective_genome_size": "mm",
            },
        }
    }


def test_parse_config_minimal_qc_block(tmp_path):
    from encode_pipeline.config.validator import _parse_config_minimal

    cfg_path, _ = _make_minimal_config(
        tmp_path,
        "qc:\n  blacklist_filter: true\n  cross_correlation: false\n",
    )
    result = _parse_config_minimal(str(cfg_path))
    assert result == {
        "qc": {
            "blacklist_filter": True,
            "cross_correlation": False,
        }
    }


def test_parse_config_minimal_tool_parameters(tmp_path):
    from encode_pipeline.config.validator import _parse_config_minimal

    cfg_path, _ = _make_minimal_config(
        tmp_path,
        "tool_parameters:\n  macs3:\n    qvalue: 0.01\n    broad_cutoff: 0.1\n",
    )
    result = _parse_config_minimal(str(cfg_path))
    assert result == {
        "tool_parameters": {
            "macs3": {
                "qvalue": "0.01",
                "broad_cutoff": "0.1",
            },
        }
    }


def test_parse_config_minimal_cuttag_seacr(tmp_path):
    from encode_pipeline.config.validator import _parse_config_minimal

    cfg_path, _ = _make_minimal_config(
        tmp_path,
        "cuttag:\n"
        "  peak_caller: macs3\n"
        "  seacr:\n"
        "    enabled: false\n"
        "    threshold: 0.01\n",
    )
    result = _parse_config_minimal(str(cfg_path))
    assert result == {
        "cuttag": {
            "peak_caller": "macs3",
            "seacr": {
                "enabled": False,
                "threshold": "0.01",
            },
        }
    }


def test_parse_config_minimal_ignores_comments_and_blanks(tmp_path):
    from encode_pipeline.config.validator import _parse_config_minimal

    cfg_path, _ = _make_minimal_config(
        tmp_path,
        "\n# comment\nthreads: 4\n\n  # indented comment\nuse_control: true\n",
    )
    result = _parse_config_minimal(str(cfg_path))
    assert result == {"threads": 4, "use_control": True}


def test_parse_config_minimal_unsupported_list_line(tmp_path):
    """Top-level list lines are not represented in the output dict.

    The current parser treats ``items:`` as a top-level key with value ``""``
    and ignores the following ``- a`` line because it has no active section.
    """
    from encode_pipeline.config.validator import _parse_config_minimal

    cfg_path, _ = _make_minimal_config(
        tmp_path,
        "items:\n  - a\n  - b\n",
    )
    result = _parse_config_minimal(str(cfg_path))
    assert result == {"items": ""}


def test_parse_config_minimal_nested_list_ignored_in_section(tmp_path):
    """List lines inside a known section are ignored by the current parser."""
    from encode_pipeline.config.validator import _parse_config_minimal

    cfg_path, _ = _make_minimal_config(
        tmp_path,
        "qc:\n  blacklist_filter: true\n  - ignored_item\n",
    )
    result = _parse_config_minimal(str(cfg_path))
    assert result == {"qc": {"blacklist_filter": True}}


def test_parse_config_minimal_malformed_indent(tmp_path):
    """Unexpected indent inside a top-level section returns the partial dict."""
    from encode_pipeline.config.validator import _parse_config_minimal

    cfg_path, _ = _make_minimal_config(
        tmp_path,
        "threads: 4\n    unexpected_indent: value\nuse_control: true\n",
    )
    result = _parse_config_minimal(str(cfg_path))
    # The indented line has no active section, so it is skipped.
    assert result == {"threads": 4, "use_control": True}


# ---------------------------------------------------------------------------
# _load_yaml subprocess tests
# ---------------------------------------------------------------------------


def _run_in_subprocess(code, extra_pythonpath=""):
    """Run Python code in a subprocess with PYTHONPATH=extra:src."""
    repo_src = str(Path(__file__).resolve().parents[2] / "src")
    pythonpath = repo_src
    if extra_pythonpath:
        pythonpath = f"{extra_pythonpath}{os.pathsep}{repo_src}"
    env = os.environ.copy()
    env["PYTHONPATH"] = pythonpath
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    proc = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        env=env,
    )
    return proc


def test_load_yaml_with_pyyaml_available(tmp_path):
    """When PyYAML is present, _load_yaml returns None for an empty file."""
    pytest.importorskip("yaml")

    cfg_path = tmp_path / "empty.yaml"
    cfg_path.write_text("", encoding="utf-8")

    code = (
        "from encode_pipeline.config.validator import _load_yaml\n"
        f"print(repr(_load_yaml({str(cfg_path)!r})))"
    )
    proc = _run_in_subprocess(code)
    assert proc.returncode == 0, proc.stderr
    # PyYAML returns None for empty files; fallback returns {}.
    assert proc.stdout.strip() == "None"


def test_load_yaml_pyyaml_parses_inline_lists(tmp_path):
    """PyYAML parses inline list syntax; fallback does not."""
    pytest.importorskip("yaml")

    cfg_path = tmp_path / "list.yaml"
    cfg_path.write_text("items: [a, b]\n", encoding="utf-8")

    code = (
        "from encode_pipeline.config.validator import _load_yaml\n"
        f"print(repr(_load_yaml({str(cfg_path)!r})))"
    )
    proc = _run_in_subprocess(code)
    assert proc.returncode == 0, proc.stderr
    assert "['a', 'b']" in proc.stdout


def test_load_yaml_fallback_without_pyyaml(tmp_path):
    """When PyYAML cannot be imported, _load_yaml falls back to the stdlib parser."""
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(
        "threads: 4\n"
        "use_control: true\n"
        "genome_resources:\n"
        "  hs:\n"
        "    effective_genome_size: hs\n",
        encoding="utf-8",
    )

    fake_yaml = tmp_path / "fake_yaml"
    fake_yaml.mkdir()
    (fake_yaml / "yaml.py").write_text(
        "raise ImportError('simulated missing PyYAML')\n",
        encoding="utf-8",
    )

    code = (
        "from encode_pipeline.config.validator import _load_yaml\n"
        f"print(repr(_load_yaml({str(cfg_path)!r})))"
    )
    proc = _run_in_subprocess(code, extra_pythonpath=str(fake_yaml))
    assert proc.returncode == 0, proc.stderr
    assert "threads" in proc.stdout
    assert "genome_resources" in proc.stdout
    assert "hs" in proc.stdout


def test_load_yaml_empty_file_fallback(tmp_path):
    """Fallback parser returns {} for an empty file."""
    cfg_path = tmp_path / "empty.yaml"
    cfg_path.write_text("", encoding="utf-8")

    fake_yaml = tmp_path / "fake_yaml"
    fake_yaml.mkdir()
    (fake_yaml / "yaml.py").write_text(
        "raise ImportError('simulated missing PyYAML')\n",
        encoding="utf-8",
    )

    code = (
        "from encode_pipeline.config.validator import _load_yaml\n"
        f"print(repr(_load_yaml({str(cfg_path)!r})))"
    )
    proc = _run_in_subprocess(code, extra_pythonpath=str(fake_yaml))
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "{}"


# ---------------------------------------------------------------------------
# config.validator.main tests
# ---------------------------------------------------------------------------


def test_validator_main_success_returns_none(tmp_path, capsys):
    from encode_pipeline.config.validator import main

    cfg_path, _ = _make_minimal_config(tmp_path, "")
    result = main(["--config", str(cfg_path)])
    assert result is None

    captured = capsys.readouterr()
    assert "OK:" in captured.out
    assert "use_control:" in captured.out


def test_validator_main_missing_config_exits_nonzero(tmp_path):
    cfg_path = tmp_path / "missing.yaml"
    code = (
        "from encode_pipeline.config.validator import main\n"
        "import sys\n"
        f"sys.exit(main(['--config', {str(cfg_path)!r}]))"
    )
    proc = _run_in_subprocess(code)
    assert proc.returncode == 1
    assert "ERROR: Config file not found" in proc.stderr


def test_validator_main_invalid_config_exits_nonzero(tmp_path):
    cfg_path, _ = _make_minimal_config(
        tmp_path,
        "threads: not_an_int\n",
    )
    code = (
        "from encode_pipeline.config.validator import main\n"
        "import sys\n"
        f"sys.exit(main(['--config', {str(cfg_path)!r}]))"
    )
    proc = _run_in_subprocess(code)
    assert proc.returncode == 1
    assert "ERROR:" in proc.stderr


# ---------------------------------------------------------------------------
# Wrapper / package CLI equivalence tests
# ---------------------------------------------------------------------------


def test_cli_validate_main_emits_structured_log(tmp_path):
    cfg_path, _ = _make_minimal_config(tmp_path, "")
    code = (
        "import sys\n"
        "from encode_pipeline.cli.validate import main\n"
        f"sys.exit(main(['--config', {str(cfg_path)!r}]))"
    )
    proc = _run_in_subprocess(code)
    assert proc.returncode == 0, proc.stderr
    combined = proc.stdout + proc.stderr
    assert "encode-pipeline" in combined
    assert "OK:" in combined
    assert "use_control:" in combined


def test_scripts_validate_samples_wrapper(tmp_path):
    cfg_path, _ = _make_minimal_config(tmp_path, "")
    repo_root = Path(__file__).resolve().parents[2]
    script = repo_root / "scripts" / "validate_samples.py"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(repo_root / "src")
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    proc = subprocess.run(
        [sys.executable, str(script), "--config", str(cfg_path)],
        capture_output=True,
        text=True,
        env=env,
    )
    assert proc.returncode == 0, proc.stderr
    combined = proc.stdout + proc.stderr
    assert "encode-pipeline" in combined
    assert "OK:" in combined


def test_encode_validate_console_script_if_installed(tmp_path):
    import shutil

    exe = shutil.which("encode-validate")
    if exe is None:
        pytest.skip("encode-validate console script not installed")

    cfg_path, _ = _make_minimal_config(tmp_path, "")
    env = os.environ.copy()
    env["PYTHONPATH"] = str(Path(__file__).resolve().parents[2] / "src")
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    proc = subprocess.run(
        [exe, "--config", str(cfg_path)],
        capture_output=True,
        text=True,
        env=env,
    )
    assert proc.returncode == 0, proc.stderr
    combined = proc.stdout + proc.stderr
    assert "encode-pipeline" in combined
    assert "OK:" in combined
