"""Public CLI contracts for the strict input-validation switch."""

from pathlib import Path

import pytest
import yaml

from encode_pipeline.cli.validate import main as validate_main


SAMPLE_HEADER = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
    "bowtie2_index\n"
)


def _placeholder_config(tmp_path: Path) -> Path:
    samples_path = tmp_path / "samples.tsv"
    samples_path.write_text(
        SAMPLE_HEADER
        + (
            f"S1\t{tmp_path / 'missing_R1.fastq'}\t{tmp_path / 'missing_R2.fastq'}"
            f"\tPE\tchipseq\tCTCF\tnarrow\ths\t{tmp_path / 'missing_index'}\n"
        ),
        encoding="utf-8",
    )
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "samples": str(samples_path),
                "use_control": False,
            },
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return config_path


def test_public_validate_cli_forwards_strict_inputs_and_rejects_missing_files(
    tmp_path, capsys
):
    config_path = _placeholder_config(tmp_path)

    with pytest.raises(SystemExit) as exit_info:
        validate_main(["--config", str(config_path), "--strict-inputs"])

    captured = capsys.readouterr()
    assert exit_info.value.code == 1
    assert "fastq_1 file not found" in captured.err
    assert str(tmp_path) not in captured.out


def test_public_validate_cli_keeps_placeholder_paths_in_non_strict_mode(
    tmp_path, capsys
):
    config_path = _placeholder_config(tmp_path)

    result = validate_main(["--config", str(config_path)])

    captured = capsys.readouterr()
    assert result is None
    assert "OK: 1 sample(s) validated" in captured.out
    assert captured.err == ""
