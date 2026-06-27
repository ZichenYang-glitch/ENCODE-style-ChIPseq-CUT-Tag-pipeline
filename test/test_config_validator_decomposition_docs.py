"""Guard that the validator decomposition spec doc stays in sync.

This test does not import or execute the validator; it only checks that the
policy document exists and references the key files and concepts it claims to
describe.
"""

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
DOCS_PATH = REPO_ROOT / "docs" / "development" / "config-validator-decomposition.md"


REQUIRED_MENTIONS = (
    "validate_config",
    "load_and_validate_samples",
    "ValidationError",
    "scripts/validate_samples.py",
    "src/encode_pipeline/config/validator.py",
    "src/encode_pipeline/config/validate.py",
    "src/encode_pipeline/samples/load.py",
    "src/encode_pipeline/samples/validate.py",
    "src/encode_pipeline/cli/validate.py",
    "test/config/test_validation.py",
    "test/test_stage8_smoke_profiles.py",
    "defaults.py",
    "Recommended PR52 first extraction",
    "must not introduce an `encode_pipeline.config.cli`",
)


def test_validator_decomposition_doc_exists():
    assert DOCS_PATH.is_file(), f"Missing validator decomposition doc: {DOCS_PATH}"


def test_validator_decomposition_doc_mentions_required_concepts():
    text = DOCS_PATH.read_text(encoding="utf-8")
    missing = [m for m in REQUIRED_MENTIONS if m not in text]
    assert not missing, f"Doc missing required references: {missing}"
