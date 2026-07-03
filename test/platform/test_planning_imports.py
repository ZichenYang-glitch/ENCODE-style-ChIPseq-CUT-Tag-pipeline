"""Import-boundary tests for planning modules."""

import sys


def test_planning_platform_does_not_import_forbidden_modules():
    before = set(sys.modules.keys())
    import encode_pipeline.platform.planning  # noqa: F401

    imported = set(sys.modules.keys()) - before
    forbidden = {"fastapi", "pydantic", "snakemake", "subprocess"}
    found = forbidden & imported
    assert not found, f"platform.planning imported forbidden modules: {found}"


def test_planning_service_does_not_import_forbidden_modules():
    before = set(sys.modules.keys())
    import encode_pipeline.services.planning  # noqa: F401

    imported = set(sys.modules.keys()) - before
    forbidden = {"fastapi", "pydantic", "snakemake", "subprocess"}
    found = forbidden & imported
    assert not found, f"services.planning imported forbidden modules: {found}"
