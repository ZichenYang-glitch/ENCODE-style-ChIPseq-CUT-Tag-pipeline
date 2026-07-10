"""Tests for the ENCODE-style workflow adapter."""

from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter


def test_encode_adapter_declares_workspace_plan_capability():
    adapter = EncodeStyleWorkflowAdapter()
    assert "workspace_plan" in adapter.capabilities.supports
