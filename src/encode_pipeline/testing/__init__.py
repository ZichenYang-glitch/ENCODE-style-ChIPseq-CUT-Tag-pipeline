"""Reusable test support for independently packaged workflow adapters."""

from encode_pipeline.testing.adapter_conformance import (
    AdapterConformanceCase,
    AdapterConformanceError,
    verify_adapter_conformance,
)

__all__ = [
    "AdapterConformanceCase",
    "AdapterConformanceError",
    "verify_adapter_conformance",
]
