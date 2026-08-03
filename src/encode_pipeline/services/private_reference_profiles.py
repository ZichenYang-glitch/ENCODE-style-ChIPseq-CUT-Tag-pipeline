"""Descriptor-safe loader for operator-private Reference Profile bindings."""

from __future__ import annotations

import json
import math
from pathlib import Path
import re
from types import MappingProxyType
from typing import Mapping

import yaml

from encode_pipeline.services.private_storage_pools import (
    _read_bounded_regular_file,
)


PRIVATE_REFERENCE_PROFILE_SCHEMA_VERSION = "helixweave-reference-profiles-v1"

_MAX_CONFIG_BYTES = 64 * 1024
_MAX_PROFILES = 256
_MAX_BINDINGS_PER_PROFILE = 32
_MAX_NODES = 4096
_MAX_DEPTH = 24
_MAX_STRING_BYTES = 16 * 1024
_SAFE_KEY = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]{0,254}$")
_WORKFLOW_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:+/-]{0,254}$")
_PUBLIC_MESSAGE = "private reference profile configuration is invalid"


class PrivateReferenceProfileConfigError(RuntimeError):
    """Stable redacted failure that never reflects paths or private payloads."""

    reason_code = "REFERENCE_PROFILE_CONFIG_INVALID"

    def __init__(self) -> None:
        super().__init__(_PUBLIC_MESSAGE)


class PrivateReferenceProfileConfig:
    """Bounded immutable mapping from config keys to adapter-private payloads."""

    __slots__ = ("_bindings",)

    def __init__(
        self,
        bindings: Mapping[str, Mapping[str, Mapping[str, object]]],
    ) -> None:
        self._bindings = MappingProxyType(
            {
                config_key: MappingProxyType(
                    {
                        workflow_id: _freeze(payload)
                        for workflow_id, payload in workflow_bindings.items()
                    }
                )
                for config_key, workflow_bindings in bindings.items()
            }
        )

    @classmethod
    def from_file(cls, config_path: Path) -> PrivateReferenceProfileConfig:
        try:
            content = _read_bounded_regular_file(config_path)
            if len(content) > _MAX_CONFIG_BYTES:
                raise ValueError("oversized")
            if config_path.suffix.lower() in {".yaml", ".yml"}:
                payload = _strict_yaml_loads(content)
            else:
                payload = _strict_json_loads(content)
            bindings = _validate_document(payload)
        except Exception:
            raise PrivateReferenceProfileConfigError() from None
        return cls(bindings)

    @property
    def config_keys(self) -> tuple[str, ...]:
        return tuple(sorted(self._bindings))

    def workflow_ids_for(self, config_key: str) -> tuple[str, ...]:
        try:
            return tuple(sorted(self._bindings[_checked_key(config_key)]))
        except (KeyError, ValueError):
            raise PrivateReferenceProfileConfigError() from None

    def binding_for(
        self,
        config_key: str,
        workflow_id: str,
    ) -> Mapping[str, object]:
        try:
            payload = self._bindings[_checked_key(config_key)][
                _checked_workflow_id(workflow_id)
            ]
        except (KeyError, ValueError):
            raise PrivateReferenceProfileConfigError() from None
        thawed = _thaw(payload)
        if not isinstance(thawed, dict):
            raise PrivateReferenceProfileConfigError() from None
        return thawed

    def __repr__(self) -> str:
        binding_count = sum(len(bindings) for bindings in self._bindings.values())
        return (
            f"{type(self).__name__}(configured_profile_count={len(self._bindings)}, "
            f"binding_count={binding_count})"
        )


def load_private_reference_profile_config(
    config_path: Path,
) -> PrivateReferenceProfileConfig:
    return PrivateReferenceProfileConfig.from_file(config_path)


def _strict_json_loads(content: bytes) -> object:
    def object_pairs(pairs: list[tuple[str, object]]) -> dict[str, object]:
        result: dict[str, object] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError("duplicate")
            result[key] = value
        return result

    def reject_constant(_value: str) -> float:
        raise ValueError("constant")

    try:
        return json.loads(
            content.decode("utf-8"),
            object_pairs_hook=object_pairs,
            parse_constant=reject_constant,
        )
    except (UnicodeError, ValueError, json.JSONDecodeError, RecursionError):
        raise ValueError("document") from None


class _UniqueSafeLoader(yaml.SafeLoader):
    pass


def _construct_unique_mapping(
    loader: _UniqueSafeLoader,
    node: yaml.MappingNode,
    deep: bool = False,
) -> dict[object, object]:
    loader.flatten_mapping(node)
    result: dict[object, object] = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if not isinstance(key, str) or key in result:
            raise ValueError("mapping key")
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


_UniqueSafeLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
    _construct_unique_mapping,
)


def _strict_yaml_loads(content: bytes) -> object:
    try:
        document = content.decode("utf-8")
        documents = list(yaml.load_all(document, Loader=_UniqueSafeLoader))
    except (UnicodeError, ValueError, yaml.YAMLError, RecursionError):
        raise ValueError("document") from None
    if len(documents) != 1:
        raise ValueError("document count")
    return documents[0]


def _validate_document(
    payload: object,
) -> dict[str, dict[str, dict[str, object]]]:
    _validate_json_value(payload, depth=0, counter=[0], ancestors=set())
    if not isinstance(payload, dict) or set(payload) != {
        "schema_version",
        "profiles",
    }:
        raise ValueError("shape")
    if payload["schema_version"] != PRIVATE_REFERENCE_PROFILE_SCHEMA_VERSION:
        raise ValueError("version")
    profiles = payload["profiles"]
    if not isinstance(profiles, dict) or not profiles or len(profiles) > _MAX_PROFILES:
        raise ValueError("profiles")
    result: dict[str, dict[str, dict[str, object]]] = {}
    for config_key, raw_profile in profiles.items():
        config_key = _checked_key(config_key)
        if not isinstance(raw_profile, dict) or set(raw_profile) != {"bindings"}:
            raise ValueError("profile")
        raw_bindings = raw_profile["bindings"]
        if (
            not isinstance(raw_bindings, dict)
            or not raw_bindings
            or len(raw_bindings) > _MAX_BINDINGS_PER_PROFILE
        ):
            raise ValueError("bindings")
        workflow_bindings: dict[str, dict[str, object]] = {}
        for workflow_id, raw_binding in raw_bindings.items():
            workflow_id = _checked_workflow_id(workflow_id)
            if not isinstance(raw_binding, dict):
                raise ValueError("binding")
            workflow_bindings[workflow_id] = dict(raw_binding)
        result[config_key] = workflow_bindings
    return result


def _validate_json_value(
    value: object,
    *,
    depth: int,
    counter: list[int],
    ancestors: set[int],
) -> None:
    counter[0] += 1
    if counter[0] > _MAX_NODES or depth > _MAX_DEPTH:
        raise ValueError("complexity")
    if value is None or isinstance(value, bool):
        return
    if isinstance(value, str):
        if len(value.encode("utf-8")) > _MAX_STRING_BYTES or any(
            ord(character) < 32 and character not in "\t\n\r" for character in value
        ):
            raise ValueError("string")
        return
    if isinstance(value, int) and not isinstance(value, bool):
        return
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("number")
        return
    if isinstance(value, (dict, list)):
        identity = id(value)
        if identity in ancestors:
            raise ValueError("cycle")
        ancestors.add(identity)
        try:
            if isinstance(value, dict):
                for key, item in value.items():
                    if not isinstance(key, str):
                        raise ValueError("key")
                    _validate_json_value(
                        key,
                        depth=depth + 1,
                        counter=counter,
                        ancestors=ancestors,
                    )
                    _validate_json_value(
                        item,
                        depth=depth + 1,
                        counter=counter,
                        ancestors=ancestors,
                    )
            else:
                for item in value:
                    _validate_json_value(
                        item,
                        depth=depth + 1,
                        counter=counter,
                        ancestors=ancestors,
                    )
        finally:
            ancestors.remove(identity)
        return
    raise ValueError("type")


def _checked_key(value: object) -> str:
    if not isinstance(value, str) or _SAFE_KEY.fullmatch(value) is None:
        raise ValueError("config key")
    return value


def _checked_workflow_id(value: object) -> str:
    if not isinstance(value, str) or _WORKFLOW_ID.fullmatch(value) is None:
        raise ValueError("workflow id")
    return value


def _freeze(value: object) -> object:
    if isinstance(value, dict):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, list):
        return tuple(_freeze(item) for item in value)
    return value


def _thaw(value: object) -> object:
    if isinstance(value, Mapping):
        return {key: _thaw(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_thaw(item) for item in value]
    return value
