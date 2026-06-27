"""Primitive config value coercion helpers."""

__all__ = ["coerce_bool", "coerce_int"]


def coerce_int(value, *, name: str, minimum: int, error_cls=ValueError) -> int:
    """Return a strictly parsed integer, rejecting bools and floats."""
    if isinstance(value, bool):
        raise error_cls(f"config {name} must be an integer, got {value!r}")
    if isinstance(value, int):
        parsed = value
    elif isinstance(value, str):
        text = value.strip()
        if not text.isdigit():
            raise error_cls(f"config {name} must be an integer, got {value!r}")
        parsed = int(text)
    else:
        raise error_cls(f"config {name} must be an integer, got {value!r}")

    if parsed < minimum:
        if minimum == 1:
            raise error_cls(f"config {name} must be positive, got {parsed}")
        raise error_cls(f"config {name} must be non-negative, got {parsed}")
    return parsed


def coerce_bool(value, name: str, error_cls=ValueError) -> bool:
    """Coerce a value to bool. Accepts bool or string bool."""
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value.lower() in ("true", "false"):
        return value.lower() == "true"
    raise error_cls(f"{name} must be true or false, got {value!r}")
