"""Make this checkout's sources win over another installed copy."""

from __future__ import annotations

import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
_SRC = str(_ROOT / "src")
_BROWSER = str(Path(__file__).resolve().parent)
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
if _BROWSER not in sys.path:
    sys.path.append(_BROWSER)
