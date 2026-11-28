"""The computation configuration: target, bindings, assumptions, dimensions."""

from __future__ import annotations

import ctypes
from enum import IntEnum, IntFlag

from . import _ffi
from .errors import check

__all__ = ["Assumption", "Config", "DimensionMode"]


class Assumption(IntFlag):
    """The fact bits `rssn_config_assume` accepts.

    `POSITIVE` is `REAL | NONZERO | NONNEGATIVE` (7) and `NEGATIVE` is
    `REAL | NONZERO | NONPOSITIVE` (11), which is how the engine defines them.
    """

    REAL = 1
    NONZERO = 2
    NONNEGATIVE = 4
    NONPOSITIVE = 8
    INTEGER = 16

    POSITIVE = REAL | NONZERO | NONNEGATIVE
    NEGATIVE = REAL | NONZERO | NONPOSITIVE


class DimensionMode(IntEnum):
    """`rssn_config_dimension_mode`."""

    OFF = 0
    ADVISORY = 1
    STRICT = 2


_ASSUMPTION_WORDS = {
    "real": Assumption.REAL,
    "nonzero": Assumption.NONZERO,
    "nonnegative": Assumption.NONNEGATIVE,
    "nonpositive": Assumption.NONPOSITIVE,
    "integer": Assumption.INTEGER,
    "positive": Assumption.POSITIVE,
    "negative": Assumption.NEGATIVE,
}


def _facts(value) -> int:
    """An `Assumption`, an integer bit set, or a set of names as an int."""
    if isinstance(value, Assumption):
        return int(value)
    if isinstance(value, int):
        return value
    if isinstance(value, str):
        return int(_ASSUMPTION_WORDS[value.lower()])
    total = Assumption(0)
    for item in value:
        word = item.lower() if isinstance(item, str) else None
        total |= _ASSUMPTION_WORDS[word] if word else Assumption(item)
    return int(total)


class Config:
    """An owned `RssnConfig`.

    Every method mutates the configuration and returns `self`, so a target can
    be written in one expression::

        cfg = Config().numeric(1e-12).bind("x", 0.5).assume("n", Assumption.INTEGER)

    A `Config` is single-threaded, like the session it is used with, and is
    released by :meth:`close`, by the context manager, or when it is collected.
    """

    __slots__ = ("_ptr",)

    def __init__(self) -> None:
        pointer = _ffi.library().rssn_config_new()
        if not pointer:
            raise MemoryError(f"rssn_config_new failed: {_ffi.last_error()}")
        self._ptr = pointer

    # -- lifetime ---------------------------------------------------------
    @property
    def ptr(self) -> int:
        if not self._ptr:
            raise RuntimeError("this Config has been closed")
        return self._ptr

    def close(self) -> None:
        if self._ptr:
            _ffi.library().rssn_config_free(ctypes.c_void_p(self._ptr))
            self._ptr = 0

    def __enter__(self) -> "Config":
        return self

    def __exit__(self, *exc_info) -> None:
        self.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:  # interpreter shutdown, missing library
            pass

    # -- target -----------------------------------------------------------
    def symbolic(self) -> "Config":
        """Ask for a closed form (the default)."""
        check(_ffi.library().rssn_config_symbolic(ctypes.c_void_p(self.ptr)))
        return self

    def numeric(self, tolerance: float) -> "Config":
        """Ask for a numeric answer within `tolerance`."""
        check(_ffi.library().rssn_config_numeric(ctypes.c_void_p(self.ptr), float(tolerance)))
        return self

    # -- environment ------------------------------------------------------
    def bind(self, name: str, value: float) -> "Config":
        """Bind the symbol `name` to `value` for numeric evaluation."""
        check(
            _ffi.library().rssn_config_bind(
                ctypes.c_void_p(self.ptr), name.encode("utf-8"), float(value)
            )
        )
        return self

    def bind_all(self, bindings: dict[str, float]) -> "Config":
        for name, value in bindings.items():
            self.bind(name, value)
        return self

    def assume(self, name: str, facts) -> "Config":
        """Declare facts about `name`: an :class:`Assumption` or names such as `"positive"`."""
        check(
            _ffi.library().rssn_config_assume(
                ctypes.c_void_p(self.ptr), name.encode("utf-8"), ctypes.c_uint8(_facts(facts))
            )
        )
        return self

    def budget(self, max_iterations: int = 0, max_nodes: int = 0) -> "Config":
        """Limit the search; zero keeps the current value."""
        check(
            _ffi.library().rssn_config_budget(
                ctypes.c_void_p(self.ptr), ctypes.c_size_t(max_iterations), ctypes.c_size_t(max_nodes)
            )
        )
        return self

    # -- dimensions -------------------------------------------------------
    def dimension_mode(self, mode) -> "Config":
        """Set :class:`DimensionMode` (`"off"`, `"advisory"`, `"strict"` are accepted)."""
        if isinstance(mode, str):
            mode = DimensionMode[mode.upper()]
        check(
            _ffi.library().rssn_config_dimension_mode(
                ctypes.c_void_p(self.ptr), ctypes.c_uint8(int(mode))
            )
        )
        return self

    def bypass_dimensions(self) -> "Config":
        """Disable hard dimensional errors even in strict mode."""
        check(_ffi.library().rssn_config_bypass_dimensions(ctypes.c_void_p(self.ptr)))
        return self

    def dimension(self, symbol: str, text: str) -> "Config":
        """Declare the dimension of `symbol` as text (`"m/s"`, `"Length/Time"`, `"any"`)."""
        check(
            _ffi.library().rssn_config_dimension(
                ctypes.c_void_p(self.ptr), symbol.encode("utf-8"), text.encode("utf-8")
            )
        )
        return self

    def define_base_dimension(self, name: str) -> "Config":
        """Register a base dimension the engine does not have (for example `"Currency"`)."""
        check(
            _ffi.library().rssn_config_define_base_dimension(
                ctypes.c_void_p(self.ptr), name.encode("utf-8")
            )
        )
        return self

    def __repr__(self) -> str:
        return f"<Config {'closed' if not self._ptr else hex(self._ptr)}>"
