"""Terms: handles into a session's expression graph."""

from __future__ import annotations

import ctypes
from fractions import Fraction
from typing import TYPE_CHECKING, Any

from . import _ffi
from .errors import RssnError, Status, check

if TYPE_CHECKING:  # pragma: no cover
    from .answer import Answer
    from .config import Config
    from .session import Session

__all__ = ["Term", "coerce"]


class Term:
    """A node of a :class:`~rssn.session.Session` graph.

    A term is a `uint32_t` handle valid for the session that minted it; the
    session owns it and releases it with itself. Terms are immutable values:
    `x * x` makes a new node, it does not change `x`.

    Operators are *requests*, not results. `rssn.diff(f, x)` returns the
    unevaluated derivative; :meth:`compute` reduces it::

        >>> import rssn
        >>> s = rssn.Session()
        >>> x = s.sym("x")
        >>> str(rssn.diff(rssn.sin(x) * rssn.cos(x), x))
        'diff(sin(x)*cos(x), x)'
        >>> s.compute(rssn.diff(rssn.sin(x) * rssn.cos(x), x)).text
        'cos(2*x)'
    """

    __slots__ = ("_session", "_handle")

    def __init__(self, session: "Session", handle: int, status: Status = Status.UNKNOWN_OPERATOR) -> None:
        # The C constructors report failure as `RSSN_NO_TERM` plus a message;
        # the code is a caller-supplied hint about what was attempted.
        if handle == _ffi.RSSN_NO_TERM:
            raise RssnError(status, _ffi.last_error() or "no term")
        self._session = session
        self._handle = int(handle)

    # -- identity ---------------------------------------------------------
    @property
    def session(self) -> "Session":
        """The session that owns this handle."""
        return self._session

    @property
    def handle(self) -> int:
        """The raw `uint32_t` handle."""
        return self._handle

    def __repr__(self) -> str:
        return f"<Term {self._handle} {self.text!r}>"

    def __str__(self) -> str:
        return self.text

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Term):
            return self._handle == other._handle and self._session is other._session
        return NotImplemented

    def __hash__(self) -> int:
        return hash((id(self._session), self._handle))

    # -- renderings -------------------------------------------------------
    @property
    def text(self) -> str:
        """The term as infix text, parseable by :meth:`Session.parse`."""
        return _ffi.take_string(_ffi.library().rssn_term_to_string(self._ptr(), self._handle))

    @property
    def latex(self) -> str:
        """The term as LaTeX math."""
        return _ffi.take_string(_ffi.library().rssn_term_to_latex(self._ptr(), self._handle))

    @property
    def typst(self) -> str:
        """The term as Typst math."""
        return _ffi.take_string(_ffi.library().rssn_term_to_typst(self._ptr(), self._handle))

    @property
    def pretty(self) -> str:
        """The term as a multi-line Unicode drawing."""
        return _ffi.take_string(_ffi.library().rssn_term_to_pretty(self._ptr(), self._handle))

    # -- values -----------------------------------------------------------
    @property
    def value(self) -> float | None:
        """The value of a literal number, or `None` for a symbolic term."""
        out = ctypes.c_double()
        code = _ffi.library().rssn_term_to_f64(self._ptr(), self._handle, ctypes.byref(out))
        if int(code) == int(Status.OK):
            return out.value
        if int(code) == int(Status.NO_VALUE):
            return None
        raise RssnError(Status(int(code)), _ffi.last_error())

    def tensor(self) -> Any:
        """A dense `f64` extraction of a list term.

        Returns a `numpy.ndarray` when numpy is installed and a nested list
        otherwise. A scalar term gives a 0-d array (or a float).
        """
        tensor = _ffi.RssnTensor()
        check(_ffi.library().rssn_term_tensor_data(self._ptr(), self._handle, ctypes.byref(tensor)))
        try:
            shape = [int(tensor.shape[index]) for index in range(tensor.rank)]
            data = [float(tensor.data[index]) for index in range(tensor.len)]
        finally:
            _ffi.library().rssn_tensor_free(ctypes.byref(tensor))
        try:  # numpy is optional
            import numpy as np  # noqa: PLC0415

            return np.array(data, dtype=float).reshape(shape) if shape else np.array(data[0])
        except ImportError:
            if not shape:
                return data[0] if data else 0.0
            return _nest(data, shape)

    # -- computation ------------------------------------------------------
    def compute(self, config: "Config | None" = None) -> "Answer":
        """Reduce this term to the target of `config` (symbolic by default)."""
        return self._session.compute(self, config)

    def simplify(self) -> "Answer":
        """Reduce this term with the default symbolic configuration."""
        return self._session.compute(self)

    def eval(self, **bindings: float) -> float:
        """Evaluate numerically without running any rule."""
        return self._session.eval(self, bindings)

    def to_float(self) -> float:
        """The numeric value of the term, computing it if necessary."""
        value = self.value
        if value is not None:
            return value
        answer = self.compute(self._session.numeric_config())
        if answer.value is None:
            raise RssnError(Status.NOT_NUMERIC, f"{self.text!r} has no numeric value")
        return answer.value

    def __float__(self) -> float:
        return self.to_float()

    # -- arithmetic, as requests ------------------------------------------
    def __add__(self, other: Any) -> "Term":
        return self._session.apply("add", self, other)

    def __radd__(self, other: Any) -> "Term":
        return self._session.apply("add", other, self)

    def __sub__(self, other: Any) -> "Term":
        return self._session.apply("sub", self, other)

    def __rsub__(self, other: Any) -> "Term":
        return self._session.apply("sub", other, self)

    def __mul__(self, other: Any) -> "Term":
        return self._session.apply("mul", self, other)

    def __rmul__(self, other: Any) -> "Term":
        return self._session.apply("mul", other, self)

    def __truediv__(self, other: Any) -> "Term":
        return self._session.apply("divide", self, other)

    def __rtruediv__(self, other: Any) -> "Term":
        return self._session.apply("divide", other, self)

    def __pow__(self, other: Any, modulo: Any = None) -> "Term":
        if modulo is not None:
            return NotImplemented
        return self._session.apply("pow", self, other)

    def __rpow__(self, other: Any) -> "Term":
        return self._session.apply("pow", other, self)

    def __neg__(self) -> "Term":
        return self._session.apply("neg", self)

    def __pos__(self) -> "Term":
        return self

    def __abs__(self) -> "Term":
        return self._session.apply("abs", self)

    # -- helpers ----------------------------------------------------------
    # `sub`, `neg` and `divide` are the engine's named forms; it rewrites them
    # to `add(a, mul(-1, b))` and `mul(a, pow(b, -1))` at its normalize step and
    # prints both spellings the same way. The arithmetic quotient is `divide`
    # because `div` is the divergence of a vector field.
    def _ptr(self) -> int:
        return self._session.ptr

    def _handle_of(self, value: Any) -> int:
        return coerce(value, self._session)


def _nest(data: list[float], shape: list[int]) -> Any:
    """Reshape a flat list into nested lists of `shape`."""
    if len(shape) == 1:
        return list(data)
    step = len(data) // shape[0]
    return [_nest(data[index * step : (index + 1) * step], shape[1:]) for index in range(shape[0])]


def coerce(value: Any, session: "Session") -> int:
    """A term handle for `value` in `session`.

    `Term` handles must belong to `session`. Numbers become exact integers or
    floats, `Fraction` becomes an exact rational, and `str` is parsed as infix
    text (so `"x"` is the symbol `x`).
    """
    if isinstance(value, Term):
        if value._session is not session:
            raise ValueError("term belongs to another session")
        return value._handle
    library = _ffi.library()
    if isinstance(value, bool):
        return int(library.rssn_term_int(session.ptr, int(value)))
    if isinstance(value, int):
        return int(library.rssn_term_int(session.ptr, value))
    if isinstance(value, Fraction):
        return int(library.rssn_term_rational(session.ptr, value.numerator, value.denominator))
    if isinstance(value, float):
        return int(library.rssn_term_float(session.ptr, value))
    if isinstance(value, str):
        handle = int(library.rssn_term_parse(session.ptr, value.encode("utf-8")))
        if handle == _ffi.RSSN_NO_TERM:
            raise RssnError(Status.PARSE, _ffi.last_error() or f"cannot parse {value!r}")
        return handle
    # numpy scalars and anything else with __float__
    if hasattr(value, "item"):
        item = value.item()
        if isinstance(item, (int, float)):
            return coerce(item, session)
    raise TypeError(f"cannot use {type(value).__name__} as a term")
