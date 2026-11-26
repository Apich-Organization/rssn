"""Host callbacks: Python callables the engine can call.

Some operators are *higher-order*: they take a function and call it, once per
quadrature node, per optimisation step, per root-finder probe. `num_bfgs`,
`num_brent`, `num_analysis_integrate_adaptive_scaled`, `num_bootstrap` and about
180 others are like this. In Rust they take a closure; in Python they take any
callable, because this module turns one into the C function pointer the engine
knows how to call.

    import rssn

    s = rssn.Session()
    # root of x**2 - 2 in [1, 2]
    answer = s.call_value(
        "num_brent", 1.0, 2.0, 1e-12, 64,
        lambda x: x * x - 2.0,
    )
    print(answer["root"])            # 1.4142135623730951

The engine hands the callable native Python objects — `float`, `int`, `bool`,
`complex`, `str`, `numpy` arrays (or nested lists without numpy) and dicts for
records — and converts whatever it returns back into a value. Returning a
`rssn.Value` or a `rssn.Term` is also allowed.

A `Callback` can be created explicitly when it should outlive one call, and it
is a context manager:

    with rssn.Callback(s, lambda x: x * x, signature="(float) -> float") as f:
        ...  # pass `f` as an argument, repeatedly
"""

from __future__ import annotations

import ctypes
import json
from typing import TYPE_CHECKING, Any, Callable

from . import _ffi
from .errors import RssnError, Status, check
from .value import Value

if TYPE_CHECKING:  # pragma: no cover
    from .session import Session
    from .term import Term

__all__ = ["Callback", "callback_signatures", "parse_signature"]

#: The status a callback returns when the Python function raised.
#: `RssnStatus::InvalidArgument` — the enclosing call raises `RssnError` with
#: the Python exception's text in the message.
FAILED = int(Status.INVALID_ARGUMENT)

#: Scalar spellings the compact signature grammar accepts.
_SCALARS = {
    "int": "int",
    "integer": "int",
    "real": "real",
    "rational": "real",
    "float": "float",
    "double": "float",
    "complex": "complex",
    "bool": "bool",
    "boolean": "bool",
    "string": "string",
    "str": "string",
    "text": "string",
    "symbol": "symbol",
}


def _split_top(text: str) -> list[str]:
    """Split `text` on commas that are not inside brackets."""
    parts: list[str] = []
    depth = 0
    current: list[str] = []
    for character in text:
        if character in "([{":
            depth += 1
        elif character in ")]}":
            depth -= 1
        if character == "," and depth == 0:
            parts.append("".join(current))
            current = []
        else:
            current.append(character)
    if "".join(current).strip():
        parts.append("".join(current))
    return [part.strip() for part in parts]


def _type_expr(token: str) -> Any:
    """One compact type (`real`, `float[]`, `int[][]`) as a `TypeExpr` JSON."""
    token = token.strip()
    rank = 0
    while token.endswith("[]"):
        rank += 1
        token = token[:-2].strip()
    base = _SCALARS.get(token.lower())
    if base is None:
        raise ValueError(
            f"unknown type {token!r} in a callback signature; use one of "
            f"{', '.join(sorted(set(_SCALARS.values())))} or append [] for a tensor"
        )
    if rank == 0:
        return base
    return {"tensor": {"rank": {"lit": rank}, "shape": None, "elem": base}}


def parse_signature(text: str) -> dict:
    """`"(float, float) -> float"` as the JSON object the C ABI reads.

    The arrow may be `->` or `=>`; the parameter list may be omitted for a
    nullary function (`"() -> bool"`). A trailing `[]` makes a tensor of that
    rank: `"float[]"` is a vector, `"float[][]"` a matrix. A string that is
    already JSON is passed through unchanged.
    """
    stripped = text.strip()
    if stripped.startswith(("{", "[")):
        return json.loads(stripped)
    arrow = "->" if "->" in stripped else "=>"
    if arrow not in stripped:
        raise ValueError(f"a callback signature needs `->`: {text!r}")
    left, right = stripped.split(arrow, 1)
    left = left.strip()
    if left.startswith("(") and left.endswith(")"):
        left = left[1:-1]
    params = [_type_expr(part) for part in _split_top(left)] if left.strip() else []
    return {"params": params, "ret": _type_expr(right)}


def _number(payload: Any) -> float:
    """A `Number` JSON (`{"Float": 1.0}`, `{"Int": 3}`, `{"Rat": [1, 3]}`)."""
    if isinstance(payload, (int, float)):
        return float(payload)
    if isinstance(payload, dict) and len(payload) == 1:
        tag, body = next(iter(payload.items()))
        if tag in ("Float", "float"):
            return float(body)
        if tag in ("Int", "int"):
            return float(body)
        if tag in ("Rat", "rat"):
            if isinstance(body, dict):
                return float(body["numer"]) / float(body["denom"])
            return float(body[0]) / float(body[1])
    raise ValueError(f"not a number: {payload!r}")


def to_python(session: "Session", cell: int) -> Any:
    """The Python object a value cell holds.

    This is what a callback receives for each of its arguments. A tensor becomes
    a `numpy` array when numpy is installed and a nested list otherwise; a
    record becomes a `dict`; a variant becomes `(tag, [payload...])`.
    """
    value = Value.borrowed(session, cell)
    kind = value.kind()
    if kind == "int":
        return int(value.float())
    if kind in ("real", "float"):
        return value.float()
    if kind == "bool":
        return bool(value.json()["bool"])
    if kind == "complex":
        body = value.json()["complex"]
        exact = body.get("Exact")
        if exact is not None:
            return complex(_number(exact["re"]), _number(exact["im"]))
        floating = body.get("Float", {})
        return complex(float(floating.get("re", 0.0)), float(floating.get("im", 0.0)))
    if kind in ("string", "symbol"):
        return value.text()
    if kind == "tensor":
        try:
            return value.tensor()
        except RssnError:
            return [to_python(session, element.cell) for element in value]
    if kind == "record":
        library = _ffi.library()
        fields: dict[str, Any] = {}
        for index in range(len(value)):
            raw = ctypes.c_void_p()
            check(
                library.rssn_value_field_name(
                    ctypes.c_void_p(session.ptr), cell, ctypes.c_size_t(index),
                    ctypes.byref(raw),
                )
            )
            fields[_ffi.take_string(raw.value)] = to_python(session, value[index].cell)
        return fields
    if kind == "variant":
        library = _ffi.library()
        raw = ctypes.c_void_p()
        check(
            library.rssn_value_tag(
                ctypes.c_void_p(session.ptr), cell, ctypes.byref(raw)
            )
        )
        tag = _ffi.take_string(raw.value)
        return tag, [to_python(session, element.cell) for element in value]
    return value.json()


def from_python(session: "Session", value: Any) -> int:
    """A cell id holding `value`, which the caller then owns."""
    from .term import Term

    library = _ffi.library()
    out = ctypes.c_uint32()
    if isinstance(value, Value):
        return value.cell
    if isinstance(value, Callback):
        return value.cell
    if isinstance(value, Term):
        if value.session is not session:
            raise ValueError("term belongs to another session")
        check(library.rssn_value_new_term(ctypes.c_void_p(session.ptr), value.handle,
                                          ctypes.byref(out)))
        return out.value
    if isinstance(value, bool):
        check(library.rssn_value_new_bool(ctypes.c_void_p(session.ptr), int(value),
                                          ctypes.byref(out)))
        return out.value
    if isinstance(value, int):
        check(library.rssn_value_new_int(ctypes.c_void_p(session.ptr), value,
                                         ctypes.byref(out)))
        return out.value
    if isinstance(value, float):
        check(library.rssn_value_new_float(ctypes.c_void_p(session.ptr), value,
                                           ctypes.byref(out)))
        return out.value
    if isinstance(value, complex):
        check(library.rssn_value_new_complex(ctypes.c_void_p(session.ptr), value.real,
                                             value.imag, ctypes.byref(out)))
        return out.value
    if isinstance(value, str):
        check(library.rssn_value_new_text(ctypes.c_void_p(session.ptr),
                                          value.encode("utf-8"), ctypes.byref(out)))
        return out.value
    # Anything else goes through `Value.of`, which knows how to build a record
    # from a dict and a tensor from a rectangular sequence. The cell it mints
    # is handed to the caller, so the temporary must not free it.
    return Value.of(session, value).detach()


class Callback:
    """A Python callable registered with a session as a host callback.

    Create one when the same function is passed to several calls, or when the
    engine should keep it: a callback registered implicitly by
    :meth:`Session.call_value` lives for that call only.

    Args:
        session: the session the callback belongs to.
        function: any callable.
        signature: the function's type, as `"(float) -> float"` or as the JSON
            object `{"params": [...], "ret": ...}`. Required unless `op` and
            `index` are given.
        op: take the signature from this operator's declared parameter instead.
        index: which argument of `op` the callback fills.
    """

    __slots__ = ("_session", "_value", "_function", "_thunk")

    def __init__(
        self,
        session: "Session",
        function: Callable[..., Any],
        *,
        signature: str | dict | None = None,
        op: str | None = None,
        index: int = 0,
    ) -> None:
        if not callable(function):
            raise TypeError("a callback must be callable")
        self._session = session
        self._function = function
        self._value: Value | None = None
        # ctypes needs the thunk object itself to stay alive for as long as the
        # engine may call it.
        self._thunk = _ffi.RssnCallbackFn(self._invoke)
        library = _ffi.library()
        cell = ctypes.c_uint32()
        if signature is not None:
            payload = signature if isinstance(signature, dict) else parse_signature(signature)
            check(
                library.rssn_callback_new(
                    ctypes.c_void_p(session.ptr),
                    json.dumps(payload).encode("utf-8"),
                    self._thunk,
                    None,
                    ctypes.byref(cell),
                )
            )
        elif op is not None:
            check(
                library.rssn_callback_new_for(
                    ctypes.c_void_p(session.ptr),
                    op.encode("utf-8"),
                    ctypes.c_size_t(index),
                    self._thunk,
                    None,
                    ctypes.byref(cell),
                )
            )
        else:
            raise ValueError("a Callback needs a signature, or an operator and an index")
        self._value = Value(session, cell.value)

    # -- lifetime ---------------------------------------------------------
    @property
    def session(self) -> "Session":
        return self._session

    @property
    def cell(self) -> int:
        """The cell id holding the `function` value."""
        if self._value is None:
            raise RuntimeError("this Callback has been closed")
        return self._value.cell

    @property
    def value(self) -> Value:
        """The callback as a value, to pass to :meth:`Session.call_value`."""
        if self._value is None:
            raise RuntimeError("this Callback has been closed")
        return self._value

    def close(self) -> None:
        """Retire the callback and release its cell."""
        if self._value is not None:
            _ffi.library().rssn_callback_free(ctypes.c_void_p(self._session.ptr),
                                              self._value.cell)
            self._value = None

    def __enter__(self) -> "Callback":
        return self

    def __exit__(self, *exc_info) -> None:
        self.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:  # interpreter shutdown, closed session
            pass

    def __repr__(self) -> str:
        state = "closed" if self._value is None else f"cell={self._value.cell}"
        return f"<Callback {self._function!r} {state}>"

    # -- the C side -------------------------------------------------------
    def _invoke(self, _user, session_ptr, args, count, out) -> int:
        try:
            values = [to_python(self._session, args[index]) for index in range(count)]
            cell = from_python(self._session, self._function(*values))
            out[0] = cell
            return 0
        except BaseException as error:  # noqa: BLE001 - the boundary must not unwind
            message = f"{type(error).__name__}: {error}"
            try:
                _ffi.library().rssn_callback_fail(message.encode("utf-8", "replace"))
            except Exception:
                pass
            return FAILED


def callback_signatures(session: "Session", op: str, index: int) -> list[dict]:
    """The signatures operator `op` declares for a callback at argument `index`.

    Usually one; an operator with an objective *and* a gradient parameter
    declares two different function types at different indices, and a few
    declare more than one at the same index. Each element is the JSON object
    :class:`Callback` accepts.
    """
    library = _ffi.library()
    raw = ctypes.c_void_p()
    check(
        library.rssn_callback_signatures(
            ctypes.c_void_p(session.ptr), op.encode("utf-8"), ctypes.c_size_t(index),
            ctypes.byref(raw),
        )
    )
    return json.loads(_ffi.take_string(raw.value) or "[]")
