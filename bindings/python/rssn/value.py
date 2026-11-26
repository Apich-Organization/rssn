"""Dynamic values: the cell ABI behind operators that take native objects.

The term API in :mod:`rssn.term` carries only values that are nodes of the
expression graph. An operator that returns a native object — an AD tape, a
compiled backend, a solver handle — returns a *value* instead, and an operator
that takes one can only be called with a value. Those are the operators the
manifest marks with `ffi = False` in :data:`rssn.OPS`; they are reached with
:meth:`rssn.Session.call_value`.

A :class:`Value` is a handle into a session-owned cell table, released when it
is collected or by :meth:`close`. Values are cheap to pass: the session keeps
them, and only the session that minted one can read it.

    >>> s = rssn.Session()
    >>> backend = s.call_value("num_backend_current")
    >>> s.call_value("num_backend_name", backend).text()
    'interpreter'

    >>> tape = s.call_value("num_ad_tape_new", 16)
    >>> x = s.call_value("num_ad_tape_var", tape, 1.5)
    >>> z = s.call_value("num_ad_var_mul", s.call_value("num_ad_var_sin", x), x)
    >>> s.call_value("num_ad_var_value", z).float()
    1.4962424799060816
"""

from __future__ import annotations

import ctypes
import json
from typing import TYPE_CHECKING, Any

from . import _ffi
from .errors import RssnError, Status, check

if TYPE_CHECKING:  # pragma: no cover
    from .session import Session
    from .term import Term

__all__ = ["KIND_NAMES", "Value"]

#: The kind codes `rssn_value_kind` reports, in declaration order.
KIND_NAMES = (
    "int",
    "real",
    "float",
    "complex",
    "bool",
    "string",
    "symbol",
    "unit",
    "quantity",
    "tensor",
    "sparse",
    "record",
    "variant",
    "term",
    "function",
    "dyn",
    "opaque",
)


def _shape_of(values) -> tuple[int, ...]:
    """The extents of a nested sequence, or `()` when it is not nested."""
    shape: list[int] = []
    node = values
    while isinstance(node, (list, tuple)):
        shape.append(len(node))
        if not node:
            break
        node = node[0]
    return tuple(shape)


def _flatten(values, shape: tuple[int, ...]) -> list[float] | None:
    """The row-major floats of a rectangular numeric sequence, else `None`."""
    flat: list[float] = []

    def walk(node, depth: int) -> bool:
        if depth == len(shape):
            if isinstance(node, (list, tuple)) or isinstance(node, bool):
                return False
            if not isinstance(node, (int, float)):
                return False
            flat.append(float(node))
            return True
        if not isinstance(node, (list, tuple)) or len(node) != shape[depth]:
            return False
        return all(walk(child, depth + 1) for child in node)

    return flat if walk(values, 0) else None


class Value:
    """A dynamic value owned by a session.

    Attributes:
        session: the session that owns the cell.
        cell: the raw `uint32_t` cell id.
    """

    __slots__ = ("_session", "_cell", "_owned")

    def __init__(self, session: "Session", cell: int, owned: bool = True) -> None:
        if not cell:
            raise RssnError(Status.NO_VALUE, _ffi.last_error() or "no value")
        self._session = session
        self._cell = int(cell)
        self._owned = owned

    @classmethod
    def borrowed(cls, session: "Session", cell: int) -> "Value":
        """A view of a cell the caller does not own.

        This is what a host callback receives for its arguments: rssn owns those
        cells and releases them when the callback returns, so closing the view
        must not free them.
        """
        return cls(session, cell, owned=False)

    # -- lifetime ---------------------------------------------------------
    @property
    def session(self) -> "Session":
        """The session that owns this value."""
        return self._session

    @property
    def cell(self) -> int:
        """The raw cell id."""
        return self._cell

    def close(self) -> None:
        """Release the cell. The session releases the rest with itself."""
        if self._cell and self._owned:
            _ffi.library().rssn_value_free(ctypes.c_void_p(self._session.ptr), self._cell)
        self._cell = 0

    def detach(self) -> int:
        """The cell id, leaving ownership to the caller.

        Used at the callback boundary: the value a Python callable returns is
        handed to the engine, which reads it and frees it itself.
        """
        cell = self._cell
        self._cell = 0
        return cell

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:  # interpreter shutdown, missing library, closed session
            pass

    def __enter__(self) -> "Value":
        return self

    def __exit__(self, *exc_info) -> None:
        self.close()

    def __repr__(self) -> str:
        try:
            return f"<Value {self.kind()} {self.text()!r}>"
        except Exception:
            return f"<Value cell={self._cell}>"

    # -- construction -----------------------------------------------------
    @classmethod
    def of(cls, session: "Session", value: Any) -> "Value":
        """`value` as a cell of `session`.

        A :class:`Value` must belong to `session`; a
        :class:`~rssn.term.Term` is copied as a term value; numbers, strings
        and sequences of them become the matching value kinds.
        """
        if isinstance(value, Value):
            if value._session is not session:
                raise ValueError("value belongs to another session")
            return value
        from .callback import Callback
        from .term import Term

        if isinstance(value, Callback):
            if value.session is not session:
                raise ValueError("callback belongs to another session")
            # The callback owns the cell; this view must not free it.
            return cls.borrowed(session, value.cell)

        library = _ffi.library()
        out = ctypes.c_uint32()
        if isinstance(value, Term):
            if value.session is not session:
                raise ValueError("term belongs to another session")
            check(library.rssn_value_new_term(ctypes.c_void_p(session.ptr), value.handle,
                                              ctypes.byref(out)))
            return cls(session, out.value)
        if isinstance(value, bool):
            check(library.rssn_value_new_bool(ctypes.c_void_p(session.ptr), int(value),
                                              ctypes.byref(out)))
            return cls(session, out.value)
        if isinstance(value, int):
            check(library.rssn_value_new_int(ctypes.c_void_p(session.ptr), value,
                                             ctypes.byref(out)))
            return cls(session, out.value)
        if isinstance(value, float):
            check(library.rssn_value_new_float(ctypes.c_void_p(session.ptr), value,
                                               ctypes.byref(out)))
            return cls(session, out.value)
        if isinstance(value, complex):
            check(library.rssn_value_new_complex(
                ctypes.c_void_p(session.ptr), value.real, value.imag, ctypes.byref(out)))
            return cls(session, out.value)
        if isinstance(value, str):
            check(library.rssn_value_new_text(ctypes.c_void_p(session.ptr),
                                              value.encode("utf-8"), ctypes.byref(out)))
            return cls(session, out.value)
        if isinstance(value, dict):
            return cls._record(session, value, out)
        # A numpy array is a nested list with a shape; anything else that
        # spells `.tolist()` is treated the same way.
        if hasattr(value, "tolist") and hasattr(value, "shape"):
            value = value.tolist()
        if isinstance(value, (list, tuple)):
            return cls._sequence(session, value, out)
        raise TypeError(f"cannot use {type(value).__name__} as a value")

    @classmethod
    def _record(cls, session: "Session", fields: dict, out) -> "Value":
        """A record value from a dict of field name -> value."""
        names = [str(name).encode("utf-8") for name in fields]
        # The field values own their cells: a list comprehension of bare cells
        # would let the temporary `Value`s be collected -- and their cells freed
        # -- before the record is built.
        members = [cls.of(session, field) for field in fields.values()]
        cells = [member.cell for member in members]
        name_buffer = (ctypes.c_char_p * max(len(names), 1))(*names)
        cell_buffer = _ffi.as_handle_buffer(cells)
        check(
            _ffi.library().rssn_value_new_record(
                ctypes.c_void_p(session.ptr), name_buffer, cell_buffer,
                ctypes.c_size_t(len(names)), ctypes.byref(out),
            )
        )
        return cls(session, out.value)

    @classmethod
    def _sequence(cls, session: "Session", values, out) -> "Value":
        """A dense tensor for a rectangular rank-2+ numeric sequence, else a list.

        A flat sequence keeps the list form: a list of numbers is already a
        rank-1 tensor to the engine, and building it element by element keeps
        each element's own kind (an integer vector stays `tensor<int>`).
        """
        shape = _shape_of(values)
        flat = _flatten(values, shape) if len(shape) >= 2 else None
        if flat is not None:
            data = (ctypes.c_double * max(len(flat), 1))(*flat)
            extents = (ctypes.c_size_t * max(len(shape), 1))(*shape)
            check(
                _ffi.library().rssn_value_new_tensor(
                    ctypes.c_void_p(session.ptr), data, ctypes.c_size_t(len(flat)),
                    extents, ctypes.c_size_t(len(shape)), ctypes.byref(out),
                )
            )
            return cls(session, out.value)
        elements = [cls.of(session, element) for element in values]
        buffer = _ffi.as_handle_buffer([element.cell for element in elements])
        check(
            _ffi.library().rssn_value_new_list(
                ctypes.c_void_p(session.ptr), buffer, ctypes.c_size_t(len(elements)),
                ctypes.byref(out),
            )
        )
        return cls(session, out.value)

    # -- inspection -------------------------------------------------------
    def kind(self) -> str:
        """The value's kind: `"float"`, `"tensor"`, `"opaque"`, ..."""
        code = ctypes.c_uint8()
        check(
            _ffi.library().rssn_value_kind(
                ctypes.c_void_p(self._session.ptr), self._cell, ctypes.byref(code)
            )
        )
        return KIND_NAMES[code.value] if code.value < len(KIND_NAMES) else f"kind{code.value}"

    def float(self) -> float:
        """The numeric value."""
        out = ctypes.c_double()
        check(
            _ffi.library().rssn_value_to_f64(
                ctypes.c_void_p(self._session.ptr), self._cell, ctypes.byref(out)
            )
        )
        return out.value

    def __float__(self) -> float:
        return self.float()

    def __int__(self) -> int:
        return int(self.float())

    def __bool__(self) -> bool:
        try:
            return bool(self.float())
        except RssnError:
            pass
        try:
            return len(self) > 0
        except RssnError:
            # A variant, a unit, an opaque handle: it exists, so it is true.
            return True

    def __array__(self, dtype=None, copy=None):  # numpy interop
        array = self.tensor()
        if hasattr(array, "astype"):
            return array.astype(dtype) if dtype is not None else array
        import numpy as np  # noqa: PLC0415

        return np.asarray(array, dtype=dtype)

    def text(self) -> str:
        """The text of a string or symbol, or any value's display form."""
        pointer = ctypes.c_void_p()
        check(
            _ffi.library().rssn_value_to_text(
                ctypes.c_void_p(self._session.ptr), self._cell, ctypes.byref(pointer)
            )
        )
        return _ffi.take_string(pointer.value)

    def json(self) -> Any:
        """The value in its `serde` JSON form, decoded."""
        pointer = ctypes.c_void_p()
        check(
            _ffi.library().rssn_value_json(
                ctypes.c_void_p(self._session.ptr), self._cell, ctypes.byref(pointer)
            )
        )
        return json.loads(_ffi.take_string(pointer.value))

    def to_term(self) -> "Term":
        """The value as a term, when it has a graph form."""
        from .term import Term

        out = ctypes.c_uint32()
        check(
            _ffi.library().rssn_value_to_term(
                ctypes.c_void_p(self._session.ptr), self._cell, ctypes.byref(out)
            )
        )
        return Term(self._session, out.value)

    def tensor(self) -> Any:
        """A dense `f64` extraction, as a numpy array or a nested list."""
        from .term import _nest

        tensor = _ffi.RssnTensor()
        check(
            _ffi.library().rssn_value_to_tensor(
                ctypes.c_void_p(self._session.ptr), self._cell, ctypes.byref(tensor)
            )
        )
        try:
            shape = [int(tensor.shape[index]) for index in range(tensor.rank)]
            data = [float(tensor.data[index]) for index in range(tensor.len)]
        finally:
            _ffi.library().rssn_tensor_free(ctypes.byref(tensor))
        try:
            import numpy as np  # noqa: PLC0415

            return np.array(data, dtype=float).reshape(shape) if shape else np.array(data[0])
        except ImportError:
            if not shape:
                return data[0] if data else 0.0
            return _nest(data, shape)

    def __len__(self) -> int:
        out = ctypes.c_size_t()
        check(
            _ffi.library().rssn_value_len(
                ctypes.c_void_p(self._session.ptr), self._cell, ctypes.byref(out)
            )
        )
        return int(out.value)

    def __getitem__(self, index: int) -> "Value":
        out = ctypes.c_uint32()
        check(
            _ffi.library().rssn_value_at(
                ctypes.c_void_p(self._session.ptr), self._cell, ctypes.c_size_t(index),
                ctypes.byref(out),
            )
        )
        return Value(self._session, out.value)

    def __iter__(self):
        for index in range(len(self)):
            yield self[index]
