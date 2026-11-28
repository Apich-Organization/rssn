"""Sessions: the owner of every term, and the entry point of the API."""

from __future__ import annotations

import ctypes
import itertools
import json
import operator
import threading
from collections.abc import Mapping
from numbers import Complex, Real
from typing import Any, Iterable, cast

from . import _ffi
from .answer import Answer
from .config import Config
from .errors import RssnError, Status, check
from .term import Term, coerce

__all__ = [
    "Session",
    "default_session",
    "set_default_session",
    "use_session",
]


class Session:
    """An owned `RssnSession`.

    A session owns every term made through it and is single-threaded: its
    store sits behind a `RefCell` in Rust, so two threads using one session at
    once is a borrow conflict. Use one session per thread (see
    :func:`default_session`), or `compute_all` for parallel requests inside
    one session.

    Sessions are released by :meth:`close`, by the context manager, or when
    they are collected::

        with rssn.Session() as s:
            x = s.sym("x")
            print(s.compute(rssn.diff(x**2, x)).text)
    """

    __slots__ = ("_ptr", "_numeric_config", "_lock")

    def __init__(self) -> None:
        pointer = _ffi.library().rssn_session_new()
        if not pointer:
            raise RssnError(Status.RULE, _ffi.last_error() or "rssn_session_new failed")
        self._ptr = pointer
        self._numeric_config: Config | None = None
        self._lock = threading.Lock()

    # -- lifetime ---------------------------------------------------------
    @property
    def ptr(self) -> int:
        """The raw `RssnSession *`."""
        if not self._ptr:
            raise RuntimeError("this Session has been closed")
        return self._ptr

    def close(self) -> None:
        """Release the session and every term made through it."""
        if self._ptr:
            if self._numeric_config is not None:
                self._numeric_config.close()
                self._numeric_config = None
            _ffi.library().rssn_session_free(ctypes.c_void_p(self._ptr))
            self._ptr = 0

    def __enter__(self) -> "Session":
        return self

    def __exit__(self, *exc_info) -> None:
        self.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:  # interpreter shutdown, missing library
            pass

    def __repr__(self) -> str:
        return f"<Session {'closed' if not self._ptr else hex(self._ptr)}>"

    # -- terms ------------------------------------------------------------
    def sym(self, name: str) -> Term:
        """The symbol `name`."""
        return Term(self, _ffi.library().rssn_term_sym(self.ptr, name.encode("utf-8")))

    def int(self, value: int) -> Term:
        """An exact integer."""
        value = int(value)
        if -(1 << 63) <= value < (1 << 63):
            return Term(self, _ffi.library().rssn_term_int(self.ptr, value))
        from .value import Value

        with Value.of(self, value) as exact:
            return exact.to_term()

    def float(self, value: float) -> Term:
        """A floating point number."""
        return Term(self, _ffi.library().rssn_term_float(self.ptr, float(value)))

    def rational(self, numerator: int, denominator: int = 1) -> Term:
        """The exact fraction `numerator / denominator`."""
        handle = _ffi.library().rssn_term_rational(self.ptr, int(numerator), int(denominator))
        if handle == _ffi.RSSN_NO_TERM:
            raise RssnError(Status.INVALID_ARGUMENT, _ffi.last_error() or "zero denominator")
        return Term(self, handle)

    def parse(self, source: str) -> Term:
        """Parse infix text such as `"diff(sin(x)^2, x)"`."""
        handle = _ffi.library().rssn_term_parse(self.ptr, source.encode("utf-8"))
        if handle == _ffi.RSSN_NO_TERM:
            raise RssnError(Status.PARSE, _ffi.last_error() or f"cannot parse {source!r}")
        return Term(self, handle)

    def term(self, value: Any) -> Term:
        """`value` as a term: a `Term`, a number, a `Fraction`, or infix text."""
        if isinstance(value, Term):
            if value.session is not self:
                raise ValueError("term belongs to another session")
            return value
        return Term(self, coerce(value, self))

    def apply(self, op: str, *args: Any) -> Term:
        """Apply the operator `op` to `args`, as an unevaluated request.

        The operator name is not validated against the manifest by the C ABI;
        an unknown name or a wrong arity is reported as
        :class:`~rssn.errors.RssnError` with `Status.UNKNOWN_OPERATOR`.
        """
        from ._ops import OPS  # generated table: name -> (min, max, ffi)

        spec = OPS.get(op)
        if spec is not None:
            minimum, maximum = spec[0], spec[1]
            if len(args) < minimum or (maximum is not None and len(args) > maximum):
                expected = f"{minimum}" if maximum == minimum else f"{minimum}..{maximum or 'n'}"
                raise RssnError(
                    Status.UNKNOWN_OPERATOR,
                    f"{op!r} takes {expected} arguments, got {len(args)}",
                )
        for arg in args:
            if callable(arg):
                raise RssnError(
                    Status.INVALID_ARGUMENT,
                    f"{op!r} takes a function argument, and a term cannot hold one: "
                    "use Session.call_value, which registers a host callback",
                )
        handles = [coerce(arg, self) for arg in args]
        buffer = _ffi.as_handle_buffer(handles)
        handle = _ffi.library().rssn_term_apply(
            self.ptr, op.encode("utf-8"), buffer, ctypes.c_size_t(len(handles))
        )
        if handle == _ffi.RSSN_NO_TERM:
            raise RssnError(Status.UNKNOWN_OPERATOR, _ffi.last_error() or f"operator {op!r} failed")
        return Term(self, handle)

    # `call` is the name the kernel ABI uses for the same operation.
    call = apply

    def call_value(self, op: str, *args: Any) -> "Value":
        """Call `op` through the value ABI and return the result as a :class:`~rssn.value.Value`.

        This is how the operators whose arguments have no term form are
        reached: native object handles, records, host callbacks. A *callable*
        argument becomes a host callback automatically, with the signature the
        operator declares for that argument::

            s.call_value("num_brent", 1.0, 2.0, 1e-12, 64, lambda x: x*x - 2)

        Arguments may be :class:`~rssn.value.Value` objects, terms, numbers,
        strings, lists of those, dicts (records), or callables.
        """
        from .callback import Callback, callback_signatures
        from .value import Value

        callables = [index for index, arg in enumerate(args) if callable(arg)]
        if not callables:
            return self._value_call(op, [Value.of(self, arg) for arg in args])

        # A few operators declare more than one function type at the same
        # position (a residual that answers a number or a vector, say). Try the
        # declared ones in order and keep the first that the engine accepts.
        options = [callback_signatures(self, op, index) or [None] for index in callables]
        first_error: Exception | None = None
        for combination in itertools.islice(itertools.product(*options), 8):
            callbacks: list[Callback] = []
            try:
                cells: list[Value] = []
                for position, arg in enumerate(args):
                    if position in callables:
                        declared = combination[callables.index(position)]
                        callback = (
                            Callback(self, arg, signature=declared)
                            if declared is not None
                            else Callback(self, arg, op=op, index=position)
                        )
                        callbacks.append(callback)
                        cells.append(callback.value)
                    else:
                        cells.append(Value.of(self, arg))
                try:
                    return self._value_call(op, cells)
                except RssnError as error:
                    first_error = first_error or error
                    continue
            finally:
                for callback in callbacks:
                    callback.close()
        assert first_error is not None
        raise first_error

    def call_key(self, key: Mapping[str, Any], *args: Any) -> "Value":
        """Call the exact overload in a key from :func:`rssn.op_info`.

        Use a registered :class:`~rssn.callback.Callback` value for a function
        parameter. Argument types and dimensions are checked by the kernel.
        """
        from .value import Value

        rule_set, name = key.get("rule_set"), key.get("name")
        if any(not isinstance(part, str) or not part or "\0" in part for part in (rule_set, name)):
            raise ValueError("operator key requires a rule set and name")
        digest = key.get("signature_hash")
        if not isinstance(digest, (list, tuple, bytes, bytearray)) or len(digest) != 16:
            raise ValueError("operator key requires 16 signature hash bytes")
        hash_buffer = (ctypes.c_uint8 * 16).from_buffer_copy(bytes(digest))
        entry = getattr(_ffi.library(), "rssn_value_call_key", None)
        if entry is None:
            raise RssnError(Status.INVALID_ARGUMENT, "this RSSN library lacks keyed calls")
        cells = [Value.of(self, arg) for arg in args]
        return self._call_prepared(entry, (rule_set.encode(), name.encode(), hash_buffer), cells)

    def _value_call(self, op: str, cells: list["Value"]) -> "Value":
        return self._call_prepared(_ffi.library().rssn_value_call, (op.encode("utf-8"),), cells)

    def _call_prepared(self, entry, identifiers: tuple, cells: list["Value"]) -> "Value":
        from .value import Value

        buffer = _ffi.as_handle_buffer([value.cell for value in cells])
        out = ctypes.c_uint32()
        check(entry(self.ptr, *identifiers, buffer, len(cells), ctypes.byref(out)))
        return Value(self, out.value)

    def available_backends(self) -> tuple[str, ...]:
        """Selectors the loaded RSSN library can use on this host."""
        with self.call_value("num_backend_available") as choices:
            names = []
            for index in range(len(choices)):
                with choices[index] as choice:
                    names.append(choice.text())
            return tuple(names)

    def current_backend(self) -> str:
        """Return the active backend's name without changing the selection."""
        try:
            with self.call_value("num_backend_current_name") as name:
                return name.text()
        except RssnError as error:
            if error.status != Status.UNKNOWN_OPERATOR:
                raise

        with self.call_value("num_backend_current") as backend:
            try:
                with self.call_value("num_backend_name", backend) as name:
                    return name.text()
            finally:
                self.release_backend(backend)

    def open_backend(self, name: str) -> "Value":
        """Open a named backend for explicit compilation in this session.

        Use :meth:`release_backend` to release the native handle when finished.
        """
        if "\0" in name:
            raise ValueError("the backend selector contains NUL")
        return self.call_value("num_backend_named", name)

    def release_backend(self, backend: "Value") -> None:
        """Release a native backend handle and invalidate its copies."""
        with self.call_value("num_backend_release", backend):
            pass

    def set_default_backend(self, name: str) -> str:
        """Select the process-wide backend and return its effective name.

        Pass `"auto"` to restore the library's automatic selection policy.
        """
        if "\0" in name:
            raise ValueError("the backend selector contains NUL")
        with self.call_value("num_backend_select_default", name) as selected:
            return selected.text()

    def encode_npy(self, shape: Iterable[int], values: Iterable[float]) -> bytes:
        """Encode real row-major samples as a NumPy float64 array in memory."""
        from .value import Value

        extents = tuple(operator.index(size) for size in shape)
        samples = list(values)
        # NumPy complex scalars can cast to float by dropping their imaginary part.
        if any(isinstance(sample, Complex) and not isinstance(sample, Real) for sample in samples):
            raise TypeError("NumPy encoding requires real samples")
        with Value.tensor_cells(self, extents, (len(extents),),
                                exemplar=None if extents else 0) as dimensions:
            with Value._dense(self, samples, (len(samples),), ctypes.c_uint32()) as data:
                with self.call_value("num_npy_encode", dimensions, data) as wire:
                    return bytes.fromhex(wire.text())

    def decode_npy(self, payload: bytes | bytearray | memoryview) -> "Value":
        """Read NumPy float32 or float64 bytes as a binary64 tensor."""
        from .value import Value

        with self.call_value("num_npy_decode", memoryview(payload).hex()) as record:
            with record["shape"] as dimensions, record["data"] as data:
                extents = []
                for index in range(len(dimensions)):
                    with dimensions[index] as extent:
                        extents.append(int(extent))
                _, samples = data._tensor_parts()
                return Value._dense(self, samples, tuple(extents), ctypes.c_uint32())

    def callback(self, function: Any, *, signature: str | dict | None = None,
                 op: str | None = None, index: int = 0):
        """`function` as a reusable :class:`~rssn.callback.Callback`.

        Pass `signature="(float) -> float"` (or the JSON form), or name the
        operator and argument position the callback is for and let the engine
        supply the type::

            f = s.callback(lambda x: x*x, op="num_brent", index=4)
        """
        from .callback import Callback

        return Callback(self, function, signature=signature, op=op, index=index)

    def call_value_json(self, op: str, *args: Any) -> Any:
        """Call `op` through the value ABI and return the result as decoded JSON.

        A callable argument is registered exactly as in :meth:`call_value`.
        """
        from .callback import Callback, callback_signatures
        from .value import Value

        def invoke(cells: list["Value"]) -> Any:
            buffer = _ffi.as_handle_buffer([value.cell for value in cells])
            pointer = ctypes.c_void_p()
            check(
                _ffi.library().rssn_value_call_json(
                    ctypes.c_void_p(self.ptr), op.encode("utf-8"), buffer,
                    ctypes.c_size_t(len(cells)), ctypes.byref(pointer),
                )
            )
            return json.loads(_ffi.take_string(pointer.value))

        callables = [index for index, arg in enumerate(args) if callable(arg)]
        if not callables:
            return invoke([Value.of(self, arg) for arg in args])
        options = [callback_signatures(self, op, index) or [None] for index in callables]
        first_error: Exception | None = None
        for combination in itertools.islice(itertools.product(*options), 8):
            callbacks: list[Callback] = []
            try:
                cells: list["Value"] = []
                for position, arg in enumerate(args):
                    if position in callables:
                        declared = combination[callables.index(position)]
                        callback = (
                            Callback(self, arg, signature=declared)
                            if declared is not None
                            else Callback(self, arg, op=op, index=position)
                        )
                        callbacks.append(callback)
                        cells.append(callback.value)
                    else:
                        cells.append(Value.of(self, arg))
                try:
                    return invoke(cells)
                except RssnError as error:
                    first_error = first_error or error
                    continue
            finally:
                for callback in callbacks:
                    callback.close()
        assert first_error is not None
        raise first_error

    # -- computation ------------------------------------------------------
    def compute(self, request: Any, config: Config | None = None) -> Answer:
        """Reduce `request` (a term or infix text) to the target of `config`."""
        return cast(Answer, self.compute_all((request,), config)[0])

    def simplify(self, source: Any) -> Answer:
        """Parse and reduce `source` with the default symbolic configuration."""
        term = self.parse(source) if isinstance(source, str) else self.term(source)
        return self.compute(term)

    def compute_all(
        self,
        requests: Iterable[Any],
        config: Config | None = None,
        *,
        errors: str = "raise",
    ) -> list[Answer | RssnError | None]:
        """Compute several requests in parallel inside the library.

        `errors` decides what happens to a request that failed:
        `"raise"` raises the first failure, `"collect"` puts the
        :class:`~rssn.errors.RssnError` in its place, and `"ignore"` puts
        `None` there. The answers equal one :meth:`compute` per request.
        """
        if errors not in ("raise", "collect", "ignore"):
            raise ValueError("errors must be 'raise', 'collect', or 'ignore'")
        terms = [self.term(request) for request in requests]
        count = len(terms)
        if count == 0:
            return []
        handles = _ffi.as_handle_buffer([term.handle for term in terms])
        answers = (_ffi.RssnAnswer * count)()
        statuses = (ctypes.c_int * count)()
        library = _ffi.library()
        with_metadata = getattr(library, "rssn_compute_all_with_metadata", None)
        arguments = (self.ptr, handles, count, config.ptr if config is not None else None,
                     answers, statuses)
        metadata = None
        out: list[Answer | RssnError | None] = []
        first_error: RssnError | None = None
        try:
            if with_metadata is None:
                check(library.rssn_compute_all(*arguments))
            else:
                metadata = (_ffi.RssnAnswerMetadata * count)()
                check(with_metadata(*arguments, metadata))
            batch_error = _ffi.last_error()
            for index, term in enumerate(terms):
                status = int(statuses[index])
                row = metadata[index] if metadata is not None else None
                if status == int(Status.OK):
                    if row is None:
                        dimension = self.dimension(term, config)
                        codes = self._dimension_codes()
                    else:
                        dimension = row.take("dimension")
                        codes = tuple(row.take("diagnostic_codes").splitlines())
                    out.append(self._answer(answers[index], dimension, codes))
                else:
                    message = batch_error if row is None else row.take("message")
                    error = RssnError(Status(status), message)
                    first_error = first_error or error
                    out.append(error if errors == "collect" else None)
        finally:
            if metadata is not None:
                for row in metadata:
                    row.close()
        if errors == "raise" and first_error is not None:
            raise first_error
        return out

    def compute_json(self, request: Any, config: Config | None = None) -> dict:
        """The compute record as a dict (`RssnComputeRecord` as JSON)."""
        term = self.term(request)
        pointer = ctypes.c_void_p()
        check(
            _ffi.library().rssn_compute_json(
                self.ptr, term.handle, config.ptr if config is not None else None,
                ctypes.byref(pointer),
            )
        )
        return json.loads(_ffi.take_string(pointer.value))

    def compute_bincode(self, request: Any, config: Config | None = None) -> bytes:
        """The compute record in the deterministic `bincode-next` wire format."""
        term = self.term(request)
        pointer = ctypes.c_void_p()
        length = ctypes.c_size_t()
        check(
            _ffi.library().rssn_compute_bincode(
                self.ptr, term.handle, config.ptr if config is not None else None,
                ctypes.byref(pointer), ctypes.byref(length),
            )
        )
        return _ffi.take_bytes(pointer.value, length.value)

    def compute_all_bincode(self, requests: Iterable[Any], config: Config | None = None) -> bytes:
        """The batch record (`RssnComputeBatchRecord`) in the `bincode-next` format."""
        terms = [self.term(request) for request in requests]
        handles = _ffi.as_handle_buffer([term.handle for term in terms])
        pointer = ctypes.c_void_p()
        length = ctypes.c_size_t()
        check(
            _ffi.library().rssn_compute_all_bincode(
                self.ptr, handles, ctypes.c_size_t(len(terms)),
                config.ptr if config is not None else None,
                ctypes.byref(pointer), ctypes.byref(length),
            )
        )
        return _ffi.take_bytes(pointer.value, length.value)

    def eval(self, request: Any, bindings: dict[str, float] | None = None) -> float:
        """Evaluate a term numerically without running any rule."""
        term = self.term(request)
        bindings = bindings or {}
        names = list(bindings)
        name_buffer = (ctypes.c_char_p * max(len(names), 1))(
            *(name.encode("utf-8") for name in names)
        )
        value_buffer = (ctypes.c_double * max(len(names), 1))(*(bindings[name] for name in names))
        out = ctypes.c_double()
        check(
            _ffi.library().rssn_term_eval(
                self.ptr, term.handle, name_buffer, value_buffer,
                ctypes.c_size_t(len(names)), ctypes.byref(out),
            )
        )
        return out.value

    def dimension(self, request: Any, config: Config | None = None) -> str:
        """The inferred dimension of `request` as text (`""` when none)."""
        term = self.term(request)
        pointer = ctypes.c_void_p()
        check(
            _ffi.library().rssn_dimension(
                self.ptr, term.handle, config.ptr if config is not None else None,
                ctypes.byref(pointer),
            )
        )
        return _ffi.take_string(pointer.value)

    def numeric_config(self, tolerance: float = 1e-12) -> Config:
        """A cached configuration asking for a numeric answer."""
        if self._numeric_config is None:
            self._numeric_config = Config().numeric(tolerance)
        return self._numeric_config

    # -- kernel metadata --------------------------------------------------
    def manifest(self) -> dict:
        """The process-wide standard kernel manifest, decoded from JSON."""
        return manifest()

    @staticmethod
    def manifest_hash() -> int:
        """The first 64 bits of the manifest hash."""
        return int(_ffi.library().rssn_manifest_hash())

    def value_count(self) -> int:
        """The number of live value cells; a diagnostic for leak hunting."""
        return int(_ffi.library().rssn_value_count(ctypes.c_void_p(self.ptr)))

    def callback_count(self) -> int:
        """The number of live host callbacks; a diagnostic for leak hunting."""
        return int(_ffi.library().rssn_callback_count(ctypes.c_void_p(self.ptr)))

    def last_error(self) -> str:
        """The message of the last failed call on this thread."""
        return _ffi.last_error()

    # -- internals --------------------------------------------------------
    @staticmethod
    def _dimension_codes() -> tuple[str, ...]:
        count = int(_ffi.library().rssn_dimension_count())
        codes: list[str] = []
        for index in range(count):
            code = _ffi.library().rssn_dimension_code(ctypes.c_size_t(index))
            codes.append(code.decode("utf-8", "replace") if code else "")
        return tuple(codes)

    def _answer(self, raw: _ffi.RssnAnswer, dimension: str, codes: tuple[str, ...]) -> Answer:
        answer_term = Term(self, raw.term)
        text = _ffi.take_string(_ffi.library().rssn_term_to_string(self.ptr, raw.term))
        return Answer(
            term=answer_term,
            text=text,
            value=float(raw.value) if raw.has_value else None,
            error=float(raw.error) if raw.has_value else None,
            reduced=bool(raw.reduced),
            iterations=int(raw.iterations),
            dimension=dimension or None,
            diagnostic_codes=codes,
        )


# --------------------------------------------------------------------------
# a default session per thread, for the module-level operator functions
# --------------------------------------------------------------------------

_local = threading.local()


def default_session() -> Session:
    """The calling thread's session, created on first use."""
    session = getattr(_local, "session", None)
    if session is None or not session._ptr:
        session = Session()
        _local.session = session
    return session


def set_default_session(session: Session | None) -> Session | None:
    """Install `session` as this thread's default and return the previous one."""
    previous = getattr(_local, "session", None)
    _local.session = session
    return previous


class use_session:  # noqa: N801 - used as a context manager
    """Run a block with `session` as the default::

        with rssn.use_session(s):
            rssn.diff(x**2, x)
    """

    __slots__ = ("_session", "_previous")

    def __init__(self, session: Session) -> None:
        self._session = session
        self._previous: Session | None = None

    def __enter__(self) -> Session:
        self._previous = set_default_session(self._session)
        return self._session

    def __exit__(self, *exc_info) -> None:
        set_default_session(self._previous)


_manifest_cache: tuple[ctypes.CDLL, dict] | None = None


def manifest() -> dict:
    """The standard kernel manifest, reused while the library stays bound."""
    global _manifest_cache
    library = _ffi.library()
    if _manifest_cache is None or _manifest_cache[0] is not library:
        pointer = library.rssn_manifest_json()
        _manifest_cache = (library, json.loads(_ffi.take_string(pointer)))
    return _manifest_cache[1]
