"""Sessions: the owner of every term, and the entry point of the API."""

from __future__ import annotations

import ctypes
import itertools
import json
import threading
from typing import Any, Iterable

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
        return Term(self, _ffi.library().rssn_term_int(self.ptr, int(value)))

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

    def _value_call(self, op: str, cells: list["Value"]) -> "Value":
        """One `rssn_value_call` with prepared cells."""
        from .value import Value

        buffer = _ffi.as_handle_buffer([value.cell for value in cells])
        out = ctypes.c_uint32()
        check(
            _ffi.library().rssn_value_call(
                ctypes.c_void_p(self.ptr), op.encode("utf-8"), buffer,
                ctypes.c_size_t(len(cells)), ctypes.byref(out),
            )
        )
        return Value(self, out.value)

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
        term = self.term(request)
        answer = _ffi.RssnAnswer()
        code = _ffi.library().rssn_compute(
            self.ptr, term.handle, config.ptr if config is not None else None, ctypes.byref(answer)
        )
        if int(code) != int(Status.OK):
            raise RssnError(Status(int(code)), _ffi.last_error())
        return self._answer(term, answer)

    def simplify(self, source: Any) -> Answer:
        """Parse and reduce `source` with the default symbolic configuration."""
        term = self.parse(source) if isinstance(source, str) else self.term(source)
        answer = _ffi.RssnAnswer()
        code = _ffi.library().rssn_compute(self.ptr, term.handle, None, ctypes.byref(answer))
        if int(code) != int(Status.OK):
            raise RssnError(Status(int(code)), _ffi.last_error())
        return self._answer(term, answer)

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
        terms = [self.term(request) for request in requests]
        count = len(terms)
        if count == 0:
            return []
        handles = _ffi.as_handle_buffer([term.handle for term in terms])
        answers = (_ffi.RssnAnswer * count)()
        statuses = (ctypes.c_int * count)()
        code = _ffi.library().rssn_compute_all(
            self.ptr,
            handles,
            ctypes.c_size_t(count),
            config.ptr if config is not None else None,
            answers,
            statuses,
        )
        if int(code) != int(Status.OK):
            raise RssnError(Status(int(code)), _ffi.last_error())
        out: list[Answer | RssnError | None] = []
        first_error: RssnError | None = None
        for index, term in enumerate(terms):
            status = int(statuses[index])
            if status == int(Status.OK):
                out.append(self._answer(term, answers[index]))
            else:
                error = RssnError(Status(status), _ffi.last_error())
                first_error = first_error or error
                out.append(error if errors == "collect" else None)
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
    def _answer(self, term: Term, raw: _ffi.RssnAnswer) -> Answer:
        """Build an :class:`Answer` from the C struct plus the dimension report."""
        # Read the thread's dimension codes first: `rssn_dimension` runs the
        # request again (a cache hit, but it re-records the findings).
        count = int(_ffi.library().rssn_dimension_count())
        codes: list[str] = []
        for index in range(count):
            code = _ffi.library().rssn_dimension_code(ctypes.c_size_t(index))
            codes.append(code.decode("utf-8", "replace") if code else "")
        answer_term = Term(self, raw.term)
        text = _ffi.take_string(_ffi.library().rssn_term_to_string(self.ptr, raw.term))
        dimension = self.dimension(answer_term)
        return Answer(
            term=answer_term,
            text=text,
            value=float(raw.value) if raw.has_value else None,
            error=float(raw.error) if raw.has_value else None,
            reduced=bool(raw.reduced),
            iterations=int(raw.iterations),
            dimension=dimension or None,
            diagnostic_codes=tuple(codes),
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


_manifest_cache: dict | None = None


def manifest() -> dict:
    """The standard kernel manifest as a dict (decoded once per process)."""
    global _manifest_cache
    if _manifest_cache is None:
        pointer = _ffi.library().rssn_manifest_json()
        _manifest_cache = json.loads(_ffi.take_string(pointer))
    return _manifest_cache
