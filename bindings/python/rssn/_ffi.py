"""ctypes binding of the rssn C ABI.

This module owns every `extern "C"` declaration the package uses. Nothing here
knows about ergonomics: handles are integers, strings are `bytes`, statuses are
integers. `session.py`, `term.py` and `config.py` build on it.

The library is located on first use:

1. `RSSN_LIB` — an explicit path to `rssn.dll` / `librssn.so` / `librssn.dylib`;
2. `rssn/_lib/` — a copy shipped inside the package;
3. the extracted distribution's root when importing its in-place SDK;
4. the Bazel FFI output or other build directories of a source checkout
   (`target/debug`, `target/ffi`, `target/release`, `target-bindings/debug`);
5. the platform loader search path (an installed `rssn` library).

Call :func:`load` to point the package at another library before creating a
session.
"""

from __future__ import annotations

import ctypes
import os
import sys
from ctypes import (
    POINTER,
    CFUNCTYPE,
    Structure,
    c_bool,
    c_char_p,
    c_double,
    c_int,
    c_int64,
    c_size_t,
    c_uint32,
    c_uint8,
    c_void_p,
)

__all__ = [
    "RSSN_NO_TERM",
    "RssnAnswer",
    "RssnAnswerMetadata",
    "RssnCallbackFn",
    "RssnTensor",
    "library",
    "library_path",
    "load",
    "raw",
]

#: The handle value that means "no term" (a failed constructor).
RSSN_NO_TERM = 0xFFFFFFFF


class RssnAnswer(Structure):
    """`struct RssnAnswer` from `rssn.h`."""

    _fields_ = [
        ("term", c_uint32),
        ("has_value", c_bool),
        ("value", c_double),
        ("error", c_double),
        ("reduced", c_bool),
        ("iterations", c_size_t),
    ]


class RssnAnswerMetadata(Structure):
    """Owned metadata strings for a batch answer."""

    _fields_ = [
        ("dimension", c_void_p),
        ("diagnostic_codes", c_void_p),
        ("message", c_void_p),
    ]

    def take(self, field: str) -> str:
        pointer = getattr(self, field)
        setattr(self, field, None)
        return take_string(pointer)

    def close(self) -> None:
        for field, _ in self._fields_:
            pointer = getattr(self, field)
            if pointer:
                setattr(self, field, None)
                library().rssn_string_free(pointer)


class RssnTensor(Structure):
    """`struct RssnTensor` from `rssn.h`; row-major `f64` data."""

    _fields_ = [
        ("rank", c_size_t),
        ("shape", POINTER(c_size_t)),
        ("len", c_size_t),
        ("data", POINTER(c_double)),
    ]


#: The C signature of a host callback (`RssnCallbackFn` in `rssn.h`).
#:
#: `(user, session, args, count, out) -> RssnStatus`. `args` holds `count` live
#: cell ids owned by rssn; `out` must be set to a cell id holding the result.
RssnCallbackFn = CFUNCTYPE(
    c_int, c_void_p, c_void_p, POINTER(c_uint32), c_size_t, POINTER(c_uint32)
)


if os.name == "nt":
    _LIB_NAMES = ("rssn.dll", "librssn.dll")
elif sys.platform == "darwin":
    _LIB_NAMES = ("librssn.dylib", "rssn.dylib")
else:
    _LIB_NAMES = ("librssn.so", "librssn.so.0", "rssn.so")
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, "..", "..", ".."))

_lib: ctypes.CDLL | None = None
_lib_path: str | None = None


def _candidate_dirs() -> list[str]:
    dirs = [os.path.join(_HERE, "_lib")]
    sdk = os.path.join(_REPO, "bindings", "python", "rssn")
    if os.path.normcase(_HERE) == os.path.normcase(sdk) and os.path.isfile(os.path.join(_REPO, "rssn.h")):
        dirs.append(_REPO)
        if os.path.isfile(os.path.join(_REPO, "MODULE.bazel")):
            dirs.extend(os.path.join(_REPO, "bazel-bin", target, "ffi")
                        for target in ("ffi_distribution", "ffi", "ffi_distribution_portable",
                                       "ffi_dev", "ffi_distribution_dev"))
    for profile in ("debug", "ffi", "release"):
        dirs.append(os.path.join(_REPO, "target", profile))
        dirs.append(os.path.join(_REPO, "target-bindings", profile))
        dirs.append(os.path.join(_REPO, "target-ffi", profile))
        dirs.append(os.path.join(_REPO, "dist", "rssn", profile))
    return dirs


def _library_in(directory: str) -> str | None:
    for name in _LIB_NAMES:
        candidate = os.path.join(directory, name)
        if os.path.isfile(candidate):
            return os.path.abspath(candidate)
    return None


def _resolve_path(path: str) -> str:
    if os.path.isdir(path):
        candidate = _library_in(path)
        if candidate is None:
            raise FileNotFoundError(f"{path!r} holds no rssn shared library")
        return candidate
    return os.path.abspath(path) if os.path.isfile(path) else path


def _bind(lib: ctypes.CDLL) -> None:
    """Declare the prototype of every entry point in `rssn.h`."""
    u32, size, double = c_uint32, c_size_t, c_double
    status = ctypes.c_int
    ptr = c_void_p

    def proto(name: str, restype, *argtypes) -> None:
        function = getattr(lib, name)
        function.restype = restype
        function.argtypes = list(argtypes)

    # session and configuration
    proto("rssn_version", c_char_p)
    proto("rssn_session_new", ptr)
    proto("rssn_session_free", None, ptr)
    proto("rssn_config_new", ptr)
    proto("rssn_config_free", None, ptr)
    proto("rssn_config_symbolic", status, ptr)
    proto("rssn_config_numeric", status, ptr, double)
    proto("rssn_config_bind", status, ptr, c_char_p, double)
    proto("rssn_config_assume", status, ptr, c_char_p, c_uint8)
    proto("rssn_config_budget", status, ptr, size, size)
    proto("rssn_config_dimension_mode", status, ptr, c_uint8)
    proto("rssn_config_bypass_dimensions", status, ptr)
    proto("rssn_config_dimension", status, ptr, c_char_p, c_char_p)
    proto("rssn_config_define_base_dimension", status, ptr, c_char_p)

    # terms
    proto("rssn_term_parse", u32, ptr, c_char_p)
    proto("rssn_term_sym", u32, ptr, c_char_p)
    proto("rssn_term_int", u32, ptr, c_int64)
    proto("rssn_term_float", u32, ptr, double)
    proto("rssn_term_rational", u32, ptr, c_int64, c_int64)
    proto("rssn_term_apply", u32, ptr, c_char_p, POINTER(c_uint32), size)
    proto("rssn_term_eval", status, ptr, u32, POINTER(c_char_p), POINTER(double), size,
          POINTER(double))
    proto("rssn_term_to_f64", status, ptr, u32, POINTER(double))
    proto("rssn_term_tensor_data", status, ptr, u32, POINTER(RssnTensor))
    proto("rssn_tensor_free", None, POINTER(RssnTensor))
    proto("rssn_term_to_string", ptr, ptr, u32)
    proto("rssn_term_to_latex", ptr, ptr, u32)
    proto("rssn_term_to_typst", ptr, ptr, u32)
    proto("rssn_term_to_pretty", ptr, ptr, u32)

    # compute
    proto("rssn_compute", status, ptr, u32, ptr, POINTER(RssnAnswer))
    proto("rssn_compute_all", status, ptr, POINTER(c_uint32), size, ptr,
          POINTER(RssnAnswer), POINTER(ctypes.c_int))
    if hasattr(lib, "rssn_compute_all_with_metadata"):
        proto("rssn_compute_all_with_metadata", status, ptr, POINTER(c_uint32), size, ptr,
              POINTER(RssnAnswer), POINTER(ctypes.c_int), POINTER(RssnAnswerMetadata))
    proto("rssn_compute_json", status, ptr, u32, ptr, POINTER(ptr))
    proto("rssn_compute_bincode", status, ptr, u32, ptr, POINTER(ptr), POINTER(size))
    proto("rssn_compute_all_bincode", status, ptr, POINTER(c_uint32), size, ptr,
          POINTER(ptr), POINTER(size))
    proto("rssn_simplify", u32, ptr, c_char_p)
    proto("rssn_dimension", status, ptr, u32, ptr, POINTER(ptr))
    proto("rssn_dimension_count", size)
    proto("rssn_dimension_code", c_char_p, size)

    # dynamic values (the value ABI: operators with native-object arguments)
    proto("rssn_value_new_int", status, ptr, c_int64, POINTER(c_uint32))
    if hasattr(lib, "rssn_value_new_int_decimal"):
        proto("rssn_value_new_int_decimal", status, ptr, c_char_p, POINTER(c_uint32))
    proto("rssn_value_new_float", status, ptr, double, POINTER(c_uint32))
    proto("rssn_value_new_text", status, ptr, c_char_p, POINTER(c_uint32))
    proto("rssn_value_new_symbol", status, ptr, c_char_p, POINTER(c_uint32))
    proto("rssn_value_new_term", status, ptr, u32, POINTER(c_uint32))
    proto("rssn_value_new_list", status, ptr, POINTER(c_uint32), size, POINTER(c_uint32))
    proto("rssn_value_new_tensor", status, ptr, POINTER(double), size, POINTER(size), size,
          POINTER(c_uint32))
    if hasattr(lib, "rssn_value_new_tensor_cells"):
        proto("rssn_value_new_tensor_cells", status, ptr, POINTER(c_uint32), size,
              POINTER(size), size, u32, POINTER(c_uint32))
    if hasattr(lib, "rssn_value_with_coherent_dimension"):
        proto("rssn_value_with_coherent_dimension", status, ptr, u32, POINTER(c_uint8), size,
              POINTER(c_uint32))
    proto("rssn_value_new_bool", status, ptr, c_int, POINTER(c_uint32))
    proto("rssn_value_new_complex", status, ptr, double, double, POINTER(c_uint32))
    proto("rssn_value_new_record", status, ptr, POINTER(c_char_p), POINTER(c_uint32), size,
          POINTER(c_uint32))
    proto("rssn_value_call", status, ptr, c_char_p, POINTER(c_uint32), size,
          POINTER(c_uint32))
    if hasattr(lib, "rssn_value_call_key"):
        proto("rssn_value_call_key", status, ptr, c_char_p, c_char_p,
              POINTER(c_uint8), POINTER(c_uint32), size, POINTER(c_uint32))
    proto("rssn_value_call_json", status, ptr, c_char_p, POINTER(c_uint32), size,
          POINTER(ptr))
    proto("rssn_value_kind", status, ptr, u32, POINTER(c_uint8))
    proto("rssn_value_to_f64", status, ptr, u32, POINTER(double))
    proto("rssn_value_to_text", status, ptr, u32, POINTER(ptr))
    proto("rssn_value_to_term", status, ptr, u32, POINTER(c_uint32))
    proto("rssn_value_to_tensor", status, ptr, u32, POINTER(RssnTensor))
    proto("rssn_value_len", status, ptr, u32, POINTER(size))
    proto("rssn_value_at", status, ptr, u32, size, POINTER(c_uint32))
    proto("rssn_value_field_name", status, ptr, u32, size, POINTER(ptr))
    proto("rssn_value_tag", status, ptr, u32, POINTER(ptr))
    proto("rssn_value_json", status, ptr, u32, POINTER(ptr))
    if hasattr(lib, "rssn_value_clone"):
        proto("rssn_value_clone", status, ptr, u32, POINTER(c_uint32))
    proto("rssn_value_free", None, ptr, u32)
    proto("rssn_value_count", size, ptr)

    # host callbacks (the value ABI's function arguments)
    proto("rssn_callback_new", status, ptr, c_char_p, RssnCallbackFn, ptr, POINTER(c_uint32))
    proto("rssn_callback_new_for", status, ptr, c_char_p, size, RssnCallbackFn, ptr,
          POINTER(c_uint32))
    proto("rssn_callback_signatures", status, ptr, c_char_p, size, POINTER(ptr))
    proto("rssn_callback_free", None, ptr, c_uint32)
    proto("rssn_callback_count", size, ptr)
    proto("rssn_callback_fail", None, c_char_p)

    # embedding kernel and simulations
    proto("rssn_manifest_json", ptr)
    proto("rssn_manifest_hash", ctypes.c_uint64)
    if hasattr(lib, "rssn_op_info_json"):
        proto("rssn_op_info_json", status, c_char_p, POINTER(ptr))
    proto("rssn_sim_run", status, c_char_p, c_char_p, POINTER(ptr))

    # ownership
    proto("rssn_string_free", None, ptr)
    proto("rssn_bytes_free", None, ptr, size)
    proto("rssn_last_error", c_char_p)


def load(path: str | os.PathLike[str] | None = None) -> ctypes.CDLL:
    """Bind to a shared-library file or directory before creating sessions."""
    global _lib, _lib_path
    requested = os.fspath(path) if path is not None else os.environ.get("RSSN_LIB")
    if requested:
        if path is None and os.path.dirname(requested) and not os.path.exists(requested):
            raise FileNotFoundError(f"RSSN_LIB={requested!r} does not exist")
        candidates = (_resolve_path(requested),)
    else:
        local = tuple(
            candidate for directory in _candidate_dirs()
            if (candidate := _library_in(directory)) is not None
        )
        candidates = local + _LIB_NAMES

    last_error = None
    for candidate in candidates:
        try:
            library = ctypes.CDLL(candidate)
        except OSError as error:
            last_error = error
            continue
        try:
            _bind(library)
        except AttributeError as error:
            raise OSError(f"{candidate!r} is missing an RSSN ABI entry point: {error}") from error
        _lib, _lib_path = library, candidate
        return library

    if requested:
        raise OSError(f"cannot load the RSSN shared library at {requested!r}: {last_error}") from last_error
    raise FileNotFoundError(
        "no rssn shared library found: set RSSN_LIB to the library file or its directory, "
        "or build one with `bazel build //:ffi_dev`"
    ) from last_error


def library() -> ctypes.CDLL:
    """The bound library, loading it on first use."""
    if _lib is None:
        load()
    assert _lib is not None
    return _lib


def library_path() -> str:
    """The path or platform-loader name used for the bound library."""
    if _lib is None:
        load()
    assert _lib_path is not None
    return _lib_path


def raw() -> ctypes.CDLL:
    """The raw `ctypes` library, for calls this package does not wrap."""
    return library()


# Strings rssn hands out are owned by the caller and must be released with
# `rssn_string_free`; the pointer is kept so the right address is freed.
def take_string(pointer: int) -> str:
    """Read a caller-owned `char *` and free it."""
    if not pointer:
        return ""
    try:
        return ctypes.string_at(pointer).decode("utf-8", "replace")
    finally:
        library().rssn_string_free(c_void_p(pointer))


def take_bytes(pointer: int, length: int) -> bytes:
    """Read a caller-owned `uint8_t *` buffer and free it."""
    if not pointer:
        return b""
    try:
        return ctypes.string_at(pointer, length)
    finally:
        library().rssn_bytes_free(c_void_p(pointer), c_size_t(length))


def last_error() -> str:
    """The message of the last failed call on this thread (may be empty)."""
    message = library().rssn_last_error()
    return message.decode("utf-8", "replace") if message else ""


#: `c_uint32 * n` builds a new ctypes class on every call, which costs more than
#: the call itself for the one- and two-argument operators that dominate. The
#: classes are immutable, so they are cached by length.
_HANDLE_TYPES: dict[int, type] = {}


def as_handle_buffer(handles: list[int]) -> "ctypes.Array[c_uint32]":
    """The handles as the `uint32_t *` the ABI reads."""
    length = len(handles) or 1
    factory = _HANDLE_TYPES.get(length)
    if factory is None:
        factory = _HANDLE_TYPES[length] = c_uint32 * length
    return factory(*handles)


#: The cell id that means "no value".
RSSN_NO_VALUE = 0

__all__ += ["RSSN_NO_VALUE", "as_handle_buffer", "last_error", "take_bytes", "take_string"]
