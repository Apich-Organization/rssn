"""ctypes binding of the rssn C ABI.

This module owns every `extern "C"` declaration the package uses. Nothing here
knows about ergonomics: handles are integers, strings are `bytes`, statuses are
integers. `session.py`, `term.py` and `config.py` build on it.

The library is located once, on import:

1. `RSSN_LIB` — an explicit path to `rssn.dll` / `librssn.so` / `librssn.dylib`;
2. `rssn/_lib/` — a copy shipped inside the package;
3. the usual build directories of a source checkout (`target/debug`,
   `target/ffi`, `target/release`, `target-bindings/debug`);
4. the platform loader search path (an installed `rssn` library).

Call :func:`load` to point the package at another library before creating a
session.
"""

from __future__ import annotations

import ctypes
import os
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


_LIB_NAMES = ("rssn.dll", "librssn.so", "librssn.dylib", "rssn.so", "rssn.dylib")
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))  # bindings/python/rssn -> repo

_lib: ctypes.CDLL | None = None
_lib_path: str | None = None


def _candidate_dirs() -> list[str]:
    dirs = [os.path.join(_HERE, "_lib")]
    for profile in ("debug", "ffi", "release"):
        dirs.append(os.path.join(_REPO, "target", profile))
        dirs.append(os.path.join(_REPO, "target-bindings", profile))
        dirs.append(os.path.join(_REPO, "target-ffi", profile))
        dirs.append(os.path.join(_REPO, "dist", "rssn", profile))
    return dirs


def _search() -> str | None:
    explicit = os.environ.get("RSSN_LIB")
    if explicit:
        if os.path.isdir(explicit):
            for name in _LIB_NAMES:
                candidate = os.path.join(explicit, name)
                if os.path.exists(candidate):
                    return candidate
            raise FileNotFoundError(f"RSSN_LIB={explicit!r} holds no rssn library")
        if not os.path.exists(explicit):
            raise FileNotFoundError(f"RSSN_LIB={explicit!r} does not exist")
        return explicit
    for directory in _candidate_dirs():
        for name in _LIB_NAMES:
            candidate = os.path.join(directory, name)
            if os.path.exists(candidate):
                return candidate
    return None


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
    proto("rssn_value_new_float", status, ptr, double, POINTER(c_uint32))
    proto("rssn_value_new_text", status, ptr, c_char_p, POINTER(c_uint32))
    proto("rssn_value_new_symbol", status, ptr, c_char_p, POINTER(c_uint32))
    proto("rssn_value_new_term", status, ptr, u32, POINTER(c_uint32))
    proto("rssn_value_new_list", status, ptr, POINTER(c_uint32), size, POINTER(c_uint32))
    proto("rssn_value_new_tensor", status, ptr, POINTER(double), size, POINTER(size), size,
          POINTER(c_uint32))
    proto("rssn_value_new_bool", status, ptr, c_int, POINTER(c_uint32))
    proto("rssn_value_new_complex", status, ptr, double, double, POINTER(c_uint32))
    proto("rssn_value_new_record", status, ptr, POINTER(c_char_p), POINTER(c_uint32), size,
          POINTER(c_uint32))
    proto("rssn_value_call", status, ptr, c_char_p, POINTER(c_uint32), size,
          POINTER(c_uint32))
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
    proto("rssn_sim_run", status, c_char_p, c_char_p, POINTER(ptr))

    # ownership
    proto("rssn_string_free", None, ptr)
    proto("rssn_bytes_free", None, ptr, size)
    proto("rssn_last_error", c_char_p)


def load(path: str | os.PathLike[str] | None = None) -> ctypes.CDLL:
    """Bind the package to a specific library, or to the one found on import."""
    global _lib, _lib_path
    if path is None:
        resolved = _search()
        if resolved is None:
            raise FileNotFoundError(
                "no rssn shared library found: set RSSN_LIB to rssn.dll / librssn.so, "
                "or build one with `bash scripts/build_ffi.sh`"
            )
    else:
        resolved = os.fspath(path)
    library = ctypes.CDLL(resolved)
    _bind(library)
    _lib, _lib_path = library, resolved
    return library


def library() -> ctypes.CDLL:
    """The bound library, loading it on first use."""
    if _lib is None:
        load()
    assert _lib is not None
    return _lib


def library_path() -> str:
    """The filesystem path of the bound library."""
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


def as_handle_buffer(handles: list[int]) -> "ctypes.Array[c_uint32]":
    return (c_uint32 * max(len(handles), 1))(*handles)


#: The cell id that means "no value".
RSSN_NO_VALUE = 0

__all__ += ["RSSN_NO_VALUE", "as_handle_buffer", "last_error", "take_bytes", "take_string"]
