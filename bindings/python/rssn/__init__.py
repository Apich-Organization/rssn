"""Python bindings for rssn, a scientific computing library.

rssn is a computer algebra system, a set of numerical kernels and a physics
simulation engine behind one small API. This package binds its C ABI
(`rssn.h`), so it needs a shared library: set `RSSN_LIB` to `rssn.dll` /
`librssn.so` / `librssn.dylib`, or build one in a source checkout with
`bazel build //:ffi_dev`.

    >>> import rssn
    >>> s = rssn.Session()
    >>> x = s.sym("x")
    >>> s.compute(rssn.diff(rssn.sin(x) * rssn.cos(x), x)).text
    'cos(2*x)'
    >>> s.compute(rssn.diff(rssn.sin(x) * rssn.cos(x), x),
    ...           rssn.Config().numeric(1e-12).bind("x", 0.0)).value
    1.0

Operators are *requests*, not results: `rssn.diff(f, x)` returns an unevaluated
term and `Session.compute` reduces it. Every registered operator has a function
here, named as the manifest names it (`_ops.OPS` is the table); names that are
Python keywords get a trailing underscore (`and_`, `not_`, `or_`).

Because the operator names include Python builtins (`abs`, `sum`, `pow`,
`round`, `list`), `import rssn` is the recommended form: `from rssn import *`
would shadow those builtins in the importing module.
"""

from __future__ import annotations

from . import _ffi, constants, scenarios, sim, units
from ._ffi import load, library_path
from ._ops import *  # noqa: F403 - the operator facade
from ._ops import BY_CATEGORY, BY_RULE_SET, NAMES, OPS
from .answer import Answer
from .callback import Callback, callback_signatures, parse_signature
from .config import Assumption, Config, DimensionMode
from .docs import op_categories, op_help, op_info, op_names, op_rule_sets, op_search
from .errors import RssnError, Status
from .session import (
    Session,
    default_session,
    manifest,
    set_default_session,
    use_session,
)
from .term import Term
from .value import KIND_NAMES, Value

__all__ = sorted(
    {
        "Answer",
        "Assumption",
        "BY_CATEGORY",
        "BY_RULE_SET",
        "Callback",
        "Config",
        "DimensionMode",
        "NAMES",
        "OPS",
        "RssnError",
        "Session",
        "Status",
        "Term",
        "Value",
        "KIND_NAMES",
        "callback_signatures",
        "constants",
        "default_session",
        "library_path",
        "load",
        "manifest",
        "manifest_hash",
        "op_categories",
        "op_help",
        "op_info",
        "op_names",
        "op_rule_sets",
        "op_search",
        "parse_signature",
        "scenarios",
        "set_default_session",
        "sim",
        "units",
        "use_session",
        "version",
    }
    | {NAMES[name] for name in OPS}
)


def version() -> str:
    """The rssn version string, e.g. `"3.0.0-beta.1"`."""
    return _ffi.library().rssn_version().decode("utf-8")


def manifest_hash() -> int:
    """The first 64 bits of the kernel manifest hash."""
    return int(_ffi.library().rssn_manifest_hash())


def __getattr__(name: str):
    if name == "__version__":
        return version()
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
