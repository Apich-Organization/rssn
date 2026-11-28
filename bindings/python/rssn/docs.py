"""The operator reference: the engine's own documentation, searchable.

Every operator's documentation travels with the package: it is generated from
the kernel manifest into the function's docstring (`_ops.py`), so `help()` and
`inspect.getdoc()` already show it. This module adds the two things a REPL user
wants on top of that -- one operator by name, and a search over all of them.

    >>> import rssn
    >>> print(rssn.op_help("diff"))
    >>> rssn.op_search("brent")            # [(name, first line of the doc), ...]
    >>> rssn.op_search(category="calculus")[:5]
"""

from __future__ import annotations

import ctypes
import inspect
import json
import re
from typing import Iterable

from . import _ffi, _ops
from .errors import RssnError, Status, check

__all__ = ["op_info", "op_help", "op_search", "op_names", "op_categories", "op_rule_sets"]


def _arity(minimum: int, maximum: int | None) -> str:
    if maximum is None:
        return f"{minimum}+ arguments"
    if minimum == maximum:
        return f"{minimum} argument" if minimum == 1 else f"{minimum} arguments"
    return f"{minimum} to {maximum} arguments"


def op_info(name: str) -> dict:
    """Return one live operator contract, including overload keys and dimensions."""
    if not isinstance(name, str) or not name or "\0" in name:
        raise ValueError("provide one operator name")
    entry = getattr(_ffi.library(), "rssn_op_info_json", None)
    if entry is None:
        raise RssnError(Status.INVALID_ARGUMENT, "this RSSN library lacks operator contract lookup")
    result = ctypes.c_void_p()
    check(entry(name.encode("utf-8"), ctypes.byref(result)))
    return json.loads(_ffi.take_string(result.value))


def op_help(name: str) -> str:
    """The documentation of one operator, as text.

    The signature is the generated function's own, so the parameter names are
    the ones the engine's reference uses.

    Args:
        name: the operator's name as the manifest spells it (`"diff"`).

    Returns:
        The signature, the arity, the category and the documentation.

    Raises:
        KeyError: when the name is not a registered operator.
    """
    spec = _ops.OPS.get(name)
    if spec is None:
        raise KeyError(f"{name!r} is not a registered operator")
    binding = _ops.NAMES[name]
    function = getattr(_ops, binding)
    try:
        signature = str(inspect.signature(function))
    except (TypeError, ValueError):  # pragma: no cover - every operator has one
        signature = "(...)"
    minimum, maximum, reachable = spec
    lines = [
        f"rssn.{binding}{signature}",
        "",
        f"  operator   : {name}",
        f"  arity      : {_arity(minimum, maximum)}",
        f"  category   : {_category(name)}",
        f"  callable   : {'term ABI' if reachable else 'not reachable'}",
        "",
        inspect.getdoc(function) or "(no documentation)",
    ]
    return "\n".join(lines)


def _category(name: str) -> str:
    for category, members in _ops.BY_CATEGORY.items():
        if _ops.NAMES[name] in members:
            return category
    return "uncategorised"


def op_names() -> list[str]:
    """Every registered operator name, sorted."""
    return sorted(_ops.OPS)


def op_categories() -> list[str]:
    """The manifest's categories, sorted."""
    return sorted(_ops.BY_CATEGORY)


def op_rule_sets() -> list[str]:
    """The engine's rule sets, sorted."""
    return sorted(_ops.BY_RULE_SET)


def op_search(
    pattern: str | None = None,
    *,
    category: str | None = None,
    rule_set: str | None = None,
    limit: int = 50,
) -> list[tuple[str, str]]:
    """Operators whose name or documentation matches `pattern`.

    Args:
        pattern: a regular expression, matched case-insensitively against the
            name and the first line of the documentation; `None` matches all.
        category: keep one manifest category (`"calculus"`, `"numerical_matrix"`).
        rule_set: keep one rule set.
        limit: the maximum number of results.

    Returns:
        `(name, first line of the documentation)` pairs, in name order.
    """
    names: Iterable[str] = sorted(_ops.OPS)
    if category is not None:
        wanted = set(_ops.BY_CATEGORY.get(category, ()))
        names = [name for name in names if _ops.NAMES[name] in wanted]
    if rule_set is not None:
        wanted = set(_ops.BY_RULE_SET.get(rule_set, ()))
        names = [name for name in names if _ops.NAMES[name] in wanted]
    compiled = re.compile(pattern, re.IGNORECASE) if pattern is not None else None
    results: list[tuple[str, str]] = []
    for name in names:
        doc = inspect.getdoc(getattr(_ops, _ops.NAMES[name])) or ""
        summary = doc.splitlines()[0] if doc else ""
        if len(summary) > 120:
            summary = summary[:119].rstrip() + "…"
        if compiled is None or compiled.search(name) or compiled.search(doc):
            results.append((name, summary))
            if len(results) >= limit:
                break
    return results
