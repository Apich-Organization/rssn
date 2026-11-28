"""The answer of a computation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    from .term import Term

__all__ = ["Answer"]


@dataclass(frozen=True)
class Answer:
    """A reduced term plus everything the engine reported about it.

    Attributes:
        term: the answer term (a closed form, or a float literal for a numeric
            target).
        text: the answer as infix text.
        value: the numeric value, when one is known.
        error: the estimated absolute error of `value`.
        reduced: whether the request was fully reduced; `False` means the term
            still contains unevaluated operators.
        iterations: outer iterations the engine performed.
        dimension: the inferred dimension as text, or `None`.
        diagnostic_codes: the stable machine-readable kind of each advisory
            finding (`"mismatch"`, `"unknown_unit"`, `"unannotated"`, ...),
            read from the thread-local report of the computation. The prose of
            each finding is in :meth:`Session.compute_json` under
            `"diagnostics"`; the C ABI only keeps the codes for a handle-based
            computation.
    """

    term: "Term"
    text: str
    value: float | None
    error: float | None
    reduced: bool
    iterations: int
    dimension: str | None = None
    diagnostic_codes: tuple[str, ...] = ()

    def __str__(self) -> str:
        return self.text

    def __repr__(self) -> str:
        parts = [f"text={self.text!r}"]
        if self.value is not None:
            parts.append(f"value={self.value!r}")
            if self.error is not None:
                parts.append(f"error={self.error!r}")
        if not self.reduced:
            parts.append("reduced=False")
        if self.dimension:
            parts.append(f"dimension={self.dimension!r}")
        if self.diagnostic_codes:
            parts.append(f"diagnostic_codes={list(self.diagnostic_codes)!r}")
        return f"<Answer {' '.join(parts)}>"

    @property
    def latex(self) -> str:
        """The answer term as LaTeX math."""
        return self.term.latex

    @property
    def typst(self) -> str:
        """The answer term as Typst math."""
        return self.term.typst

    @property
    def pretty(self) -> str:
        """The answer term as a Unicode drawing."""
        return self.term.pretty

    def tensor(self) -> Any:
        """The answer as a dense array, when it is a rectangular list."""
        return self.term.tensor()

    def float(self) -> float:
        """The numeric value, raising if the answer has none."""
        if self.value is None:
            from .errors import RssnError, Status

            raise RssnError(Status.NOT_NUMERIC, f"{self.text!r} has no numeric value")
        return self.value
