"""Higher-order operators: Python callables the engine calls.

    RSSN_LIB=/path/to/rssn.dll python bindings/examples/python/callbacks.py

Some operators take a function and call it — root finders, quadrature,
optimisers, bootstrap. In Python any callable works: the binding turns it into
the C function pointer the engine knows how to call, and converts the arguments
and the result for you.

Run from the repository root, or with `RSSN_LIB` pointing at a built library.
"""

from __future__ import annotations

import rssn


def root_of_two() -> None:
    """`sqrt(2)` as the root of `x**2 - 2`, by Brent's method."""
    session = rssn.Session()
    answer = session.call_value(
        "num_brent",
        1.0, 2.0,          # a, b: a bracketing interval
        1e-12,             # tolerance
        64,                # iteration budget
        lambda x: x * x - 2.0,
    )
    root = float(answer[0])
    print(f"num_brent:        sqrt(2) = {root!r}")
    print(f"                  evaluations = {int(answer[1])}")


def through_the_facade() -> None:
    """The generated facade routes a higher-order operator the same way."""
    session = rssn.Session()
    answer = rssn.num_brent(1.0, 2.0, 1e-12, 64, lambda x: x * x - 2.0, session=session)
    print(f"rssn.num_brent:   {float(answer[0])!r}")


def minimise() -> None:
    """A two-dimensional objective and its gradient, as plain Python functions.

    `num_bfgs(x0, tolerance, max_iterations, max_evaluations, objective, gradient)`.
    The objective receives the point as a numpy array (or a list without numpy)
    and answers a float; the gradient answers a list of the same length.
    """
    session = rssn.Session()

    def objective(xs):
        return float(sum((x - 1.0) ** 2 for x in xs))

    def gradient(xs):
        return [2.0 * (x - 1.0) for x in xs]

    answer = session.call_value(
        "num_bfgs", [3.0, -2.0], 1e-10, 100, 1000, objective, gradient
    )
    point = [float(value) for value in answer[0]]
    print(f"num_bfgs:         minimum at {point}")
    print(f"                  value = {float(answer[1])!r}, "
          f"converged = {answer[4].json()}")


def explicit_callback() -> None:
    """A callback that outlives one call, with its signature spelled out.

    `callback_signatures` says what the engine expects at an argument position,
    and `Callback` registers one explicitly — useful when the same function is
    passed to several calls, or when the type is not the declared one.
    """
    session = rssn.Session()
    print(f"declared at num_brent/4: {rssn.callback_signatures(session, 'num_brent', 4)}")

    with session.callback(lambda x: x * x - 2.0, signature="(float) -> float") as residual:
        answer = session.call_value("num_brent", 1.0, 2.0, 1e-12, 64, residual)
        print(f"explicit Callback: {float(answer[0])!r}")


def failures_cross() -> None:
    """An exception in the callback aborts the operator with its message."""
    session = rssn.Session()
    try:
        session.call_value(
            "num_brent", 1.0, 2.0, 1e-12, 64,
            lambda x: (_ for _ in ()).throw(ValueError("no thanks")),
        )
    except rssn.RssnError as error:
        print(f"callback failure: {error.status.name}: {error}")


if __name__ == "__main__":
    root_of_two()
    through_the_facade()
    minimise()
    explicit_callback()
    failures_cross()
