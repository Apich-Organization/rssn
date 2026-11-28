"""Write `cases.tsv`: the cross-language conformance table.

Run once to (re)create the file:

    python bindings/tests/make_cases.py

The expectation columns are left empty; `run_conformance.py --record` fills them
from the engine, and the Python, Julia and C++ runners then all verify the same
rows.
"""

from __future__ import annotations

import pathlib

HEADER = [
    "name", "kind", "payload", "target", "tolerance", "bindings",
    "expect_text", "expect_value", "expect_status",
]

# name, kind, payload, target, tolerance, bindings
CASES = [
    ("calculus.derivative", "compute", "diff(sin(x)*cos(x), x)", "symbolic", "", ""),
    ("calculus.derivative_numeric", "compute", "diff(sin(x)*cos(x), x)", "numeric", "1e-12", "x=0"),
    ("calculus.integral", "compute", "integral(x^2, x)", "symbolic", "", ""),
    ("calculus.definite", "compute", "defint(x^2, x, 0, 1)", "symbolic", "", ""),
    ("calculus.definite_numeric", "compute", "defint(x^2, x, 0, 1)", "numeric", "1e-12", ""),
    ("calculus.limit", "compute", "limit(sin(x)/x, x, 0)", "symbolic", "", ""),
    ("calculus.sum", "compute", "sum(1/(k*(k+1)), k, 1, n)", "symbolic", "", ""),
    ("algebra.expand", "compute", "expand((x+1)^3)", "symbolic", "", ""),
    ("algebra.factor", "compute", "factor(x^2 - 1)", "symbolic", "", ""),
    ("algebra.solve", "compute", "solve(x^2 - 4, x)", "symbolic", "", ""),
    ("number_theory.factorint", "compute", "factorint(360)", "symbolic", "", ""),
    ("number_theory.gcd", "compute", "gcd(12, 18)", "symbolic", "", ""),
    ("linear_algebra.det", "compute", "det(list(list(1, 2), list(3, 4)))", "symbolic", "", ""),
    (
        "linear_algebra.matmul",
        "compute",
        "matmul(list(list(1,2),list(3,4)), list(list(0,1),list(1,0)))",
        "symbolic",
        "",
        "",
    ),
    ("numeric.pythagorean", "compute", "sin(x)^2 + cos(x)^2", "numeric", "1e-12", "x=0.7"),
    ("numeric.sqrt", "compute", "sqrt(x)", "numeric", "1e-12", "x=2"),
    ("simplify.ratio", "simplify", "(x^2 - 1)/(x - 1)", "symbolic", "", ""),
    ("arith.sub", "apply", "sub", "", "", "5,3"),
    ("arith.neg", "apply", "neg", "", "", "5"),
    ("arith.divide", "apply", "divide", "", "", "6,3"),
    # A higher-order operator: every runner supplies the same residual
    # `f(x) = x*x - 2`, so this row checks that all three wrappers hand a host
    # callback over the same way.
    ("callbacks.brent", "callback", "num_brent 1.0,2.0,1e-12,64", "", "", ""),
    ("errors.parse", "compute", "sin(", "symbolic", "", ""),
    ("errors.unknown_operator", "apply", "definitely_not_an_operator", "", "", "1"),
    ("errors.not_numeric", "compute", "diff(x^2, x)", "numeric", "1e-12", ""),
    (
        "sim.ising",
        "sim",
        'ising {"width":4,"height":4,"temperature":2.0,"mc_steps":5,"seed":7}',
        "",
        "",
        "",
    ),
    (
        "sim.poisson_multigrid",
        "sim",
        'poisson_multigrid {"dim":2,"n":9,"sigma":0.0,"cycle":"v","tol":1e-8}',
        "",
        "",
        "",
    ),
    ("manifest.hash", "manifest", "", "", "", ""),
    ("units.table", "units", "", "", "", ""),
]

STATUS = {
    "errors.parse": "parse",
    "errors.unknown_operator": "unknown_operator",
    "errors.not_numeric": "not_numeric",
}


def main() -> None:
    path = pathlib.Path(__file__).with_name("cases.tsv")
    lines = [
        "# Cross-language conformance cases for the rssn bindings.",
        "#",
        "# One row per case, tab-separated; `#` starts a comment. The Python runner",
        "# (`bindings/tests/run_conformance.py --record`) fills the expectation",
        "# columns by running the engine; the Python, Julia and C++ runners then all",
        "# check the same file, so a wrapper that disagrees with the engine — or with",
        "# another wrapper — fails.",
        "#",
        "# columns: " + ", ".join(f"{index + 1} {name}" for index, name in enumerate(HEADER)),
        "#",
        "\t".join(HEADER),
    ]
    for name, kind, payload, target, tolerance, bindings in CASES:
        status = STATUS.get(name, "ok")
        lines.append(
            "\t".join([name, kind, payload, target, tolerance, bindings, "", "", status])
        )
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines) + "\n")
    print(f"wrote {path} ({len(CASES)} cases)")


if __name__ == "__main__":
    main()
