"""Assemble the actual RSSN DE kernels and their local dependencies."""

import argparse
import re
from pathlib import Path

RSSN = Path(__file__).resolve().parents[4]


def item(source, head):
    # RSSN's top-level bodies have unindented closing braces; inner ones do not.
    match = re.search(rf"^{head}\b.*?^}}", source, re.M | re.S)
    if not match:
        raise ValueError(f"RSSN item missing: {head}")
    return match[0]


def generate(root):
    q = (root / "src/kernels/quadrature.rs").read_text()
    l = (root / "src/kernels/special_ext.rs").read_text()
    s = (root / "src/kernels/special.rs").read_text()
    k = (root / "src/kernels/mod.rs").read_text()
    positive, count = re.subn(
        r"^    if s > 0\.0 \{\n(.*?)^    }",
        r"\1",
        re.search(r"^    if s > 0\.0 \{\n.*?^    }", item(l, "pub fn lerch"),
                  re.M | re.S)[0],
        flags=re.M | re.S,
    )
    if count != 1 or "quadrature.value / gamma_numerical(s)" not in positive:
        raise ValueError("RSSN Lerch positive-order DE branch changed")
    appell = re.search(r"^    if c > a && a > 0\.0 \{\n.*?^    }",
                       item(l, "pub fn appell_f1"), re.M | re.S)
    if not appell or appell[0].count(".powf(") != 4 or "quadrature::tanh_sinh" not in appell[0]:
        raise ValueError("RSSN Appell F1 tanh-sinh branch changed")
    # The source branch has a trailing return because it sits inside lerch().
    parts = [
        "pub mod kernels {",
        item(k, r"pub\(crate\) fn exact_eq"),
        "pub mod quadrature { use std::f64::consts::PI;",
        "#[derive(Debug, Clone, Copy, PartialEq)]",
        item(q, "pub struct Integral"),
        *[item(q, x) for x in ("fn de_driver", "pub fn tanh_sinh",
                               "pub fn tanh_sinh_with_levels", "pub fn exp_sinh")],
        "}",
        "pub mod special_ext { use statrs::function::gamma::{gamma, ln_gamma};",
        *[item(s, x) for x in ("pub fn gamma_numerical", "pub fn ln_gamma_numerical")],
        *[item(l, x) for x in ("fn lerch_denominator", "fn lerch_small_positive_order",
                               "fn lerch_negative_order")],
        "pub fn small(z: f64, s: f64, a: f64) -> f64 { lerch_small_positive_order(z,s,a) }",
        "pub fn negative(z: f64, s: f64, a: f64) -> f64 { lerch_negative_order(z,s,a) }",
        "pub fn positive(z: f64, s: f64, a: f64) -> f64 {",
        positive,
        "}",
        "pub fn appell(a: f64, b1: f64, b2: f64, c: f64, x: f64, y: f64) -> f64 {",
        appell[0],
        "f64::NAN",
        "} } }",
    ]
    return "\n\n".join(parts) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--rssn", type=Path, default=RSSN)
    args = parser.parse_args()
    args.output.write_text(generate(args.rssn))
