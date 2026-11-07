# Symbolic audit

A pass over the kernel's symbolic answers: exactness, normal forms and coverage
of standard calculus. Every request in the table below is run through
`Session::compute` with `Config::new()` (standard rules, symbolic target) by
`tests/symbolic_battery.rs`, and the answer is checked against an independent
numeric oracle:

* antiderivatives are differentiated back (central differences) and compared
  with the integrand at five sample points;
* derivatives are compared with central differences of the input;
* simplifications must agree with the input numerically, and the cases that pin
  a normal form compare the printed text;
* limits, sums and definite integrals are compared with a known value, and the
  cases that do not exist (`limit(1/x, x, 0)`, `sum(1/k)`) must stay unevaluated;
* `solve` roots are substituted back into the equation; `taylor` is compared with
  the function near 0 to the truncation order.

Regenerate the table with
`AUDIT_TABLE=docs/table.md cargo test --test symbolic_battery -- --nocapture`.
`tests/symbolic_exactness.rs` holds the targeted regression tests for the defects
below.

Status words: `ok` (correct, and the pinned text when there is one), `WRONG`
(numeric disagreement), `FLOAT` (a float literal in a symbolic answer), `FORMAT`
(correct but not the pinned normal form), `UNEVALUATED` (no closed form found),
`UNCHECKED` (the oracle could not evaluate the answer, e.g. a real logarithm of a
negative number on the sample points).

## Findings

No numerically wrong answer was found in the battery (the two `WRONG` entries of
the first run were mistakes in the expected values of the battery itself). The
defects were exactness leaks, missed closed forms and normal forms.

### Fixed

| # | Defect | Root cause | Fix |
|---|---|---|---|
| 1 | `limit((1-cos t)/t^2)` "returns 0.5" | The kernel answers `1/2` for `Target::Symbolic`. The float is the literal `0.5` that `Target::Numeric` returns by contract, and `Answer::value` (`Some(0.5)`) which is filled for symbolic answers too. A caller that asks for a numeric target, or reads `Answer::value` instead of `Answer::term`, sees a float. | None needed in the kernel; documented. |
| 1b | `solve(x^4 - x - 1)` returned `list(-0.7244…, 1.2207…)` in a symbolic run | `heuristics::numeric_roots` (the fallback of the quartic formula) wrote bisection floats into the symbolic answer. | The roots are now exact `root_of(k, c0, …, cn)` objects; the float only orders them. |
| 2 | `limit(t*ln(t), t, 0)`, `limit(x^2*ln(x), x, 0)` unevaluated | The Gruntz core needs both sides of the point; the left side of a real `ln` is undefined, so it gave up although the one-sided limit exists. | A function that is undefined (NaN) on one side takes the limit of the side where it exists. |
| 2b | `limit((1+x)^(1/x), x, 0)`, `(cos x)^(1/x^2)`, `x^(1/(1-x))` unevaluated | The continuity shortcut evaluated `1^oo` with IEEE `pow` (= 1), got a "finite" value, and the numeric check then rejected it, ending the search. | Powers that are `1^oo`, `0^0` or `oo^0` at the point are not "continuous"; the search goes on to the `exp(lim e ln b)` form. |
| 2c | `limit(ln(x)/sqrt(x), x, oo)` unevaluated | The numeric check compared the claimed limit 0 with the probe at `x = 2^26`, where the value is still 0.002. | A slowly converging limit is accepted when the distance to the claim keeps shrinking along 2^10 … 2^50. |
| 2d | `limit(sqrt(x^2+x) - x, x, oo)` unevaluated | Catastrophic cancellation blinds the probe, and the series fallback was finite-point only. | The fallback now expands `f(±1/t)` about `t = 0+`. |
| 3 | `defint(sin(t)/t, 0, 1)` = `si(1) - si(0)`; `erf(oo)`, `si(oo)`, `zeta(0)`, … | No rule for the special values. | Rules for `si(0)`, `shi(0)`, `li(0)`, parity of `si`/`shi`, `erf(oo)`, `erfc(±oo)`, `si(oo)`, `ci(oo)`, `ei(-oo)`, `besselj/besseli(n, 0)`; `zeta` at 0, the negative integers and every even integer up to 60 (Bernoulli numbers); `digamma` at half-integers. |
| 4 | `integral(sqrt(1-t^2))` = `… - 1/2*asin(-t)` | The quadratic-radical integrator wrote `-asin(u)/sqrt(-a)` with a negative slope in `u`; the engine had already stopped when the odd-function rewrite produced the alternative spelling, and the two spellings tie in cost. | The integrator writes `asin(-u)/sqrt(-a)` (positive slope), so the answer is `1/2*t*sqrt(1-t^2) + 1/2*asin(t)`. The parity rules (`sin`, `tan`, `asin`, `atan`, `sinh`, `tanh`, `asinh`, `erf`, `cos`, `cosh`, …) already existed and are pinned by tests. `antiderivative` and the algebraic integrator also return their answer in simplified form. |
| 5 | `2*(x^2 + 1)^-1*x` | Not produced by rssn: `Graph::display`, LaTeX, Typst and the Unicode printer all print `2*x/(x^2 + 1)`, `\frac{…}{…}`, `frac(…)`. The `^-1` text is rendered by `aevia-core/src/sym.rs` (`Term::render`, the `Pow` arm), which writes `b^n` for every integer `n`. | Regression tests in rssn; the fix belongs to aevia-core. |
| 6 | `1/120*x^4 - 1/5040*x^6 - 1/6*x^2 + 1` | Two things: the printer's "lead with a positive term" rule is applied to a sum whose terms are in descending order, which scrambles a four-term series; and no ordering said that a series reads upwards. | `taylor` / `laurent` results are marked as truncated series (`Graph::mark_ascending`, carried by `Graph::import`) and print lowest power first: `1 - 1/6*x^2 + 1/120*x^4 - 1/5040*x^6`. Ordinary polynomials stay descending. All printers (text, LaTeX, Typst, pretty) share `ordered_terms`. |
| 7 | `solve(exp(t) - 2)` has empty `multiplicities` | `Solutions.multiplicities` was empty whenever the solver did not distinguish the roots (and `to_value` printed an empty list). | The list always has one entry per root: the exact multiplicity where it is known (rational roots; square-free polynomials; irrational roots of a repeated factor matched to the exact square-free part with a margin), otherwise the default `1`, a lower bound. |
| 8 | `(x+y)^2 - x^2 - 2*x*y - y^2` and `(x+1)^2 - x^2 - 2x - 1` not zero | The collapse kernel only looked for products of sums, not powers of sums. | Powers `(a+b)^n`, `2 <= n <= 8`, count as hidden sums. |
| 9 | `sin(x + pi/2)`, `cos(x + pi)`, `sin(pi - x)`, `exp(ln x + ln y)` | Only shifts by `n*pi` existed. | Shifts and reflections by `pi` and `pi/2`; `exp(a + ln b) = exp(a)*b`. |
| 10 | `integral(exp(x)*sin(x)*cos(x))`, `1/(exp(x)+exp(-x))`, `exp(x)/(1+exp(2x))`, `sech(x)^2`, `x^x*(1+ln x)`, `abs(x)` | Missing methods. | Double angle when other factors are present; substitution `w = exp(b x)` for commensurable exponent slopes; reciprocal hyperbolic functions through `cosh`/`sinh` on retry; a power with variable exponent is a substitution candidate; `integral(abs(u))`. |
| 11 | `sum(k/2^k)`, `sum(k^2/2^k)`, `sum(z^k/k!)`, `sum((-1)^k/k^s)`, `sum(1/((2k-1)(2k+1)))` | Missing closed forms. | Polynomial times rational geometric (Stirling numbers), exponential series, Dirichlet eta; the half-integer digamma values close the last one. |

### Open

* `integral(sin(x^2), x)` and `integral(x^x, x)` have no elementary antiderivative
  and no Fresnel operator exists; they stay unevaluated (`reduced = false`).
* `limit(sqrt(x+1) - sqrt(x), x, oo)` needs half-integer valuations in the series
  engine and stays unevaluated.
* `ln(x*y)` for positive `x`, `y` stays `ln(x*y)`; the expanded `ln(x) + ln(y)` is
  larger. `sqrt(x^2)` stays without an assumption (correct for complex `x`).
* The spelling of a square root is chosen by cost: `sqrt(1-x^2)` when the input said
  `sqrt`, `(1-x^2)^(1/2)` when the engine derived it. Both are equal terms.
* `solve` is a real solver: `x^2 + 1` has no roots, `x^4 - 1` has `-1, 1`, and a
  periodic equation lists the principal root with `complete = false` in
  `solve_structured`.
* Antiderivatives are complex-analytic (`ln(x - 1)`, not `ln|x - 1|`); `real_integral`
  gives the absolute-value form.

## Printing conventions

* A product with negative integer powers prints as a quotient: `2*x/(x^2 + 1)`;
  LaTeX `\frac`, Typst `frac(…)`.
* A rational coefficient is written in front: `1/120*x^5` (not `x^5/120`).
* A polynomial prints in descending degree, leading with a positive term when there
  is one (`1 - x^2`); a truncated series (`taylor`, `laurent`) prints in ascending
  degree.

## Cases

| # | status | request | answer |
|---|---|---|---|
| 1 | ok | `integral(x^3 - 2*x + 5, x)` | `1/4*x^4 - x^2 + 5*x` |
| 2 | ok | `integral(1/x, x)` | `ln(x)` |
| 3 | ok | `integral(1/(x^2 + 1), x)` | `atan(x)` |
| 4 | ok | `integral(1/(x^2 + 4), x)` | `1/2*atan(1/2*x)` |
| 5 | ok | `integral(1/(x^2 - 1), x)` | `1/2*ln(x - 1) - 1/2*ln(x + 1)` |
| 6 | ok | `integral(1/(x*(x + 1)), x)` | `ln(x) - ln(x + 1)` |
| 7 | ok | `integral(1/(x^2 + 3*x + 2), x)` | `ln(x + 1) - ln(x + 2)` |
| 8 | ok | `integral(1/((x + 1)^2*(x + 2)), x)` | `ln(x + 2) - ln(x + 1) - 1/(x + 1)` |
| 9 | ok | `integral((x + 1)/(x^2 + 2*x + 5), x)` | `1/2*ln(x^2 + 2*x + 5)` |
| 10 | ok | `integral(x/(x^2 + 1), x)` | `1/2*ln(x^2 + 1)` |
| 11 | ok | `integral(x^2/(x^2 + 1), x)` | `x - atan(x)` |
| 12 | ok | `integral(1/(x^3 + 1), x)` | `atan((2*x - 1)/3^(1/2))/3^(1/2) + 1/3*ln(x + 1) - 1/6*ln(x^2 - x + 1)` |
| 13 | ok | `integral(1/(x^4 - 1), x)` | `1/4*ln(x - 1) - 1/2*atan(x) - 1/4*ln(x + 1)` |
| 14 | ok | `integral(1/(x^4 + 1), x)` | `1/8*2^(1/2)*ln(x^2 + 2^(1/2)*x + 1) - 1/8*2^(1/2)*ln(x^2 - 2^(1/2)*x + 1) + 1/2*atan(2^(1/2)*x - 1)/2^(1/2) + 1/2*atan(2^(1/2)*x + 1)/2^(1/2)` |
| 15 | ok | `integral((3*x + 2)/((x - 1)*(x + 3)), x)` | `5/4*ln(x - 1) + 7/4*ln(x + 3)` |
| 16 | ok | `integral(1/(x*(x^2 + 1)), x)` | `ln(x) - 1/2*ln(x^2 + 1)` |
| 17 | ok | `integral((x^3 + x)/(x^2 + 1)^2, x)` | `1/2*ln(x^2 + 1)` |
| 18 | ok | `integral(1/(x^2 + 2*x + 2)^2, x)` | `1/2*atan(x + 1) + (1/2*x + 1/2)/(x^2 + 2*x + 2)` |
| 19 | ok | `integral(sin(x), x)` | `-cos(x)` |
| 20 | ok | `integral(sin(x)^2, x)` | `1/2*x - 1/4*sin(2*x)` |
| 21 | ok | `integral(cos(x)^2, x)` | `1/2*x + 1/4*sin(2*x)` |
| 22 | ok | `integral(sin(x)*cos(x), x)` | `-1/2*cos(x)^2` |
| 23 | ok | `integral(sin(x)^3, x)` | `1/3*cos(x)^3 - cos(x)` |
| 24 | ok | `integral(cos(x)^3, x)` | `sin(x) - 1/3*sin(x)^3` |
| 25 | ok | `integral(sin(x)^4, x)` | `3/8*x - 1/4*sin(2*x) + 1/32*sin(4*x)` |
| 26 | ok | `integral(sin(2*x)*cos(3*x), x)` | `1/2*cos(x) - 1/10*cos(5*x)` |
| 27 | ok | `integral(sin(x)*sin(2*x), x)` | `1/2*sin(x) - 1/6*sin(3*x)` |
| 28 | ok | `integral(cos(x)*cos(3*x), x)` | `1/4*sin(2*x) + 1/8*sin(4*x)` |
| 29 | ok | `integral(tan(x), x)` | `-ln(cos(x))` |
| 30 | ok | `integral(tan(x)^2, x)` | `tan(x) - x` |
| 31 | ok | `integral(1/cos(x)^2, x)` | `tan(x)` |
| 32 | ok | `integral(1/sin(x), x)` | `ln(tan(1/2*x))` |
| 33 | ok | `integral(1/cos(x), x)` | `ln(tan(1/2*x) + 1) - ln(tan(1/2*x) - 1)` |
| 34 | ok | `integral(1/(1 + cos(x)), x)` | `tan(1/2*x)` |
| 35 | ok | `integral(1/(2 + sin(x)), x)` | `2*atan((2*tan(1/2*x) + 1)/3^(1/2))/3^(1/2)` |
| 36 | ok | `integral(exp(x), x)` | `exp(x)` |
| 37 | ok | `integral(exp(-x^2), x)` | `1/2*erf(x)*pi^(1/2)` |
| 38 | ok | `integral(exp(2*x)*sin(x), x)` | `1/5*(2*sin(x) - cos(x))*exp(2*x)` |
| 39 | ok | `integral(x*exp(x), x)` | `x*exp(x) - exp(x)` |
| 40 | ok | `integral(x^2*exp(x), x)` | `x^2*exp(x) - 2*x*exp(x) + 2*exp(x)` |
| 41 | ok | `integral(x*exp(-x), x)` | `-x*exp(-x) - exp(-x)` |
| 42 | ok | `integral(x*sin(x), x)` | `sin(x) - x*cos(x)` |
| 43 | ok | `integral(x*cos(x), x)` | `x*sin(x) + cos(x)` |
| 44 | ok | `integral(x^2*sin(x), x)` | `2*x*sin(x) - x^2*cos(x) + 2*cos(x)` |
| 45 | ok | `integral(x*ln(x), x)` | `1/2*x^2*ln(x) - 1/4*x^2` |
| 46 | ok | `integral(ln(x), x)` | `x*ln(x) - x` |
| 47 | ok | `integral(ln(x)^2, x)` | `x*ln(x)^2 - 2*x*ln(x) + 2*x` |
| 48 | ok | `integral(x^2*ln(x), x)` | `1/3*x^3*ln(x) - 1/9*x^3` |
| 49 | ok | `integral(atan(x), x)` | `x*atan(x) - 1/2*ln(x^2 + 1)` |
| 50 | ok | `integral(asin(x), x)` | `x*asin(x) + (1 - x^2)^(1/2)` |
| 51 | ok | `integral(x*atan(x), x)` | `1/2*x^2*atan(x) - 1/2*x + 1/2*atan(x)` |
| 52 | ok | `integral(exp(x)*cos(x), x)` | `1/2*(cos(x) + sin(x))*exp(x)` |
| 53 | ok | `integral(x*exp(x^2), x)` | `1/2*exp(x^2)` |
| 54 | ok | `integral(2*x*cos(x^2), x)` | `sin(x^2)` |
| 55 | ok | `integral(cos(x)/sin(x), x)` | `ln(sin(x))` |
| 56 | ok | `integral(ln(x)/x, x)` | `1/2*ln(x)^2` |
| 57 | ok | `integral(1/(x*ln(x)), x)` | `ln(ln(x))` |
| 58 | ok | `integral(exp(x)/(1 + exp(x)), x)` | `ln(exp(x) + 1)` |
| 59 | ok | `integral(1/(1 + exp(x)), x)` | `x - ln(exp(x) + 1)` |
| 60 | ok | `integral(sinh(x), x)` | `cosh(x)` |
| 61 | ok | `integral(cosh(x)^2, x)` | `1/2*x - 1/8*exp(-2*x) + 1/8*exp(2*x)` |
| 62 | ok | `integral(sqrt(x), x)` | `2/3*sqrt(x)^3` |
| 63 | ok | `integral(1/sqrt(x), x)` | `2*sqrt(x)` |
| 64 | ok | `integral(sqrt(1 - x^2), x)` | `1/2*x*sqrt(1 - x^2) + 1/2*asin(x)` |
| 65 | ok | `integral(1/sqrt(1 - x^2), x)` | `asin(x)` |
| 66 | ok | `integral(1/sqrt(x^2 + 1), x)` | `ln(2*x + 2*sqrt(x^2 + 1))` |
| 67 | ok | `integral(sqrt(x^2 + 1), x)` | `1/2*x*sqrt(x^2 + 1) + 1/2*ln(2*x + 2*sqrt(x^2 + 1))` |
| 68 | ok | `integral(x/sqrt(x^2 + 1), x)` | `sqrt(x^2 + 1)` |
| 69 | ok | `integral(1/sqrt(x^2 - 1), x)` | `ln(2*x + 2*sqrt(x^2 - 1))` |
| 70 | ok | `integral(x*sqrt(x + 1), x)` | `2/5*sqrt(x + 1)^5 - 2/3*sqrt(x + 1)^3` |
| 71 | ok | `integral(sin(x)/x, x)` | `si(x)` |
| 72 | ok | `integral(exp(x)/x, x)` | `ei(x)` |
| 73 | ok | `integral(1/ln(x), x)` | `li(x)` |
| 74 | UNEVALUATED | `integral(sin(x^2), x)` | `integral(sin(x^2), x)` |
| 75 | UNEVALUATED | `integral(x^x, x)` | `integral(x^x, x)` |
| 76 | ok | `diff(x^5, x)` | `5*x^4` |
| 77 | ok | `diff(sin(x)*cos(x), x)` | `cos(2*x)` |
| 78 | ok | `diff(exp(x^2), x)` | `2*x*exp(x^2)` |
| 79 | ok | `diff(ln(x^2 + 1), x)` | `2*x/(x^2 + 1)` |
| 80 | ok | `diff(atan(x), x)` | `1/(x^2 + 1)` |
| 81 | ok | `diff(asin(x), x)` | `1/(1 - x^2)^(1/2)` |
| 82 | ok | `diff(x^x, x)` | `ln(x)*x^x + x^x` |
| 83 | ok | `diff(sqrt(1 + x^2), x)` | `x/(x^2 + 1)^(1/2)` |
| 84 | ok | `diff(tan(x), x)` | `tan(x)^2 + 1` |
| 85 | ok | `diff(sinh(x)/x, x)` | `cosh(x)/x - sinh(x)/x^2` |
| 86 | ok | `diff(x/(x^2 + 1), x)` | `1/(x^2 + 1) - 2*x^2/(x^2 + 1)^2` |
| 87 | ok | `diff(exp(sin(x)), x)` | `cos(x)*exp(sin(x))` |
| 88 | ok | `diff(ln(sin(x)), x)` | `cot(x)` |
| 89 | ok | `diff(abs(x)^3, x)` | `3*x*abs(x)` |
| 90 | ok | `diff(atan(x)/x, x)` | `1/(x^3 + x) - atan(x)/x^2` |
| 91 | ok | `diff(acosh(x + 1), x)` | `1/((x + 1)^2 - 1)^(1/2)` |
| 92 | ok | `diff(gamma(x), x)` | `digamma(x)*gamma(x)` |
| 93 | ok | `diff(erf(x), x)` | `2*exp(-x^2)/pi^(1/2)` |
| 94 | ok | `sin(x)^2 + cos(x)^2` | `1` |
| 95 | ok | `exp(ln(x))` | `x` |
| 96 | ok | `ln(exp(x))` | `x` |
| 97 | ok | `(x^2 - 1)/(x - 1)` | `x + 1` |
| 98 | ok | `sin(-x)` | `-sin(x)` |
| 99 | ok | `cos(-x)` | `cos(x)` |
| 100 | ok | `tan(-x)` | `-tan(x)` |
| 101 | ok | `asin(-x)` | `-asin(x)` |
| 102 | ok | `atan(-x)` | `-atan(x)` |
| 103 | ok | `sinh(-x)` | `-sinh(x)` |
| 104 | ok | `cosh(-x)` | `cosh(x)` |
| 105 | ok | `tanh(-x)` | `-tanh(x)` |
| 106 | ok | `asinh(-x)` | `-asinh(x)` |
| 107 | ok | `erf(-x)` | `-erf(x)` |
| 108 | ok | `sqrt(x^2)` | `sqrt(x^2)` |
| 109 | ok | `(x + 1)^2 - (x^2 + 2*x + 1)` | `0` |
| 110 | ok | `(x^3 - 1)/(x^2 + x + 1)` | `x - 1` |
| 111 | ok | `cos(x)^2 - sin(x)^2` | `cos(2*x)` |
| 112 | ok | `2*sin(x)*cos(x)` | `sin(2*x)` |
| 113 | ok | `exp(x)*exp(2*x)` | `exp(3*x)` |
| 114 | ok | `ln(x^2)` | `ln(x^2)` |
| 115 | ok | `ln(x*y)` | `ln(x*y)` |
| 116 | ok | `x/x` | `1` |
| 117 | ok | `(x^2 - 4)/(x + 2)` | `x - 2` |
| 118 | ok | `1/(1/x)` | `x` |
| 119 | ok | `sqrt(x)*sqrt(x)` | `x` |
| 120 | ok | `exp(x)^2` | `exp(2*x)` |
| 121 | ok | `sin(x)/cos(x)` | `tan(x)` |
| 122 | ok | `cosh(x)^2 - sinh(x)^2` | `1` |
| 123 | ok | `(x + y)^2 - x^2 - 2*x*y - y^2` | `0` |
| 124 | ok | `si(0)` | `0` |
| 125 | ok | `erf(0)` | `0` |
| 126 | ok | `gamma(1)` | `1` |
| 127 | ok | `gamma(3)` | `2` |
| 128 | ok | `gamma(5)` | `24` |
| 129 | ok | `gamma(1/2)` | `pi^(1/2)` |
| 130 | ok | `zeta(2)` | `1/6*pi^2` |
| 131 | ok | `zeta(4)` | `1/90*pi^4` |
| 132 | ok | `ci(1) - ci(1)` | `0` |
| 133 | ok | `li(0)` | `0` |
| 134 | ok | `ei(-oo)` | `0` |
| 135 | ok | `erf(oo)` | `1` |
| 136 | ok | `erfc(oo)` | `0` |
| 137 | ok | `erfc(0)` | `1` |
| 138 | ok | `ln(1)` | `0` |
| 139 | ok | `asin(0)` | `0` |
| 140 | ok | `si(oo)` | `1/2*pi` |
| 141 | ok | `gamma(3)` | `2` |
| 142 | ok | `limit((1 - cos(x))/x^2, x, 0)` | `1/2` |
| 143 | ok | `limit((x^2 + 1)/(2*x^2 - x), x, oo)` | `1/2` |
| 144 | ok | `limit(x*ln(x), x, 0, plus)` | `0` |
| 145 | ok | `limit(x*ln(x), x, 0)` | `0` |
| 146 | ok | `limit(x^x, x, 0, plus)` | `1` |
| 147 | ok | `limit((1 + 1/x)^x, x, oo)` | `E` |
| 148 | ok | `limit(x*exp(-x), x, oo)` | `0` |
| 149 | ok | `limit(ln(x)/x, x, oo)` | `0` |
| 150 | ok | `limit(sin(x)/x, x, oo)` | `0` |
| 151 | ok | `limit((exp(x) - 1)/x, x, 0)` | `1` |
| 152 | ok | `limit(x^2*ln(x), x, 0, plus)` | `0` |
| 153 | ok | `limit(sin(x)/x, x, 0)` | `1` |
| 154 | ok | `limit((exp(x) - 1 - x)/x^2, x, 0)` | `1/2` |
| 155 | ok | `limit((tan(x) - x)/x^3, x, 0)` | `1/3` |
| 156 | ok | `limit(x^3*exp(-x), x, oo)` | `0` |
| 157 | ok | `limit(x^(1/x), x, oo)` | `1` |
| 158 | ok | `limit(x^(1/x), x, oo)` | `1` |
| 159 | ok | `limit((x^2 - 4)/(x - 2), x, 2)` | `4` |
| 160 | ok | `limit(exp(x)/x^10, x, oo)` | `oo` |
| 161 | ok | `limit(1/x, x, oo)` | `0` |
| 162 | ok | `limit(x*sin(1/x), x, oo)` | `1` |
| 163 | ok | `limit((1 + x)^(1/x), x, 0)` | `E` |
| 164 | ok | `limit(sqrt(x^2 + x) - x, x, oo)` | `1/2` |
| 165 | ok | `limit(atan(x)/x, x, 0)` | `1` |
| 166 | ok | `limit(atan(x), x, oo)` | `1/2*pi` |
| 167 | ok | `limit(ln(1 + x)/exp(x), x, oo)` | `0` |
| 168 | ok | `limit(sinh(2*x)/x, x, 0)` | `2` |
| 169 | ok | `limit((1 - cos(x))/(x*sin(x)), x, 0)` | `1/2` |
| 170 | ok | `limit(ln(x)/(1 - x), x, 1)` | `-1` |
| 171 | ok | `defint(2*x, x, 0, 1)` | `1` |
| 172 | ok | `defint(1/(1 + x^2), x, 0, 1)` | `1/4*pi` |
| 173 | ok | `defint(sin(x), x, 0, pi/2)` | `1` |
| 174 | ok | `defint(sin(x)/x, x, 0, 1)` | `si(1)` |
| 175 | ok | `defint(exp(-x), x, 0, oo)` | `1` |
| 176 | ok | `defint(exp(-x^2), x, 0, oo)` | `1/2*pi^(1/2)` |
| 177 | ok | `defint(1/(1 + x^2), x, 0, oo)` | `1/2*pi` |
| 178 | ok | `defint(ln(x), x, 0, 1)` | `-1` |
| 179 | ok | `defint(x*ln(x), x, 0, 1)` | `-1/4` |
| 180 | ok | `defint(sqrt(1 - x^2), x, 0, 1)` | `1/4*pi` |
| 181 | ok | `defint(x*exp(-x^2), x, 0, oo)` | `1/2` |
| 182 | ok | `defint(x*exp(-x), x, 0, oo)` | `1` |
| 183 | ok | `defint(sin(x), x, 0, pi)` | `2` |
| 184 | ok | `defint(sin(x)^2, x, 0, pi)` | `1/2*pi` |
| 185 | ok | `defint(x^2, x, 0, 1)` | `1/3` |
| 186 | ok | `defint(1/(1 + x^2), x, -oo, oo)` | `pi` |
| 187 | ok | `defint(x*exp(-x), x, 0, 1)` | `1 - 2*exp(-1)` |
| 188 | ok | `sum(1/k^2, k, 1, oo)` | `1/6*pi^2` |
| 189 | ok | `sum(1/2^k, k, 0, oo)` | `2` |
| 190 | ok | `sum(1/(k*(k + 1)), k, 1, oo)` | `1` |
| 191 | ok | `sum(k, k, 1, 10)` | `55` |
| 192 | ok | `sum(k^2, k, 1, 10)` | `385` |
| 193 | ok | `sum(k^3, k, 1, 10)` | `3025` |
| 194 | ok | `sum(2^k, k, 0, 9)` | `1023` |
| 195 | ok | `sum(1/factorial(k), k, 0, oo)` | `E` |
| 196 | ok | `sum(1/3^k, k, 0, oo)` | `3/2` |
| 197 | ok | `sum(1/k^4, k, 1, oo)` | `1/90*pi^4` |
| 198 | ok | `sum(1/((k)*(k + 2)), k, 1, oo)` | `3/4` |
| 199 | ok | `sum(k/2^k, k, 1, oo)` | `2` |
| 200 | ok | `taylor(sin(x)/x, x, 0, 6)` | `1 - 1/6*x^2 + 1/120*x^4 - 1/5040*x^6` |
| 201 | ok | `taylor(sin(x), x, 0, 7)` | `x - 1/6*x^3 + 1/120*x^5 - 1/5040*x^7` |
| 202 | ok | `taylor(exp(x), x, 0, 6)` | `1 + x + 1/2*x^2 + 1/6*x^3 + 1/24*x^4 + 1/120*x^5 + 1/720*x^6` |
| 203 | ok | `taylor(cos(x), x, 0, 6)` | `1 - 1/2*x^2 + 1/24*x^4 - 1/720*x^6` |
| 204 | ok | `taylor(1/(1 - x), x, 0, 6)` | `1 + x + x^2 + x^3 + x^4 + x^5 + x^6` |
| 205 | ok | `taylor(ln(1 + x), x, 0, 6)` | `x - 1/2*x^2 + 1/3*x^3 - 1/4*x^4 + 1/5*x^5 - 1/6*x^6` |
| 206 | ok | `taylor(atan(x), x, 0, 6)` | `x - 1/3*x^3 + 1/5*x^5` |
| 207 | ok | `taylor(1/(1 + x^2), x, 0, 6)` | `1 - x^2 + x^4 - x^6` |
| 208 | ok | `taylor(exp(x)*cos(x), x, 0, 6)` | `1 + x - 1/3*x^3 - 1/6*x^4 - 1/30*x^5` |
| 209 | ok | `taylor(sqrt(1 + x), x, 0, 6)` | `1 + 1/2*x - 1/8*x^2 + 1/16*x^3 - 5/128*x^4 + 7/256*x^5 - 21/1024*x^6` |
| 210 | ok | `taylor(tan(x), x, 0, 6)` | `x + 1/3*x^3 + 2/15*x^5` |
| 211 | ok | `solve(x^2 - 4, x)` | `list(-2, 2)` |
| 212 | ok | `solve(x^2 - 2, x)` | `list(-2^(1/2), 2^(1/2))` |
| 213 | ok | `solve(x^2 + 1, x)` | `list()` |
| 214 | ok | `solve(x^3 - 6*x^2 + 11*x - 6, x)` | `list(1, 2, 3)` |
| 215 | ok | `solve(x^3 - 2, x)` | `list(2^(1/3))` |
| 216 | ok | `solve(x^4 - 5*x^2 + 4, x)` | `list(-2, -1, 1, 2)` |
| 217 | ok | `solve(x^4 - 1, x)` | `list(-1, 1)` |
| 218 | ok | `solve(x^4 + x^2 - 2, x)` | `list(-1, 1)` |
| 219 | ok | `solve(2*x + 3, x)` | `list(-3/2)` |
| 220 | ok | `solve(x^2 - 2*x + 1, x)` | `list(1)` |
| 221 | ok | `solve(exp(x) - 2, x)` | `list(ln(2))` |
| 222 | ok | `solve(ln(x) - 1, x)` | `list(E)` |
| 223 | ok | `solve(sin(x) - 1/2, x)` | `list(1/6*pi)` |
| 224 | ok | `solve(cos(x), x)` | `list(1/2*pi)` |
| 225 | ok | `solve(2^x - 8, x)` | `list(ln(8)/ln(2))` |
| 226 | ok | `solve(sqrt(x) - 3, x)` | `list(9)` |
| 227 | ok | `solve(x^2 - 3*x + 2, x)` | `list(1, 2)` |
| 228 | ok | `solve(x^3 + x, x)` | `list(0)` |
| 229 | ok | `solve(tan(x) - 1, x)` | `list(1/4*pi)` |
| 230 | ok | `solve(exp(2*x) - 3*exp(x) + 2, x)` | `list(0, ln(2))` |
| 231 | ok | `integral(1/(x^2 + 1)^2, x)` | `1/2*atan(x) + 1/2*x/(x^2 + 1)` |
| 232 | ok | `integral(x/(x^2 + 1)^2, x)` | `-1/2/(x^2 + 1)` |
| 233 | ok | `integral(x^2/(x^2 + 1)^2, x)` | `1/2*atan(x) - 1/2*x/(x^2 + 1)` |
| 234 | ok | `integral(1/(x^3 - x), x)` | `1/2*ln(x - 1) - ln(x) + 1/2*ln(x + 1)` |
| 235 | ok | `integral((2*x + 3)/(x^2 + 4*x + 5), x)` | `ln(x^2 + 4*x + 5) - atan(x + 2)` |
| 236 | ok | `integral(1/(x^2 - 2*x + 5), x)` | `1/2*atan(1/2*x - 1/2)` |
| 237 | ok | `integral(x^3/(x^2 + 1), x)` | `1/2*x^2 - 1/2*ln(x^2 + 1)` |
| 238 | ok | `integral(1/(x^2*(x - 1)), x)` | `ln(x - 1) - ln(x) + 1/x` |
| 239 | ok | `integral(1/(x^3*(x + 1)), x)` | `ln(x) - ln(x + 1) + 1/x - 1/2/x^2` |
| 240 | ok | `integral((x + 2)/(x^2 - 4), x)` | `ln(x - 2)` |
| 241 | ok | `integral((x^2 + 1)/(x^4 + 1), x)` | `atan(2^(1/2)*x - 1)/2^(1/2) + atan(2^(1/2)*x + 1)/2^(1/2)` |
| 242 | ok | `integral(1/(x^4 + 5*x^2 + 4), x)` | `1/3*atan(x) - 1/6*atan(1/2*x)` |
| 243 | ok | `integral((x^5 + 1)/(x^2 - 1), x)` | `1/4*x^4 + 1/2*x^2 + ln(x - 1)` |
| 244 | UNCHECKED | `integral(1/((x - 1)*(x - 2)*(x - 3)), x)` | `1/2*ln(x - 3) - ln(x - 2) + 1/2*ln(x - 1)` |
| 245 | ok | `integral(1/(x^2 + x + 1), x)` | `2*atan((2*x + 1)/3^(1/2))/3^(1/2)` |
| 246 | ok | `integral(sin(x)^2*cos(x), x)` | `1/3*sin(x)^3` |
| 247 | ok | `integral(sin(x)*cos(x)^2, x)` | `-1/3*cos(x)^3` |
| 248 | ok | `integral(sin(x)^5, x)` | `2/3*cos(x)^3 - 1/5*cos(x)^5 - cos(x)` |
| 249 | ok | `integral(cos(x)^4, x)` | `3/8*x + 1/4*sin(2*x) + 1/32*sin(4*x)` |
| 250 | ok | `integral(tan(x)^3, x)` | `1/2*tan(x)^2 - 1/2*ln(tan(x)^2 + 1)` |
| 251 | ok | `integral(1/cos(x)^3, x)` | `1/2*ln(tan(1/2*x) + 1) - 1/2*ln(tan(1/2*x) - 1) + 1/2*sin(x)/cos(x)^2` |
| 252 | ok | `integral(sin(3*x), x)` | `-1/3*cos(3*x)` |
| 253 | ok | `integral(sin(x)*cos(2*x), x)` | `1/2*cos(x) - 1/6*cos(3*x)` |
| 254 | ok | `integral(x*sin(x)^2, x)` | `1/4*x^2 - 1/4*x*sin(2*x) - 1/8*cos(2*x)` |
| 255 | ok | `integral(sin(x)/cos(x)^2, x)` | `sec(x)` |
| 256 | ok | `integral(cos(x)/(1 + sin(x)^2), x)` | `atan(sin(x))` |
| 257 | ok | `integral(1/(sin(x)*cos(x)), x)` | `ln(tan(1/2*x)) - ln(tan(1/2*x) - 1) - ln(tan(1/2*x) + 1)` |
| 258 | ok | `integral(1/(1 + sin(x)), x)` | `-2/(tan(1/2*x) + 1)` |
| 259 | ok | `integral(1/(3 + 2*cos(x)), x)` | `4*atan(2*tan(1/2*x)/20^(1/2))/20^(1/2)` |
| 260 | ok | `integral(sin(x)^2*cos(x)^2, x)` | `1/8*x - 1/32*sin(4*x)` |
| 261 | ok | `integral(cot(x), x)` | `ln(sin(x))` |
| 262 | ok | `integral(1/sin(x)^2, x)` | `-cot(x)` |
| 263 | ok | `integral(sec(x)*tan(x), x)` | `sec(x)` |
| 264 | ok | `integral(exp(-x)*cos(2*x), x)` | `1/5*(2*sin(2*x) - cos(2*x))*exp(-x)` |
| 265 | ok | `integral(x*exp(2*x), x)` | `1/2*x*exp(2*x) - 1/4*exp(2*x)` |
| 266 | ok | `integral(x^3*exp(x), x)` | `x^3*exp(x) - 3*x^2*exp(x) + 6*x*exp(x) - 6*exp(x)` |
| 267 | ok | `integral(x^2*exp(-x), x)` | `-x^2*exp(-x) - 2*x*exp(-x) - 2*exp(-x)` |
| 268 | ok | `integral(sqrt(x)*ln(x), x)` | `2/3*ln(x)*sqrt(x)^3 - 4/9*sqrt(x)^3` |
| 269 | ok | `integral(x^2*atan(x), x)` | `1/3*x^3*atan(x) - 1/6*x^2 + 1/6*ln(x^2 + 1)` |
| 270 | ok | `integral(ln(x^2 + 1), x)` | `x*ln(x^2 + 1) - 2*x + 2*atan(x)` |
| 271 | ok | `integral(x*ln(x^2 + 1), x)` | `1/2*((x^2 + 1)*ln(x^2 + 1) - x^2 - 1)` |
| 272 | ok | `integral(atan(x)/(1 + x^2), x)` | `1/2*atan(x)^2` |
| 273 | ok | `integral(asin(x)/sqrt(1 - x^2), x)` | `1/2*asin(x)^2` |
| 274 | ok | `integral(ln(x)/x^2, x)` | `(-ln(x) - 1)/x` |
| 275 | ok | `integral(x*acos(x), x)` | `1/2*x^2*acos(x) - 1/4*x*(1 - x^2)^(1/2) + 1/4*asin(x)` |
| 276 | ok | `integral(exp(x)*sin(x)*cos(x), x)` | `1/10*(sin(2*x) - 2*cos(2*x))*exp(x)` |
| 277 | ok | `integral(exp(2*x)/(1 + exp(x)), x)` | `exp(x) - ln(exp(x) + 1)` |
| 278 | ok | `integral(1/(exp(x) + exp(-x)), x)` | `atan(exp(x))` |
| 279 | ok | `integral(sech(x)^2, x)` | `-2/(exp(2*x) + 1)` |
| 280 | ok | `integral(tanh(x), x)` | `ln(cosh(x))` |
| 281 | ok | `integral(x*sinh(x), x)` | `1/2*x*exp(x) + 1/2*x*exp(-x) - 1/2*exp(x) + 1/2*exp(-x)` |
| 282 | ok | `integral(exp(sqrt(x)), x)` | `2*exp(sqrt(x))*sqrt(x) - 2*exp(sqrt(x))` |
| 283 | ok | `integral(ln(ln(x))/x, x)` | `ln(x)*ln(ln(x)) - ln(x)` |
| 284 | ok | `integral(x*exp(-x^2/2), x)` | `-exp(-1/2*x^2)` |
| 285 | ok | `integral(2^x, x)` | `2^x/ln(2)` |
| 286 | ok | `integral(x*2^x, x)` | `x*2^x/ln(2) - 2^x/ln(2)^2` |
| 287 | ok | `integral(x^x*(1 + ln(x)), x)` | `x^x` |
| 288 | ok | `integral(x*sqrt(x^2 + 1), x)` | `1/3*sqrt(x^2 + 1)^3` |
| 289 | ok | `integral(x^2*sqrt(x^3 + 1), x)` | `2/9*(x^3 + 1)^(3/2)` |
| 290 | ok | `integral(cos(x)*exp(sin(x)), x)` | `exp(sin(x))` |
| 291 | ok | `integral(exp(x)/(1 + exp(2*x)), x)` | `atan(exp(x))` |
| 292 | ok | `integral(1/(x*sqrt(x^2 - 1)), x)` | `-2*atan(sqrt(x^2 - 1)/(x - 1))` |
| 293 | ok | `integral(1/(1 + sqrt(x)), x)` | `2*sqrt(x) - 2*ln(sqrt(x) + 1)` |
| 294 | ok | `integral(1/(sqrt(x) + x), x)` | `2*ln(1/2*sqrt(x) + 1/2)` |
| 295 | ok | `integral(sqrt(1 + x)/x, x)` | `ln(abs(sqrt(x + 1) - 1)) - ln(abs(sqrt(x + 1) + 1)) + 2*sqrt(x + 1)` |
| 296 | ok | `integral(1/(x*sqrt(1 - x^2)), x)` | `-x*ln(2*(1/x^2 - 1)^(1/2) + 2/x)/abs(x)` |
| 297 | ok | `integral(sqrt(x^2 - 1), x)` | `1/2*x*sqrt(x^2 - 1) - 1/2*ln(2*x + 2*sqrt(x^2 - 1))` |
| 298 | ok | `integral(1/(x^2*sqrt(x^2 + 1)), x)` | `-sqrt(x^2 + 1)/x` |
| 299 | ok | `integral(sqrt(4 - x^2), x)` | `1/2*x*sqrt(4 - x^2) + 2*asin(1/2*x)` |
| 300 | ok | `integral(sqrt(x)/(1 + x), x)` | `2*sqrt(x) - 2*atan(sqrt(x))` |
| 301 | ok | `integral(x/sqrt(1 - x^2), x)` | `-sqrt(1 - x^2)` |
| 302 | ok | `integral(1/sqrt(2*x - x^2), x)` | `asin(x - 1)` |
| 303 | ok | `integral(abs(x), x)` | `1/2*x*abs(x)` |
| 304 | ok | `diff(x^2*sin(x)*exp(x), x)` | `x^2*cos(x)*exp(x) + x^2*exp(x)*sin(x) + 2*x*exp(x)*sin(x)` |
| 305 | ok | `diff(sin(cos(x)), x)` | `-cos(cos(x))*sin(x)` |
| 306 | ok | `diff(ln(x + sqrt(x^2 + 1)), x)` | `(1 + x/(x^2 + 1)^(1/2))/(x + sqrt(x^2 + 1))` |
| 307 | ok | `diff(atan(x^2), x)` | `2*x/(x^4 + 1)` |
| 308 | ok | `diff(x^(1/3), x)` | `1/3/x^(2/3)` |
| 309 | ok | `diff((x^2 + 1)^(3/2), x)` | `3*x*(x^2 + 1)^(1/2)` |
| 310 | ok | `diff(exp(x)/(1 + x^2), x)` | `exp(x)/(x^2 + 1) - 2*x*exp(x)/(x^2 + 1)^2` |
| 311 | ok | `diff(sin(x)^cos(x), x)` | `cos(x)^2*sin(x)^(cos(x) - 1) - ln(sin(x))*sin(x)^(cos(x) + 1)` |
| 312 | ok | `diff(asinh(x), x)` | `1/(x^2 + 1)^(1/2)` |
| 313 | ok | `diff(atanh(x), x)` | `1/(1 - x^2)` |
| 314 | ok | `diff(acos(x)*x, x)` | `acos(x) - x/(1 - x^2)^(1/2)` |
| 315 | ok | `diff(ln(abs(x)), x)` | `x/abs(x)^2` |
| 316 | ok | `diff(sqrt(x)*ln(x), x)` | `1/2*ln(x)/x^(1/2) + 1/x^(1/2)` |
| 317 | ok | `diff(tanh(x)^2, x)` | `2*tanh(x) - 2*tanh(x)^3` |
| 318 | ok | `diff(sec(x), x)` | `sec(x)*tan(x)` |
| 319 | ok | `diff(1/(x^2 + 1)^2, x)` | `-4*x/(x^2 + 1)^3` |
| 320 | ok | `diff(exp(exp(x)), x)` | `exp(x + exp(x))` |
| 321 | ok | `diff(x*exp(-x^2), x)` | `exp(-x^2) - 2*x^2*exp(-x^2)` |
| 322 | ok | `diff(si(x), x)` | `sin(x)/x` |
| 323 | ok | `diff(erfc(x), x)` | `-2*exp(-x^2)/pi^(1/2)` |
| 324 | ok | `diff(besselj(0, x), x)` | `besselj(-1, x)` |
| 325 | ok | `diff(lgamma(x), x)` | `digamma(x)` |
| 326 | ok | `exp(-ln(x))` | `1/x` |
| 327 | ok | `exp(ln(x) + ln(y))` | `x*y` |
| 328 | ok | `ln(exp(2))` | `2` |
| 329 | ok | `exp(ln(2) + ln(3))` | `6` |
| 330 | ok | `sin(pi/6)` | `1/2` |
| 331 | ok | `cos(pi/3)` | `1/2` |
| 332 | ok | `sin(x + pi/2)` | `cos(x)` |
| 333 | ok | `cos(x + pi)` | `-cos(x)` |
| 334 | ok | `sin(x + 2*pi)` | `sin(x)` |
| 335 | ok | `sin(pi - x)` | `sin(x)` |
| 336 | ok | `cos(2*x) + 2*sin(x)^2` | `2*sin(x)^2 + cos(2*x)` |
| 337 | ok | `tan(x)^2 + 1` | `tan(x)^2 + 1` |
| 338 | ok | `1/cos(x)^2 - tan(x)^2` | `1/cos(x)^2 - tan(x)^2` |
| 339 | ok | `(a + b)*(a - b)` | `(a + b)*(a - b)` |
| 340 | ok | `(x + 1)^3 - x^3 - 3*x^2 - 3*x - 1` | `0` |
| 341 | ok | `(x^2 - 1)/(x^2 + 2*x + 1)` | `(x^2 - 1)/(x^2 + 2*x + 1)` |
| 342 | ok | `x/(x^2 + x)` | `1/(x + 1)` |
| 343 | ok | `1/x + 1/(x + 1)` | `1/x + 1/(x + 1)` |
| 344 | ok | `x/(x + 1) - 1/(x + 1)` | `x/(x + 1) - 1/(x + 1)` |
| 345 | ok | `sqrt(x)^2` | `x` |
| 346 | ok | `(x^(1/2))^4` | `x^2` |
| 347 | ok | `abs(x)^2` | `abs(x)^2` |
| 348 | ok | `exp(x)*exp(-x)` | `1` |
| 349 | ok | `ln(x) - ln(x)` | `0` |
| 350 | ok | `2^x*2^x` | `2^(2*x)` |
| 351 | ok | `(2^x)^2` | `2^(2*x)` |
| 352 | ok | `sinh(x) + cosh(x)` | `cosh(x) + sinh(x)` |
| 353 | ok | `(exp(x) - exp(-x))/2` | `1/2*(exp(x) - exp(-x))` |
| 354 | ok | `(exp(x) + exp(-x))/2` | `1/2*(exp(x) + exp(-x))` |
| 355 | ok | `sin(x)^4 - cos(x)^4` | `sin(x)^4 - cos(x)^4` |
| 356 | ok | `x^3/x` | `x^2` |
| 357 | ok | `x^a*x^b` | `x^(a + b)` |
| 358 | ok | `(x^a)^b` | `(x^a)^b` |
| 359 | ok | `log(2, 8)` | `3` |
| 360 | ok | `log(10, 100)` | `2` |
| 361 | ok | `atan(0)` | `0` |
| 362 | ok | `atan(1)` | `1/4*pi` |
| 363 | ok | `acos(1/2)` | `1/3*pi` |
| 364 | ok | `cosh(0)` | `1` |
| 365 | ok | `sinh(0)` | `0` |
| 366 | ok | `cos(pi)` | `-1` |
| 367 | ok | `gamma(2)` | `1` |
| 368 | ok | `gamma(4)` | `6` |
| 369 | ok | `gamma(5/2)` | `3/4*pi^(1/2)` |
| 370 | ok | `gamma(-1/2)` | `-2*pi^(1/2)` |
| 371 | ok | `zeta(0)` | `-1/2` |
| 372 | ok | `zeta(-1)` | `-1/12` |
| 373 | ok | `zeta(-3)` | `1/120` |
| 374 | ok | `zeta(-2)` | `0` |
| 375 | ok | `zeta(6)` | `1/945*pi^6` |
| 376 | ok | `shi(0)` | `0` |
| 377 | ok | `erfi(0)` | `0` |
| 378 | ok | `sinc(0)` | `1` |
| 379 | ok | `besselj(0, 0)` | `1` |
| 380 | ok | `besselj(1, 0)` | `0` |
| 381 | ok | `limit((1 + 2/x)^x, x, oo)` | `exp(2)` |
| 382 | ok | `limit(x^(1/ln(x)), x, oo)` | `E` |
| 383 | ok | `limit(ln(x)/sqrt(x), x, oo)` | `0` |
| 384 | ok | `limit(exp(1/x), x, 0, plus)` | `oo` |
| 385 | ok | `limit(exp(-1/x^2), x, 0)` | `0` |
| 386 | ok | `limit(exp(1/x), x, 0, minus)` | `0` |
| 387 | ok | `limit(x*(sqrt(x^2 + 1) - x), x, oo)` | `1/2` |
| 388 | ok | `limit(atan(x), x, -oo)` | `-1/2*pi` |
| 389 | ok | `limit(tanh(x), x, oo)` | `1` |
| 390 | ok | `limit((x - sin(x))/x^3, x, 0)` | `1/6` |
| 391 | ok | `limit((2^x - 1)/x, x, 0)` | `ln(2)` |
| 392 | ok | `limit((sin(x) - x*cos(x))/x^3, x, 0)` | `1/3` |
| 393 | ok | `limit(ln(1 + x)/x, x, 0)` | `1` |
| 394 | ok | `limit(x^2/exp(x), x, oo)` | `0` |
| 395 | ok | `limit((cos(x) - 1)/x^2, x, 0)` | `-1/2` |
| 396 | ok | `limit(x/sin(x), x, 0)` | `1` |
| 397 | ok | `limit((3*x^2 + x)/(x^2 + 5), x, oo)` | `3` |
| 398 | ok | `limit((x + 1)/(x^2 + 1), x, oo)` | `0` |
| 399 | ok | `limit((x^3 + 1)/(x^2 + 1), x, oo)` | `oo` |
| 400 | ok | `limit(-x^2, x, oo)` | `-oo` |
| 401 | ok | `limit(x*exp(-x^2), x, oo)` | `0` |
| 402 | ok | `limit(exp(x)/(exp(x) + 1), x, oo)` | `1` |
| 403 | ok | `limit(exp(x)/(exp(x) + 1), x, -oo)` | `0` |
| 404 | ok | `limit(1/(1 + exp(-x)), x, 0)` | `1/2` |
| 405 | ok | `limit(x^(x^2), x, 1)` | `1` |
| 406 | ok | `limit(sin(x)*exp(-x), x, oo)` | `0` |
| 407 | ok | `limit(sin(x)/(x - pi), x, pi)` | `-1` |
| 408 | ok | `limit((sqrt(x) - 2)/(x - 4), x, 4)` | `1/4` |
| 409 | ok | `limit(1/x, x, 0)` | `limit(1/x, x, 0)` |
| 410 | ok | `limit(abs(x)/x, x, 0)` | `limit(abs(x)/x, x, 0)` |
| 411 | ok | `limit(sin(x), x, oo)` | `limit(sin(x), x, oo)` |
| 412 | ok | `limit(sin(1/x), x, 0)` | `limit(sin(1/x), x, 0)` |
| 413 | ok | `sum((-1)^(k + 1)/k, k, 1, oo)` | `ln(2)` |
| 414 | ok | `sum(1/(k*(k + 1)*(k + 2)), k, 1, oo)` | `1/4` |
| 415 | ok | `sum(1/((2*k - 1)*(2*k + 1)), k, 1, oo)` | `1/2` |
| 416 | ok | `sum((-1)^(k + 1)/k^2, k, 1, oo)` | `1/12*pi^2` |
| 417 | ok | `sum(1/k^3, k, 1, oo)` | `zeta(3)` |
| 418 | ok | `sum(1/2^k, k, 1, oo)` | `1` |
| 419 | ok | `sum(k^2/2^k, k, 1, oo)` | `6` |
| 420 | ok | `sum((-1)^k/factorial(k), k, 0, oo)` | `exp(-1)` |
| 421 | ok | `sum(1, k, 1, 10)` | `10` |
| 422 | ok | `sum(k, k, 1, 100)` | `5050` |
| 423 | ok | `sum(1/2^k, k, 1, 10)` | `1023/1024` |
| 424 | ok | `sum(1/k, k, 1, oo)` | `sum(1/k, k, 1, oo)` |
| 425 | ok | `sum(k, k, 1, oo)` | `sum(k, k, 1, oo)` |
| 426 | ok | `taylor(1/cos(x), x, 0, 6)` | `1 + 1/2*x^2 + 5/24*x^4 + 61/720*x^6` |
| 427 | ok | `taylor(exp(-x^2), x, 0, 6)` | `1 - x^2 + 1/2*x^4 - 1/6*x^6` |
| 428 | ok | `taylor(asin(x), x, 0, 5)` | `x + 1/6*x^3 + 3/40*x^5` |
| 429 | ok | `taylor(sinh(x), x, 0, 5)` | `x + 1/6*x^3 + 1/120*x^5` |
| 430 | ok | `taylor(ln(1 + x)/(1 + x), x, 0, 5)` | `x - 3/2*x^2 + 11/6*x^3 - 25/12*x^4 + 137/60*x^5` |
| 431 | ok | `taylor((1 + x)^(1/3), x, 0, 5)` | `1 + 1/3*x - 1/9*x^2 + 5/81*x^3 - 10/243*x^4 + 22/729*x^5` |
| 432 | ok | `taylor(exp(sin(x)), x, 0, 5)` | `1 + x + 1/2*x^2 - 1/8*x^4 - 1/15*x^5` |
| 433 | ok | `taylor(x/(exp(x) - 1), x, 0, 5)` | `1 - 1/2*x + 1/12*x^2 - 1/720*x^4` |
| 434 | ok | `taylor(sqrt(1 - x^2), x, 0, 5)` | `1 - 1/2*x^2 - 1/8*x^4` |
| 435 | ok | `taylor(cos(x)^2, x, 0, 5)` | `1 - x^2 + 1/3*x^4` |
| 436 | ok | `solve(x^3 - 3*x + 1, x)` | `list(2*cos(10/9*pi), 2*cos(4/9*pi), 2*cos(2/9*pi))` |
| 437 | ok | `solve(x^4 - 10*x^2 + 1, x)` | `list(-(2^(1/2) + 3^(1/2)), 2^(1/2) - 3^(1/2), 3^(1/2) - 2^(1/2), 2^(1/2) + 3^(1/2))` |
| 438 | ok | `solve(x^3 - x^2 - 2*x + 2, x)` | `list(-2^(1/2), 1, 2^(1/2))` |
| 439 | ok | `solve(x^4 - 16, x)` | `list(-2, 2)` |
| 440 | ok | `solve(x^4 + 4, x)` | `list()` |
| 441 | ok | `solve(x^2 - 5*x + 6, x)` | `list(2, 3)` |
| 442 | ok | `solve(x^3 + 3*x^2 + 3*x + 1, x)` | `list(-1)` |
| 443 | ok | `solve(4*x^3 - 3*x - 1/2, x)` | `list(cos(11/9*pi), cos(5/9*pi), cos(1/9*pi))` |
| 444 | ok | `solve(2*cos(x) - 1, x)` | `list(1/3*pi)` |
| 445 | ok | `solve(sin(x), x)` | `list(0)` |
| 446 | ok | `solve(tan(x) - sqrt(3), x)` | `list(1/3*pi)` |
| 447 | ok | `solve(exp(x) + exp(-x) - 3, x)` | `list(ln(3/2 - 1/2*5^(1/2)), ln(3/2 + 1/2*5^(1/2)))` |
| 448 | ok | `solve(ln(x) + ln(x - 1) - ln(2), x)` | `list(2)` |
| 449 | ok | `solve(abs(x - 1) - 2, x)` | `list(-1, 3)` |
| 450 | ok | `solve(sqrt(x + 1) - x + 1, x)` | `list(3)` |
| 451 | ok | `solve(1/x + 1/(x + 1) - 1, x)` | `list(1/2 - 1/2*5^(1/2), 1/2 + 1/2*5^(1/2))` |
| 452 | ok | `solve(x^2 - 2*x - 1, x)` | `list(1 - 2^(1/2), 1 + 2^(1/2))` |
| 453 | ok | `solve(3^x - 9, x)` | `list(ln(9)/ln(3))` |
| 454 | ok | `solve(x*exp(x) - 1, x)` | `list(lambertw(1))` |
| 455 | ok | `solve(sin(x)^2 - 1/4, x)` | `list(-1/6*pi, 1/6*pi)` |
| 456 | ok | `solve(x^5 - x, x)` | `list(-1, 0, 1)` |
| 457 | ok | `solve(x^6 - 1, x)` | `list(-1, 1)` |
| 458 | ok | `solve(2*x^2 - 8, x)` | `list(-2, 2)` |
| 459 | ok | `solve(x^2 + 2*x + 1, x)` | `list(-1)` |
