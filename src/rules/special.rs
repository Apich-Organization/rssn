//! Special functions: gamma and beta, error functions, zeta, Bessel
//! functions, orthogonal polynomials and the step-like functions.
//!
//! Numeric values come from [`crate::kernels::special`] (which wraps
//! `statrs`) except `zeta`, whose only available implementation is far too
//! coarse and is replaced by a short Euler-Maclaurin sum. Exact values are
//! kernels on literal arguments; orthogonal polynomials of literal degree
//! expand to explicit polynomials with exact rational coefficients, so the
//! differentiation and integration kernels see them as ordinary
//! polynomials.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Facts;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OnReals;
use crate::graph::OpDescriptor;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::graph::op::EvalFn;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::kernels::special as num;

use super::calculus::at_infinity;
use super::calculus::calculus;
use super::calculus::integral_rule;
use super::calculus::partials;
use crate::graph::Cx as KernelCx;
use crate::rules::poly::repr::from_term;
use crate::rules::poly::repr::to_term;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;
use super::elementary::elementary;

/// Highest degree of an orthogonal polynomial that is expanded.
const MAX_DEGREE: u32 = 64;
/// Largest integer or half-integer whose gamma value is computed exactly.
const MAX_GAMMA: u32 = 1000;

/// The special-function rule set.
#[must_use]
pub fn special() -> RuleSet {
    RuleSet::new("special", install)
        .needs(elementary())
        .needs(calculus())
}

/// `(n, x)` of a two-argument evaluation with a non-negative integer
/// degree `n`.
fn degree_and_point(args: &[f64]) -> Option<(u32, f64)> {
    let (&n, &x) = (args.first()?, args.get(1)?);
    (n >= 0.0 && n <= f64::from(u32::MAX) && n.fract() == 0.0).then(|| (n.to_u32().unwrap_or(0), x))
}

/// A two-argument evaluation `f(order, x)`.
fn order_and_point(
    args: &[f64],
    f: fn(f64, f64) -> f64,
) -> f64 {
    match args {
        | [n, x] => f(*n, *x),
        | _ => f64::NAN,
    }
}

/// `ln |gamma(x)|`.
fn lgamma(x: f64) -> f64 {
    if x > 0.0 {
        num::ln_gamma_numerical(x)
    } else {
        num::gamma_numerical(x).abs().ln()
    }
}

/// Beta function on all reals where it is defined.
fn beta(
    a: f64,
    b: f64,
) -> f64 {
    if a > 0.0 && b > 0.0 {
        num::beta_numerical(a, b)
    } else {
        num::gamma_numerical(a) * num::gamma_numerical(b) / num::gamma_numerical(a + b)
    }
}

/// Bernoulli corrections `B_{2j} / (2j)!` for `j = 1..=7`.
const ZETA_CORRECTIONS: [f64; 7] = [
    1.0 / 12.0,
    -1.0 / 720.0,
    1.0 / 30_240.0,
    -1.0 / 1_209_600.0,
    1.0 / 47_900_160.0,
    -691.0 / 1_307_674_368_000.0,
    1.0 / 74_724_249_600.0,
];
/// Terms of the zeta series summed directly.
const ZETA_TERMS: u32 = 12;

/// Riemann zeta for real `s > 1` by Euler-Maclaurin summation: the first
/// terms directly, the tail by its integral plus Bernoulli corrections.
fn zeta(s: f64) -> f64 {
    // Euler–Maclaurin summation is also the analytic continuation for
    // s < 1, accurate for moderate |s|.
    if s.is_nan() || s < -20.0 {
        return f64::NAN;
    }
    if s.total_cmp(&1.0).is_eq() {
        return f64::INFINITY;
    }
    let n = f64::from(ZETA_TERMS);
    let mut sum: f64 = (1..ZETA_TERMS).map(|k| f64::from(k).powf(-s)).sum();
    sum += 0.5_f64.mul_add(n.powf(-s), n.powf(1.0 - s) / (s - 1.0));
    // s (s + 1) ... (s + 2j - 2) * n^(-s - 2j + 1)
    let mut rising = s;
    let mut power = n.powf(-s - 1.0);
    let mut k = 1.0;
    for c in ZETA_CORRECTIONS {
        sum = (c * rising).mul_add(power, sum);
        rising *= (s + k) * (s + k + 1.0);
        power /= n * n;
        k += 2.0;
    }
    sum
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    // name, arity, semantics, what is known of the value for real arguments
    let functions: [(&str, u8, EvalFn, Facts); 29] = [
        (
            "gamma",
            1,
            |a| a.first().map_or(f64::NAN, |&x| num::gamma_numerical(x)),
            Facts::REAL,
        ),
        (
            "lgamma",
            1,
            |a| a.first().map_or(f64::NAN, |&x| lgamma(x)),
            Facts::REAL,
        ),
        (
            "digamma",
            1,
            |a| a.first().map_or(f64::NAN, |&x| num::digamma_numerical(x)),
            Facts::REAL,
        ),
        (
            "beta",
            2,
            |a| match a {
                | [x, y] => beta(*x, *y),
                | _ => f64::NAN,
            },
            Facts::REAL,
        ),
        (
            "erf",
            1,
            |a| a.first().map_or(f64::NAN, |&x| num::erf_numerical(x)),
            Facts::REAL,
        ),
        (
            "erfc",
            1,
            |a| a.first().map_or(f64::NAN, |&x| num::erfc_numerical(x)),
            Facts::REAL,
        ),
        (
            "zeta",
            1,
            |a| a.first().map_or(f64::NAN, |&x| zeta(x)),
            Facts::NONE,
        ),
        (
            "besselj",
            2,
            |a| order_and_point(a, num::bessel_j),
            Facts::NONE,
        ),
        (
            "bessely",
            2,
            |a| order_and_point(a, num::bessel_y),
            Facts::NONE,
        ),
        (
            "besseli",
            2,
            |a| order_and_point(a, num::bessel_i),
            Facts::NONE,
        ),
        (
            "besselk",
            2,
            |a| order_and_point(a, num::bessel_k),
            Facts::NONE,
        ),
        (
            "erfi",
            1,
            |a| a.first().map_or(f64::NAN, |&x| num::erfi(x)),
            Facts::REAL,
        ),
        (
            "erfcinv",
            1,
            |a| a.first().map_or(f64::NAN, |&x| num::inverse_erfc(x)),
            Facts::REAL,
        ),
        (
            "polygamma",
            2,
            |a| match a {
                | [n, x] if *n >= 0.0 && n.fract() == 0.0 && *n <= f64::from(u32::MAX) => {
                    num::polygamma_numerical(n.to_u32().unwrap_or(0), *x)
                },
                | _ => f64::NAN,
            },
            Facts::REAL,
        ),
        (
            "laguerre_gen",
            3,
            |a| match a {
                | [n, alpha, x] if *n >= 0.0 && n.fract() == 0.0 && *n <= 1e6 => {
                    num::generalized_laguerre(n.to_u32().unwrap_or(0), *alpha, *x)
                },
                | _ => f64::NAN,
            },
            Facts::REAL,
        ),
        (
            "legendre",
            2,
            |a| degree_and_point(a).map_or(f64::NAN, |(n, x)| num::legendre_p(n, x)),
            Facts::NONE,
        ),
        (
            "chebyshevt",
            2,
            |a| degree_and_point(a).map_or(f64::NAN, |(n, x)| num::chebyshev_t(n, x)),
            Facts::NONE,
        ),
        (
            "chebyshevu",
            2,
            |a| degree_and_point(a).map_or(f64::NAN, |(n, x)| num::chebyshev_u(n, x)),
            Facts::NONE,
        ),
        (
            "hermite",
            2,
            |a| degree_and_point(a).map_or(f64::NAN, |(n, x)| num::hermite_h(n, x)),
            Facts::NONE,
        ),
        (
            "laguerre",
            2,
            |a| degree_and_point(a).map_or(f64::NAN, |(n, x)| num::laguerre_l(n, x)),
            Facts::NONE,
        ),
        (
            "sinc",
            1,
            |a| a.first().map_or(f64::NAN, |&x| num::sinc(x)),
            Facts::REAL,
        ),
        (
            "heaviside",
            1,
            |a| a.first().map_or(f64::NAN, |&x| step(x, 0.0, 0.5, 1.0)),
            Facts::NONNEGATIVE,
        ),
        ("hurwitz_zeta", 2, |a| match a {
            | [s, q] => num::hurwitz_zeta(*s, *q),
            | _ => f64::NAN,
        }, Facts::NONE),
        ("lambertw", 1, |a| a.first().map_or(f64::NAN, |&x| num::lambert_w(x, true)), Facts::REAL),
        ("lambertw_m1", 1, |a| a.first().map_or(f64::NAN, |&x| num::lambert_w(x, false)), Facts::NEGATIVE),
        ("floor", 1, |a| a.first().map_or(f64::NAN, |x| x.floor()), Facts::INTEGER),
        ("ceil", 1, |a| a.first().map_or(f64::NAN, |x| x.ceil()), Facts::INTEGER),
        ("round", 1, |a| a.first().map_or(f64::NAN, |x| x.round()), Facts::INTEGER),
        (
            "sign",
            1,
            |a| a.first().map_or(f64::NAN, |&x| step(x, -1.0, 0.0, 1.0)),
            Facts::REAL,
        ),
    ];
    let mut registered = Vec::with_capacity(functions.len());
    for (name, arity, eval, on_reals) in functions {
        let op = i.op(OpDescriptor::new(name, Arity::Fixed(arity)).eval(eval))?;
        if on_reals != Facts::NONE {
            i.graph().ops_mut().set_attr(op, OnReals(on_reals));
        }
        registered.push((name, op));
    }
    let op_named = |name: &str| registered.iter().find(|r| r.0 == name).map(|r| r.1);
    let pi = i
        .graph()
        .ops()
        .lookup("pi")
        .ok_or_else(|| RuleError::Invalid {
            rule: "special".to_owned(),
            reason: "`pi` is not registered",
        })?;
    let (Some(gamma), Some(zeta), Some(sign), Some(heaviside)) = (
        op_named("gamma"),
        op_named("zeta"),
        op_named("sign"),
        op_named("heaviside"),
    ) else {
        return Err(RuleError::Invalid {
            rule: "special".to_owned(),
            reason: "operator missing",
        });
    };

    partials(i, "erf", &["2/pi^(1/2) * exp(-?a^2)"])?;
    partials(i, "erfc", &["-2/pi^(1/2) * exp(-?a^2)"])?;
    partials(i, "gamma", &["gamma(?a) * digamma(?a)"])?;
    partials(i, "lgamma", &["digamma(?a)"])?;
    partials(
        i,
        "beta",
        &[
            "beta(?a, ?b) * (digamma(?a) - digamma(?a + ?b))",
            "beta(?a, ?b) * (digamma(?b) - digamma(?a + ?b))",
        ],
    )?;
    partials(i, "sinc", &["(cos(pi * ?a) - sinc(?a)) / ?a"])?;
    partials(i, "digamma", &["polygamma(1, ?a)"])?;
    partials(i, "polygamma", &["", "polygamma(?a + 1, ?b)"])?;
    partials(i, "erfi", &["2/pi^(1/2) * exp(?a^2)"])?;
    partials(i, "lambertw", &["lambertw(?a)/(?a*(1 + lambertw(?a)))"])?;
    partials(i, "lambertw_m1", &["lambertw_m1(?a)/(?a*(1 + lambertw_m1(?a)))"])?;
    partials(i, "besselj", &["", "(besselj(?a - 1, ?b) - besselj(?a + 1, ?b))/2"])?;
    partials(i, "bessely", &["", "(bessely(?a - 1, ?b) - bessely(?a + 1, ?b))/2"])?;
    partials(i, "besseli", &["", "(besseli(?a - 1, ?b) + besseli(?a + 1, ?b))/2"])?;
    partials(i, "besselk", &["", "-(besselk(?a - 1, ?b) + besselk(?a + 1, ?b))/2"])?;
    at_infinity(i, "erfi", Some("oo"), Some("-oo"))?;
    at_infinity(i, "erf", Some("1"), Some("-1"))?;
    at_infinity(i, "erfc", Some("0"), Some("2"))?;
    at_infinity(i, "heaviside", Some("1"), Some("0"))?;
    at_infinity(i, "sign", Some("1"), Some("-1"))?;
    integral_rule(i, gaussian_integral);
    integral_rule(i, erf_integral);

    let families = [
        ("legendre", Family::Legendre),
        ("chebyshevt", Family::ChebyshevT),
        ("chebyshevu", Family::ChebyshevU),
        ("hermite", Family::Hermite),
        ("laguerre", Family::Laguerre),
    ];
    let mut expandable = Vec::with_capacity(families.len());
    for (name, family) in families {
        if let Some(op) = op_named(name) {
            expandable.push((op, family));
        }
    }
    i.kernel("special/exact", Tier::Normalize, Exact { gamma, zeta, pi });
    i.kernel("special/step", Tier::Normalize, Step { sign, heaviside });
    let euler_gamma =
        i.op(OpDescriptor::new("euler_gamma", Arity::Fixed(0)).cost(3).eval(|_| 0.577_215_664_901_532_9))?;
    i.graph().ops_mut().set_attr(euler_gamma, OnReals(Facts::POSITIVE));
    if let (Some(digamma), Some(polygamma)) = (op_named("digamma"), op_named("polygamma")) {
        i.kernel(
            "special/polygamma-exact",
            Tier::Normalize,
            PolygammaExact { digamma, polygamma, zeta, euler_gamma, pi },
        );
    }
    if let (Some(floor), Some(ceil), Some(round)) = (op_named("floor"), op_named("ceil"), op_named("round")) {
        i.kernel("special/rounding", Tier::Normalize, Rounding { floor, ceil, round });
    }

    i.kernel(
        "special/orthogonal",
        Tier::Normalize,
        Orthogonal { families: expandable },
    );

    i.rewrites(
        Tier::Normalize,
        &[
            "special/gamma-1: gamma(1) => 1",
            "special/gamma-half: gamma(1/2) => pi^(1/2)",
            "special/gamma-recurrence: gamma(1 + ?x) / gamma(?x) => ?x",
            "special/lgamma-1: lgamma(1) => 0",
            "special/lgamma-2: lgamma(2) => 0",
            "special/erf-0: erf(0) => 0",
            "special/erfc-0: erfc(0) => 1",
            "special/erf-odd: erf(-?x) => -erf(?x)",
            "special/erf-erfc: erfc(?x) + erf(?x) => 1",
            "special/sinc-0: sinc(0) => 1",
            "special/sinc-even: sinc(-?x) => sinc(?x)",
            "special/sign-odd: sign(-?x) => -sign(?x)",
            "special/erfi-0: erfi(0) => 0",
            "special/erfi-odd: erfi(-?x) => -erfi(?x)",
            "special/polygamma-0: polygamma(0, ?x) => digamma(?x)",
            "special/besselj-half: besselj(1/2, ?x) => (2/(pi*?x))^(1/2) * sin(?x) if positive(?x)",
            "special/besselj-minus-half: besselj(-1/2, ?x) => (2/(pi*?x))^(1/2) * cos(?x) if positive(?x)",
            "special/bessely-half: bessely(1/2, ?x) => -(2/(pi*?x))^(1/2) * cos(?x) if positive(?x)",
            "special/besseli-half: besseli(1/2, ?x) => (2/(pi*?x))^(1/2) * sinh(?x) if positive(?x)",
            "special/besselk-half: besselk(1/2, ?x) => (pi/(2*?x))^(1/2) * exp(-?x) if positive(?x)",
            "special/besselk-minus-half: besselk(-1/2, ?x) => (pi/(2*?x))^(1/2) * exp(-?x) if positive(?x)",
            "special/besselk-even: besselk(-?n, ?x) => besselk(?n, ?x)",
            "special/besselj-0-at-0: besselj(0, 0) => 1",
            "special/besseli-0-at-0: besseli(0, 0) => 1",
            "special/besselj-negative-order: besselj(?n, ?x) => (-1)^?n * besselj(-?n, ?x) if integer(?n), negative(?n)",
            "special/bessely-negative-order: bessely(?n, ?x) => (-1)^?n * bessely(-?n, ?x) if integer(?n), negative(?n)",
            "special/besseli-negative-order: besseli(?n, ?x) => besseli(-?n, ?x) if integer(?n), negative(?n)",
            "special/erfcinv-1: erfcinv(1) => 0",
            "special/lambertw-0: lambertw(0) => 0",
            "special/hurwitz-1: hurwitz_zeta(?s, 1) => zeta(?s)",
            "special/lambertw-e: lambertw(E) => 1",
            "special/lambertw-inverse: lambertw(?x*exp(?x)) => ?x if nonnegative(?x)",
            "special/digamma-half: digamma(1/2) => -euler_gamma - 2*ln(2)",
            "special/digamma-shift: digamma(?x + 1) - digamma(?x) => 1/?x",
            "special/floor-integer: floor(?n) => ?n if integer(?n)",
            "special/ceil-integer: ceil(?n) => ?n if integer(?n)",
            "special/round-integer: round(?n) => ?n if integer(?n)",
            "special/floor-shift: floor(?x + ?n) => floor(?x) + ?n if integer(?n)",
            "special/ceil-shift: ceil(?x + ?n) => ceil(?x) + ?n if integer(?n)",
            "special/round-shift: round(?x + ?n) => round(?x) + ?n if integer(?n)",
            "special/ceil-floor: ceil(-?x) => -floor(?x)",
            "special/floor-ceil: floor(-?x) => -ceil(?x)",
        ],
    )?;
    i.rewrites(
        Tier::Explore,
        &[
            "special/beta: beta(?a, ?b) => gamma(?a) * gamma(?b) / gamma(?a + ?b)",
            "special/lgamma: lgamma(?x) => ln(gamma(?x)) if positive(?x)",
        ],
    )?;
    if let Some(op) = op_named("laguerre_gen") {
        i.kernel("special/laguerre-gen", Tier::Normalize, GeneralizedLaguerre { op });
    }
    i.define(&[
        "frac(x) := x - floor(x)",
        "ln_factorial(n) := lgamma(n + 1)",
        "bessel_k0(x) := besselk(0, x)",
        "bessel_k1(x) := besselk(1, x)",
        "inverse_erfc(x) := erfcinv(x)",
        "legendre_rodrigues(n, x) := diffn((x^2 - 1)^n, x, n) / (2^n * gamma(n + 1))",
        "hermite_rodrigues(n, x) := (-1)^n * exp(x^2) * diffn(exp(-x^2), x, n)",
        "laguerre_rodrigues(n, x) := exp(x) / gamma(n + 1) * diffn(exp(-x) * x^n, x, n)",
        "bessel_differential_equation(y, x, n) := x^2 * diff(diff(y, x), x) + x * diff(y, x) + (x^2 - n^2) * y = 0",
        "legendre_differential_equation(y, x, n) := (1 - x^2) * diff(diff(y, x), x) - 2 * x * diff(y, x) + n * (n + 1) * y = 0",
        "laguerre_differential_equation(y, x, n) := x * diff(diff(y, x), x) + (1 - x) * diff(y, x) + n * y = 0",
        "hermite_differential_equation(y, x, n) := diff(diff(y, x), x) - 2 * x * diff(y, x) + 2 * n * y = 0",
        "chebyshev_differential_equation(y, x, n) := (1 - x^2) * diff(diff(y, x), x) - x * diff(y, x) + n^2 * y = 0",
    ])
}

/// `L_n^{(α)}(x) = Σ_k (-1)^k C(n + α, n - k) x^k / k!` for a literal
/// degree `n`, with any `α`.
struct GeneralizedLaguerre {
    op: OpId,
}

impl Kernel for GeneralizedLaguerre {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let &[degree, alpha, x] = graph.children(node) else {
            return Outcome::Pass;
        };
        let Some(n) = integer(graph, degree).and_then(|d| d.to_u32()).filter(|&d| d <= MAX_DEGREE) else {
            return Outcome::Pass;
        };
        let mut terms = Vec::new();
        for k in 0..=n {
            // C(n + α, n - k) = Π_{j=1}^{n-k} (α + k + j) / (n - k)!
            let mut factors = Vec::new();
            for j in 1..=n - k {
                let shift = graph.int(i64::from(k + j));
                factors.push(graph.node(core::ADD, &[alpha, shift]));
            }
            let denominator = factorial(n - k) * factorial(k);
            let sign = if k % 2 == 0 { BigInt::one() } else { -BigInt::one() };
            let c = graph.num(Number::rat(BigRational::new(sign, denominator)));
            factors.push(c);
            if k > 0 {
                let e = graph.int(i64::from(k));
                factors.push(graph.node(core::POW, &[x, e]));
            }
            terms.push(if factors.len() == 1 { factors[0] } else { graph.node(core::MUL, &factors) });
        }
        let sum = if terms.len() == 1 { terms[0] } else { graph.node(core::ADD, &terms) };
        let polynomial = cx.simplify(sum);
        // Like the other orthogonal polynomials: the expansion is the
        // requested form, a number at a numeric point.
        if cx.graph.number_of(polynomial).is_some() { Outcome::Equal(polynomial) } else { Outcome::Pinned(polynomial) }
    }
}

/// `below`, `at` or `above` according to the sign of `x`.
fn step(
    x: f64,
    below: f64,
    at: f64,
    above: f64,
) -> f64 {
    match x.partial_cmp(&0.0) {
        | Some(std::cmp::Ordering::Less) => below,
        | Some(std::cmp::Ordering::Equal) => at,
        | Some(std::cmp::Ordering::Greater) => above,
        | None => f64::NAN,
    }
}

/// A literal integer.
fn integer(
    graph: &Graph,
    node: NodeId,
) -> Option<BigInt> {
    match graph.number_of(node)? {
        | Number::Int(v) => Some(v.clone()),
        | _ => None,
    }
}

fn factorial(n: u32) -> BigInt {
    (1..=n).fold(BigInt::one(), |acc, k| acc * k)
}

/// `coefficient * pi^exponent` as a term.
fn times_pi_power(
    graph: &mut Graph,
    pi: OpId,
    coefficient: &BigRational,
    exponent: Number,
) -> NodeId {
    let base = graph.node(pi, &[]);
    let power_of = graph.num(exponent);
    let power = graph.node(core::POW, &[base, power_of]);
    if coefficient.is_one() {
        power
    } else {
        let c = graph.num(Number::rat(coefficient.clone()));
        graph.node(core::MUL, &[c, power])
    }
}

/// `gamma(x)` as `(c, half)` meaning `c` or `c * sqrt(pi)`, for integers
/// and half-integers.
fn gamma_exact(x: &BigRational) -> Option<(BigRational, bool)> {
    let limit = BigInt::from(MAX_GAMMA);
    if x.is_integer() {
        let n = x.to_integer();
        if !n.is_positive() || n > limit {
            return None;
        }
        return Some((
            BigRational::from_integer(factorial(n.to_u32()?.saturating_sub(1))),
            false,
        ));
    }
    if *x.denom() != BigInt::from(2) {
        return None;
    }
    let m = x.numer();
    if m.abs() > BigInt::from(2 * MAX_GAMMA) {
        return None;
    }
    let value = if m.is_positive() {
        // Gamma(k + 1/2) = (2k)! / (4^k k!) sqrt(pi)
        let k = ((m - BigInt::one()) / BigInt::from(2)).to_u32()?;
        BigRational::new(factorial(2 * k), BigInt::from(4).pow(k) * factorial(k))
    } else {
        // Gamma(1/2 - k) = (-4)^k k! / (2k)! sqrt(pi)
        let k = ((BigInt::one() - m) / BigInt::from(2)).to_u32()?;
        BigRational::new(BigInt::from(-4).pow(k) * factorial(k), factorial(2 * k))
    };
    Some((value, true))
}

/// `zeta(2k)` for `k = 1..=6` as the rational `r` in `r * pi^(2k)`.
fn zeta_even(n: &BigInt) -> Option<BigRational> {
    let (numer, denom): (i64, i64) = match n.to_u32()? {
        | 2 => (1, 6),
        | 4 => (1, 90),
        | 6 => (1, 945),
        | 8 => (1, 9450),
        | 10 => (1, 93_555),
        | 12 => (691, 638_512_875),
        | _ => return None,
    };
    Some(BigRational::new(numer.into(), denom.into()))
}

/// Exact values of `gamma` at integers and half-integers and of `zeta` at
/// small even integers.
struct Exact {
    gamma: OpId,
    zeta: OpId,
    pi: OpId,
}

impl Kernel for Exact {
    fn ops(&self) -> Vec<OpId> {
        vec![self.gamma, self.zeta]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let Some(&arg) = graph.children(node).first() else {
            return Outcome::Pass;
        };
        if graph.op(node) == self.zeta {
            let Some(n) = integer(graph, arg) else {
                return Outcome::Pass;
            };
            let Some(coefficient) = zeta_even(&n) else {
                return Outcome::Pass;
            };
            return Outcome::Pinned(times_pi_power(graph, self.pi, &coefficient, Number::Int(n)));
        }
        let Some((value, half)) = graph
            .number_of(arg)
            .and_then(Number::to_rational)
            .and_then(|x| gamma_exact(&x))
        else {
            return Outcome::Pass;
        };
        if half {
            let exponent = Number::fraction(1, 2).unwrap_or_else(|| 0.into());
            Outcome::Pinned(times_pi_power(graph, self.pi, &value, exponent))
        } else {
            Outcome::Equal(graph.num(Number::rat(value)))
        }
    }

    fn revisit(&self) -> bool {
        true
    }
}

/// `sign` and `heaviside` where the sign of the argument is known.
/// `digamma(n) = H_{n-1} - γ` and
/// `polygamma(m, n) = (-1)^{m+1} m! (ζ(m+1) - Σ_{j<n} j^{-(m+1)})` at
/// positive integers `n`.
struct PolygammaExact {
    digamma: OpId,
    polygamma: OpId,
    zeta: OpId,
    euler_gamma: OpId,
    pi: OpId,
}

impl Kernel for PolygammaExact {
    fn ops(&self) -> Vec<OpId> {
        vec![self.digamma, self.polygamma]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let args = graph.children(node).to_vec();
        let (order, at) = match args.as_slice() {
            | &[x] if graph.op(node) == self.digamma => (0, x),
            | &[m, x] => match graph.number_of(m).and_then(Number::to_i64) {
                | Some(m) if (1..=60).contains(&m) => (m, x),
                | _ => return Outcome::Pass,
            },
            | _ => return Outcome::Pass,
        };
        let Some(n) = graph.number_of(at).and_then(Number::to_i64).filter(|n| (1..=10_000).contains(n)) else {
            return Outcome::Pass;
        };
        let power = order.saturating_add(1);
        let mut partial = BigRational::zero();
        for j in 1..n {
            let mut denominator = BigInt::one();
            for _ in 0..power {
                denominator *= BigInt::from(j);
            }
            partial += BigRational::new(BigInt::one(), denominator);
        }
        if order == 0 {
            let constant = graph.node(self.euler_gamma, &[]);
            let minus_one = graph.int(-1);
            let negated = graph.node(core::MUL, &[minus_one, constant]);
            let h = graph.num(Number::rat(partial));
            return Outcome::Pinned(graph.node(core::ADD, &[h, negated]));
        }
        let mut factorial = BigInt::one();
        for j in 2..=order {
            factorial *= BigInt::from(j);
        }
        let sign = if order % 2 == 1 { BigInt::one() } else { -BigInt::one() };
        let scale = BigRational::from_integer(sign * factorial);
        let leading = match zeta_even(&BigInt::from(power)) {
            | Some(c) => times_pi_power(graph, self.pi, &(&scale * c), Number::from(power)),
            | None => {
                let s = graph.int(power);
                let zeta = graph.node(self.zeta, &[s]);
                let scale_term = graph.num(Number::rat(scale.clone()));
                graph.node(core::MUL, &[scale_term, zeta])
            },
        };
        if partial.is_zero() {
            return Outcome::Pinned(leading);
        }
        let correction = graph.num(Number::rat(-(scale * partial)));
        Outcome::Pinned(graph.node(core::ADD, &[leading, correction]))
    }
}

/// `floor`, `ceil` and `round` (half away from zero) of a value: exact
/// for rational literals, and for other constant terms whenever the
/// numeric value is far enough from a rounding boundary to be certain.
struct Rounding {
    floor: OpId,
    ceil: OpId,
    round: OpId,
}

impl Kernel for Rounding {
    fn ops(&self) -> Vec<OpId> {
        vec![self.floor, self.ceil, self.round]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let Some(&arg) = graph.children(node).first() else {
            return Outcome::Pass;
        };
        let op = graph.op(node);
        if let Some(r) = graph.as_number(arg).and_then(Number::to_rational) {
            let rounded = if op == self.floor {
                r.floor()
            } else if op == self.ceil {
                r.ceil()
            } else {
                r.round()
            };
            return Outcome::Equal(graph.num(Number::rat(rounded)));
        }
        if !graph.free_symbols(graph.find(arg)).is_empty() {
            return Outcome::Pass;
        }
        let Some(value) = graph.eval(arg, &crate::graph::Env::numeric(0.0)) else {
            return Outcome::Pass;
        };
        if !value.is_finite() || value.abs() > 1e15 {
            return Outcome::Pass;
        }
        // The distance to the nearest boundary must dwarf the evaluation error.
        let boundary = if op == self.round { (value - 0.5).round() + 0.5 } else { value.round() };
        if (value - boundary).abs() <= 1e-9 * value.abs().max(1.0) {
            return Outcome::Pass;
        }
        let rounded = if op == self.floor {
            value.floor()
        } else if op == self.ceil {
            value.ceil()
        } else {
            value.round()
        };
        rounded
            .to_i64()
            .map_or(Outcome::Pass, |v| Outcome::Equal(graph.int(v)))
    }
}

struct Step {
    sign: OpId,
    heaviside: OpId,
}

impl Kernel for Step {
    fn ops(&self) -> Vec<OpId> {
        vec![self.sign, self.heaviside]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let Some(&arg) = graph.children(node).first() else {
            return Outcome::Pass;
        };
        let facts = graph.facts(arg);
        let (below, at, above) = if graph.op(node) == self.sign {
            (Number::from(-1), Number::from(0), Number::from(1))
        } else {
            (
                Number::from(0),
                Number::fraction(1, 2).unwrap_or_else(|| 0.into()),
                Number::from(1),
            )
        };
        let value = if facts.has(Facts::POSITIVE) {
            above
        } else if facts.has(Facts::NEGATIVE) {
            below
        } else if facts.has(Facts::NONNEGATIVE) && facts.has(Facts::NONPOSITIVE) {
            at
        } else {
            return Outcome::Pass;
        };
        Outcome::Equal(graph.num(value))
    }

    fn revisit(&self) -> bool {
        true
    }
}

/// An orthogonal polynomial family.
#[derive(Copy, Clone)]
enum Family {
    Legendre,
    ChebyshevT,
    ChebyshevU,
    Hermite,
    Laguerre,
}

/// Coefficients (ascending powers) of the degree-`n` member of `family`,
/// by its three-term recurrence
/// `p[k+1] = (a x + b) p[k] - c p[k-1]` in exact arithmetic.
fn coefficients(
    family: Family,
    n: u32,
) -> Vec<BigRational> {
    let q = |num: i64, den: i64| BigRational::new(num.into(), den.into());
    let mut previous = vec![q(1, 1)];
    let mut current = match family {
        | Family::Laguerre => vec![q(1, 1), q(-1, 1)],
        | Family::ChebyshevU | Family::Hermite => vec![q(0, 1), q(2, 1)],
        | Family::Legendre | Family::ChebyshevT => vec![q(0, 1), q(1, 1)],
    };
    if n == 0 {
        return previous;
    }
    for k in 1..n {
        let k = i64::from(k);
        let (a, b, c) = match family {
            | Family::Legendre => (q(2 * k + 1, k + 1), q(0, 1), q(k, k + 1)),
            | Family::ChebyshevT | Family::ChebyshevU => (q(2, 1), q(0, 1), q(1, 1)),
            | Family::Hermite => (q(2, 1), q(0, 1), q(2 * k, 1)),
            | Family::Laguerre => (q(-1, k + 1), q(2 * k + 1, k + 1), q(k, k + 1)),
        };
        let mut next = vec![q(0, 1); current.len() + 1];
        for (j, coef) in current.iter().enumerate() {
            if let Some(slot) = next.get_mut(j + 1) {
                *slot += &a * coef;
            }
            if let Some(slot) = next.get_mut(j) {
                *slot += &b * coef;
            }
        }
        for (slot, coef) in next.iter_mut().zip(&previous) {
            *slot -= &c * coef;
        }
        previous = current;
        current = next;
    }
    current
}

/// Expands orthogonal polynomials of literal degree.
struct Orthogonal {
    families: Vec<(OpId, Family)>,
}

impl Kernel for Orthogonal {
    fn ops(&self) -> Vec<OpId> {
        self.families.iter().map(|f| f.0).collect()
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let op = graph.op(node);
        let (Some(&(_, family)), &[degree, x]) = (
            self.families.iter().find(|f| f.0 == op),
            graph.children(node),
        ) else {
            return Outcome::Pass;
        };
        let Some(n) = integer(graph, degree)
            .and_then(|d| d.to_u32())
            .filter(|&d| d <= MAX_DEGREE)
        else {
            return Outcome::Pass;
        };
        let coefficients = coefficients(family, n);
        // At a literal point the value itself is the answer; a pinned
        // polynomial in a number would stay unevaluated.
        match graph.number_of(x).cloned() {
            | Some(Number::Float(v)) => {
                let value = coefficients.iter().rev().fold(0.0_f64, |acc, c| {
                    acc.mul_add(v, c.to_f64().unwrap_or(f64::NAN))
                });
                return if value.is_finite() {
                    Outcome::Equal(graph.float(value))
                } else {
                    Outcome::Pass
                };
            },
            | Some(exact) => {
                if let Some(point) = exact.to_rational() {
                    let value = coefficients
                        .iter()
                        .rev()
                        .fold(BigRational::zero(), |acc, c| acc * &point + c);
                    return Outcome::Equal(graph.num(Number::rat(value)));
                }
            },
            | None => {},
        }
        let mut terms = Vec::new();
        for (k, c) in coefficients.into_iter().enumerate() {
            if c.is_zero() {
                continue;
            }
            let power = match k {
                | 0 => None,
                | 1 => Some(x),
                | _ => {
                    let exponent = graph.int(i64::try_from(k).unwrap_or(0));
                    Some(graph.node(core::POW, &[x, exponent]))
                },
            };
            let coefficient = graph.num(Number::rat(c.clone()));
            terms.push(match power {
                | None => coefficient,
                | Some(p) if c.is_one() => p,
                | Some(p) => graph.node(core::MUL, &[coefficient, p]),
            });
        }
        // The expansion is the requested form; the compact spelling would
        // otherwise win extraction and hide the polynomial from the
        // differentiation and integration kernels.
        Outcome::Pinned(graph.node(core::ADD, &terms))
    }
}

// ----------------------------------------------------------------------
// Antiderivatives involving the error function
// ----------------------------------------------------------------------

fn depends(
    graph: &Graph,
    node: NodeId,
    x: crate::graph::SymbolId,
) -> bool {
    graph.depends_on(graph.find(node), x)
}

/// `∫ x^n exp(a x² + b x + c) dx` for a literal `n ≥ 0`: the Gaussian
/// `n = 0` through the error function, higher moments by the reduction
/// `I_n = (x^(n-1) e^q - (n-1) I_(n-2) - b I_(n-1)) / (2a)`.
fn gaussian_integral(
    cx: &mut KernelCx<'_>,
    f: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let graph = &mut *cx.graph;
    let symbol = graph.symbol_of(x)?;
    let (exp, erf, pi) = (graph.ops().lookup("exp")?, graph.ops().lookup("erf")?, graph.ops().lookup("pi")?);
    let factors = if graph.op(f) == core::MUL { graph.children(f).to_vec() } else { vec![f] };
    let mut n: u32 = 0;
    let mut exponent = None;
    for factor in factors {
        if graph.symbol_of(factor) == Some(symbol) {
            n = n.checked_add(1)?;
        } else if let (true, &[base, e]) = (graph.op(factor) == core::POW, graph.children(factor)) {
            if graph.symbol_of(base) != Some(symbol) {
                return None;
            }
            n = n.checked_add(integer(graph, e)?.to_u32()?)?;
        } else if graph.op(factor) == exp && exponent.is_none() {
            exponent = graph.children(factor).first().copied();
        } else {
            return None;
        }
    }
    let q = exponent?;
    if n > 12 {
        return None;
    }
    let mut gens = Gens::default();
    let gx = gens.index(graph, x);
    let poly = from_term(graph, &mut gens, q, Limits::default())?;
    for g in poly.support() {
        if g != gx && gens.node(g).is_some_and(|node| depends(graph, node, symbol)) {
            return None;
        }
    }
    if poly.degree_in(gx) != 2 {
        return None;
    }
    let coefficients = poly.coefficients_in(gx);
    let a = to_term(graph, &gens, coefficients.get(2)?);
    let b = to_term(graph, &gens, coefficients.get(1)?);
    let c = to_term(graph, &gens, coefficients.first()?);
    let int = |graph: &mut Graph, v: i64| graph.int(v);
    let mul = |graph: &mut Graph, f: &[NodeId]| graph.node(core::MUL, f);
    let add = |graph: &mut Graph, t: &[NodeId]| graph.node(core::ADD, t);
    let pow = |graph: &mut Graph, b: NodeId, e: NodeId| graph.node(core::POW, &[b, e]);
    let half = graph.num(Number::fraction(1, 2)?);
    let minus_one = int(graph, -1);
    // I0 = sqrt(pi) / (2 sqrt(-a)) * exp(c - b²/(4a)) * erf(sqrt(-a) (x + b/(2a))),
    // or with erfi and sqrt(a) when a is known positive.
    let growing = graph.facts(a).has(Facts::POSITIVE);
    let erf = if growing { graph.ops().lookup("erfi")? } else { erf };
    let minus_a = if growing { a } else { mul(graph, &[minus_one, a]) };
    let root = pow(graph, minus_a, half);
    let pi_node = graph.node(pi, &[]);
    let root_pi = pow(graph, pi_node, half);
    let two = int(graph, 2);
    let inv_two_root = {
        let d = mul(graph, &[two, root]);
        pow(graph, d, minus_one)
    };
    let b2 = pow(graph, b, two);
    let four_a = {
        let four = int(graph, 4);
        mul(graph, &[four, a])
    };
    let inv_four_a = pow(graph, four_a, minus_one);
    let shift_exp = {
        let t = mul(graph, &[minus_one, b2, inv_four_a]);
        let e = add(graph, &[c, t]);
        graph.node(exp, &[e])
    };
    let two_a = mul(graph, &[two, a]);
    let inv_two_a = pow(graph, two_a, minus_one);
    let shifted = {
        let t = mul(graph, &[b, inv_two_a]);
        let s = add(graph, &[x, t]);
        mul(graph, &[root, s])
    };
    let error = graph.node(erf, &[shifted]);
    let i0 = mul(graph, &[root_pi, inv_two_root, shift_exp, error]);
    let e_q = graph.node(exp, &[q]);
    let mut moments = vec![i0];
    for k in 1..=n {
        // I_k = (x^(k-1) e^q - (k-1) I_(k-2) - b I_(k-1)) / (2a)
        let k1 = int(graph, i64::from(k) - 1);
        let xk = pow(graph, x, k1);
        let mut terms = vec![mul(graph, &[xk, e_q])];
        if k >= 2 {
            let coefficient = int(graph, -(i64::from(k) - 1));
            terms.push(mul(graph, &[coefficient, moments[(k - 2) as usize]]));
        }
        terms.push(mul(graph, &[minus_one, b, moments[(k - 1) as usize]]));
        let sum = add(graph, &terms);
        moments.push(mul(graph, &[sum, inv_two_a]));
    }
    moments.pop()
}

/// `∫ erf(a x + b) dx = ((a x + b) erf(a x + b) + exp(-(a x + b)²)/√π) / a`,
/// and the same for `erfc` with the sign of the second term flipped.
fn erf_integral(
    cx: &mut KernelCx<'_>,
    f: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let graph = &mut *cx.graph;
    let symbol = graph.symbol_of(x)?;
    let (exp, erf, erfc, pi) =
        (graph.ops().lookup("exp")?, graph.ops().lookup("erf")?, graph.ops().lookup("erfc")?, graph.ops().lookup("pi")?);
    let op = graph.op(f);
    if op != erf && op != erfc {
        return None;
    }
    let &[u] = graph.children(f) else {
        return None;
    };
    let mut gens = Gens::default();
    let gx = gens.index(graph, x);
    let poly = from_term(graph, &mut gens, u, Limits::default())?;
    for g in poly.support() {
        if g != gx && gens.node(g).is_some_and(|node| depends(graph, node, symbol)) {
            return None;
        }
    }
    if poly.degree_in(gx) != 1 {
        return None;
    }
    let a = to_term(graph, &gens, poly.coefficients_in(gx).get(1)?);
    let (minus_one, two) = (graph.int(-1), graph.int(2));
    let half = graph.num(Number::fraction(-1, 2)?);
    let square = graph.node(core::POW, &[u, two]);
    let negated = graph.node(core::MUL, &[minus_one, square]);
    let gauss = graph.node(exp, &[negated]);
    let pi_node = graph.node(pi, &[]);
    let inv_root_pi = graph.node(core::POW, &[pi_node, half]);
    let sign = graph.int(if op == erf { 1 } else { -1 });
    let tail = graph.node(core::MUL, &[sign, gauss, inv_root_pi]);
    let head = graph.node(core::MUL, &[u, f]);
    let sum = graph.node(core::ADD, &[head, tail]);
    let inv_a = graph.node(core::POW, &[a, minus_one]);
    Some(graph.node(core::MUL, &[sum, inv_a]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Budget;
    use crate::graph::ClosedForm;
    use crate::graph::Engine;
    use crate::graph::Env;
    use crate::graph::Extractor;
    use crate::graph::Saturate;
    use crate::rules::testing::eval;
    use crate::rules::testing::numeric;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    #[test]
    fn polygamma_at_integers() {
        let sets = [special()];
        assert_eq!(simplify(&sets, "digamma(1)"), "-euler_gamma");
        assert_eq!(simplify(&sets, "digamma(4)"), "11/6 - euler_gamma");
        assert_eq!(simplify(&sets, "polygamma(1, 1)"), "1/6*pi^2");
        assert_eq!(simplify(&sets, "polygamma(1, 3)"), "1/6*pi^2 - 5/4");
        let value = |src: &str| numeric(&sets, src, &[], 1e-12).0;
        assert!((value("polygamma(2, 2)") + 2.0 * (1.202_056_903_159_594 - 1.0)).abs() < 1e-9);
        assert!((value("digamma(1/2)") + 1.963_510_026_021_423).abs() < 1e-9);
    }

    #[test]
    fn rounding() {
        let sets = [special()];
        for (src, want) in [
            ("floor(7/2)", "3"),
            ("ceil(7/2)", "4"),
            ("round(7/2)", "4"),
            ("round(-7/2)", "-4"),
            ("floor(-7/2)", "-4"),
            ("floor(pi)", "3"),
            ("ceil(exp(1))", "3"),
            ("round(2^(1/2)*10)", "14"),
            ("floor(x + 3)", "floor(x) + 3"),
            ("ceil(-x)", "-floor(x)"),
            ("frac(5/2)", "1/2"),
        ] {
            assert_eq!(simplify(&sets, src), want, "{src}");
        }
        let (text, reduced) = reduce_with(&[crate::rules::standard()].concat(), "floor(n)", &[("n", Facts::INTEGER)]);
        assert!(reduced && text == "n", "{text}");
    }

    #[test]
    fn bessel_functions_of_any_order() {
        let sets = [special()];
        let value = |src: &str| numeric(&sets, src, &[], 1e-12).0;
        assert!((value("besselj(2, 1)") - 0.114_903_484_9).abs() < 1e-9);
        assert!((value("besselk(0, 1)") - 0.421_024_438_2).abs() < 1e-9);
        assert!((value("bessel_k1(2)") - 0.139_865_881_8).abs() < 1e-9);
        // d/dx J_0 = -J_1 and the defining equation holds numerically.
        assert_eq!(simplify(&sets, "diff(besselj(0, x), x)"), "besselj(-1, x)");
        assert_eq!(simplify(&sets, "besselj(-1, x) + besselj(1, x)"), "0");
        let residual = simplify(&sets, "x^2*diff(diff(besselj(3, x), x), x) + x*diff(besselj(3, x), x) + (x^2 - 9)*besselj(3, x)");
        assert!(eval(&sets, &residual, &[("x", 2.7)]).abs() < 1e-10, "{residual}");
        assert!(simplify(&sets, "bessel_differential_equation(y(x), x, n)").ends_with("= 0"));
    }

    #[test]
    fn polygamma_erfi_and_relatives() {
        let sets = [special()];
        assert_eq!(simplify(&sets, "polygamma(0, x)"), "digamma(x)");
        assert_eq!(simplify(&sets, "diff(digamma(x), x)"), "polygamma(1, x)");
        let value = |src: &str| numeric(&sets, src, &[], 1e-12).0;
        // ψ'(1) = π²/6
        assert!((value("polygamma(1, 1)") - std::f64::consts::PI.powi(2) / 6.0).abs() < 1e-10);
        assert!((value("erfi(1)") - 1.650_425_758_797_542_8).abs() < 1e-12);
        assert!((value("erfc(inverse_erfc(0.3))") - 0.3).abs() < 1e-12);
        assert!((value("ln_factorial(5)") - 120.0_f64.ln()).abs() < 1e-12);
        assert_eq!(simplify(&sets, "integral(exp(x^2), x)"), "1/2*erfi(x)*pi^(1/2)");
    }

    #[test]
    fn rodrigues_formulas_and_generalized_laguerre() {
        let sets = [special()];
        assert_eq!(simplify(&sets, "legendre_rodrigues(3, x)"), simplify(&sets, "legendre(3, x)"));
        assert_eq!(simplify(&sets, "hermite_rodrigues(3, x)"), simplify(&sets, "hermite(3, x)"));
        assert_eq!(simplify(&sets, "laguerre_rodrigues(2, x)"), simplify(&sets, "laguerre(2, x)"));
        assert_eq!(simplify(&sets, "laguerre_gen(2, 1, x)"), "1/2*x^2 - 3*x + 3");
        assert_eq!(simplify(&sets, "laguerre_gen(1, a, x)"), simplify(&sets, "1 + a - x"));
        let (v, _) = numeric(&sets, "laguerre_gen(3, 0.5, 1.2)", &[], 1e-12);
        assert!((v - crate::kernels::special::generalized_laguerre(3, 0.5, 1.2)).abs() < 1e-12);
    }

    #[test]
    fn gaussian_integrals() {
        let sets = [special()];
        let run = |src: &str| simplify(&sets, src);
        assert_eq!(run("integral(exp(-x^2), x)"), "1/2*erf(x)*pi^(1/2)");
        assert_eq!(run("defint(exp(-x^2), x, -oo, oo)"), "pi^(1/2)");
        assert_eq!(run("defint(x^2*exp(-x^2), x, -oo, oo)"), "1/2*pi^(1/2)");
        assert_eq!(run("defint(x^4*exp(-x^2), x, -oo, oo)"), "3/4*pi^(1/2)");
        assert_eq!(run("integral(erf(2*x), x)"), "x*erf(2*x) + 1/2*exp(-4*x^2)/pi^(1/2)");
        assert_eq!(run("limit(erf(x), x, -oo)"), "-1");
        // A shifted, scaled Gaussian against its numeric value.
        let closed = run("defint(exp(-x^2/2 + x), x, -oo, oo)");
        let (value, _) = numeric(&sets, "defint(exp(-x^2/2 + x), x, -oo, oo)", &[], 1e-10);
        assert!((eval(&sets, &closed, &[]) - value).abs() < 1e-8, "{closed}");
    }

    fn s(src: &str) -> String {
        simplify(&[special()], src)
    }

    fn ev(
        src: &str,
        bindings: &[(&str, f64)],
    ) -> f64 {
        eval(&[special()], src, bindings)
    }

    #[track_caller]
    fn close(
        got: f64,
        want: f64,
        relative: f64,
    ) {
        let scale = want.abs().max(1e-300);
        assert!(
            (got - want).abs() <= relative * scale,
            "got {got}, want {want} (relative tolerance {relative})"
        );
    }

    /// Value of a text like `-7/12` or `24`.
    fn rational(text: &str) -> f64 {
        let value = |t: &str| t.trim().parse::<f64>().unwrap_or(f64::NAN);
        text.split_once('/')
            .map_or_else(|| value(text), |(n, d)| value(n) / value(d))
    }

    /// Agreement to `1e-12`, relative for large values and absolute near 0.
    #[track_caller]
    fn agree(
        a: f64,
        b: f64,
        tolerance: f64,
    ) {
        assert!((a - b).abs() <= tolerance * (1.0 + b.abs()), "{a} vs {b}");
    }

    #[test]
    fn gamma_family_values() {
        // Relative tolerance 1e-12: statrs' Lanczos gamma is good to ~1e-15.
        close(ev("gamma(x)", &[("x", 0.5)]), 1.772_453_850_905_516, 1e-12);
        close(
            ev("gamma(x)", &[("x", 10.5)]),
            1_133_278.388_948_785_5,
            1e-12,
        );
        close(ev("gamma(x)", &[("x", 5.0)]), 24.0, 1e-12);
        close(ev("gamma(x)", &[("x", -1.5)]), 2.363_271_801_207_355, 1e-12);
        close(
            ev("lgamma(x)", &[("x", 10.0)]),
            12.801_827_480_081_469,
            1e-12,
        );
        close(
            ev("lgamma(x)", &[("x", 0.5)]),
            0.572_364_942_924_700_1,
            1e-12,
        );
        close(
            ev("digamma(x)", &[("x", 1.0)]),
            -0.577_215_664_901_532_9,
            1e-10,
        );
        close(
            ev("digamma(x)", &[("x", 0.5)]),
            -1.963_510_026_021_423_5,
            1e-10,
        );
        close(
            ev("digamma(x)", &[("x", 10.0)]),
            2.251_752_589_066_721,
            1e-10,
        );
        close(
            ev("beta(a, b)", &[("a", 2.0), ("b", 3.0)]),
            1.0 / 12.0,
            1e-12,
        );
        close(
            ev("beta(a, b)", &[("a", 0.5), ("b", 0.5)]),
            std::f64::consts::PI,
            1e-12,
        );
        // Outside the positive quadrant through the gamma quotient.
        close(ev("beta(a, b)", &[("a", -0.5), ("b", 2.0)]), -4.0, 1e-12);
    }

    #[test]
    fn error_function_values() {
        // statrs' erf is good to about 1e-11 (erf(1) is off by 7e-12).
        close(ev("erf(x)", &[("x", 1.0)]), 0.842_700_792_949_714_9, 1e-10);
        close(ev("erf(x)", &[("x", 0.5)]), 0.520_499_877_813_046_5, 1e-10);
        close(
            ev("erf(x)", &[("x", -2.0)]),
            -0.995_322_265_018_952_7,
            1e-10,
        );
        close(
            ev("erfc(x)", &[("x", 1.0)]),
            0.157_299_207_050_285_13,
            1e-10,
        );
        close(ev("erfc(x)", &[("x", 3.0)]), 2.209_049_699_858_544e-5, 1e-9);
    }

    #[test]
    fn zeta_values() {
        // Euler-Maclaurin with 12 direct terms and 7 corrections.
        close(ev("zeta(s)", &[("s", 2.0)]), 1.644_934_066_848_226_4, 1e-14);
        close(ev("zeta(s)", &[("s", 3.0)]), 1.202_056_903_159_594_2, 1e-14);
        close(ev("zeta(s)", &[("s", 4.0)]), 1.082_323_233_711_138_2, 1e-14);
        close(ev("zeta(s)", &[("s", 1.5)]), 2.612_375_348_685_488, 1e-14);
        close(ev("zeta(s)", &[("s", 1.1)]), 10.584_448_464_950_809, 1e-13);
        close(
            ev("zeta(s)", &[("s", 10.0)]),
            1.000_994_575_127_818_1,
            1e-14,
        );
        close(ev("zeta(s)", &[("s", 50.0)]), 1.0, 1e-14);
        // The analytic continuation: ζ(1/2) and ζ(-1) = -1/12.
        close(ev("zeta(s)", &[("s", 0.5)]), -1.460_354_508_809_587, 1e-10);
        close(ev("zeta(s)", &[("s", -1.0)]), -1.0 / 12.0, 1e-10);
        assert!(ev("zeta(s)", &[("s", 1.0)]).is_infinite());
    }

    #[test]
    fn bessel_values() {
        // The available implementations are polynomial fits (Numerical
        // Recipes) good to about 1e-7, only for orders 0 and 1 (order 1 of `J` and `Y` only from |x| = 8): tolerance
        // 5e-7 relative to `1 + |value|`. Other orders are NaN.
        let at = |name: &str, n: f64, x: f64| ev(&format!("{name}(n, x)"), &[("n", n), ("x", x)]);
        agree(at("besselj", 0.0, 1.0), 0.765_197_686_557_966_6, 5e-7);
        agree(at("besselj", 1.0, 10.0), 0.043_472_746_168_861_44, 5e-7);
        agree(at("besselj", 1.0, 20.0), 0.066_833_124_175_850_05, 5e-7);
        agree(at("besselj", 0.0, 5.0), -0.177_596_771_314_338_3, 5e-7);
        agree(at("besselj", 0.0, 20.0), 0.167_024_664_340_583_1, 5e-7);
        agree(at("bessely", 0.0, 1.0), 0.088_256_964_215_676_97, 5e-7);
        agree(at("bessely", 1.0, 10.0), 0.249_015_424_206_953_9, 5e-7);
        agree(at("bessely", 0.0, 10.0), 0.055_671_167_283_599_39, 5e-7);
        agree(at("besseli", 0.0, 1.0), 1.266_065_877_752_008_4, 5e-7);
        agree(at("besseli", 1.0, 1.0), 0.565_159_103_992_485_1, 5e-7);
        agree(at("besseli", 0.0, 5.0) / 27.239_871_823_604_442, 1.0, 5e-7);
        agree(at("besseli", 1.0, 5.0), 24.335_642_142_450_52, 5e-7);
        // Any order, small arguments included.
        agree(at("besselj", 1.0, 1.0), 0.440_050_585_744_933_5, 1e-12);
        agree(at("bessely", 1.0, 1.0), -0.781_212_821_300_288_7, 1e-9);
        agree(at("besselj", 2.0, 1.0), 0.114_903_484_931_900_5, 1e-9);
        let half = (2.0 / std::f64::consts::PI).sqrt() * 1.0_f64.sin();
        agree(at("besselj", 0.5, 1.0), half, 1e-12);
    }

    #[test]
    fn step_like_values() {
        close(
            ev("sinc(x)", &[("x", 0.5)]),
            2.0 / std::f64::consts::PI,
            1e-14,
        );
        assert!((ev("sinc(x)", &[("x", 0.0)]) - 1.0).abs() < 1e-9);
        assert!(ev("sinc(x)", &[("x", 1.0)]).abs() < 1e-15);
        for (x, h, sg) in [(-2.0, 0.0, -1.0), (0.0, 0.5, 0.0), (3.0, 1.0, 1.0)] {
            assert!((ev("heaviside(x)", &[("x", x)]) - h).abs() < 1e-9);
            assert!((ev("sign(x)", &[("x", x)]) - sg).abs() < 1e-9);
        }
    }

    #[test]
    fn exact_gamma_values() {
        assert_eq!(s("gamma(1)"), "1");
        assert_eq!(s("gamma(5)"), "24");
        assert_eq!(s("gamma(20)"), "121645100408832000");
        assert_eq!(s("gamma(1/2)"), "pi^(1/2)");
        assert_eq!(s("gamma(3/2)"), "1/2*pi^(1/2)");
        assert_eq!(s("gamma(5/2)"), "3/4*pi^(1/2)");
        assert_eq!(s("gamma(-1/2)"), "-2*pi^(1/2)");
        assert_eq!(s("gamma(-3/2)"), "4/3*pi^(1/2)");
        // Poles and other rationals stay.
        assert_eq!(s("gamma(0)"), "gamma(0)");
        assert_eq!(s("gamma(-2)"), "gamma(-2)");
        assert_eq!(s("gamma(1/3)"), "gamma(1/3)");
        assert_eq!(s("gamma(x)"), "gamma(x)");
        // The exact value agrees with the numeric one.
        close(ev("gamma(7/2)", &[]), 3.323_350_970_447_843, 1e-12);
        close(ev("15/8*pi^(1/2)", &[]), 3.323_350_970_447_843, 1e-12);
    }

    #[test]
    fn exact_zeta_values() {
        assert_eq!(s("zeta(2)"), "1/6*pi^2");
        assert_eq!(s("zeta(4)"), "1/90*pi^4");
        assert_eq!(s("zeta(6)"), "1/945*pi^6");
        assert_eq!(s("zeta(8)"), "1/9450*pi^8");
        assert_eq!(s("zeta(10)"), "1/93555*pi^10");
        assert_eq!(s("zeta(12)"), "691/638512875*pi^12");
        for k in [2.0, 4.0, 6.0, 8.0, 10.0, 12.0] {
            let text = s(&format!("zeta({k})"));
            close(ev(&text, &[]), ev("zeta(s)", &[("s", k)]), 1e-13);
        }
        assert_eq!(s("zeta(3)"), "zeta(3)");
        assert_eq!(s("zeta(14)"), "zeta(14)");
    }

    #[test]
    fn rewrites_fire() {
        assert_eq!(s("gamma(x + 1) / gamma(x)"), "x");
        assert_eq!(s("2 * gamma(y + 1) / gamma(y)"), "2*y");
        assert_eq!(s("lgamma(1)"), "0");
        assert_eq!(s("lgamma(2)"), "0");
        assert_eq!(s("erf(0)"), "0");
        assert_eq!(s("erfc(0)"), "1");
        assert_eq!(s("erf(-x)"), "-erf(x)");
        assert_eq!(s("erfc(x) + erf(x)"), "1");
        assert_eq!(s("erfc(x) + 3 + erf(x)"), "4");
        assert_eq!(s("sinc(0)"), "1");
        assert_eq!(s("sinc(-x)"), "sinc(x)");
        assert_eq!(s("sign(-x)"), "-sign(x)");
        // Exploring tier: the beta function through gamma values.
        assert_eq!(s("beta(2, 3)"), "1/12");
        assert_eq!(s("beta(1/2, 1/2)"), "pi");
        assert_eq!(s("exp(lgamma(5))"), "24");
        // `lgamma` needs a positive argument to become ln(gamma).
        assert_eq!(
            reduce_with(&[special()], "exp(lgamma(x))", &[]).0,
            "exp(lgamma(x))"
        );
        assert_eq!(
            reduce_with(&[special()], "exp(lgamma(x))", &[("x", Facts::POSITIVE)]).0,
            "gamma(x)"
        );
    }

    #[test]
    fn signs_from_facts() {
        let run = |src: &str, facts: Facts| reduce_with(&[special()], src, &[("x", facts)]).0;
        assert_eq!(run("sign(x)", Facts::POSITIVE), "1");
        assert_eq!(run("sign(x)", Facts::NEGATIVE), "-1");
        assert_eq!(run("sign(x)", Facts::REAL), "sign(x)");
        assert_eq!(run("sign(x^2 + 1)", Facts::REAL), "1");
        assert_eq!(run("heaviside(x)", Facts::POSITIVE), "1");
        assert_eq!(run("heaviside(x)", Facts::NEGATIVE), "0");
        assert_eq!(s("sign(0)"), "0");
        assert_eq!(s("sign(-7/3)"), "-1");
        assert_eq!(s("heaviside(0)"), "1/2");
        assert_eq!(s("heaviside(5)"), "1");
    }

    /// Central finite difference of the operator expression `f`.
    fn finite_difference(
        f: &str,
        x: f64,
    ) -> f64 {
        let h = 1e-5;
        (ev(f, &[("x", x + h)]) - ev(f, &[("x", x - h)])) / (2.0 * h)
    }

    /// The closed form `src` reduces to, evaluated at `x` directly on the
    /// term (its printed form does not always parse back the same).
    fn symbolic_value(
        src: &str,
        x: f64,
    ) -> f64 {
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &[special()]).unwrap_or_else(|e| panic!("{e}"));
        let root = g.parse(src).unwrap_or_else(|e| panic!("{src}: {e}"));
        engine.run(
            &mut g,
            &[root],
            &Env::symbolic(),
            &Saturate,
            &Budget::default(),
        );
        let node = Extractor::new(&g, &[root], &ClosedForm).build(&mut g, root);
        let mut env = Env::numeric(0.0);
        env.bind(g.interner_mut().symbol("x"), x);
        node.and_then(|n| g.eval(n, &env)).unwrap_or(f64::NAN)
    }

    #[test]
    fn derivatives_in_both_phases() {
        let x = 0.3;
        for f in [
            "erf(x)",
            "erfc(x)",
            "gamma(x)",
            "lgamma(x)",
            "beta(x, 2)",
            "beta(3/2, x)",
            "sinc(x)",
            "erf(2*x^2)",
        ] {
            let want = finite_difference(f, x);
            // Numeric phase.
            let (numeric_value, _) =
                numeric(&[special()], &format!("diff({f}, x)"), &[("x", x)], 1e-9);
            assert!(
                (numeric_value - want).abs() < 1e-6,
                "{f}: numeric {numeric_value} vs {want}"
            );
            // Symbolic phase: the reduced derivative, evaluated.
            let symbolic = s(&format!("diff({f}, x)"));
            assert!(!symbolic.contains("diff"), "{f}: {symbolic}");
            let value = symbolic_value(&format!("diff({f}, x)"), x);
            assert!(
                (value - want).abs() < 1e-6,
                "{f}: {symbolic} = {value} vs {want}"
            );
        }
        assert_eq!(s("diff(lgamma(x), x)"), "digamma(x)");
        assert_eq!(s("diff(gamma(x), x)"), "digamma(x)*gamma(x)");
    }

    #[test]
    fn numeric_phase_uses_the_evaluations() {
        let at = |src: &str, x: f64| numeric(&[special()], src, &[("x", x)], 1e-12).0;
        close(at("gamma(x)", 4.5), 11.631_728_396_567_448, 1e-12);
        close(at("erf(x) + erfc(x)", 0.7), 1.0, 1e-12);
        close(at("zeta(x)", 2.0), 1.644_934_066_848_226_4, 1e-13);
    }

    const CLOSED_FORMS: [(&str, [&str; 6]); 5] = [
        (
            "legendre",
            [
                "1",
                "x",
                "(3*x^2 - 1)/2",
                "(5*x^3 - 3*x)/2",
                "(35*x^4 - 30*x^2 + 3)/8",
                "(63*x^5 - 70*x^3 + 15*x)/8",
            ],
        ),
        (
            "chebyshevt",
            [
                "1",
                "x",
                "2*x^2 - 1",
                "4*x^3 - 3*x",
                "8*x^4 - 8*x^2 + 1",
                "16*x^5 - 20*x^3 + 5*x",
            ],
        ),
        (
            "chebyshevu",
            [
                "1",
                "2*x",
                "4*x^2 - 1",
                "8*x^3 - 4*x",
                "16*x^4 - 12*x^2 + 1",
                "32*x^5 - 32*x^3 + 6*x",
            ],
        ),
        (
            "hermite",
            [
                "1",
                "2*x",
                "4*x^2 - 2",
                "8*x^3 - 12*x",
                "16*x^4 - 48*x^2 + 12",
                "32*x^5 - 160*x^3 + 120*x",
            ],
        ),
        (
            "laguerre",
            [
                "1",
                "1 - x",
                "(x^2 - 4*x + 2)/2",
                "(-x^3 + 9*x^2 - 18*x + 6)/6",
                "(x^4 - 16*x^3 + 72*x^2 - 96*x + 24)/24",
                "(-x^5 + 25*x^4 - 200*x^3 + 600*x^2 - 600*x + 120)/120",
            ],
        ),
    ];

    #[test]
    fn low_degree_polynomials_match_closed_forms() {
        for (name, forms) in CLOSED_FORMS {
            for (n, closed) in forms.iter().enumerate() {
                let expanded = s(&format!("{name}({n}, x)"));
                assert!(!expanded.contains(name), "{name}({n}) stayed: {expanded}");
                for x in [-1.3, -0.4, 0.0, 0.25, 0.9, 2.1] {
                    agree(ev(&expanded, &[("x", x)]), ev(closed, &[("x", x)]), 1e-12);
                }
            }
        }
        assert_eq!(s("legendre(3, x)"), "5/2*x^3 - 3/2*x");
        assert_eq!(s("chebyshevt(3, x)"), "4*x^3 - 3*x");
        assert_eq!(s("hermite(4, x)"), "16*x^4 - 48*x^2 + 12");
        assert_eq!(s("legendre(0, y)"), "1");
        assert_eq!(s("legendre(1, y)"), "y");
    }

    #[test]
    fn polynomials_of_any_degree_match_the_recurrence() {
        for name in [
            "legendre",
            "chebyshevt",
            "chebyshevu",
            "hermite",
            "laguerre",
        ] {
            for n in 0..=20 {
                let expanded = s(&format!("{name}({n}, x)"));
                let recurrence = format!("{name}({n}, x)");
                for x in [-0.9, -0.3, 0.2, 0.7] {
                    agree(
                        ev(&expanded, &[("x", x)]),
                        ev(&recurrence, &[("x", x)]),
                        1e-9,
                    );
                }
                // At an exact point the result is an exact rational.
                for point in ["1/3", "-2/5"] {
                    let exact = s(&format!("{name}({n}, {point})"));
                    let want = ev(&recurrence, &[("x", rational(point))]);
                    agree(rational(&exact), want, 1e-11);
                }
            }
        }
        // Degrees beyond the expansion limit and symbolic degrees stay.
        assert_eq!(s("legendre(65, x)"), "legendre(65, x)");
        assert_eq!(s("hermite(n, x)"), "hermite(n, x)");
        assert_eq!(s("legendre(-1, x)"), "legendre(-1, x)");
        assert_eq!(s("legendre(1/2, x)"), "legendre(1/2, x)");
        // Recurrence-defined values at special points.
        assert_eq!(s("legendre(10, 1)"), "1");
        assert_eq!(s("chebyshevt(7, 1)"), "1");
        assert_eq!(s("hermite(3, 1)"), "-4");
    }

    #[test]
    fn polynomial_arguments_may_be_expressions() {
        assert_eq!(s("chebyshevt(2, y + 1)"), "2*(y + 1)^2 - 1");
        assert_eq!(s("diff(legendre(3, x), x)"), "15/2*x^2 - 3/2");
        assert_eq!(s("diff(chebyshevt(3, 2*x), x)"), "96*x^2 - 6");
        // The float value at a float point.
        agree(ev("legendre(3, 0.5)", &[]), -0.437_5, 1e-14);
        assert_eq!(s("legendre(2, 0.5)"), format!("{}", -0.125));
    }

    /// The residual of a defining equation `lhs = 0` evaluated at `x`, for the
    /// solution `y` of the family and its order `n`.
    fn ode_residual(
        equation: &str,
        y: &str,
        n: u32,
        x: f64,
    ) -> f64 {
        let sets = [special()];
        let text = simplify(&sets, &format!("{equation}({y}, x, {n})"));
        let lhs = text.strip_suffix(" = 0").unwrap_or_else(|| panic!("not an equation `= 0`: {text}"));
        eval(&sets, lhs, &[("x", x)])
    }

    #[test]
    fn classical_orthogonal_polynomials_solve_their_equations() {
        for n in 0..=5 {
            for x in [-0.8, -0.3, 0.2, 0.9] {
                let tol = 1e-8;
                let r = ode_residual("legendre_differential_equation", &format!("legendre({n}, x)"), n, x);
                assert!(r.abs() < tol, "Legendre {n} at {x}: {r}");
                let r = ode_residual("hermite_differential_equation", &format!("hermite({n}, x)"), n, x);
                assert!(r.abs() < tol, "Hermite {n} at {x}: {r}");
                let r = ode_residual("chebyshev_differential_equation", &format!("chebyshevt({n}, x)"), n, x);
                assert!(r.abs() < tol, "Chebyshev T {n} at {x}: {r}");
                let r = ode_residual("laguerre_differential_equation", &format!("laguerre({n}, x)"), n, x + 1.0);
                assert!(r.abs() < tol, "Laguerre {n} at {x}: {r}");
            }
        }
    }

    #[test]
    fn a_wrong_degree_does_not_solve_the_equation() {
        // P_3 solves the degree-3 equation but not the degree-2 one, and likewise for the others.
        assert!(ode_residual("legendre_differential_equation", "legendre(3, x)", 2, 0.4).abs() > 1e-3);
        assert!(ode_residual("hermite_differential_equation", "hermite(3, x)", 2, 0.4).abs() > 1e-3);
        assert!(ode_residual("chebyshev_differential_equation", "chebyshevt(3, x)", 2, 0.4).abs() > 1e-3);
        assert!(ode_residual("laguerre_differential_equation", "laguerre(3, x)", 2, 0.4).abs() > 1e-3);
        // U_3 is not a solution of Chebyshev's first-kind equation.
        assert!(ode_residual("chebyshev_differential_equation", "chebyshevu(3, x)", 3, 0.4).abs() > 1e-3);
    }
}
