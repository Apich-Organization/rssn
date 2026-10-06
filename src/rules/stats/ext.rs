//! More distributions and probability calculus: further families, derived
//! properties of distributions, expectations of transformations, sums and
//! affine maps of independent variables, conjugate Bayesian updates and
//! exact small-sample tests.
//!
//! # Families
//!
//! `cauchy(x0, gamma)`, `laplace_dist(mu, b)`, `lognormal(mu, sigma)`,
//! `weibull(k, lambda)` (shape, scale), `geometric(p)` (trials up to the
//! first success, support `1, 2, ...`), `negative_binomial(r, p)`
//! (failures before the `r`-th success), `f_dist(d1, d2)`,
//! `pareto(xm, alpha)`, `rayleigh(sigma)`, `logistic(mu, s)` and
//! `gumbel(mu, beta)` are inert terms like the families of the parent
//! module and answer the same requests `pdf`, `cdf`, `expectation`,
//! `variance_of`, `mgf` and `entropy_of` (a request without a closed form,
//! such as the mean of the Cauchy distribution, stays unreduced; moments
//! that do not exist for the given numeric parameters likewise).
//!
//! # Properties of a distribution `d`
//!
//! | operator | value |
//! |---|---|
//! | `charfn(d, t)` | the characteristic function `E[exp(i t X)]` |
//! | `skewness_of(d)`, `kurtosis_of(d)` | standardised third moment; *excess* kurtosis |
//! | `median_of(d)`, `mode_of(d)` | where the closed form exists |
//! | `quantile_of(d, p)` | the inverse distribution function (through `erfcinv`, `ln`, `tan`, ...) |
//! | `std_of(d)` | the square root of the variance |
//! | `moment_of(d, k)`, `central_moment_of(d, k)` | raw and central moments, `k <= 8`, from the derivatives of the moment generating function (the uniform and heavy-tailed families by their own formulas where needed) |
//! | `expect_of(g, x, d)` | `E[g(X)]` as the integral `∫ g(x) f(x) dx` over the support (a sum for discrete families), left to the integrator and the summation rules |
//!
//! # Independent variables
//!
//! | operator | value |
//! |---|---|
//! | `sum_dist(d1, d2)` | the distribution of `X1 + X2` for independent variables, when the family is closed under convolution (normal, Cauchy, Poisson, chi-squared, binomial and Bernoulli with a common probability, gamma and exponential with a common rate, negative binomial with a common probability) |
//! | `iid_sum(d, n)`, `iid_mean(d, n)` | the sum and the mean of `n` independent copies |
//! | `scale_dist(c, d)`, `shift_dist(s, d)` | the distribution of `c X` (`c > 0` for the families that need it) and of `X + s` |
//!
//! # Bayesian updates and exact tests
//!
//! | operator | value |
//! |---|---|
//! | `bayes_update(prior, model, data)` | the conjugate posterior: `beta_dist` prior with the model `bernoulli`, `binomial_trials(m)` or `geometric`; `gamma_dist` prior with `poisson` or `exponential`; `normal(m0, s0)` prior with `normal_known_sd(s)`; `data` is the list of observations |
//! | `binomial_test(k, n, p0)` | the exact two-sided p-value (rational) |
//! | `fisher_exact(a, b, c, d)` | the exact two-sided p-value of the 2×2 table `[[a, b], [c, d]]` (rational) |
//! | `wilson_interval(k, n, z)` | `list(lower, upper)`, the Wilson score interval, in closed form |
//! | `mean_interval(mean, s, n, z)` | `list(mean - z s/√n, mean + z s/√n)` |
//! | `proportion_z_test(k, n, p0)` | `list(z, p_value)` with the p-value through the normal `cdf` |
//! | `r_squared(xs, ys)`, `residuals(xs, ys)` | the coefficient of determination of the simple regression and its residuals (exact for exact data) |

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;

use super::call;
use super::difference;
use super::f64_of;
use super::items_of;
use super::list_node;
use super::negate;
use super::sum_node;
use super::Family;
use super::FAMILIES as CORE_FAMILIES;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::Pat;
use crate::graph::RuleError;
use crate::graph::Tier;
use crate::graph::VarNames;
use crate::rules::calculus::derivative;
use crate::rules::combinatorics::choose;
use crate::rules::poly::best;

/// Families added to the ones of the parent module. An empty template
/// means that the quantity does not exist.
pub(super) const FAMILIES: [Family; 11] = [
    Family {
        name: "cauchy",
        params: 2,
        pdf: "1 / (pi * ?b * (1 + ((?x - ?a) / ?b)^2))",
        cdf: "1/2 + atan((?x - ?a) / ?b) / pi",
        mean: "",
        variance: "",
        mgf: None,
        entropy: Some("ln(4 * pi * ?b)"),
    },
    Family {
        name: "laplace_dist",
        params: 2,
        pdf: "exp(-abs(?x - ?a) / ?b) / (2 * ?b)",
        cdf: "1/2 + sign(?x - ?a) / 2 * (1 - exp(-abs(?x - ?a) / ?b))",
        mean: "?a",
        variance: "2 * ?b^2",
        mgf: Some("exp(?a * ?x) / (1 - ?b^2 * ?x^2)"),
        entropy: Some("1 + ln(2 * ?b)"),
    },
    Family {
        name: "lognormal",
        params: 2,
        pdf: "exp(-(ln(abs(?x)) - ?a)^2 / (2 * ?b^2)) / (abs(?x) * ?b * (2 * pi)^(1/2)) * heaviside(?x)",
        cdf: "1/2 * (1 + erf((ln(abs(?x)) - ?a) / (?b * 2^(1/2)))) * heaviside(?x)",
        mean: "exp(?a + ?b^2 / 2)",
        variance: "(exp(?b^2) - 1) * exp(2 * ?a + ?b^2)",
        mgf: None,
        entropy: Some("?a + 1/2 + ln(?b * (2 * pi)^(1/2))"),
    },
    Family {
        name: "weibull",
        params: 2,
        pdf: "?a / ?b * (abs(?x) / ?b)^(?a - 1) * exp(-(abs(?x) / ?b)^?a) * heaviside(?x)",
        cdf: "(1 - exp(-(abs(?x) / ?b)^?a)) * heaviside(?x)",
        mean: "?b * gamma(1 + 1 / ?a)",
        variance: "?b^2 * (gamma(1 + 2 / ?a) - gamma(1 + 1 / ?a)^2)",
        mgf: None,
        entropy: Some("euler_gamma * (1 - 1 / ?a) + ln(?b / ?a) + 1"),
    },
    Family {
        name: "geometric",
        params: 1,
        pdf: "(1 - ?a)^(?x - 1) * ?a * heaviside(?x - 1/2)",
        cdf: "(1 - (1 - ?a)^floor(?x)) * heaviside(?x - 1/2)",
        mean: "1 / ?a",
        variance: "(1 - ?a) / ?a^2",
        mgf: Some("?a * exp(?x) / (1 - (1 - ?a) * exp(?x))"),
        entropy: Some("(-(1 - ?a) * ln(1 - ?a) - ?a * ln(?a)) / ?a"),
    },
    Family {
        name: "negative_binomial",
        params: 2,
        pdf: "gamma(?x + ?a) / (gamma(?x + 1) * gamma(?a)) * ?b^?a * (1 - ?b)^?x * heaviside(?x + 1/2)",
        cdf: "beta_reg(?a, floor(?x) + 1, ?b) * heaviside(?x + 1/2)",
        mean: "?a * (1 - ?b) / ?b",
        variance: "?a * (1 - ?b) / ?b^2",
        mgf: Some("(?b / (1 - (1 - ?b) * exp(?x)))^?a"),
        entropy: None,
    },
    Family {
        name: "f_dist",
        params: 2,
        pdf: "((?a * abs(?x))^?a * ?b^?b / (?a * abs(?x) + ?b)^(?a + ?b))^(1/2) / (abs(?x) * beta(?a / 2, ?b / 2)) * heaviside(?x)",
        cdf: "beta_reg(?a / 2, ?b / 2, ?a * ?x / (?a * ?x + ?b)) * heaviside(?x)",
        mean: "?b / (?b - 2)",
        variance: "2 * ?b^2 * (?a + ?b - 2) / (?a * (?b - 2)^2 * (?b - 4))",
        mgf: None,
        entropy: None,
    },
    Family {
        name: "pareto",
        params: 2,
        pdf: "?b * ?a^?b / abs(?x)^(?b + 1) * heaviside(?x - ?a)",
        cdf: "(1 - (?a / abs(?x))^?b) * heaviside(?x - ?a)",
        mean: "?b * ?a / (?b - 1)",
        variance: "?a^2 * ?b / ((?b - 1)^2 * (?b - 2))",
        mgf: None,
        entropy: Some("ln(?a / ?b) + 1 / ?b + 1"),
    },
    Family {
        name: "rayleigh",
        params: 1,
        pdf: "?x / ?a^2 * exp(-?x^2 / (2 * ?a^2)) * heaviside(?x)",
        cdf: "(1 - exp(-?x^2 / (2 * ?a^2))) * heaviside(?x)",
        mean: "?a * (pi / 2)^(1/2)",
        variance: "(4 - pi) / 2 * ?a^2",
        mgf: None,
        entropy: Some("1 + ln(?a / 2^(1/2)) + euler_gamma / 2"),
    },
    Family {
        name: "logistic",
        params: 2,
        pdf: "exp(-(?x - ?a) / ?b) / (?b * (1 + exp(-(?x - ?a) / ?b))^2)",
        cdf: "1 / (1 + exp(-(?x - ?a) / ?b))",
        mean: "?a",
        variance: "pi^2 * ?b^2 / 3",
        mgf: None,
        entropy: Some("ln(?b) + 2"),
    },
    Family {
        name: "gumbel",
        params: 2,
        pdf: "exp(-(?x - ?a) / ?b - exp(-(?x - ?a) / ?b)) / ?b",
        cdf: "exp(-exp(-(?x - ?a) / ?b))",
        mean: "?a + ?b * euler_gamma",
        variance: "pi^2 * ?b^2 / 6",
        mgf: Some("gamma(1 - ?b * ?x) * exp(?a * ?x)"),
        entropy: Some("ln(?b) + euler_gamma + 1"),
    },
];

/// Parameter-space check for the families of this module (`None` for a
/// family of the parent module).
pub(super) fn admissible(
    name: &str,
    params: &[Option<f64>],
) -> Option<bool> {
    let p = |i: usize| params.get(i).copied().flatten();
    let positive = |i: usize| p(i).is_none_or(|v| v > 0.0);
    let probability = |i: usize| p(i).is_none_or(|v| v > 0.0 && v <= 1.0);
    Some(match name {
        | "cauchy" | "laplace_dist" | "lognormal" | "logistic" | "gumbel" => positive(1),
        | "weibull" | "f_dist" | "pareto" => positive(0) && positive(1),
        | "geometric" => probability(0),
        | "negative_binomial" => positive(0) && probability(1),
        | "rayleigh" => positive(0),
        | _ => return None,
    })
}

/// Whether the moment of the given order (1: mean, 2: variance) exists for
/// the numeric parameters.
pub(super) fn moment_exists(
    name: &str,
    order: u8,
    params: &[Option<f64>],
) -> bool {
    let p = |i: usize| params.get(i).copied().flatten();
    match (name, order) {
        | ("cauchy", _) => false,
        | ("f_dist", 1) => p(1).is_none_or(|d| d > 2.0),
        | ("f_dist", _) => p(1).is_none_or(|d| d > 4.0),
        | ("pareto", 1) => p(1).is_none_or(|a| a > 1.0),
        | ("pareto", _) => p(1).is_none_or(|a| a > 2.0),
        | _ => true,
    }
}

/// Properties by family name: characteristic function, skewness, excess
/// kurtosis, median, mode and quantile function (of `?x`).
struct Extra {
    name: &'static str,
    charfn: Option<&'static str>,
    skewness: Option<&'static str>,
    kurtosis: Option<&'static str>,
    median: Option<&'static str>,
    mode: Option<&'static str>,
    quantile: Option<&'static str>,
}

const fn extra(
    name: &'static str,
    charfn: Option<&'static str>,
    skewness: Option<&'static str>,
    kurtosis: Option<&'static str>,
    median: Option<&'static str>,
    mode: Option<&'static str>,
    quantile: Option<&'static str>,
) -> Extra {
    Extra { name, charfn, skewness, kurtosis, median, mode, quantile }
}

const EXTRAS: [Extra; 21] = [
    extra(
        "normal",
        Some("exp(I * ?a * ?x - ?b^2 * ?x^2 / 2)"),
        Some("0"),
        Some("0"),
        Some("?a"),
        Some("?a"),
        Some("?a - ?b * 2^(1/2) * erfcinv(2 * ?x)"),
    ),
    extra(
        "uniform",
        Some("(exp(I * ?b * ?x) - exp(I * ?a * ?x)) / (I * ?x * (?b - ?a))"),
        Some("0"),
        Some("-6/5"),
        Some("(?a + ?b) / 2"),
        None,
        Some("?a + ?x * (?b - ?a)"),
    ),
    extra(
        "exponential",
        Some("?a / (?a - I * ?x)"),
        Some("2"),
        Some("6"),
        Some("ln(2) / ?a"),
        Some("0"),
        Some("ln(1 / (1 - ?x)) / ?a"),
    ),
    extra(
        "bernoulli",
        Some("1 - ?a + ?a * exp(I * ?x)"),
        Some("(1 - 2 * ?a) / (?a * (1 - ?a))^(1/2)"),
        Some("(1 - 6 * ?a * (1 - ?a)) / (?a * (1 - ?a))"),
        None,
        None,
        None,
    ),
    extra(
        "binomial_dist",
        Some("(1 - ?b + ?b * exp(I * ?x))^?a"),
        Some("(1 - 2 * ?b) / (?a * ?b * (1 - ?b))^(1/2)"),
        Some("(1 - 6 * ?b * (1 - ?b)) / (?a * ?b * (1 - ?b))"),
        None,
        Some("floor((?a + 1) * ?b)"),
        None,
    ),
    extra(
        "poisson",
        Some("exp(?a * (exp(I * ?x) - 1))"),
        Some("1 / ?a^(1/2)"),
        Some("1 / ?a"),
        None,
        Some("floor(?a)"),
        None,
    ),
    extra(
        "gamma_dist",
        Some("(1 - I * ?x / ?b)^(-?a)"),
        Some("2 / ?a^(1/2)"),
        Some("6 / ?a"),
        None,
        Some("(?a - 1) / ?b"),
        None,
    ),
    extra(
        "beta_dist",
        None,
        Some("2 * (?b - ?a) * (?a + ?b + 1)^(1/2) / ((?a + ?b + 2) * (?a * ?b)^(1/2))"),
        Some("6 * ((?a - ?b)^2 * (?a + ?b + 1) - ?a * ?b * (?a + ?b + 2)) / (?a * ?b * (?a + ?b + 2) * (?a + ?b + 3))"),
        None,
        Some("(?a - 1) / (?a + ?b - 2)"),
        None,
    ),
    extra("student_t", None, Some("0"), Some("6 / (?a - 4)"), Some("0"), Some("0"), None),
    extra(
        "chi_squared",
        Some("(1 - 2 * I * ?x)^(-?a / 2)"),
        Some("(8 / ?a)^(1/2)"),
        Some("12 / ?a"),
        None,
        Some("?a - 2"),
        None,
    ),
    extra(
        "cauchy",
        Some("exp(I * ?a * ?x - ?b * abs(?x))"),
        None,
        None,
        Some("?a"),
        Some("?a"),
        Some("?a + ?b * tan(pi * (?x - 1/2))"),
    ),
    extra(
        "laplace_dist",
        Some("exp(I * ?a * ?x) / (1 + ?b^2 * ?x^2)"),
        Some("0"),
        Some("3"),
        Some("?a"),
        Some("?a"),
        Some("?a - ?b * sign(?x - 1/2) * ln(1 - 2 * abs(?x - 1/2))"),
    ),
    extra(
        "lognormal",
        None,
        Some("(exp(?b^2) + 2) * (exp(?b^2) - 1)^(1/2)"),
        Some("exp(4 * ?b^2) + 2 * exp(3 * ?b^2) + 3 * exp(2 * ?b^2) - 6"),
        Some("exp(?a)"),
        Some("exp(?a - ?b^2)"),
        Some("exp(?a - ?b * 2^(1/2) * erfcinv(2 * ?x))"),
    ),
    extra(
        "weibull",
        None,
        None,
        None,
        Some("?b * ln(2)^(1 / ?a)"),
        Some("?b * ((?a - 1) / ?a)^(1 / ?a)"),
        Some("?b * (-ln(1 - ?x))^(1 / ?a)"),
    ),
    extra(
        "geometric",
        Some("?a * exp(I * ?x) / (1 - (1 - ?a) * exp(I * ?x))"),
        Some("(2 - ?a) / (1 - ?a)^(1/2)"),
        Some("6 + ?a^2 / (1 - ?a)"),
        Some("ceil(-ln(2) / ln(1 - ?a))"),
        Some("1"),
        Some("ceil(ln(1 - ?x) / ln(1 - ?a))"),
    ),
    extra(
        "negative_binomial",
        Some("(?b / (1 - (1 - ?b) * exp(I * ?x)))^?a"),
        Some("(2 - ?b) / (?a * (1 - ?b))^(1/2)"),
        Some("6 / ?a + ?b^2 / (?a * (1 - ?b))"),
        None,
        None,
        None,
    ),
    extra("f_dist", None, None, None, None, None, None),
    extra(
        "pareto",
        None,
        None,
        None,
        Some("?a * 2^(1 / ?b)"),
        Some("?a"),
        Some("?a * (1 - ?x)^(-1 / ?b)"),
    ),
    extra(
        "rayleigh",
        None,
        Some("2 * pi^(1/2) * (pi - 3) / (4 - pi)^(3/2)"),
        Some("-(6 * pi^2 - 24 * pi + 16) / (4 - pi)^2"),
        Some("?a * (2 * ln(2))^(1/2)"),
        Some("?a"),
        Some("?a * (-2 * ln(1 - ?x))^(1/2)"),
    ),
    extra(
        "logistic",
        Some("exp(I * ?a * ?x) * pi * ?b * ?x / sinh(pi * ?b * ?x)"),
        Some("0"),
        Some("6/5"),
        Some("?a"),
        Some("?a"),
        Some("?a + ?b * ln(?x / (1 - ?x))"),
    ),
    extra(
        "gumbel",
        Some("gamma(1 - I * ?b * ?x) * exp(I * ?a * ?x)"),
        Some("12 * 6^(1/2) * zeta(3) / pi^3"),
        Some("12/5"),
        Some("?a - ?b * ln(ln(2))"),
        Some("?a"),
        Some("?a - ?b * ln(-ln(?x))"),
    ),
];

/// The parameter count of a family of either table.
fn parameter_count(name: &str) -> Option<usize> {
    CORE_FAMILIES
        .iter()
        .chain(FAMILIES.iter())
        .find(|f| f.name == name)
        .map(|f| usize::from(f.params))
}

fn family_of(name: &str) -> Option<&'static Family> {
    CORE_FAMILIES.iter().chain(FAMILIES.iter()).find(|f| f.name == name)
}

/// Instantiates pattern text over named variables.
fn instantiate(
    graph: &mut Graph,
    text: &str,
    names: &[&str],
    args: &[NodeId],
) -> Option<NodeId> {
    let mut vars = VarNames::default();
    for name in names {
        vars.index(name);
    }
    let pat = Pat::parse(text, graph, &mut vars).ok()?;
    pat.instantiate(graph, args)
}

/// A distribution term: the family name and its parameters.
fn distribution(
    graph: &mut Graph,
    node: NodeId,
) -> Option<(String, Vec<NodeId>)> {
    let node = best(graph, node)?;
    let name = graph.ops().get(graph.op(node)).name.to_string();
    let params = graph.children(node).to_vec();
    (parameter_count(&name)? == params.len()).then_some((name, params))
}

/// `a, b, x` instantiation for a family template.
fn family_template(
    graph: &mut Graph,
    text: &str,
    params: &[NodeId],
    x: NodeId,
) -> Option<NodeId> {
    let a = *params.first()?;
    let b = params.get(1).copied().unwrap_or(a);
    instantiate(graph, text, &["a", "b", "x"], &[a, b, x])
}

fn heaviside_free(text: &str) -> String {
    let mut out = text.to_owned();
    while let Some(start) = out.find("heaviside(") {
        let bytes = out.as_bytes();
        let mut depth = 0_i32;
        let mut end = start;
        for (k, &byte) in bytes.iter().enumerate().skip(start) {
            match byte {
                | b'(' => depth += 1,
                | b')' => {
                    depth -= 1;
                    if depth == 0 {
                        end = k + 1;
                        break;
                    }
                },
                | _ => {},
            }
        }
        let before = out[..start].trim_end_matches(' ');
        let after = &out[end..];
        let mut head = before.strip_suffix('*').map_or(before, str::trim_end).to_owned();
        let mut tail = after.to_owned();
        if head.is_empty() {
            if let Some(rest) = tail.trim_start().strip_prefix("* ") {
                tail = rest.to_owned();
            } else if tail.trim_start().starts_with('/') {
                "1".clone_into(&mut head);
            }
        }
        let joined = if tail.starts_with(' ') || tail.is_empty() || head.is_empty() { format!("{head}{tail}") } else { format!("{head} {tail}") };
        out = joined;
    }
    out
}

fn is_equal(
    cx: &mut Cx<'_>,
    a: NodeId,
    b: NodeId,
) -> bool {
    let d = difference(cx.graph, a, b);
    cx.is_zero(d)
}

#[derive(Copy, Clone)]
enum Kind {
    CharFn,
    Skewness,
    Kurtosis,
    Median,
    Mode,
    Quantile,
    Std,
    Moment,
    CentralMoment,
    ExpectOf,
    SumDist,
    IidSum,
    IidMean,
    ScaleDist,
    ShiftDist,
    Bayes,
    BinomialTest,
    FisherExact,
    Wilson,
    MeanInterval,
    ProportionZ,
    RSquared,
    Residuals,
}

struct Extension {
    op: OpId,
    kind: Kind,
}

fn extras_of(name: &str) -> Option<&'static Extra> {
    EXTRAS.iter().find(|e| e.name == name)
}

fn property(
    graph: &mut Graph,
    args: &[NodeId],
    pick: fn(&Extra) -> Option<&'static str>,
    with_argument: bool,
) -> Option<NodeId> {
    let (name, params) = distribution(graph, *args.first()?)?;
    let text = pick(extras_of(&name)?)?;
    let x = if with_argument { *args.get(1)? } else { *params.first()? };
    if text.contains('I') && graph.ops().lookup("I").is_none() {
        return None;
    }
    family_template(graph, text, &params, x)
}

/// The raw moment `E[X^k]` from derivatives of the moment generating
/// function at zero.
fn raw_moment(
    cx: &mut Cx<'_>,
    dist: NodeId,
    k: u32,
) -> Option<NodeId> {
    if k == 0 {
        return Some(cx.graph.int(1));
    }
    let (name, params) = distribution(cx.graph, dist)?;
    if name == "uniform" {
        let (a, b) = (*params.first()?, *params.get(1)?);
        let text = "(?b^(?k + 1) - ?a^(?k + 1)) / ((?k + 1) * (?b - ?a))";
        let kn = cx.graph.int(i64::from(k));
        let term = instantiate(cx.graph, text, &["a", "b", "k"], &[a, b, kn])?;
        return Some(cx.simplify(term));
    }
    let family = family_of(&name)?;
    let mgf_text = family.mgf?;
    let t = {
        let s = cx.graph.interner_mut().fresh_symbol("t");
        cx.graph.symbol_node(s)
    };
    let mgf = family_template(cx.graph, mgf_text, &params, t)?;
    let mut d = cx.simplify(mgf);
    for _ in 0..k {
        let next = derivative(cx.graph, d, t)?;
        d = cx.simplify(next);
    }
    let zero = cx.graph.int(0);
    let at_zero = cx.graph.substitute(d, t, zero);
    Some(cx.simplify(at_zero))
}

fn central_moment(
    cx: &mut Cx<'_>,
    dist: NodeId,
    k: u32,
) -> Option<NodeId> {
    let mean = raw_moment(cx, dist, 1)?;
    let mut terms = Vec::new();
    let mut binomial = BigInt::one();
    for j in 0..=k {
        if j > 0 {
            binomial = binomial * (k - j + 1) / j;
        }
        let m = raw_moment(cx, dist, j)?;
        let exponent = cx.graph.int(i64::from(k - j));
        let minus_mean = negate(cx.graph, mean);
        let power = cx.graph.node(core::POW, &[minus_mean, exponent]);
        let coefficient = cx.graph.num(Number::Int(binomial.clone()));
        terms.push(cx.graph.node(core::MUL, &[coefficient, m, power]));
    }
    let total = cx.graph.node(core::ADD, &terms);
    Some(cx.simplify(total))
}

fn support(
    graph: &mut Graph,
    name: &str,
    params: &[NodeId],
) -> Option<(NodeId, NodeId, bool)> {
    let infinity = call(graph, "oo", &[])?;
    let minus_infinity = negate(graph, infinity);
    let zero = graph.int(0);
    let one = graph.int(1);
    Some(match name {
        | "normal" | "cauchy" | "laplace_dist" | "logistic" | "gumbel" | "student_t" => (minus_infinity, infinity, false),
        | "uniform" => (*params.first()?, *params.get(1)?, false),
        | "exponential" | "gamma_dist" | "chi_squared" | "weibull" | "rayleigh" | "f_dist" | "lognormal" => {
            (zero, infinity, false)
        },
        | "beta_dist" => (zero, one, false),
        | "pareto" => (*params.first()?, infinity, false),
        | "bernoulli" => (zero, one, true),
        | "binomial_dist" => (zero, *params.first()?, true),
        | "poisson" | "negative_binomial" => (zero, infinity, true),
        | "geometric" => (one, infinity, true),
        | _ => return None,
    })
}

fn expectation_of(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[g, x, dist] = args else { return None };
    let (name, params) = distribution(cx.graph, dist)?;
    let family = family_of(&name)?;
    let (lo, hi, discrete) = support(cx.graph, &name, &params)?;
    let density = family_template(cx.graph, &heaviside_free(family.pdf), &params, x)?;
    let integrand = cx.graph.node(core::MUL, &[g, density]);
    let head = if discrete { "sum" } else { "defint" };
    call(cx.graph, head, &[integrand, x, lo, hi])
}

/// Replaces the parameter pattern of `text` with the parameters of two
/// distributions.
fn pair_template(
    graph: &mut Graph,
    text: &str,
    first: &[NodeId],
    second: &[NodeId],
) -> Option<NodeId> {
    let a = *first.first()?;
    let b = first.get(1).copied().unwrap_or(a);
    let c = *second.first()?;
    let d = second.get(1).copied().unwrap_or(c);
    instantiate(graph, text, &["a", "b", "c", "d"], &[a, b, c, d])
}

fn sum_dist(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[first, second] = args else { return None };
    let ((n1, p1), (n2, p2)) = (distribution(cx.graph, first)?, distribution(cx.graph, second)?);
    let same = |cx: &mut Cx<'_>, a: NodeId, b: NodeId| is_equal(cx, a, b);
    let text = match (n1.as_str(), n2.as_str()) {
        | ("normal", "normal") => "normal(?a + ?c, (?b^2 + ?d^2)^(1/2))",
        | ("poisson", "poisson") => "poisson(?a + ?c)",
        | ("chi_squared", "chi_squared") => "chi_squared(?a + ?c)",
        | ("cauchy", "cauchy") => "cauchy(?a + ?c, ?b + ?d)",
        | ("exponential", "exponential") if same(cx, p1[0], p2[0]) => "gamma_dist(2, ?a)",
        | ("gamma_dist", "gamma_dist") if same(cx, p1[1], p2[1]) => "gamma_dist(?a + ?c, ?b)",
        | ("exponential", "gamma_dist") if same(cx, p1[0], p2[1]) => "gamma_dist(?c + 1, ?a)",
        | ("gamma_dist", "exponential") if same(cx, p1[1], p2[0]) => "gamma_dist(?a + 1, ?b)",
        | ("binomial_dist", "binomial_dist") if same(cx, p1[1], p2[1]) => "binomial_dist(?a + ?c, ?b)",
        | ("bernoulli", "bernoulli") if same(cx, p1[0], p2[0]) => "binomial_dist(2, ?a)",
        | ("bernoulli", "binomial_dist") if same(cx, p1[0], p2[1]) => "binomial_dist(?c + 1, ?a)",
        | ("binomial_dist", "bernoulli") if same(cx, p1[1], p2[0]) => "binomial_dist(?a + 1, ?b)",
        | ("negative_binomial", "negative_binomial") if same(cx, p1[1], p2[1]) => "negative_binomial(?a + ?c, ?b)",
        | _ => return None,
    };
    pair_template(cx.graph, text, &p1, &p2)
}

fn iid(
    cx: &mut Cx<'_>,
    args: &[NodeId],
    mean: bool,
) -> Option<NodeId> {
    let &[dist, n] = args else { return None };
    let (name, params) = distribution(cx.graph, dist)?;
    let text = if mean {
        match name.as_str() {
            | "normal" => "normal(?a, ?b / ?x^(1/2))",
            | "cauchy" => "cauchy(?a, ?b)",
            | "gamma_dist" => "gamma_dist(?x * ?a, ?x * ?b)",
            | "exponential" => "gamma_dist(?x, ?x * ?a)",
            | _ => return None,
        }
    } else {
        match name.as_str() {
            | "normal" => "normal(?x * ?a, ?x^(1/2) * ?b)",
            | "gamma_dist" => "gamma_dist(?x * ?a, ?b)",
            | "exponential" => "gamma_dist(?x, ?a)",
            | "poisson" => "poisson(?x * ?a)",
            | "binomial_dist" => "binomial_dist(?x * ?a, ?b)",
            | "bernoulli" => "binomial_dist(?x, ?a)",
            | "chi_squared" => "chi_squared(?x * ?a)",
            | "cauchy" => "cauchy(?x * ?a, ?x * ?b)",
            | "negative_binomial" => "negative_binomial(?x * ?a, ?b)",
            | _ => return None,
        }
    };
    if f64_of(cx.graph, n).is_some_and(|v| v <= 0.0) {
        return None;
    }
    family_template(cx.graph, text, &params, n)
}

fn affine(
    cx: &mut Cx<'_>,
    args: &[NodeId],
    scale: bool,
) -> Option<NodeId> {
    let &[amount, dist] = args else { return None };
    let (name, params) = distribution(cx.graph, dist)?;
    let text = if scale {
        if f64_of(cx.graph, amount).is_some_and(|c| c <= 0.0) {
            return None;
        }
        match name.as_str() {
            | "normal" => "normal(?x * ?a, abs(?x) * ?b)",
            | "uniform" => "uniform(?x * ?a, ?x * ?b)",
            | "exponential" => "exponential(?a / ?x)",
            | "gamma_dist" => "gamma_dist(?a, ?b / ?x)",
            | "cauchy" => "cauchy(?x * ?a, abs(?x) * ?b)",
            | "laplace_dist" => "laplace_dist(?x * ?a, abs(?x) * ?b)",
            | "lognormal" => "lognormal(?a + ln(?x), ?b)",
            | "weibull" => "weibull(?a, ?x * ?b)",
            | "pareto" => "pareto(?x * ?a, ?b)",
            | "rayleigh" => "rayleigh(?x * ?a)",
            | "logistic" => "logistic(?x * ?a, abs(?x) * ?b)",
            | "gumbel" => "gumbel(?x * ?a, ?x * ?b)",
            | _ => return None,
        }
    } else {
        match name.as_str() {
            | "normal" => "normal(?a + ?x, ?b)",
            | "uniform" => "uniform(?a + ?x, ?b + ?x)",
            | "cauchy" => "cauchy(?a + ?x, ?b)",
            | "laplace_dist" => "laplace_dist(?a + ?x, ?b)",
            | "logistic" => "logistic(?a + ?x, ?b)",
            | "gumbel" => "gumbel(?a + ?x, ?b)",
            | _ => return None,
        }
    };
    family_template(cx.graph, text, &params, amount)
}

fn bayes(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[prior, model, data] = args else { return None };
    let (name, params) = distribution(cx.graph, prior)?;
    let model = best(cx.graph, model)?;
    let model_name = cx.graph.ops().get(cx.graph.op(model)).name.to_string();
    let model_name = if cx.graph.children(model).is_empty() {
        cx.graph.display(model)
    } else {
        model_name
    };
    let model_args = cx.graph.children(model).to_vec();
    let data = best(cx.graph, data)?;
    let items = items_of(cx.graph, data)?;
    if items.is_empty() {
        return None;
    }
    let n = cx.graph.int(i64::try_from(items.len()).ok()?);
    let total = sum_node(cx.graph, &items);
    let numbers: Option<Vec<f64>> = items.iter().map(|&v| f64_of(cx.graph, v)).collect();
    let text = match (name.as_str(), model_name.as_str()) {
        | ("beta_dist", "bernoulli") => {
            if numbers.is_some_and(|v| v.iter().any(|&x| x.abs() > 1e-12 && (x - 1.0).abs() > 1e-12)) {
                return None;
            }
            "beta_dist(?a + ?s, ?b + ?n - ?s)"
        },
        | ("beta_dist", "binomial_trials") => "beta_dist(?a + ?s, ?b + ?n * ?m - ?s)",
        | ("beta_dist", "geometric") => "beta_dist(?a + ?n, ?b + ?s - ?n)",
        | ("gamma_dist", "poisson") => "gamma_dist(?a + ?s, ?b + ?n)",
        | ("gamma_dist", "exponential") => "gamma_dist(?a + ?n, ?b + ?s)",
        | ("normal", "normal_known_sd") => {
            "normal((?a / ?b^2 + ?s / ?m^2) / (1 / ?b^2 + ?n / ?m^2), (1 / ?b^2 + ?n / ?m^2)^(-1/2))"
        },
        | _ => return None,
    };
    let a = *params.first()?;
    let b = params.get(1).copied().unwrap_or(a);
    let m = model_args.first().copied().unwrap_or(a);
    instantiate(cx.graph, text, &["a", "b", "s", "n", "m"], &[a, b, total, n, m])
}

// ----------------------------------------------------------------------
// Exact tests
// ----------------------------------------------------------------------

fn integer(
    graph: &Graph,
    node: NodeId,
    cap: i64,
) -> Option<u64> {
    let v = graph.number_of(node).and_then(Number::to_i64).filter(|&v| (0..=cap).contains(&v))?;
    u64::try_from(v).ok()
}

fn binomial_test(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[k, n, p0] = args else { return None };
    let (k, n) = (integer(cx.graph, k, 2000)?, integer(cx.graph, n, 2000)?);
    if k > n {
        return None;
    }
    let p = cx.graph.number_of(p0).and_then(Number::to_rational)?;
    if p <= BigRational::zero() || p >= BigRational::one() {
        return None;
    }
    let q = BigRational::one() - &p;
    let pmf = |i: u64| -> Option<BigRational> {
        let c = choose(&BigInt::from(n), i)?;
        Some(BigRational::from_integer(c) * num_traits::pow::pow(p.clone(), usize::try_from(i).ok()?) * num_traits::pow::pow(q.clone(), usize::try_from(n - i).ok()?))
    };
    let observed = pmf(k)?;
    let mut total = BigRational::zero();
    for i in 0..=n {
        let value = pmf(i)?;
        if value <= observed {
            total += value;
        }
    }
    Some(cx.graph.num(Number::rat(total)))
}

fn fisher_exact(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[a, b, c, d] = args else { return None };
    let (a, b, c, d) = (integer(cx.graph, a, 1000)?, integer(cx.graph, b, 1000)?, integer(cx.graph, c, 1000)?, integer(cx.graph, d, 1000)?);
    let (row1, row2, col1) = (a + b, c + d, a + c);
    let total = row1 + row2;
    let weight = |x: u64| -> Option<BigInt> { Some(choose(&BigInt::from(row1), x)? * choose(&BigInt::from(row2), col1 - x)?) };
    let observed = weight(a)?;
    let mut sum = BigInt::zero();
    for x in col1.saturating_sub(row2)..=col1.min(row1) {
        let w = weight(x)?;
        if w <= observed {
            sum += w;
        }
    }
    let denominator = choose(&BigInt::from(total), col1)?;
    Some(cx.graph.num(Number::rat(BigRational::new(sum, denominator))))
}

fn template_over(
    graph: &mut Graph,
    text: &str,
    names: &[&str],
    args: &[NodeId],
) -> Option<NodeId> {
    instantiate(graph, text, names, args)
}

fn regression_coefficients(
    cx: &mut Cx<'_>,
    xs: NodeId,
    ys: NodeId,
) -> Option<(NodeId, NodeId)> {
    let request = call(cx.graph, "linear_regression", &[xs, ys])?;
    let result = cx.simplify(request);
    let items = items_of(cx.graph, result)?;
    let &[intercept, slope] = items.as_slice() else { return None };
    Some((intercept, slope))
}

impl Extension {
    #[allow(clippy::too_many_lines)]
    fn compute(
        &self,
        cx: &mut Cx<'_>,
        args: &[NodeId],
    ) -> Option<NodeId> {
        match self.kind {
            | Kind::CharFn => property(cx.graph, args, |e| e.charfn, true),
            | Kind::Skewness => property(cx.graph, args, |e| e.skewness, false),
            | Kind::Kurtosis => property(cx.graph, args, |e| e.kurtosis, false),
            | Kind::Median => property(cx.graph, args, |e| e.median, false),
            | Kind::Mode => property(cx.graph, args, |e| e.mode, false),
            | Kind::Quantile => property(cx.graph, args, |e| e.quantile, true),
            | Kind::Std => {
                let variance = call(cx.graph, "variance_of", &[*args.first()?])?;
                let half = cx.graph.num(Number::fraction(1, 2)?);
                Some(cx.graph.node(core::POW, &[variance, half]))
            },
            | Kind::Moment | Kind::CentralMoment => {
                let &[dist, k] = args else { return None };
                let k = u32::try_from(integer(cx.graph, k, 8)?).ok()?;
                if matches!(self.kind, Kind::Moment) {
                    raw_moment(cx, dist, k)
                } else {
                    central_moment(cx, dist, k)
                }
            },
            | Kind::ExpectOf => expectation_of(cx, args),
            | Kind::SumDist => sum_dist(cx, args),
            | Kind::IidSum => iid(cx, args, false),
            | Kind::IidMean => iid(cx, args, true),
            | Kind::ScaleDist => affine(cx, args, true),
            | Kind::ShiftDist => affine(cx, args, false),
            | Kind::Bayes => bayes(cx, args),
            | Kind::BinomialTest => binomial_test(cx, args),
            | Kind::FisherExact => fisher_exact(cx, args),
            | Kind::Wilson => {
                let &[k, n, z] = args else { return None };
                let names = ["s", "n", "z"];
                let lower = template_over(
                    cx.graph,
                    "(?s + ?z^2 / 2 - ?z * (?s * (?n - ?s) / ?n + ?z^2 / 4)^(1/2)) / (?n + ?z^2)",
                    &names,
                    &[k, n, z],
                )?;
                let upper = template_over(
                    cx.graph,
                    "(?s + ?z^2 / 2 + ?z * (?s * (?n - ?s) / ?n + ?z^2 / 4)^(1/2)) / (?n + ?z^2)",
                    &names,
                    &[k, n, z],
                )?;
                Some(list_node(cx.graph, &[lower, upper]))
            },
            | Kind::MeanInterval => {
                let &[m, s, n, z] = args else { return None };
                let names = ["m", "s", "n", "z"];
                let lower = template_over(cx.graph, "?m - ?z * ?s / ?n^(1/2)", &names, &[m, s, n, z])?;
                let upper = template_over(cx.graph, "?m + ?z * ?s / ?n^(1/2)", &names, &[m, s, n, z])?;
                Some(list_node(cx.graph, &[lower, upper]))
            },
            | Kind::ProportionZ => {
                let &[k, n, p0] = args else { return None };
                let names = ["s", "n", "m"];
                let z = "(?s / ?n - ?m) / (?m * (1 - ?m) / ?n)^(1/2)";
                let zt = template_over(cx.graph, z, &names, &[k, n, p0])?;
                let p = template_over(
                    cx.graph,
                    &format!("2 * (1 - cdf(normal(0, 1), abs({z})))"),
                    &names,
                    &[k, n, p0],
                )?;
                Some(list_node(cx.graph, &[zt, p]))
            },
            | Kind::RSquared => {
                let &[xs, ys] = args else { return None };
                template_over(cx.graph, "correlation(?a, ?b)^2", &["a", "b"], &[xs, ys])
            },
            | Kind::Residuals => {
                let &[xs, ys] = args else { return None };
                let (intercept, slope) = regression_coefficients(cx, xs, ys)?;
                let (bx, by) = (best(cx.graph, xs)?, best(cx.graph, ys)?);
                let (xs, ys) = (items_of(cx.graph, bx)?, items_of(cx.graph, by)?);
                let mut out = Vec::with_capacity(xs.len());
                for (&x, &y) in xs.iter().zip(&ys) {
                    let fitted = template_over(cx.graph, "?a + ?b * ?x", &["a", "b", "x"], &[intercept, slope, x])?;
                    let residual = difference(cx.graph, y, fitted);
                    out.push(cx.simplify(residual));
                }
                Some(list_node(cx.graph, &out))
            },
        }
    }
}

impl Kernel for Extension {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        self.compute(cx, &args).map_or(Outcome::Pass, Outcome::Equal)
    }

    fn revisit(&self) -> bool {
        true
    }
}

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for family in &FAMILIES {
        i.op(OpDescriptor::new(family.name, Arity::Fixed(family.params)))?;
    }
    // Inert model descriptors for `bayes_update`.
    i.op(OpDescriptor::new("binomial_trials", Arity::Fixed(1)))?;
    i.op(OpDescriptor::new("normal_known_sd", Arity::Fixed(1)))?;
    for (name, arity, kind) in [
        ("charfn", 2, Kind::CharFn),
        ("skewness_of", 1, Kind::Skewness),
        ("kurtosis_of", 1, Kind::Kurtosis),
        ("median_of", 1, Kind::Median),
        ("mode_of", 1, Kind::Mode),
        ("quantile_of", 2, Kind::Quantile),
        ("std_of", 1, Kind::Std),
        ("moment_of", 2, Kind::Moment),
        ("central_moment_of", 2, Kind::CentralMoment),
        ("expect_of", 3, Kind::ExpectOf),
        ("sum_dist", 2, Kind::SumDist),
        ("iid_sum", 2, Kind::IidSum),
        ("iid_mean", 2, Kind::IidMean),
        ("scale_dist", 2, Kind::ScaleDist),
        ("shift_dist", 2, Kind::ShiftDist),
        ("bayes_update", 3, Kind::Bayes),
        ("binomial_test", 3, Kind::BinomialTest),
        ("fisher_exact", 4, Kind::FisherExact),
        ("wilson_interval", 3, Kind::Wilson),
        ("mean_interval", 4, Kind::MeanInterval),
        ("proportion_z_test", 3, Kind::ProportionZ),
        ("r_squared", 2, Kind::RSquared),
        ("residuals", 2, Kind::Residuals),
    ] {
        let op = i.op(OpDescriptor::new(name, Arity::Fixed(arity)).flags(OpFlags::HEAVY).cost(100))?;
        i.kernel(&format!("stats/{name}"), Tier::Reduce, Extension { op, kind });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
        use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[crate::rules::stats::stats()], src)
    }

    fn value(
        src: &str,
        bindings: &[(&str, f64)],
    ) -> f64 {
        crate::rules::testing::numeric(&[crate::rules::stats::stats()], src, bindings, 1e-12).0
    }

    fn close(
        a: f64,
        b: f64,
        what: &str,
    ) {
        assert!((a - b).abs() < 1e-7 * (1.0 + b.abs()), "{what}: {a} vs {b}");
    }

    #[test]
    fn new_families_have_closed_forms() {
        assert_eq!(run("cdf(cauchy(0, 1), 1)"), "3/4");
        assert_eq!(run("expectation(laplace_dist(2, 3))"), "2");
        assert_eq!(run("variance_of(laplace_dist(2, 3))"), "18");
        assert_eq!(run("expectation(geometric(1/4))"), "4");
        assert_eq!(run("variance_of(geometric(1/2))"), "2");
        assert_eq!(run("expectation(negative_binomial(3, 1/2))"), "3");
        assert_eq!(run("expectation(logistic(1, 2))"), "1");
        assert_eq!(run("variance_of(rayleigh(1))"), "2 - 1/2*pi");
        assert_eq!(run("expectation(pareto(1, 3))"), "3/2");
        assert_eq!(run("expectation(f_dist(3, 6))"), "3/2");
        // The Cauchy mean and a Pareto variance with alpha <= 2 do not exist.
        let (text, reduced) = crate::rules::testing::reduce_with(&[crate::rules::stats::stats()], "expectation(cauchy(0, 1))", &[]);
        assert!(!reduced, "{text}");
        let (text, reduced) = crate::rules::testing::reduce_with(&[crate::rules::stats::stats()], "variance_of(pareto(1, 2))", &[]);
        assert!(!reduced, "{text}");
    }

    #[test]
    fn densities_integrate_to_the_cdf() {
        // Weibull and lognormal cdf against the pdf at a point via a symmetric difference quotient.
        for (dist, x) in [("weibull(2, 3)", 2.0), ("lognormal(0, 1/2)", 1.5), ("gumbel(1, 2)", 0.5), ("logistic(0, 1)", 0.7), ("f_dist(3, 5)", 1.2), ("pareto(1, 2)", 2.0), ("rayleigh(2)", 1.0), ("laplace_dist(0, 2)", 0.3)] {
            let h = 1e-5;
            let upper = value(&format!("cdf({dist}, x)"), &[("x", x + h)]);
            let lower = value(&format!("cdf({dist}, x)"), &[("x", x - h)]);
            let density = value(&format!("pdf({dist}, x)"), &[("x", x)]);
            close((upper - lower) / (2.0 * h), density, dist);
        }
    }

    #[test]
    fn quantiles_invert_the_cdf() {
        for dist in ["normal(1, 2)", "exponential(3)", "weibull(2, 3)", "logistic(0, 1)", "gumbel(0, 1)", "laplace_dist(0, 1)", "rayleigh(2)", "pareto(1, 2)", "lognormal(0, 1)", "cauchy(0, 1)", "uniform(1, 4)"] {
            let p = 0.3;
            let q = value(&format!("quantile_of({dist}, p)"), &[("p", p)]);
            let back = value(&format!("cdf({dist}, x)"), &[("x", q)]);
            close(back, p, dist);
        }
    }

    #[test]
    fn derived_properties() {
        assert_eq!(run("skewness_of(exponential(2))"), "2");
        assert_eq!(run("kurtosis_of(normal(0, 1))"), "0");
        assert_eq!(run("kurtosis_of(uniform(0, 1))"), "-6/5");
        assert_eq!(run("median_of(exponential(1))"), "ln(2)");
        assert_eq!(run("mode_of(gamma_dist(3, 2))"), "1");
        assert_eq!(run("std_of(normal(0, 3))"), "3");
        assert_eq!(run("quantile_of(uniform(0, 4), 1/4)"), "1");
        assert_eq!(run("quantile_of(exponential(1), 1/2)"), "ln(2)");
    }

    #[test]
    fn characteristic_functions() {
        let sets = [crate::rules::stats::stats(), crate::rules::complex::complex()];
        let text = crate::rules::testing::simplify(&sets, "charfn(exponential(2), t)");
        assert!(text.contains('I'), "{text}");
        let text = crate::rules::testing::simplify(&sets, "charfn(cauchy(0, 1), t)");
        assert!(text.contains("abs(t)"), "{text}");
    }

    #[test]
    fn moments_from_the_mgf() {
        assert_eq!(run("moment_of(exponential(2), 3)"), "3/4");
        assert_eq!(run("moment_of(normal(0, 1), 4)"), "3");
        assert_eq!(run("moment_of(poisson(2), 2)"), "6");
        assert_eq!(run("moment_of(uniform(0, 1), 3)"), "1/4");
        assert_eq!(run("central_moment_of(normal(1, 2), 2)"), "4");
        assert_eq!(run("central_moment_of(exponential(1), 3)"), "2");
        assert_eq!(run("central_moment_of(poisson(3), 3)"), "3");
    }

    #[test]
    fn sums_and_affine_maps() {
        assert_eq!(run("sum_dist(normal(1, 3), normal(2, 4))"), "normal(3, 5)");
        assert_eq!(run("sum_dist(poisson(2), poisson(3))"), "poisson(5)");
        assert_eq!(run("sum_dist(exponential(2), exponential(2))"), "gamma_dist(2, 2)");
        assert_eq!(run("sum_dist(binomial_dist(3, 1/2), binomial_dist(4, 1/2))"), "binomial_dist(7, 1/2)");
        assert_eq!(run("sum_dist(bernoulli(1/3), bernoulli(1/3))"), "binomial_dist(2, 1/3)");
        assert_eq!(run("iid_sum(normal(1, 2), n)"), "normal(n, 2*n^(1/2))");
        assert_eq!(run("iid_mean(normal(1, 2), 16)"), "normal(1, 1/2)");
        assert_eq!(run("iid_sum(exponential(2), 5)"), "gamma_dist(5, 2)");
        assert_eq!(run("scale_dist(2, uniform(0, 1))"), "uniform(0, 2)");
        assert_eq!(run("scale_dist(2, exponential(1))"), "exponential(1/2)");
        assert_eq!(run("shift_dist(5, normal(0, 1))"), "normal(5, 1)");
        // The sum of two independent normals has the variance of the parts added.
        assert_eq!(run("variance_of(sum_dist(normal(0, 1), normal(0, 2)))"), "5");
    }

    #[test]
    fn conjugate_updates() {
        assert_eq!(run("bayes_update(beta_dist(1, 1), bernoulli, list(1, 0, 1, 1))"), "beta_dist(4, 2)");
        assert_eq!(run("bayes_update(beta_dist(2, 3), binomial_trials(10), list(4, 6))"), "beta_dist(12, 13)");
        assert_eq!(run("bayes_update(gamma_dist(2, 1), poisson, list(3, 1, 2))"), "gamma_dist(8, 4)");
        assert_eq!(run("bayes_update(gamma_dist(1, 1), exponential, list(2, 3))"), "gamma_dist(3, 6)");
        assert_eq!(run("bayes_update(normal(0, 1), normal_known_sd(1), list(2, 4))"), "normal(2, 1/3^(1/2))");
        assert_eq!(run("expectation(bayes_update(beta_dist(1, 1), bernoulli, list(1, 0, 1, 1)))"), "2/3");
    }

    #[test]
    fn exact_tests_and_intervals() {
        // 3 successes in 3 trials, p0 = 1/2: p = 2/8.
        assert_eq!(run("binomial_test(3, 3, 1/2)"), "1/4");
        assert_eq!(run("binomial_test(5, 10, 1/2)"), "1");
        // The tea-tasting table: p = 17/35 for [[3, 1], [1, 3]].
        assert_eq!(run("fisher_exact(3, 1, 1, 3)"), "17/35");
        assert_eq!(run("fisher_exact(10, 0, 0, 10)"), "1/92378");
        assert_eq!(run("mean_interval(10, 2, 4, 2)"), "list(8, 12)");
        let text = run("wilson_interval(5, 10, 2)");
        assert!(text.starts_with("list("), "{text}");
        assert_eq!(run("r_squared(list(1, 2, 3, 4), list(2, 4, 6, 8))"), "1");
        assert_eq!(run("residuals(list(0, 1, 2), list(1, 2, 2))"), "list(-1/6, 1/3, -1/6)");
    }
}
