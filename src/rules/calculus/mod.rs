//! Calculus: the operators of continuous analysis and their reduction
//! kernels.
//!
//! `diff(f, x)` is an *identity transformation request*: a heavy node equal
//! to the derivative of `f`. Kernels compete to reduce it — structurally
//! when every operator in `f` has known partial derivatives, numerically by
//! finite differences when a numeric answer is wanted and the structural
//! route is blocked.

mod diff;
mod gosper;
mod integrate;
mod limits;
mod powerseries;
pub(crate) use powerseries::laurent_expansion;
mod series;
#[cfg(test)]
mod tests;

use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Facts;
use crate::graph::OnReals;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::Pat;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::graph::VarNames;
use crate::graph::rule::Installer;

use super::elementary::elementary;
use super::poly::poly;

pub use diff::Partials;
pub use integrate::IntegralTable;
pub use integrate::TableFn;

/// Operator attribute: the limits of a one-argument function as its
/// argument tends to `+∞` and `-∞`, as nullary patterns. The limit kernel
/// consults it for functions it has no built-in knowledge of.
#[derive(Clone, Debug)]
pub struct AtInfinity {
    /// Limit at `+∞`.
    pub plus: Option<Pat>,
    /// Limit at `-∞`.
    pub minus: Option<Pat>,
}

/// Teaches the limit kernel the values of `name` at `±∞` (`None` where
/// there is no finite limit or it is not known).
///
/// # Errors
/// Fails when the operator is unknown or a value does not parse.
pub(crate) fn at_infinity(
    i: &mut Installer<'_>,
    name: &str,
    plus: Option<&str>,
    minus: Option<&str>,
) -> Result<(), RuleError> {
    let op = i.graph().ops().lookup(name).ok_or_else(|| RuleError::Invalid {
        rule: format!("limits of {name}"),
        reason: "unknown operator",
    })?;
    let mut parse = |text: Option<&str>| -> Result<Option<Pat>, RuleError> {
        text.map(|t| {
            Pat::parse(t, i.graph(), &mut VarNames::default())
                .map_err(|error| RuleError::Parse { rule: format!("limits of {name}: {t}"), error })
        })
        .transpose()
    };
    let plus = parse(plus)?;
    let minus = parse(minus)?;
    i.graph().ops_mut().set_attr(op, AtInfinity { plus, minus });
    Ok(())
}

/// Teaches the integrator an antiderivative rule.
pub(crate) fn integral_rule(
    i: &mut Installer<'_>,
    rule: TableFn,
) {
    if let Some(op) = i.graph().ops().lookup("integral") {
        let mut table = i.graph().ops().attr::<IntegralTable>(op).cloned().unwrap_or_default();
        table.0.push(rule);
        i.graph().ops_mut().set_attr(op, table);
    }
}

/// The calculus rule set.
#[must_use]
pub fn calculus() -> RuleSet {
    RuleSet::new("calculus", install).needs(elementary()).needs(poly())
}

/// Attaches the partial derivatives of operator `name`, one pattern per
/// argument, written over the variables `?a`, `?b`, `?c`, `?d`.
///
/// Other rule sets call this to teach differentiation about their own
/// operators.
///
/// # Errors
/// Fails when the operator is unknown or a pattern does not parse.
pub(crate) fn partials(
    i: &mut Installer<'_>,
    name: &str,
    texts: &[&str],
) -> Result<(), RuleError> {
    let invalid = |reason: &'static str| RuleError::Invalid {
        rule: format!("d/{name}"),
        reason,
    };
    let op = i
        .graph()
        .ops()
        .lookup(name)
        .ok_or_else(|| invalid("unknown operator"))?;
    let mut patterns = Vec::with_capacity(texts.len());
    for text in texts {
        if text.is_empty() {
            patterns.push(None);
            continue;
        }
        let mut vars = VarNames::default();
        for v in ["a", "b", "c", "d"] {
            vars.index(v);
        }
        let pat = Pat::parse(text, i.graph(), &mut vars).map_err(|error| RuleError::Parse {
            rule: format!("d/{name}: {text}"),
            error,
        })?;
        if vars.len() > 4 {
            return Err(invalid("partial derivatives may only use ?a, ?b, ?c, ?d"));
        }
        patterns.push(Some(pat));
    }
    i.graph().ops_mut().set_attr(op, Partials(patterns));
    Ok(())
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let diff_op = i.op(OpDescriptor::new("diff", Arity::Fixed(2))
        .flags(OpFlags::HEAVY.with(OpFlags::OPAQUE_ON_APPLY))
        .cost(100))?;

    partials(i, "pow", &["?b * ?a^(?b - 1)", "?a^?b * ln(?a)"])?;
    partials(i, "exp", &["exp(?a)"])?;
    partials(i, "ln", &["1/?a"])?;
    partials(i, "sin", &["cos(?a)"])?;
    partials(i, "cos", &["-sin(?a)"])?;
    partials(i, "tan", &["1 + tan(?a)^2"])?;
    partials(i, "asin", &["(1 - ?a^2)^(-1/2)"])?;
    partials(i, "acos", &["-(1 - ?a^2)^(-1/2)"])?;
    partials(i, "atan", &["1/(1 + ?a^2)"])?;
    partials(i, "sinh", &["cosh(?a)"])?;
    partials(i, "cosh", &["sinh(?a)"])?;
    partials(i, "tanh", &["1 - tanh(?a)^2"])?;
    partials(i, "sqrt", &["1/(2*sqrt(?a))"])?;
    partials(i, "abs", &["?a/abs(?a)"])?;
    partials(i, "cot", &["-(1 + cot(?a)^2)"])?;
    partials(i, "sec", &["sec(?a) * tan(?a)"])?;
    partials(i, "csc", &["-csc(?a) * cot(?a)"])?;
    partials(i, "acot", &["-1/(1 + ?a^2)"])?;
    partials(i, "asec", &["1/(?a^2 * (1 - ?a^(-2))^(1/2))"])?;
    partials(i, "acsc", &["-1/(?a^2 * (1 - ?a^(-2))^(1/2))"])?;
    partials(i, "coth", &["1 - coth(?a)^2"])?;
    partials(i, "sech", &["-sech(?a) * tanh(?a)"])?;
    partials(i, "csch", &["-csch(?a) * coth(?a)"])?;
    partials(i, "asinh", &["(?a^2 + 1)^(-1/2)"])?;
    partials(i, "acosh", &["(?a^2 - 1)^(-1/2)"])?;
    partials(i, "atanh", &["1/(1 - ?a^2)"])?;
    partials(i, "acoth", &["1/(1 - ?a^2)"])?;
    partials(i, "asech", &["-1/(?a * (1 - ?a^2)^(1/2))"])?;
    partials(i, "acsch", &["-1/(?a^2 * (1 + ?a^(-2))^(1/2))"])?;
    partials(i, "atan2", &["?b/(?a^2 + ?b^2)", "-?a/(?a^2 + ?b^2)"])?;
    partials(i, "log", &["-ln(?b)/(?a * ln(?a)^2)", "1/(?b * ln(?a))"])?;

    i.kernel(
        "calculus/diff",
        Tier::Reduce,
        diff::Differentiate { diff: diff_op },
    );
    // diffn(f, x, n): the n-th derivative for a literal n.
    let diffn = i.op(OpDescriptor::new("diffn", Arity::Fixed(3)).flags(OpFlags::HEAVY).cost(100))?;
    i.kernel("calculus/diffn", Tier::Reduce, diff::Repeated { diffn, diff: diff_op });
    i.kernel(
        "calculus/diff-numeric",
        Tier::Reduce,
        diff::FiniteDifference { diff: diff_op },
    );

    // Integration. `integral(f, x)` is an antiderivative (x stays free);
    // `defint(f, x, a, b)` binds x inside f.
    let integral = i.op(OpDescriptor::new("integral", Arity::Fixed(2)).flags(OpFlags::HEAVY).cost(100))?;
    let defint =
        i.op(OpDescriptor::new("defint", Arity::Fixed(4)).flags(OpFlags::HEAVY).cost(100).binder(1, 0b1))?;
    let infinity = i.op(OpDescriptor::new("oo", Arity::Fixed(0)).eval(|_| f64::INFINITY))?;
    i.graph().ops_mut().set_attr(infinity, OnReals(Facts::POSITIVE));
    let lookup = |i: &mut Installer<'_>, name: &'static str| {
        i.graph().ops().lookup(name).ok_or_else(|| RuleError::Invalid {
            rule: format!("calculus needs `{name}`"),
            reason: "operator of the elementary rule set is missing",
        })
    };
    let functions = integrate::Functions {
        diff: diff_op,
        exp: lookup(i, "exp")?,
        ln: lookup(i, "ln")?,
        sin: lookup(i, "sin")?,
        cos: lookup(i, "cos")?,
        tan: lookup(i, "tan")?,
        asin: lookup(i, "asin")?,
        acos: lookup(i, "acos")?,
        atan: lookup(i, "atan")?,
        sinh: lookup(i, "sinh")?,
        cosh: lookup(i, "cosh")?,
        tanh: lookup(i, "tanh")?,
        sqrt: lookup(i, "sqrt")?,
    };
    i.graph().ops_mut().set_attr(integral, functions);
    i.kernel("calculus/integral", Tier::Reduce, integrate::Antiderivative { integral, functions });
    i.kernel("calculus/defint", Tier::Reduce, integrate::Definite { defint, functions });
    i.kernel("calculus/quadrature", Tier::Reduce, integrate::Quadrature { defint });
    // limit(f, x, a) and limit(f, x, a, plus|minus); x is bound in f.
    let limit =
        i.op(OpDescriptor::new("limit", Arity::Variadic).flags(OpFlags::HEAVY).cost(100).binder(1, 0b1))?;
    i.kernel("calculus/limit", Tier::Reduce, limits::SymbolicLimit { limit, infinity, functions });
    i.kernel("calculus/limit-numeric", Tier::Reduce, limits::NumericLimit { limit });
    // Series, sums and products. The index of a sum or product is bound.
    let request = |name: &str, arity: u8| OpDescriptor::new(name, Arity::Fixed(arity)).flags(OpFlags::HEAVY).cost(100);
    let ops = series::SeriesOps {
        functions,
        infinity,
        defint,
        taylor: i.op(request("taylor", 4))?,
        laurent: i.op(request("laurent", 4))?,
        fourier: i.op(request("fourier_series", 4))?,
        sum: i.op(request("sum", 4).binder(1, 0b1))?,
        product: i.op(request("product", 4).binder(1, 0b1))?,
        converges: i.op(request("converges", 2).binder(1, 0b1))?,
    };
    i.kernel("calculus/series", Tier::Reduce, series::SeriesKernel { ops });
    // antidifference(t, k): T with T(k+1) - T(k) = t(k).
    let antidifference = i.op(request("antidifference", 2))?;
    i.kernel("calculus/antidifference", Tier::Reduce, Antidifference { op: antidifference });
    i.kernel("calculus/sum-numeric", Tier::Reduce, series::NumericSum { sum: ops.sum, product: ops.product });
    i.rewrites(
        Tier::Reduce,
        &[
            // The fundamental theorem, in the direction that removes both
            // requests at once.
            "calculus/ftc: diff(integral(?f, ?x), ?x) => ?f",
            "calculus/defint-empty: defint(?f, ?x, ?a, ?a) => 0",
        ],
    )?;
    Ok(())
}

/// A verified antiderivative of `f` with respect to `x`, for other rule
/// sets (differential equations, transforms). `None` if none is found or
/// the calculus rule set is not installed in this graph.
pub fn antiderivative(
    cx: &mut Cx<'_>,
    f: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let integral = cx.graph.ops().lookup("integral")?;
    let functions = *cx.graph.ops().attr::<integrate::Functions>(integral)?;
    integrate::antiderivative(cx, functions, f, x)
}

/// The derivative of the concrete term `f` with respect to the symbol
/// `x`, unsimplified. Parts that cannot be differentiated stay as `diff`
/// requests.
pub fn derivative(
    graph: &mut Graph,
    f: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let diff = graph.ops().lookup("diff")?;
    let symbol = graph.symbol_of(x)?;
    Some(diff::Differentiate { diff }.derive(graph, f, symbol, x))
}

/// Kernel for `antidifference(t, k)`, by Gosper's algorithm.
struct Antidifference {
    op: crate::graph::OpId,
}

impl crate::graph::Kernel for Antidifference {
    fn ops(&self) -> Vec<crate::graph::OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> crate::graph::Outcome {
        let &[t, k] = cx.graph.children(node) else {
            return crate::graph::Outcome::Pass;
        };
        gosper::antidifference(cx, t, k, true).map_or(crate::graph::Outcome::Pass, crate::graph::Outcome::Equal)
    }
}
