#![allow(deprecated)]


use std::convert::AsRef;
use std::fmt::Debug;
use std::fmt::Write;
use std::fmt::{
    self,
};

use num_traits::ToPrimitive;

use super::dag_mgr::DagOp;
use super::expr::Expr;

impl Clone for Expr {
    fn clone(&self) -> Self {

        match self {
            | Self::Constant(c) => Self::Constant(*c),
            | Self::BigInt(i) => Self::BigInt(i.clone()),
            | Self::Rational(r) => Self::Rational(r.clone()),
            | Self::Boolean(b) => Self::Boolean(*b),
            | Self::Variable(s) => Self::Variable(s.clone()),
            | Self::Pattern(s) => Self::Pattern(s.clone()),
            | Self::Add(a, b) => Self::Add(a.clone(), b.clone()),
            | Self::AddList(list) => Self::AddList(list.clone()),
            | Self::Sub(a, b) => Self::Sub(a.clone(), b.clone()),
            | Self::Mul(a, b) => Self::Mul(a.clone(), b.clone()),
            | Self::MulList(list) => Self::MulList(list.clone()),
            | Self::Div(a, b) => Self::Div(a.clone(), b.clone()),
            | Self::Power(a, b) => Self::Power(a.clone(), b.clone()),
            | Self::Sin(a) => Self::Sin(a.clone()),
            | Self::Cos(a) => Self::Cos(a.clone()),
            | Self::Tan(a) => Self::Tan(a.clone()),
            | Self::Exp(a) => Self::Exp(a.clone()),
            | Self::Log(a) => Self::Log(a.clone()),
            | Self::Neg(a) => Self::Neg(a.clone()),
            | Self::Eq(a, b) => Self::Eq(a.clone(), b.clone()),
            | Self::Matrix(m) => Self::Matrix(m.clone()),
            | Self::Vector(v) => Self::Vector(v.clone()),
            | Self::Complex(re, im) => {
                Self::Complex(
                    re.clone(),
                    im.clone(),
                )
            },
            | Self::Derivative(e, s) => Self::Derivative(e.clone(), s.clone()),
            | Self::Sum {
                body,
                var,
                from,
                to,
            } => {
                Self::Sum {
                    body : body.clone(),
                    var : var.clone(),
                    from : from.clone(),
                    to : to.clone(),
                }
            },
            | Self::Integral {
                integrand,
                var,
                lower_bound,
                upper_bound,
            } => {
                Self::Integral {
                    integrand : integrand.clone(),
                    var : var.clone(),
                    lower_bound : lower_bound.clone(),
                    upper_bound : upper_bound.clone(),
                }
            },
            | Self::Path(pt, p1, p2) => {
                Self::Path(
                    pt.clone(),
                    p1.clone(),
                    p2.clone(),
                )
            },
            | Self::Abs(a) => Self::Abs(a.clone()),
            | Self::Sqrt(a) => Self::Sqrt(a.clone()),
            | Self::Sec(a) => Self::Sec(a.clone()),
            | Self::Csc(a) => Self::Csc(a.clone()),
            | Self::Cot(a) => Self::Cot(a.clone()),
            | Self::ArcSin(a) => Self::ArcSin(a.clone()),
            | Self::ArcCos(a) => Self::ArcCos(a.clone()),
            | Self::ArcTan(a) => Self::ArcTan(a.clone()),
            | Self::ArcSec(a) => Self::ArcSec(a.clone()),
            | Self::ArcCsc(a) => Self::ArcCsc(a.clone()),
            | Self::ArcCot(a) => Self::ArcCot(a.clone()),
            | Self::Sinh(a) => Self::Sinh(a.clone()),
            | Self::Cosh(a) => Self::Cosh(a.clone()),
            | Self::Tanh(a) => Self::Tanh(a.clone()),
            | Self::Sech(a) => Self::Sech(a.clone()),
            | Self::Csch(a) => Self::Csch(a.clone()),
            | Self::Coth(a) => Self::Coth(a.clone()),
            | Self::ArcSinh(a) => Self::ArcSinh(a.clone()),
            | Self::ArcCosh(a) => Self::ArcCosh(a.clone()),
            | Self::ArcTanh(a) => Self::ArcTanh(a.clone()),
            | Self::ArcSech(a) => Self::ArcSech(a.clone()),
            | Self::ArcCsch(a) => Self::ArcCsch(a.clone()),
            | Self::ArcCoth(a) => Self::ArcCoth(a.clone()),
            | Self::LogBase(b, a) => Self::LogBase(b.clone(), a.clone()),
            | Self::Atan2(y, x) => Self::Atan2(y.clone(), x.clone()),
            | Self::Binomial(n, k) => Self::Binomial(n.clone(), k.clone()),
            | Self::Boundary(e) => Self::Boundary(e.clone()),
            | Self::Domain(s) => Self::Domain(s.clone()),
            | Self::VolumeIntegral {
                scalar_field,
                volume,
            } => {
                Self::VolumeIntegral {
                    scalar_field : scalar_field.clone(),
                    volume : volume.clone(),
                }
            },
            | Self::SurfaceIntegral {
                vector_field,
                surface,
            } => {
                Self::SurfaceIntegral {
                    vector_field : vector_field.clone(),
                    surface : surface.clone(),
                }
            },
            | Self::Pi => Self::Pi,
            | Self::E => Self::E,
            | Self::Infinity => Self::Infinity,
            | Self::NegativeInfinity => Self::NegativeInfinity,
            | Self::Apply(a, b) => Self::Apply(a.clone(), b.clone()),
            | Self::Tuple(v) => Self::Tuple(v.clone()),
            | Self::Gamma(a) => Self::Gamma(a.clone()),
            | Self::Beta(a, b) => Self::Beta(a.clone(), b.clone()),
            | Self::Erf(a) => Self::Erf(a.clone()),
            | Self::Erfc(a) => Self::Erfc(a.clone()),
            | Self::Erfi(a) => Self::Erfi(a.clone()),
            | Self::Zeta(a) => Self::Zeta(a.clone()),
            | Self::BesselJ(a, b) => Self::BesselJ(a.clone(), b.clone()),
            | Self::BesselY(a, b) => Self::BesselY(a.clone(), b.clone()),
            | Self::LegendreP(a, b) => Self::LegendreP(a.clone(), b.clone()),
            | Self::LaguerreL(a, b) => Self::LaguerreL(a.clone(), b.clone()),
            | Self::HermiteH(a, b) => Self::HermiteH(a.clone(), b.clone()),
            | Self::Digamma(a) => Self::Digamma(a.clone()),
            | Self::KroneckerDelta(a, b) => Self::KroneckerDelta(a.clone(), b.clone()),
            | Self::DerivativeN(e, s, n) => {
                Self::DerivativeN(
                    e.clone(),
                    s.clone(),
                    n.clone(),
                )
            },
            | Self::Series(a, b, c, d) => {
                Self::Series(
                    a.clone(),
                    b.clone(),
                    c.clone(),
                    d.clone(),
                )
            },
            | Self::Summation(a, b, c, d) => {
                Self::Summation(
                    a.clone(),
                    b.clone(),
                    c.clone(),
                    d.clone(),
                )
            },
            | Self::Product(a, b, c, d) => {
                Self::Product(
                    a.clone(),
                    b.clone(),
                    c.clone(),
                    d.clone(),
                )
            },
            | Self::ConvergenceAnalysis(e, s) => Self::ConvergenceAnalysis(e.clone(), s.clone()),
            | Self::AsymptoticExpansion(a, b, c, d) => {
                Self::AsymptoticExpansion(
                    a.clone(),
                    b.clone(),
                    c.clone(),
                    d.clone(),
                )
            },
            | Self::Lt(a, b) => Self::Lt(a.clone(), b.clone()),
            | Self::Gt(a, b) => Self::Gt(a.clone(), b.clone()),
            | Self::Le(a, b) => Self::Le(a.clone(), b.clone()),
            | Self::Ge(a, b) => Self::Ge(a.clone(), b.clone()),
            | Self::Union(v) => Self::Union(v.clone()),
            | Self::Interval(a, b, c, d) => {
                Self::Interval(
                    a.clone(),
                    b.clone(),
                    *c,
                    *d,
                )
            },
            | Self::Solve(e, s) => Self::Solve(e.clone(), s.clone()),
            | Self::Substitute(a, b, c) => {
                Self::Substitute(
                    a.clone(),
                    b.clone(),
                    c.clone(),
                )
            },
            | Self::Limit(a, b, c) => {
                Self::Limit(
                    a.clone(),
                    b.clone(),
                    c.clone(),
                )
            },
            | Self::InfiniteSolutions => Self::InfiniteSolutions,
            | Self::NoSolution => Self::NoSolution,
            | Self::Dag(n) => Self::Dag(n.clone()),
            | Self::Factorial(a) => Self::Factorial(a.clone()),
            | Self::Permutation(a, b) => Self::Permutation(a.clone(), b.clone()),
            | Self::Combination(a, b) => Self::Combination(a.clone(), b.clone()),
            | Self::FallingFactorial(a, b) => Self::FallingFactorial(a.clone(), b.clone()),
            | Self::RisingFactorial(a, b) => Self::RisingFactorial(a.clone(), b.clone()),
            | Self::Ode {
                equation,
                func,
                var,
            } => {
                Self::Ode {
                    equation : equation.clone(),
                    func : func.clone(),
                    var : var.clone(),
                }
            },
            | Self::Pde {
                equation,
                func,
                vars,
            } => {
                Self::Pde {
                    equation : equation.clone(),
                    func : func.clone(),
                    vars : vars.clone(),
                }
            },
            | Self::GeneralSolution(e) => Self::GeneralSolution(e.clone()),
            | Self::ParticularSolution(e) => Self::ParticularSolution(e.clone()),
            | Self::Fredholm(a, b, c, d) => {
                Self::Fredholm(
                    a.clone(),
                    b.clone(),
                    c.clone(),
                    d.clone(),
                )
            },
            | Self::Volterra(a, b, c, d) => {
                Self::Volterra(
                    a.clone(),
                    b.clone(),
                    c.clone(),
                    d.clone(),
                )
            },
            | Self::And(v) => Self::And(v.clone()),
            | Self::Or(v) => Self::Or(v.clone()),
            | Self::Not(a) => Self::Not(a.clone()),
            | Self::Xor(a, b) => Self::Xor(a.clone(), b.clone()),
            | Self::Implies(a, b) => Self::Implies(a.clone(), b.clone()),
            | Self::Equivalent(a, b) => Self::Equivalent(a.clone(), b.clone()),
            | Self::Predicate {
                name,
                args,
            } => {
                Self::Predicate {
                    name : name.clone(),
                    args : args.clone(),
                }
            },
            | Self::ForAll(s, e) => Self::ForAll(s.clone(), e.clone()),
            | Self::Exists(s, e) => Self::Exists(s.clone(), e.clone()),
            | Self::Polynomial(c) => Self::Polynomial(c.clone()),
            | Self::SparsePolynomial(p) => Self::SparsePolynomial(p.clone()),
            | Self::Floor(a) => Self::Floor(a.clone()),
            | Self::IsPrime(a) => Self::IsPrime(a.clone()),
            | Self::Gcd(a, b) => Self::Gcd(a.clone(), b.clone()),
            | Self::Distribution(d) => Self::Distribution(d.clone()),
            | Self::Mod(a, b) => Self::Mod(a.clone(), b.clone()),
            | Self::Max(a, b) => Self::Max(a.clone(), b.clone()),
            | Self::Quantity(q) => Self::Quantity(q.clone()),
            | Self::QuantityWithValue(v, u) => Self::QuantityWithValue(v.clone(), u.clone()),
            | Self::Transpose(a) => Self::Transpose(a.clone()),
            | Self::MatrixMul(a, b) => Self::MatrixMul(a.clone(), b.clone()),
            | Self::MatrixVecMul(a, b) => Self::MatrixVecMul(a.clone(), b.clone()),
            | Self::Inverse(a) => Self::Inverse(a.clone()),
            | Self::System(v) => Self::System(v.clone()),
            | Self::Solutions(v) => Self::Solutions(v.clone()),
            | Self::ParametricSolution {
                x,
                y,
            } => {
                Self::ParametricSolution {
                    x : x.clone(),
                    y : y.clone(),
                }
            },
            | Self::RootOf {
                poly,
                index,
            } => {
                Self::RootOf {
                    poly : poly.clone(),
                    index : *index,
                }
            },

            | Self::CustomZero => Self::CustomZero,
            | Self::CustomString(a) => Self::CustomString(a.clone()),
            | Self::CustomArcOne(a) => Self::CustomArcOne(a.clone()),
            | Self::CustomArcTwo(a, b) => Self::CustomArcTwo(a.clone(), b.clone()),
            | Self::CustomArcThree(a, b, c) => {
                Self::CustomArcThree(
                    a.clone(),
                    b.clone(),
                    c.clone(),
                )
            },
            | Self::CustomArcFour(a, b, c, d) => {
                Self::CustomArcFour(
                    a.clone(),
                    b.clone(),
                    c.clone(),
                    d.clone(),
                )
            },
            | Self::CustomArcFive(a, b, c, d, e) => {
                Self::CustomArcFive(
                    a.clone(),
                    b.clone(),
                    c.clone(),
                    d.clone(),
                    e.clone(),
                )
            },
            | Self::CustomVecOne(a) => Self::CustomVecOne(a.clone()),
            | Self::CustomVecTwo(a, b) => Self::CustomVecTwo(a.clone(), b.clone()),
            | Self::CustomVecThree(a, b, c) => {
                Self::CustomVecThree(
                    a.clone(),
                    b.clone(),
                    c.clone(),
                )
            },
            | Self::CustomVecFour(a, b, c, d) => {
                Self::CustomVecFour(
                    a.clone(),
                    b.clone(),
                    c.clone(),
                    d.clone(),
                )
            },
            | Self::CustomVecFive(a, b, c, d, e) => {
                Self::CustomVecFive(
                    a.clone(),
                    b.clone(),
                    c.clone(),
                    d.clone(),
                    e.clone(),
                )
            },
            | Self::UnaryList(s, a) => Self::UnaryList(s.clone(), a.clone()),
            | Self::BinaryList(s, a, b) => {
                Self::BinaryList(
                    s.clone(),
                    a.clone(),
                    b.clone(),
                )
            },
            | Self::NaryList(s, v) => Self::NaryList(s.clone(), v.clone()),
        }
    }
}

impl Debug for Expr {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {

        // Use Display for a more compact representation in debug outputs
        write!(f, "{self}")
    }
}

impl fmt::Display for Expr {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {

        match self {
            | Self::Dag(node) => {
                match node.to_expr() {
                    | Ok(expr) => write!(f, "{expr}"),
                    | Err(e) => {

                        write!(
                            f,
                            "<Error converting DAG to Expr: {e}>"
                        )
                    },
                }
            },
            | Self::Constant(c) => write!(f, "{c}"),
            | Self::BigInt(i) => write!(f, "{i}"),
            | Self::Rational(r) => write!(f, "{r}"),
            | Self::Boolean(b) => write!(f, "{b}"),
            | Self::Variable(s) => write!(f, "{s}"),
            | Self::Pattern(s) => write!(f, "{s}"),
            | Self::Add(a, b) => write!(f, "({a} + {b})"),
            | Self::AddList(list) => {

                write!(f, "(")?;

                for (i, item) in list
                    .iter()
                    .enumerate()
                {

                    if i > 0 {

                        write!(f, " + ")?;
                    }

                    write!(f, "{item}")?;
                }

                write!(f, ")")
            },
            | Self::Sub(a, b) => write!(f, "({a} - {b})"),
            | Self::Mul(a, b) => write!(f, "({a} * {b})"),
            | Self::MulList(list) => {

                write!(f, "(")?;

                for (i, item) in list
                    .iter()
                    .enumerate()
                {

                    if i > 0 {

                        write!(f, " * ")?;
                    }

                    write!(f, "{item}")?;
                }

                write!(f, ")")
            },
            | Self::Div(a, b) => write!(f, "({a} / {b})"),
            | Self::Power(a, b) => write!(f, "({a}^({b}))"),
            | Self::Sin(a) => write!(f, "sin({a})"),
            | Self::Cos(a) => write!(f, "cos({a})"),
            | Self::Tan(a) => write!(f, "tan({a})"),
            | Self::Exp(a) => write!(f, "exp({a})"),
            | Self::Log(a) => write!(f, "ln({a})"),
            | Self::Neg(a) => write!(f, "-({a})"),
            | Self::Eq(a, b) => write!(f, "{a} = {b}"),
            | Self::Matrix(m) => {

                let mut s = String::new();

                s.push('[');

                for (i, row) in m.iter().enumerate() {

                    if i > 0 {

                        s.push_str("; ");
                    }

                    s.push('[');

                    for (j, val) in row
                        .iter()
                        .enumerate()
                    {

                        if j > 0 {

                            s.push_str(", ");
                        }

                        s.push_str(&val.to_string());
                    }

                    s.push(']');
                }

                s.push(']');

                write!(f, "{s}")
            },
            | Self::Vector(v) => {

                let mut s = String::new();

                s.push('[');

                for (i, val) in v.iter().enumerate() {

                    if i > 0 {

                        s.push_str(", ");
                    }

                    s.push_str(&val.to_string());
                }

                s.push(']');

                write!(f, "{s}")
            },
            | Self::Complex(re, im) => write!(f, "({re} + {im}i)"),
            | Self::Derivative(expr, var) => {

                write!(
                    f,
                    "d/d{var}({expr})"
                )
            },
            | Self::Integral {
                integrand,
                var,
                lower_bound,
                upper_bound,
            } => {

                write!(
                    f,
                    "integral({integrand}, {var}, {lower_bound}, {upper_bound})"
                )
            },
            | Self::Sum {
                body,
                var,
                from,
                to,
            } => {

                write!(
                    f,
                    "sum({body}, {var}, {from}, {to})"
                )
            },
            | Self::Path(path_type, p1, p2) => {

                write!(
                    f,
                    "path({path_type:?}, {p1}, {p2})"
                )
            },
            | Self::Abs(a) => write!(f, "|{a}|"),
            | Self::Sqrt(a) => write!(f, "sqrt({a})"),
            | Self::Sec(a) => write!(f, "sec({a})"),
            | Self::Csc(a) => write!(f, "csc({a})"),
            | Self::Cot(a) => write!(f, "cot({a})"),
            | Self::ArcSin(a) => write!(f, "asin({a})"),
            | Self::ArcCos(a) => write!(f, "acos({a})"),
            | Self::ArcTan(a) => write!(f, "atan({a})"),
            | Self::ArcSec(a) => write!(f, "asec({a})"),
            | Self::ArcCsc(a) => write!(f, "acsc({a})"),
            | Self::ArcCot(a) => write!(f, "acot({a})"),
            | Self::Sinh(a) => write!(f, "sinh({a})"),
            | Self::Cosh(a) => write!(f, "cosh({a})"),
            | Self::Tanh(a) => write!(f, "tanh({a})"),
            | Self::Sech(a) => write!(f, "sech({a})"),
            | Self::Csch(a) => write!(f, "csch({a})"),
            | Self::Coth(a) => write!(f, "coth({a})"),
            | Self::ArcSinh(a) => write!(f, "asinh({a})"),
            | Self::ArcCosh(a) => write!(f, "acosh({a})"),
            | Self::ArcTanh(a) => write!(f, "atanh({a})"),
            | Self::ArcSech(a) => write!(f, "asech({a})"),
            | Self::ArcCsch(a) => write!(f, "acsch({a})"),
            | Self::ArcCoth(a) => write!(f, "acoth({a})"),
            | Self::LogBase(b, a) => {

                write!(
                    f,
                    "log_base({b}, {a})"
                )
            },
            | Self::Atan2(y, x) => write!(f, "atan2({y}, {x})"),
            | Self::Pi => write!(f, "Pi"),
            | Self::E => write!(f, "E"),
            | Self::Infinity => write!(f, "Infinity"),
            | Self::NegativeInfinity => write!(f, "-Infinity"),
            | Self::Ode {
                equation,
                func,
                var,
            } => {

                write!(
                    f,
                    "ode({equation}, {func}, {var})"
                )
            },
            | Self::Pde {
                equation,
                func,
                vars,
            } => {

                write!(
                    f,
                    "pde({equation}, {func}, {vars:?})"
                )
            },
            | Self::Fredholm(a, b, c, d) => {

                write!(
                    f,
                    "fredholm({a}, {b}, {c}, {d})"
                )
            },
            | Self::Volterra(a, b, c, d) => {

                write!(
                    f,
                    "volterra({a}, {b}, {c}, {d})"
                )
            },
            | Self::And(v) => {

                write!(
                    f,
                    "({})",
                    v.iter()
                        .map(std::string::ToString::to_string)
                        .collect::<Vec<String>>()
                        .join(" && ")
                )
            },
            | Self::Or(v) => {

                write!(
                    f,
                    "({})",
                    v.iter()
                        .map(std::string::ToString::to_string)
                        .collect::<Vec<String>>()
                        .join(" || ")
                )
            },
            | Self::Not(a) => write!(f, "!({a})"),
            | Self::Xor(a, b) => write!(f, "({a} ^ {b})"),
            | Self::Implies(a, b) => write!(f, "({a} => {b})"),
            | Self::Equivalent(a, b) => write!(f, "({a} <=> {b})"),
            | Self::Predicate {
                name,
                args,
            } => {

                let args_str = args
                    .iter()
                    .map(std::string::ToString::to_string)
                    .collect::<Vec<_>>()
                    .join(", ");

                write!(
                    f,
                    "{name}({args_str})"
                )
            },
            | Self::ForAll(s, e) => write!(f, "∀{s}. ({e})"),
            | Self::Exists(s, e) => write!(f, "∃{s}. ({e})"),
            | Self::Polynomial(coeffs) => {

                let mut s = String::new();

                for (i, coeff) in coeffs
                    .iter()
                    .enumerate()
                    .rev()
                {

                    if !s.is_empty() {

                        s.push_str(" + ");
                    }

                    let _ = write!(s, "{coeff}*x^{i}");
                }

                write!(f, "{s}")
            },
            | Self::SparsePolynomial(p) => {

                let mut s = String::new();

                for (monomial, coeff) in &p.terms {

                    if !s.is_empty() {

                        s.push_str(" + ");
                    }

                    let _ = write!(s, "({coeff})");

                    for (var, exp) in &monomial.0 {

                        let _ = write!(s, "*{var}^{exp}");
                    }
                }

                write!(f, "{s}")
            },
            | Self::Floor(a) => write!(f, "floor({a})"),
            | Self::IsPrime(a) => write!(f, "is_prime({a})"),
            | Self::Gcd(a, b) => write!(f, "gcd({a}, {b})"),
            | Self::Factorial(a) => write!(f, "factorial({a})"),
            | Self::Distribution(d) => write!(f, "{d:?}"),
            | Self::Mod(a, b) => write!(f, "({a} mod {b})"),
            | Self::Max(a, b) => write!(f, "max({a}, {b})"),
            | Self::System(v) => write!(f, "system({v:?})"),
            | Self::Solutions(v) => {

                write!(
                    f,
                    "solutions({v:?})"
                )
            },
            | Self::ParametricSolution {
                x,
                y,
            } => {

                write!(
                    f,
                    "parametric_solution({x}, {y})"
                )
            },
            | Self::RootOf {
                poly,
                index,
            } => {

                write!(
                    f,
                    "root_of({poly}, {index})"
                )
            },
            | Self::Erfc(a) => write!(f, "erfc({a})"),
            | Self::Erfi(a) => write!(f, "erfi({a})"),
            | Self::Zeta(a) => write!(f, "zeta({a})"),

            | Self::CustomZero => write!(f, "CustomZero"),
            | Self::CustomString(s) => {

                write!(
                    f,
                    "CustomString({s})"
                )
            },
            | Self::CustomArcOne(a) => {

                write!(
                    f,
                    "CustomArcOne({a})"
                )
            },
            | Self::CustomArcTwo(a, b) => {

                write!(
                    f,
                    "CustomArcTwo({a}, {b})"
                )
            },
            | Self::CustomArcThree(a, b, c) => {

                write!(
                    f,
                    "CustomArcThree({a}, {b}, {c})"
                )
            },
            | Self::CustomArcFour(a, b, c, d) => {

                write!(
                    f,
                    "CustomArcFour({a}, {b}, {c}, {d})"
                )
            },
            | Self::CustomArcFive(a, b, c, d, e) => {

                write!(
                    f,
                    "CustomArcFive({a}, {b}, {c}, {d}, {e})"
                )
            },
            | Self::CustomVecOne(v) => {

                write!(
                    f,
                    "CustomVecOne({v:?})"
                )
            },
            | Self::CustomVecTwo(v1, v2) => {

                write!(
                    f,
                    "CustomVecTwo({v1:?}, {v2:?})"
                )
            },
            | Self::CustomVecThree(v1, v2, v3) => {

                write!(
                    f,
                    "CustomVecThree({v1:?}, {v2:?}, {v3:?})"
                )
            },
            | Self::CustomVecFour(v1, v2, v3, v4) => {

                write!(
                    f,
                    "CustomVecFour({v1:?}, {v2:?}, {v3:?}, {v4:?})"
                )
            },
            | Self::CustomVecFive(v1, v2, v3, v4, v5) => {

                write!(
                    f,
                    "CustomVecFive({v1:?}, {v2:?}, {v3:?}, {v4:?}, {v5:?})"
                )
            },

            | Self::UnaryList(s, a) => write!(f, "{s}({a})"),
            | Self::BinaryList(s, a, b) => write!(f, "{s}({a}, {b})"),
            | Self::NaryList(s, v) => {

                write!(f, "{s}(")?;

                for (i, item) in v.iter().enumerate() {

                    if i > 0 {

                        write!(f, ", ")?;
                    }

                    write!(f, "{item}")?;
                }

                write!(f, ")")
            },

            | Self::GeneralSolution(e) => {

                write!(
                    f,
                    "general_solution({e})"
                )
            },
            | Self::ParticularSolution(e) => {

                write!(
                    f,
                    "particular_solution({e})"
                )
            },
            | Self::InfiniteSolutions => {

                write!(
                    f,
                    "InfiniteSolutions"
                )
            },
            | Self::NoSolution => write!(f, "NoSolution"),
            | Self::KroneckerDelta(a, b) => {

                write!(
                    f,
                    "kronecker_delta({a}, {b})"
                )
            },
            | Self::MatrixMul(a, b) => {

                write!(
                    f,
                    "matrix_mul({a}, {b})"
                )
            },
            | Self::MatrixVecMul(a, b) => {

                write!(
                    f,
                    "matrix_vec_mul({a}, {b})"
                )
            },
            | Self::QuantityWithValue(v, u) => {

                write!(
                    f,
                    "quantity_with_value({v}, \"{u}\")"
                )
            },
            | Self::Tuple(v) => write!(f, "tuple({v:?})"),
            | Self::Interval(a, b, c, d) => {

                write!(
                    f,
                    "interval({a}, {b}, {c}, {d})"
                )
            },
            | Self::Domain(s) => write!(f, "domain({s})"),
            | Self::AsymptoticExpansion(a, b, c, d) => {

                write!(
                    f,
                    "asymptotic_expansion({a}, {b}, {c}, {d})"
                )
            },
            | Self::ConvergenceAnalysis(e, s) => {

                write!(
                    f,
                    "convergence_analysis({e}, {s})"
                )
            },
            | Self::DerivativeN(e, s, n) => {

                write!(
                    f,
                    "derivative_n({e}, {s}, {n})"
                )
            },
            | Self::FallingFactorial(a, b) => {

                write!(
                    f,
                    "falling_factorial({a}, {b})"
                )
            },
            | Self::RisingFactorial(a, b) => {

                write!(
                    f,
                    "rising_factorial({a}, {b})"
                )
            },
            | Self::Product(a, b, c, d) => {

                write!(
                    f,
                    "product({a}, {b}, {c}, {d})"
                )
            },
            | Self::Series(a, b, c, d) => {

                write!(
                    f,
                    "series({a}, {b}, {c}, {d})"
                )
            },
            | Self::Summation(a, b, c, d) => {

                write!(
                    f,
                    "summation({a}, {b}, {c}, {d})"
                )
            },
            | Self::Lt(a, b) => write!(f, "({a} < {b})"),
            | Self::Gt(a, b) => write!(f, "({a} > {b})"),
            | Self::Le(a, b) => write!(f, "({a} <= {b})"),
            | Self::Ge(a, b) => write!(f, "({a} >= {b})"),
            | Self::Transpose(a) => write!(f, "transpose({a})"),
            | Self::VolumeIntegral {
                scalar_field,
                volume,
            } => {

                write!(
                    f,
                    "volume_integral({scalar_field}, {volume})"
                )
            },
            | Self::SurfaceIntegral {
                vector_field,
                surface,
            } => {

                write!(
                    f,
                    "surface_integral({vector_field}, {surface})"
                )
            },
            | Self::Union(v) => write!(f, "union({v:?})"),
            | Self::Solve(e, s) => write!(f, "solve({e}, {s})"),
            | Self::Apply(a, b) => write!(f, "apply({a}, {b})"),
            | Self::Quantity(q) => write!(f, "{q:?}"),
            | Self::Inverse(a) => write!(f, "inverse({a})"),
            | Self::Limit(a, b, c) => {

                write!(
                    f,
                    "limit({a}, {b}, {c})"
                )
            },
            | Self::Binomial(a, b) => {

                write!(
                    f,
                    "binomial({a}, {b})"
                )
            },
            | Self::Permutation(a, b) => {

                write!(
                    f,
                    "permutation({a}, {b})"
                )
            },
            | Self::Combination(a, b) => {

                write!(
                    f,
                    "combination({a}, {b})"
                )
            },
            | Self::Boundary(a) => write!(f, "boundary({a})"),
            | Self::Gamma(a) => write!(f, "gamma({a})"),
            | Self::Beta(a, b) => write!(f, "beta({a}, {b})"),
            | Self::Erf(a) => write!(f, "erf({a})"),
            | Self::BesselJ(a, b) => {

                write!(
                    f,
                    "BesselJ({a}, {b})"
                )
            },
            | Self::BesselY(a, b) => {

                write!(
                    f,
                    "BesselY({a}, {b})"
                )
            },
            | Self::LegendreP(a, b) => {

                write!(
                    f,
                    "LegendreP({a}, {b})"
                )
            },
            | Self::LaguerreL(a, b) => {

                write!(
                    f,
                    "LaguerreL({a}, {b})"
                )
            },
            | Self::HermiteH(a, b) => {

                write!(
                    f,
                    "HermiteH({a}, {b})"
                )
            },
            | Self::Digamma(a) => write!(f, "Digamma({a})"),
            | Self::Substitute(a, b, c) => {

                write!(
                    f,
                    "substitute({a}, {b}, {c})"
                )
            },
        }
    }
}

impl Expr {
    /// Returns the real part of the expression.
    #[must_use]
    #[inline]

    pub fn re(&self) -> Self {

        if let Self::Complex(re, _) =
            self
        {

            re.as_ref().clone()
        } else {

            self.clone()
        }
    }

    /// Returns the imaginary part of the expression.
    #[must_use]
    #[inline]

    pub fn im(&self) -> Self {

        if let Self::Complex(_, im) =
            self
        {

            im.as_ref().clone()
        } else {

            Self::Constant(0.0)
        }
    }

    /// Attempts to convert the expression to a 64-bit float.
    #[inline]
    #[must_use]

    pub fn to_f64(
        &self
    ) -> Option<f64> {

        match self {
            | Self::Constant(val) => {
                Some(*val)
            },
            | Self::BigInt(val) => {
                val.to_f64()
            },
            | Self::Rational(val) => {
                val.to_f64()
            },
            | Self::Pi => Some(
                std::f64::consts::PI,
            ),
            | Self::E => {
                Some(
                    std::f64::consts::E,
                )
            },
            | Self::Dag(node) => {
                node.to_expr()
                    .ok()?
                    .to_f64()
            },
            | _ => None,
        }
    }

    /// Returns the operation type of the expression in a unified way.
    ///
    /// This method handles both regular expressions and DAG nodes, returning
    /// the operation regardless of internal representation.
    ///
    /// # Returns
    /// * `DagOp` - The operation type corresponding to this expression
    #[must_use]

    pub fn op(&self) -> DagOp {

        match self {
            | Self::Dag(node) => node.op.clone(),
            | _ => {
                self.to_dag_op_internal()
                    .expect(
                        "Failed to convert Expr to DagOp; this should be impossible for any valid \
                         Expr",
                    )
            },
        }
    }

    /// Returns the children of the expression in a unified way.
    ///
    /// This method handles both regular expressions and DAG nodes, returning
    /// the direct child expressions regardless of internal representation.
    ///
    /// # Returns
    /// * `Vec<Expr>` - A vector containing the direct children of this expression
    #[must_use]

    pub fn children(
        &self
    ) -> Vec<Self> {

        match self {
            | Self::Dag(node) => {
                node.children
                    .iter()
                    .map(|n| {

                        Self::Dag(
                            n.clone(),
                        )
                    })
                    .collect()
            },
            | _ => self
                .get_children_internal(
                ),
        }
    }

    #[allow(dead_code)]

    pub(crate) const fn variant_order(
        &self
    ) -> i32 {

        match self {
            | Self::Constant(_) => 0,
            | Self::BigInt(_) => 1,
            | Self::Rational(_) => 2,
            | Self::Boolean(_) => 3,
            | Self::Variable(_) => 4,
            | Self::Pattern(_) => 5,
            | Self::Add(_, _) => 6,
            | Self::AddList(_) => 6, // Same order as Add
            | Self::Sub(_, _) => 7,
            | Self::Mul(_, _) => 8,
            | Self::MulList(_) => 8, // Same order as Mul
            | Self::Div(_, _) => 9,
            | Self::Power(_, _) => 10,
            | Self::Sin(_) => 11,
            | Self::Cos(_) => 12,
            | Self::Tan(_) => 13,
            | Self::Exp(_) => 14,
            | Self::Log(_) => 15,
            | Self::Neg(_) => 16,
            | Self::Eq(_, _) => 17,
            | Self::Matrix(_) => 18,
            | Self::Vector(_) => 19,
            | Self::Complex(_, _) => 20,
            | Self::Derivative(_, _) => 21,
            | Self::Integral {
                ..
            } => 22,
            | Self::Sum {
                ..
            } => 22, // Assign same order as Integral for now
            | Self::Path(_, _, _) => 23,
            | Self::Abs(_) => 24,
            | Self::Sqrt(_) => 25,
            | Self::Sec(_) => 26,
            | Self::Csc(_) => 27,
            | Self::Cot(_) => 28,
            | Self::ArcSin(_) => 29,
            | Self::ArcCos(_) => 30,
            | Self::ArcTan(_) => 31,
            | Self::ArcSec(_) => 32,
            | Self::ArcCsc(_) => 33,
            | Self::ArcCot(_) => 34,
            | Self::Sinh(_) => 35,
            | Self::Cosh(_) => 36,
            | Self::Tanh(_) => 37,
            | Self::Sech(_) => 38,
            | Self::Csch(_) => 39,
            | Self::Coth(_) => 40,
            | Self::ArcSinh(_) => 41,
            | Self::ArcCosh(_) => 42,
            | Self::ArcTanh(_) => 43,
            | Self::ArcSech(_) => 44,
            | Self::ArcCsch(_) => 45,
            | Self::ArcCoth(_) => 46,
            | Self::LogBase(_, _) => 47,
            | Self::Atan2(_, _) => 48,
            | Self::Binomial(_, _) => 49,
            | Self::Boundary(_) => 50,
            | Self::Domain(_) => 51,
            | Self::VolumeIntegral {
                ..
            } => 52,
            | Self::SurfaceIntegral {
                ..
            } => 53,
            | Self::Pi => 54,
            | Self::E => 55,
            | Self::Infinity => 56,
            | Self::NegativeInfinity => 57,
            | Self::Apply(_, _) => 58,
            | Self::Tuple(_) => 59,
            | Self::Gamma(_) => 60,
            | Self::Beta(_, _) => 61,
            | Self::Erf(_) => 62,
            | Self::Erfc(_) => 63,
            | Self::Erfi(_) => 64,
            | Self::Zeta(_) => 65,
            | Self::BesselJ(_, _) => 66,
            | Self::BesselY(_, _) => 67,
            | Self::LegendreP(_, _) => 68,
            | Self::LaguerreL(_, _) => 69,
            | Self::HermiteH(_, _) => 70,
            | Self::Digamma(_) => 71,
            | Self::KroneckerDelta(_, _) => 72,
            | Self::DerivativeN(_, _, _) => 73,
            | Self::Series(_, _, _, _) => 74,
            | Self::Summation(_, _, _, _) => 75,
            | Self::Product(_, _, _, _) => 76,
            | Self::ConvergenceAnalysis(_, _) => 77,
            | Self::AsymptoticExpansion(_, _, _, _) => 78,
            | Self::Lt(_, _) => 79,
            | Self::Gt(_, _) => 80,
            | Self::Le(_, _) => 81,
            | Self::Ge(_, _) => 82,
            | Self::Union(_) => 83,
            | Self::Interval(_, _, _, _) => 84,
            | Self::Solve(_, _) => 85,
            | Self::Substitute(_, _, _) => 86,
            | Self::Limit(_, _, _) => 87,
            | Self::InfiniteSolutions => 88,
            | Self::NoSolution => 89,
            | Self::Dag(_) => 90,
            | Self::Factorial(_) => 91,
            | Self::Permutation(_, _) => 92,
            | Self::Combination(_, _) => 93,
            | Self::FallingFactorial(_, _) => 94,
            | Self::RisingFactorial(_, _) => 95,
            | Self::Ode {
                ..
            } => 96,
            | Self::Pde {
                ..
            } => 97,
            | Self::GeneralSolution(_) => 98,
            | Self::ParticularSolution(_) => 99,
            | Self::Fredholm(_, _, _, _) => 100,
            | Self::Volterra(_, _, _, _) => 101,
            | Self::And(_) => 102,
            | Self::Or(_) => 103,
            | Self::Not(_) => 104,
            | Self::Xor(_, _) => 105,
            | Self::Implies(_, _) => 106,
            | Self::Equivalent(_, _) => 107,
            | Self::Predicate {
                ..
            } => 108,
            | Self::ForAll(_, _) => 109,
            | Self::Exists(_, _) => 110,
            | Self::Polynomial(_) => 111,
            | Self::SparsePolynomial(_) => 112,
            | Self::Floor(_) => 113,
            | Self::IsPrime(_) => 114,
            | Self::Gcd(_, _) => 115,
            | Self::Distribution(_) => 116,
            | Self::Mod(_, _) => 117,
            | Self::Max(_, _) => 118,
            | Self::Transpose(_) => 119,
            | Self::MatrixMul(_, _) => 120,
            | Self::MatrixVecMul(_, _) => 121,
            | Self::Inverse(_) => 122,
            | Self::System(_) => 123,
            | Self::Solutions(_) => 124,
            | Self::ParametricSolution {
                ..
            } => 125,
            | Self::RootOf {
                ..
            } => 126,
            | Self::Quantity(_) => 127,
            | Self::QuantityWithValue(_, _) => 128,
            | Self::CustomZero => 129,
            | Self::CustomString(_) => 130,
            | Self::CustomArcOne(_) => 131,
            | Self::CustomArcTwo(_, _) => 132,
            | Self::CustomArcThree(_, _, _) => 133,
            | Self::CustomArcFour(_, _, _, _) => 134,
            | Self::CustomArcFive(_, _, _, _, _) => 135,
            | Self::CustomVecOne(_) => 136,
            | Self::CustomVecTwo(_, _) => 137,
            | Self::CustomVecThree(_, _, _) => 138,
            | Self::CustomVecFour(_, _, _, _) => 139,
            | Self::CustomVecFive(_, _, _, _, _) => 140,
            | Self::UnaryList(_, _) => 141,
            | Self::BinaryList(_, _, _) => 142,
            | Self::NaryList(_, _) => 143,
        }
    }
}
