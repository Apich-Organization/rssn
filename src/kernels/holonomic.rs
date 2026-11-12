//! D-finite (holonomic) functions of one variable.
//!
//! A function is *D-finite* when it satisfies a linear differential equation
//! with polynomial coefficients. The class contains the rational functions,
//! `exp`, `log`, the circular and hyperbolic functions and their inverses,
//! algebraic powers, and it is closed under sum, product and derivative.
//! That closure is effective, which gives three things a rewriting engine
//! cannot offer:
//!
//! * [`annihilator`]: a differential equation for any expression built from
//!   those pieces, computed by linear algebra over `Q(x)` rather than by
//!   pattern matching;
//! * [`recurrence`]: the linear recurrence its Taylor coefficients satisfy;
//! * [`is_zero`]: a **decision procedure** for identities. An expression is
//!   identically zero exactly when finitely many of its Taylor coefficients
//!   vanish, and the equation says how many: every coefficient that the
//!   recurrence does not determine from earlier ones. `sin(x)^2 + cos(x)^2 -
//!   1` is recognised as zero by computing three rational numbers, with no
//!   trigonometric identity involved.
//!
//! Expressions are [`Expr`] trees over atoms that are analytic at the origin
//! with rational Taylor coefficients (the conditions are stated on each
//! [`Atom`]); the `holonomic` rule set converts terms to them and leaves
//! anything else unreduced.
//!
//! Everything is exact rational arithmetic. Sizes are bounded by
//! [`Limits`]; exceeding one returns `None`, never a wrong answer.

#![allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::many_single_char_names,
    clippy::module_name_repetitions,
    clippy::single_call_fn,
    clippy::missing_docs_in_private_items,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    reason = "dense polynomial arithmetic over vectors whose lengths are computed alongside"
)]

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

type Q = BigRational;

fn qi(n: i64) -> Q {
    Q::from(BigInt::from(n))
}

/// A polynomial over the rationals, lowest degree first, without trailing
/// zeros (the zero polynomial is empty).
pub type Poly = Vec<Q>;

fn trim(mut p: Poly) -> Poly {
    while p.last().is_some_and(Zero::is_zero) {
        p.pop();
    }
    p
}

fn p_add(
    a: &[Q],
    b: &[Q],
) -> Poly {
    let n = a.len().max(b.len());
    trim(
        (0..n)
            .map(|k| {
                a.get(k).cloned().unwrap_or_else(Q::zero)
                    + b.get(k).cloned().unwrap_or_else(Q::zero)
            })
            .collect(),
    )
}

fn p_neg(a: &[Q]) -> Poly {
    a.iter().map(|c| -c).collect()
}

fn p_sub(
    a: &[Q],
    b: &[Q],
) -> Poly {
    p_add(a, &p_neg(b))
}

fn p_mul(
    a: &[Q],
    b: &[Q],
) -> Poly {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![Q::zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        if x.is_zero() {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    trim(out)
}

fn p_scale(
    a: &[Q],
    c: &Q,
) -> Poly {
    trim(a.iter().map(|x| x * c).collect())
}

fn p_deriv(a: &[Q]) -> Poly {
    trim(
        a.iter()
            .enumerate()
            .skip(1)
            .map(|(k, c)| c * qi(k as i64))
            .collect(),
    )
}

fn p_divrem(
    a: &[Q],
    b: &[Q],
) -> (Poly, Poly) {
    let mut r = a.to_vec();
    let Some(lead) = b.last() else {
        return (Vec::new(), r);
    };
    if r.len() < b.len() {
        return (Vec::new(), trim(r));
    }
    let mut q = vec![Q::zero(); r.len() - b.len() + 1];
    for k in (0..q.len()).rev() {
        let factor = &r[k + b.len() - 1] / lead;
        if !factor.is_zero() {
            for (j, c) in b.iter().enumerate() {
                r[k + j] -= &factor * c;
            }
        }
        q[k] = factor;
    }
    r.truncate(b.len() - 1);
    (trim(q), trim(r))
}

/// The monic greatest common divisor (the zero polynomial for two zeros).
fn p_gcd(
    a: &[Q],
    b: &[Q],
) -> Poly {
    let (mut x, mut y) = (a.to_vec(), b.to_vec());
    while !y.is_empty() {
        let (_, r) = p_divrem(&x, &y);
        x = y;
        y = r;
    }
    match x.last().cloned() {
        | Some(lead) => x.iter().map(|c| c / &lead).collect(),
        | None => x,
    }
}

fn p_eval(
    a: &[Q],
    x: &Q,
) -> Q {
    a.iter().rev().fold(Q::zero(), |acc, c| acc * x + c)
}

/// A rational function in lowest terms with a monic denominator.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RatFn {
    /// Numerator.
    pub num: Poly,
    /// Denominator (monic, never zero).
    pub den: Poly,
}

impl RatFn {
    /// `num / den` in lowest terms, or `None` for a zero denominator.
    #[must_use]
    pub fn new(
        num: Poly,
        den: Poly,
    ) -> Option<Self> {
        let (num, den) = (trim(num), trim(den));
        if den.is_empty() {
            return None;
        }
        if num.is_empty() {
            return Some(Self {
                num,
                den: vec![Q::one()],
            });
        }
        let g = p_gcd(&num, &den);
        let (n, _) = p_divrem(&num, &g);
        let (d, _) = p_divrem(&den, &g);
        let lead = d.last().cloned().unwrap_or_else(Q::one);
        Some(Self {
            num: n.iter().map(|c| c / &lead).collect(),
            den: d.iter().map(|c| c / &lead).collect(),
        })
    }

    /// A polynomial as a rational function.
    #[must_use]
    pub fn poly(p: Poly) -> Self {
        Self {
            num: trim(p),
            den: vec![Q::one()],
        }
    }

    /// The constant `c`.
    #[must_use]
    pub fn constant(c: Q) -> Self {
        Self::poly(vec![c])
    }

    /// Zero.
    #[must_use]
    pub fn zero() -> Self {
        Self::poly(Vec::new())
    }

    /// Whether this is zero.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.num.is_empty()
    }

    /// The sum.
    #[must_use]
    pub fn add(
        &self,
        other: &Self,
    ) -> Self {
        let num = p_add(&p_mul(&self.num, &other.den), &p_mul(&other.num, &self.den));
        Self::new(num, p_mul(&self.den, &other.den)).unwrap_or_else(Self::zero)
    }

    /// The negation.
    #[must_use]
    pub fn neg(&self) -> Self {
        Self {
            num: p_neg(&self.num),
            den: self.den.clone(),
        }
    }

    /// The difference.
    #[must_use]
    pub fn sub(
        &self,
        other: &Self,
    ) -> Self {
        self.add(&other.neg())
    }

    /// The product.
    #[must_use]
    pub fn mul(
        &self,
        other: &Self,
    ) -> Self {
        Self::new(p_mul(&self.num, &other.num), p_mul(&self.den, &other.den))
            .unwrap_or_else(Self::zero)
    }

    /// The quotient, or `None` when dividing by zero.
    #[must_use]
    pub fn div(
        &self,
        other: &Self,
    ) -> Option<Self> {
        Self::new(p_mul(&self.num, &other.den), p_mul(&self.den, &other.num))
    }

    /// The derivative.
    #[must_use]
    pub fn deriv(&self) -> Self {
        let num = p_sub(
            &p_mul(&p_deriv(&self.num), &self.den),
            &p_mul(&self.num, &p_deriv(&self.den)),
        );
        Self::new(num, p_mul(&self.den, &self.den)).unwrap_or_else(Self::zero)
    }

    /// The first `n` Taylor coefficients at the origin, or `None` when the
    /// function has a pole there.
    #[must_use]
    pub fn series(
        &self,
        n: usize,
    ) -> Option<Vec<Q>> {
        let d0 = self.den.first().filter(|c| !c.is_zero())?;
        let mut out: Vec<Q> = Vec::with_capacity(n);
        for k in 0..n {
            let mut acc = self.num.get(k).cloned().unwrap_or_else(Q::zero);
            for j in 1..=k.min(self.den.len() - 1) {
                acc -= &self.den[j] * &out[k - j];
            }
            out.push(acc / d0);
        }
        Some(out)
    }
}

/// A function analytic at the origin with rational Taylor coefficients.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Atom {
    /// A rational function without a pole at the origin.
    Rational(RatFn),
    /// `exp(g)` for a rational `g` with `g(0) = 0`.
    Exp(RatFn),
    /// `u^q` for a rational `u` with `u(0) = 1` and a rational exponent `q`.
    Power(RatFn, Q),
    /// `log(u)` for a rational `u` with `u(0) = 1`.
    Log(RatFn),
    /// `sin(g)` for a polynomial `g` with `g(0) = 0`.
    Sin(Poly),
    /// `cos(g)`, likewise.
    Cos(Poly),
    /// `sinh(g)`, likewise.
    Sinh(Poly),
    /// `cosh(g)`, likewise.
    Cosh(Poly),
    /// `atan(g)` for a rational `g` with `g(0) = 0`.
    Atan(RatFn),
    /// `atanh(g)`, likewise.
    Atanh(RatFn),
    /// `asin(g)` for a rational `g` with `g(0) = 0`.
    Asin(RatFn),
}

/// An expression in the D-finite closure of the atoms.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Expr {
    /// An atom.
    Atom(Atom),
    /// A sum.
    Add(Box<Self>, Box<Self>),
    /// A product.
    Mul(Box<Self>, Box<Self>),
    /// A rational multiple.
    Scale(Q, Box<Self>),
    /// A non-negative integer power.
    Pow(Box<Self>, u32),
    /// The derivative.
    Diff(Box<Self>),
}

/// Size limits; exceeding one makes a computation return `None`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Limits {
    /// Highest order of a differential equation.
    pub max_order: usize,
    /// Highest degree of a coefficient polynomial.
    pub max_degree: usize,
    /// Most Taylor coefficients computed for a zero test.
    pub max_terms: usize,
}

impl Default for Limits {
    fn default() -> Self {
        Self {
            max_order: 12,
            max_degree: 400,
            max_terms: 4000,
        }
    }
}

/// A linear differential operator `p_0 + p_1 D + ... + p_r D^r` with
/// polynomial coefficients and `p_r != 0`. The zero function has the
/// operator `[1]` (order 0).
pub type Operator = Vec<Poly>;

fn valid_at_origin(
    r: &RatFn,
    value: i64,
) -> bool {
    r.series(1).is_some_and(|s| s.first() == Some(&qi(value)))
}

impl Atom {
    /// Whether the atom meets the conditions stated on its variant.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        match self {
            | Self::Rational(r) => r.series(1).is_some(),
            | Self::Exp(g) | Self::Atan(g) | Self::Atanh(g) | Self::Asin(g) => {
                valid_at_origin(g, 0)
            },
            | Self::Power(u, _) | Self::Log(u) => valid_at_origin(u, 1),
            | Self::Sin(g) | Self::Cos(g) | Self::Sinh(g) | Self::Cosh(g) => {
                g.first().is_none_or(Zero::is_zero)
            },
        }
    }

    /// The first `n` Taylor coefficients, or `None` for an invalid atom.
    #[must_use]
    pub fn series(
        &self,
        n: usize,
    ) -> Option<Vec<Q>> {
        if !self.is_valid() {
            return None;
        }
        match self {
            | Self::Rational(r) => r.series(n),
            | Self::Exp(g) => solve_log_derivative(&g.deriv().series(n)?, n),
            | Self::Power(u, q) => {
                let s = RatFn::constant(q.clone()).mul(&u.deriv()).div(u)?;
                solve_log_derivative(&s.series(n)?, n)
            },
            | Self::Log(u) => Some(integrate(&u.deriv().div(u)?.series(n)?, n)),
            | Self::Atan(g) => {
                let w = g.deriv().div(&RatFn::constant(Q::one()).add(&g.mul(g)))?;
                Some(integrate(&w.series(n)?, n))
            },
            | Self::Atanh(g) => {
                let w = g.deriv().div(&RatFn::constant(Q::one()).sub(&g.mul(g)))?;
                Some(integrate(&w.series(n)?, n))
            },
            | Self::Asin(g) => {
                // (1 - g^2)^(-1/2) has logarithmic derivative g g' / (1 - g^2).
                let one_minus = RatFn::constant(Q::one()).sub(&g.mul(g));
                let s = g.mul(&g.deriv()).div(&one_minus)?;
                let root = solve_log_derivative(&s.series(n)?, n)?;
                let w = series_mul(&g.deriv().series(n)?, &root, n);
                Some(integrate(&w, n))
            },
            | Self::Sin(g) => Some(circular(g, n, false).0),
            | Self::Cos(g) => Some(circular(g, n, false).1),
            | Self::Sinh(g) => Some(circular(g, n, true).0),
            | Self::Cosh(g) => Some(circular(g, n, true).1),
        }
    }

    /// A differential operator annihilating the atom.
    #[must_use]
    pub fn annihilator(&self) -> Option<Operator> {
        if !self.is_valid() {
            return None;
        }
        // f' = s f  gives  den(s) D - num(s).
        let first_order = |s: &RatFn| -> Operator { vec![p_neg(&s.num), s.den.clone()] };
        // f' = w with w' = s w  gives  (den(s) D - num(s)) D.
        let integral_of =
            |s: &RatFn| -> Operator { vec![Vec::new(), p_neg(&s.num), s.den.clone()] };
        let log_derivative = |w: &RatFn| -> Option<RatFn> { w.deriv().div(w) };
        let one = RatFn::constant(Q::one());
        Some(match self {
            | Self::Rational(r) => {
                if r.is_zero() {
                    vec![vec![Q::one()]]
                } else {
                    first_order(&log_derivative(r)?)
                }
            },
            | Self::Exp(g) => first_order(&g.deriv()),
            | Self::Power(u, q) => first_order(&RatFn::constant(q.clone()).mul(&u.deriv()).div(u)?),
            | Self::Log(u) => {
                let w = u.deriv().div(u)?;
                if w.is_zero() {
                    vec![Vec::new(), vec![Q::one()]]
                } else {
                    integral_of(&log_derivative(&w)?)
                }
            },
            | Self::Atan(g) | Self::Atanh(g) => {
                let g2 = g.mul(g);
                let denominator = if matches!(self, Self::Atan(_)) {
                    one.add(&g2)
                } else {
                    one.sub(&g2)
                };
                let w = g.deriv().div(&denominator)?;
                if w.is_zero() {
                    vec![Vec::new(), vec![Q::one()]]
                } else {
                    integral_of(&log_derivative(&w)?)
                }
            },
            | Self::Asin(g) => {
                let dg = g.deriv();
                if dg.is_zero() {
                    vec![Vec::new(), vec![Q::one()]]
                } else {
                    // w = g' (1 - g^2)^(-1/2): w'/w = g''/g' + g g'/(1 - g^2).
                    let s = dg
                        .deriv()
                        .div(&dg)?
                        .add(&g.mul(&dg).div(&one.sub(&g.mul(g)))?);
                    integral_of(&s)
                }
            },
            | Self::Sin(g) | Self::Cos(g) | Self::Sinh(g) | Self::Cosh(g) => {
                let d1 = p_deriv(g);
                if d1.is_empty() {
                    // A constant: sin(0) = 0 or cos(0) = 1.
                    return Some(vec![Vec::new(), vec![Q::one()]]);
                }
                // g' f'' - g'' f' +- g'^3 f = 0.
                let cube = p_mul(&d1, &p_mul(&d1, &d1));
                let p0 = if matches!(self, Self::Sin(_) | Self::Cos(_)) {
                    cube
                } else {
                    p_neg(&cube)
                };
                vec![p0, p_neg(&p_deriv(&d1)), d1]
            },
        })
    }
}

/// The series of `f` with `f' = s f`, `f(0) = 1`.
fn solve_log_derivative(
    s: &[Q],
    n: usize,
) -> Option<Vec<Q>> {
    let mut f: Vec<Q> = Vec::with_capacity(n);
    for k in 0..n {
        if k == 0 {
            f.push(Q::one());
            continue;
        }
        // k f_k = sum_{j=0}^{k-1} s_j f_{k-1-j}.
        let mut acc = Q::zero();
        for j in 0..k {
            acc += s.get(j)? * &f[k - 1 - j];
        }
        f.push(acc / qi(k as i64));
    }
    Some(f)
}

/// The series of the antiderivative vanishing at the origin.
fn integrate(
    w: &[Q],
    n: usize,
) -> Vec<Q> {
    (0..n)
        .map(|k| {
            if k == 0 {
                Q::zero()
            } else {
                w.get(k - 1).cloned().unwrap_or_else(Q::zero) / qi(k as i64)
            }
        })
        .collect()
}

fn series_mul(
    a: &[Q],
    b: &[Q],
    n: usize,
) -> Vec<Q> {
    let mut out = vec![Q::zero(); n];
    for (i, x) in a.iter().enumerate().take(n) {
        if x.is_zero() {
            continue;
        }
        for (j, y) in b.iter().enumerate().take(n - i) {
            out[i + j] += x * y;
        }
    }
    out
}

fn series_deriv(
    a: &[Q],
    n: usize,
) -> Vec<Q> {
    (0..n)
        .map(|k| a.get(k + 1).map_or_else(Q::zero, |c| c * qi(k as i64 + 1)))
        .collect()
}

/// The series of `(sin g, cos g)`, or of `(sinh g, cosh g)`: the solution of
/// `s' = g' c`, `c' = -+ g' s` with `s(0) = 0`, `c(0) = 1`.
fn circular(
    g: &[Q],
    n: usize,
    hyperbolic: bool,
) -> (Vec<Q>, Vec<Q>) {
    let dg = p_deriv(g);
    let mut s: Vec<Q> = Vec::with_capacity(n);
    let mut c: Vec<Q> = Vec::with_capacity(n);
    for k in 0..n {
        if k == 0 {
            s.push(Q::zero());
            c.push(Q::one());
            continue;
        }
        let (mut sa, mut ca) = (Q::zero(), Q::zero());
        for j in 0..k.min(dg.len()) {
            sa += &dg[j] * &c[k - 1 - j];
            ca += &dg[j] * &s[k - 1 - j];
        }
        s.push(sa / qi(k as i64));
        c.push(if hyperbolic { ca } else { -ca } / qi(k as i64));
    }
    (s, c)
}

impl Expr {
    /// The first `n` Taylor coefficients at the origin.
    #[must_use]
    pub fn series(
        &self,
        n: usize,
    ) -> Option<Vec<Q>> {
        match self {
            | Self::Atom(atom) => atom.series(n),
            | Self::Add(a, b) => {
                let (x, y) = (a.series(n)?, b.series(n)?);
                Some(x.iter().zip(&y).map(|(p, q)| p + q).collect())
            },
            | Self::Mul(a, b) => Some(series_mul(&a.series(n)?, &b.series(n)?, n)),
            | Self::Scale(c, a) => Some(a.series(n)?.iter().map(|x| x * c).collect()),
            | Self::Pow(a, e) => {
                let base = a.series(n)?;
                let mut out = vec![Q::zero(); n];
                if let Some(first) = out.first_mut() {
                    *first = Q::one();
                }
                for _ in 0..*e {
                    out = series_mul(&out, &base, n);
                }
                Some(out)
            },
            | Self::Diff(a) => Some(series_deriv(&a.series(n + 1)?, n)),
        }
    }
}

// ---------------------------------------------------------------------------
// Closure: annihilators of sums, products and derivatives
// ---------------------------------------------------------------------------

/// Removes the polynomial content of an operator and clears denominators.
fn primitive(op: Vec<RatFn>) -> Operator {
    // Common denominator.
    let mut den: Poly = vec![Q::one()];
    for c in &op {
        let g = p_gcd(&den, &c.den);
        let (q, _) = p_divrem(&c.den, &g);
        den = p_mul(&den, &q);
    }
    let mut polys: Vec<Poly> = op
        .iter()
        .map(|c| {
            let (q, _) = p_divrem(&den, &c.den);
            p_mul(&c.num, &q)
        })
        .collect();
    let mut content: Poly = Vec::new();
    for p in &polys {
        content = p_gcd(&content, p);
    }
    if !content.is_empty() {
        for p in &mut polys {
            *p = p_divrem(p, &content).0;
        }
    }
    // Integer coefficients with a positive leading one, for a canonical form.
    let mut scale = BigInt::one();
    for c in polys.iter().flatten() {
        scale = num_integer::Integer::lcm(&scale, c.denom());
    }
    let mut factor = Q::from(scale);
    if polys
        .last()
        .and_then(|p| p.last())
        .is_some_and(Signed::is_negative)
    {
        factor = -factor;
    }
    polys.iter().map(|p| p_scale(p, &factor)).collect()
}

/// The derivative of a coordinate vector `sum c_i f^(i)` (`i < r`) modulo
/// the operator `op` of order `r` that annihilates `f`.
fn derive_coordinates(
    c: &[RatFn],
    op: &[Poly],
) -> Option<Vec<RatFn>> {
    let r = op.len() - 1;
    let lead = RatFn::poly(op[r].clone());
    let mut out: Vec<RatFn> = c.iter().map(RatFn::deriv).collect();
    for i in 1..r {
        out[i] = out[i].add(&c[i - 1]);
    }
    if r > 0 {
        // c_{r-1} f^(r) = -c_{r-1} sum_{i<r} (p_i / p_r) f^(i).
        let top = &c[r - 1];
        for i in 0..r {
            let ratio = RatFn::poly(op[i].clone()).div(&lead)?;
            out[i] = out[i].sub(&top.mul(&ratio));
        }
    }
    Some(out)
}

/// Solves `sum_k a_k v_k = target` over `Q(x)` for the given vectors, or
/// returns `None` when `target` is not in their span.
fn express(
    vectors: &[Vec<RatFn>],
    target: &[RatFn],
) -> Option<Vec<RatFn>> {
    let rows = target.len();
    let cols = vectors.len();
    // Augmented matrix, one row per coordinate.
    let mut m: Vec<Vec<RatFn>> = (0..rows)
        .map(|i| {
            let mut row: Vec<RatFn> = vectors.iter().map(|v| v[i].clone()).collect();
            row.push(target[i].clone());
            row
        })
        .collect();
    let mut pivot_cols = Vec::new();
    let mut at = 0;
    for col in 0..cols {
        let Some(found) = (at..rows).find(|&i| !m[i][col].is_zero()) else {
            continue;
        };
        m.swap(at, found);
        let pivot = m[at][col].clone();
        for j in col..=cols {
            m[at][j] = m[at][j].div(&pivot)?;
        }
        for i in 0..rows {
            if i != at && !m[i][col].is_zero() {
                let factor = m[i][col].clone();
                for j in col..=cols {
                    let term = factor.mul(&m[at][j]);
                    m[i][j] = m[i][j].sub(&term);
                }
            }
        }
        pivot_cols.push(col);
        at += 1;
        if at == rows {
            break;
        }
    }
    if (at..rows).any(|i| !m[i][cols].is_zero()) {
        return None;
    }
    let mut solution = vec![RatFn::zero(); cols];
    for (row, &col) in pivot_cols.iter().enumerate() {
        solution[col] = m[row][cols].clone();
    }
    Some(solution)
}

/// Finds the operator of least order annihilating the element whose
/// successive derivatives (as coordinate vectors) `next` produces, starting
/// from `start`.
fn minimal_operator(
    start: Vec<RatFn>,
    mut next: impl FnMut(&[RatFn]) -> Option<Vec<RatFn>>,
    limits: &Limits,
) -> Option<Operator> {
    if start.iter().all(RatFn::is_zero) {
        return Some(vec![vec![Q::one()]]);
    }
    let mut derivatives = vec![start];
    loop {
        let order = derivatives.len();
        if order > limits.max_order {
            return None;
        }
        let candidate = next(&derivatives[order - 1])?;
        if let Some(combination) = express(&derivatives, &candidate) {
            // D^order - sum a_k D^k.
            let mut op: Vec<RatFn> = combination.iter().map(RatFn::neg).collect();
            op.push(RatFn::constant(Q::one()));
            let op = primitive(op);
            if op.iter().any(|p| p.len() > limits.max_degree + 1) {
                return None;
            }
            return Some(op);
        }
        derivatives.push(candidate);
    }
}

/// An operator annihilating `f + g`, given operators for `f` and `g`.
fn sum_operator(
    a: &[Poly],
    b: &[Poly],
    limits: &Limits,
) -> Option<Operator> {
    let (ra, rb) = (a.len() - 1, b.len() - 1);
    // Coordinates: f, f', ..., then g, g', ...
    let mut start = vec![RatFn::zero(); ra + rb];
    if ra > 0 {
        start[0] = RatFn::constant(Q::one());
    }
    if rb > 0 {
        start[ra] = RatFn::constant(Q::one());
    }
    minimal_operator(
        start,
        |v| {
            let mut left = derive_coordinates(&v[..ra], a)?;
            left.extend(derive_coordinates(&v[ra..], b)?);
            Some(left)
        },
        limits,
    )
}

/// An operator annihilating `f g`.
fn product_operator(
    a: &[Poly],
    b: &[Poly],
    limits: &Limits,
) -> Option<Operator> {
    let (ra, rb) = (a.len() - 1, b.len() - 1);
    if ra == 0 || rb == 0 {
        // One factor is zero.
        return Some(vec![vec![Q::one()]]);
    }
    // Coordinates: f^(i) g^(j) at index i * rb + j.
    let mut start = vec![RatFn::zero(); ra * rb];
    start[0] = RatFn::constant(Q::one());
    minimal_operator(
        start,
        |v| {
            let mut out = vec![RatFn::zero(); ra * rb];
            // Leibniz: differentiate the coefficient, then f, then g.
            for i in 0..ra {
                // Row i as coordinates in g.
                let row: Vec<RatFn> = v[i * rb..(i + 1) * rb].to_vec();
                let dg = derive_coordinates(&row, b)?;
                for j in 0..rb {
                    out[i * rb + j] = out[i * rb + j].add(&dg[j]);
                }
            }
            for j in 0..rb {
                // Column j as coordinates in f; the coefficient derivative
                // was counted above, so remove it here.
                let column: Vec<RatFn> = (0..ra).map(|i| v[i * rb + j].clone()).collect();
                let df = derive_coordinates(&column, a)?;
                for i in 0..ra {
                    let pure = df[i].sub(&column[i].deriv());
                    out[i * rb + j] = out[i * rb + j].add(&pure);
                }
            }
            Some(out)
        },
        limits,
    )
}

/// An operator annihilating `f'`.
fn derivative_operator(
    a: &[Poly],
    limits: &Limits,
) -> Option<Operator> {
    let r = a.len() - 1;
    if r == 0 {
        return Some(vec![vec![Q::one()]]);
    }
    // f' in the coordinates of f.
    let mut unit = vec![RatFn::zero(); r];
    unit[0] = RatFn::constant(Q::one());
    let start = derive_coordinates(&unit, a)?;
    minimal_operator(start, |v| derive_coordinates(v, a), limits)
}

/// A differential operator annihilating `expr`.
///
/// It is the operator of least order found in the vector space that the
/// closure construction provides, which is usually but not always the
/// minimal operator of the function. Returns `None` when an atom is invalid
/// or a limit is exceeded.
#[must_use]
pub fn annihilator(
    expr: &Expr,
    limits: &Limits,
) -> Option<Operator> {
    match expr {
        | Expr::Atom(atom) => atom.annihilator().map(|op| primitive_polys(&op)),
        | Expr::Add(a, b) => {
            sum_operator(&annihilator(a, limits)?, &annihilator(b, limits)?, limits)
        },
        | Expr::Mul(a, b) => {
            product_operator(&annihilator(a, limits)?, &annihilator(b, limits)?, limits)
        },
        | Expr::Scale(c, a) => {
            if c.is_zero() {
                Some(vec![vec![Q::one()]])
            } else {
                annihilator(a, limits)
            }
        },
        | Expr::Pow(a, e) => {
            let base = annihilator(a, limits)?;
            // The constant 1 for a zeroth power.
            let mut acc: Operator = vec![Vec::new(), vec![Q::one()]];
            for _ in 0..*e {
                acc = product_operator(&acc, &base, limits)?;
            }
            Some(acc)
        },
        | Expr::Diff(a) => derivative_operator(&annihilator(a, limits)?, limits),
    }
}

fn primitive_polys(op: &[Poly]) -> Operator {
    primitive(op.iter().cloned().map(RatFn::poly).collect())
}

// ---------------------------------------------------------------------------
// Recurrences and the zero test
// ---------------------------------------------------------------------------

/// `(m + s)(m + s - 1)...(m + s - i + 1)` as a polynomial in `m`.
fn falling(
    s: i64,
    i: usize,
) -> Poly {
    let mut out: Poly = vec![Q::one()];
    for t in 0..i {
        out = p_mul(&out, &[qi(s - t as i64), Q::one()]);
    }
    out
}

/// The recurrence of the Taylor coefficients `a_n` of a solution of `op`:
/// pairs `(s, c_s)` with `sum_s c_s(m) a_{m+s} = 0` for every integer `m`
/// (coefficients with a negative index are zero), shifts ascending.
#[must_use]
pub fn recurrence(op: &[Poly]) -> Vec<(i64, Poly)> {
    let mut by_shift: std::collections::BTreeMap<i64, Poly> = std::collections::BTreeMap::new();
    for (i, p) in op.iter().enumerate() {
        for (j, c) in p.iter().enumerate() {
            if c.is_zero() {
                continue;
            }
            // x^j D^i contributes c * (m + s)^(falling i) a_{m+s}, s = i - j.
            let s = i as i64 - j as i64;
            let term = p_scale(&falling(s, i), c);
            let slot = by_shift.entry(s).or_default();
            *slot = p_add(slot, &term);
        }
    }
    by_shift
        .into_iter()
        .filter(|(_, p)| !p.is_empty())
        .collect()
}

/// How many leading Taylor coefficients decide whether a solution of `op`
/// is zero: every coefficient the recurrence does not determine from earlier
/// ones lies below this index. `None` when the bound exceeds
/// `limits.max_terms`.
#[must_use]
pub fn deciding_terms(
    op: &[Poly],
    limits: &Limits,
) -> Option<usize> {
    let rec = recurrence(op);
    let Some((top, lead)) = rec.last() else {
        return Some(0);
    };
    // a_{m+top} is determined whenever lead(m) != 0 and m >= 0; undetermined
    // are the indices below `top` and those at non-negative integer roots.
    let mut need = (*top).max(0);
    // Integer roots are bounded by Cauchy's bound.
    let lead_coefficient = lead.last()?.abs();
    let bound = lead
        .iter()
        .map(|c| (c.abs() / &lead_coefficient).ceil().to_integer())
        .max()?
        + BigInt::one();
    let bound = bound.to_i64().filter(|b| *b <= limits.max_terms as i64)?;
    for m in 0..=bound {
        if p_eval(lead, &qi(m)).is_zero() {
            need = need.max(m + top + 1);
        }
    }
    let need = usize::try_from(need.max(0)).ok()?;
    (need <= limits.max_terms).then_some(need)
}

/// Whether `expr` is identically zero. `None` when it cannot be decided
/// within the limits (or an atom is invalid); `Some` is a proof either way.
#[must_use]
pub fn is_zero(
    expr: &Expr,
    limits: &Limits,
) -> Option<bool> {
    let op = annihilator(expr, limits)?;
    let terms = deciding_terms(&op, limits)?;
    let series = expr.series(terms)?;
    Some(series.iter().all(Zero::is_zero))
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

    use super::*;

    fn x() -> Poly {
        vec![qi(0), qi(1)]
    }

    fn rat(p: Poly) -> RatFn {
        RatFn::poly(p)
    }

    fn atom(a: Atom) -> Expr {
        Expr::Atom(a)
    }

    fn add(
        a: Expr,
        b: Expr,
    ) -> Expr {
        Expr::Add(Box::new(a), Box::new(b))
    }

    fn mul(
        a: Expr,
        b: Expr,
    ) -> Expr {
        Expr::Mul(Box::new(a), Box::new(b))
    }

    fn scale(
        c: i64,
        a: Expr,
    ) -> Expr {
        Expr::Scale(qi(c), Box::new(a))
    }

    fn sq(a: Expr) -> Expr {
        Expr::Pow(Box::new(a), 2)
    }

    fn constant(c: i64) -> Expr {
        atom(Atom::Rational(RatFn::constant(qi(c))))
    }

    fn sin() -> Expr {
        atom(Atom::Sin(x()))
    }

    fn cos() -> Expr {
        atom(Atom::Cos(x()))
    }

    /// The operator applied to a truncated series vanishes to the order the
    /// truncation allows.
    fn annihilates(
        op: &[Poly],
        expr: &Expr,
    ) -> bool {
        let n = 24;
        let mut d = expr.series(n + op.len()).unwrap();
        let mut total = vec![Q::zero(); n];
        for p in op {
            let term = series_mul(p, &d, n);
            for (t, v) in total.iter_mut().zip(&term) {
                *t += v;
            }
            d = series_deriv(&d, d.len() - 1);
        }
        total.iter().all(Zero::is_zero)
    }

    #[test]
    fn series_of_the_atoms() {
        let q = |n: i64, d: i64| Q::new(BigInt::from(n), BigInt::from(d));
        assert_eq!(
            Atom::Exp(rat(x())).series(5).unwrap(),
            vec![qi(1), qi(1), q(1, 2), q(1, 6), q(1, 24)]
        );
        assert_eq!(
            Atom::Sin(x()).series(6).unwrap(),
            vec![qi(0), qi(1), qi(0), q(-1, 6), qi(0), q(1, 120)]
        );
        assert_eq!(
            Atom::Cosh(x()).series(5).unwrap(),
            vec![qi(1), qi(0), q(1, 2), qi(0), q(1, 24)]
        );
        let one_plus_x = rat(vec![qi(1), qi(1)]);
        assert_eq!(
            Atom::Log(one_plus_x.clone()).series(5).unwrap(),
            vec![qi(0), qi(1), q(-1, 2), q(1, 3), q(-1, 4)]
        );
        assert_eq!(
            Atom::Power(one_plus_x.clone(), q(1, 2)).series(4).unwrap(),
            vec![qi(1), q(1, 2), q(-1, 8), q(1, 16)]
        );
        assert_eq!(
            Atom::Atan(rat(x())).series(6).unwrap(),
            vec![qi(0), qi(1), qi(0), q(-1, 3), qi(0), q(1, 5)]
        );
        assert_eq!(
            Atom::Asin(rat(x())).series(6).unwrap(),
            vec![qi(0), qi(1), qi(0), q(1, 6), qi(0), q(3, 40)]
        );
        // 1 / (1 - x - x^2): the Fibonacci numbers.
        let fib = RatFn::new(vec![qi(1)], vec![qi(1), qi(-1), qi(-1)]).unwrap();
        assert_eq!(
            fib.series(7).unwrap(),
            [1, 1, 2, 3, 5, 8, 13].map(qi).to_vec()
        );
        // Conditions on the atoms are enforced.
        assert!(Atom::Exp(rat(vec![qi(1)])).series(3).is_none());
        assert!(Atom::Log(rat(x())).series(3).is_none());
        assert!(
            Atom::Rational(RatFn::new(vec![qi(1)], x()).unwrap())
                .series(3)
                .is_none()
        );
    }

    #[test]
    fn annihilators_annihilate() {
        let limits = Limits::default();
        let one_plus_x = rat(vec![qi(1), qi(1)]);
        let x2 = vec![qi(0), qi(0), qi(1)];
        let atoms = vec![
            Atom::Rational(RatFn::new(vec![qi(1), qi(2)], vec![qi(1), qi(0), qi(1)]).unwrap()),
            Atom::Exp(rat(x2.clone())),
            Atom::Power(
                one_plus_x.clone(),
                Q::new(BigInt::from(-3), BigInt::from(2)),
            ),
            Atom::Log(one_plus_x),
            Atom::Sin(x2.clone()),
            Atom::Cos(vec![qi(0), qi(3), qi(1)]),
            Atom::Sinh(x()),
            Atom::Cosh(x2),
            Atom::Atan(rat(x())),
            Atom::Atanh(rat(x())),
            Atom::Asin(rat(x())),
        ];
        for a in &atoms {
            let e = atom(a.clone());
            let op = annihilator(&e, &limits).unwrap();
            assert!(annihilates(&op, &e), "{a:?}: {op:?}");
        }
        // Closure: a sum, a product, a power and a derivative.
        let e = add(mul(sin(), atom(atoms[1].clone())), sq(cos()));
        let op = annihilator(&e, &limits).unwrap();
        assert!(annihilates(&op, &e));
        let d = Expr::Diff(Box::new(mul(atom(atoms[8].clone()), sin())));
        let op = annihilator(&d, &limits).unwrap();
        assert!(annihilates(&op, &d));
    }

    #[test]
    fn known_equations_are_found() {
        let limits = Limits::default();
        // sin: f'' + f = 0.
        assert_eq!(
            annihilator(&sin(), &limits).unwrap(),
            vec![vec![qi(1)], vec![], vec![qi(1)]]
        );
        // sin + cos satisfies the same equation, not one of order four.
        assert_eq!(
            annihilator(&add(sin(), cos()), &limits).unwrap(),
            vec![vec![qi(1)], vec![], vec![qi(1)]]
        );
        // The product of two *arbitrary* solutions of f'' + f = 0 satisfies
        // a third-order equation, f''' + 4 f' = 0, and nothing of lower
        // order; that is what the closure finds, although sin * cos =
        // sin(2x)/2 happens to satisfy f'' + 4 f = 0 as well (see the caveat
        // on `annihilator`).
        assert_eq!(
            annihilator(&mul(sin(), cos()), &limits).unwrap(),
            vec![vec![], vec![qi(4)], vec![], vec![qi(1)]]
        );
        // exp(x) / (1 - x): (x - 1) f' + (2 - x) f = 0.
        let e = mul(
            atom(Atom::Exp(rat(x()))),
            atom(Atom::Rational(
                RatFn::new(vec![qi(1)], vec![qi(1), qi(-1)]).unwrap(),
            )),
        );
        assert_eq!(
            annihilator(&e, &limits).unwrap(),
            vec![vec![qi(2), qi(-1)], vec![qi(-1), qi(1)]]
        );
    }

    #[test]
    fn recurrences_of_taylor_coefficients() {
        // exp: f' - f = 0 gives (m + 1) a_{m+1} - a_m = 0.
        let rec = recurrence(&[vec![qi(-1)], vec![qi(1)]]);
        assert_eq!(rec, vec![(0, vec![qi(-1)]), (1, vec![qi(1), qi(1)])]);
        // sin: (m + 2)(m + 1) a_{m+2} + a_m = 0.
        let rec = recurrence(&[vec![qi(1)], vec![], vec![qi(1)]]);
        assert_eq!(rec[1], (2, vec![qi(2), qi(3), qi(1)]));
        assert_eq!(
            deciding_terms(&[vec![qi(1)], vec![], vec![qi(1)]], &Limits::default()),
            Some(2)
        );
        // x f' - 3 f = 0 has the solution x^3: the coefficient a_3 is free.
        assert_eq!(
            deciding_terms(&[vec![qi(-3)], vec![qi(0), qi(1)]], &Limits::default()),
            Some(4)
        );
    }

    #[test]
    fn identities_are_decided() {
        let limits = Limits::default();
        let one_plus_x = rat(vec![qi(1), qi(1)]);
        let exp = |sign: i64| atom(Atom::Exp(rat(vec![qi(0), qi(sign)])));
        let cases: Vec<(&str, Expr, bool)> = vec![
            (
                "sin^2 + cos^2 - 1",
                add(add(sq(sin()), sq(cos())), constant(-1)),
                true,
            ),
            (
                "sin^2 + cos^2 - 2",
                add(add(sq(sin()), sq(cos())), constant(-2)),
                false,
            ),
            (
                "sin(2x) - 2 sin cos",
                add(
                    atom(Atom::Sin(vec![qi(0), qi(2)])),
                    scale(-2, mul(sin(), cos())),
                ),
                true,
            ),
            (
                "cosh^2 - sinh^2 - 1",
                add(
                    add(
                        sq(atom(Atom::Cosh(x()))),
                        scale(-1, sq(atom(Atom::Sinh(x())))),
                    ),
                    constant(-1),
                ),
                true,
            ),
            (
                "exp(x) exp(-x) - 1",
                add(mul(exp(1), exp(-1)), constant(-1)),
                true,
            ),
            (
                "sqrt(1+x)^2 - (1+x)",
                add(
                    sq(atom(Atom::Power(
                        one_plus_x.clone(),
                        Q::new(BigInt::one(), BigInt::from(2)),
                    ))),
                    atom(Atom::Rational(one_plus_x.clone().neg())),
                ),
                true,
            ),
            (
                "d/dx atan - 1/(1+x^2)",
                add(
                    Expr::Diff(Box::new(atom(Atom::Atan(rat(x()))))),
                    atom(Atom::Rational(
                        RatFn::new(vec![qi(-1)], vec![qi(1), qi(0), qi(1)]).unwrap(),
                    )),
                ),
                true,
            ),
            (
                "2 sinh - exp(x) + exp(-x)",
                add(
                    scale(2, atom(Atom::Sinh(x()))),
                    add(scale(-1, exp(1)), exp(-1)),
                ),
                true,
            ),
            (
                "log(1+x) - x agrees to first order only",
                add(
                    atom(Atom::Log(one_plus_x)),
                    atom(Atom::Rational(rat(x()).neg())),
                ),
                false,
            ),
            (
                "sin - x + x^3/6 agrees to fourth order only",
                add(
                    sin(),
                    atom(Atom::Rational(rat(vec![
                        qi(0),
                        qi(-1),
                        qi(0),
                        Q::new(BigInt::one(), BigInt::from(6)),
                    ]))),
                ),
                false,
            ),
        ];
        for (name, expr, want) in cases {
            assert_eq!(is_zero(&expr, &limits), Some(want), "{name}");
        }
        // A limit that is too small gives no answer rather than a wrong one.
        let tight = Limits {
            max_order: 1,
            ..limits
        };
        assert_eq!(
            is_zero(&add(add(sq(sin()), sq(cos())), constant(-1)), &tight),
            None
        );
    }
}
