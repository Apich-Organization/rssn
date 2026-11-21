//! Forward-mode dual numbers and reverse-mode tape automatic differentiation.
//!
//! Automatic Differentiation (AD): forward-mode dual numbers ([`Dual`]) and
//! tape-based reverse-mode AD ([`Tape`], [`Var`]).

use std::cell::RefCell;
use std::ops::Add;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Neg;
use std::ops::Sub;
use std::rc::Rc;

/// Error returned by automatic differentiation helpers.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AdError {
    /// Input vector or output vector was empty.
    Empty,
    /// The result belongs to a tape other than the input variables.
    ForeignTape,
    /// Closure output length changed between directional passes.
    ShapeMismatch {
        /// Expected output length.
        expected: usize,
        /// Actual output length.
        got: usize,
    },
}

impl std::fmt::Display for AdError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            | Self::Empty => f.write_str("autodiff input or output is empty"),
            | Self::ForeignTape => f.write_str("autodiff result belongs to a different tape"),
            | Self::ShapeMismatch { expected, got } => {
                write!(f, "autodiff expected {expected} outputs, got {got}")
            },
        }
    }
}

impl std::error::Error for AdError {}

/// A forward-mode dual number `val + deriv * eps` (`eps^2 = 0`).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Dual {
    /// Primal value.
    pub val: f64,
    /// Tangent (derivative) value.
    pub deriv: f64,
}

impl Dual {
    /// Creates a constant dual number with derivative `0`.
    #[must_use]
    pub const fn constant(val: f64) -> Self {
        Self { val, deriv: 0.0 }
    }

    /// Creates an independent variable dual number with derivative `1`.
    #[must_use]
    pub const fn var(val: f64) -> Self {
        Self { val, deriv: 1.0 }
    }

    /// Creates a dual number with explicit value and derivative.
    #[must_use]
    pub const fn new(
        val: f64,
        deriv: f64,
    ) -> Self {
        Self { val, deriv }
    }

    /// Sine of `self`.
    #[must_use]
    pub fn sin(self) -> Self {
        Self::new(self.val.sin(), self.deriv * self.val.cos())
    }

    /// Cosine of `self`.
    #[must_use]
    pub fn cos(self) -> Self {
        Self::new(self.val.cos(), -self.deriv * self.val.sin())
    }

    /// Tangent of `self`.
    #[must_use]
    pub fn tan(self) -> Self {
        let c = self.val.cos();
        Self::new(self.val.tan(), self.deriv / (c * c))
    }

    /// Natural exponential `exp(self)`.
    #[must_use]
    pub fn exp(self) -> Self {
        let e = self.val.exp();
        Self::new(e, self.deriv * e)
    }

    /// Natural logarithm `ln(self)`.
    #[must_use]
    pub fn ln(self) -> Self {
        Self::new(self.val.ln(), self.deriv / self.val)
    }

    /// Square root `sqrt(self)`.
    #[must_use]
    pub fn sqrt(self) -> Self {
        let s = self.val.sqrt();
        Self::new(s, self.deriv / (2.0 * s))
    }

    /// Floating-point power `self^p`.
    #[must_use]
    pub fn powf(
        self,
        p: f64,
    ) -> Self {
        if p == 0.0 {
            return Self::constant(1.0);
        }
        Self::new(self.val.powf(p), self.deriv * p * self.val.powf(p - 1.0))
    }

    /// Integer power `self^n`.
    ///
    /// For a finite nonzero base and finite tangent, a representable tangent
    /// is retained even when the primal or an unscaled power overflows or
    /// underflows.
    #[must_use]
    pub fn powi(
        self,
        n: i32,
    ) -> Self {
        self.integer_power::<false>(n).0
    }

    /// Evaluates `self^n` and both derivatives along `self.deriv`.
    ///
    /// Treats the input as `self.val + t * self.deriv` and differentiates
    /// at `t = 0`. The returned dual holds the value and first derivative;
    /// the second component is the second derivative with respect to `t`.
    /// The powers zero and one have zero second derivative.
    #[must_use]
    pub fn powi_with_second_derivative(
        self,
        n: i32,
    ) -> (Self, f64) {
        self.integer_power::<true>(n)
    }

    fn integer_power<const SECOND: bool>(
        self,
        n: i32,
    ) -> (Self, f64) {
        match n {
            | 0 => return (Self::constant(1.0), 0.0),
            | 1 => return (self, 0.0),
            | _ => {},
        }
        let n = i64::from(n);
        let offset = if SECOND { 2 } else { 1 };
        // Within this box every nonzero intermediate is normal. Otherwise
        // check the actual products before using the unscaled result.
        let moderate = (-64..=64).contains(&n)
            && moderate_power_factor(self.val)
            && (self.deriv == 0.0 || moderate_power_factor(self.deriv));
        let tangent = self.deriv * n as f64;
        let squared_tangent = if SECOND {
            tangent * self.deriv
        } else {
            0.0
        };
        let anchor = floating_integer_power(self.val, n - offset);
        let curvature = if SECOND {
            squared_tangent * (n - 1) as f64 * anchor
        } else {
            0.0
        };

        // Keep the IEEE behavior at poles and for nonfinite inputs.
        if !moderate && (self.val == 0.0 || !self.val.is_finite() || !self.deriv.is_finite()) {
            let value = floating_integer_power(self.val, n);
            let slope = if SECOND {
                floating_integer_power(self.val, n - 1)
            } else {
                anchor
            };
            return (Self::new(value, tangent * slope), curvature);
        }

        // Reuse one power to evaluate the primal and all requested derivatives.
        let slope = if SECOND {
            anchor * self.val
        } else {
            anchor
        };
        let value = slope * self.val;
        let derivative = tangent * slope;
        if moderate
            || (anchor.is_normal()
                && slope.is_normal()
                && value.is_normal()
                && (self.deriv == 0.0
                    || (tangent.is_normal()
                        && derivative.is_normal()
                        && (!SECOND || (squared_tangent.is_normal() && curvature.is_normal())))))
        {
            return (Self::new(value, derivative), curvature);
        }
        self.scaled_integer_power::<SECOND>(n)
    }

    fn scaled_integer_power<const SECOND: bool>(
        self,
        n: i64,
    ) -> (Self, f64) {
        let offset = if SECOND { 2 } else { 1 };
        let base = PowerFactor::positive(self.val.abs());
        let anchor = PowerFactor::power(self.val.abs(), n - offset);
        let slope = if SECOND {
            anchor.product(base)
        } else {
            anchor
        };
        let sign = if self.val.is_sign_negative() && n % 2 != 0 {
            -1.0
        } else {
            1.0
        };
        let slope_sign = if self.val.is_sign_negative() {
            -sign
        } else {
            sign
        };
        let tangent_sign = slope_sign * (n as f64).signum() * self.deriv.signum();
        let value = slope.product(base).round(sign);
        if self.deriv == 0.0 {
            let curvature = if SECOND {
                0.0_f64.copysign(sign)
            } else {
                0.0
            };
            return (Self::new(value, 0.0_f64.copysign(tangent_sign)), curvature);
        }
        let seed = PowerFactor::positive(self.deriv.abs());
        let weighted_seed = seed.product(PowerFactor::positive((n as f64).abs()));
        let derivative = slope.product(weighted_seed).round(tangent_sign);
        let curvature = if SECOND {
            anchor
                .product(weighted_seed)
                .product(seed)
                .product(PowerFactor::positive(((n - 1) as f64).abs()))
                .round(sign)
        } else {
            0.0
        };
        (Self::new(value, derivative), curvature)
    }

    /// Absolute value `|self|`.
    #[must_use]
    pub fn abs(self) -> Self {
        Self::new(self.val.abs(), self.deriv * self.val.signum())
    }

    /// Hyperbolic sine `sinh(self)`.
    #[must_use]
    pub fn sinh(self) -> Self {
        Self::new(self.val.sinh(), self.deriv * self.val.cosh())
    }

    /// Hyperbolic cosine `cosh(self)`.
    #[must_use]
    pub fn cosh(self) -> Self {
        Self::new(self.val.cosh(), self.deriv * self.val.sinh())
    }

    /// Hyperbolic tangent `tanh(self)`.
    #[must_use]
    pub fn tanh(self) -> Self {
        let c = self.val.cosh();
        Self::new(self.val.tanh(), self.deriv / (c * c))
    }

    /// Inverse sine `asin(self)`.
    #[must_use]
    pub fn asin(self) -> Self {
        Self::new(
            self.val.asin(),
            self.deriv / (1.0 - self.val * self.val).sqrt(),
        )
    }

    /// Inverse cosine `acos(self)`.
    #[must_use]
    pub fn acos(self) -> Self {
        Self::new(
            self.val.acos(),
            -self.deriv / (1.0 - self.val * self.val).sqrt(),
        )
    }

    /// Inverse tangent `atan(self)`.
    #[must_use]
    pub fn atan(self) -> Self {
        Self::new(self.val.atan(), self.deriv / (1.0 + self.val * self.val))
    }
}

fn floating_integer_power(
    base: f64,
    power: i64,
) -> f64 {
    match power {
        | 0 => 1.0,
        | 1 => base,
        | 2 => base * base,
        | -64..=64 => base.powi(power as i32),
        | _ => {
            // Derivative exponents can extend below i32::MIN. Keep integer
            // parity even when a wider exponent requires powf.
            let magnitude = base.abs().powf(power as f64);
            if base.is_sign_negative() && power % 2 != 0 {
                -magnitude
            } else {
                magnitude
            }
        },
    }
}

const fn moderate_power_factor(value: f64) -> bool {
    let exponent = (value.to_bits() >> 52) & 2047;
    exponent.wrapping_sub(1015) <= 16 // Binary exponent -8 through 8.
}

/// A positive factor with a normalized significand and unrestricted exponent.
#[derive(Clone, Copy)]
struct PowerFactor {
    significand: f64,
    exponent: i64,
}

impl PowerFactor {
    fn positive(value: f64) -> Self {
        let (normalized, offset) = if value.is_subnormal() {
            (value * 18_014_398_509_481_984.0, -54)
        } else {
            (value, 0)
        };
        let bits = normalized.to_bits();
        Self {
            significand: f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | (1023_u64 << 52)),
            exponent: ((bits >> 52) & 2047) as i64 - 1023 + offset,
        }
    }

    #[allow(
        clippy::single_call_fn,
        reason = "isolates range reduction for a powered factor"
    )]
    fn power(
        base: f64,
        power: i64,
    ) -> Self {
        let mut exponent = power as f64;
        let mut value = base.powf(exponent);
        if value.is_normal() {
            return Self::positive(value);
        }
        let estimate = base.log2() * exponent;
        // At most two binary64 factors and two i32 factors are applied to
        // this power. Beyond this bound they cannot restore a finite value.
        if estimate.abs() > 8192.0 {
            return Self {
                significand: 1.0,
                exponent: if estimate > 0.0 {
                    8192
                } else {
                    -8192
                },
            };
        }
        let mut squarings = 0;
        while !value.is_normal() {
            exponent *= 0.5;
            squarings += 1;
            value = base.powf(exponent);
        }
        let mut factor = Self::positive(value);
        for _ in 0..squarings {
            factor = factor.product(factor);
        }
        factor
    }

    fn product(
        self,
        rhs: Self,
    ) -> Self {
        let significand = self.significand * rhs.significand;
        let exponent = self.exponent + rhs.exponent;
        if significand >= 2.0 {
            Self {
                significand: significand * 0.5,
                exponent: exponent + 1,
            }
        } else {
            Self {
                significand,
                exponent,
            }
        }
    }

    fn round(
        self,
        sign: f64,
    ) -> f64 {
        let value = match self.exponent {
            | 1024.. => f64::INFINITY,
            | e @ -1022..=1023 => {
                let fraction = self.significand.to_bits() & 0x000f_ffff_ffff_ffff;
                f64::from_bits(fraction | ((e + 1023) as u64) << 52)
            },
            | e @ -1074..=-1023 => self.significand * f64::from_bits(1 << (e + 1074)),
            | -1075 => (self.significand * 0.5) * f64::from_bits(1),
            | _ => 0.0,
        };
        value.copysign(sign)
    }
}

impl Add for Dual {
    type Output = Self;

    fn add(
        self,
        rhs: Self,
    ) -> Self {
        Self::new(self.val + rhs.val, self.deriv + rhs.deriv)
    }
}

impl Add<f64> for Dual {
    type Output = Self;

    fn add(
        self,
        rhs: f64,
    ) -> Self {
        Self::new(self.val + rhs, self.deriv)
    }
}

impl Sub for Dual {
    type Output = Self;

    fn sub(
        self,
        rhs: Self,
    ) -> Self {
        Self::new(self.val - rhs.val, self.deriv - rhs.deriv)
    }
}

impl Sub<f64> for Dual {
    type Output = Self;

    fn sub(
        self,
        rhs: f64,
    ) -> Self {
        Self::new(self.val - rhs, self.deriv)
    }
}

impl Mul for Dual {
    type Output = Self;

    fn mul(
        self,
        rhs: Self,
    ) -> Self {
        Self::new(
            self.val * rhs.val,
            self.deriv * rhs.val + self.val * rhs.deriv,
        )
    }
}

impl Mul<f64> for Dual {
    type Output = Self;

    fn mul(
        self,
        rhs: f64,
    ) -> Self {
        Self::new(self.val * rhs, self.deriv * rhs)
    }
}

impl Div for Dual {
    type Output = Self;

    fn div(
        self,
        rhs: Self,
    ) -> Self {
        let value = self.val / rhs.val;
        // Form the quotient first: squaring the denominator can overflow or
        // underflow even when both the value and tangent are representable.
        let correction = value * rhs.deriv;
        let direct = (self.deriv - correction) / rhs.val;
        let deriv =
            if !direct.is_finite() || (correction == 0.0 && value != 0.0 && rhs.deriv != 0.0) {
                let scaled = self.deriv / rhs.val - value * (rhs.deriv / rhs.val);
                if scaled.is_finite() {
                    scaled
                } else {
                    direct
                }
            } else {
                direct
            };
        Self::new(value, deriv)
    }
}

impl Div<f64> for Dual {
    type Output = Self;

    fn div(
        self,
        rhs: f64,
    ) -> Self {
        Self::new(self.val / rhs, self.deriv / rhs)
    }
}

impl Neg for Dual {
    type Output = Self;

    fn neg(self) -> Self {
        Self::new(-self.val, -self.deriv)
    }
}

/// Exact first derivative `(f(x), f'(x))` of `f(Dual) -> Dual` via forward-mode AD.
pub fn forward_derivative<F>(
    f: &mut F,
    x: f64,
) -> (f64, f64)
where
    F: FnMut(Dual) -> Dual,
{
    let out = f(Dual::var(x));
    (out.val, out.deriv)
}

/// Exact gradient of `f(&[Dual]) -> Dual` at `x` via forward-mode AD.
///
/// # Errors
///
/// [`AdError::Empty`] when `x` is empty.
pub fn forward_grad<F>(
    f: &mut F,
    x: &[f64],
) -> Result<Vec<f64>, AdError>
where
    F: FnMut(&[Dual]) -> Dual,
{
    if x.is_empty() {
        return Err(AdError::Empty);
    }
    let n = x.len();
    let mut g = Vec::with_capacity(n);
    let mut d_in: Vec<Dual> = x.iter().map(|&v| Dual::constant(v)).collect();
    for i in 0..n {
        if let Some(slot) = d_in.get_mut(i) {
            slot.deriv = 1.0;
        }
        g.push(f(&d_in).deriv);
        if let Some(slot) = d_in.get_mut(i) {
            slot.deriv = 0.0;
        }
    }
    Ok(g)
}

/// Exact Jacobian of `f(&[Dual]) -> Vec<Dual>` at `x` via forward-mode AD.
///
/// Returns a row-major `m × n` Jacobian.
///
/// # Errors
///
/// [`AdError::Empty`] when `x` is empty or `f` returns an empty vector, and
/// [`AdError::ShapeMismatch`] if the output length varies across passes.
pub fn forward_jacobian<F>(
    f: &mut F,
    x: &[f64],
) -> Result<Vec<f64>, AdError>
where
    F: FnMut(&[Dual]) -> Vec<Dual>,
{
    if x.is_empty() {
        return Err(AdError::Empty);
    }
    let n = x.len();
    let mut d_in: Vec<Dual> = x.iter().map(|&v| Dual::constant(v)).collect();
    if let Some(slot) = d_in.get_mut(0) {
        slot.deriv = 1.0;
    }
    let first = f(&d_in);
    let m = first.len();
    if m == 0 {
        return Err(AdError::Empty);
    }
    let total = m.saturating_mul(n);
    let mut jac = vec![0.0; total];
    for (r, val) in first.iter().enumerate() {
        if let Some(cell) = jac.get_mut(r.saturating_mul(n)) {
            *cell = val.deriv;
        }
    }
    if let Some(slot) = d_in.get_mut(0) {
        slot.deriv = 0.0;
    }
    for c in 1..n {
        if let Some(slot) = d_in.get_mut(c) {
            slot.deriv = 1.0;
        }
        let col_out = f(&d_in);
        if col_out.len() != m {
            return Err(AdError::ShapeMismatch {
                expected: m,
                got: col_out.len(),
            });
        }
        for (r, val) in col_out.iter().enumerate() {
            if let Some(cell) = jac.get_mut(r.saturating_mul(n).saturating_add(c)) {
                *cell = val.deriv;
            }
        }
        if let Some(slot) = d_in.get_mut(c) {
            slot.deriv = 0.0;
        }
    }
    Ok(jac)
}

// --- Reverse-mode AD (Wengert tape) ---

#[derive(Clone, Debug)]
struct TapeNode {
    weights: [f64; 2],
    deps: [usize; 2],
    arity: u8,
}

/// Wengert tape for reverse-mode automatic differentiation.
#[derive(Clone, Debug, Default)]
pub struct Tape {
    nodes: Rc<RefCell<Vec<TapeNode>>>,
}

impl Tape {
    /// Creates an empty AD tape.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of recorded nodes, for bounded host adapters.
    #[must_use]
    pub fn len(&self) -> usize {
        self.nodes.borrow().len()
    }

    /// Whether the tape contains no variables or operations.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.nodes.borrow().is_empty()
    }

    /// Creates a differentiable input variable on the tape.
    #[must_use]
    pub fn var(
        &self,
        val: f64,
    ) -> Var {
        let mut nodes = self.nodes.borrow_mut();
        let idx = nodes.len();
        nodes.push(TapeNode {
            weights: [0.0, 0.0],
            deps: [0, 0],
            arity: 0,
        });
        Var {
            val,
            idx,
            tape: Rc::clone(&self.nodes),
        }
    }

    fn push_unary(
        &self,
        val: f64,
        dep: usize,
        w: f64,
    ) -> Var {
        let mut nodes = self.nodes.borrow_mut();
        let idx = nodes.len();
        nodes.push(TapeNode {
            weights: [w, 0.0],
            deps: [dep, 0],
            arity: 1,
        });
        Var {
            val,
            idx,
            tape: Rc::clone(&self.nodes),
        }
    }

    fn push_binary(
        &self,
        val: f64,
        dep0: usize,
        w0: f64,
        dep1: usize,
        w1: f64,
    ) -> Var {
        let mut nodes = self.nodes.borrow_mut();
        let idx = nodes.len();
        nodes.push(TapeNode {
            weights: [w0, w1],
            deps: [dep0, dep1],
            arity: 2,
        });
        Var {
            val,
            idx,
            tape: Rc::clone(&self.nodes),
        }
    }
}

/// A reverse-mode AD variable tracked on a [`Tape`].
#[derive(Clone, Debug)]
pub struct Var {
    /// Primal value.
    pub val: f64,
    idx: usize,
    tape: Rc<RefCell<Vec<TapeNode>>>,
}

impl Var {
    /// Whether two variables belong to exactly the same differentiation tape.
    #[must_use]
    pub fn shares_tape(
        &self,
        other: &Self,
    ) -> bool {
        Rc::ptr_eq(&self.tape, &other.tape)
    }

    /// Number of nodes retained by this variable's tape.
    #[must_use]
    pub fn node_count(&self) -> usize {
        self.tape.borrow().len()
    }

    // Only checked adapters use this immediately after a failed single-node
    // operation. Never removes nodes referenced by an existing handle.
    pub(crate) fn discard_newest(&self) {
        let mut nodes = self.tape.borrow_mut();
        if self.idx.checked_add(1) == Some(nodes.len()) {
            nodes.pop();
        }
    }

    /// Runs reverse-mode backpropagation from `self` and returns the adjoint
    /// vector indexed by [`Var::grad_in`].
    #[must_use]
    pub fn backward(&self) -> Vec<f64> {
        let nodes = self.tape.borrow();
        let mut grads = vec![0.0; nodes.len()];
        if let Some(g) = grads.get_mut(self.idx) {
            *g = 1.0;
        }
        for i in (0..=self.idx).rev() {
            let g = grads.get(i).copied().unwrap_or(0.0);
            if g == 0.0 {
                continue;
            }
            let Some(node) = nodes.get(i) else {
                continue;
            };
            if node.arity >= 1
                && let Some(d0) = grads.get_mut(node.deps[0])
            {
                *d0 += g * node.weights[0];
            }
            if node.arity == 2
                && let Some(d1) = grads.get_mut(node.deps[1])
            {
                *d1 += g * node.weights[1];
            }
        }
        grads
    }

    /// Reads this variable's derivative from a gradient vector returned by [`Var::backward`].
    #[must_use]
    pub fn grad_in(
        &self,
        grads: &[f64],
    ) -> f64 {
        grads.get(self.idx).copied().unwrap_or(0.0)
    }

    fn tape(&self) -> Tape {
        Tape {
            nodes: Rc::clone(&self.tape),
        }
    }

    /// Sine of `self`.
    #[must_use]
    pub fn sin(&self) -> Self {
        self.tape()
            .push_unary(self.val.sin(), self.idx, self.val.cos())
    }

    /// Cosine of `self`.
    #[must_use]
    pub fn cos(&self) -> Self {
        self.tape()
            .push_unary(self.val.cos(), self.idx, -self.val.sin())
    }

    /// Natural exponential `exp(self)`.
    #[must_use]
    pub fn exp(&self) -> Self {
        let e = self.val.exp();
        self.tape().push_unary(e, self.idx, e)
    }

    /// Natural logarithm `ln(self)`.
    #[must_use]
    pub fn ln(&self) -> Self {
        self.tape()
            .push_unary(self.val.ln(), self.idx, 1.0 / self.val)
    }

    /// Square root `sqrt(self)`.
    #[must_use]
    pub fn sqrt(&self) -> Self {
        let s = self.val.sqrt();
        self.tape().push_unary(s, self.idx, 0.5 / s)
    }

    /// Floating-point power `self^p`.
    #[must_use]
    pub fn powf(
        &self,
        p: f64,
    ) -> Self {
        if p == 0.0 {
            return self.tape().push_unary(1.0, self.idx, 0.0);
        }
        self.tape()
            .push_unary(self.val.powf(p), self.idx, p * self.val.powf(p - 1.0))
    }

    /// Hyperbolic tangent `tanh(self)`.
    #[must_use]
    pub fn tanh(&self) -> Self {
        let c = self.val.cosh();
        self.tape()
            .push_unary(self.val.tanh(), self.idx, 1.0 / (c * c))
    }
}

impl Add for &Var {
    type Output = Var;

    fn add(
        self,
        rhs: &Var,
    ) -> Var {
        assert!(
            Rc::ptr_eq(&self.tape, &rhs.tape),
            "AD variables must belong to the same tape"
        );
        self.tape()
            .push_binary(self.val + rhs.val, self.idx, 1.0, rhs.idx, 1.0)
    }
}

impl Sub for &Var {
    type Output = Var;

    fn sub(
        self,
        rhs: &Var,
    ) -> Var {
        assert!(
            Rc::ptr_eq(&self.tape, &rhs.tape),
            "AD variables must belong to the same tape"
        );
        self.tape()
            .push_binary(self.val - rhs.val, self.idx, 1.0, rhs.idx, -1.0)
    }
}

impl Mul for &Var {
    type Output = Var;

    fn mul(
        self,
        rhs: &Var,
    ) -> Var {
        assert!(
            Rc::ptr_eq(&self.tape, &rhs.tape),
            "AD variables must belong to the same tape"
        );
        self.tape()
            .push_binary(self.val * rhs.val, self.idx, rhs.val, rhs.idx, self.val)
    }
}

impl Div for &Var {
    type Output = Var;

    fn div(
        self,
        rhs: &Var,
    ) -> Var {
        assert!(
            Rc::ptr_eq(&self.tape, &rhs.tape),
            "AD variables must belong to the same tape"
        );
        self.tape().push_binary(
            self.val / rhs.val,
            self.idx,
            1.0 / rhs.val,
            rhs.idx,
            -(self.val / rhs.val) / rhs.val,
        )
    }
}

impl Neg for &Var {
    type Output = Var;

    fn neg(self) -> Var {
        self.tape().push_unary(-self.val, self.idx, -1.0)
    }
}

/// Exact value and gradient `(f(x), grad f(x))` of `f(&[Var]) -> Var` at `x`
/// in a single reverse-mode pass.
///
/// # Errors
///
/// [`AdError::Empty`] when `x` is empty.
pub fn reverse_grad<F>(
    f: &mut F,
    x: &[f64],
) -> Result<(f64, Vec<f64>), AdError>
where
    F: FnMut(&[Var]) -> Var,
{
    if x.is_empty() {
        return Err(AdError::Empty);
    }
    let tape = Tape::new();
    let vars: Vec<Var> = x.iter().map(|&v| tape.var(v)).collect();
    let out = f(&vars);
    if !Rc::ptr_eq(&out.tape, &tape.nodes) {
        return Err(AdError::ForeignTape);
    }
    let grads = out.backward();
    let g = vars.iter().map(|v| v.grad_in(&grads)).collect();
    Ok((out.val, g))
}
