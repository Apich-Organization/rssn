//! Complex balls with outward-rounded interval arithmetic (Krawczyk proofs).
//!
//! # Complex balls with outward-rounded real interval arithmetic
//!
//! A [`Cball`] is a rectangle `[re.lo, re.hi] + i [im.lo, im.hi]` of the
//! complex plane. Its arithmetic is the rectangle arithmetic of
//! [`Interval`], which rounds every
//! elementary operation outward by one ulp, so the computed rectangle
//! encloses the exact image of every point of the operands. That is what
//! makes the Krawczyk test of [`crate::kernels::homotopy_cert`] a proof
//! rather than a numerical check.
//!
//! A system is a list of [`CTerm`]s (a complex coefficient and an exponent
//! vector per term). [`real_terms`] lifts the real `PolyTerm` systems of the
//! homotopy kernel; [`eval_terms`], [`eval_jacobian_terms`] and
//! [`values_at_terms`] evaluate them over balls and points.
//! [`eval_value_and_jacobian_terms`] and [`values_and_jacobian_at_terms`]
//! share the work when both outputs are needed.

use crate::backend::tape::Interval;
use crate::backend::tape::Scalar;
use crate::kernels::rootfind::PolyTerm;
use num_complex::Complex64 as Cx;
use std::collections::HashMap;
use std::ops::Add;
use std::ops::Mul;
use std::ops::Neg;
use std::ops::Sub;

/// A complex rectangle `[re] + i [im]`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cball {
    /// Real part interval.
    pub re: Interval,
    /// Imaginary part interval.
    pub im: Interval,
}

impl Cball {
    /// The point `re + i im`.
    #[must_use]
    pub const fn point(
        re: f64,
        im: f64,
    ) -> Self {
        Self {
            re: Interval::point(re),
            im: Interval::point(im),
        }
    }

    /// A [`Cx`] as a point ball.
    #[must_use]
    pub const fn of(value: Cx) -> Self {
        Self::point(value.re, value.im)
    }

    /// The square `|z - (re + i im)|_inf <= radius` (one radius for both
    /// parts), widened outward.
    #[must_use]
    pub fn square(
        re: f64,
        im: f64,
        radius: f64,
    ) -> Self {
        let (lo, hi) = ((re - radius).next_down(), (re + radius).next_up());
        let (ilo, ihi) = ((im - radius).next_down(), (im + radius).next_up());
        Self {
            re: Interval { lo, hi },
            im: Interval { lo: ilo, hi: ihi },
        }
    }

    /// A ball with separate end points.
    #[must_use]
    pub const fn rectangle(
        re: (f64, f64),
        im: (f64, f64),
    ) -> Self {
        Self {
            re: Interval { lo: re.0, hi: re.1 },
            im: Interval { lo: im.0, hi: im.1 },
        }
    }

    /// `0`.
    #[must_use]
    pub const fn zero() -> Self {
        Self {
            re: Interval::point(0.0),
            im: Interval::point(0.0),
        }
    }

    /// `1`.
    #[must_use]
    pub const fn one() -> Self {
        Self {
            re: Interval::point(1.0),
            im: Interval::point(0.0),
        }
    }

    /// The centre.
    #[must_use]
    pub const fn centre(&self) -> Cx {
        Cx::new(
            f64::midpoint(self.re.lo, self.re.hi),
            f64::midpoint(self.im.lo, self.im.hi),
        )
    }

    /// The larger of the two half-widths.
    #[must_use]
    pub fn half_width(&self) -> f64 {
        interval_half_width(self.re).max(interval_half_width(self.im))
    }

    /// Whether both parts are finite.
    #[must_use]
    pub const fn is_finite(&self) -> bool {
        self.re.hi.is_finite()
            && self.re.lo.is_finite()
            && self.im.hi.is_finite()
            && self.im.lo.is_finite()
    }

    /// Whether `z` lies in the rectangle.
    #[must_use]
    pub fn contains(
        &self,
        z: Cx,
    ) -> bool {
        self.re.contains(z.re) && self.im.contains(z.im)
    }

    /// A lower bound of `|z|` for `z` in the rectangle.
    #[must_use]
    pub fn modulus_lower(&self) -> f64 {
        if [self.re.lo, self.re.hi, self.im.lo, self.im.hi]
            .into_iter()
            .any(f64::is_nan)
        {
            return 0.0;
        }
        let re = if self.re.lo <= 0.0 && self.re.hi >= 0.0 {
            0.0
        } else {
            self.re.lo.abs().min(self.re.hi.abs())
        };
        let im = if self.im.lo <= 0.0 && self.im.hi >= 0.0 {
            0.0
        } else {
            self.im.lo.abs().min(self.im.hi.abs())
        };
        modulus_bound(re, im, false)
    }

    /// An upper bound of `|z|` over the rectangle.
    #[must_use]
    pub fn modulus_upper(&self) -> f64 {
        if [self.re.lo, self.re.hi, self.im.lo, self.im.hi]
            .into_iter()
            .any(f64::is_nan)
        {
            return f64::INFINITY;
        }
        modulus_bound(interval_abs_max(self.re), interval_abs_max(self.im), true)
    }

    /// The intersection, when the rectangles meet.
    #[must_use]
    pub fn intersect(
        &self,
        other: &Self,
    ) -> Option<Self> {
        let (rlo, rhi) = (self.re.lo.max(other.re.lo), self.re.hi.min(other.re.hi));
        let (ilo, ihi) = (self.im.lo.max(other.im.lo), self.im.hi.min(other.im.hi));
        (rlo <= rhi && ilo <= ihi).then_some(Self {
            re: Interval { lo: rlo, hi: rhi },
            im: Interval { lo: ilo, hi: ihi },
        })
    }

    /// Whether this rectangle lies in the interior of `outer` (strictly in
    /// both coordinates, which is the hypothesis of the Krawczyk theorem).
    #[must_use]
    pub fn inside_strictly(
        &self,
        outer: &Self,
    ) -> bool {
        outer.re.lo < self.re.lo
            && self.re.hi < outer.re.hi
            && outer.im.lo < self.im.lo
            && self.im.hi < outer.im.hi
    }

    /// A rectangular enclosure of `z^n` by binary exponentiation.
    #[must_use]
    pub fn powi(
        &self,
        n: u32,
    ) -> Self {
        if n != 0 && is_zero_interval(self.re) && is_zero_interval(self.im) {
            return *self;
        }
        let mut power = n;
        let mut factor = *self;
        let mut result = None;
        while power != 0 {
            if power & 1 != 0 {
                result = Some(result.map_or(factor, |prefix| prefix * factor));
            }
            power >>= 1;
            if power != 0 {
                factor = factor.square_value();
            }
        }
        result.unwrap_or_else(Self::one)
    }

    /// `self / other`. Returns `None` when the denominator enclosure
    /// cannot be separated from zero.
    #[must_use]
    pub fn div(
        &self,
        other: &Self,
    ) -> Option<Self> {
        let maximum = interval_abs_max(other.re).max(interval_abs_max(other.im));
        let (numerator, divisor) = if maximum > 0.0
            && maximum.is_finite()
            && !(f64::from_bits(523_u64 << 52)..=f64::from_bits(1523_u64 << 52)).contains(&maximum)
        {
            let exponent = ((maximum.to_bits() >> 52) & 0x7ff) as i32 - 1023;
            let shift = (-exponent).clamp(-1022, 1023);
            let scale = f64::from_bits(((shift + 1023) as u64) << 52);
            (self.times_real(scale), other.times_real(scale))
        } else {
            (*self, *other)
        };
        let denominator = interval_square(divisor.re).add(interval_square(divisor.im));
        if denominator.lo <= 0.0 {
            return None;
        }
        let inverse = Self {
            re: divisor.re.div(denominator),
            im: divisor.im.neg().div(denominator),
        };
        Some(numerator * inverse)
    }

    fn square_value(self) -> Self {
        if is_zero_interval(self.im) {
            return Self {
                re: interval_square(self.re),
                im: self.im,
            };
        }
        if is_zero_interval(self.re) {
            return Self {
                re: interval_square(self.im).neg(),
                im: self.re,
            };
        }
        Self {
            re: interval_square(self.re).sub(interval_square(self.im)),
            im: self.re.mul(self.im).mul(Interval::point(2.0)),
        }
    }

    fn times_real(
        self,
        factor: f64,
    ) -> Self {
        let interval = Interval::point(factor);
        let scaled = |value: Interval| {
            if is_zero_interval(value) {
                value
            } else {
                value.mul(interval)
            }
        };
        Self {
            re: scaled(self.re),
            im: scaled(self.im),
        }
    }
}

fn interval_half_width(interval: Interval) -> f64 {
    let width = interval.hi - interval.lo;
    if width.is_finite() {
        0.5 * width
    } else {
        0.5 * interval.hi - 0.5 * interval.lo
    }
}

fn is_zero_interval(interval: Interval) -> bool {
    interval.lo == 0.0 && interval.hi == 0.0
}

const fn interval_abs_max(interval: Interval) -> f64 {
    interval.lo.abs().max(interval.hi.abs())
}

fn interval_square(interval: Interval) -> Interval {
    if interval.lo.is_nan() || interval.hi.is_nan() {
        return Interval::ENTIRE;
    }
    if is_zero_interval(interval) {
        return Interval::point(0.0);
    }
    let far = interval_abs_max(interval);
    let near = if interval.contains(0.0) {
        0.0
    } else {
        interval.lo.abs().min(interval.hi.abs())
    };
    Interval {
        lo: (near * near).next_down().max(0.0),
        hi: (far * far).next_up(),
    }
}

/// Bound the modulus as max(a,b) * sqrt(1 + (min(a,b)/max(a,b))²).
/// Round each monotone operation in the requested direction.
fn modulus_bound(
    a: f64,
    b: f64,
    upper: bool,
) -> f64 {
    let far = a.max(b);
    let near = a.min(b);
    if near == 0.0 || !far.is_finite() {
        return far;
    }
    if let Some(exact) = exact_modulus(far, near) {
        return exact;
    }
    let round = |value: f64| {
        if upper {
            value.next_up()
        } else {
            value.next_down()
        }
    };
    let ratio = round(near / far).max(0.0);
    let square = round(ratio * ratio).max(0.0);
    let distance = round(1.0 + square);
    let bound = round(far * round(distance.sqrt()));
    bound.max(far)
}

/// Recognize an exact Pythagorean norm when all its residuals fit in f64.
#[allow(clippy::float_cmp, clippy::single_call_fn)]
fn exact_modulus(
    far: f64,
    near: f64,
) -> Option<f64> {
    // Squared products above this input threshold have no residual bits
    // below the f64 subnormal range.
    const MIN: f64 = f64::from_bits(543_u64 << 52); // 2^-480
    if near < MIN {
        return None;
    }
    let far_square = far * far;
    let near_square = near * near;
    let square = far_square + near_square;
    if !square.is_finite()
        || square - far_square != near_square
        || far.mul_add(far, -far_square) != 0.0
        || near.mul_add(near, -near_square) != 0.0
    {
        return None;
    }
    let root = square.sqrt();
    (root.mul_add(root, -square) == 0.0).then_some(root)
}

impl Add for Cball {
    type Output = Self;

    fn add(
        self,
        other: Self,
    ) -> Self {
        Self {
            re: self.re.add(other.re),
            im: self.im.add(other.im),
        }
    }
}

impl Sub for Cball {
    type Output = Self;

    fn sub(
        self,
        other: Self,
    ) -> Self {
        Self {
            re: self.re.sub(other.re),
            im: self.im.sub(other.im),
        }
    }
}

impl Mul for Cball {
    type Output = Self;

    fn mul(
        self,
        other: Self,
    ) -> Self {
        // (a + bi)(c + di) = (ac - bd) + (ad + bc)i, every product an
        // outward-rounded interval product.
        let re = self.re.mul(other.re).sub(self.im.mul(other.im));
        let im = self.re.mul(other.im).add(self.im.mul(other.re));
        Self { re, im }
    }
}

impl Neg for Cball {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            re: self.re.neg(),
            im: self.im.neg(),
        }
    }
}

/// One term of a complex polynomial: `coef * z_0^e_0 * z_1^e_1 * ...`.
#[derive(Clone, Debug, PartialEq)]
pub struct CTerm {
    /// Complex coefficient.
    pub coef: Cx,
    /// Exponents, one per variable.
    pub exps: Vec<u32>,
}

impl CTerm {
    /// A real-coefficient term.
    #[must_use]
    pub const fn real(
        coef: f64,
        exps: Vec<u32>,
    ) -> Self {
        Self {
            coef: Cx::new(coef, 0.0),
            exps,
        }
    }
}

/// Lifts a real system to complex terms.
#[must_use]
pub fn real_terms(system: &[Vec<PolyTerm>]) -> Vec<Vec<CTerm>> {
    system
        .iter()
        .map(|equation| {
            equation
                .iter()
                .map(|term| {
                    CTerm {
                        coef: Cx::new(term.coef, 0.0),
                        exps: term.exps.clone(),
                    }
                })
                .collect()
        })
        .collect()
}

/// Arithmetic used by the box and point polynomial evaluators.
trait PolynomialValue: Copy + Add<Output = Self> + Mul<Output = Self> {
    fn zero() -> Self;
    fn coefficient(value: Cx) -> Self;
    fn power(
        self,
        exponent: u32,
    ) -> Self;
    fn scale(
        self,
        factor: f64,
    ) -> Self;
}

impl PolynomialValue for Cball {
    fn zero() -> Self {
        Self::zero()
    }

    fn coefficient(value: Cx) -> Self {
        Self::of(value)
    }

    fn power(
        self,
        exponent: u32,
    ) -> Self {
        self.powi(exponent)
    }

    fn scale(
        self,
        factor: f64,
    ) -> Self {
        self.times_real(factor)
    }
}

impl PolynomialValue for Cx {
    fn zero() -> Self {
        Self::new(0.0, 0.0)
    }

    fn coefficient(value: Cx) -> Self {
        value
    }

    fn power(
        self,
        exponent: u32,
    ) -> Self {
        self.powu(exponent)
    }

    fn scale(
        self,
        factor: f64,
    ) -> Self {
        self * factor
    }
}

const INLINE_POWER_COLUMNS: usize = 8;
const INLINE_POWER_COUNT: usize = 15;

/// Powers zero and one need no storage. Small exponents have direct slots
/// and arbitrary exponents use a sparse map instead of an exponent-sized table.
struct PowerCache<'a, T> {
    point: &'a [T],
    enabled: bool,
    small: Vec<[Option<T>; INLINE_POWER_COUNT]>,
    large: Option<HashMap<(usize, u32), T>>,
}

impl<'a, T: PolynomialValue> PowerCache<'a, T> {
    const fn new(
        point: &'a [T],
        enabled: bool,
    ) -> Self {
        Self {
            point,
            enabled,
            small: Vec::new(),
            large: None,
        }
    }

    fn get(
        &mut self,
        column: usize,
        power: u32,
    ) -> T {
        let value = self.point[column];
        if power == 1 {
            return value;
        }
        if !self.enabled {
            return value.power(power);
        }
        if column < INLINE_POWER_COLUMNS && power <= 1 + INLINE_POWER_COUNT as u32 {
            if self.small.len() <= column {
                self.small
                    .resize_with(self.point.len().min(INLINE_POWER_COLUMNS), || {
                        [None; INLINE_POWER_COUNT]
                    });
            }
            return *self.small[column][(power - 2) as usize]
                .get_or_insert_with(|| value.power(power));
        }
        *self
            .large
            .get_or_insert_with(HashMap::new)
            .entry((column, power))
            .or_insert_with(|| value.power(power))
    }
}

#[derive(Clone, Copy)]
struct TermFactor<T> {
    column: usize,
    exponent: u32,
    prefix: T,
    value: T,
}

/// Evaluates a polynomial system with a power cache shared across its rows.
/// Prefix products and a backwards suffix sweep evaluate a term and all
/// of its partial derivatives without dividing by any input coordinate.
fn evaluate_system<T: PolynomialValue, const VALUES: bool, const JACOBIAN: bool>(
    system: &[Vec<CTerm>],
    point: &[T],
) -> (Vec<T>, Vec<Vec<T>>) {
    let mut values = Vec::with_capacity(if VALUES { system.len() } else { 0 });
    let mut jacobian = Vec::with_capacity(if JACOBIAN { system.len() } else { 0 });
    let use_cache = system.len() > 2 || system.iter().any(|row| row.len() > 2);
    let mut powers = PowerCache::new(point, use_cache);
    let mut inline = [TermFactor {
        column: 0,
        exponent: 0,
        prefix: T::zero(),
        value: T::zero(),
    }; INLINE_POWER_COLUMNS];
    let mut spilled = Vec::new();
    let short_term = point.len() <= INLINE_POWER_COLUMNS;

    for equation in system {
        let mut value = T::zero();
        let mut row = if JACOBIAN {
            vec![T::zero(); point.len()]
        } else {
            Vec::new()
        };
        for term in equation {
            if term.coef.re == 0.0 && term.coef.im == 0.0 {
                continue;
            }
            if JACOBIAN && !short_term {
                spilled.clear();
            }
            let mut count = 0;
            let mut prefix = T::coefficient(term.coef);
            for (column, &exponent) in term.exps.iter().take(point.len()).enumerate() {
                if exponent == 0 {
                    continue;
                }
                let factor = powers.get(column, exponent);
                if JACOBIAN {
                    let part = TermFactor {
                        column,
                        exponent,
                        prefix,
                        value: factor,
                    };
                    if short_term {
                        inline[count] = part;
                        count += 1;
                    } else {
                        spilled.push(part);
                    }
                }
                prefix = prefix * factor;
            }
            if VALUES {
                value = value + prefix;
            }
            if JACOBIAN {
                let factors = if short_term {
                    &inline[..count]
                } else {
                    spilled.as_slice()
                };
                let mut suffix = None;
                for factor in factors.iter().rev() {
                    let left = if factor.exponent == 1 {
                        factor.prefix
                    } else {
                        factor.prefix.scale(f64::from(factor.exponent))
                            * powers.get(factor.column, factor.exponent - 1)
                    };
                    let derivative = suffix.map_or(left, |right| left * right);
                    row[factor.column] = row[factor.column] + derivative;
                    suffix = Some(suffix.map_or(factor.value, |right| factor.value * right));
                }
            }
        }
        if VALUES {
            values.push(value);
        }
        if JACOBIAN {
            jacobian.push(row);
        }
    }
    (values, jacobian)
}

/// `f_i(X)` for every equation of the system (an enclosure of the range
/// over the box).
#[must_use]
pub fn eval_terms(
    system: &[Vec<CTerm>],
    x: &[Cball],
) -> Vec<Cball> {
    evaluate_system::<_, true, false>(system, x).0
}

/// The complex Jacobian `∂f_i/∂z_j` over the box.
///
/// Each entry is an enclosure. A missing or zero exponent gives no
/// contribution to that column.
#[must_use]
pub fn eval_jacobian_terms(
    system: &[Vec<CTerm>],
    x: &[Cball],
) -> Vec<Vec<Cball>> {
    evaluate_system::<_, false, true>(system, x).1
}

/// Encloses both the values and the complex Jacobian of a system over a box.
///
/// The evaluator can reuse powers across rows and shares each term's products
/// across the value and partial derivatives. Missing exponents are zero, and
/// exponents beyond the number of coordinates are ignored.
#[must_use]
pub fn eval_value_and_jacobian_terms(
    system: &[Vec<CTerm>],
    x: &[Cball],
) -> (Vec<Cball>, Vec<Vec<Cball>>) {
    evaluate_system::<_, true, true>(system, x)
}

/// The values `f_i(x)` at a point, in plain `f64` arithmetic.
#[must_use]
pub fn values_at_terms(
    system: &[Vec<CTerm>],
    x: &[Cx],
) -> Vec<Cx> {
    evaluate_system::<_, true, false>(system, x).0
}

/// The complex Jacobian at a point.
#[must_use]
pub fn jacobian_at_terms(
    system: &[Vec<CTerm>],
    x: &[Cx],
) -> Vec<Vec<Cx>> {
    evaluate_system::<_, false, true>(system, x).1
}

/// Evaluates the values and complex Jacobian at a point using one power cache.
///
/// The arithmetic is ordinary complex `f64` arithmetic. The exponent
/// convention matches [`eval_value_and_jacobian_terms`].
#[must_use]
pub fn values_and_jacobian_at_terms(
    system: &[Vec<CTerm>],
    x: &[Cx],
) -> (Vec<Cx>, Vec<Vec<Cx>>) {
    evaluate_system::<_, true, true>(system, x)
}

/// `f_i(X)` of a real system over a box.
#[must_use]
pub fn eval_polys(
    system: &[Vec<PolyTerm>],
    x: &[Cball],
) -> Vec<Cball> {
    eval_terms(&real_terms(system), x)
}

/// The Jacobian of a real system over a box.
#[must_use]
pub fn eval_jacobian(
    system: &[Vec<PolyTerm>],
    x: &[Cball],
) -> Vec<Vec<Cball>> {
    eval_jacobian_terms(&real_terms(system), x)
}

/// The values of a real system at a point.
#[must_use]
pub fn values_at(
    system: &[Vec<PolyTerm>],
    x: &[Cx],
) -> Vec<Cx> {
    values_at_terms(&real_terms(system), x)
}

/// The Jacobian of a real system at a point.
#[must_use]
pub fn jacobian_at(
    system: &[Vec<PolyTerm>],
    x: &[Cx],
) -> Vec<Vec<Cx>> {
    jacobian_at_terms(&real_terms(system), x)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ball(
        re: (f64, f64),
        im: (f64, f64),
    ) -> Cball {
        Cball::rectangle(re, im)
    }

    #[test]
    fn arithmetic_encloses_the_pointwise_image() {
        let a = ball((1.0, 2.0), (-1.0, 0.5));
        let b = ball((-3.0, 0.5), (0.25, 2.0));
        let sum = a + b;
        let difference = a - b;
        let product = a * b;
        for &xr in &[1.0_f64, 1.5, 2.0] {
            for &xi in &[-1.0_f64, -0.5, 0.5] {
                for &yr in &[-3.0_f64, -1.0, 0.5] {
                    for &yi in &[0.25_f64, 1.0, 2.0] {
                        let (x, y) = (Cx::new(xr, xi), Cx::new(yr, yi));
                        assert!(sum.contains(x + y));
                        assert!(difference.contains(x - y));
                        assert!(product.contains(x * y));
                    }
                }
            }
        }
        assert!(sum.contains(Cx::new(2.0, -0.5)));
        assert!(!sum.contains(Cx::new(3.0, 1.0)));
    }

    #[test]
    fn division_and_powers() {
        let z = ball((1.0, 2.0), (0.5, 1.0));
        let w = ball((3.0, 4.0), (-2.0, -1.0));
        let q = z.div(&w).unwrap();
        for &zr in &[1.0_f64, 2.0] {
            for &zi in &[0.5_f64, 1.0] {
                for &wr in &[3.0_f64, 4.0] {
                    for &wi in &[-2.0_f64, -1.0] {
                        let value = Cx::new(zr, zi) / Cx::new(wr, wi);
                        assert!(q.contains(value), "{value} not in {q:?}");
                    }
                }
            }
        }
        // a denominator that may vanish gives no information
        assert!(z.div(&Cball::square(0.0, 0.0, 1.0)).is_none());
        let cube = ball((1.0, 2.0), (0.0, 1.0)).powi(3);
        for &re in &[1.0_f64, 2.0] {
            for &im in &[0.0_f64, 1.0] {
                assert!(cube.contains(Cx::new(re, im).powi(3)));
            }
        }
        assert!(z.powi(0).contains(Cx::new(1.0, 0.0)));
    }

    #[test]
    fn rectangles_and_norms() {
        let a = ball((1.0, 2.0), (3.0, 4.0));
        let b = ball((1.5, 3.0), (3.5, 9.0));
        let c = a.intersect(&b).unwrap();
        assert!((c.re.lo - 1.5).abs() < 1e-12 && (c.im.hi - 4.0).abs() < 1e-12);
        assert!(a.intersect(&ball((5.0, 6.0), (0.0, 1.0))).is_none());
        assert!(a.inside_strictly(&ball((0.0, 3.0), (2.0, 5.0))));
        assert!(!a.inside_strictly(&a));
        let d = ball((3.0, 4.0), (0.0, 0.0));
        assert!((d.modulus_lower() - 3.0).abs() < 1e-12 && (d.modulus_upper() - 4.0).abs() < 1e-12);
        let through = ball((-1.0, 1.0), (2.0, 3.0));
        assert!((through.modulus_lower() - 2.0).abs() < 1e-12);
        let s = Cball::square(0.5, -0.25, 0.1);
        assert!(s.contains(Cx::new(0.55, -0.2)));
        assert!(!s.contains(Cx::new(0.7, -0.2)));
        assert!((s.centre() - Cx::new(0.5, -0.25)).norm() < 1e-15);
        assert!((s.half_width() - 0.1).abs() < 1e-12);
    }

    #[test]
    fn polynomial_evaluation_encloses_values() {
        // f(x, y) = x^2 y - 3 y + 1, g(x, y) = x + y^2.
        let system = vec![
            vec![
                PolyTerm {
                    coef: 1.0,
                    exps: vec![2, 1],
                },
                PolyTerm {
                    coef: -3.0,
                    exps: vec![0, 1],
                },
                PolyTerm {
                    coef: 1.0,
                    exps: vec![0, 0],
                },
            ],
            vec![
                PolyTerm {
                    coef: 1.0,
                    exps: vec![1, 0],
                },
                PolyTerm {
                    coef: 1.0,
                    exps: vec![0, 2],
                },
            ],
        ];
        let x = ball((0.5, 1.5), (-0.5, 0.5));
        let y = ball((1.0, 2.0), (0.0, 1.0));
        let values = eval_polys(&system, &[x, y]);
        let jacobian = eval_jacobian(&system, &[x, y]);
        for &xr in &[0.5_f64, 1.5] {
            for &xi in &[-0.5_f64, 0.5] {
                for &yr in &[1.0_f64, 2.0] {
                    for &yi in &[0.0_f64, 1.0] {
                        let (p, q) = (Cx::new(xr, xi), Cx::new(yr, yi));
                        assert!(values[0].contains(p * p * q - 3.0 * q + 1.0));
                        assert!(values[1].contains(p + q * q));
                        assert!(jacobian[0][0].contains(2.0 * p * q));
                        assert!(jacobian[0][1].contains(p * p - 3.0));
                        assert!(jacobian[1][0].contains(Cx::new(1.0, 0.0)));
                        assert!(jacobian[1][1].contains(2.0 * q));
                    }
                }
            }
        }
        let point = [Cx::new(0.7, -0.2), Cx::new(1.3, 0.4)];
        let values = values_at(&system, &point);
        assert!(
            (values[0] - (point[0] * point[0] * point[1] - 3.0 * point[1] + 1.0)).norm() < 1e-14
        );
        let jacobian = jacobian_at(&system, &point);
        assert!((jacobian[0][0] - 2.0 * point[0] * point[1]).norm() < 1e-14);
        assert!((jacobian[1][1] - 2.0 * point[1]).norm() < 1e-14);
    }

    #[test]
    fn complex_coefficients_and_the_real_lift_agree() {
        let terms = vec![vec![
            CTerm {
                coef: Cx::new(1.0, 1.0),
                exps: vec![2],
            },
            CTerm {
                coef: Cx::new(-2.0, 0.5),
                exps: vec![0],
            },
        ]];
        let x = [ball((0.4, 0.6), (-0.1, 0.1))];
        let values = eval_terms(&terms, &x);
        let jacobian = eval_jacobian_terms(&terms, &x);
        for &re in &[0.4_f64, 0.6] {
            for &im in &[-0.1_f64, 0.1] {
                let z = Cx::new(re, im);
                assert!(values[0].contains(Cx::new(1.0, 1.0) * z * z + Cx::new(-2.0, 0.5)));
                assert!(jacobian[0][0].contains(Cx::new(2.0, 2.0) * z));
            }
        }
        let real = vec![vec![PolyTerm {
            coef: 3.0,
            exps: vec![1],
        }]];
        let lifted = real_terms(&real);
        assert_eq!(lifted.len(), 1);
        assert_eq!(lifted[0][0].coef, Cx::new(3.0, 0.0));
        let at = values_at_terms(&lifted, &[Cx::new(2.0, 1.0)]);
        assert!((at[0] - Cx::new(6.0, 3.0)).norm() < 1e-15);
    }
}
