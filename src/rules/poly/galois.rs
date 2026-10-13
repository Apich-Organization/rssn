//! Galois groups of irreducible rational quintics and the data for their
//! solution in radicals.
//!
//! **Decision.** Dummit's theorem: an irreducible quintic is solvable by
//! radicals exactly when its sextic resolvent has a rational root. The
//! resolvent is the polynomial `Π (w - θ_i)` over the six values
//! `θ = (Σ x_k x_{k+1} - Σ x_k x_{k+2})²` taken over the cyclic orderings
//! of the roots. Here the roots are computed to about 190 digits
//! (Durand–Kerner in `f64`, then Newton in fixed-point integer arithmetic),
//! the coefficients of the resolvent are rounded to integers and verified
//! to be integers to 90 digits, and the rational roots are found exactly
//! with integer arithmetic. The group is then `F20` (the discriminant is
//! not a square), `C5` or `D5` (square discriminant), `A5` or `S5`.
//!
//! **Radicals.** For a solvable quintic with one real root `r`, order the
//! roots so that the Galois group lies in `⟨k -> k+1, k -> 2k⟩` and complex
//! conjugation is `k -> -k` (the ordering with `θ` rational). With the
//! depressed monic integer form `z^5 + P z^3 + Q z^2 + R z + S` the
//! Lagrange resolvents `t_j = Σ ζ^{jk} z_k` are then real and
//! `r = (t_1 + t_2 + t_3 + t_4)/5`. The numbers `R_j = t_j^5` are the roots
//! of a quartic `Φ` with coefficients in `Q(√5)`; they are recognised from
//! their numerical values by lattice reduction (LLL) on `(c, 1, √5)`. The
//! real root is then `r = Σ_{Φ(R)=0} R^{1/5} / 5` with real fifth roots.
//! Solvable quintics with five real roots are classified but their roots
//! are left as `root_of` objects.

use num_bigint::BigInt;
use num_integer::Integer;
use num_integer::Roots;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

/// Precision of the fixed-point arithmetic, in bits.
const BITS: usize = 640;

/// The Galois group of an irreducible rational quintic.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum QuinticGroup {
    /// Cyclic of order 5 or dihedral of order 10 (solvable, square discriminant).
    C5OrD5,
    /// The Frobenius group of order 20 (solvable).
    F20,
    /// `A5` (not solvable by radicals).
    A5,
    /// `S5` (not solvable by radicals).
    S5,
}

impl QuinticGroup {
    /// Whether the equation is solvable by radicals.
    #[must_use]
    pub const fn solvable(self) -> bool {
        matches!(self, Self::C5OrD5 | Self::F20)
    }
}

/// Data for writing the real root of a solvable quintic with one real
/// root in radicals.
#[derive(Clone, Debug)]
pub struct Radical {
    /// `-c4 / (5 c5)`: the shift of the depression.
    pub shift: BigRational,
    /// The scale `λ` of the monic integer form (`z = λ y`).
    pub scale: BigInt,
    /// The coefficients `c0 … c3` of `Φ(w) = w^4 + c3 w^3 + … + c0` as
    /// pairs `(u, v)` meaning `u + v √5`.
    pub phi: Vec<(BigRational, BigRational)>,
}

/// The analysis of an irreducible quintic.
#[derive(Clone, Debug)]
pub struct Quintic {
    /// The Galois group.
    pub group: QuinticGroup,
    /// Radical data, when the group is solvable and there is one real root.
    pub radical: Option<Radical>,
}

#[derive(Clone, Debug)]
struct Fx {
    re: BigInt,
    im: BigInt,
}

fn one() -> BigInt {
    BigInt::one() << BITS
}

impl Fx {
    fn zero() -> Self {
        Self { re: BigInt::zero(), im: BigInt::zero() }
    }

    fn real(re: BigInt) -> Self {
        Self { re, im: BigInt::zero() }
    }

    fn integer(n: &BigInt) -> Self {
        Self::real(n << BITS)
    }

    fn add(&self, o: &Self) -> Self {
        Self { re: &self.re + &o.re, im: &self.im + &o.im }
    }

    fn sub(&self, o: &Self) -> Self {
        Self { re: &self.re - &o.re, im: &self.im - &o.im }
    }

    fn mul(&self, o: &Self) -> Self {
        Self { re: (&self.re * &o.re - &self.im * &o.im) >> BITS, im: (&self.re * &o.im + &self.im * &o.re) >> BITS }
    }

    fn conj(&self) -> Self {
        Self { re: self.re.clone(), im: -self.im.clone() }
    }

    fn div(&self, o: &Self) -> Option<Self> {
        let d = &o.re * &o.re + &o.im * &o.im;
        if d.is_zero() {
            return None;
        }
        let re = (&self.re * &o.re + &self.im * &o.im) << BITS;
        let im = (&self.im * &o.re - &self.re * &o.im) << BITS;
        Some(Self { re: re / &d, im: im / &d })
    }

    fn close(&self, o: &Self, exponent: usize) -> bool {
        let bound = BigInt::one() << (BITS - exponent);
        (&self.re - &o.re).abs() < bound && (&self.im - &o.im).abs() < bound
    }

    fn from_f64(re: f64, im: f64) -> Option<Self> {
        let scale = BigRational::from_integer(one());
        let to = |v: f64| BigRational::from_float(v).map(|r| (r * &scale).to_integer());
        Some(Self { re: to(re)?, im: to(im)? })
    }
}

fn pow5(z: &Fx) -> Fx {
    let z2 = z.mul(z);
    z2.mul(&z2).mul(z)
}

/// Value of the integer polynomial `g` (ascending) at `z`.
fn horner(g: &[BigInt], z: &Fx) -> Fx {
    g.iter().rev().fold(Fx::zero(), |acc, c| acc.mul(z).add(&Fx::integer(c)))
}

fn derivative_of(g: &[BigInt]) -> Vec<BigInt> {
    g.iter().enumerate().skip(1).map(|(i, c)| c * BigInt::from(i)).collect()
}

/// `p(y + h)` for ascending coefficients.
fn taylor_shift(
    p: &[BigRational],
    h: &BigRational,
) -> Vec<BigRational> {
    let mut work = p.to_vec();
    let n = work.len();
    for i in 0..n {
        for j in (i..n.saturating_sub(1)).rev() {
            let next = work.get(j + 1).cloned().unwrap_or_else(BigRational::zero);
            if let Some(slot) = work.get_mut(j) {
                *slot = &*slot + h * next;
            }
        }
    }
    work
}

/// Complex roots of the monic integer polynomial by Durand–Kerner.
fn approximate_roots(g: &[BigInt]) -> Option<Vec<(f64, f64)>> {
    let coeffs: Vec<f64> = g.iter().map(ToPrimitive::to_f64).collect::<Option<_>>()?;
    durand_kerner(&coeffs)
}

/// Complex roots of a polynomial (ascending `f64` coefficients, any leading
/// coefficient) by the Durand–Kerner iteration.
pub(crate) fn durand_kerner(raw: &[f64]) -> Option<Vec<(f64, f64)>> {
    let n = raw.len().checked_sub(1).filter(|&n| n >= 1)?;
    let lead = *raw.last()?;
    if lead == 0.0 || !lead.is_finite() {
        return None;
    }
    let coeffs: Vec<f64> = raw.iter().map(|c| c / lead).collect();
    let bound = 1.0
        + (0..n)
            .map(|i| coeffs[i].abs().powf(1.0 / (n - i) as f64))
            .fold(0.0_f64, f64::max)
            * 2.0;
    let mut z: Vec<(f64, f64)> = (0..n)
        .map(|k| {
            let angle = 2.0 * std::f64::consts::PI * k as f64 / n as f64 + 0.4;
            (bound * angle.cos(), bound * angle.sin())
        })
        .collect();
    let mul = |a: (f64, f64), b: (f64, f64)| (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0);
    let div = |a: (f64, f64), b: (f64, f64)| {
        let d = b.0 * b.0 + b.1 * b.1;
        ((a.0 * b.0 + a.1 * b.1) / d, (a.1 * b.0 - a.0 * b.1) / d)
    };
    for _ in 0..2000 {
        let mut moved = 0.0_f64;
        for i in 0..n {
            let zi = z[i];
            let value = coeffs.iter().rev().fold((0.0, 0.0), |acc, &c| {
                let m = mul(acc, zi);
                (m.0 + c, m.1)
            });
            let mut denominator = (1.0, 0.0);
            for (j, &zj) in z.iter().enumerate() {
                if j != i {
                    denominator = mul(denominator, (zi.0 - zj.0, zi.1 - zj.1));
                }
            }
            if denominator == (0.0, 0.0) {
                return None;
            }
            let step = div(value, denominator);
            z[i] = (zi.0 - step.0, zi.1 - step.1);
            moved = moved.max(step.0.abs() + step.1.abs());
        }
        if moved < 1e-14 * bound {
            break;
        }
    }
    z.iter().all(|r| r.0.is_finite() && r.1.is_finite()).then_some(z)
}

/// Roots refined to `BITS` bits.
fn refined_roots(g: &[BigInt]) -> Option<Vec<Fx>> {
    let derivative = derivative_of(g);
    let mut out = Vec::new();
    for (re, im) in approximate_roots(g)? {
        let mut z = Fx::from_f64(re, im)?;
        for _ in 0..12 {
            let step = horner(g, &z).div(&horner(&derivative, &z))?;
            z = z.sub(&step);
        }
        out.push(z);
    }
    Some(out)
}

fn permutations(items: &[usize]) -> Vec<Vec<usize>> {
    if items.len() <= 1 {
        return vec![items.to_vec()];
    }
    let mut out = Vec::new();
    for (i, &first) in items.iter().enumerate() {
        let mut rest = items.to_vec();
        rest.remove(i);
        for mut tail in permutations(&rest) {
            tail.insert(0, first);
            out.push(tail);
        }
    }
    out
}

/// `θ` of the cyclic ordering `x`.
fn theta(x: &[Fx]) -> Option<Fx> {
    let at = |i: usize| x.get(i % 5);
    let mut q = Fx::zero();
    for k in 0..5 {
        let edge = at(k)?.mul(at(k + 1)?);
        let diagonal = at(k)?.mul(at(k + 2)?);
        q = q.add(&edge).sub(&diagonal);
    }
    Some(q.mul(&q))
}

fn is_square(n: &BigInt) -> bool {
    !n.is_negative() && {
        let r = n.sqrt();
        &r * &r == *n
    }
}

/// Integer value of the real fixed-point `v` when it is within `2^-90` of one.
fn round_integer(v: &Fx) -> Option<BigInt> {
    let half = one() >> 1;
    let n = (&v.re + half) >> BITS;
    let near = Fx::integer(&n);
    v.close(&near, 300).then_some(n)
}

fn identify(
    c: &Fx,
    sqrt5: &BigInt,
) -> Option<(BigRational, BigRational)> {
    let rows = vec![
        vec![BigInt::one(), BigInt::zero(), BigInt::zero(), c.re.clone()],
        vec![BigInt::zero(), BigInt::one(), BigInt::zero(), one()],
        vec![BigInt::zero(), BigInt::zero(), BigInt::one(), sqrt5.clone()],
    ];
    let reduced = lll(rows)?;
    let best = reduced.into_iter().find(|r| r.first().is_some_and(|m| !m.is_zero()) && r.get(3).is_some_and(|e| e.abs() < BigInt::one() << (BITS / 2)))?;
    let (m0, m1, m2) = (best.first()?.clone(), best.get(1)?.clone(), best.get(2)?.clone());
    let u = BigRational::new(-m1, m0.clone());
    let v = BigRational::new(-m2, m0);
    // Verify.
    let lhs = BigRational::from_integer(c.re.clone());
    let rhs = &u * BigRational::from_integer(one()) + &v * BigRational::from_integer(sqrt5.clone());
    ((lhs - rhs).abs() < BigRational::from_integer(BigInt::one() << (BITS / 2))).then_some((u, v))
}

/// LLL reduction (δ = 3/4) of integer row vectors.
fn lll(mut b: Vec<Vec<BigInt>>) -> Option<Vec<Vec<BigInt>>> {
    let n = b.len();
    let dot = |x: &[BigRational], y: &[BigRational]| -> BigRational { x.iter().zip(y).map(|(p, q)| p * q).sum() };
    let to_q = |v: &[BigInt]| -> Vec<BigRational> { v.iter().cloned().map(BigRational::from_integer).collect() };
    let delta = BigRational::new(BigInt::from(3), BigInt::from(4));
    let gram_schmidt = |b: &[Vec<BigInt>]| -> (Vec<Vec<BigRational>>, Vec<Vec<BigRational>>) {
        let mut star: Vec<Vec<BigRational>> = Vec::new();
        let mut mu = vec![vec![BigRational::zero(); b.len()]; b.len()];
        for i in 0..b.len() {
            let mut v = to_q(&b[i]);
            for j in 0..i {
                let denominator = dot(&star[j], &star[j]);
                let m = if denominator.is_zero() { BigRational::zero() } else { dot(&to_q(&b[i]), &star[j]) / denominator };
                for (a, s) in v.iter_mut().zip(&star[j]) {
                    *a -= &m * s;
                }
                mu[i][j] = m;
            }
            star.push(v);
        }
        (star, mu)
    };
    let mut k = 1;
    let mut guard = 0;
    while k < n {
        guard += 1;
        if guard > 100_000 {
            return None;
        }
        let (_, mu) = gram_schmidt(&b);
        for j in (0..k).rev() {
            let (_, mu_now) = gram_schmidt(&b);
            let q = (&mu_now[k][j] + BigRational::new(BigInt::one(), BigInt::from(2))).floor().to_integer();
            if !q.is_zero() {
                let row = b[j].clone();
                for (a, r) in b[k].iter_mut().zip(&row) {
                    *a -= &q * r;
                }
            }
        }
        let _ = mu;
        let (star, mu) = gram_schmidt(&b);
        let lhs = dot(&star[k], &star[k]);
        let rhs = (&delta - &mu[k][k - 1] * &mu[k][k - 1]) * dot(&star[k - 1], &star[k - 1]);
        if lhs >= rhs {
            k += 1;
        } else {
            b.swap(k, k - 1);
            k = k.saturating_sub(1).max(1);
        }
    }
    Some(b)
}

/// Analyses the irreducible quintic with rational coefficients `c`
/// (ascending, degree five).
#[must_use]
pub fn analyse_quintic(c: &[BigRational]) -> Option<Quintic> {
    if c.len() != 6 {
        return None;
    }
    let lead = c.last()?.clone();
    let monic: Vec<BigRational> = c.iter().map(|a| a / &lead).collect();
    let shift = -monic.get(4)? / BigRational::from_integer(BigInt::from(5));
    let depressed = taylor_shift(&monic, &shift);
    // z = λ y, integer coefficients.
    let mut scale = BigInt::one();
    for a in depressed.iter().take(4) {
        scale = scale.lcm(a.denom());
    }
    let power = |k: usize| -> BigInt { (0..k).fold(BigInt::one(), |acc, _| acc * &scale) };
    let g: Vec<BigInt> = depressed
        .iter()
        .enumerate()
        .map(|(i, a)| {
            let degree_gap = 5 - i;
            (a * BigRational::from_integer(power(degree_gap))).to_integer()
        })
        .collect();
    // Check integrality: the value must be exact.
    for (i, a) in depressed.iter().enumerate() {
        if (a * BigRational::from_integer(power(5 - i))).denom() != &BigInt::one() {
            return None;
        }
    }
    let roots = refined_roots(&g)?;
    // Resolvent θ values.
    let mut values: Vec<(Fx, Vec<usize>)> = Vec::new();
    for perm in permutations(&[1, 2, 3, 4]) {
        let mut order = vec![0_usize];
        order.extend(perm);
        let x: Vec<Fx> = order.iter().map(|&i| roots.get(i).cloned()).collect::<Option<_>>()?;
        let t = theta(&x)?;
        match values.iter_mut().find(|(v, _)| v.close(&t, 300)) {
            | Some((_, orders)) => orders.push(order.iter().fold(0, |acc, &i| acc * 5 + i)),
            | None => values.push((t, vec![order.iter().fold(0, |acc, &i| acc * 5 + i)])),
        }
    }
    if values.len() != 6 {
        return None;
    }
    // Sextic with integer coefficients.
    let mut sextic = vec![Fx::integer(&BigInt::one())];
    for (t, _) in &values {
        let mut next = vec![Fx::zero(); sextic.len() + 1];
        for (i, a) in sextic.iter().enumerate() {
            next[i + 1] = next[i + 1].add(a);
            next[i] = next[i].sub(&a.mul(t));
        }
        sextic = next;
    }
    let integer_sextic: Vec<BigInt> = sextic.iter().map(round_integer).collect::<Option<_>>()?;
    let mut rational_theta: Option<Fx> = None;
    for (t, _) in &values {
        if let Some(n) = round_integer(t) {
            let value = integer_sextic.iter().rev().fold(BigInt::zero(), |acc, a| acc * &n + a);
            if value.is_zero() {
                rational_theta = Some(t.clone());
                break;
            }
        }
    }
    // Discriminant: product of squared root differences.
    let mut disc = Fx::integer(&BigInt::one());
    for i in 0..5 {
        for j in i + 1..5 {
            let d = roots.get(i)?.sub(roots.get(j)?);
            disc = disc.mul(&d.mul(&d));
        }
    }
    let disc = round_integer(&disc)?;
    let square = is_square(&disc);
    let Some(rational_theta) = rational_theta else {
        return Some(Quintic { group: if square { QuinticGroup::A5 } else { QuinticGroup::S5 }, radical: None });
    };
    let group = if square { QuinticGroup::C5OrD5 } else { QuinticGroup::F20 };
    // One real root?
    let real: Vec<usize> = (0..5).filter(|&i| roots.get(i).is_some_and(|z| z.im.abs() < BigInt::one() << (BITS - 100))).collect();
    let [r] = real.as_slice() else {
        return Some(Quintic { group, radical: None });
    };
    let upper: Vec<usize> = (0..5).filter(|&i| roots.get(i).is_some_and(|z| z.im.is_positive()) && i != *r).collect();
    let [a, b] = upper.as_slice() else {
        return Some(Quintic { group, radical: None });
    };
    let (root_r, root_a, root_b) = (roots.get(*r)?.clone(), roots.get(*a)?.clone(), roots.get(*b)?.clone());
    let mut chosen: Option<Vec<Fx>> = None;
    'search: for (alpha, other) in [(&root_a, &root_b), (&root_b, &root_a)] {
        for alpha in [alpha.clone(), alpha.conj()] {
            for beta in [other.clone(), other.conj()] {
                let order = vec![root_r.clone(), alpha.clone(), beta.clone(), beta.conj(), alpha.conj()];
                if theta(&order).is_some_and(|t| t.close(&rational_theta, 300)) {
                    chosen = Some(order);
                    break 'search;
                }
            }
        }
    }
    let Some(order) = chosen else {
        return Some(Quintic { group, radical: None });
    };
    // Lagrange resolvents.
    let sqrt5 = ((BigInt::from(5)) << (2 * BITS)).sqrt();
    let cos = Fx::real((&sqrt5 - one()) >> 2);
    let sin_sq: BigInt = (BigInt::from(10) << BITS) + (&sqrt5 << 1);
    let sin = Fx::real(Roots::sqrt(&(sin_sq << BITS)) >> 2);
    let zeta = Fx { re: cos.re, im: sin.re };
    let mut zeta_powers = vec![Fx::integer(&BigInt::one())];
    for _ in 1..5 {
        let last = zeta_powers.last()?.clone();
        zeta_powers.push(last.mul(&zeta));
    }
    let mut resolvents = Vec::new();
    for j in 1..5_usize {
        let mut t = Fx::zero();
        for (k, x) in order.iter().enumerate() {
            t = t.add(&zeta_powers.get((j * k) % 5)?.mul(x));
        }
        // Real by construction.
        if t.im.abs() > BigInt::one() << (BITS - 200) {
            return Some(Quintic { group, radical: None });
        }
        resolvents.push(pow5(&t));
    }
    // Φ(w) = Π (w - R_j): elementary symmetric functions.
    let mut phi = vec![Fx::integer(&BigInt::one())];
    for r5 in &resolvents {
        let mut next = vec![Fx::zero(); phi.len() + 1];
        for (i, a) in phi.iter().enumerate() {
            next[i + 1] = next[i + 1].add(a);
            next[i] = next[i].sub(&a.mul(r5));
        }
        phi = next;
    }
    // phi is ascending: [c0, c1, c2, c3, 1].
    let mut pairs = Vec::new();
    for coefficient in phi.iter().take(4) {
        pairs.push(identify(coefficient, &sqrt5)?);
    }
    Some(Quintic { group, radical: Some(Radical { shift, scale, phi: pairs }) })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn quintic(c: &[i64]) -> Vec<BigRational> {
        c.iter().map(|&v| BigRational::from_integer(BigInt::from(v))).collect()
    }

    #[test]
    fn de_moivre_numeric_quintic() {
        let q = analyse_quintic(&quintic(&[3, 5, 0, 5, 0, 1]));
        eprintln!("{q:?}");
        assert!(q.is_some_and(|q| q.group.solvable() && q.radical.is_some()));
    }

    #[test]
    fn groups_of_known_quintics() {
        // x^5 - x - 1: S5.
        assert_eq!(analyse_quintic(&quintic(&[-1, -1, 0, 0, 0, 1])).map(|q| q.group), Some(QuinticGroup::S5));
        // x^5 - 5x + 12 (D5), x^5 + 15x + 12 (F20): solvable.
        let a = analyse_quintic(&quintic(&[12, -5, 0, 0, 0, 1]));
        assert!(a.as_ref().is_some_and(|q| q.group.solvable() && q.radical.is_some()), "{a:?}");
        let b = analyse_quintic(&quintic(&[12, 15, 0, 0, 0, 1]));
        assert!(b.as_ref().is_some_and(|q| q.group.solvable() && q.radical.is_some()), "{b:?}");
        // x^5 - 4x + 2 and x^5 - 5x + 3... are S5 / unsolvable.
        assert_eq!(analyse_quintic(&quintic(&[2, -4, 0, 0, 0, 1])).map(|q| q.group), Some(QuinticGroup::S5));
        // The cyclic quintic with three... x^5 + x^4 - 4x^3 - 3x^2 + 3x + 1 (C5, five real roots).
        let c = analyse_quintic(&quintic(&[1, 3, -3, -4, 1, 1]));
        assert_eq!(c.map(|q| q.group), Some(QuinticGroup::C5OrD5));
    }
}
