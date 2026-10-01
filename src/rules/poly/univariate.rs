//! Dense univariate polynomials over the rationals and their factorisation.
//!
//! Coefficient vectors are in ascending degree with no trailing zeros; the
//! empty vector is the zero polynomial. Factorisation over `Q` follows the
//! classical route: make the polynomial primitive over `Z`, split off
//! repeated factors (Yun), and factor each square-free part with
//! Berlekamp–Zassenhaus — factor modulo a small prime (distinct-degree and
//! Cantor–Zassenhaus equal-degree splitting), Hensel-lift the modular
//! factors past the Landau–Mignotte bound, and recombine them into true
//! factors by trial division.

use num_bigint::BigInt;
use num_bigint::BigUint;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

/// Coefficients over `Q`, ascending degree.
pub type QPoly = Vec<BigRational>;
/// Coefficients over `Z`, ascending degree.
pub type ZPoly = Vec<BigInt>;

fn trim<T: Zero>(mut p: Vec<T>) -> Vec<T> {
    while p.last().is_some_and(Zero::is_zero) {
        p.pop();
    }
    p
}

/// Degree, or `None` for the zero polynomial.
#[must_use]
pub const fn degree<T>(p: &[T]) -> Option<usize> {
    p.len().checked_sub(1)
}

/// Sum of two polynomials over `Q`.
#[must_use]
pub fn add(
    a: &[BigRational],
    b: &[BigRational],
) -> QPoly {
    let mut out: QPoly = a.to_vec();
    out.resize(a.len().max(b.len()), BigRational::zero());
    for (slot, x) in out.iter_mut().zip(b) {
        *slot += x;
    }
    trim(out)
}

/// `a - b` over `Q`.
#[must_use]
pub fn sub(
    a: &[BigRational],
    b: &[BigRational],
) -> QPoly {
    let negated: QPoly = b.iter().map(|x| -x).collect();
    add(a, &negated)
}

/// Product over `Q`.
#[must_use]
pub fn mul(
    a: &[BigRational],
    b: &[BigRational],
) -> QPoly {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![BigRational::zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    trim(out)
}

/// Quotient and remainder of `a` by `b` over `Q`. `None` when `b` is zero.
#[must_use]
pub fn divrem(
    a: &[BigRational],
    b: &[BigRational],
) -> Option<(QPoly, QPoly)> {
    let lead = b.last()?;
    let mut rem: QPoly = a.to_vec();
    if a.len() < b.len() {
        return Some((Vec::new(), trim(rem)));
    }
    let mut quo = vec![BigRational::zero(); a.len() - b.len() + 1];
    for k in (0..quo.len()).rev() {
        let factor = &rem[k + b.len() - 1] / lead;
        if !factor.is_zero() {
            for (j, y) in b.iter().enumerate() {
                rem[k + j] -= &factor * y;
            }
        }
        quo[k] = factor;
    }
    Some((trim(quo), trim(rem)))
}

/// Divides by the leading coefficient.
#[must_use]
pub fn monic(p: &[BigRational]) -> QPoly {
    match p.last() {
        | Some(lead) => p.iter().map(|c| c / lead).collect(),
        | None => Vec::new(),
    }
}

/// Monic greatest common divisor over `Q` (Euclid).
#[must_use]
pub fn gcd(
    a: &[BigRational],
    b: &[BigRational],
) -> QPoly {
    let (mut x, mut y) = (trim(a.to_vec()), trim(b.to_vec()));
    while !y.is_empty() {
        let Some((_, r)) = divrem(&x, &y) else {
            break;
        };
        x = y;
        y = r;
    }
    monic(&x)
}

/// Derivative.
#[must_use]
pub fn derivative(p: &[BigRational]) -> QPoly {
    p.iter()
        .enumerate()
        .skip(1)
        .map(|(i, c)| c * BigRational::from_integer(BigInt::from(i)))
        .collect()
}

/// Value at `x` (Horner).
#[must_use]
pub fn eval(
    p: &[BigRational],
    x: &BigRational,
) -> BigRational {
    p.iter().rev().fold(BigRational::zero(), |acc, c| acc * x + c)
}

/// Square-free decomposition (Yun): `p = c * prod f_i^i` with the `f_i`
/// monic, square-free and pairwise coprime. Returns `(f_i, i)` for the
/// non-trivial `f_i`.
#[must_use]
pub fn square_free(p: &[BigRational]) -> Vec<(QPoly, u32)> {
    let mut out = Vec::new();
    let f = monic(p);
    if degree(&f).is_none_or(|d| d == 0) {
        return out;
    }
    let df = derivative(&f);
    let g = gcd(&f, &df);
    let exact = |a: &[BigRational], b: &[BigRational]| divrem(a, b).map(|(q, _)| q).unwrap_or_default();
    let mut c = exact(&f, &g);
    let mut d = sub(&exact(&df, &g), &derivative(&c));
    let mut multiplicity = 1_u32;
    while degree(&c).is_some_and(|deg| deg > 0) {
        let part = gcd(&c, &d);
        if degree(&part).is_some_and(|deg| deg > 0) {
            out.push((part.clone(), multiplicity));
        }
        c = exact(&c, &part);
        d = sub(&exact(&d, &part), &derivative(&c));
        multiplicity += 1;
    }
    out
}

/// Splits a rational polynomial into its content and primitive integer
/// part with positive leading coefficient: `p = content * primitive`.
#[must_use]
pub fn primitive(p: &[BigRational]) -> (BigRational, ZPoly) {
    if p.is_empty() {
        return (BigRational::zero(), Vec::new());
    }
    let denominators = p.iter().fold(BigInt::one(), |acc, c| acc.lcm(c.denom()));
    let integers: ZPoly = p.iter().map(|c| (c * &denominators).to_integer()).collect();
    let mut content = integers.iter().fold(BigInt::zero(), |acc, c| acc.gcd(c));
    if integers.last().is_some_and(Signed::is_negative) {
        content = -content;
    }
    let primitive = integers.iter().map(|c| c / &content).collect();
    (BigRational::new(content, denominators), primitive)
}

/// Exact quotient `a / b` over `Z`, or `None` if `b` does not divide `a`.
#[must_use]
pub fn exact_div(
    a: &[BigInt],
    b: &[BigInt],
) -> Option<ZPoly> {
    let lead = b.last()?;
    if a.len() < b.len() {
        return a.is_empty().then(Vec::new);
    }
    let mut rem: ZPoly = a.to_vec();
    let mut quo = vec![BigInt::zero(); a.len() - b.len() + 1];
    for k in (0..quo.len()).rev() {
        let (factor, leftover) = rem[k + b.len() - 1].div_rem(lead);
        if !leftover.is_zero() {
            return None;
        }
        if !factor.is_zero() {
            for (j, y) in b.iter().enumerate() {
                rem[k + j] -= &factor * y;
            }
        }
        quo[k] = factor;
    }
    rem.iter().all(Zero::is_zero).then(|| trim(quo))
}

fn z_mul(
    a: &[BigInt],
    b: &[BigInt],
) -> ZPoly {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![BigInt::zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    out
}

// ----------------------------------------------------------------------
// Arithmetic modulo a small prime
// ----------------------------------------------------------------------

/// Polynomials over `F_p` for a prime `p < 2^31`, ascending degree.
type PPoly = Vec<u64>;

#[derive(Copy, Clone)]
struct Fp(u64);

impl Fp {
    const fn pow(
        self,
        mut base: u64,
        mut exp: u64,
    ) -> u64 {
        let mut acc = 1_u64;
        base %= self.0;
        while exp > 0 {
            if exp & 1 == 1 {
                acc = acc * base % self.0;
            }
            base = base * base % self.0;
            exp >>= 1;
        }
        acc
    }

    const fn inv(
        self,
        x: u64,
    ) -> u64 {
        self.pow(x, self.0 - 2)
    }

    fn reduce(
        self,
        p: &[BigInt],
    ) -> PPoly {
        let modulus = BigInt::from(self.0);
        trim(p.iter().map(|c| c.mod_floor(&modulus).to_u64().unwrap_or(0)).collect())
    }

    fn sub(
        self,
        a: &[u64],
        b: &[u64],
    ) -> PPoly {
        let mut out = a.to_vec();
        out.resize(a.len().max(b.len()), 0);
        for (slot, &y) in out.iter_mut().zip(b) {
            *slot = (*slot + self.0 - y) % self.0;
        }
        trim(out)
    }

    fn mul(
        self,
        a: &[u64],
        b: &[u64],
    ) -> PPoly {
        if a.is_empty() || b.is_empty() {
            return Vec::new();
        }
        let mut out = vec![0_u64; a.len() + b.len() - 1];
        for (i, &x) in a.iter().enumerate() {
            if x == 0 {
                continue;
            }
            for (j, &y) in b.iter().enumerate() {
                out[i + j] = (out[i + j] + x * y) % self.0;
            }
        }
        trim(out)
    }

    /// Quotient and remainder; `b` must be non-zero.
    fn divrem(
        self,
        a: &[u64],
        b: &[u64],
    ) -> (PPoly, PPoly) {
        let Some(&lead) = b.last() else {
            return (Vec::new(), a.to_vec());
        };
        let mut rem = a.to_vec();
        if a.len() < b.len() {
            return (Vec::new(), trim(rem));
        }
        let inv = self.inv(lead);
        let mut quo = vec![0_u64; a.len() - b.len() + 1];
        for k in (0..quo.len()).rev() {
            let factor = rem[k + b.len() - 1] * inv % self.0;
            if factor != 0 {
                for (j, &y) in b.iter().enumerate() {
                    rem[k + j] = (rem[k + j] + self.0 - factor * y % self.0) % self.0;
                }
            }
            quo[k] = factor;
        }
        (trim(quo), trim(rem))
    }

    fn rem(
        self,
        a: &[u64],
        b: &[u64],
    ) -> PPoly {
        self.divrem(a, b).1
    }

    fn monic(
        self,
        p: &[u64],
    ) -> PPoly {
        match p.last() {
            | Some(&lead) => {
                let inv = self.inv(lead);
                p.iter().map(|&c| c * inv % self.0).collect()
            },
            | None => Vec::new(),
        }
    }

    fn gcd(
        self,
        a: &[u64],
        b: &[u64],
    ) -> PPoly {
        let (mut x, mut y) = (a.to_vec(), b.to_vec());
        while !y.is_empty() {
            let r = self.rem(&x, &y);
            x = y;
            y = r;
        }
        self.monic(&x)
    }

    /// `(s, t)` with `s*a + t*b = 1`, for coprime `a`, `b`.
    fn bezout(
        self,
        a: &[u64],
        b: &[u64],
    ) -> Option<(PPoly, PPoly)> {
        let (mut r0, mut r1) = (a.to_vec(), b.to_vec());
        let (mut s0, mut s1): (PPoly, PPoly) = (vec![1], Vec::new());
        let (mut t0, mut t1): (PPoly, PPoly) = (Vec::new(), vec![1]);
        while !r1.is_empty() {
            let (q, r) = self.divrem(&r0, &r1);
            let s2 = self.sub(&s0, &self.mul(&q, &s1));
            let t2 = self.sub(&t0, &self.mul(&q, &t1));
            (r0, r1) = (r1, r);
            (s0, s1) = (s1, s2);
            (t0, t1) = (t1, t2);
        }
        let &[unit] = r0.as_slice() else {
            return None;
        };
        let inv = self.inv(unit);
        let scale = |p: PPoly| -> PPoly { p.into_iter().map(|c| c * inv % self.0).collect() };
        Some((scale(s0), scale(t0)))
    }

    fn derivative(
        self,
        p: &[u64],
    ) -> PPoly {
        trim(p.iter().enumerate().skip(1).map(|(i, &c)| c * (i as u64 % self.0) % self.0).collect())
    }

    /// `base^exp mod modulus`.
    fn pow_mod(
        self,
        base: &[u64],
        exp: &BigUint,
        modulus: &[u64],
    ) -> PPoly {
        let mut acc: PPoly = vec![1];
        let base = self.rem(base, modulus);
        for bit in (0..exp.bits()).rev() {
            acc = self.rem(&self.mul(&acc, &acc), modulus);
            if exp.bit(bit) {
                acc = self.rem(&self.mul(&acc, &base), modulus);
            }
        }
        acc
    }

    /// Splits a monic square-free polynomial into products of irreducibles
    /// of equal degree: returns `(product, degree)` pairs.
    fn distinct_degree(
        self,
        f: &[u64],
    ) -> Vec<(PPoly, usize)> {
        let mut out = Vec::new();
        let mut f = f.to_vec();
        let x: PPoly = vec![0, 1];
        let mut h = x.clone();
        let mut d = 0_usize;
        let p = BigUint::from(self.0);
        while degree(&f).is_some_and(|deg| deg >= 2 * (d + 1)) {
            d += 1;
            h = self.pow_mod(&h, &p, &f);
            let g = self.gcd(&self.sub(&h, &x), &f);
            if degree(&g).is_some_and(|deg| deg > 0) {
                f = self.divrem(&f, &g).0;
                h = self.rem(&h, &f);
                out.push((g, d));
            }
        }
        if let Some(deg) = degree(&f).filter(|&deg| deg > 0) {
            out.push((f, deg));
        }
        out
    }

    /// Splits a product of irreducibles of degree `d` (Cantor–Zassenhaus,
    /// odd `p`).
    fn equal_degree(
        self,
        f: &[u64],
        d: usize,
        seed: &mut u64,
        out: &mut Vec<PPoly>,
    ) {
        let n = degree(f).unwrap_or(0);
        if n <= d {
            if n > 0 {
                out.push(self.monic(f));
            }
            return;
        }
        let exponent = (BigUint::from(self.0).pow(u32::try_from(d).unwrap_or(u32::MAX)) - 1_u32) >> 1;
        loop {
            let trial: PPoly = trim(
                (0..n)
                    .map(|_| {
                        *seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
                        (*seed >> 33) % self.0
                    })
                    .collect(),
            );
            if degree(&trial).is_none_or(|deg| deg == 0) {
                continue;
            }
            let shared = self.gcd(&trial, f);
            let split = if degree(&shared).is_some_and(|deg| deg > 0) {
                shared
            } else {
                let power = self.pow_mod(&trial, &exponent, f);
                self.gcd(&self.sub(&power, &[1]), f)
            };
            if degree(&split).is_some_and(|deg| deg > 0 && deg < n) {
                let other = self.divrem(f, &split).0;
                self.equal_degree(&split, d, seed, out);
                self.equal_degree(&other, d, seed, out);
                return;
            }
        }
    }

    /// Monic irreducible factors of a monic square-free polynomial.
    fn factor(
        self,
        f: &[u64],
    ) -> Vec<PPoly> {
        let mut out = Vec::new();
        let mut seed = 0x9e37_79b9_7f4a_7c15_u64 ^ self.0;
        for (product, d) in self.distinct_degree(f) {
            self.equal_degree(&product, d, &mut seed, &mut out);
        }
        out
    }
}

// ----------------------------------------------------------------------
// Hensel lifting and recombination
// ----------------------------------------------------------------------

fn symmetric(
    c: &BigInt,
    modulus: &BigInt,
) -> BigInt {
    let r = c.mod_floor(modulus);
    if &r * 2 > *modulus { r - modulus } else { r }
}

fn lift_coeffs(p: &[u64]) -> ZPoly {
    p.iter().map(|&c| BigInt::from(c)).collect()
}

/// Lifts `f = g*h (mod p)` to a factorisation modulo `p^k >= target`,
/// with `g` monic and `h` carrying the leading coefficient of `f`.
fn hensel_pair(
    fp: Fp,
    f: &[BigInt],
    g: &[u64],
    h: &[u64],
    target: &BigInt,
) -> Option<(ZPoly, ZPoly)> {
    let (_, t) = fp.bezout(g, h)?;
    let p = BigInt::from(fp.0);
    let mut big_g = lift_coeffs(g);
    let mut big_h = lift_coeffs(h);
    // Pin the leading coefficient of h to that of f so that f - g*h has
    // lower degree than f at every step.
    if let (Some(slot), Some(lead)) = (big_h.last_mut(), f.last()) {
        slot.clone_from(lead);
    }
    let mut modulus = p.clone();
    while modulus < *target {
        let product = z_mul(&big_g, &big_h);
        let mut error: ZPoly = f.to_vec();
        error.resize(error.len().max(product.len()), BigInt::zero());
        for (slot, y) in error.iter_mut().zip(&product) {
            *slot -= y;
        }
        let e = fp.reduce(&error.iter().map(|c| c / &modulus).collect::<Vec<_>>());
        if !e.is_empty() {
            // Solve dg*h + dh*g = e (mod p) with deg dg < deg g.
            let dg = fp.rem(&fp.mul(&t, &e), g);
            let (dh, leftover) = fp.divrem(&fp.sub(&e, &fp.mul(&dg, h)), g);
            if !leftover.is_empty() {
                return None;
            }
            for (slot, &c) in big_g.iter_mut().zip(&dg) {
                *slot += &modulus * c;
            }
            big_h.resize(big_h.len().max(dh.len()), BigInt::zero());
            for (slot, &c) in big_h.iter_mut().zip(&dh) {
                *slot += &modulus * c;
            }
        }
        modulus *= &p;
    }
    Some((big_g, big_h))
}

/// Lifts all modular factors of `f` to monic factors modulo a power of `p`
/// that is at least `target`. Returns the lifted factors and that modulus.
fn hensel_all(
    fp: Fp,
    f: &[BigInt],
    factors: &[PPoly],
    target: &BigInt,
) -> Option<(Vec<ZPoly>, BigInt)> {
    let p = BigInt::from(fp.0);
    let mut modulus = p.clone();
    while modulus < *target {
        modulus *= &p;
    }
    let mut lifted = Vec::with_capacity(factors.len());
    let mut remaining: ZPoly = f.to_vec();
    for (index, g) in factors.iter().enumerate() {
        let rest = factors.get(index + 1..).unwrap_or(&[]);
        if rest.is_empty() {
            // The last factor is what is left, made monic.
            let lead = remaining.last()?.clone();
            let inverse = lead.extended_gcd(&modulus).x.mod_floor(&modulus);
            lifted.push(remaining.iter().map(|c| symmetric(&(c * &inverse), &modulus)).collect());
            break;
        }
        let lead = fp.reduce(&[remaining.last()?.clone()]);
        let mut h = lead;
        for other in rest {
            h = fp.mul(&h, other);
        }
        let (big_g, big_h) = hensel_pair(fp, &remaining, g, &h, target)?;
        lifted.push(big_g.iter().map(|c| symmetric(c, &modulus)).collect());
        remaining = big_h.iter().map(|c| symmetric(c, &modulus)).collect();
    }
    Some((lifted, modulus))
}

fn z_primitive(p: &[BigInt]) -> ZPoly {
    let mut content = p.iter().fold(BigInt::zero(), |acc, c| acc.gcd(c));
    if content.is_zero() {
        return Vec::new();
    }
    if p.last().is_some_and(Signed::is_negative) {
        content = -content;
    }
    p.iter().map(|c| c / &content).collect()
}

const PRIMES: [u64; 24] =
    [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97];

/// Irreducible factors over `Z` of a primitive square-free polynomial with
/// positive leading coefficient.
fn factor_square_free(f: &[BigInt]) -> Vec<ZPoly> {
    let n = degree(f).unwrap_or(0);
    if n <= 1 {
        return if n == 1 { vec![f.to_vec()] } else { Vec::new() };
    }
    let Some(lead) = f.last() else {
        return Vec::new();
    };
    // Choose the prime giving the fewest modular factors among a few
    // admissible ones: fewer factors, fewer recombinations.
    let mut best: Option<(Fp, Vec<PPoly>)> = None;
    let mut tried = 0;
    for &p in &PRIMES {
        let fp = Fp(p);
        if (lead % p).is_zero() {
            continue;
        }
        let reduced = fp.reduce(f);
        let shared = fp.gcd(&reduced, &fp.derivative(&reduced));
        if degree(&shared) != Some(0) {
            continue;
        }
        let factors = fp.factor(&fp.monic(&reduced));
        if best.as_ref().is_none_or(|(_, old)| factors.len() < old.len()) {
            best = Some((fp, factors));
        }
        tried += 1;
        if tried == 5 || best.as_ref().is_some_and(|(_, fs)| fs.len() == 1) {
            break;
        }
    }
    let Some((fp, modular)) = best else {
        return vec![f.to_vec()];
    };
    if modular.len() == 1 {
        return vec![f.to_vec()];
    }
    // Landau–Mignotte: a factor's coefficients are bounded by 2^n * |f|_2;
    // candidates are scaled by the leading coefficient.
    let norm_squared: BigInt = f.iter().map(|c| c * c).sum();
    let bound = (norm_squared.sqrt() + 1) * (BigInt::one() << n) * lead.abs() * 2;
    let Some((lifted, modulus)) = hensel_all(fp, f, &modular, &bound) else {
        return vec![f.to_vec()];
    };

    let mut out = Vec::new();
    let mut work: ZPoly = f.to_vec();
    let mut pool: Vec<ZPoly> = lifted;
    let mut size = 1;
    while 2 * size <= pool.len() {
        let mut found = None;
        let mut choice: Vec<usize> = (0..size).collect();
        'subsets: loop {
            let work_lead = work.last().cloned().unwrap_or_else(BigInt::one);
            let mut candidate: ZPoly = vec![work_lead.clone()];
            for &i in &choice {
                candidate = z_mul(&candidate, &pool[i]).iter().map(|c| symmetric(c, &modulus)).collect();
            }
            let candidate = z_primitive(&candidate);
            let scaled: ZPoly = work.iter().map(|c| c * &work_lead).collect();
            if !candidate.is_empty() && exact_div(&scaled, &candidate).is_some() {
                found = Some((choice.clone(), candidate));
                break 'subsets;
            }
            // Next subset of `size` indices in lexicographic order.
            let mut i = size;
            loop {
                if i == 0 {
                    break 'subsets;
                }
                i -= 1;
                if choice[i] != i + pool.len() - size {
                    break;
                }
            }
            choice[i] += 1;
            for j in i + 1..size {
                choice[j] = choice[j - 1] + 1;
            }
        }
        match found {
            | Some((choice, factor)) => {
                if let Some(quotient) = exact_div(&work, &factor) {
                    work = quotient;
                } else {
                    // The candidate divides lead*work; the primitive part of
                    // the quotient is what remains.
                    let lead_now = work.last().cloned().unwrap_or_else(BigInt::one);
                    let scaled: ZPoly = work.iter().map(|c| c * &lead_now).collect();
                    work = z_primitive(&exact_div(&scaled, &factor).unwrap_or_default());
                }
                out.push(factor);
                for &i in choice.iter().rev() {
                    pool.remove(i);
                }
            },
            | None => size += 1,
        }
    }
    if degree(&work).is_some_and(|d| d > 0) {
        out.push(z_primitive(&work));
    }
    out
}

/// Complete factorisation over `Q`.
///
/// `p = content * prod factor_i^e_i` with
/// each factor primitive over `Z`, irreducible, and with positive leading
/// coefficient. Factors are ordered by degree, then by coefficients.
#[must_use]
pub fn factor(p: &[BigRational]) -> (BigRational, Vec<(ZPoly, u32)>) {
    let (content, _) = primitive(p);
    let mut out = Vec::new();
    for (part, multiplicity) in square_free(p) {
        let (_, integer_part) = primitive(&part);
        for factor in factor_square_free(&integer_part) {
            out.push((factor, multiplicity));
        }
    }
    out.sort_by(|(a, ea), (b, eb)| a.len().cmp(&b.len()).then_with(|| a.cmp(b)).then_with(|| ea.cmp(eb)));
    // The content of the input relative to the product of primitive
    // factors: compare leading coefficients.
    let lead_of_product = out.iter().fold(BigInt::one(), |acc, (f, e)| {
        acc * num_traits::pow(f.last().cloned().unwrap_or_else(BigInt::one), *e as usize)
    });
    let unit = match p.last() {
        | Some(lead) => lead / BigRational::from_integer(lead_of_product),
        | None => content,
    };
    (unit, out)
}

/// The rational roots of `p`, each once, ascending.
#[must_use]
pub fn rational_roots(p: &[BigRational]) -> Vec<BigRational> {
    let mut roots: Vec<BigRational> = factor(p)
        .1
        .into_iter()
        .filter_map(|(f, _)| match f.as_slice() {
            | [c0, c1] => Some(BigRational::new(-c0.clone(), c1.clone())),
            | _ => None,
        })
        .collect();
    roots.sort();
    roots.dedup();
    roots
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(coeffs: &[i64]) -> QPoly {
        coeffs.iter().map(|&c| BigRational::from_integer(BigInt::from(c))).collect()
    }

    fn z(coeffs: &[i64]) -> ZPoly {
        coeffs.iter().map(|&c| BigInt::from(c)).collect()
    }

    fn product(
        unit: &BigRational,
        factors: &[(ZPoly, u32)],
    ) -> QPoly {
        let mut out: QPoly = vec![unit.clone()];
        for (f, e) in factors {
            let fq: QPoly = f.iter().map(|c| BigRational::from_integer(c.clone())).collect();
            for _ in 0..*e {
                out = mul(&out, &fq);
            }
        }
        out
    }

    #[test]
    fn division_and_gcd() {
        // (x^2 - 1) / (x - 1) = x + 1
        assert_eq!(divrem(&q(&[-1, 0, 1]), &q(&[-1, 1])), Some((q(&[1, 1]), Vec::new())));
        assert_eq!(divrem(&q(&[1, 0, 1]), &q(&[-1, 1])), Some((q(&[1, 1]), q(&[2]))));
        assert_eq!(divrem(&q(&[1]), &[]), None);
        // gcd((x-1)(x+2), (x-1)(x+3)) = x - 1
        assert_eq!(gcd(&mul(&q(&[-1, 1]), &q(&[2, 1])), &mul(&q(&[-1, 1]), &q(&[3, 1]))), q(&[-1, 1]));
        assert_eq!(gcd(&q(&[1, 1]), &q(&[2, 1])), q(&[1]));
    }

    #[test]
    fn square_free_decomposition() {
        // (x - 1)^2 * (x + 2)^3 * (x + 5)
        let a = q(&[-1, 1]);
        let b = q(&[2, 1]);
        let c = q(&[5, 1]);
        let p = mul(&mul(&mul(&a, &a), &mul(&mul(&b, &b), &b)), &c);
        assert_eq!(square_free(&p), vec![(c, 1), (a, 2), (b, 3)]);
        assert!(square_free(&q(&[7])).is_empty());
    }

    #[test]
    fn factorisation_reproduces_the_input() {
        let cases: [&[i64]; 9] = [
            &[-1, 0, 1],                    // x^2 - 1
            &[1, 0, 1],                     // x^2 + 1 (irreducible)
            &[-1, 0, 0, 0, 1],              // x^4 - 1
            &[1, 0, 0, 0, 1],               // x^4 + 1 (irreducible, splits mod every p)
            &[-2, 0, 1],                    // x^2 - 2
            &[6, -5, -2, 1],                // (x-1)(x+2)(x-3)
            &[1, 1, 1, 1, 1, 1, 1],         // cyclotomic 7
            &[-1, 0, 0, 0, 0, 0, 1],        // x^6 - 1
            &[4, 8, 5, 1],                  // (x+1)(x+2)^2
        ];
        for coeffs in cases {
            let p = q(coeffs);
            let (unit, factors) = factor(&p);
            assert_eq!(product(&unit, &factors), p, "factors of {coeffs:?}: {factors:?}");
            for (f, _) in &factors {
                assert!(f.last().is_some_and(Signed::is_positive));
            }
        }
    }

    #[test]
    fn known_factorisations() {
        let degrees = |coeffs: &[i64]| -> Vec<(usize, u32)> {
            factor(&q(coeffs)).1.iter().map(|(f, e)| (f.len() - 1, *e)).collect()
        };
        assert_eq!(degrees(&[-1, 0, 1]), vec![(1, 1), (1, 1)]);
        assert_eq!(degrees(&[1, 0, 1]), vec![(2, 1)]);
        assert_eq!(degrees(&[1, 0, 0, 0, 1]), vec![(4, 1)]);
        assert_eq!(degrees(&[-1, 0, 0, 0, 1]), vec![(1, 1), (1, 1), (2, 1)]);
        assert_eq!(degrees(&[-1, 0, 0, 0, 0, 0, 1]), vec![(1, 1), (1, 1), (2, 1), (2, 1)]);
        assert_eq!(degrees(&[4, 8, 5, 1]), vec![(1, 1), (1, 2)]);
        // x^8 + x^4 + 1 = (x^2 + x + 1)(x^2 - x + 1)(x^4 - x^2 + 1)
        assert_eq!(degrees(&[1, 0, 0, 0, 1, 0, 0, 0, 1]), vec![(2, 1), (2, 1), (4, 1)]);
    }

    #[test]
    fn non_monic_and_rational_inputs() {
        // 6x^2 + 5x + 1 = (2x + 1)(3x + 1)
        let (unit, factors) = factor(&q(&[1, 5, 6]));
        assert_eq!(unit, BigRational::one());
        assert_eq!(factors, vec![(z(&[1, 2]), 1), (z(&[1, 3]), 1)]);
        // x^2/2 - 1/2 = 1/2 (x - 1)(x + 1)
        let half = BigRational::new(BigInt::from(1), BigInt::from(2));
        let p: QPoly = vec![-half.clone(), BigRational::zero(), half.clone()];
        let (unit, factors) = factor(&p);
        assert_eq!(unit, half);
        assert_eq!(product(&unit, &factors), p);
        // -x^2 + 1
        let (unit, factors) = factor(&q(&[1, 0, -1]));
        assert_eq!(unit, -BigRational::one());
        assert_eq!(product(&unit, &factors), q(&[1, 0, -1]));
    }

    #[test]
    fn products_of_random_factors_are_recovered() {
        // Multiply irreducibles together and check the factor multiset.
        let pieces: [&[i64]; 5] = [&[3, 1], &[-7, 2], &[1, 1, 1], &[2, 0, 1], &[-3, 0, 0, 1]];
        let mut seed = 12_345_u64;
        for _ in 0..20 {
            let mut p = q(&[1]);
            let mut expected: Vec<ZPoly> = Vec::new();
            for piece in pieces {
                seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
                if (seed >> 40).is_multiple_of(2) {
                    p = mul(&p, &q(piece));
                    expected.push(z(piece));
                }
            }
            let (_, factors) = factor(&p);
            let mut got: Vec<ZPoly> = factors.into_iter().map(|(f, _)| f).collect();
            got.sort();
            expected.sort();
            assert_eq!(got, expected);
        }
    }

    #[test]
    fn roots() {
        // (2x - 1)(x + 3)(x^2 + 1)
        let p = mul(&mul(&q(&[-1, 2]), &q(&[3, 1])), &q(&[1, 0, 1]));
        let roots = rational_roots(&p);
        assert_eq!(roots, vec![
            BigRational::from_integer(BigInt::from(-3)),
            BigRational::new(BigInt::from(1), BigInt::from(2))
        ]);
        assert!(roots.iter().all(|r| eval(&p, r).is_zero()));
    }
}
