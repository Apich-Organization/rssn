//! Factorisation of multivariate polynomials over `Q`.
//!
//! Kronecker's reduction: with `d_i` the degree of `f` in `x_i`, the
//! substitution `x_i -> t^(w_i)`, `w_1 = 1`, `w_(i+1) = w_i (d_i + 1)`,
//! maps distinct monomials of degree at most `d_i` in each variable to
//! distinct powers of `t`, so a factor `g` of `f` maps to a factor of the
//! univariate image, and `g` is recovered from its image by reading the
//! exponents in the mixed radix `(d_1 + 1, d_2 + 1, ...)`. The image is
//! factored over `Z` (Zassenhaus, [`univariate::factor`]); products of
//! subsets of its irreducible factors, by increasing size, are mapped back
//! and tried as divisors. The first divisor found is irreducible (any
//! proper factor of it would have come from a smaller subset), so it is
//! divided out as often as it divides and the process repeats on the
//! cofactor.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::Signed;
use num_traits::Zero;

use super::repr::Poly;
use super::univariate;
use crate::graph::Number;

const CAP: usize = 200_000;
/// The largest univariate image degree attempted.
const MAX_IMAGE_DEGREE: u64 = 3000;
/// The largest number of univariate factors recombined.
const MAX_PIECES: usize = 18;
/// The largest number of candidate subsets tried per split.
const MAX_SUBSETS: usize = 100_000;

/// Mixed-radix weights for the variables `vars` of `f`.
fn weights(
    f: &Poly,
    vars: &[u32],
) -> Option<Vec<u64>> {
    let mut out = Vec::with_capacity(vars.len());
    let mut w: u64 = 1;
    for &v in vars {
        out.push(w);
        w = w.checked_mul(u64::from(f.degree_in(v)) + 1)?;
    }
    (w <= MAX_IMAGE_DEGREE + 1).then_some(out)
}

/// `f(t^(w_1), ..., t^(w_n))` as dense rational coefficients.
fn image(
    f: &Poly,
    vars: &[u32],
    w: &[u64],
) -> Option<Vec<BigRational>> {
    let mut dense: Vec<BigRational> = Vec::new();
    for (mono, c) in f.terms() {
        let mut e: u64 = 0;
        for &(g, k) in mono {
            let i = vars.iter().position(|&v| v == g)?;
            e += u64::from(k) * w[i];
        }
        let e = usize::try_from(e).ok()?;
        if dense.len() <= e {
            dense.resize(e + 1, BigRational::zero());
        }
        dense[e] += c.to_rational()?;
    }
    Some(dense)
}

/// The polynomial whose image is `z`, read in the mixed radix `w`.
fn preimage(
    z: &[BigInt],
    vars: &[u32],
    w: &[u64],
) -> Poly {
    let mut out = Poly::zero();
    for (e, c) in z.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let mut rest = e as u64;
        let mut mono = Vec::new();
        for i in (0..vars.len()).rev() {
            let k = rest / w[i];
            rest %= w[i];
            if k > 0 {
                mono.push((vars[i], u32::try_from(k).unwrap_or(u32::MAX)));
            }
        }
        mono.sort_unstable();
        out = out.add(&Poly::monomial(mono, Number::rat(BigRational::from_integer(c.clone()))));
    }
    out
}

fn z_mul(
    a: &[BigInt],
    b: &[BigInt],
) -> Vec<BigInt> {
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

/// `g` with a positive leading coefficient.
fn normalise(g: Poly) -> Poly {
    match g.leading() {
        | Some((_, c)) if c.to_rational().is_some_and(|r| r.is_negative()) => g.neg(),
        | _ => g,
    }
}

/// Advances `c` to the next `c.len()`-subset of `0..n` in lexicographic
/// order; `false` after the last one.
fn next_combination(
    c: &mut [usize],
    n: usize,
) -> bool {
    let k = c.len();
    let mut i = k;
    while i > 0 {
        i -= 1;
        if c[i] < n - k + i {
            c[i] += 1;
            for j in i + 1..k {
                c[j] = c[j - 1] + 1;
            }
            return true;
        }
    }
    false
}

/// One irreducible factor of the non-constant `f`, or `None` when `f` is
/// irreducible (or too large to split).
fn split_off(
    f: &Poly,
    vars: &[u32],
) -> Option<Poly> {
    let w = weights(f, vars)?;
    let dense = image(f, vars, &w)?;
    let (_, factors) = univariate::factor(&dense);
    let pieces: Vec<&Vec<BigInt>> =
        factors.iter().flat_map(|(z, m)| std::iter::repeat_n(z, *m as usize)).collect();
    if pieces.len() <= 1 || pieces.len() > MAX_PIECES {
        return None;
    }
    let n = pieces.len();
    let mut tried = 0;
    let mut seen: Vec<Vec<BigInt>> = Vec::new();
    for size in 1..=n / 2 {
        let mut choice: Vec<usize> = (0..size).collect();
        loop {
            let product = choice.iter().fold(vec![BigInt::from(1)], |acc, &i| z_mul(&acc, pieces[i]));
            if !seen.contains(&product) {
                let g = preimage(&product, vars, &w);
                if g.as_constant().is_none() && f.div_exact(&g, CAP).is_some() {
                    return Some(normalise(g));
                }
                seen.push(product);
            }
            tried += 1;
            if tried > MAX_SUBSETS {
                return None;
            }
            if !next_combination(&mut choice, n) {
                break;
            }
        }
    }
    None
}

/// `f = unit * prod g_i^m_i` with irreducible `g_i` of positive leading
/// coefficient, for a polynomial with exact coefficients in the
/// generators `vars`.
///
/// `None` when the polynomial is too large for the
/// reduction (the caller keeps it unfactored).
#[must_use]
pub fn factor(
    f: &Poly,
    vars: &[u32],
) -> Option<(Number, Vec<(Poly, u32)>)> {
    let mut work = f.clone();
    let mut out: Vec<(Poly, u32)> = Vec::new();
    for _ in 0..64 {
        if let Some(c) = work.as_constant() {
            out.sort_by_key(|(g, _)| (g.total_degree(), g.len()));
            return Some((c, out));
        }
        weights(&work, vars)?;
        let g = split_off(&work, vars).unwrap_or_else(|| normalise(work.clone()));
        let mut m = 0;
        while let Some(q) = work.div_exact(&g, CAP) {
            work = q;
            m += 1;
            if work.as_constant().is_some() {
                break;
            }
        }
        if m == 0 {
            return None;
        }
        out.push((g, m));
    }
    None
}
