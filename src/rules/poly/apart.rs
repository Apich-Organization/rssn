//! Partial fraction decomposition over the rationals (and rational function
//! fields with symbolic parameters).
//!
//! `n / d = q + sum_i sum_{j <= m_i} a_ij / f_i^j`, where `d = c prod f_i^m_i`
//! is the factorisation of the denominator into irreducible factors over `Q`
//! and every numerator `a_ij` has degree below `deg f_i`.
//!
//! To eliminate rational `BigRational` GCD swell inside polynomial divisions
//! and Bézout inversions:
//! * Polynomial division and cofactor extraction clear denominators into
//!   primitive integer polynomials over `Z[x]`.
//! * Simple linear factors `(a x + b)` use closed-form homogeneous integer
//!   evaluation `Z_r(-b/a) / Z_d'(-b/a)`.
//! * Repeated linear factors `(a x + b)^m` use shifted Horner Taylor
//!   expansion in `u = a x + b` and fraction-free power-series division.
//! * Higher-degree factors compute `B^(-1) mod f_i^m_i` via Half-GCD modular
//!   Bézout over `F_p[x]` (`multifactor::ext_gcd_fp` with pooled scratch),
//!   CRT reconstruction, early Wang rational reconstruction, and exact integer
//!   certification, followed by balanced/sequential `f_i`-adic expansion.
//! * Multi-factor decompositions execute factor blocks concurrently via
//!   `crate::kernels::par` when `crate::execution::should_split_work` allows,
//!   preserving exact deterministic factor and power ordering.
//!
//! `split_real_quadratics` goes one step further and rewrites such a
//! decomposition over `R`, adjoining `sqrt(discriminant)` for each real
//! irreducible quadratic factor.

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::Zero;

use super::multifactor;
use super::repr::Mono;
use super::repr::Poly;
use super::univariate;
use super::univariate::QPoly;
use super::univariate::ZPoly;
use crate::execution;
use crate::graph::Number;
use crate::kernels::par;

const POLY_CAP: usize = 200_000;

/// One term `numerator / factor^power` of a decomposition.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Piece {
    /// An irreducible factor of the denominator, primitive over `Z`.
    pub factor: QPoly,
    /// Its power in this term.
    pub power: u32,
    /// A numerator of degree below `deg factor` (never zero).
    pub numerator: QPoly,
}

/// A complete decomposition of a rational function.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Apart {
    /// The polynomial part.
    pub quotient: QPoly,
    /// The proper fractions.
    pub pieces: Vec<Piece>,
}

/// A term `(real + radical*sqrt(radicand)) / (x - pole)^power` of a
/// decomposition taken over a real quadratic extension of `Q`.
///
/// The pole is `offset + scale*sqrt(radicand)`; for a term with a rational
/// numerator and a rational pole `radicand`, `radical` and `scale` are zero.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RadicalPiece {
    /// Rational part of the numerator.
    pub real: BigRational,
    /// Coefficient of `sqrt(radicand)` in the numerator.
    pub radical: BigRational,
    /// Positive rational whose square root is adjoined; zero for a term
    /// that is entirely rational.
    pub radicand: BigRational,
    /// Rational part of the pole.
    pub offset: BigRational,
    /// Coefficient of `sqrt(radicand)` in the pole.
    pub scale: BigRational,
    /// Power of the linear denominator.
    pub power: u32,
}

/// One term `(numer / denom) / factor^power` of a partial fraction
/// decomposition with symbolic parameters.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct PolyPiece {
    /// Factor polynomial in the decomposition variable (and possibly parameters).
    pub(crate) factor: Poly,
    /// Power of `factor` in this term.
    pub(crate) power: u32,
    /// Numerator polynomial (degree in `var` below `deg factor`).
    pub(crate) numer: Poly,
    /// Parameter-only denominator (free of `var`, primitive with positive leading coeff).
    pub(crate) denom: Poly,
}

/// A complete partial fraction decomposition of a multivariate rational
/// function with respect to a chosen variable `var`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct PolyApart {
    /// Numerator of the polynomial part in `var`.
    pub(crate) quotient_numer: Poly,
    /// Parameter-only denominator of the polynomial part (free of `var`).
    pub(crate) quotient_denom: Poly,
    /// Proper partial fraction terms.
    pub(crate) pieces: Vec<PolyPiece>,
}

fn trim_q(mut p: QPoly) -> QPoly {
    while p.last().is_some_and(Zero::is_zero) {
        p.pop();
    }
    p
}

fn trim_z(p: &mut ZPoly) {
    while p.last().is_some_and(Zero::is_zero) {
        p.pop();
    }
}

fn z_content_and_primitive(p: &[BigInt]) -> (BigInt, ZPoly) {
    let mut trimmed = p.to_vec();
    trim_z(&mut trimmed);
    if trimmed.is_empty() {
        return (BigInt::zero(), Vec::new());
    }
    let mut c = trimmed.iter().fold(BigInt::zero(), |acc, x| acc.gcd(x));
    if trimmed.last().is_some_and(Signed::is_negative) {
        c = -c;
    }
    if c.is_one() {
        return (c, trimmed);
    }
    let prim = trimmed.iter().map(|x| x / &c).collect();
    (c, prim)
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
        if x.is_zero() {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            if !y.is_zero() {
                out[i + j] += x * y;
            }
        }
    }
    trim_z(&mut out);
    out
}

fn z_pow(
    base: &[BigInt],
    mut exponent: u32,
) -> ZPoly {
    if exponent == 0 {
        return vec![BigInt::one()];
    }
    if base.is_empty() {
        return Vec::new();
    }
    let mut result = vec![BigInt::one()];
    let mut cur = base.to_vec();
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = z_mul(&result, &cur);
        }
        exponent >>= 1;
        if exponent > 0 {
            cur = z_mul(&cur, &cur);
        }
    }
    result
}

fn z_derivative(p: &[BigInt]) -> ZPoly {
    let mut out: ZPoly = p
        .iter()
        .enumerate()
        .skip(1)
        .map(|(i, c)| c * BigInt::from(i))
        .collect();
    trim_z(&mut out);
    out
}

fn z_sub_const(
    mut p: ZPoly,
    c: &BigInt,
) -> ZPoly {
    if p.is_empty() {
        if !c.is_zero() {
            p.push(-c);
        }
        return p;
    }
    p[0] -= c;
    trim_z(&mut p);
    p
}

/// Evaluates `a^d * P(-b / a)` for `P` of degree `d` via homogeneous Horner
/// without any rational arithmetic.
fn eval_homog_z(
    p: &[BigInt],
    num: &BigInt,
    den: &BigInt,
) -> BigInt {
    if p.is_empty() {
        return BigInt::zero();
    }
    let mut acc = BigInt::zero();
    let mut den_pow = BigInt::one();
    for (i, c) in p.iter().rev().enumerate() {
        if i == 0 {
            acc.clone_from(c);
        } else {
            den_pow *= den;
            acc = acc * num + c * &den_pow;
        }
    }
    acc
}

/// Fraction-free integer polynomial division: returns `(quo, rem, den)` with
/// `den > 0` such that `den * a == quo * b + rem` and `deg(rem) < deg(b)`.
fn divrem_z(
    a: &[BigInt],
    b: &[BigInt],
) -> Option<(ZPoly, ZPoly, BigInt)> {
    let lead = b.last()?;
    if lead.is_zero() {
        return None;
    }
    let mut a_len = a.len();
    while a_len > 0 && a[a_len - 1].is_zero() {
        a_len -= 1;
    }
    if a_len < b.len() {
        return Some((Vec::new(), a[..a_len].to_vec(), BigInt::one()));
    }
    let abs_lead = lead.abs();
    let lead_neg = lead.is_negative();
    let q_len = a_len - b.len() + 1;
    let mut rem: ZPoly = a[..a_len].to_vec();
    let mut quo = vec![BigInt::zero(); q_len];
    let mut den = BigInt::one();
    for k in (0..q_len).rev() {
        let top = rem[k + b.len() - 1].clone();
        if top.is_zero() {
            continue;
        }
        let g = top.gcd(&abs_lead);
        let step = &abs_lead / &g;
        let mut factor = &top / &g;
        if lead_neg {
            factor = -factor;
        }
        if !step.is_one() {
            for r in &mut rem[..k + b.len() - 1] {
                *r *= &step;
            }
            for q_slot in &mut quo[k + 1..] {
                *q_slot *= &step;
            }
            den *= &step;
        }
        rem[k + b.len() - 1] = BigInt::zero();
        for j in 0..(b.len() - 1) {
            if !b[j].is_zero() {
                rem[k + j] -= &factor * &b[j];
            }
        }
        quo[k] = factor;
    }
    trim_z(&mut quo);
    rem.truncate(b.len() - 1);
    trim_z(&mut rem);
    Some((quo, rem, den))
}

fn z_to_q_scaled(
    z: &[BigInt],
    scale: &BigRational,
) -> QPoly {
    if z.is_empty() || scale.is_zero() {
        return Vec::new();
    }
    let mut out: QPoly = z
        .iter()
        .map(|c| scale * BigRational::from_integer(c.clone()))
        .collect();
    while out.last().is_some_and(Zero::is_zero) {
        out.pop();
    }
    out
}

/// Wang's rational number reconstruction modulo `modulus`: recovers `(num, den)`
/// with `|num|, den <= limit` and `num / den == u (mod modulus)` when it exists.
fn rat_recon(
    u: &BigInt,
    modulus: &BigInt,
    limit: &BigInt,
) -> Option<(BigInt, BigInt)> {
    if u.is_zero() {
        return Some((BigInt::zero(), BigInt::one()));
    }
    let mut r0 = modulus.clone();
    let mut r1 = u.clone();
    let mut t0 = BigInt::zero();
    let mut t1 = BigInt::one();
    while &r1 > limit {
        let (q, r2) = r0.div_rem(&r1);
        let t2 = &t0 - &q * &t1;
        r0 = r1;
        r1 = r2;
        t0 = t1;
        t1 = t2;
    }
    if t1.is_zero() {
        return None;
    }
    let (num, den) = if t1.is_negative() {
        (-r1, -t1)
    } else {
        (r1, t1)
    };
    if &den > limit || !num.gcd(&den).is_one() {
        return None;
    }
    Some((num, den))
}

/// Attempts Wang rational reconstruction on `crt_s mod crt_mod` and certifies
/// `cand_num * a == cand_den (mod m)` via exact integer division.
fn try_rat_recon_inverse(
    crt_s: &[BigInt],
    crt_mod: &BigInt,
    a: &[BigInt],
    m: &[BigInt],
) -> Option<(ZPoly, BigInt)> {
    let halved: BigInt = crt_mod >> 1;
    let limit = halved.sqrt();
    if limit.is_zero() {
        return None;
    }
    let mut pairs = Vec::with_capacity(crt_s.len());
    let mut lcm_den = BigInt::one();
    for c in crt_s {
        let (num, den) = rat_recon(c, crt_mod, &limit)?;
        lcm_den = lcm_den.lcm(&den);
        pairs.push((num, den));
    }
    let mut cand_num: ZPoly = pairs
        .into_iter()
        .map(|(num, den)| num * (&lcm_den / &den))
        .collect();
    trim_z(&mut cand_num);
    if cand_num.is_empty() {
        return None;
    }
    let g = cand_num.iter().fold(lcm_den.clone(), |acc, x| acc.gcd(x));
    if g > BigInt::one() {
        for c in &mut cand_num {
            *c /= &g;
        }
        lcm_den /= &g;
    }
    let diff = z_sub_const(z_mul(&cand_num, a), &lcm_den);
    univariate::exact_div(&diff, m).map(|_| (cand_num, lcm_den))
}

/// Simplifies `(num_poly, den_scalar)` by dividing out common content and ensuring
/// `den_scalar > 0`.
fn simplify_z_inv(
    mut num: ZPoly,
    mut den: BigInt,
) -> Option<(ZPoly, BigInt)> {
    trim_z(&mut num);
    if den.is_zero() || num.is_empty() {
        return None;
    }
    if den.is_negative() {
        den = -den;
        for c in &mut num {
            *c = -&*c;
        }
    }
    let g = num.iter().fold(den.clone(), |acc, c| acc.gcd(c));
    if g > BigInt::one() {
        for c in &mut num {
            *c /= &g;
        }
        den /= &g;
    }
    Some((num, den))
}

/// Computes `(inv_num, inv_den)` with `inv_den > 0` such that
/// `inv_num * a == inv_den (mod m)`, where `m` is primitive over `Z` with
/// positive leading coefficient.
///
/// Uses Half-GCD extended GCD over `F_p[x]` (`multifactor::ext_gcd_fp`),
/// incremental CRT, early Wang rational reconstruction, resultant-scaled CRT
/// fallback, and exact integer verification — completely eliminating rational
/// Euclid coefficient swell.
fn modular_inverse_mod_z(
    a: &[BigInt],
    m: &[BigInt],
) -> Option<(ZPoly, BigInt)> {
    if m.len() <= 1 || a.is_empty() {
        return None;
    }
    let (c_a, a_prim) = z_content_and_primitive(a);
    if a_prim.is_empty() {
        return None;
    }
    let (_, a_rem, a_rem_den) = divrem_z(&a_prim, m)?;
    let (c_rem, a_inner) = z_content_and_primitive(&a_rem);
    if a_inner.is_empty() {
        return None;
    }
    if a_inner.len() == 1 {
        return simplify_z_inv(vec![a_rem_den], c_a * c_rem * &a_inner[0]);
    }

    let lc_m = m.last()?;
    let lc_a = a_inner.last()?;
    let bound_bits = multifactor::hadamard_resultant_bound_bits(m, &a_inner) + 32;

    let mut crt_s: Vec<BigInt> = Vec::new();
    let mut crt_mod = BigInt::zero();
    let mut crt_r: Vec<BigInt> = Vec::new();
    let mut crt_r_mod = BigInt::zero();
    let mut prime_idx = 0_usize;
    let mut good_primes = 0_usize;
    let mut bad_primes = 0_usize;

    loop {
        let p = multifactor::nth_word_prime(prime_idx);
        prime_idx += 1;
        let p_big = BigInt::from(p);
        if lc_m.mod_floor(&p_big).is_zero() || lc_a.mod_floor(&p_big).is_zero() {
            continue;
        }
        let mp = multifactor::fp_reduce_z(m, p);
        let ap = multifactor::fp_reduce_z(&a_inner, p);
        let (gp, sp, _) = multifactor::ext_gcd_fp(&ap, &mp, p);
        if gp.len() != 1 {
            bad_primes += 1;
            if bad_primes >= 2 && multifactor::gcd_z(m, &a_inner).len() > 1 {
                return None;
            }
            continue;
        }
        bad_primes = 0;

        let sp_rem = multifactor::fp_poly_rem(&sp, &mp, p);
        let rp = multifactor::resultant_fp(&mp, &ap, p);
        if rp == 0 {
            continue;
        }

        let changed_s = multifactor::crt_step_poly(&mut crt_s, &mut crt_mod, &sp_rem, p);
        let changed_r = multifactor::crt_step_poly(&mut crt_r, &mut crt_r_mod, &[rp], p);
        good_primes += 1;

        if (good_primes <= 4 || good_primes.is_power_of_two())
            && let Some((u_inner, d_inner)) = try_rat_recon_inverse(&crt_s, &crt_mod, &a_inner, m)
        {
            let num: ZPoly = u_inner.into_iter().map(|c| c * &a_rem_den).collect();
            return simplify_z_inv(num, c_a * c_rem * d_inner);
        }

        if (!changed_s && !changed_r) || (crt_mod.bits() as usize) > bound_bits {
            let half = &crt_mod >> 1;
            let r_res = multifactor::symmetric_mod(crt_r.first()?, &crt_mod, &half);
            if !r_res.is_zero() {
                let mut u_inner: ZPoly = crt_s
                    .iter()
                    .map(|s| {
                        let prod = (s * &r_res).mod_floor(&crt_mod);
                        multifactor::symmetric_mod(&prod, &crt_mod, &half)
                    })
                    .collect();
                trim_z(&mut u_inner);
                let diff = z_sub_const(z_mul(&u_inner, &a_inner), &r_res);
                if univariate::exact_div(&diff, m).is_some() {
                    let num: ZPoly = u_inner.into_iter().map(|c| c * &a_rem_den).collect();
                    return simplify_z_inv(num, c_a * c_rem * r_res);
                }
            }
        }

        if (crt_mod.bits() as usize) > bound_bits + 256 {
            return None;
        }
    }
}

/// Computes `a^{deg P} * P(x) mod (a*x + b)^m` expressed in powers of
/// `u = a*x + b` (ascending degrees `0..m`).
fn shift_linear_z(
    p: &[BigInt],
    a: &BigInt,
    b: &BigInt,
    m: usize,
) -> ZPoly {
    if p.is_empty() || m == 0 {
        return Vec::new();
    }
    let d_p = p.len() - 1;
    let mut acc = vec![BigInt::zero(); m];
    acc[0].clone_from(&p[d_p]);
    let mut cur_len = 1_usize;
    let mut a_pow = BigInt::one();
    let neg_b = -b;
    for c in p[..d_p].iter().rev() {
        let new_len = (cur_len + 1).min(m);
        for j in (1..new_len).rev() {
            let prev = acc[j - 1].clone();
            acc[j] = prev + &neg_b * &acc[j];
        }
        acc[0] *= &neg_b;
        cur_len = new_len;
        a_pow *= a;
        if !c.is_zero() {
            acc[0] += c * &a_pow;
        }
    }
    acc
}

/// Linear factor block `(a*x + b)^multiplicity` using homogeneous evaluation for
/// `m == 1` and shifted Horner Taylor power-series expansion for `m >= 2`.
fn linear_factor_block(
    z_rem: &[BigInt],
    z_denom: &[BigInt],
    z_denom_deriv: &[BigInt],
    rem_scale: &BigRational,
    spec: &FactorSpec,
) -> Option<Vec<Piece>> {
    let [b, a] = spec.f_z.as_slice() else {
        return None;
    };
    let neg_b = -b;
    if spec.multiplicity == 1 {
        if !eval_homog_z(z_denom, &neg_b, a).is_zero() {
            return None;
        }
        let v_deriv = eval_homog_z(z_denom_deriv, &neg_b, a);
        if v_deriv.is_zero() {
            return None;
        }
        if z_rem.is_empty() {
            return Some(Vec::new());
        }
        let v_rem = eval_homog_z(z_rem, &neg_b, a);
        if v_rem.is_zero() {
            return Some(Vec::new());
        }
        let d = z_denom.len() - 1;
        let d_r = z_rem.len() - 1;
        let exp_a = u32::try_from(d.checked_sub(d_r)?).ok()?;
        let coeff = rem_scale
            * &spec.f_orig_scale
            * BigRational::new(v_rem * num_traits::pow(a.clone(), exp_a as usize), v_deriv);
        if coeff.is_zero() {
            return Some(Vec::new());
        }
        return Some(vec![Piece {
            factor: spec.factor_q.clone(),
            power: 1,
            numerator: vec![coeff],
        }]);
    }

    let f_power_z = z_pow(&spec.f_z, spec.multiplicity);
    let cofactor_z = univariate::exact_div(z_denom, &f_power_z)?;
    if eval_homog_z(&cofactor_z, &neg_b, a).is_zero() {
        return None;
    }
    if z_rem.is_empty() {
        return Some(Vec::new());
    }
    let m = usize::try_from(spec.multiplicity).ok()?;
    let d_r = z_rem.len() - 1;
    let d_b = cofactor_z.len() - 1;
    let a_u = shift_linear_z(z_rem, a, b, m);
    let b_u = shift_linear_z(&cofactor_z, a, b, m);
    let b0 = b_u.first()?.clone();
    if b0.is_zero() {
        return None;
    }

    // `q_nums[k] = b0^(k + 1) * Q_u[k]` where `Q_u = A_u / B_u mod u^m`.
    let mut b0_pows = Vec::with_capacity(m + 1);
    b0_pows.push(BigInt::one());
    for _ in 0..m {
        b0_pows.push(b0_pows.last()? * &b0);
    }
    let mut q_nums: Vec<BigInt> = Vec::with_capacity(m);
    for k in 0..m {
        let mut acc = &a_u[k] * &b0_pows[k];
        for i in 1..=k {
            if !b_u[i].is_zero() && !q_nums[k - i].is_zero() {
                acc -= &b_u[i] * &b0_pows[i - 1] * &q_nums[k - i];
            }
        }
        q_nums.push(acc);
    }

    let a_factor_scale = if d_b >= d_r {
        BigRational::from_integer(num_traits::pow(a.clone(), d_b - d_r))
    } else {
        BigRational::new(BigInt::one(), num_traits::pow(a.clone(), d_r - d_b))
    };
    let base_scale = rem_scale * &a_factor_scale;
    let mut pieces = Vec::new();
    for k in (0..m).rev() {
        if q_nums[k].is_zero() {
            continue;
        }
        let power = u32::try_from(m - k).ok()?;
        let c_f_pow = if spec.f_orig_scale.is_one() {
            BigRational::one()
        } else {
            num_traits::pow(spec.f_orig_scale.clone(), power as usize)
        };
        let coeff =
            &base_scale * c_f_pow * BigRational::new(q_nums[k].clone(), b0_pows[k + 1].clone());
        if !coeff.is_zero() {
            pieces.push(Piece {
                factor: spec.factor_q.clone(),
                power,
                numerator: vec![coeff],
            });
        }
    }
    Some(pieces)
}

/// Expands `scale * local_z / spec.f_z^multiplicity` into pieces with `power`
/// in `power_offset + 1 ..= power_offset + multiplicity`, ordered ascending.
fn f_adic_expand_z(
    local_z: &[BigInt],
    local_scale: &BigRational,
    spec: &FactorSpec,
    multiplicity: u32,
    power_offset: u32,
    out: &mut Vec<Piece>,
) -> Option<()> {
    if local_z.is_empty() || local_scale.is_zero() || multiplicity == 0 {
        return Some(());
    }
    if multiplicity == 1 {
        let power = power_offset + 1;
        let c_f_pow = if spec.f_orig_scale.is_one() {
            BigRational::one()
        } else {
            num_traits::pow(spec.f_orig_scale.clone(), power as usize)
        };
        let numerator = z_to_q_scaled(local_z, &(local_scale * c_f_pow));
        if !numerator.is_empty() {
            out.push(Piece {
                factor: spec.factor_q.clone(),
                power,
                numerator,
            });
        }
        return Some(());
    }

    if multiplicity >= 8 {
        let low_m = multiplicity / 2;
        let high_m = multiplicity - low_m;
        let f_pow_low = z_pow(&spec.f_z, low_m);
        let (quo_z, rem_z, step_den) = divrem_z(local_z, &f_pow_low)?;
        let (c_q, quo_prim) = z_content_and_primitive(&quo_z);
        let (c_r, rem_prim) = z_content_and_primitive(&rem_z);
        let quo_scale = local_scale * BigRational::new(c_q, step_den.clone());
        let rem_scale = local_scale * BigRational::new(c_r, step_den);
        f_adic_expand_z(&quo_prim, &quo_scale, spec, high_m, power_offset, out)?;
        f_adic_expand_z(
            &rem_prim,
            &rem_scale,
            spec,
            low_m,
            power_offset + high_m,
            out,
        )?;
        return Some(());
    }

    let mut cur_z = local_z.to_vec();
    let mut cur_scale = local_scale.clone();
    let mut rev_pieces = Vec::new();
    for exponent in (1..=multiplicity).rev() {
        if cur_z.is_empty() || cur_scale.is_zero() {
            break;
        }
        let (quo_z, rem_z, step_den) = divrem_z(&cur_z, &spec.f_z)?;
        if !rem_z.is_empty() {
            let power = power_offset + exponent;
            let c_f_pow = if spec.f_orig_scale.is_one() {
                BigRational::one()
            } else {
                num_traits::pow(spec.f_orig_scale.clone(), power as usize)
            };
            let rem_scale =
                &cur_scale * BigRational::new(BigInt::one(), step_den.clone()) * c_f_pow;
            let numerator = z_to_q_scaled(&rem_z, &rem_scale);
            if !numerator.is_empty() {
                rev_pieces.push(Piece {
                    factor: spec.factor_q.clone(),
                    power,
                    numerator,
                });
            }
        }
        let (c_q, q_prim) = z_content_and_primitive(&quo_z);
        cur_z = q_prim;
        cur_scale *= BigRational::new(c_q, step_den);
    }
    rev_pieces.reverse();
    out.extend(rev_pieces);
    Some(())
}

/// Whether `parts` recombines to `numer / denom` exactly.
///
/// A partial fraction decomposition is a *rewrite*: if the pieces do not add
/// up to the input, the engine would replace a request with a term that is not
/// equal to it. The check costs one polynomial product per piece — far less
/// than the decomposition itself — and turns a wrong decomposition into an
/// unreduced request instead of a wrong answer.
fn recombines(
    parts: &Apart,
    numer: &[BigRational],
    denom: &[BigRational],
) -> bool {
    let mut lhs_numer = parts.quotient.clone();
    let mut lhs_denom: QPoly = vec![BigRational::one()];
    for piece in &parts.pieces {
        let mut f_pow: QPoly = vec![BigRational::one()];
        for _ in 0..piece.power {
            f_pow = univariate::mul(&f_pow, &piece.factor);
        }
        lhs_numer = univariate::add(
            &univariate::mul(&lhs_numer, &f_pow),
            &univariate::mul(&piece.numerator, &lhs_denom),
        );
        lhs_denom = univariate::mul(&lhs_denom, &f_pow);
    }
    univariate::mul(&lhs_numer, denom) == univariate::mul(numer, &lhs_denom)
}

/// Decomposes `numer / denom` over the rationals.
///
/// Input coefficients are in
/// ascending degree. The result gives the polynomial quotient and nonzero
/// fraction numerators, in factor order and increasing denominator power.
/// Returns `None` for a zero denominator, an incomplete factorisation, or a
/// decomposition that does not recombine.
#[must_use]
pub fn apart(
    numer: &[BigRational],
    denom: &[BigRational],
) -> Option<Apart> {
    let numer = trim_q(numer.to_vec());
    let denom = trim_q(denom.to_vec());
    if denom.is_empty() {
        return None;
    }
    if numer.is_empty() {
        return Some(Apart {
            quotient: Vec::new(),
            pieces: Vec::new(),
        });
    }
    if denom.len() == 1 {
        let c = &denom[0];
        let quotient = numer.into_iter().map(|a| a / c).collect();
        return Some(Apart {
            quotient,
            pieces: Vec::new(),
        });
    }

    let (c_n, z_n) = univariate::primitive(&numer);
    let (c_d, z_d) = univariate::primitive(&denom);
    let scale = c_n / &c_d;
    let (quo_z, rem_z, div_den) = divrem_z(&z_n, &z_d)?;
    let quo_scale = &scale * BigRational::new(BigInt::one(), div_den.clone());
    let quotient = z_to_q_scaled(&quo_z, &quo_scale);
    if rem_z.is_empty() {
        return Some(Apart {
            quotient,
            pieces: Vec::new(),
        });
    }

    let (_, factors) = univariate::factor(&denom);
    let (c_r, z_r) = z_content_and_primitive(&rem_z);
    let rem_scale = scale * BigRational::new(c_r, div_den);
    let factor_infos: Vec<FactorSpec> = factors
        .into_iter()
        .map(|(f_z, m)| {
            let f_q: QPoly = f_z.iter().cloned().map(BigRational::from_integer).collect();
            FactorSpec {
                factor_q: f_q,
                f_z,
                f_orig_scale: BigRational::one(),
                multiplicity: m,
            }
        })
        .collect();
    let pieces = decompose_proper_z(&z_r, &rem_scale, &z_d, &factor_infos, false)?;
    let result = Apart { quotient, pieces };
    recombines(&result, &numer, &denom).then_some(result)
}

struct FactorSpec {
    factor_q: QPoly,
    f_z: ZPoly,
    f_orig_scale: BigRational,
    multiplicity: u32,
}

fn product_tree_z(polys: &[ZPoly]) -> ZPoly {
    match polys {
        | [] => vec![BigInt::one()],
        | [only] => only.clone(),
        | _ => {
            let mid = polys.len() / 2;
            let left = product_tree_z(&polys[..mid]);
            let right = product_tree_z(&polys[mid..]);
            z_mul(&left, &right)
        },
    }
}

#[allow(dead_code)] // Used directly by unit tests for custom factor lists.
fn decompose_proper(
    remainder: &[BigRational],
    denominator: &[BigRational],
    factors: &[(QPoly, u32)],
) -> Option<Vec<Piece>> {
    let denom_trim = trim_q(denominator.to_vec());
    let degree = denom_trim.len().checked_sub(1)?;
    if degree == 0 {
        return None;
    }
    let rem_trim = trim_q(remainder.to_vec());
    let (c_d, z_d) = univariate::primitive(&denom_trim);
    let (c_r, z_r) = univariate::primitive(&rem_trim);
    let rem_scale = if c_r.is_zero() {
        BigRational::zero()
    } else {
        c_r / c_d
    };

    let mut factored_degree = 0_usize;
    let mut specs = Vec::with_capacity(factors.len());
    for (factor, multiplicity) in factors {
        let trimmed_f = trim_q(factor.clone());
        let factor_degree = trimmed_f.len().checked_sub(1)?;
        if factor_degree == 0 || *multiplicity == 0 {
            return None;
        }
        factored_degree = factored_degree
            .checked_add(factor_degree.checked_mul(usize::try_from(*multiplicity).ok()?)?)?;
        let (f_orig_scale, f_z) = univariate::primitive(&trimmed_f);
        specs.push(FactorSpec {
            factor_q: factor.clone(),
            f_z,
            f_orig_scale,
            multiplicity: *multiplicity,
        });
    }
    if factored_degree != degree {
        return None;
    }

    decompose_proper_z(&z_r, &rem_scale, &z_d, &specs, true)
}

fn decompose_proper_z(
    z_rem: &[BigInt],
    rem_scale: &BigRational,
    z_denom: &[BigInt],
    factors: &[FactorSpec],
    verify_product: bool,
) -> Option<Vec<Piece>> {
    let degree = z_denom.len().checked_sub(1)?;
    if degree == 0 {
        return None;
    }
    if verify_product {
        let powers: Vec<ZPoly> = factors
            .iter()
            .map(|spec| z_pow(&spec.f_z, spec.multiplicity))
            .collect();
        if product_tree_z(&powers) != z_denom {
            return None;
        }
    }

    let has_simple_linear = factors
        .iter()
        .any(|spec| spec.f_z.len() == 2 && spec.multiplicity == 1);
    let z_denom_deriv = if has_simple_linear {
        z_derivative(z_denom)
    } else {
        Vec::new()
    };

    let run_block = |spec: &FactorSpec| -> Option<Vec<Piece>> {
        if spec.f_z.len() == 2 {
            return linear_factor_block(z_rem, z_denom, &z_denom_deriv, rem_scale, spec);
        }
        let f_power_z = z_pow(&spec.f_z, spec.multiplicity);
        let cofactor_z = univariate::exact_div(z_denom, &f_power_z)?;
        let (inv_num_z, inv_den_z) = modular_inverse_mod_z(&cofactor_z, &f_power_z)?;
        if z_rem.is_empty() || rem_scale.is_zero() {
            return Some(Vec::new());
        }
        let prod_z = z_mul(z_rem, &inv_num_z);
        let (_, local_z, local_den_z) = divrem_z(&prod_z, &f_power_z)?;
        if local_z.is_empty() {
            return Some(Vec::new());
        }
        let (c_loc, local_prim) = z_content_and_primitive(&local_z);
        let local_scale = rem_scale * BigRational::new(c_loc, inv_den_z * local_den_z);
        let mut pieces = Vec::new();
        f_adic_expand_z(
            &local_prim,
            &local_scale,
            spec,
            spec.multiplicity,
            0,
            &mut pieces,
        )?;
        Some(pieces)
    };

    let work_per_factor = degree.saturating_mul(degree).saturating_mul(16);
    if execution::should_split_work(factors.len(), work_per_factor) {
        let blocks = par::try_map_indexed(factors.len(), 2, |index| {
            let spec = factors.get(index).ok_or(())?;
            run_block(spec).ok_or(())
        })
        .ok()?;
        return Some(blocks.into_iter().flatten().collect());
    }

    let mut pieces = Vec::new();
    for spec in factors {
        pieces.extend(run_block(spec)?);
    }
    Some(pieces)
}

// ---------------------------------------------------------------------------
// Symbolic-parameter multivariate partial fractions in one variable
// ---------------------------------------------------------------------------

type UPoly = Vec<Poly>;

fn trim_upoly(u: &mut UPoly) {
    while u.last().is_some_and(Poly::is_zero) {
        u.pop();
    }
}

fn poly_to_upoly(
    p: &Poly,
    var: u32,
) -> UPoly {
    let mut coeffs = p.coefficients_in(var);
    trim_upoly(&mut coeffs);
    coeffs
}

fn upoly_to_poly(
    u: &[Poly],
    var: u32,
) -> Option<Poly> {
    let mut acc = Poly::zero();
    let x = Poly::generator(var);
    let mut x_pow = Poly::constant(Number::from(1));
    for (k, c) in u.iter().enumerate() {
        if k > 0 {
            x_pow = x_pow.mul(&x, POLY_CAP)?;
        }
        if !c.is_zero() {
            acc = acc.add(&c.mul(&x_pow, POLY_CAP)?);
            if acc.len() > POLY_CAP {
                return None;
            }
        }
    }
    Some(acc)
}

fn upoly_mul(
    a: &[Poly],
    b: &[Poly],
) -> Option<UPoly> {
    if a.is_empty() || b.is_empty() {
        return Some(Vec::new());
    }
    let mut out = vec![Poly::zero(); a.len() + b.len() - 1];
    for (i, ca) in a.iter().enumerate() {
        if ca.is_zero() {
            continue;
        }
        for (j, cb) in b.iter().enumerate() {
            if !cb.is_zero() {
                out[i + j] = out[i + j].add(&ca.mul(cb, POLY_CAP)?);
            }
        }
    }
    trim_upoly(&mut out);
    Some(out)
}

fn upoly_sub(
    a: &[Poly],
    b: &[Poly],
) -> UPoly {
    let mut out = a.to_vec();
    out.resize(a.len().max(b.len()), Poly::zero());
    for (slot, cb) in out.iter_mut().zip(b) {
        if !cb.is_zero() {
            *slot = slot.sub(cb);
        }
    }
    trim_upoly(&mut out);
    out
}

fn upoly_scale(
    a: &[Poly],
    s: &Poly,
) -> Option<UPoly> {
    if s.is_zero() || a.is_empty() {
        return Some(Vec::new());
    }
    if s.as_constant().is_some_and(|n| n.is_one()) {
        return Some(a.to_vec());
    }
    let mut out: UPoly = a
        .iter()
        .map(|c| c.mul(s, POLY_CAP))
        .collect::<Option<Vec<_>>>()?;
    trim_upoly(&mut out);
    Some(out)
}

/// Fraction-free polynomial division over `Q[params][var]`: returns
/// `(quo, rem, scale)` with `scale` free of `var` such that
/// `scale * a == quo * b + rem`.
fn upoly_pseudo_divrem(
    a: &[Poly],
    b: &[Poly],
) -> Option<(UPoly, UPoly, Poly)> {
    let lead = b.last()?;
    if lead.is_zero() {
        return None;
    }
    let mut a_len = a.len();
    while a_len > 0 && a[a_len - 1].is_zero() {
        a_len -= 1;
    }
    if a_len < b.len() {
        return Some((
            Vec::new(),
            a[..a_len].to_vec(),
            Poly::constant(Number::from(1)),
        ));
    }
    let one_const = |p: &Poly| p.as_constant().is_some_and(|n| n.is_one());
    let q_len = a_len - b.len() + 1;
    let mut rem: UPoly = a[..a_len].to_vec();
    let mut quo = vec![Poly::zero(); q_len];
    let mut scale = Poly::constant(Number::from(1));
    for k in (0..q_len).rev() {
        let top = rem[k + b.len() - 1].clone();
        if top.is_zero() {
            continue;
        }
        let g = if one_const(lead) || one_const(&top) {
            Poly::constant(Number::from(1))
        } else {
            multifactor::poly_gcd(&top, lead)?
        };
        let (step, factor) = if one_const(&g) {
            (lead.clone(), top)
        } else {
            (lead.div_exact(&g, POLY_CAP)?, top.div_exact(&g, POLY_CAP)?)
        };
        if !one_const(&step) {
            for r in &mut rem[..k + b.len() - 1] {
                if !r.is_zero() {
                    *r = r.mul(&step, POLY_CAP)?;
                }
            }
            for q_slot in &mut quo[k + 1..] {
                if !q_slot.is_zero() {
                    *q_slot = q_slot.mul(&step, POLY_CAP)?;
                }
            }
            scale = scale.mul(&step, POLY_CAP)?;
        }
        rem[k + b.len() - 1] = Poly::zero();
        for j in 0..(b.len() - 1) {
            if !b[j].is_zero() {
                rem[k + j] = rem[k + j].sub(&factor.mul(&b[j], POLY_CAP)?);
            }
        }
        quo[k] = factor;
    }
    trim_upoly(&mut quo);
    rem.truncate(b.len() - 1);
    trim_upoly(&mut rem);
    Some((quo, rem, scale))
}

fn upoly_content(u: &[Poly]) -> Option<Poly> {
    let mut g = Poly::zero();
    for c in u.iter().filter(|c| !c.is_zero()) {
        if g.is_zero() {
            g = c.clone();
        } else {
            g = multifactor::poly_gcd(&g, c)?;
        }
        if g.as_constant().is_some_and(|n| n.is_one()) {
            return Some(g);
        }
    }
    if g.is_zero() {
        Some(Poly::constant(Number::from(1)))
    } else {
        let (_, prim) = multifactor::poly_primitive_z(&g)?;
        Some(prim)
    }
}

fn simplify_poly_pair(
    num: Poly,
    den: Poly,
) -> Option<(Poly, Poly)> {
    if den.is_zero() {
        return None;
    }
    if num.is_zero() {
        return Some((Poly::zero(), Poly::constant(Number::from(1))));
    }
    if let Some(c) = den.as_constant() {
        let inv = c.recip()?;
        return Some((num.scale(&inv), Poly::constant(Number::from(1))));
    }
    let g = multifactor::poly_gcd(&num, &den)?;
    let (num, den) = if g.as_constant().is_none() {
        (num.div_exact(&g, POLY_CAP)?, den.div_exact(&g, POLY_CAP)?)
    } else {
        (num, den)
    };
    let (c_d, den_prim) = multifactor::poly_primitive_z(&den)?;
    let inv_c = c_d.recip();
    Some((num.scale(&Number::rat(inv_c)), den_prim))
}

/// Computes `(s, scale)` with `scale` free of `var` and nonzero such that
/// `s * a == scale (mod m)` in `Frac(Q[params])[var]`.
fn upoly_inverse_mod(
    a: &[Poly],
    m: &[Poly],
) -> Option<(UPoly, Poly)> {
    if m.len() <= 1 || a.is_empty() {
        return None;
    }
    let mut r0 = m.to_vec();
    let mut r1 = a.to_vec();
    let mut t0: UPoly = Vec::new();
    let mut t1: UPoly = vec![Poly::constant(Number::from(1))];
    while r1.len() > 1 {
        let (q, mut r2, step_scale) = upoly_pseudo_divrem(&r0, &r1)?;
        let mut t2 = upoly_sub(&upoly_scale(&t0, &step_scale)?, &upoly_mul(&q, &t1)?);
        let c_r = upoly_content(&r2)?;
        let c_t = upoly_content(&t2)?;
        let common = multifactor::poly_gcd(&c_r, &c_t)?;
        if common.as_constant().is_none() {
            for slot in &mut r2 {
                if !slot.is_zero() {
                    *slot = slot.div_exact(&common, POLY_CAP)?;
                }
            }
            for slot in &mut t2 {
                if !slot.is_zero() {
                    *slot = slot.div_exact(&common, POLY_CAP)?;
                }
            }
        }
        r0 = r1;
        r1 = r2;
        t0 = t1;
        t1 = t2;
    }
    let [unit] = r1.as_slice() else {
        return None;
    };
    if unit.is_zero() {
        return None;
    }
    let (_, s_rem, rem_scale) = upoly_pseudo_divrem(&t1, m)?;
    let mut scale = unit.mul(&rem_scale, POLY_CAP)?;
    let mut s_out = s_rem;
    let c_s = upoly_content(&s_out)?;
    let common = multifactor::poly_gcd(&c_s, &scale)?;
    if common.as_constant().is_none() {
        for slot in &mut s_out {
            if !slot.is_zero() {
                *slot = slot.div_exact(&common, POLY_CAP)?;
            }
        }
        scale = scale.div_exact(&common, POLY_CAP)?;
    }
    Some((s_out, scale))
}

/// The [`recombines`] check for the multivariate path: the pieces stand for
/// `numer / (factor^power * denom)`, so the sum is accumulated over a common
/// denominator and cross-multiplied against the input.
fn recombines_poly(
    parts: &PolyApart,
    numer: &Poly,
    denom: &Poly,
) -> bool {
    let mut lhs_numer = parts.quotient_numer.clone();
    let mut lhs_denom = parts.quotient_denom.clone();
    for piece in &parts.pieces {
        let Some(f_pow) = piece.factor.pow(piece.power, POLY_CAP) else {
            return false;
        };
        let Some(piece_den) = f_pow.mul(&piece.denom, POLY_CAP) else {
            return false;
        };
        let (Some(left), Some(right)) = (
            lhs_numer.mul(&piece_den, POLY_CAP),
            piece.numer.mul(&lhs_denom, POLY_CAP),
        ) else {
            return false;
        };
        lhs_numer = left.add(&right);
        let Some(next) = lhs_denom.mul(&piece_den, POLY_CAP) else {
            return false;
        };
        lhs_denom = next;
    }
    match (
        lhs_numer.mul(denom, POLY_CAP),
        numer.mul(&lhs_denom, POLY_CAP),
    ) {
        | (Some(left), Some(right)) => left == right,
        | _ => false,
    }
}

/// Decomposes `numer / denom` with respect to `var` over `Frac(Q[params])`.
#[must_use]
pub(crate) fn apart_poly(
    numer: &Poly,
    denom: &Poly,
    var: u32,
) -> Option<PolyApart> {
    if denom.is_zero() || !numer.is_exact() || !denom.is_exact() {
        return None;
    }
    // Fast path when both numerator and denominator involve only `var`.
    if let (Some(n_u), Some(d_u)) = (numer.univariate_in(var), denom.univariate_in(var)) {
        let n_q: QPoly = n_u.iter().map(Number::to_rational).collect::<Option<_>>()?;
        let d_q: QPoly = d_u.iter().map(Number::to_rational).collect::<Option<_>>()?;
        let parts = apart(&n_q, &d_q)?;
        let q_to_poly = |q: &[BigRational]| -> Poly {
            let nums: Vec<Number> = q.iter().cloned().map(Number::rat).collect();
            Poly::from_univariate(var, &nums)
        };
        let one = Poly::constant(Number::from(1));
        return Some(PolyApart {
            quotient_numer: q_to_poly(&parts.quotient),
            quotient_denom: one.clone(),
            pieces: parts
                .pieces
                .into_iter()
                .map(|p| {
                    PolyPiece {
                        factor: q_to_poly(&p.factor),
                        power: p.power,
                        numer: q_to_poly(&p.numerator),
                        denom: one.clone(),
                    }
                })
                .collect(),
        });
    }

    let (denom_rat, denom_prim) = multifactor::poly_primitive_z(denom)?;
    let denom_mono = denom_prim.monomial_content();
    let denom_core = if denom_mono.is_empty() {
        denom_prim
    } else {
        denom_prim.div_monomial(&denom_mono)
    };
    let (unit_core, raw_factors) = if let Some(c) = denom_core.as_constant() {
        (c, Vec::new())
    } else {
        multifactor::factor(&denom_core, &denom_core.support())?
    };

    let mut const_mono: Mono = Vec::new();
    let mut var_mono_power = 0_u32;
    for &(g, e) in &denom_mono {
        if g == var {
            var_mono_power = e;
        } else {
            const_mono.push((g, e));
        }
    }
    let scale_const = Number::rat(denom_rat * unit_core.to_rational()?);
    let mut const_denom = Poly::monomial(const_mono, scale_const);
    let mut var_factors: Vec<(Poly, u32)> = Vec::new();
    if var_mono_power > 0 {
        var_factors.push((Poly::generator(var), var_mono_power));
    }
    for (g, m) in raw_factors {
        if g.degree_in(var) == 0 {
            const_denom = const_denom.mul(&g.pow(m, POLY_CAP)?, POLY_CAP)?;
        } else {
            var_factors.push((g, m));
        }
    }

    if var_factors.is_empty() {
        let (quotient_numer, quotient_denom) = simplify_poly_pair(numer.clone(), const_denom)?;
        return Some(PolyApart {
            quotient_numer,
            quotient_denom,
            pieces: Vec::new(),
        });
    }

    let mut denom_x = Poly::constant(Number::from(1));
    for (f, m) in &var_factors {
        denom_x = denom_x.mul(&f.pow(*m, POLY_CAP)?, POLY_CAP)?;
    }

    let numer_u = poly_to_upoly(numer, var);
    let denom_x_u = poly_to_upoly(&denom_x, var);
    let (q_u, r_u, q_scale) = upoly_pseudo_divrem(&numer_u, &denom_x_u)?;
    let (quotient_numer, quotient_denom) = simplify_poly_pair(
        upoly_to_poly(&q_u, var)?,
        q_scale.mul(&const_denom, POLY_CAP)?,
    )?;

    if r_u.is_empty() {
        return Some(PolyApart {
            quotient_numer,
            quotient_denom,
            pieces: Vec::new(),
        });
    }

    let base_den = q_scale.mul(&const_denom, POLY_CAP)?;
    let mut pieces = Vec::new();
    for (factor, multiplicity) in var_factors {
        let f_pow = factor.pow(multiplicity, POLY_CAP)?;
        let cofactor = denom_x.div_exact(&f_pow, POLY_CAP)?;
        let f_u = poly_to_upoly(&factor, var);
        let f_pow_u = poly_to_upoly(&f_pow, var);
        let b_u = poly_to_upoly(&cofactor, var);
        let (s_u, inv_scale) = upoly_inverse_mod(&b_u, &f_pow_u)?;
        let prod_u = upoly_mul(&r_u, &s_u)?;
        let (_, mut cur_u, loc_scale) = upoly_pseudo_divrem(&prod_u, &f_pow_u)?;
        if cur_u.is_empty() {
            continue;
        }
        let mut cur_den = base_den
            .mul(&inv_scale, POLY_CAP)?
            .mul(&loc_scale, POLY_CAP)?;
        let mut block_pieces = Vec::new();
        for power in (1..=multiplicity).rev() {
            if cur_u.is_empty() {
                break;
            }
            let (quo_u, rem_u, step_den) = upoly_pseudo_divrem(&cur_u, &f_u)?;
            let step_total_den = cur_den.mul(&step_den, POLY_CAP)?;
            if !rem_u.is_empty() {
                let rem_poly = upoly_to_poly(&rem_u, var)?;
                let (piece_num, piece_den) = simplify_poly_pair(rem_poly, step_total_den.clone())?;
                if !piece_num.is_zero() {
                    block_pieces.push(PolyPiece {
                        factor: factor.clone(),
                        power,
                        numer: piece_num,
                        denom: piece_den,
                    });
                }
            }
            cur_u = quo_u;
            cur_den = step_total_den;
        }
        block_pieces.reverse();
        pieces.extend(block_pieces);
    }

    let result = PolyApart {
        quotient_numer,
        quotient_denom,
        pieces,
    };
    recombines_poly(&result, numer, denom).then_some(result)
}

/// Rewrites the proper fractions of `parts` as terms with linear
/// denominators, adjoining `sqrt(discriminant)` for every real irreducible
/// quadratic factor.
///
/// `1/(x^2 - 2)` decomposes as `1/(x^2 - 2)` over `Q` and here becomes
/// `(sqrt(2)/4)/(x - sqrt(2)) - (sqrt(2)/4)/(x + sqrt(2))`, carried exactly as
/// `real + radical*sqrt(radicand)` with `radicand = 8`.
///
/// Returns `None`, rather than a wrong answer, when a factor has degree three
/// or more, when a quadratic factor has no real roots (it is irreducible over
/// `R` as well), or when a real quadratic factor occurs with power two or more
/// — a higher-order pole over the extension needs the full expansion there,
/// which this does not attempt. The polynomial part `parts.quotient` is
/// unchanged and is not repeated in the result.
#[must_use]
pub fn split_real_quadratics(parts: &Apart) -> Option<Vec<RadicalPiece>> {
    let mut out = Vec::new();
    for piece in &parts.pieces {
        match piece.factor.len() {
            | 2 => {
                let lead = piece.factor.get(1)?.clone();
                let constant = piece.factor.first()?.clone();
                out.push(RadicalPiece {
                    real: piece
                        .numerator
                        .first()
                        .cloned()
                        .unwrap_or_else(BigRational::zero),
                    radical: BigRational::zero(),
                    radicand: BigRational::zero(),
                    offset: -constant / lead,
                    scale: BigRational::zero(),
                    power: piece.power,
                });
            },
            | 3 => {
                if piece.power != 1 {
                    return None;
                }
                let a = piece.factor.get(2)?.clone();
                let b = piece.factor.get(1)?.clone();
                let c = piece.factor.first()?.clone();
                let two_a = BigRational::from_integer(BigInt::from(2)) * &a;
                let disc = &b * &b - BigRational::from_integer(BigInt::from(4)) * &a * &c;
                if disc <= BigRational::zero() {
                    return None;
                }
                let alpha = piece
                    .numerator
                    .get(1)
                    .cloned()
                    .unwrap_or_else(BigRational::zero);
                let beta = piece
                    .numerator
                    .first()
                    .cloned()
                    .unwrap_or_else(BigRational::zero);
                // Residues at the two roots are conjugates: `p +- q sqrt(disc)`
                // with `p = alpha / (2a)` and `q = (2a beta - alpha b)/(2a disc)`.
                let real = &alpha / &two_a;
                let radical = (&two_a * &beta - &alpha * &b) / (&two_a * &disc);
                let offset = -b / &two_a;
                let scale = BigRational::one() / &two_a;
                out.push(RadicalPiece {
                    real: real.clone(),
                    radical: radical.clone(),
                    radicand: disc.clone(),
                    offset: offset.clone(),
                    scale: scale.clone(),
                    power: 1,
                });
                out.push(RadicalPiece {
                    real,
                    radical: -radical,
                    radicand: disc,
                    offset,
                    scale: -scale,
                    power: 1,
                });
            },
            | _ => return None,
        }
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;

    use super::*;

    fn q(values: &[i64]) -> QPoly {
        values
            .iter()
            .map(|&v| BigRational::from_integer(BigInt::from(v)))
            .collect()
    }

    fn frac(
        numer: i64,
        denom: i64,
    ) -> BigRational {
        BigRational::new(BigInt::from(numer), BigInt::from(denom))
    }

    fn polynomial_power(
        factor: &[BigRational],
        exponent: u32,
    ) -> QPoly {
        match exponent {
            | 0 => vec![BigRational::one()],
            | 1 => factor.to_vec(),
            | _ => {
                let half = polynomial_power(factor, exponent / 2);
                let squared = univariate::mul(&half, &half);
                if exponent.is_multiple_of(2) {
                    squared
                } else {
                    univariate::mul(&squared, factor)
                }
            },
        }
    }

    /// Recombines a decomposition over the original denominator.
    fn recombine(
        parts: &Apart,
        denom: &[BigRational],
    ) -> QPoly {
        let mut total = univariate::mul(&parts.quotient, denom);
        for piece in &parts.pieces {
            let power = polynomial_power(&piece.factor, piece.power);
            let (cofactor, rest) = univariate::divrem(denom, &power).unwrap_or_default();
            assert_eq!(rest, []);
            total = univariate::add(&total, &univariate::mul(&cofactor, &piece.numerator));
        }
        total
    }

    /// `recombine - numer` after cross-multiplication; empty exactly when the
    /// decomposition is the original rational function.
    fn recombination_defect(
        parts: &Apart,
        numer: &[BigRational],
        denom: &[BigRational],
    ) -> QPoly {
        univariate::sub(&recombine(parts, denom), numer)
    }

    /// Taylor coefficients of `p(r + t)` up to `t^(n - 1)`, by repeated
    /// division by `x - r`.
    fn taylor(
        p: &[BigRational],
        r: &BigRational,
        n: usize,
    ) -> Vec<BigRational> {
        let divisor = vec![-r.clone(), BigRational::one()];
        let mut rest: QPoly = p.to_vec();
        let mut out = Vec::with_capacity(n);
        for _ in 0..n {
            match univariate::divrem(&rest, &divisor) {
                | Some((quotient, remainder)) => {
                    out.push(remainder.first().cloned().unwrap_or_else(BigRational::zero));
                    rest = quotient;
                },
                | None => out.push(BigRational::zero()),
            }
        }
        out
    }

    /// The first `n` coefficients of `a(t) / b(t)` at `t = 0`, `b(0) != 0`.
    fn series(
        a: &[BigRational],
        b: &[BigRational],
        n: usize,
    ) -> Vec<BigRational> {
        let lead = b.first().cloned().unwrap_or_else(BigRational::one);
        let mut out: Vec<BigRational> = Vec::with_capacity(n);
        for k in 0..n {
            let mut acc = a.get(k).cloned().unwrap_or_else(BigRational::zero);
            for i in 1..=k {
                let bi = b.get(i).cloned().unwrap_or_else(BigRational::zero);
                let si = out.get(k - i).cloned().unwrap_or_else(BigRational::zero);
                acc -= bi * si;
            }
            out.push(acc / &lead);
        }
        out
    }

    /// Multiplicity of `root` in `denom`, and `denom / (x - root)^m`.
    fn local(
        denom: &[BigRational],
        root: &BigRational,
    ) -> (u32, QPoly) {
        let divisor = vec![-root.clone(), BigRational::one()];
        let mut rest: QPoly = denom.to_vec();
        let mut power = 0;
        while rest.len() >= 2 {
            match univariate::divrem(&rest, &divisor) {
                | Some((quotient, remainder)) if remainder.is_empty() => {
                    rest = quotient;
                    power += 1;
                },
                | _ => break,
            }
        }
        (power, rest)
    }

    type Ext = (BigRational, BigRational);

    fn ext_mul(
        a: &Ext,
        b: &Ext,
        d: &BigRational,
    ) -> Ext {
        (&a.0 * &b.0 + &a.1 * &b.1 * d, &a.0 * &b.1 + &a.1 * &b.0)
    }

    fn ext_add(
        a: &Ext,
        b: &Ext,
    ) -> Ext {
        (&a.0 + &b.0, &a.1 + &b.1)
    }

    /// Asserts `c0/(x - r0) + c1/(x - r1) == numer / factor` exactly, for
    /// `factor = a (x - r0)(x - r1)`, by cross-multiplying in `Q(sqrt d)`:
    /// the sum of the conjugate pair is rational, and its cleared numerator
    /// must be the original numerator with the leading coefficient removed.
    fn quadratic_pair(
        split: &[RadicalPiece],
        numer: &[BigRational],
        factor: &[BigRational],
    ) {
        assert_eq!(split.len(), 2, "a quadratic factor splits into two terms");
        let d = split[0].radicand.clone();
        let a = factor.last().cloned().unwrap_or_else(BigRational::one);
        let coef: Vec<Ext> = split
            .iter()
            .map(|p| (p.real.clone(), p.radical.clone()))
            .collect();
        let pole: Vec<Ext> = split
            .iter()
            .map(|p| (p.offset.clone(), p.scale.clone()))
            .collect();
        // N = c0 (x - r1) + c1 (x - r0), D = (x - r0)(x - r1).
        let n_x = ext_add(&coef[0], &coef[1]);
        let n_0 = {
            let sum = ext_add(
                &ext_mul(&coef[0], &pole[1], &d),
                &ext_mul(&coef[1], &pole[0], &d),
            );
            (-sum.0, -sum.1)
        };
        let d_x = (-(&pole[0].0 + &pole[1].0), -(&pole[0].1 + &pole[1].1));
        let d_0 = ext_mul(&pole[0], &pole[1], &d);
        let zero = BigRational::zero();
        assert_eq!(n_x.1, zero, "numerator is not rational");
        assert_eq!(n_0.1, zero, "numerator is not rational");
        assert_eq!(d_x.1, zero, "denominator is not rational");
        assert_eq!(d_0.1, zero, "denominator is not rational");
        assert_eq!(
            &a * &n_x.0,
            numer.get(1).cloned().unwrap_or_else(BigRational::zero)
        );
        assert_eq!(
            &a * &n_0.0,
            numer.first().cloned().unwrap_or_else(BigRational::zero)
        );
        assert_eq!(
            d_x.0,
            &factor.get(1).cloned().unwrap_or_else(BigRational::zero) / &a
        );
        assert_eq!(
            d_0.0,
            &factor.first().cloned().unwrap_or_else(BigRational::zero) / &a
        );
    }

    #[test]
    fn decompositions_recombine() {
        let cases: [(&[i64], &[i64]); 5] = [
            (&[1], &[-1, 0, 1]),            // 1/(x^2 - 1)
            (&[0, 0, 0, 1], &[-1, 0, 1]),   // x^3/(x^2 - 1)
            (&[1, 2], &[0, 0, 1, 1]),       // (2x + 1)/(x^2 (x + 1))
            (&[3, 0, 1], &[1, 0, 2, 0, 1]), // (x^2 + 3)/(x^2 + 1)^2
            (&[1], &[-1, 0, 0, 1]),         // 1/(x^3 - 1)
        ];
        for (n, d) in cases {
            let (n, d) = (q(n), q(d));
            let parts = apart(&n, &d).unwrap_or_else(|| panic!("no decomposition"));
            assert_eq!(recombine(&parts, &d), n, "{parts:?}");
            for piece in &parts.pieces {
                assert!(piece.numerator.len() < piece.factor.len());
            }
        }
    }

    /// The cross-multiplied numerator is identically zero for every case,
    /// including repeated factors, irreducible quadratics, an irreducible
    /// cubic, a cancelling pair and an improper fraction.
    #[test]
    fn recombination_is_exactly_zero() {
        let cases: [(&[i64], &[i64]); 9] = [
            (&[1], &[-1, 0, 1]),          // 1/(x^2 - 1)
            (&[0, 0, 0, 1], &[-1, 0, 1]), // x^3/(x^2 - 1)
            (&[1, 2], &[0, 0, 1, 1]),     // (2x + 1)/(x^2 (x + 1))
            (&[1], &[0, 0, 1, 1]),        // 1/(x^2 (x + 1))
            (&[1], &[1, 2, 1]),           // 1/(x + 1)^2
            (&[1], &[-1, 0, 0, 1]),       // 1/(x^3 - 1)
            (&[1], &[1, 0, 0, 1]),        // 1/(x^3 + 1)
            (&[2], &[4, 0, 0, 0, 1]),     // 2/(x^4 + 4)
            (&[1, 1], &[2, 3, 1]),        // (x + 1)/((x + 1)(x + 2))
        ];
        for (n, d) in cases {
            let (n, d) = (q(n), q(d));
            let parts = apart(&n, &d).unwrap_or_else(|| panic!("no decomposition for {n:?}/{d:?}"));
            let defect = recombination_defect(&parts, &n, &d);
            assert!(defect.is_empty(), "{n:?}/{d:?} leaves {defect:?}");
            for piece in &parts.pieces {
                assert_ne!(piece.numerator, []);
                assert!(piece.numerator.len() < piece.factor.len());
                assert!(piece.power >= 1);
            }
        }
    }

    /// A repeated irreducible quadratic keeps both powers, each with a
    /// constant numerator: `(x^2 + 3)/(x^2 + 1)^2 = 1/(x^2 + 1) + 2/(x^2 + 1)^2`.
    #[test]
    fn repeated_irreducible_quadratic() {
        let parts = apart(&q(&[3, 0, 1]), &q(&[1, 0, 2, 0, 1])).unwrap_or_else(|| panic!("none"));
        assert_eq!(
            parts.pieces,
            vec![
                Piece {
                    factor: q(&[1, 0, 1]),
                    power: 1,
                    numerator: q(&[1]),
                },
                Piece {
                    factor: q(&[1, 0, 1]),
                    power: 2,
                    numerator: q(&[2]),
                },
            ]
        );
    }

    #[test]
    fn rational_unit_and_non_monic_repeated_factors_keep_order_and_omit_zero_pieces() {
        let linear = q(&[1, 2]);
        let quadratic = q(&[1, 0, 1]);
        let denominator = univariate::mul(
            &[frac(-3, 7)],
            &univariate::mul(
                &polynomial_power(&linear, 3),
                &polynomial_power(&quadratic, 2),
            ),
        );
        let expected = Apart {
            quotient: q(&[-2, 3]),
            pieces: vec![
                Piece {
                    factor: linear.clone(),
                    power: 1,
                    numerator: vec![frac(1, 3)],
                },
                Piece {
                    factor: linear,
                    power: 3,
                    numerator: vec![frac(5, 2)],
                },
                Piece {
                    factor: quadratic.clone(),
                    power: 1,
                    numerator: q(&[-3, 1]),
                },
                Piece {
                    factor: quadratic,
                    power: 2,
                    numerator: q(&[4, -2]),
                },
            ],
        };
        let mut numerator = recombine(&expected, &denominator);
        numerator.push(BigRational::zero());
        let mut padded_denominator = denominator;
        padded_denominator.push(BigRational::zero());
        assert_eq!(apart(&numerator, &padded_denominator), Some(expected));
    }

    #[test]
    fn independent_factor_blocks_preserve_factor_and_power_order_across_pools() {
        let factors: Vec<_> = (1..=12)
            .map(|k| (q(&[k, 1]), u32::try_from(k % 2 + 1).unwrap()))
            .collect();
        let denominator = factors
            .iter()
            .fold(q(&[1]), |product, (factor, multiplicity)| {
                univariate::mul(&product, &polynomial_power(factor, *multiplicity))
            });
        let expected = Apart {
            quotient: Vec::new(),
            pieces: factors
                .iter()
                .enumerate()
                .flat_map(|(index, (factor, multiplicity))| {
                    (1..=*multiplicity).map(move |power| {
                        Piece {
                            factor: factor.clone(),
                            power,
                            numerator: vec![frac(
                                i64::try_from(index).unwrap() + 1,
                                i64::from(power),
                            )],
                        }
                    })
                })
                .collect(),
        };
        let remainder = recombine(&expected, &denominator);
        for workers in [1, 2, 4] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap();
            assert_eq!(
                pool.install(|| decompose_proper(&remainder, &denominator, &factors)),
                Some(expected.pieces.clone())
            );
        }
    }

    #[test]
    fn invalid_factor_products_have_no_decomposition() {
        let denominator = q(&[1, 2, 1]);
        let remainder = q(&[1]);
        for factors in [
            vec![(q(&[1, 1]), 1)],
            vec![(q(&[1, 1]), 1), (q(&[1, 1]), 1)],
            vec![(q(&[0, 1]), 1), (q(&[1, 1]), 1)],
            vec![(q(&[1, 1]), 2), (q(&[2, 1]), 0)],
        ] {
            assert_eq!(decompose_proper(&remainder, &denominator, &factors), None);
        }
        // Duplicate / associate linear factors for simple denominator (x + 1)(x + 2).
        let denom2 = q(&[2, 3, 1]);
        assert_eq!(
            decompose_proper(&remainder, &denom2, &[(q(&[1, 1]), 1), (q(&[2, 2]), 1)]),
            None
        );
    }

    /// `1/(x^2 (x + 1)) = -1/x + 1/x^2 + 1/(x + 1)`, asserted as an exact
    /// term list.
    #[test]
    fn repeated_factor_terms_are_exact() {
        let parts = apart(&q(&[1]), &q(&[0, 0, 1, 1])).unwrap_or_else(|| panic!("none"));
        assert_eq!(parts.quotient, []);
        assert_eq!(
            parts.pieces,
            vec![
                Piece {
                    factor: q(&[0, 1]),
                    power: 1,
                    numerator: q(&[-1]),
                },
                Piece {
                    factor: q(&[0, 1]),
                    power: 2,
                    numerator: q(&[1]),
                },
                Piece {
                    factor: q(&[1, 1]),
                    power: 1,
                    numerator: q(&[1]),
                },
            ]
        );
    }

    /// `1/(x^3 + 1) = 1/3/(x + 1) + (2/3 - x/3)/(x^2 - x + 1)`: a linear term
    /// plus a proper term over a genuinely irreducible quadratic.
    #[test]
    fn irreducible_quadratic_terms_are_exact() {
        let parts = apart(&q(&[1]), &q(&[1, 0, 0, 1])).unwrap_or_else(|| panic!("none"));
        assert_eq!(parts.quotient, []);
        assert_eq!(
            parts.pieces,
            vec![
                Piece {
                    factor: q(&[1, 1]),
                    power: 1,
                    numerator: vec![frac(1, 3)],
                },
                Piece {
                    factor: q(&[1, -1, 1]),
                    power: 1,
                    numerator: vec![frac(2, 3), frac(-1, 3)],
                },
            ]
        );
        let quadratic = &parts.pieces[1].factor;
        assert!(
            univariate::rational_roots(quadratic).is_empty(),
            "x^2 - x + 1 must have no rational root"
        );
        assert_eq!(quadratic.len(), 3);
        assert_eq!(parts.pieces[1].numerator.len(), 2);
    }

    /// At every rational pole, the coefficient of `(x - r)^-k` equals the
    /// standard local limit: the Taylor coefficient of order `m - k` of
    /// `numer(r + t) / (denom / (x - r)^m)` evaluated at `t = 0`.
    #[test]
    fn pole_coefficients_are_the_local_limits() {
        let cases: [(&[i64], &[i64]); 4] = [
            (&[1], &[-1, 0, 1]),      // 1/(x^2 - 1)
            (&[1], &[0, 0, 1, 1]),    // 1/(x^2 (x + 1))
            (&[1, 2], &[0, 0, 1, 1]), // (2x + 1)/(x^2 (x + 1))
            (&[1], &[-1, 0, 0, 1]),   // 1/(x^3 - 1)
        ];
        for (n, d) in cases {
            let (numer, denom) = (q(n), q(d));
            let parts = apart(&numer, &denom).unwrap_or_else(|| panic!("none"));
            let roots = univariate::rational_roots(&denom);
            assert!(!roots.is_empty(), "{d:?} has rational poles");
            for root in roots {
                let (power, cofactor) = local(&denom, &root);
                assert!(power >= 1);
                let order = power as usize;
                let expanded = series(
                    &taylor(&numer, &root, order),
                    &taylor(&cofactor, &root, order),
                    order,
                );
                for k in 1..=power {
                    let expected = expanded[order - k as usize].clone();
                    let got = parts
                        .pieces
                        .iter()
                        .find(|p| {
                            p.power == k
                                && p.factor.len() == 2
                                && univariate::eval(&p.factor, &root).is_zero()
                        })
                        .and_then(|p| p.numerator.first().cloned());
                    assert_eq!(got, Some(expected), "pole {root} order {k} of {n:?}/{d:?}");
                }
            }
        }
    }

    /// `1/(x^2 - 2) = (sqrt 2/4)/(x - sqrt 2) - (sqrt 2/4)/(x + sqrt 2)`,
    /// exactly, with the recombination cross-multiplied in `Q(sqrt 8)`.
    #[test]
    fn real_quadratic_splits_over_the_extension() {
        let numer = q(&[1]);
        let denom = q(&[-2, 0, 1]);
        let parts = apart(&numer, &denom).unwrap_or_else(|| panic!("none"));
        assert_eq!(
            parts.pieces,
            vec![Piece {
                factor: denom.clone(),
                power: 1,
                numerator: q(&[1]),
            }]
        );
        let split = split_real_quadratics(&parts).unwrap_or_else(|| panic!("no split"));
        assert_eq!(split.len(), 2);
        assert_eq!(split[0].radicand, frac(8, 1));
        assert_eq!(split[0].real, BigRational::zero());
        assert_eq!(split[0].radical, frac(1, 8));
        assert_eq!(split[0].offset, BigRational::zero());
        assert_eq!(split[0].scale, frac(1, 2));
        assert_eq!(split[1].radical, frac(-1, 8));
        assert_eq!(split[1].scale, frac(-1, 2));
        quadratic_pair(&split, &numer, &denom);
    }

    /// The same for a quadratic with a non-zero linear term.
    #[test]
    fn shifted_real_quadratic_splits() {
        let numer = q(&[1]);
        let denom = q(&[-1, 1, 1]); // x^2 + x - 1
        let parts = apart(&numer, &denom).unwrap_or_else(|| panic!("none"));
        let split = split_real_quadratics(&parts).unwrap_or_else(|| panic!("no split"));
        assert_eq!(split.len(), 2);
        assert_eq!(split[0].radicand, frac(5, 1));
        assert_eq!(split[0].radical, frac(1, 5));
        assert_eq!(split[0].offset, frac(-1, 2));
        assert_eq!(split[0].scale, frac(1, 2));
        assert_eq!(split[1].radical, frac(-1, 5));
        assert_eq!(split[1].scale, frac(-1, 2));
        quadratic_pair(&split, &numer, &denom);
        assert_eq!(univariate::rational_roots(&denom), []);
    }

    /// A rational factor keeps its exact residue when a real quadratic sits
    /// beside it: `(x + 1)/((x + 2)(x^2 - 3))` has residue `-1` at `x = -2`.
    #[test]
    fn linear_and_quadratic_factors_mix() {
        let numer = q(&[1, 1]);
        let denom = univariate::mul(&q(&[2, 1]), &q(&[-3, 0, 1]));
        let parts = apart(&numer, &denom).unwrap_or_else(|| panic!("none"));
        let split = split_real_quadratics(&parts).unwrap_or_else(|| panic!("no split"));
        assert_eq!(split.len(), 3);
        assert_eq!(split[0].offset, frac(-2, 1));
        assert_eq!(split[0].scale, BigRational::zero());
        assert_eq!(split[0].radicand, BigRational::zero());
        assert_eq!(split[0].real, frac(-1, 1));
        assert_eq!(split[0].power, 1);
        assert_eq!(split[1].radicand, frac(12, 1));
        assert_eq!(split[1].scale, frac(1, 2));
        assert_eq!(split[2].scale, frac(-1, 2));
        quadratic_pair(&split[1..], &q(&[-1, 1]), &q(&[-3, 0, 1]));
    }

    /// Refusal is explicit: a cubic that is irreducible over `Q`, a repeated
    /// real quadratic, and a zero denominator are all reported as `None`
    /// rather than answered with a wrong structure.
    #[test]
    fn refusal_is_honest() {
        assert_eq!(apart(&q(&[1]), &[]), None);
        assert_eq!(apart(&q(&[1]), &q(&[0])), None);
        // x^3 - 2 has no rational root, so it stays one cubic factor.
        let cubic = apart(&q(&[1]), &q(&[-2, 0, 0, 1])).unwrap_or_else(|| panic!("none"));
        assert_eq!(cubic.pieces.len(), 1);
        assert_eq!(cubic.pieces[0].factor, q(&[-2, 0, 0, 1]));
        assert_eq!(split_real_quadratics(&cubic), None);
        // (x^2 - 2)^2 is a double pole over Q(sqrt 2), which is not attempted.
        let repeated = apart(&q(&[1]), &q(&[4, 0, -4, 0, 1])).unwrap_or_else(|| panic!("none"));
        assert_eq!(
            repeated.pieces,
            vec![Piece {
                factor: q(&[-2, 0, 1]),
                power: 2,
                numerator: q(&[1]),
            }]
        );
        assert_eq!(split_real_quadratics(&repeated), None);
        // A quadratic with complex roots is already irreducible over R.
        let complex = apart(&q(&[1]), &q(&[1, 0, 1])).unwrap_or_else(|| panic!("none"));
        assert_eq!(split_real_quadratics(&complex), None);
    }

    #[test]
    fn high_multiplicity_linear_and_quadratic_factors_recombine_exactly() {
        // Non-monic linear factor (2x - 3)^10 alongside (x + 1)^8 and (2x^2 - x + 3)^5.
        let lin1 = q(&[-3, 2]);
        let lin2 = q(&[1, 1]);
        let quad = q(&[3, -1, 2]);
        let denom = univariate::mul(
            &polynomial_power(&lin1, 10),
            &univariate::mul(&polynomial_power(&lin2, 8), &polynomial_power(&quad, 5)),
        );
        // Improper numerator exceeding denominator degree.
        let numer: QPoly = (0..=32)
            .map(|i| frac((i * 7 - 5) % 13 + 1, (i % 3) + 1))
            .collect();
        let parts = apart(&numer, &denom).expect("high multiplicity decomposition");
        assert_eq!(
            recombination_defect(&parts, &numer, &denom),
            [] as [num_rational::Ratio<num_bigint::BigInt>; 0]
        );
        for piece in &parts.pieces {
            assert_ne!(
                piece.numerator,
                [] as [num_rational::Ratio<num_bigint::BigInt>; 0]
            );
            assert!(piece.numerator.len() < piece.factor.len());
        }
        // Verify binary-splitting f-adic path (`multiplicity >= 8`) directly on (x^2 + 1)^9.
        let q_fac = q(&[1, 0, 1]);
        let d_hi = univariate::mul(&polynomial_power(&q_fac, 9), &q(&[2, 3]));
        let n_hi: QPoly = (0..=20).map(|i| frac(i - 4, 5)).collect();
        let parts_hi = apart(&n_hi, &d_hi).expect("multiplicity-9 quadratic decomposition");
        assert_eq!(
            recombination_defect(&parts_hi, &n_hi, &d_hi),
            [] as [num_rational::Ratio<num_bigint::BigInt>; 0]
        );
    }

    #[test]
    fn irreducible_higher_degree_denominators_and_large_crt_inverses() {
        // Product of irreducible cubic (x^3 - 2), quartic (x^4 + x + 1),
        // and cyclotomic sextic (x^6 + x^3 + 1).
        let cub = q(&[-2, 0, 0, 1]);
        let quart = q(&[1, 1, 0, 0, 1]);
        let sext = q(&[1, 0, 0, 1, 0, 0, 1]);
        let denom = univariate::mul(&cub, &univariate::mul(&quart, &sext));
        let numer: QPoly = (0..=15).map(|i| frac(i * 11 - 17, 7)).collect();
        let parts = apart(&numer, &denom).expect("higher-degree decomposition");
        assert_eq!(parts.pieces.len(), 3);
        assert_eq!(
            recombination_defect(&parts, &numer, &denom),
            [] as [num_rational::Ratio<num_bigint::BigInt>; 0]
        );

        // Large-coefficient irreducible quartics forcing multi-prime CRT inverse beyond 31 bits.
        let f1: ZPoly = [-100_003_i64, 17, -23, 5, 19]
            .into_iter()
            .map(BigInt::from)
            .collect();
        let f2: ZPoly = [99_991_i64, -31, 41, 13, 29]
            .into_iter()
            .map(BigInt::from)
            .collect();
        let (inv_num, inv_den) = modular_inverse_mod_z(&f2, &f1).expect("modular inverse");
        let check = z_sub_const(z_mul(&inv_num, &f2), &inv_den);
        assert!(univariate::exact_div(&check, &f1).is_some());
    }

    #[test]
    fn edge_cases_and_parallel_determinism() {
        // Zero numerator, constant denominator, exact divisibility.
        assert_eq!(
            apart(&[], &q(&[1, 0, 1])),
            Some(Apart {
                quotient: Vec::new(),
                pieces: Vec::new(),
            })
        );
        assert_eq!(
            apart(&q(&[0, 0]), &q(&[1, 1])),
            Some(Apart {
                quotient: Vec::new(),
                pieces: Vec::new(),
            })
        );
        assert_eq!(
            apart(&q(&[2, 4, 6]), &q(&[2])),
            Some(Apart {
                quotient: q(&[1, 2, 3]),
                pieces: Vec::new(),
            })
        );
        assert_eq!(
            apart(&q(&[-1, 0, 1]), &q(&[-1, 1])),
            Some(Apart {
                quotient: q(&[1, 1]),
                pieces: Vec::new(),
            })
        );

        // Large multi-factor decomposition deterministic under serial and parallel execution.
        let mut denom = q(&[1]);
        for k in 1_i64..=14 {
            denom = univariate::mul(&denom, &q(&[-k, 2]));
            if k <= 4 {
                denom = univariate::mul(&denom, &q(&[k, 0, 1]));
            }
        }
        let numer: QPoly = (0..=25).map(|i| frac(i + 1, 3)).collect();
        let seq = execution::with_scoped_pool(Some(1), || apart(&numer, &denom));
        let par = execution::with_scoped_pool(Some(4), || apart(&numer, &denom));
        assert!(seq.is_some());
        assert_eq!(seq, par);
        let parts = seq.unwrap();
        assert_eq!(
            recombination_defect(&parts, &numer, &denom),
            [] as [num_rational::Ratio<num_bigint::BigInt>; 0]
        );
    }

    #[test]
    fn mixed_symbolic_parameters_recombine_exactly() {
        let x = Poly::generator(0);
        let a = Poly::generator(1);
        let b = Poly::generator(2);
        let c = |n: i64| Poly::constant(Number::from(n));

        // Recombination compares the two sides as rational functions. A piece
        // stands for `numer / (factor^power * denom)`, so the sum is
        // accumulated over a common denominator and cross-multiplied against
        // the original ratio; no exact division is involved, which a piece
        // with a parameter denominator would not admit.
        let recombine_poly = |parts: &PolyApart| -> (Poly, Poly) {
            let mut numerator = parts.quotient_numer.clone();
            let mut denominator = parts.quotient_denom.clone();
            for piece in &parts.pieces {
                let f_pow = piece.factor.pow(piece.power, POLY_CAP).unwrap();
                let piece_den = f_pow.mul(&piece.denom, POLY_CAP).unwrap();
                // n/d + numer/piece_den
                let left = numerator.mul(&piece_den, POLY_CAP).unwrap();
                let right = piece.numer.mul(&denominator, POLY_CAP).unwrap();
                numerator = left.add(&right);
                denominator = denominator.mul(&piece_den, POLY_CAP).unwrap();
            }
            (numerator, denominator)
        };
        let same_ratio = |parts: &PolyApart, denom: &Poly, numer: &Poly| {
            let (n, d) = recombine_poly(parts);
            n.mul(denom, POLY_CAP).unwrap() == numer.mul(&d, POLY_CAP).unwrap()
        };

        // Case 1: (a*x + b) / (x^2 - 1)
        let n1 = a.mul(&x, POLY_CAP).unwrap().add(&b);
        let d1 = x.pow(2, POLY_CAP).unwrap().sub(&c(1));
        let p1 = apart_poly(&n1, &d1, 0).expect("symbolic numerator");
        assert_eq!(p1.pieces.len(), 2);
        assert!(same_ratio(&p1, &d1, &n1), "case 1 does not recombine");

        // Case 2: 1 / (x^2 - a^2)
        let n2 = c(1);
        let d2 = x
            .pow(2, POLY_CAP)
            .unwrap()
            .sub(&a.pow(2, POLY_CAP).unwrap());
        let p2 = apart_poly(&n2, &d2, 0).expect("symbolic difference of squares");
        assert_eq!(p2.pieces.len(), 2);
        assert!(same_ratio(&p2, &d2, &n2), "case 2 does not recombine");

        // Case 3: (x^3 + b) / ((x - a)^2 * (x + b))
        let n3 = x.pow(3, POLY_CAP).unwrap().add(&b);
        let d3 = x
            .sub(&a)
            .pow(2, POLY_CAP)
            .unwrap()
            .mul(&x.add(&b), POLY_CAP)
            .unwrap();
        let p3 = apart_poly(&n3, &d3, 0).expect("symbolic repeated pole");
        assert!(same_ratio(&p3, &d3, &n3), "case 3 does not recombine");

        // Case 4: (a*x + b) / (a * (x^2 + 1)) -> parameter denominator only.
        let n4 = a.mul(&x, POLY_CAP).unwrap().add(&b);
        let d4 = a
            .mul(&x.pow(2, POLY_CAP).unwrap().add(&c(1)), POLY_CAP)
            .unwrap();
        let p4 = apart_poly(&n4, &d4, 0).expect("parameter factor in denominator");
        assert!(same_ratio(&p4, &d4, &n4), "case 4 does not recombine");
    }
}
