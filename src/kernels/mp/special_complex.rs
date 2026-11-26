//! Arbitrary-precision complex special functions on [`MpComplex`].

use num_bigint::BigInt;

use super::ball::Ball;
use super::complex::MpComplex;
use super::consts;
use super::funcs;
use super::real::Real;
use super::real::Round;

/// Exact rational Bernoulli numbers `B_{2k}` for `k = 1..=8`: `(num, den)`.
const BERNOULLI_2K: [(i64, i64); 8] = [
    (1, 6),
    (-1, 30),
    (1, 42),
    (-1, 30),
    (5, 66),
    (-691, 2730),
    (7, 6),
    (-3617, 510),
];

fn real_rat(
    num: i64,
    den: i64,
    prec: u32,
) -> Result<Ball, String> {
    let n = Real::from_i64(num, prec)?;
    let d = Real::from_i64(den, prec)?;
    Ball::point(n).div(&Ball::point(d), prec)
}

fn mp_from_i64(
    n: i64,
    prec: u32,
) -> Result<MpComplex, String> {
    Ok(MpComplex::from_real(Ball::point(Real::from_i64(n, prec)?)))
}

fn mp_rat(
    num: i64,
    den: i64,
    prec: u32,
) -> Result<MpComplex, String> {
    Ok(MpComplex::from_real(real_rat(num, den, prec)?))
}

fn pi_ball(prec: u32) -> Result<Ball, String> {
    let inner = prec.saturating_add(16);
    let p = consts::pi(inner)?;
    let mid = p.with_precision(prec, Round::NearestEven)?;
    let exp = mid.floor_log2().map_or(-i64::from(prec), |t| {
        t.saturating_sub(i64::from(prec)).saturating_add(1)
    });
    let rad = Real::rounded(BigInt::from(1), exp, 1, Round::AwayFromZero)?;
    Ball::new(mid, rad)
}

impl MpComplex {
    /// Complex error function `erf(z)` at `prec` bits.
    ///
    /// # Errors
    /// Fails if series does not converge or exponents overflow.
    pub fn erf(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        let inner = prec.saturating_add(48);
        let z2 = self.mul(self, inner)?;
        let mut term = self.clone();
        let mut sum = term.clone();
        let max_terms = 400_usize;
        for k in 1..=max_terms {
            let kf = i64::try_from(k).unwrap_or(i64::MAX);
            let denom = mp_from_i64(kf, inner)?;
            term = term.mul(&z2, inner)?.neg().div(&denom, inner)?;
            let factor = mp_from_i64(2 * kf + 1, inner)?;
            let add = term.div(&factor, inner)?;
            sum = sum.add(&add, inner)?;
            let mag_re = add.re.mid.floor_log2().unwrap_or(i64::MIN);
            let mag_im = add.im.mid.floor_log2().unwrap_or(i64::MIN);
            if mag_re < -i64::from(prec + 16) && mag_im < -i64::from(prec + 16) {
                let pi_b = pi_ball(inner)?;
                let sqrt_pi = pi_b.map_monotone(&|v, p, r| v.sqrt(p, r), true, inner)?;
                let two_over_sqrt_pi =
                    Ball::point(Real::from_i64(2, inner)?).div(&sqrt_pi, inner)?;
                return sum.mul(&Self::from_real(two_over_sqrt_pi), prec);
            }
        }
        Err("mp_complex erf series did not converge".to_string())
    }

    /// Complex complementary error function `erfc(z) = 1 - erf(z)` at `prec` bits.
    ///
    /// # Errors
    /// Fails if `erf` fails.
    pub fn erfc(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        let inner = prec.saturating_add(24);
        let e = self.erf(inner)?;
        Self::one(inner).sub(&e, prec)
    }

    /// Principal complex log-gamma `ln(Gamma(z))` at `prec` bits (away from `(-inf, 0]`).
    ///
    /// # Errors
    /// Fails if `z` meets `(-inf, 0]` or exponents overflow.
    pub fn ln_gamma(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        let inner = prec.saturating_add(48);
        let shift_target = (prec / 8).clamp(12, 48);
        let mut z = self.clone();
        let mut shift_acc = Self::zero(inner);
        for _ in 0..128 {
            let re_hi = z.re.mid.add(&z.re.rad, inner, Round::TowardPositive)?;
            let target_r = Real::from_i64(i64::from(shift_target), inner)?;
            if re_hi.cmp_value(&target_r) >= 0 {
                break;
            }
            let lz = z.ln(inner)?;
            shift_acc = shift_acc.sub(&lz, inner)?;
            z = z.add(&Self::one(inner), inner)?;
        }
        let half = mp_rat(1, 2, inner)?;
        let z_minus_half = z.sub(&half, inner)?;
        let ln_z = z.ln(inner)?;
        let two_pi = Self::from_real(pi_ball(inner)?).mul(&mp_from_i64(2, inner)?, inner)?;
        let half_ln_2pi = two_pi.ln(inner)?.mul(&half, inner)?;
        let mut res = z_minus_half
            .mul(&ln_z, inner)?
            .sub(&z, inner)?
            .add(&half_ln_2pi, inner)?;

        let inv_z = Self::one(inner).div(&z, inner)?;
        let inv_z2 = inv_z.mul(&inv_z, inner)?;
        let mut pow_term = inv_z;
        for (idx, &(b_num, b_den)) in BERNOULLI_2K.iter().enumerate() {
            let k = i64::try_from(idx + 1).unwrap_or(1);
            let coeff = mp_rat(b_num, b_den * (2 * k) * (2 * k - 1), inner)?;
            res = res.add(&coeff.mul(&pow_term, inner)?, inner)?;
            pow_term = pow_term.mul(&inv_z2, inner)?;
        }
        res.add(&shift_acc, prec)
    }

    /// Complex gamma function `Gamma(z) = exp(ln_gamma(z))` at `prec` bits.
    ///
    /// # Errors
    /// Fails if `z` meets `(-inf, 0]` or exponents overflow.
    pub fn gamma(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        let inner = prec.saturating_add(32);
        self.ln_gamma(inner)?.exp(prec)
    }

    /// Complex digamma function `psi(z) = Gamma'(z) / Gamma(z)` at `prec` bits.
    ///
    /// # Errors
    /// Fails if `z` meets `(-inf, 0]` or exponents overflow.
    pub fn digamma(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        let inner = prec.saturating_add(48);
        let shift_target = (prec / 8).clamp(12, 48);
        let mut z = self.clone();
        let mut shift_acc = Self::zero(inner);
        for _ in 0..128 {
            let re_hi = z.re.mid.add(&z.re.rad, inner, Round::TowardPositive)?;
            let target_r = Real::from_i64(i64::from(shift_target), inner)?;
            if re_hi.cmp_value(&target_r) >= 0 {
                break;
            }
            shift_acc = shift_acc.sub(&Self::one(inner).div(&z, inner)?, inner)?;
            z = z.add(&Self::one(inner), inner)?;
        }
        let inv_z = Self::one(inner).div(&z, inner)?;
        let inv_z2 = inv_z.mul(&inv_z, inner)?;
        let mut res = z
            .ln(inner)?
            .sub(&inv_z.mul(&mp_rat(1, 2, inner)?, inner)?, inner)?;
        let mut pow_term = inv_z2.clone();
        for (idx, &(b_num, b_den)) in BERNOULLI_2K.iter().enumerate() {
            let k = i64::try_from(idx + 1).unwrap_or(1);
            let coeff = mp_rat(b_num, b_den * (2 * k), inner)?;
            res = res.sub(&coeff.mul(&pow_term, inner)?, inner)?;
            pow_term = pow_term.mul(&inv_z2, inner)?;
        }
        res.add(&shift_acc, prec)
    }

    /// Complex Bessel function `J_nu(z)` at `prec` bits for `nu >= 0`.
    ///
    /// # Errors
    /// Fails if `z` meets `(-inf, 0]` for non-integer `nu` or series fails to converge.
    pub fn bessel_j(
        &self,
        nu: &Self,
        prec: u32,
    ) -> Result<Self, String> {
        let inner = prec.saturating_add(48);
        let half = mp_rat(1, 2, inner)?;
        let hz = self.mul(&half, inner)?;
        let hz2 = hz.mul(&hz, inner)?;
        let nu_plus_1 = nu.add(&Self::one(inner), inner)?;
        let mut term = Self::one(inner).div(&nu_plus_1.gamma(inner)?, inner)?;
        let mut sum = term.clone();
        for k in 1..=300 {
            let kf = i64::try_from(k).unwrap_or(i64::MAX);
            let denom =
                mp_from_i64(kf, inner)?.mul(&nu.add(&mp_from_i64(kf, inner)?, inner)?, inner)?;
            term = term.mul(&hz2, inner)?.neg().div(&denom, inner)?;
            sum = sum.add(&term, inner)?;
            let mag_re = term.re.mid.floor_log2().unwrap_or(i64::MIN);
            let mag_im = term.im.mid.floor_log2().unwrap_or(i64::MIN);
            if mag_re < -i64::from(prec + 16) && mag_im < -i64::from(prec + 16) {
                let pref = hz.pow(nu, inner)?;
                return pref.mul(&sum, prec);
            }
        }
        Err("mp_complex bessel_j series did not converge".to_string())
    }

    /// Complex Riemann zeta function `zeta(s)` at `prec` bits for `Re(s) > 1` or `s != 1`.
    ///
    /// # Errors
    /// Fails at `s = 1` or on exponent overflow.
    pub fn riemann_zeta(
        &self,
        prec: u32,
    ) -> Result<Self, String> {
        let inner = prec.saturating_add(48);
        let one = Self::one(inner);
        let n_steps = (prec / 4).clamp(16, 128);
        let mut sum = Self::zero(inner);
        let neg_s = self.neg();
        for k in 1..n_steps {
            let kf = mp_from_i64(i64::from(k), inner)?;
            sum = sum.add(&kf.pow(&neg_s, inner)?, inner)?;
        }
        let n_mp = mp_from_i64(i64::from(n_steps), inner)?;
        let n_ms = n_mp.pow(&neg_s, inner)?;
        let int_term = n_mp
            .pow(&one.sub(self, inner)?, inner)?
            .div(&self.sub(&one, inner)?, inner)?;
        let half_term = n_ms.mul(&mp_rat(1, 2, inner)?, inner)?;

        let inv_n = one.div(&n_mp, inner)?;
        let inv_n2 = inv_n.mul(&inv_n, inner)?;
        let mut poch = self.clone();
        let mut pow_term = inv_n;
        let mut em = Self::zero(inner);
        // B_{2k}/(2k)! for k=1..=6
        const B_FACT: [(i64, i64); 6] = [
            (1, 12),
            (-1, 720),
            (1, 30_240),
            (-1, 1_209_600),
            (1, 47_900_160),
            (-691, 1_307_674_368_000),
        ];
        for (idx, &(num, den)) in B_FACT.iter().enumerate() {
            if idx > 0 {
                let two_k = i64::try_from(2 * idx).unwrap_or(2);
                poch = poch
                    .mul(&self.add(&mp_from_i64(two_k - 1, inner)?, inner)?, inner)?
                    .mul(&self.add(&mp_from_i64(two_k, inner)?, inner)?, inner)?;
                pow_term = pow_term.mul(&inv_n2, inner)?;
            }
            let c = mp_rat(num, den, inner)?;
            em = em.add(
                &c.mul(&poch, inner)?
                    .mul(&n_ms, inner)?
                    .mul(&pow_term, inner)?,
                inner,
            )?;
        }
        let _ = funcs::exp;
        sum.add(&int_term, inner)?
            .add(&half_term, inner)?
            .add(&em, prec)
    }

    /// Complex polylogarithm `Li_s(z) = sum_{k>=1} z^k / k^s` at `prec` bits for `|z| < 0.85`.
    ///
    /// # Errors
    /// Fails if `|z| >= 0.85` or series fails to converge.
    pub fn polylog(
        &self,
        s: &Self,
        prec: u32,
    ) -> Result<Self, String> {
        let inner = prec.saturating_add(48);
        let neg_s = s.neg();
        let mut z_pow = self.clone();
        let mut sum = Self::zero(inner);
        for k in 1..=512 {
            let kf = mp_from_i64(i64::from(k), inner)?;
            let term = z_pow.mul(&kf.pow(&neg_s, inner)?, inner)?;
            sum = sum.add(&term, inner)?;
            let mag_re = term.re.mid.floor_log2().unwrap_or(i64::MIN);
            let mag_im = term.im.mid.floor_log2().unwrap_or(i64::MIN);
            if k > 4 && mag_re < -i64::from(prec + 16) && mag_im < -i64::from(prec + 16) {
                return sum.mul(&Self::one(prec), prec);
            }
            z_pow = z_pow.mul(self, inner)?;
        }
        Err("mp_complex polylog series did not converge (requires |z| < 0.85)".to_string())
    }
}
