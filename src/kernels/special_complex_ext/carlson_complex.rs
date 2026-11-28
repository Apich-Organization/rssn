//! Shared duplication and homogeneous scaling for the complex Carlson forms.

use std::f64::consts::FRAC_PI_4;

use num_complex::Complex64;

const ZERO: Complex64 = Complex64::new(0.0, 0.0);
const ONE: Complex64 = Complex64::new(1.0, 0.0);
pub(super) const NAN: Complex64 = Complex64::new(f64::NAN, f64::NAN);
const MAX_ITERATIONS: usize = 80;
// The seventh-order expansion is well within binary64 error at this radius.
const NEAR_CENTRE_SQUARED: f64 = 4.9e-5;

fn normalized<const N: usize>(arguments: [Complex64; N]) -> Option<([Complex64; N], f64)> {
    if arguments.iter().any(|value| !value.is_finite()) {
        return None;
    }
    let scale = arguments.iter().fold(0.0_f64, |scale, value| {
        scale.max(value.re.abs()).max(value.im.abs())
    });
    // Leave ordinary arguments and their rounding unscaled. This also avoids
    // dividing every input on the common path.
    if scale == 0.0 || (1.0e-60..=1.0e60).contains(&scale) {
        return Some((arguments, 1.0));
    }
    let scaled = arguments.map(|argument| argument / scale);
    if arguments.iter().zip(scaled).any(|(original, scaled)| {
        (original.re != 0.0 && scaled.re == 0.0) || (original.im != 0.0 && scaled.im == 0.0)
    }) {
        return None;
    }
    Some((scaled, scale))
}

fn stable_root(z: Complex64) -> Complex64 {
    let root = crate::kernels::complex_scalar::Complex::new(z.re, z.im).sqrt();
    Complex64::new(root.re, root.im)
}

fn duplicate([x, y, z]: [Complex64; 3]) -> ([Complex64; 3], [Complex64; 3]) {
    let roots @ [sx, sy, sz] = [x.sqrt(), y.sqrt(), z.sqrt()];
    // DLMF 19.36.4. Form x + lambda as (sx + sy)(sx + sz). The
    // factorization avoids subtracting large nearly equal values near a cut.
    let xy = sx + sy;
    let xz = sx + sz;
    let yz = sy + sz;
    ([0.25 * xy * xz, 0.25 * xy * yz, 0.25 * xz * yz], roots)
}

fn deviations<const N: usize>(
    values: [Complex64; N],
    centre: Complex64,
) -> Option<[Complex64; N]> {
    // A single Taylor expansion cannot represent arguments on opposite
    // sides of the cut, even when the distance between them is small.
    if centre.re < 0.0
        && values
            .iter()
            .any(|value| value.im.is_sign_negative() != centre.im.is_sign_negative())
    {
        return None;
    }
    let deviations = values.map(|value| (centre - value) / centre);
    deviations
        .iter()
        .all(|deviation| deviation.norm_sqr() < NEAR_CENTRE_SQUARED)
        .then_some(deviations)
}

fn rf(mut values: [Complex64; 3]) -> Option<Complex64> {
    for _ in 0..MAX_ITERATIONS {
        let [x, y, z] = values;
        let centre = (x + y + z) / 3.0;
        if let Some([dx, dy, dz]) = deviations(values, centre) {
            // DLMF 19.36.1
            let xy = dx * dy;
            let e2 = xy - dz * dz;
            let e3 = xy * dz;
            let e2_squared = e2 * e2;
            let polynomial = ONE - e2 / 10.0 + e3 / 14.0 + e2_squared / 24.0
                - 3.0 * e2 * e3 / 44.0
                - 5.0 * e2_squared * e2 / 208.0
                + 3.0 * e3 * e3 / 104.0
                + e2_squared * e3 / 16.0;
            return Some(polynomial / centre.sqrt());
        }
        (values, _) = duplicate(values);
    }
    None
}

fn third_kind_series(
    [x, y, z, p]: [Complex64; 4],
    centre: Complex64,
) -> Complex64 {
    // The two extra copies of p turn the four distinct variables into the
    // five variables in the DLMF 19.36.2 series. For R_D, p is z.
    let sum_two = x * (y + z) + y * z;
    let product = x * y * z;
    let p_squared = p * p;
    let e2 = sum_two - 3.0 * p_squared;
    let e3 = product + 2.0 * p * (sum_two - p_squared);
    let e4 = 2.0 * p * product + p_squared * sum_two;
    let e5 = p_squared * product;
    let e2_squared = e2 * e2;
    let polynomial = ONE - 3.0 * e2 / 14.0 + e3 / 6.0 + 9.0 * e2_squared / 88.0
        - 3.0 * e4 / 22.0
        - 9.0 * e2 * e3 / 52.0
        + 3.0 * e5 / 26.0
        - e2_squared * e2 / 16.0
        + 3.0 * e3 * e3 / 40.0
        + 3.0 * e2 * e4 / 20.0
        + 45.0 * e2_squared * e3 / 272.0
        - 9.0 * (e3 * e4 + e2 * e5) / 68.0;
    polynomial / (centre * centre.sqrt())
}

fn rd(mut values: [Complex64; 3]) -> Option<Complex64> {
    let mut sum = ZERO;
    let mut factor = 1.0;
    for _ in 0..MAX_ITERATIONS {
        let [x, y, z] = values;
        let centre = (x + y + 3.0 * z) / 5.0;
        if let Some([dx, dy, dz]) = deviations(values, centre) {
            return Some(3.0 * sum + factor * third_kind_series([dx, dy, dz, dz], centre));
        }
        let (next, [_, _, sz]) = duplicate(values);
        sum += (0.25 * factor * ONE) / (sz * next[2]);
        factor *= 0.25;
        values = next;
    }
    None
}

fn rj([x, y, z, mut p]: [Complex64; 4]) -> Option<Complex64> {
    if super::carlson_same_branch(p, x) {
        return rd([y, z, x]);
    }
    if super::carlson_same_branch(p, y) {
        return rd([x, z, y]);
    }
    if super::carlson_same_branch(p, z) {
        return rd([x, y, z]);
    }
    let mut values = [x, y, z];
    let delta = (p - x) * (p - y) * (p - z);
    let mut sum = ZERO;
    let mut factor = 1.0;
    for _ in 0..MAX_ITERATIONS {
        let [x, y, z] = values;
        let centre = (x + y + z + 2.0 * p) / 5.0;
        if let Some(delta) = deviations([x, y, z, p], centre) {
            return Some(6.0 * sum + factor * third_kind_series(delta, centre));
        }
        let (next, [sx, sy, sz]) = duplicate(values);
        let sp = p.sqrt();
        let d = (sp + sx) * (sp + sy) * (sp + sz);
        let e = (factor * factor * factor * delta) / (d * d);
        sum += factor * rf([ONE, ONE + e, ONE + e])? / d;
        factor *= 0.25;
        values = next;
        p = next[2] + 0.25 * (p - z);
    }
    None
}
fn rg(
    x: Complex64,
    y: Complex64,
    z: Complex64,
) -> Option<Complex64> {
    let (x, mut y, mut z) = if y == ZERO {
        (y, x, z)
    } else if z == ZERO {
        (z, x, y)
    } else {
        (x, y, z)
    };
    if x == ZERO {
        // DLMF 19.20.4. Two zeros are allowed for R_G but not R_F or R_D.
        if y == ZERO {
            return Some(0.5 * stable_root(z));
        }
        if z == ZERO {
            return Some(0.5 * stable_root(y));
        }
        // The equal-argument limit also retains the selected signed side of
        // the negative real axis through the principal square root.
        if y == z && y.im.is_sign_negative() == z.im.is_sign_negative() {
            return Some(FRAC_PI_4 * stable_root(z));
        }
        // For one zero, keep the larger nonzero magnitude in the R_D
        // denominator and preserve each argument's original branch side.
        if y.norm() > z.norm() {
            std::mem::swap(&mut y, &mut z);
        }
    }
    // DLMF 19.21.10 uses the product of the principal square roots.
    let roots = if x == ZERO {
        ZERO
    } else {
        stable_root(x) * stable_root(y) / stable_root(z)
    };
    Some(0.5 * (z * rf([x, y, z])? - (x - z) * (y - z) * rd([x, y, z])? / 3.0 + roots))
}

pub(super) fn first_kind(arguments: [Complex64; 3]) -> Option<Complex64> {
    let (arguments, scale) = normalized(arguments)?;
    if arguments.iter().filter(|value| **value == ZERO).count() > 1 {
        return None;
    }
    let result = rf(arguments)? / scale.sqrt();
    result.is_finite().then_some(result)
}

pub(super) fn third_kind_degenerate(arguments: [Complex64; 3]) -> Option<Complex64> {
    let (arguments, scale) = normalized(arguments)?;
    let [x, y, z] = arguments;
    if z == ZERO || x == ZERO && y == ZERO {
        return None;
    }
    let result = (rd(arguments)? / scale.sqrt()) / scale;
    result.is_finite().then_some(result)
}

pub(super) fn second_kind(arguments: [Complex64; 3]) -> Option<Complex64> {
    let (arguments, scale) = normalized(arguments)?;
    let [x, y, z] = arguments;
    let result = rg(x, y, z)? * scale.sqrt();
    result.is_finite().then_some(result)
}

pub(super) fn third_kind(arguments: [Complex64; 4]) -> Option<Complex64> {
    let (arguments, scale) = normalized(arguments)?;
    let [x, y, z, p] = arguments;
    if p == ZERO || [x, y, z].iter().filter(|value| **value == ZERO).count() > 1 {
        return None;
    }
    let result = (rj(arguments)? / scale.sqrt()) / scale;
    result.is_finite().then_some(result)
}
