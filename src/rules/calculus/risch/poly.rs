//! Generic dense polynomials over a coefficient field.
#![allow(
    // The tower arithmetic indexes coefficient vectors and shifts degrees by
    // hand; every index comes from a length that the surrounding loop just
    // checked, and a tower of two or three levels keeps the vectors short.
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::pedantic,
    clippy::nursery
)]

use super::elem::Coeff;

pub fn trim<F: Coeff>(mut p: Vec<F>) -> Vec<F> {
    while p.last().is_some_and(Coeff::is_zero) {
        p.pop();
    }
    p
}

pub fn trim_ref<F: Coeff>(p: &[F]) -> &[F] {
    let mut n = p.len();
    while n > 0 && p[n - 1].is_zero() {
        n -= 1;
    }
    &p[..n]
}

pub fn is_zero<F: Coeff>(p: &[F]) -> bool {
    trim_ref(p).is_empty()
}

pub fn degree<F: Coeff>(p: &[F]) -> Option<usize> {
    trim_ref(p).len().checked_sub(1)
}

pub fn add<F: Coeff>(
    a: &[F],
    b: &[F],
) -> Vec<F> {
    let mut out: Vec<F> = a.to_vec();
    out.resize(a.len().max(b.len()), F::zero());
    for (slot, x) in out.iter_mut().zip(b) {
        *slot = slot.add(x);
    }
    trim(out)
}

pub fn neg<F: Coeff>(a: &[F]) -> Vec<F> {
    a.iter().map(Coeff::neg).collect()
}

pub fn sub<F: Coeff>(
    a: &[F],
    b: &[F],
) -> Vec<F> {
    add(a, &neg(b))
}

pub fn scale<F: Coeff>(
    a: &[F],
    c: &F,
) -> Vec<F> {
    if c.is_zero() {
        return Vec::new();
    }
    if c.is_one() {
        return trim_ref(a).to_vec();
    }
    trim(a.iter().map(|x| x.mul(c)).collect())
}

pub fn mul<F: Coeff>(
    a: &[F],
    b: &[F],
) -> Vec<F> {
    let a = trim_ref(a);
    let b = trim_ref(b);
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    if a.len() == 1 {
        return scale(b, &a[0]);
    }
    if b.len() == 1 {
        return scale(a, &b[0]);
    }
    let mut out = vec![F::zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        if x.is_zero() {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            out[i + j] = out[i + j].add(&x.mul(y));
        }
    }
    trim(out)
}

pub fn divrem<F: Coeff>(
    a: &[F],
    b: &[F],
) -> Option<(Vec<F>, Vec<F>)> {
    let b = trim_ref(b);
    if b.is_empty() {
        return None;
    }
    let mut r = trim_ref(a).to_vec();
    if r.len() < b.len() {
        return Some((Vec::new(), r));
    }
    let lead = b.last()?.clone();
    let inv = lead.inv()?;
    let mut quo = vec![F::zero(); r.len() - b.len() + 1];
    while r.len() >= b.len() {
        let d = r.len() - b.len();
        let c = r.last()?.mul(&inv);
        quo[d] = c.clone();
        for (i, bc) in b.iter().enumerate() {
            r[d + i] = r[d + i].sub(&c.mul(bc));
        }
        while r.last().is_some_and(Coeff::is_zero) {
            r.pop();
        }
        if d == 0 {
            break;
        }
    }
    Some((trim(quo), trim(r)))
}

pub fn exact_div<F: Coeff>(
    a: &[F],
    b: &[F],
) -> Option<Vec<F>> {
    let (quo, rem) = divrem(a, b)?;
    if rem.is_empty() {
        Some(quo)
    } else {
        None
    }
}

pub fn monic<F: Coeff>(p: &[F]) -> Vec<F> {
    let p = trim_ref(p);
    match p.last() {
        | Some(lead) => {
            if lead.is_one() {
                return p.to_vec();
            }
            match lead.inv() {
                | Some(inv) => p.iter().map(|c| c.mul(&inv)).collect(),
                | None => p.to_vec(),
            }
        },
        | None => Vec::new(),
    }
}

/// Monic gcd over the coefficient field.
pub fn gcd<F: Coeff>(
    a: &[F],
    b: &[F],
) -> Vec<F> {
    let (mut x, mut y) = (trim_ref(a).to_vec(), trim_ref(b).to_vec());
    if x.is_empty() {
        return monic(&y);
    }
    if y.is_empty() {
        return monic(&x);
    }
    if x.len() == 1 || y.len() == 1 {
        return vec![F::one()];
    }
    let mut guard = 0;
    while !y.is_empty() {
        guard += 1;
        if guard > 2048 {
            return vec![F::one()];
        }
        let Some((_, r)) = divrem(&x, &y) else {
            return vec![F::one()];
        };
        x = y;
        y = r;
    }
    monic(&x)
}

/// Formal derivative with respect to the polynomial variable.
pub fn formal_derivative<F: Coeff>(p: &[F]) -> Vec<F> {
    p.iter()
        .enumerate()
        .skip(1)
        .map(|(i, c)| {
            let k = F::from_q(&num_rational::BigRational::from_integer(
                num_bigint::BigInt::from(i),
            ));
            c.mul(&k)
        })
        .collect()
}

/// `d/dx` of a polynomial in `theta` with `theta' = slope*theta + constant`.
pub fn derivative_x<F: Coeff>(
    p: &[F],
    slope: &F,
    constant: &F,
    tower: &super::tower::Tower,
) -> Vec<F> {
    let mut out = Vec::with_capacity(p.len());
    for (i, c) in p.iter().enumerate() {
        let mut term = c.deriv(tower);
        if i > 0 {
            let k = F::from_q(&num_rational::BigRational::from_integer(
                num_bigint::BigInt::from(i),
            ));
            term = term.add(&k.mul(&slope.mul(c)));
        }
        if let Some(next) = p.get(i + 1) {
            let k = F::from_q(&num_rational::BigRational::from_integer(
                num_bigint::BigInt::from(i + 1),
            ));
            term = term.add(&k.mul(&constant.mul(next)));
        }
        out.push(term);
    }
    trim(out)
}

/// `(D_sqf, D_2)` with `D = D_sqf * D_2`, `D_sqf` squarefree and `D_2 = gcd(D, D')`.
pub fn squarefree_parts<F: Coeff>(d: &[F]) -> (Vec<F>, Vec<F>) {
    let g = gcd(d, &formal_derivative(d));
    let sqf = exact_div(d, &g).unwrap_or_else(|| d.to_vec());
    (monic(&sqf), g)
}

/// Sylvester resultant over the coefficient field.
pub fn resultant<F: Coeff>(
    p: &[F],
    r: &[F],
) -> Option<F> {
    let (mut left, mut right) = (trim_ref(p).to_vec(), trim_ref(r).to_vec());
    if left.is_empty() || right.is_empty() {
        return None;
    }
    let mut value = F::one();
    let mut guard = 0;
    while right.len() > 1 {
        guard += 1;
        if guard > 2048 {
            return None;
        }
        let m = left.len() - 1;
        let n = right.len() - 1;
        if m % 2 == 1 && n % 2 == 1 {
            value = value.neg();
        }
        if m < n {
            std::mem::swap(&mut left, &mut right);
            continue;
        }
        let (_, remainder) = divrem(&left, &right)?;
        if remainder.is_empty() {
            return Some(F::zero());
        }
        let rr = remainder.len() - 1;
        let lead = right.last()?.clone();
        let e = i64::try_from(m - rr).ok()?;
        value = value.mul(&lead.pow_i(e)?);
        (left, right) = (right, remainder);
    }
    let e = i64::try_from(left.len() - 1).ok()?;
    Some(value.mul(&right.first()?.pow_i(e)?))
}
