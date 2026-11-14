//! Exact rational dense decompositions and perfect-square arithmetic.
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::Zero;

/// Failure of an exact dense factorization.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ExactError {
    /// Input dimensions or output storage overflow do not match the data.
    Shape,
    /// A nonzero Gram-Schmidt norm has no rational square root.
    IrrationalNorm,
}
impl std::fmt::Display for ExactError {
    fn fmt(
        &self,
        formatter: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        formatter.write_str(match self {
            | Self::Shape => "incompatible exact matrix dimensions",
            | Self::IrrationalNorm => "an exact QR of this matrix needs an irrational square root",
        })
    }
}
impl std::error::Error for ExactError {}

/// Exact partial-pivot LU with row-major factors satisfying `P A = L U`.
/// Singular matrices retain valid factors with zero pivots.
///
/// # Errors
/// Returns [`ExactError::Shape`] for inconsistent dimensions or size overflow.
#[allow(clippy::type_complexity)] // Three named matrix factors preserve the established API.
pub fn lu_exact(
    n: usize,
    a: &[BigRational],
) -> Result<(Vec<BigRational>, Vec<BigRational>, Vec<BigRational>), ExactError> {
    let size = n.checked_mul(n).ok_or(ExactError::Shape)?;
    if size != a.len() {
        return Err(ExactError::Shape);
    }
    let mut u: Vec<BigRational> = a.to_vec();
    let mut l = vec![BigRational::zero(); size];
    let mut p = vec![BigRational::zero(); size];
    for index in 0..n {
        l[index * n + index] = BigRational::one();
        p[index * n + index] = BigRational::one();
    }
    for column in 0..n {
        let mut pivot = column;
        while pivot < n && u[pivot * n + column].is_zero() {
            pivot += 1;
        }
        if pivot == n {
            // Singular: no pivot in this column. Continue with the identity,
            // which is what an exact `L U` of a singular matrix is.
            continue;
        }
        if pivot != column {
            for inner in 0..n {
                u.swap(pivot * n + inner, column * n + inner);
                p.swap(pivot * n + inner, column * n + inner);
            }
            for inner in 0..column {
                l.swap(pivot * n + inner, column * n + inner);
            }
        }
        let lead = u[column * n + column].clone();
        for row in column + 1..n {
            let factor = u[row * n + column].clone() / &lead;
            l[row * n + column] = factor.clone();
            for inner in column..n {
                let term = factor.clone() * u[column * n + inner].clone();
                u[row * n + inner] -= term;
            }
        }
    }
    Ok((p, l, u))
}

/// An exact QR decomposition by modified Gram-Schmidt.
///
/// # Errors
///
/// Invalid matrix dimensions, size overflow, or an irrational column norm.
pub fn qr_exact(
    rows: usize,
    cols: usize,
    a: &[BigRational],
) -> Result<(Vec<BigRational>, Vec<BigRational>), ExactError> {
    let size = rows.checked_mul(cols).ok_or(ExactError::Shape)?;
    let square = rows.checked_mul(rows).ok_or(ExactError::Shape)?;
    if size != a.len() {
        return Err(ExactError::Shape);
    }
    let mut q = vec![BigRational::zero(); square];
    let mut r = vec![BigRational::zero(); size];
    let mut rank = 0;
    for column in 0..cols {
        let mut v: Vec<BigRational> = (0..rows)
            .map(|row| a[row * cols + column].clone())
            .collect();
        for inner in 0..rank {
            let mut dot = BigRational::zero();
            for row in 0..rows {
                dot += q[row * rows + inner].clone() * v[row].clone();
            }
            r[inner * cols + column] = dot.clone();
            for row in 0..rows {
                let term = dot.clone() * q[row * rows + inner].clone();
                v[row] -= term;
            }
        }
        let norm = v.iter().fold(BigRational::zero(), |acc, value| {
            acc + value.clone() * value.clone()
        });
        if norm.is_zero() {
            continue;
        }
        let Some(length) = exact_sqrt(&norm) else {
            return Err(ExactError::IrrationalNorm);
        };
        r[rank * cols + column] = length.clone();
        for row in 0..rows {
            q[row * rows + rank] = v[row].clone() / &length;
        }
        rank += 1;
    }
    Ok((q, r))
}

/// The exact square root of a non-negative rational, when it has one.
#[must_use]
pub fn exact_sqrt(value: &BigRational) -> Option<BigRational> {
    if value.is_negative() {
        return None;
    }
    let numerator = integer_sqrt(value.numer())?;
    let denominator = integer_sqrt(value.denom())?;
    Some(BigRational::new(numerator, denominator))
}

/// The integer square root, when the argument is a perfect square.
fn integer_sqrt(value: &num_bigint::BigInt) -> Option<num_bigint::BigInt> {
    if value.is_zero() {
        return Some(BigInt::zero());
    }
    if value.is_negative() {
        return None;
    }
    // Newton's method on integers: exact, and it never needs a float that
    // cannot hold the value.
    let mut guess: BigInt = {
        let bits = value.bits();
        BigInt::one() << ((bits as usize).div_ceil(2))
    };
    loop {
        let next = (&guess + value / &guess) >> 1;
        if next >= guess {
            break;
        }
        guess = next;
    }
    if &guess * &guess == *value {
        Some(guess)
    } else {
        None
    }
}
