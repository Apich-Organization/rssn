// Checked binary64 polynomial coefficients are in descending power order.
use crate::kernels::polynomial::Polynomial as RealPolynomial;

const MAX_BARYCENTRIC_NODES: usize = 2048;
const MAX_POLYNOMIAL_COEFFICIENTS: usize = 4097;
const MAX_GCD_COEFFICIENTS: usize = 513;

fn normalized_polynomial(mut coefficients: Vec<f64>) -> RealPolynomial {
    let scale = coefficients
        .iter()
        .copied()
        .map(f64::abs)
        .fold(0.0, f64::max);
    if scale > 0.0 {
        for coefficient in &mut coefficients {
            *coefficient /= scale;
        }
    }
    RealPolynomial::new(coefficients)
}

binding!(
    polynomial_barycentric_eval,
    "num_real_polynomial_barycentric_eval",
    "numerical_special",
    "Interpolates finite samples with up to 2048 distinct nodes. The nodes and query share abscissa units, and the result uses ordinate units.",
    "crate::kernels::polynomial::Polynomial::interpolate_barycentric_eval",
    DimensionPolicy::Interpolation,
    (nodes: Vec<f64>, values: Vec<f64>, query: f64) -> Result<f64, OpError>,
    |nodes: Vec<f64>, values: Vec<f64>, query: f64| {
        if nodes.len() > super::MAX_BARYCENTRIC_NODES {
            return Err(OpError::Domain("barycentric interpolation exceeds the node budget".into()));
        }
        super::RealPolynomial::interpolate_barycentric_eval(&nodes, &values, query).map_err(OpError::Domain)
    }
);

binding!(
    polynomial_power,
    "num_real_polynomial_power",
    "numerical_special",
    "Raises a polynomial with dimensionless descending coefficients to a nonnegative integer power. The result contains at most 4097 coefficients.",
    "crate::kernels::polynomial::Polynomial::pow",
    DimensionPolicy::Dimensionless,
    (coefficients: Vec<f64>, exponent: usize) -> Result<Vec<f64>, OpError>,
    |coefficients: Vec<f64>, exponent: usize| {
        if coefficients.is_empty() || coefficients.len() > super::MAX_POLYNOMIAL_COEFFICIENTS {
            return Err(OpError::Domain("polynomial coefficients exceed the supported shape".into()));
        }
        let polynomial = super::RealPolynomial::new(coefficients);
        if polynomial.degree().saturating_mul(exponent) >= super::MAX_POLYNOMIAL_COEFFICIENTS {
            return Err(OpError::Domain("polynomial power exceeds the coefficient budget".into()));
        }
        Ok(polynomial.pow(exponent).coeffs)
    }
);

binding!(
    polynomial_gcd,
    "num_real_polynomial_gcd",
    "numerical_special",
    "Returns the monic numerical greatest common divisor of two dimensionless coefficient vectors. Tolerance is relative to normalized coefficients and must lie in [0, 1).",
    "crate::kernels::polynomial::Polynomial::gcd",
    DimensionPolicy::Dimensionless,
    (left: Vec<f64>, right: Vec<f64>, relative_tolerance: f64) -> Result<Vec<f64>, OpError>,
    |left: Vec<f64>, right: Vec<f64>, relative_tolerance: f64| {
        if left.is_empty() || right.is_empty()
            || left.len().max(right.len()) > super::MAX_GCD_COEFFICIENTS
            || !(0.0..1.0).contains(&relative_tolerance)
        {
            return Err(OpError::Domain("polynomial GCD requires bounded nonempty coefficient vectors and a relative tolerance in [0, 1)".into()));
        }
        let left = super::normalized_polynomial(left);
        let right = super::normalized_polynomial(right);
        Ok(left.gcd(&right, relative_tolerance).coeffs)
    }
);
