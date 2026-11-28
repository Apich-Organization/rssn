// Exact integer contracts for Bézout coefficients, symbols and prime groups.
struct BezoutIdentity {
    gcd: u64,
    x: i128,
    y: i128,
}
record!(BezoutIdentity {
    gcd: u64,
    x: i128,
    y: i128
});

binding!(
    extended_gcd_u64,
    "num_number_theory_extended_gcd_u64",
    "numerical_special",
    "Given unsigned 64-bit integers a and b, return the exact record {gcd, x, y} satisfying a*x + b*y = gcd(a,b). Coefficients x and y are signed.",
    "crate::kernels::modular::extended_gcd_u64",
    DimensionPolicy::Dimensionless,
    (a: u64, b: u64) -> super::BezoutIdentity,
    |a: u64, b: u64| {
        let (gcd, x, y) = crate::kernels::modular::extended_gcd_u64(a, b);
        super::BezoutIdentity { gcd, x, y }
    }
);
binding!(
    modular_inverse_u64,
    "num_number_theory_mod_inv_u64",
    "numerical_special",
    "Canonical inverse of an unsigned 64-bit integer modulo m >= 2. Returns a domain error when a and m are not coprime.",
    "crate::kernels::modular::mod_inv_u64",
    DimensionPolicy::Dimensionless,
    (a: u64, m: u64) -> Result<u64, OpError>,
    |a: u64, m: u64| {
        super::require_residue_modulus(m)?;
        crate::kernels::modular::mod_inv_u64(a, m)
            .ok_or_else(|| OpError::Domain("value and modulus must be coprime".into()))
    }
);
binding!(
    jacobi_i128,
    "num_number_theory_jacobi_symbol",
    "numerical_special",
    "Jacobi symbol of a signed 128-bit numerator and a positive odd unsigned 64-bit denominator. Returns -1, 0, or 1.",
    "crate::kernels::number_theory::jacobi_symbol",
    DimensionPolicy::Dimensionless,
    (a: i128, n: u64) -> Result<i32, OpError>,
    |a: i128, n: u64| {
        if n.is_multiple_of(2) {
            return Err(OpError::Domain("Jacobi denominator must be positive and odd".into()));
        }
        Ok(super::number_theory::jacobi_symbol(a, n))
    }
);
binding!(
    legendre_u64,
    "num_number_theory_legendre_symbol",
    "numerical_special",
    "Legendre symbol of an unsigned 64-bit numerator and an odd prime denominator. Returns zero when the prime divides the numerator.",
    "crate::kernels::number_theory::legendre_symbol",
    DimensionPolicy::Dimensionless,
    (a: u64, p: u64) -> Result<i32, OpError>,
    |a: u64, p: u64| super::number_theory::legendre_symbol_checked(a, p)
        .ok_or_else(|| OpError::Domain("Legendre denominator must be an odd prime".into()))
);
binding!(
    primitive_root_u64,
    "num_number_theory_primitive_root",
    "numerical_special",
    "Smallest positive primitive root modulo an unsigned 64-bit integer. Returns a domain error if the multiplicative group has no primitive root.",
    "crate::kernels::number_theory::primitive_root",
    DimensionPolicy::Dimensionless,
    (modulus: u64) -> Result<u64, OpError>,
    |modulus: u64| super::number_theory::primitive_root(modulus)
        .ok_or_else(|| OpError::Domain("modulus has no positive primitive root".into()))
);
binding!(
    finite_field_prime_generator,
    "num_finite_field_primitive_root_prime",
    "numerical_vector",
    "Smallest positive primitive root modulo a prime. Accepts 2 and the odd primes representable in u64.",
    "crate::kernels::finite_field::primitive_root_prime",
    DimensionPolicy::Dimensionless,
    (modulus: u64) -> Result<u64, OpError>,
    |modulus: u64| super::finite_field::primitive_root_prime(modulus)
        .ok_or_else(|| OpError::Domain("primitive-root modulus must be prime".into()))
);
binding!(
    finite_field_poly_mul_ntt,
    "num_finite_field_poly_mul_ntt_mod",
    "numerical_vector",
    "Exact coefficient convolution modulo an odd prime. Inputs use the same coefficient order. For a nonempty product the next power of two of the output length must divide p-1. An empty input gives an empty product. Output uses the native tensor size limit.",
    "crate::kernels::finite_field::poly_mul_ntt_mod",
    DimensionPolicy::Dimensionless,
    (a: Vec<u64>, b: Vec<u64>, modulus: u64) -> Result<Vec<u64>, OpError>,
    |a: Vec<u64>, b: Vec<u64>, modulus: u64| super::checked_ntt_product(&a, &b, modulus)
);

fn checked_ntt_product(
    a: &[u64],
    b: &[u64],
    modulus: u64,
) -> Result<Vec<u64>, OpError> {
    let invalid_prime = || OpError::Domain("NTT modulus must be an odd prime".into());
    if modulus < 3 || modulus.is_multiple_of(2) {
        return Err(invalid_prime());
    }
    if a.is_empty() || b.is_empty() {
        return if finite_field::is_prime_u64(modulus) {
            Ok(Vec::new())
        } else {
            Err(invalid_prime())
        };
    }
    let out_len = a.len().saturating_add(b.len()).saturating_sub(1);
    if out_len >= super::native_object_api::VALUE_LIMIT {
        return Err(OpError::Domain(
            "NTT product exceeds the native tensor size limit".into(),
        ));
    }
    let transform_len = out_len.next_power_of_two() as u64;
    if !(modulus - 1).is_multiple_of(transform_len) {
        return Err(OpError::Domain(format!(
            "NTT transform length {transform_len} must divide modulus minus one"
        )));
    }
    finite_field::poly_mul_ntt_mod(a, b, modulus).ok_or_else(invalid_prime)
}
