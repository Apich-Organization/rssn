// Strongly typed calls over the existing dimensionless Weierstrass kernels.
use crate::kernels::special_complex_ext::elliptic_modular as weierstrass;
use num_complex::Complex64;

struct WeierstrassInvariants {
    g2: Complex64,
    g3: Complex64,
}
record!(WeierstrassInvariants {
    g2: Complex64,
    g3: Complex64
});

fn weierstrass_lattice_ratio(
    omega1: Complex64,
    omega2: Complex64,
) -> Result<Complex64, OpError> {
    // Complex arguments already passed the ABI's shared finite-value check.
    // Keep the supplied orientation for reducing positions. The native
    // kernels normalize the orientation for the theta nome themselves.
    let ratio = omega2 / omega1;
    if !ratio.is_finite() || ratio.im == 0.0 {
        return Err(OpError::Domain(
            "Weierstrass half-periods must be nonzero and non-collinear with a finite ratio".into(),
        ));
    }
    Ok(ratio)
}

fn weierstrass_centered_point(
    z: Complex64,
    omega1: Complex64,
    omega2: Complex64,
    ratio: Complex64,
) -> Result<Complex64, OpError> {
    let relative = (0.5 * z) / omega1;
    let second = (relative.im / ratio.im).round();
    let first = (-second).mul_add(ratio.re, relative.re).round();
    let center = 2.0 * (first * omega1 + second * omega2);
    let reduced = z - center;
    if !reduced.is_finite() {
        return Err(OpError::Domain(
            "Weierstrass lattice coordinates exceed binary64 range".into(),
        ));
    }
    if reduced == Complex64::new(0.0, 0.0) {
        return Err(OpError::Domain(
            "the Weierstrass derivative has a pole at this lattice point".into(),
        ));
    }
    Ok(reduced)
}

binding!(
    weierstrass_p_prime,
    "num_weierstrass_p_prime_complex",
    "numerical_special",
    "Numerically evaluates the Weierstrass derivative at z. The arguments omega1 and omega2 are nonzero, non-collinear half-periods, so the basic periods are 2*omega1 and 2*omega2. All arguments and the result are dimensionless. A detected lattice pole, an invalid lattice or a nonfinite result produces a domain error.",
    "crate::kernels::special_complex_ext::elliptic_modular::weierstrass_p_prime_complex",
    DimensionPolicy::Dimensionless,
    (z: num_complex::Complex64, omega1: num_complex::Complex64, omega2: num_complex::Complex64)
        -> Result<num_complex::Complex64, OpError>,
    |z: num_complex::Complex64, omega1: num_complex::Complex64, omega2: num_complex::Complex64| {
        let ratio = super::weierstrass_lattice_ratio(omega1, omega2)?;
        let point = super::weierstrass_centered_point(z, omega1, omega2, ratio)?;
        Ok(super::weierstrass::weierstrass_p_prime_complex(point, omega1, omega2))
    }
);
binding!(
    weierstrass_invariants,
    "num_weierstrass_invariants_complex",
    "numerical_special",
    "Computes the Weierstrass invariants g2 and g3 as a named complex record. The dimensionless arguments omega1 and omega2 are nonzero, non-collinear half-periods. An invalid lattice or a nonfinite result produces a domain error.",
    "crate::kernels::special_complex_ext::elliptic_modular::weierstrass_invariants_complex",
    DimensionPolicy::Dimensionless,
    (omega1: num_complex::Complex64, omega2: num_complex::Complex64)
        -> Result<super::WeierstrassInvariants, OpError>,
    |omega1: num_complex::Complex64, omega2: num_complex::Complex64| {
        super::weierstrass_lattice_ratio(omega1, omega2)?;
        let (g2, g3) = super::weierstrass::weierstrass_invariants_complex(omega1, omega2);
        Ok(super::WeierstrassInvariants {g2, g3})
    }
);
