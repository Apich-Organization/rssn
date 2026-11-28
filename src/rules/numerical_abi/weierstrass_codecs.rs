// Strongly typed calls over the existing dimensionless Weierstrass kernels.
use crate::kernels::special_complex_ext::elliptic_modular as weierstrass;

struct WeierstrassInvariants {
    g2: Complex64,
    g3: Complex64,
}
record!(WeierstrassInvariants {
    g2: Complex64,
    g3: Complex64
});

#[allow(
    clippy::arithmetic_side_effects,
    reason = "calculates the complex lattice ratio"
)]
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

#[allow(
    clippy::arithmetic_side_effects,
    clippy::single_call_fn,
    reason = "isolates the lattice pole check"
)]
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

binding!(
    weierstrass_sigma,
    "num_weierstrass_sigma_complex",
    "numerical_special",
    "Weierstrass sigma for dimensionless complex z and two nonzero, non-collinear half-periods. Sigma vanishes at lattice points and its derivative at zero is one.",
    "crate::kernels::special_complex_ext::elliptic_modular::weierstrass_sigma_complex",
    DimensionPolicy::Dimensionless,
    (z: num_complex::Complex64, omega1: num_complex::Complex64, omega2: num_complex::Complex64) -> num_complex::Complex64,
    |z: num_complex::Complex64, omega1: num_complex::Complex64, omega2: num_complex::Complex64| {
        super::weierstrass::weierstrass_sigma_complex(z, omega1, omega2)
    }
);

binding!(
    weierstrass_zeta,
    "num_weierstrass_zeta_complex",
    "numerical_special",
    "Weierstrass zeta for dimensionless complex z and two nonzero, non-collinear half-periods. Lattice poles return a domain error.",
    "crate::kernels::special_complex_ext::elliptic_modular::weierstrass_zeta_complex",
    DimensionPolicy::Dimensionless,
    (z: num_complex::Complex64, omega1: num_complex::Complex64, omega2: num_complex::Complex64) -> num_complex::Complex64,
    |z: num_complex::Complex64, omega1: num_complex::Complex64, omega2: num_complex::Complex64| {
        super::weierstrass::weierstrass_zeta_complex(z, omega1, omega2)
    }
);

binding!(
    complex_modular_discriminant,
    "num_weierstrass_modular_discriminant_complex",
    "numerical_special",
    "Discriminant (2*pi)^12 * eta(tau)^24 for a dimensionless tau in the open upper half-plane.",
    "crate::kernels::special_complex_ext::elliptic_modular::modular_discriminant_complex",
    DimensionPolicy::Dimensionless,
    (tau: num_complex::Complex64) -> num_complex::Complex64,
    |tau: num_complex::Complex64| { super::weierstrass::modular_discriminant_complex(tau) }
);
