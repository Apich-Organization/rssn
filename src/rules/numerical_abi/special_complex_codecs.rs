use crate::kernels::special_complex_ext::orthopoly_hyper as orthogonal_complex;

fn complex_polynomial_degree(n: u32) -> Result<u32, OpError> {
    if n > 4096 {
        return Err(OpError::Domain(
            "the complex polynomial ABI supports degrees up to 4096".into(),
        ));
    }
    Ok(n)
}

macro_rules! complex_polynomial {
    ($id:ident, $name:literal, $doc:literal, $path:literal, $kernel:ident) => {
        binding!(
            $id, $name, "numerical_special", $doc, $path, DimensionPolicy::Dimensionless,
            (n: u32, z: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
            |n: u32, z: num_complex::Complex64| {
                Ok(super::orthogonal_complex::$kernel(super::complex_polynomial_degree(n)?, z))
            }
        );
    };
}

complex_polynomial!(
    special_complex_legendre_p,
    "num_special_complex_legendre_p",
    "Legendre polynomial P_n(z) for a dimensionless complex point and integer degree 0 through 4096.",
    "crate::kernels::special_complex_ext::orthopoly_hyper::legendre_p_complex",
    legendre_p_complex
);
complex_polynomial!(
    special_complex_chebyshev_t,
    "num_special_complex_chebyshev_t",
    "First-kind Chebyshev polynomial T_n(z) for a dimensionless complex point and integer degree 0 through 4096.",
    "crate::kernels::special_complex_ext::orthopoly_hyper::chebyshev_t_complex",
    chebyshev_t_complex
);
complex_polynomial!(
    special_complex_chebyshev_u,
    "num_special_complex_chebyshev_u",
    "Second-kind Chebyshev polynomial U_n(z) for a dimensionless complex point and integer degree 0 through 4096.",
    "crate::kernels::special_complex_ext::orthopoly_hyper::chebyshev_u_complex",
    chebyshev_u_complex
);
complex_polynomial!(
    special_complex_hermite_h,
    "num_special_complex_hermite_h",
    "Physicists’ Hermite polynomial H_n(z) for a dimensionless complex point and integer degree 0 through 4096.",
    "crate::kernels::special_complex_ext::orthopoly_hyper::hermite_h_complex",
    hermite_h_complex
);

binding!(
    special_complex_associated_legendre_p,
    "num_special_complex_associated_legendre_p",
    "numerical_special",
    "Associated Legendre P_l^m(z) with the Condon–Shortley phase and principal sqrt(1-z²). Degree 0 through 4096. Orders m > l return zero.",
    "crate::kernels::special_complex_ext::orthopoly_hyper::assoc_legendre_p_complex",
    DimensionPolicy::Dimensionless,
    (l: u32, m: u32, z: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |l: u32, m: u32, z: num_complex::Complex64| {
        Ok(super::orthogonal_complex::assoc_legendre_p_complex(super::complex_polynomial_degree(l)?, m, z))
    }
);
binding!(
    special_complex_laguerre_l,
    "num_special_complex_laguerre_l",
    "numerical_special",
    "Generalized Laguerre L_n^(alpha)(z) with dimensionless complex parameter and point. Degree 0 through 4096.",
    "crate::kernels::special_complex_ext::orthopoly_hyper::laguerre_l_complex",
    DimensionPolicy::Dimensionless,
    (n: u32, alpha: num_complex::Complex64, z: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |n: u32, alpha: num_complex::Complex64, z: num_complex::Complex64| {
        Ok(super::orthogonal_complex::laguerre_l_complex(super::complex_polynomial_degree(n)?, alpha, z))
    }
);
binding!(
    special_complex_jacobi_p,
    "num_special_complex_jacobi_p",
    "numerical_special",
    "Jacobi P_n^(alpha,beta)(z) with dimensionless complex parameters and point. Degree 0 through 4096.",
    "crate::kernels::special_complex_ext::orthopoly_hyper::jacobi_p_complex",
    DimensionPolicy::Dimensionless,
    (n: u32, alpha: num_complex::Complex64, beta: num_complex::Complex64, z: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |n: u32, alpha: num_complex::Complex64, beta: num_complex::Complex64, z: num_complex::Complex64| {
        Ok(super::orthogonal_complex::jacobi_p_complex(super::complex_polynomial_degree(n)?, alpha, beta, z))
    }
);
binding!(
    special_complex_gegenbauer_c,
    "num_special_complex_gegenbauer_c",
    "numerical_special",
    "Gegenbauer C_n^(lambda)(z) with dimensionless complex parameter and point. Degree 0 through 4096.",
    "crate::kernels::special_complex_ext::orthopoly_hyper::gegenbauer_c_complex",
    DimensionPolicy::Dimensionless,
    (n: u32, lambda: num_complex::Complex64, z: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |n: u32, lambda: num_complex::Complex64, z: num_complex::Complex64| {
        Ok(super::orthogonal_complex::gegenbauer_c_complex(super::complex_polynomial_degree(n)?, lambda, z))
    }
);
binding!(
    special_complex_spherical_harmonic,
    "num_special_complex_spherical_harmonic",
    "numerical_special",
    "Normalized spherical harmonic Y_l^m(theta, phi) with Condon–Shortley phase. Theta is colatitude, phi is azimuth, in radians. Degree 0 through 4096. Orders |m| > l return zero.",
    "crate::kernels::special_complex_ext::orthopoly_hyper::spherical_harmonic_complex",
    DimensionPolicy::Dimensionless,
    (l: u32, m: i32, theta: num_complex::Complex64, phi: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |l: u32, m: i32, theta: num_complex::Complex64, phi: num_complex::Complex64| {
        Ok(super::orthogonal_complex::spherical_harmonic_complex(super::complex_polynomial_degree(l)?, m, theta, phi))
    }
);

use crate::kernels::special_complex_ext::erf_integrals;
use crate::kernels::special_complex_ext::gamma_beta;
use crate::kernels::special_complex_ext::zeta_polylog;

fn bernoulli_family_degree(n: u32) -> Result<u32, OpError> {
    if n > 128 {
        return Err(OpError::Domain(
            "the complex Bernoulli and Euler polynomial ABI supports degrees up to 128".into(),
        ));
    }
    Ok(n)
}

binding!(
    special_complex_bernoulli,
    "num_special_complex_bernoulli",
    "numerical_special",
    "Bernoulli polynomial B_n(z), with B_1(0) = -1/2. Degree 0 through 128.",
    "crate::kernels::special_complex_ext::zeta_polylog::bernoulli_poly_complex",
    DimensionPolicy::Dimensionless,
    (n: u32, z: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |n: u32, z: num_complex::Complex64| {
        Ok(super::zeta_polylog::bernoulli_poly_complex(super::bernoulli_family_degree(n)?, z))
    }
);

binding!(
    special_complex_euler,
    "num_special_complex_euler",
    "numerical_special",
    "Euler polynomial E_n(z), with E_1(z) = z - 1/2. Degree 0 through 128.",
    "crate::kernels::special_complex_ext::zeta_polylog::euler_poly_complex",
    DimensionPolicy::Dimensionless,
    (n: u32, z: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |n: u32, z: num_complex::Complex64| {
        Ok(super::zeta_polylog::euler_poly_complex(super::bernoulli_family_degree(n)?, z))
    }
);

binding!(
    special_complex_clausen,
    "num_special_complex_clausen",
    "numerical_special",
    "Complex Clausen function using the antisymmetric polylog combination for even order and the symmetric combination for odd order. Integer order 0 through 4096.",
    "crate::kernels::special_complex_ext::zeta_polylog::clausen_complex",
    DimensionPolicy::Dimensionless,
    (n: u32, z: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |n: u32, z: num_complex::Complex64| {
        Ok(super::zeta_polylog::clausen_complex(super::complex_polynomial_degree(n)?, z))
    }
);

binding!(
    special_complex_log_barnes_g,
    "num_special_complex_log_barnes_g",
    "numerical_special",
    "Analytic log Barnes G, real on the positive axis. The signed imaginary zero selects a side of the non-positive real cut. Zeros of G return a domain error.",
    "crate::kernels::special_complex_ext::gamma_beta::ln_barnes_g_complex",
    DimensionPolicy::Dimensionless,
    (z: num_complex::Complex64) -> num_complex::Complex64,
    |z: num_complex::Complex64| { super::gamma_beta::ln_barnes_g_complex(z) }
);

struct ComplexFresnelAuxiliary {
    f: Complex64,
    g: Complex64,
}

record!(ComplexFresnelAuxiliary {
    f: Complex64,
    g: Complex64
});

binding!(
    special_complex_fresnel_auxiliary,
    "num_special_complex_fresnel_auxiliary",
    "numerical_special",
    "Fresnel auxiliary functions f and g, returned as the named components of a record.",
    "crate::kernels::special_complex_ext::erf_integrals::fresnel_fg_complex",
    DimensionPolicy::Dimensionless,
    (z: num_complex::Complex64) -> super::ComplexFresnelAuxiliary,
    |z: num_complex::Complex64| {
        let (f, g) = super::erf_integrals::fresnel_fg_complex(z);
        super::ComplexFresnelAuxiliary { f, g }
    }
);

use crate::kernels::special_complex_ext::bessel_airy as cylindrical_complex;

fn finite_special(value: num_complex::Complex64) -> Result<num_complex::Complex64, OpError> {
    value.is_finite().then_some(value).ok_or_else(|| {
        OpError::NotReduced("complex special function is singular, nonfinite or unconverged".into())
    })
}

binding!(
    special_complex_bessel_j_order,
    "num_special_complex_bessel_j_order",
    "numerical_special",
    "Bessel J with complex order and argument on the principal branch.",
    "crate::kernels::special_complex_ext::bessel_airy::bessel_j_complex_order",
    DimensionPolicy::Dimensionless,
    (order: num_complex::Complex64, z: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |order, z| super::finite_special(super::cylindrical_complex::bessel_j_complex_order(order, z))
);

binding!(
    special_complex_bessel_i_order,
    "num_special_complex_bessel_i_order",
    "numerical_special",
    "Modified Bessel I with complex order and argument on the principal branch.",
    "crate::kernels::special_complex_ext::bessel_airy::bessel_i_complex_order",
    DimensionPolicy::Dimensionless,
    (order: num_complex::Complex64, z: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |order, z| super::finite_special(super::cylindrical_complex::bessel_i_complex_order(order, z))
);

binding!(
    special_complex_anger_j,
    "num_special_complex_anger_j",
    "numerical_special",
    "Anger J for real order and complex argument using its finite integral.",
    "crate::kernels::special_complex_ext::bessel_airy::anger_j_complex",
    DimensionPolicy::Dimensionless,
    (order: f64, z: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |order, z| super::finite_special(super::cylindrical_complex::anger_j_complex(order, z))
);

binding!(
    special_complex_weber_e,
    "num_special_complex_weber_e",
    "numerical_special",
    "Weber E for real order and complex argument using its finite integral.",
    "crate::kernels::special_complex_ext::bessel_airy::weber_e_complex",
    DimensionPolicy::Dimensionless,
    (order: f64, z: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |order, z| super::finite_special(super::cylindrical_complex::weber_e_complex(order, z))
);

struct ComplexKelvin {
    ber: num_complex::Complex64,
    bei: num_complex::Complex64,
    ker: num_complex::Complex64,
    kei: num_complex::Complex64,
}
record!(ComplexKelvin {
    ber: num_complex::Complex64,
    bei: num_complex::Complex64,
    ker: num_complex::Complex64,
    kei: num_complex::Complex64
});
binding!(
    special_complex_kelvin,
    "num_special_complex_kelvin",
    "numerical_special",
    "Kelvin ber, bei, ker and kei for real order and complex argument as named components.",
    "crate::kernels::special_complex_ext::bessel_airy::kelvin_complex",
    DimensionPolicy::Dimensionless,
    (order: f64, z: num_complex::Complex64) -> Result<super::ComplexKelvin, OpError>,
    |order, z| {
        let (ber, bei, ker, kei) = super::cylindrical_complex::kelvin_complex(order, z);
        Ok(super::ComplexKelvin {
            ber: super::finite_special(ber)?,
            bei: super::finite_special(bei)?,
            ker: super::finite_special(ker)?,
            kei: super::finite_special(kei)?,
        })
    }
);

binding!(
    special_complex_exp_integral_en,
    "num_special_complex_exp_integral_en",
    "numerical_special",
    "Generalized exponential integral E_n for nonnegative integer orders up to 4096. At the negative real cut, the sign of the imaginary zero selects the boundary value.",
    "crate::kernels::special_complex_ext::erf_integrals::exp_integral_en_complex",
    DimensionPolicy::Dimensionless,
    (order: u32, z: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |order, z| super::finite_special(super::erf_integrals::exp_integral_en_complex(super::complex_polynomial_degree(order)?, z))
);

binding!(
    special_complex_legendre_q,
    "num_special_complex_legendre_q",
    "numerical_special",
    "Principal Legendre Q of the second kind for nonnegative integer degrees up to 4096. The sign of the imaginary zero selects the side of the cut from -1 to 1.",
    "crate::kernels::special_complex_ext::orthopoly_hyper::legendre_q_complex",
    DimensionPolicy::Dimensionless,
    (degree: u32, z: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |degree, z| super::finite_special(super::orthogonal_complex::legendre_q_complex(super::complex_polynomial_degree(degree)?, z))
);
