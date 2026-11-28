//! Independent Laurent, degeneration, and quasiperiod checks for sigma and zeta.
use std::f64::consts::PI;

use num_complex::Complex64;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::value::ComplexValue;
use rssn::kernels::special_complex_ext::elliptic_modular as elliptic;

const ZERO: Complex64 = Complex64::new(0.0, 0.0);
const ONE: Complex64 = Complex64::new(1.0, 0.0);
const I: Complex64 = Complex64::new(0.0, 1.0);

fn close(
    actual: Complex64,
    expected: Complex64,
    tolerance: f64,
) {
    assert!(
        (actual - expected).norm() <= tolerance * expected.norm().max(1.0),
        "actual {actual}, expected {expected}"
    );
}

fn value(point: Complex64) -> Value {
    Value::Complex(ComplexValue::Float {
        re: point.re,
        im: point.im,
    })
}

#[test]
fn origin_normalization_and_tall_lattice_limit() {
    let z = Complex64::new(0.0002, -0.0003);
    let a = 0.5 * ONE;
    let b = 0.5 * I;
    close(elliptic::weierstrass_sigma_complex(z, a, b) / z, ONE, 1e-13);
    close(z * elliptic::weierstrass_zeta_complex(z, a, b), ONE, 1e-12);
    assert_eq!(elliptic::weierstrass_sigma_complex(ZERO, a, b), ZERO);
    assert!(!elliptic::weierstrass_zeta_complex(ZERO, a, b).is_finite());

    let tall = 1e4 * I;
    let point = Complex64::new(0.21, 0.19);
    let w = PI * point;
    close(
        elliptic::weierstrass_sigma_complex(point, a, tall),
        (w * w / 6.0).exp() * w.sin() / PI,
        2e-13,
    );
    close(
        elliptic::weierstrass_zeta_complex(point, a, tall),
        PI * (w / 3.0 + w.cos() / w.sin()),
        2e-13,
    );
}

#[test]
fn square_lattice_quasiperiods_and_basis_invariance() {
    let z = Complex64::new(0.21, -0.13);
    let a = 0.5 * ONE;
    let b = 0.5 * I;
    let zeta = elliptic::weierstrass_zeta_complex(z, a, b);
    let sigma = elliptic::weierstrass_sigma_complex(z, a, b);

    close(
        elliptic::weierstrass_zeta_complex(z + ONE, a, b),
        zeta + PI,
        2e-13,
    );
    close(
        elliptic::weierstrass_zeta_complex(z + I, a, b),
        zeta - PI * I,
        2e-13,
    );
    close(
        elliptic::weierstrass_sigma_complex(z + ONE, a, b),
        -(PI * (z + a)).exp() * sigma,
        2e-13,
    );
    close(
        elliptic::weierstrass_sigma_complex(z + I, a, b),
        -(-PI * I * (z + b)).exp() * sigma,
        2e-13,
    );
    for (first, second) in [(b, -a), (a, b + 73.0 * a), (a + 31.0 * b, b), (a, -b)] {
        close(
            elliptic::weierstrass_zeta_complex(z, first, second),
            zeta,
            2e-12,
        );
        close(
            elliptic::weierstrass_sigma_complex(z, first, second),
            sigma,
            2e-12,
        );
    }
    let dilation = Complex64::new(0.75, 0.44);
    close(
        elliptic::weierstrass_sigma_complex(z * dilation, a * dilation, b * dilation),
        sigma * dilation,
        2e-13,
    );
    close(
        elliptic::weierstrass_zeta_complex(z * dilation, a * dilation, b * dilation),
        zeta / dilation,
        2e-13,
    );
}

#[test]
fn typed_weierstrass_singularities_and_modular_transform() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let a = 0.5 * ONE;
    let b = 0.5 * I;
    let at_lattice = [value(ONE + I), value(a), value(b)];
    assert_eq!(
        session
            .call_value("num_weierstrass_sigma_complex", &at_lattice, &ctx)
            .unwrap(),
        value(ZERO)
    );
    assert!(
        session
            .call_value("num_weierstrass_zeta_complex", &at_lattice, &ctx)
            .is_err()
    );

    let tau = Complex64::new(0.23, 0.81);
    let discriminant = elliptic::modular_discriminant_complex(tau);
    close(
        elliptic::modular_discriminant_complex(tau + ONE),
        discriminant,
        4e-13,
    );
    close(
        elliptic::modular_discriminant_complex(-ONE / tau),
        tau.powi(12) * discriminant,
        4e-13,
    );
    assert!(
        session
            .call_value(
                "num_weierstrass_modular_discriminant_complex",
                &[value(-I)],
                &ctx
            )
            .is_err()
    );
}
