//! Accuracy, scaling and ABI regressions for the complex Carlson forms.

use num_complex::Complex64;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::Value;
use rssn::kernels::special_complex_ext::elliptic_modular as kernels;

const ZERO: Complex64 = Complex64::new(0.0, 0.0);
const ONE: Complex64 = Complex64::new(1.0, 0.0);

fn near(
    got: Complex64,
    expected: Complex64,
    relative: f64,
) {
    let residual = (got - expected).norm();
    assert!(
        residual <= relative * expected.norm(),
        "{got} vs {expected}, absolute error {residual}"
    );
}

fn call(
    session: &Session,
    name: &str,
    input: &[Complex64],
) -> Result<Complex64, String> {
    let input: Vec<_> = input
        .iter()
        .map(|value| {
            Value::Complex(ComplexValue::Float {
                re: value.re,
                im: value.im,
            })
        })
        .collect();
    let result = session
        .call_value(
            &format!("num_special_complex_carlson_{name}"),
            &input,
            &CallCtx::new(),
        )
        .map_err(|error| error.to_string())?;
    let Value::Complex(value) = result else {
        panic!("expected a Complex result, received {result:?}");
    };
    let (re, im) = value.to_f64();
    Ok(Complex64::new(re, im))
}

#[test]
fn complex_rj_allows_a_zero_argument_and_counts_the_full_correction() {
    // 65-digit independent defining integral, rounded to binary64.
    let x = Complex64::new(1.0, 2.0);
    let y = Complex64::new(2.0, -2.0);
    let p = Complex64::new(3.0, -3.0);
    let expected = Complex64::new(0.377_872_701_513_904_8, 0.211_743_082_815_471_4);
    near(kernels::carlson_rj_complex(x, y, ZERO, p), expected, 4e-15);
    let session = Session::new();
    near(
        call(&session, "rj", &[x, y, ZERO, p]).unwrap(),
        expected,
        4e-15,
    );
    assert!(call(&session, "rj", &[x, y, ZERO, Complex64::new(-4.0, 0.75)]).is_err());
    for p in [x, y] {
        let other = if p == x { y } else { x };
        assert_eq!(
            call(&session, "rj", &[x, y, ZERO, p]).unwrap(),
            call(&session, "rd", &[other, ZERO, p]).unwrap(),
        );
    }
}

#[test]
fn the_five_routes_use_their_principal_diagonals_and_report_poles() {
    let session = Session::new();
    let four = Complex64::new(4.0, 0.0);
    let half = Complex64::new(0.5, 0.0);
    let eighth = Complex64::new(0.125, 0.0);
    assert_eq!(call(&session, "rf", &[four; 3]).unwrap(), half);
    assert_eq!(call(&session, "rc", &[four; 2]).unwrap(), half);
    assert_eq!(call(&session, "rd", &[four; 3]).unwrap(), eighth);
    assert_eq!(call(&session, "rj", &[four; 4]).unwrap(), eighth);
    assert_eq!(call(&session, "rg", &[four, ZERO, ZERO]).unwrap(), ONE);
    for (name, args) in [
        ("rc", vec![ONE, ZERO]),
        ("rf", vec![ONE, ZERO, ZERO]),
        ("rd", vec![ONE, ONE, ZERO]),
        ("rj", vec![ONE, ONE, ONE, ZERO]),
    ] {
        assert!(call(&session, name, &args).is_err());
    }
}

#[test]
fn complex_carlson_homogeneity_preserves_extreme_results() {
    let x = Complex64::new(0.3, 0.5);
    let y = Complex64::new(0.7, -0.08);
    let z = Complex64::new(0.3, 0.9);
    let rf = kernels::carlson_rf_complex(x, y, z);
    let rd = kernels::carlson_rd_complex(x, y, z);
    let rg = kernels::carlson_rg_complex(x, y, z);
    let large = 1.0e190_f64;
    let tiny = 1.0e-200_f64;
    near(
        kernels::carlson_rf_complex(large * x, large * y, large * z) * large.sqrt(),
        rf,
        4e-15,
    );
    near(
        kernels::carlson_rd_complex(large * x, large * y, large * z) * large * large.sqrt(),
        rd,
        4e-15,
    );
    near(
        kernels::carlson_rd_complex(tiny * x, tiny * y, tiny * z) * tiny * tiny.sqrt(),
        rd,
        4e-15,
    );
    near(
        kernels::carlson_rg_complex(large * x, large * y, large * z) / large.sqrt(),
        rg,
        4e-15,
    );
}

#[test]
fn duplication_avoids_cancelling_nearly_opposite_square_roots() {
    let x = Complex64::new(-1.0, 1e-4);
    let y = x.conj();
    let z = Complex64::new(-1.0, 3e-4);
    let expected = Complex64::new(131.102_827_714_607_32, -130.102_877_722_356);
    near(kernels::carlson_rf_complex(x, y, z), expected, 1e-12);
    near(
        kernels::carlson_rf_complex(x.conj(), y.conj(), z.conj()),
        expected.conj(),
        1e-12,
    );
}
