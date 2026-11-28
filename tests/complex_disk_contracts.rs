//! Portable measured disks and outward bounds for their checked arithmetic.

use num_complex::Complex64;
use num_rational::BigRational;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::overload::matches_type;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
use rssn::kernels::complex_ball::ComplexBall;

fn measured(
    value: Value,
    dim: DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(value),
        dim,
    })
}

fn disk(
    coords: [f64; 3],
    dim: DimExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![3],
        elem: TypeExpr::Float.with_dim(dim),
        data: coords.into_iter().map(Value::Float).collect(),
    })
}

fn real(value: &Value) -> f64 {
    if let Value::Dyn(value) = value {
        return real(&value.value);
    }
    match value {
        | Value::Float(n) => *n,
        | _ => panic!("expected a real result"),
    }
}

fn coords(
    value: &Value,
    dim: DimExpr,
) -> [f64; 3] {
    let Value::Tensor(tensor) = value else {
        panic!("expected a disk")
    };
    assert_eq!(tensor.shape, [3]);
    assert!(matches_type(&TypeExpr::Float.with_dim(dim), &tensor.elem).accepted());
    std::array::from_fn(|index| real(&tensor.data[index]))
}

#[test]
fn portable_disks_have_exact_shape_and_physical_dimensions() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let new = session
        .call_value(
            "num_complex_disk_point",
            &[
                measured(Value::Float(3.0), length.clone()),
                measured(Value::Float(4.0), length.clone()),
            ],
            &ctx,
        )
        .unwrap();
    assert_eq!(coords(&new, length.clone()), [3.0, 4.0, 0.0]);
    let a = disk([3.0, 4.0, 0.125], length.clone());
    let b = disk([2.0, 0.0, 0.125], time.clone());
    assert!(
        session
            .call_value("num_complex_disk_add", &[a.clone(), b.clone()], &ctx)
            .is_err()
    );
    let product = session
        .call_value("num_complex_disk_mul", &[a.clone(), b.clone()], &ctx)
        .unwrap();
    let product_coords = coords(&product, DimExpr::Mul(vec![length.clone(), time.clone()]));
    assert_eq!(product_coords[..2], [6.0, 8.0]);
    assert!(product_coords[2] >= 0.890_625);
    let quotient = session
        .call_value("num_complex_disk_div", &[a.clone(), b], &ctx)
        .unwrap();
    let Value::Variant(quotient) = quotient else {
        panic!("expected Option")
    };
    assert_eq!(quotient.tag, "Some");
    let quotient_coords = coords(
        &quotient.payload[0],
        DimExpr::Mul(vec![length.clone(), time.pow(-1)]),
    );
    assert!((quotient_coords[0] - 1.5).abs() < 1e-13);
    let zero = disk([0.0, 0.0, 1.0], length.clone());
    let unresolved = session
        .call_value("num_complex_disk_div", &[a.clone(), zero], &ctx)
        .unwrap();
    assert!(matches!(unresolved, Value::Variant(v) if v.tag=="None" && v.payload.is_empty()));
    let negative = disk([3.0, 4.0, -0.125], length);
    assert!(
        session
            .call_value("num_complex_disk_radius", &[negative], &ctx)
            .is_err()
    );
    let malformed = Value::Tensor(TensorValue {
        shape: vec![2],
        elem: TypeExpr::Float,
        data: vec![Value::Float(0.0); 2],
    });
    assert!(
        session
            .call_value("num_complex_disk_radius", &[malformed], &ctx)
            .is_err()
    );
}

fn rational(x: f64) -> BigRational {
    BigRational::from_float(x).unwrap()
}

fn contains_exact(
    ball: ComplexBall,
    re: BigRational,
    im: BigRational,
) {
    assert!(ball.is_finite());
    let re = re - rational(ball.mid.re);
    let im = im - rational(ball.mid.im);
    assert!(&re * &re + &im * &im <= rational(ball.rad).pow(2));
}

#[test]
fn disk_multiplication_encloses_exact_rounded_and_canceled_products() {
    let lhs = ComplexBall::point(1e140, (1e140_f64).next_up());
    let rhs = ComplexBall::point((1e140_f64).next_up(), 1e140);
    let (a, b) = (rational(lhs.mid.re), rational(lhs.mid.im));
    let (c, d) = (rational(rhs.mid.re), rational(rhs.mid.im));
    contains_exact(lhs.mul_ball(&rhs), &a * &c - &b * &d, &a * &d + &b * &c);
    let lhs = ComplexBall::new(Complex64::new(2.0, -3.0), 0.125).unwrap();
    let rhs = ComplexBall::new(Complex64::new(1.0, 2.0), 0.25).unwrap();
    let enclosure = lhs.mul_ball(&rhs);
    for n in [Complex64::new(0.125, 0.0), Complex64::new(0.0, -0.125)] {
        for m in [Complex64::new(-0.25, 0.0), Complex64::new(0.0, 0.25)] {
            let (a, b) = (rational(lhs.mid.re + n.re), rational(lhs.mid.im + n.im));
            let (c, d) = (rational(rhs.mid.re + m.re), rational(rhs.mid.im + m.im));
            contains_exact(enclosure, &a * &c - &b * &d, &a * &d + &b * &c);
        }
    }
}

#[test]
fn disk_reciprocal_encloses_extreme_points_and_rejects_zero() {
    for scale in [1e-260, 1e-150, 1.0, 1e150, 1e260] {
        let point = ComplexBall::point(scale, -0.5 * scale);
        let (re, im) = (rational(point.mid.re), rational(point.mid.im));
        let norm_sq = &re * &re + &im * &im;
        contains_exact(point.recip().unwrap(), re / &norm_sq, -im / norm_sq);
    }
    let ball = ComplexBall::new(Complex64::new(2.0, 0.0), 0.25).unwrap();
    let enclosure = ball.recip().unwrap();
    for x in [1.75, 2.25] {
        contains_exact(enclosure, rational(1.0) / rational(x), rational(0.0));
    }
    assert!(
        ComplexBall::new(Complex64::new(1.0, 0.0), 1.0)
            .unwrap()
            .recip()
            .is_none()
    );
    assert!(ball.contains(Complex64::new(2.25, 0.0)));
    assert!(!ball.contains(Complex64::new(2.26, 0.0)));
    assert!(ComplexBall::point(1.75, 0.0).overlaps(&ball));
    assert!(!ComplexBall::point(1.75, 0.0).inside_strictly(&ball));
}

#[test]
fn disk_jets_enclose_polynomial_values_and_derivatives() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let z = disk([2.0, 0.0, 0.125], DimExpr::Dimless);
    let variable = session
        .call_value("num_complex_disk_jet_variable", &[z], &ctx)
        .unwrap();
    let result = session
        .call_value(
            "num_complex_disk_jet_mul",
            &[variable.clone(), variable.clone()],
            &ctx,
        )
        .unwrap();
    let Value::Record(record) = result else {
        panic!("expected a jet")
    };
    for (label, expected, sample) in [("value", 4.0, 2.125_f64.powi(2)), ("derivative", 4.0, 4.25)]
    {
        let field = &record
            .fields
            .iter()
            .find(|(name, _)| name == label)
            .unwrap()
            .1;
        let [center, imag, radius] = coords(field, DimExpr::Dimless);
        assert_eq!([center, imag], [expected, 0.0]);
        assert!((sample - center).abs() <= radius);
    }
    let out = session
        .call_value(
            "num_complex_disk_jet_powi",
            &[variable, Value::int(3)],
            &ctx,
        )
        .unwrap();
    assert!(matches!(out, Value::Variant(v) if v.tag=="Some"));
}

fn disk_tensor(
    shape: &[usize],
    data: &[[f64; 3]],
    dim: DimExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: shape.to_vec(),
        elem: TypeExpr::Float.with_dim(dim),
        data: data.iter().copied().flatten().map(Value::Float).collect(),
    })
}

#[test]
fn disk_matrix_calls_preserve_shapes_and_solution_units() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let len = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let coeff = disk_tensor(
        &[2, 2, 3],
        &[[2.0, 0.0, 0.0], [0.0; 3], [0.0; 3], [4.0, 0.0, 0.0]],
        time.clone(),
    );
    let rhs = disk_tensor(&[2, 3], &[[2.0, 0.0, 0.01], [8.0, 0.0, 0.02]], len.clone());
    let output = session
        .call_value(
            "num_complex_disk_matvec",
            &[coeff.clone(), rhs.clone()],
            &ctx,
        )
        .unwrap();
    let Value::Tensor(output) = output else {
        panic!("expected disk vector")
    };
    assert_eq!(output.shape, [2, 3]);
    let mixed = DimExpr::Mul(vec![len.clone(), time.clone()]);
    assert!(matches_type(&TypeExpr::Float.with_dim(mixed), &output.elem).accepted());
    assert_eq!(real(&output.data[3]), 32.0);
    let product = session
        .call_value(
            "num_complex_disk_matmul",
            &[coeff.clone(), coeff.clone()],
            &ctx,
        )
        .unwrap();
    let Value::Tensor(product) = product else {
        panic!("expected disk matrix")
    };
    assert_eq!(product.shape, [2, 2, 3]);
    assert_eq!(real(&product.data[9]), 16.0);
    let result = session
        .call_value("num_complex_disk_solve", &[coeff, rhs.clone()], &ctx)
        .unwrap();
    let Value::Variant(solution) = result else {
        panic!("expected optional solution")
    };
    assert_eq!(solution.tag, "Some");
    let Value::Tensor(solution) = &solution.payload[0] else {
        panic!("expected disk vector")
    };
    let velocity = DimExpr::Mul(vec![len, time.clone().pow(-1)]);
    assert!(matches_type(&TypeExpr::Float.with_dim(velocity), &solution.elem).accepted());
    assert_eq!(
        [real(&solution.data[0]), real(&solution.data[3])],
        [1.0, 2.0]
    );
    let singular = disk_tensor(&[2, 2, 3], &[[1.0, 0.0, 0.0]; 4], time);
    let result = session
        .call_value("num_complex_disk_solve", &[singular, rhs], &ctx)
        .unwrap();
    assert!(matches!(result, Value::Variant(v) if v.tag=="None"));
}

#[test]
fn disk_matrix_calls_report_invalid_shape_and_work_limits() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let d = DimExpr::Dimless;
    let matrix = disk_tensor(&[1, 2, 3], &[[0.0; 3]; 2], d.clone());
    let rhs = disk_tensor(&[1, 3], &[[0.0; 3]], d.clone());
    let err = session
        .call_value("num_complex_disk_matvec", &[matrix, rhs], &ctx)
        .unwrap_err();
    assert!(err.to_string().contains("dimensions do not agree"));

    let n = 65;
    let matrix = disk_tensor(&[n, n, 3], &vec![[0.0; 3]; n * n], d.clone());
    let rhs = disk_tensor(&[n, 3], &vec![[0.0; 3]; n], d);
    let err = session
        .call_value("num_complex_disk_solve", &[matrix, rhs], &ctx)
        .unwrap_err();
    assert!(err.to_string().contains("work limit"));
}

#[test]
fn acosh_jet_encloses_the_principal_derivative_on_both_sides() {
    use rssn::kernels::complex_ball::ComplexBallJet;
    let one = Complex64::new(1.0, 0.0);
    for re in [-2.0, 2.0] {
        for im in [-1.0, 1.0] {
            let input = ComplexBall::new(Complex64::new(re, im), 0.002).unwrap();
            let output = ComplexBallJet::variable(input).acosh().unwrap();
            for sample in [input.mid, input.mid + Complex64::new(0.002, 0.0)] {
                let expect = one / ((sample - one).sqrt() * (sample + one).sqrt());
                assert!(output.val.contains(sample.acosh()));
                assert!(output.deriv.contains(expect));
                assert!(!output.deriv.contains(-expect));
            }
        }
    }
    for z in [Complex64::new(-2.0, 0.005), Complex64::new(1.0, 0.005)] {
        let disk = ComplexBall::new(z, 0.01).unwrap();
        assert!(ComplexBallJet::variable(disk).acosh().is_none());
    }
}

#[test]
fn elementary_jets_and_reciprocal_disks_enclose_samples() {
    use rssn::kernels::complex_ball::ComplexBallJet;
    let center = Complex64::new(0.3, 0.2);
    let input = ComplexBall::new(center, 0.001).unwrap();
    let x = ComplexBallJet::variable(input);
    let one = Complex64::new(1.0, 0.0);
    for z in [center, center + Complex64::new(0.0, 0.001)] {
        for (result, value, derivative) in [
            (x.exp(), z.exp(), z.exp()),
            (x.sin(), z.sin(), z.cos()),
            (x.cos(), z.cos(), -z.sin()),
            (x.sinh(), z.sinh(), z.cosh()),
            (x.cosh(), z.cosh(), z.sinh()),
            (x.tan().unwrap(), z.tan(), one / z.cos().powi(2)),
            (x.tanh().unwrap(), z.tanh(), one / z.cosh().powi(2)),
            (x.ln().unwrap(), z.ln(), one / z),
            (x.sqrt().unwrap(), z.sqrt(), one / (2.0 * z.sqrt())),
            (x.asin().unwrap(), z.asin(), one / (one - z * z).sqrt()),
            (x.acos().unwrap(), z.acos(), -one / (one - z * z).sqrt()),
            (x.atan().unwrap(), z.atan(), one / (one + z * z)),
            (x.asinh().unwrap(), z.asinh(), one / (one + z * z).sqrt()),
            (x.atanh().unwrap(), z.atanh(), one / (one - z * z)),
            (x.powc(&x).unwrap(), z.powc(z), z.powc(z) * (z.ln() + one)),
        ] {
            assert!(result.val.is_finite() && result.deriv.is_finite());
            assert!(result.val.contains(value), "{result:?} lost {value}");
            assert!(
                result.deriv.contains(derivative),
                "{result:?} lost {derivative}"
            );
        }
        for (result, expected) in [
            (input.sec(), one / z.cos()),
            (input.csc(), one / z.sin()),
            (input.cot(), z.cos() / z.sin()),
            (input.sech(), one / z.cosh()),
            (input.csch(), one / z.sinh()),
            (input.coth(), z.cosh() / z.sinh()),
        ] {
            assert!(result.unwrap().contains(expected));
        }
    }
    assert!(ComplexBall::zero().csc().is_none());
    assert!(ComplexBall::zero().csch().is_none());
}

#[test]
fn analytic_jet_abi_checks_units_and_distinguishes_rejected_cuts() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let z = disk([-2.0, 1.0, 0.0], DimExpr::Dimless);
    let variable = session
        .call_value("num_complex_disk_jet_variable", &[z], &ctx)
        .unwrap();
    let result = session
        .call_value("num_complex_disk_jet_acosh", &[variable], &ctx)
        .unwrap();
    let Value::Variant(result) = result else {
        panic!("expected option")
    };
    assert_eq!(result.tag, "Some");
    let Value::Record(jet) = &result.payload[0] else {
        panic!("expected jet")
    };
    let derivative = &jet
        .fields
        .iter()
        .find(|(field, _)| field == "derivative")
        .unwrap()
        .1;
    assert!(coords(derivative, DimExpr::Dimless)[0] < 0.0);
    let cut = disk([-2.0, 0.0, 0.0], DimExpr::Dimless);
    let variable = session
        .call_value("num_complex_disk_jet_variable", &[cut], &ctx)
        .unwrap();
    let result = session
        .call_value("num_complex_disk_jet_acosh", &[variable], &ctx)
        .unwrap();
    assert!(matches!(result, Value::Variant(v) if v.tag=="None"));
    let measured = disk([1.0, 0.0, 0.0], DimExpr::base("Time"));
    assert!(
        session
            .call_value("num_complex_disk_sec", &[measured], &ctx)
            .is_err()
    );
}

#[test]
fn principal_log_jet_bounds_the_real_part_lost_by_rounded_hypot() {
    use rssn::kernels::complex_ball::ComplexBallJet;
    let imag = 2.0_f64.powi(-28);
    let z = Complex64::new(1.0, imag);
    let reference = Complex64::new(0.5 * (imag * imag).ln_1p(), imag.atan());
    let output = ComplexBallJet::variable(ComplexBall::of(z)).ln().unwrap();
    assert!(reference.re > 0.0);
    assert!(output.val.contains(reference));
    assert!(output.val.rad < 1e-9);
    assert!(output.deriv.contains(Complex64::new(1.0, 0.0) / z));
}
