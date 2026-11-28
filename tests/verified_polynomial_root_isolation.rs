use num_complex::Complex64 as C;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::TensorValue;
use rssn::kernels::complex_roots::isolate_roots_verified;

fn coefficients(ascending: &[f64]) -> Value {
    Value::Tensor(TensorValue::of(
        vec![ascending.len()],
        ascending
            .iter()
            .map(|&re| Value::Complex(ComplexValue::Float { re, im: 0.0 }))
            .collect(),
        TypeExpr::Complex,
    ))
}

#[test]
fn simple_roots_certify_on_normal_and_extreme_coefficient_scales() {
    for scale in [1.0, f64::MAX, f64::from_bits(1)] {
        let mut coeffs = vec![C::new(0.0, 0.0); 13];
        coeffs[0] = C::new(-scale, 0.0);
        coeffs[12] = C::new(scale, 0.0);
        let disks = isolate_roots_verified(&coeffs, 1e-9);
        assert_eq!(disks.len(), 12);
        for (i, disk) in disks.iter().enumerate() {
            assert!(disk.radius() < 1e-8);
            assert!(disks[..i].iter().all(|other| !other.overlaps(disk)));
        }
        for k in 0..12 {
            let point = C::from_polar(1.0, std::f64::consts::TAU * f64::from(k) / 12.0);
            assert!(disks.iter().any(|disk| disk.contains(point)));
        }
    }
}

#[test]
fn isolated_roots_keep_canonical_disk_shape_and_report_partial_sets() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let call = |c: &[f64]| {
        session.call_value(
            "num_complex_polynomial_isolate_verified",
            &[coefficients(c), Value::Float(1e-8)],
            &ctx,
        )
    };
    for (c, degree, complete, roots) in [
        (&[-1.0, 0.0, 1.0][..], 2, true, 2),
        (&[2.0, -3.0, 0.0, 1.0][..], 3, false, 1),
        (&[1.0][..], 0, true, 0),
    ] {
        let Value::Record(report) = call(c).unwrap() else {
            panic!("missing certificate record")
        };
        assert_eq!(report.get("complete"), Some(&Value::Bool(complete)));
        assert_eq!(report.get("degree"), Some(&Value::int(degree)));
        let Some(Value::Tensor(disks)) = report.get("disks") else {
            panic!("missing disk array")
        };
        assert_eq!(disks.shape, [roots, 3]);
        assert_eq!(disks.elem, TypeExpr::Float);
        if !complete {
            let center = disks.data[0].as_f64().unwrap();
            let radius = disks.data[2].as_f64().unwrap();
            assert!((center + 2.0).abs() <= radius);
        }
    }
    assert!(call(&[0.0]).is_err());
    assert!(call(&vec![1.0; 130]).is_err());
}
