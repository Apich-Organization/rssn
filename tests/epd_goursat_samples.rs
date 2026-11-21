//! Sampled Euler–Poisson–Darboux characteristics and their dimensioned ABI.

use rayon::ThreadPoolBuilder;
use rssn::kernels::numeric_common::NumericError;
use rssn::rules::pde::riemann::epd_goursat;
use rssn::rules::pde::riemann::epd_goursat_samples;

const RHO: [f64; 4] = [1.0, 1.4, 2.0, 3.0];
const SIGMA: [f64; 4] = [-0.5, -0.15, 0.3, 0.95];

#[test]
fn affine_analytic_solution_covers_signed_and_extreme_coordinate_scales() {
    let phi = RHO.map(|rho| rho + SIGMA[0]);
    let psi = SIGMA.map(|sigma| RHO[0] + sigma);
    let targets = [
        [2.4, 0.6],
        [2.89, 0.91],
        [1.0, 0.2],
        [2.1, -0.5],
        [1.0, -0.5],
    ];
    for order in [-0.5, 0.0, 1.0, 1.5, 2.5] {
        for scale in [1e-110, 1.0, 1e110] {
            let rho = RHO.map(|p| p * scale);
            let sigma = SIGMA.map(|p| p * scale);
            let scaled = targets.map(|[r, s]| [r * scale, s * scale]);
            let actual = epd_goursat_samples(order, &rho, &phi, &sigma, &psi, &scaled).unwrap();
            for (value, [r, s]) in actual.iter().zip(targets) {
                assert!(
                    (value - (r + s)).abs() < 2e-10,
                    "order {order}, scale {scale}: {value}"
                );
            }
        }
    }
}

#[test]
fn irregular_knots_and_worker_dispatch_preserve_output_order_and_bits() {
    let phi = [0.5, 0.9, 1.8, 2.2];
    let psi = [0.5, 1.2, 1.1, 3.0];
    let targets: Vec<_> = (0..96)
        .map(|i| {
            [
                1.3 + f64::from(i % 13) / 8.0,
                -0.4 + f64::from(i % 17) / 15.0,
            ]
        })
        .collect();
    let execute = |workers| {
        ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap()
            .install(|| {
                epd_goursat_samples(1.5, &RHO, &phi, &SIGMA, &psi, &targets)
                    .unwrap()
                    .into_iter()
                    .map(f64::to_bits)
                    .collect::<Vec<_>>()
            })
    };
    assert_eq!(execute(1), execute(4));
    let [r, s] = [1.7, 0.5];
    let actual = epd_goursat_samples(0.0, &RHO, &phi, &SIGMA, &psi, &[[r, s]]).unwrap()[0];
    let expected = (0.9 + 0.5 * (1.8 - 0.9)) + (1.1 + (s - 0.3) / (0.95 - 0.3) * (3.0 - 1.1)) - 0.5;
    assert!((actual - expected).abs() < 1e-14);
}

#[test]
fn validates_shapes_corner_geometry_and_nonfinite_callbacks() {
    let rho = [1.0, 2.0];
    let sigma = [-1.0, 0.0];
    let data = [3.0, 5.0];
    assert!(matches!(
        epd_goursat_samples(1.5, &rho, &[3.0], &sigma, &data, &[]),
        Err(NumericError::ShapeMismatch(_))
    ));
    for (r, phi, s, psi, target) in [
        (rho, data, sigma, [3.5, 5.0], [1.4, -0.3]),
        ([1.0, 1.0], data, sigma, data, [1.0, -0.3]),
        (rho, [3.0, f64::INFINITY], sigma, data, [1.4, -0.3]),
        (rho, data, [-1.0, 1.1], data, [1.4, 1.0]),
        (rho, data, sigma, data, [1.4, f64::NAN]),
        (rho, data, sigma, data, [2.1, -0.3]),
    ] {
        assert!(matches!(
            epd_goursat_samples(1.5, &r, &phi, &s, &psi, &[target]),
            Err(NumericError::Domain(_))
        ));
    }
    assert_eq!(
        epd_goursat_samples(1.5, &[1.0], &[3.0], &sigma, &data, &[[1.0, -0.5]]).unwrap(),
        [4.0]
    );
    assert!(epd_goursat(1.5, (1.0, -1.0), &|_| f64::NAN, &|_| 1.0, (1.5, 0.0)).is_none());
}

#[test]
fn native_signature_and_session_preserve_measures_and_reject_wrong_domains() {
    use rssn::api::Session;
    use rssn::kernel::CallCtx;
    use rssn::kernel::Value;
    use rssn::kernel::signature::DimExpr;
    use rssn::kernel::signature::TypeExpr;
    use rssn::kernel::value::DynValue;
    use rssn::kernel::value::TensorValue;

    let measured = |shape: &[usize], data: &[f64], kind: &str| {
        let dimension = DimExpr::base(kind);
        Value::Tensor(TensorValue {
            shape: shape.to_vec(),
            data: data
                .iter()
                .copied()
                .map(|x| {
                    Value::Dyn(DynValue {
                        dim: dimension.clone(),
                        value: Box::new(Value::Float(x)),
                    })
                })
                .collect(),
            elem: TypeExpr::Float.with_dim(dimension),
        })
    };
    let call = |session: &Session, args: &[Value]| {
        session.call_value("num_native_epd_goursat_samples", args, &CallCtx::new())
    };
    let session = Session::new();
    let manifest = session.manifest();
    let sig = &manifest
        .op("num_native_epd_goursat_samples")
        .unwrap()
        .overloads[0];
    let r = TypeExpr::Float.with_dim(DimExpr::Var("Y".into()));
    assert!(matches!(&sig.ret, TypeExpr::Tensor { elem, .. } if elem.as_ref() == &r));

    let phi = RHO.map(|rho| rho + SIGMA[0]);
    let psi = SIGMA.map(|sigma| RHO[0] + sigma);
    let mut args = [
        Value::Float(1.5),
        measured(&[4], &RHO, "Length"),
        measured(&[4], &phi, "Time"),
        measured(&[4], &SIGMA, "Length"),
        measured(&[4], &psi, "Time"),
        measured(&[2, 2], &[2.4, 0.6, 1.0, 0.2], "Length"),
    ];
    let Value::Tensor(out) = call(&session, &args).unwrap() else {
        panic!("expected a measured solution tensor");
    };
    assert_eq!(out.shape, [2]);
    assert_eq!(out.elem, TypeExpr::Float.with_dim(DimExpr::base("Time")));
    for (value, want) in out.data.into_iter().zip([3.0, 1.2]) {
        let Value::Dyn(value) = value else {
            panic!("expected a measured sample");
        };
        assert!((value.value.as_f64().unwrap() - want).abs() < 1e-10);
    }
    args[4] = measured(&[4], &psi, "Length");
    assert!(call(&session, &args).is_err());
    args[4] = measured(&[4], &psi, "Time");
    args[5] = measured(&[1, 3], &[2.4, 0.6, 0.0], "Length");
    assert!(call(&session, &args).is_err());
    args[5] = measured(&[257, 2], &[0.0; 514], "Length");
    assert!(call(&session, &args).is_err());
    args[5] = measured(&[0, 2], &[], "Length");
    let Value::Tensor(empty) = call(&session, &args).unwrap() else {
        panic!("expected empty measured solution tensor");
    };
    assert_eq!(empty.shape, [0]);
    assert_eq!(empty.elem, TypeExpr::Float.with_dim(DimExpr::base("Time")));
}
