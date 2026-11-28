//! Checked degree and total path count before homotopy construction.

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
use rssn::kernels::homotopy::HomotopyOptions;
use rssn::kernels::homotopy::homotopy_solve;
use rssn::kernels::rootfind::PolyTerm;
use rssn::kernels::rootfind::RootError;
use rssn::kernels::rootfind::homotopy_system_with_gamma;

fn diagonal_system(degrees: &[u32]) -> Vec<Vec<PolyTerm>> {
    degrees
        .iter()
        .enumerate()
        .map(|(i, &degree)| {
            let mut exps = vec![0; degrees.len()];
            exps[i] = degree;
            vec![
                PolyTerm { coef: 1.0, exps },
                PolyTerm {
                    coef: -1.0,
                    exps: vec![0; degrees.len()],
                },
            ]
        })
        .collect()
}

#[test]
fn too_many_paths_are_rejected_even_when_the_machine_product_wraps() {
    let opts = HomotopyOptions::default();
    for degrees in [&[100_001][..], &[317, 316], &[256; 8]] {
        assert_eq!(
            homotopy_solve(&diagonal_system(degrees), &opts),
            Err(RootError::Invalid)
        );
    }
    // The exact product for the third case is 256^8 = 2^64.
    assert_eq!(256_u128.pow(8), 1_u128 << 64);
}

#[test]
fn a_term_degree_above_u32_is_not_a_small_positive_degree() {
    let mut system = diagonal_system(&[2, 2]);
    system[0].push(PolyTerm {
        coef: 1.0,
        exps: vec![u32::MAX, 2],
    });
    assert_eq!(u64::from(u32::MAX) + 2, (1_u64 << 32) + 1);
    assert_eq!(
        homotopy_solve(&system, &HomotopyOptions::default()),
        Err(RootError::Invalid)
    );
}

#[test]
fn empty_constant_and_misshapen_systems_keep_their_invalid_error() {
    let opts = HomotopyOptions::default();
    assert_eq!(homotopy_solve(&[], &opts), Err(RootError::Invalid));
    for system in [
        vec![Vec::new()],
        vec![vec![PolyTerm {
            coef: 1.0,
            exps: vec![0],
        }]],
        vec![vec![PolyTerm {
            coef: 1.0,
            exps: vec![1, 1],
        }]],
    ] {
        assert_eq!(homotopy_solve(&system, &opts), Err(RootError::Invalid));
    }
}

#[test]
fn accepted_mixed_degrees_keep_start_path_order() {
    let system = diagonal_system(&[2, 3]);
    let opts = HomotopyOptions::default();
    let actual = homotopy_solve(&system, &opts).unwrap();
    let root_im = 3.0_f64.sqrt() * 0.5;
    let expected = [(1.0, 0.0), (-0.5, root_im), (-0.5, -root_im)];
    assert_eq!(actual.len(), 6);
    for (i, path) in actual.iter().enumerate() {
        let (x, y) = (path.x[0], path.x[1]);
        let rx = if i % 2 == 0 { 1.0 } else { -1.0 };
        let (ry, iy) = expected[i / 2];
        assert!((x.0 - rx).hypot(x.1) < 1e-12);
        assert!((y.0 - ry).hypot(y.1 - iy) < 1e-12);
        assert!(path.converged && !path.suspect && !path.at_infinity);
    }
}

#[test]
fn affine_homotopy_validates_degrees_coefficients_and_gamma_before_tracking() {
    for degrees in [&[200_001][..], &[500, 401], &[256; 8]] {
        assert_eq!(
            homotopy_system_with_gamma(&diagonal_system(degrees), (1.0, 0.0)),
            Err(RootError::Invalid)
        );
    }
    let mut system = diagonal_system(&[2, 2]);
    system[0].push(PolyTerm {
        coef: 1.0,
        exps: vec![u32::MAX, 2],
    });
    assert_eq!(
        homotopy_system_with_gamma(&system, (1.0, 0.0)),
        Err(RootError::Invalid)
    );
    let mut system = diagonal_system(&[2]);
    for gamma in [(0.0, -0.0), (f64::NAN, 1.0), (1.0, f64::INFINITY)] {
        assert_eq!(
            homotopy_system_with_gamma(&system, gamma),
            Err(RootError::Invalid)
        );
    }
    system[0][1].coef = f64::INFINITY;
    assert_eq!(
        homotopy_system_with_gamma(&system, (1.0, 0.0)),
        Err(RootError::Invalid)
    );
    assert_eq!(
        homotopy_solve(&system, &HomotopyOptions::default()),
        Err(RootError::Invalid)
    );
}

#[test]
fn affine_path_results_are_ordered_and_bitwise_equal_with_one_or_two_workers() {
    let system = diagonal_system(&[2, 3]);
    let solve = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| homotopy_system_with_gamma(&system, (1.0, 0.0)).unwrap())
    };
    let serial = solve(1);
    let parallel = solve(2);
    assert_eq!(serial.len(), 6);
    assert_eq!(parallel.len(), serial.len());
    let root_im = 3.0_f64.sqrt() * 0.5;
    let y = [(1.0, 0.0), (-0.5, root_im), (-0.5, -root_im)];
    for (index, (a, b)) in serial.iter().zip(&parallel).enumerate() {
        assert!(a.finite && b.finite);
        assert_eq!(a.residual.to_bits(), b.residual.to_bits());
        assert_eq!(a.x.len(), b.x.len());
        for ((ar, ai), (br, bi)) in a.x.iter().zip(&b.x) {
            assert_eq!((ar.to_bits(), ai.to_bits()), (br.to_bits(), bi.to_bits()));
        }
        let x = if index % 2 == 0 { 1.0 } else { -1.0 };
        assert!((a.x[0].0 - x).hypot(a.x[0].1) < 1e-12);
        assert!((a.x[1].0 - y[index / 2].0).hypot(a.x[1].1 - y[index / 2].1) < 1e-12);
    }
}

#[test]
fn native_affine_limits_preserve_the_small_measured_system() {
    let session = Session::new();
    let length = DimExpr::Base("Length".into());
    let solve = |degrees: &[u32]| {
        let variables = degrees.len();
        let len = variables * 2;
        let coefficients = (0..variables)
            .flat_map(|_| [Value::Float(1.0), Value::Float(-1.0)])
            .collect();
        let exponents = diagonal_system(degrees)
            .into_iter()
            .flatten()
            .flat_map(|term| term.exps.into_iter().map(Value::int))
            .collect();
        let offsets = (0..=variables).map(|i| Value::int(i * 2)).collect();
        let tensor = |shape, elem, data| Value::Tensor(TensorValue { shape, elem, data });
        session.call_value(
            "num_analysis_homotopy_system_scaled",
            &[
                tensor(vec![len], TypeExpr::Float, coefficients),
                tensor(vec![len, variables], TypeExpr::Int, exponents),
                tensor(vec![variables + 1], TypeExpr::Int, offsets),
                Value::Dyn(DynValue {
                    dim: length.clone(),
                    value: Box::new(Value::Float(2.0)),
                }),
                Value::Dyn(DynValue {
                    dim: DimExpr::Base("Time".into()),
                    value: Box::new(Value::Float(3.0)),
                }),
            ],
            &CallCtx::new(),
        )
    };
    let Value::Record(solution) = solve(&[1]).unwrap() else {
        panic!("expected a homotopy result");
    };
    let Some(Value::Tensor(roots)) = solution.get("roots") else {
        panic!("expected measured roots");
    };
    assert_eq!(roots.shape, [1]);
    assert_eq!(roots.elem, TypeExpr::Complex.with_dim(length.clone()));
    let Some(Value::Tensor(finite)) = solution.get("finite") else {
        panic!("expected finite flags");
    };
    assert_eq!(finite.data, [Value::Bool(true)]);
    for degrees in [&[17][..], &[16, 16, 2], &[u32::MAX]] {
        assert!(solve(degrees).is_err());
    }
}
