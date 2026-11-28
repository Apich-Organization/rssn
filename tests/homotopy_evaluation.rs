//! Analytic polynomial systems for the projective path evaluator.

use rssn::kernels::homotopy::HomotopyOptions;
use rssn::kernels::homotopy::Precision;
use rssn::kernels::homotopy::homotopy_solve;
use rssn::kernels::rootfind::PolyTerm;

fn term(
    coef: f64,
    exps: &[u32],
) -> PolyTerm {
    PolyTerm {
        coef,
        exps: exps.to_vec(),
    }
}

#[test]
fn analytic_roots_include_zero_coordinates_with_nonzero_derivatives() {
    // x^3 - x = 0 and y - x = 0 give the three finite roots x = y = -1, 0, 1.
    let system = vec![
        vec![term(1.0, &[3, 0]), term(-1.0, &[1, 0])],
        vec![term(1.0, &[0, 1]), term(-1.0, &[1, 0])],
    ];
    for precision in [Precision::Double, Precision::DoubleDouble] {
        let opts = HomotopyOptions {
            precision,
            ..HomotopyOptions::default()
        };
        let mut ends = homotopy_solve(&system, &opts).unwrap();
        ends.sort_by(|a, b| a.x[0].0.total_cmp(&b.x[0].0));
        assert_eq!(ends.len(), 3);
        for (path, expected) in ends.iter().zip([-1.0, 0.0, 1.0]) {
            assert!(path.converged && !path.at_infinity && !path.singular && !path.suspect);
            assert_eq!(path.winding, 1);
            assert!(path.residual < 1e-12);
            for (&(re, im), &(re_lo, im_lo)) in path.x.iter().zip(&path.x_lo) {
                assert!(((re + re_lo) - expected).hypot(im + im_lo) < 1e-12);
            }
        }
    }
}

#[test]
fn paths_keep_start_order_across_worker_counts_and_precision_modes() {
    // The target equals the start system. Start order is the roots of unity
    // in increasing argument, even while projective coordinates use a patch.
    let system = vec![vec![term(1.0, &[4]), term(-1.0, &[0])]];
    let serial = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let parallel = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    for precision in [
        Precision::Double,
        Precision::DoubleDouble,
        Precision::Adaptive,
    ] {
        let opts = HomotopyOptions {
            precision,
            ..HomotopyOptions::default()
        };
        let expected = serial.install(|| homotopy_solve(&system, &opts)).unwrap();
        let actual = parallel.install(|| homotopy_solve(&system, &opts)).unwrap();
        assert_eq!(actual, expected);
        assert_eq!(actual.len(), 4);
        for (path, (rx, ry)) in
            actual
                .iter()
                .zip([(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)])
        {
            let (real, imag) = path.x[0];
            assert!((real - rx).hypot(imag - ry) < 1e-12, "{path:?}");
            assert!(path.converged && !path.singular && !path.suspect && !path.at_infinity);
        }
    }
}
