//! Analytic checks for prepared simplex residual indicators.

use rssn::sim::boundary::MarkerBcs;
use rssn::sim::fem_simplex::FemProblem;
use rssn::sim::fem_simplex::FemProblemVar;
use rssn::sim::fem_simplex::FemSpace;
use rssn::sim::fem_simplex::SimplexMesh;
use rssn::sim::fem_simplex::SimplexResidualEstimator;
use rssn::sim::fem_simplex::error_indicators;
use rssn::sim::fem_simplex::error_indicators_var;

fn custom_space(
    dim: usize,
    degree: usize,
    nodes: Vec<[f64; 3]>,
    cells: Vec<Vec<usize>>,
) -> FemSpace {
    let level = vec![0; cells.len()];
    FemSpace::new(
        SimplexMesh {
            dim,
            nodes,
            cells,
            facets: Vec::new(),
            level,
        },
        degree,
    )
    .unwrap_or_else(|error| panic!("{error}"))
}

fn assert_relative(
    actual: f64,
    expected: f64,
    tolerance: f64,
) {
    assert!(
        actual.is_finite() && (actual / expected - 1.0).abs() <= tolerance,
        "got {actual:e}, expected {expected:e}"
    );
}

#[test]
fn p2_face_integral_is_invariant_under_large_translation() {
    let bcs = MarkerBcs::new();
    let problem = FemProblem {
        k: 1.0,
        c: 0.0,
        f: &|_, _| 0.0,
        bcs: &bcs,
    };
    // The upper triangle is u = s y, while the lower triangle is zero.
    // The common unit edge has jump s. Each cell receives half the
    // integral of s² on the edge, and the volume residual is zero.
    let expected = (1.0_f64 / 6.0).sqrt();
    for offset in [0.0, 2.0_f64.powi(48), -2.0_f64.powi(48)] {
        let space = custom_space(
            2,
            2,
            vec![
                [offset, 0.0, 0.0],
                [offset + 1.0, 0.0, 0.0],
                [offset, 1.0, 0.0],
                [offset, -1.0, 0.0],
            ],
            vec![vec![0, 1, 2], vec![0, 3, 1]],
        );
        let u = space.interpolate(|x| (x[0] - offset) * x[1].max(0.0));
        let prepared = SimplexResidualEstimator::new(&space).unwrap_or_else(|e| panic!("{e}"));
        let eta = prepared
            .indicators(&u, &problem, 0.0)
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(eta, error_indicators(&space, &u, &problem, 0.0));
        assert_eq!(eta.len(), 2);
        for value in eta {
            assert_relative(value, expected, 8e-15);
        }
    }
}

#[test]
fn tetrahedron_residual_keeps_canceling_geometric_and_source_scales() {
    let bcs = MarkerBcs::new();
    for scale in [1e-170_f64, 1e-150, 1e-100, 1e100, 1e150, 1e170] {
        let space = custom_space(
            3,
            1,
            vec![
                [0.0; 3],
                [scale, 0.0, 0.0],
                [0.0, scale, 0.0],
                [0.0, 0.0, scale],
            ],
            vec![vec![0, 1, 2, 3]],
        );
        let source = |_, _| 1.0 / scale;
        let problem = FemProblem {
            k: 1.0,
            c: 0.0,
            f: &source,
            bcs: &bcs,
        };
        let eta = error_indicators(&space, &[0.0; 4], &problem, 0.0);
        // h = sqrt(2)s, V = s³/6 and f = 1/s.
        let expected = scale * scale.sqrt() / 3.0_f64.sqrt();
        assert_relative(eta[0], expected, 1e-14);
    }
}

#[test]
fn interior_tetrahedron_jump_keeps_the_unit_normal_at_extreme_scales() {
    let bcs = MarkerBcs::new();
    let problem = FemProblem {
        k: 1.0,
        c: 0.0,
        f: &|_, _| 0.0,
        bcs: &bcs,
    };
    for scale in [1e-160_f64, 1e-100, 1.0, 1e100, 1e160] {
        let space = custom_space(
            3,
            1,
            vec![
                [0.0; 3],
                [scale, 0.0, 0.0],
                [0.0, scale, 0.0],
                [0.0, 0.0, scale],
                [0.0, 0.0, -scale],
            ],
            vec![vec![0, 1, 2, 3], vec![0, 2, 1, 4]],
        );
        let u = space.interpolate(|x| x[2].max(0.0));
        let eta = error_indicators(&space, &u, &problem, 0.0);
        // Constant jump 1, face area s²/2 and face diameter sqrt(2)s.
        let expected = scale * scale.sqrt() * (std::f64::consts::SQRT_2 * 0.25).sqrt();
        assert_relative(eta[0], expected, 8e-15);
        assert_eq!(eta[0], eta[1]);
    }
}

#[test]
fn variable_residual_uses_distinct_representable_difference_points() {
    let offset = 2.0_f64.powi(48);
    let space = custom_space(
        1,
        1,
        vec![[offset, 0.0, 0.0], [offset + 1.0, 0.0, 0.0]],
        vec![vec![0, 1]],
    );
    let bcs = MarkerBcs::new();
    let k = |x: [f64; 3]| 1.0 + (x[0] - offset);
    let problem = FemProblemVar {
        k: &k,
        c: &|_| 0.0,
        f: &|_, _| -1.0,
        bcs: &bcs,
    };
    // u_x = 1 and (k u_x)_x = 1, despite the physical-coordinate ulp being
    // larger than the nominal difference step on this translated mesh.
    let u = [0.0, 1.0];
    let eta = error_indicators_var(&space, &u, &problem, 0.0);
    assert!(eta[0].abs() < 1e-14, "{eta:?}");
    let prepared = SimplexResidualEstimator::new(&space).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(prepared.indicators_var(&u, &problem, 0.0), Ok(eta));
}

#[test]
fn prepared_estimator_owns_the_geometry_and_validates_shapes() {
    let mut space = custom_space(1, 1, vec![[0.0; 3], [1.0, 0.0, 0.0]], vec![vec![0, 1]]);
    let prepared = SimplexResidualEstimator::new(&space).unwrap_or_else(|e| panic!("{e}"));
    let bcs = MarkerBcs::new();
    let problem = FemProblem {
        k: 1.0,
        c: 0.0,
        f: &|t, _| t,
        bcs: &bcs,
    };
    assert!(prepared.indicators(&[0.0], &problem, 0.0).is_err());
    for time in [0.0, 1.0, 2.0] {
        let eta = prepared
            .indicators(&[0.0, 0.0], &problem, time)
            .unwrap_or_else(|e| panic!("{e}"));
        assert!((eta[0] - time).abs() < 2e-14);
    }

    space.mesh.nodes[1][0] = 0.0;
    assert!(SimplexResidualEstimator::new(&space).is_err());
    assert!(error_indicators(&space, &[0.0, 0.0], &problem, 1.0)[0].is_nan());
    let eta = prepared
        .indicators(&[0.0, 0.0], &problem, 2.0)
        .unwrap_or_else(|e| panic!("{e}"));
    assert!((eta[0] - 2.0).abs() < 2e-14);
}
