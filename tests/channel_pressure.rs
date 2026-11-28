//! Projection tests for the collocated channel model.

use ndarray::Array2;
use ndarray::ShapeBuilder;
use rssn::sim::models_navier_stokes_fluid::ChannelFlowConfig;
use rssn::sim::models_navier_stokes_fluid::InflowProfile;
use rssn::sim::models_navier_stokes_fluid::WallCondition;
use rssn::sim::models_navier_stokes_fluid::run_channel_flow_with;

fn cfg(mask: &Array2<bool>) -> ChannelFlowConfig<'_> {
    ChannelFlowConfig {
        n: mask.nrows(),
        re: 100.0,
        dt: 1e-3,
        n_iter: 1,
        obstacle_mask: mask,
        inflow: InflowProfile::Uniform(1.0),
        walls: WallCondition::Slip,
    }
}

#[test]
fn one_cell_pressure_and_velocity_match_the_analytic_projection() {
    let mask = Array2::from_elem((3, 3), false);
    for wall in [WallCondition::Slip, WallCondition::NoSlip] {
        for (re, inflow) in [(100.0, 1.0), (20.0, -2.0)] {
            let cfg = ChannelFlowConfig {
                re,
                inflow: InflowProfile::Uniform(inflow),
                walls: wall,
                ..cfg(&mask)
            };
            let h = 0.5;
            let intermediate = cfg.dt * inflow / (re * h * h);
            let expected_p = (inflow - intermediate) * h / (2.0 * cfg.dt);
            let expected_u = (3.0 * intermediate + inflow) / 4.0;
            let (u, v, p) = run_channel_flow_with(&cfg).unwrap();
            assert!((p[[1, 1]] - expected_p).abs() < 2e-13);
            assert!((u[[1, 1]] - expected_u).abs() < 1e-15);
            assert_eq!(p[[1, 0]].to_bits(), p[[1, 1]].to_bits());
            assert_eq!(p[[1, 2]].to_bits(), 0.0_f64.to_bits());
            assert!(v.iter().all(|&x| x.abs() < 1e-15));
        }
    }
}

#[test]
fn a_closed_driven_region_reports_incompatible_flux() {
    for (n, blocked_column) in [(7, 3), (9, 8), (3, 2)] {
        let mask = Array2::from_shape_fn((n, n), |(_, i)| i == blocked_column);
        let error = run_channel_flow_with(&cfg(&mask)).unwrap_err();
        assert!(error.contains("incompatible net flux"), "{error}");
    }
}

#[test]
fn separate_quiet_pockets_have_zero_pressure_gauges_and_zero_velocity() {
    let n = 13;
    let mut mask = Array2::from_elem((n, n), false);
    for (y, x) in [(4, 4), (8, 7)] {
        for j in y - 1..=y + 1 {
            for i in x - 1..=x + 2 {
                let interior = j == y && (i == x || i == x + 1);
                mask[[j, i]] = !interior;
            }
        }
    }
    let cfg = ChannelFlowConfig {
        n_iter: 10,
        ..cfg(&mask)
    };
    let (u, v, p) = run_channel_flow_with(&cfg).unwrap();
    for node in [[4, 4], [4, 5], [8, 7], [8, 8]] {
        assert!(u[node].abs() < 1e-20 && v[node].abs() < 1e-20);
        assert!(p[node].abs() < 1e-20);
    }
    assert!(p[[1, 1]].abs() > 0.1);
    for (node, &solid) in mask.indexed_iter() {
        if solid {
            assert_eq!(p[node].to_bits(), 0.0_f64.to_bits());
        }
    }
}

#[test]
fn closed_quiescent_and_all_solid_fields_take_the_exact_zero_solution() {
    let n = 7;
    for mask in [
        Array2::from_elem((n, n), true),
        Array2::from_shape_fn((n, n), |(_, i)| i == n - 1),
    ] {
        let cfg = ChannelFlowConfig {
            inflow: InflowProfile::Uniform(0.0),
            n_iter: 5,
            ..cfg(&mask)
        };
        let (u, v, p) = run_channel_flow_with(&cfg).unwrap();
        assert!(u.iter().chain(v.iter()).chain(p.iter()).all(|&x| x == 0.0));
    }
}

#[test]
fn invalid_material_step_and_inflow_are_checked_before_the_first_update() {
    let mask = Array2::from_elem((5, 5), false);
    let solver = cfg(&mask);
    for invalid in [
        ChannelFlowConfig {
            dt: 0.0,
            ..solver.clone()
        },
        ChannelFlowConfig {
            dt: f64::NAN,
            ..solver.clone()
        },
        ChannelFlowConfig {
            dt: -1.0,
            ..solver.clone()
        },
        ChannelFlowConfig {
            re: f64::INFINITY,
            ..solver.clone()
        },
        ChannelFlowConfig {
            re: f64::from_bits(1),
            ..solver.clone()
        },
        ChannelFlowConfig {
            re: -10.0,
            ..solver.clone()
        },
        ChannelFlowConfig {
            inflow: InflowProfile::Uniform(f64::NAN),
            ..solver.clone()
        },
        ChannelFlowConfig {
            inflow: InflowProfile::Parabolic(f64::INFINITY),
            ..solver.clone()
        },
        ChannelFlowConfig {
            inflow: InflowProfile::Uniform(f64::MAX),
            ..solver.clone()
        },
    ] {
        assert!(run_channel_flow_with(&invalid).is_err());
    }

    let zero = ChannelFlowConfig {
        dt: 0.0,
        n_iter: 0,
        ..solver
    };
    let (u, v, p) = run_channel_flow_with(&zero).unwrap();
    assert_eq!(u[[2, 0]].to_bits(), 1.0_f64.to_bits());
    assert_eq!(u[[2, 1]].to_bits(), 0.0_f64.to_bits());
    assert!(v.iter().chain(&p).all(|x| x.to_bits() == 0));
}

#[test]
fn regular_obstacle_flow_is_identical_across_active_pools() {
    let n = 45;
    let mask = Array2::from_shape_fn((n, n), |(j, i)| {
        (15..24).contains(&j) && (16..19).contains(&i)
    });
    let cfg = ChannelFlowConfig {
        n_iter: 6,
        dt: 0.001,
        ..cfg(&mask)
    };
    let run = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| run_channel_flow_with(&cfg))
            .unwrap()
    };
    let (u1, v1, p1) = run(1);
    let (u2, v2, p2) = run(2);
    for (a, b) in [(&u1, &u2), (&v1, &v2), (&p1, &p2)] {
        assert!(a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits()));
    }
}

#[test]
fn the_pressure_graph_uses_mask_coordinates_for_column_major_input() {
    let n = 9;
    let blocked = |(j, i)| (3..5).contains(&j) && (4..6).contains(&i);
    let row_major = Array2::from_shape_fn((n, n), blocked);
    let column_major = Array2::from_shape_fn((n, n).f(), blocked);
    let (u1, v1, p1) = run_channel_flow_with(&cfg(&row_major)).unwrap();
    let (u2, v2, p2) = run_channel_flow_with(&cfg(&column_major)).unwrap();
    for (a, b) in [(&u1, &u2), (&v1, &v2), (&p1, &p2)] {
        assert!(a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits()));
    }
}
