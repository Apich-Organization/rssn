//! Reusable triangular finite-element geometry and assembly.

use std::cell::RefCell;

use rssn::sim::fem_triangles::Degree;
use rssn::sim::fem_triangles::PreparedTriangles;
use rssn::sim::fem_triangles::Problem;
use rssn::sim::fem_triangles::TriMesh;
use rssn::sim::fem_triangles::solve;

#[test]
fn prepared_mesh_owns_geometry_and_evaluates_new_callbacks_in_order() {
    for degree in [Degree::P1, Degree::P2] {
        let mesh = TriMesh::rectangle(0.0, 1.0, 0.0, 1.0, 3, 4);
        let prepared = PreparedTriangles::new(&mesh, degree).unwrap();
        let calls = RefCell::new(Vec::new());
        let record =
            |kind, x: f64, y: f64| calls.borrow_mut().push((kind, x.to_bits(), y.to_bits()));
        let u = |x, y| 1.0 + x + 2.0 * y;
        let c = |_, y| 0.4 + 0.2 * y;
        let problem = Problem {
            k: &|x, y| {
                record(0, x, y);
                1.0 + 0.3 * x
            },
            c: &|x, y| {
                record(1, x, y);
                c(x, y)
            },
            f: &|x, y| {
                record(2, x, y);
                c(x, y) * u(x, y) - 0.3
            },
            g: &|x, y| {
                record(3, x, y);
                u(x, y)
            },
        };
        let first = solve(&mesh, degree, &problem).unwrap();
        let first_calls = calls.take();
        let quadrature_calls = 3 * 7 * mesh.triangles.len();
        let boundary_points = if degree == Degree::P1 {
            14
        } else {
            28
        };
        assert_eq!(first_calls.len(), quadrature_calls + boundary_points);
        for row in first_calls[..quadrature_calls].as_chunks::<3>().0 {
            let (_, x, y) = row[0];
            assert_eq!(*row, [(0, x, y), (1, x, y), (2, x, y)]);
        }
        assert!(
            first_calls[quadrature_calls..]
                .iter()
                .all(|call| call.0 == 3)
        );
        drop(mesh);
        let second = prepared.solve(&problem).unwrap();
        assert_eq!(first_calls, calls.take());
        assert_eq!(first.points, second.points);
        assert_eq!(first.elements, second.elements);
        assert_eq!(prepared.degrees_of_freedom(), second.values.len());
        assert!(prepared.retained_bytes() >= 16 * second.points.len());
        assert!(first.max_nodal_error(&u) < 2e-10);
        assert!(second.max_nodal_error(&u) < 2e-10);

        let v = |x, y| 2.0 - 3.0 * x + y;
        let another = Problem {
            k: &|_, y| 0.5 + y,
            c: &|_, _| 0.75,
            f: &|x, y| 0.75 * v(x, y) - 1.0,
            g: &v,
        };
        assert!(prepared.solve(&another).unwrap().max_nodal_error(&v) < 2e-10);
    }
}

#[test]
fn affine_manufactured_solution_at_extreme_geometry_scales() {
    for (scale, reaction) in [(1e-155, 1e307), (1e155, 1e-310)] {
        let mesh = TriMesh::rectangle(0.0, scale, 0.0, scale, 3, 4);
        let u = |x, y| 1.0 + x / scale + 2.0 * (y / scale);
        for degree in [Degree::P1, Degree::P2] {
            let prepared = PreparedTriangles::new(&mesh, degree).unwrap();
            for reaction in [0.0, reaction] {
                let solution = prepared
                    .solve(&Problem {
                        k: &|_, _| 1.0,
                        c: &|_, _| reaction,
                        f: &|x, y| reaction * u(x, y),
                        g: &u,
                    })
                    .unwrap();
                assert!(solution.values.iter().all(|x| x.is_finite()));
                let error = solution.max_nodal_error(&u);
                assert!(error < 2e-10, "{scale:e} {degree:?} {reaction:e} {error:e}");
            }
        }
    }
}

#[test]
fn p2_midpoints_stay_finite_when_the_endpoints_sum_overflows() {
    let mesh = TriMesh {
        nodes: vec![[9e307, 9e307], [1.2e308, 9e307], [9e307, 1.2e308]],
        triangles: vec![[2, 1, 0]],
    };
    let expected = |x, y| (x / 1e308) + (y / 1e308);
    let solution = solve(
        &mesh,
        Degree::P2,
        &Problem {
            k: &|_, _| 1.0,
            c: &|_, _| 0.0,
            f: &|_, _| 0.0,
            g: &expected,
        },
    )
    .unwrap();
    assert_eq!(solution.points.len(), 6);
    assert!(solution.points.iter().flatten().all(|x| x.is_finite()));
    assert!(solution.values.iter().all(|x| x.is_finite()));
    assert!(solution.max_nodal_error(&expected) < 1e-13);
}

#[test]
fn extreme_loads_keep_relative_residuals_and_nonzero_solutions() {
    let mesh = TriMesh::rectangle(0.0, 1.0, 0.0, 1.0, 2, 2);
    for degree in [Degree::P1, Degree::P2] {
        let plan = PreparedTriangles::new(&mesh, degree).unwrap();
        let solve_for = |load| {
            plan.solve(&Problem {
                k: &|_, _| 1.0,
                c: &|_, _| 0.0,
                f: &|_, _| load,
                g: &|_, _| 0.0,
            })
            .unwrap()
        };
        let reference = solve_for(1.0);
        if degree == Degree::P1 {
            assert!((reference.values[4] - 0.0625).abs() < 1e-15);
        }
        for load in [1e-300, 1e300] {
            let result = solve_for(load);
            assert!(result.iterations > 0);
            for (&actual, &expected) in result.values.iter().zip(&reference.values) {
                assert!(
                    (actual / load - expected).abs() < 2e-13,
                    "{degree:?} {load:e} {actual:e} {expected:e}"
                );
            }
        }
    }
}

#[test]
fn an_unused_free_node_does_not_masquerade_as_a_solved_equation() {
    let mut mesh = TriMesh::rectangle(0.0, 1.0, 0.0, 1.0, 2, 2);
    mesh.nodes.push([10.0, 10.0]);
    let problem = Problem {
        k: &|_, _| 1.0,
        c: &|_, _| 0.0,
        f: &|_, _| 0.0,
        g: &|_, _| 0.0,
    };
    let plan = PreparedTriangles::new(&mesh, Degree::P1).unwrap();
    assert!(plan.solve(&problem).is_err());
}

#[test]
fn malformed_geometry_and_coefficients_return_errors() {
    let mut mesh = TriMesh::rectangle(0.0, 1.0, 0.0, 1.0, 2, 2);
    let prepared = PreparedTriangles::new(&mesh, Degree::P1).unwrap();
    let valid = Problem {
        k: &|_, _| 1.0,
        c: &|_, _| 0.0,
        f: &|_, _| 0.0,
        g: &|_, _| 0.0,
    };
    assert!(
        prepared
            .solve(&Problem {
                k: &|_, _| 0.0,
                ..valid
            })
            .is_err()
    );
    assert!(
        prepared
            .solve(&Problem {
                c: &|_, _| -1.0,
                ..valid
            })
            .is_err()
    );
    assert!(
        prepared
            .solve(&Problem {
                f: &|_, _| f64::NAN,
                ..valid
            })
            .is_err()
    );
    assert!(
        prepared
            .solve(&Problem {
                g: &|_, _| f64::INFINITY,
                ..valid
            })
            .is_err()
    );
    assert!(prepared.solve(&valid).is_ok());

    let saved = mesh.nodes[0];
    mesh.nodes[0][0] = f64::NAN;
    assert!(PreparedTriangles::new(&mesh, Degree::P2).is_err());
    mesh.nodes[0] = saved;
    mesh.triangles[0][0] = mesh.nodes.len();
    assert!(solve(&mesh, Degree::P2, &valid).is_err());

    let collinear = TriMesh {
        nodes: vec![[0.0, 0.0], [1e-155, 1e-155], [2e-155, 2e-155]],
        triangles: vec![[0, 1, 2]],
    };
    assert!(PreparedTriangles::new(&collinear, Degree::P1).is_err());
}
