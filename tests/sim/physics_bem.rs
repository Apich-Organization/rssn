//! Boundary element method for the 2D Laplace equation (ported from
//! `physics_bem_test.rs`).

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::sim::physics_bem::*;

fn cfg() -> ProptestConfig {
    ProptestConfig {
        cases: 32,
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

/// Unit square, counter-clockwise from the origin, with `n` nodes per side.
/// Left side u = 100, right side u = 0, top and bottom insulated.
fn conduction_problem(n: usize) -> (Vec<(f64, f64)>, Vec<BoundaryCondition<f64>>) {
    let mut points = Vec::new();
    let mut bcs = Vec::new();
    for i in 0..n {
        points.push((i as f64 / n as f64, 0.0));
        bcs.push(BoundaryCondition::Flux(0.0));
    }
    for i in 0..n {
        points.push((1.0, i as f64 / n as f64));
        bcs.push(BoundaryCondition::Potential(0.0));
    }
    for i in 0..n {
        points.push((1.0 - i as f64 / n as f64, 1.0));
        bcs.push(BoundaryCondition::Flux(0.0));
    }
    for i in 0..n {
        points.push((0.0, 1.0 - i as f64 / n as f64));
        bcs.push(BoundaryCondition::Potential(100.0));
    }
    (points, bcs)
}

fn elements_of(points: &[(f64, f64)]) -> Vec<Element2D> {
    let n = points.len();
    (0..n)
        .map(|i| {
            Element2D::new(
                Vector2D::new(points[i].0, points[i].1),
                Vector2D::new(points[(i + 1) % n].0, points[(i + 1) % n].1),
            )
        })
        .collect()
}

#[test]
fn vector_helpers_and_element_geometry() {
    let v = Vector2D::new(3.0, 4.0);
    assert_eq!(v.norm(), 5.0);
    let w = (v + Vector2D::new(1.0, 1.0)) - Vector2D::new(0.0, 2.0);
    assert_eq!((w.x, w.y), (4.0, 3.0));
    assert_eq!((v * 2.0).x, 6.0);
    let e = Element2D::new(Vector2D::new(0.0, 0.0), Vector2D::new(2.0, 0.0));
    assert_eq!(e.length, 2.0);
    assert_eq!((e.midpoint.x, e.midpoint.y), (1.0, 0.0));
    // Counter-clockwise traversal => the normal points to the right of the direction of travel, i.e. outward.
    assert_eq!((e.normal.x, e.normal.y), (0.0, -1.0));
    let d = Element2D::new(Vector2D::new(1.0, 1.0), Vector2D::new(1.0, 3.0));
    assert_eq!((d.normal.x, d.normal.y), (1.0, 0.0));
    assert!((Vector3D::new(2.0, 3.0, 6.0).norm() - 7.0).abs() < 1e-12);
}

#[test]
fn rectangle_conduction_matches_the_linear_profile() {
    let n = 10;
    let (points, bcs) = conduction_problem(n);
    let (u, q) = solve_laplace_bem_2d(&points, &bcs).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!((u.len(), q.len()), (4 * n, 4 * n));
    // Prescribed values are returned untouched.
    assert_approx_eq!(u[n + n / 2], 0.0, 1e-12);
    assert_approx_eq!(u[3 * n + n / 2], 100.0, 1e-12);
    assert_approx_eq!(q[n / 2], 0.0, 1e-12);
    // Exact solution u = 100 (1 - x): the insulated sides read 100 (1 - x_mid) at their elements.
    for i in 0..n {
        let x_mid = (i as f64 + 0.5) / n as f64;
        assert!(
            (u[i] - 100.0 * (1.0 - x_mid)).abs() < 5.0,
            "bottom node {i}: {} vs {}",
            u[i],
            100.0 * (1.0 - x_mid)
        );
    }
    // The flux out of the cold side is u_x = -100 (pointing outward: q = du/dn = -du/dx = 100 at x = 1 ... sign per convention).
    let cold: f64 = q[n..2 * n].iter().sum::<f64>() / n as f64;
    assert!(
        cold.abs() > 50.0 && cold.abs() < 150.0,
        "mean flux on the cold side {cold}"
    );
    // Global balance: net flux through the boundary vanishes (harmonic function).
    let net: f64 = q
        .iter()
        .enumerate()
        .map(|(i, v)| v * elements_of(&points)[i].length)
        .sum();
    assert!(net.abs() < 5.0, "net flux {net}");
}

#[test]
fn evaluated_interior_potential_follows_the_linear_profile() {
    let (points, bcs) = conduction_problem(10);
    let (u, q) = solve_laplace_bem_2d(&points, &bcs).unwrap_or_else(|e| panic!("{e}"));
    let elements = elements_of(&points);
    for (x, want) in [(0.25, 75.0), (0.5, 50.0), (0.75, 25.0)] {
        let pot = evaluate_potential_2d((x, 0.5), &elements, &u, &q);
        assert!(
            (pot - want).abs() < 6.0,
            "u({x}, 0.5) = {pot}, expected {want}"
        );
    }
}

#[test]
fn cylinder_scenario_reproduces_the_harmonic_function_x() {
    // u = x on the unit circle has the exact interior solution u = x and boundary flux du/dn = x.
    let (u, q) = simulate_2d_cylinder_scenario().unwrap_or_else(|e| panic!("{e}"));
    assert_eq!((u.len(), q.len()), (40, 40));
    for i in 0..40 {
        let a = 2.0 * std::f64::consts::PI * i as f64 / 40.0;
        assert!((u[i] - a.cos()).abs() < 1e-12);
        // The flux is defined on the element between node i and i + 1: its midpoint angle is a + pi/40.
        let mid = a + std::f64::consts::PI / 40.0;
        assert!(
            (q[i] - mid.cos()).abs() < 0.1,
            "q[{i}] = {} vs {}",
            q[i],
            mid.cos()
        );
    }
    // The flux integrates to zero over the closed boundary.
    assert!(q.iter().sum::<f64>().abs() < 1e-6);
}

#[test]
fn mismatched_inputs_and_singular_systems_are_errors() {
    assert!(
        solve_laplace_bem_2d(
            &[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)],
            &[BoundaryCondition::Potential(1.0)]
        )
        .is_err()
    );
    // Pure Neumann data leaves the potential undetermined up to a constant.
    let square = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
    let neumann = [BoundaryCondition::Flux(0.0); 4];
    assert!(solve_laplace_bem_2d(&square, &neumann).is_err());
}

#[test]
fn three_dimensional_bem_self_check_passes() {
    assert_eq!(solve_laplace_bem_3d(), Ok(()));
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_constant_potential_has_zero_flux(u_val in 0.0..1000.0f64) {
        let points = vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        let bcs = vec![BoundaryCondition::Potential(u_val); 4];
        let (u, q) = solve_laplace_bem_2d(&points, &bcs).map_err(TestCaseError::fail)?;
        for i in 0..4 {
            prop_assert!((u[i] - u_val).abs() < 1e-10);
            prop_assert!(q[i].abs() < 1e-5);
        }
    }

    #[test]
    fn prop_linear_boundary_data_is_reproduced_inside(a in -5.0..5.0f64, b in -5.0..5.0f64) {
        // u = a x + b y is harmonic; prescribe it on a circle and test the interior value at the centre offset.
        let n = 48;
        let (mut points, mut bcs) = (Vec::new(), Vec::new());
        for i in 0..n {
            let t = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
            let (x, y) = (t.cos(), t.sin());
            points.push((x, y));
            bcs.push(BoundaryCondition::Potential(a * x + b * y));
        }
        let (u, q) = solve_laplace_bem_2d(&points, &bcs).map_err(TestCaseError::fail)?;
        let elements = elements_of(&points);
        let pot = evaluate_potential_2d((0.2, -0.3), &elements, &u, &q);
        let want = a * 0.2 - b * 0.3;
        prop_assert!((pot - want).abs() < 0.08 * (a.abs() + b.abs() + 1.0), "{pot} vs {want}");
    }
}

mod bem_3d {
    use rssn::sim::physics_bem::{
        BoundaryCondition, SurfaceMesh3D, evaluate_potential_3d, solve_laplace_bem_3d_mesh,
    };

    fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }

    #[test]
    fn generated_meshes_are_valid_closed_surfaces() {
        for sub in 0..3 {
            let m = SurfaceMesh3D::icosphere(1.0, [0.0; 3], sub);
            assert_eq!(m.triangles.len(), 20 * 4usize.pow(sub as u32));
            m.validate().unwrap();
            // Euler characteristic 2: V - E + F with E = 3F/2
            let f = m.triangles.len() as i64;
            assert_eq!(m.vertices.len() as i64 - f / 2, 2);
            // inscribed polyhedron approaches the sphere
            let exact_area = 4.0 * std::f64::consts::PI;
            assert!(m.total_area() < exact_area && m.total_area() > 0.7 * exact_area);
            assert!(m.volume() > 0.0 && m.volume() < exact_area / 3.0);
        }
        let fine = SurfaceMesh3D::icosphere(2.0, [1.0, 0.0, 0.0], 4);
        assert!((fine.volume() - 4.0 / 3.0 * std::f64::consts::PI * 8.0).abs() < 0.2);
        for n in 1..4 {
            let c = SurfaceMesh3D::cube(2.0, [0.5, 0.5, 0.5], n);
            assert_eq!(c.triangles.len(), 12 * n * n);
            c.validate().unwrap();
            assert!((c.volume() - 8.0).abs() < 1e-12);
            assert!((c.total_area() - 24.0).abs() < 1e-12);
            // outward normals: centroid offset from the centre has positive projection
            for t in 0..c.triangles.len() {
                let (cen, nor) = (c.centroid(t).unwrap(), c.normal(t).unwrap());
                let rel = [cen[0] - 0.5, cen[1] - 0.5, cen[2] - 0.5];
                assert!(dot(rel, nor) > 0.0);
            }
        }
    }

    #[test]
    fn validation_rejects_bad_meshes() {
        let mut m = SurfaceMesh3D::icosphere(1.0, [0.0; 3], 0);
        let ok = m.clone();
        m.triangles.pop();
        assert!(m.validate().is_err()); // open surface
        let mut flipped = ok.clone();
        for t in &mut flipped.triangles {
            t.swap(1, 2);
        }
        assert!(flipped.validate().is_err()); // inward
        let mut bad = ok.clone();
        bad.triangles[0][0] = 99;
        assert!(bad.validate().is_err());
        let mut degenerate = ok;
        degenerate.triangles[0] = [0, 0, 1];
        assert!(degenerate.validate().is_err());
    }

    fn solve_dirichlet(mesh: &SurfaceMesh3D, f: impl Fn([f64; 3]) -> f64) -> Vec<f64> {
        let bcs: Vec<_> =
            (0..mesh.triangles.len()).map(|t| BoundaryCondition::Potential(f(mesh.centroid(t).unwrap()))).collect();
        solve_laplace_bem_3d_mesh(mesh, &bcs).unwrap().q
    }

    fn flux_error(mesh: &SurfaceMesh3D, q: &[f64], grad: impl Fn([f64; 3]) -> [f64; 3]) -> f64 {
        // relative RMS error of q against grad(u) . n at the centroids
        let (mut num, mut den) = (0.0, 0.0);
        for (t, &qt) in q.iter().enumerate() {
            let exact = dot(grad(mesh.centroid(t).unwrap()), mesh.normal(t).unwrap());
            num += (qt - exact).powi(2);
            den += exact * exact;
        }
        (num / den).sqrt()
    }

    #[test]
    fn linear_harmonic_function_on_a_sphere_converges() {
        let mut errs = Vec::new();
        for sub in 1..4 {
            let mesh = SurfaceMesh3D::icosphere(1.0, [0.0; 3], sub);
            let q = solve_dirichlet(&mesh, |p| p[0]);
            errs.push(flux_error(&mesh, &q, |_| [1.0, 0.0, 0.0]));
        }
        assert!(errs[2] < 0.02, "{errs:?}");
        assert!(errs[1] < errs[0] && errs[2] < errs[1], "{errs:?}");
        assert!(errs[0] / errs[2] > 3.0, "{errs:?}");
    }

    #[test]
    fn point_source_outside_the_sphere() {
        let x0 = [2.0, 0.3, -0.4];
        let u = |p: [f64; 3]| {
            let d = [p[0] - x0[0], p[1] - x0[1], p[2] - x0[2]];
            1.0 / dot(d, d).sqrt()
        };
        let grad = |p: [f64; 3]| {
            let d = [p[0] - x0[0], p[1] - x0[1], p[2] - x0[2]];
            let r3 = dot(d, d).powf(1.5);
            [-d[0] / r3, -d[1] / r3, -d[2] / r3]
        };
        let mesh = SurfaceMesh3D::icosphere(1.0, [0.0; 3], 3);
        let bcs: Vec<_> =
            (0..mesh.triangles.len()).map(|t| BoundaryCondition::Potential(u(mesh.centroid(t).unwrap()))).collect();
        let sol = solve_laplace_bem_3d_mesh(&mesh, &bcs).unwrap();
        assert!(flux_error(&mesh, &sol.q, grad) < 0.03);
        // interior evaluation against the exact harmonic function
        for p in [[0.0, 0.0, 0.0], [0.3, -0.2, 0.1], [-0.4, 0.4, 0.3], [0.5, 0.0, 0.0]] {
            let v = evaluate_potential_3d(p, &mesh, &sol).unwrap();
            assert!((v - u(p)).abs() < 0.01 * u(p).abs() + 2e-3, "u({p:?}) = {v}, exact {}", u(p));
        }
    }

    #[test]
    fn constant_potential_has_zero_flux_and_center_value() {
        let mesh = SurfaceMesh3D::icosphere(1.5, [0.2, 0.0, 0.0], 2);
        let q = solve_dirichlet(&mesh, |_| 3.0);
        assert!(q.iter().all(|v| v.abs() < 1e-9), "{q:?}");
    }

    #[test]
    fn mixed_dirichlet_neumann_data_recovers_the_harmonic_function() {
        // u = x z: Dirichlet on the upper hemisphere, Neumann below.
        let u = |p: [f64; 3]| p[0] * p[2];
        let grad = |p: [f64; 3]| [p[2], 0.0, p[0]];
        let mesh = SurfaceMesh3D::icosphere(1.0, [0.0; 3], 3);
        let bcs: Vec<_> = (0..mesh.triangles.len())
            .map(|t| {
                let c = mesh.centroid(t).unwrap();
                if c[1] >= 0.0 {
                    BoundaryCondition::Potential(u(c))
                } else {
                    BoundaryCondition::Flux(dot(grad(c), mesh.normal(t).unwrap()))
                }
            })
            .collect();
        let sol = solve_laplace_bem_3d_mesh(&mesh, &bcs).unwrap();
        let (mut eu, mut nu) = (0.0, 0.0);
        let (mut eq, mut nq) = (0.0, 0.0);
        for t in 0..mesh.triangles.len() {
            let c = mesh.centroid(t).unwrap();
            let qx = dot(grad(c), mesh.normal(t).unwrap());
            eu += (sol.u[t] - u(c)).powi(2);
            nu += u(c).powi(2);
            eq += (sol.q[t] - qx).powi(2);
            nq += qx * qx;
        }
        assert!((eu / nu).sqrt() < 0.05, "u error {}", (eu / nu).sqrt());
        assert!((eq / nq).sqrt() < 0.05, "q error {}", (eq / nq).sqrt());
    }

    #[test]
    fn pure_neumann_problem_is_solved_up_to_a_constant() {
        let mesh = SurfaceMesh3D::icosphere(1.0, [0.0; 3], 3);
        let bcs: Vec<_> =
            (0..mesh.triangles.len()).map(|t| BoundaryCondition::Flux(mesh.normal(t).unwrap()[0])).collect();
        let sol = solve_laplace_bem_3d_mesh(&mesh, &bcs).unwrap();
        // u = x + c with area-weighted mean zero: x itself up to symmetry error
        let (mut num, mut den) = (0.0, 0.0);
        for t in 0..mesh.triangles.len() {
            num += (sol.u[t] - mesh.centroid(t).unwrap()[0]).powi(2);
            den += mesh.centroid(t).unwrap()[0].powi(2);
        }
        assert!((num / den).sqrt() < 0.05, "{}", (num / den).sqrt());
        // incompatible data is rejected
        let bad = vec![BoundaryCondition::Flux(1.0); mesh.triangles.len()];
        assert!(solve_laplace_bem_3d_mesh(&mesh, &bad).is_err());
    }

    #[test]
    fn cube_with_a_linear_potential() {
        let mesh = SurfaceMesh3D::cube(2.0, [0.0; 3], 4);
        let q = solve_dirichlet(&mesh, |p| p[0]);
        // faces away from the edges: flux +-1 on x faces, 0 on y/z faces
        let mut checked = 0;
        for (t, &qt) in q.iter().enumerate() {
            let (c, n) = (mesh.centroid(t).unwrap(), mesh.normal(t).unwrap());
            if c.iter().all(|v| v.abs() < 0.6 || v.abs() > 0.99) && c.iter().filter(|v| v.abs() > 0.99).count() == 1 {
                assert!((qt - n[0]).abs() < 0.12, "face centre flux {qt} vs {}", n[0]);
                checked += 1;
            }
        }
        assert!(checked > 10);
        let sol = solve_laplace_bem_3d_mesh(
            &mesh,
            &(0..mesh.triangles.len())
                .map(|t| BoundaryCondition::Potential(mesh.centroid(t).unwrap()[0]))
                .collect::<Vec<_>>(),
        )
        .unwrap();
        let v = evaluate_potential_3d([0.3, 0.1, -0.2], &mesh, &sol).unwrap();
        assert!((v - 0.3).abs() < 0.03, "{v}");
    }

    #[test]
    fn input_errors_are_reported() {
        let mesh = SurfaceMesh3D::icosphere(1.0, [0.0; 3], 0);
        assert!(solve_laplace_bem_3d_mesh(&mesh, &[]).is_err());
        let nan = vec![BoundaryCondition::Potential(f64::NAN); mesh.triangles.len()];
        assert!(solve_laplace_bem_3d_mesh(&mesh, &nan).is_err());
        let sol = solve_laplace_bem_3d_mesh(&mesh, &vec![BoundaryCondition::Potential(1.0); 20]).unwrap();
        let other = SurfaceMesh3D::icosphere(1.0, [0.0; 3], 1);
        assert!(evaluate_potential_3d([0.0; 3], &other, &sol).is_err());
    }
}
