//! Coverage round 2: the triangular finite-element module (mesh construction, boundary
//! detection, refinement, solves with manufactured solutions, error estimation, adaptivity).

use std::collections::HashMap;
use std::f64::consts::PI;

use rssn::sim::fem_triangles::Degree;
use rssn::sim::fem_triangles::Problem;
use rssn::sim::fem_triangles::TriMesh;
use rssn::sim::fem_triangles::error_indicators;
use rssn::sim::fem_triangles::solve;
use rssn::sim::fem_triangles::solve_adaptive;

fn signed_area(
    m: &TriMesh,
    t: &[usize; 3],
) -> f64 {
    let [a, b, c] = t.map(|v| m.nodes[v]);
    0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]))
}

fn total_area(m: &TriMesh) -> f64 {
    m.triangles.iter().map(|t| signed_area(m, t).abs()).sum()
}

/// Every undirected edge is shared by one or two triangles, and no node lies in the interior of a
/// foreign edge (no hanging nodes): checked through the Euler characteristic of a disc.
fn assert_conforming_disc(m: &TriMesh) {
    let mut edges: HashMap<(usize, usize), usize> = HashMap::new();
    for t in &m.triangles {
        for i in 0..3 {
            let (a, b) = (t[i], t[(i + 1) % 3]);
            *edges.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    assert!(
        edges.values().all(|&c| c == 1 || c == 2),
        "non-manifold edge"
    );
    let used: std::collections::HashSet<usize> = m.triangles.iter().flatten().copied().collect();
    let chi = used.len() as i64 - edges.len() as i64 + m.triangles.len() as i64;
    assert_eq!(
        chi, 1,
        "V - E + F must be 1 for a simply connected conforming mesh"
    );
}

#[test]
fn rectangle_mesh_has_the_expected_counts_orientation_and_boundary() {
    let (nx, ny) = (5, 3);
    let m = TriMesh::rectangle(-1.0, 2.0, 0.5, 1.5, nx, ny);
    assert_eq!(m.nodes.len(), (nx + 1) * (ny + 1));
    assert_eq!(m.triangles.len(), 2 * nx * ny);
    assert!((total_area(&m) - 3.0).abs() < 1e-13);
    for t in &m.triangles {
        assert!(signed_area(&m, t) > 0.0, "counter-clockwise triangles");
    }
    let on = m.boundary_vertices();
    assert_eq!(on.len(), m.nodes.len());
    assert_eq!(on.iter().filter(|&&b| b).count(), 2 * (nx + ny));
    for (p, &b) in m.nodes.iter().zip(&on) {
        let at_wall = (p[0] + 1.0).abs() < 1e-12
            || (p[0] - 2.0).abs() < 1e-12
            || (p[1] - 0.5).abs() < 1e-12
            || (p[1] - 1.5).abs() < 1e-12;
        assert_eq!(b, at_wall, "node {p:?}");
    }
    assert_conforming_disc(&m);
    // an empty mesh has no boundary
    assert!(TriMesh::default().boundary_vertices().is_empty());
}

#[test]
fn l_shape_has_three_quarters_of_the_area_and_a_reentrant_corner() {
    for n in [1usize, 2, 4] {
        let m = TriMesh::l_shape(n);
        assert!((total_area(&m) - 3.0).abs() < 1e-12, "n = {n}");
        assert_eq!(m.triangles.len(), 6 * n * n, "n = {n}");
        // every stored node is used: the mesh is compacted
        let used: std::collections::HashSet<usize> =
            m.triangles.iter().flatten().copied().collect();
        assert_eq!(used.len(), m.nodes.len());
        // 8 unit-length sides, each cut into n segments
        let on = m.boundary_vertices();
        assert_eq!(on.iter().filter(|&&b| b).count(), 8 * n, "n = {n}");
        // the forbidden quadrant (x > 0, y < 0) has no node in its open interior
        assert!(m.nodes.iter().all(|p| !(p[0] > 1e-12 && p[1] < -1e-12)));
        // the reentrant corner is a boundary node
        let o = m
            .nodes
            .iter()
            .position(|p| p[0].abs() < 1e-12 && p[1].abs() < 1e-12)
            .unwrap();
        assert!(on[o]);
        assert_conforming_disc(&m);
    }
}

#[test]
fn compaction_renumbers_without_changing_the_geometry() {
    let m = TriMesh {
        nodes: vec![
            [9.0, 9.0],
            [0.0, 0.0],
            [7.0, 7.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [5.0, 5.0],
        ],
        triangles: vec![[1, 3, 4]],
    };
    let c = m.clone().compacted();
    assert_eq!(c.nodes.len(), 3);
    assert_eq!(c.triangles, vec![[0, 1, 2]]);
    assert!((total_area(&c) - 0.5).abs() < 1e-15);
    assert_eq!(c.nodes[0], [0.0, 0.0]);
    // compacting twice changes nothing
    let cc = c.clone().compacted();
    assert_eq!(cc.nodes, c.nodes);
    assert_eq!(cc.triangles, c.triangles);
}

#[test]
fn refinement_is_conforming_conserves_area_and_keeps_shape_regularity() {
    let mut m = TriMesh::rectangle(0.0, 1.0, 0.0, 1.0, 2, 2);
    // nothing marked -> nothing changes
    let same = m.refine(&vec![false; m.triangles.len()]);
    assert_eq!(same.triangles.len(), m.triangles.len());
    assert_eq!(same.nodes.len(), m.nodes.len());
    let mut count = m.triangles.len();
    for round in 0..4 {
        // mark only the triangles touching the corner, producing a graded mesh
        let marked: Vec<bool> = m
            .triangles
            .iter()
            .map(|t| {
                t.iter()
                    .any(|&v| m.nodes[v][0] < 0.3 && m.nodes[v][1] < 0.3)
            })
            .collect();
        m = m.refine(&marked);
        assert!(m.triangles.len() > count, "round {round}");
        count = m.triangles.len();
        assert!((total_area(&m) - 1.0).abs() < 1e-12);
        assert_conforming_disc(&m);
        for t in &m.triangles {
            assert!(
                signed_area(&m, t) > 0.0,
                "orientation preserved by bisection"
            );
        }
        // newest-vertex bisection keeps the angles bounded away from zero
        let mut min_angle = f64::MAX;
        for t in &m.triangles {
            let p = t.map(|v| m.nodes[v]);
            for i in 0..3 {
                let (a, b, c) = (p[i], p[(i + 1) % 3], p[(i + 2) % 3]);
                let u = [b[0] - a[0], b[1] - a[1]];
                let v = [c[0] - a[0], c[1] - a[1]];
                let ang = (u[0] * v[1] - u[1] * v[0])
                    .abs()
                    .atan2(u[0] * v[0] + u[1] * v[1]);
                min_angle = min_angle.min(ang);
            }
        }
        assert!(min_angle > 0.1, "min angle {min_angle} after round {round}");
    }
    // boundary stays exactly the square's boundary
    let on = m.boundary_vertices();
    for (p, &b) in m.nodes.iter().zip(&on) {
        let wall = p[0].abs() < 1e-12
            || (p[0] - 1.0).abs() < 1e-12
            || p[1].abs() < 1e-12
            || (p[1] - 1.0).abs() < 1e-12;
        assert_eq!(b, wall, "{p:?}");
    }
}

fn poisson_l2(
    n: usize,
    degree: Degree,
) -> (f64, f64) {
    let exact = |x: f64, y: f64| (PI * x).sin() * (PI * y).sin();
    let source = |x: f64, y: f64| 2.0 * PI * PI * exact(x, y);
    let (one, zero) = (|_: f64, _: f64| 1.0, |_: f64, _: f64| 0.0);
    let problem = Problem {
        k: &one,
        c: &zero,
        f: &source,
        g: &zero,
    };
    let mesh = TriMesh::rectangle(0.0, 1.0, 0.0, 1.0, n, n);
    let sol = solve(&mesh, degree, &problem).unwrap();
    (sol.l2_error(&mesh, &exact), sol.max_nodal_error(&exact))
}

#[test]
fn p1_and_p2_converge_at_their_optimal_l2_orders() {
    let (e4, _) = poisson_l2(4, Degree::P1);
    let (e8, _) = poisson_l2(8, Degree::P1);
    let (e16, n16) = poisson_l2(16, Degree::P1);
    let r1 = (e4 / e8).log2();
    let r2 = (e8 / e16).log2();
    assert!(r1 > 1.8 && r2 > 1.9, "P1 rates {r1} {r2}");
    assert!(n16 < 1e-2);
    let (p4, _) = poisson_l2(4, Degree::P2);
    let (p8, _) = poisson_l2(8, Degree::P2);
    let rp = (p4 / p8).log2();
    assert!(rp > 2.7, "P2 rate {rp}");
    assert!(p8 < e8, "P2 is more accurate than P1 on the same mesh");
    // P2 carries edge degrees of freedom
    let mesh = TriMesh::rectangle(0.0, 1.0, 0.0, 1.0, 3, 3);
    let (one, zero) = (|_: f64, _: f64| 1.0, |_: f64, _: f64| 0.0);
    let problem = Problem {
        k: &one,
        c: &zero,
        f: &one,
        g: &zero,
    };
    let s1 = solve(&mesh, Degree::P1, &problem).unwrap();
    let s2 = solve(&mesh, Degree::P2, &problem).unwrap();
    assert_eq!(s1.values.len(), mesh.nodes.len());
    assert!(s2.values.len() > s1.values.len());
    assert_eq!(s2.points.len(), s2.values.len());
}

#[test]
fn solution_respects_dirichlet_data_the_maximum_principle_and_linearity() {
    let mesh = TriMesh::rectangle(0.0, 1.0, 0.0, 1.0, 8, 8);
    let (one, zero) = (|_: f64, _: f64| 1.0, |_: f64, _: f64| 0.0);
    let g = |x: f64, y: f64| 1.0 + x + 2.0 * y;
    // harmonic data (affine) with f = 0 is reproduced exactly by P1
    let problem = Problem {
        k: &one,
        c: &zero,
        f: &zero,
        g: &g,
    };
    let sol = solve(&mesh, Degree::P1, &problem).unwrap();
    assert!(
        sol.max_nodal_error(&g) < 1e-10,
        "affine solutions are in the P1 space"
    );
    // maximum principle: f >= 0 and g = 0 -> 0 <= u
    let f = |x: f64, y: f64| 1.0 + x * y;
    let problem = Problem {
        k: &one,
        c: &zero,
        f: &f,
        g: &zero,
    };
    let s = solve(&mesh, Degree::P1, &problem).unwrap();
    assert!(s.values.iter().all(|&v| v > -1e-12));
    // linearity in the source: u(2 f) = 2 u(f)
    let f2 = |x: f64, y: f64| 2.0 * f(x, y);
    let problem2 = Problem {
        k: &one,
        c: &zero,
        f: &f2,
        g: &zero,
    };
    let s2 = solve(&mesh, Degree::P1, &problem2).unwrap();
    for (a, b) in s.values.iter().zip(&s2.values) {
        assert!((2.0 * a - b).abs() < 1e-9 * (1.0 + b.abs()));
    }
    // scaling the conductivity by 4 divides the solution by 4
    let k4 = |_: f64, _: f64| 4.0;
    let problem4 = Problem {
        k: &k4,
        c: &zero,
        f: &f,
        g: &zero,
    };
    let s4 = solve(&mesh, Degree::P1, &problem4).unwrap();
    for (a, b) in s.values.iter().zip(&s4.values) {
        assert!((a - 4.0 * b).abs() < 1e-9 * (1.0 + a.abs()));
    }
}

#[test]
fn degenerate_triangles_are_reported_as_errors() {
    let mesh = TriMesh {
        nodes: vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [0.0, 1.0]],
        triangles: vec![[0, 1, 2], [0, 1, 3]],
    };
    let (one, zero) = (|_: f64, _: f64| 1.0, |_: f64, _: f64| 0.0);
    let problem = Problem {
        k: &one,
        c: &zero,
        f: &one,
        g: &zero,
    };
    assert!(solve(&mesh, Degree::P1, &problem).is_err());
    assert!(solve(&mesh, Degree::P2, &problem).is_err());
}

#[test]
fn error_estimator_is_zero_for_exactly_represented_solutions_and_adaptivity_converges() {
    let (one, zero) = (|_: f64, _: f64| 1.0, |_: f64, _: f64| 0.0);
    let g = |x: f64, y: f64| 3.0 * x - y;
    let mesh = TriMesh::rectangle(0.0, 1.0, 0.0, 1.0, 4, 4);
    let problem = Problem {
        k: &one,
        c: &zero,
        f: &zero,
        g: &g,
    };
    let sol = solve(&mesh, Degree::P1, &problem).unwrap();
    let eta = error_indicators(&mesh, &sol, &problem);
    assert_eq!(eta.len(), mesh.triangles.len());
    assert!(
        eta.iter().all(|&e| e < 1e-9),
        "affine solution has no residual or jump"
    );

    // Laplace on the L-shape with the corner singularity r^{2/3} sin(2 theta / 3)
    let u = |x: f64, y: f64| {
        let r = x.hypot(y);
        let th = y.atan2(x);
        let th = if th < 0.0 {
            th + 2.0 * PI
        } else {
            th
        };
        r.powf(2.0 / 3.0) * (2.0 * th / 3.0).sin()
    };
    let problem = Problem {
        k: &one,
        c: &zero,
        f: &zero,
        g: &u,
    };
    let start = TriMesh::l_shape(2);
    let (final_mesh, sol, history) = solve_adaptive(&start, &problem, 5, 0.5).unwrap();
    assert_eq!(history.len(), 6);
    assert!(history[5].0 > history[0].0, "degrees of freedom grow");
    assert!(
        history[5].1 < 0.6 * history[0].1,
        "estimator decreases: {history:?}"
    );
    assert!(final_mesh.triangles.len() > start.triangles.len());
    assert!((total_area(&final_mesh) - 3.0).abs() < 1e-10);
    // the true error also decreases and refinement concentrates at the corner
    let e0 = {
        let s = solve(&start, Degree::P1, &problem).unwrap();
        s.l2_error(&start, &u)
    };
    let e1 = sol.l2_error(&final_mesh, &u);
    assert!(e1 < 0.5 * e0, "adaptive L2 error {e1} vs initial {e0}");
    let near = |m: &TriMesh| {
        m.triangles
            .iter()
            .filter(|t| t.iter().any(|&v| m.nodes[v][0].hypot(m.nodes[v][1]) < 0.3))
            .count() as f64
            / m.triangles.len() as f64
    };
    assert!(
        near(&final_mesh) > near(&start),
        "refinement clusters at the reentrant corner"
    );
}
