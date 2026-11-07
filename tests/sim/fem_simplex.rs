//! P1/P2 simplex finite elements in 1D, 2D and 3D.

use rssn::sim::boundary::BcValue;
use rssn::sim::boundary::BoundaryCondition as Bc;
use rssn::sim::boundary::MarkerBcs;
use rssn::sim::fem_simplex::FemProblem;
use rssn::sim::fem_simplex::FemSpace;
use rssn::sim::fem_simplex::SimplexMesh;
use rssn::sim::fem_simplex::error_indicators;
use rssn::sim::fem_simplex::marker;
use rssn::sim::fem_simplex::simplex_rule;
use rssn::sim::fem_simplex::solve_adaptive;
use rssn::sim::fem_simplex::solve_heat;
use rssn::sim::fem_simplex::solve_static;
use rssn::sim::structured::LinearSolver;
use std::collections::HashMap;
use std::f64::consts::PI;

#[test]
fn grundmann_moeller_rule_integrates_monomials_exactly() {
    let fact = |m: usize| (1..=m).map(|x| x as f64).product::<f64>();
    for n in 0..=3usize {
        for s in 0..=3usize {
            let (pts, w) = simplex_rule(n, s);
            assert!((w.iter().sum::<f64>() - 1.0).abs() < 1e-12);
            // all monomials of total degree <= 2s+1
            for a0 in 0..=(2 * s + 1) {
                for a1 in 0..=(2 * s + 1 - a0) {
                    for a2 in 0..=(2 * s + 1 - a0 - a1) {
                        let a = [a0, a1, a2];
                        if n < 2 && a2 > 0 || n < 1 && a1 > 0 {
                            continue;
                        }
                        let deg: usize = a.iter().take(n + 1).sum();
                        // monomial over barycentric coordinates 0..=n
                        let q: f64 = pts
                            .iter()
                            .zip(&w)
                            .map(|(p, wq)| {
                                wq * p
                                    .iter()
                                    .zip(a.iter())
                                    .map(|(l, &e)| l.powi(e as i32))
                                    .product::<f64>()
                            })
                            .sum();
                        let exact = fact(n)
                            * a.iter().take(n + 1).map(|&e| fact(e)).product::<f64>()
                            / fact(n + deg);
                        assert!(
                            (q - exact).abs() < 1e-12,
                            "n {n} s {s} a {a:?}: {q} vs {exact}"
                        );
                    }
                }
            }
        }
    }
}

fn exact(
    dim: usize,
    x: [f64; 3],
) -> f64 {
    (0..dim).map(|a| x[a]).sum::<f64>().exp()
}

fn bcs_for(dim: usize) -> MarkerBcs<Bc> {
    let u = move |x: [f64; 3]| exact(dim, x);
    let dir = Bc::Dirichlet(BcValue::space(u));
    let mut b = MarkerBcs::new();
    if dim == 1 {
        return b.set(marker::XMIN, dir).set(
            marker::XMAX,
            Bc::Robin {
                alpha: 2.0,
                beta: 1.0,
                g: BcValue::space(move |x| 3.0 * u(x)),
            },
        );
    }
    b = b
        .set(marker::XMIN, dir.clone())
        .set(marker::XMAX, Bc::Neumann(BcValue::space(u)))
        .set(
            marker::YMIN,
            Bc::Robin {
                alpha: 2.0,
                beta: 1.0,
                g: BcValue::space(u),
            },
        )
        .set(marker::YMAX, dir);
    if dim == 3 {
        b = b
            .set(marker::ZMIN, Bc::Neumann(BcValue::space(move |x| -u(x))))
            .set(marker::ZMAX, Bc::Neumann(BcValue::space(u)));
    }
    b
}

fn mesh(
    dim: usize,
    n: usize,
) -> SimplexMesh {
    match dim {
        | 1 => SimplexMesh::interval(n, 0.0, 1.0),
        | 2 => SimplexMesh::rectangle(n, n, [0.0; 2], [1.0; 2]),
        | _ => SimplexMesh::boxed([n; 3], [0.0; 3], [1.0; 3]),
    }
}

fn poisson_l2(
    dim: usize,
    degree: usize,
    n: usize,
) -> f64 {
    let space = FemSpace::new(mesh(dim, n), degree).unwrap();
    let bcs = bcs_for(dim);
    // -Δu + c u = f with u = exp(sum x): f = (c - dim) u
    let c = 1.5;
    let f = move |_t: f64, x: [f64; 3]| (c - dim as f64) * exact(dim, x);
    let prob = FemProblem {
        k: 1.0,
        c,
        f: &f,
        bcs: &bcs,
    };
    let u = solve_static(&space, &prob, 0.0, LinearSolver::Auto).unwrap();
    space.l2_error(&u, |x| exact(dim, x))
}

#[test]
fn mixed_marker_bcs_converge_with_optimal_order_p1_and_p2_in_1d_2d_3d() {
    let cases: [(usize, usize, [usize; 2]); 6] = [
        (1, 1, [8, 16]),
        (1, 2, [4, 8]),
        (2, 1, [8, 16]),
        (2, 2, [4, 8]),
        (3, 1, [4, 8]),
        (3, 2, [2, 4]),
    ];
    for (dim, deg, ns) in cases {
        let e1 = poisson_l2(dim, deg, ns[0]);
        let e2 = poisson_l2(dim, deg, ns[1]);
        let order = (e1 / e2).log2();
        let min = deg as f64 + 0.7;
        assert!(order > min, "dim {dim} P{deg}: order {order} ({e1}, {e2})");
    }
}

#[test]
fn bad_inputs_are_rejected() {
    assert!(FemSpace::new(SimplexMesh::interval(4, 0.0, 1.0), 3).is_err());
    let space = FemSpace::new(SimplexMesh::interval(4, 0.0, 1.0), 1).unwrap();
    let bcs = MarkerBcs::with_default(Bc::Periodic);
    let f = |_t: f64, _x: [f64; 3]| 0.0;
    let prob = FemProblem {
        k: 1.0,
        c: 1.0,
        f: &f,
        bcs: &bcs,
    };
    assert!(solve_static(&space, &prob, 0.0, LinearSolver::Auto).is_err());
    let bcs = MarkerBcs::with_default(Bc::Robin {
        alpha: 1.0,
        beta: 0.0,
        g: BcValue::constant(0.0),
    });
    let prob = FemProblem {
        k: 1.0,
        c: 1.0,
        f: &f,
        bcs: &bcs,
    };
    assert!(solve_static(&space, &prob, 0.0, LinearSolver::Auto).is_err());
}

#[test]
fn heat_equation_with_time_dependent_dirichlet_data_2d_p2_and_3d_p1() {
    // u = exp(-t) sin(pi x) + t y : u_t - Δu = -exp(-t) sin(pi x) + y + pi^2 exp(-t) sin(pi x) ... computed below
    let u = |t: f64, x: [f64; 3]| (-t).exp() * (PI * x[0]).sin() * (PI * x[1]).sin() + t * x[0];
    let f = move |t: f64, x: [f64; 3]| {
        // u_t = -e sin sin + x ; -Δu = 2 pi^2 e sin sin
        let e = (-t).exp() * (PI * x[0]).sin() * (PI * x[1]).sin();
        -e + x[0] + 2.0 * PI * PI * e
    };
    let bcs = MarkerBcs::with_default(Bc::Dirichlet(BcValue::spacetime(u)));
    let run = |n: usize, steps: usize| {
        let space = FemSpace::new(SimplexMesh::rectangle(n, n, [0.0; 2], [1.0; 2]), 2).unwrap();
        let prob = FemProblem {
            k: 1.0,
            c: 0.0,
            f: &f,
            bcs: &bcs,
        };
        let v = solve_heat(
            &space,
            &prob,
            |x| u(0.0, x),
            0.2,
            steps,
            0.5,
            LinearSolver::Auto,
        )
        .unwrap();
        space.l2_error(&v, |x| u(0.2, x))
    };
    let (e1, e2) = (run(4, 8), run(8, 16));
    assert!((e1 / e2).log2() > 1.8 && e2 < 1e-3, "{e1} {e2}");

    // 3D: u = exp(-t) cos(pi x) with Neumann sides (zero flux in x, linear in y,z handled by symmetry)
    let u3 = |t: f64, x: [f64; 3]| (-PI * PI * t).exp() * (PI * x[0]).cos();
    let bcs3 = MarkerBcs::with_default(Bc::neumann(0.0));
    let f3 = |_t: f64, _x: [f64; 3]| 0.0;
    let space = FemSpace::new(SimplexMesh::boxed([16, 2, 2], [0.0; 3], [1.0; 3]), 1).unwrap();
    let prob = FemProblem {
        k: 1.0,
        c: 0.0,
        f: &f3,
        bcs: &bcs3,
    };
    let v = solve_heat(
        &space,
        &prob,
        |x| u3(0.0, x),
        0.05,
        20,
        0.5,
        LinearSolver::Auto,
    )
    .unwrap();
    assert!(space.l2_error(&v, |x| u3(0.05, x)) < 6e-3);
    assert!(space.eval(&v, [0.5, 0.5, 0.5]).is_some());
}

fn centroid(
    mesh: &SimplexMesh,
    cell: &[usize],
) -> [f64; 3] {
    let mut x = [0.0; 3];
    for &v in cell {
        for a in 0..3 {
            x[a] += mesh.nodes[v][a] / cell.len() as f64;
        }
    }
    x
}

fn face_vertices(
    mesh: &SimplexMesh,
    cell: usize,
    face: usize,
) -> Vec<usize> {
    let mut v: Vec<usize> = mesh.cells[cell]
        .iter()
        .enumerate()
        .filter(|(j, _)| *j != face)
        .map(|(_, &q)| q)
        .collect();
    v.sort_unstable();
    v
}

fn grad_sine(x: [f64; 3]) -> [f64; 3] {
    [
        PI * (PI * x[0]).cos() * (PI * x[1]).sin() * (PI * x[2]).sin(),
        PI * (PI * x[0]).sin() * (PI * x[1]).cos() * (PI * x[2]).sin(),
        PI * (PI * x[0]).sin() * (PI * x[1]).sin() * (PI * x[2]).cos(),
    ]
}

/// P1 `H¹` seminorm error against an exact gradient.
fn h1_error_p1(
    space: &FemSpace,
    u: &[f64],
    grad_exact: impl Fn([f64; 3]) -> [f64; 3],
) -> f64 {
    let dim = space.mesh.dim;
    let (pts, w) = simplex_rule(dim, 3);
    let mut s = 0.0;
    for ci in 0..space.mesh.cells.len() {
        let vol = space.mesh.measure(ci);
        let g = space.mesh.bary_gradients(ci);
        let cd = &space.cell_dofs[ci];
        let mut grad = [0.0; 3];
        for m in 0..=dim {
            for a in 0..3 {
                grad[a] += u[cd[m]] * g[m][a];
            }
        }
        for (lam, wq) in pts.iter().zip(&w) {
            let mut x = [0.0; 3];
            for (m, &l) in lam.iter().enumerate() {
                for a in 0..3 {
                    x[a] += l * space.mesh.nodes[space.mesh.cells[ci][m]][a];
                }
            }
            let ge = grad_exact(x);
            s += vol
                * wq
                * ((grad[0] - ge[0]).powi(2)
                    + (grad[1] - ge[1]).powi(2)
                    + (grad[2] - ge[2]).powi(2));
        }
    }
    s.sqrt()
}

fn face_centroid(
    mesh: &SimplexMesh,
    cell: usize,
    face: usize,
) -> [f64; 3] {
    let key = face_vertices(mesh, cell, face);
    let mut c = [0.0; 3];
    for &v in &key {
        for a in 0..3 {
            c[a] += mesh.nodes[v][a] / key.len() as f64;
        }
    }
    c
}

/// True when `x` lies in the closure of the local face `face` of `cell`.
fn on_face(
    mesh: &SimplexMesh,
    cell: usize,
    face: usize,
    x: [f64; 3],
) -> bool {
    let g = mesh.bary_gradients(cell);
    let x0 = mesh.nodes[mesh.cells[cell][0]];
    (0..=mesh.dim).all(|m| {
        let l: f64 =
            (0..3).map(|a| g[m][a] * (x[a] - x0[a])).sum::<f64>() + if m == 0 { 1.0 } else { 0.0 };
        if m == face {
            l.abs() < 1e-10
        } else {
            l >= -1e-10
        }
    })
}

/// Signed volume (positive for correctly oriented tetrahedra).
fn signed_volume(
    mesh: &SimplexMesh,
    cell: usize,
) -> f64 {
    let c = &mesh.cells[cell];
    let p = |i: usize| mesh.nodes[c[i]];
    let e1 = [p(1)[0] - p(0)[0], p(1)[1] - p(0)[1], p(1)[2] - p(0)[2]];
    let e2 = [p(2)[0] - p(0)[0], p(2)[1] - p(0)[1], p(2)[2] - p(0)[2]];
    let e3 = [p(3)[0] - p(0)[0], p(3)[1] - p(0)[1], p(3)[2] - p(0)[2]];
    let cr = [
        e1[1] * e2[2] - e1[2] * e2[1],
        e1[2] * e2[0] - e1[0] * e2[2],
        e1[0] * e2[1] - e1[1] * e2[0],
    ];
    (cr[0] * e3[0] + cr[1] * e3[1] + cr[2] * e3[2]) / 6.0
}

/// Volume, positivity, conformity and marker inheritance of a refined mesh.
///
/// `coarse` is the mesh the refinement started from: every face used by a
/// single cell must lie on its boundary, otherwise it is a hanging node.
fn assert_conforming(
    mesh: &SimplexMesh,
    coarse: &SimplexMesh,
    expected_volume: f64,
) -> usize {
    assert_eq!(mesh.level.len(), mesh.cells.len());
    let volume: f64 = (0..mesh.cells.len()).map(|c| mesh.measure(c)).sum();
    assert!(
        (volume - expected_volume).abs() < 1e-12,
        "volume {volume} vs {expected_volume}"
    );
    for c in 0..mesh.cells.len() {
        assert!(
            signed_volume(mesh, c) > 0.0,
            "cell {c} is not positively oriented"
        );
    }
    let mut faces: HashMap<Vec<usize>, usize> = HashMap::new();
    for ci in 0..mesh.cells.len() {
        for f in 0..=mesh.dim {
            *faces.entry(face_vertices(mesh, ci, f)).or_insert(0) += 1;
        }
    }
    let mut interior = 0;
    for (key, &n) in &faces {
        assert!(n <= 2, "face {key:?} is shared by {n} cells");
        if n == 2 {
            interior += 1;
        } else {
            let mut cen = [0.0; 3];
            for &v in key {
                for a in 0..3 {
                    cen[a] += mesh.nodes[v][a] / key.len() as f64;
                }
            }
            assert!(
                coarse
                    .facets
                    .iter()
                    .any(|&(ci, f, _)| on_face(coarse, ci, f, cen)),
                "face {key:?} at {cen:?} is used once and is not on the coarse boundary (hanging node)"
            );
        }
    }
    assert_eq!(interior + mesh.facets.len(), faces.len());
    for &(ci, f, mk) in &mesh.facets {
        let cen = face_centroid(mesh, ci, f);
        assert!(
            coarse
                .facets
                .iter()
                .any(|&(cci, cf, cmk)| cmk == mk && on_face(coarse, cci, cf, cen)),
            "refined facet at {cen:?} does not inherit marker {mk}"
        );
    }
    interior
}

/// A fresh mesh with the given cells, compacted, all boundary facets `marker`.
fn rebuild(
    coarse: &SimplexMesh,
    cells: &[Vec<usize>],
    marker: i32,
) -> SimplexMesh {
    let mut nodes: Vec<[f64; 3]> = Vec::new();
    let mut map = vec![usize::MAX; coarse.nodes.len()];
    let mut remapped = Vec::with_capacity(cells.len());
    for c in cells {
        let mut nc = Vec::with_capacity(c.len());
        for &v in c {
            if map[v] == usize::MAX {
                map[v] = nodes.len();
                nodes.push(coarse.nodes[v]);
            }
            nc.push(map[v]);
        }
        remapped.push(nc);
    }
    let mut faces: HashMap<Vec<usize>, Vec<(usize, usize)>> = HashMap::new();
    for (ci, c) in remapped.iter().enumerate() {
        for f in 0..=coarse.dim {
            let mut key: Vec<usize> = c
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != f)
                .map(|(_, &v)| v)
                .collect();
            key.sort_unstable();
            faces.entry(key).or_default().push((ci, f));
        }
    }
    let mut facets = Vec::new();
    for list in faces.values() {
        if let [(ci, f)] = list.as_slice() {
            facets.push((*ci, *f, marker));
        }
    }
    facets.sort_unstable();
    let level = vec![0; remapped.len()];
    SimplexMesh {
        dim: coarse.dim,
        nodes,
        cells: remapped,
        facets,
        level,
    }
}

/// L-shaped prism `([0, 2]² \ [1, 2]²) × [0, 1]`: its re-entrant edge
/// `x = y = 1`, `z ∈ [0, 1]` carries a `r^(2/3)` singularity.
fn l_prism(n: usize) -> SimplexMesh {
    let base = SimplexMesh::boxed([n, n, 1], [0.0; 3], [2.0, 2.0, 1.0]);
    let cells: Vec<Vec<usize>> = base
        .cells
        .iter()
        .filter(|c| {
            let x = centroid(&base, c.as_slice());
            !(x[0] > 1.0 + 1e-12 && x[1] > 1.0 + 1e-12)
        })
        .cloned()
        .collect();
    rebuild(&base, &cells, 0)
}

/// `u = r^(2/3) sin(2φ/3)` in the polar coordinates of the re-entrant edge
/// (`φ = 0` and `φ = 3π/2` are the two re-entrant faces, where `u` vanishes).
fn edge_singular(x: [f64; 3]) -> f64 {
    let (dx, dy) = (x[0] - 1.0, x[1] - 1.0);
    let r = dx.hypot(dy);
    if r < 1e-14 {
        return 0.0;
    }
    let mut theta = dy.atan2(dx);
    if theta < PI / 2.0 {
        theta += 2.0 * PI;
    }
    let phi = theta - PI / 2.0;
    r.powf(2.0 / 3.0) * (2.0 * phi / 3.0).sin()
}

/// Gradient of [`edge_singular`] (z-component zero).
fn grad_edge(x: [f64; 3]) -> [f64; 3] {
    let (dx, dy) = (x[0] - 1.0, x[1] - 1.0);
    let r = dx.hypot(dy);
    if r < 1e-13 {
        return [0.0; 3];
    }
    let mut theta = dy.atan2(dx);
    if theta < PI / 2.0 {
        theta += 2.0 * PI;
    }
    let phi = theta - PI / 2.0;
    let c = (2.0 / 3.0) * r.powf(-1.0 / 3.0);
    let (s, co) = (2.0 * phi / 3.0).sin_cos();
    let (st, ct) = theta.sin_cos();
    [c * (s * ct - co * st), c * (s * st + co * ct), 0.0]
}

#[allow(clippy::type_complexity)] // a tuple of closures is the natural return value here
fn sine_problem() -> (
    impl Fn([f64; 3]) -> f64,
    impl Fn(f64, [f64; 3]) -> f64,
    MarkerBcs<Bc>,
) {
    let exact = |x: [f64; 3]| (PI * x[0]).sin() * (PI * x[1]).sin() * (PI * x[2]).sin();
    let source = move |_t: f64, x: [f64; 3]| 3.0 * PI * PI * exact(x);
    let bcs = MarkerBcs::with_default(Bc::Dirichlet(BcValue::space(exact)));
    (exact, source, bcs)
}

#[test]
fn marked_refinement_is_conforming_and_conserves_volume() {
    let mesh = SimplexMesh::boxed([2, 2, 2], [0.0; 3], [1.0; 3]);
    let marked: Vec<bool> = mesh
        .cells
        .iter()
        .map(|c| {
            let x = centroid(&mesh, c);
            x[0] < 0.5 && x[1] < 0.5
        })
        .collect();
    assert!(marked.iter().any(|&m| m) && marked.iter().any(|&m| !m));
    let refined = mesh.refine_marked(&marked);
    assert!(refined.cells.len() > mesh.cells.len());
    assert!(refined.level.iter().all(|&l| l <= 1));
    assert!(refined.level.contains(&1));

    let interior = assert_conforming(&refined, &mesh, 1.0);
    assert!(interior > 0);

    // Boundary facets keep the marker of the box face they lie on.
    let mut seen = [false; 6];
    for &(ci, f, mk) in &refined.facets {
        assert!((0..6).contains(&mk), "unexpected marker {mk}");
        let cen = face_centroid(&refined, ci, f);
        let axis = (mk / 2) as usize;
        assert!(
            (cen[axis] - f64::from(mk % 2)).abs() < 1e-12,
            "marker {mk} on {cen:?}"
        );
        seen[mk as usize] = true;
    }
    assert!(seen.iter().all(|&s| s));

    // Every marked parent is covered by exactly eight children.
    let bary = |p: usize, x: [f64; 3]| -> [f64; 4] {
        let g = mesh.bary_gradients(p);
        let x0 = mesh.nodes[mesh.cells[p][0]];
        core::array::from_fn(|m| {
            (0..3).map(|a| g[m][a] * (x[a] - x0[a])).sum::<f64>() + if m == 0 { 1.0 } else { 0.0 }
        })
    };
    let mut checked = 0;
    for (p, &m) in marked.iter().enumerate() {
        if !m {
            continue;
        }
        let parent_volume = mesh.measure(p);
        let mut covered = 0.0;
        let mut children = 0;
        for ci in 0..refined.cells.len() {
            let x = centroid(&refined, &refined.cells[ci]);
            if bary(p, x).iter().all(|&l| l >= -1e-12) {
                children += 1;
                let vol = refined.measure(ci);
                assert!(
                    vol <= parent_volume + 1e-15,
                    "child of {p} larger than its parent"
                );
                covered += vol;
            }
        }
        assert_eq!(children, 8, "parent {p} has {children} children");
        assert!(
            (covered - parent_volume).abs() < 1e-12,
            "parent {p}: {covered} vs {parent_volume}"
        );
        checked += 1;
    }
    assert!(checked > 0);
}

#[test]
fn repeated_refinement_at_the_boundary_remains_conforming() {
    let mesh = SimplexMesh::boxed([2, 2, 2], [0.0; 3], [1.0; 3]);
    let marked: Vec<bool> = mesh
        .cells
        .iter()
        .map(|c| centroid(&mesh, c)[2] < 0.5)
        .collect();
    let first = mesh.refine_marked(&marked);
    assert_conforming(&first, &mesh, 1.0);
    assert!(first.level.contains(&1));
    assert!(first.level.contains(&0));

    // refine every second cell of the first mesh, boundary cells included
    let marked2: Vec<bool> = (0..first.cells.len()).map(|i| i % 2 == 0).collect();
    let second = first.refine_marked(&marked2);
    assert_conforming(&second, &first, 1.0);
    assert!(second.level.contains(&2));
    assert!(second.level.iter().all(|&l| l <= 2));

    // a uniform pass on a mesh produced by earlier refinements
    let third = second.refine_marked(&vec![true; second.cells.len()]);
    assert_conforming(&third, &second, 1.0);
    assert_eq!(third.level.iter().copied().max(), Some(3));

    // no marks is the identity
    let unchanged = second.refine_marked(&vec![false; second.cells.len()]);
    assert_eq!(unchanged.cells.len(), second.cells.len());
    assert_eq!(unchanged.level, second.level);
}

#[test]
fn boundary_markers_are_inherited_from_the_parent_facet() {
    let mut mesh = SimplexMesh::boxed([2, 2, 1], [0.0; 3], [1.0; 3]);
    for f in &mut mesh.facets {
        if f.2 == marker::ZMAX {
            f.2 = 42;
        } else if f.2 == marker::XMIN {
            f.2 = 7;
        }
    }
    let marked: Vec<bool> = mesh
        .cells
        .iter()
        .map(|c| centroid(&mesh, c)[0] < 0.5)
        .collect();
    let refined = mesh.refine_marked(&marked);
    let mut top = false;
    let mut left = false;
    for &(ci, f, mk) in &refined.facets {
        let cen = face_centroid(&refined, ci, f);
        if mk == 42 {
            assert!((cen[2] - 1.0).abs() < 1e-12, "marker 42 at {cen:?}");
            top = true;
        } else if mk == 7 {
            assert!(cen[0].abs() < 1e-12, "marker 7 at {cen:?}");
            left = true;
        }
    }
    assert!(top && left);
    assert_conforming(&refined, &mesh, 1.0);
}

#[test]
fn p1_poisson_sine_on_the_unit_cube_is_second_order() {
    let (exact, source, bcs) = sine_problem();
    let error = |n: usize| {
        let space = FemSpace::new(SimplexMesh::boxed([n; 3], [0.0; 3], [1.0; 3]), 1).unwrap();
        let prob = FemProblem {
            k: 1.0,
            c: 0.0,
            f: &source,
            bcs: &bcs,
        };
        let u = solve_static(&space, &prob, 0.0, LinearSolver::BiCgStab(1e-11)).unwrap();
        (
            space.l2_error(&u, &exact),
            h1_error_p1(&space, &u, grad_sine),
        )
    };
    let ((l2a, h1a), (l2b, h1b)) = (error(8), error(16));
    let (order_l2, order_h1) = ((l2a / l2b).log2(), (h1a / h1b).log2());
    eprintln!(
        "P1 cube: L2 {l2a:.3e} -> {l2b:.3e} (order {order_l2:.2}), H1 {h1a:.3e} -> {h1b:.3e} (order {order_h1:.2})"
    );
    assert!(order_l2 > 1.9, "P1 L2 order {order_l2} ({l2a}, {l2b})");
    assert!(order_h1 > 0.9, "P1 H1 order {order_h1} ({h1a}, {h1b})");
}

#[test]
fn p2_smoke_test_on_the_unit_cube() {
    let (exact, source, bcs) = sine_problem();
    let space = FemSpace::new(SimplexMesh::boxed([4, 4, 4], [0.0; 3], [1.0; 3]), 2).unwrap();
    let prob = FemProblem {
        k: 1.0,
        c: 0.0,
        f: &source,
        bcs: &bcs,
    };
    let u = solve_static(&space, &prob, 0.0, LinearSolver::Auto).unwrap();
    let err = space.l2_error(&u, exact);
    assert!(err < 1e-2, "P2 L2 error {err}");
    let eta = error_indicators(&space, &u, &prob, 0.0);
    assert_eq!(eta.len(), space.mesh.cells.len());
    assert!(eta.iter().all(|e| e.is_finite() && *e >= 0.0));
    assert!(eta.iter().map(|e| e * e).sum::<f64>() > 0.0);
}

#[test]
fn residual_estimator_tracks_the_true_error() {
    let (exact, source, bcs) = sine_problem();
    let run = |n: usize| {
        let space = FemSpace::new(SimplexMesh::boxed([n; 3], [0.0; 3], [1.0; 3]), 1).unwrap();
        let prob = FemProblem {
            k: 1.0,
            c: 0.0,
            f: &source,
            bcs: &bcs,
        };
        let u = solve_static(&space, &prob, 0.0, LinearSolver::Auto).unwrap();
        let eta = error_indicators(&space, &u, &prob, 0.0);
        let total: f64 = eta.iter().map(|e| e * e).sum::<f64>().sqrt();
        (
            total,
            h1_error_p1(&space, &u, grad_sine),
            space.l2_error(&u, &exact),
        )
    };
    let (coarse, mid, fine) = (run(2), run(4), run(8));
    eprintln!("estimator: {:?} {:?} {:?}", coarse, mid, fine);
    assert!(
        mid.0 < coarse.0 && fine.0 < mid.0,
        "estimator not decreasing: {} {} {}",
        coarse.0,
        mid.0,
        fine.0
    );
    assert!(fine.1 < mid.1 && mid.1 < coarse.1);
    // The indicator is an energy-type quantity: its ratio to the H1 error
    // must stay bounded when the mesh is refined.
    let r_coarse = coarse.0 / coarse.1;
    let r_fine = fine.0 / fine.1;
    assert!(
        (0.1..10.0).contains(&r_coarse) && (0.1..10.0).contains(&r_fine),
        "ratio {r_coarse} {r_fine}"
    );
    assert!(
        r_fine / r_coarse > 0.5 && r_fine / r_coarse < 2.0,
        "estimator strays from the H1 error: {r_coarse} -> {r_fine}"
    );
}

#[test]
fn adaptive_refinement_beats_uniform_for_an_edge_singularity() {
    let (exact, source, bcs) = (
        edge_singular,
        |_t: f64, _x: [f64; 3]| -> f64 { 0.0 },
        MarkerBcs::with_default(Bc::Dirichlet(BcValue::space(edge_singular))),
    );
    let prob = FemProblem {
        k: 1.0,
        c: 0.0,
        f: &source,
        bcs: &bcs,
    };
    let start = l_prism(8);
    let solve = |mesh: SimplexMesh| {
        let space = FemSpace::new(mesh, 1).unwrap();
        let u = solve_static(&space, &prob, 0.0, LinearSolver::Auto).unwrap();
        (
            space.l2_error(&u, exact),
            h1_error_p1(&space, &u, grad_edge),
            space.size(),
        )
    };
    // one uniform red refinement of the same start mesh
    let uniform = start.refine_marked(&vec![true; start.cells.len()]);
    let (err_u, h1_u, dofs_u) = solve(uniform);
    // adaptive refinement grades the mesh towards the re-entrant edge
    let (space_a, u_a, history) =
        solve_adaptive(&start, &prob, 1, 0.0, 0.5, 0.0, 3, LinearSolver::Auto).unwrap();
    let err_a = space_a.l2_error(&u_a, exact);
    let h1_a = h1_error_p1(&space_a, &u_a, grad_edge);
    let dofs_a = space_a.size();
    eprintln!(
        "adaptive: L2 {err_a:.3e} H1 {h1_a:.3e} at {dofs_a} dofs; uniform: L2 {err_u:.3e} H1 {h1_u:.3e} at {dofs_u} dofs"
    );
    eprintln!("history: {history:?}");
    assert!(
        dofs_a as f64 <= 1.3 * dofs_u as f64,
        "adaptive {dofs_a} vs uniform {dofs_u} dofs"
    );
    assert!(
        err_a < 0.75 * err_u,
        "adaptive L2 {err_a} at {dofs_a} dofs vs uniform {err_u} at {dofs_u}"
    );
    assert!(
        h1_a < 0.85 * h1_u,
        "adaptive H1 {h1_a} at {dofs_a} dofs vs uniform {h1_u} at {dofs_u}"
    );
    assert!(
        history.windows(2).all(|w| w[1].0 > w[0].0),
        "dofs did not grow"
    );
    assert!(
        history.last().unwrap().1 < history[0].1,
        "estimator did not decrease"
    );
    // repeated adaptive refinement stays conforming in the singular region
    assert_conforming(&space_a.mesh, &start, 3.0);
}

mod variable_coefficients {
    use super::*;
    use rssn::sim::fem_simplex::FemProblemVar;
    use rssn::sim::fem_simplex::solve_heat_var;
    use rssn::sim::fem_simplex::solve_static_var;

    fn all_dirichlet_zero() -> MarkerBcs<Bc> {
        let mut b = MarkerBcs::new();
        for m in [marker::XMIN, marker::XMAX, marker::YMIN, marker::YMAX] {
            b = b.set(m, Bc::Dirichlet(BcValue::constant(0.0)));
        }
        b
    }

    fn k_of(x: [f64; 3]) -> f64 {
        1.0 + x[0] * x[1] + 0.5 * x[1] * x[1]
    }

    fn c_of(x: [f64; 3]) -> f64 {
        1.0 + x[0] + x[1] * x[1]
    }

    /// `u = sin(pi x) sin(pi y)`, `f = -div(k grad u) + c u`.
    fn manufactured_error(
        degree: usize,
        n: usize,
    ) -> f64 {
        let space =
            FemSpace::new(SimplexMesh::rectangle(n, n, [0.0; 2], [1.0; 2]), degree).unwrap();
        let u = |x: [f64; 3]| (PI * x[0]).sin() * (PI * x[1]).sin();
        let f = |_t: f64, x: [f64; 3]| {
            let (sx, cx) = ((PI * x[0]).sin(), (PI * x[0]).cos());
            let (sy, cy) = ((PI * x[1]).sin(), (PI * x[1]).cos());
            let (ux, uy) = (PI * cx * sy, PI * sx * cy);
            let lap = -2.0 * PI * PI * sx * sy;
            // grad k = (y, x + y)
            let div = x[1] * ux + (x[0] + x[1]) * uy + k_of(x) * lap;
            -div + c_of(x) * u(x)
        };
        let bcs = all_dirichlet_zero();
        let prob = FemProblemVar {
            k: &k_of,
            c: &c_of,
            f: &f,
            bcs: &bcs,
        };
        let sol = solve_static_var(&space, &prob, 0.0, LinearSolver::Auto).unwrap();
        space.l2_error(&sol, u)
    }

    #[test]
    fn manufactured_solution_converges_with_optimal_order_for_p1_and_p2() {
        for (degree, ns, min) in [(1, [8, 16, 32], 1.85), (2, [4, 8, 16], 2.8)] {
            let e: Vec<f64> = ns.iter().map(|&n| manufactured_error(degree, n)).collect();
            let rates = [(e[0] / e[1]).log2(), (e[1] / e[2]).log2()];
            println!("P{degree} variable k, c: errors {e:?} rates {rates:?}");
            assert!(rates[1] > min, "P{degree}: {rates:?}");
            assert!(rates[0] > min - 0.3, "P{degree}: {rates:?}");
        }
    }

    #[test]
    fn constant_closures_reproduce_the_constant_api_bit_for_bit() {
        for degree in [1, 2] {
            let space =
                FemSpace::new(SimplexMesh::rectangle(6, 5, [0.0; 2], [1.0; 2]), degree).unwrap();
            let f = |_t: f64, x: [f64; 3]| 1.0 + x[0] * x[1];
            let mut bcs = all_dirichlet_zero();
            bcs = bcs.set(
                marker::XMAX,
                Bc::Robin {
                    alpha: 2.0,
                    beta: 1.0,
                    g: BcValue::constant(0.5),
                },
            );
            let constant = FemProblem {
                k: 2.5,
                c: 0.75,
                f: &f,
                bcs: &bcs,
            };
            let (k, c) = (|_: [f64; 3]| 2.5, |_: [f64; 3]| 0.75);
            let var = FemProblemVar {
                k: &k,
                c: &c,
                f: &f,
                bcs: &bcs,
            };
            let a = solve_static(&space, &constant, 0.0, LinearSolver::Auto).unwrap();
            let b = solve_static_var(&space, &var, 0.0, LinearSolver::Auto).unwrap();
            assert_eq!(a, b);
            let h1 =
                solve_heat(&space, &constant, |x| x[0], 0.1, 4, 0.5, LinearSolver::Auto).unwrap();
            let h2 =
                solve_heat_var(&space, &var, |x| x[0], 0.1, 4, 0.5, LinearSolver::Auto).unwrap();
            assert_eq!(h1, h2);
        }
    }

    #[test]
    fn layered_conductivity_gives_the_flux_continuous_solution() {
        // -(k u')' = 0 on (0,1), k = 1 on x < 1/2 and 4 beyond, u(0) = 0, u(1) = 1:
        // the flux q is constant, q (1/2 + 1/8) = 1, so u(1/2) = q/2 = 0.8.
        let space = FemSpace::new(SimplexMesh::interval(64, 0.0, 1.0), 1).unwrap();
        let bcs = MarkerBcs::new()
            .set(marker::XMIN, Bc::Dirichlet(BcValue::constant(0.0)))
            .set(marker::XMAX, Bc::Dirichlet(BcValue::constant(1.0)));
        let f = |_t: f64, _x: [f64; 3]| 0.0;
        let k = |x: [f64; 3]| if x[0] < 0.5 { 1.0 } else { 4.0 };
        let c = |_x: [f64; 3]| 0.0;
        let prob = FemProblemVar {
            k: &k,
            c: &c,
            f: &f,
            bcs: &bcs,
        };
        let u = solve_static_var(&space, &prob, 0.0, LinearSolver::Auto).unwrap();
        let mid = space.eval(&u, [0.5, 0.0, 0.0]).unwrap();
        assert!((mid - 0.8).abs() < 1e-10, "u(1/2) = {mid}");
    }
}

mod variable_estimator {
    use rssn::sim::fem_simplex::FemProblemVar;
    use rssn::sim::fem_simplex::error_indicators_var;
    use rssn::sim::fem_simplex::solve_adaptive_var;
    use rssn::sim::fem_simplex::solve_static_var;

    use super::*;

    /// Re-entrant corner-free test problem with a sharp interior layer in the diffusivity:
    /// `k = 1 + 40 exp(-80 ((x - .5)^2 + (y - .5)^2))` and the source of the exact solution
    /// `u = sin(pi x) sin(pi y)`.
    fn k_of(x: [f64; 3]) -> f64 {
        1.0 + 40.0 * (-80.0 * ((x[0] - 0.5).powi(2) + (x[1] - 0.5).powi(2))).exp()
    }

    fn dk_of(x: [f64; 3]) -> [f64; 2] {
        let e = -80.0 * ((x[0] - 0.5).powi(2) + (x[1] - 0.5).powi(2));
        let g = 40.0 * e.exp() * -160.0;
        [g * (x[0] - 0.5), g * (x[1] - 0.5)]
    }

    fn c_of(_x: [f64; 3]) -> f64 {
        0.5
    }

    fn exact(x: [f64; 3]) -> f64 {
        (PI * x[0]).sin() * (PI * x[1]).sin()
    }

    fn source(
        _t: f64,
        x: [f64; 3],
    ) -> f64 {
        let (sx, cx) = ((PI * x[0]).sin(), (PI * x[0]).cos());
        let (sy, cy) = ((PI * x[1]).sin(), (PI * x[1]).cos());
        let (ux, uy) = (PI * cx * sy, PI * sx * cy);
        let lap = -2.0 * PI * PI * sx * sy;
        let dk = dk_of(x);
        -(dk[0] * ux + dk[1] * uy + k_of(x) * lap) + c_of(x) * exact(x)
    }

    fn zero_boundary() -> MarkerBcs<Bc> {
        let mut b = MarkerBcs::new();
        for m in [marker::XMIN, marker::XMAX, marker::YMIN, marker::YMAX] {
            b = b.set(m, Bc::Dirichlet(BcValue::constant(0.0)));
        }
        b
    }

    #[test]
    fn the_estimator_matches_the_constant_one_and_sees_the_layer() {
        let bcs = zero_boundary();
        // constant closures: bit-identical to the constant estimator
        let space = FemSpace::new(SimplexMesh::rectangle(6, 6, [0.0; 2], [1.0; 2]), 2).unwrap();
        let f = |_t: f64, x: [f64; 3]| 1.0 + x[0] * x[1];
        let constant = FemProblem {
            k: 2.5,
            c: 0.75,
            f: &f,
            bcs: &bcs,
        };
        let (k, c) = (|_: [f64; 3]| 2.5, |_: [f64; 3]| 0.75);
        let var = FemProblemVar {
            k: &k,
            c: &c,
            f: &f,
            bcs: &bcs,
        };
        let u = solve_static(&space, &constant, 0.0, LinearSolver::Auto).unwrap();
        assert_eq!(
            error_indicators(&space, &u, &constant, 0.0),
            error_indicators_var(&space, &u, &var, 0.0)
        );
        // the variable problem: the indicators are largest where k varies fast
        let (kf, cf) = (&k_of, &c_of);
        let layer = FemProblemVar {
            k: kf,
            c: cf,
            f: &source,
            bcs: &bcs,
        };
        let space = FemSpace::new(SimplexMesh::rectangle(16, 16, [0.0; 2], [1.0; 2]), 1).unwrap();
        let u = solve_static_var(&space, &layer, 0.0, LinearSolver::Auto).unwrap();
        let eta = error_indicators_var(&space, &u, &layer, 0.0);
        let total: f64 = eta.iter().map(|e| e * e).sum::<f64>().sqrt();
        let (mut near, mut far, mut nn, mut nf) = (0.0, 0.0, 0, 0);
        for (ci, e) in eta.iter().enumerate() {
            let cell = &space.mesh.cells[ci];
            let mut m = [0.0; 3];
            for &v in cell {
                for a in 0..3 {
                    m[a] += space.mesh.nodes[v][a] / cell.len() as f64;
                }
            }
            let r = ((m[0] - 0.5).powi(2) + (m[1] - 0.5).powi(2)).sqrt();
            if r < 0.3 {
                near += e * e;
                nn += 1;
            } else {
                far += e * e;
                nf += 1;
            }
        }
        println!(
            "estimator: total {total:.3e}, mean eta^2 near the layer {:.2e}, far {:.2e}",
            near / nn as f64,
            far / nf as f64
        );
        assert!(near / nn as f64 > 2.0 * far / nf as f64);
    }

    #[test]
    fn adaptive_refinement_converges_for_a_variable_diffusivity() {
        let bcs = zero_boundary();
        let (kf, cf) = (&k_of, &c_of);
        let prob = FemProblemVar {
            k: kf,
            c: cf,
            f: &source,
            bcs: &bcs,
        };
        let mesh = SimplexMesh::rectangle(4, 4, [0.0; 2], [1.0; 2]);
        let (space, u, history) =
            solve_adaptive_var(&mesh, &prob, 1, 0.0, 0.5, 0.0, 6, LinearSolver::Auto).unwrap();
        println!("history (dofs, estimator): {history:?}");
        // the estimator and the true error both decrease
        assert!(history.windows(2).all(|w| w[1].1 < w[0].1), "{history:?}");
        let err = space.l2_error(&u, exact);
        let coarse = {
            let s0 = FemSpace::new(mesh.clone(), 1).unwrap();
            let u0 = solve_static_var(&s0, &prob, 0.0, LinearSolver::Auto).unwrap();
            s0.l2_error(&u0, exact)
        };
        assert!(err < 0.5 * coarse, "L2 error {err} vs coarse {coarse}");
        // the estimator tracks the error: the ratio stays within a factor ~20 from first to last
        let ratio = |e: f64, est: f64| est / e;
        let r_last = ratio(err, history.last().unwrap().1);
        let r_first = ratio(coarse, history[0].1);
        assert!(
            r_last / r_first < 20.0 && r_first / r_last < 20.0,
            "{r_first} vs {r_last}"
        );
        // invalid controls are rejected
        assert!(solve_adaptive_var(&mesh, &prob, 1, 0.0, 0.5, 0.0, 0, LinearSolver::Auto).is_err());
        assert!(solve_adaptive_var(&mesh, &prob, 1, 0.0, 0.0, 0.0, 3, LinearSolver::Auto).is_err());
    }
}
