//! Tests for the AMD and nested dissection orderings.

use rssn::kernels::krylov::{Csr, SparseOrdering};
use rssn::kernels::ordering::{amd_order, cholesky_fill, nested_dissection_order};
use rssn::kernels::ordering_ml::{
    hybrid_nested_dissection_order, multilevel_nested_dissection_order, vertex_separator_bisection, NdOptions,
};

fn grid_adj(nx: usize, ny: usize) -> Vec<Vec<usize>> {
    let mut adj = vec![Vec::new(); nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let v = j * nx + i;
            if i + 1 < nx {
                adj[v].push(v + 1);
                adj[v + 1].push(v);
            }
            if j + 1 < ny {
                adj[v].push(v + nx);
                adj[v + nx].push(v);
            }
        }
    }
    adj
}

fn is_perm(p: &[usize], n: usize) -> bool {
    let mut seen = vec![false; n];
    p.len() == n
        && p.iter().all(|&v| {
            let ok = v < n && !seen[v];
            if ok {
                seen[v] = true;
            }
            ok
        })
}

fn grid_csr(nx: usize, ny: usize) -> Csr {
    let mut t = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            let v = j * nx + i;
            t.push((v, v, 4.0));
            if i + 1 < nx {
                t.push((v, v + 1, -1.0));
                t.push((v + 1, v, -1.0));
            }
            if j + 1 < ny {
                t.push((v, v + nx, -1.0));
                t.push((v + nx, v, -1.0));
            }
        }
    }
    Csr::from_triplets(nx * ny, nx * ny, &t)
}

#[test]
fn amd_reduces_fill_on_grids() {
    let adj = grid_adj(60, 60);
    let n = adj.len();
    let nat: Vec<usize> = (0..n).collect();
    let amd = amd_order(&adj);
    assert!(is_perm(&amd, n));
    let (f_nat, f_amd) = (cholesky_fill(&adj, &nat), cholesky_fill(&adj, &amd));
    assert!(f_amd * 2 < f_nat, "amd {f_amd} vs natural {f_nat}");
    // comparable with exact minimum degree
    let a = grid_csr(60, 60);
    let md = a.min_degree_ordering();
    let f_md = cholesky_fill(&adj, &md);
    assert!((f_amd as f64) < 1.35 * f_md as f64, "amd {f_amd} vs exact md {f_md}");
}

#[test]
fn nested_dissection_on_grids() {
    let adj = grid_adj(90, 90);
    let n = adj.len();
    let nat: Vec<usize> = (0..n).collect();
    let nd = nested_dissection_order(&adj);
    assert!(is_perm(&nd, n));
    let (f_nat, f_nd) = (cholesky_fill(&adj, &nat), cholesky_fill(&adj, &nd));
    let f_amd = cholesky_fill(&adj, &amd_order(&adj));
    assert!(f_nd * 2 < f_nat, "nd {f_nd} vs natural {f_nat}");
    assert!((f_nd as f64) < 2.5 * f_amd as f64, "nd {f_nd} vs amd {f_amd}");
}

#[test]
fn amd_structured_graphs() {
    // arrow matrix: the dense row/column must go last, fill-free
    let n = 400;
    let mut adj = vec![Vec::new(); n];
    for i in 1..n {
        adj[0].push(i);
        adj[i].push(0);
    }
    let p = amd_order(&adj);
    assert!(is_perm(&p, n));
    assert_eq!(cholesky_fill(&adj, &p), 2 * n - 1);
    // disjoint cliques (supervariables): no fill beyond the cliques
    let (k, m) = (7usize, 9usize);
    let mut adj = vec![Vec::new(); k * m];
    for b in 0..k {
        for i in 0..m {
            for j in 0..m {
                if i != j {
                    adj[b * m + i].push(b * m + j);
                }
            }
        }
    }
    let p = amd_order(&adj);
    assert!(is_perm(&p, k * m));
    assert_eq!(cholesky_fill(&adj, &p), k * m * (m + 1) / 2);
    // trees admit a perfect elimination ordering
    let n = 1000;
    let mut adj = vec![Vec::new(); n];
    for i in 1..n {
        let par = (i - 1) / 3;
        adj[i].push(par);
        adj[par].push(i);
    }
    let p = amd_order(&adj);
    assert_eq!(cholesky_fill(&adj, &p), 2 * n - 1);
    assert!(amd_order(&[]).is_empty());
    assert_eq!(amd_order(&[Vec::new(), Vec::new()]).len(), 2);
}

#[test]
fn amd_scales_to_millions_of_nodes() {
    // a path of two million vertices and a long thin grid
    let n = 2_000_000;
    let mut adj = vec![Vec::new(); n];
    for i in 1..n {
        adj[i].push(i - 1);
        adj[i - 1].push(i);
    }
    let p = amd_order(&adj);
    assert!(is_perm(&p, n));
    let adj = grid_adj(5, 100_000);
    let p = amd_order(&adj);
    assert!(is_perm(&p, adj.len()));
}

#[test]
fn sparse_lu_with_new_orderings() {
    let a = grid_csr(30, 30);
    let n = 900;
    let xs: Vec<f64> = (0..n).map(|i| (i as f64 * 0.37).sin()).collect();
    let mut b = vec![0.0; n];
    a.matvec(&xs, &mut b);
    let mut nnz = Vec::new();
    for ord in [SparseOrdering::Natural, SparseOrdering::Amd, SparseOrdering::NestedDissection] {
        let lu = a.sparse_lu(ord).unwrap();
        let x = lu.solve(&b);
        let err = x.iter().zip(&xs).map(|(u, v)| (u - v).abs()).fold(0.0, f64::max);
        assert!(err < 1e-10, "{ord:?}: {err}");
        nnz.push(lu.nnz());
    }
    assert!(nnz[1] < nnz[0] && nnz[2] < nnz[0], "{nnz:?}");
}

fn grid3_adj(n: usize) -> Vec<Vec<usize>> {
    let idx = |i: usize, j: usize, k: usize| (k * n + j) * n + i;
    let mut adj = vec![Vec::new(); n * n * n];
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let v = idx(i, j, k);
                for (di, dj, dk) in [(1, 0, 0), (0, 1, 0), (0, 0, 1)] {
                    let (a, b, c) = (i + di, j + dj, k + dk);
                    if a < n && b < n && c < n {
                        let u = idx(a, b, c);
                        adj[v].push(u);
                        adj[u].push(v);
                    }
                }
            }
        }
    }
    adj
}

#[test]
fn multilevel_nd_fill_matches_amd_on_grids() {
    let adj = grid_adj(150, 150);
    let n = adj.len();
    let ml = multilevel_nested_dissection_order(&adj);
    assert!(is_perm(&ml, n));
    let f_ml = cholesky_fill(&adj, &ml);
    let f_amd = cholesky_fill(&adj, &amd_order(&adj));
    let f_nd = cholesky_fill(&adj, &nested_dissection_order(&adj));
    eprintln!("2D 150x150: multilevel {f_ml}, amd {f_amd}, level-structure nd {f_nd}");
    assert!((f_ml as f64) < 1.15 * f_amd as f64, "multilevel {f_ml} vs amd {f_amd}");
    assert!(f_ml * 2 < f_nd, "multilevel {f_ml} vs level-structure {f_nd}");
    let f_h = cholesky_fill(&adj, &hybrid_nested_dissection_order(&adj));
    assert!(f_h <= f_amd);
    // optimal top-level separator of a k x k mesh
    let (l, r, sp) = vertex_separator_bisection(&adj, &NdOptions::default());
    assert_eq!(sp.len(), 150, "separator");
    assert!(l.len().abs_diff(r.len()) < 1500);

    let adj = grid3_adj(22);
    let n = adj.len();
    let ml = multilevel_nested_dissection_order(&adj);
    assert!(is_perm(&ml, n));
    let f_ml = cholesky_fill(&adj, &ml);
    let f_amd = cholesky_fill(&adj, &amd_order(&adj));
    eprintln!("3D 22^3: multilevel {f_ml}, amd {f_amd}");
    assert!((f_ml as f64) < 1.15 * f_amd as f64, "multilevel {f_ml} vs amd {f_amd}");
    let f_h = cholesky_fill(&adj, &hybrid_nested_dissection_order(&adj));
    assert!(f_h <= f_amd);
    let (_, _, sp) = vertex_separator_bisection(&adj, &NdOptions::default());
    assert_eq!(sp.len(), 22 * 22, "separator");
}

#[test]
fn multilevel_nd_handles_irregular_graphs() {
    // disconnected union of a grid, a path and isolated vertices
    let g = grid_adj(30, 30);
    let n0 = g.len();
    let mut adj = g;
    let path = 200;
    for i in 0..path {
        adj.push(Vec::new());
        let v = n0 + i;
        if i > 0 {
            adj[v].push(v - 1);
            adj[v - 1].push(v);
        }
    }
    for _ in 0..7 {
        adj.push(Vec::new());
    }
    let n = adj.len();
    let p = multilevel_nested_dissection_order(&adj);
    assert!(is_perm(&p, n));
    let f = cholesky_fill(&adj, &p);
    let f_nat = cholesky_fill(&adj, &(0..n).collect::<Vec<_>>());
    assert!(f < f_nat);
    // tiny and empty inputs
    assert!(multilevel_nested_dissection_order(&[]).is_empty());
    assert_eq!(multilevel_nested_dissection_order(&[Vec::new()]), vec![0]);
    // a star (hub must be last)
    let m = 500;
    let mut star = vec![Vec::new(); m];
    for i in 1..m {
        star[0].push(i);
        star[i].push(0);
    }
    let p = multilevel_nested_dissection_order(&star);
    assert!(is_perm(&p, m));
    assert_eq!(cholesky_fill(&star, &p), 2 * m - 1);
}
