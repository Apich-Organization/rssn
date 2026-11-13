//! Sparse and preconditioned linear solves, fill-reducing ordering and
//! polynomial systems.
//!
//! | operator | value |
//! |---|---|
//! | `sparse_solve(entries, b, tol, method)` | the solution of `A x = b` for a sparse `A` given as the triplets `list(list(i, j, v), ..)` (0-based, duplicates summed, order `len(b)`), by a Krylov `method` with the Jacobi (diagonal) preconditioner: `minres`, `fgmres`, `bicgstab_l` or `idr_s` |
//! | `amg_solve(entries, b, tol, kind)` | the same system by general algebraic multigrid (`ruge_stuben` or `smoothed_aggregation`); the diagonal must be positive |
//! | `nested_dissection_order(adjacency)` | a fill-reducing elimination order `perm` (`perm[new] = old`, 0-based) of the symmetric pattern given as one list of neighbour indices per vertex, by multilevel nested dissection |
//! | `polynomial_system_roots(system)` | the finite roots of a square polynomial system by total-degree homotopy continuation, as `list(list(list(re, im), ..), ..)`, one entry per root and coordinate |
//!
//! A polynomial `system` is a list of equations, each a list of terms
//! `list(coefficient, e_1, .., e_n)` for the monomial `x_1^e_1 .. x_n^e_n`;
//! the number of equations is the number of variables. The roots are those
//! whose path converged to a finite endpoint: roots at infinity are dropped,
//! and a multiple root is listed once per path that reached it.
//!
//! The solvers return an answer only when the iteration reached `tol`; the
//! sizes are capped so one call stays bounded.

use super::bridge::Data;
use super::bridge::Entry;
use super::bridge::name;
use super::bridge::rows;
use super::bridge::scalar;
use super::bridge::vec_out;
use super::bridge::vector;
use crate::kernels::amg_general::AmgGeneral;
use crate::kernels::amg_general::AmgGeneralOptions;
use crate::kernels::amg_general::AmgKind;
use crate::kernels::homotopy;
use crate::kernels::krylov::Csr;
use crate::kernels::krylov_ext;
use crate::kernels::ordering_ml;
use crate::kernels::rootfind::PolyTerm;

/// Most rows of a sparse system.
const MAX_ORDER: usize = 20_000;
/// Most stored entries of a sparse matrix.
const MAX_ENTRIES: usize = 200_000;
/// Most iterations (or cycles) of a sparse solve.
const MAX_ITERATIONS: usize = 5_000;

/// The matrix `entries` of order `b.len()` and its diagonal.
fn sparse_matrix(
    a: &[Data],
    n: usize,
) -> Option<(Csr, Vec<f64>)> {
    let triplets = rows(a.first()?)?;
    (!triplets.is_empty() && triplets.len() <= MAX_ENTRIES).then_some(())?;
    let mut entries = Vec::with_capacity(triplets.len());
    let mut diagonal = vec![0.0; n];
    for t in &triplets {
        let &[i, j, v] = t.as_slice() else {
            return None;
        };
        let index = |x: f64| (x >= 0.0 && x.fract() == 0.0 && x < n as f64).then_some(x as usize);
        let (i, j) = (index(i)?, index(j)?);
        if i == j {
            diagonal[i] += v;
        }
        entries.push((i, j, v));
    }
    Some((Csr::from_triplets(n, n, &entries), diagonal))
}

/// The right-hand side and the tolerance shared by the sparse solvers.
fn right_hand_side(a: &[Data]) -> Option<(&[f64], f64)> {
    let b = vector(a, 1)?;
    (!b.is_empty() && b.len() <= MAX_ORDER).then_some(())?;
    let tol = scalar(a, 2)?;
    (tol > 0.0 && tol < 1.0).then_some((b, tol))
}

fn finished(
    x: Vec<f64>,
    converged: bool,
) -> Option<Data> {
    (converged && x.iter().all(|v| v.is_finite())).then(|| Data::vector(x))
}

fn sparse_solve(a: &[Data]) -> Option<Data> {
    let (b, tol) = right_hand_side(a)?;
    let method = name(a, 3)?;
    let n = b.len();
    let (matrix, diagonal) = sparse_matrix(a, n)?;
    let apply = |x: &[f64], y: &mut [f64]| matrix.matvec(x, y);
    // The Jacobi preconditioner needs a nonzero diagonal, and a positive one
    // for MINRES, which wants a positive definite preconditioner.
    let usable = diagonal.iter().all(|&d| d.is_finite() && d != 0.0)
        && (method != "minres" || diagonal.iter().all(|&d| d > 0.0));
    let precondition = |r: &[f64], z: &mut [f64]| {
        for ((zi, &ri), &d) in z.iter_mut().zip(r).zip(&diagonal) {
            *zi = if usable { ri / d } else { ri };
        }
    };
    let max_iter = (10 * n + 100).min(MAX_ITERATIONS);
    let result = match method {
        | "minres" => krylov_ext::minres(apply, precondition, b, None, tol, max_iter),
        | "fgmres" => krylov_ext::fgmres(apply, precondition, b, None, n.min(50), tol, max_iter),
        | "bicgstab_l" => krylov_ext::bicgstab_l(apply, precondition, 2, b, None, tol, max_iter),
        | "idr_s" => {
            krylov_ext::idr_s_with(
                apply,
                precondition,
                b,
                None,
                &krylov_ext::IdrOptions {
                    s: 4,
                    tol,
                    max_iter,
                    seed: 1,
                },
            )
        },
        | _ => return None,
    };
    finished(result.x, result.converged)
}

fn amg_solve(a: &[Data]) -> Option<Data> {
    let (b, tol) = right_hand_side(a)?;
    let kind = match name(a, 3)? {
        | "ruge_stuben" => AmgKind::RugeStuben,
        | "smoothed_aggregation" => AmgKind::SmoothedAggregation,
        | _ => return None,
    };
    let (matrix, diagonal) = sparse_matrix(a, b.len())?;
    diagonal
        .iter()
        .all(|&d| d.is_finite() && d > 0.0)
        .then_some(())?;
    let options = AmgGeneralOptions {
        kind,
        ..AmgGeneralOptions::default()
    };
    let hierarchy = AmgGeneral::build(&matrix, &options).ok()?;
    let solved = hierarchy.solve(b, None, tol, MAX_ITERATIONS.min(500));
    finished(solved.x, solved.converged)
}

/// Most vertices of an ordered graph and most adjacency entries.
const MAX_VERTICES: usize = 50_000;
const MAX_NEIGHBOURS: usize = 500_000;

fn nested_dissection_order(a: &[Data]) -> Option<Data> {
    let lists = rows(a.first()?)?;
    let n = lists.len();
    ((1..=MAX_VERTICES).contains(&n)).then_some(())?;
    let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut stored = 0_usize;
    for (v, list) in lists.iter().enumerate() {
        stored += list.len();
        (stored <= MAX_NEIGHBOURS).then_some(())?;
        for &x in list {
            let w = (x >= 0.0 && x.fract() == 0.0 && x < n as f64).then_some(x as usize)?;
            if w != v {
                adjacency[v].push(w);
                adjacency[w].push(v);
            }
        }
    }
    for list in &mut adjacency {
        list.sort_unstable();
        list.dedup();
    }
    let perm = ordering_ml::multilevel_nested_dissection_order(&adjacency);
    let mut seen = vec![false; n];
    let valid = perm.len() == n
        && perm
            .iter()
            .all(|&v| v < n && !std::mem::replace(&mut seen[v], true));
    valid.then_some(())?;
    vec_out(perm.iter().map(|&v| v as f64).collect())
}

/// Most variables, terms per equation, exponent, and paths of a polynomial system.
const MAX_VARIABLES: usize = 6;
const MAX_TERMS: usize = 40;
const MAX_EXPONENT: f64 = 8.0;
const MAX_PATHS: u64 = 128;

fn polynomial_system_roots(a: &[Data]) -> Option<Data> {
    let equations = super::bridge::parts_of(a.first()?)?;
    let n = equations.len();
    ((1..=MAX_VARIABLES).contains(&n)).then_some(())?;
    let mut system: Vec<Vec<PolyTerm>> = Vec::with_capacity(n);
    let mut paths = 1_u64;
    for equation in &equations {
        let terms = rows(equation)?;
        (!terms.is_empty() && terms.len() <= MAX_TERMS).then_some(())?;
        let mut parsed = Vec::with_capacity(terms.len());
        let mut degree = 0_u32;
        for term in &terms {
            let (&coef, exponents) = term.split_first()?;
            (exponents.len() == n && coef.is_finite()).then_some(())?;
            let exps: Option<Vec<u32>> = exponents
                .iter()
                .map(|&e| (e >= 0.0 && e.fract() == 0.0 && e <= MAX_EXPONENT).then_some(e as u32))
                .collect();
            let exps = exps?;
            degree = degree.max(exps.iter().sum());
            parsed.push(PolyTerm { coef, exps });
        }
        (degree >= 1).then_some(())?;
        paths = paths.checked_mul(u64::from(degree))?;
        (paths <= MAX_PATHS).then_some(())?;
        system.push(parsed);
    }
    let ends = homotopy::homotopy_solve(&system, &homotopy::HomotopyOptions::default()).ok()?;
    let roots: Vec<&homotopy::PathEnd> = ends
        .iter()
        .filter(|end| {
            end.converged
                && !end.at_infinity
                && end.x.len() == n
                && end
                    .x
                    .iter()
                    .all(|&(re, im)| re.is_finite() && im.is_finite())
        })
        .collect();
    let data = roots
        .iter()
        .flat_map(|end| end.x.iter().flat_map(|&(re, im)| [re, im]))
        .collect();
    Some(Data::Array {
        shape: vec![roots.len(), n, 2],
        data,
    })
}

pub(super) const ENTRIES: &[Entry] = &[
    Entry {
        name: "sparse_solve",
        params: "nvsw",
        eval: sparse_solve,
    },
    Entry {
        name: "amg_solve",
        params: "nvsw",
        eval: amg_solve,
    },
    Entry {
        name: "nested_dissection_order",
        params: "n",
        eval: nested_dissection_order,
    },
    Entry {
        name: "polynomial_system_roots",
        params: "n",
        eval: polynomial_system_roots,
    },
];

#[cfg(test)]
mod tests {
    use crate::kernels::ordering;
    use crate::rules::numeric::numeric;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[numeric()], src)
    }

    fn flat(text: &str) -> Vec<f64> {
        text.replace(['(', ')', ','], " ")
            .split_whitespace()
            .filter(|word| *word != "list")
            .map(|word| word.parse().unwrap_or(f64::NAN))
            .collect()
    }

    /// A matrix as triplets and as a dense nested list.
    fn matrices(dense: &[Vec<f64>]) -> (String, String) {
        let mut triplets = Vec::new();
        for (i, row) in dense.iter().enumerate() {
            for (j, v) in row.iter().enumerate() {
                if *v != 0.0 {
                    triplets.push(format!("list({i}, {j}, {v})"));
                }
            }
        }
        let rows: Vec<String> = dense
            .iter()
            .map(|row| {
                let cells: Vec<String> = row.iter().map(|v| format!("{v}")).collect();
                format!("list({})", cells.join(", "))
            })
            .collect();
        (
            format!("list({})", triplets.join(", ")),
            format!("list({})", rows.join(", ")),
        )
    }

    fn laplacian(n: usize) -> Vec<Vec<f64>> {
        (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        if i == j {
                            2.0
                        } else if i.abs_diff(j) == 1 {
                            -1.0
                        } else {
                            0.0
                        }
                    })
                    .collect()
            })
            .collect()
    }

    fn ones(n: usize) -> String {
        format!("list({})", vec!["1"; n].join(", "))
    }

    #[test]
    fn sparse_krylov_solves_agree_with_the_dense_ones() {
        let n = 12;
        let (sparse, dense) = matrices(&laplacian(n));
        let b = ones(n);
        let want = flat(&run(&format!("minres_solve({dense}, {b}, 1e-12)")));
        // The exact solution of -u'' = 1 on the grid: u_i = (i + 1)(n - i) / 2.
        for (i, w) in want.iter().enumerate() {
            assert!(
                (w - (i + 1) as f64 * (n - i) as f64 / 2.0).abs() < 1e-6,
                "{want:?}"
            );
        }
        for method in ["minres", "fgmres", "bicgstab_l", "idr_s"] {
            let got = flat(&run(&format!(
                "sparse_solve({sparse}, {b}, 1e-10, {method})"
            )));
            assert_eq!(got.len(), n, "{method}");
            for (g, w) in got.iter().zip(&want) {
                assert!((g - w).abs() < 1e-6, "{method}: {got:?} vs {want:?}");
            }
        }
        // A nonsymmetric, diagonally dominant matrix against the dense FGMRES.
        let skew: Vec<Vec<f64>> = vec![
            vec![4.0, 1.0, 0.0, 0.0],
            vec![-1.0, 4.0, 2.0, 0.0],
            vec![0.0, 0.5, 5.0, 1.0],
            vec![0.0, 0.0, -2.0, 6.0],
        ];
        let (sparse, dense) = matrices(&skew);
        let b = "list(1, 2, 3, 4)";
        let want = flat(&run(&format!("fgmres_solve({dense}, {b}, 1e-12)")));
        for method in ["fgmres", "bicgstab_l", "idr_s"] {
            let got = flat(&run(&format!(
                "sparse_solve({sparse}, {b}, 1e-10, {method})"
            )));
            assert_eq!(got.len(), 4, "{method}");
            for (g, w) in got.iter().zip(&want) {
                assert!((g - w).abs() < 1e-7, "{method}: {got:?} vs {want:?}");
            }
        }
        // Duplicated triplets are summed: two halves of the diagonal.
        let got = flat(&run(
            "sparse_solve(list(list(0, 0, 1), list(0, 0, 1)), list(4), 1e-10, fgmres)",
        ));
        assert!((got[0] - 2.0).abs() < 1e-9, "{got:?}");
    }

    #[test]
    fn multigrid_solves_a_poisson_matrix() {
        let n = 120;
        let (sparse, _) = matrices(&laplacian(n));
        let b = ones(n);
        for kind in ["ruge_stuben", "smoothed_aggregation"] {
            let got = flat(&run(&format!("amg_solve({sparse}, {b}, 1e-8, {kind})")));
            assert_eq!(got.len(), n, "{kind}");
            for (i, g) in got.iter().enumerate() {
                let exact = (i + 1) as f64 * (n - i) as f64 / 2.0;
                assert!(
                    (g - exact).abs() < 1e-4 * exact.max(1.0),
                    "{kind}: {i} {g} vs {exact}"
                );
            }
        }
    }

    #[test]
    fn nested_dissection_is_a_permutation_with_less_fill_than_the_natural_order() {
        // A 10 x 10 grid graph, vertex (r, c) = 10 r + c, listed with ragged neighbour lists.
        let side = 10_usize;
        let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); side * side];
        for r in 0..side {
            for c in 0..side {
                let v = r * side + c;
                if c + 1 < side {
                    adjacency[v].push(v + 1);
                }
                if r + 1 < side {
                    adjacency[v].push(v + side);
                }
            }
        }
        let text: Vec<String> = adjacency
            .iter()
            .map(|list| {
                let items: Vec<String> = list.iter().map(ToString::to_string).collect();
                format!("list({})", items.join(", "))
            })
            .collect();
        let perm = flat(&run(&format!(
            "nested_dissection_order(list({}))",
            text.join(", ")
        )));
        let mut sorted: Vec<usize> = perm.iter().map(|&v| v as usize).collect();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..side * side).collect::<Vec<_>>());
        // Symmetrised pattern for the fill count.
        let mut full = adjacency.clone();
        for (v, list) in adjacency.iter().enumerate() {
            for &w in list {
                full[w].push(v);
            }
        }
        let natural: Vec<usize> = (0..side * side).collect();
        let order: Vec<usize> = perm.iter().map(|&v| v as usize).collect();
        assert!(
            ordering::cholesky_fill(&full, &order) < ordering::cholesky_fill(&full, &natural),
            "ordering did not reduce the fill"
        );
        // An isolated vertex (an empty neighbour list) is still ordered.
        assert_eq!(
            flat(&run(
                "nested_dissection_order(list(list(1), list(0), list()))"
            ))
            .len(),
            3
        );
    }

    #[test]
    fn polynomial_systems_are_solved_by_homotopy_continuation() {
        // x^2 - 3x + 2 = 0: the roots 1 and 2.
        let out = flat(&run(
            "polynomial_system_roots(list(list(list(1, 2), list(-3, 1), list(2, 0))))",
        ));
        assert_eq!(out.len(), 4, "{out:?}");
        let mut real: Vec<f64> = vec![out[0], out[2]];
        real.sort_by(f64::total_cmp);
        assert!(
            (real[0] - 1.0).abs() < 1e-8 && (real[1] - 2.0).abs() < 1e-8,
            "{out:?}"
        );
        assert!(out[1].abs() < 1e-8 && out[3].abs() < 1e-8);
        // x^2 + 1 = 0: the roots are +-i.
        let out = flat(&run(
            "polynomial_system_roots(list(list(list(1, 2), list(1, 0))))",
        ));
        assert_eq!(out.len(), 4, "{out:?}");
        assert!(
            out[0].abs() < 1e-8 && out[2].abs() < 1e-8 && (out[1] + out[3]).abs() < 1e-8,
            "{out:?}"
        );
        // x^2 + y^2 = 1, x y = 1/4: four real roots, checked by substitution.
        let out = flat(&run(
            "polynomial_system_roots(list(list(list(1, 2, 0), list(1, 0, 2), list(-1, 0, 0)), list(list(1, 1, 1), list(-0.25, 0, 0))))",
        ));
        assert_eq!(out.len(), 4 * 2 * 2, "{out:?}");
        for root in out.chunks(4) {
            let (x, y) = (
                num_complex::Complex64::new(root[0], root[1]),
                num_complex::Complex64::new(root[2], root[3]),
            );
            assert!((x * x + y * y - 1.0).norm() < 1e-8, "{root:?}");
            assert!((x * y - 0.25).norm() < 1e-8, "{root:?}");
        }
    }

    #[test]
    fn invalid_systems_are_left_as_requests() {
        for src in [
            "sparse_solve(list(list(0, 0, 1)), list(1, 2), 1e-8, minres)",
            "sparse_solve(list(list(0, 5, 1)), list(1), 1e-8, minres)",
            "sparse_solve(list(list(0, 0, 1)), list(1), 1e-8, nonesuch)",
            "sparse_solve(list(list(0, 0, 1)), list(1), 2, fgmres)",
            "sparse_solve(list(list(0, 0, a)), list(1), 1e-8, fgmres)",
            "sparse_solve(list(), list(1), 1e-8, fgmres)",
            "amg_solve(list(list(0, 0, -1)), list(1), 1e-8, ruge_stuben)",
            "amg_solve(list(list(0, 0, 1)), list(1), 1e-8, nonesuch)",
            "nested_dissection_order(list(list(7), list(0)))",
            "nested_dissection_order(list(1, 2))",
            "polynomial_system_roots(list(list(list(1, 2, 3))))",
            "polynomial_system_roots(list(list(list(1, 20))))",
            "polynomial_system_roots(list(list(list(1, 0))))",
            // 6 x 6 x 6 = 216 paths is over the cap.
            "polynomial_system_roots(list(list(list(1, 6, 0, 0), list(1, 0, 6, 0), list(1, 0, 0, 6)), list(list(1, 6, 0, 0), list(1, 0, 6, 0), list(1, 0, 0, 6)), list(list(1, 6, 0, 0), list(1, 0, 6, 0), list(1, 0, 0, 6))))",
        ] {
            let (text, reduced) = reduce_with(&[numeric()], src, &[]);
            assert!(!reduced, "{text}");
        }
    }
}
