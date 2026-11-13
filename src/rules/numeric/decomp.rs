//! Rank-revealing and multilinear decompositions.
//!
//! | operator | value |
//! |---|---|
//! | `qr_pivoted_rank(A, tol)` | the numerical rank of `A`: the diagonal entries of `R` in the column-pivoted QR above `tol` times the largest |
//! | `qr_pivoted_permutation(A)` | the column order of that factorisation, 0-based: column `j` of `A P` is column `perm[j]` of `A` |
//! | `qr_pivoted_diagonal(A)` | the absolute diagonal of `R`, non-increasing, whose decay reveals the rank |
//! | `hosvd(T, ranks)` | `list(core, list(U_1, .., U_N), list(s_1, .., s_N))`: the truncated higher-order SVD of the tensor `T`, with the singular values of each mode unfolding (a rank of `0` keeps the full rank of its mode; a rank above the numerical one is clamped) |
//! | `hooi(T, ranks, sweeps, tol)` | `list(core, list(U_1, .., U_N), error)`: the Tucker approximation refined by higher-order orthogonal iteration, with its relative Frobenius error |
//! | `cp_als(T, rank, sweeps, seed)` | `list(weights, list(A_1, .., A_N), residual)`: the CP decomposition `sum_r w_r a_r o b_r o ..` by alternating least squares from a seeded start, with the relative residual; a diverged run is not an answer |
//! | `tt_svd(T, tol, max_rank)` | `list(list(G_1, .., G_N), error)`: the tensor-train cores (`G_k[left][index][right]`) of `T` at relative tolerance `tol` (`max_rank` of `0` is no cap) and the truncation error |
//! | `tt_dense(cores)` | the dense tensor a list of tensor-train cores represents |
//!
//! The tensors are nested lists of any depth from 2 to 6; the number of
//! elements of each operand is capped so one call stays bounded. CP factors
//! are determined only up to permutation and scaling, so compare
//! reconstructions rather than factors.

use ndarray::ArrayD;

use super::bridge::Data;
use super::bridge::Entry;
use super::bridge::count;
use super::bridge::from_tensor;
use super::bridge::matrix;
use super::bridge::num;
use super::bridge::parts_of;
use super::bridge::scalar;
use super::bridge::tensor_arg;
use super::bridge::vec_out;
use super::bridge::vector;
use crate::kernels::dense::Mat;
use crate::kernels::lowrank;
use crate::kernels::tensor;
use crate::kernels::tensor_train;

/// Most elements of a decomposed matrix or tensor.
const MAX_ELEMENTS: usize = 200_000;
/// Most modes of a decomposed tensor.
const MAX_MODES: usize = 6;

/// The matrix argument as a kernel matrix.
fn kernel_matrix(
    a: &[Data],
    at: usize,
) -> Option<Mat> {
    let (rows, cols, data) = matrix(a, at)?;
    (rows * cols <= MAX_ELEMENTS).then(|| {
        Mat {
            rows,
            cols,
            data: data.to_vec(),
        }
    })
}

fn pivoted(a: &[Data]) -> Option<lowrank::PivotedQr> {
    lowrank::qr_pivoted(&kernel_matrix(a, 0)?).ok()
}

fn qr_pivoted_rank(a: &[Data]) -> Option<Data> {
    let tol = scalar(a, 1)?;
    (0.0..1.0).contains(&tol).then_some(())?;
    num(pivoted(a)?.rank(Some(tol)) as f64)
}

fn qr_pivoted_permutation(a: &[Data]) -> Option<Data> {
    vec_out(pivoted(a)?.perm.iter().map(|&j| j as f64).collect())
}

fn qr_pivoted_diagonal(a: &[Data]) -> Option<Data> {
    let r = pivoted(a)?.r;
    vec_out((0..r.rows.min(r.cols)).map(|i| r.at(i, i).abs()).collect())
}

/// A tensor operand of 2 to 6 modes with every extent positive.
fn dense_tensor(
    a: &[Data],
    at: usize,
) -> Option<ArrayD<f64>> {
    let t = tensor_arg(a, at)?;
    (t.ndim() >= 2
        && t.ndim() <= MAX_MODES
        && t.len() <= MAX_ELEMENTS
        && t.shape().iter().all(|&n| n >= 1))
    .then_some(t)
}

/// Multilinear ranks, one per mode; `0` stands for the full extent.
fn multilinear_ranks(
    a: &[Data],
    at: usize,
    shape: &[usize],
) -> Option<Vec<usize>> {
    let ranks = vector(a, at)?;
    (ranks.len() == shape.len()).then_some(())?;
    ranks
        .iter()
        .zip(shape)
        .map(|(&r, &n)| {
            (r >= 0.0 && r.fract() == 0.0 && r <= n as f64).then_some(if r == 0.0 {
                n
            } else {
                r as usize
            })
        })
        .collect()
}

fn matrix_data(m: &Mat) -> Data {
    Data::Array {
        shape: vec![m.rows, m.cols],
        data: m.data.clone(),
    }
}

fn factors_data(factors: &[Mat]) -> Data {
    Data::Tree(factors.iter().map(matrix_data).collect())
}

fn hosvd(a: &[Data]) -> Option<Data> {
    let t = dense_tensor(a, 0)?;
    let ranks = multilinear_ranks(a, 1, t.shape())?;
    let h = tensor::hosvd(&t, Some(&ranks)).ok()?;
    Some(Data::Tree(vec![
        from_tensor(&h.core)?,
        factors_data(&h.factors),
        Data::Tree(h.singular_values.into_iter().map(Data::vector).collect()),
    ]))
}

fn hooi(a: &[Data]) -> Option<Data> {
    let t = dense_tensor(a, 0)?;
    let ranks = multilinear_ranks(a, 1, t.shape())?;
    let sweeps = count(a, 2, 200)?;
    let tol = scalar(a, 3)?;
    (sweeps >= 1 && (0.0..1.0).contains(&tol)).then_some(())?;
    let fit = tensor::hooi(&t, &ranks, sweeps, tol).ok()?;
    Some(Data::Tree(vec![
        from_tensor(&fit.core)?,
        factors_data(&fit.factors),
        Data::Scalar(fit.relative_error),
    ]))
}

fn cp_als(a: &[Data]) -> Option<Data> {
    let t = dense_tensor(a, 0)?;
    let rank = count(a, 1, 24)?;
    let sweeps = count(a, 2, 1_000)?;
    let seed = scalar(a, 3)?;
    (rank >= 1 && sweeps >= 1 && seed >= 0.0 && seed.fract() == 0.0 && seed <= 9.0e15)
        .then_some(())?;
    let mut options = tensor::CpOptions::new(rank);
    options.max_sweeps = sweeps;
    options.seed = seed as u64;
    let fit = tensor::cp_als(&t, &options).ok()?;
    // A diverged run is the signature of a degenerate problem, not a fit.
    (!fit.diverged).then_some(())?;
    Some(Data::Tree(vec![
        Data::vector(fit.weights),
        factors_data(&fit.factors),
        Data::Scalar(fit.relative_residual),
    ]))
}

fn tt_svd(a: &[Data]) -> Option<Data> {
    let t = dense_tensor(a, 0)?;
    let tol = scalar(a, 1)?;
    let cap = count(a, 2, 1_000)?;
    (0.0..1.0).contains(&tol).then_some(())?;
    let data: Vec<f64> = t.iter().copied().collect();
    let train = tensor_train::tt_svd(&data, t.shape(), tol, (cap > 0).then_some(cap)).ok()?;
    let cores = train
        .cores()
        .iter()
        .map(|core| {
            Data::Array {
                shape: core.shape().to_vec(),
                data: core.data().to_vec(),
            }
        })
        .collect();
    Some(Data::Tree(vec![
        Data::Tree(cores),
        Data::Scalar(train.truncation_error()),
    ]))
}

fn tt_dense(a: &[Data]) -> Option<Data> {
    let parts = parts_of(a.first()?)?;
    (!parts.is_empty() && parts.len() <= MAX_MODES * 4).then_some(())?;
    let mut cores = Vec::with_capacity(parts.len());
    let mut stored = 0_usize;
    for part in &parts {
        let Data::Array { shape, data } = part else {
            return None;
        };
        let &[left, mode, right] = shape.as_slice() else {
            return None;
        };
        stored += data.len();
        (stored <= MAX_ELEMENTS).then_some(())?;
        cores.push(tensor_train::TrainCore::new([left, mode, right], data.clone()).ok()?);
    }
    let train = tensor_train::TensorTrain::from_cores(cores).ok()?;
    let dense = train.to_dense_with_limit(MAX_ELEMENTS).ok()?;
    Some(Data::Array {
        shape: train.shape().to_vec(),
        data: dense,
    })
}

pub(super) const ENTRIES: &[Entry] = &[
    Entry {
        name: "qr_pivoted_rank",
        params: "vs",
        eval: qr_pivoted_rank,
    },
    Entry {
        name: "qr_pivoted_permutation",
        params: "v",
        eval: qr_pivoted_permutation,
    },
    Entry {
        name: "qr_pivoted_diagonal",
        params: "v",
        eval: qr_pivoted_diagonal,
    },
    Entry {
        name: "hosvd",
        params: "nv",
        eval: hosvd,
    },
    Entry {
        name: "hooi",
        params: "nvss",
        eval: hooi,
    },
    Entry {
        name: "cp_als",
        params: "nsss",
        eval: cp_als,
    },
    Entry {
        name: "tt_svd",
        params: "nss",
        eval: tt_svd,
    },
    Entry {
        name: "tt_dense",
        params: "n",
        eval: tt_dense,
    },
];

#[cfg(test)]
mod tests {
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

    /// The first top-level element of a printed `list(...)`.
    fn first_item(text: &str) -> String {
        let inner = &text["list(".len()..];
        let mut depth = 0_i32;
        for (at, ch) in inner.char_indices() {
            match ch {
                | '(' => depth += 1,
                | ')' | ',' if depth == 0 => return inner[..at].to_owned(),
                | ')' => depth -= 1,
                | _ => {},
            }
        }
        inner.to_owned()
    }

    /// The outer product `a o b o c` of three vectors as a nested list.
    fn outer3(
        a: &[f64],
        b: &[f64],
        c: &[f64],
    ) -> String {
        let layer = |x: f64| {
            let rows: Vec<String> = b
                .iter()
                .map(|y| {
                    let cells: Vec<String> = c.iter().map(|z| format!("{}", x * y * z)).collect();
                    format!("list({})", cells.join(", "))
                })
                .collect();
            format!("list({})", rows.join(", "))
        };
        let layers: Vec<String> = a.iter().map(|&x| layer(x)).collect();
        format!("list({})", layers.join(", "))
    }

    #[test]
    fn pivoted_qr_reveals_the_rank_and_orders_the_columns() {
        // The second column is twice the first: rank one, and the longer
        // column is pivoted first.
        let m = "list(list(1, 2), list(2, 4), list(3, 6))";
        assert_eq!(run(&format!("qr_pivoted_rank({m}, 1e-10)")), "1");
        assert_eq!(
            flat(&run(&format!("qr_pivoted_permutation({m})"))),
            vec![1.0, 0.0]
        );
        let diagonal = flat(&run(&format!("qr_pivoted_diagonal({m})")));
        assert!(
            (diagonal[0] - 56.0_f64.sqrt()).abs() < 1e-12,
            "{diagonal:?}"
        );
        assert!(diagonal[1] < 1e-12, "{diagonal:?}");
        assert_eq!(
            run("qr_pivoted_rank(list(list(1, 0, 0), list(0, 1, 0), list(0, 0, 1)), 1e-12)"),
            "3"
        );
        // Orthogonal columns of growing length come out in reverse order.
        let growing = "list(list(1, 0, 0), list(0, 2, 0), list(0, 0, 3))";
        assert_eq!(
            flat(&run(&format!("qr_pivoted_permutation({growing})"))),
            vec![2.0, 1.0, 0.0]
        );
    }

    #[test]
    fn hosvd_of_a_rank_one_tensor_is_one_core_entry() {
        let t = outer3(&[1.0, 2.0], &[3.0, 0.0, 1.0], &[2.0, -1.0]);
        let norm = flat(&run(&format!("tensor_norm({t})")))[0];
        let out = flat(&run(&format!("hosvd({t}, list(0, 0, 0))")));
        // core (1) + factors (2 + 3 + 2) + singular values (1 + 1 + 1).
        assert_eq!(out.len(), 1 + 7 + 3, "{out:?}");
        assert!((out[0].abs() - norm).abs() < 1e-10, "{out:?} vs {norm}");
        assert!(out[8..].iter().all(|s| (s - norm).abs() < 1e-10), "{out:?}");
        // A full-rank tensor keeps its energy in the core: the factors are orthonormal.
        let full = "list(list(list(1, 2), list(3, 4)), list(list(5, 6), list(7, 9)))";
        let norm = flat(&run(&format!("tensor_norm({full})")))[0];
        let out = flat(&run(&format!("hosvd({full}, list(2, 2, 2))")));
        assert_eq!(out.len(), 8 + 12 + 6, "{out:?}");
        let energy: f64 = out[..8].iter().map(|v| v * v).sum();
        assert!((energy.sqrt() - norm).abs() < 1e-9, "{energy} vs {norm}");
        // Truncating a mode loses energy but never gains it.
        let cut = flat(&run(&format!("hosvd({full}, list(1, 2, 2))")));
        let kept: f64 = cut[..4].iter().map(|v| v * v).sum();
        assert!(kept <= energy + 1e-9);
    }

    #[test]
    fn hooi_fits_a_low_multilinear_rank_tensor_exactly() {
        let t = outer3(&[1.0, -1.0, 2.0], &[1.0, 3.0], &[0.5, 1.0, 2.0, -1.0]);
        let out = flat(&run(&format!("hooi({t}, list(1, 1, 1), 10, 1e-12)")));
        assert!(out.last().is_some_and(|e| e.abs() < 1e-9), "{out:?}");
    }

    #[test]
    fn cp_als_recovers_a_rank_two_tensor() {
        // 3 e1 o e1 o e1 + 2 e2 o e2 o e2.
        let t = "list(list(list(3, 0), list(0, 0)), list(list(0, 0), list(0, 2)))";
        let out = flat(&run(&format!("cp_als({t}, 2, 500, 1)")));
        // weights (2) + factors (3 x 4) + residual (1).
        assert_eq!(out.len(), 2 + 12 + 1, "{out:?}");
        let mut weights = out[..2].to_vec();
        weights.sort_by(f64::total_cmp);
        assert!(
            (weights[0] - 2.0).abs() < 1e-6 && (weights[1] - 3.0).abs() < 1e-6,
            "{weights:?}"
        );
        assert!(out[14] < 1e-8, "{out:?}");
        // The same seed gives the same run.
        assert_eq!(
            run(&format!("cp_als({t}, 2, 500, 1)")),
            run(&format!("cp_als({t}, 2, 500, 1)"))
        );
    }

    #[test]
    fn tensor_trains_round_trip_through_their_cores() {
        let t = outer3(&[1.0, 2.0], &[1.0, -1.0, 3.0], &[2.0, 1.0]);
        let compressed = run(&format!("tt_svd({t}, 1e-12, 0)"));
        let error = *flat(&compressed).last().unwrap_or(&f64::NAN);
        assert!(error.abs() < 1e-10, "{compressed}");
        let cores = first_item(&compressed);
        let rebuilt = flat(&run(&format!("tt_dense({cores})")));
        let want = flat(&t);
        assert_eq!(rebuilt.len(), want.len());
        for (g, w) in rebuilt.iter().zip(&want) {
            assert!((g - w).abs() < 1e-10, "{rebuilt:?} vs {want:?}");
        }
        // A rank cap of one is exact for a rank-one tensor, and the cores are then 1 x n x 1.
        let capped = run(&format!("tt_svd({t}, 0, 1)"));
        assert_eq!(flat(&capped).len(), 2 + 3 + 2 + 1);
    }

    #[test]
    fn invalid_decompositions_are_left_as_requests() {
        for src in [
            "qr_pivoted_rank(list(1, 2, 3), 0.1)",
            "qr_pivoted_rank(list(list(1, 2), list(3, 4)), 2)",
            "qr_pivoted_rank(list(list(1, a), list(3, 4)), 0.1)",
            "hosvd(list(list(1, 2), list(3, 4)), list(1))",
            "hosvd(list(list(1, 2), list(3, 4)), list(3, 1))",
            "hosvd(list(1, 2, 3), list(1))",
            "hosvd(list(list(0, 0), list(0, 0)), list(1, 1))",
            "hooi(list(list(1, 2), list(3, 4)), list(1, 1), 0, 0.1)",
            "cp_als(list(list(1, 2), list(3, 4)), 0, 10, 1)",
            "cp_als(list(1, 2, 3), 1, 10, 1)",
            "tt_svd(list(list(1, 2), list(3, 4)), 1.5, 0)",
            "tt_dense(list(list(1, 2)))",
            "tt_dense(list(list(list(list(1, 2))), list(list(list(1)))))",
        ] {
            let (text, reduced) = reduce_with(&[numeric()], src, &[]);
            assert!(!reduced, "{text}");
        }
    }
}
