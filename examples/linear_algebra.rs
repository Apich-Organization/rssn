//! Matrices two ways: exactly, through the symbolic `linalg` operators, and
//! in floating point, through `kernels::dense`.
//!
//! Run with `cargo run --release --example linear_algebra`.

use rssn::api::*;
use rssn::kernels::dense::Mat;
use rssn::kernels::dense::{
    self,
};

/// Reduces one request with the standard rules and prints the answer.
fn ask(
    s: &Session,
    request: &str,
) -> Result<String, Box<dyn std::error::Error>> {
    let term = s.parse(request)?;
    Ok(s.compute(term, &Config::new())?.term.to_string())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let s = Session::new();

    // Exact linear algebra: entries may be symbols, and the answer is a
    // closed form in them.
    println!(
        "det [[a, b], [c, d]]     = {}",
        ask(&s, "det(list(list(a, b), list(c, d)))")?
    );
    println!(
        "inverse [[1, 2], [3, 4]]  = {}",
        ask(&s, "inverse(list(list(1, 2), list(3, 4)))")?
    );
    // A^-1 A is the identity, which is the cheapest proof that the inverse
    // is right — and it is a symbolic identity, not a numerical one.
    println!(
        "inverse(A) A             = {}",
        ask(
            &s,
            "matmul(inverse(list(list(2, 1, 1), list(1, 3, 2), list(1, 0, 0))), list(list(2, 1, 1), list(1, 3, 2), list(1, 0, 0)))"
        )?
    );
    // A linear solve keeps rationals exact: 2x + y = 3, x + 3y = 5.
    println!(
        "linsolve [[2,1],[1,3]]   = {}",
        ask(&s, "linsolve(list(list(2, 1), list(1, 3)), list(3, 5))")?
    );
    println!(
        "eigenvals [[2,1],[1,2]]  = {}",
        ask(&s, "eigenvals(list(list(2, 1), list(1, 2)))")?
    );
    // A singular system still has an answer: a line, with a free parameter.
    println!(
        "singular system          = {}",
        ask(&s, "linsolve(list(list(1, 1, 1)), list(3))")?
    );

    // Floating point: the same problems through the dense kernels.
    let a = Mat::from_rows(&[
        vec![2.0, -1.0, 0.0],
        vec![-1.0, 2.0, -1.0],
        vec![0.0, -1.0, 2.0],
    ])?;
    let b = [1.0, 0.0, 3.0];

    let x = dense::solve(&a, &b)?;
    let residual = (0..3)
        .map(|i| (0..3).map(|j| a.at(i, j) * x[j]).sum::<f64>() - b[i])
        .fold(0.0_f64, |m, r| m.max(r.abs()));
    println!("\ndense solve: x = [{:.6}, {:.6}, {:.6}]", x[0], x[1], x[2]);
    println!("  max |A x - b| = {residual:.3e}");

    println!("  det (by LU)   = {:.6}", dense::lu_factor(&a)?.det());
    // The exact determinant of this matrix is 4; eigenvalues 2 - sqrt 2,
    // 2 and 2 + sqrt 2, whose product is 4.
    let (values, _vectors) = dense::eigen_symmetric(&a)?;
    println!(
        "  eigenvalues   = [{:.6}, {:.6}, {:.6}]",
        values[0], values[1], values[2]
    );
    println!(
        "  predicted     = [{:.6}, {:.6}, {:.6}]",
        2.0 - 2.0_f64.sqrt(),
        2.0,
        2.0 + 2.0_f64.sqrt()
    );
    Ok(())
}
