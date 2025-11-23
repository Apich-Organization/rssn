use rssn::symbolic::core::Expr;
use rssn::symbolic::matrix::*;
use rssn::symbolic::simplify_dag::simplify;

fn main() {
    println!("Symbolic Matrix Demo");
    println!("====================");

    // Create matrices
    let m1 = Expr::Matrix(vec![
        vec![Expr::new_constant(1.0), Expr::new_constant(2.0)],
        vec![Expr::new_constant(3.0), Expr::new_constant(4.0)],
    ]);
    let m2 = Expr::Matrix(vec![
        vec![Expr::new_constant(2.0), Expr::new_constant(0.0)],
        vec![Expr::new_constant(1.0), Expr::new_constant(2.0)],
    ]);

    println!("Matrix 1:");
    println!("{:?}", m1);
    println!("Matrix 2:");
    println!("{:?}", m2);

    // Addition
    let sum = add_matrices(&m1, &m2);
    println!("Sum (M1 + M2):");
    println!("{:?}", sum);

    // Multiplication
    let prod = mul_matrices(&m1, &m2);
    println!("Product (M1 * M2):");
    println!("{:?}", prod);

    // Transpose
    let t = transpose_matrix(&m1);
    println!("Transpose (M1^T):");
    println!("{:?}", t);

    // Determinant
    let det = determinant(&m1);
    println!("Determinant (det(M1)):");
    println!("{:?}", det);

    // Inverse
    let inv = inverse_matrix(&m1);
    println!("Inverse (M1^-1):");
    println!("{:?}", inv);

    // Solve linear system Ax = b
    // 4x + 7y = 5
    // 2x + 6y = -2
    let a = Expr::Matrix(vec![
        vec![Expr::new_constant(4.0), Expr::new_constant(7.0)],
        vec![Expr::new_constant(2.0), Expr::new_constant(6.0)],
    ]);
    let b = Expr::Matrix(vec![
        vec![Expr::new_constant(5.0)],
        vec![Expr::new_constant(-2.0)],
    ]);
    println!("Solving Ax = b:");
    println!("A = {:?}", a);
    println!("b = {:?}", b);
    match solve_linear_system(&a, &b) {
        Ok(sol) => println!("Solution x = {:?}", sol),
        Err(e) => println!("Error: {}", e),
    }
}
