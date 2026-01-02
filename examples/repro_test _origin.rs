use rssn::prelude::numerical::*;

fn main() {
    let rows = 4;
    let cols = 4;
    let seed = 650u64;
    
    let mut data = Vec::with_capacity(rows * cols);
    let mut rng = seed;
    for _ in 0..(rows*cols) {
        rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
        let val = (rng % 100) as f64 / 10.0 - 5.0; // -5.0 to 5.0
        data.push(val);
    }

    let a = numerical_Matrix::new(rows, cols, data.clone());

    // Generate random solution x
    let mut x = Vec::with_capacity(rows);
    for _ in 0..rows {
        rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
        x.push((rng % 100) as f64 / 10.0 - 5.0);
    }

    // Compute b = A * x
    let mut b = vec![0.0; rows];
    for i in 0..rows {
        let mut sum = 0.0;
        for j in 0..cols {
            sum += a.get(i, j) * x[j];
        }
        b[i] = sum;
    }

    println!("A = {:?}", data);
    println!("b = {:?}", b);

    // Solve Ax = b
    let result = numerical_solve_linear_system(&a, &b).unwrap();

    match result {
        numerical_LinearSolution::Unique(sol) => {
            println!("Solution found: {:?}", sol);
            for i in 0..rows {
                let mut sum = 0.0;
                for j in 0..cols {
                    sum += a.get(i, j) * sol[j];
                }
                println!("Row {}: A*sol = {}, b = {}, diff = {}", i, sum, b[i], (sum - b[i]).abs());
            }
        },
        _ => println!("Non-unique solution"),
    }
}
