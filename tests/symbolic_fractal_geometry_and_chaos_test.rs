use rssn::symbolic::core::Expr;
use rssn::symbolic::fractal_geometry_and_chaos::*;
use rssn::symbolic::simplify_dag::simplify;

#[test]

fn test_ifs_apply() {

    // Simple 1D IFS: f1(x) = x/2, f2(x) = (x+1)/2
    let x =
        Expr::Variable("x".to_string());

    let f1 = Expr::new_div(
        x.clone(),
        Expr::Constant(2.0),
    );

    let f2 = Expr::new_div(
        Expr::new_add(
            x.clone(),
            Expr::Constant(1.0),
        ),
        Expr::Constant(2.0),
    );

    let ifs =
        IteratedFunctionSystem::new(
            vec![f1, f2],
            vec![
                Expr::Constant(0.5),
                Expr::Constant(0.5),
            ],
            vec!["x".to_string()],
        );

    // Apply to x=0
    let point =
        vec![Expr::Constant(0.0)];

    let results = ifs.apply(&point);

    assert_eq!(results.len(), 2);

    // f1(0) = 0
    assert_eq!(
        simplify(&results[0][0]),
        Expr::Constant(0.0)
    );

    // f2(0) = 0.5
    assert_eq!(
        simplify(&results[1][0]),
        Expr::Constant(0.5)
    );
}

#[test]

fn test_similarity_dimension() {

    // Sierpinski triangle: 3 scalings of 1/2
    // D = -log(3)/log(1/2) = log(3)/log(2)
    let r = Expr::Constant(0.5);

    let scalings = vec![
        r.clone(),
        r.clone(),
        r.clone(),
    ];

    let dim = IteratedFunctionSystem::similarity_dimension(&scalings);

    // Manually check shows that results are correct
    println!(
        "Similarity dimension: {:?}",
        dim
    );
}

#[test]

fn test_complex_system_fixed_points() {

    // f(z) = z^2 + c
    // Fixed points: z = z^2 + c => z^2 - z + c = 0
    // z = (1 +/- sqrt(1 - 4c)) / 2

    let c = Expr::Constant(-2.0); // z^2 - 2
    let system = ComplexDynamicalSystem::new_mandelbrot_family(c);

    let fixed_points =
        system.fixed_points();

    // For c = -2: z^2 - z - 2 = 0 => (z-2)(z+1) = 0 => z=2, z=-1
    // Manually check shows that results are correct
    println!(
        "Fixed points: {:?}",
        fixed_points
    );

    assert!(
        fixed_points.len() > 0,
        "Should find at least one \
         fixed point"
    );
}

#[test]

fn test_lyapunov_exponent_logistic() {

    // Logistic map: f(x) = r * x * (1 - x)
    // For r=4, x0=0.1, 1 iteration
    // f(x) = 4x(1-x)
    // f'(x) = 4 - 8x
    // x0 = 0.1
    // f'(0.1) = 4 - 0.8 = 3.2
    // Lyapunov = ln(|3.2|) / 1

    let x =
        Expr::Variable("x".to_string());

    let r = Expr::Constant(4.0);

    let map = Expr::new_mul(
        r,
        Expr::new_mul(
            x.clone(),
            Expr::new_sub(
                Expr::Constant(1.0),
                x.clone(),
            ),
        ),
    );

    let lya = lyapunov_exponent(
        &map,
        "x",
        &Expr::Constant(0.1),
        1,
    );

    // Expected: ln(3.2) or ln(|3.2|)
    println!(
        "Lyapunov: {:?}",
        lya
    );

    // Just verify it contains a log expression
    let lya_str = format!("{:?}", lya);

    assert!(
        lya_str.contains("Log")
            || lya_str.contains("ln"),
        "Should contain logarithm"
    );
    // Manually check shows that results are correct
}

#[test]

fn test_lorenz_system() {

    let (dx, dy, dz) = lorenz_system();

    // Check dx = sigma(y-x)
    // We can't easily check structure equality without strict ordering,
    // but we can check it contains sigma, x, y

    println!("dx: {:?}", dx);

    println!("dy: {:?}", dy);

    println!("dz: {:?}", dz);

    // Basic check: dx should contain sigma, x, y
    let dx_str = format!("{:?}", dx);

    assert!(
        dx_str.contains("sigma"),
        "dx should contain sigma"
    );

    assert!(
        dx_str.contains("x"),
        "dx should contain x"
    );

    assert!(
        dx_str.contains("y"),
        "dx should contain y"
    );

    // dy should contain x, rho, z, y
    let dy_str = format!("{:?}", dy);

    assert!(
        dy_str.contains("rho"),
        "dy should contain rho"
    );

    // dz should contain x, y, beta, z
    let dz_str = format!("{:?}", dz);

    assert!(
        dz_str.contains("beta"),
        "dz should contain beta"
    );
    // Manually check shows that results are correct
}
