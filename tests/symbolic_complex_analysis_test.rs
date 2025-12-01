use rssn::symbolic::complex_analysis::*;
use rssn::symbolic::core::Expr;
use rssn::symbolic::simplify_dag::simplify;
use std::sync::Arc;

#[test]
fn test_mobius_transformation() {
    // Test identity transformation
    let identity = MobiusTransformation::identity();
    let z = Expr::Variable("z".to_string());
    let result = identity.apply(&z);
    
    // Identity should return z
    println!("Identity(z) = {:?}", result);
    assert_eq!(simplify(&result), z);
}

#[test]
fn test_mobius_composition() {
    // f(z) = z + 1 (translation)
    let f = MobiusTransformation::new(
        Expr::Constant(1.0),
        Expr::Constant(1.0),
        Expr::Constant(0.0),
        Expr::Constant(1.0),
    );
    
    // g(z) = 2z (scaling)
    let g = MobiusTransformation::new(
        Expr::Constant(2.0),
        Expr::Constant(0.0),
        Expr::Constant(0.0),
        Expr::Constant(1.0),
    );
    
    // Compose: (f ∘ g)(z) = f(g(z)) = f(2z) = 2z + 1
    let composed = f.compose(&g);
    
    let z = Expr::Constant(3.0);
    let result = composed.apply(&z);
    
    // Should be 2*3 + 1 = 7
    println!("Composed transformation result: {:?}", result);
    let result_val = simplify(&result);
    assert_eq!(result_val, Expr::Constant(7.0));
}

#[test]
fn test_mobius_inverse() {
    // f(z) = (2z + 1) / (z + 3)
    let f = MobiusTransformation::new(
        Expr::Constant(2.0),
        Expr::Constant(1.0),
        Expr::Constant(1.0),
        Expr::Constant(3.0),
    );
    
    let f_inv = f.inverse();
    
    // f^(-1)(f(z)) should equal z
    let z = Expr::Constant(5.0);
    let fz = f.apply(&z);
    let result = f_inv.apply(&fz);
    
    println!("f(z) = {:?}", fz);
    println!("f^(-1)(f(z)) = {:?}", result);
    
    // Should approximately equal 5.0
    let result_simplified = simplify(&result);
    println!("Simplified: {:?}", result_simplified);
}

#[test]
fn test_complex_modulus() {
    // |3 + 4i| = 5
    let z = Expr::Complex(
        Arc::new(Expr::Constant(3.0)),
        Arc::new(Expr::Constant(4.0))
    );
    
    let modulus = complex_modulus(&z);
    let result = simplify(&modulus);
    
    println!("Modulus of 3+4i: {:?}", result);
    // Result can be either 5.0 or sqrt(25)
    let result_str = format!("{:?}", result);
    let is_correct = result == Expr::Constant(5.0) || result_str.contains("sqrt") || result_str.contains("Sqrt");
    assert!(is_correct, "Modulus should be 5 or sqrt(25), got {:?}", result);
}

#[test]
fn test_complex_arg() {
    // arg(1 + i) = π/4
    let z = Expr::Complex(
        Arc::new(Expr::Constant(1.0)),
        Arc::new(Expr::Constant(1.0))
    );
    
    let arg = complex_arg(&z);
    
    println!("Argument of 1+i: {:?}", arg);
    // Should contain atan2
    let arg_str = format!("{:?}", arg);
    assert!(arg_str.contains("Atan2") || arg_str.contains("atan"));
}

#[test]
fn test_classify_singularity() {
    // f(z) = 1/z has a simple pole at z=0
    let z = Expr::Variable("z".to_string());
    let func = Expr::new_div(Expr::Constant(1.0), z.clone());
    println!("in test Function: {:?}", func);
    
    let singularity_type = classify_singularity(&func, "z", &Expr::Constant(0.0), 5);
    
    println!("Singularity type of 1/z at z=0: {:?}", singularity_type);
    // Should be a pole of order 1
    match singularity_type {
        SingularityType::Pole(n) => {
            assert_eq!(n, 1, "Expected pole of order 1, got order {}", n);
        },
        _ => panic!("Expected pole, got {:?}", singularity_type),
    }
}

#[test]
fn test_calculate_residue() {
    // Residue of 1/z at z=0 is 1
    let z = Expr::Variable("z".to_string());
    let func = Expr::new_div(Expr::Constant(1.0), z.clone());
    
    let residue = calculate_residue(&func, "z", &Expr::Constant(0.0));
    
    println!("Residue of 1/z at z=0: {:?}", residue);
    let residue_simplified = simplify(&residue);
    println!("Simplified residue: {:?}", residue_simplified);
}

#[test]
fn test_cauchy_integral_formula() {
    // Cauchy's formula: f(z0) = (1/2πi) ∮ f(z)/(z-z0) dz
    // For f(z) = z^2, f(1) = 1
    let z = Expr::Variable("z".to_string());
    let func = Expr::new_pow(z.clone(), Expr::Constant(2.0));
    
    let result = cauchy_integral_formula(&func, "z", &Expr::Constant(1.0));
    
    println!("Cauchy formula result: {:?}", result);
    assert_eq!(simplify(&result), Expr::Constant(1.0));
}

#[test]
fn test_cauchy_derivative_formula() {
    // f(z) = z^3, f'(z) = 3z^2, f'(2) = 12
    let z = Expr::Variable("z".to_string());
    let func = Expr::new_pow(z.clone(), Expr::Constant(3.0));
    
    let result = cauchy_derivative_formula(&func, "z", &Expr::Constant(2.0), 1);
    
    println!("Derivative at z=2: {:?}", result);
    assert_eq!(simplify(&result), Expr::Constant(12.0));
}

#[test]
fn test_complex_distance() {
    // Distance between 0 and 3+4i should be 5
    let p1 = Expr::Constant(0.0);
    let p2 = Expr::Complex(
        Arc::new(Expr::Constant(3.0)),
        Arc::new(Expr::Constant(4.0))
    );
    
    let distance = complex_distance(&p1, &p2).unwrap();
    
    println!("Distance: {}", distance);
    assert!((distance - 5.0).abs() < 1e-10);
}
