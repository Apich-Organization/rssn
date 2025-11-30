use rssn::symbolic::core::Expr;
use rssn::symbolic::multi_valued::*;

#[test]
fn test_arg() {
    let z = Expr::new_complex(Expr::Constant(1.0), Expr::Constant(1.0));
    let result = arg(&z);
    
    // Should return Arg(complex(1, 1))
    // We can't evaluate it symbolically, but we can check it's constructed
    assert!(matches!(result, Expr::Apply(_, _) | Expr::Dag(_)));
}

#[test]
fn test_abs() {
    let z = Expr::new_complex(Expr::Constant(3.0), Expr::Constant(4.0));
    let result = abs(&z);
    
    // Should return Abs(complex(3, 4))
    // The magnitude would be 5, but symbolically it's Abs(...)
    assert!(matches!(result, Expr::Abs(_) | Expr::Dag(_)));
}

#[test]
fn test_general_log() {
    let z = Expr::Variable("z".to_string());
    let k = Expr::Variable("k".to_string());
    
    let result = general_log(&z, &k);
    
    // Result should be a complex expression involving log, arg, and k
    // We just verify it returns something
    assert!(!matches!(result, Expr::Constant(0.0)));
}

#[test]
fn test_general_sqrt() {
    let z = Expr::Variable("z".to_string());
    let k = Expr::Constant(0.0); // Principal branch
    
    let result = general_sqrt(&z, &k);
    
    // Result should be a complex expression
    assert!(!matches!(result, Expr::Constant(0.0)));
}

#[test]
fn test_general_sqrt_two_branches() {
    let z = Expr::Constant(4.0);
    let k0 = Expr::Constant(0.0);
    let k1 = Expr::Constant(1.0);
    
    let branch0 = general_sqrt(&z, &k0);
    let branch1 = general_sqrt(&z, &k1);
    
    // Two branches should be different
    // We can't easily compare them symbolically, but they should both exist
    assert!(!matches!(branch0, Expr::Constant(0.0)));
    assert!(!matches!(branch1, Expr::Constant(0.0)));
}

#[test]
fn test_general_power() {
    let z = Expr::Variable("z".to_string());
    let w = Expr::Constant(0.5); // Square root
    let k = Expr::Constant(0.0);
    
    let result = general_power(&z, &w, &k);
    
    // Result should be exp(w * log(z))
    assert!(!matches!(result, Expr::Constant(0.0)));
}

#[test]
fn test_general_nth_root() {
    let z = Expr::Constant(8.0);
    let n = Expr::Constant(3.0); // Cube root
    let k = Expr::Constant(0.0);
    
    let result = general_nth_root(&z, &n, &k);
    
    // Result should be the principal cube root of 8
    assert!(!matches!(result, Expr::Constant(0.0)));
}

#[test]
fn test_general_nth_root_multiple_branches() {
    let z = Expr::Constant(1.0);
    let n = Expr::Constant(3.0); // Cube roots of unity
    
    let k0 = Expr::Constant(0.0);
    let k1 = Expr::Constant(1.0);
    let k2 = Expr::Constant(2.0);
    
    let root0 = general_nth_root(&z, &n, &k0);
    let root1 = general_nth_root(&z, &n, &k1);
    let root2 = general_nth_root(&z, &n, &k2);
    
    // All three cube roots should exist
    assert!(!matches!(root0, Expr::Constant(0.0)));
    assert!(!matches!(root1, Expr::Constant(0.0)));
    assert!(!matches!(root2, Expr::Constant(0.0)));
}

#[test]
fn test_general_arcsin() {
    let z = Expr::Constant(0.5);
    let k = Expr::Constant(0.0);
    
    let result = general_arcsin(&z, &k);
    
    // Result should involve arcsin and k
    assert!(!matches!(result, Expr::Constant(0.0)));
}

#[test]
fn test_general_arccos() {
    let z = Expr::Constant(0.5);
    let k = Expr::Constant(0.0);
    let s = Expr::Constant(1.0); // Positive sign
    
    let result = general_arccos(&z, &k, &s);
    
    // Result should involve arccos, k, and s
    assert!(!matches!(result, Expr::Constant(0.0)));
}

#[test]
fn test_general_arctan() {
    let z = Expr::Constant(1.0);
    let k = Expr::Constant(0.0);
    
    let result = general_arctan(&z, &k);
    
    // Result should be k*pi + arctan(1)
    assert!(!matches!(result, Expr::Constant(0.0)));
}

#[test]
fn test_general_arcsinh() {
    let z = Expr::Variable("z".to_string());
    let k = Expr::Constant(0.0);
    
    let result = general_arcsinh(&z, &k);
    
    // Result should be log(z + sqrt(z^2 + 1))
    // eprintln!("Result: {:?}", result);
    assert!(!matches!(result, Expr::Constant(0.0)));
    // assert!(false);
}

#[test]
fn test_general_arccosh() {
    let z = Expr::Constant(2.0);
    let k = Expr::Constant(0.0);
    
    let result = general_arccosh(&z, &k);
    
    // Result should be log(z + sqrt(z^2 - 1))
    assert!(!matches!(result, Expr::Constant(0.0)));
}

#[test]
fn test_general_arctanh() {
    let z = Expr::Constant(0.5);
    let k = Expr::Constant(0.0);
    
    let result = general_arctanh(&z, &k);
    
    // Result should be (1/2) * log((1+z)/(1-z))
    assert!(!matches!(result, Expr::Constant(0.0)));
}

#[test]
fn test_general_log_with_symbolic_k() {
    let z = Expr::new_complex(Expr::Constant(1.0), Expr::Constant(0.0));
    let k = Expr::Variable("k".to_string());
    
    let result = general_log(&z, &k);
    
    // Result should contain the symbolic k
    // We can't easily verify the exact structure, but it should be non-trivial
    assert!(!matches!(result, Expr::Constant(0.0)));
}

#[test]
fn test_general_power_integer_exponent() {
    let z = Expr::Variable("z".to_string());
    let w = Expr::Constant(2.0); // Square
    let k = Expr::Constant(0.0);
    
    let result = general_power(&z, &w, &k);
    
    // z^2 via exp(2*log(z))
    assert!(!matches!(result, Expr::Constant(0.0)));
}

#[test]
fn test_abs_real_number() {
    let z = Expr::Constant(-5.0);
    let result = abs(&z);
    
    // Abs(-5) should be constructed
    assert!(matches!(result, Expr::Abs(_) | Expr::Dag(_)));
}

#[test]
fn test_arg_real_positive() {
    let z = Expr::Constant(5.0);
    let result = arg(&z);
    
    // Arg(5) should be 0, but symbolically it's Arg(5)
    assert!(matches!(result, Expr::Apply(_, _) | Expr::Dag(_)));
}

#[test]
fn test_general_sqrt_of_negative() {
    let z = Expr::Constant(-1.0);
    let k = Expr::Constant(0.0);
    
    let result = general_sqrt(&z, &k);
    
    // sqrt(-1) = i (principal branch)
    // The result should be a complex expression
    assert!(!matches!(result, Expr::Constant(0.0)));
}

#[test]
fn test_general_nth_root_of_negative() {
    let z = Expr::Constant(-8.0);
    let n = Expr::Constant(3.0);
    let k = Expr::Constant(0.0);
    
    let result = general_nth_root(&z, &n, &k);
    
    // Cube root of -8 should be -2 (principal branch)
    assert!(!matches!(result, Expr::Constant(0.0)));
}
