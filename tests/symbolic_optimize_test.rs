use rssn::symbolic::core::Expr;
use rssn::symbolic::optimize::*;
use rssn::symbolic::simplify_dag::simplify;

#[test]
fn test_find_extrema_1d() {
    // f(x) = x^2
    let x = Expr::Variable("x".to_string());
    let f = Expr::new_pow(x.clone(), Expr::Constant(2.0));
    
    let extrema = find_extrema(&f, &["x"]).unwrap();
    
    assert_eq!(extrema.len(), 1);
    let point = &extrema[0];
    assert_eq!(point.point_type, ExtremumType::LocalMin);
    
    // Check if x is 0
    let x_val = point.point.get(&x).unwrap();
    let simplified_x = simplify(x_val);
    assert!(rssn::symbolic::simplify::is_zero(&simplified_x));
}

#[test]
fn test_find_extrema_2d_min() {
    // f(x, y) = x^2 + y^2
    let x = Expr::Variable("x".to_string());
    let y = Expr::Variable("y".to_string());
    let f = Expr::new_add(
        Expr::new_pow(x.clone(), Expr::Constant(2.0)),
        Expr::new_pow(y.clone(), Expr::Constant(2.0))
    );
    
    let extrema = find_extrema(&f, &["x", "y"]).unwrap();
    
    assert_eq!(extrema.len(), 1);
    let point = &extrema[0];
    assert_eq!(point.point_type, ExtremumType::LocalMin);
    
    // Just check that we got solutions for both variables
    assert!(point.point.contains_key(&x));
    assert!(point.point.contains_key(&y));
}

#[test]
fn test_find_extrema_2d_saddle() {
    // f(x, y) = x^2 - y^2
    let x = Expr::Variable("x".to_string());
    let y = Expr::Variable("y".to_string());
    let f = Expr::new_sub(
        Expr::new_pow(x.clone(), Expr::Constant(2.0)),
        Expr::new_pow(y.clone(), Expr::Constant(2.0))
    );
    
    let extrema = find_extrema(&f, &["x", "y"]).unwrap();
    
    assert_eq!(extrema.len(), 1);
    let point = &extrema[0];
    assert_eq!(point.point_type, ExtremumType::SaddlePoint);
    
    // Just check that we got solutions for both variables
    assert!(point.point.contains_key(&x));
    assert!(point.point.contains_key(&y));
}

#[test]
fn test_hessian_matrix() {
    // f(x, y) = x^2 + y^2
    let x = Expr::Variable("x".to_string());
    let y = Expr::Variable("y".to_string());
    let f = Expr::new_add(
        Expr::new_pow(x.clone(), Expr::Constant(2.0)),
        Expr::new_pow(y.clone(), Expr::Constant(2.0))
    );
    
    let hessian = hessian_matrix(&f, &["x", "y"]);
    
    // Check that we got a matrix
    if let Expr::Matrix(rows) = hessian {
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].len(), 2);
    } else {
        panic!("Expected matrix");
    }
}
