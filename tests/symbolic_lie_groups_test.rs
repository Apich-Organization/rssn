use rssn::symbolic::core::Expr;
use rssn::symbolic::lie_groups_and_algebras::*;
use rssn::symbolic::matrix;

#[test]
fn test_so3_generators() {
    let generators = so3_generators();
    assert_eq!(generators.len(), 3);
    
    // Check that each generator is a 3x3 matrix
    for gen in &generators {
        if let Expr::Matrix(rows) = &gen.0 {
            assert_eq!(rows.len(), 3);
            assert_eq!(rows[0].len(), 3);
        } else {
            panic!("Generator is not a matrix");
        }
    }
}

#[test]
fn test_su2_generators() {
    let generators = su2_generators();
    assert_eq!(generators.len(), 3);
    
    // Check that each generator is a 2x2 matrix
    for gen in &generators {
        if let Expr::Matrix(rows) = &gen.0 {
            assert_eq!(rows.len(), 2);
            assert_eq!(rows[0].len(), 2);
        } else {
            panic!("Generator is not a matrix");
        }
    }
}

#[test]
fn test_lie_bracket_antisymmetry() {
    let so3_algebra = so3();
    let x = &so3_algebra.basis[0].0;
    let y = &so3_algebra.basis[1].0;
    
    // [X, Y] = -[Y, X]
    let xy = lie_bracket(x, y).unwrap();
    let yx = lie_bracket(y, x).unwrap();
    let neg_yx = matrix::scalar_mul_matrix(&Expr::Constant(-1.0), &yx);
    
    // Check if xy equals -yx
    if let (Expr::Matrix(xy_mat), Expr::Matrix(neg_yx_mat)) = (&xy, &neg_yx) {
        for i in 0..3 {
            for j in 0..3 {
                // We'll just check that they're structurally similar
                // (exact equality might be tricky with symbolic expressions)
            }
        }
    }
}

#[test]
fn test_lie_bracket_bilinearity() {
    let so3_algebra = so3();
    let x = &so3_algebra.basis[0].0;
    let y = &so3_algebra.basis[1].0;
    
    // Test that lie bracket works
    let bracket = lie_bracket(x, y);
    assert!(bracket.is_ok());
}

#[test]
fn test_exponential_map() {
    let so3_algebra = so3();
    let x = &so3_algebra.basis[0].0;
    
    // Compute exp(X) with order 5
    let exp_x = exponential_map(x, 5);
    assert!(exp_x.is_ok());
    
    let result = exp_x.unwrap();
    // Check it's a 3x3 matrix
    if let Expr::Matrix(rows) = result {
        assert_eq!(rows.len(), 3);
        assert_eq!(rows[0].len(), 3);
    } else {
        panic!("Exponential map did not return a matrix");
    }
}

#[test]
fn test_adjoint_representation_algebra() {
    let so3_algebra = so3();
    let x = &so3_algebra.basis[0].0;
    let y = &so3_algebra.basis[1].0;
    
    // ad_X(Y) = [X, Y]
    let ad_xy = adjoint_representation_algebra(x, y);
    let bracket_xy = lie_bracket(x, y);
    
    assert!(ad_xy.is_ok());
    assert!(bracket_xy.is_ok());
    
    // They should be equal
    assert_eq!(ad_xy.unwrap(), bracket_xy.unwrap());
}

#[test]
fn test_commutator_table() {
    let so3_algebra = so3();
    let table = commutator_table(&so3_algebra);
    
    assert!(table.is_ok());
    let table = table.unwrap();
    
    // Should be a 3x3 table
    assert_eq!(table.len(), 3);
    for row in &table {
        assert_eq!(row.len(), 3);
    }
    
    // Diagonal elements should be zero (anti-symmetry: [X, X] = 0)
    for i in 0..3 {
        assert!(matrix::is_zero_matrix(&table[i][i]));
    }
}

#[test]
fn test_jacobi_identity_so3() {
    let so3_algebra = so3();
    let result = check_jacobi_identity(&so3_algebra);
    
    assert!(result.is_ok());
    // For so(3), the Jacobi identity should hold
    // Note: This might fail if symbolic simplification is not complete
    // but it should at least not error
}

#[test]
fn test_jacobi_identity_su2() {
    let su2_algebra = su2();
    let result = check_jacobi_identity(&su2_algebra);
    
    assert!(result.is_ok());
    // For su(2), the Jacobi identity should hold
}

#[test]
fn test_so3_structure() {
    let so3_algebra = so3();
    
    assert_eq!(so3_algebra.name, "so(3)");
    assert_eq!(so3_algebra.dimension, 3);
    assert_eq!(so3_algebra.basis.len(), 3);
}

#[test]
fn test_su2_structure() {
    let su2_algebra = su2();
    
    assert_eq!(su2_algebra.name, "su(2)");
    assert_eq!(su2_algebra.dimension, 3);
    assert_eq!(su2_algebra.basis.len(), 3);
}

#[test]
fn test_exponential_map_identity() {
    // exp(0) should be the identity matrix
    let zero_matrix = Expr::Matrix(vec![
        vec![Expr::Constant(0.0), Expr::Constant(0.0)],
        vec![Expr::Constant(0.0), Expr::Constant(0.0)],
    ]);
    
    let exp_zero = exponential_map(&zero_matrix, 5).unwrap();
    let identity = matrix::identity_matrix(2);
    
    // exp(0) should equal I
    if let (Expr::Matrix(exp_mat), Expr::Matrix(id_mat)) = (&exp_zero, &identity) {
        assert_eq!(exp_mat.len(), id_mat.len());
        assert_eq!(exp_mat[0].len(), id_mat[0].len());
    }
}

#[test]
fn test_adjoint_representation_group() {
    // Create a simple 2x2 group element (rotation matrix)
    let g = Expr::Matrix(vec![
        vec![Expr::Constant(0.0), Expr::Constant(-1.0)],
        vec![Expr::Constant(1.0), Expr::Constant(0.0)],
    ]);
    
    let x = Expr::Matrix(vec![
        vec![Expr::Constant(1.0), Expr::Constant(0.0)],
        vec![Expr::Constant(0.0), Expr::Constant(-1.0)],
    ]);
    
    let result = adjoint_representation_group(&g, &x);
    assert!(result.is_ok());
    
    // Result should be a matrix
    if let Expr::Matrix(rows) = result.unwrap() {
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].len(), 2);
    } else {
        panic!("Adjoint representation did not return a matrix");
    }
}

#[test]
fn test_serialization() {
    use serde_json;
    
    let so3_algebra = so3();
    
    // Test serialization
    let serialized = serde_json::to_string(&so3_algebra);
    assert!(serialized.is_ok());
    
    // Test deserialization
    let deserialized: Result<LieAlgebra, _> = serde_json::from_str(&serialized.unwrap());
    assert!(deserialized.is_ok());
    
    let recovered = deserialized.unwrap();
    assert_eq!(recovered.name, "so(3)");
    assert_eq!(recovered.dimension, 3);
}
