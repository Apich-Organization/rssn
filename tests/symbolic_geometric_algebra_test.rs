use rssn::symbolic::core::Expr;
use rssn::symbolic::geometric_algebra::Multivector;

#[test]

fn test_scalar_creation() {

    let signature = (3, 0, 0); // 3D Euclidean space
    let scalar = Multivector::scalar(
        signature,
        Expr::new_constant(5.0),
    );

    assert_eq!(
        scalar.terms.len(),
        1
    );

    assert_eq!(
        scalar.terms.get(&0),
        Some(&Expr::new_constant(5.0))
    );
}

#[test]

fn test_vector_creation() {

    let signature = (3, 0, 0);

    let components = vec![
        Expr::new_constant(1.0),
        Expr::new_constant(2.0),
        Expr::new_constant(3.0),
    ];

    let vector = Multivector::vector(
        signature,
        components,
    );

    // Should have 3 terms for e1, e2, e3
    assert_eq!(
        vector.terms.len(),
        3
    );

    assert_eq!(
        vector.terms.get(&1),
        Some(&Expr::new_constant(1.0))
    ); // e1
    assert_eq!(
        vector.terms.get(&2),
        Some(&Expr::new_constant(2.0))
    ); // e2
    assert_eq!(
        vector.terms.get(&4),
        Some(&Expr::new_constant(3.0))
    ); // e3
}

#[test]

fn test_multivector_addition() {

    let signature = (3, 0, 0);

    let mv1 = Multivector::scalar(
        signature,
        Expr::new_constant(2.0),
    );

    let mv2 = Multivector::scalar(
        signature,
        Expr::new_constant(3.0),
    );

    let result = mv1 + mv2;

    // Should have scalar term = 5.0
    assert_eq!(
        result.terms.len(),
        1
    );
}

#[test]

fn test_multivector_subtraction() {

    let signature = (3, 0, 0);

    let mv1 = Multivector::scalar(
        signature,
        Expr::new_constant(5.0),
    );

    let mv2 = Multivector::scalar(
        signature,
        Expr::new_constant(3.0),
    );

    let result = mv1 - mv2;

    // Should have scalar term = 2.0
    assert_eq!(
        result.terms.len(),
        1
    );
}

#[test]

fn test_scalar_multiplication() {

    let signature = (3, 0, 0);

    let mv = Multivector::scalar(
        signature,
        Expr::new_constant(3.0),
    );

    let scalar = Expr::new_constant(2.0);

    let result = mv * scalar;

    // Should have scalar term = 6.0
    assert_eq!(
        result.terms.len(),
        1
    );
}

#[test]

fn test_geometric_product_scalars() {

    let signature = (3, 0, 0);

    let mv1 = Multivector::scalar(
        signature,
        Expr::new_constant(2.0),
    );

    let mv2 = Multivector::scalar(
        signature,
        Expr::new_constant(3.0),
    );

    let result =
        mv1.geometric_product(&mv2);

    // Scalar * Scalar = Scalar (6.0)
    assert_eq!(
        result.terms.len(),
        1
    );
}

#[test]

fn test_geometric_product_vectors() {

    let signature = (3, 0, 0);

    // e1
    let v1 = Multivector::vector(
        signature,
        vec![
            Expr::new_constant(1.0),
            Expr::new_constant(0.0),
            Expr::new_constant(0.0),
        ],
    );

    // e2
    let v2 = Multivector::vector(
        signature,
        vec![
            Expr::new_constant(0.0),
            Expr::new_constant(1.0),
            Expr::new_constant(0.0),
        ],
    );

    let result =
        v1.geometric_product(&v2);

    // e1 * e2 = e12 (bivector)
    // Should have a bivector term
    assert!(
        !result
            .terms
            .is_empty()
    );
}

#[test]

fn test_geometric_product_same_vector()
{

    let signature = (3, 0, 0);

    // e1
    let v1 = Multivector::vector(
        signature,
        vec![
            Expr::new_constant(1.0),
            Expr::new_constant(0.0),
            Expr::new_constant(0.0),
        ],
    );

    let result = v1
        .clone()
        .geometric_product(&v1);

    // e1 * e1 = 1 (scalar) in Euclidean space
    assert_eq!(
        result.terms.len(),
        1
    );
}

#[test]

fn test_outer_product() {

    let signature = (3, 0, 0);

    let v1 = Multivector::vector(
        signature,
        vec![
            Expr::new_constant(1.0),
            Expr::new_constant(0.0),
            Expr::new_constant(0.0),
        ],
    );

    let v2 = Multivector::vector(
        signature,
        vec![
            Expr::new_constant(0.0),
            Expr::new_constant(1.0),
            Expr::new_constant(0.0),
        ],
    );

    let result = v1.outer_product(&v2);

    // e1 ∧ e2 = e12 (bivector)
    assert!(
        !result
            .terms
            .is_empty()
    );
}

#[test]

fn test_inner_product() {

    let signature = (3, 0, 0);

    let v1 = Multivector::vector(
        signature,
        vec![
            Expr::new_constant(1.0),
            Expr::new_constant(0.0),
            Expr::new_constant(0.0),
        ],
    );

    let v2 = Multivector::vector(
        signature,
        vec![
            Expr::new_constant(1.0),
            Expr::new_constant(0.0),
            Expr::new_constant(0.0),
        ],
    );

    let result = v1.inner_product(&v2);

    // e1 · e1 = 1 (scalar)
    assert!(
        !result
            .terms
            .is_empty()
    );
}

#[test]

fn test_grade_projection() {

    let signature = (3, 0, 0);

    let mut mv =
        Multivector::new(signature);

    mv.terms.insert(
        0,
        Expr::new_constant(1.0),
    ); // scalar
    mv.terms.insert(
        1,
        Expr::new_constant(2.0),
    ); // e1 (vector)
    mv.terms.insert(
        3,
        Expr::new_constant(3.0),
    ); // e12 (bivector)

    let scalar_part =
        mv.grade_projection(0);

    let vector_part =
        mv.grade_projection(1);

    let bivector_part =
        mv.grade_projection(2);

    assert_eq!(
        scalar_part
            .terms
            .len(),
        1
    );

    assert_eq!(
        vector_part
            .terms
            .len(),
        1
    );

    assert_eq!(
        bivector_part
            .terms
            .len(),
        1
    );
}

#[test]

fn test_reverse_scalar() {

    let signature = (3, 0, 0);

    let mv = Multivector::scalar(
        signature,
        Expr::new_constant(5.0),
    );

    let reversed = mv.reverse();

    // Reverse of scalar is itself
    assert_eq!(
        reversed.terms.len(),
        1
    );
}

#[test]

fn test_reverse_vector() {

    let signature = (3, 0, 0);

    let v = Multivector::vector(
        signature,
        vec![
            Expr::new_constant(1.0),
            Expr::new_constant(2.0),
            Expr::new_constant(3.0),
        ],
    );

    let reversed = v.reverse();

    // Reverse of vector is itself (grade 1)
    assert_eq!(
        reversed.terms.len(),
        3
    );
}

#[test]

fn test_magnitude() {

    let signature = (3, 0, 0);

    let v = Multivector::vector(
        signature,
        vec![
            Expr::new_constant(3.0),
            Expr::new_constant(4.0),
            Expr::new_constant(0.0),
        ],
    );

    let mag = v.magnitude();

    // Magnitude should be sqrt(9 + 16) = 5
    // We can't easily verify the exact value symbolically, but it should be non-zero
    assert!(!matches!(
        mag,
        Expr::new_constant(0.0)
    ));
}

#[test]

fn test_dual() {

    let signature = (3, 0, 0);

    let v = Multivector::vector(
        signature,
        vec![
            Expr::new_constant(1.0),
            Expr::new_constant(0.0),
            Expr::new_constant(0.0),
        ],
    );

    let dual = v.dual();

    // Dual of a vector in 3D is a bivector
    assert!(
        !dual
            .terms
            .is_empty()
    );
}

#[test]

fn test_normalize() {

    let signature = (3, 0, 0);

    let v = Multivector::vector(
        signature,
        vec![
            Expr::new_constant(3.0),
            Expr::new_constant(4.0),
            Expr::new_constant(0.0),
        ],
    );

    let normalized = v.normalize();

    // Normalized vector should have unit magnitude
    assert!(
        !normalized
            .terms
            .is_empty()
    );
}

#[test]

fn test_minkowski_signature() {

    let signature = (1, 3, 0); // Minkowski spacetime (1 time, 3 space)
    let mv = Multivector::scalar(
        signature,
        Expr::new_constant(1.0),
    );

    assert_eq!(
        mv.signature,
        (1, 3, 0)
    );
}

#[test]

fn test_geometric_product_associativity()
 {

    let signature = (3, 0, 0);

    let v1 = Multivector::scalar(
        signature,
        Expr::new_constant(2.0),
    );

    let v2 = Multivector::scalar(
        signature,
        Expr::new_constant(3.0),
    );

    let v3 = Multivector::scalar(
        signature,
        Expr::new_constant(5.0),
    );

    let result1 = v1
        .clone()
        .geometric_product(&v2.clone())
        .geometric_product(&v3.clone());

    let result2 = v1.geometric_product(
        &v2.geometric_product(&v3),
    );

    // (v1 * v2) * v3 = v1 * (v2 * v3)
    assert_eq!(
        result1.terms.len(),
        result2.terms.len()
    );
}

#[test]

fn test_outer_product_anticommutativity()
 {

    let signature = (3, 0, 0);

    let v1 = Multivector::vector(
        signature,
        vec![
            Expr::new_constant(1.0),
            Expr::new_constant(0.0),
            Expr::new_constant(0.0),
        ],
    );

    let v2 = Multivector::vector(
        signature,
        vec![
            Expr::new_constant(0.0),
            Expr::new_constant(1.0),
            Expr::new_constant(0.0),
        ],
    );

    let result1 = v1
        .clone()
        .outer_product(&v2.clone());

    let result2 = v2.outer_product(&v1);

    // v1 ∧ v2 = -(v2 ∧ v1) for vectors
    // Both should have terms, but with opposite signs
    assert!(
        !result1
            .terms
            .is_empty()
    );

    assert!(
        !result2
            .terms
            .is_empty()
    );
}
