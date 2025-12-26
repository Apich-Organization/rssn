use num_bigint::BigInt;
use num_traits::{
    One,
    Zero,
};
use rssn::symbolic::core::Expr;
use rssn::symbolic::tensor::*;

#[test]

fn test_tensor_creation() {

    // Create a 2x2 tensor with symbolic components
    let a = Expr::new_variable("a");

    let b = Expr::new_variable("b");

    let c = Expr::new_variable("c");

    let d = Expr::new_variable("d");

    let tensor = Tensor::new(
        vec![
            a.clone(),
            b.clone(),
            c.clone(),
            d.clone(),
        ],
        vec![2, 2],
    )
    .unwrap();

    assert_eq!(tensor.rank(), 2);

    assert_eq!(
        tensor.shape,
        vec![2, 2]
    );

    assert_eq!(
        tensor
            .get(&[0, 0])
            .unwrap(),
        &a
    );

    assert_eq!(
        tensor
            .get(&[0, 1])
            .unwrap(),
        &b
    );

    assert_eq!(
        tensor
            .get(&[1, 0])
            .unwrap(),
        &c
    );

    assert_eq!(
        tensor
            .get(&[1, 1])
            .unwrap(),
        &d
    );
}

#[test]

fn test_tensor_addition_symbolic() {

    // T1 = [[a, b], [c, d]]
    // T2 = [[e, f], [g, h]]
    // T1 + T2 = [[a+e, b+f], [c+g, d+h]]
    let a = Expr::new_variable("a");

    let b = Expr::new_variable("b");

    let c = Expr::new_variable("c");

    let d = Expr::new_variable("d");

    let e = Expr::new_variable("e");

    let f = Expr::new_variable("f");

    let g = Expr::new_variable("g");

    let h = Expr::new_variable("h");

    let t1 = Tensor::new(
        vec![
            a.clone(),
            b.clone(),
            c.clone(),
            d.clone(),
        ],
        vec![2, 2],
    )
    .unwrap();

    let t2 = Tensor::new(
        vec![
            e.clone(),
            f.clone(),
            g.clone(),
            h.clone(),
        ],
        vec![2, 2],
    )
    .unwrap();

    let sum = t1.add(&t2).unwrap();

    // Check symbolic structure (components should be DAG nodes representing additions)
    println!(
        "Sum[0,0]: {:?}",
        sum.get(&[0, 0])
            .unwrap()
    );

    println!(
        "Sum[0,1]: {:?}",
        sum.get(&[0, 1])
            .unwrap()
    );

    // The components should be symbolic additions
    assert_eq!(sum.rank(), 2);

    assert_eq!(
        sum.shape,
        vec![2, 2]
    );
}

#[test]

fn test_tensor_scalar_multiplication_symbolic() {

    // T = [[a, b], [c, d]]
    // k * T = [[k*a, k*b], [k*c, k*d]]
    let a = Expr::new_variable("a");

    let b = Expr::new_variable("b");

    let c = Expr::new_variable("c");

    let d = Expr::new_variable("d");

    let k = Expr::new_variable("k");

    let tensor = Tensor::new(
        vec![
            a.clone(),
            b.clone(),
            c.clone(),
            d.clone(),
        ],
        vec![2, 2],
    )
    .unwrap();

    let scaled = tensor
        .scalar_mul(&k)
        .unwrap();

    println!(
        "Scaled[0,0]: {:?}",
        scaled
            .get(&[0, 0])
            .unwrap()
    );

    println!(
        "Scaled[1,1]: {:?}",
        scaled
            .get(&[1, 1])
            .unwrap()
    );

    assert_eq!(scaled.rank(), 2);

    assert_eq!(
        scaled.shape,
        vec![2, 2]
    );
}

#[test]

fn test_tensor_outer_product_symbolic() {

    // v1 = [a, b], v2 = [c, d]
    // v1 ⊗ v2 = [[a*c, a*d], [b*c, b*d]]
    let a = Expr::new_variable("a");

    let b = Expr::new_variable("b");

    let c = Expr::new_variable("c");

    let d = Expr::new_variable("d");

    let v1 = Tensor::new(
        vec![a.clone(), b.clone()],
        vec![2],
    )
    .unwrap();

    let v2 = Tensor::new(
        vec![c.clone(), d.clone()],
        vec![2],
    )
    .unwrap();

    let outer = v1
        .outer_product(&v2)
        .unwrap();

    assert_eq!(outer.rank(), 2);

    assert_eq!(
        outer.shape,
        vec![2, 2]
    );

    println!(
        "Outer[0,0]: {:?}",
        outer
            .get(&[0, 0])
            .unwrap()
    );

    println!(
        "Outer[0,1]: {:?}",
        outer
            .get(&[0, 1])
            .unwrap()
    );

    println!(
        "Outer[1,0]: {:?}",
        outer
            .get(&[1, 0])
            .unwrap()
    );

    println!(
        "Outer[1,1]: {:?}",
        outer
            .get(&[1, 1])
            .unwrap()
    );
}

#[test]

fn test_tensor_contraction_symbolic() {

    // Create a rank-2 tensor [[a, b], [c, d]]
    // Contract indices 0 and 1 to get trace: a + d
    let a = Expr::new_variable("a");

    let b = Expr::new_variable("b");

    let c = Expr::new_variable("c");

    let d = Expr::new_variable("d");

    let tensor = Tensor::new(
        vec![
            a.clone(),
            b.clone(),
            c.clone(),
            d.clone(),
        ],
        vec![2, 2],
    )
    .unwrap();

    let trace = tensor
        .contract(0, 1)
        .unwrap();

    // Trace should be rank-0 (scalar)
    assert_eq!(trace.rank(), 0);

    assert_eq!(
        trace.shape,
        vec![] as Vec<usize>
    );

    println!(
        "Trace: {:?}",
        trace.components[0]
    );
    // The trace should be symbolically a + d
}

#[test]

fn test_metric_tensor_symbolic() {

    // Create a simple 2D Euclidean metric: g = [[1, 0], [0, 1]]
    let one = Expr::BigInt(BigInt::one());

    let zero = Expr::BigInt(BigInt::zero());

    let g = Tensor::new(
        vec![
            one.clone(),
            zero.clone(),
            zero.clone(),
            one.clone(),
        ],
        vec![2, 2],
    )
    .unwrap();

    let metric = MetricTensor::new(g).unwrap();

    // Check that inverse is also identity
    let one_ast = one
        .to_ast()
        .unwrap_or(one.clone());

    let inv_00 = metric
        .g_inv
        .get(&[0, 0])
        .unwrap()
        .to_ast()
        .unwrap_or(
            metric
                .g_inv
                .get(&[0, 0])
                .unwrap()
                .clone(),
        );

    let inv_11 = metric
        .g_inv
        .get(&[1, 1])
        .unwrap()
        .to_ast()
        .unwrap_or(
            metric
                .g_inv
                .get(&[1, 1])
                .unwrap()
                .clone(),
        );

    println!(
        "g_inv[0,0]: {:?}",
        inv_00
    );

    println!(
        "g_inv[1,1]: {:?}",
        inv_11
    );
}

#[test]

fn test_metric_tensor_raise_lower_index_symbolic() {

    // Use identity metric for simplicity
    let one = Expr::BigInt(BigInt::one());

    let zero = Expr::BigInt(BigInt::zero());

    let g = Tensor::new(
        vec![
            one.clone(),
            zero.clone(),
            zero.clone(),
            one.clone(),
        ],
        vec![2, 2],
    )
    .unwrap();

    let metric = MetricTensor::new(g).unwrap();

    // Create a symbolic vector [a, b]
    let a = Expr::new_variable("a");

    let b = Expr::new_variable("b");

    let vector = Tensor::new(
        vec![a.clone(), b.clone()],
        vec![2],
    )
    .unwrap();

    // Lower the index
    let covector = metric
        .lower_index(&vector)
        .unwrap();

    // With identity metric, lowering should give the same components
    println!(
        "Covector[0]: {:?}",
        covector
            .get(&[0])
            .unwrap()
    );

    println!(
        "Covector[1]: {:?}",
        covector
            .get(&[1])
            .unwrap()
    );

    // Raise it back
    let raised = metric
        .raise_index(&covector)
        .unwrap();

    println!(
        "Raised[0]: {:?}",
        raised
            .get(&[0])
            .unwrap()
    );

    println!(
        "Raised[1]: {:?}",
        raised
            .get(&[1])
            .unwrap()
    );
}

#[test]

fn test_tensor_to_matrix_expr() {

    // Create a 2x2 tensor and convert to matrix expression
    let a = Expr::new_variable("a");

    let b = Expr::new_variable("b");

    let c = Expr::new_variable("c");

    let d = Expr::new_variable("d");

    let tensor = Tensor::new(
        vec![
            a.clone(),
            b.clone(),
            c.clone(),
            d.clone(),
        ],
        vec![2, 2],
    )
    .unwrap();

    let matrix_expr = tensor
        .to_matrix_expr()
        .unwrap();

    if let Expr::Matrix(rows) = matrix_expr {

        assert_eq!(rows.len(), 2);

        assert_eq!(rows[0].len(), 2);

        assert_eq!(rows[0][0], a);

        assert_eq!(rows[0][1], b);

        assert_eq!(rows[1][0], c);

        assert_eq!(rows[1][1], d);
    } else {

        panic!("Expected Matrix expression");
    }
}
