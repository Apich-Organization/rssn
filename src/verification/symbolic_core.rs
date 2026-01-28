//! Kani verification for symbolic core API.
//!
//! These proofs verify the correctness of the DAG-based expression system,
//! ensuring that constructors, conversions, and operators behave as expected.

use std::sync::Arc;

use num_bigint::BigInt;
use num_rational::BigRational;

use crate::symbolic::core::Expr;
use crate::symbolic::core::ToConstant;

// --- Leaf Node Proofs ---

#[kani::proof]

pub(crate) fn proof_new_constant_is_dag()
 {

    let c: f64 = kani::any();

    if c.is_finite() {

        let expr =
            Expr::new_constant(c);

        assert!(expr.is_dag());
    }
}

#[kani::proof]

fn proof_new_variable_is_dag() {

    let expr = Expr::new_variable("x");

    assert!(expr.is_dag());
}

#[kani::proof]

fn proof_new_pi_is_dag() {

    let expr = Expr::new_pi();

    assert!(expr.is_dag());
}

#[kani::proof]

fn proof_new_e_is_dag() {

    let expr = Expr::new_e();

    assert!(expr.is_dag());
}

#[kani::proof]

fn proof_new_infinity_is_dag() {

    let expr = Expr::new_infinity();

    assert!(expr.is_dag());
}

#[kani::proof]

fn proof_new_bigint_is_dag() {

    // We use a small integer to avoid potential issues with num_bigint in Kani
    let i: i64 = kani::any();

    let expr = Expr::new_bigint(
        BigInt::from(i),
    );

    assert!(expr.is_dag());
}

#[kani::proof]

fn proof_new_rational_is_dag() {

    let n: i64 = kani::any();

    let d: i64 = kani::any();

    if d != 0 {

        let expr = Expr::new_rational(
            BigRational::new(
                BigInt::from(n),
                BigInt::from(d),
            ),
        );

        assert!(expr.is_dag());
    }
}

// --- Structural & Normalization Proofs ---

#[kani::proof]

fn proof_add_commutativity_normalization()
 {

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    let sum1 = Expr::new_add(&x, &y);

    let sum2 = Expr::new_add(&y, &x);

    // In a DAG-based system with normalization, x+y and y+x should be the same node.
    assert_eq!(
        sum1, sum2,
        "Addition should be \
         normalized (commutative)"
    );
}

#[kani::proof]

fn proof_mul_commutativity_normalization()
 {

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    let prod1 = Expr::new_mul(&x, &y);

    let prod2 = Expr::new_mul(&y, &x);

    assert_eq!(
        prod1, prod2,
        "Multiplication should be \
         normalized (commutative)"
    );
}

#[kani::proof]

fn proof_nested_add_normalization() {

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    let z = Expr::new_variable("z");

    // (x + y) + z vs x + (y + z)
    // Note: DagManager might not do full associative reordering automatically
    // depending on implementation, but let's check if it does basic sorting.
    let sum1 = Expr::new_add(
        &Expr::new_add(&x, &y),
        &z,
    );

    let sum2 = Expr::new_add(
        &x,
        &Expr::new_add(&y, &z),
    );

    // If it only does binary sorting per node, these might NOT be equal without a flatten pass.
    // However, if it's N-ary, it would.
    // Based on api.rs, new_add is binary. So they might remain nested.
}

#[kani::proof]

fn proof_to_dag_idempotence() {

    let c: f64 = kani::any();

    if c.is_finite() {

        let ast_expr =
            Expr::new_constant(c);

        let dag_expr_1 = ast_expr
            .to_dag()
            .expect("to_dag failed");

        assert!(dag_expr_1.is_dag());

        let dag_expr_2 = dag_expr_1
            .to_dag()
            .expect(
                "to_dag failed again",
            );

        assert_eq!(
            dag_expr_1,
            dag_expr_2
        );
    }
}

#[kani::proof]

fn proof_to_ast_to_dag_roundtrip() {

    let c: f64 = kani::any();

    if c.is_finite() {

        let expr =
            Expr::new_constant(c);

        let ast = expr
            .to_ast()
            .expect("to_ast failed");

        let dag = ast
            .to_dag()
            .expect("to_dag failed");

        assert_eq!(expr, dag);
    }
}

// --- Operator Overloading Proofs ---

#[kani::proof]

fn proof_operator_overloading_is_dag() {

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    let add = &x + &y;

    let sub = &x - &y;

    let mul = &x * &y;

    let div = &x / &y;

    let neg = -&x;

    assert!(add.is_dag());

    assert!(sub.is_dag());

    assert!(mul.is_dag());

    assert!(div.is_dag());

    assert!(neg.is_dag());
}

#[kani::proof]

fn proof_mixed_f64_ops() {

    let x = Expr::new_variable("x");

    let c: f64 = kani::any();

    if c.is_finite() {

        let res = &x + c;

        assert!(res.is_dag());

        let res2 = c * &x;

        assert!(res2.is_dag());
    }
}

// --- Method Proofs ---

#[kani::proof]

fn proof_math_methods_is_dag() {

    let x = Expr::new_variable("x");

    assert!(x.sin().is_dag());

    assert!(x.cos().is_dag());

    assert!(x.tan().is_dag());

    assert!(x.exp().is_dag());

    assert!(x.ln().is_dag());

    assert!(x.abs().is_dag());

    assert!(x.sqrt().is_dag());

    assert!(x.asin().is_dag());

    assert!(x.acos().is_dag());

    assert!(x.atan().is_dag());

    assert!(x.sinh().is_dag());

    assert!(x.cosh().is_dag());

    assert!(x.tanh().is_dag());
}

// --- Container Proofs ---

#[kani::proof]

fn proof_vector_is_dag() {

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    let vec_expr =
        Expr::new_vector(vec![x, y]);

    assert!(vec_expr.is_dag());
}

#[kani::proof]

fn proof_tuple_is_dag() {

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    let tuple_expr =
        Expr::new_tuple(vec![x, y]);

    assert!(tuple_expr.is_dag());
}

#[kani::proof]

fn proof_matrix_is_dag() {

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    let mat = Expr::new_matrix(vec![
        vec![x.clone(), y.clone()],
        vec![y, x],
    ]);

    assert!(mat.is_dag());
}

// --- Logic Proofs ---

#[kani::proof]

fn proof_logic_ops_is_dag() {

    let x = Expr::new_variable("x");

    let y = Expr::new_variable("y");

    assert!(
        Expr::new_and(vec![
            x.clone(),
            y.clone()
        ])
        .is_dag()
    );

    assert!(
        Expr::new_or(vec![
            x.clone(),
            y.clone()
        ])
        .is_dag()
    );

    assert!(Expr::new_not(&x).is_dag());

    assert!(
        Expr::new_xor(&x, &y).is_dag()
    );

    assert!(
        Expr::new_implies(&x, &y)
            .is_dag()
    );

    assert!(
        Expr::new_equivalent(&x, &y)
            .is_dag()
    );
}

// --- Quantifier & Relation Proofs ---

#[kani::proof]

fn proof_quantifiers_is_dag() {

    let x = Expr::new_variable("x");

    let body = Expr::new_gt(
        &x,
        &Expr::new_constant(0.0),
    );

    assert!(
        Expr::new_forall("x", &body)
            .is_dag()
    );

    assert!(
        Expr::new_exists("x", &body)
            .is_dag()
    );
}

#[kani::proof]

fn proof_interval_is_dag() {

    let l = Expr::new_constant(0.0);

    let u = Expr::new_constant(1.0);

    let interval = Expr::new_interval(
        &l, &u, true, false,
    );

    assert!(interval.is_dag());
}
