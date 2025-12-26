use rssn::symbolic::core::Expr;
use rssn::symbolic::unit_unification::{
    unify_expression,
    SupportedQuantity,
    UnitQuantity,
};
use std::sync::Arc;

#[test]

fn test_parse_and_unify_length() {

    // 5 meters
    let expr = Expr::QuantityWithValue(
        Arc::new(Expr::new_constant(
            5.0,
        )),
        "meter".to_string(),
    );

    let unified = unify_expression(&expr).unwrap();

    if let Expr::Quantity(q) = unified {

        match &q.0 {
            SupportedQuantity::Length(l) => assert_eq!(l.value, 5.0),
            _ => panic!("Expected Length"),
        }
    } else {

        panic!("Expected Quantity");
    }
}

#[test]

fn test_add_same_units() {

    // 5m + 3m = 8m
    let q1 = Expr::QuantityWithValue(
        Arc::new(Expr::new_constant(
            5.0,
        )),
        "m".to_string(),
    );

    let q2 = Expr::QuantityWithValue(
        Arc::new(Expr::new_constant(
            3.0,
        )),
        "m".to_string(),
    );

    let expr = Expr::new_add(q1, q2);

    let unified = unify_expression(&expr).unwrap();

    if let Expr::Quantity(q) = unified {

        match &q.0 {
            SupportedQuantity::Length(l) => assert_eq!(l.value, 8.0),
            _ => panic!("Expected Length"),
        }
    } else {

        panic!("Expected Quantity");
    }
}

#[test]

fn test_multiply_units() {

    // 2m * 3m = 6m^2 (Area)
    let q1 = Expr::QuantityWithValue(
        Arc::new(Expr::new_constant(
            2.0,
        )),
        "m".to_string(),
    );

    let q2 = Expr::QuantityWithValue(
        Arc::new(Expr::new_constant(
            3.0,
        )),
        "m".to_string(),
    );

    let expr = Expr::new_mul(q1, q2);

    let unified = unify_expression(&expr).unwrap();

    if let Expr::Quantity(q) = unified {

        match &q.0 {
            SupportedQuantity::Area(a) => assert_eq!(a.value, 6.0),
            _ => panic!("Expected Area"),
        }
    } else {

        panic!("Expected Quantity");
    }
}

#[test]

fn test_divide_units() {

    // 10m / 2s = 5m/s (Velocity)
    let q1 = Expr::QuantityWithValue(
        Arc::new(Expr::new_constant(
            10.0,
        )),
        "m".to_string(),
    );

    let q2 = Expr::QuantityWithValue(
        Arc::new(Expr::new_constant(
            2.0,
        )),
        "s".to_string(),
    );

    let expr = Expr::new_div(q1, q2);

    let unified = unify_expression(&expr).unwrap();

    if let Expr::Quantity(q) = unified {

        match &q.0 {
            SupportedQuantity::Velocity(v) => assert_eq!(v.value, 5.0),
            _ => panic!("Expected Velocity"),
        }
    } else {

        panic!("Expected Quantity");
    }
}

#[test]

fn test_scalar_multiplication() {

    // 3 * 4kg = 12kg
    let scalar = Expr::new_constant(3.0);

    let q = Expr::QuantityWithValue(
        Arc::new(Expr::new_constant(
            4.0,
        )),
        "kg".to_string(),
    );

    let expr = Expr::new_mul(scalar, q);

    let unified = unify_expression(&expr).unwrap();

    if let Expr::Quantity(q) = unified {

        match &q.0 {
            SupportedQuantity::Mass(m) => assert_eq!(m.value, 12.0),
            _ => panic!("Expected Mass"),
        }
    } else {

        panic!("Expected Quantity");
    }
}
