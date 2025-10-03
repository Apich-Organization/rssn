use rssn::symbolic::core::Expr;
use rssn::symbolic::unit_unification::{unify_expression, UnitQuantity, SupportedQuantity};
use uom::si::f64::Length;
use uom::si::length::meter;

#[test]
fn test_add_meters_and_centimeters() {
    let ten_meters = Expr::QuantityWithValue(Box::new(Expr::Constant(10.0)), "m".to_string());
    let five_cm = Expr::QuantityWithValue(Box::new(Expr::Constant(5.0)), "cm".to_string());

    let addition_expr = Expr::Add(Box::new(ten_meters), Box::new(five_cm));

    let result = unify_expression(&addition_expr).unwrap();

    if let Expr::Quantity(unit_quantity) = result {
        if let SupportedQuantity::Length(l) = unit_quantity.0 {
            // uom defaults to the base unit, which is meter for length.
            assert!((l.value - 10.05).abs() < 1e-9);
        } else {
            panic!("Expected the result to be a Length quantity");
        }
    } else {
        panic!("Expected the result to be a Quantity, but got {:?}", result);
    }
}
