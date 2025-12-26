use num_complex::Complex;
use rssn::symbolic::core::Expr;
use serde::Serialize;
use std::collections::HashMap;

#[derive(Serialize)]

struct EvalInput {
    expr: Expr,
    vars: HashMap<String, Complex<f64>>,
}

#[test]

fn see_json() {

    let mut vars = HashMap::new();

    vars.insert("z".to_string(), Complex::new(0.0, 1.0));

    let input = EvalInput {
        expr: Expr::Variable("z".to_string()),
        vars,
    };

    println!("{}", serde_json::to_string_pretty(&input).unwrap());
}
