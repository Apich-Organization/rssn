//! Symbolic indefinite sum (anti-difference) examples.
//!
//! Demonstrates:
//!   1. Constructing `IndefiniteSum` / `IndefiniteProduct` AST nodes
//!   2. Simplifying them to closed-form expressions via `rssn::simplify`
//!   3. Evaluating the simplified expression and verifying the recurrence
//!      F(x+1) - F(x) = f(x)
//!   4. Linearity: Δ⁻¹(a·f + b·g) = a·Δ⁻¹f + b·Δ⁻¹g
//!   5. Non-unit step sizes
//!   6. Symbolic indefinite product P(x)/P(x-1) = f(x-1)

use rssn::numerical::elementary::eval_expr_single;
use rssn::numerical::indefinite_sum::eval_antidiff;
use rssn::numerical::indefinite_sum::try_closed_form_sum;
use rssn::simplify;
use rssn::symbolic::core::Expr;
use std::sync::Arc;

fn main() {
    println!("=== Symbolic Indefinite Sum Examples ===\n");

    // ------------------------------------------------------------------
    // 1. IndefiniteSum(3, x, 1)  →  simplified to  3·x
    // ------------------------------------------------------------------
    println!("--- 1. Constant: IndefiniteSum(3, x, 1) ---");
    let body = Expr::new_constant(3.0);
    let step = Expr::new_constant(1.0);
    let expr = Expr::new_indefinite_sum(body, "x".to_string(), step);
    println!("  Before simplify: {expr}");
    let simplified = simplify(&expr);
    println!("  After  simplify: {simplified}");
    // Evaluate: Δ⁻¹ 3 = 3x  →  F(5) = 15
    let val = eval_expr_single(&simplified, "x", 5.0).unwrap();
    println!("  F(5) = {val}  (expected 3·5 = 15)\n");

    // ------------------------------------------------------------------
    // 2. IndefiniteSum(eˣ, x, 1)  →  simplified to  eˣ / (e - 1)
    // ------------------------------------------------------------------
    println!("--- 2. Exponential: IndefiniteSum(e^x, x, 1) ---");
    let body = Expr::Exp(Arc::new(Expr::new_variable("x")));
    let step = Expr::new_constant(1.0);
    let sum_expr = Expr::new_indefinite_sum(body, "x".to_string(), step);
    println!("  Before simplify: {sum_expr}");
    let simplified = simplify(&sum_expr);
    println!("  After  simplify: {simplified}");
    // Verify recurrence F(x+1) - F(x) = e^x
    for x in [1.0_f64, 2.0, 3.0] {
        let fxp1 = eval_expr_single(&simplified, "x", x + 1.0).unwrap();
        let fx = eval_expr_single(&simplified, "x", x).unwrap();
        println!(
            "  F({x:.0}+1) - F({x:.0}) = {:.8}  (expected e^{x:.0} = {:.8})",
            fxp1 - fx,
            x.exp()
        );
    }
    println!();

    // ------------------------------------------------------------------
    // 3. IndefiniteSum(sin(0.5·x), x, 1)
    // ------------------------------------------------------------------
    println!("--- 3. Trigonometric: IndefiniteSum(sin(0.5·x), x, 1) ---");
    let a = 0.5_f64;
    let ax = Expr::new_mul(Expr::new_constant(a), Expr::new_variable("x"));
    let body = Expr::Sin(Arc::new(ax));
    let step = Expr::new_constant(1.0);
    let sum_expr = Expr::new_indefinite_sum(body, "x".to_string(), step);
    let simplified = simplify(&sum_expr);
    println!("  Simplified: {simplified}");
    let x = 3.0_f64;
    let fxp1 = eval_expr_single(&simplified, "x", x + 1.0).unwrap();
    let fx = eval_expr_single(&simplified, "x", x).unwrap();
    println!(
        "  F(4) - F(3) = {:.8}  (expected sin(0.5·3) = {:.8})\n",
        fxp1 - fx,
        (a * x).sin()
    );

    // ------------------------------------------------------------------
    // 4. IndefiniteSum(x^3, x, 1)  via Hurwitz zeta  ζ(−3,1) − ζ(−3,x)
    // ------------------------------------------------------------------
    println!("--- 4. Power: IndefiniteSum(x^3, x, 1) ---");
    let body = Expr::Power(
        Arc::new(Expr::new_variable("x")),
        Arc::new(Expr::new_constant(3.0)),
    );
    let step = Expr::new_constant(1.0);
    let sum_expr = Expr::new_indefinite_sum(body.clone(), "x".to_string(), step);
    let simplified = simplify(&sum_expr);
    println!("  Simplified form: {simplified}");
    // eval_expr_single may not handle BinaryList("hurwitz_zeta_antidiff"); use eval_antidiff
    let anti = try_closed_form_sum(&body, "x").unwrap();
    for x in [2.0_f64, 4.0, 6.0] {
        let fxp1 = eval_antidiff(&anti, "x", x + 1.0).unwrap();
        let fx = eval_antidiff(&anti, "x", x).unwrap();
        println!(
            "  F({x:.0}+1) - F({x:.0}) = {:.2}  (expected {x:.0}^3 = {:.0})",
            fxp1 - fx,
            x.powi(3)
        );
    }
    println!();

    // ------------------------------------------------------------------
    // 5. Linearity: Δ⁻¹(2·sin(x) + 3·eˣ) = 2·Δ⁻¹sin(x) + 3·Δ⁻¹eˣ
    // ------------------------------------------------------------------
    println!("--- 5. Linearity: IndefiniteSum(2·sin(x) + 3·e^x, x, 1) ---");
    let sin_x = Expr::Sin(Arc::new(Expr::new_variable("x")));
    let exp_x = Expr::Exp(Arc::new(Expr::new_variable("x")));
    let body = Expr::new_add(
        Expr::new_mul(Expr::new_constant(2.0), sin_x),
        Expr::new_mul(Expr::new_constant(3.0), exp_x),
    );
    let sum_expr = Expr::new_indefinite_sum(body, "x".to_string(), Expr::new_constant(1.0));
    let simplified = simplify(&sum_expr);
    let x = 2.5_f64;
    if let (Ok(fxp1), Ok(fx)) = (
        eval_expr_single(&simplified, "x", x + 1.0),
        eval_expr_single(&simplified, "x", x),
    ) {
        let expected = 2.0 * x.sin() + 3.0 * x.exp();
        println!(
            "  F({x}+1) - F({x}) = {:.8}  (expected {:.8})\n",
            fxp1 - fx,
            expected
        );
    } else {
        println!("  (simplified form requires numeric fallback)\n");
    }

    // ------------------------------------------------------------------
    // 6. Classic result: 1 + 2 + ... + n = n(n+1)/2
    //    IndefiniteSum(x, x, 1) → Δ⁻¹ x = x(x-1)/2
    //    Σ_{k=1}^n k = F(n+1) - F(1)
    // ------------------------------------------------------------------
    println!("--- 6. Definite sum: 1 + 2 + ... + n = n(n+1)/2 ---");
    // Use plain Variable("x") so IndefiniteSum(x, x, 1) simplifies via the Variable rule.
    let body_x = Expr::Variable("x".to_string());
    let sum_of_k =
        Expr::new_indefinite_sum(body_x.clone(), "x".to_string(), Expr::new_constant(1.0));
    let anti_x = simplify(&sum_of_k);
    println!("  IndefiniteSum(x, x, 1) simplified: {anti_x}");
    let anti = try_closed_form_sum(&body_x, "x").unwrap();
    let f1 = eval_antidiff(&anti, "x", 1.0).unwrap(); // F(1) = 0
    for n in [5u32, 10, 50, 100] {
        let fn1 = eval_antidiff(&anti, "x", n as f64 + 1.0).unwrap();
        let sum = fn1 - f1;
        let expected = (n * (n + 1) / 2) as f64;
        println!("  Σ k, k=1..{n:3} = {sum:8.1}  (n(n+1)/2 = {expected:8.0})");
    }
    println!();

    // ------------------------------------------------------------------
    // 8. Non-unit step: IndefiniteSum(eˣ, x, 2)
    //    Recurrence: F(x+2) - F(x) = eˣ
    // ------------------------------------------------------------------
    println!("--- 8. Non-unit step: IndefiniteSum(e^x, x, step=2) ---");
    let body = Expr::Exp(Arc::new(Expr::new_variable("x")));
    let sum_expr = Expr::new_indefinite_sum(body, "x".to_string(), Expr::new_constant(2.0));
    let simplified = simplify(&sum_expr);
    println!("  Simplified: {simplified}");
    // Verify recurrence using the simplified symbolic expression directly.
    // simplify correctly produces e^x/(e^2-1) for step=2.
    for x in [1.0_f64, 2.0, 3.0] {
        let fxp2 = eval_expr_single(&simplified, "x", x + 2.0).unwrap();
        let fx_val = eval_expr_single(&simplified, "x", x).unwrap();
        println!(
            "  F({x:.0}+2) - F({x:.0}) = {:.8}  (expected e^{x:.0} = {:.8})",
            fxp2 - fx_val,
            x.exp()
        );
    }
    println!();

    // ------------------------------------------------------------------
    // 9. IndefiniteProduct: symbolic construction
    //    ∏ f(k) = exp(Δ⁻¹ ln f(x))  (log-sum reduction)
    // ------------------------------------------------------------------
    println!("--- 9. Symbolic IndefiniteProduct: IndefiniteProduct(2^x, x, 1) ---");
    // f(x) = 2^x  →  ln(2^x) = x·ln2  →  Δ⁻¹(x·ln2) = ln2·x(x-1)/2
    // →  ∏_{k=1}^{x-1} 2^k = exp(ln2·x(x-1)/2) = 2^(x(x-1)/2)
    let base = Expr::new_constant(2.0);
    let f_power = Expr::Power(Arc::new(base), Arc::new(Expr::new_variable("x")));
    let prod_expr = Expr::new_indefinite_product(f_power, "x".to_string(), Expr::new_constant(1.0));
    println!("  Before simplify: {prod_expr}");
    let simplified = simplify(&prod_expr);
    println!("  After  simplify: {simplified}");
    // Verify numerically: ∏_{k=1}^{x-1} 2^k  vs  2^(x(x-1)/2)
    use rssn::numerical::indefinite_sum::eval_indefinite_product_numerical;
    let f_2k = |k: f64| -> Result<f64, String> { Ok(2.0_f64.powf(k)) };
    for x in [3.0_f64, 4.0, 5.0] {
        let product = eval_indefinite_product_numerical(x, &f_2k, 1.0, 1.0).unwrap();
        let formula = 2.0_f64.powf(x * (x - 1.0) / 2.0);
        println!("  prod_{{k=1..{x:.0}-1}} 2^k = {product:.2}  (2^(x(x-1)/2) = {formula:.2})");
    }
    println!();

    // ------------------------------------------------------------------
    // 10. Display and inspection of IndefiniteSum
    // ------------------------------------------------------------------
    println!("--- 10. Display of symbolic expression ---");
    let cos_x = Expr::Cos(Arc::new(Expr::new_variable("x")));
    let sum_cos = Expr::new_indefinite_sum(cos_x, "x".to_string(), Expr::new_constant(1.0));
    println!("  Unsimplified: {sum_cos}");
    let simplified_cos = simplify(&sum_cos);
    println!("  Simplified:   {simplified_cos}");
    let x = 1.0_f64;
    let fxp1 = eval_expr_single(&simplified_cos, "x", x + 1.0).unwrap();
    let fx = eval_expr_single(&simplified_cos, "x", x).unwrap();
    println!(
        "  F(2) - F(1) = {:.8}  (expected cos(1) = {:.8})",
        fxp1 - fx,
        1.0_f64.cos()
    );
    println!();

    // ------------------------------------------------------------------
    // 11. Classic result: 1^3 + 2^3 + ... + n^3 = n^2(n+1)^2/4
    //    IndefiniteSum(x^3, x, 1) → Δ⁻¹ x^3 = x(x-1)/2
    //    Σ_{k=1}^n k^3 = F(n+1) - F(1)
    // ------------------------------------------------------------------
    println!("--- 11. Definite sum: 1^3 + 2^3 + ... + n^3 = n^2(n+1)^2/4 ---");
    let body_x = Expr::Power(
        Arc::new(Expr::new_variable("x")),
        Arc::new(Expr::new_constant(3.0)),
    );
    let sum_of_k =
        Expr::new_indefinite_sum(body_x.clone(), "x".to_string(), Expr::new_constant(1.0));
    let anti_x = simplify(&sum_of_k);
    println!("  IndefiniteSum(x^3, x, 1) simplified: {anti_x}");
    // Evaluate via eval_antidiff since the simplified form is a BinaryList tag
    let anti = try_closed_form_sum(&body_x, "x").unwrap();
    let f1 = eval_antidiff(&anti, "x", 1.0).unwrap(); // F(1) = 0
    for n in [5u32, 10, 50, 100] {
        let fn1 = eval_antidiff(&anti, "x", n as f64 + 1.0).unwrap();
        let sum = fn1 - f1;
        let expected = (n.pow(2) * (n + 1).pow(2) / 4) as f64;
        println!("  Σ k^3, k=1..{n:3} = {sum:8.1}  (n^2(n+1)^2/4 = {expected:8.0})");
    }
    println!();
}
