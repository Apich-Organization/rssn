//! Numerical indefinite sum (anti-difference) examples.
//!
//! Demonstrates the three evaluation strategies available in
//! `rssn::numerical::indefinite_sum`:
//!
//!   1. Closed-form rules (constant, exponential, trigonometric, power)
//!   2. Unified pipeline (`eval_indefinite_sum_numerical`)
//!   3. Series-expansion fallback (`series_antidiff`)
//!   4. Indefinite product via log-sum reduction

use rssn::numerical::indefinite_sum::IndefiniteSumConfig;
use rssn::numerical::indefinite_sum::compute_taylor_coeffs_numerical;
use rssn::numerical::indefinite_sum::eval_antidiff;
use rssn::numerical::indefinite_sum::eval_indefinite_product_numerical;
use rssn::numerical::indefinite_sum::eval_indefinite_sum_numerical;
use rssn::numerical::indefinite_sum::series_antidiff;
use rssn::numerical::indefinite_sum::try_closed_form_sum;
use rssn::numerical::special::hurwitz_zeta;
use rssn::symbolic::core::Expr;
use std::sync::Arc;

fn main() {
    println!("=== Numerical Indefinite Sum Examples ===\n");

    // ------------------------------------------------------------------
    // 1. Closed-form: constant  Δ⁻¹ c = c·x
    // ------------------------------------------------------------------
    println!("--- 1. Closed-form: constant ---");
    let f_const = Expr::new_constant(5.0);
    let anti_const = try_closed_form_sum(&f_const, "x").unwrap();
    // F(x+1) - F(x) should equal 5
    let f4 = eval_antidiff(&anti_const, "x", 4.0).unwrap();
    let f3 = eval_antidiff(&anti_const, "x", 3.0).unwrap();
    println!("  f(x) = 5");
    println!("  F(4) - F(3) = {} (expected 5)", f4 - f3);
    println!();

    // ------------------------------------------------------------------
    // 2. Closed-form: exponential  Δ⁻¹ eˣ = eˣ / (e - 1)
    // ------------------------------------------------------------------
    println!("--- 2. Closed-form: exponential ---");
    let f_exp = Expr::Exp(Arc::new(Expr::new_variable("x")));
    let anti_exp = try_closed_form_sum(&f_exp, "x").unwrap();
    for x in [1.0_f64, 2.0, 3.0] {
        let fxp1 = eval_antidiff(&anti_exp, "x", x + 1.0).unwrap();
        let fx = eval_antidiff(&anti_exp, "x", x).unwrap();
        println!(
            "  F({:.0}+1) - F({:.0}) = {:.6}  (expected e^{:.0} = {:.6})",
            x,
            x,
            fxp1 - fx,
            x,
            x.exp()
        );
    }
    println!();

    // ------------------------------------------------------------------
    // 3. Closed-form: trigonometric  Δ⁻¹ sin(ax)
    // ------------------------------------------------------------------
    println!("--- 3. Closed-form: sin(ax) ---");
    let a = 0.8_f64;
    let f_sin = Expr::Sin(Arc::new(Expr::new_mul(
        Expr::new_constant(a),
        Expr::new_variable("x"),
    )));
    let anti_sin = try_closed_form_sum(&f_sin, "x").unwrap();
    let x = 2.5_f64;
    let fxp1 = eval_antidiff(&anti_sin, "x", x + 1.0).unwrap();
    let fx = eval_antidiff(&anti_sin, "x", x).unwrap();
    println!("  f(x) = sin({a}·x),  x = {x}",);
    println!(
        "  F(x+1) - F(x) = {:.8}  (expected sin({:.1}·{x}) = {:.8})",
        fxp1 - fx,
        a,
        (a * x).sin()
    );
    println!();

    // ------------------------------------------------------------------
    // 4. Closed-form: power  Δ⁻¹ xⁿ = ζ(−n, 1) − ζ(−n, x)
    // ------------------------------------------------------------------
    println!("--- 4. Closed-form: x^n (Hurwitz zeta) ---");
    let f_pow = Expr::Power(
        Arc::new(Expr::new_variable("x")),
        Arc::new(Expr::new_constant(3.0)),
    );
    let anti_pow = try_closed_form_sum(&f_pow, "x").unwrap();
    for x in [2.0_f64, 4.0, 6.0] {
        let fxp1 = eval_antidiff(&anti_pow, "x", x + 1.0).unwrap();
        let fx = eval_antidiff(&anti_pow, "x", x).unwrap();
        println!(
            "  F({:.0}+1) - F({:.0}) = {:.2}  (expected {:.0}^3 = {:.0})",
            x,
            x,
            fxp1 - fx,
            x,
            x.powi(3)
        );
    }
    println!();

    // ------------------------------------------------------------------
    // 5. Unified pipeline: eval_indefinite_sum_numerical
    //    Tries closed-form → Abel-Plana → series in order
    // ------------------------------------------------------------------
    println!("--- 5. Unified pipeline: eval_indefinite_sum_numerical ---");
    let config = IndefiniteSumConfig::default(); // h = 0, step = 1
    let f_exp2 = Expr::Exp(Arc::new(Expr::new_variable("x")));
    for x in [1.0_f64, 2.0, 5.0] {
        let fxp1 = eval_indefinite_sum_numerical(&f_exp2, "x", x + 1.0, &config).unwrap();
        let fx = eval_indefinite_sum_numerical(&f_exp2, "x", x, &config).unwrap();
        println!(
            "  [e^x] F({:.0}+1) - F({:.0}) = {:.6}  (expected {:.6})",
            x,
            x,
            fxp1 - fx,
            x.exp()
        );
    }
    println!();

    // ------------------------------------------------------------------
    // 6. Taylor coefficients via forward differences
    // ------------------------------------------------------------------
    println!("--- 6. Taylor coefficients (forward differences) ---");
    let f_quad = Expr::Power(
        Arc::new(Expr::new_variable("x")),
        Arc::new(Expr::new_constant(2.0)),
    );
    let coeffs = compute_taylor_coeffs_numerical(&f_quad, "x", 1.0, 4);
    println!("  f(x) = x^2 expanded around x = 1:");
    println!("    c_0 = {:.6}  (expected 1.0)", coeffs[0]);
    println!("    c_1 = {:.6}  (expected 2.0)", coeffs[1]);
    println!("    c_2 = {:.6}  (expected 1.0)", coeffs[2]);
    println!();

    // ------------------------------------------------------------------
    // 7. Series anti-difference with exact Taylor coefficients
    //    Backward convention: G(z+1) - G(z) = f(z+1)
    // ------------------------------------------------------------------
    println!("--- 7. Series anti-difference (exact coefficients for e^x) ---");
    let mut fact = 1.0_f64;
    let exact_coeffs: Vec<f64> = (0..6)
        .map(|m| {
            if m > 0 {
                fact *= m as f64;
            }
            1.0 / fact
        })
        .collect();
    let h = 0.0_f64;
    println!(
        "  Normalization: G(h=0) = {:.2e}",
        series_antidiff(&exact_coeffs, h, h)
    );
    let gz0 = series_antidiff(&exact_coeffs, h, 0.0);
    let gz1 = series_antidiff(&exact_coeffs, h, 1.0);
    println!(
        "  G(1) - G(0) = {:.6}  (expected e^1 = {:.6})",
        gz1 - gz0,
        std::f64::consts::E
    );
    println!();

    // ------------------------------------------------------------------
    // 8. Hurwitz zeta raw recurrence identity (mathematical check)
    // ------------------------------------------------------------------
    println!("--- 8. Hurwitz zeta identity: −ζ(−2, x+1) ---");
    println!("  This formula satisfies F(x+1) - F(x) = (x+1)^2");
    for x in [1.0_f64, 3.0, 5.0] {
        let fx = -hurwitz_zeta(-2.0, x + 1.0);
        let fxp1 = -hurwitz_zeta(-2.0, x + 2.0);
        println!(
            "  x={:.0}: F(x+1)-F(x) = {:.2}  (expected (x+1)^2 = {:.0})",
            x,
            fxp1 - fx,
            (x + 1.0).powi(2)
        );
    }
    println!();

    // ------------------------------------------------------------------
    // 9. Classic result: 1 + 2 + ... + n = n(n+1)/2
    //    Δ⁻¹ x = x(x-1)/2  →  Σ_{k=1}^n k = F(n+1) - F(1)
    // ------------------------------------------------------------------
    println!("--- 9. Definite sum: 1 + 2 + ... + n = n(n+1)/2 ---");
    let f_linear = Expr::Power(
        Arc::new(Expr::new_variable("x")),
        Arc::new(Expr::new_constant(1.0)),
    );
    let anti_linear = try_closed_form_sum(&f_linear, "x").unwrap();
    println!("  Δ⁻¹ x = {anti_linear}");
    for n in [5u32, 10, 100] {
        // Σ_{k=1}^n k = F(n+1) - F(1)
        let fn1 = eval_antidiff(&anti_linear, "x", n as f64 + 1.0).unwrap();
        let f1 = eval_antidiff(&anti_linear, "x", 1.0).unwrap();
        let sum = fn1 - f1;
        let expected = (n * (n + 1) / 2) as f64;
        println!("  Σ k, k=1..{n:3} = {sum:8.1}  (formula n(n+1)/2 = {expected:8.0})");
    }
    println!();

    // ------------------------------------------------------------------
    // 10. Indefinite product  ∏ f(k) via log-sum reduction
    // ------------------------------------------------------------------
    println!("--- 10. Indefinite product: prod_{{k=h}}^{{x-1}} f(k) ---");
    // f(k) = 2  →  product from k=h to x-1 is 2^(x-h)
    let f_prod = |_k: f64| -> Result<f64, String> { Ok(2.0_f64) };
    let h_prod = 1.0_f64;
    for x in [2.0_f64, 4.0, 6.0] {
        let px = eval_indefinite_product_numerical(x, &f_prod, h_prod, 1.0).unwrap();
        println!(
            "  P({:.0}) = {:.4}  (expected 2^{:.0} = {:.4})",
            x,
            px,
            x - h_prod,
            2.0_f64.powf(x - h_prod)
        );
    }
    println!();

    // ------------------------------------------------------------------
    // 11. Classic result: 1^3 + 2^3 + ... + n^3 = n^2(n+1)^2/4
    //    Δ⁻¹ x^3 = x(x-1)/2  →  Σ_{k=1}^n k = F(n+1) - F(1)
    // ------------------------------------------------------------------
    println!("--- 11. Definite sum: 1^3 + 2^3 + ... + n^3 = n^2(n+1)^2/4 ---");
    // let f_linear = Expr::Power(
    //     Arc::new(Expr::new_variable("x")),
    //     Arc::new(Expr::new_constant(3.0)),
    // );
    let f_linear = Expr::new_pow(Expr::new_variable("x"), Expr::new_constant(3.0));
    let anti_linear = try_closed_form_sum(&f_linear, "x").unwrap();
    println!("  Δ⁻¹ x^3 = {anti_linear}");
    for n in [5u32, 10, 100] {
        // Σ_{k=1}^n k^3 = F(n+1) - F(1)
        let fn1 = eval_antidiff(&anti_linear, "x", n as f64 + 1.0).unwrap();
        let f1 = eval_antidiff(&anti_linear, "x", 1.0).unwrap();
        let sum = fn1 - f1;
        let expected = (n.pow(2) * (n + 1).pow(2) / 4) as f64;
        println!("  Σ k^3, k=1..{n:3} = {sum:8.1}  (formula n^2(n+1)^2/4 = {expected:8.0})");
    }
    println!();
}
