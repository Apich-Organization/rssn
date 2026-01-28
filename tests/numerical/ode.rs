#[cfg(test)]
mod tests {
    use rssn::numerical::ode::solve_ode_system_rk4;
    use rssn::symbolic::core::Expr;
    use std::sync::Arc;

    #[test]
    fn test_rk4_overflow() {
        let funcs = vec![Expr::Power(Arc::new(Expr::Variable("y0".to_string())), Arc::new(Expr::new_constant(2.0)))];
        let y0 = &[1.0];
        let x_range = (0.0, 2.0); // Go past the singularity at x=1
        let num_steps = 100;
        let result = solve_ode_system_rk4(&funcs, y0, x_range, num_steps);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "Overflow or invalid value encountered during ODE solving.");
    }
}
