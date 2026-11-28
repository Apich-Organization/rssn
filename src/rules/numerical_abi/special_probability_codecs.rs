binding!(
    harmonic,
    "num_harmonic",
    "numerical_special",
    "`num_harmonic(n)` returns the floating-point harmonic number H_n for a nonnegative u64 integer. H_0 is zero.",
    "crate::kernels::special::harmonic",
    DimensionPolicy::Dimensionless,
    (n: u64) -> f64,
    special::harmonic
);
binding!(
    log_binomial,
    "num_ln_binomial",
    "numerical_special",
    "`num_ln_binomial(n, k)` returns the natural logarithm of the binomial coefficient on 0 <= k <= n. The result stays finite when the coefficient would overflow binary64.",
    "crate::kernels::special::ln_binomial",
    DimensionPolicy::Dimensionless,
    (n: u64, k: u64) -> Result<f64, OpError>,
    |n: u64, k: u64| {
        if k > n {
            return Err(OpError::Domain("a finite log binomial requires k <= n".into()));
        }
        Ok(special::ln_binomial(n, k))
    }
);
binding!(
    checked_multinomial,
    "num_multinomial",
    "numerical_special",
    "`num_multinomial(n, parts)` returns the floating-point multinomial coefficient. Parts are nonnegative integers and must sum exactly to n. The empty partition of zero gives one. A coefficient that exceeds binary64 produces a domain error.",
    "crate::kernels::special::checked_multinomial",
    DimensionPolicy::Dimensionless,
    (n: u64, parts: Vec<u64>) -> Result<f64, OpError>,
    |n: u64, parts: Vec<u64>| {
        special::checked_multinomial(n, &parts)
            .ok_or_else(|| OpError::Domain("multinomial parts must sum to n".into()))
    }
);
binding!(
    inverse_beta_reg,
    "num_inverse_beta_reg",
    "numerical_special",
    "`num_inverse_beta_reg(a, b, p)` returns the regularized beta quantile at a probability in [0, 1]. Shape parameters must be finite and positive. The general inversion uses a bounded search and reports nonrepresentable forward evaluations as domain errors.",
    "crate::kernels::special::inverse_beta_reg",
    DimensionPolicy::Dimensionless,
    (a: f64, b: f64, p: f64) -> Result<f64, OpError>,
    |a: f64, b: f64, p: f64| special::inverse_beta_reg(a, b, p).map_err(OpError::Domain)
);
