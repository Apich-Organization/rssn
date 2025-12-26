//! # Fractal Geometry and Chaos Theory
//!
//! This module provides symbolic tools for exploring concepts in fractal geometry
//! and chaos theory. It includes representations for Iterated Function Systems (IFS),
//! complex dynamical systems (Mandelbrot/Julia sets), and tools for analyzing chaotic
//! behavior such as fixed points, stability, and Lyapunov exponents.

use crate::symbolic::calculus::{
    differentiate,
    substitute,
};
use crate::symbolic::core::Expr;
use crate::symbolic::simplify_dag::simplify;
use crate::symbolic::solve::solve;
use serde::{
    Deserialize,
    Serialize,
};
use std::sync::Arc;

// ============================================================================
// Iterated Function Systems (IFS)
// ============================================================================

/// Represents an Iterated Function System (IFS).
///
/// An IFS is a finite set of contraction mappings on a complete metric space.
/// It is often used to construct fractals (e.g., Sierpinski triangle, Barnsley fern).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]

pub struct IteratedFunctionSystem {
    /// The set of contraction mappings. Each expression should be a function of the coordinates.
    /// For 2D, these are typically affine transformations.
    pub functions: Vec<Expr>,
    /// The probabilities associated with each function (for randomized algorithms).
    /// These should sum to 1.
    pub probabilities: Vec<Expr>,
    /// The variables involved in the functions (e.g., ["x", "y"]).
    pub variables: Vec<String>,
}

impl IteratedFunctionSystem {
    /// Creates a new Iterated Function System.
    #[must_use]

    pub const fn new(
        functions: Vec<Expr>,
        probabilities: Vec<Expr>,
        variables: Vec<String>,
    ) -> Self {

        Self {
            functions,
            probabilities,
            variables,
        }
    }

    /// Applies the IFS to a point (symbolically) to generate the set of possible next points.
    ///
    /// # Arguments
    /// * `point` - A vector of expressions representing the coordinates of the point.
    ///
    /// # Returns
    /// A vector of points (vectors of expressions), one for each function in the IFS.
    #[must_use]

    pub fn apply(
        &self,
        point: &[Expr],
    ) -> Vec<Vec<Expr>> {

        if point.len() != self.variables.len() {

            // In a real scenario, return Result. For now, panic or return empty.
            return vec![];
        }

        let mut results = Vec::new();

        for func in &self.functions {

            // Assuming func is a vector/list expression [f_1, f_2, ...] or we handle
            // the transformation logic differently.
            // Actually, usually an IFS in 2D is defined by matrices or a list of equations.
            // Let's assume `func` is a List/Vector expression containing the new coordinates.

            let mut new_point = Vec::new();

            if let Expr::Vector(coords) = func {

                for coord_expr in coords {

                    let mut substituted = coord_expr.clone();

                    for (i, var) in self
                        .variables
                        .iter()
                        .enumerate()
                    {

                        substituted = substitute(
                            &substituted,
                            var,
                            &point[i],
                        );
                    }

                    new_point.push(simplify(
                        &substituted,
                    ));
                }
            } else {

                // If func is not a list, maybe it's a 1D IFS?
                let mut substituted = func.clone();

                for (i, var) in self
                    .variables
                    .iter()
                    .enumerate()
                {

                    substituted = substitute(
                        &substituted,
                        var,
                        &point[i],
                    );
                }

                new_point.push(simplify(
                    &substituted,
                ));
            }

            results.push(new_point);
        }

        results
    }

    /// Calculates the similarity dimension for a self-similar IFS.
    ///
    /// If the IFS consists of `N` similarities with scaling ratios `r_1, ..., r_N`,
    /// the similarity dimension `D` is the unique solution to `sum(r_i^D) = 1`.
    ///
    /// This function assumes the `probabilities` field actually holds the scaling factors `r_i`
    /// if they are provided, or it requires the user to provide scaling factors.
    ///
    /// For simplicity here, we accept a list of scaling factors.
    #[must_use]

    pub fn similarity_dimension(scaling_factors: &[Expr]) -> Expr {

        // We need to solve sum(r_i^D) = 1 for D.
        // This is generally transcendental.
        // If all r_i are equal to r, then N * r^D = 1 => D = - log(N) / log(r).

        // Check if all scaling factors are the same
        let first = &scaling_factors[0];

        let all_same = scaling_factors
            .iter()
            .all(|r| r == first);

        if all_same {

            let n = Expr::Constant(scaling_factors.len() as f64);

            let r = first.clone();

            // D = log(N) / log(1/r) = - log(N) / log(r)
            let num = Expr::new_neg(Expr::new_log(n));

            let den = Expr::new_log(r);

            simplify(&Expr::new_div(
                num, den,
            ))
        } else {

            // Return the equation sum(r_i^D) = 1
            let d = Expr::Variable("D".to_string());

            let mut sum = Expr::Constant(0.0);

            for r in scaling_factors {

                sum = Expr::new_add(
                    sum,
                    Expr::new_pow(r.clone(), d.clone()),
                );
            }

            Expr::Eq(
                Arc::new(sum),
                Arc::new(Expr::Constant(1.0)),
            )
        }
    }
}

// ============================================================================
// Complex Dynamical Systems (Mandelbrot / Julia)
// ============================================================================

/// Represents a complex dynamical system defined by z_{n+1} = `f(z_n)` + c.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]

pub struct ComplexDynamicalSystem {
    /// The function f(z) (e.g., z^2).
    pub function: Expr,
    /// The parameter c.
    pub c: Expr,
}

impl ComplexDynamicalSystem {
    /// Creates a new Mandelbrot/Julia system z -> z^2 + c.
    #[must_use]

    pub fn new_mandelbrot_family(c: Expr) -> Self {

        // f(z) = z^2
        let z = Expr::Variable("z".to_string());

        let f = Expr::new_pow(
            z,
            Expr::Constant(2.0),
        );

        Self { function: f, c }
    }

    /// Iterates the system once: z_{n+1} = `f(z_n)` + c.
    #[must_use]

    pub fn iterate(
        &self,
        z_n: &Expr,
    ) -> Expr {

        let f_z = substitute(
            &self.function,
            "z",
            z_n,
        );

        simplify(&Expr::new_add(
            f_z,
            self.c.clone(),
        ))
    }

    /// Computes the orbit of a point up to n iterations.
    #[must_use]

    pub fn orbit(
        &self,
        start_z: Expr,
        n: usize,
    ) -> Vec<Expr> {

        let mut orbit = Vec::with_capacity(n + 1);

        orbit.push(start_z.clone());

        let mut current = start_z;

        for _ in 0..n {

            current = self.iterate(&current);

            orbit.push(current.clone());
        }

        orbit
    }

    /// Finds fixed points of the system: z = f(z) + c.
    #[must_use]

    pub fn fixed_points(&self) -> Vec<Expr> {

        let z = Expr::Variable("z".to_string());

        // Solve z = f(z) + c  =>  f(z) + c - z = 0
        let rhs = Expr::new_add(
            self.function
                .clone(),
            self.c.clone(),
        );

        let eq = Expr::new_sub(rhs, z);

        // Use the solver
        solve(&eq, "z")
    }

    /// Checks the stability of a fixed point z*.
    /// Stability is determined by |(f(z) + c)'| at z*.
    /// If modulus < 1, stable (attracting).
    /// If modulus > 1, unstable (repelling).
    /// Returns the symbolic magnitude of the derivative.
    #[must_use]

    pub fn stability_index(
        &self,
        fixed_point: &Expr,
    ) -> Expr {

        let map = Expr::new_add(
            self.function
                .clone(),
            self.c.clone(),
        );

        let deriv = differentiate(&map, "z");

        let val = substitute(
            &deriv,
            "z",
            fixed_point,
        );

        simplify(&Expr::new_abs(val))
    }
}

// ============================================================================
// Chaos Theory Tools
// ============================================================================

/// Calculates the fixed points of a 1D map f(x).
/// Solves f(x) = x.
#[must_use]

pub fn find_fixed_points(
    map_function: &Expr,
    var: &str,
) -> Vec<Expr> {

    let x = Expr::Variable(var.to_string());

    // f(x) - x = 0
    let eq = Expr::new_sub(
        map_function.clone(),
        x,
    );

    solve(&eq, var)
}

/// Analyzes the stability of a fixed point for a 1D map f(x).
/// Returns the derivative evaluated at the fixed point: f'(x*).
/// |f'(x*)| < 1 => Stable.
#[must_use]

pub fn analyze_stability(
    map_function: &Expr,
    var: &str,
    fixed_point: &Expr,
) -> Expr {

    let deriv = differentiate(map_function, var);

    let val = substitute(
        &deriv,
        var,
        fixed_point,
    );

    simplify(&val)
}

/// Calculates the symbolic Lyapunov exponent for a 1D chaotic map `x_{n+1} = f(x_n)`.
///
/// The Lyapunov exponent `λ` is defined as:
/// `λ = lim (n->inf) (1/n) * sum_{i=0}^{n-1} ln(|f'(x_i)|)`
///
/// # Arguments
/// * `map_function` - The symbolic expression for `f(x)`.
/// * `var` - The variable name (e.g., "x").
/// * `initial_x` - The initial value `x_0`.
/// * `n_iterations` - The number of iterations to sum.
///
/// # Returns
/// An `Expr` representing the approximate Lyapunov exponent after `n` iterations.
#[must_use]

pub fn lyapunov_exponent(
    map_function: &Expr,
    var: &str,
    initial_x: &Expr,
    n_iterations: usize,
) -> Expr {

    let mut current_x = initial_x.clone();

    let mut sum_log_derivs = Expr::Constant(0.0);

    // Pre-calculate derivative function
    let deriv_func = differentiate(map_function, var);

    for _ in 0..n_iterations {

        // Evaluate f'(current_x)
        let deriv_val = substitute(
            &deriv_func,
            var,
            &current_x,
        );

        // ln(|f'(x)|)
        let log_abs = Expr::new_log(Expr::new_abs(
            deriv_val,
        ));

        sum_log_derivs = Expr::new_add(
            sum_log_derivs,
            log_abs,
        );

        // Update x: x_{n+1} = f(x_n)
        current_x = substitute(
            map_function,
            var,
            &current_x,
        );

        // Simplify occasionally to prevent expression explosion?
        // For purely symbolic, this might grow huge.
        // We'll simplify at each step.
        current_x = simplify(&current_x);

        sum_log_derivs = simplify(&sum_log_derivs);
    }

    let n = Expr::Constant(n_iterations as f64);

    simplify(&Expr::new_div(
        sum_log_derivs,
        n,
    ))
}

/// Returns the Lorenz System equations.
///
/// dx/dt = sigma * (y - x)
/// dy/dt = x * (rho - z) - y
/// dz/dt = x * y - beta * z
#[must_use]

pub fn lorenz_system() -> (Expr, Expr, Expr) {

    let x = Expr::Variable("x".to_string());

    let y = Expr::Variable("y".to_string());

    let z = Expr::Variable("z".to_string());

    let sigma = Expr::Variable("sigma".to_string());

    let rho = Expr::Variable("rho".to_string());

    let beta = Expr::Variable("beta".to_string());

    let dx = Expr::new_mul(
        sigma,
        Expr::new_sub(y.clone(), x.clone()),
    );

    let dy_term1 = Expr::new_mul(
        x.clone(),
        Expr::new_sub(rho, z.clone()),
    );

    let dy = Expr::new_sub(dy_term1, y.clone());

    let dz = Expr::new_sub(
        Expr::new_mul(x, y),
        Expr::new_mul(beta, z),
    );

    (dx, dy, dz)
}
