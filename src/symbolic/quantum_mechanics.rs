//! # Quantum Mechanics
//!
//! This module provides symbolic tools for quantum mechanics, including representations
//! of quantum states (Bra-Ket notation), operators, and fundamental equations.
//!
//! ## Key Concepts
//! - **Bra-Ket Notation**: `|ψ>` (Ket) and `<ψ|` (Bra).
//! - **Operators**: Symbolic representations of physical observables.
//! - **Schrödinger Equation**: Time-independent and time-dependent versions.
//! - **Commutation Relations**: `[A, B] = AB - BA`.
//! - **Expectation Values and Uncertainty**.
//! - **Relativistic Quantum Mechanics**: Dirac and Klein-Gordon equations.

use std::sync::Arc;

use serde::Deserialize;
use serde::Serialize;

use crate::symbolic::calculus::differentiate;
use crate::symbolic::core::Expr;
use crate::symbolic::simplify_dag::simplify;
use crate::symbolic::solve::solve;

/// Represents a quantum state using Dirac notation (Ket).
///
/// Symbolically, a ket is represented as `|state>`.
#[derive(
    Clone,
    Debug,
    Serialize,
    Deserialize,
    PartialEq,
    Eq,
)]

pub struct Ket {
    pub state: Expr,
}

/// Represents a quantum state using Dirac notation (Bra).
///
/// Symbolically, a bra is represented as `<state|`.
#[derive(
    Clone,
    Debug,
    Serialize,
    Deserialize,
    PartialEq,
    Eq,
)]

pub struct Bra {
    pub state: Expr,
}

/// Computes the inner product of a Bra and a Ket, `<Bra|Ket>`.
///
/// This is a symbolic representation of the inner product over all space,
/// typically defined as `∫ ψ*(x)φ(x) dx`.
///
/// # Returns
/// An `Expr` representing `∫ bra.state * ket.state dx`.
#[must_use]

pub fn bra_ket(
    bra: &Bra,
    ket: &Ket,
) -> Expr {

    let integrand = Expr::new_mul(
        bra.state.clone(),
        ket.state.clone(),
    );

    simplify(&Expr::Integral {
        integrand: Arc::new(integrand),
        var: Arc::new(Expr::Variable(
            "x".to_string(),
        )),
        lower_bound: Arc::new(
            Expr::NegativeInfinity,
        ),
        upper_bound: Arc::new(
            Expr::Infinity,
        ),
    })
}

/// Represents a quantum operator.
///
/// Symbolically, an operator `A` acts on a state `|ψ>` as `A|ψ>`.
#[derive(
    Clone,
    Debug,
    Serialize,
    Deserialize,
    PartialEq,
    Eq,
)]

pub struct Operator {
    pub op: Expr,
}

impl Operator {
    /// Creates a new operator from an expression.
    #[must_use]

    pub const fn new(op: Expr) -> Self {

        Self {
            op,
        }
    }

    /// Applies an operator to a Ket, `O|Ket>`.
    ///
    /// # Returns
    /// A new `Ket` with state `O * ket.state`.
    #[must_use]

    pub fn apply(
        &self,
        ket: &Ket,
    ) -> Ket {

        Ket {
            state: simplify(
                &Expr::new_mul(
                    self.op.clone(),
                    ket.state.clone(),
                ),
            ),
        }
    }
}

/// Computes the commutator of two operators: `[A, B] = AB - BA`.
///
/// When applied to a state `|ψ>`, it returns `A(B|ψ>) - B(A|ψ>)`.
#[must_use]

pub fn commutator(
    a: &Operator,
    b: &Operator,
    ket: &Ket,
) -> Expr {

    let ab_psi = a.apply(&b.apply(ket));

    let ba_psi = b.apply(&a.apply(ket));

    simplify(&Expr::new_sub(
        ab_psi.state,
        ba_psi.state,
    ))
}

/// Computes the expectation value of an operator: `<A> = <ψ|A|ψ> / <ψ|ψ>`.
#[must_use]

pub fn expectation_value(
    op: &Operator,
    psi: &Ket,
) -> Expr {

    let bra = Bra {
        state: psi.state.clone(),
    };

    let numerator =
        bra_ket(&bra, &op.apply(psi));

    let denominator =
        bra_ket(&bra, psi);

    simplify(&Expr::new_div(
        numerator,
        denominator,
    ))
}

/// Computes the uncertainty (standard deviation) of an operator: `ΔA = sqrt(<A^2> - <A>^2)`.
#[must_use]

pub fn uncertainty(
    op: &Operator,
    psi: &Ket,
) -> Expr {

    let op_sq = Operator {
        op: Expr::new_pow(
            op.op.clone(),
            Expr::Constant(2.0),
        ),
    };

    let exp_a_sq =
        expectation_value(&op_sq, psi);

    let exp_a =
        expectation_value(op, psi);

    let exp_a_whole_sq = Expr::new_pow(
        exp_a,
        Expr::Constant(2.0),
    );

    simplify(&Expr::new_sqrt(
        Expr::new_sub(
            exp_a_sq,
            exp_a_whole_sq,
        ),
    ))
}

/// Computes the probability density at a point: `ρ(x) = |ψ(x)|^2`.
#[must_use]

pub fn probability_density(
    psi: &Ket
) -> Expr {

    simplify(&Expr::new_pow(
        Expr::new_abs(
            psi.state.clone(),
        ),
        Expr::Constant(2.0),
    ))
}

/// Hamiltonian for a free particle: `H = -ħ² / (2m) * ∇²`.
#[must_use]

pub fn hamiltonian_free_particle(
    m: &Expr
) -> Operator {

    let hbar =
        Expr::new_variable("hbar");

    let hbar_sq = Expr::new_pow(
        hbar,
        Expr::Constant(2.0),
    );

    let two_m = Expr::new_mul(
        Expr::Constant(2.0),
        m.clone(),
    );

    let coeff = Expr::new_neg(
        Expr::new_div(hbar_sq, two_m),
    );

    let laplacian =
        Expr::new_variable("d2_dx2");

    Operator {
        op: simplify(&Expr::new_mul(
            coeff,
            laplacian,
        )),
    }
}

/// Hamiltonian for a harmonic oscillator: `H = -ħ² / (2m) * ∇² + 1/2 * m * ω² * x²`.
#[must_use]

pub fn hamiltonian_harmonic_oscillator(
    m: &Expr,
    omega: &Expr,
) -> Operator {

    let free_h =
        hamiltonian_free_particle(m);

    let half = Expr::Constant(0.5);

    let omega_sq = Expr::new_pow(
        omega.clone(),
        Expr::Constant(2.0),
    );

    let x = Expr::new_variable("x");

    let x_sq = Expr::new_pow(
        x,
        Expr::Constant(2.0),
    );

    let potential = Expr::new_mul(
        half,
        Expr::new_mul(
            m.clone(),
            Expr::new_mul(
                omega_sq,
                x_sq,
            ),
        ),
    );

    Operator {
        op: simplify(&Expr::new_add(
            free_h.op,
            potential,
        )),
    }
}

/// Angular momentum operator `L_z`: `L_z = -i * ħ * ∂/∂φ`.
#[must_use]

pub fn angular_momentum_z() -> Operator
{

    let i = Expr::new_complex(
        Expr::Constant(0.0),
        Expr::Constant(1.0),
    );

    let hbar =
        Expr::new_variable("hbar");

    let i_hbar = Expr::new_mul(i, hbar);

    let d_dphi =
        Expr::new_variable("d_dphi");

    Operator {
        op: simplify(&Expr::new_neg(
            Expr::new_mul(
                i_hbar,
                d_dphi,
            ),
        )),
    }
}

/// Returns the Pauli matrices: `σ_x, σ_y, σ_z`.
#[must_use]

pub fn pauli_matrices(
) -> (Expr, Expr, Expr) {

    let zero = Expr::Constant(0.0);

    let one = Expr::Constant(1.0);

    let i = Expr::new_complex(
        Expr::Constant(0.0),
        Expr::Constant(1.0),
    );

    let neg_i =
        Expr::new_neg(i.clone());

    let sigma_x = Expr::Matrix(vec![
        vec![
            zero.clone(),
            one.clone(),
        ],
        vec![
            one.clone(),
            zero.clone(),
        ],
    ]);

    let sigma_y = Expr::Matrix(vec![
        vec![zero.clone(), neg_i],
        vec![i, zero.clone()],
    ]);

    let sigma_z = Expr::Matrix(vec![
        vec![one, zero.clone()],
        vec![
            zero,
            Expr::Constant(-1.0),
        ],
    ]);

    (
        sigma_x,
        sigma_y,
        sigma_z,
    )
}

/// Spin operator: `S = ħ/2 * σ`.
#[must_use]

pub fn spin_operator(
    pauli: &Expr
) -> Expr {

    let hbar =
        Expr::new_variable("hbar");

    let half_hbar = Expr::new_mul(
        Expr::Constant(0.5),
        hbar,
    );

    simplify(&Expr::new_mul(
        half_hbar,
        pauli.clone(),
    ))
}

/// Solves the time-independent Schrödinger equation: `H|ψ> = E|ψ>`.
#[must_use]

pub fn solve_time_independent_schrodinger(
    hamiltonian: &Operator,
    wave_function: &Ket,
) -> (Vec<Expr>, Vec<Ket>) {

    let h_psi = hamiltonian
        .apply(wave_function);

    let e =
        Expr::Variable("E".to_string());

    let e_psi = Expr::new_mul(
        e,
        wave_function
            .state
            .clone(),
    );

    let equation = Expr::new_sub(
        h_psi.state,
        e_psi,
    );

    let solutions =
        solve(&equation, "E");

    let eigenfunctions = solutions
        .iter()
        .map(|_sol| {

            wave_function.clone()
        })
        .collect();

    (
        solutions,
        eigenfunctions,
    )
}

/// Time-dependent Schrödinger equation: `iħ ∂/∂t |ψ> = H|ψ>`.
#[must_use]

pub fn time_dependent_schrodinger_equation(
    hamiltonian: &Operator,
    wave_function: &Ket,
) -> Expr {

    let i = Expr::new_complex(
        Expr::Constant(0.0),
        Expr::Constant(1.0),
    );

    let hbar = Expr::Variable(
        "hbar".to_string(),
    );

    let i_hbar = Expr::new_mul(i, hbar);

    let d_psi_dt = differentiate(
        &wave_function.state,
        "t",
    );

    let lhs =
        Expr::new_mul(i_hbar, d_psi_dt);

    let rhs = hamiltonian
        .apply(wave_function)
        .state;

    simplify(&Expr::new_sub(
        lhs, rhs,
    ))
}

/// Dirac equation for a free particle: `(iħ γ^μ ∂_μ - mc)ψ = 0`.
#[must_use]

pub fn dirac_equation(
    psi: &Expr,
    m: &Expr,
) -> Expr {

    let i = Expr::new_complex(
        Expr::Constant(0.0),
        Expr::Constant(1.0),
    );

    let hbar =
        Expr::new_variable("hbar");

    let c = Expr::new_variable("c");

    let gamma_mu =
        Expr::new_variable("gamma_mu");

    let d_mu = Expr::new_variable(
        "partial_mu",
    );

    let term1 = Expr::new_mul(
        i,
        Expr::new_mul(
            hbar,
            Expr::new_mul(
                gamma_mu,
                d_mu,
            ),
        ),
    );

    let term2 =
        Expr::new_mul(m.clone(), c);

    let dirac_op =
        Expr::new_sub(term1, term2);

    simplify(&Expr::new_mul(
        dirac_op,
        psi.clone(),
    ))
}

/// Klein-Gordon equation: `(∂^μ ∂_μ + (mc/ħ)²)ψ = 0`.
#[must_use]

pub fn klein_gordon_equation(
    psi: &Expr,
    m: &Expr,
) -> Expr {

    let hbar =
        Expr::new_variable("hbar");

    let c = Expr::new_variable("c");

    let dalembertian =
        Expr::new_variable("d_mu_d_mu");

    let mc_hbar = Expr::new_div(
        Expr::new_mul(m.clone(), c),
        hbar,
    );

    let mass_term = Expr::new_pow(
        mc_hbar,
        Expr::Constant(2.0),
    );

    let kg_op = Expr::new_add(
        dalembertian,
        mass_term,
    );

    simplify(&Expr::new_mul(
        kg_op,
        psi.clone(),
    ))
}

/// Computes the first-order energy correction in perturbation theory: `E^(1) = <ψ^(0)|H'|ψ^(0)>`.
#[must_use]

pub fn first_order_energy_correction(
    perturbation: &Operator,
    unperturbed_state: &Ket,
) -> Expr {

    bra_ket(
        &Bra {
            state: unperturbed_state
                .state
                .clone(),
        },
        &perturbation
            .apply(unperturbed_state),
    )
}

/// Scattering amplitude in quantum mechanics: `f(θ, φ) ∝ <φ|V|ψ>`.
#[must_use]

pub fn scattering_amplitude(
    initial_state: &Ket,
    final_state: &Ket,
    potential: &Operator,
) -> Expr {

    let term =
        potential.apply(initial_state);

    bra_ket(
        &Bra {
            state: final_state
                .state
                .clone(),
        },
        &term,
    )
}
