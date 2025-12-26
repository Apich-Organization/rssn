//! # Thermodynamics Module
//!
//! This module provides symbolic representations of fundamental laws,
//! potentials, and distributions in thermodynamics and statistical mechanics.
//!
//! ## Key Concepts:
//! - **The Laws of Thermodynamics**: Conservation of energy, entropy, and absolute zero.
//! - **Thermodynamical Potentials**: Enthalpy, Helmholtz free energy, and Gibbs free energy.
//! - **Ideal Gas Law**: $PV = nRT$.
//! - **Statistical Distributions**: Boltzmann, Fermi-Dirac, and Bose-Einstein.

use crate::symbolic::calculus::differentiate;
use crate::symbolic::core::Expr;
use crate::symbolic::simplify_dag::simplify;
use std::sync::Arc;

/// Represents the First Law of Thermodynamics: $dU = dQ - dW$.
///
/// States that the change in internal energy $dU$ of a closed system is equal to the
/// heat $dQ$ added to the system minus the work $dW$ done by the system.
#[must_use]

pub fn first_law_thermodynamics(
    internal_energy_change: &Expr,
    heat_added: &Expr,
    work_done: &Expr,
) -> Expr {

    simplify(&Expr::new_sub(
        internal_energy_change.clone(),
        Expr::new_sub(
            heat_added.clone(),
            work_done.clone(),
        ),
    ))
}

/// Represents the Ideal Gas Law: $PV = nRT$.
///
/// This function can solve for any one of the variables ($P, V, n, R, T$) if the others are provided.
/// Here we return the expression $(PV) - (nRT)$ which should equal zero.
#[must_use]

pub fn ideal_gas_law(
    p: &Expr,
    v: &Expr,
    n: &Expr,
    r: &Expr,
    t: &Expr,
) -> Expr {

    simplify(&Expr::new_sub(
        Expr::new_mul(p.clone(), v.clone()),
        Expr::new_mul(
            n.clone(),
            Expr::new_mul(r.clone(), t.clone()),
        ),
    ))
}

/// Calculates Enthalpy: $H = U + PV$.
#[must_use]

pub fn enthalpy(
    internal_energy: &Expr,
    pressure: &Expr,
    volume: &Expr,
) -> Expr {

    simplify(&Expr::new_add(
        internal_energy.clone(),
        Expr::new_mul(
            pressure.clone(),
            volume.clone(),
        ),
    ))
}

/// Calculates Helmholtz Free Energy: $A = U - TS$.
#[must_use]

pub fn helmholtz_free_energy(
    internal_energy: &Expr,
    temperature: &Expr,
    entropy: &Expr,
) -> Expr {

    simplify(&Expr::new_sub(
        internal_energy.clone(),
        Expr::new_mul(
            temperature.clone(),
            entropy.clone(),
        ),
    ))
}

/// Calculates Gibbs Free Energy: $G = H - TS$.
#[must_use]

pub fn gibbs_free_energy(
    enthalpy: &Expr,
    temperature: &Expr,
    entropy: &Expr,
) -> Expr {

    simplify(&Expr::new_sub(
        enthalpy.clone(),
        Expr::new_mul(
            temperature.clone(),
            entropy.clone(),
        ),
    ))
}

/// Calculates Entropy via Boltzmann's formula: $S = `k_B` \ln(\Omega)$.
#[must_use]

pub fn boltzmann_entropy(omega: &Expr) -> Expr {

    let k_b = Expr::new_variable("k_B");

    simplify(&Expr::new_mul(
        k_b,
        Expr::new_log(omega.clone()),
    ))
}

/// Calculates the efficiency of a Carnot engine: $\eta = 1 - \`frac{T_c}{T_h`}$.
#[must_use]

pub fn carnot_efficiency(
    t_cold: &Expr,
    t_hot: &Expr,
) -> Expr {

    simplify(&Expr::new_sub(
        Expr::Constant(1.0),
        Expr::new_div(
            t_cold.clone(),
            t_hot.clone(),
        ),
    ))
}

/// Represents the Boltzmann Distribution: $`P_i` = \frac{e^{-E_i / (`k_B` T)}}{Z}$.
#[must_use]

pub fn boltzmann_distribution(
    energy: &Expr,
    temperature: &Expr,
    partition_function: &Expr,
) -> Expr {

    let k_b = Expr::new_variable("k_B");

    let exponent = Expr::new_neg(Arc::new(
        Expr::new_div(
            energy.clone(),
            Expr::new_mul(
                k_b,
                temperature.clone(),
            ),
        ),
    ));

    simplify(&Expr::new_div(
        Expr::new_exp(exponent),
        partition_function.clone(),
    ))
}

/// Calculates the Partition Function: $Z = \`sum_i` e^{-E_i / (`k_B` T)}$.
#[must_use]

pub fn partition_function(
    energies: &[Expr],
    temperature: &Expr,
) -> Expr {

    let k_b = Expr::new_variable("k_B");

    let mut z = Expr::Constant(0.0);

    for e in energies {

        let exponent = Expr::new_neg(Arc::new(
            Expr::new_div(
                e.clone(),
                Expr::new_mul(
                    k_b.clone(),
                    temperature.clone(),
                ),
            ),
        ));

        z = Expr::new_add(
            z,
            Expr::new_exp(exponent),
        );
    }

    simplify(&z)
}

/// Represents the Fermi-Dirac Distribution for fermions.
///
/// Formula: $f(E) = \frac{1}{e^{(E-\mu)/(k_B T)} + 1}$
#[must_use]

pub fn fermi_dirac_distribution(
    energy: &Expr,
    chemical_potential: &Expr,
    temperature: &Expr,
) -> Expr {

    let k_b = Expr::new_variable("k_B");

    let exponent = Expr::new_div(
        Expr::new_sub(
            energy.clone(),
            chemical_potential.clone(),
        ),
        Expr::new_mul(
            k_b,
            temperature.clone(),
        ),
    );

    simplify(&Expr::new_div(
        Expr::Constant(1.0),
        Expr::new_add(
            Expr::new_exp(exponent),
            Expr::Constant(1.0),
        ),
    ))
}

/// Represents the Bose-Einstein Distribution for bosons.
///
/// Formula: $f(E) = \frac{1}{e^{(E-\mu)/(k_B T)} - 1}$
#[must_use]

pub fn bose_einstein_distribution(
    energy: &Expr,
    chemical_potential: &Expr,
    temperature: &Expr,
) -> Expr {

    let k_b = Expr::new_variable("k_B");

    let exponent = Expr::new_div(
        Expr::new_sub(
            energy.clone(),
            chemical_potential.clone(),
        ),
        Expr::new_mul(
            k_b,
            temperature.clone(),
        ),
    );

    simplify(&Expr::new_div(
        Expr::Constant(1.0),
        Expr::new_sub(
            Expr::new_exp(exponent),
            Expr::Constant(1.0),
        ),
    ))
}

/// Calculates the work done during an isothermal expansion: $W = nRT \`ln(V_2/V_1)`$.
#[must_use]

pub fn work_isothermal_expansion(
    n: &Expr,
    r: &Expr,
    t: &Expr,
    v1: &Expr,
    v2: &Expr,
) -> Expr {

    simplify(&Expr::new_mul(
        Expr::new_mul(
            n.clone(),
            Expr::new_mul(r.clone(), t.clone()),
        ),
        Expr::new_log(Expr::new_div(
            v2.clone(),
            v1.clone(),
        )),
    ))
}

/// Maxwell Relation examples: $(\partial S / \partial V)_T = (\partial P / \partial T)_V$.
///
/// This function takes a thermodynamic potential (e.g., Helmholtz $A(T, V)$)
/// and verifies the Maxwell relation by taking mixed partial derivatives.
#[must_use]

pub fn verify_maxwell_relation_helmholtz(
    a: &Expr,
    t_var: &str,
    v_var: &str,
) -> Expr {

    // d2A / (dT dV) should equal d2A / (dV dT)
    let da_dt = differentiate(a, t_var);

    let d2a_dt_dv = differentiate(&da_dt, v_var);

    let da_dv = differentiate(a, v_var);

    let d2a_dv_dt = differentiate(&da_dv, t_var);

    simplify(&Expr::new_sub(
        d2a_dt_dv, d2a_dv_dt,
    ))
}
