//! # Symbolic Thermodynamics
//!
//! This module provides symbolic representations of fundamental thermodynamic laws
//! and distributions. It includes the First Law of Thermodynamics, definitions of
//! Helmholtz and Gibbs free energies, and statistical distributions like Boltzmann,
//! Fermi-Dirac, and Bose-Einstein.

use crate::symbolic::core::Expr;

/// Represents the First Law of Thermodynamics: `dU = dQ - dW`.
///
/// The First Law of Thermodynamics states that energy cannot be created or destroyed,
/// only transferred or changed from one form to another. It relates the change in
/// internal energy (`dU`) to the heat added to the system (`dQ`) and the work done
/// by the system (`dW`).
///
/// # Arguments
/// * `internal_energy_change` - The change in internal energy `dU`.
/// * `heat_added` - The heat added to the system `dQ`.
/// * `work_done` - The work done by the system `dW`.
///
/// # Returns
/// An `Expr` representing the symbolic First Law of Thermodynamics.
pub fn first_law_thermodynamics(
    internal_energy_change: Expr,
    heat_added: Expr,
    work_done: Expr,
) -> Expr {
    Expr::Sub(
        Box::new(internal_energy_change),
        Box::new(Expr::Sub(Box::new(heat_added), Box::new(work_done))),
    )
}

/// Represents the Helmholtz Free Energy: `A = U - T*S`.
///
/// Helmholtz free energy is a thermodynamic potential that measures the "useful" or
/// process-initiating work obtainable from a closed thermodynamic system at a constant
/// temperature and volume. `U` is internal energy, `T` is temperature, and `S` is entropy.
///
/// # Arguments
/// * `internal_energy` - The internal energy `U`.
/// * `temperature` - The temperature `T`.
/// * `entropy` - The entropy `S`.
///
/// # Returns
/// An `Expr` representing the symbolic Helmholtz Free Energy.
pub fn helmholtz_free_energy(internal_energy: Expr, temperature: Expr, entropy: Expr) -> Expr {
    Expr::Sub(
        Box::new(internal_energy),
        Box::new(Expr::Mul(Box::new(temperature), Box::new(entropy))),
    )
}

/// Represents the Gibbs Free Energy: `G = H - T*S = U + P*V - T*S`.
///
/// Gibbs free energy is a thermodynamic potential that measures the "useful" or
/// process-initiating work obtainable from an isothermal, isobaric thermodynamic system.
/// `H` is enthalpy, `U` is internal energy, `P` is pressure, `V` is volume, `T` is temperature,
/// and `S` is entropy.
///
/// # Arguments
/// * `internal_energy` - The internal energy `U`.
/// * `pressure` - The pressure `P`.
/// * `volume` - The volume `V`.
/// * `temperature` - The temperature `T`.
/// * `entropy` - The entropy `S`.
///
/// # Returns
/// An `Expr` representing the symbolic Gibbs Free Energy.
pub fn gibbs_free_energy(
    internal_energy: Expr,
    pressure: Expr,
    volume: Expr,
    temperature: Expr,
    entropy: Expr,
) -> Expr {
    let pv_term = Expr::Mul(Box::new(pressure), Box::new(volume));
    let ts_term = Expr::Mul(Box::new(temperature), Box::new(entropy));
    Expr::Sub(
        Box::new(Expr::Add(Box::new(internal_energy), Box::new(pv_term))),
        Box::new(ts_term),
    )
}

/// Represents the Boltzmann Distribution: `P_i = exp(-E_i / (k*T)) / Z`.
///
/// The Boltzmann distribution describes the probability of a system being in a certain
/// state `i` with energy `E_i` at a given temperature `T`. `k` is the Boltzmann constant,
/// and `Z` is the partition function.
///
/// # Arguments
/// * `energy` - The energy `E_i` of the state.
/// * `temperature` - The temperature `T`.
/// * `partition_function` - The partition function `Z`.
///
/// # Returns
/// An `Expr` representing the symbolic Boltzmann Distribution.
pub fn boltzmann_distribution(energy: Expr, temperature: Expr, partition_function: Expr) -> Expr {
    let k = Expr::Variable("k".to_string());
    let kt_term = Expr::Mul(Box::new(k), Box::new(temperature));
    let exponent = Expr::Div(
        Box::new(Expr::Mul(Box::new(Expr::Constant(-1.0)), Box::new(energy))),
        Box::new(kt_term),
    );
    let numerator = Expr::Exp(Box::new(exponent));
    Expr::Div(Box::new(numerator), Box::new(partition_function))
}

/// Represents the Partition Function: `Z = sum(exp(-E_i / (k*T)))`.
///
/// The partition function `Z` is a fundamental quantity in statistical mechanics
/// that encodes the statistical properties of a system in thermodynamic equilibrium.
/// It is a sum over all possible states `i` of the system.
///
/// # Arguments
/// * `energies` - A vector of `Expr` representing the energies `E_i` of the states.
/// * `temperature` - The temperature `T`.
///
/// # Returns
/// An `Expr` representing the symbolic Partition Function.
pub fn partition_function(energies: Vec<Expr>, temperature: Expr) -> Expr {
    let k = Expr::Variable("k".to_string());
    let kt_term = Expr::Mul(Box::new(k), Box::new(temperature.clone()));

    let terms = energies.into_iter().map(|energy| {
        let exponent = Expr::Div(
            Box::new(Expr::Mul(Box::new(Expr::Constant(-1.0)), Box::new(energy))),
            Box::new(kt_term.clone()),
        );
        Expr::Exp(Box::new(exponent))
    });

    terms
        .reduce(|acc, term| Expr::Add(Box::new(acc), Box::new(term)))
        .unwrap_or(Expr::Constant(0.0))
}

/// Represents the Fermi-Dirac Distribution.
///
/// The Fermi-Dirac distribution describes the probability that a fermion will occupy
/// a given quantum state at a given temperature. It is crucial for understanding
/// the behavior of electrons in solids.
///
/// # Arguments
/// * `energy` - The energy `E` of the state.
/// * `fermi_level` - The Fermi level `μ`.
/// * `temperature` - The temperature `T`.
///
/// # Returns
/// An `Expr` representing the symbolic Fermi-Dirac Distribution.
pub fn fermi_dirac_distribution(energy: Expr, fermi_level: Expr, temperature: Expr) -> Expr {
    let k = Expr::Variable("k".to_string());
    let kt_term = Expr::Mul(Box::new(k), Box::new(temperature));
    let energy_diff = Expr::Sub(Box::new(energy), Box::new(fermi_level));
    let exponent = Expr::Div(Box::new(energy_diff), Box::new(kt_term));
    let exp_term = Expr::Exp(Box::new(exponent));
    let denominator = Expr::Add(Box::new(exp_term), Box::new(Expr::Constant(1.0)));
    Expr::Div(Box::new(Expr::Constant(1.0)), Box::new(denominator))
}

/// Represents the Bose-Einstein Distribution.
///
/// The Bose-Einstein distribution describes the probability that a boson will occupy
/// a given quantum state at a given temperature. It is used for systems of identical,
/// indistinguishable bosons.
///
/// # Arguments
/// * `energy` - The energy `E` of the state.
/// * `chemical_potential` - The chemical potential `μ`.
/// * `temperature` - The temperature `T`.
///
/// # Returns
/// An `Expr` representing the symbolic Bose-Einstein Distribution.
pub fn bose_einstein_distribution(
    energy: Expr,
    chemical_potential: Expr,
    temperature: Expr,
) -> Expr {
    let k = Expr::Variable("k".to_string());
    let kt_term = Expr::Mul(Box::new(k), Box::new(temperature));
    let energy_diff = Expr::Sub(Box::new(energy), Box::new(chemical_potential));
    let exponent = Expr::Div(Box::new(energy_diff), Box::new(kt_term));
    let exp_term = Expr::Exp(Box::new(exponent));
    let denominator = Expr::Sub(Box::new(exp_term), Box::new(Expr::Constant(1.0)));
    Expr::Div(Box::new(Expr::Constant(1.0)), Box::new(denominator))
}
