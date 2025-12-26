//! # Solid-State Physics
//!
//! This module provides symbolic tools for solid-state physics, including representations
//! of crystal lattices, Bloch's theorem for electron wave functions in periodic potentials,
//! energy band models, and various physical properties like the density of states and
//! Fermi energy.

use crate::symbolic::core::Expr;
use crate::symbolic::simplify_dag::simplify;
use crate::symbolic::vector::Vector;
use serde::{Deserialize, Serialize};

/// Represents a crystal lattice with basis vectors.
#[derive(Clone, Debug, Serialize, Deserialize)]

pub struct CrystalLattice {
    pub a1: Vector,
    pub a2: Vector,
    pub a3: Vector,
}

impl CrystalLattice {
    /// Creates a new crystal lattice.
    #[must_use]

    pub const fn new(
        a1: Vector,
        a2: Vector,
        a3: Vector,
    ) -> Self {

        Self { a1, a2, a3 }
    }

    /// Computes the volume of the unit cell: V = |a1 . (a2 x a3)|.
    ///
    /// # Examples
    ///
    /// ```
    /// 
    /// use rssn::symbolic::core::Expr;
    /// use rssn::symbolic::solid_state_physics::CrystalLattice;
    /// use rssn::symbolic::vector::Vector;
    ///
    /// let a1 = Vector::new(
    ///     Expr::Constant(1.0),
    ///     Expr::Constant(0.0),
    ///     Expr::Constant(0.0),
    /// );
    ///
    /// let a2 = Vector::new(
    ///     Expr::Constant(0.0),
    ///     Expr::Constant(1.0),
    ///     Expr::Constant(0.0),
    /// );
    ///
    /// let a3 = Vector::new(
    ///     Expr::Constant(0.0),
    ///     Expr::Constant(0.0),
    ///     Expr::Constant(1.0),
    /// );
    ///
    /// let lattice = CrystalLattice::new(a1, a2, a3);
    ///
    /// assert_eq!(
    ///     lattice.volume(),
    ///     Expr::Constant(1.0)
    /// );
    /// ```
    #[must_use]

    pub fn volume(&self) -> Expr {

        let a2_cross_a3 = self
            .a2
            .cross(&self.a3);

        simplify(
            &self
                .a1
                .dot(&a2_cross_a3),
        )
    }

    /// Computes the reciprocal lattice vectors:
    /// b1 = 2π * (a2 x a3) / V
    /// b2 = 2π * (a3 x a1) / V
    /// b3 = 2π * (a1 x a2) / V
    ///
    /// # Examples
    ///
    /// ```
    /// 
    /// use rssn::symbolic::core::Expr;
    /// use rssn::symbolic::solid_state_physics::CrystalLattice;
    /// use rssn::symbolic::vector::Vector;
    ///
    /// let a1 = Vector::new(
    ///     Expr::Constant(1.0),
    ///     Expr::Constant(0.0),
    ///     Expr::Constant(0.0),
    /// );
    ///
    /// let a2 = Vector::new(
    ///     Expr::Constant(0.0),
    ///     Expr::Constant(1.0),
    ///     Expr::Constant(0.0),
    /// );
    ///
    /// let a3 = Vector::new(
    ///     Expr::Constant(0.0),
    ///     Expr::Constant(0.0),
    ///     Expr::Constant(1.0),
    /// );
    ///
    /// let lattice = CrystalLattice::new(a1, a2, a3);
    ///
    /// let (b1, b2, b3) = lattice.reciprocal_lattice_vectors();
    /// ```
    #[must_use]

    pub fn reciprocal_lattice_vectors(
        &self
    ) -> (
        Vector,
        Vector,
        Vector,
    ) {

        let v = self.volume();

        let two_pi = Expr::new_mul(
            Expr::Constant(2.0),
            Expr::new_variable("pi"),
        );

        let b1 = self
            .a2
            .cross(&self.a3)
            .scalar_mul(&two_pi)
            .scalar_mul(&Expr::new_div(
                Expr::Constant(1.0),
                v.clone(),
            ));

        let b2 = self
            .a3
            .cross(&self.a1)
            .scalar_mul(&two_pi)
            .scalar_mul(&Expr::new_div(
                Expr::Constant(1.0),
                v.clone(),
            ));

        let b3 = self
            .a1
            .cross(&self.a2)
            .scalar_mul(&two_pi)
            .scalar_mul(&Expr::new_div(
                Expr::Constant(1.0),
                v,
            ));

        (
            Vector::new(
                simplify(&b1.x),
                simplify(&b1.y),
                simplify(&b1.z),
            ),
            Vector::new(
                simplify(&b2.x),
                simplify(&b2.y),
                simplify(&b2.z),
            ),
            Vector::new(
                simplify(&b3.x),
                simplify(&b3.y),
                simplify(&b3.z),
            ),
        )
    }
}

/// Represents Bloch's Theorem: `ψ_k(r) = exp(i*k*r) * u_k(r)`.
///
/// Bloch's theorem states that the wave function of an electron in a periodic potential
/// can be expressed as a product of a plane wave `exp(i*k*r)` and a periodic function `u_k(r)`.
#[must_use]

pub fn bloch_theorem(
    k_vector: &Vector,
    r_vector: &Vector,
    periodic_function: &Expr,
) -> Expr {

    let i = Expr::new_complex(
        Expr::Constant(0.0),
        Expr::Constant(1.0),
    );

    let k_dot_r = k_vector.dot(r_vector);

    let ikr = Expr::new_mul(i, k_dot_r);

    let exp_term = Expr::new_exp(ikr);

    simplify(&Expr::new_mul(
        exp_term,
        periodic_function.clone(),
    ))
}

/// Represents a simple energy band model using the parabolic band approximation.
///
/// `E(k) = E_c + (hbar^2 * k^2) / (2 * m*)`
#[must_use]

pub fn energy_band(
    k_magnitude: &Expr,
    effective_mass: &Expr,
    band_edge: &Expr,
) -> Expr {

    let hbar = Expr::new_variable("hbar");

    let hbar_sq = Expr::new_pow(
        hbar,
        Expr::Constant(2.0),
    );

    let k_sq = Expr::new_pow(
        k_magnitude.clone(),
        Expr::Constant(2.0),
    );

    let kinetic_term = Expr::new_div(
        Expr::new_mul(hbar_sq, k_sq),
        Expr::new_mul(
            Expr::Constant(2.0),
            effective_mass.clone(),
        ),
    );

    simplify(&Expr::new_add(
        band_edge.clone(),
        kinetic_term,
    ))
}

/// Computes the Density of States (DOS) for a 3D electron gas.
///
/// `D(E) = (V / 2π^2) * (2m* / hbar^2)^(3/2) * sqrt(E)`
#[must_use]

pub fn density_of_states_3d(
    energy: &Expr,
    effective_mass: &Expr,
    volume: &Expr,
) -> Expr {

    let hbar = Expr::new_variable("hbar");

    let pi = Expr::new_variable("pi");

    let factor1 = Expr::new_div(
        volume.clone(),
        Expr::new_mul(
            Expr::Constant(2.0),
            Expr::new_pow(
                pi,
                Expr::Constant(2.0),
            ),
        ),
    );

    let factor2 = Expr::new_pow(
        Expr::new_div(
            Expr::new_mul(
                Expr::Constant(2.0),
                effective_mass.clone(),
            ),
            Expr::new_pow(
                hbar,
                Expr::Constant(2.0),
            ),
        ),
        Expr::new_div(
            Expr::Constant(3.0),
            Expr::Constant(2.0),
        ),
    );

    let sqrt_e = Expr::new_pow(
        energy.clone(),
        Expr::new_div(
            Expr::Constant(1.0),
            Expr::Constant(2.0),
        ),
    );

    simplify(&Expr::new_mul(
        factor1,
        Expr::new_mul(factor2, sqrt_e),
    ))
}

/// Fermi Energy for a 3D electron gas: `E_F = (hbar^2 / 2m*) * (3π^2 * n)^(2/3)`
/// where `n = N / V` is the electron concentration.
#[must_use]

pub fn fermi_energy_3d(
    electron_concentration: &Expr,
    effective_mass: &Expr,
) -> Expr {

    let hbar = Expr::new_variable("hbar");

    let pi = Expr::new_variable("pi");

    let term1 = Expr::new_div(
        Expr::new_pow(
            hbar,
            Expr::Constant(2.0),
        ),
        Expr::new_mul(
            Expr::Constant(2.0),
            effective_mass.clone(),
        ),
    );

    let term2 = Expr::new_pow(
        Expr::new_mul(
            Expr::new_mul(
                Expr::Constant(3.0),
                Expr::new_pow(
                    pi,
                    Expr::Constant(2.0),
                ),
            ),
            electron_concentration.clone(),
        ),
        Expr::new_div(
            Expr::Constant(2.0),
            Expr::Constant(3.0),
        ),
    );

    simplify(&Expr::new_mul(
        term1, term2,
    ))
}

/// Drude model electrical conductivity: `σ = (n * e^2 * τ) / m*`
#[must_use]

pub fn drude_conductivity(
    n: &Expr,
    e_charge: &Expr,
    relaxation_time: &Expr,
    effective_mass: &Expr,
) -> Expr {

    simplify(&Expr::new_div(
        Expr::new_mul(
            n.clone(),
            Expr::new_mul(
                Expr::new_pow(
                    e_charge.clone(),
                    Expr::Constant(2.0),
                ),
                relaxation_time.clone(),
            ),
        ),
        effective_mass.clone(),
    ))
}

/// Hall Coefficient: `R_H = 1 / (n * q)`
#[must_use]

pub fn hall_coefficient(
    carrier_concentration: &Expr,
    carrier_charge: &Expr,
) -> Expr {

    simplify(&Expr::new_div(
        Expr::Constant(1.0),
        Expr::new_mul(
            carrier_concentration.clone(),
            carrier_charge.clone(),
        ),
    ))
}

/// Debye Frequency: `ω_D = v_s * (6π^2 * n)^(1/3)`
/// where `v_s` is the speed of sound and `n` is the number density of atoms.
#[must_use]

pub fn debye_frequency(
    sound_velocity: &Expr,
    atom_density: &Expr,
) -> Expr {

    let pi = Expr::new_variable("pi");

    let inner = Expr::new_mul(
        Expr::new_mul(
            Expr::Constant(6.0),
            Expr::new_pow(
                pi,
                Expr::Constant(2.0),
            ),
        ),
        atom_density.clone(),
    );

    simplify(&Expr::new_mul(
        sound_velocity.clone(),
        Expr::new_pow(
            inner,
            Expr::new_div(
                Expr::Constant(1.0),
                Expr::Constant(3.0),
            ),
        ),
    ))
}

/// Einstein Heat Capacity: `C_v = 3Nk_B * (Θ_E / T)^2 * exp(Θ_E / T) / (exp(Θ_E / T) - 1)^2`
#[must_use]

pub fn einstein_heat_capacity(
    n_atoms: &Expr,
    einstein_temp: &Expr,
    temperature: &Expr,
) -> Expr {

    let k_b = Expr::new_variable("k_B");

    let x = Expr::new_div(
        einstein_temp.clone(),
        temperature.clone(),
    );

    let exp_x = Expr::new_exp(x.clone());

    let numerator = Expr::new_mul(
        Expr::new_mul(
            Expr::Constant(3.0),
            Expr::new_mul(n_atoms.clone(), k_b),
        ),
        Expr::new_mul(
            Expr::new_pow(
                x.clone(),
                Expr::Constant(2.0),
            ),
            exp_x,
        ),
    );

    let denominator = Expr::new_pow(
        Expr::new_sub(
            Expr::new_exp(x),
            Expr::Constant(1.0),
        ),
        Expr::Constant(2.0),
    );

    simplify(&Expr::new_div(
        numerator,
        denominator,
    ))
}

/// Plasma Frequency: `ω_p = sqrt((n * e^2) / (ε_0 * m*))`
#[must_use]

pub fn plasma_frequency(
    n: &Expr,
    e_charge: &Expr,
    epsilon_0: &Expr,
    effective_mass: &Expr,
) -> Expr {

    let numerator = Expr::new_mul(
        n.clone(),
        Expr::new_pow(
            e_charge.clone(),
            Expr::Constant(2.0),
        ),
    );

    let denominator = Expr::new_mul(
        epsilon_0.clone(),
        effective_mass.clone(),
    );

    simplify(&Expr::new_pow(
        Expr::new_div(
            numerator,
            denominator,
        ),
        Expr::new_div(
            Expr::Constant(1.0),
            Expr::Constant(2.0),
        ),
    ))
}

/// London penetration depth: `λ_L = sqrt(m / (μ_0 * n_s * q^2))`
#[must_use]

pub fn london_penetration_depth(
    mass: &Expr,
    mu_0: &Expr,
    supercarrier_density: &Expr,
    charge: &Expr,
) -> Expr {

    let numerator = mass.clone();

    let denominator = Expr::new_mul(
        mu_0.clone(),
        Expr::new_mul(
            supercarrier_density.clone(),
            Expr::new_pow(
                charge.clone(),
                Expr::Constant(2.0),
            ),
        ),
    );

    simplify(&Expr::new_pow(
        Expr::new_div(
            numerator,
            denominator,
        ),
        Expr::new_div(
            Expr::Constant(1.0),
            Expr::Constant(2.0),
        ),
    ))
}
