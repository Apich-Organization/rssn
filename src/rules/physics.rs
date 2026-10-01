//! Physics: the formula libraries of classical mechanics,
//! electromagnetism, relativity, quantum mechanics, quantum field theory,
//! solid-state physics and thermodynamics.
//!
//! Every formula is an ordinary operator, defined by its right-hand side
//! (see [`Installer::define`](crate::graph::Installer::define)) or, where
//! it needs a loop or a case analysis, by a kernel. Arguments may be
//! numbers, symbols or whole expressions, and vectors are lists, so the
//! formulas compose with each other and with the rest of the engine:
//! `kinetic_energy(m, norm(diff(r, t)))` differentiates, takes the norm
//! and squares in one request.
//!
//! # Constants
//!
//! Physical constants are nullary operators that stay symbolic in exact
//! answers and evaluate to their CODATA values in numeric ones:
//! `c_0`, `hbar`, `h_planck`, `k_B`, `epsilon_0`, `mu_0`, `G_N`, `q_e`
//! (elementary charge), `m_e`, `m_p`, `N_A`, `R_gas`, `sigma_SB`, `a_0`.
//!
//! # Quantum operators
//!
//! An operator acting on wave functions is a term built from
//! `op_mul(f)` (multiplication by `f`), `op_d(x)` (`∂/∂x`),
//! `op_laplacian(list(x, ...))`, `op_identity`, `op_add(A, B, ...)`,
//! `op_scale(c, A)`, `op_compose(A, B, ...)` (`B` acts first) and
//! `op_power(A, n)`; a matrix acts on a spinor by `matmul`.
//! `qm_apply(A, psi)` applies one. Ready-made operators:
//! `position_operator(x)`, `momentum_operator(x)`, `kinetic_operator(m, x)`,
//! `hamiltonian_operator(m, V, x)`, `hamiltonian_3d(m, V, vars)`,
//! `hamiltonian_free_particle(m, x)`, `hamiltonian_harmonic_oscillator(m, w, x)`,
//! `angular_momentum_z(phi)`, the Pauli matrices `pauli_x`, `pauli_y`,
//! `pauli_z` and the Dirac matrices `dirac_gamma(mu)` (Dirac
//! representation, `mu = 0..3`, `5` for `γ⁵`).
//!
//! Inner products integrate over the whole line: `braket(phi, psi, x)` is
//! `∫ conj(phi) psi dx`; `braket_on(phi, psi, x, a, b)` over `[a, b]`.
//!
//! # Conventions
//!
//! Equations are returned as residuals (`lhs - rhs`, zero when the
//! equation holds). Four-vectors are `list(p0, p1, p2, p3)` with
//! contravariant components and the metric `diag(1, -1, -1, -1)`; field
//! theory uses natural units. The Schwarzschild metric uses the signature
//! `(-, +, +, +)` and coordinates `list(t, r, theta, phi)`.

use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Facts;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::OnReals;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::powi;
use crate::rules::complex::build::sub;
use crate::rules::poly::best;
use crate::rules::poly::repr::from_term;
use crate::rules::poly::repr::to_term;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;

use super::calculus::derivative;
use super::complex::complex;
use super::geometry::geometry;
use super::special::special;
use super::variational::variational;

/// The physics rule set.
#[must_use]
pub fn physics() -> RuleSet {
    RuleSet::new("physics", install)
        .needs(variational())
        .needs(geometry())
        .needs(complex())
        .needs(special())
}

/// `(name, value)` of the physical constants, SI units.
const CONSTANTS: [(&str, f64); 14] = [
    ("c_0", crate::constant::SPEED_OF_LIGHT),
    ("hbar", crate::constant::REDUCED_PLANCK_CONSTANT),
    ("h_planck", crate::constant::PLANCK_CONSTANT),
    ("k_B", crate::constant::BOLTZMANN_CONSTANT),
    ("epsilon_0", crate::constant::VACUUM_ELECTRIC_PERMITTIVITY),
    ("mu_0", crate::constant::VACUUM_MAGNETIC_PERMEABILITY),
    ("G_N", crate::constant::GRAVITATIONAL_CONSTANT),
    ("q_e", crate::constant::ELEMENTARY_CHARGE),
    ("m_e", crate::constant::ELECTRON_MASS),
    ("m_p", crate::constant::PROTON_MASS_KG),
    ("N_A", crate::constant::AVOGADRO_CONSTANT),
    ("R_gas", crate::constant::MOLAR_GAS_CONSTANT),
    ("sigma_SB", crate::constant::STEFAN_BOLTZMANN_CONSTANT),
    ("a_0", crate::constant::BOHR_RADIUS),
];

/// The value of a physical constant by name.
fn constant_value(name: &str) -> f64 {
    CONSTANTS.iter().find(|c| c.0 == name).map_or(f64::NAN, |c| c.1)
}

macro_rules! constant_eval {
    ($($name:literal),*) => {
        [$(($name, (|_: &[f64]| constant_value($name)) as crate::graph::op::EvalFn)),*]
    };
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Request {
    Apply,
    Eigenvalue,
    PoissonBracket,
    GeodesicAcceleration,
    Schwarzschild,
    DiracGamma,
    DiracEquation,
    KleinGordon,
    Slash,
    DiracAdjoint,
    MinkowskiDot,
    ScalarLagrangian,
    FermionPropagator,
    PartitionFunction,
}

/// The operator-algebra constructors.
#[derive(Copy, Clone, Debug)]
struct OpAlgebra {
    mul: OpId,
    d: OpId,
    laplacian: OpId,
    identity: OpId,
    add: OpId,
    scale: OpId,
    compose: OpId,
    power: OpId,
    matmul: OpId,
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, eval) in constant_eval!(
        "c_0", "hbar", "h_planck", "k_B", "epsilon_0", "mu_0", "G_N", "q_e", "m_e", "m_p", "N_A", "R_gas", "sigma_SB", "a_0"
    ) {
        let op = i.op(OpDescriptor::new(name, Arity::Fixed(0)).eval(eval))?;
        i.graph().ops_mut().set_attr(op, OnReals(Facts::POSITIVE));
    }
    let plain = |name: &str, arity: Arity| OpDescriptor::new(name, arity);
    let algebra = OpAlgebra {
        mul: i.op(plain("op_mul", Arity::Fixed(1)))?,
        d: i.op(plain("op_d", Arity::Fixed(1)))?,
        laplacian: i.op(plain("op_laplacian", Arity::Fixed(1)))?,
        identity: i.op(plain("op_identity", Arity::Fixed(0)))?,
        add: i.op(plain("op_add", Arity::Variadic))?,
        scale: i.op(plain("op_scale", Arity::Fixed(2)))?,
        compose: i.op(plain("op_compose", Arity::Variadic))?,
        power: i.op(plain("op_power", Arity::Fixed(2)))?,
        matmul: i.graph().ops().lookup("matmul").ok_or(RuleError::Invalid {
            rule: "physics".to_owned(),
            reason: "needs the linear algebra rule set",
        })?,
    };
    let heavy = |name: &str, arity: u8| OpDescriptor::new(name, Arity::Fixed(arity)).flags(OpFlags::HEAVY).cost(100);
    for (name, arity, request) in [
        ("qm_apply", 2, Request::Apply),
        ("energy_eigenvalue", 3, Request::Eigenvalue),
        ("poisson_bracket", 4, Request::PoissonBracket),
        ("geodesic_acceleration", 2, Request::GeodesicAcceleration),
        ("schwarzschild_metric", 2, Request::Schwarzschild),
        ("dirac_gamma", 1, Request::DiracGamma),
        ("dirac_equation", 3, Request::DiracEquation),
        ("klein_gordon_equation", 3, Request::KleinGordon),
        ("feynman_slash", 1, Request::Slash),
        ("dirac_adjoint", 1, Request::DiracAdjoint),
        ("minkowski_dot", 2, Request::MinkowskiDot),
        ("scalar_field_lagrangian", 3, Request::ScalarLagrangian),
        ("fermion_propagator", 2, Request::FermionPropagator),
        ("partition_function", 2, Request::PartitionFunction),
    ] {
        let op = i.op(heavy(name, arity))?;
        i.kernel(&format!("physics/{name}"), Tier::Reduce, Physics { op, request, algebra });
    }
    i.define(&CLASSICAL)?;
    i.define(&ELECTROMAGNETISM)?;
    i.define(&RELATIVITY)?;
    i.define(&QUANTUM)?;
    i.define(&FIELD_THEORY)?;
    i.define(&SOLID_STATE)?;
    i.define(&THERMODYNAMICS)
}

const CLASSICAL: [&str; 21] = [
    "kinematics(r, t) := list(r, diff(r, t), diff(diff(r, t), t))",
    "newtons_second_law(m, a) := m * a",
    "momentum(m, v) := m * v",
    "kinetic_energy(m, v) := m * v^2 / 2",
    "potential_energy_gravity_uniform(m, h, g) := m * g * h",
    "potential_energy_gravity_universal(m1, m2, r, G) := -G * m1 * m2 / r",
    "potential_energy_spring(k, x) := k * x^2 / 2",
    "work_constant_force(F, d) := dot(F, d)",
    "work_line_integral(F, curve, t, a, b) := line_integral_vec(F, curve, t, a, b)",
    "power(F, v) := dot(F, v)",
    "torque(r, F) := cross(r, F)",
    "angular_momentum(r, p) := cross(r, p)",
    "centripetal_acceleration(v, r) := v^2 / r",
    "moment_of_inertia_point_mass(m, r) := m * r^2",
    "rotational_kinetic_energy(inertia, w) := inertia * w^2 / 2",
    "lagrangian(T, V) := T - V",
    "hamiltonian(T, V) := T + V",
    "euler_lagrange_equation(L, q, t) := euler_lagrange(L, q, t)",
    "hamilton_equations(H, q, p) := list(diff(H, p), -diff(H, q))",
    "gravitational_force(m1, m2, r) := G_N * m1 * m2 / r^2",
    "harmonic_oscillator_frequency(k, m) := (k / m)^(1/2)",
];

const ELECTROMAGNETISM: [&str; 10] = [
    "maxwell_equations(E, B, rho, J, vars, t) := list(div(E, vars) - rho / epsilon_0, div(B, vars), curl(E, vars) + diff(B, t), curl(B, vars) - mu_0 * (J + epsilon_0 * diff(E, t)))",
    "lorentz_force(q, E, v, B) := q * (E + cross(v, B))",
    "electric_field_from_potentials(V, A, vars, t) := -(grad(V, vars) + diff(A, t))",
    "electric_field_from_potential(V, vars) := -grad(V, vars)",
    "magnetic_field_from_vector_potential(A, vars) := curl(A, vars)",
    "poynting_vector(E, B) := cross(E, B) / mu_0",
    "em_energy_density(E, B) := (epsilon_0 * dot(E, E) + dot(B, B) / mu_0) / 2",
    "coulombs_law(q, r) := q * r / (4 * pi * epsilon_0 * norm(r)^3)",
    "coulomb_force(q1, q2, r) := q1 * q2 / (4 * pi * epsilon_0 * r^2)",
    "point_charge_potential(q, r) := q / (4 * pi * epsilon_0 * r)",
];

const RELATIVITY: [&str; 12] = [
    "lorentz_factor(v) := (1 - v^2 / c_0^2)^(-1/2)",
    "lorentz_transformation_x(x, t, v) := list(lorentz_factor(v) * (x - v * t), lorentz_factor(v) * (t - v * x / c_0^2))",
    "velocity_addition(v, u) := (v + u) / (1 + v * u / c_0^2)",
    "mass_energy_equivalence(m) := m * c_0^2",
    "relativistic_momentum(m, v) := lorentz_factor(v) * m * v",
    "relativistic_energy(m, p) := ((p * c_0)^2 + (m * c_0^2)^2)^(1/2)",
    "doppler_effect(f, v) := f * ((1 - v / c_0) / (1 + v / c_0))^(1/2)",
    "schwarzschild_radius(M) := 2 * G_N * M / c_0^2",
    "gravitational_time_dilation(t, r, M) := t * (1 - schwarzschild_radius(M) / r)^(1/2)",
    "einstein_tensor_from(Ric, R, g) := Ric - R * g / 2",
    "einstein_field_equations(Ric, R, g, T) := Ric - R * g / 2 - 8 * pi * G_N * T / c_0^4",
    "minkowski_metric := list(list(1, 0, 0, 0), list(0, -1, 0, 0), list(0, 0, -1, 0), list(0, 0, 0, -1))",
];

const QUANTUM: [&str; 27] = [
    "position_operator(x) := op_mul(x)",
    "momentum_operator(x) := op_scale(-I * hbar, op_d(x))",
    "kinetic_operator(m, x) := op_scale(-hbar^2 / (2 * m), op_compose(op_d(x), op_d(x)))",
    "hamiltonian_operator(m, V, x) := op_add(kinetic_operator(m, x), op_mul(V))",
    "hamiltonian_3d(m, V, vars) := op_add(op_scale(-hbar^2 / (2 * m), op_laplacian(vars)), op_mul(V))",
    "hamiltonian_free_particle(m, x) := kinetic_operator(m, x)",
    "hamiltonian_harmonic_oscillator(m, w, x) := hamiltonian_operator(m, m * w^2 * x^2 / 2, x)",
    "angular_momentum_z(phi) := op_scale(-I * hbar, op_d(phi))",
    "commutator(A, B, psi) := qm_apply(A, qm_apply(B, psi)) - qm_apply(B, qm_apply(A, psi))",
    "matrix_commutator(A, B) := matmul(A, B) - matmul(B, A)",
    "matrix_anticommutator(A, B) := matmul(A, B) + matmul(B, A)",
    "braket(phi, psi, x) := defint(conj(phi) * psi, x, -oo, oo)",
    "braket_on(phi, psi, x, a, b) := defint(conj(phi) * psi, x, a, b)",
    "expectation_value(A, psi, x) := braket(psi, qm_apply(A, psi), x) / braket(psi, psi, x)",
    "uncertainty(A, psi, x) := (expectation_value(op_compose(A, A), psi, x) - expectation_value(A, psi, x)^2)^(1/2)",
    "probability_density(psi) := conj(psi) * psi",
    "normalize_wavefunction(psi, x) := psi / braket(psi, psi, x)^(1/2)",
    "time_dependent_schrodinger_equation(H, psi, t) := I * hbar * diff(psi, t) - qm_apply(H, psi)",
    "time_independent_schrodinger_equation(H, psi, E) := qm_apply(H, psi) - E * psi",
    "solve_time_independent_schrodinger(H, psi, x) := list(energy_eigenvalue(H, psi, x), psi)",
    "first_order_energy_correction(V, psi, x) := expectation_value(V, psi, x)",
    "scattering_amplitude(psi_i, psi_f, V, x) := braket(psi_f, qm_apply(V, psi_i), x)",
    "pauli_x := list(list(0, 1), list(1, 0))",
    "pauli_y := list(list(0, -I), list(I, 0))",
    "pauli_z := list(list(1, 0), list(0, -1))",
    "pauli_matrices := list(pauli_x, pauli_y, pauli_z)",
    "spin_operator(sigma) := hbar * sigma / 2",
];

const FIELD_THEORY: [&str; 5] = [
    "qed_lagrangian(psibar, psi, A, m, e) := psibar * (I * gamma_mu * partial_mu - e * gamma_mu * A - m) * psi - F_mu_nu^2 / 4",
    "qcd_lagrangian(psibar, psi, G, m, gs) := psibar * (I * gamma_mu * partial_mu + gs * gamma_mu * G - m) * psi - G_mu_nu_a^2 / 4",
    "propagator(p, m) := I / (p^2 - m^2 + I * epsilon)",
    "scattering_cross_section(M, flux, dphi) := abs(M)^2 / flux * dphi",
    "feynman_propagator_position_space(x, y, m) := defint(propagator(p, m) * exp(-I * p * (x - y)), p, -oo, oo)",
];

const SOLID_STATE: [&str; 14] = [
    "lattice_volume(a1, a2, a3) := dot(a1, cross(a2, a3))",
    "reciprocal_lattice_vectors(a1, a2, a3) := list(2 * pi * cross(a2, a3) / lattice_volume(a1, a2, a3), 2 * pi * cross(a3, a1) / lattice_volume(a1, a2, a3), 2 * pi * cross(a1, a2) / lattice_volume(a1, a2, a3))",
    "bloch_wave(k, r, u) := exp(I * dot(k, r)) * u",
    "energy_band(k, m, E0) := E0 + hbar^2 * k^2 / (2 * m)",
    "density_of_states_3d(E, m, V) := V / (2 * pi^2) * (2 * m / hbar^2)^(3/2) * E^(1/2)",
    "fermi_energy_3d(n, m) := hbar^2 / (2 * m) * (3 * pi^2 * n)^(2/3)",
    "drude_conductivity(n, e, tau, m) := n * e^2 * tau / m",
    "hall_coefficient(n, q) := 1 / (n * q)",
    "debye_frequency(v, n) := v * (6 * pi^2 * n)^(1/3)",
    "einstein_heat_capacity(N, TE, T) := 3 * N * k_B * (TE / T)^2 * exp(TE / T) / (exp(TE / T) - 1)^2",
    "debye_heat_capacity(N, TD, T) := 9 * N * k_B * (T / TD)^3 * defint(u^4 * exp(u) / (exp(u) - 1)^2, u, 0, TD / T)",
    "plasma_frequency(n, e, eps, m) := (n * e^2 / (eps * m))^(1/2)",
    "london_penetration_depth(m, mu, ns, q) := (m / (mu * ns * q^2))^(1/2)",
    "fermi_wavevector_3d(n) := (3 * pi^2 * n)^(1/3)",
];

const THERMODYNAMICS: [&str; 21] = [
    "first_law_thermodynamics(dU, Q, W) := dU - (Q - W)",
    "ideal_gas_law(P, V, n, R, T) := P * V - n * R * T",
    "enthalpy(U, P, V) := U + P * V",
    "helmholtz_free_energy(U, T, S) := U - T * S",
    "gibbs_free_energy(H, T, S) := H - T * S",
    "boltzmann_entropy(Omega) := k_B * ln(Omega)",
    "carnot_efficiency(Tc, Th) := 1 - Tc / Th",
    "boltzmann_distribution(E, T, Z) := exp(-E / (k_B * T)) / Z",
    "fermi_dirac_distribution(E, mu, T) := 1 / (exp((E - mu) / (k_B * T)) + 1)",
    "bose_einstein_distribution(E, mu, T) := 1 / (exp((E - mu) / (k_B * T)) - 1)",
    "work_isothermal_expansion(n, R, T, V1, V2) := n * R * T * ln(V2 / V1)",
    "verify_maxwell_relation_helmholtz(A, T, V) := diff(diff(A, T), V) - diff(diff(A, V), T)",
    "helmholtz_from_partition(Z, T) := -k_B * T * ln(Z)",
    "internal_energy_from_partition(Z, T) := k_B * T^2 * diff(ln(Z), T)",
    "entropy_from_helmholtz(A, T) := -diff(A, T)",
    "pressure_from_helmholtz(A, V) := -diff(A, V)",
    "heat_capacity(U, T) := diff(U, T)",
    "planck_law(nu, T) := 2 * h_planck * nu^3 / c_0^2 / (exp(h_planck * nu / (k_B * T)) - 1)",
    "stefan_boltzmann_law(T) := sigma_SB * T^4",
    "wien_peak_wavelength(T) := 2.897771955e-3 / T",
    "ideal_gas_pressure(n, T, V) := n * R_gas * T / V",
];

struct Physics {
    op: OpId,
    request: Request,
    algebra: OpAlgebra,
}

impl Kernel for Physics {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        let result = match self.request {
            | Request::Apply => args.get(..2).and_then(|a| apply(cx, self.algebra, a[0], a[1], 0)),
            | Request::Eigenvalue => eigenvalue(cx, self.algebra, &args),
            | Request::PoissonBracket => poisson_bracket(cx, &args),
            | Request::GeodesicAcceleration => geodesic_acceleration(cx, &args),
            | Request::Schwarzschild => schwarzschild(cx, &args),
            | Request::DiracGamma => {
                let mu = cx.graph.number_of(*args.first().unwrap_or(&NodeId::NONE)).and_then(crate::graph::Number::to_i64);
                mu.and_then(|mu| gamma_matrix(cx.graph, mu)).map(|g| matrix_term(cx.graph, &g))
            },
            | Request::DiracEquation => dirac_equation(cx, &args),
            | Request::KleinGordon => klein_gordon(cx, &args),
            | Request::Slash => vector(cx.graph, args[0]).and_then(|p| slash(cx.graph, &p)).map(|m| matrix_term(cx.graph, &m)),
            | Request::DiracAdjoint => dirac_adjoint(cx, args[0]),
            | Request::MinkowskiDot => minkowski_dot(cx.graph, args[0], args[1]),
            | Request::ScalarLagrangian => scalar_lagrangian(cx, &args),
            | Request::FermionPropagator => fermion_propagator(cx, &args),
            | Request::PartitionFunction => partition_function(cx, &args),
        };
        result.map_or(Outcome::Pass, Outcome::Equal)
    }
}

// ----------------------------------------------------------------------
// Helpers
// ----------------------------------------------------------------------

fn constant(
    graph: &mut Graph,
    name: &str,
) -> Option<NodeId> {
    let op = graph.ops().lookup(name)?;
    Some(graph.node(op, &[]))
}

fn imaginary_unit(graph: &mut Graph) -> Option<NodeId> {
    constant(graph, "I")
}

/// The entries of a list, best forms.
fn vector(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Vec<NodeId>> {
    let term = best(graph, node)?;
    (graph.op(term) == core::LIST).then(|| graph.children(term).to_vec())
}

fn list(
    graph: &mut Graph,
    items: &[NodeId],
) -> NodeId {
    graph.node(core::LIST, items)
}

fn matrix_term(
    graph: &mut Graph,
    rows: &[Vec<NodeId>],
) -> NodeId {
    let rows: Vec<NodeId> = rows.iter().map(|r| list(graph, r)).collect();
    list(graph, &rows)
}

fn is_heavy(
    graph: &Graph,
    node: NodeId,
) -> bool {
    let mut stack = vec![node];
    let mut seen = std::collections::HashSet::new();
    while let Some(n) = stack.pop() {
        if !seen.insert(n) {
            continue;
        }
        if graph.ops().get(graph.op(n)).flags.has(OpFlags::HEAVY) {
            return true;
        }
        stack.extend_from_slice(graph.children(n));
    }
    false
}

// ----------------------------------------------------------------------
// Quantum operators
// ----------------------------------------------------------------------

/// `A psi` for an operator term `A`.
fn apply(
    cx: &mut Cx<'_>,
    algebra: OpAlgebra,
    a: NodeId,
    psi: NodeId,
    depth: usize,
) -> Option<NodeId> {
    if depth > 32 {
        return None;
    }
    let a = best(cx.graph, a)?;
    let psi = best(cx.graph, psi)?;
    let op = cx.graph.op(a);
    let args = cx.graph.children(a).to_vec();
    let graph = &mut *cx.graph;
    if op == algebra.mul {
        return Some(mul(graph, &[args[0], psi]));
    }
    if op == algebra.d {
        graph.symbol_of(args[0])?;
        return derivative(graph, psi, args[0]);
    }
    if op == algebra.identity {
        return Some(psi);
    }
    if op == algebra.laplacian {
        let vars = vector(graph, args[0])?;
        let mut terms = Vec::with_capacity(vars.len());
        for x in vars {
            graph.symbol_of(x)?;
            let first = derivative(graph, psi, x)?;
            terms.push(derivative(graph, first, x)?);
        }
        return Some(add(graph, &terms));
    }
    if op == algebra.add {
        let mut terms = Vec::with_capacity(args.len());
        for part in args {
            terms.push(apply(cx, algebra, part, psi, depth + 1)?);
        }
        return Some(add(cx.graph, &terms));
    }
    if op == algebra.scale {
        let inner = apply(cx, algebra, args[1], psi, depth + 1)?;
        return Some(mul(cx.graph, &[args[0], inner]));
    }
    if op == algebra.compose {
        let mut state = psi;
        for &part in args.iter().rev() {
            state = apply(cx, algebra, part, state, depth + 1)?;
        }
        return Some(state);
    }
    if op == algebra.power {
        let n = graph.number_of(args[1]).and_then(crate::graph::Number::to_i64).filter(|n| (0..=16).contains(n))?;
        let mut state = psi;
        for _ in 0..n {
            state = apply(cx, algebra, args[0], state, depth + 1)?;
        }
        return Some(state);
    }
    if op == core::LIST {
        return Some(graph.node(algebra.matmul, &[a, psi]));
    }
    // Any other term multiplies, once it is a closed form (an operator
    // definition that has not been expanded yet is not a multiplier).
    if is_heavy(graph, a) {
        return None;
    }
    Some(mul(graph, &[a, psi]))
}

/// `E` with `H psi = E psi`, when `H psi / psi` is free of the position
/// variables.
fn eigenvalue(
    cx: &mut Cx<'_>,
    algebra: OpAlgebra,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[h, psi, vars] = args else {
        return None;
    };
    let h_psi = apply(cx, algebra, h, psi, 0)?;
    let h_psi = cx.simplify(h_psi);
    let psi = best(cx.graph, psi)?;
    // Divide term by term after expanding, so that the wave function
    // cancels from every term on its own.
    let mut gens = Gens::default();
    let terms = match from_term(cx.graph, &mut gens, h_psi, Limits { terms: 256, exponent: 8 }) {
        | Some(p) => {
            let expanded = to_term(cx.graph, &gens, &p);
            if cx.graph.op(expanded) == core::ADD { cx.graph.children(expanded).to_vec() } else { vec![expanded] }
        },
        | None => vec![h_psi],
    };
    let inverse = powi(cx.graph, psi, -1);
    let mut quotients = Vec::with_capacity(terms.len());
    for term in terms {
        let q = mul(cx.graph, &[term, inverse]);
        quotients.push(cx.simplify(q));
    }
    let ratio = add(cx.graph, &quotients);
    let ratio = cx.simplify(ratio);
    let vars = vector(cx.graph, vars).unwrap_or_else(|| vec![vars]);
    for x in vars {
        let symbol = cx.graph.symbol_of(x)?;
        if cx.graph.depends_on(cx.graph.find(ratio), symbol) {
            return None;
        }
    }
    Some(ratio)
}

// ----------------------------------------------------------------------
// Classical mechanics and relativity
// ----------------------------------------------------------------------

/// `{f, g} = Σ ∂f/∂q ∂g/∂p - ∂f/∂p ∂g/∂q`.
fn poisson_bracket(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[f, g, q, p] = args else {
        return None;
    };
    let f = best(cx.graph, f)?;
    let g = best(cx.graph, g)?;
    let qs = vector(cx.graph, q).unwrap_or_else(|| vec![q]);
    let ps = vector(cx.graph, p).unwrap_or_else(|| vec![p]);
    if qs.len() != ps.len() {
        return None;
    }
    let mut terms = Vec::with_capacity(2 * qs.len());
    for (&q, &p) in qs.iter().zip(&ps) {
        cx.graph.symbol_of(q)?;
        cx.graph.symbol_of(p)?;
        let (fq, gp) = (derivative(cx.graph, f, q)?, derivative(cx.graph, g, p)?);
        let (fp, gq) = (derivative(cx.graph, f, p)?, derivative(cx.graph, g, q)?);
        let first = mul(cx.graph, &[fq, gp]);
        let second = mul(cx.graph, &[fp, gq]);
        terms.push(sub(cx.graph, first, second));
    }
    let sum = add(cx.graph, &terms);
    Some(cx.simplify(sum))
}

/// `a^k = -Γ^k_ij u^i u^j` for `Γ` given as `[k][i][j]`.
fn geodesic_acceleration(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[gamma, u] = args else {
        return None;
    };
    let u = vector(cx.graph, u)?;
    let n = u.len();
    let planes = vector(cx.graph, gamma)?;
    if planes.len() != n {
        return None;
    }
    let mut out = Vec::with_capacity(n);
    for plane in planes {
        let rows = vector(cx.graph, plane)?;
        if rows.len() != n {
            return None;
        }
        let mut terms = Vec::new();
        for (i, row) in rows.into_iter().enumerate() {
            let entries = vector(cx.graph, row)?;
            if entries.len() != n {
                return None;
            }
            for (j, entry) in entries.into_iter().enumerate() {
                terms.push(mul(cx.graph, &[entry, u[i], u[j]]));
            }
        }
        let sum = add(cx.graph, &terms);
        let value = neg(cx.graph, sum);
        out.push(cx.simplify(value));
    }
    Some(list(cx.graph, &out))
}

/// The Schwarzschild metric in `list(t, r, theta, phi)`, signature
/// `(-, +, +, +)`.
fn schwarzschild(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[m, vars] = args else {
        return None;
    };
    let vars = vector(cx.graph, vars)?;
    let &[_, r, theta, _] = vars.as_slice() else {
        return None;
    };
    let graph = &mut *cx.graph;
    let (c, g) = (constant(graph, "c_0")?, constant(graph, "G_N")?);
    let sin = graph.ops().lookup("sin")?;
    // 1 - 2 G M / (c² r)
    let two = graph.int(2);
    let c2 = powi(graph, c, 2);
    let denominator = mul(graph, &[c2, r]);
    let inverse = powi(graph, denominator, -1);
    let rs_over_r = mul(graph, &[two, g, m, inverse]);
    let one = graph.int(1);
    let factor = sub(graph, one, rs_over_r);
    let zero = graph.int(0);
    let gtt = {
        let product = mul(graph, &[factor, c2]);
        neg(graph, product)
    };
    let grr = powi(graph, factor, -1);
    let r2 = powi(graph, r, 2);
    let s = graph.node(sin, &[theta]);
    let s2 = powi(graph, s, 2);
    let gpp = mul(graph, &[r2, s2]);
    let rows = vec![
        vec![gtt, zero, zero, zero],
        vec![zero, grr, zero, zero],
        vec![zero, zero, r2, zero],
        vec![zero, zero, zero, gpp],
    ];
    Some(matrix_term(graph, &rows))
}

// ----------------------------------------------------------------------
// Dirac algebra
// ----------------------------------------------------------------------

/// `γ^mu` in the Dirac representation (`mu = 5` for `γ⁵`).
fn gamma_matrix(
    graph: &mut Graph,
    mu: i64,
) -> Option<Vec<Vec<NodeId>>> {
    let i = imaginary_unit(graph)?;
    let (zero, one, minus_one) = (graph.int(0), graph.int(1), graph.int(-1));
    let minus_i = neg(graph, i);
    // Pauli blocks.
    let sigma: [[NodeId; 4]; 3] = [[zero, one, one, zero], [zero, minus_i, i, zero], [one, zero, zero, minus_one]];
    let negate = |graph: &mut Graph, n: NodeId| {
        if n == zero {
            zero
        } else {
            let negated = neg(graph, n);
            best(graph, negated).unwrap_or(negated)
        }
    };
    let rows = match mu {
        | 0 => vec![
            vec![one, zero, zero, zero],
            vec![zero, one, zero, zero],
            vec![zero, zero, minus_one, zero],
            vec![zero, zero, zero, minus_one],
        ],
        | 1..=3 => {
            let s = sigma[usize::try_from(mu - 1).ok()?];
            let (n0, n1, n2, n3) = (negate(graph, s[0]), negate(graph, s[1]), negate(graph, s[2]), negate(graph, s[3]));
            vec![
                vec![zero, zero, s[0], s[1]],
                vec![zero, zero, s[2], s[3]],
                vec![n0, n1, zero, zero],
                vec![n2, n3, zero, zero],
            ]
        },
        | 5 => vec![
            vec![zero, zero, one, zero],
            vec![zero, zero, zero, one],
            vec![one, zero, zero, zero],
            vec![zero, one, zero, zero],
        ],
        | _ => return None,
    };
    Some(rows)
}

fn mat_vec(
    graph: &mut Graph,
    m: &[Vec<NodeId>],
    v: &[NodeId],
) -> Vec<NodeId> {
    m.iter()
        .map(|row| {
            let nonzero: Vec<(NodeId, NodeId)> = row
                .iter()
                .zip(v)
                .filter(|&(&e, _)| graph.number_of(e).is_none_or(|n| !n.is_zero()))
                .map(|(&e, &x)| (e, x))
                .collect();
            let terms: Vec<NodeId> = nonzero.into_iter().map(|(e, x)| mul(graph, &[e, x])).collect();
            add(graph, &terms)
        })
        .collect()
}

/// `p̸ = γ^0 p^0 - γ^1 p^1 - γ^2 p^2 - γ^3 p^3`.
fn slash(
    graph: &mut Graph,
    p: &[NodeId],
) -> Option<Vec<Vec<NodeId>>> {
    if p.len() != 4 {
        return None;
    }
    let mut sum: Vec<Vec<Vec<NodeId>>> = Vec::with_capacity(4);
    for (mu, &component) in p.iter().enumerate() {
        let gamma = gamma_matrix(graph, i64::try_from(mu).ok()?)?;
        let sign = graph.int(if mu == 0 { 1 } else { -1 });
        sum.push(gamma.iter().map(|row| row.iter().map(|&e| mul(graph, &[sign, component, e])).collect()).collect());
    }
    let mut out = vec![vec![NodeId::NONE; 4]; 4];
    for (a, row) in out.iter_mut().enumerate() {
        for (b, entry) in row.iter_mut().enumerate() {
            let terms: Vec<NodeId> = sum.iter().map(|m| m[a][b]).collect();
            *entry = add(graph, &terms);
        }
    }
    Some(out)
}

/// `(iħ γ^μ ∂_μ - m c) ψ` with `∂_0 = ∂_t / c`; `vars = list(t, x, ...)`.
fn dirac_equation(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[psi, m, vars] = args else {
        return None;
    };
    let psi = vector(cx.graph, psi)?;
    let vars = vector(cx.graph, vars)?;
    if psi.len() != 4 || vars.is_empty() || vars.len() > 4 {
        return None;
    }
    let graph = &mut *cx.graph;
    let (hbar, c, i) = (constant(graph, "hbar")?, constant(graph, "c_0")?, imaginary_unit(graph)?);
    let mut total: Vec<NodeId> = Vec::new();
    for (mu, &x) in vars.iter().enumerate() {
        graph.symbol_of(x)?;
        let gamma = gamma_matrix(graph, i64::try_from(mu).ok()?)?;
        let mut d_psi = Vec::with_capacity(4);
        for &component in &psi {
            d_psi.push(derivative(graph, component, x)?);
        }
        let mut term = mat_vec(graph, &gamma, &d_psi);
        let factor = if mu == 0 {
            let inverse = powi(graph, c, -1);
            mul(graph, &[i, hbar, inverse])
        } else {
            mul(graph, &[i, hbar])
        };
        for t in &mut term {
            *t = mul(graph, &[factor, *t]);
        }
        if total.is_empty() {
            total = term;
        } else {
            for (acc, t) in total.iter_mut().zip(term) {
                *acc = add(graph, &[*acc, t]);
            }
        }
    }
    let mut out = Vec::with_capacity(4);
    for (acc, &component) in total.into_iter().zip(&psi) {
        let mass_term = mul(cx.graph, &[m, c, component]);
        let value = sub(cx.graph, acc, mass_term);
        out.push(cx.simplify(value));
    }
    Some(list(cx.graph, &out))
}

/// `(1/c²) ∂²φ/∂t² - ∇²φ + (m c / ħ)² φ` with `vars = list(t, x, ...)`.
fn klein_gordon(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[phi, m, vars] = args else {
        return None;
    };
    let phi = best(cx.graph, phi)?;
    let vars = vector(cx.graph, vars)?;
    let (&t, space) = vars.split_first()?;
    let graph = &mut *cx.graph;
    let (hbar, c) = (constant(graph, "hbar")?, constant(graph, "c_0")?);
    let second = |graph: &mut Graph, x: NodeId| -> Option<NodeId> {
        graph.symbol_of(x)?;
        let first = derivative(graph, phi, x)?;
        derivative(graph, first, x)
    };
    let mut terms = Vec::new();
    let dtt = second(graph, t)?;
    let inverse_c2 = powi(graph, c, -2);
    terms.push(mul(graph, &[inverse_c2, dtt]));
    for &x in space {
        let dxx = second(graph, x)?;
        terms.push(neg(graph, dxx));
    }
    let inverse_hbar = powi(graph, hbar, -1);
    let ratio = mul(graph, &[m, c, inverse_hbar]);
    let squared = powi(graph, ratio, 2);
    terms.push(mul(graph, &[squared, phi]));
    let sum = add(graph, &terms);
    Some(cx.simplify(sum))
}

/// `ψ̄ = ψ† γ⁰` as a row (list).
fn dirac_adjoint(
    cx: &mut Cx<'_>,
    psi: NodeId,
) -> Option<NodeId> {
    let psi = vector(cx.graph, psi)?;
    if psi.len() != 4 {
        return None;
    }
    let conj = cx.graph.ops().lookup("conj")?;
    let mut out = Vec::with_capacity(4);
    for (k, &component) in psi.iter().enumerate() {
        let c = cx.graph.node(conj, &[component]);
        let value = if k < 2 { c } else { neg(cx.graph, c) };
        out.push(cx.simplify(value));
    }
    Some(list(cx.graph, &out))
}

/// `a·b = a⁰b⁰ - a¹b¹ - a²b² - a³b³`.
fn minkowski_dot(
    graph: &mut Graph,
    a: NodeId,
    b: NodeId,
) -> Option<NodeId> {
    let (a, b) = (vector(graph, a)?, vector(graph, b)?);
    if a.len() != b.len() || a.is_empty() {
        return None;
    }
    let terms: Vec<NodeId> = a
        .iter()
        .zip(&b)
        .enumerate()
        .map(|(k, (&x, &y))| {
            let product = mul(graph, &[x, y]);
            if k == 0 { product } else { neg(graph, product) }
        })
        .collect();
    Some(add(graph, &terms))
}

/// `½ ∂_μφ ∂^μφ - ½ m² φ²` (natural units) with `vars = list(t, x, ...)`.
fn scalar_lagrangian(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[phi, m, vars] = args else {
        return None;
    };
    let phi = best(cx.graph, phi)?;
    let vars = vector(cx.graph, vars)?;
    let graph = &mut *cx.graph;
    let mut terms = Vec::new();
    for (k, &x) in vars.iter().enumerate() {
        graph.symbol_of(x)?;
        let d = derivative(graph, phi, x)?;
        let square = powi(graph, d, 2);
        terms.push(if k == 0 { square } else { neg(graph, square) });
    }
    let m2 = powi(graph, m, 2);
    let phi2 = powi(graph, phi, 2);
    let mass = mul(graph, &[m2, phi2]);
    terms.push(neg(graph, mass));
    let sum = add(graph, &terms);
    let half = graph.num(crate::graph::Number::fraction(1, 2)?);
    let value = mul(graph, &[half, sum]);
    Some(cx.simplify(value))
}

/// `i (p̸ + m) / (p·p - m²)`.
fn fermion_propagator(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[p, m] = args else {
        return None;
    };
    let components = vector(cx.graph, p)?;
    let graph = &mut *cx.graph;
    let mut matrix = slash(graph, &components)?;
    for (k, row) in matrix.iter_mut().enumerate() {
        row[k] = add(graph, &[row[k], m]);
    }
    let square = minkowski_dot(graph, p, p)?;
    let m2 = powi(graph, m, 2);
    let denominator = sub(graph, square, m2);
    let inverse = powi(graph, denominator, -1);
    let i = imaginary_unit(graph)?;
    let factor = mul(graph, &[i, inverse]);
    let mut rows = Vec::with_capacity(4);
    for row in matrix {
        let mut out = Vec::with_capacity(4);
        for e in row {
            let value = mul(cx.graph, &[factor, e]);
            out.push(cx.simplify(value));
        }
        rows.push(out);
    }
    Some(matrix_term(cx.graph, &rows))
}

/// `Z = Σ exp(-E_i / (k_B T))`.
fn partition_function(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[energies, t] = args else {
        return None;
    };
    let energies = vector(cx.graph, energies)?;
    let graph = &mut *cx.graph;
    let exp = graph.ops().lookup("exp")?;
    let k = constant(graph, "k_B")?;
    let kt = mul(graph, &[k, t]);
    let inverse = powi(graph, kt, -1);
    let terms: Vec<NodeId> = energies
        .into_iter()
        .map(|e| {
            let ratio = mul(graph, &[e, inverse]);
            let exponent = neg(graph, ratio);
            graph.node(exp, &[exponent])
        })
        .collect();
    Some(add(graph, &terms))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Facts;
    use crate::rules::testing::eval;
    use crate::rules::testing::numeric;
    use crate::rules::testing::reduce_with;

    /// Symbols that are real (or positive) in every test.
    const ASSUME: [(&str, Facts); 10] = [
        ("x", Facts::REAL),
        ("y", Facts::REAL),
        ("z", Facts::REAL),
        ("t", Facts::REAL),
        ("k", Facts::REAL),
        ("a", Facts::POSITIVE),
        ("m", Facts::POSITIVE),
        ("w", Facts::POSITIVE),
        ("r", Facts::POSITIVE),
        ("q", Facts::REAL),
    ];

    fn run(src: &str) -> String {
        let (text, reduced) = reduce_with(&[physics()], src, &ASSUME);
        assert!(reduced, "`{src}` was not fully reduced: {text}");
        text
    }

    fn same(
        a: &str,
        b: &str,
    ) {
        assert_eq!(run(a), run(b), "{a} != {b}");
    }

    fn value(src: &str) -> f64 {
        numeric(&[physics()], src, &[], 1e-12).0
    }

    #[test]
    fn constants_are_symbolic_and_numeric() {
        assert_eq!(run("mass_energy_equivalence(m)"), "m*c_0^2");
        let e = value("mass_energy_equivalence(1)");
        assert!((e - 299_792_458.0_f64.powi(2)).abs() / e < 1e-12);
        assert!((eval(&[physics()], "hbar * 2 * pi / h_planck", &[]) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn classical_mechanics() {
        same("kinetic_energy(m, v)", "m*v^2/2");
        same("potential_energy_gravity_universal(m1, m2, r, G_N)", "-G_N*m1*m2/r");
        same("lagrangian(kinetic_energy(m, v), potential_energy_spring(k, x))", "m*v^2/2 - k*x^2/2");
        same("torque(list(1, 0, 0), list(0, 1, 0))", "list(0, 0, 1)");
        same("work_constant_force(list(F, 0, 0), list(d, 0, 0))", "F*d");
        same("kinematics(list(t^2, 3*t, 0), t)", "list(list(t^2, 3*t, 0), list(2*t, 3, 0), list(2, 0, 0))");
        same(
            "euler_lagrange_equation(lagrangian(kinetic_energy(m, diff(q(t), t)), potential_energy_spring(k, q(t))), q(t), t)",
            "k*q(t) + m*diff(diff(q(t), t), t)",
        );
        same("poisson_bracket(x, p, x, p)", "1");
        same("poisson_bracket(x^2*p, p^2, x, p)", "4*x*p^2");
        same("poisson_bracket(x1*p2, x2*p1, list(x1, x2), list(p1, p2))", "x2*p2 - x1*p1");
        same("hamilton_equations(p^2/(2*m) + k*x^2/2, x, p)", "list(p/m, -k*x)");
    }

    #[test]
    fn electromagnetism() {
        // A uniform field: E = -grad(-E0 x) = (E0, 0, 0).
        same("electric_field_from_potential(-E0*x, list(x, y, z))", "list(E0, 0, 0)");
        same("magnetic_field_from_vector_potential(list(-y/2, x/2, 0), list(x, y, z))", "list(0, 0, 1)");
        same("lorentz_force(q, list(0, 0, 0), list(v, 0, 0), list(0, 0, B))", "list(0, -B*q*v, 0)");
        // A plane wave in vacuum satisfies the source-free equations
        // when its speed is 1/sqrt(epsilon_0 mu_0).
        let maxwell = run(
            "maxwell_equations(list(0, cos(z - w*t), 0), list(-cos(z - w*t)/w, 0, 0), 0, list(0, 0, 0), list(x, y, z), t)",
        );
        assert!(maxwell.starts_with("list(0, 0, list(0, 0, 0), list(0, "), "{maxwell}");
        same("poynting_vector(list(Ex, 0, 0), list(0, B, 0))", "list(0, 0, B*Ex/mu_0)");
        same("coulomb_force(q, q, r)", "q^2/(4*pi*epsilon_0*r^2)");
    }

    #[test]
    fn relativity() {
        same("lorentz_factor(0)", "1");
        same("velocity_addition(c_0, v)", "c_0");
        let gamma = value("lorentz_factor(0.6*c_0)");
        assert!((gamma - 1.25).abs() < 1e-12);
        same("lorentz_transformation_x(x, t, 0)", "list(x, t)");
        same("schwarzschild_radius(M)", "2*G_N*M/c_0^2");
        // The flat limit of the Schwarzschild metric.
        let g = run("schwarzschild_metric(0, list(t, r, theta, phi))");
        assert_eq!(g, "list(list(-c_0^2, 0, 0, 0), list(0, 1, 0, 0), list(0, 0, r^2, 0), list(0, 0, 0, r^2*sin(theta)^2))");
        // Geodesics of the plane in polar coordinates: a_r = r θ'^2.
        same(
            "geodesic_acceleration(christoffel(list(list(1, 0), list(0, r^2)), list(r, theta)), list(u, v))",
            "list(r*v^2, -2*u*v/r)",
        );
    }

    #[test]
    fn quantum_operators() {
        same("qm_apply(momentum_operator(x), exp(I*k*x))", "hbar*k*exp(k*x*I)");
        same("commutator(position_operator(x), momentum_operator(x), f(x))", "hbar*f(x)*I");
        same("matrix_commutator(pauli_x, pauli_y)", "list(list(2*I, 0), list(0, -2*I))");
        same("probability_density(exp(I*k*x))", "1");
        // The ground state of the harmonic oscillator.
        let psi = "exp(-m*w*x^2/(2*hbar))";
        same(&format!("energy_eigenvalue(hamiltonian_harmonic_oscillator(m, w, x), {psi}, x)"), "hbar*w/2");
        same(&format!("time_independent_schrodinger_equation(hamiltonian_harmonic_oscillator(m, w, x), {psi}, hbar*w/2)"), "0");
        // Expectation values with a normalised Gaussian.
        same("expectation_value(position_operator(x^2), exp(-a*x^2), x)", "1/(4*a)");
        same("expectation_value(momentum_operator(x), exp(-a*x^2), x)", "0");
        same("uncertainty(position_operator(x), exp(-a*x^2), x)", "1/(2*a^(1/2))");
        same("braket(exp(-x^2/2), exp(-x^2/2), x)", "pi^(1/2)");
    }

    #[test]
    fn spin_and_dirac_algebra() {
        same("spin_operator(pauli_z)", "list(list(hbar/2, 0), list(0, -hbar/2))");
        // γ⁰γ⁰ = 1, {γ¹, γ¹} = -2.
        same("matmul(dirac_gamma(0), dirac_gamma(0))", "identity(4)");
        same("matrix_anticommutator(dirac_gamma(1), dirac_gamma(1))", "smul(-2, identity(4))");
        same("matrix_anticommutator(dirac_gamma(0), dirac_gamma(2))", "zeros(4, 4)");
        same("feynman_slash(list(en, 0, 0, 0))", "smul(en, dirac_gamma(0))");
        same("minkowski_dot(list(en, p, 0, 0), list(en, p, 0, 0))", "en^2 - p^2");
        // A spinor at rest: (iħ/c ∂_t - m c) e^(-i m c² t / ħ) (1, 0, 0, 0) = 0.
        same("dirac_equation(list(exp(-I*m*c_0^2*t/hbar), 0, 0, 0), m, list(t, x))", "list(0, 0, 0, 0)");
        // A plane wave with the relativistic dispersion relation.
        same("klein_gordon_equation(exp(I*(k*x - (c_0^2*k^2 + m^2*c_0^4/hbar^2)^(1/2)*t)), m, list(t, x))", "0");
        same("dirac_adjoint(list(a1, a2, a3, a4))", "list(conj(a1), conj(a2), -conj(a3), -conj(a4))");
        same("scalar_field_lagrangian(f(t, x), m, list(t, x))", "(diff(f(t, x), t)^2 - diff(f(t, x), x)^2 - m^2*f(t, x)^2)/2");
    }

    #[test]
    fn solid_state() {
        same("lattice_volume(list(a, 0, 0), list(0, a, 0), list(0, 0, a))", "a^3");
        same(
            "reciprocal_lattice_vectors(list(a, 0, 0), list(0, a, 0), list(0, 0, a))",
            "list(list(2*pi/a, 0, 0), list(0, 2*pi/a, 0), list(0, 0, 2*pi/a))",
        );
        same("hall_coefficient(n, -q_e)", "1/(-n*q_e)");
        same("energy_band(0, m, E0)", "E0");
        // Debye's law at high temperature tends to 3 N k_B (Dulong–Petit).
        let c = value("debye_heat_capacity(1, 1, 1000) / k_B");
        assert!((c - 3.0).abs() < 1e-4, "{c}");
        let e = value("einstein_heat_capacity(1, 1, 1000) / k_B");
        assert!((e - 3.0).abs() < 1e-4, "{e}");
    }

    #[test]
    fn thermodynamics() {
        same("carnot_efficiency(300, 600)", "1/2");
        same("gibbs_free_energy(enthalpy(U, P, V), T, S)", "U + P*V - T*S");
        same("verify_maxwell_relation_helmholtz(T^2*ln(V), T, V)", "0");
        same("partition_function(list(0, en), T)", "1 + exp(-en/(k_B*T))");
        // Two-level system: U = k T² ∂ ln Z / ∂T.
        let u = run("internal_energy_from_partition(1 + exp(-en/(k_B*T)), T)");
        let at = eval(&[physics()], &u, &[("en", 2.0), ("T", 3.0e23)]);
        let expected = 2.0 / (1.0 + (2.0 / (crate::constant::BOLTZMANN_CONSTANT * 3.0e23)).exp());
        assert!((at - expected).abs() < 1e-12, "{u}: {at} vs {expected}");
        same("entropy_from_helmholtz(-k_B*T*ln(T), T)", "k_B + k_B*ln(T)");
        let ratio = value("fermi_dirac_distribution(1, 1, 300)");
        assert!((ratio - 0.5).abs() < 1e-15);
    }
}
