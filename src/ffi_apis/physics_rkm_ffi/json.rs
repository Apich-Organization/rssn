//! JSON-based FFI API for physics RKM (Runge-Kutta Methods) functions.

use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_rkm::BogackiShampine23;
use crate::physics::physics_rkm::CashKarp45;
use crate::physics::physics_rkm::DormandPrince54;
use crate::physics::physics_rkm::{
    self,
};

#[derive(Deserialize)]

struct LorenzInput {
    sigma: f64,
    rho: f64,
    beta: f64,
    y0: Vec<f64>,
    t_span: (f64, f64),
    dt_initial: f64,
    tol: (f64, f64),
}

#[derive(Deserialize)]

struct DampedOscillatorInput {
    omega: f64,
    zeta: f64,
    y0: Vec<f64>,
    t_span: (f64, f64),
    dt: f64,
}

#[derive(Deserialize)]

struct VanDerPolInput {
    mu: f64,
    y0: Vec<f64>,
    t_span: (f64, f64),
    dt_initial: f64,
    tol: (f64, f64),
}

#[derive(Deserialize)]

struct LotkaVolterraInput {
    alpha: f64,
    beta: f64,
    delta: f64,
    gamma: f64,
    y0: Vec<f64>,
    t_span: (f64, f64),
    dt_initial: f64,
    tol: (f64, f64),
}

#[derive(Serialize)]

struct OdeResult {
    time: Vec<f64>,
    states: Vec<Vec<f64>>,
}

/// Solves the Lorenz system using adaptive Dormand-Prince RK5(4) method via JSON serialization.
///
/// The Lorenz system is a chaotic dynamical system defined by:
/// dx/dt = σ(y - x), dy/dt = x(ρ - z) - y, dz/dt = xy - βz.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `sigma`: Prandtl number σ
///   - `rho`: Rayleigh number ρ
///   - `beta`: Geometric parameter β
///   - `y0`: Initial state [x₀, y₀, z₀]
///   - `t_span`: Time interval [`t_start`, `t_end`]
///   - `dt_initial`: Initial time step size
///   - `tol`: Error tolerances [absolute, relative]
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<OdeResult, String>` with
/// `time` and `states` arrays.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_physics_rkm_lorenz_json(
    input: *const c_char
) -> *mut c_char {

    let input : LorenzInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    OdeResult,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        },
    };

    let system =
        physics_rkm::LorenzSystem {
            sigma: input.sigma,
            rho: input.rho,
            beta: input.beta,
        };

    let solver = DormandPrince54::new();

    let results = solver.solve(
        &system,
        &input.y0,
        input.t_span,
        input.dt_initial,
        input.tol,
    );

    let mut time = Vec::with_capacity(
        results.len(),
    );

    let mut states = Vec::with_capacity(
        results.len(),
    );

    for (t, y) in results {

        time.push(t);

        states.push(y);
    }

    to_c_string(
        serde_json::to_string(
            &FfiResult::<
                OdeResult,
                String,
            >::ok(
                OdeResult {
                    time,
                    states,
                },
            ),
        )
        .unwrap(),
    )
}

/// Solves the damped oscillator system using RK4 method via JSON serialization.
///
/// The damped harmonic oscillator is defined by:
/// d²x/dt² + 2ζωdx/dt + ω²x = 0.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `omega`: Natural frequency ω
///   - `zeta`: Damping ratio ζ
///   - `y0`: Initial state [x₀, v₀]
///   - `t_span`: Time interval [`t_start`, `t_end`]
///   - `dt`: Time step size
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<OdeResult, String>` with
/// `time` and `states` arrays.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_physics_rkm_damped_oscillator_json(
    input: *const c_char
) -> *mut c_char {

    let input : DampedOscillatorInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    OdeResult,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        },
    };

    let system = physics_rkm::DampedOscillatorSystem {
        omega : input.omega,
        zeta : input.zeta,
    };

    let results =
        physics_rkm::solve_rk4(
            &system,
            &input.y0,
            input.t_span,
            input.dt,
        );

    let mut time = Vec::with_capacity(
        results.len(),
    );

    let mut states = Vec::with_capacity(
        results.len(),
    );

    for (t, y) in results {

        time.push(t);

        states.push(y);
    }

    to_c_string(
        serde_json::to_string(
            &FfiResult::<
                OdeResult,
                String,
            >::ok(
                OdeResult {
                    time,
                    states,
                },
            ),
        )
        .unwrap(),
    )
}

/// Solves the Van der Pol oscillator using adaptive Cash-Karp RK4(5) method via JSON serialization.
///
/// The Van der Pol equation models nonlinear oscillations with self-excitation:
/// d²x/dt² - μ(1 - x²)dx/dt + x = 0.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `mu`: Nonlinearity parameter μ
///   - `y0`: Initial state [x₀, v₀]
///   - `t_span`: Time interval [`t_start`, `t_end`]
///   - `dt_initial`: Initial time step size
///   - `tol`: Error tolerances [absolute, relative]
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<OdeResult, String>` with
/// `time` and `states` arrays.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_physics_rkm_vanderpol_json(
    input: *const c_char
) -> *mut c_char {

    let input : VanDerPolInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    OdeResult,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        },
    };

    let system =
        physics_rkm::VanDerPolSystem {
            mu: input.mu,
        };

    let solver = CashKarp45::default();

    let results = solver.solve(
        &system,
        &input.y0,
        input.t_span,
        input.dt_initial,
        input.tol,
    );

    let mut time = Vec::with_capacity(
        results.len(),
    );

    let mut states = Vec::with_capacity(
        results.len(),
    );

    for (t, y) in results {

        time.push(t);

        states.push(y);
    }

    to_c_string(
        serde_json::to_string(
            &FfiResult::<
                OdeResult,
                String,
            >::ok(
                OdeResult {
                    time,
                    states,
                },
            ),
        )
        .unwrap(),
    )
}

/// Solves the Lotka-Volterra predator-prey system using Bogacki-Shampine RK2(3) via JSON serialization.
///
/// The Lotka-Volterra equations model population dynamics:
/// dx/dt = αx - βxy, dy/dt = δxy - γy.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `alpha`: Prey growth rate α
///   - `beta`: Predation rate β
///   - `delta`: Predator efficiency δ
///   - `gamma`: Predator death rate γ
///   - `y0`: Initial state [prey₀, predator₀]
///   - `t_span`: Time interval [`t_start`, `t_end`]
///   - `dt_initial`: Initial time step size
///   - `tol`: Error tolerances [absolute, relative]
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<OdeResult, String>` with
/// `time` and `states` arrays.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_physics_rkm_lotka_volterra_json(
    input: *const c_char
) -> *mut c_char {

    let input : LotkaVolterraInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    OdeResult,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        },
    };

    let system = physics_rkm::LotkaVolterraSystem {
        alpha : input.alpha,
        beta : input.beta,
        delta : input.delta,
        gamma : input.gamma,
    };

    let solver =
        BogackiShampine23::default();

    let results = solver.solve(
        &system,
        &input.y0,
        input.t_span,
        input.dt_initial,
        input.tol,
    );

    let mut time = Vec::with_capacity(
        results.len(),
    );

    let mut states = Vec::with_capacity(
        results.len(),
    );

    for (t, y) in results {

        time.push(t);

        states.push(y);
    }

    to_c_string(
        serde_json::to_string(
            &FfiResult::<
                OdeResult,
                String,
            >::ok(
                OdeResult {
                    time,
                    states,
                },
            ),
        )
        .unwrap(),
    )
}
