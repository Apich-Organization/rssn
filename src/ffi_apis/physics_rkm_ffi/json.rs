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
use crate::physics::physics_rkm::OdeSystem;
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

#[no_mangle]

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

#[no_mangle]

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

#[no_mangle]

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

#[no_mangle]

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
