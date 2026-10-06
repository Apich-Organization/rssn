//! Scenarios by name, with JSON parameters and results.
//!
//! The unified simulation entry point: a scenario name plus JSON
//! parameters in, JSON results out. This is what the C interface and other
//! language bindings call; Rust code can use the typed functions of
//! [`models`](super::models) directly.

use ndarray::Array2;
use num_complex::Complex;
use serde::Deserialize;
use serde::de::DeserializeOwned;
use serde_json::Value;
use serde_json::json;

use super::models::fdtd_electrodynamics;
use super::models::geodesic_relativity;
use super::models::gpe_superfluidity;
use super::models::ising_statistical;
use super::models::linear_elasticity;
use super::models::navier_stokes_fluid;
use super::models::schrodinger_quantum;

/// The scenario names [`run`] accepts.
pub const SCENARIOS: &[&str] = &[
    "fdtd",
    "geodesic",
    "gpe",
    "ising",
    "elasticity",
    "lid_driven_cavity",
    "channel_flow",
    "schrodinger",
    "ising_wolff",
    "sod_shock_tube",
    "kepler_symplectic",
];

fn parse<T: DeserializeOwned>(params: &str) -> Result<T, String> {
    serde_json::from_str(params).map_err(|e| format!("invalid parameters: {e}"))
}

fn encode(value: &impl serde::Serialize) -> Result<Value, String> {
    serde_json::to_value(value).map_err(|e| e.to_string())
}

#[derive(Deserialize)]
struct IsingRequest {
    #[serde(flatten)]
    params: ising_statistical::IsingParameters,
    seed: Option<u64>,
}

#[derive(Deserialize)]
struct WolffRequest {
    width: usize,
    height: usize,
    temperature: f64,
    thermalisation: usize,
    sweeps: usize,
    seed: Option<u64>,
}

#[derive(Deserialize)]
struct SodRequest {
    cells: usize,
    t: f64,
    limiter: Option<String>,
}

#[derive(Deserialize)]
struct KeplerRequest {
    q: Vec<f64>,
    p: Vec<f64>,
    dt: f64,
    steps: usize,
    order: Option<u32>,
}

#[derive(Deserialize)]
struct ChannelRequest {
    nx: usize,
    ny: usize,
    re: f64,
    dt: f64,
    n_iter: usize,
    obstacle_mask: Array2<bool>,
}

#[derive(Deserialize)]
struct SchrodingerRequest {
    params: schrodinger_quantum::SchrodingerParameters,
    /// The initial wave function, flattened row-major, as `[re, im]` pairs.
    initial_psi: Vec<Complex<f64>>,
}

fn flow(output: navier_stokes_fluid::NavierStokesOutput) -> Result<Value, String> {
    let (u, v, p) = output?;
    Ok(json!({ "u": encode(&u)?, "v": encode(&v)?, "p": encode(&p)? }))
}

/// Runs the scenario `name` with parameters given as a JSON object and
/// returns its results as JSON.
///
/// Arrays are encoded the way `ndarray` serialises them:
/// `{"v": 1, "dim": [rows, cols], "data": [...]}`.
///
/// # Errors
/// An unknown scenario, parameters that do not parse, or a failure of the
/// simulation itself.
pub fn run(
    name: &str,
    params: &str,
) -> Result<String, String> {
    let value = match name {
        | "fdtd" => {
            let frames = fdtd_electrodynamics::run_fdtd_simulation(&parse(params)?);
            json!({ "frames": encode(&frames)? })
        },
        | "geodesic" => {
            let orbit = geodesic_relativity::run_geodesic_simulation(&parse(params)?);
            json!({ "orbit": encode(&orbit)? })
        },
        | "gpe" => {
            let density = gpe_superfluidity::run_gpe_ground_state_finder(&parse(params)?)?;
            json!({ "density": encode(&density)? })
        },
        | "ising" => {
            let request: IsingRequest = parse(params)?;
            let (spins, magnetization) = match request.seed {
                | Some(seed) => ising_statistical::run_ising_simulation_seeded(&request.params, seed),
                | None => ising_statistical::run_ising_simulation(&request.params),
            };
            json!({ "spins": spins, "magnetization": magnetization })
        },
        | "elasticity" => {
            let displacements = linear_elasticity::run_elasticity_simulation(&parse(params)?)?;
            json!({ "displacements": displacements })
        },
        | "lid_driven_cavity" => flow(navier_stokes_fluid::run_lid_driven_cavity(&parse(params)?))?,
        | "channel_flow" => {
            let r: ChannelRequest = parse(params)?;
            flow(navier_stokes_fluid::run_channel_flow(r.nx, r.ny, r.re, r.dt, r.n_iter, &r.obstacle_mask))?
        },
        | "ising_wolff" => {
            let r: WolffRequest = parse(params)?;
            let stats = ising_statistical::run_wolff_simulation(r.width, r.height, r.temperature, r.thermalisation, r.sweeps, r.seed.unwrap_or(ising_statistical::DEFAULT_ISING_SEED));
            encode(&stats)?
        },
        | "sod_shock_tube" => {
            let r: SodRequest = parse(params)?;
            let limiter = match r.limiter.as_deref() {
                | None | Some("van_leer") => super::gas_dynamics::Limiter::VanLeer,
                | Some("minmod") => super::gas_dynamics::Limiter::Minmod,
                | Some("none") => super::gas_dynamics::Limiter::None,
                | Some(other) => return Err(format!("unknown limiter `{other}`")),
            };
            let (x, density, exact) = super::gas_dynamics::sod_shock_tube(r.cells, r.t, limiter);
            json!({ "x": x, "density": density, "exact": exact })
        },
        | "kepler_symplectic" => {
            let r: KeplerRequest = parse(params)?;
            if r.q.len() != r.p.len() {
                return Err("q and p must have the same length".to_owned());
            }
            let path = super::integrators::symplectic_trajectory(&r.q, &r.p, 1.0, r.dt, r.steps, r.order.unwrap_or(4), |q: &[f64]| {
                let r3 = q.iter().map(|x| x * x).sum::<f64>().powf(1.5);
                q.iter().map(|x| -x / r3).collect()
            });
            json!({ "trajectory": path })
        },
        | "schrodinger" => {
            let mut r: SchrodingerRequest = parse(params)?;
            if r.initial_psi.len() != r.params.nx.saturating_mul(r.params.ny) {
                return Err("initial_psi must have nx * ny entries".to_owned());
            }
            let frames = schrodinger_quantum::run_schrodinger_simulation(&r.params, &mut r.initial_psi)?;
            json!({ "frames": encode(&frames)? })
        },
        | _ => {
            return Err(format!(
                "unknown scenario `{name}`; known: {}",
                SCENARIOS.join(", ")
            ));
        },
    };
    serde_json::to_string(&value).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scenarios_run_from_json() {
        let out = run("ising", r#"{"width": 8, "height": 8, "temperature": 1.5, "mc_steps": 20, "seed": 7}"#).unwrap();
        let v: Value = serde_json::from_str(&out).unwrap();
        assert_eq!(v["spins"].as_array().unwrap().len(), 64);
        assert_eq!(out, run("ising", r#"{"width": 8, "height": 8, "temperature": 1.5, "mc_steps": 20, "seed": 7}"#).unwrap());

        let out = run("ising_wolff", r#"{"width": 16, "height": 16, "temperature": 1.5, "thermalisation": 50, "sweeps": 200, "seed": 3}"#).unwrap();
        let v: Value = serde_json::from_str(&out).unwrap();
        assert!(v["magnetization"].as_f64().unwrap() > 0.9, "{out}");
        let out = run("sod_shock_tube", r#"{"cells": 50, "t": 0.2, "limiter": "minmod"}"#).unwrap();
        let v: Value = serde_json::from_str(&out).unwrap();
        assert_eq!(v["density"].as_array().unwrap().len(), 50);
        let out = run("kepler_symplectic", r#"{"q": [1.0, 0.0], "p": [0.0, 1.0], "dt": 0.01, "steps": 10}"#).unwrap();
        let v: Value = serde_json::from_str(&out).unwrap();
        assert_eq!(v["trajectory"].as_array().unwrap().len(), 11);

        let out = run(
            "geodesic",
            r#"{"black_hole_mass": 1.0, "initial_state": [10.0, 0.0, 0.0, 0.035], "proper_time_end": 50.0, "initial_dt": 0.1}"#,
        )
        .unwrap();
        let v: Value = serde_json::from_str(&out).unwrap();
        assert!(!v["orbit"].as_array().unwrap().is_empty());

        let out = run("lid_driven_cavity", r#"{"nx": 8, "ny": 8, "re": 10.0, "dt": 0.001, "n_iter": 5, "lid_velocity": 1.0}"#).unwrap();
        let v: Value = serde_json::from_str(&out).unwrap();
        assert_eq!(v["u"]["dim"], json!([8, 8]));
    }

    #[test]
    fn bad_requests_are_errors() {
        assert!(run("warp_drive", "{}").unwrap_err().contains("unknown scenario"));
        assert!(run("ising", "{\"width\": 1}").unwrap_err().contains("invalid parameters"));
    }
}
