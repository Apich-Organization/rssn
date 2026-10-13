//! The continuum and particle engines by scenario name with JSON parameters.

use rssn::sim::scenario::{ENGINE_SCENARIOS, SCENARIOS, run};
use serde_json::Value;

fn go(name: &str, params: &str) -> Value {
    let out = run(name, params).unwrap_or_else(|e| panic!("{name}: {e}"));
    serde_json::from_str(&out).unwrap()
}

fn f(v: &Value, key: &str) -> f64 {
    v[key].as_f64().unwrap_or_else(|| panic!("missing number `{key}` in {v}"))
}

#[test]
fn every_engine_scenario_is_listed_and_defaults_are_valid_requests() {
    for name in ENGINE_SCENARIOS {
        assert!(SCENARIOS.contains(name), "{name} missing from SCENARIOS");
    }
    assert_eq!(ENGINE_SCENARIOS.len(), 18);
    // an empty object selects the documented defaults for the cheap scenarios
    for name in ["poisson_multigrid", "fem_simplex", "wave_nd"] {
        let v = go(name, "{}");
        assert!(v.is_object());
    }
    assert!(run("heat_nd", "").is_ok(), "an empty string means {{}}");
}

#[test]
fn continuum_scenarios_return_accurate_results() {
    let v = go("poisson_multigrid", r#"{"dim": 3, "n": 17, "cycle": "fmg"}"#);
    assert_eq!(v["converged"], true);
    assert!(f(&v, "max_error") < 0.02 && f(&v, "convergence_factor") < 0.3, "{v}");

    for bc in ["dirichlet", "neumann", "periodic"] {
        let v = go("heat_nd", &format!(r#"{{"dim": 1, "n": 40, "steps": 40, "bc": "{bc}"}}"#));
        assert!(f(&v, "max_error") < 2e-3, "{bc}: {v}");
    }
    let v = go("heat_nd", r#"{"dim": 2, "n": 25, "adaptive_tol": 1e-6}"#);
    assert!(f(&v, "max_error") < 5e-3 && v["stats"]["accepted"].as_u64().unwrap() > 0, "{v}");

    let v = go("wave_nd", r#"{"dim": 3, "n": 13, "t_end": 0.2}"#);
    assert!(f(&v, "max_error") < 0.02, "{v}");

    for dim in 1..=3 {
        let v = go("euler_fv", &format!(r#"{{"dim": {dim}, "n": 64, "recon": "muscl_vanleer", "riemann": "roe"}}"#));
        assert!(f(&v, "l1_error") < 0.02, "dim {dim}: {}", f(&v, "l1_error"));
    }
    let v = go("dg_conservation", r#"{"law": "advection", "ne": 10, "p": 3}"#);
    assert!(f(&v, "l2_error") < 1e-3, "{v}");
    let v = go("dg_conservation", r#"{"law": "euler", "ne": 40, "p": 2, "limiter": true, "t_end": 0.2}"#);
    assert!(f(&v, "l2_error") < 0.1, "{v}");
    assert!(run("dg_conservation", r#"{"law": "magic"}"#).is_err());

    let v = go("acoustics_pml", r#"{"dim": 1, "n": 160, "pml_cells": 16}"#);
    assert!(f(&v, "energy_ratio") < 1e-3, "{v}");
    let v = go("acoustics_pml", r#"{"dim": 1, "n": 160, "absorbing": false}"#);
    assert!(f(&v, "energy_ratio") > 0.9, "{v}");

    let v = go("yee_fdtd", r#"{"dim": 1, "n": 100, "steps": 300}"#);
    assert!(f(&v, "energy") < 1e-3 * f(&v, "peak_energy"), "{v}");
}

#[test]
fn fluid_quantum_and_structural_scenarios() {
    let v = go("navier_stokes_mac", r#"{"case": "taylor_green", "n": 16, "t_end": 0.5, "dt": 0.05, "re": 10.0}"#);
    assert!(f(&v, "max_divergence") < 1e-8, "{v}");
    let v = go("navier_stokes_mac", r#"{"case": "cavity", "n": 16, "t_end": 1.0, "dt": 0.02}"#);
    assert!(f(&v, "kinetic_energy") > 0.0);
    let v = go("navier_stokes_mac", r#"{"case": "poiseuille", "n": 16, "re": 10.0, "t_end": 20.0, "dt": 0.1}"#);
    let line = v["u_line"].as_array().unwrap();
    assert!(line[8].as_f64().unwrap() > 1.1, "centreline speed {}", line[8]);
    assert!(run("navier_stokes_mac", r#"{"case": "vortex"}"#).is_err());

    for dim in 1..=3 {
        let n = if dim == 3 { 16 } else { 64 };
        let v = go("split_step", &format!(r#"{{"dim": {dim}, "n": {n}, "len": 12.0, "dt": 0.05, "steps": 1000}}"#));
        let exact = 0.5 * dim as f64;
        assert!((f(&v, "energy") - exact).abs() < 0.01 * exact, "dim {dim}: {v}");
        assert!((f(&v, "norm") - 1.0).abs() < 1e-9);
    }
    let v = go("split_step", r#"{"mode": "evolve", "dt": 0.02, "steps": 100, "g": 1.0}"#);
    assert!((f(&v, "norm") - go_norm_ref()).abs() < 1e-6, "{v}");

    for (dim, degree) in [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1)] {
        let n = if dim == 3 { 4 } else { 8 };
        let v = go("fem_simplex", &format!(r#"{{"dim": {dim}, "degree": {degree}, "n": {n}}}"#));
        assert!(f(&v, "l2_error") < if dim == 3 { 0.3 } else { 0.03 }, "dim {dim} P{degree}: {v}");
    }
    let v = go("elasticity_fem", r#"{"dim": 2, "degree": 2, "length_cells": 20}"#);
    let tip = v["tip_displacement"][1].as_f64().unwrap();
    let beam = f(&v, "euler_bernoulli_deflection");
    assert!(tip < 0.0 && (-tip / beam - 1.0).abs() < 0.12, "{v}");
}

fn go_norm_ref() -> f64 {
    // the evolve scenario starts from a unit Gaussian displaced by one length unit: norm is pi^{d/2}
    std::f64::consts::PI.sqrt()
}

#[test]
fn bem_and_particle_scenarios() {
    for (problem, method) in [("interior", "aca"), ("exterior", "aca"), ("interior", "p1"), ("exterior", "p1")] {
        let v = go("bem_laplace", &format!(r#"{{"subdivisions": 2, "problem": "{problem}", "method": "{method}"}}"#));
        assert!(f(&v, "max_flux_error") < 0.2, "{problem} {method}: {v}");
    }
    let v = go("sph", r#"{"t_end": 0.05, "grad_h": true}"#);
    assert!(f(&v, "energy_drift").abs() < 1e-3, "{v}");
    let v = go("dpd", r#"{"dim": 3, "box_len": 4.5, "steps": 150}"#);
    assert!((f(&v, "temperature") - 1.0).abs() < 0.1, "{v}");
    assert!(v["momentum"].as_array().unwrap().iter().all(|m| m.as_f64().unwrap().abs() < 1e-8));
    let v = go("md_lj", r#"{"dim": 3, "cells": 6, "steps": 60}"#);
    assert!(f(&v, "energy_drift").abs() < 0.02, "{v}");
    let v = go("md_lj", r#"{"dim": 2, "cells": 8, "steps": 60, "thermostat": "langevin"}"#);
    assert!(f(&v, "temperature") > 0.1);
    // spin models in 1-3D: hot lattices are disordered, cold ones ordered
    for dim in 1..=3 {
        let hot = go("spin_model", &format!(r#"{{"dim": {dim}, "n": 8, "beta": 0.05, "sweeps": 50}}"#));
        assert!(f(&hot, "abs_magnetisation_per_site") < 0.5, "{hot}");
    }
    let cold = go("spin_model", r#"{"dim": 2, "n": 12, "beta": 1.0, "sweeps": 100, "update": "wolff"}"#);
    assert!(f(&cold, "abs_magnetisation_per_site") > 0.95, "{cold}");
    let xy = go("spin_model", r#"{"dim": 2, "n": 8, "model": "xy", "update": "metropolis", "beta": 2.0, "sweeps": 50}"#);
    assert!(f(&xy, "energy_per_site") < 0.0);
    assert!(run("spin_model", r#"{"model": "potts7"}"#).is_err());

    let v = go("nbody_gravity", r#"{"positions": [[0,0,0],[2,0,0]], "masses": [1.0, 3.0], "softening": 0.0, "theta": 0.1}"#);
    let a = v["accelerations"].as_array().unwrap();
    let (a0, a1) = (a[0][0].as_f64().unwrap(), a[1][0].as_f64().unwrap());
    assert!((a0 - 3.0 / 4.0).abs() < 1e-9 && (a1 + 1.0 / 4.0).abs() < 1e-9, "{v}");

    let v = go("kerr_geodesic", r#"{"radius": 10.0, "steps": 200}"#);
    assert!(f(&v, "constraint").abs() < 1e-8 && (f(&v, "r") - 10.0).abs() < 1e-3, "{v}");
    assert!(run("kerr_geodesic", r#"{"radius": 1.0, "spin": 0.5}"#).is_err());
}
