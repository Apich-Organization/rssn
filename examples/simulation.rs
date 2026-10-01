//! Running simulations by name with JSON parameters — the same entry point
//! the C interface exposes as `rssn_sim_run`.
//!
//! Run with `cargo run --release --example simulation`.

use rssn::sim::scenario;

fn main() -> Result<(), String> {
    println!("scenarios: {}", scenario::SCENARIOS.join(", "));
    let ising = scenario::run(
        "ising",
        r#"{"width": 32, "height": 32, "temperature": 1.8, "mc_steps": 200, "seed": 42}"#,
    )?;
    let value: serde_json::Value = serde_json::from_str(&ising).map_err(|e| e.to_string())?;
    println!("Ising 32x32 at T = 1.8: magnetization {}", value["magnetization"]);

    let orbit = scenario::run(
        "geodesic",
        r#"{"black_hole_mass": 1.0, "initial_state": [10.0, 0.0, 0.0, 0.038], "proper_time_end": 500.0, "initial_dt": 0.1}"#,
    )?;
    let value: serde_json::Value = serde_json::from_str(&orbit).map_err(|e| e.to_string())?;
    let points = value["orbit"].as_array().map_or(0, Vec::len);
    println!("Schwarzschild geodesic: {points} points");
    Ok(())
}
