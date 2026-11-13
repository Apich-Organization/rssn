//! Prints the typed scenario catalog (`rssn::sim::schema`) as JSON.
//!
//! `cargo run --profile tools --example scenario_catalog -- --bless` re-runs
//! every scenario on its defaults or example and rewrites
//! `src/sim/scenario_results.json`, the witnessed result shapes.

fn main() {
    if std::env::args().any(|a| a == "--bless") {
        let witnessed = rssn::sim::schema::witness().unwrap_or_else(|e| panic!("{e}"));
        let text = serde_json::to_string_pretty(&witnessed).expect("serialises");
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src/sim/scenario_results.json");
        std::fs::write(&path, format!("{text}\n")).expect("writes");
        println!("wrote {}", path.display());
        return;
    }
    println!("{}", rssn::sim::schema::catalog_json());
}
