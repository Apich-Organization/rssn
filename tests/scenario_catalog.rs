//! Catalog coverage and deterministic witness generation for native scenarios.

use rssn::sim::scenario::ENGINE_SCENARIOS;
use rssn::sim::scenario::SCENARIOS;
use rssn::sim::scenario::run;
use rssn::sim::schema::catalog;
use rssn::sim::schema::witness;
use serde_json::Value;
use std::collections::BTreeSet;

#[test]
fn catalog_lists_every_dispatch_name_once() {
    let listed = catalog();
    let actual: Vec<_> = listed.iter().map(|entry| entry.name.as_str()).collect();
    assert_eq!(actual, SCENARIOS);
    assert!(actual.ends_with(ENGINE_SCENARIOS));
    assert_eq!(actual.iter().collect::<BTreeSet<_>>().len(), actual.len());

    for name in SCENARIOS {
        let error = run(name, "{").expect_err("malformed JSON was accepted");
        assert!(error.contains("invalid parameters"), "{name}: {error}");
    }
    assert!(run("nbody_gravity", "{}").is_err());
    let nbody = listed
        .iter()
        .find(|entry| entry.name == "nbody_gravity")
        .unwrap();
    assert!(!nbody.all_defaults);
    assert!(run(&nbody.name, &nbody.example).is_ok());
}

#[test]
fn witness_in_a_parallel_pool_matches_the_checked_in_schemas() {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(2)
        .build()
        .unwrap();
    let observed = pool.install(witness).unwrap();
    let expected: Value =
        serde_json::from_str(include_str!("../src/sim/scenario_results.json")).unwrap();
    assert_eq!(observed, expected);
}
