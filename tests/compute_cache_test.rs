use rssn::compute::cache::{
    ComputationResultCache,
    ParsingCache,
};
use rssn::symbolic::core::Expr;
use std::sync::Arc;

#[test]

fn test_parsing_cache() {

    let cache = ParsingCache::new();

    let input = "x + 1";

    let expr = Arc::new(
        Expr::new_variable("x"),
    ); // Dummy expr for test

    assert!(cache
        .get(input)
        .is_none());

    cache.set(
        input.to_string(),
        expr.clone(),
    );

    let cached = cache.get(input);

    assert!(cached.is_some());

    // Note: We can't easily compare Arc<Expr> for equality without dereferencing,
    // but here we just check if it's cached.

    cache.clear();

    assert!(cache
        .get(input)
        .is_none());
}

#[test]

fn test_computation_result_cache() {

    let cache =
        ComputationResultCache::new();

    let expr = Arc::new(
        Expr::new_variable("x"),
    );

    let value = "42.0".to_string();

    assert!(cache
        .get(&expr)
        .is_none());

    cache.set(
        expr.clone(),
        value.clone(),
    );

    let cached = cache.get(&expr);

    assert!(cached.is_some());

    if let Some(v) = cached {

        assert_eq!(v, "42.0");
    } else {

        panic!("Unexpected value type");
    }

    cache.clear();

    assert!(cache
        .get(&expr)
        .is_none());
}
