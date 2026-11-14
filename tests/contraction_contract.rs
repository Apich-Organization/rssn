use rssn::kernels::contractions::ContractionError;
use rssn::kernels::contractions::contract_rectangular;
use std::cell::Cell;

#[test]
fn scalar_failure_stops_the_contraction_before_any_further_host_work() {
    let calls = Cell::new(0);
    let result = contract_rectangular(
        2,
        3,
        2,
        || 0i32,
        |_, _, _, _| {
            calls.set(calls.get() + 1);
            Err("scalar failure")
        },
    );
    assert_eq!(result, Err(ContractionError::Element("scalar failure")));
    assert_eq!(calls.get(), 1);
}

#[test]
fn overflow_is_rejected_before_callbacks_and_empty_inner_uses_the_host_identity() {
    let calls = Cell::new(0);
    let bad = contract_rectangular(
        usize::MAX,
        1,
        2,
        || 0i32,
        |_, _, _, _| {
            calls.set(calls.get() + 1);
            Ok::<_, ()>(1)
        },
    );
    assert_eq!(bad, Err(ContractionError::ShapeOverflow));
    assert_eq!(calls.get(), 0);
    let zeros = contract_rectangular(
        2,
        0,
        2,
        || 42i32,
        |_, _, _, _| {
            calls.set(calls.get() + 1);
            Ok::<_, ()>(1)
        },
    )
    .unwrap();
    assert_eq!(zeros, vec![42; 4]);
    assert_eq!(calls.get(), 0);
}
