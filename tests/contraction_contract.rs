use rssn::kernels::contractions::ContractionError;
use rssn::kernels::contractions::contract_rectangular;
use rssn::kernels::contractions::contract_rectangular_append;
use rssn::kernels::contractions::contract_rectangular_for_each;
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

#[test]
fn streaming_preserves_row_major_first_product_order_for_host_scalars() {
    #[derive(Debug, PartialEq, Eq)]
    struct Physical {
        value: usize,
        unit: &'static str,
    }

    let mut calls = Vec::new();
    let mut cells = Vec::new();
    contract_rectangular_for_each(
        2,
        3,
        2,
        || panic!("nonempty products have no seeded identity"),
        |row, index, column, previous: Option<Physical>| {
            calls.push((row, column, index, previous.is_some()));
            let value = previous.map_or(0, |value| {
                assert_eq!(value.unit, "m s");
                value.value * 10
            });
            Ok::<_, ()>(Physical {
                value: value + 4 * row + column + index,
                unit: "m s",
            })
        },
        |row, column, value| {
            cells.push((row, column, value.value, value.unit));
            Ok(())
        },
    )
    .unwrap();
    assert_eq!(
        cells,
        [
            (0, 0, 12, "m s"),
            (0, 1, 123, "m s"),
            (1, 0, 456, "m s"),
            (1, 1, 567, "m s")
        ]
    );
    let expected: Vec<_> = (0..2)
        .flat_map(|row| {
            (0..2).flat_map(move |col| (0..3).map(move |index| (row, col, index, index != 0)))
        })
        .collect();
    assert_eq!(calls, expected);
}

#[test]
fn empty_axes_skip_all_callbacks_and_streaming_uses_per_cell_identities() {
    for (rows, cols) in [(usize::MAX, 0), (0, usize::MAX)] {
        let product = contract_rectangular::<i32, ()>(
            rows,
            usize::MAX,
            cols,
            || panic!(),
            |_, _, _, _| panic!(),
        );
        assert_eq!(product, Ok(vec![]));
        contract_rectangular_for_each::<i32, ()>(
            rows,
            usize::MAX,
            cols,
            || panic!(),
            |_, _, _, _| panic!(),
            |_, _, _| panic!(),
        )
        .unwrap();
    }
    let mut identity = 0;
    let mut cells = Vec::new();
    contract_rectangular_for_each(
        2,
        0,
        2,
        || {
            identity += 1;
            identity
        },
        |_, _, _, _| panic!("empty reductions have no scalar steps"),
        |row, column, value| {
            cells.push((row, column, value));
            Ok::<_, ()>(())
        },
    )
    .unwrap();
    assert_eq!(cells, [(0, 0, 1), (0, 1, 2), (1, 0, 3), (1, 1, 4)]);
}

#[test]
fn streaming_stops_at_consumer_error_even_when_the_whole_output_cannot_fit() {
    let mut scalar_calls = 0;
    let mut consumer_calls = 0;
    let result = contract_rectangular_for_each(
        usize::MAX,
        1,
        1,
        || 0u64,
        |row, index, column, previous| {
            assert_eq!((row, index, column, previous), (0, 0, 0, None));
            scalar_calls += 1;
            Ok(42)
        },
        |row, column, value| {
            assert_eq!((row, column, value), (0, 0, 42));
            consumer_calls += 1;
            Err("consumer stopped")
        },
    );
    assert_eq!(result, Err(ContractionError::Element("consumer stopped")));
    assert_eq!((scalar_calls, consumer_calls), (1, 1));
}

#[test]
fn a_scalar_error_preserves_only_already_consumed_cells() {
    let mut cells = Vec::new();
    let result = contract_rectangular_for_each(
        1,
        1,
        3,
        || 0,
        |_, _, column, _| {
            if column == 1 {
                Err("scalar failed")
            } else {
                Ok(column)
            }
        },
        |row, column, value| {
            cells.push((row, column, value));
            Ok(())
        },
    );
    assert_eq!(result, Err(ContractionError::Element("scalar failed")));
    assert_eq!(cells, [(0, 0, 0)]);
}

#[test]
fn appending_rolls_back_a_partial_result_and_reuses_capacity_on_retry() {
    let mut output = Vec::with_capacity(5);
    output.push("prefix".to_owned());
    let allocation = output.as_ptr();
    let failed =
        contract_rectangular_append(2, 1, 2, &mut output, String::new, |row, _, column, _| {
            if row == 1 {
                Err("scalar failure")
            } else {
                Ok(format!("{row}{column}"))
            }
        });
    assert_eq!(failed, Err(ContractionError::Element("scalar failure")));
    assert_eq!(output, ["prefix"]);
    contract_rectangular_append(2, 1, 2, &mut output, String::new, |row, _, column, _| {
        Ok::<_, ()>(format!("{row}{column}"))
    })
    .unwrap();
    assert_eq!(output, ["prefix", "00", "01", "10", "11"]);
    assert_eq!(output.as_ptr(), allocation);
}

#[test]
fn append_guard_restores_output_when_a_host_callback_unwinds() {
    let mut output = vec![7];
    let failure = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ = contract_rectangular_append::<i32, ()>(
            1,
            1,
            2,
            &mut output,
            || 0,
            |_, _, column, _| {
                assert_eq!(column, 0, "host panic");
                Ok(10)
            },
        );
    }));
    assert!(failure.is_err());
    assert_eq!(output, [7]);
}

#[test]
fn allocation_extents_are_rejected_before_scalar_work() {
    let too_many_u64 = isize::MAX as usize / size_of::<u64>() + 1;
    let result =
        contract_rectangular::<u64, ()>(too_many_u64, 1, 1, || panic!(), |_, _, _, _| panic!());
    assert_eq!(result, Err(ContractionError::ShapeOverflow));
    let mut unit = vec![()];
    let result = contract_rectangular_append::<(), ()>(
        usize::MAX,
        1,
        1,
        &mut unit,
        || panic!(),
        |_, _, _, _| panic!(),
    );
    assert_eq!(result, Err(ContractionError::ShapeOverflow));
    assert_eq!(unit, [()]);
    let result = contract_rectangular_for_each::<u64, ()>(
        usize::MAX,
        1,
        2,
        || panic!(),
        |_, _, _, _| panic!(),
        |_, _, _| panic!(),
    );
    assert_eq!(result, Err(ContractionError::ShapeOverflow));
}
