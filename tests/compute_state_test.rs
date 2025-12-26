use rssn::compute::state::State;

#[test]

fn test_state_new() {

    let state = State::new();

    assert_eq!(
        state.intermediate_value,
        ""
    );
}

#[test]

fn test_state_default() {

    let state = State::default();

    assert_eq!(
        state.intermediate_value,
        ""
    );
}
