use rssn::compute::computable::Computable;
use rssn::compute::computation::ComputationProgress;
use rssn::compute::state::State;

struct DummyComputable;

impl Computable for DummyComputable {
    fn compute(
        &self,
        _state: &mut State,
        progress: &mut ComputationProgress,
    ) -> Result<(), String> {

        progress.percentage = 100.0;

        progress.description =
            "Done".to_string();

        Ok(())
    }
}

#[test]

fn test_computable() {

    let computable = DummyComputable;

    let mut state = State::new();

    let mut progress =
        ComputationProgress {
            percentage: 0.0,
            description: "Starting"
                .to_string(),
        };

    let result = computable.compute(
        &mut state,
        &mut progress,
    );

    assert!(result.is_ok());

    assert_eq!(
        progress.percentage,
        100.0
    );

    assert_eq!(
        progress.description,
        "Done"
    );
}
