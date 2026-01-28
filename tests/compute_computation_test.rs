use std::sync::Arc;
use std::sync::Condvar;
use std::sync::Mutex;
use std::sync::atomic::AtomicBool;

use rssn::compute::computation::Computation;
use rssn::compute::computation::ComputationProgress;
use rssn::compute::computation::ComputationStatus;
use rssn::compute::state::State;
use rssn::symbolic::core::Expr;

#[test]

fn test_computation_creation() {

    let expr =
        Arc::new(Expr::new_constant(1.0));

    let computation = Computation {
        id: "test_id".to_string(),
        expr: expr.clone(),
        status:
            ComputationStatus::Pending,
        progress: ComputationProgress {
            percentage: 0.0,
            description: "Init"
                .to_string(),
        },
        result: None,
        state: State::new(),
        pause: Arc::new((
            Mutex::new(false),
            Condvar::new(),
        )),
        cancel_signal: Arc::new(
            AtomicBool::new(false),
        ),
    };

    assert_eq!(
        computation.id,
        "test_id"
    );

    assert_eq!(
        computation.status,
        ComputationStatus::Pending
    );

    assert_eq!(
        computation
            .progress
            .percentage,
        0.0
    );
}

#[test]

fn test_computation_serialization() {

    let expr =
        Arc::new(Expr::new_constant(1.0));

    let computation = Computation {
        id: "test_id".to_string(),
        expr: expr.clone(),
        status:
            ComputationStatus::Completed,
        progress: ComputationProgress {
            percentage: 100.0,
            description: "Done"
                .to_string(),
        },
        result: Some("42".to_string()),
        state: State::new(),
        pause: Arc::new((
            Mutex::new(false),
            Condvar::new(),
        )),
        cancel_signal: Arc::new(
            AtomicBool::new(false),
        ),
    };

    let serialized =
        serde_json::to_string(
            &computation,
        )
        .unwrap();

    let deserialized: Computation =
        serde_json::from_str(
            &serialized,
        )
        .unwrap();

    assert_eq!(
        deserialized.id,
        computation.id
    );

    assert_eq!(
        deserialized.status,
        computation.status
    );

    assert_eq!(
        deserialized.result,
        computation.result
    );
    // Note: pause and cancel_signal are skipped in serialization, so they will be default/missing?
    // Wait, Computation doesn't implement Default.
    // serde will use Default::default() for skipped fields if they implement Default, or error if not?
    // Arc<(Mutex<bool>, Condvar)> does not implement Default.
    // Ah, `#[serde(skip)]` requires the field to implement `Default` OR use `default` attribute.
    // If it doesn't implement Default, deserialization will fail unless we provide a default.
    // Let's check if `Computation` deserialization works.
    // If it fails, I need to add `#[serde(default)]` or implement Default for the skipped fields wrapper or use a function.
}
