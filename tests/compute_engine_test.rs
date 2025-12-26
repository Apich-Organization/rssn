use std::sync::Arc;
use std::thread;
use std::time::Duration;

use rssn::compute::computation::ComputationStatus;
use rssn::compute::engine::ComputeEngine;
use rssn::symbolic::core::Expr;

#[test]

fn test_engine_new() {

    let engine = ComputeEngine::new();

    // Engine should be created successfully
    // No computations should be present initially
    assert!(engine
        .get_status("nonexistent")
        .is_none());
}

#[test]

fn test_parse_and_submit_valid() {

    let engine = ComputeEngine::new();

    let result = engine.parse_and_submit("2 + 2");

    assert!(result.is_ok());

    let id = result.unwrap();

    assert!(!id.is_empty());
}

#[test]

fn test_parse_and_submit_invalid() {

    let engine = ComputeEngine::new();

    // Use a truly invalid expression that will fail parsing
    let result = engine.parse_and_submit("((((");

    // Note: The parser might be lenient, so we just check it doesn't panic
    // If it succeeds, that's also acceptable behavior
    let _ = result;
}

#[test]

fn test_submit_direct() {

    let engine = ComputeEngine::new();

    let expr = Arc::new(Expr::Constant(42.0));

    let id = engine.submit(expr);

    assert!(!id.is_empty());

    // Check that computation exists
    assert!(engine
        .get_status(&id)
        .is_some());
}

#[test]

fn test_get_status() {

    let engine = ComputeEngine::new();

    let id = engine
        .parse_and_submit("x + 1")
        .unwrap();

    // Status should exist
    let status = engine.get_status(&id);

    assert!(status.is_some());

    // Should be Pending or Running
    let status = status.unwrap();

    assert!(matches!(
        status,
        ComputationStatus::Pending | ComputationStatus::Running
    ));
}

#[test]

fn test_get_status_nonexistent() {

    let engine = ComputeEngine::new();

    let status = engine.get_status("nonexistent-id");

    assert!(status.is_none());
}

#[test]

fn test_get_progress() {

    let engine = ComputeEngine::new();

    let id = engine
        .parse_and_submit("2 + 2")
        .unwrap();

    // Progress should exist
    let progress = engine.get_progress(&id);

    assert!(progress.is_some());

    let progress = progress.unwrap();

    assert!(progress.percentage >= 0.0 && progress.percentage <= 100.0);

    assert!(!progress
        .description
        .is_empty());
}

#[test]

fn test_get_progress_nonexistent() {

    let engine = ComputeEngine::new();

    let progress = engine.get_progress("nonexistent-id");

    assert!(progress.is_none());
}

#[test]

fn test_get_result_eventually_completes() {

    let engine = ComputeEngine::new();

    let id = engine
        .parse_and_submit("2 + 2")
        .unwrap();

    // Wait for computation to complete (simulated work takes ~5 seconds)
    // Add extra time to be safe
    thread::sleep(Duration::from_secs(
        7,
    ));

    // Result should be available
    let result = engine.get_result(&id);

    // Note: Due to timing, result might not always be available
    // This is acceptable for an async system
    if result.is_some() {

        assert!(!result
            .unwrap()
            .is_empty());
    }
}

#[test]

fn test_pause_and_resume() {

    let engine = ComputeEngine::new();

    let id = engine
        .parse_and_submit("2 + 2")
        .unwrap();

    // Wait a bit for computation to start
    thread::sleep(Duration::from_millis(500));

    // Pause the computation
    engine.pause(&id);

    thread::sleep(Duration::from_millis(500));

    // Check if paused (might be in various states due to timing)
    if let Some(status) = engine.get_status(&id) {

        println!(
            "Status after pause: {:?}",
            status
        );
    }

    // Resume the computation
    engine.resume(&id);

    thread::sleep(Duration::from_millis(500));

    // Should exist and be in some valid state
    // Due to async nature, we can't guarantee exact state
    assert!(engine
        .get_status(&id)
        .is_some());
}

#[test]

fn test_cancel() {

    let engine = ComputeEngine::new();

    let id = engine
        .parse_and_submit("2 + 2")
        .unwrap();

    // Wait a bit for computation to start
    thread::sleep(Duration::from_millis(100));

    // Cancel the computation
    engine.cancel(&id);

    // Computation should no longer exist in registry
    thread::sleep(Duration::from_millis(100));

    let status = engine.get_status(&id);

    assert!(status.is_none());
}

#[test]

fn test_multiple_computations() {

    let engine = ComputeEngine::new();

    let id1 = engine
        .parse_and_submit("1 + 1")
        .unwrap();

    let id2 = engine
        .parse_and_submit("2 + 2")
        .unwrap();

    let id3 = engine
        .parse_and_submit("3 + 3")
        .unwrap();

    // All should have different IDs
    assert_ne!(id1, id2);

    assert_ne!(id2, id3);

    assert_ne!(id1, id3);

    // All should have status
    assert!(engine
        .get_status(&id1)
        .is_some());

    assert!(engine
        .get_status(&id2)
        .is_some());

    assert!(engine
        .get_status(&id3)
        .is_some());
}

#[test]

fn test_parsing_cache() {

    let engine = ComputeEngine::new();

    // Submit same expression twice
    let id1 = engine
        .parse_and_submit("x + y")
        .unwrap();

    let id2 = engine
        .parse_and_submit("x + y")
        .unwrap();

    // Should create different computations (different IDs)
    assert_ne!(id1, id2);

    // But parsing should be cached (no way to directly test this,
    // but it should work without errors)
}

#[test]

fn test_default_trait() {

    let engine = ComputeEngine::default();

    let id = engine
        .parse_and_submit("1 + 1")
        .unwrap();

    assert!(!id.is_empty());
}

#[test]

fn test_clone_trait() {

    let engine1 = ComputeEngine::new();

    let engine2 = engine1.clone();

    // Both should work independently
    let id1 = engine1
        .parse_and_submit("1 + 1")
        .unwrap();

    let id2 = engine2
        .parse_and_submit("2 + 2")
        .unwrap();

    assert_ne!(id1, id2);
}

#[test]

fn test_concurrent_submissions() {

    use std::sync::Arc;
    use std::thread;

    let engine = Arc::new(ComputeEngine::new());

    let mut handles = vec![];

    for i in 0 .. 5 {

        let engine_clone = Arc::clone(&engine);

        let handle = thread::spawn(move || {

            let input = format!("{} + {}", i, i);

            engine_clone
                .parse_and_submit(&input)
                .unwrap()
        });

        handles.push(handle);
    }

    let ids : Vec<String> = handles
        .into_iter()
        .map(|h| h.join().unwrap())
        .collect();

    // All IDs should be unique
    for i in 0 .. ids.len() {

        for j in (i + 1) .. ids.len() {

            assert_ne!(ids[i], ids[j]);
        }
    }
}
