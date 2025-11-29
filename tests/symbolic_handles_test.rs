use rssn::symbolic::core::Expr;
use rssn::symbolic::handles::HANDLE_MANAGER;

#[test]
fn test_handle_insert_and_get() {
    HANDLE_MANAGER.clear();
    
    let expr = Expr::new_variable("x");
    let handle = HANDLE_MANAGER.insert(expr.clone());
    
    assert!(handle > 0, "Handle should be non-zero");
    
    let retrieved = HANDLE_MANAGER.get(handle);
    assert!(retrieved.is_some(), "Should retrieve inserted expression");
    
    if let Some(arc_expr) = retrieved {
        assert_eq!(format!("{}", arc_expr), "x");
    }
}

#[test]
fn test_handle_exists() {
    HANDLE_MANAGER.clear();
    
    let expr = Expr::new_constant(42.0);
    let handle = HANDLE_MANAGER.insert(expr);
    
    assert!(HANDLE_MANAGER.exists(handle), "Handle should exist");
    assert!(!HANDLE_MANAGER.exists(99999), "Non-existent handle should not exist");
}

#[test]
fn test_handle_free() {
    HANDLE_MANAGER.clear();
    
    let expr = Expr::new_variable("y");
    let handle = HANDLE_MANAGER.insert(expr);
    
    assert!(HANDLE_MANAGER.exists(handle), "Handle should exist before free");
    
    let freed = HANDLE_MANAGER.free(handle);
    assert!(freed.is_some(), "Free should return the expression");
    
    assert!(!HANDLE_MANAGER.exists(handle), "Handle should not exist after free");
    assert!(HANDLE_MANAGER.get(handle).is_none(), "Get should return None after free");
}

#[test]
fn test_handle_count() {
    HANDLE_MANAGER.clear();
    
    let initial_count = HANDLE_MANAGER.count();
    assert_eq!(initial_count, 0, "Should start with 0 handles after clear");
    
    let h1 = HANDLE_MANAGER.insert(Expr::new_constant(1.0));
    assert_eq!(HANDLE_MANAGER.count(), initial_count + 1);
    
    let h2 = HANDLE_MANAGER.insert(Expr::new_constant(2.0));
    assert_eq!(HANDLE_MANAGER.count(), initial_count + 2);
    
    HANDLE_MANAGER.free(h1);
    assert_eq!(HANDLE_MANAGER.count(), initial_count + 1);
    
    HANDLE_MANAGER.free(h2);
    assert_eq!(HANDLE_MANAGER.count(), initial_count);
}

#[test]
fn test_handle_clear() {
    HANDLE_MANAGER.clear();
    
    let h1 = HANDLE_MANAGER.insert(Expr::new_constant(1.0));
    let h2 = HANDLE_MANAGER.insert(Expr::new_constant(2.0));
    let h3 = HANDLE_MANAGER.insert(Expr::new_constant(3.0));
    
    assert_eq!(HANDLE_MANAGER.count(), 3);
    
    HANDLE_MANAGER.clear();
    
    assert_eq!(HANDLE_MANAGER.count(), 0);
    assert!(!HANDLE_MANAGER.exists(h1));
    assert!(!HANDLE_MANAGER.exists(h2));
    assert!(!HANDLE_MANAGER.exists(h3));
}

#[test]
fn test_handle_clone_expr() {
    HANDLE_MANAGER.clear();
    
    let expr = Expr::new_add(Expr::new_variable("x"), Expr::new_constant(5.0));
    let handle = HANDLE_MANAGER.insert(expr);
    
    let cloned = HANDLE_MANAGER.clone_expr(handle);
    assert!(cloned.is_some(), "Should clone expression");
    
    if let Some(cloned_expr) = cloned {
        let original = HANDLE_MANAGER.get(handle).unwrap();
        assert_eq!(format!("{}", cloned_expr), format!("{}", original));
    }
}

#[test]
fn test_handle_get_all_handles() {
    HANDLE_MANAGER.clear();
    
    let h1 = HANDLE_MANAGER.insert(Expr::new_constant(1.0));
    let h2 = HANDLE_MANAGER.insert(Expr::new_constant(2.0));
    let h3 = HANDLE_MANAGER.insert(Expr::new_constant(3.0));
    
    let all_handles = HANDLE_MANAGER.get_all_handles();
    
    assert_eq!(all_handles.len(), 3);
    assert!(all_handles.contains(&h1));
    assert!(all_handles.contains(&h2));
    assert!(all_handles.contains(&h3));
}

#[test]
fn test_handle_unique_ids() {
    HANDLE_MANAGER.clear();
    
    let h1 = HANDLE_MANAGER.insert(Expr::new_constant(1.0));
    let h2 = HANDLE_MANAGER.insert(Expr::new_constant(1.0)); // Same expression
    
    assert_ne!(h1, h2, "Different handles should be generated even for same expression");
}

#[test]
fn test_handle_thread_safety() {
    use std::sync::Arc;
    use std::thread;
    
    HANDLE_MANAGER.clear();
    
    let handles = Arc::new(std::sync::Mutex::new(Vec::new()));
    let mut threads = vec![];
    
    // Spawn multiple threads inserting expressions
    for i in 0..10 {
        let handles_clone = Arc::clone(&handles);
        let thread = thread::spawn(move || {
            let expr = Expr::new_constant(i as f64);
            let handle = HANDLE_MANAGER.insert(expr);
            handles_clone.lock().unwrap().push(handle);
        });
        threads.push(thread);
    }
    
    // Wait for all threads to complete
    for thread in threads {
        thread.join().unwrap();
    }
    
    let handles = handles.lock().unwrap();
    assert_eq!(handles.len(), 10, "All threads should have inserted expressions");
    
    // Verify all handles are unique
    let mut sorted_handles = handles.clone();
    sorted_handles.sort();
    sorted_handles.dedup();
    assert_eq!(sorted_handles.len(), 10, "All handles should be unique");
    
    // Verify all handles still exist
    for &handle in handles.iter() {
        assert!(HANDLE_MANAGER.exists(handle), "Handle {} should exist", handle);
    }
}

#[test]
fn test_handle_complex_expression() {
    HANDLE_MANAGER.clear();
    
    // Create a complex expression: (x^2 + y) * sin(z)
    let x = Expr::new_variable("x");
    let y = Expr::new_variable("y");
    let z = Expr::new_variable("z");
    
    let x_squared = Expr::new_pow(x, Expr::new_constant(2.0));
    let sum = Expr::new_add(x_squared, y);
    let sin_z = Expr::new_sin(z);
    let complex_expr = Expr::new_mul(sum, sin_z);
    
    let handle = HANDLE_MANAGER.insert(complex_expr);
    
    let retrieved = HANDLE_MANAGER.get(handle);
    assert!(retrieved.is_some());
    
    if let Some(arc_expr) = retrieved {
        let expr_str = format!("{}", arc_expr);
        assert!(expr_str.contains("x"));
        assert!(expr_str.contains("y"));
        assert!(expr_str.contains("z"));
        assert!(expr_str.contains("sin"));
    }
}

#[test]
fn test_handle_persistence_across_operations() {
    HANDLE_MANAGER.clear();
    
    let expr1 = Expr::new_variable("a");
    let expr2 = Expr::new_variable("b");
    let expr3 = Expr::new_variable("c");
    
    let h1 = HANDLE_MANAGER.insert(expr1);
    let h2 = HANDLE_MANAGER.insert(expr2);
    let h3 = HANDLE_MANAGER.insert(expr3);
    
    // Free middle handle
    HANDLE_MANAGER.free(h2);
    
    // Other handles should still work
    assert!(HANDLE_MANAGER.exists(h1), "Handle h1 should exist");
    assert!(!HANDLE_MANAGER.exists(h2), "Handle h2 should not exist after free");
    assert!(HANDLE_MANAGER.exists(h3), "Handle h3 should exist");
    
    // Verify we can still get the expressions
    let retrieved_h1 = HANDLE_MANAGER.get(h1);
    assert!(retrieved_h1.is_some(), "Should be able to get h1");
    if let Some(expr) = retrieved_h1 {
        assert_eq!(format!("{}", expr), "a");
    }
    
    let retrieved_h3 = HANDLE_MANAGER.get(h3);
    assert!(retrieved_h3.is_some(), "Should be able to get h3");
    if let Some(expr) = retrieved_h3 {
        assert_eq!(format!("{}", expr), "c");
    }
}
