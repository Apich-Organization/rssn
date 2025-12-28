//! Computation engine and related infrastructure.
//!
//! This module provides a comprehensive framework for managing and executing
//! computational tasks asynchronously. It includes:
//!
//! - **Engine**: The main computation engine for submitting and managing tasks
//! - **Computation**: Task representation with status, progress, and results
//! - **State**: Computation state management
//! - **Cache**: Parsing and result caching for performance
//! - **Computable**: Trait for types that can be computed
//!
//! # Architecture
//!
//! The compute module is designed around an asynchronous, event-driven architecture:
//!
//! ```text
//! ┌─────────────────┐
//! │  ComputeEngine  │  ← Main entry point
//! └────────┬────────┘
//!          │
//!          ├─→ ParsingCache        (String → Expr)
//!          ├─→ ResultCache         (Expr → Value)
//!          └─→ Computations        (ID → Computation)
//!                    │
//!                    ├─→ Status    (Pending/Running/Completed/Failed)
//!                    ├─→ Progress  (Percentage + Description)
//!                    ├─→ State     (Intermediate values)
//!                    └─→ Result    (Final value)
//! ```
//!
//! # Usage
//!
//! ## Basic Computation
//!
//! ```
//! 
//! use rssn::compute::engine::ComputeEngine;
//!
//! // Create an engine
//! let engine = ComputeEngine::new();
//!
//! // Submit a computation
//! let id = engine
//!     .parse_and_submit("2 + 2")
//!     .unwrap();
//!
//! // Check status
//! if let Some(status) = engine.get_status(&id) {
//!
//!     println!(
//!         "Status: {:?}",
//!         status
//!     );
//! }
//!
//! // Wait for result
//! std::thread::sleep(std::time::Duration::from_secs(6));
//!
//! if let Some(result) = engine.get_result(&id) {
//!
//!     println!("Result: {}", result);
//! }
//! ```
//!
//! ## Advanced: Pause/Resume/Cancel
//!
//! ```
//! 
//! use rssn::compute::engine::ComputeEngine;
//!
//! let engine = ComputeEngine::new();
//!
//! let id = engine
//!     .parse_and_submit("complex_calculation")
//!     .unwrap();
//!
//! // Pause the computation
//! engine.pause(&id);
//!
//! // Do something else...
//!
//! // Resume when ready
//! engine.resume(&id);
//!
//! // Or cancel if no longer needed
//! engine.cancel(&id);
//! ```
//!
//! ## Using Caches
//!
//! ```
//! 
//! use std::sync::Arc;
//!
//! use rssn::compute::cache::ComputationResultCache;
//! use rssn::compute::cache::ParsingCache;
//! use rssn::symbolic::core::Expr;
//!
//! // Parsing cache
//! let parsing_cache = ParsingCache::new();
//!
//! let expr = Arc::new(Expr::Constant(42.0));
//!
//! parsing_cache.set(
//!     "my_expr".to_string(),
//!     expr.clone(),
//! );
//!
//! // Later...
//! if let Some(cached_expr) = parsing_cache.get("my_expr")
//! {
//!
//!     println!("Found cached expression");
//! }
//!
//! // Result cache
//! let result_cache = ComputationResultCache::new();
//!
//! result_cache.set(
//!     expr.clone(),
//!     "42".to_string(),
//! );
//! ```
//!
//! # Thread Safety
//!
//! All types in this module are thread-safe and can be shared across threads:
//! - `ComputeEngine` uses `Arc<RwLock<...>>` for internal state
//! - `ParsingCache` and `ComputationResultCache` use `Mutex` for synchronization
//! - `Computation` uses `Arc<Mutex<...>>` for safe concurrent access
//!
//! # Performance
//!
//! The module is optimized for:
//! - **Caching**: Avoid re-parsing and re-computing identical expressions
//! - **Async execution**: Non-blocking computation using Rayon thread pool
//! - **Minimal locking**: Fine-grained locks to reduce contention
//!
//! # Examples
//!
//! See individual module documentation for more examples:
//! - [`engine`] - Computation engine
//! - [`computation`] - Computation task representation
//! - [`cache`] - Caching infrastructure
//! - [`state`] - State management
//! - [`computable`] - Computable trait

/// Caching for parsing and computation results.
pub mod cache;
/// Trait for computable mathematical objects.
pub mod computable;
/// Task representation and tracking.
pub mod computation;
pub mod engine;
/// State management for computations.
pub mod state;
