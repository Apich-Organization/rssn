//! Computation engine for managing and executing asynchronous computations.
//!
//! The `ComputeEngine` provides a high-level interface for submitting, tracking,
//! and managing computational tasks. It handles:
//! - Expression parsing and caching
//! - Asynchronous computation execution
//! - Progress tracking and status monitoring
//! - Pause/resume/cancel operations
//! - Result caching
//!
//! # Examples
//!
//! ```
//! 
//! use rssn::compute::engine::ComputeEngine;
//!
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
//! // Get result when complete
//! std::thread::sleep(std::time::Duration::from_secs(1));
//!
//! if let Some(result) = engine.get_result(&id) {
//!
//!     println!("Result: {}", result);
//! }
//! ```

#![allow(unused_imports)]

use std::collections::HashMap;
use std::io::prelude::*;
use std::sync::atomic::AtomicBool;
use std::sync::Arc;
use std::sync::Condvar;
use std::sync::Mutex;
use std::sync::RwLock;

/// Development in place.
use rayon::prelude::*;
use uuid::Uuid;

use crate::compute::cache::ComputationResultCache;
use crate::compute::cache::ParsingCache;
use crate::compute::computation::Computation;
use crate::compute::computation::ComputationProgress;
use crate::compute::computation::ComputationStatus;
use crate::compute::computation::Value;
use crate::compute::state::State;
use crate::symbolic::core::Expr;

/// A computation engine for managing asynchronous computations.
///
/// `ComputeEngine` maintains a registry of active computations and provides
/// methods for submitting new computations, querying their status, and
/// controlling their execution (pause/resume/cancel).
///
/// # Thread Safety
///
/// `ComputeEngine` is thread-safe and can be shared across multiple threads.
/// All internal state is protected by appropriate synchronization primitives.
///
/// # Caching
///
/// The engine maintains two caches:
/// - **Parsing cache**: Stores parsed expressions to avoid re-parsing
/// - **Result cache**: Stores computation results for reuse
#[allow(dead_code)]
#[derive(Clone)]

pub struct ComputeEngine {
    /// Registry of active computations, indexed by computation ID.
    computations: Arc<
        RwLock<
            HashMap<
                String,
                Arc<Mutex<Computation>>,
            >,
        >,
    >,
    /// Cache for parsed expressions.
    parsing_cache: Arc<ParsingCache>,
    /// Cache for computation results.
    result_cache:
        Arc<ComputationResultCache>,
}

impl ComputeEngine {
    /// Creates a new `ComputeEngine`.
    ///
    /// # Examples
    ///
    /// ```
    /// 
    /// use rssn::compute::engine::ComputeEngine;
    ///
    /// let engine = ComputeEngine::new();
    /// ```
    #[must_use]

    pub fn new() -> Self {

        Self {
            computations: Arc::new(RwLock::new(
                HashMap::new(),
            )),
            parsing_cache: Arc::new(ParsingCache::new()),
            result_cache: Arc::new(
                ComputationResultCache::new(),
            ),
        }
    }

    /// Parses an input string and submits it as a computation.
    ///
    /// This method first attempts to retrieve the parsed expression from the
    /// parsing cache. If not found, it parses the input and caches the result.
    /// The parsed expression is then submitted for computation.
    ///
    /// # Arguments
    ///
    /// * `input` - The input string to parse and compute
    ///
    /// # Returns
    ///
    /// * `Ok(String)` - The computation ID on success
    /// * `Err(String)` - An error message if parsing fails
    ///
    /// # Examples
    ///
    /// ```
    /// 
    /// use rssn::compute::engine::ComputeEngine;
    ///
    /// let engine = ComputeEngine::new();
    ///
    /// match engine.parse_and_submit("x + 1") {
    ///     | Ok(id) => {
    ///
    ///         println!(
    ///             "Computation ID: {}",
    ///             id
    ///         )
    ///     },
    ///     | Err(e) => eprintln!("Parse error: {}", e),
    /// }
    /// ```

    pub fn parse_and_submit(
        &self,
        input: &str,
    ) -> Result<String, String> {

        let expr = match self
            .parsing_cache
            .get(input)
        {
            | Some(expr) => expr,
            | None => {
                match crate::input::parser::parse_expr(
                    input,
                ) {
                    | Ok((_, expr)) => {

                        let expr = Arc::new(expr);

                        self.parsing_cache
                            .set(
                                input.to_string(),
                                expr.clone(),
                            );

                        expr
                    },
                    | Err(e) => return Err(e.to_string()),
                }
            },
        };

        Ok(self.submit(expr))
    }

    /// Gets the current status of a computation.
    ///
    /// # Arguments
    ///
    /// * `id` - The computation ID
    ///
    /// # Returns
    ///
    /// * `Some(ComputationStatus)` - The current status if the computation exists
    /// * `None` - If the computation ID is not found
    ///
    /// # Examples
    ///
    /// ```
    /// 
    /// use rssn::compute::engine::ComputeEngine;
    ///
    /// let engine = ComputeEngine::new();
    ///
    /// let id = engine
    ///     .parse_and_submit("2 + 2")
    ///     .unwrap();
    ///
    /// if let Some(status) = engine.get_status(&id) {
    ///
    ///     println!(
    ///         "Status: {:?}",
    ///         status
    ///     );
    /// }
    /// ```
    #[must_use]

    pub fn get_status(
        &self,
        id: &str,
    ) -> Option<ComputationStatus> {

        let computations = self
            .computations
            .read()
            .expect(
                "ComputeEngine \
                 computations lock \
                 poisoned",
            );

        computations
            .get(id)
            .map(|comp| {

                comp.lock()
                    .expect(
                        "Computation \
                         lock poisoned",
                    )
                    .status
                    .clone()
            })
    }

    /// Gets the current progress of a computation.
    ///
    /// # Arguments
    ///
    /// * `id` - The computation ID
    ///
    /// # Returns
    ///
    /// * `Some(ComputationProgress)` - The current progress if the computation exists
    /// * `None` - If the computation ID is not found
    ///
    /// # Examples
    ///
    /// ```
    /// 
    /// use rssn::compute::engine::ComputeEngine;
    ///
    /// let engine = ComputeEngine::new();
    ///
    /// let id = engine
    ///     .parse_and_submit("2 + 2")
    ///     .unwrap();
    ///
    /// if let Some(progress) = engine.get_progress(&id) {
    ///
    ///     println!(
    ///         "Progress: {}%",
    ///         progress.percentage
    ///     );
    /// }
    /// ```
    #[must_use]

    pub fn get_progress(
        &self,
        id: &str,
    ) -> Option<ComputationProgress>
    {

        let computations = self
            .computations
            .read()
            .expect(
                "ComputeEngine \
                 computations lock \
                 poisoned",
            );

        computations
            .get(id)
            .map(|comp| {

                comp.lock()
                    .expect(
                        "Computation \
                         lock poisoned",
                    )
                    .progress
                    .clone()
            })
    }

    /// Gets the result of a completed computation.
    ///
    /// # Arguments
    ///
    /// * `id` - The computation ID
    ///
    /// # Returns
    ///
    /// * `Some(Value)` - The result if the computation is complete
    /// * `None` - If the computation is not found or not yet complete
    ///
    /// # Examples
    ///
    /// ```
    /// 
    /// use rssn::compute::engine::ComputeEngine;
    ///
    /// let engine = ComputeEngine::new();
    ///
    /// let id = engine
    ///     .parse_and_submit("2 + 2")
    ///     .unwrap();
    ///
    /// // Wait for completion
    /// std::thread::sleep(std::time::Duration::from_secs(6));
    ///
    /// if let Some(result) = engine.get_result(&id) {
    ///
    ///     println!("Result: {}", result);
    /// }
    /// ```
    #[must_use]

    pub fn get_result(
        &self,
        id: &str,
    ) -> Option<Value> {

        let computations = self
            .computations
            .read()
            .expect(
                "ComputeEngine \
                 computations lock \
                 poisoned",
            );

        computations
            .get(id)
            .and_then(|comp| {

                comp.lock()
                    .expect(
                        "Computation \
                         lock poisoned",
                    )
                    .result
                    .clone()
            })
    }

    /// Submits an expression for asynchronous computation.
    ///
    /// This method creates a new computation task and executes it asynchronously
    /// using Rayon's thread pool. The computation can be monitored, paused,
    /// resumed, or cancelled using the returned ID.
    ///
    /// # Arguments
    ///
    /// * `expr` - The expression to compute
    ///
    /// # Returns
    ///
    /// A unique computation ID (UUID) as a string
    ///
    /// # Examples
    ///
    /// ```
    /// 
    /// use rssn::compute::engine::ComputeEngine;
    /// use rssn::symbolic::core::Expr;
    /// use std::sync::Arc;
    ///
    /// let engine = ComputeEngine::new();
    ///
    /// let expr = Arc::new(Expr::Constant(42.0));
    ///
    /// let id = engine.submit(expr);
    ///
    /// println!(
    ///     "Submitted computation: {}",
    ///     id
    /// );
    /// ```
    #[must_use]

    pub fn submit(
        &self,
        expr: Arc<Expr>,
    ) -> String {

        let id =
            Uuid::new_v4().to_string();

        let pause = Arc::new((
            Mutex::new(false),
            Condvar::new(),
        ));

        let computation = Arc::new(Mutex::new(
            Computation {
                id: id.clone(),
                expr,
                status: ComputationStatus::Pending,
                progress: ComputationProgress {
                    percentage: 0.0,
                    description: "Pending".to_string(),
                },
                result: None,
                cancel_signal: Arc::new(AtomicBool::new(
                    false,
                )),
                state: State {
                    intermediate_value: String::new(),
                },
                pause: pause.clone(),
            },
        ));

        {

            let mut computations = self
                .computations
                .write()
                .expect(
                    "ComputeEngine \
                     computations \
                     lock poisoned",
                );

            computations.insert(
                id.clone(),
                computation.clone(),
            );
        }

        let _engine = self.clone();

        let result_cache = self
            .result_cache
            .clone();

        rayon::spawn(move || {

            let (lock, cvar) = &*pause;

            let mut comp_guard =
                computation
                    .lock()
                    .expect(
                        "Computation \
                         lock poisoned",
                    );

            comp_guard.status = ComputationStatus::Running;

            // Simulate work
            for i in 0..100 {

                let mut paused =
                    lock.lock().expect(
                        "Pause lock \
                         poisoned",
                    );

                while *paused {

                    comp_guard.status =
                        ComputationStatus::Paused;

                    println!(
                        "Computation \
                         {} paused.",
                        comp_guard.id
                    );

                    paused = cvar
                        .wait(paused)
                        .expect(
                            "Condition variable wait \
                             failed",
                        );
                }

                comp_guard.status =
                    ComputationStatus::Running;

                if comp_guard.status
                    == ComputationStatus::Failed(
                        "Cancelled".to_string(),
                    )
                {

                    println!(
                        "Computation {} cancelled.",
                        comp_guard.id
                    );

                    return;
                }

                std::thread::sleep(
                    std::time::Duration::from_millis(50),
                );

                comp_guard
                    .progress
                    .percentage =
                    i as f32;

                comp_guard
                    .progress
                    .description = format!(
                    "{i}% complete"
                );
            }

            comp_guard.status =
                ComputationStatus::Completed;

            comp_guard
                .progress
                .percentage = 100.0;

            comp_guard
                .progress
                .description =
                "Completed".to_string();

            let result =
                "Result of the \
                 computation"
                    .to_string();

            comp_guard.result =
                Some(result.clone());

            result_cache.set(
                comp_guard
                    .expr
                    .clone(),
                result,
            );
        });

        id
    }

    /// Pauses a running computation.
    ///
    /// The computation will pause at the next checkpoint and can be resumed
    /// using the `resume` method.
    ///
    /// # Arguments
    ///
    /// * `id` - The computation ID
    ///
    /// # Examples
    ///
    /// ```
    /// 
    /// use rssn::compute::engine::ComputeEngine;
    ///
    /// let engine = ComputeEngine::new();
    ///
    /// let id = engine
    ///     .parse_and_submit("2 + 2")
    ///     .unwrap();
    ///
    /// // Pause the computation
    /// engine.pause(&id);
    /// ```

    pub fn pause(
        &self,
        id: &str,
    ) {

        if let Some(computation) = self
            .computations
            .read()
            .expect(
                "ComputeEngine \
                 computations lock \
                 poisoned",
            )
            .get(id)
        {

            let comp = computation
                .lock()
                .expect(
                    "Computation lock \
                     poisoned",
                );

            let (lock, cvar) =
                &*comp.pause;

            let mut paused =
                lock.lock().expect(
                    "Pause lock \
                     poisoned",
                );

            *paused = true;

            cvar.notify_one();
        }
    }

    /// Resumes a paused computation.
    ///
    /// # Arguments
    ///
    /// * `id` - The computation ID
    ///
    /// # Examples
    ///
    /// ```
    /// 
    /// use rssn::compute::engine::ComputeEngine;
    ///
    /// let engine = ComputeEngine::new();
    ///
    /// let id = engine
    ///     .parse_and_submit("2 + 2")
    ///     .unwrap();
    ///
    /// engine.pause(&id);
    ///
    /// // ... later ...
    /// engine.resume(&id);
    /// ```

    pub fn resume(
        &self,
        id: &str,
    ) {

        if let Some(computation) = self
            .computations
            .read()
            .expect(
                "ComputeEngine \
                 computations lock \
                 poisoned",
            )
            .get(id)
        {

            let comp = computation
                .lock()
                .expect(
                    "Computation lock \
                     poisoned",
                );

            let (lock, cvar) =
                &*comp.pause;

            let mut paused =
                lock.lock().expect(
                    "Pause lock \
                     poisoned",
                );

            *paused = false;

            cvar.notify_one();
        }
    }

    /// Cancels a computation and removes it from the registry.
    ///
    /// The computation will be marked as failed with status "Cancelled" and
    /// removed from the active computations.
    ///
    /// # Arguments
    ///
    /// * `id` - The computation ID
    ///
    /// # Examples
    ///
    /// ```
    /// 
    /// use rssn::compute::engine::ComputeEngine;
    ///
    /// let engine = ComputeEngine::new();
    ///
    /// let id = engine
    ///     .parse_and_submit("2 + 2")
    ///     .unwrap();
    ///
    /// // Cancel the computation
    /// engine.cancel(&id);
    /// ```

    pub fn cancel(
        &self,
        id: &str,
    ) {

        if let Some(computation) = self
            .computations
            .read()
            .expect(
                "ComputeEngine \
                 computations lock \
                 poisoned",
            )
            .get(id)
        {

            let mut comp = computation
                .lock()
                .expect(
                    "Computation lock \
                     poisoned",
                );

            comp.status = ComputationStatus::Failed(
                "Cancelled".to_string(),
            );

            let (lock, cvar) =
                &*comp.pause;

            let mut paused =
                lock.lock().expect(
                    "Pause lock \
                     poisoned",
                );

            *paused = false;

            cvar.notify_one();
        }

        self.computations
            .write()
            .expect(
                "ComputeEngine \
                 computations lock \
                 poisoned",
            )
            .remove(id);
    }
}

impl Default for ComputeEngine {
    fn default() -> Self {

        Self::new()
    }
}
