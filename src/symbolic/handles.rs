//! # Thread-Safe Handle Manager for Expressions
//!
//! Provides a global, concurrent manager for `Expr` handles used across the FFI boundary.
//! This module enables safe sharing of symbolic expressions between Rust and foreign code
//! (C, Python, etc.) through opaque integer handles.
//!
//! ## Features
//!
//! - **Thread-safe**: Uses `DashMap` for lock-free concurrent access
//! - **Unique handles**: Atomic counter ensures no handle collisions
//! - **Memory safe**: Automatic cleanup when handles are freed
//! - **Reference counted**: Expressions are stored in `Arc` for efficient cloning
//!
//! ## Example
//!
//! ```rust
//! 
//! use rssn::symbolic::core::Expr;
//! use rssn::symbolic::handles::HANDLE_MANAGER;
//!
//! // Insert an expression and get a handle
//! let expr = Expr::new_variable("x");
//!
//! let handle = HANDLE_MANAGER.insert(expr);
//!
//! // Retrieve the expression
//! if let Some(expr_arc) = HANDLE_MANAGER.get(handle) {
//!
//!     println!(
//!         "Expression: {}",
//!         expr_arc
//!     );
//! }
//!
//! // Free the handle when done
//! HANDLE_MANAGER.free(handle);
//! ```

use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;
use std::sync::Arc;

use dashmap::DashMap;

use crate::symbolic::core::Expr;

/// A thread-safe, global manager for expression handles.
///
/// The `HandleManager` provides a centralized way to manage symbolic expressions
/// across FFI boundaries using integer handles. It uses a `DashMap` for concurrent
/// access and an `AtomicUsize` to generate unique handles.
///
/// Handles start at 1 (0 is reserved for null/invalid handles in FFI contexts).
#[derive(Debug)]

pub struct HandleManager {
    expressions :
        DashMap<usize, Arc<Expr>>,
    next_handle : AtomicUsize,
}

impl HandleManager {
    /// Creates a new, empty `HandleManager`.
    ///
    /// This is typically not called directly; use the global `HANDLE_MANAGER` instead.

    pub(crate) fn new() -> Self {

        Self {
            expressions : DashMap::new(
            ),
            next_handle:
                AtomicUsize::new(1), /* Start at 1, reserve 0 for null */
        }
    }

    /// Inserts a new expression into the manager and returns a unique handle for it.
    ///
    /// # Arguments
    ///
    /// * `expr` - The expression to store
    ///
    /// # Returns
    ///
    /// A unique handle (usize) that can be used to retrieve or free the expression later.
    ///
    /// # Example
    ///
    /// ```rust
    /// 
    /// use rssn::symbolic::core::Expr;
    /// use rssn::symbolic::handles::HANDLE_MANAGER;
    ///
    /// let expr = Expr::new_constant(42.0);
    ///
    /// let handle = HANDLE_MANAGER.insert(expr);
    ///
    /// assert!(handle > 0);
    /// ```

    pub fn insert(
        &self,
        expr : Expr,
    ) -> usize {

        let handle = self
            .next_handle
            .fetch_add(
                1,
                Ordering::SeqCst,
            );

        self.expressions
            .insert(
                handle,
                Arc::new(expr),
            );

        handle
    }

    /// Retrieves a clone of the `Arc<Expr>` associated with a given handle.
    ///
    /// Returns `None` if the handle is not found. The returned `Arc` is a cheap clone
    /// that shares ownership of the underlying expression.
    ///
    /// # Arguments
    ///
    /// * `handle` - The handle to look up
    ///
    /// # Returns
    ///
    /// `Some(Arc<Expr>)` if the handle exists, `None` otherwise.
    ///
    /// # Example
    ///
    /// ```rust
    /// 
    /// use rssn::symbolic::core::Expr;
    /// use rssn::symbolic::handles::HANDLE_MANAGER;
    ///
    /// let expr = Expr::new_variable("x");
    ///
    /// let handle = HANDLE_MANAGER.insert(expr);
    ///
    /// if let Some(retrieved) = HANDLE_MANAGER.get(handle) {
    ///
    ///     println!(
    ///         "Found: {}",
    ///         retrieved
    ///     );
    /// }
    /// ```

    pub fn get(
        &self,
        handle : usize,
    ) -> Option<Arc<Expr>> {

        self.expressions
            .get(&handle)
            .map(|arc_expr| {
                arc_expr.clone()
            })
    }

    /// Retrieves a deep clone of the expression (not the Arc).
    ///
    /// This is useful when you need an owned copy of the expression rather than
    /// a shared reference.
    ///
    /// # Arguments
    ///
    /// * `handle` - The handle to look up
    ///
    /// # Returns
    ///
    /// `Some(Expr)` if the handle exists, `None` otherwise.

    pub fn clone_expr(
        &self,
        handle : usize,
    ) -> Option<Expr> {

        self.expressions
            .get(&handle)
            .map(|arc_expr| {
                (**arc_expr).clone()
            })
    }

    /// Removes an expression from the manager, freeing its memory if this was the last reference.
    ///
    /// Returns the `Arc<Expr>` if it was found, otherwise `None`.
    ///
    /// # Arguments
    ///
    /// * `handle` - The handle to free
    ///
    /// # Returns
    ///
    /// `Some(Arc<Expr>)` if the handle existed, `None` otherwise.
    ///
    /// # Example
    ///
    /// ```rust
    /// 
    /// use rssn::symbolic::core::Expr;
    /// use rssn::symbolic::handles::HANDLE_MANAGER;
    ///
    /// let expr = Expr::new_constant(3.14);
    ///
    /// let handle = HANDLE_MANAGER.insert(expr);
    ///
    /// // Free the handle
    /// let freed = HANDLE_MANAGER.free(handle);
    ///
    /// assert!(freed.is_some());
    ///
    /// // Trying to get it again returns None
    /// assert!(HANDLE_MANAGER
    ///     .get(handle)
    ///     .is_none());
    /// ```

    pub fn free(
        &self,
        handle : usize,
    ) -> Option<Arc<Expr>> {

        self.expressions
            .remove(&handle)
            .map(|(_, arc_expr)| {
                arc_expr
            })
    }

    /// Checks if a handle exists in the manager.
    ///
    /// # Arguments
    ///
    /// * `handle` - The handle to check
    ///
    /// # Returns
    ///
    /// `true` if the handle exists, `false` otherwise.

    pub fn exists(
        &self,
        handle : usize,
    ) -> bool {

        self.expressions
            .contains_key(&handle)
    }

    /// Returns the number of expressions currently managed.
    ///
    /// # Returns
    ///
    /// The count of active handles.
    ///
    /// # Example
    ///
    /// ```rust
    /// 
    /// use rssn::symbolic::core::Expr;
    /// use rssn::symbolic::handles::HANDLE_MANAGER;
    ///
    /// let initial_count = HANDLE_MANAGER.count();
    ///
    /// let handle = HANDLE_MANAGER.insert(Expr::new_constant(
    ///     1.0,
    /// ));
    ///
    /// assert_eq!(
    ///     HANDLE_MANAGER.count(),
    ///     initial_count + 1
    /// );
    ///
    /// HANDLE_MANAGER.free(handle);
    ///
    /// assert_eq!(
    ///     HANDLE_MANAGER.count(),
    ///     initial_count
    /// );
    /// ```

    pub fn count(&self) -> usize {

        self.expressions
            .len()
    }

    /// Clears all expressions from the manager.
    ///
    /// **Warning**: This will invalidate all existing handles. Use with caution.
    ///
    /// # Example
    ///
    /// ```rust
    /// 
    /// use rssn::symbolic::core::Expr;
    /// use rssn::symbolic::handles::HANDLE_MANAGER;
    ///
    /// let h1 = HANDLE_MANAGER.insert(Expr::new_constant(
    ///     1.0,
    /// ));
    ///
    /// let h2 = HANDLE_MANAGER.insert(Expr::new_constant(
    ///     2.0,
    /// ));
    ///
    /// HANDLE_MANAGER.clear();
    ///
    /// assert_eq!(
    ///     HANDLE_MANAGER.count(),
    ///     0
    /// );
    ///
    /// assert!(HANDLE_MANAGER
    ///     .get(h1)
    ///     .is_none());
    ///
    /// assert!(HANDLE_MANAGER
    ///     .get(h2)
    ///     .is_none());
    /// ```

    pub fn clear(&self) {

        self.expressions
            .clear();
    }

    /// Returns a vector of all active handles.
    ///
    /// The handles are returned in arbitrary order.
    ///
    /// # Returns
    ///
    /// A `Vec<usize>` containing all active handles.
    ///
    /// # Example
    ///
    /// ```rust
    /// 
    /// use rssn::symbolic::core::Expr;
    /// use rssn::symbolic::handles::HANDLE_MANAGER;
    ///
    /// HANDLE_MANAGER.clear(); // Start fresh
    /// let h1 = HANDLE_MANAGER.insert(Expr::new_constant(
    ///     1.0,
    /// ));
    ///
    /// let h2 = HANDLE_MANAGER.insert(Expr::new_constant(
    ///     2.0,
    /// ));
    ///
    /// let handles = HANDLE_MANAGER.get_all_handles();
    ///
    /// assert_eq!(handles.len(), 2);
    ///
    /// assert!(handles.contains(&h1));
    ///
    /// assert!(handles.contains(&h2));
    /// ```

    pub fn get_all_handles(
        &self
    ) -> Vec<usize> {

        self.expressions
            .iter()
            .map(|entry| *entry.key())
            .collect()
    }
}

impl Default for HandleManager {
    fn default() -> Self {

        Self::new()
    }
}

/// The global singleton instance of the `HandleManager`.
///
/// This is the primary way to interact with the handle manager. It's initialized
/// lazily on first use and is thread-safe.
///
/// # Example
///
/// ```rust
/// 
/// use rssn::symbolic::core::Expr;
/// use rssn::symbolic::handles::HANDLE_MANAGER;
///
/// // Use the global instance
/// let handle = HANDLE_MANAGER.insert(Expr::new_variable(
///     "x",
/// ));
///
/// let expr = HANDLE_MANAGER.get(handle);
///
/// HANDLE_MANAGER.free(handle);
/// ```

pub static HANDLE_MANAGER :
    std::sync::LazyLock<HandleManager> =
    std::sync::LazyLock::new(
        HandleManager::new,
    );
