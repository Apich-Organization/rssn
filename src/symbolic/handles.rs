// src/symbolic/handles.rs

//! # Thread-Safe Handle Manager for Expressions
//! Provides a global, concurrent manager for `Expr` handles used across the FFI boundary.

use crate::symbolic::core::Expr;
use dashmap::DashMap;
use once_cell::sync::Lazy;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

/// A thread-safe, global manager for expression handles.
///
/// It uses a `DashMap` for concurrent access and an `AtomicUsize` to generate unique handles.
pub struct HandleManager {
    expressions: DashMap<usize, Arc<Expr>>,
    next_handle: AtomicUsize,
}

impl HandleManager {
    /// Creates a new, empty `HandleManager`.
    pub(crate) fn new() -> Self {
        HandleManager {
            expressions: DashMap::new(),
            // Start handles from 1, as 0 can be used as a null/invalid handle.
            next_handle: AtomicUsize::new(1),
        }
    }

    /// Inserts a new expression into the manager and returns a unique handle for it.
    pub fn insert(&self, expr: Expr) -> usize {
        let handle = self.next_handle.fetch_add(1, Ordering::SeqCst);
        self.expressions.insert(handle, Arc::new(expr));
        handle
    }

    /// Retrieves a clone of the `Arc<Expr>` associated with a given handle.
    ///
    /// Returns `None` if the handle is not found.
    pub fn get(&self, handle: usize) -> Option<Arc<Expr>> {
        self.expressions
            .get(&handle)
            .map(|arc_expr| arc_expr.clone())
    }

    /// Removes an expression from the manager, freeing its memory if the handle was the last reference.
    ///
    /// Returns the `Arc<Expr>` if it was found, otherwise `None`.
    pub fn free(&self, handle: usize) -> Option<Arc<Expr>> {
        self.expressions
            .remove(&handle)
            .map(|(_, arc_expr)| arc_expr)
    }
}

/// The global singleton instance of the `HandleManager`.
pub static HANDLE_MANAGER: Lazy<HandleManager> = Lazy::new(HandleManager::new);
