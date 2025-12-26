use crate::compute::computation::Value;
use crate::symbolic::core::Expr;
use std::collections::HashMap;
use std::sync::{Arc, Mutex};

/// A thread-safe cache for parsed expressions.
///
/// This cache stores the mapping from input strings to parsed `Expr` objects.
/// It avoids re-parsing the same string multiple times.

pub struct ParsingCache {
    cache: Mutex<HashMap<String, Arc<Expr>>>,
}

impl ParsingCache {
    /// Creates a new, empty `ParsingCache`.
    #[must_use]

    pub fn new() -> Self {

        Self {
            cache: Mutex::new(HashMap::new()),
        }
    }

    /// Retrieves a parsed expression from the cache, if available.
    ///
    /// # Arguments
    /// * `input` - The input string to look up.
    ///
    /// # Returns
    /// * `Option<Arc<Expr>>` - The cached expression, or `None` if not found.

    pub fn get(&self, input: &str) -> Option<Arc<Expr>> {

        let cache = self
            .cache
            .lock()
            .expect("ParsingCache lock poisoned");

        cache
            .get(input)
            .cloned()
    }

    /// Stores a parsed expression in the cache.
    ///
    /// # Arguments
    /// * `input` - The input string.
    /// * `expr` - The parsed expression.

    pub fn set(&self, input: String, expr: Arc<Expr>) {

        let mut cache = self
            .cache
            .lock()
            .expect("ParsingCache lock poisoned");

        cache.insert(input, expr);
    }

    /// Clears the cache.

    pub fn clear(&self) {

        let mut cache = self
            .cache
            .lock()
            .expect("ParsingCache lock poisoned");

        cache.clear();
    }
}

impl Default for ParsingCache {
    fn default() -> Self {

        Self::new()
    }
}

/// A thread-safe cache for computation results.
///
/// This cache stores the mapping from expressions to their computed values.
/// It avoids re-computing the value of the same expression multiple times.

pub struct ComputationResultCache {
    cache: Mutex<HashMap<Arc<Expr>, Value>>,
}

impl ComputationResultCache {
    /// Creates a new, empty `ComputationResultCache`.
    #[must_use]

    pub fn new() -> Self {

        Self {
            cache: Mutex::new(HashMap::new()),
        }
    }

    /// Retrieves a computed value from the cache, if available.
    ///
    /// # Arguments
    /// * `expr` - The expression to look up.
    ///
    /// # Returns
    /// * `Option<Value>` - The cached value, or `None` if not found.

    pub fn get(&self, expr: &Arc<Expr>) -> Option<Value> {

        let cache = self
            .cache
            .lock()
            .expect("ComputationResultCache lock poisoned");

        cache
            .get(expr)
            .cloned()
    }

    /// Stores a computed value in the cache.
    ///
    /// # Arguments
    /// * `expr` - The expression.
    /// * `value` - The computed value.

    pub fn set(&self, expr: Arc<Expr>, value: Value) {

        let mut cache = self
            .cache
            .lock()
            .expect("ComputationResultCache lock poisoned");

        cache.insert(expr, value);
    }

    /// Clears the cache.

    pub fn clear(&self) {

        let mut cache = self
            .cache
            .lock()
            .expect("ComputationResultCache lock poisoned");

        cache.clear();
    }
}

impl Default for ComputationResultCache {
    fn default() -> Self {

        Self::new()
    }
}
