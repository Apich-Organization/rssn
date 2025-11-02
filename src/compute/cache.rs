use crate::compute::computation::Value;
use crate::symbolic::core::Expr;
use std::collections::HashMap;
use std::sync::{Arc, Mutex};

pub struct ParsingCache {
    cache: Mutex<HashMap<String, Arc<Expr>>>,
}

impl ParsingCache {
    pub fn new() -> Self {
        Self {
            cache: Mutex::new(HashMap::new()),
        }
    }

    pub fn get(&self, input: &str) -> Option<Arc<Expr>> {
        let cache = self.cache.lock().expect("ParsingCache lock poisoned");
        cache.get(input).cloned()
    }

    pub fn set(&self, input: String, expr: Arc<Expr>) {
        let mut cache = self.cache.lock().expect("ParsingCache lock poisoned");
        cache.insert(input, expr);
    }
}

impl Default for ParsingCache {
    fn default() -> Self {
        Self::new()
    }
}

pub struct ComputationResultCache {
    cache: Mutex<HashMap<Arc<Expr>, Value>>,
}

impl ComputationResultCache {
    pub fn new() -> Self {
        Self {
            cache: Mutex::new(HashMap::new()),
        }
    }

    pub fn get(&self, expr: &Arc<Expr>) -> Option<Value> {
        let cache = self
            .cache
            .lock()
            .expect("ComputationResultCache lock poisoned");
        cache.get(expr).cloned()
    }

    pub fn set(&self, expr: Arc<Expr>, value: Value) {
        let mut cache = self
            .cache
            .lock()
            .expect("ComputationResultCache lock poisoned");
        cache.insert(expr, value);
    }
}

impl Default for ComputationResultCache {
    fn default() -> Self {
        Self::new()
    }
}
