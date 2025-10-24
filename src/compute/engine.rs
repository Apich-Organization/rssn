#![allow(unused_imports)]
use crate::compute::cache::{ComputationResultCache, ParsingCache};
use crate::compute::computation::{Computation, ComputationProgress, ComputationStatus, Value};
use crate::compute::state::State;
use crate::symbolic::core::Expr;
/// Development in place.
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::prelude::*;
use std::sync::atomic::AtomicBool;
use std::sync::Condvar;
use std::sync::{Arc, Mutex, RwLock};
use uuid::Uuid;

#[allow(dead_code)]
#[derive(Clone)]
pub struct ComputeEngine {
    computations: Arc<RwLock<HashMap<String, Arc<Mutex<Computation>>>>>,
    parsing_cache: Arc<ParsingCache>,
    result_cache: Arc<ComputationResultCache>,
}

impl ComputeEngine {
    pub fn new() -> Self {
        Self {
            computations: Arc::new(RwLock::new(HashMap::new())),
            parsing_cache: Arc::new(ParsingCache::new()),
            result_cache: Arc::new(ComputationResultCache::new()),
        }
    }

    pub fn parse_and_submit(&self, input: &str) -> Result<String, String> {
        let expr = match self.parsing_cache.get(input) {
            Some(expr) => expr,
            None => match crate::input::parser::parse_expr(input) {
                Ok((_, expr)) => {
                    let expr = Arc::new(expr);
                    self.parsing_cache.set(input.to_string(), expr.clone());
                    expr
                }
                Err(e) => return Err(e.to_string()),
            },
        };
        Ok(self.submit(expr))
    }

    pub fn get_status(&self, id: &str) -> Option<ComputationStatus> {
        let computations = self.computations.read().unwrap();
        computations
            .get(id)
            .map(|comp| comp.lock().unwrap().status.clone())
    }

    pub fn get_progress(&self, id: &str) -> Option<ComputationProgress> {
        let computations = self.computations.read().unwrap();
        computations
            .get(id)
            .map(|comp| comp.lock().unwrap().progress.clone())
    }

    pub fn get_result(&self, id: &str) -> Option<Value> {
        let computations = self.computations.read().unwrap();
        computations
            .get(id)
            .and_then(|comp| comp.lock().unwrap().result.clone())
    }

    pub fn submit(&self, expr: Arc<Expr>) -> String {
        let id = Uuid::new_v4().to_string();
        let pause = Arc::new((Mutex::new(false), Condvar::new()));
        let computation = Arc::new(Mutex::new(Computation {
            id: id.clone(),
            expr,
            status: ComputationStatus::Pending,
            progress: ComputationProgress {
                percentage: 0.0,
                description: "Pending".to_string(),
            },
            result: None,
            cancel_signal: Arc::new(AtomicBool::new(false)),
            state: State {
                intermediate_value: "".to_string(),
            },
            pause: pause.clone(),
        }));

        {
            let mut computations = self.computations.write().unwrap();
            computations.insert(id.clone(), computation.clone());
        }

        let _engine = self.clone();
        let result_cache = self.result_cache.clone();
        rayon::spawn(move || {
            let (lock, cvar) = &*pause;
            let mut comp_guard = computation.lock().unwrap();
            comp_guard.status = ComputationStatus::Running;

            // Simulate work
            for i in 0..100 {
                let mut paused = lock.lock().unwrap();
                while *paused {
                    comp_guard.status = ComputationStatus::Paused;
                    println!("Computation {} paused.", comp_guard.id);
                    paused = cvar.wait(paused).unwrap();
                }
                comp_guard.status = ComputationStatus::Running;

                if comp_guard.status == ComputationStatus::Failed("Cancelled".to_string()) {
                    println!("Computation {} cancelled.", comp_guard.id);
                    return;
                }

                std::thread::sleep(std::time::Duration::from_millis(50));
                comp_guard.progress.percentage = i as f32;
                comp_guard.progress.description = format!("{}% complete", i);
            }

            comp_guard.status = ComputationStatus::Completed;
            comp_guard.progress.percentage = 100.0;
            comp_guard.progress.description = "Completed".to_string();
            let result = "Result of the computation".to_string();
            comp_guard.result = Some(result.clone());
            result_cache.set(comp_guard.expr.clone(), result);
        });

        id
    }

    pub fn pause(&self, id: &str) {
        if let Some(computation) = self.computations.read().unwrap().get(id) {
            let comp = computation.lock().unwrap();
            let (lock, cvar) = &*comp.pause;
            let mut paused = lock.lock().unwrap();
            *paused = true;
            cvar.notify_one();
        }
    }

    pub fn resume(&self, id: &str) {
        if let Some(computation) = self.computations.read().unwrap().get(id) {
            let comp = computation.lock().unwrap();
            let (lock, cvar) = &*comp.pause;
            let mut paused = lock.lock().unwrap();
            *paused = false;
            cvar.notify_one();
        }
    }

    pub fn cancel(&self, id: &str) {
        if let Some(computation) = self.computations.read().unwrap().get(id) {
            let mut comp = computation.lock().unwrap();
            comp.status = ComputationStatus::Failed("Cancelled".to_string());
            let (lock, cvar) = &*comp.pause;
            let mut paused = lock.lock().unwrap();
            *paused = false;
            cvar.notify_one();
        }
        self.computations.write().unwrap().remove(id);
    }
}
