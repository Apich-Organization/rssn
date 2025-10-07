// src/plugins/mod.rs

//! # RSSN Plugin System
//! Defines the core traits and data structures for the `rssn` plugin architecture.

pub mod manager;

use crate::symbolic::core::Expr;
use std::error::Error;
use std::fmt;

/// Represents the health status of a plugin, for use in heartbeat checks.
#[derive(Debug)]
pub enum PluginHealth {
    /// The plugin is operating correctly.
    Ok,
    /// The plugin is in a degraded state but is still functional.
    Degraded(String),
    /// The plugin has encountered a critical error.
    Error(String),
}

/// A specialized error type for plugin-related failures.
#[derive(Debug)]
pub struct PluginError {
    message: String,
}

impl PluginError {
    pub fn new(msg: &str) -> Self {
        PluginError {
            message: msg.to_string(),
        }
    }
}

impl fmt::Display for PluginError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Plugin Error: {}", self.message)
    }
}

impl Error for PluginError {}

/// The central trait that all `rssn` plugins must implement.
///
/// This trait defines the contract between the main library and a dynamically loaded plugin,
/// covering identity, lifecycle, execution, and health monitoring.
pub trait Plugin: Send + Sync {
    /// Returns the unique, machine-readable name of the plugin (e.g., "fortran_solver").
    fn name(&self) -> &'static str;

    /// Returns the semantic version of the RSSN API the plugin was built against.
    /// The PluginManager will use this to ensure compatibility.
    /// Example: `"0.1.0"`
    fn api_version(&self) -> &'static str;

    /// Called once when the plugin is loaded by the `PluginManager`.
    /// Use this for any necessary setup or initialization.
    fn on_load(&self) -> Result<(), PluginError>;

    /// The primary entry point for executing plugin functionality.
    ///
    /// # Arguments
    /// * `command` - A string identifier for the specific function to execute within the plugin.
    /// * `args` - A symbolic `Expr` tree passed as an argument.
    ///
    /// # Returns
    /// A `Result` containing either a resulting `Expr` on success or a `PluginError` on failure.
    fn execute(&self, command: &str, args: &Expr) -> Result<Expr, PluginError>;

    /// Performs a health check on the plugin.
    /// This is used by the `PluginManager` for heartbeat monitoring.
    fn health_check(&self) -> PluginHealth;
}
