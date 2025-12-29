//! # RSSN Plugin System
//! Defines the core traits and data structures for the `rssn` plugin architecture.
#![allow(unsafe_code)]
#![allow(clippy::indexing_slicing)]
#![allow(
    clippy::no_mangle_with_rust_abi
)]

use std::collections::HashMap;
use std::error::Error;
use std::fmt;

/// Represents the health status of a plugin, for use in heartbeat checks.
use abi_stable::std_types::RString;

use crate::symbolic::core::Expr;

#[repr(C)]
#[derive(
    abi_stable::StableAbi,
    Debug,
    Clone,
    PartialEq,
    Eq,
)]
/// Represents the health status of a plugin.

pub enum PluginHealth {
    /// The plugin is operating correctly.
    Ok,
    /// The plugin is in a degraded state but is still functional.
    Warning(RString),
    /// The plugin has encountered a critical error.
    Error(RString),
    /// The plugin needs to be reinitialized.
    RequiresReinitialization,
    /// The plugin has been terminated.
    Terminated,
    /// The plugin's health is unknown.
    Unknown,
}

/// Represents the kind of error that a plugin can encounter.
#[derive(
    Debug, Clone, PartialEq, Eq,
)]

pub enum PluginErrorKind {
    /// The plugin was not found.
    NotFound,
    /// The plugin's API version is incompatible with the host.
    IncompatibleVersion,
    /// The plugin failed to load.
    LoadFailed,
    /// The plugin failed to execute a command.
    ExecutionFailed,
    /// An error occurred during serialization or deserialization.
    SerializationError,
    /// An internal error occurred within the plugin.
    InternalError,
}

/// A specialized error type for plugin-related failures.
#[derive(Debug, Clone)]

pub struct PluginError {
    /// The kind of error.
    pub kind: PluginErrorKind,
    /// A descriptive error message.
    pub message: String,
}

impl PluginError {
    /// Creates a new `PluginError`.

    #[must_use]

    pub fn new(
        kind: PluginErrorKind,
        msg: &str,
    ) -> Self {

        Self {
            kind,
            message: msg.to_string(),
        }
    }
}

impl fmt::Display for PluginError {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {

        write!(
            f,
            "Plugin Error ({:?}): {}",
            self.kind, self.message
        )
    }
}

impl Error for PluginError {
}

/// The central trait that all `rssn` plugins must implement.
///
/// This trait defines the contract between the main library and a dynamically loaded plugin,
/// covering identity, lifecycle, execution, and health monitoring.

pub trait Plugin: Send + Sync {
    /// Returns the unique, machine-readable name of the plugin (e.g., "`fortran_solver`").

    fn name(&self) -> &'static str;

    /// Returns a human-readable description of the plugin.

    fn description(
        &self
    ) -> &'static str {

        ""
    }

    /// Returns the semantic version of the RSSN API the plugin was built against.
    /// The `PluginManager` will use this to ensure compatibility.
    /// Example: `"0.1.0"`

    fn api_version(
        &self
    ) -> &'static str;

    /// Called once when the plugin is loaded by the `PluginManager`.
    /// Use this for any necessary setup or initialization.
    ///
    /// # Errors
    ///
    /// This function will return a `PluginError` if the plugin fails to initialize
    /// or encounters an unrecoverable error during setup.

    fn on_load(
        &self
    ) -> Result<(), PluginError>;

    /// The primary entry point for executing plugin functionality.
    ///
    /// # Arguments
    /// * `command` - A string identifier for the specific function to execute within the plugin.
    /// * `args` - A symbolic `Expr` tree passed as an argument.
    ///
    /// # Returns
    /// A `Result` containing either a resulting `Expr` on success or a `PluginError` on failure.
    ///
    /// # Errors
    ///
    /// This function will return a `PluginError` if the specified `command` is not
    /// recognized or if an error occurs during the execution of the plugin's logic.

    fn execute(
        &self,
        command: &str,
        args: &Expr,
    ) -> Result<Expr, PluginError>;

    /// Performs a health check on the plugin.
    /// This is used by the `PluginManager` for heartbeat monitoring.

    fn health_check(
        &self
    ) -> PluginHealth;

    /// Returns a map of metadata for the plugin.

    fn metadata(
        &self
    ) -> HashMap<String, String> {

        HashMap::new()
    }
}
