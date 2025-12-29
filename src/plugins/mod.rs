//! # FFI Plugin System
//!
//! The `plugins` module and its submodules are for testing the FFI plugin system.
//!
//! # Usage
//!
//! ```rust
//! use rssn::plugins::manager::PluginManager;
/// Plugin manager for loading and managing plugins.
pub mod manager;
/// C-style plugin interface.
pub mod plugin_c;
/// Stable ABI plugin interface.
pub mod stable_abi;
