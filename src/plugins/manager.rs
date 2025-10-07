// src/plugins/manager.rs

//! # Plugin Manager
//! Responsible for loading, managing, and interacting with plugins.

use crate::plugins::{Plugin, PluginHealth};
use libloading::{Library, Symbol};
use std::collections::HashMap;
use std::error::Error;
use std::ffi::c_void;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, RwLock};
use std::thread;
use std::time::Duration;

/// The expected signature for the function that instantiates the plugin.
///
/// Each plugin library must export a C-compatible function named `_plugin_create`
/// with this signature, returning a pointer to a Box containing the plugin trait object.
type PluginCreate = unsafe extern "C" fn() -> *mut Box<dyn Plugin>;

/// Holds the plugin instance and its current health state.
pub struct ManagedPlugin {
    pub plugin: Box<dyn Plugin>,
    pub health: RwLock<PluginHealth>,
}

/// Manages the lifecycle of all loaded plugins.
pub struct PluginManager {
    plugins: Arc<RwLock<HashMap<String, ManagedPlugin>>>,
    health_check_thread: Option<thread::JoinHandle<()>>,
    stop_signal: Arc<AtomicBool>,
    _libraries: Vec<Library>,
}

impl PluginManager {
    /// Creates a new `PluginManager` and loads plugins from a specified directory.
    pub fn new(plugin_dir: &str) -> Result<Self, Box<dyn Error>> {
        let mut manager = PluginManager {
            plugins: Arc::new(RwLock::new(HashMap::new())),
            health_check_thread: None,
            stop_signal: Arc::new(AtomicBool::new(false)),
            _libraries: Vec::new(),
        };

        manager.load_plugins(plugin_dir)?;
        manager.start_health_checks(Duration::from_secs(30));

        Ok(manager)
    }

    /// Scans a directory for dynamic libraries and attempts to load them as plugins.
    pub(crate) fn load_plugins(&mut self, directory: &str) -> Result<(), Box<dyn Error>> {
        for entry in std::fs::read_dir(directory)? {
            let entry = entry?;
            let path = entry.path();

            if path.is_file()
                && path
                    .extension()
                    .map_or(false, |ext| ext == std::env::consts::DLL_EXTENSION)
            {
                println!("Attempting to load plugin: {:?}", path);
                unsafe {
                    self.load_plugin(&path)?;
                }
            }
        }
        Ok(())
    }

    /// Loads a single plugin from a dynamic library file.
    unsafe fn load_plugin(&mut self, library_path: &std::path::Path) -> Result<(), Box<dyn Error>> {
        let library = Library::new(library_path)?;
        self._libraries.push(library);
        let library = self._libraries.last().unwrap();

        let constructor: Symbol<'_, PluginCreate> = library.get(b"_plugin_create")?;
        let plugin_box_ptr = constructor();
        // Reconstitute the outer box, which gives us ownership of the inner Box<dyn Plugin>.
        let plugin_box = Box::from_raw(plugin_box_ptr);
        // Unbox it to get the `Box<dyn Plugin>`.
        let plugin = *plugin_box;

        let plugin_api_version = plugin.api_version();
        let crate_version = env!("CARGO_PKG_VERSION");

        let (p_major, p_minor) = parse_version(plugin_api_version)?;
        let (c_major, c_minor) = parse_version(crate_version)?;

        if p_major != c_major || p_minor != c_minor {
            return Err(format!(
                "Plugin '{}' has incompatible API version {}. Expected a version compatible with {}.",
                plugin.name(), plugin_api_version, crate_version
            ).into());
        }

        println!(
            "Successfully loaded plugin: {} (API version: {})",
            plugin.name(),
            plugin_api_version
        );

        plugin.on_load()?;

        let managed_plugin = ManagedPlugin {
            plugin,
            health: RwLock::new(PluginHealth::Ok),
        };

        self.plugins
            .write()
            .unwrap()
            .insert(managed_plugin.plugin.name().to_string(), managed_plugin);

        Ok(())
    }

    /// Spawns a background thread to periodically perform health checks on all plugins.
    pub(crate) fn start_health_checks(&mut self, interval: Duration) {
        let plugins_clone = Arc::clone(&self.plugins);
        let stop_signal_clone = Arc::clone(&self.stop_signal);

        let handle = thread::spawn(move || {
            while !stop_signal_clone.load(Ordering::SeqCst) {
                println!("Running plugin health checks...");
                let plugins_map = plugins_clone.read().unwrap();
                for (name, managed_plugin) in plugins_map.iter() {
                    let new_health = managed_plugin.plugin.health_check();
                    let mut health_writer = managed_plugin.health.write().unwrap();
                    *health_writer = new_health;
                }
                drop(plugins_map); // Release read lock before sleeping
                thread::sleep(interval);
            }
            println!("Plugin health check thread shutting down.");
        });

        self.health_check_thread = Some(handle);
    }
}

impl Drop for PluginManager {
    fn drop(&mut self) {
        println!("Shutting down plugin manager...");
        self.stop_signal.store(true, Ordering::SeqCst);
        if let Some(handle) = self.health_check_thread.take() {
            handle.join().expect("Failed to join health check thread");
        }
        println!("Plugin manager shut down.");
    }
}

/// A helper to parse a semantic version string into (major, minor) components.
pub(crate) fn parse_version(version: &str) -> Result<(u32, u32), Box<dyn Error>> {
    let mut parts = version.split('.');
    let major = parts
        .next()
        .ok_or("Invalid version string")?
        .parse::<u32>()?;
    let minor = parts
        .next()
        .ok_or("Invalid version string")?
        .parse::<u32>()?;
    Ok((major, minor))
}
