//! # Plugin Manager
//! Responsible for loading, managing, and interacting with plugins.
#![allow(unsafe_code)]
#![allow(clippy::indexing_slicing)]
#![allow(clippy::no_mangle_with_rust_abi)]
use crate::plugins::plugin_c::{Plugin, PluginError, PluginHealth};
use crate::symbolic::core::Expr;
use libloading::{Library, Symbol};
use std::collections::HashMap;
use std::error::Error;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, RwLock};
use std::thread;
use std::time::Duration;
use bincode::config;
use bincode::serde;
/// The expected signature for the function that instantiates the plugin.
///
/// Each plugin library must export a C-compatible function named `_plugin_create`
/// with this signature, returning a pointer to a Box containing the plugin trait object.
type PluginCreate = unsafe extern "C" fn() -> *mut Box<dyn Plugin>;
/// Holds the plugin instance and its current health state.
use crate::plugins::stable_abi::{StablePlugin_TO, StablePluginModule};
use abi_stable::std_types::RBox;

pub struct ManagedPlugin {
    pub plugin: Box<dyn Plugin>,
    pub health: RwLock<PluginHealth>,
}

pub struct ManagedStablePlugin {
    pub plugin: StablePlugin_TO<'static,RBox<()>>,
    pub health: RwLock<PluginHealth>,
}

/// Manages the lifecycle of all loaded plugins.
pub struct PluginManager {
    pub plugins: Arc<RwLock<HashMap<String, ManagedPlugin>>>,
    pub stable_plugins: Arc<RwLock<HashMap<String, ManagedStablePlugin>>>,
    pub health_check_thread: Option<thread::JoinHandle<()>>,
    pub stop_signal: Arc<AtomicBool>,
    pub libraries: Vec<Library>,
}
impl PluginManager {
    /// Creates a new `PluginManager` and loads plugins from a specified directory.
    pub fn new(plugin_dir: &str) -> Result<Self, Box<dyn Error>> {
        let mut manager = PluginManager {
            plugins: Arc::new(RwLock::new(HashMap::new())),
            stable_plugins: Arc::new(RwLock::new(HashMap::new())),
            health_check_thread: None,
            stop_signal: Arc::new(AtomicBool::new(false)),
            libraries: Vec::new(),
        };
        manager.load_plugins(plugin_dir)?;
        manager.start_health_checks(Duration::from_secs(30));
        Ok(manager)
    }
    /// Executes a command on a specific plugin.
    ///
    /// # Arguments
    /// * `plugin_name` - The name of the target plugin.
    /// * `command` - The command to execute on the plugin.
    /// * `args` - The symbolic expression to pass as an argument.
    ///
    /// # Returns
    /// A `Result` containing the output `Expr` or a `PluginError`.
    pub fn execute_plugin(
        &self,
        plugin_name: &str,
        command: &str,
        args: &Expr,
    ) -> Result<Expr, PluginError> {
        let stable_plugins_map = self.stable_plugins.read().expect("No data was found.");
        if let Some(managed_plugin) = stable_plugins_map.get(plugin_name) {
		let config = config::standard(); // Get the configuration object
		let args_vec = serde::encode_to_vec(args, config).map_err(|e| PluginError::new(&e.to_string()))?;
		let result_vec = managed_plugin.plugin.execute(command.into(), args_vec.into()).into_result().map_err(|e| PluginError::new(&e.to_string()))?;
		let config = config::standard(); // Get the configuration object
		let (result_expr, _) = serde::decode_from_slice(&result_vec, config).map_err(|e| PluginError::new(&e.to_string()))?;
		return Ok(result_expr);
        }

        let plugins_map = self.plugins.read().expect("No data was found.");
        if let Some(managed_plugin) = plugins_map.get(plugin_name) {
            return managed_plugin.plugin.execute(command, args);
        }

        Err(PluginError::new(&format!(
            "Plugin '{}' not found.",
            plugin_name
        )))
    }
    /// Scans a directory for dynamic libraries and attempts to load them as plugins.
    pub(crate) fn load_plugins(&mut self, directory: &str) -> Result<(), Box<dyn Error>> {
        for entry in std::fs::read_dir(directory)? {
            let entry = entry?;
            let path = entry.path();
            if path.is_file()
                && path
                    .extension()
                    .is_some_and(|ext| ext == std::env::consts::DLL_EXTENSION)
            {
                println!("Attempting to load plugin: {:?}", path.display());
                unsafe {
                    if self.load_stable_plugin(&path).is_err() {
                        self.load_plugin(&path)?;
                    }
                }
            }
        }
        Ok(())
    }

    unsafe fn load_stable_plugin(&mut self, library_path: &std::path::Path) -> Result<(), Box<dyn Error>> {
        let library = Library::new(library_path)?;
        self.libraries.push(library);
        let library = self.libraries.last().expect("Library invalid");

        let module: Symbol<'_, *const StablePluginModule> = library.get(b"STABLE_PLUGIN_MODULE")?;
        let module = &**module;

        let plugin = (module.new)();
        let plugin_name = (module.name)();
        let plugin_api_version = (module.version)();

        let crate_version = env!("CARGO_PKG_VERSION");
        let (p_major, p_minor) = parse_version(&plugin_api_version)?;
        let (c_major, c_minor) = parse_version(crate_version)?;

        if p_major != c_major || p_minor != c_minor {
            return Err(
                format!(
                    "Plugin '{}' has incompatible API version {}. Expected a version compatible with {}.",
                    plugin_name, plugin_api_version, crate_version
                )
                    .into(),
            );
        }

        println!(
            "Successfully loaded stable plugin: {} (API version: {})",
            plugin_name,
            plugin_api_version
        );

        plugin.on_load().into_result().map_err(|e| e.to_string())?;

        let managed_plugin = ManagedStablePlugin {
            plugin,
            health: RwLock::new(PluginHealth::Ok),
        };

        self.stable_plugins
            .write()
            .expect("Unexpected Error")
            .insert(plugin_name.to_string(), managed_plugin);

        Ok(())
    }
    /// Loads a single plugin from a dynamic library file.
    unsafe fn load_plugin(&mut self, library_path: &std::path::Path) -> Result<(), Box<dyn Error>> {
        let library = Library::new(library_path)?;
        self.libraries.push(library);
        let library = self.libraries.last().expect("Library invalid");
        let constructor: Symbol<'_, PluginCreate> = library.get(b"_plugin_create")?;
        let plugin_box_ptr = constructor();
        let plugin_box = Box::from_raw(plugin_box_ptr);
        let plugin = *plugin_box;
        let plugin_api_version = plugin.api_version();
        let crate_version = env!("CARGO_PKG_VERSION");
        let (p_major, p_minor) = parse_version(plugin_api_version)?;
        let (c_major, c_minor) = parse_version(crate_version)?;
        if p_major != c_major || p_minor != c_minor {
            return Err(
                format!(
                    "Plugin '{}' has incompatible API version {}. Expected a version compatible with {}.",
                    plugin.name(), plugin_api_version, crate_version
                )
                    .into(),
            );
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
            .expect("Unexpected Error")
            .insert(managed_plugin.plugin.name().to_string(), managed_plugin);
        Ok(())
    }
    /// Spawns a background thread to periodically perform health checks on all plugins.
    pub(crate) fn start_health_checks(&mut self, interval: Duration) {
        let plugins_clone = Arc::clone(&self.plugins);
        let stable_plugins_clone = Arc::clone(&self.stable_plugins);
        let stop_signal_clone = Arc::clone(&self.stop_signal);
        let handle = thread::spawn(move || {
            while !stop_signal_clone.load(Ordering::SeqCst) {
                println!("Running plugin health checks...");
                let plugins_map = plugins_clone.read().expect("No data was found.");
                for (_name, managed_plugin) in plugins_map.iter() {
                    let new_health = managed_plugin.plugin.health_check();
                    let mut health_writer = managed_plugin
                        .health
                        .write()
                        .expect("Plugin healthy data not found.");
                    *health_writer = new_health;
                }
                drop(plugins_map);

                let stable_plugins_map = stable_plugins_clone.read().expect("No data was found.");
                for (_name, managed_plugin) in stable_plugins_map.iter() {
                    let new_health = managed_plugin.plugin.health_check().into_result().unwrap_or(PluginHealth::Error("Health check failed".to_string().into()));
                    let mut health_writer = managed_plugin
                        .health
                        .write()
                        .expect("Plugin healthy data not found.");
                    *health_writer = new_health;
                }
                drop(stable_plugins_map);

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
