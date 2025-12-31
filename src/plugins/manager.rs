//! # Plugin Manager
//! Responsible for loading, managing, and interacting with plugins.
#![allow(unsafe_code)]
#![allow(clippy::indexing_slicing)]
#![allow(
    clippy::no_mangle_with_rust_abi
)]

use std::collections::HashMap;
use std::error::Error;
use std::sync::Arc;
use std::sync::LazyLock;
use std::sync::RwLock;
use std::sync::atomic::AtomicBool;
use std::sync::atomic::Ordering;
use std::thread;
use std::time::Duration;

use bincode_next::config;
use bincode_next::serde;
use libloading::Library;
use libloading::Symbol;

/// Global instance of the PluginManager for FFI access.

pub static GLOBAL_PLUGIN_MANAGER:
    LazyLock<RwLock<PluginManager>> =
    LazyLock::new(|| {

        RwLock::new(
            PluginManager::empty(),
        )
    });

use crate::plugins::plugin_c::Plugin;
use crate::plugins::plugin_c::PluginError;
use crate::plugins::plugin_c::PluginErrorKind;
use crate::plugins::plugin_c::PluginHealth;
use crate::symbolic::core::Expr;

/// The expected signature for the function that instantiates the plugin.
///
/// Each plugin library must export a C-compatible function named `_plugin_create`
/// with this signature, returning a pointer to a Box containing the plugin trait object.

type PluginCreate =
    unsafe extern "C" fn() -> *mut Box<
        dyn Plugin,
    >;

use abi_stable::std_types::RBox;

/// Holds the plugin instance and its current health state.
use crate::plugins::stable_abi::StablePlugin_TO;
/// Holds the plugin instance and its current health state.
use crate::plugins::stable_abi::StablePluginModule;

/// Holds the plugin instance and its current health state.

pub struct ManagedPlugin {
    /// The plugin instance
    pub plugin: Box<dyn Plugin>,
    /// The plugin health state
    pub health: RwLock<PluginHealth>,
}

/// A plugin that is managed by the plugin manager.

pub struct ManagedStablePlugin {
    /// The plugin instance.
    pub plugin: StablePlugin_TO<
        'static,
        RBox<()>,
    >,
    /// The plugin's health status.
    pub health: RwLock<PluginHealth>,
}

/// Manages the lifecycle of all loaded plugins.

pub struct PluginManager {
    /// A map of plugin names to their managed instances.
    pub plugins: Arc<
        RwLock<
            HashMap<
                String,
                ManagedPlugin,
            >,
        >,
    >,
    /// A map of stable plugin names to their managed instances.
    pub stable_plugins: Arc<
        RwLock<
            HashMap<
                String,
                ManagedStablePlugin,
            >,
        >,
    >,
    /// The handle to the health check thread.
    pub health_check_thread:
        Option<thread::JoinHandle<()>>,
    /// A signal to stop the health check thread.
    pub stop_signal: Arc<AtomicBool>,
    /// A list of loaded libraries.
    pub libraries: Vec<Library>,
}

impl PluginManager {
    /// Creates a new `PluginManager` without loading any plugins.

    #[must_use]

    pub fn empty() -> Self {

        Self {
            plugins: Arc::new(
                RwLock::new(
                    HashMap::new(),
                ),
            ),
            stable_plugins: Arc::new(
                RwLock::new(
                    HashMap::new(),
                ),
            ),
            health_check_thread: None,
            stop_signal: Arc::new(
                AtomicBool::new(false),
            ),
            libraries: Vec::new(),
        }
    }

    /// Creates a new `PluginManager` and loads plugins from a specified directory.
    ///
    /// # Errors
    ///
    /// This function will return an error if it fails to load plugins from the
    /// specified directory.

    pub fn new(
        plugin_dir: &str
    ) -> Result<Self, Box<dyn Error>>
    {

        let mut manager = Self {
            plugins: Arc::new(
                RwLock::new(
                    HashMap::new(),
                ),
            ),
            stable_plugins: Arc::new(
                RwLock::new(
                    HashMap::new(),
                ),
            ),
            health_check_thread: None,
            stop_signal: Arc::new(
                AtomicBool::new(false),
            ),
            libraries: Vec::new(),
        };

        manager
            .load_plugins(plugin_dir)?;

        manager.start_health_checks(
            Duration::from_secs(30),
        );

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
    ///
    /// # Errors
    ///
    /// This function will return a `PluginError` if:
    /// - The plugin with `plugin_name` is not found.
    /// - Serialization of arguments or deserialization of results fails.
    /// - The plugin's execution of the command fails.
    ///
    /// # Panics
    ///
    /// This function may panic if the internal `RwLock` protecting the plugin maps
    /// is poisoned (e.g., a thread holding the lock has panicked).

    pub fn execute_plugin(
        &self,
        plugin_name: &str,
        command: &str,
        args: &Expr,
    ) -> Result<Expr, PluginError> {

        if let Some(managed_plugin) =
            self.stable_plugins
                .read()
                .expect("Lock poisoned")
                .get(plugin_name)
        {

            let config =
                config::standard();

            let args_vec = serde::encode_to_vec(args, config).map_err(|e| {
                PluginError::new(PluginErrorKind::SerializationError, &e.to_string())
            })?;

            let result_vec = managed_plugin
                .plugin
                .execute(command.into(), args_vec.into())
                .into_result()
                .map_err(|e| PluginError::new(PluginErrorKind::ExecutionFailed, &e.to_string()))?;

            let (result_expr, _) = serde::decode_from_slice(&result_vec, config)
                .map_err(|e| PluginError::new(PluginErrorKind::SerializationError, &e.to_string()))?;

            return Ok(result_expr);
        }

        if let Some(managed_plugin) =
            self.plugins
                .read()
                .expect("Lock poisoned")
                .get(plugin_name)
        {

            return managed_plugin
                .plugin
                .execute(
                    command,
                    args,
                );
        }

        Err(PluginError::new(
            PluginErrorKind::NotFound,
            &format!("Plugin '{plugin_name}' not found"),
        ))
    }

    /// Registers a plugin manually. Useful for testing or internal plugins.
    ///
    /// # Panics
    ///
    /// Panics if the internal `RwLock` for plugins is poisoned.

    pub fn register_plugin(
        &self,
        plugin: Box<dyn Plugin>,
    ) {

        let name = plugin
            .name()
            .to_string();

        let mut plugins_map = self
            .plugins
            .write()
            .expect("Lock poisoned");

        plugins_map.insert(
            name,
            ManagedPlugin {
                plugin,
                health: RwLock::new(
                    PluginHealth::Ok,
                ),
            },
        );
    }

    /// Returns a list of all loaded plugin names.
    ///
    /// # Panics
    ///
    /// Panics if the internal `RwLock` for plugins or stable plugins is poisoned.
    #[must_use]

    pub fn get_loaded_plugin_names(
        &self
    ) -> Vec<String> {

        let mut names = Vec::new();

        names.extend(
            self.plugins
                .read()
                .expect("Lock poisoned")
                .keys()
                .cloned(),
        );

        names.extend(
            self.stable_plugins
                .read()
                .expect("Lock poisoned")
                .keys()
                .cloned(),
        );

        names
    }

    /// Unloads a plugin by name.
    ///
    /// # Panics
    ///
    /// Panics if the internal `RwLock` for plugins or stable plugins is poisoned.
    #[must_use]

    pub fn unload_plugin(
        &self,
        name: &str,
    ) -> bool {

        if self
            .plugins
            .write()
            .expect("Lock poisoned")
            .remove(name)
            .is_some()
        {

            return true;
        }

        self.stable_plugins
            .write()
            .expect("Lock poisoned")
            .remove(name)
            .is_some()
    }

    /// Gets metadata for all plugins.
    ///
    /// # Panics
    ///
    /// Panics if the internal `RwLock` for plugins is poisoned.
    #[must_use]

    pub fn get_plugin_metadata(
        &self
    ) -> HashMap<
        String,
        HashMap<String, String>,
    > {

        let mut all_metadata =
            HashMap::new();

        {

            let plugins_map = self
                .plugins
                .read()
                .expect(
                    "Lock poisoned",
                );

            for (name, managed) in
                plugins_map.iter()
            {

                all_metadata.insert(
                    name.clone(),
                    managed
                        .plugin
                        .metadata(),
                );
            }
        }

        all_metadata
    }

    /// Scans a directory for dynamic libraries and attempts to load them as plugins.

    pub(crate) fn load_plugins(
        &mut self,
        directory: &str,
    ) -> Result<(), Box<dyn Error>>
    {

        for entry in std::fs::read_dir(
            directory,
        )? {

            let entry = entry?;

            let path = entry.path();

            if path.is_file()
                && path
                    .extension()
                    .is_some_and(|ext| ext == std::env::consts::DLL_EXTENSION)
            {

                println!(
                    "Attempting to load plugin: {:?}",
                    path.display()
                );

                unsafe {

                    if self
                        .load_stable_plugin(&path)
                        .is_err()
                    {

                        self.load_plugin(&path)?;
                    }
                }
            }
        }

        Ok(())
    }

    /// Loads a stable plugin from a dynamic library file.
    ///
    /// # Errors
    ///
    /// This function will return an error if:
    /// - The library cannot be loaded.
    /// - The `STABLE_PLUGIN_MODULE` symbol cannot be found.
    /// - The plugin's API version is incompatible with the current crate version.
    /// - The plugin's `on_load` method returns an error.
    ///
    /// # Panics
    ///
    /// This function will panic if:
    /// - `self.libraries.last()` returns `None`, indicating an internal logic error where a library was added but not found.
    /// - The `RwLock` for `self.stable_plugins` is poisoned.

    unsafe fn load_stable_plugin(
        &mut self,
        library_path: &std::path::Path,
    ) -> Result<(), Box<dyn Error>>
    {

        unsafe {

            let library = Library::new(
                library_path,
            )?;

            self.libraries
                .push(library);

            let library = self
                .libraries
                .last()
                .expect(
                    "Library invalid",
                );

            let module: Symbol<
            '_,
            *const StablePluginModule,
        > = library.get(
            b"STABLE_PLUGIN_MODULE",
        )?;

            let module = &**module;

            let plugin = (module.new)();

            let plugin_name =
                (module.name)();

            let plugin_api_version =
                (module.version)();

            let crate_version = env!(
                "CARGO_PKG_VERSION"
            );

            let (p_major, p_minor) =
                parse_version(
                    &plugin_api_version,
                )?;

            let (c_major, c_minor) =
                parse_version(
                    crate_version,
                )?;

            if p_major != c_major
                || p_minor != c_minor
            {

                return Err(format!(
                "Plugin '{plugin_name}' has \
                 incompatible API \
                 version {plugin_api_version}. Expected \
                 a version compatible \
                 with {crate_version}."
            )
            .into());
            }

            println!(
            "Successfully loaded \
             stable plugin: {plugin_name} (API \
             version: {plugin_api_version})"
        );

            plugin
                .on_load()
                .into_result()
                .map_err(|e| {

                    e.to_string()
                })?;

            let managed_plugin =
            ManagedStablePlugin {
                plugin,
                health: RwLock::new(
                    PluginHealth::Ok,
                ),
            };

            self.stable_plugins
                .write()
                .expect(
                    "Unexpected Error",
                )
                .insert(
                    plugin_name
                        .to_string(),
                    managed_plugin,
                );

            Ok(())
        }
    }

    /// Loads a single plugin from a dynamic library file.

    /// Loads a standard plugin from a dynamic library file.
    ///
    /// # Errors
    ///
    /// This function will return an error if:
    /// - The library cannot be loaded.
    /// - The `_plugin_create` symbol cannot be found.
    /// - The plugin's API version is incompatible with the current crate version.
    /// - The plugin's `on_load` method returns an error.
    ///
    /// # Panics
    ///
    /// This function will panic if:
    /// - `self.libraries.last()` returns `None`, indicating an internal logic error where a library was added but not found.
    /// - The `RwLock` for `self.plugins` is poisoned.

    unsafe fn load_plugin(
        &mut self,
        library_path: &std::path::Path,
    ) -> Result<(), Box<dyn Error>>
    {

        unsafe {

            let library = Library::new(
                library_path,
            )?;

            self.libraries
                .push(library);

            let library = self
                .libraries
                .last()
                .expect(
                    "Library invalid",
                );

            let constructor: Symbol<
                '_,
                PluginCreate,
            > = library.get(
                b"_plugin_create",
            )?;

            let plugin_box_ptr =
                constructor();

            let plugin_box =
                Box::from_raw(
                    plugin_box_ptr,
                );

            let plugin = *plugin_box;

            let plugin_api_version =
                plugin.api_version();

            let crate_version = env!(
                "CARGO_PKG_VERSION"
            );

            let (p_major, p_minor) =
                parse_version(
                    plugin_api_version,
                )?;

            let (c_major, c_minor) =
                parse_version(
                    crate_version,
                )?;

            if p_major != c_major
                || p_minor != c_minor
            {

                return Err(format!(
                "Plugin '{}' has \
                 incompatible API \
                 version {}. Expected \
                 a version compatible \
                 with {}.",
                plugin.name(),
                plugin_api_version,
                crate_version
            )
                .into());
            }

            println!(
                "Successfully loaded \
                 plugin: {} (API \
                 version: {})",
                plugin.name(),
                plugin_api_version
            );

            plugin.on_load()?;

            let managed_plugin =
            ManagedPlugin {
                plugin,
                health: RwLock::new(
                    PluginHealth::Ok,
                ),
            };

            self.plugins
                .write()
                .expect(
                    "Unexpected Error",
                )
                .insert(
                    managed_plugin
                        .plugin
                        .name()
                        .to_string(),
                    managed_plugin,
                );

            Ok(())
        }
    }

    /// Spawns a background thread to periodically perform health checks on all plugins.

    pub(crate) fn start_health_checks(
        &mut self,
        interval: Duration,
    ) {

        let plugins_clone =
            Arc::clone(&self.plugins);

        let stable_plugins_clone =
            Arc::clone(
                &self.stable_plugins,
            );

        let stop_signal_clone =
            Arc::clone(
                &self.stop_signal,
            );

        let handle = thread::spawn(
            move || {

                while !stop_signal_clone.load(Ordering::SeqCst) {

                println!("Running plugin health checks...");

                let plugins_map = plugins_clone
                    .read()
                    .expect("No data was found.");

                for (_name, managed_plugin) in plugins_map.iter() {

                    let new_health = managed_plugin
                        .plugin
                        .health_check();

                    let mut health_writer = managed_plugin
                        .health
                        .write()
                        .expect("Plugin healthy data not found.");

                    *health_writer = new_health;
                }

                drop(plugins_map);

                let stable_plugins_map = stable_plugins_clone
                    .read()
                    .expect("No data was found.");

                for (_name, managed_plugin) in stable_plugins_map.iter() {

                    let new_health = managed_plugin
                        .plugin
                        .health_check()
                        .into_result()
                        .unwrap_or(PluginHealth::Error(
                            "Health check failed"
                                .to_string()
                                .into(),
                        ));

                    let mut health_writer = managed_plugin
                        .health
                        .write()
                        .expect("Plugin healthy data not found.");

                    *health_writer = new_health;
                }

                drop(stable_plugins_map);

                thread::sleep(interval);
            }

                println!(
                    "Plugin health \
                     check thread \
                     shutting down."
                );
            },
        );

        self.health_check_thread =
            Some(handle);
    }
}

impl Drop for PluginManager {
    fn drop(&mut self) {

        println!(
            "Shutting down plugin \
             manager..."
        );

        self.stop_signal
            .store(
                true,
                Ordering::SeqCst,
            );

        if let Some(handle) = self
            .health_check_thread
            .take()
        {

            handle
                .join()
                .expect(
                    "Failed to join \
                     health check \
                     thread",
                );
        }

        println!(
            "Plugin manager shut down."
        );
    }
}

/// A helper to parse a semantic version string into (major, minor) components.

pub(crate) fn parse_version(
    version: &str
) -> Result<(u32, u32), Box<dyn Error>>
{

    let mut parts = version.split('.');

    let major = parts
        .next()
        .ok_or(
            "Invalid version string",
        )?
        .parse::<u32>()?;

    let minor = parts
        .next()
        .ok_or(
            "Invalid version string",
        )?
        .parse::<u32>()?;

    Ok((major, minor))
}
