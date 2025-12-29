//! Handle-based FFI API for the Plugin Manager.

use std::os::raw::c_char;

use crate::ffi_apis::common::{c_str_to_str, to_json_string};
use crate::plugins::manager::GLOBAL_PLUGIN_MANAGER;
use crate::symbolic::handles::HANDLE_MANAGER;

/// Loads plugins from a specified directory.
///
/// # Arguments
/// * `path` - Path to the plugin directory.
///
/// # Returns
/// True if successful, false otherwise.
#[no_mangle]

pub unsafe extern "C" fn rssn_plugins_load(
    path: *const c_char
) -> bool {

    if let Some(path_str) =
        c_str_to_str(path)
    {

        let mut manager =
            match GLOBAL_PLUGIN_MANAGER
                .write()
            {
                | Ok(m) => m,
                | Err(_) => {
                    return false
                },
            };

        match manager
            .load_plugins(path_str)
        {
            | Ok(()) => true,
            | Err(e) => {

                // Ideally log error
                eprintln!(
                    "Failed to load \
                     plugins: {e}"
                );

                false
            },
        }
    } else {

        false
    }
}

/// Returns a JSON array of loaded plugin names.
///
/// The caller must free the string using `rssn_free_string`.
#[no_mangle]

pub extern "C" fn rssn_plugins_get_loaded(
) -> *mut c_char {

    let manager = match GLOBAL_PLUGIN_MANAGER.read() {
        Ok(m) => m,
        Err(_) => return std::ptr::null_mut(),
    };

    let names = manager
        .get_loaded_plugin_names();

    to_json_string(&names)
}

/// Unloads a plugin by name.
#[no_mangle]

pub unsafe extern "C" fn rssn_plugins_unload(
    name: *const c_char
) -> bool {

    if let Some(name_str) =
        c_str_to_str(name)
    {

        let manager =
            match GLOBAL_PLUGIN_MANAGER
                .read()
            {
                | Ok(m) => m,
                | Err(_) => {
                    return false
                },
            };

        manager.unload_plugin(name_str)
    } else {

        false
    }
}

/// Executes a plugin command.
///
/// # Arguments
/// * `name` - Plugin name.
/// * `command` - Command string.
/// * `args_handle` - Handle to the argument expression.
///
/// # Returns
/// Handle to the result expression, or 0 on error.
#[no_mangle]

pub unsafe extern "C" fn rssn_plugins_execute(
    name: *const c_char,
    command: *const c_char,
    args_handle: usize,
) -> usize {

    let name_str =
        match c_str_to_str(name) {
            | Some(s) => s,
            | None => return 0,
        };

    let command_str =
        match c_str_to_str(command) {
            | Some(s) => s,
            | None => return 0,
        };

    let args_expr = match HANDLE_MANAGER
        .get(args_handle)
    {
        | Some(expr) => expr,
        | None => return 0,
    };

    let manager =
        match GLOBAL_PLUGIN_MANAGER
            .read()
        {
            | Ok(m) => m,
            | Err(_) => return 0,
        };

    // We need to clone args to pass potentially? execute takes &Expr.
    // PluginManager::execute_plugin takes &Expr.

    match manager.execute_plugin(
        name_str,
        command_str,
        &args_expr,
    ) {
        | Ok(result_expr) => {
            HANDLE_MANAGER
                .insert(result_expr)
        },
        | Err(e) => {

            eprintln!(
                "Plugin execution \
                 failed: {e}"
            );

            0
        },
    }
}
