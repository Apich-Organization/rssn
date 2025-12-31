//! Handle-based FFI API for the Plugin Manager.

use std::os::raw::c_char;

use crate::ffi_apis::common::c_str_to_str;
use crate::ffi_apis::common::to_json_string;
use crate::plugins::manager::GLOBAL_PLUGIN_MANAGER;
use crate::symbolic::handles::HANDLE_MANAGER;

/// Loads plugins from a specified directory.
///
/// # Arguments
/// * `path` - Path to the plugin directory.
///
/// # Returns
/// True if successful, false otherwise.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_plugins_load(
    path: *const c_char
) -> bool {

    unsafe {

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
}

/// Returns a JSON array of loaded plugin names.
///
/// The caller must free the string using `rssn_free_string`.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_plugins_get_loaded()
-> *mut c_char {

    let names = match GLOBAL_PLUGIN_MANAGER.read() {
        Ok(m) => m.get_loaded_plugin_names(),
        Err(_) => return std::ptr::null_mut(),
    };

    to_json_string(&names)
}

/// Unloads a plugin by name.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_plugins_unload(
    name: *const c_char
) -> bool {

    unsafe {

        if let Some(name_str) =
            c_str_to_str(name)
        {

            match GLOBAL_PLUGIN_MANAGER
                .read()
            {
                | Ok(m) => {
                    m.unload_plugin(
                        name_str,
                    )
                },
                | Err(_) => false,
            }
        } else {

            false
        }
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
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_plugins_execute(
    name: *const c_char,
    command: *const c_char,
    args_handle: usize,
) -> usize {

    unsafe {

        let name_str =
            match c_str_to_str(name) {
                | Some(s) => s,
                | None => return 0,
            };

        let command_str =
            match c_str_to_str(command)
            {
                | Some(s) => s,
                | None => return 0,
            };

        let args_expr =
            match HANDLE_MANAGER
                .get(args_handle)
            {
                | Some(expr) => expr,
                | None => return 0,
            };

        let result =
            match GLOBAL_PLUGIN_MANAGER
                .read()
            {
                | Ok(m) => {
                    m.execute_plugin(
                        name_str,
                        command_str,
                        &args_expr,
                    )
                },
                | Err(_) => return 0,
            };

        match result {
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
}
