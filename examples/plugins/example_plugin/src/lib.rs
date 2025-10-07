use rssn::plugins::{Plugin, PluginError, PluginHealth};
use rssn::symbolic::core::Expr;

/// A simple example plugin that returns whatever expression it receives.
struct IdentityPlugin;

impl Plugin for IdentityPlugin {
    fn name(&self) -> &'static str {
        "identity"
    }

    fn api_version(&self) -> &'static str {
        // This must match the version of the main rssn crate it was compiled against.
        env!("CARGO_PKG_VERSION")
    }

    fn on_load(&self) -> Result<(), PluginError> {
        println!("IdentityPlugin loaded!");
        Ok(())
    }

    fn execute(&self, command: &str, args: &Expr) -> Result<Expr, PluginError> {
        println!(
            "IdentityPlugin executing command '{}' with args: {}",
            command,
            args.to_string()
        );
        // This simple plugin just returns the arguments it received, regardless of the command.
        Ok(args.clone())
    }

    fn health_check(&self) -> PluginHealth {
        PluginHealth::Ok
    }
}

/// The mandatory entry point for the plugin.
///
/// This function is called by the `PluginManager` to create an instance of the plugin.
#[no_mangle]
pub extern "C" fn _plugin_create() -> *mut Box<dyn Plugin> {
    // Create the plugin instance and put it into a Box.
    let plugin_box: Box<dyn Plugin> = Box::new(IdentityPlugin);
    // Wrap the trait object Box in another Box to create a thin pointer,
    // then return it as a raw pointer.
    Box::into_raw(Box::new(plugin_box))
}
