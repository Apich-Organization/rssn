#![allow(unsafe_code)]
#![allow(clippy::indexing_slicing)]
#![allow(
    clippy::no_mangle_with_rust_abi
)]
#![allow(non_local_definitions)]

use abi_stable::sabi_trait;
use abi_stable::std_types::RBox;
use abi_stable::std_types::RResult;
use abi_stable::std_types::RString;
use abi_stable::std_types::RVec;
use abi_stable::StableAbi;

use crate::plugins::plugin_c::PluginHealth;

#[allow(non_local_definitions)]
#[sabi_trait]
/// This trait defines the stable ABI for plugins.
pub trait StablePlugin:
    Send + Sync
{
    /// Returns the name of the plugin.
    fn name(&self) -> RString;

    /// Returns the API version of the plugin.
    fn api_version(&self) -> RString;

    /// Called when the plugin is loaded.
    fn on_load(
        &self
    ) -> RResult<(), RString>;

    /// Executes a command on the plugin.
    fn execute(
        &self,
        command: RString,
        args: RVec<u8>,
    ) -> RResult<RVec<u8>, RString>;

    /// Performs a health check on the plugin.
    fn health_check(
        &self
    ) -> RResult<PluginHealth, RString>;
}

#[repr(C)]
#[derive(StableAbi)]
#[sabi(kind(Prefix(prefix_ref = RoVtable)))]
#[sabi(missing_field(panic))]
/// This struct defines the stable ABI for a plugin module.
pub struct StablePluginModule {
    /// Returns the name of the plugin.
    pub name : extern "C" fn() -> RString,
    /// Returns the version of the plugin.
    pub version : extern "C" fn() -> RString,
    /// Creates a new instance of the plugin.
    pub new : extern "C" fn() -> StablePlugin_TO<'static, RBox<()>>,
}
