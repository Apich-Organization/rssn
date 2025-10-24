#![allow(unsafe_code)]
#![allow(clippy::indexing_slicing)]
#![allow(clippy::no_mangle_with_rust_abi)]
#![allow(non_local_definitions)]
use crate::plugins::plugin_c::PluginHealth;
use abi_stable::sabi_trait;
use abi_stable::std_types::{RBox, RResult, RString, RVec};
use abi_stable::StableAbi;

#[allow(non_local_definitions)]
#[sabi_trait]
pub trait StablePlugin: Send + Sync {
    fn name(&self) -> RString;

    fn api_version(&self) -> RString;

    fn on_load(&self) -> RResult<(), RString>;

    fn execute(&self, command: RString, args: RVec<u8>) -> RResult<RVec<u8>, RString>;

    fn health_check(&self) -> RResult<PluginHealth, RString>;
}

#[repr(C)]
#[derive(StableAbi)]
#[sabi(kind(Prefix(prefix_ref = RoVtable)))]
#[sabi(missing_field(panic))]
pub struct StablePluginModule {
    pub name: extern "C" fn() -> RString,
    pub version: extern "C" fn() -> RString,
    pub new: extern "C" fn() -> StablePlugin_TO<'static, RBox<()>>,
}
