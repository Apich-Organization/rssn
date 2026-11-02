//! Loading a plugin from a shared library.
//!
//! # The ABI, and why it is bytes
//!
//! Everything crosses the boundary as one `bincode-next` buffer in and one out.
//! A `#[repr(C)]` mirror of [`Value`] would have to track every variant, every
//! field and every enum discriminant across two independently compiled
//! artifacts, and would break silently the day a variant is added. The
//! deterministic encoding is already the stable description of a value — it is
//! what the manifest is hashed with — so the ABI is a byte buffer, and the
//! types on both sides are `Vec<Value>` and `Result<Value, String>`.
//!
//! ```c
//! typedef struct { uint8_t *ptr; size_t len; } RssnBuf;
//!
//! typedef struct {
//!     uint32_t abi;                            /* RSSN_PLUGIN_ABI */
//!     void (*log)(uint32_t level, const uint8_t *msg, size_t len);
//! } RssnHost;
//!
//! /* bincode(Vec<Value>) -> bincode(Result<Value, String>) */
//! typedef RssnBuf (*RssnCallFn)(void *state, const uint8_t *args, size_t len);
//!
//! typedef struct {
//!     const char    *name;
//!     uint32_t       arity_min;
//!     uint32_t       arity_max;      /* 0xFFFFFFFF means variadic */
//!     const uint8_t *signatures;     /* bincode(Vec<Signature>) */
//!     size_t         signatures_len;
//!     RssnCallFn     call;
//! } RssnOpEntry;
//!
//! typedef struct {
//!     uint32_t             abi;
//!     const char          *name;
//!     const char          *version;
//!     const char          *description;
//!     const RssnOpEntry   *ops;
//!     size_t               op_count;
//!     void *(*init)(const RssnHost *host);
//!     void  (*free_buf)(RssnBuf);
//! } RssnPluginDesc;
//!
//! const RssnPluginDesc *rssn_plugin_v1(const RssnHost *host);
//! ```
//!
//! # Ownership
//!
//! The plugin allocates the buffer it returns and frees it through its own
//! `free_buf`; rssn never frees memory it did not allocate. rssn owns the
//! argument buffer it passes in and frees it after the call. Both sides
//! therefore only ever release their own memory, which is the one rule that
//! makes a cross-allocator boundary safe.
//!
//! # Why this is behind a feature
//!
//! Loading a shared library runs arbitrary code in the host process: there is
//! no sandbox, no capability check beyond the ABI number, and no way to recover
//! from a plugin that corrupts memory. It is therefore opt-in
//! (`plugin-cdylib`) and reported in [`Features::plugin_kinds`], so a host can
//! refuse a cdylib plugin with a precise diagnostic instead of failing
//! obscurely.
//!
//! [`Features::plugin_kinds`]: crate::kernel::manifest::Features::plugin_kinds

// The ABI is a C boundary: it dereferences pointers the plugin hands over.
#![allow(unsafe_code)]

use std::ffi::c_void;

use crate::kernel::signature::PLUGIN_ABI;

/// The plugin ABI this module implements.
///
/// A plugin compares it with the `abi` it was compiled against; a mismatch is
/// refused by the loader rather than producing a call through the wrong
/// function-pointer layout.
pub const ABI: u32 = PLUGIN_ABI;
#[cfg(feature = "plugin-cdylib")]
use std::path::Path;

#[cfg(feature = "plugin-cdylib")]
use crate::graph::Arity;
#[cfg(feature = "plugin-cdylib")]
use crate::graph::OpDescriptor;
#[cfg(feature = "plugin-cdylib")]
use crate::graph::op::EvalDyn;
#[cfg(feature = "plugin-cdylib")]
use crate::kernel::call::CallCtx;
#[cfg(feature = "plugin-cdylib")]
use crate::kernel::call::OpError;
#[cfg(feature = "plugin-cdylib")]
use crate::kernel::codec;
#[cfg(feature = "plugin-cdylib")]
use crate::kernel::plugin::Plugin;
#[cfg(feature = "plugin-cdylib")]
use crate::kernel::plugin::PluginError;
#[cfg(feature = "plugin-cdylib")]
use crate::kernel::plugin::PluginManifest;
#[cfg(feature = "plugin-cdylib")]
use crate::kernel::plugin::Registry;
#[cfg(feature = "plugin-cdylib")]
use crate::kernel::signature::Signature;
#[cfg(feature = "plugin-cdylib")]
use crate::kernel::value::Value;

/// A byte buffer owned by whoever produced it.
#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct RssnBuf {
    /// The bytes, or null when `len` is zero.
    pub ptr: *mut u8,
    /// How many bytes there are.
    pub len: usize,
}

#[cfg(feature = "plugin-cdylib")]
impl RssnBuf {
    /// The buffer as a slice.
    ///
    /// # Safety
    /// `ptr` must point to `len` initialised bytes that stay alive for `'a`.
    #[must_use]
    pub const unsafe fn as_slice<'a>(self) -> &'a [u8] {
        if self.ptr.is_null() || self.len == 0 {
            return &[];
        }
        // SAFETY: the caller guarantees the pointer and length.
        unsafe { std::slice::from_raw_parts(self.ptr, self.len) }
    }
}

/// What the host offers the plugin.
#[repr(C)]
pub struct RssnHost {
    /// [`PLUGIN_ABI`] as the host compiled it.
    pub abi: u32,
    /// Writes a message to the host log; `level` is 0 for debug, 1 for warn.
    pub log: Option<extern "C" fn(level: u32, message: *const u8, len: usize)>,
}

/// One operator the plugin provides.
#[repr(C)]
pub struct RssnOpEntry {
    /// NUL-terminated operator name.
    pub name: *const std::os::raw::c_char,
    /// Fewest accepted arguments.
    pub arity_min: u32,
    /// Most accepted arguments; `u32::MAX` means variadic.
    pub arity_max: u32,
    /// `bincode(Vec<Signature>)`.
    pub signatures: *const u8,
    /// Length of `signatures`.
    pub signatures_len: usize,
    /// The call itself.
    pub call: Option<RssnCallFn>,
}

/// `bincode(Vec<Value>)` in, `bincode(Result<Value, String>)` out.
pub type RssnCallFn = extern "C" fn(state: *mut c_void, args: *const u8, len: usize) -> RssnBuf;

/// What the plugin declares about itself.
#[repr(C)]
pub struct RssnPluginDesc {
    /// [`PLUGIN_ABI`] as the plugin compiled it.
    pub abi: u32,
    /// NUL-terminated `vendor/name`.
    pub name: *const std::os::raw::c_char,
    /// NUL-terminated version.
    pub version: *const std::os::raw::c_char,
    /// NUL-terminated one-line description.
    pub description: *const std::os::raw::c_char,
    /// The operator table.
    pub ops: *const RssnOpEntry,
    /// How many entries `ops` has.
    pub op_count: usize,
    /// Called once after loading; the returned pointer is passed to every
    /// `call` as `state`.
    pub init: Option<extern "C" fn(host: *const RssnHost) -> *mut c_void>,
    /// Frees a buffer this plugin returned.
    pub free_buf: Option<extern "C" fn(buffer: RssnBuf)>,
}

/// The symbol a plugin library must export.
pub const ENTRY_SYMBOL: &[u8] = b"rssn_plugin_v1\0";

/// A plugin loaded from a shared library.
#[cfg(feature = "plugin-cdylib")]
pub struct CdylibPlugin {
    /// Keeps the library mapped for as long as the plugin lives; dropping it
    /// would unmap the code the function pointers point into.
    _library: Library,
    /// What the plugin declared.
    manifest: PluginManifest,
    /// The operators, in declaration order.
    ops: Vec<LoadedOp>,
    /// The plugin's own per-instance state.
    state: *mut c_void,
    /// Frees a buffer the plugin returned.
    free_buf: Option<extern "C" fn(RssnBuf)>,
}

/// One operator read out of a plugin's descriptor.
#[cfg(feature = "plugin-cdylib")]
struct LoadedOp {
    name: String,
    arity: Arity,
    signatures: Vec<Signature>,
    call: RssnCallFn,
}

// The plugin promises `Send + Sync` by being a `Plugin`; the raw pointers are
// only ever used behind the trait's own bounds.
#[cfg(feature = "plugin-cdylib")]
unsafe impl Send for CdylibPlugin {}
#[cfg(feature = "plugin-cdylib")]
unsafe impl Sync for CdylibPlugin {}

#[cfg(feature = "plugin-cdylib")]
impl std::fmt::Debug for CdylibPlugin {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        f.debug_struct("CdylibPlugin")
            .field("name", &self.manifest.name)
            .field("ops", &self.ops.len())
            .finish_non_exhaustive()
    }
}

#[cfg(feature = "plugin-cdylib")]
impl CdylibPlugin {
    /// The operators this plugin provides.
    #[must_use]
    pub fn op_names(&self) -> Vec<&str> {
        self.ops.iter().map(|op| op.name.as_str()).collect()
    }
}

/// Reads a NUL-terminated C string, or the empty string.
#[cfg(feature = "plugin-cdylib")]
/// # Safety
/// `text` must be null or a valid NUL-terminated string.
unsafe fn c_str(text: *const std::os::raw::c_char) -> String {
    if text.is_null() {
        return String::new();
    }
    // SAFETY: the caller guarantees a NUL-terminated string.
    unsafe { std::ffi::CStr::from_ptr(text) }
        .to_str()
        .unwrap_or_default()
        .to_owned()
}

#[cfg(feature = "plugin-cdylib")]
impl Plugin for CdylibPlugin {
    fn manifest(&self) -> PluginManifest {
        self.manifest.clone()
    }

    fn register(
        &self,
        reg: &mut Registry<'_>,
    ) -> Result<(), PluginError> {
        for op in &self.ops {
            let mut descriptor = OpDescriptor::new(&op.name, op.arity);
            for signature in &op.signatures {
                descriptor = descriptor.signature(signature.clone());
            }
            let evaluator = CdylibEval {
                call: op.call,
                state: self.state,
                free_buf: self.free_buf,
                signature: op.signatures.first().cloned(),
            };
            descriptor = descriptor.eval_dyn(std::sync::Arc::new(evaluator));
            reg.register_op(descriptor)?;
        }
        Ok(())
    }
}

/// Forwards one call across the ABI.
#[cfg(feature = "plugin-cdylib")]
struct CdylibEval {
    call: RssnCallFn,
    state: *mut c_void,
    free_buf: Option<extern "C" fn(RssnBuf)>,
    signature: Option<Signature>,
}

#[cfg(feature = "plugin-cdylib")]
unsafe impl Send for CdylibEval {}
#[cfg(feature = "plugin-cdylib")]
unsafe impl Sync for CdylibEval {}

#[cfg(feature = "plugin-cdylib")]
impl EvalDyn for CdylibEval {
    fn eval(
        &self,
        args: &[Value],
        _ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let encoded = codec::encode_deterministic(&args.to_vec()).map_err(OpError::Codec)?;
        // SAFETY: the plugin promised this entry point follows the ABI in the
        // module documentation; it may not keep `encoded` after returning.
        let out = (self.call)(self.state, encoded.as_ptr(), encoded.len());
        // SAFETY: the plugin returned a buffer it allocated, and it stays valid
        // until `free_buf` is called below.
        let answer: Result<Value, String> = {
            let bytes = unsafe { out.as_slice() };
            codec::decode_deterministic(bytes).map_err(OpError::Codec)?
        };
        if let Some(free) = self.free_buf {
            // SAFETY: `out` came from this plugin and is freed exactly once.
            free(out);
        }
        answer.map_err(|message| OpError::Plugin(PluginError::Rejected(message)))
    }

    fn signature(&self) -> Signature {
        self.signature.clone().unwrap_or_else(|| {
            Signature::new(
                vec![crate::kernel::infer::term()],
                crate::kernel::infer::term(),
            )
        })
    }
}

/// Loads the plugin in `path`.
///
/// # Errors
/// Fails when the library cannot be opened, does not export
/// [`ENTRY_SYMBOL`], reports a different [`PLUGIN_ABI`], or declares an
/// operator without a signature.
#[cfg(feature = "plugin-cdylib")]
pub fn load(path: &Path) -> Result<CdylibPlugin, PluginError> {
    // SAFETY: opening an arbitrary shared library runs its initialisers. This
    // is why the whole entry point is behind a feature the host opts into.
    let library = unsafe { Library::new(path) }.map_err(|error| {
        PluginError::Registration(format!("cannot open {}: {error}", path.display()))
    })?;
    // SAFETY: the symbol is read as the entry point the ABI declares.
    let entry: Symbol<'_, EntryFn> = unsafe { library.get(ENTRY_SYMBOL) }.map_err(|error| {
        PluginError::Registration(format!("{} has no rssn_plugin_v1: {error}", path.display()))
    })?;
    let host = RssnHost {
        abi: PLUGIN_ABI,
        log: None,
    };
    let descriptor = (*entry)(&raw const host);
    if descriptor.is_null() {
        return Err(PluginError::Registration(String::from(
            "the plugin returned a null descriptor",
        )));
    }
    // SAFETY: the plugin returned a non-null descriptor that stays valid for
    // as long as the library is loaded.
    let descriptor = unsafe { &*descriptor };
    if descriptor.abi != PLUGIN_ABI {
        return Err(PluginError::AbiMismatch {
            expected: PLUGIN_ABI,
            got: descriptor.abi,
        });
    }
    // SAFETY: every pointer below belongs to the descriptor the plugin owns.
    let manifest = PluginManifest {
        name: unsafe { c_str(descriptor.name) },
        version: unsafe { c_str(descriptor.version) },
        abi: descriptor.abi,
        kind: crate::kernel::plugin::PluginKind::Cdylib,
        capabilities: vec![String::from("ops"), String::from("cdylib")],
        description: unsafe { c_str(descriptor.description) },
    };
    if manifest.name.is_empty() {
        return Err(PluginError::Registration(String::from(
            "the plugin has no name",
        )));
    }
    let mut ops = Vec::with_capacity(descriptor.op_count);
    for index in 0..descriptor.op_count {
        // SAFETY: `ops` has `op_count` entries.
        let Some(entry) = (unsafe { descriptor.ops.add(index).as_ref() }) else {
            break;
        };
        let Some(call) = entry.call else {
            return Err(PluginError::Registration(String::from(
                "an operator has no call function",
            )));
        };
        let signatures: Vec<Signature> = if entry.signatures.is_null() || entry.signatures_len == 0
        {
            Vec::new()
        } else {
            // SAFETY: the plugin declared this many bytes at this address.
            let bytes =
                unsafe { std::slice::from_raw_parts(entry.signatures, entry.signatures_len) };
            codec::decode_deterministic(bytes).map_err(|error| {
                PluginError::Registration(format!("cannot read the signatures: {error}"))
            })?
        };
        if signatures.is_empty() {
            return Err(PluginError::Registration(format!(
                "operator `{}` declares no signature",
                unsafe { c_str(entry.name) }
            )));
        }
        let arity = if entry.arity_max == u32::MAX {
            Arity::Variadic
        } else {
            let max = u8::try_from(entry.arity_max).unwrap_or(u8::MAX);
            Arity::Fixed(max)
        };
        ops.push(LoadedOp {
            name: unsafe { c_str(entry.name) },
            arity,
            signatures,
            call,
        });
    }
    let state = descriptor
        .init
        .map_or(std::ptr::null_mut(), |init| init(&raw const host));
    Ok(CdylibPlugin {
        _library: library,
        manifest,
        ops,
        state,
        free_buf: descriptor.free_buf,
    })
}

/// The entry point signature a plugin exports.
#[cfg(feature = "plugin-cdylib")]
type EntryFn = extern "C" fn(host: *const RssnHost) -> *const RssnPluginDesc;

/// The dynamic loader, re-exported so a host can see what it is using.
#[cfg(feature = "plugin-cdylib")]
pub use libloading::Library;
#[cfg(feature = "plugin-cdylib")]
use libloading::Symbol;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[cfg(feature = "plugin-cdylib")]
    fn the_buffer_reads_as_empty_when_it_is_null() {
        let buffer = RssnBuf {
            ptr: std::ptr::null_mut(),
            len: 0,
        };
        // SAFETY: a null buffer with zero length is the empty slice.
        assert_eq!(unsafe { buffer.as_slice() }.len(), 0);
    }

    #[test]
    fn a_missing_library_is_reported_with_its_path() {
        #[cfg(feature = "plugin-cdylib")]
        {
            let error = load(Path::new("/nonexistent/plugin.so"));
            assert!(error.is_err(), "a missing library must not load");
        }
    }

    #[test]
    fn the_entry_symbol_is_nul_terminated() {
        assert_eq!(ENTRY_SYMBOL.last(), Some(&0_u8));
        assert_eq!(&ENTRY_SYMBOL[..ENTRY_SYMBOL.len() - 1], b"rssn_plugin_v1");
    }
}
