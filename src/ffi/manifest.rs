//! The manifest, over the C ABI.
//!
//! A host that only speaks C still needs to know what this build of rssn
//! contains: which operators exist, what they accept, which rule sets are
//! installed, and a hash it can cache against. Two calls cover that.
//!
//! The manifest is built once per process from the *standard* rule sets and
//! cached: it is a property of the build, not of a session, so there is
//! nothing session-specific to report. A host that installed extra plugins
//! reads the richer, per-session manifest through the Rust API
//! ([`Session::manifest`](crate::api::Session::manifest)).
//!
//! ```c
//! char *json = rssn_manifest_json();
//! uint64_t hash = rssn_manifest_hash();
//! /* ... */
//! rssn_string_free(json);
//! ```

use std::sync::OnceLock;

use super::to_c;
use crate::kernel::manifest::KernelManifest;

/// The process-wide manifest of the standard build.
fn standard() -> &'static KernelManifest {
    static MANIFEST: OnceLock<KernelManifest> = OnceLock::new();
    MANIFEST.get_or_init(|| {
        let sets = crate::rules::standard();
        let mut graph = crate::graph::Graph::new();
        // A build whose own rule sets do not install is broken, not
        // recoverable; the manifest is then the empty one, which hashes
        // consistently and says so.
        if crate::graph::Engine::install(&mut graph, &sets).is_err() {
            return KernelManifest::default();
        }
        KernelManifest::build(&graph, &sets, &[])
    })
}

/// The manifest as JSON, for a host that wants to read it.
///
/// Returns a NUL-terminated string the caller owns and must release with
/// [`rssn_string_free`](super::rssn_string_free). Never returns null.
///
/// JSON is used here and nowhere else: it is a human-facing debugging
/// format. Anything rssn sends to another rssn goes through
/// [`KernelManifest::to_cbor`], which is deterministic.
#[unsafe(no_mangle)]
pub extern "C" fn rssn_manifest_json() -> *mut std::os::raw::c_char {
    to_c(standard().to_json())
}

/// The first 64 bits of the manifest hash, as a cache key.
///
/// The full 256-bit hash is in the manifest itself and in
/// [`KernelManifest::hash_u64`]'s sibling `manifest_hash`. This entry point
/// exists because a `uint64_t` is the widest integer every C host has.
#[unsafe(no_mangle)]
pub extern "C" fn rssn_manifest_hash() -> u64 {
    standard().hash_u64()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_manifest_crosses_the_c_boundary() {
        let json = rssn_manifest_json();
        assert!(!json.is_null());
        // SAFETY: `rssn_manifest_json` returns a live `CString`.
        let text = unsafe { std::ffi::CStr::from_ptr(json) }
            .to_str()
            .unwrap_or_default()
            .to_owned();
        // SAFETY: the pointer came from `CString::into_raw` and is freed once.
        unsafe { super::super::rssn_string_free(json) };
        let manifest: serde_json::Value = serde_json::from_str(&text).expect("manifest JSON");
        assert_eq!(manifest["api_version"].as_u64(), Some(u64::from(crate::kernel::API_VERSION)));
        assert!(text.contains("\"ops\""));
        assert_eq!(rssn_manifest_hash(), standard().hash_u64());
    }

    #[test]
    fn the_hash_is_stable_within_a_process() {
        assert_eq!(rssn_manifest_hash(), rssn_manifest_hash());
    }
}
