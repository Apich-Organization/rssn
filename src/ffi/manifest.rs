//! The manifest, over the C ABI.
//!
//! A host that only speaks C still needs to know what this build of rssn
//! contains: which operators exist, what they accept, which rule sets are
//! installed, and a hash it can cache against.
//!
//! The first C manifest request captures the shared standard manifest at
//! the current process-wide unit table revision. The hash and the JSON use
//! that same immutable snapshot. A host that installed extra plugins reads
//! the richer, per-session manifest through the Rust API
//! ([`Session::manifest`](crate::api::Session::manifest)).
//!
//! ```c
//! char *json = rssn_manifest_json();
//! uint64_t hash = rssn_manifest_hash();
//! /* ... */
//! rssn_string_free(json);
//! ```

use super::RssnStatus;
use std::ffi::CString;
use std::ffi::c_char;
use std::sync::Arc;
use std::sync::OnceLock;

use crate::kernel::manifest::KernelManifest;
use crate::kernel::manifest::OpInfo;

/// The process-wide manifest of the standard build.
fn standard() -> &'static KernelManifest {
    static MANIFEST: OnceLock<Arc<KernelManifest>> = OnceLock::new();
    MANIFEST.get_or_init(crate::api::standard_manifest)
}

/// The manifest as JSON, for a host that wants to read it.
///
/// Returns a NUL-terminated string the caller owns and must release with
/// [`rssn_string_free`](super::rssn_string_free). Never returns null.
///
/// This JSON view is intended for inspecting build metadata. Use
/// [`KernelManifest::to_cbor`] for deterministic, self-describing exchange.
#[unsafe(no_mangle)]
pub extern "C" fn rssn_manifest_json() -> *mut std::os::raw::c_char {
    static JSON: OnceLock<CString> = OnceLock::new();
    // JSON escapes interior NULs. Each caller owns its copy of the encoding.
    JSON.get_or_init(|| CString::new(standard().to_json()).unwrap_or_default())
        .clone()
        .into_raw()
}

fn lookup_operator(name: &str) -> Result<&'static OpInfo, RssnStatus> {
    standard().op(name).ok_or_else(|| {
        super::fail(
            RssnStatus::UnknownOperator,
            &format!("unknown operator {name}"),
        )
    })
}

/// Writes one operator's live manifest entry as a NUL-terminated JSON string.
///
/// Accepts a registered name or alias. The record includes the canonical name,
/// overload keys, type signatures, their sources, and dimensional metadata.
/// The caller releases the successful result with
/// [`rssn_string_free`](super::rssn_string_free).
///
/// # Safety
/// `name` must be null or a readable NUL-terminated string. `out` must be
/// null or valid for writing a string pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_op_info_json(
    name: *const c_char,
    out: *mut *mut c_char,
) -> RssnStatus {
    super::guard(|| {
        if out.is_null() {
            return super::fail(RssnStatus::NullArgument, "null output pointer");
        }
        // SAFETY: the caller supplies a readable NUL-terminated string or null.
        let name = match unsafe { super::text(name) } {
            | Ok(name) => name,
            | Err(status) => return status,
        };
        let operator = match lookup_operator(name) {
            | Ok(operator) => operator,
            | Err(status) => return status,
        };
        let json = match serde_json::to_string(operator) {
            | Ok(json) => json,
            | Err(error) => return super::fail(RssnStatus::Panic, &error.to_string()),
        };
        let pointer = super::to_c(json);
        if pointer.is_null() {
            return super::fail(RssnStatus::Panic, "the JSON is not representable as C text");
        }
        // SAFETY: `out` is non-null and writable by the caller's contract.
        unsafe { *out = pointer };
        RssnStatus::Ok
    })
}

/// Reads one live operator key without parsing JSON.
///
/// `name` may be a registered name or alias. `index` is the zero-based offset
/// in that operator's manifest `keys` array. Free the two successful string
/// results with [`rssn_string_free`](super::rssn_string_free). Outputs are
/// written only on success.
///
/// # Safety
/// `name` must be null or a readable NUL-terminated string. The output
/// pointers must be null or valid for writing two distinct string pointers
/// and 16 hash bytes, respectively. Output storage must not overlap.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_op_key_at(
    name: *const c_char,
    index: usize,
    out_rule_set: *mut *mut c_char,
    out_name: *mut *mut c_char,
    out_signature_hash: *mut u8,
) -> RssnStatus {
    super::guard(|| {
        if out_rule_set.is_null() || out_name.is_null() || out_signature_hash.is_null() {
            return super::fail(RssnStatus::NullArgument, "null operator key output pointer");
        }
        if std::ptr::eq(out_rule_set, out_name) {
            return super::fail(
                RssnStatus::InvalidArgument,
                "operator key outputs must be distinct",
            );
        }
        // SAFETY: the caller supplies a readable NUL-terminated string or null.
        let name = match unsafe { super::text(name) } {
            | Ok(name) => name,
            | Err(status) => return status,
        };
        let operator = match lookup_operator(name) {
            | Ok(operator) => operator,
            | Err(status) => return status,
        };
        let Some(key) = operator.keys.get(index) else {
            return super::fail(
                RssnStatus::InvalidArgument,
                &format!(
                    "operator {} has no signature key at index {index}",
                    operator.name
                ),
            );
        };
        let (Ok(rule_set), Ok(canonical_name)) = (
            CString::new(key.rule_set.as_str()),
            CString::new(key.name.as_str()),
        ) else {
            return super::fail(
                RssnStatus::Panic,
                "operator key is not representable as C text",
            );
        };
        // SAFETY: all three non-null outputs are writable by the caller's contract.
        unsafe {
            *out_rule_set = rule_set.into_raw();
            *out_name = canonical_name.into_raw();
            out_signature_hash
                .cast::<[u8; 16]>()
                .write(key.signature_hash);
        }
        RssnStatus::Ok
    })
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
        assert_eq!(
            manifest["api_version"].as_u64(),
            Some(u64::from(crate::kernel::API_VERSION))
        );
        assert!(text.contains("\"ops\""));
        assert_eq!(rssn_manifest_hash(), standard().hash_u64());
    }

    #[test]
    fn the_hash_is_stable_within_a_process() {
        assert_eq!(rssn_manifest_hash(), rssn_manifest_hash());
    }
}
