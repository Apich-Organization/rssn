//! C callers and ordinary Rust sessions share the canonical manifest.

use std::ffi::CStr;

use rssn::api::Session;
use rssn::ffi::rssn_manifest_hash;
use rssn::ffi::rssn_manifest_json;
use rssn::ffi::rssn_string_free;
use rssn::kernel::manifest::KernelManifest;

#[test]
fn c_manifest_keeps_the_legacy_bytes_and_has_independent_owned_copies() {
    // Rebuild through the original C path to catch a difference in the
    // standard Rust session's registry or unit metadata.
    let sets = rssn::rules::standard();
    let mut graph = rssn::graph::Graph::new();
    rssn::graph::Engine::install(&mut graph, &sets).expect("standard rules install");
    let legacy = KernelManifest::build(&graph, &sets, &[]);
    let expected = legacy.to_json();
    assert_eq!(rssn_manifest_hash(), legacy.hash_u64());
    assert_eq!(Session::new().manifest_json(), expected);

    std::thread::scope(|scope| {
        for _ in 0..4 {
            let expected = &expected;
            scope.spawn(move || {
                let first = rssn_manifest_json();
                let second = rssn_manifest_json();
                assert!(!first.is_null() && !second.is_null());
                assert_ne!(first, second);
                // SAFETY: Each function result is a live NUL-terminated string
                // owned by this caller and each is freed once.
                unsafe {
                    assert_eq!(CStr::from_ptr(first).to_bytes(), expected.as_bytes());
                    rssn_string_free(first);
                    assert_eq!(CStr::from_ptr(second).to_bytes(), expected.as_bytes());
                    rssn_string_free(second);
                }
            });
        }
    });
}
