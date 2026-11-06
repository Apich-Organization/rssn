//! The source overlay must survive embedding with every reviewed overload intact.
use std::collections::BTreeMap;

use rssn::kernel::Signature;
use rssn::kernel::SignatureSource;
use serde::Deserialize;

#[derive(Deserialize)]
struct Overlay {
    ops: BTreeMap<String, Vec<Signature>>,
}

#[test]
fn source_overlay_matches_the_embedded_manifest() {
    let source: Overlay = serde_json::from_str(include_str!("../signatures.overlay.json"))
        .expect("source overlay must contain valid signatures");
    assert!(
        !source.ops.is_empty(),
        "the reviewed overlay must not disappear"
    );
    let session = rssn::Session::new();
    let manifest = session.manifest();
    for (name, signatures) in source.ops {
        assert!(!signatures.is_empty(), "empty overlay for {name}");
        let op = manifest
            .ops
            .iter()
            .find(|op| op.name == name)
            .unwrap_or_else(|| panic!("overlay operator {name} is not registered"));
        assert_eq!(op.signature_source, SignatureSource::Overlay, "{name}");
        assert_eq!(op.overloads, signatures, "{name}");
    }
}
