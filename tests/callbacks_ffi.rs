//! C boundary checks stay separate from the safe host callback contract suite.
#[test]
fn api_v2_manifest_crosses_the_c_string_boundary() {
    let pointer = rssn::ffi::rssn_manifest_json();
    assert!(!pointer.is_null());
    // SAFETY: this pointer is the live owned CString returned by RSSN.
    let text = unsafe { std::ffi::CStr::from_ptr(pointer) }
        .to_str()
        .unwrap()
        .to_owned();
    // SAFETY: the returned CString is released exactly once.
    unsafe { rssn::ffi::rssn_string_free(pointer) };
    let manifest: serde_json::Value = serde_json::from_str(&text).unwrap();
    assert_eq!(manifest["api_version"].as_u64(), Some(2));
    assert!(manifest["ops"].as_array().unwrap().iter().any(|op| {
        op["name"] == "num_integrate"
            && op["rust_path"] == "crate::kernels::quadrature::integrate_gk"
    }));
}
