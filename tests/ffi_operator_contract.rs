//! Operator discovery and overload selection over the C ABI.
use rssn::ffi::*;
use rssn::kernel::manifest::OpInfo;
use std::ffi::CStr;
use std::ffi::c_char;
use std::ptr;

struct Session(*mut RssnSession);

impl Session {
    fn new() -> Self {
        let raw = rssn_session_new();
        assert!(!raw.is_null());
        Self(raw)
    }
}

impl Drop for Session {
    fn drop(&mut self) {
        // SAFETY: the guard owns the session and frees it once.
        unsafe { rssn_session_free(self.0) };
    }
}

fn take_string(raw: *mut c_char) -> String {
    assert!(!raw.is_null());
    // SAFETY: all callers pass successful RSSN string results and transfer ownership.
    let text = unsafe { CStr::from_ptr(raw) }
        .to_string_lossy()
        .into_owned();
    // SAFETY: the owned RSSN allocation is no longer borrowed.
    unsafe { rssn_string_free(raw) };
    text
}

fn info(name: &CStr) -> OpInfo {
    let mut json = ptr::null_mut();
    // SAFETY: the C string is NUL-terminated and the output is writable.
    assert_eq!(
        unsafe { rssn_op_info_json(name.as_ptr(), &mut json) },
        RssnStatus::Ok
    );
    serde_json::from_str(&take_string(json)).expect("operator contract JSON")
}

#[test]
fn one_operator_contract_resolves_aliases_and_matches_the_manifest() {
    let alias = info(c"airy_ai");
    assert_eq!(alias.name, "airyai");
    assert_eq!(alias.keys, info(c"airyai").keys);
    let mesh = info(c"num_meshgen_tetra_volumes");
    let manifest: serde_json::Value =
        serde_json::from_str(&take_string(rssn_manifest_json())).unwrap();
    let entry = manifest["ops"]
        .as_array()
        .unwrap()
        .iter()
        .find(|entry| entry["name"].as_str() == Some(mesh.name.as_str()))
        .unwrap();
    assert_eq!(&serde_json::to_value(mesh).unwrap(), entry);
    let mut missing = ptr::null_mut();
    // SAFETY: the name and output follow the API's pointer contract.
    assert_eq!(
        unsafe {
            rssn_op_info_json(
                c"no_such_rssn_operation_for_contract_lookup".as_ptr(),
                &mut missing,
            )
        },
        RssnStatus::UnknownOperator
    );
    assert!(missing.is_null());
}

#[test]
fn pinned_calls_preserve_wide_cells_and_reject_stale_keys() {
    let session = Session::new();
    let a = c"79228162514264337593543950336";
    let b = c"237684487542793012780631851008";
    let operator = info(c"num_number_theory_gcd_u128");
    let key = operator.keys.first().expect("declared signature key");
    assert_eq!(key.rule_set, "numerical_special");
    assert_eq!(key.name, "num_number_theory_gcd_u128");
    let mut cells = [0, 0];
    // SAFETY: the session is owned, the two literals are C strings and the targets exist.
    unsafe {
        assert_eq!(
            rssn_value_new_int_decimal(session.0, a.as_ptr(), &mut cells[0]),
            RssnStatus::Ok
        );
        assert_eq!(
            rssn_value_new_int_decimal(session.0, b.as_ptr(), &mut cells[1]),
            RssnStatus::Ok
        );
    }
    let call = |hash: &[u8; 16], out: &mut u32| {
        // SAFETY: the session, key, argument ids, and output are live for this call.
        unsafe {
            rssn_value_call_key(
                session.0,
                c"numerical_special".as_ptr(),
                c"num_number_theory_gcd_u128".as_ptr(),
                hash.as_ptr(),
                cells.as_ptr(),
                cells.len(),
                out,
            )
        }
    };
    let mut result = 0;
    assert_eq!(call(&key.signature_hash, &mut result), RssnStatus::Ok);
    let mut text = ptr::null_mut();
    // SAFETY: the session owns the result cell and text points to writable storage.
    assert_eq!(
        unsafe { rssn_value_to_text(session.0, result, &mut text) },
        RssnStatus::Ok
    );
    assert_eq!(take_string(text), a.to_str().unwrap());

    let mut bad_hash = key.signature_hash;
    bad_hash[0] ^= 0x80;
    let mut untouched = 0;
    assert_eq!(call(&bad_hash, &mut untouched), RssnStatus::UnknownOperator);
    assert_eq!(untouched, 0);
    // SAFETY: the guard still owns this live session.
    assert_eq!(unsafe { rssn_value_count(session.0) }, 3);
}

#[test]
fn selected_key_matches_the_manifest_alias_and_rejects_out_of_range_indices() {
    let expected = info(c"airy_ai");
    let expected = expected.keys.first().expect("manifest key");
    let mut rule_set = ptr::null_mut();
    let mut name = ptr::null_mut();
    let mut hash = [0; 16];
    // SAFETY: the C name and the three distinct output targets follow the API contract.
    assert_eq!(
        unsafe {
            rssn_op_key_at(
                c"airy_ai".as_ptr(),
                0,
                &mut rule_set,
                &mut name,
                hash.as_mut_ptr(),
            )
        },
        RssnStatus::Ok
    );
    assert_eq!(take_string(rule_set), expected.rule_set);
    assert_eq!(take_string(name), expected.name);
    assert_eq!(hash, expected.signature_hash);

    let mut rule_set = ptr::null_mut();
    let mut name = ptr::null_mut();
    let mut hash = [0xaa; 16];
    // SAFETY: all outputs are valid even when the key index is invalid.
    assert_eq!(
        unsafe {
            rssn_op_key_at(
                c"airy_ai".as_ptr(),
                usize::MAX,
                &mut rule_set,
                &mut name,
                hash.as_mut_ptr(),
            )
        },
        RssnStatus::InvalidArgument
    );
    assert!(rule_set.is_null() && name.is_null());
    assert_eq!(hash, [0xaa; 16]);
}
