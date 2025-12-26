use std::ffi::CStr;

use rssn::ffi_apis::common::*;
use rssn::ffi_apis::constant_ffi::bincode_api::*;
use rssn::ffi_apis::constant_ffi::handle::*;
use rssn::ffi_apis::constant_ffi::json::*;

#[test]

fn test_handle_api_build_date() {

    unsafe {

        let ptr = rssn_get_build_date();

        assert!(!ptr.is_null());

        let c_str = CStr::from_ptr(ptr);

        let date = c_str
            .to_str()
            .unwrap();

        assert!(!date.is_empty());

        rssn_free_string(ptr);
    }
}

#[test]

fn test_handle_api_commit_sha() {

    unsafe {

        let ptr = rssn_get_commit_sha();

        assert!(!ptr.is_null());

        let c_str = CStr::from_ptr(ptr);

        let sha = c_str
            .to_str()
            .unwrap();

        assert!(!sha.is_empty());

        rssn_free_string(ptr);
    }
}

#[test]

fn test_json_api_build_info() {

    unsafe {

        let ptr = rssn_get_build_info_json();

        assert!(!ptr.is_null());

        let c_str = CStr::from_ptr(ptr);

        let json = c_str
            .to_str()
            .unwrap();

        // Parse JSON
        let info : serde_json::Value = serde_json::from_str(json).unwrap();

        assert!(info["build_date"].is_string());

        assert!(info["commit_sha"].is_string());

        assert!(info["rustc_version"].is_string());

        rssn_free_string(ptr);
    }
}

#[test]

fn test_json_api_build_date() {

    unsafe {

        let ptr = rssn_get_build_date_json();

        assert!(!ptr.is_null());

        let c_str = CStr::from_ptr(ptr);

        let json = c_str
            .to_str()
            .unwrap();

        // Should be a JSON string
        let date : String = serde_json::from_str(json).unwrap();

        assert!(!date.is_empty());

        rssn_free_string(ptr);
    }
}

#[test]

fn test_bincode_api_build_info() {

    let buffer = rssn_get_build_info_bincode();

    assert!(!buffer.is_null());

    unsafe {

        let slice = buffer.as_slice();

        let (info, _) : (BuildInfo, usize) = bincode::serde::decode_from_slice(
            slice,
            bincode::config::standard(),
        )
        .unwrap();

        assert!(!info
            .build_date
            .is_empty());

        assert!(!info
            .commit_sha
            .is_empty());

        assert!(!info
            .rustc_version
            .is_empty());
    }

    rssn_free_bincode_buffer(buffer);
}

#[test]

fn test_bincode_api_build_date() {

    let buffer = rssn_get_build_date_bincode();

    assert!(!buffer.is_null());

    unsafe {

        let slice = buffer.as_slice();

        let (date, _) : (String, usize) = bincode::serde::decode_from_slice(
            slice,
            bincode::config::standard(),
        )
        .unwrap();

        assert!(!date.is_empty());
    }

    rssn_free_bincode_buffer(buffer);
}

#[test]

fn test_all_three_apis_consistency() {

    // Get build date from all three APIs
    let handle_date = unsafe {

        let ptr = rssn_get_build_date();

        let c_str = CStr::from_ptr(ptr);

        let date = c_str
            .to_str()
            .unwrap()
            .to_string();

        rssn_free_string(ptr);

        date
    };

    let json_date = unsafe {

        let ptr = rssn_get_build_date_json();

        let c_str = CStr::from_ptr(ptr);

        let json = c_str
            .to_str()
            .unwrap();

        let date : String = serde_json::from_str(json).unwrap();

        rssn_free_string(ptr);

        date
    };

    let bincode_date = {

        let buffer = rssn_get_build_date_bincode();

        let date = unsafe {

            let slice = buffer.as_slice();

            let (d, _) : (String, usize) = bincode::serde::decode_from_slice(
                slice,
                bincode::config::standard(),
            )
            .unwrap();

            d
        };

        rssn_free_bincode_buffer(buffer);

        date
    };

    // All three should return the same value
    assert_eq!(
        handle_date,
        json_date
    );

    assert_eq!(
        json_date,
        bincode_date
    );
}
