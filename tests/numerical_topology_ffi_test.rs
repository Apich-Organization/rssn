use rssn::ffi_apis::common::{rssn_free_bincode_buffer, rssn_free_string};
use rssn::ffi_apis::numerical_topology_ffi::{bincode_api, handle, json};
use rssn::numerical::graph::Graph;
use std::ffi::{CStr, CString};

#[test]

fn test_topology_handle_ffi() {

    unsafe {

        let mut g = Graph::new(3);

        g.add_edge(0, 1, 1.0);

        let ptr = handle::rssn_num_topology_find_connected_components(&g);

        assert!(!ptr.is_null());

        let comps = &*ptr;

        assert_eq!(comps.len(), 2);

        let _ = Box::from_raw(ptr);

        // Euclidean distance
        let p1 = [0.0, 0.0];

        let p2 = [3.0, 4.0];

        let dist = handle::rssn_num_topology_euclidean_distance(p1.as_ptr(), p2.as_ptr(), 2);

        assert!((dist - 5.0).abs() < 1e-9);
    }
}

#[test]

fn test_topology_json_ffi() {

    unsafe {

        let json_input = r#"{
            "points": [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
            "epsilon": 1.1,
            "max_dim": 1
        }"#;

        let c_json = CString::new(json_input).unwrap();

        let res_ptr = json::rssn_num_topology_betti_numbers_json(c_json.as_ptr());

        assert!(!res_ptr.is_null());

        let res_str = CStr::from_ptr(res_ptr).to_str().unwrap();

        let v: serde_json::Value = serde_json::from_str(res_str).unwrap();

        let betti = v["ok"].as_array().unwrap();

        assert_eq!(betti[0].as_u64().unwrap(), 1);

        assert_eq!(betti[1].as_u64().unwrap(), 1);

        rssn_free_string(res_ptr);
    }
}

#[test]

fn test_topology_bincode_ffi() {

    unsafe {

        use rssn::ffi_apis::common::{from_bincode_buffer, to_bincode_buffer};
        use serde::{Deserialize, Serialize};

        #[derive(Serialize)]

        struct BettiInput {
            points: Vec<Vec<f64>>,
            epsilon: f64,
            max_dim: usize,
        }

        let input = BettiInput {
            points: vec![
                vec![0.0, 0.0],
                vec![1.0, 0.0],
            ],
            epsilon: 1.5,
            max_dim: 1,
        };

        let buffer = to_bincode_buffer(&input);

        let res_buffer = bincode_api::rssn_num_topology_betti_numbers_bincode(buffer);

        assert!(!res_buffer.is_null());

        #[derive(Deserialize)]

        struct FfiResult<T, E> {
            ok: Option<T>,
            #[allow(dead_code)]
            err: Option<E>,
        }

        let res: FfiResult<Vec<usize>, String> = from_bincode_buffer(&res_buffer).unwrap();

        assert_eq!(res.ok.unwrap()[0], 1);

        rssn_free_bincode_buffer(res_buffer);

        rssn_free_bincode_buffer(buffer);
    }
}
