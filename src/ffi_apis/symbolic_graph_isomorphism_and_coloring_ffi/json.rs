use crate::ffi_apis::common::*;
use crate::symbolic::graph::Graph;
use crate::symbolic::graph_isomorphism_and_coloring::*;

/// Checks if two graphs are isomorphic.
/// Input: {"g1": Graph, "g2": Graph}
/// Output: bool
#[no_mangle]

pub unsafe extern "C" fn rssn_json_are_isomorphic_heuristic(
    json: *const std::os::raw::c_char
) -> *mut std::os::raw::c_char {

    #[derive(serde::Deserialize)]

    struct Input {
        g1: Graph<String>,
        g2: Graph<String>,
    }

    let input: Input = match from_json_string(json) {
        Some(i) => i,
        None => return std::ptr::null_mut(),
    };

    let result = are_isomorphic_heuristic(&input.g1, &input.g2);

    to_json_string(&result)
}

/// Greedy coloring.
/// Input: Graph
/// Output: {node_id: color_id}
#[no_mangle]

pub unsafe extern "C" fn rssn_json_greedy_coloring(
    json: *const std::os::raw::c_char
) -> *mut std::os::raw::c_char {

    let graph: Graph<String> = match from_json_string(json) {
        Some(g) => g,
        None => return std::ptr::null_mut(),
    };

    let result = greedy_coloring(&graph);

    to_json_string(&result)
}

/// Exact chromatic number.
/// Input: Graph
/// Output: usize
#[no_mangle]

pub unsafe extern "C" fn rssn_json_chromatic_number_exact(
    json: *const std::os::raw::c_char
) -> *mut std::os::raw::c_char {

    let graph: Graph<String> = match from_json_string(json) {
        Some(g) => g,
        None => return std::ptr::null_mut(),
    };

    let result = chromatic_number_exact(&graph);

    to_json_string(&result)
}
