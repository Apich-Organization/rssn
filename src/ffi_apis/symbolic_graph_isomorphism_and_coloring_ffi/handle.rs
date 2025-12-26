use crate::ffi_apis::symbolic_graph_ffi::handle::RssnGraph;
use crate::symbolic::graph::Graph;
use crate::symbolic::graph_isomorphism_and_coloring::*;
use std::os::raw::{
    c_char,
    c_int,
};

/// Checks if two graphs are potentially isomorphic using WL test.
#[no_mangle]

pub extern "C" fn rssn_are_isomorphic_heuristic(
    g1 : *const RssnGraph,
    g2 : *const RssnGraph,
) -> c_int {

    if g1.is_null() || g2.is_null() {

        return 0;
    }

    unsafe {

        let graph1 = &*(g1
            as *const Graph<String>);

        let graph2 = &*(g2
            as *const Graph<String>);

        if are_isomorphic_heuristic(
            graph1,
            graph2,
        ) {

            1
        } else {

            0
        }
    }
}

/// Finds a valid vertex coloring using greedy heuristic.
/// Returns a JSON object mapping node IDs to colors.
#[no_mangle]

pub extern "C" fn rssn_greedy_coloring(
    graph : *const RssnGraph
) -> *mut c_char {

    if graph.is_null() {

        return std::ptr::null_mut();
    }

    unsafe {

        let g = &*(graph
            as *const Graph<String>);

        let colors = greedy_coloring(g);

        let json =
            serde_json::to_string(
                &colors,
            )
            .unwrap_or_default();

        std::ffi::CString::new(json)
            .unwrap()
            .into_raw()
    }
}

/// Finds the chromatic number exactly (NP-hard).
#[no_mangle]

pub extern "C" fn rssn_chromatic_number_exact(
    graph : *const RssnGraph
) -> usize {

    if graph.is_null() {

        return 0;
    }

    unsafe {

        let g = &*(graph
            as *const Graph<String>);

        chromatic_number_exact(g)
    }
}
