use std::os::raw::c_char;

use crate::ffi_apis::common::*;
use crate::ffi_apis::symbolic_graph_operations_ffi::handle::convert_expr_graph_to_string_graph;
use crate::symbolic::graph::Graph;
use crate::symbolic::graph_operations::*;

/// Creates an induced subgraph.
/// Input JSON: {"graph": <graph>, "nodes": ["label1", "label2"]}
#[no_mangle]

pub extern "C" fn rssn_json_graph_induced_subgraph(
    json : *const c_char
) -> *mut c_char {

    #[derive(serde::Deserialize)]

    struct Input {
        graph : Graph<String>,
        nodes : Vec<String>,
    }

    let input : Input = match from_json_string(json) {
        | Some(i) => i,
        | None => return std::ptr::null_mut(),
    };

    let result = induced_subgraph(
        &input.graph,
        &input.nodes,
    );

    to_json_string(&result)
}

/// Computes the union of two graphs.
/// Input JSON: {"g1": <graph>, "g2": <graph>}
#[no_mangle]

pub extern "C" fn rssn_json_graph_union(
    json : *const c_char
) -> *mut c_char {

    #[derive(serde::Deserialize)]

    struct Input {
        g1 : Graph<String>,
        g2 : Graph<String>,
    }

    let input : Input = match from_json_string(json) {
        | Some(i) => i,
        | None => return std::ptr::null_mut(),
    };

    let result =
        union(&input.g1, &input.g2);

    to_json_string(&result)
}

/// Computes the intersection of two graphs.
#[no_mangle]

pub extern "C" fn rssn_json_graph_intersection(
    json : *const c_char
) -> *mut c_char {

    #[derive(serde::Deserialize)]

    struct Input {
        g1 : Graph<String>,
        g2 : Graph<String>,
    }

    let input : Input = match from_json_string(json) {
        | Some(i) => i,
        | None => return std::ptr::null_mut(),
    };

    let result = intersection(
        &input.g1,
        &input.g2,
    );

    to_json_string(&result)
}

/// Computes the Cartesian product of two graphs.
#[no_mangle]

pub extern "C" fn rssn_json_graph_cartesian_product(
    json : *const c_char
) -> *mut c_char {

    #[derive(serde::Deserialize)]

    struct Input {
        g1 : Graph<String>,
        g2 : Graph<String>,
    }

    let input : Input = match from_json_string(json) {
        | Some(i) => i,
        | None => return std::ptr::null_mut(),
    };

    let result_expr = cartesian_product(
        &input.g1,
        &input.g2,
    );

    let result = convert_expr_graph_to_string_graph(result_expr);

    to_json_string(&result)
}

/// Computes the Tensor product of two graphs.
#[no_mangle]

pub extern "C" fn rssn_json_graph_tensor_product(
    json : *const c_char
) -> *mut c_char {

    #[derive(serde::Deserialize)]

    struct Input {
        g1 : Graph<String>,
        g2 : Graph<String>,
    }

    let input : Input = match from_json_string(json) {
        | Some(i) => i,
        | None => return std::ptr::null_mut(),
    };

    let result_expr = tensor_product(
        &input.g1,
        &input.g2,
    );

    let result = convert_expr_graph_to_string_graph(result_expr);

    to_json_string(&result)
}

/// Computes the complement of a graph.
#[no_mangle]

pub extern "C" fn rssn_json_graph_complement(
    json : *const c_char
) -> *mut c_char {

    let graph : Graph<String> = match from_json_string(json) {
        | Some(g) => g,
        | None => return std::ptr::null_mut(),
    };

    let result = complement(&graph);

    to_json_string(&result)
}

/// Computes the disjoint union of two graphs.
#[no_mangle]

pub extern "C" fn rssn_json_graph_disjoint_union(
    json : *const c_char
) -> *mut c_char {

    #[derive(serde::Deserialize)]

    struct Input {
        g1 : Graph<String>,
        g2 : Graph<String>,
    }

    let input : Input = match from_json_string(json) {
        | Some(i) => i,
        | None => return std::ptr::null_mut(),
    };

    let result_expr = disjoint_union(
        &input.g1,
        &input.g2,
    );

    let result = convert_expr_graph_to_string_graph(result_expr);

    to_json_string(&result)
}

/// Computes the join of two graphs.
#[no_mangle]

pub extern "C" fn rssn_json_graph_join(
    json : *const c_char
) -> *mut c_char {

    #[derive(serde::Deserialize)]

    struct Input {
        g1 : Graph<String>,
        g2 : Graph<String>,
    }

    let input : Input = match from_json_string(json) {
        | Some(i) => i,
        | None => return std::ptr::null_mut(),
    };

    let result_expr =
        join(&input.g1, &input.g2);

    let result = convert_expr_graph_to_string_graph(result_expr);

    to_json_string(&result)
}
