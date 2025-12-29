use std::os::raw::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_json_string;
use crate::symbolic::core::Expr;
use crate::symbolic::graph::Graph;
use crate::symbolic::graph_algorithms::bfs;
use crate::symbolic::graph_algorithms::connected_components;
use crate::symbolic::graph_algorithms::dfs;
use crate::symbolic::graph_algorithms::edmonds_karp_max_flow;
use crate::symbolic::graph_algorithms::has_cycle;
use crate::symbolic::graph_algorithms::is_bipartite;
use crate::symbolic::graph_algorithms::kruskal_mst;

/// Creates a new graph from JSON specification.
/// JSON format: {"`is_directed"`: true/false}
#[no_mangle]

pub extern "C" fn rssn_json_graph_new(
    json: *const c_char
) -> *mut c_char {

    #[derive(serde::Deserialize)]

    struct GraphSpec {
        is_directed: bool,
    }

    let spec : GraphSpec = match from_json_string(json) {
        | Some(s) => s,
        | None => return std::ptr::null_mut(),
    };

    let graph: Graph<String> =
        Graph::new(spec.is_directed);

    to_json_string(&graph)
}

/// Adds a node to the graph.
/// Input JSON: {"graph": <graph>, "label": "`node_label`"}
/// Returns updated graph as JSON.
#[no_mangle]

pub extern "C" fn rssn_json_graph_add_node(
    json: *const c_char
) -> *mut c_char {

    #[derive(serde::Deserialize)]

    struct Input {
        graph: Graph<String>,
        label: String,
    }

    let mut input : Input = match from_json_string(json) {
        | Some(i) => i,
        | None => return std::ptr::null_mut(),
    };

    input
        .graph
        .add_node(input.label);

    to_json_string(&input.graph)
}

/// Adds an edge to the graph.
/// Input JSON: {"graph": <graph>, "from": "label1", "to": "label2", "weight": <expr>}
#[no_mangle]

pub extern "C" fn rssn_json_graph_add_edge(
    json: *const c_char
) -> *mut c_char {

    #[derive(serde::Deserialize)]

    struct Input {
        graph: Graph<String>,
        from: String,
        to: String,
        weight: Expr,
    }

    let mut input : Input = match from_json_string(json) {
        | Some(i) => i,
        | None => return std::ptr::null_mut(),
    };

    input
        .graph
        .add_edge(
            &input.from,
            &input.to,
            input.weight,
        );

    to_json_string(&input.graph)
}

/// Gets the adjacency matrix of the graph.
/// Input JSON: <graph>
/// Returns Expr (matrix) as JSON.
#[no_mangle]

pub extern "C" fn rssn_json_graph_adjacency_matrix(
    json: *const c_char
) -> *mut c_char {

    let graph : Graph<String> = match from_json_string(json) {
        | Some(g) => g,
        | None => return std::ptr::null_mut(),
    };

    let matrix =
        graph.to_adjacency_matrix();

    to_json_string(&matrix)
}

/// Gets the Laplacian matrix of the graph.
#[no_mangle]

pub extern "C" fn rssn_json_graph_laplacian_matrix(
    json: *const c_char
) -> *mut c_char {

    let graph : Graph<String> = match from_json_string(json) {
        | Some(g) => g,
        | None => return std::ptr::null_mut(),
    };

    let matrix =
        graph.to_laplacian_matrix();

    to_json_string(&matrix)
}

/// Performs BFS traversal.
/// Input JSON: {"graph": <graph>, "`start_node"`: <index>}
#[no_mangle]

pub extern "C" fn rssn_json_graph_bfs(
    json: *const c_char
) -> *mut c_char {

    #[derive(serde::Deserialize)]

    struct Input {
        graph: Graph<String>,
        start_node: usize,
    }

    let input : Input = match from_json_string(json) {
        | Some(i) => i,
        | None => return std::ptr::null_mut(),
    };

    let result = bfs(
        &input.graph,
        input.start_node,
    );

    to_json_string(&result)
}

/// Performs DFS traversal.
#[no_mangle]

pub extern "C" fn rssn_json_graph_dfs(
    json: *const c_char
) -> *mut c_char {

    #[derive(serde::Deserialize)]

    struct Input {
        graph: Graph<String>,
        start_node: usize,
    }

    let input : Input = match from_json_string(json) {
        | Some(i) => i,
        | None => return std::ptr::null_mut(),
    };

    let result = dfs(
        &input.graph,
        input.start_node,
    );

    to_json_string(&result)
}

/// Finds connected components.
#[no_mangle]

pub extern "C" fn rssn_json_graph_connected_components(
    json: *const c_char
) -> *mut c_char {

    let graph : Graph<String> = match from_json_string(json) {
        | Some(g) => g,
        | None => return std::ptr::null_mut(),
    };

    let result =
        connected_components(&graph);

    to_json_string(&result)
}

/// Computes maximum flow.
/// Input JSON: {"graph": <graph>, "source": <index>, "sink": <index>}
#[no_mangle]

pub extern "C" fn rssn_json_graph_max_flow(
    json: *const c_char
) -> *mut c_char {

    #[derive(serde::Deserialize)]

    struct Input {
        graph: Graph<String>,
        source: usize,
        sink: usize,
    }

    let input : Input = match from_json_string(json) {
        | Some(i) => i,
        | None => return std::ptr::null_mut(),
    };

    let flow = edmonds_karp_max_flow(
        &input.graph,
        input.source,
        input.sink,
    );

    to_json_string(&flow)
}

/// Computes MST using Kruskal's algorithm.
#[no_mangle]

pub extern "C" fn rssn_json_graph_kruskal_mst(
    json: *const c_char
) -> *mut c_char {

    let graph : Graph<String> = match from_json_string(json) {
        | Some(g) => g,
        | None => return std::ptr::null_mut(),
    };

    let mst = kruskal_mst(&graph);

    to_json_string(&mst)
}

/// Checks if graph has a cycle.
#[no_mangle]

pub extern "C" fn rssn_json_graph_has_cycle(
    json: *const c_char
) -> *mut c_char {

    let graph : Graph<String> = match from_json_string(json) {
        | Some(g) => g,
        | None => return std::ptr::null_mut(),
    };

    let result = has_cycle(&graph);

    to_json_string(&result)
}

/// Checks if graph is bipartite.
#[no_mangle]

pub extern "C" fn rssn_json_graph_is_bipartite(
    json: *const c_char
) -> *mut c_char {

    let graph : Graph<String> = match from_json_string(json) {
        | Some(g) => g,
        | None => return std::ptr::null_mut(),
    };

    let result = is_bipartite(&graph);

    to_json_string(&result.is_some())
}
