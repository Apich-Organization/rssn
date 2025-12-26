use crate::symbolic::core::Expr;
use crate::symbolic::graph::Graph;
use crate::symbolic::graph_algorithms::*;
use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int};

/// Opaque type for Graph<String> to work with cbindgen
#[repr(C)]

pub struct RssnGraph {
    _private: [u8; 0],
}

/// Creates a new graph.
#[no_mangle]

pub extern "C" fn rssn_graph_new(is_directed: c_int) -> *mut RssnGraph {

    let graph = Graph::<String>::new(is_directed != 0);

    Box::into_raw(Box::new(graph)) as *mut RssnGraph
}

/// Frees a graph.
#[no_mangle]

pub extern "C" fn rssn_graph_free(ptr: *mut RssnGraph) {

    if !ptr.is_null() {

        unsafe {

            drop(Box::from_raw(ptr as *mut Graph<String>))
        };
    }
}

/// Adds a node to the graph.
#[no_mangle]

pub extern "C" fn rssn_graph_add_node(ptr: *mut RssnGraph, label: *const c_char) -> usize {

    if ptr.is_null() || label.is_null() {

        return usize::MAX;
    }

    let graph = unsafe {

        &mut *(ptr as *mut Graph<String>)
    };

    let label_str = unsafe {

        CStr::from_ptr(label)
            .to_string_lossy()
            .into_owned()
    };

    graph.add_node(label_str)
}

/// Adds an edge to the graph.
#[no_mangle]

pub extern "C" fn rssn_graph_add_edge(
    ptr: *mut RssnGraph,
    from_label: *const c_char,
    to_label: *const c_char,
    weight: *const Expr,
) {

    if ptr.is_null() || from_label.is_null() || to_label.is_null() || weight.is_null() {

        return;
    }

    let graph = unsafe {

        &mut *(ptr as *mut Graph<String>)
    };

    let from = unsafe {

        CStr::from_ptr(from_label)
            .to_string_lossy()
            .into_owned()
    };

    let to = unsafe {

        CStr::from_ptr(to_label)
            .to_string_lossy()
            .into_owned()
    };

    let weight_expr = unsafe {

        (*weight).clone()
    };

    graph.add_edge(&from, &to, weight_expr);
}

/// Gets the number of nodes in the graph.
#[no_mangle]

pub extern "C" fn rssn_graph_node_count(ptr: *const RssnGraph) -> usize {

    if ptr.is_null() {

        return 0;
    }

    let graph = unsafe {

        &*(ptr as *const Graph<String>)
    };

    graph.node_count()
}

/// Gets the adjacency matrix of the graph.
#[no_mangle]

pub extern "C" fn rssn_graph_adjacency_matrix(ptr: *const RssnGraph) -> *mut Expr {

    if ptr.is_null() {

        return std::ptr::null_mut();
    }

    let graph = unsafe {

        &*(ptr as *const Graph<String>)
    };

    Box::into_raw(Box::new(graph.to_adjacency_matrix()))
}

/// Gets the incidence matrix of the graph.
#[no_mangle]

pub extern "C" fn rssn_graph_incidence_matrix(ptr: *const RssnGraph) -> *mut Expr {

    if ptr.is_null() {

        return std::ptr::null_mut();
    }

    let graph = unsafe {

        &*(ptr as *const Graph<String>)
    };

    Box::into_raw(Box::new(graph.to_incidence_matrix()))
}

/// Gets the Laplacian matrix of the graph.
#[no_mangle]

pub extern "C" fn rssn_graph_laplacian_matrix(ptr: *const RssnGraph) -> *mut Expr {

    if ptr.is_null() {

        return std::ptr::null_mut();
    }

    let graph = unsafe {

        &*(ptr as *const Graph<String>)
    };

    Box::into_raw(Box::new(graph.to_laplacian_matrix()))
}

/// Performs BFS traversal from a start node.
/// Returns a JSON string containing the node IDs in visit order.
#[no_mangle]

pub extern "C" fn rssn_graph_bfs(ptr: *const RssnGraph, start_node: usize) -> *mut c_char {

    if ptr.is_null() {

        return std::ptr::null_mut();
    }

    let graph = unsafe {

        &*(ptr as *const Graph<String>)
    };

    let result = bfs(graph, start_node);

    match serde_json::to_string(&result) {
        Ok(json) => CString::new(json)
            .unwrap()
            .into_raw(),
        Err(_) => std::ptr::null_mut(),
    }
}

/// Performs DFS traversal from a start node.
/// Returns a JSON string containing the node IDs in visit order.
#[no_mangle]

pub extern "C" fn rssn_graph_dfs(ptr: *const RssnGraph, start_node: usize) -> *mut c_char {

    if ptr.is_null() {

        return std::ptr::null_mut();
    }

    let graph = unsafe {

        &*(ptr as *const Graph<String>)
    };

    let result = dfs(graph, start_node);

    match serde_json::to_string(&result) {
        Ok(json) => CString::new(json)
            .unwrap()
            .into_raw(),
        Err(_) => std::ptr::null_mut(),
    }
}

/// Finds connected components.
/// Returns a JSON string containing the components.
#[no_mangle]

pub extern "C" fn rssn_graph_connected_components(ptr: *const RssnGraph) -> *mut c_char {

    if ptr.is_null() {

        return std::ptr::null_mut();
    }

    let graph = unsafe {

        &*(ptr as *const Graph<String>)
    };

    let result = connected_components(graph);

    match serde_json::to_string(&result) {
        Ok(json) => CString::new(json)
            .unwrap()
            .into_raw(),
        Err(_) => std::ptr::null_mut(),
    }
}

/// Computes maximum flow using Edmonds-Karp algorithm.
#[no_mangle]

pub extern "C" fn rssn_graph_max_flow(ptr: *const RssnGraph, source: usize, sink: usize) -> f64 {

    if ptr.is_null() {

        return 0.0;
    }

    let graph = unsafe {

        &*(ptr as *const Graph<String>)
    };

    edmonds_karp_max_flow(graph, source, sink)
}

/// Computes minimum spanning tree using Kruskal's algorithm.
/// Returns a JSON string containing the MST edges.
#[no_mangle]

pub extern "C" fn rssn_graph_kruskal_mst(ptr: *const RssnGraph) -> *mut c_char {

    if ptr.is_null() {

        return std::ptr::null_mut();
    }

    let graph = unsafe {

        &*(ptr as *const Graph<String>)
    };

    let mst = kruskal_mst(graph);

    // Convert to a simpler format for JSON
    let edges: Vec<(usize, usize, String)> = mst
        .iter()
        .map(|(u, v, w)| (*u, *v, format!("{}", w)))
        .collect();

    match serde_json::to_string(&edges) {
        Ok(json) => CString::new(json)
            .unwrap()
            .into_raw(),
        Err(_) => std::ptr::null_mut(),
    }
}

/// Checks if the graph has a cycle.
#[no_mangle]

pub extern "C" fn rssn_graph_has_cycle(ptr: *const RssnGraph) -> c_int {

    if ptr.is_null() {

        return 0;
    }

    let graph = unsafe {

        &*(ptr as *const Graph<String>)
    };

    if has_cycle(graph) {

        1
    } else {

        0
    }
}

/// Checks if the graph is bipartite.
/// Returns 1 if bipartite, 0 otherwise.
#[no_mangle]

pub extern "C" fn rssn_graph_is_bipartite(ptr: *const RssnGraph) -> c_int {

    if ptr.is_null() {

        return 0;
    }

    let graph = unsafe {

        &*(ptr as *const Graph<String>)
    };

    if is_bipartite(graph).is_some() {

        1
    } else {

        0
    }
}
