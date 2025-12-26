use crate::ffi_apis::symbolic_graph_ffi::handle::RssnGraph;
use crate::symbolic::graph::Graph;
use crate::symbolic::graph_algorithms::*;
use std::ffi::CStr;
use std::os::raw::{
    c_char,
    c_int,
};

/// Performs DFS traversal starting from a given node.
/// Returns a JSON array of node indices in visit order.
#[no_mangle]

pub extern "C" fn rssn_graph_dfs_api(
    graph : *const RssnGraph,
    start_node : usize,
) -> *mut c_char {

    if graph.is_null() {

        return std::ptr::null_mut();
    }

    unsafe {

        let g = &*(graph
            as *const Graph<String>);

        let result = dfs(g, start_node);

        let json =
            serde_json::to_string(
                &result,
            )
            .unwrap_or_default();

        std::ffi::CString::new(json)
            .unwrap()
            .into_raw()
    }
}

/// Performs BFS traversal starting from a given node.
/// Returns a JSON array of node indices in visit order.
#[no_mangle]

pub extern "C" fn rssn_graph_bfs_api(
    graph : *const RssnGraph,
    start_node : usize,
) -> *mut c_char {

    if graph.is_null() {

        return std::ptr::null_mut();
    }

    unsafe {

        let g = &*(graph
            as *const Graph<String>);

        let result = bfs(g, start_node);

        let json =
            serde_json::to_string(
                &result,
            )
            .unwrap_or_default();

        std::ffi::CString::new(json)
            .unwrap()
            .into_raw()
    }
}

/// Finds all connected components in an undirected graph.
/// Returns a JSON array of arrays, where each inner array is a component.
#[no_mangle]

pub extern "C" fn rssn_graph_connected_components_api(
    graph : *const RssnGraph
) -> *mut c_char {

    if graph.is_null() {

        return std::ptr::null_mut();
    }

    unsafe {

        let g = &*(graph
            as *const Graph<String>);

        let result =
            connected_components(g);

        let json =
            serde_json::to_string(
                &result,
            )
            .unwrap_or_default();

        std::ffi::CString::new(json)
            .unwrap()
            .into_raw()
    }
}

/// Checks if the graph is connected.
#[no_mangle]

pub extern "C" fn rssn_graph_is_connected(
    graph : *const RssnGraph
) -> c_int {

    if graph.is_null() {

        return 0;
    }

    unsafe {

        let g = &*(graph
            as *const Graph<String>);

        if is_connected(g) {

            1
        } else {

            0
        }
    }
}

/// Finds all strongly connected components in a directed graph.
/// Returns a JSON array of arrays.
#[no_mangle]

pub extern "C" fn rssn_graph_strongly_connected_components(
    graph : *const RssnGraph
) -> *mut c_char {

    if graph.is_null() {

        return std::ptr::null_mut();
    }

    unsafe {

        let g = &*(graph
            as *const Graph<String>);

        let result = strongly_connected_components(g);

        let json =
            serde_json::to_string(
                &result,
            )
            .unwrap_or_default();

        std::ffi::CString::new(json)
            .unwrap()
            .into_raw()
    }
}

/// Checks if the graph has a cycle.
#[no_mangle]

pub extern "C" fn rssn_graph_has_cycle_api(
    graph : *const RssnGraph
) -> c_int {

    if graph.is_null() {

        return 0;
    }

    unsafe {

        let g = &*(graph
            as *const Graph<String>);

        if has_cycle(g) {

            1
        } else {

            0
        }
    }
}

/// Finds bridges and articulation points.
/// Returns a JSON object with "bridges" and "articulation_points" fields.
#[no_mangle]

pub extern "C" fn rssn_graph_bridges_and_articulation_points_api(
    graph : *const RssnGraph
) -> *mut c_char {

    if graph.is_null() {

        return std::ptr::null_mut();
    }

    unsafe {

        let g = &*(graph
            as *const Graph<String>);

        let (bridges, aps) =
            find_bridges_and_articulation_points(g);

        #[derive(serde::Serialize)]

        struct Result {
            bridges :
                Vec<(usize, usize)>,
            articulation_points :
                Vec<usize>,
        }

        let result = Result {
            bridges,
            articulation_points : aps,
        };

        let json =
            serde_json::to_string(
                &result,
            )
            .unwrap_or_default();

        std::ffi::CString::new(json)
            .unwrap()
            .into_raw()
    }
}

/// Computes the minimum spanning tree using Kruskal's algorithm.
/// Returns a new graph containing only the MST edges.
#[no_mangle]

pub extern "C" fn rssn_graph_kruskal_mst_api(
    graph : *const RssnGraph
) -> *mut RssnGraph {

    if graph.is_null() {

        return std::ptr::null_mut();
    }

    unsafe {

        let g = &*(graph
            as *const Graph<String>);

        let mst_edges = kruskal_mst(g);

        // Create a new graph with the MST edges
        let mut mst_graph =
            Graph::new(g.is_directed);

        for node in &g.nodes {

            mst_graph
                .add_node(node.clone());
        }

        for (u, v, weight) in mst_edges
        {

            mst_graph.add_edge(
                &g.nodes[u],
                &g.nodes[v],
                weight,
            );
        }

        Box::into_raw(Box::new(
            mst_graph,
        )) as *mut RssnGraph
    }
}

/// Computes maximum flow using Edmonds-Karp algorithm.
#[no_mangle]

pub extern "C" fn rssn_graph_edmonds_karp_max_flow(
    graph : *const RssnGraph,
    source : usize,
    sink : usize,
) -> f64 {

    if graph.is_null() {

        return 0.0;
    }

    unsafe {

        let g = &*(graph
            as *const Graph<String>);

        edmonds_karp_max_flow(
            g,
            source,
            sink,
        )
    }
}

/// Computes maximum flow using Dinic's algorithm.
#[no_mangle]

pub extern "C" fn rssn_graph_dinic_max_flow(
    graph : *const RssnGraph,
    source : usize,
    sink : usize,
) -> f64 {

    if graph.is_null() {

        return 0.0;
    }

    unsafe {

        let g = &*(graph
            as *const Graph<String>);

        dinic_max_flow(g, source, sink)
    }
}

/// Checks if a graph is bipartite.
/// Returns a JSON array of partition assignments (0 or 1 for each node), or null if not bipartite.
#[no_mangle]

pub extern "C" fn rssn_graph_is_bipartite_api(
    graph : *const RssnGraph
) -> *mut c_char {

    if graph.is_null() {

        return std::ptr::null_mut();
    }

    unsafe {

        let g = &*(graph
            as *const Graph<String>);

        match is_bipartite(g) {
            | Some(partition) => {

                let json =
                    serde_json::to_string(&partition)
                        .unwrap_or_default();

                std::ffi::CString::new(
                    json,
                )
                .unwrap()
                .into_raw()
            },
            | None => {
                std::ptr::null_mut()
            },
        }
    }
}

/// Finds maximum matching in a bipartite graph.
/// partition_json should be a JSON array of 0s and 1s indicating the partition.
/// Returns a JSON array of [u, v] pairs representing the matching.
#[no_mangle]

pub extern "C" fn rssn_graph_bipartite_maximum_matching(
    graph : *const RssnGraph,
    partition_json : *const c_char,
) -> *mut c_char {

    if graph.is_null()
        || partition_json.is_null()
    {

        return std::ptr::null_mut();
    }

    unsafe {

        let g = &*(graph
            as *const Graph<String>);

        let partition_str =
            CStr::from_ptr(
                partition_json,
            )
            .to_str()
            .unwrap_or("");

        let partition: Vec<i8> =
            match serde_json::from_str(partition_str) {
                | Ok(p) => p,
                | Err(_) => return std::ptr::null_mut(),
            };

        let matching =
            bipartite_maximum_matching(
                g,
                &partition,
            );

        let json =
            serde_json::to_string(
                &matching,
            )
            .unwrap_or_default();

        std::ffi::CString::new(json)
            .unwrap()
            .into_raw()
    }
}

/// Performs topological sort on a DAG.
/// Returns a JSON array of node indices in topological order, or null if the graph has a cycle.
#[no_mangle]

pub extern "C" fn rssn_graph_topological_sort(
    graph : *const RssnGraph
) -> *mut c_char {

    if graph.is_null() {

        return std::ptr::null_mut();
    }

    unsafe {

        let g = &*(graph
            as *const Graph<String>);

        match topological_sort(g) {
            | Some(order) => {

                let json = serde_json::to_string(&order)
                    .unwrap_or_default();

                std::ffi::CString::new(
                    json,
                )
                .unwrap()
                .into_raw()
            },
            | None => {
                std::ptr::null_mut()
            },
        }
    }
}

/// Frees a C string returned by other functions.
#[no_mangle]

pub extern "C" fn rssn_free_string_api(
    ptr : *mut c_char
) {

    if !ptr.is_null() {

        unsafe {

            let _ = std::ffi::CString::from_raw(ptr);
        }
    }
}
