//! Handle-based FFI API for numerical graph algorithms.

use std::slice;

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::graph::bfs;
use crate::numerical::graph::dijkstra;
use crate::numerical::graph::floyd_warshall;
use crate::numerical::graph::page_rank;
use crate::numerical::graph::Graph;

/// Creates a new graph.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_graph_create(
    num_nodes: usize
) -> *mut Graph {

    Box::into_raw(Box::new(
        Graph::new(num_nodes),
    ))
}

/// Frees a graph.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_graph_free(
    graph: *mut Graph
) {

    if !graph.is_null() {

        let _ = Box::from_raw(graph);
    }
}

/// Adds a directed edge.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_graph_add_edge(
    graph: *mut Graph,
    u: usize,
    v: usize,
    weight: f64,
) {

    if let Some(g) = graph.as_mut() {

        g.add_edge(u, v, weight);
    }
}

/// Computes Dijkstra's algorithm.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_graph_dijkstra(
    graph: *mut Graph,
    start_node: usize,
    dist: *mut f64,
    prev: *mut isize,
) -> i32 {

    if graph.is_null()
        || dist.is_null()
        || prev.is_null()
    {

        update_last_error(
            "Null pointer passed to \
             rssn_num_graph_dijkstra"
                .to_string(),
        );

        return -1;
    }

    let g = &*graph;

    let n = g.num_nodes();

    let (d, p) =
        dijkstra(g, start_node);

    let dist_slice =
        slice::from_raw_parts_mut(
            dist, n,
        );

    let prev_slice =
        slice::from_raw_parts_mut(
            prev, n,
        );

    dist_slice.copy_from_slice(&d);

    for i in 0 .. n {

        prev_slice[i] = match p[i] {
            | Some(node) => {
                node as isize
            },
            | None => -1,
        };
    }

    0
}

/// Computes BFS.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_graph_bfs(
    graph: *mut Graph,
    start_node: usize,
    dist: *mut usize,
) -> i32 {

    if graph.is_null() || dist.is_null()
    {

        update_last_error(
            "Null pointer passed to \
             rssn_num_graph_bfs"
                .to_string(),
        );

        return -1;
    }

    let g = &*graph;

    let n = g.num_nodes();

    let d = bfs(g, start_node);

    let dist_slice =
        slice::from_raw_parts_mut(
            dist, n,
        );

    dist_slice.copy_from_slice(&d);

    0
}

/// Computes `PageRank`.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_graph_page_rank(
    graph: *mut Graph,
    damping_factor: f64,
    tolerance: f64,
    max_iter: usize,
    scores: *mut f64,
) -> i32 {

    if graph.is_null()
        || scores.is_null()
    {

        update_last_error(
            "Null pointer passed to \
             rssn_num_graph_page_rank"
                .to_string(),
        );

        return -1;
    }

    let g = &*graph;

    let n = g.num_nodes();

    let s = page_rank(
        g,
        damping_factor,
        tolerance,
        max_iter,
    );

    let scores_slice =
        slice::from_raw_parts_mut(
            scores,
            n,
        );

    scores_slice.copy_from_slice(&s);

    0
}

/// Computes Floyd-Warshall.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_graph_floyd_warshall(
    graph: *mut Graph,
    dist_matrix: *mut f64,
) -> i32 {

    if graph.is_null()
        || dist_matrix.is_null()
    {

        update_last_error("Null pointer passed to rssn_num_graph_floyd_warshall".to_string());

        return -1;
    }

    let g = &*graph;

    let n = g.num_nodes();

    let mat = floyd_warshall(g);

    let mat_slice =
        slice::from_raw_parts_mut(
            dist_matrix,
            n * n,
        );

    mat_slice.copy_from_slice(&mat);

    0
}

/// Computes Connected Components.
/// Result array `components` must be allocated by caller with size `num_nodes`.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_graph_connected_components(
    graph: *mut Graph,
    components: *mut usize,
) -> i32 {

    if graph.is_null()
        || components.is_null()
    {

        update_last_error("Null pointer passed to rssn_num_graph_connected_components".to_string());

        return -1;
    }

    let g = &*graph;

    let n = g.num_nodes();

    let comp = crate::numerical::graph::connected_components(g);

    let comp_slice =
        slice::from_raw_parts_mut(
            components,
            n,
        );

    comp_slice.copy_from_slice(&comp);

    0
}

/// Computes Minimum Spanning Tree (MST).
/// Returns a new Graph handle.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_graph_minimum_spanning_tree(
    graph: *mut Graph
) -> *mut Graph {

    if graph.is_null() {

        update_last_error(
            "Null pointer passed to rssn_num_graph_minimum_spanning_tree".to_string(),
        );

        return std::ptr::null_mut();
    }

    let g = &*graph;

    let mst = crate::numerical::graph::minimum_spanning_tree(g);

    Box::into_raw(Box::new(mst))
}
