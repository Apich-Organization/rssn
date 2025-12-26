//! JSON-based FFI API for numerical graph algorithms.

use crate::ffi_apis::common::{from_json_string, to_c_string};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::graph::{bfs, dijkstra, floyd_warshall, page_rank, Graph};
use serde::{Deserialize, Serialize};
use std::os::raw::c_char;

#[derive(Deserialize)]

struct Edge {
    u: usize,
    v: usize,
    w: f64,
}

#[derive(Deserialize)]

struct GraphDef {
    num_nodes: usize,
    edges: Vec<Edge>,
}

impl GraphDef {
    fn to_graph(&self) -> Graph {

        let mut g = Graph::new(self.num_nodes);

        for edge in &self.edges {

            g.add_edge(
                edge.u, edge.v, edge.w,
            );
        }

        g
    }
}

#[derive(Deserialize)]

struct DijkstraInput {
    graph: GraphDef,
    start_node: usize,
}

#[derive(Serialize)]

struct DijkstraOutput {
    dist: Vec<f64>,
    prev: Vec<Option<usize>>,
}

#[derive(Deserialize)]

struct PageRankInput {
    graph: GraphDef,
    damping_factor: f64,
    tolerance: f64,
    max_iter: usize,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_graph_dijkstra_json(input_json: *const c_char) -> *mut c_char {

    let input: DijkstraInput = match from_json_string(input_json) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<DijkstraOutput, String> {
                        ok: None,
                        err: Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
            )
        }
    };

    let g = input
        .graph
        .to_graph();

    let (dist, prev) = dijkstra(&g, input.start_node);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(DijkstraOutput { dist, prev }),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_graph_bfs_json(input_json: *const c_char) -> *mut c_char {

    let input: DijkstraInput = match from_json_string(input_json) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<usize>, String> {
                        ok: None,
                        err: Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
            )
        }
    };

    let g = input
        .graph
        .to_graph();

    let dist = bfs(&g, input.start_node);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(dist),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_graph_page_rank_json(input_json: *const c_char) -> *mut c_char {

    let input: PageRankInput = match from_json_string(input_json) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<f64>, String> {
                        ok: None,
                        err: Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
            )
        }
    };

    let g = input
        .graph
        .to_graph();

    let scores = page_rank(
        &g,
        input.damping_factor,
        input.tolerance,
        input.max_iter,
    );

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(scores),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_graph_floyd_warshall_json(
    input_json: *const c_char
) -> *mut c_char {

    let input: GraphDef = match from_json_string(input_json) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<f64>, String> {
                        ok: None,
                        err: Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
            )
        }
    };

    let g = input.to_graph();

    let mat = floyd_warshall(&g);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(mat),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[derive(Serialize)]

struct EdgeOut {
    u: usize,
    v: usize,
    w: f64,
}

#[derive(Serialize)]

struct GraphDefOut {
    num_nodes: usize,
    edges: Vec<EdgeOut>,
}

impl GraphDefOut {
    fn from_graph(graph: &Graph) -> Self {

        let num_nodes = graph.num_nodes();

        let mut edges = Vec::new();

        for u in 0..num_nodes {

            for &(v, w) in graph.adj(u) {

                edges.push(EdgeOut { u, v, w });
            }
        }

        Self { num_nodes, edges }
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_graph_connected_components_json(
    input_json: *const c_char
) -> *mut c_char {

    let input: GraphDef = match from_json_string(input_json) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<usize>, String> {
                        ok: None,
                        err: Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
            )
        }
    };

    let g = input.to_graph();

    let comp = crate::numerical::graph::connected_components(&g);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(comp),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_graph_minimum_spanning_tree_json(
    input_json: *const c_char
) -> *mut c_char {

    let input: GraphDef = match from_json_string(input_json) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<GraphDefOut, String> {
                        ok: None,
                        err: Some("Invalid JSON input".to_string()),
                    },
                )
                .unwrap(),
            )
        }
    };

    let g = input.to_graph();

    let mst = crate::numerical::graph::minimum_spanning_tree(&g);

    let mst_def = GraphDefOut::from_graph(&mst);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(mst_def),
            err: None::<String>,
        })
        .unwrap(),
    )
}
