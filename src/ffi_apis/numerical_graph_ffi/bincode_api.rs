//! Bincode-based FFI API for numerical graph algorithms.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::graph::bfs;
use crate::numerical::graph::dijkstra;
use crate::numerical::graph::floyd_warshall;
use crate::numerical::graph::page_rank;
use crate::numerical::graph::Graph;

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

        let mut g =
            Graph::new(self.num_nodes);

        for edge in &self.edges {

            g.add_edge(
                edge.u,
                edge.v,
                edge.w,
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

/// Computes Dijkstra's shortest path algorithm on a graph using bincode for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_graph_dijkstra_bincode(

    buffer: BincodeBuffer

) -> BincodeBuffer {

    let input : DijkstraInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<DijkstraOutput, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
    };

    let g = input
        .graph
        .to_graph();

    let (dist, prev) =
        dijkstra(&g, input.start_node);

    to_bincode_buffer(&FfiResult {
        ok: Some(DijkstraOutput {
            dist,
            prev,
        }),
        err: None::<String>,
    })
}

/// Computes Breadth-First Search (BFS) on a graph using bincode for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_graph_bfs_bincode(

    buffer: BincodeBuffer

) -> BincodeBuffer {

    let input : DijkstraInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<usize>, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
    };

    let g = input
        .graph
        .to_graph();

    let dist =
        bfs(&g, input.start_node);

    to_bincode_buffer(&FfiResult {
        ok: Some(dist),
        err: None::<String>,
    })
}

/// Computes the PageRank scores for nodes in a graph using bincode for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_graph_page_rank_bincode(

    buffer: BincodeBuffer

) -> BincodeBuffer {

    let input : PageRankInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<f64>, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
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

    to_bincode_buffer(&FfiResult {
        ok: Some(scores),
        err: None::<String>,
    })
}

/// Computes the Floyd-Warshall all-pairs shortest path algorithm on a graph using bincode for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_graph_floyd_warshall_bincode(

    buffer: BincodeBuffer

) -> BincodeBuffer {

    let input : GraphDef = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<f64>, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
    };

    let g = input.to_graph();

    let mat = floyd_warshall(&g);

    to_bincode_buffer(&FfiResult {
        ok: Some(mat),
        err: None::<String>,
    })
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
    fn from_graph(
        graph: &Graph
    ) -> Self {

        let num_nodes =
            graph.num_nodes();

        let mut edges = Vec::new();

        for u in 0 .. num_nodes {

            for &(v, w) in graph.adj(u)
            {

                edges.push(EdgeOut {
                    u,
                    v,
                    w,
                });
            }
        }

        Self {
            num_nodes,
            edges,
        }
    }
}

/// Computes the connected components of a graph using bincode for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_graph_connected_components_bincode(

    buffer: BincodeBuffer

) -> BincodeBuffer {

    let input : GraphDef = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<usize>, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
    };

    let g = input.to_graph();

    let comp = crate::numerical::graph::connected_components(&g);

    to_bincode_buffer(&FfiResult {
        ok: Some(comp),
        err: None::<String>,
    })
}

/// Computes the Minimum Spanning Tree (MST) of a graph using bincode for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_graph_minimum_spanning_tree_bincode(

    buffer: BincodeBuffer

) -> BincodeBuffer {

    let input : GraphDef = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<GraphDefOut, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
    };

    let g = input.to_graph();

    let mst = crate::numerical::graph::minimum_spanning_tree(&g);

    let mst_def =
        GraphDefOut::from_graph(&mst);

    to_bincode_buffer(&FfiResult {
        ok: Some(mst_def),
        err: None::<String>,
    })
}
