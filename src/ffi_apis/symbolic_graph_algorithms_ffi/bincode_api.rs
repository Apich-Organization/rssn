use crate::ffi_apis::common::*;
use crate::symbolic::graph::Graph;
use crate::symbolic::graph_algorithms::*;

/// Performs DFS traversal.
#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_graph_dfs_api(input_buf: BincodeBuffer) -> BincodeBuffer {

    #[derive(serde::Deserialize)]

    struct Input {
        graph: Graph<String>,
        start_node: usize,
    }

    let input: Input = match from_bincode_buffer(&input_buf) {
        Some(i) => i,
        None => return BincodeBuffer::empty(),
    };

    let result = dfs(&input.graph, input.start_node);

    to_bincode_buffer(&result)
}

/// Performs BFS traversal.
#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_graph_bfs_api(input_buf: BincodeBuffer) -> BincodeBuffer {

    #[derive(serde::Deserialize)]

    struct Input {
        graph: Graph<String>,
        start_node: usize,
    }

    let input: Input = match from_bincode_buffer(&input_buf) {
        Some(i) => i,
        None => return BincodeBuffer::empty(),
    };

    let result = bfs(&input.graph, input.start_node);

    to_bincode_buffer(&result)
}

/// Finds connected components.
#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_graph_connected_components_api(
    graph_buf: BincodeBuffer
) -> BincodeBuffer {

    let graph: Graph<String> = match from_bincode_buffer(&graph_buf) {
        Some(g) => g,
        None => return BincodeBuffer::empty(),
    };

    let result = connected_components(&graph);

    to_bincode_buffer(&result)
}

/// Checks if graph is connected.
#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_graph_is_connected(
    graph_buf: BincodeBuffer
) -> BincodeBuffer {

    let graph: Graph<String> = match from_bincode_buffer(&graph_buf) {
        Some(g) => g,
        None => return BincodeBuffer::empty(),
    };

    let result = is_connected(&graph);

    to_bincode_buffer(&result)
}

/// Finds strongly connected components.
#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_graph_strongly_connected_components(
    graph_buf: BincodeBuffer
) -> BincodeBuffer {

    let graph: Graph<String> = match from_bincode_buffer(&graph_buf) {
        Some(g) => g,
        None => return BincodeBuffer::empty(),
    };

    let result = strongly_connected_components(&graph);

    to_bincode_buffer(&result)
}

/// Checks if graph has a cycle.
#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_graph_has_cycle_api(
    graph_buf: BincodeBuffer
) -> BincodeBuffer {

    let graph: Graph<String> = match from_bincode_buffer(&graph_buf) {
        Some(g) => g,
        None => return BincodeBuffer::empty(),
    };

    let result = has_cycle(&graph);

    to_bincode_buffer(&result)
}

/// Finds bridges and articulation points.
#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_graph_bridges_and_articulation_points(
    graph_buf: BincodeBuffer
) -> BincodeBuffer {

    let graph: Graph<String> = match from_bincode_buffer(&graph_buf) {
        Some(g) => g,
        None => return BincodeBuffer::empty(),
    };

    let (bridges, aps) = find_bridges_and_articulation_points(&graph);

    #[derive(serde::Serialize)]

    struct Result {
        bridges: Vec<(usize, usize)>,
        articulation_points: Vec<usize>,
    }

    let result = Result {
        bridges,
        articulation_points: aps,
    };

    to_bincode_buffer(&result)
}

/// Computes MST using Kruskal's algorithm.
#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_graph_kruskal_mst_api(
    graph_buf: BincodeBuffer
) -> BincodeBuffer {

    let graph: Graph<String> = match from_bincode_buffer(&graph_buf) {
        Some(g) => g,
        None => return BincodeBuffer::empty(),
    };

    let mst_edges = kruskal_mst(&graph);

    let mut mst_graph = Graph::new(graph.is_directed);

    for node in &graph.nodes {

        mst_graph.add_node(node.clone());
    }

    for (u, v, weight) in mst_edges {

        mst_graph.add_edge(&graph.nodes[u], &graph.nodes[v], weight);
    }

    to_bincode_buffer(&mst_graph)
}

/// Computes maximum flow using Edmonds-Karp.
#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_graph_edmonds_karp_max_flow(
    input_buf: BincodeBuffer
) -> BincodeBuffer {

    #[derive(serde::Deserialize)]

    struct Input {
        graph: Graph<String>,
        source: usize,
        sink: usize,
    }

    let input: Input = match from_bincode_buffer(&input_buf) {
        Some(i) => i,
        None => return BincodeBuffer::empty(),
    };

    let result = edmonds_karp_max_flow(&input.graph, input.source, input.sink);

    to_bincode_buffer(&result)
}

/// Computes maximum flow using Dinic's algorithm.
#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_graph_dinic_max_flow(
    input_buf: BincodeBuffer
) -> BincodeBuffer {

    #[derive(serde::Deserialize)]

    struct Input {
        graph: Graph<String>,
        source: usize,
        sink: usize,
    }

    let input: Input = match from_bincode_buffer(&input_buf) {
        Some(i) => i,
        None => return BincodeBuffer::empty(),
    };

    let result = dinic_max_flow(&input.graph, input.source, input.sink);

    to_bincode_buffer(&result)
}

/// Checks if graph is bipartite.
#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_graph_is_bipartite_api(
    graph_buf: BincodeBuffer
) -> BincodeBuffer {

    let graph: Graph<String> = match from_bincode_buffer(&graph_buf) {
        Some(g) => g,
        None => return BincodeBuffer::empty(),
    };

    let result = is_bipartite(&graph);

    to_bincode_buffer(&result)
}

/// Finds maximum matching in bipartite graph.
#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_graph_bipartite_maximum_matching(
    input_buf: BincodeBuffer
) -> BincodeBuffer {

    #[derive(serde::Deserialize)]

    struct Input {
        graph: Graph<String>,
        partition: Vec<i8>,
    }

    let input: Input = match from_bincode_buffer(&input_buf) {
        Some(i) => i,
        None => return BincodeBuffer::empty(),
    };

    let result = bipartite_maximum_matching(&input.graph, &input.partition);

    to_bincode_buffer(&result)
}

/// Performs topological sort.
#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_graph_topological_sort(
    graph_buf: BincodeBuffer
) -> BincodeBuffer {

    let graph: Graph<String> = match from_bincode_buffer(&graph_buf) {
        Some(g) => g,
        None => return BincodeBuffer::empty(),
    };

    let result = topological_sort(&graph);

    to_bincode_buffer(&result)
}
