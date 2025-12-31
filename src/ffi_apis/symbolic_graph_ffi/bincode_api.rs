use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::symbolic::core::Expr;
use crate::symbolic::graph::Graph;
use crate::symbolic::graph_algorithms::bfs;
use crate::symbolic::graph_algorithms::connected_components;
use crate::symbolic::graph_algorithms::dfs;
use crate::symbolic::graph_algorithms::edmonds_karp_max_flow;
use crate::symbolic::graph_algorithms::has_cycle;
use crate::symbolic::graph_algorithms::is_bipartite;
use crate::symbolic::graph_algorithms::kruskal_mst;

/// Creates a new graph from bincode specification.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_graph_new(
    spec_buf: BincodeBuffer
) -> BincodeBuffer {

    #[derive(serde::Deserialize)]

    struct GraphSpec {
        is_directed: bool,
    }

    let spec : GraphSpec = match from_bincode_buffer(&spec_buf) {
        | Some(s) => s,
        | None => return BincodeBuffer::empty(),
    };

    let graph: Graph<String> =
        Graph::new(spec.is_directed);

    to_bincode_buffer(&graph)
}

/// Adds a node to the graph.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_graph_add_node(
    input_buf: BincodeBuffer
) -> BincodeBuffer {

    #[derive(serde::Deserialize)]

    struct Input {
        graph: Graph<String>,
        label: String,
    }

    let mut input : Input = match from_bincode_buffer(&input_buf) {
        | Some(i) => i,
        | None => return BincodeBuffer::empty(),
    };

    input
        .graph
        .add_node(input.label);

    to_bincode_buffer(&input.graph)
}

/// Adds an edge to the graph.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_graph_add_edge(
    input_buf: BincodeBuffer
) -> BincodeBuffer {

    #[derive(serde::Deserialize)]

    struct Input {
        graph: Graph<String>,
        from: String,
        to: String,
        weight: Expr,
    }

    let mut input : Input = match from_bincode_buffer(&input_buf) {
        | Some(i) => i,
        | None => return BincodeBuffer::empty(),
    };

    input
        .graph
        .add_edge(
            &input.from,
            &input.to,
            input.weight,
        );

    to_bincode_buffer(&input.graph)
}

/// Gets the adjacency matrix.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_graph_adjacency_matrix(
    graph_buf: BincodeBuffer
) -> BincodeBuffer {

    let graph : Graph<String> = match from_bincode_buffer(&graph_buf) {
        | Some(g) => g,
        | None => return BincodeBuffer::empty(),
    };

    let matrix =
        graph.to_adjacency_matrix();

    to_bincode_buffer(&matrix)
}

/// Gets the Laplacian matrix.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_graph_laplacian_matrix(
    graph_buf: BincodeBuffer
) -> BincodeBuffer {

    let graph : Graph<String> = match from_bincode_buffer(&graph_buf) {
        | Some(g) => g,
        | None => return BincodeBuffer::empty(),
    };

    let matrix =
        graph.to_laplacian_matrix();

    to_bincode_buffer(&matrix)
}

/// Performs BFS traversal.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_graph_bfs(
    input_buf: BincodeBuffer
) -> BincodeBuffer {

    #[derive(serde::Deserialize)]

    struct Input {
        graph: Graph<String>,
        start_node: usize,
    }

    let input : Input = match from_bincode_buffer(&input_buf) {
        | Some(i) => i,
        | None => return BincodeBuffer::empty(),
    };

    let result = bfs(
        &input.graph,
        input.start_node,
    );

    to_bincode_buffer(&result)
}

/// Performs DFS traversal.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_graph_dfs(
    input_buf: BincodeBuffer
) -> BincodeBuffer {

    #[derive(serde::Deserialize)]

    struct Input {
        graph: Graph<String>,
        start_node: usize,
    }

    let input : Input = match from_bincode_buffer(&input_buf) {
        | Some(i) => i,
        | None => return BincodeBuffer::empty(),
    };

    let result = dfs(
        &input.graph,
        input.start_node,
    );

    to_bincode_buffer(&result)
}

/// Finds connected components.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_graph_connected_components(
    graph_buf: BincodeBuffer
) -> BincodeBuffer {

    let graph : Graph<String> = match from_bincode_buffer(&graph_buf) {
        | Some(g) => g,
        | None => return BincodeBuffer::empty(),
    };

    let result =
        connected_components(&graph);

    to_bincode_buffer(&result)
}

/// Computes maximum flow.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_graph_max_flow(
    input_buf: BincodeBuffer
) -> BincodeBuffer {

    #[derive(serde::Deserialize)]

    struct Input {
        graph: Graph<String>,
        source: usize,
        sink: usize,
    }

    let input : Input = match from_bincode_buffer(&input_buf) {
        | Some(i) => i,
        | None => return BincodeBuffer::empty(),
    };

    let flow = edmonds_karp_max_flow(
        &input.graph,
        input.source,
        input.sink,
    );

    to_bincode_buffer(&flow)
}

/// Computes MST using Kruskal's algorithm.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_graph_kruskal_mst(
    graph_buf: BincodeBuffer
) -> BincodeBuffer {

    let graph : Graph<String> = match from_bincode_buffer(&graph_buf) {
        | Some(g) => g,
        | None => return BincodeBuffer::empty(),
    };

    let mst = kruskal_mst(&graph);

    to_bincode_buffer(&mst)
}

/// Checks if graph has a cycle.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_graph_has_cycle(
    graph_buf: BincodeBuffer
) -> bool {

    let graph: Graph<String> =
        match from_bincode_buffer(
            &graph_buf,
        ) {
            | Some(g) => g,
            | None => return false,
        };

    has_cycle(&graph)
}

/// Checks if graph is bipartite.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_graph_is_bipartite(
    graph_buf: BincodeBuffer
) -> bool {

    let graph: Graph<String> =
        match from_bincode_buffer(
            &graph_buf,
        ) {
            | Some(g) => g,
            | None => return false,
        };

    is_bipartite(&graph).is_some()
}
