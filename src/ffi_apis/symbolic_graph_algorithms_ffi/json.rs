use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_json_string;
use crate::symbolic::graph::Graph;
use crate::symbolic::graph_algorithms::bfs;
use crate::symbolic::graph_algorithms::bipartite_maximum_matching;
use crate::symbolic::graph_algorithms::connected_components;
use crate::symbolic::graph_algorithms::dfs;
use crate::symbolic::graph_algorithms::dinic_max_flow;
use crate::symbolic::graph_algorithms::edmonds_karp_max_flow;
use crate::symbolic::graph_algorithms::find_bridges_and_articulation_points;
use crate::symbolic::graph_algorithms::has_cycle;
use crate::symbolic::graph_algorithms::is_bipartite;
use crate::symbolic::graph_algorithms::is_connected;
use crate::symbolic::graph_algorithms::kruskal_mst;
use crate::symbolic::graph_algorithms::strongly_connected_components;
use crate::symbolic::graph_algorithms::topological_sort;

/// Performs DFS traversal.
/// Input: {"graph": Graph, "`start_node"`: usize}
/// Output: [usize] (array of node indices)
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_graph_dfs_api(
    json: *const std::os::raw::c_char
) -> *mut std::os::raw::c_char {

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

/// Performs BFS traversal.
/// Input: {"graph": Graph, "`start_node"`: usize}
/// Output: [usize]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_graph_bfs_api(
    json: *const std::os::raw::c_char
) -> *mut std::os::raw::c_char {

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

/// Finds connected components.
/// Input: Graph
/// Output: [[usize]] (array of arrays)
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_graph_connected_components_api(
    json: *const std::os::raw::c_char
) -> *mut std::os::raw::c_char {

    let graph : Graph<String> = match from_json_string(json) {
        | Some(g) => g,
        | None => return std::ptr::null_mut(),
    };

    let result =
        connected_components(&graph);

    to_json_string(&result)
}

/// Checks if graph is connected.
/// Input: Graph
/// Output: bool
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_graph_is_connected(
    json: *const std::os::raw::c_char
) -> *mut std::os::raw::c_char {

    let graph : Graph<String> = match from_json_string(json) {
        | Some(g) => g,
        | None => return std::ptr::null_mut(),
    };

    let result = is_connected(&graph);

    to_json_string(&result)
}

/// Finds strongly connected components.
/// Input: Graph
/// Output: [[usize]]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_graph_strongly_connected_components(
    json: *const std::os::raw::c_char
) -> *mut std::os::raw::c_char {

    let graph : Graph<String> = match from_json_string(json) {
        | Some(g) => g,
        | None => return std::ptr::null_mut(),
    };

    let result =
        strongly_connected_components(
            &graph,
        );

    to_json_string(&result)
}

/// Checks if graph has a cycle.
/// Input: Graph
/// Output: bool
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_graph_has_cycle_api(
    json: *const std::os::raw::c_char
) -> *mut std::os::raw::c_char {

    let graph : Graph<String> = match from_json_string(json) {
        | Some(g) => g,
        | None => return std::ptr::null_mut(),
    };

    let result = has_cycle(&graph);

    to_json_string(&result)
}

/// Finds bridges and articulation points.
/// Input: Graph
/// Output: {"bridges": [(usize, usize)], "`articulation_points"`: [usize]}
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_graph_bridges_and_articulation_points(
    json: *const std::os::raw::c_char
) -> *mut std::os::raw::c_char {

    #[derive(serde::Serialize)]

    struct BridgesResult {
        bridges: Vec<(usize, usize)>,
        articulation_points: Vec<usize>,
    }

    let graph : Graph<String> = match from_json_string(json) {
        | Some(g) => g,
        | None => return std::ptr::null_mut(),
    };

    let (bridges, aps) = find_bridges_and_articulation_points(&graph);

    let result = BridgesResult {
        bridges,
        articulation_points: aps,
    };

    to_json_string(&result)
}

/// Computes MST using Kruskal's algorithm.
/// Input: Graph
/// Output: Graph (MST)
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_graph_kruskal_mst_api(
    json: *const std::os::raw::c_char
) -> *mut std::os::raw::c_char {

    let graph : Graph<String> = match from_json_string(json) {
        | Some(g) => g,
        | None => return std::ptr::null_mut(),
    };

    let mst_edges = kruskal_mst(&graph);

    let mut mst_graph =
        Graph::new(graph.is_directed);

    for node in &graph.nodes {

        mst_graph
            .add_node(node.clone());
    }

    for (u, v, weight) in mst_edges {

        mst_graph.add_edge(
            &graph.nodes[u],
            &graph.nodes[v],
            weight,
        );
    }

    to_json_string(&mst_graph)
}

/// Computes maximum flow using Edmonds-Karp.
/// Input: {"graph": Graph, "source": usize, "sink": usize}
/// Output: f64
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_graph_edmonds_karp_max_flow(
    json: *const std::os::raw::c_char
) -> *mut std::os::raw::c_char {

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

    let result = edmonds_karp_max_flow(
        &input.graph,
        input.source,
        input.sink,
    );

    to_json_string(&result)
}

/// Computes maximum flow using Dinic's algorithm.
/// Input: {"graph": Graph, "source": usize, "sink": usize}
/// Output: f64
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_graph_dinic_max_flow(
    json: *const std::os::raw::c_char
) -> *mut std::os::raw::c_char {

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

    let result = dinic_max_flow(
        &input.graph,
        input.source,
        input.sink,
    );

    to_json_string(&result)
}

/// Checks if graph is bipartite.
/// Input: Graph
/// Output: [i8] or null
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_graph_is_bipartite_api(
    json: *const std::os::raw::c_char
) -> *mut std::os::raw::c_char {

    let graph : Graph<String> = match from_json_string(json) {
        | Some(g) => g,
        | None => return std::ptr::null_mut(),
    };

    let result = is_bipartite(&graph);

    to_json_string(&result)
}

/// Finds maximum matching in bipartite graph.
/// Input: {"graph": Graph, "partition": [i8]}
/// Output: [(usize, usize)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_graph_bipartite_maximum_matching(
    json: *const std::os::raw::c_char
) -> *mut std::os::raw::c_char {

    #[derive(serde::Deserialize)]

    struct Input {
        graph: Graph<String>,
        partition: Vec<i8>,
    }

    let input : Input = match from_json_string(json) {
        | Some(i) => i,
        | None => return std::ptr::null_mut(),
    };

    let result =
        bipartite_maximum_matching(
            &input.graph,
            &input.partition,
        );

    to_json_string(&result)
}

/// Performs topological sort.
/// Input: Graph
/// Output: [usize] or null
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_graph_topological_sort(
    json: *const std::os::raw::c_char
) -> *mut std::os::raw::c_char {

    let graph : Graph<String> = match from_json_string(json) {
        | Some(g) => g,
        | None => return std::ptr::null_mut(),
    };

    let result =
        topological_sort(&graph);

    to_json_string(&result)
}
