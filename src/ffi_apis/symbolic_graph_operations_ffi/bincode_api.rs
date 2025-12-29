use crate::ffi_apis::common::{BincodeBuffer, from_bincode_buffer, to_bincode_buffer};
use crate::ffi_apis::symbolic_graph_operations_ffi::handle::convert_expr_graph_to_string_graph;
use crate::symbolic::graph::Graph;
use crate::symbolic::graph_operations::{induced_subgraph, union, intersection, cartesian_product, tensor_product, complement, disjoint_union, join};

/// Creates an induced subgraph.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_graph_induced_subgraph(
    input_buf: BincodeBuffer
) -> BincodeBuffer {

    #[derive(serde::Deserialize)]

    struct Input {
        graph: Graph<String>,
        nodes: Vec<String>,
    }

    let input : Input = match from_bincode_buffer(&input_buf) {
        | Some(i) => i,
        | None => return BincodeBuffer::empty(),
    };

    let result = induced_subgraph(
        &input.graph,
        &input.nodes,
    );

    to_bincode_buffer(&result)
}

/// Computes the union of two graphs.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_graph_union(
    input_buf: BincodeBuffer
) -> BincodeBuffer {

    #[derive(serde::Deserialize)]

    struct Input {
        g1: Graph<String>,
        g2: Graph<String>,
    }

    let input : Input = match from_bincode_buffer(&input_buf) {
        | Some(i) => i,
        | None => return BincodeBuffer::empty(),
    };

    let result =
        union(&input.g1, &input.g2);

    to_bincode_buffer(&result)
}

/// Computes the intersection of two graphs.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_graph_intersection(
    input_buf: BincodeBuffer
) -> BincodeBuffer {

    #[derive(serde::Deserialize)]

    struct Input {
        g1: Graph<String>,
        g2: Graph<String>,
    }

    let input : Input = match from_bincode_buffer(&input_buf) {
        | Some(i) => i,
        | None => return BincodeBuffer::empty(),
    };

    let result = intersection(
        &input.g1,
        &input.g2,
    );

    to_bincode_buffer(&result)
}

/// Computes the Cartesian product of two graphs.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_graph_cartesian_product(
    input_buf: BincodeBuffer
) -> BincodeBuffer {

    #[derive(serde::Deserialize)]

    struct Input {
        g1: Graph<String>,
        g2: Graph<String>,
    }

    let input : Input = match from_bincode_buffer(&input_buf) {
        | Some(i) => i,
        | None => return BincodeBuffer::empty(),
    };

    let result_expr = cartesian_product(
        &input.g1,
        &input.g2,
    );

    let result = convert_expr_graph_to_string_graph(result_expr);

    to_bincode_buffer(&result)
}

/// Computes the Tensor product of two graphs.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_graph_tensor_product(
    input_buf: BincodeBuffer
) -> BincodeBuffer {

    #[derive(serde::Deserialize)]

    struct Input {
        g1: Graph<String>,
        g2: Graph<String>,
    }

    let input : Input = match from_bincode_buffer(&input_buf) {
        | Some(i) => i,
        | None => return BincodeBuffer::empty(),
    };

    let result_expr = tensor_product(
        &input.g1,
        &input.g2,
    );

    let result = convert_expr_graph_to_string_graph(result_expr);

    to_bincode_buffer(&result)
}

/// Computes the complement of a graph.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_graph_complement(
    graph_buf: BincodeBuffer
) -> BincodeBuffer {

    let graph : Graph<String> = match from_bincode_buffer(&graph_buf) {
        | Some(g) => g,
        | None => return BincodeBuffer::empty(),
    };

    let result = complement(&graph);

    to_bincode_buffer(&result)
}

/// Computes the disjoint union of two graphs.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_graph_disjoint_union(
    input_buf: BincodeBuffer
) -> BincodeBuffer {

    #[derive(serde::Deserialize)]

    struct Input {
        g1: Graph<String>,
        g2: Graph<String>,
    }

    let input : Input = match from_bincode_buffer(&input_buf) {
        | Some(i) => i,
        | None => return BincodeBuffer::empty(),
    };

    let result_expr = disjoint_union(
        &input.g1,
        &input.g2,
    );

    let result = convert_expr_graph_to_string_graph(result_expr);

    to_bincode_buffer(&result)
}

/// Computes the join of two graphs.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_graph_join(
    input_buf: BincodeBuffer
) -> BincodeBuffer {

    #[derive(serde::Deserialize)]

    struct Input {
        g1: Graph<String>,
        g2: Graph<String>,
    }

    let input : Input = match from_bincode_buffer(&input_buf) {
        | Some(i) => i,
        | None => return BincodeBuffer::empty(),
    };

    let result_expr =
        join(&input.g1, &input.g2);

    let result = convert_expr_graph_to_string_graph(result_expr);

    to_bincode_buffer(&result)
}
