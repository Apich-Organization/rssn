use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::symbolic::graph::Graph;
use crate::symbolic::graph_isomorphism_and_coloring::are_isomorphic_heuristic;
use crate::symbolic::graph_isomorphism_and_coloring::chromatic_number_exact;
use crate::symbolic::graph_isomorphism_and_coloring::greedy_coloring;

/// Checks if two graphs are isomorphic.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_bincode_are_isomorphic_heuristic(
    input_buf: BincodeBuffer
) -> BincodeBuffer {
    #[derive(serde::Deserialize)]
    struct Input {
        g1: Graph<String>,
        g2: Graph<String>,
    }

    let input: Input = match from_bincode_buffer(&input_buf) {
        | Some(i) => i,
        | None => return BincodeBuffer::empty(),
    };

    let result = are_isomorphic_heuristic(&input.g1, &input.g2);

    to_bincode_buffer(&result)
}

/// Greedy coloring.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_bincode_greedy_coloring(graph_buf: BincodeBuffer) -> BincodeBuffer {
    let graph: Graph<String> = match from_bincode_buffer(&graph_buf) {
        | Some(g) => g,
        | None => return BincodeBuffer::empty(),
    };

    let result = greedy_coloring(&graph);

    to_bincode_buffer(&result)
}

/// Exact chromatic number.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_bincode_chromatic_number_exact(
    graph_buf: BincodeBuffer
) -> BincodeBuffer {
    let graph: Graph<String> = match from_bincode_buffer(&graph_buf) {
        | Some(g) => g,
        | None => return BincodeBuffer::empty(),
    };

    let result = chromatic_number_exact(&graph);

    to_bincode_buffer(&result)
}
