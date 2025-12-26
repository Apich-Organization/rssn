//! Handle-based FFI API for numerical topology.

use std::ptr;

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::graph::Graph;
use crate::numerical::topology;

/// Finds the connected components of a graph.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_topology_find_connected_components(
    graph_ptr : *const Graph
) -> *mut Vec<Vec<usize>> {

    if graph_ptr.is_null() {

        return ptr::null_mut();
    }

    let graph = &*graph_ptr;

    let components = topology::find_connected_components(graph);

    Box::into_raw(Box::new(components))
}

/// Computes the Betti numbers for a point cloud.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_topology_betti_numbers(
    points : *const *const f64,
    n_points : usize,
    dim : usize,
    epsilon : f64,
    max_dim : usize,
    result : *mut usize,
) -> i32 {

    if points.is_null()
        || result.is_null()
    {

        return -1;
    }

    let mut pt_slices =
        Vec::with_capacity(n_points);

    for i in 0 .. n_points {

        pt_slices.push(
            std::slice::from_raw_parts(
                *points.add(i),
                dim,
            ),
        );
    }

    let betti = topology::betti_numbers_at_radius(
        &pt_slices,
        epsilon,
        max_dim,
    );

    let n_betti = betti.len();

    std::ptr::copy_nonoverlapping(
        betti.as_ptr(),
        result,
        n_betti,
    );

    0
}

/// Computes the Euclidean distance between two points.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_topology_euclidean_distance(
    p1 : *const f64,
    p2 : *const f64,
    dim : usize,
) -> f64 {

    if p1.is_null() || p2.is_null() {

        return f64::NAN;
    }

    let s1 = std::slice::from_raw_parts(
        p1, dim,
    );

    let s2 = std::slice::from_raw_parts(
        p2, dim,
    );

    topology::euclidean_distance(s1, s2)
}
