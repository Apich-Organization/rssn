use crate::ffi_apis::symbolic_graph_ffi::handle::RssnGraph;
use crate::symbolic::core::Expr;
use crate::symbolic::graph::Graph;
use crate::symbolic::graph_operations::*;
use std::ffi::CStr;
use std::os::raw::c_char;

pub(crate) fn convert_expr_graph_to_string_graph(g: Graph<Expr>) -> Graph<String> {

    let mut new_graph = Graph::new(g.is_directed);

    let mut id_map = std::collections::HashMap::new();

    for (i, node_label) in g
        .nodes
        .iter()
        .enumerate()
    {

        let label_str = format!("{:?}", node_label);

        let new_id = new_graph.add_node(label_str);

        id_map.insert(i, new_id);
    }

    for (u, v, weight) in g.get_edges() {

        let u_new = id_map[&u];

        let v_new = id_map[&v];

        new_graph.add_edge(
            &new_graph.nodes[u_new].clone(),
            &new_graph.nodes[v_new].clone(),
            weight,
        );
    }

    new_graph
}

/// Creates an induced subgraph.
#[no_mangle]

pub extern "C" fn rssn_graph_induced_subgraph(
    ptr: *const RssnGraph,
    node_labels: *const *const c_char,
    count: usize,
) -> *mut RssnGraph {

    if ptr.is_null() || node_labels.is_null() {

        return std::ptr::null_mut();
    }

    let graph = unsafe {

        &*(ptr as *const Graph<String>)
    };

    let mut labels = Vec::with_capacity(count);

    for i in 0..count {

        let c_str = unsafe {

            *node_labels.add(i)
        };

        if c_str.is_null() {

            continue;
        }

        let label = unsafe {

            CStr::from_ptr(c_str)
                .to_string_lossy()
                .into_owned()
        };

        labels.push(label);
    }

    let result = induced_subgraph(graph, &labels);

    Box::into_raw(Box::new(result)) as *mut RssnGraph
}

/// Computes the union of two graphs.
#[no_mangle]

pub extern "C" fn rssn_graph_union(
    ptr1: *const RssnGraph,
    ptr2: *const RssnGraph,
) -> *mut RssnGraph {

    if ptr1.is_null() || ptr2.is_null() {

        return std::ptr::null_mut();
    }

    let g1 = unsafe {

        &*(ptr1 as *const Graph<String>)
    };

    let g2 = unsafe {

        &*(ptr2 as *const Graph<String>)
    };

    let result = union(g1, g2);

    Box::into_raw(Box::new(result)) as *mut RssnGraph
}

/// Computes the intersection of two graphs.
#[no_mangle]

pub extern "C" fn rssn_graph_intersection(
    ptr1: *const RssnGraph,
    ptr2: *const RssnGraph,
) -> *mut RssnGraph {

    if ptr1.is_null() || ptr2.is_null() {

        return std::ptr::null_mut();
    }

    let g1 = unsafe {

        &*(ptr1 as *const Graph<String>)
    };

    let g2 = unsafe {

        &*(ptr2 as *const Graph<String>)
    };

    let result = intersection(g1, g2);

    Box::into_raw(Box::new(result)) as *mut RssnGraph
}

/// Computes the Cartesian product of two graphs.
#[no_mangle]

pub extern "C" fn rssn_graph_cartesian_product(
    ptr1: *const RssnGraph,
    ptr2: *const RssnGraph,
) -> *mut RssnGraph {

    if ptr1.is_null() || ptr2.is_null() {

        return std::ptr::null_mut();
    }

    let g1 = unsafe {

        &*(ptr1 as *const Graph<String>)
    };

    let g2 = unsafe {

        &*(ptr2 as *const Graph<String>)
    };

    let result_expr = cartesian_product(g1, g2);

    let result = convert_expr_graph_to_string_graph(result_expr);

    Box::into_raw(Box::new(result)) as *mut RssnGraph
}

/// Computes the Tensor product of two graphs.
#[no_mangle]

pub extern "C" fn rssn_graph_tensor_product(
    ptr1: *const RssnGraph,
    ptr2: *const RssnGraph,
) -> *mut RssnGraph {

    if ptr1.is_null() || ptr2.is_null() {

        return std::ptr::null_mut();
    }

    let g1 = unsafe {

        &*(ptr1 as *const Graph<String>)
    };

    let g2 = unsafe {

        &*(ptr2 as *const Graph<String>)
    };

    let result_expr = tensor_product(g1, g2);

    let result = convert_expr_graph_to_string_graph(result_expr);

    Box::into_raw(Box::new(result)) as *mut RssnGraph
}

/// Computes the complement of a graph.
#[no_mangle]

pub extern "C" fn rssn_graph_complement(ptr: *const RssnGraph) -> *mut RssnGraph {

    if ptr.is_null() {

        return std::ptr::null_mut();
    }

    let graph = unsafe {

        &*(ptr as *const Graph<String>)
    };

    let result = complement(graph);

    Box::into_raw(Box::new(result)) as *mut RssnGraph
}

/// Computes the disjoint union of two graphs.
#[no_mangle]

pub extern "C" fn rssn_graph_disjoint_union(
    ptr1: *const RssnGraph,
    ptr2: *const RssnGraph,
) -> *mut RssnGraph {

    if ptr1.is_null() || ptr2.is_null() {

        return std::ptr::null_mut();
    }

    let g1 = unsafe {

        &*(ptr1 as *const Graph<String>)
    };

    let g2 = unsafe {

        &*(ptr2 as *const Graph<String>)
    };

    let result_expr = disjoint_union(g1, g2);

    let result = convert_expr_graph_to_string_graph(result_expr);

    Box::into_raw(Box::new(result)) as *mut RssnGraph
}

/// Computes the join of two graphs.
#[no_mangle]

pub extern "C" fn rssn_graph_join(
    ptr1: *const RssnGraph,
    ptr2: *const RssnGraph,
) -> *mut RssnGraph {

    if ptr1.is_null() || ptr2.is_null() {

        return std::ptr::null_mut();
    }

    let g1 = unsafe {

        &*(ptr1 as *const Graph<String>)
    };

    let g2 = unsafe {

        &*(ptr2 as *const Graph<String>)
    };

    let result_expr = join(g1, g2);

    let result = convert_expr_graph_to_string_graph(result_expr);

    Box::into_raw(Box::new(result)) as *mut RssnGraph
}
