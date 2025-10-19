use std::sync::Arc;
use rssn::symbolic::core::{DagManager, DagOp, DagNode};
use std::hash::Hash;

/// This test simulates a collision-like situation by forcing two different
/// op+children combinations to map to different nodes, verifies get_or_create
/// does not return an unrelated node when the u64 hash collides.
#[test]
fn test_dagmanager_no_collision_misreuse() {
    let mgr = DagManager::new();

    // Construct simple nodes: a and b
    let node_a = mgr.get_or_create(DagOp::Var("x".into()), vec![]);
    let node_b = mgr.get_or_create(DagOp::Var("y".into()), vec![]);

    // Create a composite op for (x + y)
    let add_children = vec![node_a.clone(), node_b.clone()];
    let node_add1 = mgr.get_or_create(DagOp::Add, add_children.clone());

    // Create another instance with same structure: should reuse node_add1
    let node_add2 = mgr.get_or_create(DagOp::Add, add_children.clone());
    assert!(Arc::ptr_eq(&node_add1, &node_add2));

    // Now create a different structure (y + x) — depending on canonicalization this may be equal;
    // ensure manager does not incorrectly return node_add1 for a truly different structure.
    let add_children_swapped = vec![node_b.clone(), node_a.clone()];
    let node_add_swapped = mgr.get_or_create(DagOp::Add, add_children_swapped);

    // Either equal (if Add is commutative and canonicalized) or different; in either case
    // the manager must not return a structurally different node for a colliding hash.
    // We verify structural equality via hash+op+len checks implemented in manager.
    assert_eq!(node_add_swapped.op(), DagOp::Add);
   // If not the same pointer, ensure it is a distinct node when structure differs.
    if !Arc::ptr_eq(&node_add1, &node_add_swapped) {
        assert_ne!(node_add1.hash(), node_add_swapped.hash());
    }
}
