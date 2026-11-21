//! Reverse Cuthill-McKee for sorted, deduplicated symmetric adjacency.

/// Applies RCM to a symmetric graph whose rows list distinct neighbours in
/// ascending order. The final reversal preserves the established ordering
/// across connected components.
pub(in crate::kernels) fn reverse_cuthill_mckee(adj: &[Vec<usize>]) -> Vec<usize> {
    let n = adj.len();
    if adj.iter().all(Vec::is_empty) {
        return (0..n).rev().collect();
    }
    let degree: Vec<_> = adj.iter().map(Vec::len).collect();
    let mut numbered = vec![false; n];
    let mut probe_seen = vec![false; n];
    let mut probe_queue = Vec::new();
    let mut order = Vec::with_capacity(n);
    for seed in 0..n {
        if numbered[seed] {
            continue;
        }
        if degree[seed] == 0 {
            order.push(seed);
            continue;
        }

        let mut start = seed;
        let (mut depth, mut last_layer) = probe(adj, start, &mut probe_seen, &mut probe_queue);
        for _ in 0..16 {
            // A tie goes to the first minimum-degree vertex in BFS order.
            let candidate = probe_queue[last_layer..]
                .iter()
                .copied()
                .min_by_key(|&v| degree[v])
                .unwrap_or(start);
            let (candidate_depth, candidate_layer) =
                probe(adj, candidate, &mut probe_seen, &mut probe_queue);
            if candidate_depth <= depth {
                break;
            }
            start = candidate;
            depth = candidate_depth;
            last_layer = candidate_layer;
        }

        // The growing numbering doubles as the Cuthill-McKee queue.
        let mut head = order.len();
        order.push(start);
        numbered[start] = true;
        while head < order.len() {
            let vertex = order[head];
            head += 1;
            let tail = order.len();
            for &neighbor in &adj[vertex] {
                if !numbered[neighbor] {
                    numbered[neighbor] = true;
                    order.push(neighbor);
                }
            }
            // Sorting by the neighbour ID breaks ties the same way as a
            // stable degree sort of this row's ascending adjacency.
            order[tail..].sort_unstable_by_key(|&v| (degree[v], v));
        }
    }
    order.reverse();
    order
}

/// Reuses the previous BFS queue and returns the depth and starting position
/// of the new last layer. Each probe clears marks only in its previous queue.
fn probe(
    adj: &[Vec<usize>],
    root: usize,
    seen: &mut [bool],
    queue: &mut Vec<usize>,
) -> (usize, usize) {
    for &v in queue.iter() {
        seen[v] = false;
    }
    queue.clear();
    queue.push(root);
    seen[root] = true;

    let mut depth = 0;
    let mut begin = 0;
    loop {
        let end = queue.len();
        for i in begin..end {
            let vertex = queue[i];
            for &neighbor in &adj[vertex] {
                if !seen[neighbor] {
                    seen[neighbor] = true;
                    queue.push(neighbor);
                }
            }
        }
        if queue.len() == end {
            return (depth, begin);
        }
        begin = end;
        depth += 1;
    }
}
