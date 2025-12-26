use std::ffi::CStr;
use std::ffi::CString;

use rssn::ffi_apis::common::from_bincode_buffer;
use rssn::ffi_apis::common::rssn_free_bincode_buffer;
use rssn::ffi_apis::common::rssn_free_string;
use rssn::ffi_apis::common::to_bincode_buffer;
use rssn::ffi_apis::common::BincodeBuffer;
use rssn::ffi_apis::ffi_api::FfiResult;
use rssn::ffi_apis::numerical_graph_ffi::bincode_api;
use rssn::ffi_apis::numerical_graph_ffi::handle;
use rssn::ffi_apis::numerical_graph_ffi::json;
use serde::Deserialize;
use serde::Serialize;

#[test]

fn test_graph_handle_ffi() {

    unsafe {

        let n = 3;

        let graph = handle::rssn_num_graph_create(n);

        assert!(!graph.is_null());

        handle::rssn_num_graph_add_edge(graph, 0, 1, 1.0);

        handle::rssn_num_graph_add_edge(graph, 1, 2, 2.0);

        let mut dist = vec![0.0; n];

        let mut prev = vec![0isize; n];

        // Dijkstra from 0
        handle::rssn_num_graph_dijkstra(
            graph,
            0,
            dist.as_mut_ptr(),
            prev.as_mut_ptr(),
        );

        assert_eq!(dist[0], 0.0);

        assert_eq!(dist[1], 1.0);

        assert_eq!(dist[2], 3.0);

        assert_eq!(prev[0], -1);

        assert_eq!(prev[1], 0);

        assert_eq!(prev[2], 1);

        // BFS from 0
        let mut bfs_dist = vec![0usize; n];

        handle::rssn_num_graph_bfs(
            graph,
            0,
            bfs_dist.as_mut_ptr(),
        );

        assert_eq!(bfs_dist[0], 0);

        assert_eq!(bfs_dist[1], 1);

        assert_eq!(bfs_dist[2], 2);

        // PageRank
        let mut scores = vec![0.0; n];

        handle::rssn_num_graph_page_rank(
            graph,
            0.85,
            1e-6,
            100,
            scores.as_mut_ptr(),
        );

        assert!(
            scores
                .iter()
                .sum::<f64>()
                - 1.0
                < 1e-6
        );

        // Floyd-Warshall
        let mut fw_dist = vec![0.0; n * n];

        handle::rssn_num_graph_floyd_warshall(
            graph,
            fw_dist.as_mut_ptr(),
        );

        assert_eq!(
            fw_dist[0 * n + 2],
            3.0
        );

        // Connected Components
        let mut components = vec![0usize; n];

        handle::rssn_num_graph_connected_components(
            graph,
            components.as_mut_ptr(),
        );

        // 0->1->2 connected
        assert_eq!(
            components[0],
            components[1]
        );

        assert_eq!(
            components[1],
            components[2]
        );

        // MST
        // Current graph 0->1(1), 1->2(2). It is already a tree if undirected?
        // But function treats directed graph. Prim's usually for undirected.
        // Let's make it undirected-ish by adding reverse edges for MST test context or assume it works on what we gave.
        // For Prim's on directed 0->1, 1->2 starting at 0:
        // 0 in tree. edge 0->1(1). pick 1. edge 1->2(2). pick 2.
        // MST has 0->1 and 1->2. (and reverse in code logic? No, code adds p->u and u->p)

        let mst_handle = handle::rssn_num_graph_minimum_spanning_tree(graph);

        assert!(!mst_handle.is_null());

        // Since it's opaque handle, we can't inspect easily without accessor or assumption.
        // We can just free it.
        handle::rssn_num_graph_free(mst_handle);

        handle::rssn_num_graph_free(graph);
    }
}

#[derive(Serialize)]

struct Edge {
    u : usize,
    v : usize,
    w : f64,
}

#[derive(Serialize)]

struct GraphDef {
    num_nodes : usize,
    edges : Vec<Edge>,
}

#[derive(Serialize)]

struct DijkstraInput {
    graph : GraphDef,
    start_node : usize,
}

#[derive(Serialize)]

struct PageRankInput {
    graph : GraphDef,
    damping_factor : f64,
    tolerance : f64,
    max_iter : usize,
}

#[test]

fn test_graph_json_ffi() {

    unsafe {

        let graph = GraphDef {
            num_nodes : 3,
            edges : vec![
                Edge {
                    u : 0,
                    v : 1,
                    w : 1.0,
                },
                Edge {
                    u : 1,
                    v : 2,
                    w : 2.0,
                },
            ],
        };

        let input = DijkstraInput {
            graph,
            start_node : 0,
        };

        let json_str = serde_json::to_string(&input).unwrap();

        let c_json = CString::new(json_str).unwrap();

        // Dijkstra
        let res_ptr = json::rssn_num_graph_dijkstra_json(c_json.as_ptr());

        let res_str = CStr::from_ptr(res_ptr)
            .to_str()
            .unwrap();

        let v : serde_json::Value = serde_json::from_str(res_str).unwrap();

        let dist = v["ok"]["dist"]
            .as_array()
            .unwrap();

        assert_eq!(
            dist[2]
                .as_f64()
                .unwrap(),
            3.0
        );

        rssn_free_string(res_ptr);

        // BFS
        let res_ptr = json::rssn_num_graph_bfs_json(c_json.as_ptr());

        let res_str = CStr::from_ptr(res_ptr)
            .to_str()
            .unwrap();

        let v : serde_json::Value = serde_json::from_str(res_str).unwrap();

        let bfs_dist = v["ok"]
            .as_array()
            .unwrap();

        assert_eq!(
            bfs_dist[2]
                .as_u64()
                .unwrap(),
            2
        );

        rssn_free_string(res_ptr);

        // Connected Components - needs GraphDef input, not DijkstraInput
        let graph_only = GraphDef {
            num_nodes : 3,
            edges : vec![
                Edge {
                    u : 0,
                    v : 1,
                    w : 1.0,
                },
                Edge {
                    u : 1,
                    v : 2,
                    w : 2.0,
                },
            ],
        };

        let graph_json_str = serde_json::to_string(&graph_only).unwrap();

        let c_graph_json = CString::new(graph_json_str).unwrap();

        let res_ptr = json::rssn_num_graph_connected_components_json(c_graph_json.as_ptr());

        let res_str = CStr::from_ptr(res_ptr)
            .to_str()
            .unwrap();

        let v : serde_json::Value = serde_json::from_str(res_str).unwrap();

        let comp = v["ok"]
            .as_array()
            .unwrap();

        // 0->1->2 connected if we consider connectivity loosely or bfs from all unvisited
        assert_eq!(comp[0], comp[1]);

        assert_eq!(comp[1], comp[2]);

        rssn_free_string(res_ptr);

        // MST
        let res_ptr = json::rssn_num_graph_minimum_spanning_tree_json(c_graph_json.as_ptr());

        let res_str = CStr::from_ptr(res_ptr)
            .to_str()
            .unwrap();

        let v : serde_json::Value = serde_json::from_str(res_str).unwrap();

        let edges = v["ok"]["edges"]
            .as_array()
            .unwrap();

        // 0->1, 1->2 plus reverse 1->0, 2->1 = 4 edges
        assert_eq!(edges.len(), 4);

        rssn_free_string(res_ptr);
    }
}

#[test]

fn test_graph_bincode_ffi() {

    unsafe {

        let graph = GraphDef {
            num_nodes : 3,
            edges : vec![
                Edge {
                    u : 0,
                    v : 1,
                    w : 1.0,
                },
                Edge {
                    u : 1,
                    v : 2,
                    w : 2.0,
                },
            ],
        };

        let input = DijkstraInput {
            graph,
            start_node : 0,
        };

        let buffer = to_bincode_buffer(&input);

        // Dijkstra
        let res_buffer = bincode_api::rssn_num_graph_dijkstra_bincode(buffer);

        #[derive(Deserialize)]

        struct DijkstraOutput {
            dist : Vec<f64>,
            prev : Vec<Option<usize>>,
        }

        let res : FfiResult<DijkstraOutput, String> = from_bincode_buffer(&res_buffer).unwrap();

        let out = res.ok.unwrap();

        assert_eq!(out.dist[2], 3.0);

        rssn_free_bincode_buffer(res_buffer);

        rssn_free_bincode_buffer(buffer);

        // Connected Components
        #[derive(Serialize)]

        struct GraphDefIn {
            // Need separate struct if reusing fields not optimal, but GraphDef is fine
            num_nodes : usize,
            edges : Vec<Edge>,
        }

        let input_graph = GraphDefIn {
            num_nodes : 3,
            edges : vec![
                Edge {
                    u : 0,
                    v : 1,
                    w : 1.0,
                },
                Edge {
                    u : 1,
                    v : 2,
                    w : 2.0,
                },
            ],
        };

        let buffer_graph = to_bincode_buffer(&input_graph);

        let res_buffer = bincode_api::rssn_num_graph_connected_components_bincode(buffer_graph);

        let res : FfiResult<Vec<usize>, String> = from_bincode_buffer(&res_buffer).unwrap();

        let comp = res.ok.unwrap();

        assert_eq!(comp[0], comp[1]);

        rssn_free_bincode_buffer(res_buffer);

        // MST
        let res_buffer = bincode_api::rssn_num_graph_minimum_spanning_tree_bincode(buffer_graph);

        #[derive(Deserialize)]

        struct EdgeOut {
            u : usize,
            v : usize,
            w : f64,
        }

        #[derive(Deserialize)]

        struct GraphDefOut {
            num_nodes : usize,
            edges : Vec<EdgeOut>,
        }

        let res : FfiResult<GraphDefOut, String> = from_bincode_buffer(&res_buffer).unwrap();

        let mst = res.ok.unwrap();

        assert_eq!(mst.edges.len(), 4);

        rssn_free_bincode_buffer(res_buffer);

        rssn_free_bincode_buffer(buffer_graph);
    }
}
