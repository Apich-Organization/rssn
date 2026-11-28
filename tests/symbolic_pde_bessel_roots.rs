//! Bessel eigenspectra used by the symbolic PDE solvers.

use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Graph;
use std::f64::consts::PI;

fn graph() -> Graph {
    let mut graph = Graph::new();
    Engine::install(&mut graph, &[rssn::rules::pde()]).unwrap();
    graph
}

fn root<const N: usize>(
    graph: &mut Graph,
    operator: &str,
    args: [f64; N],
) -> f64 {
    let op = graph.ops().lookup(operator).unwrap();
    let nodes = args.map(|value| graph.float(value));
    let term = graph.node(op, &nodes);
    graph.eval(term, &Env::numeric(0.0)).unwrap()
}

#[test]
fn bessel_spectra_support_out_of_order_requests() {
    let mut graph = graph();
    let j0 = [
        2.404_825_557_695_773,
        5.520_078_110_286_311,
        8.653_727_912_911_013,
        11.791_534_439_014_28,
    ];
    let mut results = Vec::new();
    for i in [3_usize, 0, 1, 2, 3, 1] {
        let rank = f64::from(u32::try_from(i + 1).unwrap());
        let ordinary = root(&mut graph, "bessel_zero", [0.0, rank]);
        assert!((ordinary - j0[i]).abs() < 3e-10, "rank {rank}: {ordinary}");
        // J_{1/2}(x) is proportional to sin(x)/sqrt(x).
        let half_order = root(&mut graph, "bessel_zero", [0.5, rank]);
        assert!((half_order - rank * PI).abs() < 3e-10);
        // For ν=1/2 and h=-1/2, the root condition is x cos(x)=sin(x).
        let robin = root(&mut graph, "bessel_root", [0.5, -0.5, rank]);
        assert!(robin > rank * PI && robin < (rank + 0.5) * PI);
        assert!((robin * robin.cos() - robin.sin()).abs() < 3e-9);
        let bits = [ordinary.to_bits(), half_order.to_bits(), robin.to_bits()];
        if let Some((_, old)) = results.iter().find(|(position, _)| *position == i) {
            assert_eq!(&bits, old);
        }
        results.push((i, bits));
    }
}

#[test]
fn bessel_rank_validation_and_robin_zero_mode() {
    let mut graph = graph();
    for args in [[0.0, 0.0], [0.0, -1.0], [0.0, 2.5], [-1.0, 1.0]] {
        assert!(root(&mut graph, "bessel_zero", args).is_nan());
    }
    for args in [[0.0, 1.0, -1.0], [0.0, 1.0, 1.5]] {
        assert!(root(&mut graph, "bessel_root", args).is_nan());
    }
    assert_eq!(root(&mut graph, "bessel_root", [0.0, 0.5, 0.0]), 0.0);
}

#[test]
fn large_rank_bessel_zeros_remain_ordered() {
    let mut graph = graph();
    let mut previous = 0.0;
    for rank in [255.0, 256.0, 257.0, 258.0] {
        let zero = root(&mut graph, "bessel_zero", [0.0, rank]);
        assert!(zero > previous);
        assert!((zero - (rank - 0.25) * PI).abs() < 2e-4);
        previous = zero;
    }
}
