//! Equilibrium checks for fixed-work directed-loop sampling.

use rssn::sim::sse::SseResult;
use rssn::sim::sse::SseXxz;
use rssn::sim::sse::chain_bonds;
use rssn::sim::sse::square_bonds;

const SEEDS: [u64; 8] = [17, 386, 755, 1124, 1493, 1862, 2231, 2600];

fn observables(result: SseResult) -> [f64; 4] {
    [
        result.energy,
        result.specific_heat,
        result.magnetisation_sq,
        result.staggered_sq,
    ]
}

fn check_chains(
    n: usize,
    bonds: &[(usize, usize)],
    jz: f64,
    beta: f64,
    expected: [f64; 4],
) {
    let samples = SEEDS.map(|seed| {
        let mut sim = SseXxz::new(n, bonds.to_vec(), jz, beta, seed).unwrap();
        observables(sim.run(3_000, 65_000))
    });
    let absolute = [0.01, 0.045, 0.02, 0.04];
    let n_chains = SEEDS.len() as f64;
    for (column, target) in expected.into_iter().enumerate() {
        let mean = samples.iter().map(|row| row[column]).sum::<f64>() / n_chains;
        let variance = samples
            .iter()
            .map(|row| (row[column] - mean).powi(2))
            .sum::<f64>()
            / (n_chains - 1.0);
        let se = (variance / n_chains).sqrt();
        assert!(
            (mean - target).abs() <= absolute[column] + 6.0 * se,
            "n={n}, jz={jz}, beta={beta}, observable {column}: {mean} ± {se}, expected {target}"
        );
    }
}

#[test]
fn finite_temperature_dimer_matches_the_gibbs_distribution_in_multiple_chains() {
    let jz = 0.8;
    let energies = [0.2, 0.2, 0.3, -0.7];
    for beta in [0.4, 1.2, 2.0] {
        let weights = energies.map(|e: f64| (-beta * e).exp());
        let z = weights.iter().sum::<f64>();
        let energy = weights
            .iter()
            .zip(energies)
            .map(|(w, e)| w * e)
            .sum::<f64>()
            / z;
        let var = weights
            .iter()
            .zip(energies)
            .map(|(w, e)| w * (e - energy).powi(2))
            .sum::<f64>()
            / z;
        let m2 = (weights[0] + weights[1]) / z;
        check_chains(
            2,
            &[(0, 1)],
            jz,
            beta,
            [energy, beta * beta * var, m2, 1.0 - m2],
        );
    }
}

#[test]
fn six_site_lattices_match_exact_thermal_traces_for_independent_seeds() {
    // Dense eigenvalues and eigenvectors of the six-spin XXZ Hamiltonian.
    // Columns give energy, heat capacity, uniform M^2 and staggered M^2.
    let chain = chain_bonds(6, true);
    let rectangle = square_bonds(3, 2, false);
    for (bonds, jz, beta, exact) in [
        (
            &chain,
            1.0,
            2.0,
            [
                -2.101_699_844_213_744,
                2.367_814_436_975_792,
                0.417_083_474_872_750_3,
                3.624_354_410_224_42,
            ],
        ),
        (
            &chain,
            0.5,
            1.5,
            [
                -1.298_406_318_969_7,
                1.611_285_254_882_49,
                0.819_803_874_261_455_7,
                2.398_614_862_899_668,
            ],
        ),
        (
            &chain,
            0.0,
            2.0,
            [
                -1.247_089_631_483_617,
                1.755_730_417_127_727,
                0.996_472_933_543_880_1,
                2.003_527_066_456_121,
            ],
        ),
        (
            &rectangle,
            1.0,
            2.0,
            [
                -2.493_800_320_119_937,
                2.576_616_472_621_974,
                0.352_848_846_126_570,
                3.999_753_729_527_121,
            ],
        ),
    ] {
        check_chains(6, bonds, jz, beta, exact);
    }
}

#[test]
fn empty_operator_string_samples_free_spins() {
    let mut sim = SseXxz::new(2, vec![(0, 1)], 1.0, 1e-9, 13).unwrap();
    let results = sim.run(0, 8_192);
    assert_eq!(results.mean_ops.to_bits(), 0);
    assert!((results.magnetisation_sq - 0.5).abs() < 0.025);
    assert!((results.staggered_sq - 0.5).abs() < 0.025);
}

#[test]
fn independent_instances_reproduce_the_same_seed_and_split_runs() {
    for seed in SEEDS {
        let make = || SseXxz::new(6, chain_bonds(6, true), 0.8, 1.3, seed).unwrap();
        let (mut a, mut b) = (make(), make());
        for (therm, sweeps) in [(0, 0), (10, 1), (0, 39), (7, 404), (197, 700)] {
            let lhs = observables(a.run(therm, sweeps)).map(f64::to_bits);
            let rhs = observables(b.run(therm, sweeps)).map(f64::to_bits);
            assert_eq!(lhs, rhs);
            let (lhs, rhs) = (a.sample(), b.sample());
            assert_eq!(lhs.n_ops, rhs.n_ops);
            assert_eq!(lhs.magnetisation.to_bits(), rhs.magnetisation.to_bits());
            assert_eq!(lhs.staggered.to_bits(), rhs.staggered.to_bits());
        }
    }
}
