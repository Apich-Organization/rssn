//! Seeded and thermodynamic checks for the Potts heat-bath conditional.

use rssn::sim::spin_models::Potts;
use rssn::sim::spin_models::PottsHeatBathWorkspace;

fn legacy_sweep<const D: usize>(model: &mut Potts<D>) {
    let mut weights = vec![0.0; model.q];
    for site in 0..model.spins.len() {
        let neighbors = model.lattice.neighbours(site);
        let mut total = 0.0;
        for (color, weight) in weights.iter_mut().enumerate() {
            let count = neighbors
                .iter()
                .filter(|&&n| n != usize::MAX && model.spins[n] == color)
                .count();
            *weight = (model.beta * model.j * count as f64).exp();
            total += *weight;
        }
        let mut u = model.rng.uniform() * total;
        let mut color = model.q - 1;
        for (candidate, &weight) in weights.iter().enumerate() {
            if u < weight {
                color = candidate;
                break;
            }
            u -= weight;
        }
        model.spins[site] = color;
    }
}

fn check_legacy_trace<const D: usize>(
    workspace: &mut PottsHeatBathWorkspace,
    dims: [usize; D],
    periodic: bool,
    q: usize,
) {
    let mut model = Potts::new(dims, periodic, q, 0.9, 523).unwrap();
    let mut reference = model.clone();
    for sweep in 0..30 {
        // The same workspace also serves models whose parameters are changed
        // through the public fields between sweeps.
        let (beta, j) = match sweep {
            | 0..=4 => (0.9, 0.75),
            | 5..=10 => (0.9, -1.0),
            | 11..=19 => (1.7, -1.0),
            | _ => (1.7, 0.0),
        };
        model.beta = beta;
        model.j = j;
        reference.beta = beta;
        reference.j = j;
        legacy_sweep(&mut reference);
        if sweep % 5 == 0 {
            model.heat_bath_sweep();
        } else {
            model.heat_bath_sweep_with(workspace);
        }
        assert_eq!(model.spins, reference.spins);
        assert_eq!(
            model.rng.clone().next_u64(),
            reference.rng.clone().next_u64()
        );
        assert_eq!(model.energy().to_bits(), reference.energy().to_bits());
    }
}

#[test]
fn workspace_keeps_seeded_color_order_across_couplings_and_geometries() {
    let mut workspace = PottsHeatBathWorkspace::default();
    check_legacy_trace(&mut workspace, [17], false, 2);
    check_legacy_trace(&mut workspace, [10, 13], true, 3);
    check_legacy_trace(&mut workspace, [7, 11, 5], false, 41);
    check_legacy_trace(&mut workspace, [5, 7], true, 127);
    check_legacy_trace(&mut workspace, [1; 9], false, 7);
    check_legacy_trace(&mut workspace, [21], true, 2);
}

#[test]
fn isolated_site_keeps_uniform_color_order_at_all_defined_couplings() {
    let mut workspace = PottsHeatBathWorkspace::default();
    for q in [2, 3, 31, 64] {
        let mut finite = Potts::new([1], false, q, 1.3, 78).unwrap();
        let mut infinite = finite.clone();
        infinite.beta = f64::INFINITY;
        let mut reference = finite.clone();
        for _ in 0..128 {
            legacy_sweep(&mut reference);
            finite.heat_bath_sweep();
            infinite.heat_bath_sweep_with(&mut workspace);
            assert_eq!(finite.spins, reference.spins);
            assert_eq!(infinite.spins, reference.spins);
        }
        let next = reference.rng.next_u64();
        assert_eq!(finite.rng.next_u64(), next);
        assert_eq!(infinite.rng.next_u64(), next);
    }
}

fn conditional_reference<const D: usize>(
    model: &Potts<D>,
    site: usize,
) -> Vec<f64> {
    // Enumerate the total bond energy for each candidate color independently
    // of the conditional's neighbor histogram.
    let mut sample = model.clone();
    let energies: Vec<_> = (0..model.q)
        .map(|color| {
            sample.spins[site] = color;
            sample.energy()
        })
        .collect();
    let minimum = energies.iter().copied().fold(f64::INFINITY, f64::min);
    let weights: Vec<_> = energies
        .iter()
        .map(|&energy| {
            let delta = energy - minimum;
            if delta == 0.0 {
                1.0
            } else {
                (-model.beta * delta).exp()
            }
        })
        .collect();
    let total: f64 = weights.iter().sum();
    weights.into_iter().map(|w| w / total).collect()
}

#[test]
fn conditional_probabilities_agree_with_independent_bond_energy() {
    for periodic in [false, true] {
        for (beta, j) in [
            (0.0, 2.0),
            (0.9, 1.3),
            (0.9, -1.3),
            (1e300, 1.0),
            (1e300, -1.0),
            (f64::INFINITY, 1.0),
            (f64::INFINITY, -1.0),
            (f64::INFINITY, 0.0),
        ] {
            let mut model = Potts::new([3, 4], periodic, 5, beta, 135).unwrap();
            model.j = j;
            let next = model.rng.clone().next_u64();
            for site in [0, 1, 4, 8, 11] {
                let mut p = [0.0; 5];
                model.heat_bath_probabilities(site, &mut p).unwrap();
                let reference = conditional_reference(&model, site);
                for (got, expected) in p.into_iter().zip(reference) {
                    assert!(
                        (got - expected).abs() < 4e-15,
                        "site {site}: {got} vs {expected}"
                    );
                }
                assert!((p.iter().sum::<f64>() - 1.0).abs() < 2e-15);
            }
            assert_eq!(next, model.rng.clone().next_u64());
        }
    }
}

#[test]
fn extreme_sweeps_keep_allowed_colors_and_one_draw_per_site() {
    // The ordinary sum is already infinite at this finite coupling,
    // although each of the two nonzero-count colors has a finite weight.
    let mut overflow = Potts::new([3], false, 3, 709.5, 2).unwrap();
    overflow.spins = vec![0, 2, 1];
    let mut p = [0.0; 3];
    overflow.heat_bath_probabilities(1, &mut p).unwrap();
    assert_eq!(p[0].to_bits(), 0.5_f64.to_bits());
    assert_eq!(p[1].to_bits(), 0.5_f64.to_bits());
    assert!(p[2] > 0.0 && p[2].is_subnormal());
    overflow.beta = f64::INFINITY;
    overflow.heat_bath_probabilities(1, &mut p).unwrap();
    assert_eq!(p.map(f64::to_bits), [0.5, 0.5, 0.0].map(f64::to_bits));

    let mut ferromagnetic = Potts::new([2], false, 3, 800.0, 9).unwrap();
    ferromagnetic.spins = vec![2, 0];
    let mut rng = ferromagnetic.rng.clone();
    let _ = rng.uniform();
    let _ = rng.uniform();
    ferromagnetic.heat_bath_sweep();
    assert_eq!(ferromagnetic.spins, [0, 0]);
    assert_eq!(ferromagnetic.rng.next_u64(), rng.next_u64());

    let mut antiferromagnetic = Potts::new([3], false, 2, 800.0, 18).unwrap();
    antiferromagnetic.j = -1.0;
    antiferromagnetic.spins = vec![0, 1, 1];
    let mut tie = [0.0; 2];
    antiferromagnetic
        .heat_bath_probabilities(1, &mut tie)
        .unwrap();
    assert_eq!(tie.map(f64::to_bits), [0.5, 0.5].map(f64::to_bits));
    // At the two open end sites the least-matching color is unique.
    antiferromagnetic
        .heat_bath_probabilities(0, &mut tie)
        .unwrap();
    assert_eq!(tie.map(f64::to_bits), [1.0, 0.0].map(f64::to_bits));

    let initial = antiferromagnetic.spins.clone();
    let mut workspace = PottsHeatBathWorkspace::default();
    let mut observed_zero = false;
    for seed in 0..20 {
        antiferromagnetic.spins.copy_from_slice(&initial);
        let mut queue = Potts::new([3], false, 2, 0.0, seed).unwrap().rng;
        antiferromagnetic.rng = queue.clone();
        let _ = queue.uniform();
        let middle = usize::from(queue.uniform() >= 0.5);
        let _ = queue.uniform();
        observed_zero |= middle == 0;
        antiferromagnetic.heat_bath_sweep_with(&mut workspace);
        assert_eq!(antiferromagnetic.spins, [0, middle, 1 - middle]);
        assert_eq!(antiferromagnetic.rng.next_u64(), queue.next_u64());
    }
    assert!(observed_zero);
}

#[test]
fn invalid_conditionals_leave_the_caller_output_alone() {
    let mut model = Potts::new([5], true, 3, 0.7, 11).unwrap();
    let mut output = [8.0; 3];
    assert!(model.heat_bath_probabilities(5, &mut output).is_err());
    assert!(model.heat_bath_probabilities(2, &mut output[..2]).is_err());
    model.beta = f64::NAN;
    assert!(model.heat_bath_probabilities(2, &mut output).is_err());
    model.beta = 0.7;
    model.spins[3] = 3;
    assert!(model.heat_bath_probabilities(2, &mut output).is_err());
    assert_eq!(output.map(f64::to_bits), [8.0; 3].map(f64::to_bits));
}
