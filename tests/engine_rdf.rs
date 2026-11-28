//! Radial distribution histograms from multiple particle and volume ensembles.

use rssn::sim::md_engine::MdSystem;
use rssn::sim::md_engine::Rdf;
use rssn::sim::md_engine::SimBox;

fn line(
    length: f64,
    points: &[f64],
) -> MdSystem<1> {
    MdSystem::new(
        points.iter().map(|&x| [x]).collect(),
        SimBox::periodic([length]),
    )
}

fn bits(
    actual: f64,
    expected: f64,
) {
    assert_eq!(actual.to_bits(), expected.to_bits());
}

#[test]
fn strict_cutoff_keeps_the_inside_edge_in_the_last_bin() {
    let cutoff = 0.1_f64;
    let mut rdf = Rdf::<1>::new(3, cutoff);
    let below = cutoff.next_down();
    assert_eq!((below / (cutoff / 3.0)).to_bits(), 3.0_f64.to_bits());
    rdf.add_frame(&line(4.0, &[0.0, below]));
    rdf.add_frame(&line(4.0, &[0.0, cutoff]));
    rdf.add_frame(&line(4.0, &[0.0, cutoff.next_up()]));
    assert_eq!(rdf.hist, [0.0, 0.0, 2.0]);

    let mut empty_bins = Rdf::<1>::new(0, cutoff);
    empty_bins.add_frame(&line(4.0, &[0.0, below]));
    assert!(empty_bins.result().is_empty());
}

#[test]
fn frame_normalization_uses_its_own_population_and_volume() {
    let a = line(4.0, &[0.0, 0.25]);
    let b = line(8.0, &[0.0, 0.25, 3.5]);
    let c = line(8.0, &[0.0, 0.25, 2.5, 5.6]);
    let empty = line(5.0, &[]);
    // Each frame has one pair inside the shell of length 2. Its histogram
    // contributes 2 and its reference count is 2 * N² / length.
    for (x, y, expected) in [(&a, &b, 16.0 / 17.0), (&a, &c, 2.0 / 3.0)] {
        for (first, second) in [(x, y), (y, x)] {
            let mut hist = Rdf::<1>::new(1, 1.0);
            hist.add_frame(&empty);
            assert!(hist.result()[0].1.is_nan());
            hist.add_frame(first);
            hist.add_frame(second);
            hist.add_frame(&empty);
            assert_eq!(hist.hist, [4.0]);
            bits(hist.result()[0].1, expected);
        }
    }
}

#[test]
fn enumerated_ideal_ensembles_have_the_expected_finite_population_factor() {
    // On integer sites, the shell of length one at radius 0.5 contains
    // exactly the coincident site. Thus its pair probability is 1 / V.
    fn populate(
        rdf: &mut Rdf<1>,
        sites: usize,
        particles: usize,
    ) {
        for mut code in 0..sites.pow(particles as u32) {
            let points: Vec<_> = (0..particles)
                .map(|_| {
                    let site = code % sites;
                    code /= sites;
                    [site as f64]
                })
                .collect();
            rdf.add_frame(&MdSystem::new(points, SimBox::periodic([sites as f64])));
        }
    }
    let mut two = Rdf::<1>::new(1, 0.5);
    let mut three = two.clone();
    populate(&mut two, 4, 2);
    populate(&mut three, 2, 3);
    assert_eq!(two.hist, [8.0]);
    assert_eq!(three.hist, [24.0]);
    bits(two.result()[0].1, 0.5);
    bits(three.result()[0].1, 2.0 / 3.0);
    two.merge(&three).unwrap();
    bits(two.result()[0].1, 8.0 / 13.0);
}

#[test]
fn merging_fixed_ensembles_keeps_the_original_arithmetic() {
    let frame = line(10.0, &[0.0, 0.125, 0.75, 4.0]);
    let mut independent = Rdf::<1>::new(2, 1.0);
    let mut sequential = independent.clone();
    let mut other = independent.clone();
    for i in 0..7 {
        sequential.add_frame(&frame);
        if i < 3 {
            independent.add_frame(&frame);
        } else {
            other.add_frame(&frame);
        }
    }
    independent.merge(&other).unwrap();
    let denominator = 7.0 * 4.0 * (4.0 / 10.0);
    for rdf in [independent, sequential] {
        assert_eq!(rdf.hist, [14.0, 28.0]);
        let ratios = rdf.result();
        bits(ratios[0].0, 0.25);
        bits(ratios[1].0, 0.75);
        bits(ratios[0].1, 14.0 / denominator);
        bits(ratios[1].1, 28.0 / denominator);
    }
}

#[test]
fn merging_variable_ensembles_keeps_both_completed_and_open_weights() {
    let a = line(4.0, &[0.0, 0.25]);
    let b = line(8.0, &[0.0, 0.25, 3.5]);
    let c = line(8.0, &[0.0, 0.25, 2.5, 5.6]);
    let mut first = Rdf::<1>::new(1, 1.0);
    first.add_frame(&a);
    first.add_frame(&b);
    let mut second = Rdf::<1>::new(1, 1.0);
    second.add_frame(&c);
    second.add_frame(&a);
    first.merge(&second).unwrap();
    first.merge(&Rdf::<1>::new(1, 1.0)).unwrap();
    first.add_frame(&b);
    assert_eq!(first.hist, [10.0]);
    bits(first.result()[0].1, 4.0 / 5.0);
}

#[test]
fn merge_errors_leave_the_existing_frame_counts_and_weights_in_place() {
    let mut left = Rdf::<1>::new(2, 1.0);
    left.add_frame(&line(4.0, &[0.0, 0.25]));
    for other in [Rdf::<1>::new(3, 1.0), Rdf::<1>::new(2, 0.9)] {
        assert!(left.merge(&other).is_err());
        assert_eq!(left.hist, [2.0, 0.0]);
        bits(left.result()[0].1, 2.0);
    }
    left.hist[0] = f64::MAX;
    let huge = left.clone();
    assert!(left.merge(&huge).is_err());
    assert_eq!(left.hist, huge.hist);
    assert_eq!(left.result(), huge.result());
}
