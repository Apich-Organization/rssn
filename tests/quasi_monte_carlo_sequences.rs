//! Independent bit pattern checks for Sobol seek and fixed-block quadrature.
use rssn::kernels::quadrature::LowDiscrepancy;
use rssn::kernels::quadrature::Sobol;
use rssn::kernels::quadrature::halton_point;
use rssn::kernels::quadrature::quasi_monte_carlo;
use rssn::kernels::quadrature::quasi_monte_carlo_par;

fn first_coordinate(index: u32) -> f64 {
    f64::from((index ^ (index >> 1)).reverse_bits()) * (1.0 / 4_294_967_296.0)
}

#[test]
fn sobol_random_access_and_wrap_match_gray_code() {
    let mut sequence = Sobol::new(2);
    let first_points = [[0.5, 0.5], [0.75, 0.25], [0.25, 0.75], [0.375, 0.375]];
    for (index, expected) in (1..).zip(first_points) {
        assert_eq!(sequence.point_at(index), expected);
        assert_eq!(sequence.next_point(), expected);
    }
    for index in [0, 2048, (1 << 31) - 1, 1 << 31, u32::MAX - 1, u32::MAX] {
        assert_eq!(sequence.point_at(index)[0], first_coordinate(index));
        sequence.seek(index);
        let next = index.wrapping_add(1);
        assert_eq!(sequence.next_point(), sequence.point_at(next));
        assert_eq!(sequence.point_at(next)[0], first_coordinate(next));
    }
    sequence.seek(u32::MAX);
    assert_eq!(sequence.next_point(), [0.0; 2]);
    assert_eq!(sequence.next_point(), first_points[0]);
}

#[test]
fn parallel_sobol_keeps_fixed_block_bits_and_exact_sample_order() {
    let one = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let two = rayon::ThreadPoolBuilder::new()
        .num_threads(2)
        .build()
        .unwrap();
    for samples in [1, 4095, 4096, 4097, 8192, 8193] {
        let expected_terms: Vec<_> = (1..=samples as u32)
            .map(|index| (2.0 + 3.0 * first_coordinate(index)).sin())
            .collect();
        let expected_parallel = expected_terms
            .chunks(4096)
            .map(|chunk| chunk.iter().sum::<f64>())
            .sum::<f64>();
        let expected_half = expected_terms[..samples / 2]
            .chunks(4096)
            .map(|chunk| chunk.iter().sum::<f64>())
            .sum::<f64>();
        let expected_serial = expected_terms.iter().sum::<f64>();
        let integrate = || {
            quasi_monte_carlo_par(
                |point| point[0].sin(),
                &[2.0],
                &[5.0],
                samples,
                LowDiscrepancy::Sobol,
            )
        };
        let result_one = one.install(integrate);
        let result_two = two.install(integrate);
        let serial = quasi_monte_carlo(
            |point| point[0].sin(),
            &[2.0],
            &[5.0],
            samples,
            LowDiscrepancy::Sobol,
        );
        let expected_value = 3.0 * expected_parallel / samples as f64;
        let expected_error =
            (expected_value - 3.0 * expected_half / (samples / 2).max(1) as f64).abs();
        assert_eq!(result_one.value.to_bits(), expected_value.to_bits());
        assert_eq!(result_one.error.to_bits(), expected_error.to_bits());
        assert_eq!(result_one, result_two);
        assert_eq!(
            serial.value.to_bits(),
            (3.0 * expected_serial / samples as f64).to_bits()
        );
    }
}

#[test]
fn invalid_quasi_monte_carlo_domains_do_not_call_the_integrand() {
    for (lo, hi, count, sequence) in [
        (&[0.0][..], &[][..], 1, LowDiscrepancy::Sobol),
        (&[0.0; 9][..], &[1.0; 9][..], 1, LowDiscrepancy::Sobol),
        (&[0.0; 17][..], &[1.0; 17][..], 1, LowDiscrepancy::Halton),
        (&[0.0][..], &[f64::INFINITY][..], 1, LowDiscrepancy::Halton),
        (&[0.0][..], &[1.0][..], 0, LowDiscrepancy::Sobol),
    ] {
        let untouched = |_: &[f64]| -> f64 { panic!("invalid domain called the integrand") };
        let serial = quasi_monte_carlo(untouched, lo, hi, count, sequence);
        let parallel = quasi_monte_carlo_par(untouched, lo, hi, count, sequence);
        for result in [serial, parallel] {
            assert_eq!(result.evaluations, 0);
            assert!(!result.converged);
            assert!(result.value.is_nan() && result.error.is_nan());
        }
    }
    let lo = [2.0, 3.0];
    let hi = [-1.0, 4.0];
    let result = quasi_monte_carlo_par(|_| 1.0, &lo, &hi, 8192, LowDiscrepancy::Halton);
    assert_eq!(result.value, -3.0);
    assert_eq!(result.error, 0.0);
}

#[test]
fn halton_fixed_bases_match_the_general_radical_inverse() {
    let primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53];
    for index in [0, 1, 32_768, (1_u64 << 47) - 1, u64::MAX - 1, u64::MAX] {
        let reference = primes
            .iter()
            .map(|&prime| {
                let mut absolute = index;
                let reciprocal = 1.0 / prime as f64;
                let mut place = reciprocal;
                let mut point = 0.0;
                while absolute != 0 {
                    point += (absolute % prime) as f64 * place;
                    absolute /= prime;
                    place *= reciprocal;
                }
                point.to_bits()
            })
            .collect::<Vec<_>>();
        let actual = halton_point(index, primes.len() + 1)
            .into_iter()
            .map(f64::to_bits)
            .collect::<Vec<_>>();
        assert_eq!(actual, reference);
    }
    assert!(halton_point(1, 0).is_empty());
}

#[test]
fn native_sobol_index_validates_before_changing_the_stream() {
    use rssn::Session;
    use rssn::kernel::CallCtx;
    use rssn::kernel::Value;
    let session = Session::new();
    let ctx = CallCtx::new();
    let call = |name, args: &[Value]| session.call_value(name, args, &ctx);
    let sequence = call("num_analysis_sobol_new", &[Value::int(2)]).unwrap();
    let point = "num_analysis_sobol_point_at";
    let seek = "num_analysis_sobol_seek";
    let next = "num_analysis_sobol_next";
    let at_one = call(point, &[sequence.clone(), Value::int(1)]).unwrap();
    let at_origin = call(point, &[sequence.clone(), Value::int(0)]).unwrap();
    assert!(
        call(
            seek,
            &[sequence.clone(), Value::int(u64::from(u32::MAX) + 1)]
        )
        .is_err()
    );
    assert_eq!(call(next, std::slice::from_ref(&sequence)).unwrap(), at_one);
    call(seek, &[sequence.clone(), Value::int(u32::MAX)]).unwrap();
    assert_eq!(
        call(next, std::slice::from_ref(&sequence)).unwrap(),
        at_origin
    );
    assert_eq!(call(next, &[sequence]).unwrap(), at_one);
}
