use proptest::prelude::*;
use rssn::numerical::signal::*;
use rustfft::num_complex::Complex;

#[test]

fn test_fft_basic() {

    let mut input = vec![
        Complex::new(1.0, 0.0),
        Complex::new(1.0, 0.0),
        Complex::new(1.0, 0.0),
        Complex::new(1.0, 0.0),
    ];

    let output = fft(&mut input);

    assert!(
        (output[0].re - 4.0).abs()
            < 1e-9
    );

    assert!(
        (output[1].re).abs() < 1e-9
    );
}

#[test]

fn test_convolve_basic() {

    let a = vec![1.0, 2.0, 3.0];

    let v = vec![0.0, 1.0, 0.5];

    let result = convolve(&a, &v);

    assert_eq!(
        result,
        vec![0.0, 1.0, 2.5, 4.0, 1.5]
    );
}

#[test]

fn test_cross_correlation_basic() {

    let a = vec![1.0, 2.0, 3.0];

    let v = vec![0.0, 1.0, 0.5];

    let result =
        cross_correlation(&a, &v);

    // v_rev = [0.5, 1.0, 0.0]
    // a * v_rev = [1, 2, 3] * [0.5, 1, 0]
    // k=0: 1*0.5 = 0.5
    // k=1: 1*1 + 2*0.5 = 2.0
    // k=2: 1*0 + 2*1 + 3*0.5 = 3.5
    // k=3: 2*0 + 3*1 = 3.0
    // k=4: 3*0 = 0.0
    assert_eq!(
        result,
        vec![0.5, 2.0, 3.5, 3.0, 0.0]
    );
}

#[test]

fn test_windows() {

    let hann = hann_window(5);

    assert!((hann[0]).abs() < 1e-9);

    assert!(
        (hann[2] - 1.0).abs() < 1e-9
    );

    assert!((hann[4]).abs() < 1e-9);

    let hamming = hamming_window(5);

    assert!(
        (hamming[0] - 0.08).abs()
            < 1e-9
    );

    assert!(
        (hamming[2] - 1.0).abs() < 1e-9
    );

    assert!(
        (hamming[4] - 0.08).abs()
            < 1e-9
    );
}

proptest! {
    #[test]
    fn proptest_convolve_identity(a in prop::collection::vec(-10.0..10.0f64, 1..10)) {
        let v = vec![1.0];
        let res = convolve(&a, &v);
        prop_assert_eq!(res, a);
    }
}
