use num_complex::Complex;
use proptest::prelude::*;
use rssn::numerical::transforms::*;

#[test]

fn test_fft_basic() {

    let mut data = vec![
        Complex::new(1.0, 0.0),
        Complex::new(1.0, 0.0),
        Complex::new(0.0, 0.0),
        Complex::new(0.0, 0.0),
    ];

    fft(&mut data);

    // Expected: [2, 1-i, 0, 1+i]
    assert!((data[0].re - 2.0).abs() < 1e-9);

    assert!((data[1].re - 1.0).abs() < 1e-9);

    assert!((data[1].im + 1.0).abs() < 1e-9);

    assert!((data[2].re - 0.0).abs() < 1e-9);

    assert!((data[3].re - 1.0).abs() < 1e-9);

    assert!((data[3].im - 1.0).abs() < 1e-9);
}

#[test]

fn test_ifft_basic() {

    let mut data = vec![
        Complex::new(2.0, 0.0),
        Complex::new(1.0, -1.0),
        Complex::new(0.0, 0.0),
        Complex::new(1.0, 1.0),
    ];

    ifft(&mut data);

    assert!((data[0].re - 1.0).abs() < 1e-9);

    assert!((data[1].re - 1.0).abs() < 1e-9);

    assert!((data[2].re - 0.0).abs() < 1e-9);

    assert!((data[3].re - 0.0).abs() < 1e-9);
}

#[test]

fn test_fft_ifft_roundtrip() {

    let mut data = vec![
        Complex::new(1.0, 2.0),
        Complex::new(3.0, 4.0),
        Complex::new(5.0, 6.0),
        Complex::new(7.0, 8.0),
    ];

    let original = data.clone();

    fft(&mut data);

    ifft(&mut data);

    for i in 0 .. 4 {

        assert!((data[i].re - original[i].re).abs() < 1e-9);

        assert!((data[i].im - original[i].im).abs() < 1e-9);
    }
}

proptest! {
    #[test]
    fn proptest_fft_ifft_roundtrip(re in prop::collection::vec(-100.0..100.0f64, 16), im in prop::collection::vec(-100.0..100.0f64, 16)) {
        let mut data: Vec<Complex<f64>> = re.iter().zip(im.iter()).map(|(&r, &i)| Complex::new(r, i)).collect();
        let original = data.clone();
        fft(&mut data);
        ifft(&mut data);
        for i in 0..16 {
            prop_assert!((data[i].re - original[i].re).abs() < 1e-7);
            prop_assert!((data[i].im - original[i].im).abs() < 1e-7);
        }
    }
}
