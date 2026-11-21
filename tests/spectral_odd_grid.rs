//! Native legacy spectral transport gives every odd/even transform coefficient a bin.
use rssn::kernels::fft::Complex64;
use rssn::sim::physics_sm::solve_advection_diffusion_1d;
#[test]
fn odd_five_mode_matches_native_forward_euler_amplification() {
    let n = 5;
    let dx = std::f64::consts::TAU / n as f64;
    let initial = (0..n).map(|i| (i as f64 * dx).cos()).collect::<Vec<_>>();
    let c = 0.3;
    let d = 0.1;
    let dt = 0.01;
    let steps = 7;
    let out = solve_advection_diffusion_1d(&initial, dx, c, d, dt, steps);
    let amplitude = Complex64::new(1.0 - d * dt, -c * dt).powu(steps as u32);
    for (i, value) in out.iter().enumerate() {
        let expected = (amplitude * Complex64::from_polar(1.0, i as f64 * dx)).re;
        assert!((value - expected).abs() < 1e-12);
    }
}
#[test]
fn one_point_constant_has_only_a_zero_wavenumber() {
    let out = solve_advection_diffusion_1d(&[2.5], 1.0, 10.0, 1.0, 0.5, 4);
    assert_eq!(out, vec![2.5]);
}
