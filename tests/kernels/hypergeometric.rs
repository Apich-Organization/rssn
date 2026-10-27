//! Tests for complex-parameter pFq, 2F1, 1F1 and Tricomi U.

use num_complex::Complex64;
use rssn::kernels::hypergeometric::{
    gamma_c, hyp1f1_cplx, hyp2f1_cplx, hyp_pfq, hyp_u_cplx, ln_gamma_c,
};
use rssn::kernels::special::bessel_k;
use rssn::kernels::special_ext::{bessel_jn, hyp2f1_complex};

type C = Complex64;

fn c(re: f64, im: f64) -> C {
    C::new(re, im)
}

fn close(got: C, want: C, tol: f64, what: &str) {
    let d = (got - want).norm();
    assert!(d <= tol * want.norm().max(1e-300) || d < 1e-15, "{what}: {got} vs {want} (rel {})", d / want.norm());
}

#[test]
fn complex_gamma() {
    close(gamma_c(c(1.0, 1.0)), c(0.498_015_668_118_356, -0.154_949_828_301_810_6), 1e-13, "Gamma(1+i)");
    for y in [0.3, 1.7, 4.0] {
        let g = gamma_c(c(0.5, y));
        let want = std::f64::consts::PI / (std::f64::consts::PI * y).cosh();
        assert!((g.norm_sqr() - want).abs() < 1e-13 * want, "{y}");
    }
    // reflection and recurrence at points with negative real part
    let z = c(-2.3, 0.7);
    let refl = gamma_c(z) * gamma_c(c(1.0, 0.0) - z);
    close(refl, std::f64::consts::PI / (z * std::f64::consts::PI).sin(), 1e-13, "reflection");
    close(gamma_c(z + 1.0), z * gamma_c(z), 1e-13, "recurrence");
    assert!(gamma_c(c(-3.0, 0.0)).is_nan());
    close(ln_gamma_c(c(10.0, 0.0)), c(362_880.0_f64.ln(), 0.0), 1e-14, "ln 9!");
}

#[test]
fn hyp2f1_complex_parameters_closed_forms() {
    let a = c(0.7, 0.4);
    let b = c(-0.3, 1.1);
    for z in [c(0.3, 0.2), c(0.9, 0.3), c(-2.5, 1.5), c(2.0, 0.5), c(-30.0, -4.0), c(1.2, -0.01), c(0.6, -0.8)] {
        // 2F1(a, b; b; z) = (1 - z)^(-a)
        let got = hyp2f1_cplx(a, b, b, z);
        close(got, (c(1.0, 0.0) - z).powc(-a), 1e-11, &format!("(1-z)^-a at {z}"));
        // 2F1(a, 1; 2; z) = (1 - (1-z)^(1-a)) / ((1-a) z)
        let got = hyp2f1_cplx(a, c(1.0, 0.0), c(2.0, 0.0), z);
        let want = (c(1.0, 0.0) - (c(1.0, 0.0) - z).powc(c(1.0, 0.0) - a)) / ((c(1.0, 0.0) - a) * z);
        close(got, want, 1e-11, &format!("2F1(a,1;2) at {z}"));
        // 2F1(a, a + 1/2; 1/2; z^2) = ((1+z)^(-2a) + (1-z)^(-2a)) / 2 for |Re z| < 1
        if z.re.abs() < 1.0 {
            let w = z * z;
            let got = hyp2f1_cplx(a, a + 0.5, c(0.5, 0.0), w);
            let want = ((c(1.0, 0.0) + z).powc(-2.0 * a) + (c(1.0, 0.0) - z).powc(-2.0 * a)) * 0.5;
            close(got, want, 1e-10, &format!("quadratic form at {z}"));
        }
    }
    // on the cut (upper side): (1 - z)^(-a) with 1 - z = -2 - i0
    let z = c(3.0, 0.0);
    let got = hyp2f1_cplx(a, b, b, z);
    let want = (c(2.0, 0.0)).powc(-a) * (c(0.0, std::f64::consts::PI) * a).exp();
    close(got, want, 1e-11, "cut value");
}

#[test]
fn hyp2f1_matches_real_parameter_implementation() {
    let params = [(0.3, 0.8, 1.7), (1.5, -0.4, 2.2), (0.25, 0.5, 0.75), (2.0, 3.0, 4.5)];
    for (a, b, cc) in params {
        for z in [c(0.4, 0.4), c(-4.0, 0.0), c(3.0, 0.0), c(0.95, 0.1), c(1.5, -2.0), c(-0.99, 0.0), c(0.5, 1e-7)] {
            let got = hyp2f1_cplx(c(a, 0.0), c(b, 0.0), c(cc, 0.0), z);
            let want = hyp2f1_complex(a, b, cc, z);
            close(got, want, 1e-10, &format!("{a},{b},{cc} at {z}"));
        }
    }
}

#[test]
fn gauss_theorem_and_continuity_at_one() {
    let (a, b, cc) = (c(0.3, 0.2), c(0.4, -0.1), c(2.5, 0.3));
    let g = gamma_c(cc) * gamma_c(cc - a - b) / (gamma_c(cc - a) * gamma_c(cc - b));
    close(hyp2f1_cplx(a, b, cc, c(1.0, 0.0)), g, 1e-13, "Gauss");
    close(hyp2f1_cplx(a, b, cc, c(1.0 - 1e-7, 0.0)), g, 1e-5, "approach 1");
    close(hyp2f1_cplx(a, b, cc, c(1.0 - 1e-5, 1e-5)), g, 1e-3, "approach 1 (complex)");
}

#[test]
fn entire_cases_with_large_arguments() {
    // 0F1(; 1; -x^2/4) = J_0(x)
    for x in [0.5, 3.0, 12.0, 25.0, 40.0] {
        let got = hyp_pfq(&[], &[c(1.0, 0.0)], c(-x * x / 4.0, 0.0));
        let want = bessel_jn(0, x);
        assert!((got.re - want).abs() < 1e-11, "J0({x}): {got} vs {want}");
        assert!(got.im.abs() < 1e-11);
    }
    // 1F1(1; 2; z) = (e^z - 1) / z
    for z in [c(60.0, 0.0), c(-60.0, 0.0), c(30.0, 40.0), c(-25.0, 10.0), c(0.0, 50.0), c(3.0, -2.0)] {
        let got = hyp1f1_cplx(c(1.0, 0.0), c(2.0, 0.0), z);
        close(got, (z.exp() - 1.0) / z, 1e-11, &format!("M(1,2,{z})"));
    }
    // M(a, a, z) = e^z for complex a
    let a = c(1.3, -0.8);
    for z in [c(5.0, 3.0), c(-12.0, 4.0), c(0.2, 0.1)] {
        close(hyp1f1_cplx(a, a, z), z.exp(), 1e-11, "M(a,a,z)");
    }
    // M(1/2, 1, 2 i y) = e^{i y} J_0(y)
    for y in [0.7, 5.0, 30.0] {
        let got = hyp1f1_cplx(c(0.5, 0.0), c(1.0, 0.0), c(0.0, 2.0 * y));
        let want = c(0.0, y).exp() * bessel_jn(0, y);
        close(got, want, 1e-10, &format!("M(1/2,1,2i{y})"));
    }
    // sinh(z)/z = 0F1(;3/2; z^2/4)
    let z = c(7.0, 5.0);
    close(hyp_pfq(&[], &[c(1.5, 0.0)], z * z / 4.0), z.sinh() / z, 1e-11, "sinh");
}

fn li2_series(z: C) -> C {
    let mut s = c(0.0, 0.0);
    let mut zp = c(1.0, 0.0);
    for k in 1..6000 {
        zp *= z;
        s += zp / (f64::from(k) * f64::from(k));
    }
    s
}

#[test]
fn polylogarithms_by_pfq() {
    // Li2(z) / z = 3F2(1, 1, 1; 2, 2; z)
    let one = c(1.0, 0.0);
    let two = c(2.0, 0.0);
    let f = |z: C| hyp_pfq(&[one, one, one], &[two, two], z) * z;
    for z in [c(0.3, 0.1), c(0.9, 0.4), c(-0.7, 0.2), c(0.0, 0.9)] {
        close(f(z), li2_series(z), 1e-11, &format!("Li2 at {z}"));
    }
    // inversion for z outside the unit disc: Li2(z) = -pi^2/6 - ln^2(-z)/2 - Li2(1/z)
    for z in [c(-3.0, 0.0), c(-2.0, 1.5), c(1.5, 2.5)] {
        let want = -c(std::f64::consts::PI.powi(2) / 6.0, 0.0) - (-z).ln().powi(2) * 0.5 - li2_series(one / z);
        close(f(z), want, 1e-10, &format!("Li2 inversion at {z}"));
    }
    // the cut: Li2(2 + i0) = pi^2/4 + i pi ln 2
    close(f(c(2.0, 0.0)), c(std::f64::consts::PI.powi(2) / 4.0, std::f64::consts::PI * 2.0_f64.ln()), 1e-11, "Li2(2)");
    // unit circle: Li2(-1) = -pi^2/12, Li2(i) = -pi^2/48 + i G
    close(f(c(-1.0, 0.0)), c(-std::f64::consts::PI.powi(2) / 12.0, 0.0), 1e-11, "Li2(-1)");
    close(f(c(0.0, 1.0)), c(-std::f64::consts::PI.powi(2) / 48.0, 0.915_965_594_177_219), 1e-11, "Li2(i)");
    // Li3(z)/z = 4F3(1,1,1,1; 2,2,2; z): Li3(-1) = -3 zeta(3)/4
    let g = hyp_pfq(&[one; 4], &[two; 3], c(-1.0, 0.0)) * c(-1.0, 0.0);
    close(g, c(-0.75 * 1.202_056_903_159_594_3, 0.0), 1e-11, "Li3(-1)");
}

#[test]
fn convergence_acceleration_at_unit_argument() {
    let one = c(1.0, 0.0);
    let two = c(2.0, 0.0);
    // 3F2(1,1,1; 2,2; 1) = zeta(2); 4F3(1,1,1,1; 2,2,2; 1) = zeta(3)
    let got = hyp_pfq(&[one, one, one], &[two, two], one);
    close(got, c(std::f64::consts::PI.powi(2) / 6.0, 0.0), 1e-8, "zeta(2)");
    let got = hyp_pfq(&[one; 4], &[two; 3], one);
    close(got, c(1.202_056_903_159_594_3, 0.0), 1e-8, "zeta(3)");
    // divergent: Re(sum b - sum a) <= 0
    assert!(hyp_pfq(&[one, one, one], &[one, one], one).re.is_infinite());
    // p > q + 1 diverges
    assert!(hyp_pfq(&[one, one, one], &[two], c(0.3, 0.0)).is_nan());
    // terminating series: 2F1(-3, b; c; z) is a polynomial valid at any z
    let z = c(7.0, -2.0);
    let (b, cc) = (c(0.4, 0.5), c(1.3, -0.2));
    let want = {
        let mut t = c(1.0, 0.0);
        let mut s = t;
        for k in 0..3 {
            let kf = f64::from(k);
            t *= (c(-3.0, 0.0) + kf) * (b + kf) / ((cc + kf) * (kf + 1.0)) * z;
            s += t;
        }
        s
    };
    close(hyp2f1_cplx(c(-3.0, 0.0), b, cc, z), want, 1e-13, "terminating");
    assert!(hyp2f1_cplx(b, b, c(-2.0, 0.0), c(0.3, 0.0)).is_nan());
}

fn e1_series(z: C) -> C {
    // E1(z) = -gamma - ln z - sum_{k>=1} (-z)^k / (k k!)
    let mut s = c(0.0, 0.0);
    let mut t = c(1.0, 0.0);
    for k in 1..200 {
        t *= -z / f64::from(k);
        s += t / f64::from(k);
    }
    -c(0.577_215_664_901_532_9, 0.0) - z.ln() - s
}

#[test]
fn tricomi_u_closed_forms() {
    // U(a, a + 1, z) = z^(-a)
    let a = c(0.6, -0.9);
    for z in [c(0.05, 0.02), c(1.0, 1.0), c(8.0, -3.0), c(60.0, 10.0), c(-2.0, 1.0), c(-15.0, -4.0), c(0.0, 6.0)] {
        close(hyp_u_cplx(a, a + 1.0, z), z.powc(-a), 1e-10, &format!("z^-a at {z}"));
    }
    // U(1, 1, z) = e^z E1(z), including b = 1 (integer) with Re z < 0
    for z in [c(1.0, 0.0), c(2.0, 0.0), c(0.001, 0.0), c(5.0, 3.0), c(-2.0, 0.5), c(-0.5, -1.0), c(0.3, 7.0)] {
        let want = z.exp() * e1_series(z);
        close(hyp_u_cplx(c(1.0, 0.0), c(1.0, 0.0), z), want, 1e-9, &format!("e^z E1 at {z}"));
    }
    close(hyp_u_cplx(c(1.0, 0.0), c(1.0, 0.0), c(1.0, 0.0)), c(0.596_347_362_323_194_1, 0.0), 1e-12, "e E1(1)");
    // U(1/2, 1, 2x) = e^x K_0(x) / sqrt(pi); U(a, 2a, 2x) = e^x (2x)^(1/2 - a) K_(a-1/2)(x) / sqrt(pi)
    for x in [0.3, 2.0, 12.0, 35.0] {
        // the quadrature reference for K_nu is only good to ~1e-6 at x = 35
        let tol = if x > 30.0 { 1e-5 } else { 1e-9 };
        let got = hyp_u_cplx(c(0.5, 0.0), c(1.0, 0.0), c(2.0 * x, 0.0));
        let want = x.exp() * bessel_k(0.0, x) / std::f64::consts::PI.sqrt();
        assert!((got.re - want).abs() < tol * want, "K0 x={x}: {got} vs {want}");
        let aa = 1.3;
        let got = hyp_u_cplx(c(aa, 0.0), c(2.0 * aa, 0.0), c(2.0 * x, 0.0));
        let want = x.exp() * (2.0 * x).powf(0.5 - aa) * bessel_k(aa - 0.5, x) / std::f64::consts::PI.sqrt();
        assert!((got.re - want).abs() < tol * want, "K_(a-1/2) x={x}: {got} vs {want}");
    }
}

#[test]
fn tricomi_u_recurrences_and_reflection() {
    let a = c(0.8, 0.3);
    let b = c(1.4, -0.5);
    for z in [c(0.5, 0.2), c(3.0, 2.0), c(-1.0, 2.0), c(-6.0, -1.0), c(25.0, 25.0)] {
        // (b - a - 1) U(a, b-1) + (1 - b - z) U(a, b) + z U(a, b+1) = 0
        let r = (b - a - 1.0) * hyp_u_cplx(a, b - 1.0, z)
            + (c(1.0, 0.0) - b - z) * hyp_u_cplx(a, b, z)
            + z * hyp_u_cplx(a, b + 1.0, z);
        let scale = (z * hyp_u_cplx(a, b + 1.0, z)).norm().max((z * hyp_u_cplx(a, b, z)).norm());
        assert!(r.norm() < 1e-9 * scale, "b-recurrence at {z}: {r}");
        // U(a - 1, b) + (b - 2a - z) U(a, b) + a (a - b + 1) U(a + 1, b) = 0
        let r = hyp_u_cplx(a - 1.0, b, z)
            + (b - 2.0 * a - z) * hyp_u_cplx(a, b, z)
            + a * (a - b + 1.0) * hyp_u_cplx(a + 1.0, b, z);
        let scale = hyp_u_cplx(a - 1.0, b, z).norm().max(((b - 2.0 * a - z) * hyp_u_cplx(a, b, z)).norm());
        assert!(r.norm() < 1e-9 * scale, "a-recurrence at {z}: {r}");
        // U(a, b, z) = z^(1-b) U(a - b + 1, 2 - b, z)
        let lhs = hyp_u_cplx(a, b, z);
        let rhs = z.powc(c(1.0, 0.0) - b) * hyp_u_cplx(a - b + 1.0, c(2.0, 0.0) - b, z);
        close(lhs, rhs, 1e-9, &format!("Kummer transformation at {z}"));
    }
    // polynomial cases: U(-2, b, z) = z^2 - (b + 2... ) via the general formula (b)_2 M(-2, b, z)
    let z = c(1.7, -0.4);
    let got = hyp_u_cplx(c(-2.0, 0.0), b, z);
    let want = z * z - 2.0 * (b + 1.0) * z + b * (b + 1.0);
    close(got, want, 1e-12, "U(-2,b,z)");
}
