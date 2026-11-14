//! Coverage tests: numerical complex analysis, the session compute cache,
//! the term printers and a parser/evaluator differential test. Public API
//! only.

use std::f64::consts::PI;

use num_complex::Complex64;
use rssn::api::CachePolicy;
use rssn::api::Config;
use rssn::api::Session;
use rssn::kernels::complex::Mobius;
use rssn::kernels::complex::cauchy_derivative;
use rssn::kernels::complex::circle_integral;
use rssn::kernels::complex::complex_derivative;
use rssn::kernels::complex::contour_integral;
use rssn::kernels::complex::count_zeros_poles;
use rssn::kernels::complex::newton;
use rssn::kernels::complex::residue;

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }

    fn unit(&mut self) -> f64 {
        self.next() as f64 / (1u64 << 31) as f64
    }

    fn below(
        &mut self,
        n: u64,
    ) -> u64 {
        self.next() % n
    }
}

fn c(
    re: f64,
    im: f64,
) -> Complex64 {
    Complex64::new(re, im)
}

fn near(
    a: Complex64,
    b: Complex64,
    tol: f64,
) -> bool {
    (a - b).norm() <= tol * (1.0 + b.norm())
}

// --------------------------------------------------------------- contours

#[test]
fn circle_integrals_of_laurent_monomials() {
    let center = c(0.3, -0.2);
    for radius in [0.5, 1.0, 2.5] {
        for k in -4..=4i32 {
            let got = circle_integral(|z| (z - center).powi(k), center, radius, 64);
            let want = if k == -1 {
                c(0.0, 2.0 * PI)
            } else {
                c(0.0, 0.0)
            };
            assert!(near(got, want, 1e-11), "k = {k}, r = {radius}: {got}");
        }
    }
    // too few nodes are raised to the minimum instead of failing
    let few = circle_integral(|z| 1.0 / z, c(0.0, 0.0), 1.0, 0);
    assert!(near(few, c(0.0, 2.0 * PI), 1e-12));
    // geometric convergence of the trapezoid rule for exp(1/z) on |z| = 1
    let exact = c(0.0, 2.0 * PI);
    let err = |n| (circle_integral(|z| (1.0 / z).exp(), c(0.0, 0.0), 1.0, n) - exact).norm();
    assert!(err(32) < 1e-10 && err(16) > err(32));
}

#[test]
fn polygonal_contour_integrals() {
    // path independence for an entire function: int_0^{1+i} z^2 dz = (1+i)^3 / 3
    let want = c(1.0, 1.0).powi(3) / 3.0;
    let straight = contour_integral(|z| z * z, &[c(0.0, 0.0), c(1.0, 1.0)], 4);
    let bent = contour_integral(|z| z * z, &[c(0.0, 0.0), c(1.0, 0.0), c(1.0, 1.0)], 4);
    assert!(near(straight, want, 1e-13) && near(bent, want, 1e-13));
    // closed polygon, analytic integrand: zero
    let square = [
        c(0.0, 0.0),
        c(1.0, 0.0),
        c(1.0, 1.0),
        c(0.0, 1.0),
        c(0.0, 0.0),
    ];
    assert!(contour_integral(|z| z.exp(), &square, 8).norm() < 1e-6);
    // non-analytic integrand conj(z): the integral is 2 i * area
    assert!(near(
        contour_integral(|z| z.conj(), &square, 8),
        c(0.0, 2.0),
        1e-13
    ));
    // pole outside the square gives zero, inside gives 2 pi i
    assert!(contour_integral(|z| 1.0 / (z - c(3.0, 3.0)), &square, 32).norm() < 1e-9);
    assert!(near(
        contour_integral(|z| 1.0 / (z - c(0.5, 0.5)), &square, 64),
        c(0.0, 2.0 * PI),
        1e-7
    ));
    // degenerate paths
    assert_eq!(contour_integral(|z| z, &[], 4), c(0.0, 0.0));
    assert_eq!(contour_integral(|z| z, &[c(1.0, 1.0)], 4), c(0.0, 0.0));
    let zero_subdivision = contour_integral(|z| z, &[c(0.0, 0.0), c(2.0, 0.0)], 0);
    assert!(near(zero_subdivision, c(2.0, 0.0), 1e-13));
}

#[test]
fn residues_and_the_argument_principle() {
    let a = c(1.0, 0.5);
    let b = c(-0.5, 0.25);
    let f = |z: Complex64| 1.0 / ((z - a) * (z - b));
    assert!(near(residue(f, a, 0.1, 64), 1.0 / (a - b), 1e-12));
    assert!(near(residue(f, b, 0.1, 64), 1.0 / (b - a), 1e-12));
    // exp(z)/z^3 has residue 1/2 at the origin; z^2 has none
    assert!(near(
        residue(|z| z.exp() / z.powi(3), c(0.0, 0.0), 1.0, 64),
        c(0.5, 0.0),
        1e-12
    ));
    assert!(residue(|z| z * z, c(0.0, 0.0), 1.0, 64).norm() < 1e-14);
    // zeros minus poles
    let cube = |z: Complex64| z.powi(3) - 1.0;
    assert_eq!(count_zeros_poles(cube, c(0.0, 0.0), 2.0, 256), Some(3));
    assert_eq!(count_zeros_poles(cube, c(0.0, 0.0), 0.5, 256), Some(0));
    assert_eq!(
        count_zeros_poles(|z| 1.0 / z, c(0.0, 0.0), 1.0, 128),
        Some(-1)
    );
    assert_eq!(
        count_zeros_poles(|z| z * z / (z - c(0.5, 0.0)), c(0.0, 0.0), 1.0, 512),
        Some(1)
    );
    assert_eq!(
        count_zeros_poles(|z| z.exp(), c(0.0, 0.0), 3.0, 256),
        Some(0)
    );
    // a zero on the contour, a non-finite value, a hopelessly undersampled contour
    assert_eq!(count_zeros_poles(|z| z - 1.0, c(0.0, 0.0), 1.0, 128), None);
    assert_eq!(
        count_zeros_poles(|_| c(f64::NAN, 0.0), c(0.0, 0.0), 1.0, 128),
        None
    );
    assert_eq!(
        count_zeros_poles(|z| z.powi(40), c(0.0, 0.0), 1.0, 64),
        None
    );
}

#[test]
fn complex_derivatives() {
    let z0 = c(0.4, 0.7);
    assert!(near(
        complex_derivative(&|z: Complex64| z.sin(), z0),
        z0.cos(),
        1e-8
    ));
    assert!(near(
        complex_derivative(&|z: Complex64| z.exp(), z0),
        z0.exp(),
        1e-8
    ));
    // Cauchy's formula for every order of exp: f^(n) = exp
    for order in 0..8 {
        let got = cauchy_derivative(|z| z.exp(), z0, order, 0.8, 64);
        assert!(near(got, z0.exp(), 1e-11), "order {order}");
    }
    // geometric series: d^n/dz^n 1/(1-z) at 0 is n!
    let mut fact = 1.0;
    for order in 0..8u32 {
        if order > 0 {
            fact *= f64::from(order);
        }
        let got = cauchy_derivative(|z| 1.0 / (1.0 - z), c(0.0, 0.0), order, 0.5, 128);
        assert!(near(got, c(fact, 0.0), 1e-10), "order {order}: {got}");
    }
}

#[test]
fn newton_iteration_in_the_complex_plane() {
    let f = |z: Complex64| z.powi(3) - 1.0;
    let df = |z: Complex64| 3.0 * z * z;
    for start in [c(0.9, 0.1), c(-0.6, 0.8), c(-0.6, -0.8), c(2.0, 0.0)] {
        let root = newton(f, df, start, 1e-14, 100).expect("converges");
        assert!(f(root).norm() < 1e-12, "{start} -> {root}");
        assert!((root.norm() - 1.0).abs() < 1e-12);
    }
    let i = newton(|z| z * z + 1.0, |z| 2.0 * z, c(1.0, 1.0), 1e-14, 50).unwrap();
    assert!(near(i, c(0.0, 1.0), 1e-12));
    // vanishing derivative, no iterations allowed, non-finite derivative
    assert!(newton(f, df, c(0.0, 0.0), 1e-12, 50).is_none());
    assert!(newton(f, df, c(2.0, 0.0), 1e-12, 0).is_none());
    assert!(newton(f, |_| c(f64::NAN, 0.0), c(2.0, 0.0), 1e-12, 5).is_none());
    // no real root: Newton from the real axis on x^2 + 1 stays real and fails to settle quickly
    let _ = newton(|z| z * z + 1.0, |z| 2.0 * z, c(0.5, 0.0), 1e-12, 30);
}

#[test]
fn mobius_group_laws() {
    let m = Mobius::new(c(1.0, 1.0), c(0.0, 2.0), c(1.0, 0.0), c(3.0, -1.0));
    let n = Mobius::new(c(2.0, 0.0), c(-1.0, 1.0), c(0.5, 0.0), c(1.0, 1.0));
    let id = Mobius::identity();
    let mut rng = Lcg(5);
    for _ in 0..50 {
        let z = c(rng.unit() * 4.0 - 2.0, rng.unit() * 4.0 - 2.0);
        assert!(near(id.apply(z), z, 1e-15));
        // composition acts as function composition
        assert!(near(m.compose(&n).apply(z), m.apply(n.apply(z)), 1e-10));
        // inverse (projective): m^-1(m(z)) = z
        assert!(near(m.inverse().apply(m.apply(z)), z, 1e-9));
        assert!(near(m.apply(m.inverse().apply(z)), z, 1e-9));
        // composition with the inverse is the identity map
        assert!(near(m.compose(&m.inverse()).apply(z), z, 1e-9));
    }
    // associativity of the matrix product
    let k = Mobius::new(c(0.0, 1.0), c(1.0, 0.0), c(1.0, 0.0), c(0.0, 0.0));
    let (l, r) = (m.compose(&n).compose(&k), m.compose(&n.compose(&k)));
    for (x, y) in [(l.a, r.a), (l.b, r.b), (l.c, r.c), (l.d, r.d)] {
        assert!(near(x, y, 1e-12));
    }
    // cross-ratio is invariant
    let cr = |a: Complex64, b: Complex64, c_: Complex64, d: Complex64| {
        (a - c_) * (b - d) / ((a - d) * (b - c_))
    };
    let pts = [c(0.1, 0.2), c(1.3, -0.4), c(-0.8, 0.9), c(0.5, 1.5)];
    let image: Vec<Complex64> = pts.iter().map(|&p| m.apply(p)).collect();
    assert!(near(
        cr(pts[0], pts[1], pts[2], pts[3]),
        cr(image[0], image[1], image[2], image[3]),
        1e-10
    ));
    // Cayley transform maps the real line to the unit circle
    let cayley = Mobius::new(c(1.0, 0.0), c(0.0, -1.0), c(1.0, 0.0), c(0.0, 1.0));
    for t in [-50.0, -2.0, -0.3, 0.0, 0.7, 4.0, 1e6] {
        assert!(
            (cayley.apply(c(t, 0.0)).norm() - 1.0).abs() < 1e-12,
            "t = {t}"
        );
    }
    // at the pole the image is not finite, but nothing panics
    assert!(!cayley.apply(c(0.0, -1.0)).is_finite());
}

// ------------------------------------------------------------ session cache

const REQUEST: &str = "diff(x^3 + sin(x), x)";

fn compute<'s>(
    s: &'s Session,
    text: &str,
    config: &Config,
) -> rssn::api::Answer<'s> {
    let t = s.parse(text).unwrap();
    s.compute(t, config).unwrap()
}

#[test]
fn cache_defaults_and_hits() {
    let s = Session::new();
    assert_eq!(s.cache_policy(), CachePolicy::default());
    assert_eq!(
        s.cache_policy(),
        CachePolicy {
            answers: true,
            lemmas: false,
            capacity: 4096
        }
    );
    let stats = s.cache_stats();
    assert_eq!(
        (
            stats.answer_hits,
            stats.answer_misses,
            stats.answers,
            stats.lemmas
        ),
        (0, 0, 0, 0)
    );
    let cfg = Config::new();
    let first = compute(&s, REQUEST, &cfg);
    let st = s.cache_stats();
    assert_eq!((st.answer_hits, st.answer_misses, st.answers), (0, 1, 1));
    let second = compute(&s, REQUEST, &cfg);
    let st = s.cache_stats();
    assert_eq!((st.answer_hits, st.answer_misses, st.answers), (1, 1, 1));
    // a hit is indistinguishable from the fresh answer
    assert_eq!(first.term.to_latex(), second.term.to_latex());
    assert_eq!(first.value, second.value);
    assert_eq!(first.reduced, second.reduced);
    // hash-consing: the same text is the same node, so it also hits
    let _ = compute(&s, REQUEST, &cfg);
    assert_eq!(s.cache_stats().answer_hits, 2);
    // another configuration is another key
    let _ = compute(&s, REQUEST, &Config::new().bind("y", 1.0));
    let st = s.cache_stats();
    assert_eq!((st.answer_misses, st.answers), (2, 2));
    // another request is another key
    let _ = compute(&s, "diff(x^4, x)", &cfg);
    assert_eq!(s.cache_stats().answers, 3);
}

#[test]
fn cache_policy_off_clear_and_capacity() {
    let cfg = Config::new();
    // off: nothing is stored or reused, answers are the same
    let off = Session::new();
    off.set_cache_policy(CachePolicy::off());
    let a = compute(&off, REQUEST, &cfg);
    let b = compute(&off, REQUEST, &cfg);
    let st = off.cache_stats();
    assert_eq!((st.answer_hits, st.answers), (0, 0));
    assert_eq!(a.term.to_latex(), b.term.to_latex());
    let warm = Session::new();
    let w1 = compute(&warm, REQUEST, &cfg);
    assert_eq!(
        w1.term.to_latex(),
        a.term.to_latex(),
        "cache on/off give the same answer"
    );

    // clear drops everything and a later request misses again
    let s = Session::new();
    let _ = compute(&s, REQUEST, &cfg);
    assert_eq!(s.cache_stats().answers, 1);
    s.clear_cache();
    assert_eq!(s.cache_stats().answers, 0);
    assert!(s.cache_stats().invalidations >= 1);
    let _ = compute(&s, REQUEST, &cfg);
    assert_eq!(s.cache_stats().answer_misses, 2);

    // disabling the answer table drops entries; growing the capacity keeps them
    let _ = compute(&s, REQUEST, &cfg);
    s.set_cache_policy(CachePolicy {
        answers: true,
        lemmas: false,
        capacity: 8192,
    });
    assert_eq!(s.cache_stats().answers, 1);
    s.set_cache_policy(CachePolicy::off());
    assert_eq!(s.cache_stats().answers, 0);

    // capacity 2 evicts first-in first-out
    let small = Session::new();
    small.set_cache_policy(CachePolicy {
        answers: true,
        lemmas: false,
        capacity: 2,
    });
    for k in 2..7 {
        let _ = compute(&small, &format!("diff(x^{k}, x)"), &cfg);
    }
    let st = small.cache_stats();
    assert_eq!(st.answers, 2);
    assert_eq!(st.evictions, 3);
    // the oldest was evicted (miss), the newest is still there (hit)
    let hits = st.answer_hits;
    let _ = compute(&small, "diff(x^6, x)", &cfg);
    assert_eq!(small.cache_stats().answer_hits, hits + 1);
    let misses = small.cache_stats().answer_misses;
    let _ = compute(&small, "diff(x^2, x)", &cfg);
    assert_eq!(small.cache_stats().answer_misses, misses + 1);
}

#[test]
fn cache_numeric_answers_and_lemmas_do_not_change_results() {
    let cfg = Config::new().numeric(1e-10);
    let s = Session::new();
    let a = compute(&s, "defint(x^2, x, 0, 1)", &cfg);
    let v = a.value.expect("numeric value");
    assert!((v - 1.0 / 3.0).abs() < 1e-8, "{v}");
    let b = compute(&s, "defint(x^2, x, 0, 1)", &cfg);
    assert_eq!(b.value.map(f64::to_bits), a.value.map(f64::to_bits));
    assert_eq!(s.cache_stats().answer_hits, 1);
    // lemma reuse may pick a different (equal) closed form: compare numerically
    let cold = Session::new();
    cold.set_cache_policy(CachePolicy::off());
    let lem = Session::new();
    lem.set_cache_policy(CachePolicy::full());
    let sym = Config::new();
    let _ = compute(&lem, "integral(x*exp(x), x)", &sym);
    let warm = compute(&lem, "diff(integral(x*exp(x), x), x)", &sym);
    let fresh = compute(&cold, "diff(integral(x*exp(x), x), x)", &sym);
    for x in [-1.3, 0.0, 0.7, 2.1] {
        let (w, f) = (warm.term.eval(&[("x", x)]), fresh.term.eval(&[("x", x)]));
        match (w, f) {
            | (Some(w), Some(f)) => {
                assert!(
                    (w - f).abs() < 1e-9 * (1.0 + f.abs()),
                    "x = {x}: {w} vs {f}"
                )
            },
            | (None, None) => {},
            | _ => panic!("one side evaluates, the other does not at x = {x}"),
        }
    }
    assert_eq!(warm.reduced, fresh.reduced);
}

// ---------------------------------------------------- printers and parser

#[test]
fn printers_cover_the_basic_constructs() {
    let s = Session::new();
    let p = |text: &str| s.parse(text).unwrap();
    let t = p("alpha_1 + beta*x^2 - 1/2");
    let (latex, typst, pretty) = (t.to_latex(), t.to_typst(), t.to_pretty());
    assert!(
        latex.contains(r"\alpha") && latex.contains(r"\beta"),
        "{latex}"
    );
    assert!(latex.contains(r"\frac"), "{latex}");
    assert!(latex.contains('^'), "{latex}");
    assert!(typst.contains("alpha") && typst.contains("beta"), "{typst}");
    assert!(typst.contains("frac("), "{typst}");
    assert!(pretty.contains('α') && pretty.contains('β'), "{pretty}");
    // function applications and roots
    let f = p("sqrt(x) + sin(x)/cos(x)");
    assert!(f.to_latex().contains(r"\sqrt"), "{}", f.to_latex());
    assert!(
        f.to_latex().contains(r"\sin") && f.to_latex().contains(r"\cos"),
        "{}",
        f.to_latex()
    );
    assert!(f.to_typst().contains("sqrt"), "{}", f.to_typst());
    // derivative notation and binders
    let d = p("diff(f(x), x)");
    assert!(
        d.to_latex().contains(r"\frac{d}{d")
            || d.to_latex().contains(r"\frac{\mathrm{d}")
            || d.to_latex().contains("frac"),
        "{}",
        d.to_latex()
    );
    let i = p("integral(x^2, x)");
    assert!(i.to_latex().contains(r"\int"), "{}", i.to_latex());
    // deterministic: rendering twice gives the same text
    assert_eq!(t.to_latex(), t.to_latex());
    // parenthesisation: (a+b)*c keeps the group, a+b*c does not need one
    let grouped = p("(a + b)*c").to_latex();
    assert!(
        grouped.contains('(') || grouped.contains("left"),
        "{grouped}"
    );
    let flat = p("a + b*c").to_latex();
    assert!(!flat.contains('(') && !flat.contains("left"), "{flat}");
}

/// Random fully-parenthesised expression over `x`, `y` and small integers.
fn gen_expr(
    rng: &mut Lcg,
    depth: u32,
    x: f64,
    y: f64,
) -> (String, f64) {
    if depth == 0 || rng.below(5) == 0 {
        return match rng.below(4) {
            | 0 => ("x".into(), x),
            | 1 => ("y".into(), y),
            | _ => {
                let k = 1 + rng.below(9);
                (k.to_string(), k as f64)
            },
        };
    }
    let (a, va) = gen_expr(rng, depth - 1, x, y);
    match rng.below(7) {
        | 0 => (format!("({a} + {})", "1"), va + 1.0),
        | 1 => {
            let (b, vb) = gen_expr(rng, depth - 1, x, y);
            (format!("({a} + {b})"), va + vb)
        },
        | 2 => {
            let (b, vb) = gen_expr(rng, depth - 1, x, y);
            (format!("({a} - {b})"), va - vb)
        },
        | 3 => {
            let (b, vb) = gen_expr(rng, depth - 1, x, y);
            (format!("({a} * {b})"), va * vb)
        },
        | 4 => {
            let (b, vb) = gen_expr(rng, depth - 1, x, y);
            (format!("({a} / {b})"), va / vb)
        },
        | 5 => (format!("({a})^2"), va * va),
        | _ => (format!("(-{a})"), -va),
    }
}

#[test]
fn parser_and_evaluator_agree_with_a_reference_evaluation() {
    let s = Session::new();
    let mut rng = Lcg(2024);
    let mut checked = 0;
    for case in 0..300 {
        let (x, y) = (rng.unit() * 6.0 - 3.0, rng.unit() * 6.0 - 3.0);
        let (text, want) = gen_expr(&mut rng, 4, x, y);
        let Ok(term) = s.parse(&text) else {
            panic!("case {case}: cannot parse `{text}`")
        };
        // printers never panic and keep braces and parentheses balanced
        let latex = term.to_latex();
        let balance = |open: char, close: char, s: &str| {
            s.chars().fold(0i64, |d, ch| {
                d + i64::from(ch == open) - i64::from(ch == close)
            })
        };
        assert_eq!(balance('{', '}', &latex), 0, "`{text}` -> {latex}");
        assert_eq!(
            latex.matches(r"\left(").count(),
            latex.matches(r"\right)").count(),
            "`{text}` -> {latex}"
        );
        assert!(!term.to_typst().is_empty() && !term.to_pretty().is_empty());
        if !want.is_finite() || want.abs() > 1e9 {
            continue;
        }
        let Some(got) = term.eval(&[("x", x), ("y", y)]) else {
            continue;
        };
        if got.is_finite() {
            assert!(
                (got - want).abs() <= 1e-9 * (1.0 + want.abs()),
                "case {case}: `{text}` at x = {x}, y = {y}: {got} vs {want}"
            );
            checked += 1;
        }
    }
    assert!(checked > 150, "only {checked} cases were comparable");
}

#[test]
fn operator_precedence_and_associativity() {
    let s = Session::new();
    let ev = |text: &str, x: f64| s.parse(text).unwrap().eval(&[("x", x)]).unwrap();
    assert_eq!(ev("2 + 3*4^2", 0.0), 50.0);
    assert_eq!(ev("(1 + 2)*3", 0.0), 9.0);
    assert_eq!(ev("10 - 4 - 3", 0.0), 3.0);
    assert_eq!(ev("8/4/2", 0.0), 1.0);
    assert_eq!(ev("2^3^2", 0.0), 512.0);
    assert_eq!(ev("-x^2", 3.0), -9.0);
    assert_eq!(ev("(-x)^2", 3.0), 9.0);
    assert_eq!(ev("2*x + 1", 4.0), 9.0);
    // malformed input is an error, not a panic
    for bad in ["", "1 +", "(1 + 2", "1 + 2)", "* 3", "foo(", "1 2 3 +"] {
        let _ = s.parse(bad);
    }
    assert!(s.parse("(1 + 2").is_err());
    assert!(s.parse("1 +").is_err());
}
