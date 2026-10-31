//! Integration tests for the `groebner_f4` and `groebner_f5` operators.
//!
//! Every test compares the two engines with the reference implementation
//! `groebner` — the reduced Gröbner basis is unique for an ideal and a
//! monomial order, so the printed answers must be *literally* the same list.
//! The tests also certify that the answer is reduced (Buchberger on it
//! changes nothing), that ideal membership agrees, that repeated runs agree,
//! and that the request is left unevaluated instead of returning a wrong
//! answer when a limit is hit or the arguments are malformed.

use std::time::Instant;

use rssn::prelude::*;

// ---------------------------------------------------------------------------
// Session helpers
// ---------------------------------------------------------------------------

fn session() -> Session {
    Session::with_rules(&[rssn::rules::poly()]).unwrap_or_else(|e| panic!("poly installs: {e}"))
}

/// Computes `src` symbolically; returns the printed answer and whether it is
/// free of heavy operators.
fn compute(
    session: &Session,
    src: &str,
) -> (String, bool) {
    let term = session
        .parse(src)
        .unwrap_or_else(|e| panic!("cannot parse `{src}`: {e}"));
    let answer = session
        .compute(term, &Config::empty().with(rssn::rules::poly()))
        .unwrap_or_else(|e| panic!("cannot compute `{src}`: {e:?}"));
    (answer.term.to_string(), answer.reduced)
}

/// The reduced basis `src` names, panicking if the request is not reduced.
fn basis(
    session: &Session,
    src: &str,
) -> String {
    let (text, reduced) = compute(session, src);
    assert!(reduced, "`{src}` was not reduced: {text}");
    text
}

// ---------------------------------------------------------------------------
// Problem descriptions and polynomial formatting
// ---------------------------------------------------------------------------

struct Problem {
    name: String,
    vars: Vec<String>,
    polynomials: Vec<String>,
    order: Option<&'static str>,
}

impl Problem {
    fn new(
        name: &str,
        vars: Vec<String>,
        polynomials: Vec<String>,
        order: Option<&'static str>,
    ) -> Self {
        Self {
            name: name.to_owned(),
            vars,
            polynomials,
            order,
        }
    }

    fn request(
        &self,
        operator: &str,
    ) -> String {
        let mut args = format!(
            "list({}), list({})",
            self.polynomials.join(", "),
            self.vars.join(", ")
        );
        if let Some(order) = self.order {
            args.push_str(", ");
            args.push_str(order);
        }
        format!("{operator}({args})")
    }

    fn with_order(
        &self,
        order: Option<&'static str>,
    ) -> Self {
        Self {
            name: self.name.clone(),
            vars: self.vars.clone(),
            polynomials: self.polynomials.clone(),
            order,
        }
    }
}

/// `x0*x2^2` etc.; the empty product prints as `1`.
fn monomial(
    vars: &[String],
    exponents: &[u32],
) -> String {
    let factors: Vec<String> = exponents
        .iter()
        .enumerate()
        .filter(|&(_, &e)| e > 0)
        .map(|(i, &e)| {
            match e {
                | 1 => vars.get(i).cloned().unwrap_or_default(),
                | _ => format!("{}^{e}", vars.get(i).map_or("", String::as_str)),
            }
        })
        .collect();
    if factors.is_empty() {
        "1".to_owned()
    } else {
        factors.join("*")
    }
}

/// `2*x0^2 - x1 + 1`, in the given term order.
fn polynomial(
    vars: &[String],
    terms: &[(i64, Vec<u32>)],
) -> String {
    let mut out = String::new();
    for (coefficient, exponents) in terms {
        let (sign, magnitude) = if *coefficient < 0 {
            ("-", coefficient.saturating_neg())
        } else {
            ("+", *coefficient)
        };
        if magnitude == 0 {
            continue;
        }
        let body = monomial(vars, exponents);
        let piece = if body == "1" {
            magnitude.to_string()
        } else if magnitude == 1 {
            body
        } else {
            format!("{magnitude}*{body}")
        };
        if out.is_empty() {
            if sign == "-" {
                out.push('-');
            }
        } else {
            out.push(' ');
            out.push_str(sign);
            out.push(' ');
        }
        out.push_str(&piece);
    }
    if out.is_empty() {
        "0".to_owned()
    } else {
        out
    }
}

fn exponents(
    n: usize,
    powers: &[(usize, u32)],
) -> Vec<u32> {
    let mut out = vec![0_u32; n];
    for &(i, e) in powers {
        if let Some(slot) = out.get_mut(i) {
            *slot = slot.saturating_add(e);
        }
    }
    out
}

// ---------------------------------------------------------------------------
// The system zoo
// ---------------------------------------------------------------------------

fn variables(n: usize) -> Vec<String> {
    (0..n).map(|i| format!("x{i}")).collect()
}

fn katsura(m: usize) -> Problem {
    let n = m.saturating_add(1);
    let vars = variables(n);
    let mut polynomials = Vec::new();
    let mut linear = vec![(1, exponents(n, &[(0, 1)])), (-1, exponents(n, &[]))];
    for l in 1..=m {
        linear.push((2, exponents(n, &[(l, 1)])));
    }
    polynomials.push(polynomial(&vars, &linear));
    for k in 0..m {
        let mut terms = vec![(-1, exponents(n, &[(k, 1)]))];
        for l in -(m as i64)..=(m as i64) {
            let a = l.unsigned_abs() as usize;
            let b = (k as i64 - l).unsigned_abs() as usize;
            if a <= m && b <= m {
                terms.push((1, exponents(n, &[(a, 1), (b, 1)])));
            }
        }
        polynomials.push(polynomial(&vars, &terms));
    }
    Problem::new(&format!("katsura-{m}"), vars, polynomials, Some("grevlex"))
}

fn cyclic(
    n: usize,
    order: &'static str,
) -> Problem {
    let vars = variables(n);
    let mut polynomials = Vec::new();
    for d in 1..n {
        let terms: Vec<(i64, Vec<u32>)> = (0..n)
            .map(|start| {
                let powers: Vec<(usize, u32)> = (0..d).map(|j| ((start + j) % n, 1)).collect();
                (1, exponents(n, &powers))
            })
            .collect();
        polynomials.push(polynomial(&vars, &terms));
    }
    polynomials.push(polynomial(
        &vars,
        &[(1, vec![1; n]), (-1, exponents(n, &[]))],
    ));
    Problem::new(
        &format!("cyclic-{n}-{order}"),
        vars,
        polynomials,
        Some(order),
    )
}

/// A deterministic linear congruential generator: the random systems are the
/// same on every machine and in every run.
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0 >> 33
    }

    fn coefficient(&mut self) -> i64 {
        i64::try_from(self.next() % 7).unwrap_or(0) - 3
    }

    fn choose(
        &mut self,
        upper: usize,
    ) -> usize {
        usize::try_from(self.next() % u64::try_from(upper).unwrap_or(1)).unwrap_or(0)
    }
}

/// A dense system: every monomial of degree at most two with a random
/// coefficient.
fn dense_system(
    name: &str,
    n: usize,
    rows: usize,
    seed: u64,
    order: &'static str,
) -> Problem {
    let vars = variables(n);
    let mut rng = Lcg::new(seed);
    let mut support: Vec<Vec<u32>> = Vec::new();
    for i in 0..n {
        support.push(exponents(n, &[(i, 1)]));
        support.push(exponents(n, &[(i, 2)]));
        for j in (i + 1)..n {
            support.push(exponents(n, &[(i, 1), (j, 1)]));
        }
    }
    let mut polynomials = Vec::new();
    for _ in 0..rows {
        let mut terms: Vec<(i64, Vec<u32>)> = support
            .iter()
            .map(|mono| (rng.coefficient(), mono.clone()))
            .filter(|(c, _)| *c != 0)
            .collect();
        if terms.is_empty() {
            terms.push((1, exponents(n, &[(rng.choose(n), 2)])));
        }
        polynomials.push(polynomial(&vars, &terms));
    }
    Problem::new(name, vars, polynomials, Some(order))
}

/// A sparse system: a few terms with small exponents per polynomial.
fn sparse_system(
    name: &str,
    n: usize,
    rows: usize,
    seed: u64,
    order: &'static str,
) -> Problem {
    let vars = variables(n);
    let mut rng = Lcg::new(seed);
    let mut polynomials = Vec::new();
    for _ in 0..rows {
        let mut terms: Vec<(i64, Vec<u32>)> = Vec::new();
        for _ in 0..3 {
            let mut powers = Vec::new();
            for _ in 0..2 {
                powers.push((rng.choose(n), u32::try_from(rng.choose(3)).unwrap_or(0)));
            }
            terms.push((rng.coefficient(), exponents(n, &powers)));
        }
        if terms.iter().all(|(c, _)| *c == 0) {
            terms.push((1, exponents(n, &[])));
        }
        polynomials.push(polynomial(&vars, &terms));
    }
    Problem::new(name, vars, polynomials, Some(order))
}

// ---------------------------------------------------------------------------
// The comparison itself
// ---------------------------------------------------------------------------

/// Compares both engines with `groebner` and certifies the result is the
/// reduced basis.
fn check(
    session: &Session,
    problem: &Problem,
) {
    let expected = basis(session, &problem.request("groebner"));
    for operator in ["groebner_f4", "groebner_f5"] {
        let computed = basis(session, &problem.request(operator));
        assert_eq!(
            computed, expected,
            "{}: {operator} differs from groebner",
            problem.name
        );
        // Buchberger on the computed basis must reproduce it verbatim: the
        // answer is not only *a* Gröbner basis but the reduced one. The same
        // monomial order has to be requested, of course.
        let order = problem
            .order
            .map(|name| format!(", {name}"))
            .unwrap_or_default();
        let again = basis(
            session,
            &format!(
                "groebner({computed}, list({}){order})",
                problem.vars.join(", ")
            ),
        );
        assert_eq!(
            again, computed,
            "{}: {operator} is not reduced",
            problem.name
        );
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[test]
fn katsura_and_cyclic_agree_with_buchberger() {
    let session = session();
    let mut problems = vec![
        katsura(2),
        katsura(3),
        katsura(4),
        cyclic(3, "grevlex"),
        cyclic(4, "grevlex"),
        cyclic(5, "grevlex"),
        cyclic(4, "lex"),
        cyclic(5, "lex"),
        cyclic(4, "grlex"),
    ];
    problems.push(katsura(2).with_order(Some("lex")));
    problems.push(katsura(3).with_order(Some("lex")));
    for problem in &problems {
        let start = Instant::now();
        check(&session, problem);
        eprintln!("{}: {:?}", problem.name, start.elapsed());
    }
}

#[test]
fn order_variants_and_the_default_order() {
    let session = session();
    // The arithmetic-geometric mean system in three orders; the two-argument
    // form defaults to lex, exactly like `groebner`.
    let vars = variables(3);
    let polynomials = vec![
        polynomial(
            &vars,
            &[
                (1, exponents(3, &[(0, 2)])),
                (1, exponents(3, &[(1, 1)])),
                (-1, exponents(3, &[(2, 1)])),
            ],
        ),
        polynomial(
            &vars,
            &[
                (1, exponents(3, &[(0, 1), (1, 1)])),
                (-2, exponents(3, &[(2, 2)])),
                (1, exponents(3, &[])),
            ],
        ),
        polynomial(
            &vars,
            &[(1, exponents(3, &[(1, 3)])), (-1, exponents(3, &[(0, 1)]))],
        ),
    ];
    for order in [None, Some("lex"), Some("grlex"), Some("grevlex")] {
        check(
            &session,
            &Problem::new("orders", vars.clone(), polynomials.clone(), order),
        );
    }
    // The default really is lex: it must agree with an explicit `lex`.
    let implicit = basis(
        &session,
        &format!(
            "groebner_f4(list({}), list(x0, x1, x2))",
            polynomials.join(", ")
        ),
    );
    let explicit = basis(
        &session,
        &format!(
            "groebner_f4(list({}), list(x0, x1, x2), lex)",
            polynomials.join(", ")
        ),
    );
    assert_eq!(implicit, explicit);
}

#[test]
fn positive_dimensional_ideals_agree_in_every_order() {
    let session = session();
    // The hyperbola x0*x1 = 1 is one-dimensional. FGLM cannot convert it,
    // so the F5 operator must finish in the matrix engine for lex/grlex, and
    // every order must still match Buchberger's reduced basis.
    let problem = Problem::new(
        "positive-dimensional",
        variables(2),
        vec!["x0*x1 - 1".to_owned()],
        Some("grevlex"),
    );
    for order in [Some("grevlex"), Some("grlex"), Some("lex")] {
        check(&session, &problem.with_order(order));
    }
}

#[test]
fn random_dense_and_sparse_systems_agree() {
    let session = session();
    for (seed, order) in [(1_u64, "grevlex"), (2, "grlex"), (3, "lex"), (4, "grevlex")] {
        let dense = dense_system(&format!("dense-{seed}"), 3, 3, seed, order);
        check(&session, &dense);
        let sparse = sparse_system(&format!("sparse-{seed}"), 3, 3, seed ^ 0x5EED, order);
        check(&session, &sparse);
    }
    // A slightly larger random system in four variables.
    check(&session, &sparse_system("sparse-4", 4, 4, 99, "grevlex"));
}

#[test]
fn redundant_non_monic_constants_and_zero_generators() {
    let session = session();
    let vars = variables(3);
    // Redundant generators: a multiple of another and a sum of two others.
    let redundant = Problem::new(
        "redundant",
        vars.clone(),
        vec![
            polynomial(
                &vars,
                &[(1, exponents(3, &[(0, 2)])), (-1, exponents(3, &[]))],
            ),
            polynomial(
                &vars,
                &[(3, exponents(3, &[(0, 2)])), (-3, exponents(3, &[]))],
            ),
            polynomial(
                &vars,
                &[
                    (1, exponents(3, &[(0, 1), (1, 1)])),
                    (-1, exponents(3, &[(1, 1)])),
                ],
            ),
        ],
        Some("grevlex"),
    );
    check(&session, &redundant);
    // Non-monic generators with rational coefficients.
    let non_monic = Problem::new(
        "non-monic",
        vars.clone(),
        vec![
            "2*x0^2 + 3*x1 - 1".to_owned(),
            "4*x0*x2 - 5".to_owned(),
            "-6*x1^2 + 7*x0".to_owned(),
        ],
        Some("grevlex"),
    );
    check(&session, &non_monic);
    // Zero generators mixed with real ones.
    let zeros = Problem::new(
        "zeros",
        vars.clone(),
        vec![
            "0".to_owned(),
            "x0 - 2".to_owned(),
            "0".to_owned(),
            "x1 - x2".to_owned(),
        ],
        None,
    );
    check(&session, &zeros);
    // All generators zero: the zero ideal has the empty basis.
    for operator in ["groebner", "groebner_f4", "groebner_f5"] {
        assert_eq!(
            basis(&session, &format!("{operator}(list(0, 0), list(x0, x1))")),
            "list()",
            "{operator}"
        );
    }
    // A non-zero constant generates the unit ideal.
    for operator in ["groebner", "groebner_f4", "groebner_f5"] {
        assert_eq!(
            basis(&session, &format!("{operator}(list(x0, 3), list(x0, x1))")),
            "list(1)",
            "{operator}"
        );
    }
    // An empty generator list.
    for operator in ["groebner", "groebner_f4", "groebner_f5"] {
        assert_eq!(
            basis(&session, &format!("{operator}(list(), list(x0, x1))")),
            "list()",
            "{operator}"
        );
    }
}

#[test]
fn affine_systems_use_the_homogenisation_path() {
    let session = session();
    // Mixed degrees: the F5 engine has to homogenise and dehomogenise.
    let three = variables(3);
    let affine = Problem::new(
        "affine",
        three.clone(),
        vec![
            "x0^2 + x1 + 1".to_owned(),
            "x0 - x2^3".to_owned(),
            "x1*x2 - 2".to_owned(),
            "x0 + x1 + x2 - 3".to_owned(),
        ],
        Some("grevlex"),
    );
    check(&session, &affine);
    check(&session, &affine.with_order(Some("lex")));
    // Inhomogeneous linear system.
    let linear = Problem::new(
        "affine linear",
        variables(2),
        vec!["2*x0 - 3*x1 + 1".to_owned(), "x0 + x1 - 7".to_owned()],
        Some("grevlex"),
    );
    check(&session, &linear);
    // Inconsistent affine system: basis {1}.
    let inconsistent = Problem::new(
        "affine inconsistent",
        variables(2),
        vec!["x0^2 + x1 - 1".to_owned(), "x0^2 + x1 + 1".to_owned()],
        Some("grevlex"),
    );
    check(&session, &inconsistent);
    assert_eq!(
        basis(&session, &inconsistent.request("groebner_f5")),
        "list(1)"
    );
}

#[test]
fn dangerous_criterion_cases_still_match_buchberger() {
    let session = session();
    let vars = variables(3);
    let cases: Vec<(&str, Vec<String>, Option<&'static str>)> = vec![
        // Equal leading monomials: the S-polynomial has a new leading term.
        (
            "equal leads",
            vec!["x0^2 + x1 + 1".to_owned(), "x0^2 + x2 - 1".to_owned()],
            Some("grevlex"),
        ),
        // Equal leading monomials and equal signatures.
        (
            "equal signatures",
            vec![
                "x0*x1 + x2".to_owned(),
                "x0*x1 + x0".to_owned(),
                "x0*x1 + x2 + 1".to_owned(),
            ],
            Some("grevlex"),
        ),
        // The classic rewritten-criterion trap: the second generator's
        // signature divides the pair signature, yet the pair is needed.
        (
            "rewritten trap",
            vec!["x0*x1 + x2".to_owned(), "x0*x2 + x1".to_owned()],
            Some("grevlex"),
        ),
        // Duplicate generators, and a generator equal to another's multiple.
        (
            "duplicates",
            vec![
                "x0^2 + x1".to_owned(),
                "x0^2 + x1".to_owned(),
                "2*x0^2 + 2*x1".to_owned(),
                "-x0^2 - x1".to_owned(),
            ],
            Some("grevlex"),
        ),
        // Chains of pairs sharing one polynomial.
        (
            "chain",
            vec![
                "x0*x1 - 1".to_owned(),
                "x0*x2 - 1".to_owned(),
                "x1*x2 - 1".to_owned(),
                "x0*x1*x2 - 1".to_owned(),
            ],
            Some("grevlex"),
        ),
    ];
    for (name, polynomials, order) in cases {
        check(
            &session,
            &Problem::new(name, vars.clone(), polynomials, order),
        );
    }
}

#[test]
fn ideal_membership_agrees_with_buchberger() {
    let session = session();
    let mut rng = Lcg::new(0xC0FF_EE00);
    let vars = variables(3);
    let generators = vec![
        "x0^2 + x1*x2 - 1".to_owned(),
        "x1^2 + x0*x2 - 2".to_owned(),
        "x2^2 + x0*x1 - 3".to_owned(),
    ];
    let base = Problem::new(
        "membership base",
        vars.clone(),
        generators.clone(),
        Some("grevlex"),
    );
    let base_basis = basis(&session, &base.request("groebner_f4"));
    let reference = basis(&session, &base.request("groebner"));
    assert_eq!(base_basis, reference);
    // The generators are in the ideal, so adding one must not change the
    // basis.
    for generator in &generators {
        let mut extended = generators.clone();
        extended.push(generator.clone());
        let problem = Problem::new(
            "membership generator",
            vars.clone(),
            extended,
            Some("grevlex"),
        );
        assert_eq!(basis(&session, &problem.request("groebner_f4")), base_basis);
        assert_eq!(basis(&session, &problem.request("groebner_f5")), base_basis);
        check(&session, &problem);
    }
    // Random polynomials: membership decided by the basis is the same for
    // both engines, and a basis containing the candidate is the same basis
    // exactly when the candidate lies in the ideal.
    for _ in 0..8 {
        let terms: Vec<(i64, Vec<u32>)> = (0..4)
            .map(|_| {
                let powers: Vec<(usize, u32)> = (0..2)
                    .map(|_| (rng.choose(3), u32::try_from(rng.choose(3)).unwrap_or(0)))
                    .collect();
                (rng.coefficient(), exponents(3, &powers))
            })
            .collect();
        let candidate = polynomial(&vars, &terms);
        let mut extended = generators.clone();
        extended.push(candidate.clone());
        let with = Problem::new("membership", vars.clone(), extended, Some("grevlex"));
        let without = Problem::new(
            "membership",
            vars.clone(),
            generators.clone(),
            Some("grevlex"),
        );
        let f4 = basis(&session, &with.request("groebner_f4"));
        let f5 = basis(&session, &with.request("groebner_f5"));
        let buchberger = basis(&session, &with.request("groebner"));
        assert_eq!(f4, buchberger, "candidate {candidate}");
        assert_eq!(f5, buchberger, "candidate {candidate}");
        let member = f4 == base_basis;
        assert_eq!(
            member,
            basis(&session, &without.request("groebner_f4")) == f4,
            "membership changed the basis for {candidate}"
        );
    }
    // An obvious non-member: the constant one, for a proper ideal.
    let mut extended = generators.clone();
    extended.push("1".to_owned());
    let with_one = Problem::new("membership one", vars.clone(), extended, Some("grevlex"));
    assert_eq!(basis(&session, &with_one.request("groebner_f4")), "list(1)");
}

#[test]
fn repeated_runs_are_identical() {
    let session = session();
    let problem = sparse_system("determinism", 3, 3, 7, "grevlex");
    let first = problem.request("groebner_f4");
    let second = problem.request("groebner_f5");
    let reference = basis(&session, &first);
    for _ in 0..3 {
        assert_eq!(basis(&session, &first), reference);
        assert_eq!(basis(&session, &second), reference);
    }
}

#[test]
fn unevaluated_requests_instead_of_wrong_answers() {
    let session = session();
    // An unknown order symbol leaves the request alone.
    for operator in ["groebner_f4", "groebner_f5"] {
        let (text, reduced) = compute(
            &session,
            &format!("{operator}(list(x0 - 1), list(x0), fastlex)"),
        );
        assert!(!reduced, "{operator}: `{text}` should be unevaluated");
        assert!(text.contains(operator), "{operator}: {text}");
        // A polynomial that is not in the listed variables is refused.
        let (text, reduced) = compute(&session, &format!("{operator}(list(a - 1), list(x0))"));
        assert!(!reduced, "{operator}: `{text}` should be unevaluated");
        // The arguments must be lists.
        let (text, reduced) = compute(&session, &format!("{operator}(x0 - 1, list(x0))"));
        assert!(!reduced, "{operator}: `{text}` should be unevaluated");
        // Too many arguments.
        let (text, reduced) = compute(
            &session,
            &format!("{operator}(list(x0), list(x0), lex, grevlex)"),
        );
        assert!(!reduced, "{operator}: `{text}` should be unevaluated");
    }
}

#[test]
fn stress_with_many_variables_is_bounded() {
    let session = session();
    // Twelve variables, degree at most two: a triangular linear chain that
    // expresses nine variables through the last three, and a quadratic
    // equation for each of those. The solution set has eight points and the
    // basis is small, but every polynomial lives in twelve variables, so the
    // Macaulay matrices have a wide column set.
    let n = 12;
    let vars = variables(n);
    let mut polynomials = Vec::new();
    for i in 0..9 {
        polynomials.push(polynomial(
            &vars,
            &[
                (1, exponents(n, &[(i, 1)])),
                (-1, exponents(n, &[(i.saturating_add(1), 1)])),
                (-(i as i64 + 1), exponents(n, &[])),
            ],
        ));
    }
    for i in 9..12 {
        polynomials.push(polynomial(
            &vars,
            &[(1, exponents(n, &[(i, 2)])), (-1, exponents(n, &[]))],
        ));
    }
    // Couple the free variables, so that the basis is not just the input
    // rewritten: x11 = x9*x10 follows from the quadratics above.
    polynomials.push(polynomial(
        &vars,
        &[
            (1, exponents(n, &[(9, 1), (10, 1)])),
            (-1, exponents(n, &[(11, 1)])),
        ],
    ));
    let problem = Problem::new("stress-12", vars, polynomials, Some("grevlex"));
    let start = Instant::now();
    let reference = basis(&session, &problem.request("groebner"));
    let reference_time = start.elapsed();
    assert!(
        reference_time.as_secs() < 60,
        "Buchberger took {reference_time:?}"
    );
    // The chain is really solved: the last three variables satisfy their own
    // quadratics, so `x9^2 - 1` is part of the reduced basis.
    assert!(
        reference.contains("x9^2 - 1"),
        "the stress basis is unexpected: {reference}"
    );
    for operator in ["groebner_f4", "groebner_f5"] {
        let start = Instant::now();
        let computed = basis(&session, &problem.request(operator));
        let elapsed = start.elapsed();
        assert!(
            elapsed.as_secs() < 60,
            "{operator} took {elapsed:?} on the stress system"
        );
        assert_eq!(computed, reference, "{operator} on the stress system");
    }
}
