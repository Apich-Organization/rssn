//! Quantum mechanics: angular momentum algebra, exact eigenfunctions and
//! approximation methods.
//!
//! | operator | value |
//! |---|---|
//! | `clebsch_gordan(j1, m1, j2, m2, J, M)` | `⟨j1 m1; j2 m2 | J M⟩` exactly (Racah's formula; half-integers allowed), as `± √q` |
//! | `wigner_3j(j1, j2, j3, m1, m2, m3)` | the Wigner 3j symbol, from the Clebsch–Gordan coefficient |
//! | `spherical_harmonic(l, m, theta, phi)` | `Y_l^m` with the Condon–Shortley phase, as an explicit trigonometric expression |
//! | `hydrogen_radial(n, l, r)` | the normalised radial function `R_nl(r)` (Bohr radius `a_0`) |
//! | `hydrogen_wavefunction(n, l, m, r, theta, phi)` | `R_nl(r) Y_l^m(θ, φ)` |
//! | `hydrogen_energy(n)` | `-ħ²/(2 m_e a_0² n²)` |
//! | `second_order_energy_correction(V, states, energies, n, x)` | `Σ_(m≠n) |⟨m|V|n⟩|² / (E_n - E_m)` over a basis `states` (index `n` from 1) |
//! | `variational_ground_state(H, trial, x, a)` | `list(a*, E(a*))`: the Rayleigh–Ritz bound `E(a) = ⟨ψ_a|H|ψ_a⟩/⟨ψ_a|ψ_a⟩` minimised over `a > 0` |
//! | `wkb_energy(c, nu, m, n)` | Bohr–Sommerfeld level `n` of `V = c |x|^ν`: `∮ p dx = 2πħ(n + 1/2)`, in closed form with the beta function |

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::Facts;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::graph::op::core;
use crate::graph::rule::Installer;

/// The quantum-mechanics rule set.
#[must_use]
pub fn physics_quantum() -> RuleSet {
    RuleSet::new("physics_quantum", install).needs(super::physics::physics())
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Request {
    ClebschGordan,
    Wigner3j,
    SphericalHarmonic,
    HydrogenRadial,
    SecondOrder,
    Variational,
    Wkb,
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, arity, request) in [
        ("clebsch_gordan", 6, Request::ClebschGordan),
        ("wigner_3j", 6, Request::Wigner3j),
        ("spherical_harmonic", 4, Request::SphericalHarmonic),
        ("hydrogen_radial", 3, Request::HydrogenRadial),
        ("second_order_energy_correction", 5, Request::SecondOrder),
        ("variational_ground_state", 4, Request::Variational),
        ("wkb_energy", 4, Request::Wkb),
    ] {
        let op = i.op(OpDescriptor::new(name, Arity::Fixed(arity)).flags(OpFlags::HEAVY).cost(100))?;
        i.kernel(&format!("physics_quantum/{name}"), Tier::Reduce, Quantum { op, request });
    }
    i.define(&[
        "hydrogen_wavefunction(n, l, m, r, theta, phi) := hydrogen_radial(n, l, r) * spherical_harmonic(l, m, theta, phi)",
        "hydrogen_energy(n) := -hbar^2/(2*m_e*a_0^2*n^2)",
    ])?;
    Ok(())
}

struct Quantum {
    op: OpId,
    request: Request,
}

impl Kernel for Quantum {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        let result = match self.request {
            | Request::ClebschGordan => rationals(cx.graph, &args).and_then(|v| {
                let c = clebsch_gordan(&v[0], &v[1], &v[2], &v[3], &v[4], &v[5])?;
                Some(signed_root(cx.graph, &c))
            }),
            | Request::Wigner3j => rationals(cx.graph, &args).and_then(|v| {
                let c = wigner_3j(&v[0], &v[1], &v[2], &v[3], &v[4], &v[5])?;
                Some(signed_root(cx.graph, &c))
            }),
            | Request::SphericalHarmonic => spherical_harmonic(cx, &args),
            | Request::HydrogenRadial => hydrogen_radial(cx, &args),
            | Request::SecondOrder => second_order(cx, &args),
            | Request::Variational => variational(cx, &args),
            | Request::Wkb => wkb(cx, &args),
        };
        result.map_or(Outcome::Pass, |r| {
            let r = cx.simplify(r);
            Outcome::Equal(r)
        })
    }
}

fn rationals(
    graph: &Graph,
    args: &[NodeId],
) -> Option<Vec<BigRational>> {
    args.iter().map(|&a| graph.number_of(a)?.to_rational()).collect()
}

/// A value `s √q` stored as `(s, q)` with `s ∈ {-1, 0, 1}`.
struct SignedRoot {
    sign: i8,
    square: BigRational,
}

fn signed_root(
    graph: &mut Graph,
    v: &SignedRoot,
) -> NodeId {
    if v.sign == 0 {
        return graph.int(0);
    }
    let q = graph.num(Number::rat(v.square.clone()));
    let half = graph.num(Number::fraction(1, 2).unwrap_or_else(|| Number::from(0)));
    let root = graph.node(core::POW, &[q, half]);
    if v.sign < 0 {
        let minus = graph.int(-1);
        graph.node(core::MUL, &[minus, root])
    } else {
        root
    }
}

fn factorial(n: &BigRational) -> Option<BigRational> {
    if !n.is_integer() || n.is_negative() {
        return None;
    }
    let n = n.to_integer().to_u64()?;
    if n > 500 {
        return None;
    }
    Some(BigRational::from_integer((1..=n).fold(BigInt::one(), |acc, k| acc * k)))
}

fn is_half_integer(v: &BigRational) -> bool {
    (v * BigRational::from_integer(BigInt::from(2))).is_integer()
}

/// `⟨j1 m1; j2 m2 | J M⟩` by Racah's formula.
fn clebsch_gordan(
    j1: &BigRational,
    m1: &BigRational,
    j2: &BigRational,
    m2: &BigRational,
    j: &BigRational,
    m: &BigRational,
) -> Option<SignedRoot> {
    if ![j1, m1, j2, m2, j, m].iter().all(|v| is_half_integer(v)) {
        return None;
    }
    let zero = SignedRoot { sign: 0, square: BigRational::zero() };
    // Selection rules.
    if m1 + m2 != *m || m1.abs() > *j1 || m2.abs() > *j2 || m.abs() > *j || !(j1 - m1).is_integer() || !(j2 - m2).is_integer() || !(j - m).is_integer() {
        return Some(zero);
    }
    if *j > j1 + j2 || *j < (j1 - j2).abs() || !(j1 + j2 + j).is_integer() {
        return Some(zero);
    }
    let one = BigRational::one();
    let two = BigRational::from_integer(BigInt::from(2));
    let f = |v: BigRational| factorial(&v);
    let prefactor = (&two * j + &one) * f(j + j1 - j2)? * f(j - j1 + j2)? * f(j1 + j2 - j)? / f(j1 + j2 + j + &one)?;
    let product = f(j + m)? * f(j - m)? * f(j1 - m1)? * f(j1 + m1)? * f(j2 - m2)? * f(j2 + m2)?;
    let mut sum = BigRational::zero();
    let mut k = BigRational::zero();
    loop {
        let parts = [
            k.clone(),
            j1 + j2 - j - &k,
            j1 - m1 - &k,
            j2 + m2 - &k,
            j - j2 + m1 + &k,
            j - j1 - m2 + &k,
        ];
        if parts[1].is_negative() || parts[2].is_negative() || parts[3].is_negative() {
            break;
        }
        if !parts[4].is_negative() && !parts[5].is_negative() {
            let mut denominator = BigRational::one();
            for p in &parts {
                denominator *= f(p.clone())?;
            }
            let term = BigRational::one() / denominator;
            if k.to_integer().to_i64()? % 2 == 0 { sum += term } else { sum -= term }
        }
        k += &one;
    }
    if sum.is_zero() {
        return Some(zero);
    }
    let sign = if sum.is_negative() { -1 } else { 1 };
    Some(SignedRoot { sign, square: &sum * &sum * prefactor * product })
}

/// The Wigner 3j symbol `(j1 j2 j3; m1 m2 m3) = (-1)^(j1-j2-m3) ⟨j1 m1; j2
/// m2 | j3 -m3⟩ / √(2 j3 + 1)`.
fn wigner_3j(
    j1: &BigRational,
    j2: &BigRational,
    j3: &BigRational,
    m1: &BigRational,
    m2: &BigRational,
    m3: &BigRational,
) -> Option<SignedRoot> {
    let c = clebsch_gordan(j1, m1, j2, m2, j3, &-m3)?;
    let exponent = j1 - j2 - m3;
    if !exponent.is_integer() {
        return None;
    }
    let flip = exponent.to_integer().to_i64()?.rem_euclid(2) == 1;
    let scale = BigRational::from_integer(BigInt::from(2)) * j3 + BigRational::one();
    Some(SignedRoot { sign: if flip { -c.sign } else { c.sign }, square: c.square / scale })
}

/// Coefficients (ascending) of the Legendre polynomial `P_l` by Rodrigues'
/// formula.
fn legendre(l: usize) -> Vec<BigRational> {
    // (x² - 1)^l
    let mut p = vec![BigRational::one()];
    for _ in 0..l {
        let mut next = vec![BigRational::zero(); p.len() + 2];
        for (k, c) in p.iter().enumerate() {
            next[k + 2] += c;
            next[k] -= c;
        }
        p = next;
    }
    for _ in 0..l {
        p = derivative(&p);
    }
    let scale = BigRational::from_integer(BigInt::from(2).pow(u32::try_from(l).unwrap_or(0)) * (1..=l).fold(BigInt::one(), |a, k| a * k));
    p.into_iter().map(|c| c / &scale).collect()
}

fn derivative(p: &[BigRational]) -> Vec<BigRational> {
    p.iter().enumerate().skip(1).map(|(k, c)| c * BigRational::from_integer(BigInt::from(k))).collect()
}

fn polynomial_in(
    graph: &mut Graph,
    coefficients: &[BigRational],
    x: NodeId,
) -> NodeId {
    let mut terms = Vec::new();
    for (k, c) in coefficients.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let c = graph.num(Number::rat(c.clone()));
        let e = graph.int(i64::try_from(k).unwrap_or(0));
        let power = graph.node(core::POW, &[x, e]);
        terms.push(graph.node(core::MUL, &[c, power]));
    }
    if terms.is_empty() { graph.int(0) } else { graph.node(core::ADD, &terms) }
}

fn call(
    graph: &mut Graph,
    name: &str,
    args: &[NodeId],
) -> Option<NodeId> {
    let op = graph.ops().lookup(name)?;
    graph.try_node(op, args)
}

/// `Y_l^m(θ, φ)` with the Condon–Shortley phase.
fn spherical_harmonic(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let graph = &mut *cx.graph;
    let &[l, m, theta, phi] = args else {
        return None;
    };
    let (l, m) = (graph.number_of(l)?.to_i64()?, graph.number_of(m)?.to_i64()?);
    if l < 0 || l > 40 {
        return None;
    }
    if m.abs() > l {
        return Some(graph.int(0));
    }
    let (lu, mu) = (usize::try_from(l).ok()?, usize::try_from(m.unsigned_abs()).ok()?);
    let mut p = legendre(lu);
    for _ in 0..mu {
        p = derivative(&p);
    }
    let fact = |n: usize| (1..=n).fold(BigInt::one(), |a, k| a * k);
    // N² = (2l + 1)/(4π) (l - |m|)!/(l + |m|)!
    let n2 = BigRational::new(BigInt::from(2 * l + 1) * fact(lu - mu), BigInt::from(4) * fact(lu + mu));
    let n2 = graph.num(Number::rat(n2));
    let pi = call(graph, "pi", &[])?;
    let minus_one = graph.int(-1);
    let inv_pi = graph.node(core::POW, &[pi, minus_one]);
    let half = graph.num(Number::fraction(1, 2)?);
    let norm_sq = graph.node(core::MUL, &[n2, inv_pi]);
    let norm = graph.node(core::POW, &[norm_sq, half]);
    let cos_t = call(graph, "cos", &[theta])?;
    let sin_t = call(graph, "sin", &[theta])?;
    let poly = polynomial_in(graph, &p, cos_t);
    let e = graph.int(i64::try_from(mu).ok()?);
    let sin_power = graph.node(core::POW, &[sin_t, e]);
    let i_unit = call(graph, "I", &[])?;
    let m_node = graph.int(m);
    let phase_arg = graph.node(core::MUL, &[i_unit, m_node, phi]);
    let phase = call(graph, "exp", &[phase_arg])?;
    // Condon–Shortley (-1)^m for m > 0; Y_l^(-m) = (-1)^m conj(Y_l^m)
    // carries no extra sign.
    let sign = graph.int(if m > 0 && m % 2 == 1 { -1 } else { 1 });
    let mut factors = vec![sign, norm, poly, sin_power];
    if m != 0 {
        factors.push(phase);
    }
    Some(graph.node(core::MUL, &factors))
}

/// `R_nl(r) = √((2/(n a))³ (n-l-1)!/(2n (n+l)!)) e^(-ρ/2) ρ^l L_(n-l-1)^(2l+1)(ρ)`,
/// `ρ = 2r/(n a)`.
fn hydrogen_radial(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let graph = &mut *cx.graph;
    let &[n, l, r] = args else {
        return None;
    };
    let (n, l) = (graph.number_of(n)?.to_i64()?, graph.number_of(l)?.to_i64()?);
    if n < 1 || l < 0 || l >= n || n > 40 {
        return None;
    }
    let fact = |k: i64| (1..=k).fold(BigInt::one(), |a, j| a * j);
    let binomial = |a: i64, b: i64| if b < 0 || b > a { BigInt::zero() } else { fact(a) / (fact(b) * fact(a - b)) };
    let (k, alpha) = (n - l - 1, 2 * l + 1);
    // L_k^α(ρ) = Σ_i (-1)^i C(k + α, k - i) ρ^i / i!
    let laguerre: Vec<BigRational> = (0..=k)
        .map(|i| {
            let c = BigRational::new(binomial(k + alpha, k - i), fact(i));
            if i % 2 == 1 { -c } else { c }
        })
        .collect();
    let a0 = call(graph, "a_0", &[])?;
    let n_node = graph.int(n);
    let two = graph.int(2);
    let na = graph.node(core::MUL, &[n_node, a0]);
    let minus_one = graph.int(-1);
    let inv_na = graph.node(core::POW, &[na, minus_one]);
    let rho = graph.node(core::MUL, &[two, r, inv_na]);
    // norm² = 8 (n-l-1)! / (n³ 2n (n+l)!) · a^-3
    let norm_sq = BigRational::new(BigInt::from(8) * fact(k), BigInt::from(n).pow(3) * BigInt::from(2 * n) * fact(n + l));
    let norm_sq = graph.num(Number::rat(norm_sq));
    let half = graph.num(Number::fraction(1, 2)?);
    let norm = graph.node(core::POW, &[norm_sq, half]);
    let minus_three_halves = graph.num(Number::fraction(-3, 2)?);
    let a_power = graph.node(core::POW, &[a0, minus_three_halves]);
    let minus_half = graph.num(Number::fraction(-1, 2)?);
    let decay_arg = graph.node(core::MUL, &[minus_half, rho]);
    let decay = call(graph, "exp", &[decay_arg])?;
    let l_node = graph.int(l);
    let rho_l = graph.node(core::POW, &[rho, l_node]);
    let poly = polynomial_in(graph, &laguerre, rho);
    Some(graph.node(core::MUL, &[norm, a_power, decay, rho_l, poly]))
}

/// `Σ_(m≠n) |⟨m|V|n⟩|² / (E_n - E_m)`.
fn second_order(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[v, states, energies, n, x] = args else {
        return None;
    };
    if cx.graph.op(states) != core::LIST || cx.graph.op(energies) != core::LIST {
        return None;
    }
    let states = cx.graph.children(states).to_vec();
    let energies = cx.graph.children(energies).to_vec();
    let n = usize::try_from(cx.graph.number_of(n)?.to_i64()?).ok()?.checked_sub(1)?;
    if states.len() != energies.len() || n >= states.len() {
        return None;
    }
    let mut terms = Vec::new();
    for (m, (&psi_m, &e_m)) in states.iter().zip(&energies).enumerate() {
        if m == n {
            continue;
        }
        let applied = call(cx.graph, "qm_apply", &[v, states[n]])?;
        let element = call(cx.graph, "braket", &[psi_m, applied, x])?;
        let element = cx.simplify(element);
        if cx.graph.number_of(element).is_some_and(Number::is_zero) {
            continue;
        }
        let conj = call(cx.graph, "conj", &[element])?;
        let numerator = cx.graph.node(core::MUL, &[element, conj]);
        let minus = cx.graph.int(-1);
        let neg_m = cx.graph.node(core::MUL, &[minus, e_m]);
        let gap = cx.graph.node(core::ADD, &[energies[n], neg_m]);
        let inv = cx.graph.node(core::POW, &[gap, minus]);
        terms.push(cx.graph.node(core::MUL, &[numerator, inv]));
    }
    Some(if terms.is_empty() { cx.graph.int(0) } else { cx.graph.node(core::ADD, &terms) })
}

/// Rayleigh–Ritz: minimise `E(a)` over `a > 0`.
fn variational(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[h, trial, x, a] = args else {
        return None;
    };
    let a_symbol = cx.graph.symbol_of(a)?;
    cx.graph.assume(a_symbol, Facts::POSITIVE);
    // E(a) = ∫ ψ Hψ / ∫ ψ² for a real trial function.
    let oo = call(cx.graph, "oo", &[])?;
    let minus_one = cx.graph.int(-1);
    let minus_oo = cx.graph.node(core::MUL, &[minus_one, oo]);
    let applied = call(cx.graph, "qm_apply", &[h, trial])?;
    let applied = cx.simplify(applied);
    let numerator = cx.graph.node(core::MUL, &[trial, applied]);
    let numerator = cx.simplify(numerator);
    let numerator = crate::rules::poly::expand_form(cx.graph, numerator).unwrap_or(numerator);
    let numerator = cx.simplify(numerator);
    let numerator = call(cx.graph, "defint", &[numerator, x, minus_oo, oo])?;
    let two = cx.graph.int(2);
    let square = cx.graph.node(core::POW, &[trial, two]);
    let denominator = call(cx.graph, "defint", &[square, x, minus_oo, oo])?;
    let (numerator, denominator) = (cx.simplify(numerator), cx.simplify(denominator));
    let inv = cx.graph.node(core::POW, &[denominator, minus_one]);
    let energy = cx.graph.node(core::MUL, &[numerator, inv]);
    let energy = cx.simplify(energy);
    if contains_op(cx.graph, energy, "defint") || contains_op(cx.graph, energy, "qm_apply") {
        return None;
    }
    let d = crate::rules::calculus::derivative(cx.graph, energy, a)?;
    let d = cx.simplify(d);
    let zero = cx.graph.int(0);
    let equation = cx.graph.node(core::EQ, &[d, zero]);
    let request = call(cx.graph, "solve", &[equation, a])?;
    let roots = cx.simplify(request);
    if cx.graph.op(roots) != core::LIST {
        return None;
    }
    let mut best: Option<(f64, NodeId, NodeId)> = None;
    for root in cx.graph.children(roots).to_vec() {
        let mut env = Env::numeric(0.0);
        for s in cx.graph.free_symbols(cx.graph.find(root)).to_vec() {
            env.bind(s, 1.0);
        }
        let Some(value) = cx.graph.eval(root, &env).filter(|v| *v > 0.0) else {
            continue;
        };
        let _ = value;
        let at = cx.graph.substitute(energy, a, root);
        let at = cx.simplify(at);
        let mut env = Env::numeric(0.0);
        for s in cx.graph.free_symbols(cx.graph.find(at)).to_vec() {
            env.bind(s, 1.0);
        }
        let e_value = cx.graph.eval(at, &env).unwrap_or(f64::INFINITY);
        if best.as_ref().is_none_or(|b| e_value < b.0) {
            best = Some((e_value, root, at));
        }
    }
    let (_, root, at) = best?;
    Some(cx.graph.node(core::LIST, &[root, at]))
}

fn contains_op(
    graph: &Graph,
    node: NodeId,
    name: &str,
) -> bool {
    let Some(op) = graph.ops().lookup(name) else {
        return false;
    };
    let mut stack = vec![node];
    while let Some(n) = stack.pop() {
        if graph.op(n) == op {
            return true;
        }
        stack.extend_from_slice(graph.children(n));
    }
    false
}

/// Bohr–Sommerfeld levels of `V = c |x|^ν`:
/// `∫ √(2m(E - c|x|^ν)) dx = 2 √(2m) E^(1/2 + 1/ν) c^(-1/ν) B(1/ν, 3/2)/ν`
/// over the classically allowed interval, set equal to `πħ(n + 1/2)`.
fn wkb(
    cx: &mut Cx<'_>,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[c, nu, m, n] = args else {
        return None;
    };
    let graph = &mut *cx.graph;
    let hbar = call(graph, "hbar", &[])?;
    let pi = call(graph, "pi", &[])?;
    let one = graph.int(1);
    let half = graph.num(Number::fraction(1, 2)?);
    let three_halves = graph.num(Number::fraction(3, 2)?);
    let minus_one = graph.int(-1);
    let inv_nu = graph.node(core::POW, &[nu, minus_one]);
    let beta = call(graph, "beta", &[inv_nu, three_halves])?;
    let n_half = graph.node(core::ADD, &[n, half]);
    // E^((ν+2)/(2ν)) = π ħ (n + 1/2) ν c^(1/ν) / (2 √(2m) B)
    let two = graph.int(2);
    let two_m = graph.node(core::MUL, &[two, m]);
    let root_2m = graph.node(core::POW, &[two_m, half]);
    let c_power = graph.node(core::POW, &[c, inv_nu]);
    let denominator = graph.node(core::MUL, &[two, root_2m, beta]);
    let inv_den = graph.node(core::POW, &[denominator, minus_one]);
    let rhs = graph.node(core::MUL, &[pi, hbar, n_half, nu, c_power, inv_den]);
    // exponent 2ν/(ν + 2)
    let nu_plus_2 = graph.node(core::ADD, &[nu, two]);
    let inv = graph.node(core::POW, &[nu_plus_2, minus_one]);
    let exponent = graph.node(core::MUL, &[two, nu, inv]);
    let _ = one;
    Some(graph.node(core::POW, &[rhs, exponent]))
}

#[cfg(test)]
mod tests {
    use crate::rules::testing::eval;
    use crate::rules::testing::simplify;

    fn rules() -> Vec<crate::graph::RuleSet> {
        crate::rules::standard()
    }

    #[test]
    fn clebsch_gordan_and_3j() {
        let rules = rules();
        let run = |s: &str| simplify(&rules, s);
        // Two spin-1/2: |1 0> = (|↑↓> + |↓↑>)/√2, singlet antisymmetric.
        assert_eq!(run("clebsch_gordan(1/2, 1/2, 1/2, -1/2, 1, 0)"), "1/2^(1/2)");
        assert_eq!(run("clebsch_gordan(1/2, -1/2, 1/2, 1/2, 0, 0)"), "-1/2^(1/2)");
        assert_eq!(run("clebsch_gordan(1, 1, 1/2, -1/2, 3/2, 1/2)"), "1/3^(1/2)");
        assert_eq!(run("clebsch_gordan(1, 0, 1, 0, 1, 0)"), "0");
        assert_eq!(run("clebsch_gordan(2, 1, 1, -1, 2, 0)"), "1/2^(1/2)");
        // Orthonormality: Σ_M1 m2 |C|² over a column is 1.
        let total: f64 = [(-1.0, 1.0), (0.0, 0.0), (1.0, -1.0)]
            .iter()
            .map(|(m1, m2)| {
                let v = simplify(&rules, &format!("clebsch_gordan(1, {m1}, 1, {m2}, 2, 0)^2"));
                eval(&rules, &v, &[])
            })
            .sum();
        assert!((total - 1.0).abs() < 1e-12);
        let three_j = simplify(&rules, "wigner_3j(1, 1, 0, 0, 0, 0)");
        assert!((eval(&rules, &three_j, &[]) + 1.0 / 3.0_f64.sqrt()).abs() < 1e-12, "{three_j}");
    }

    #[test]
    fn spherical_harmonics_and_hydrogen() {
        let rules = rules();
        let y10 = simplify(&rules, "spherical_harmonic(1, 0, theta, phi)");
        let want = (3.0 / (4.0 * std::f64::consts::PI)).sqrt() * 0.7_f64.cos();
        assert!((eval(&rules, &y10, &[("theta", 0.7), ("phi", 0.3)]) - want).abs() < 1e-12, "{y10}");
        // Y_2^0 integrates to 1 in |Y|² over the sphere (azimuth trivial).
        let y20 = simplify(&rules, "spherical_harmonic(2, 0, theta, phi)");
        let norm = simplify(&rules, &format!("2*pi*defint(({y20})^2*sin(theta), theta, 0, pi)"));
        assert!((eval(&rules, &norm, &[]) - 1.0).abs() < 1e-10, "{norm}");
        // R_10 = 2 a^(-3/2) e^(-r/a); ∫ R_21² r² dr = 1.
        let r10 = simplify(&rules, "hydrogen_radial(1, 0, r)");
        assert!((eval(&rules, &r10, &[("r", 0.4)]) - 2.0 * (-0.4 / 5.291_772_109e-11_f64).exp() * 5.291_772_109e-11_f64.powf(-1.5)).abs() / 1e15 < 1e-3);
        let norm = simplify(&rules, "defint(hydrogen_radial(2, 1, r)^2*r^2, r, 0, oo)");
        assert_eq!(norm, "1");
        let e = simplify(&rules, "hydrogen_energy(2)/hydrogen_energy(1)");
        assert_eq!(e, "1/4");
    }

    #[test]
    fn approximation_methods() {
        let rules = rules();
        // Gaussian trial for the harmonic oscillator is exact: E = ħω/2.
        let ground = simplify(&rules, "variational_ground_state(hamiltonian_harmonic_oscillator(m, w, x), exp(-a*x^2), x, a)");
        // The second item of `list(a*, E*)`, split at the top-level comma.
        let inner = ground.trim_start_matches("list(").trim_end_matches(')');
        let mut depth = 0;
        let split = inner
            .char_indices()
            .find(|&(_, c)| {
                match c {
                    | '(' => depth += 1,
                    | ')' => depth -= 1,
                    | _ => {},
                }
                c == ',' && depth == 0
            })
            .map(|(i, _)| i)
            .unwrap();
        let (a_star, e_star) = (&inner[..split], inner[split + 1..].trim());
        // ħ is the CODATA constant.
        let hbar = crate::constant::REDUCED_PLANCK_CONSTANT;
        let at = [("m", 1.3), ("w", 0.7)];
        assert!((eval(&rules, e_star, &at) / (0.5 * hbar * 0.7) - 1.0).abs() < 1e-10, "{ground}");
        assert!((eval(&rules, a_star, &at) / (1.3 * 0.7 / (2.0 * hbar)) - 1.0).abs() < 1e-10, "{ground}");
        // WKB for the oscillator is exact: E_n = ħω(n + 1/2) with c = m ω²/2.
        let wkb = simplify(&rules, "wkb_energy(m*w^2/2, 2, m, n)/(hbar*w*(n + 1/2))");
        assert!((eval(&rules, &wkb, &[("m", 1.3), ("w", 0.7), ("n", 2.0)]) - 1.0).abs() < 1e-10, "{wkb}");
        // Second-order correction for x perturbing the oscillator ground
        // state with a two-state basis: -|⟨1|x|0⟩|²/ħω.
        let corr = simplify(
            &rules,
            "second_order_energy_correction(position_operator(x), list(exp(-x^2/2)/pi^(1/4), 2^(1/2)*x*exp(-x^2/2)/pi^(1/4)), list(1/2, 3/2), 1, x)",
        );
        assert!((eval(&rules, &corr, &[]) + 0.5).abs() < 1e-10, "{corr}");
    }
}
