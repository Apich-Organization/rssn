//! Further integration methods (a child of the integrator, sharing its
//! state and verification).
//!
//! * **Reciprocal trigonometric functions** `sec`, `csc`, `cot` and their
//!   squares.
//! * **Rational functions of `e^{a x + b}`**: `u = e^{a x + b}` turns them
//!   into rational functions of `u`.
//! * **Rational functions of `sin`, `cos`, `tan` of one linear argument**:
//!   the Weierstrass substitution `t = tan(u/2)`.
//! * **`R(x, √Q)` for a quadratic `Q`** (Euler's algebraic method):
//!   `R = A(x) + B(x) √Q` after rationalising the denominator; `∫ B √Q =
//!   ∫ T/√Q` with `T = B Q` split by partial fractions into a polynomial
//!   part — `∫ P/√Q = S √Q + λ ∫ dx/√Q` with `S` from one linear system —
//!   and pole parts `∫ dx/((x - α)^m √Q)`, which `t = 1/(x - α)` turns into
//!   the polynomial case.
//! * **A Risch–Norman ansatz** for one transcendental extension
//!   `θ = e^{g(x)}` (`g` polynomial) or `θ = ln(x)`: `F = N(x, θ)/H(x)` (plus
//!   logarithms of the factors of the denominator) with `H` the part of the
//!   denominator that survives integration; the coefficients of `N` come
//!   from the linear identity `F' = f` in `x` and `θ`.
//! * **Biquadratic denominators** `x⁴ + c x² + e` irreducible over `Q`:
//!   the real factorisation `(x² + s x + t)(x² - s x + t)` with `t = √e`,
//!   `s = √(2t - c)`, closed-form partial fractions, and the parametric
//!   quadratic integrals.
//! * **Definite integrals by residues**: `∫_{-∞}^{∞} R(x) {1, cos(a x),
//!   sin(a x)} dx` for rational `R` whose poles are simple and come from
//!   quadratic factors over `Q`, by Jordan's lemma and the residues in the
//!   upper half-plane; and Dirichlet's `∫_0^∞ sin(a x)/x dx = π/2 sign(a)`.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::Zero;

use super::Integrator;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::op::core;
use crate::rules::poly::ratio;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;
use crate::rules::poly::repr::Poly;
use crate::rules::poly::repr::to_term;
use crate::rules::poly::univariate;
use crate::rules::poly::univariate::QPoly;

/// Solves a square or overdetermined linear system over `Q` (columns =
/// unknowns), `None` if inconsistent; free unknowns are set to zero.
fn solve_q(
    mut rows: Vec<Vec<BigRational>>,
    mut rhs: Vec<BigRational>,
    n: usize,
) -> Option<Vec<BigRational>> {
    let m = rows.len();
    let mut pivots = Vec::new();
    let mut r = 0;
    for c in 0..n {
        let Some(p) = (r..m).find(|&i| !rows[i][c].is_zero()) else {
            continue;
        };
        rows.swap(r, p);
        rhs.swap(r, p);
        let lead = rows[r][c].clone();
        for v in &mut rows[r] {
            *v /= &lead;
        }
        rhs[r] /= &lead;
        for i in 0..m {
            if i != r && !rows[i][c].is_zero() {
                let factor = rows[i][c].clone();
                let pivot_row = rows[r].clone();
                for (v, pv) in rows[i].iter_mut().zip(&pivot_row) {
                    *v -= &factor * pv;
                }
                let pr = rhs[r].clone();
                rhs[i] -= &factor * pr;
            }
        }
        pivots.push(c);
        r += 1;
    }
    if rhs[r..].iter().any(|v| !v.is_zero()) {
        return None;
    }
    let mut x = vec![BigRational::zero(); n];
    for (i, &c) in pivots.iter().enumerate() {
        x[c] = rhs[i].clone();
    }
    Some(x)
}

impl Integrator<'_, '_> {
    fn q_term(
        &mut self,
        q: &[BigRational],
    ) -> NodeId {
        self.polynomial(q)
    }

    fn q_poly_of(
        &mut self,
        f: NodeId,
    ) -> Option<QPoly> {
        let mut gens = Gens::default();
        let x = self.x;
        let gx = gens.index(self.cx.graph, x);
        let p = crate::rules::poly::repr::from_term(self.cx.graph, &mut gens, f, Limits::default())?;
        if p.support().iter().any(|&g| g != gx) {
            return None;
        }
        let mut q: QPoly = p.univariate_in(gx)?.iter().map(Number::to_rational).collect::<Option<_>>()?;
        while q.last().is_some_and(Zero::is_zero) {
            q.pop();
        }
        Some(q)
    }

    /// `sec`, `csc`, `cot` and their squares of a linear argument.
    pub(super) fn reciprocal_trig(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        let graph = &*self.cx.graph;
        let (sec, csc, cot) = (graph.ops().lookup("sec"), graph.ops().lookup("csc"), graph.ops().lookup("cot"));
        let op = graph.op(f);
        let children = graph.children(f).to_vec();
        let (inner_op, squared, u) = match children.as_slice() {
            | &[u] => (op, false, u),
            | &[base, e] if op == core::POW && self.number(e).is_some_and(|n| n == Number::from(2)) => {
                let &[u] = self.cx.graph.children(base) else { return None };
                (self.cx.graph.op(base), true, u)
            },
            | _ => return None,
        };
        if Some(inner_op) != sec && Some(inner_op) != csc && Some(inner_op) != cot {
            return None;
        }
        let a = self.linear(u)?;
        let fs = self.f;
        let tan = self.call(fs.tan, u);
        let antiderivative = if Some(inner_op) == sec {
            if squared {
                tan
            } else {
                // ln(sec u + tan u)
                let s = self.cx.graph.node(sec?, &[u]);
                let sum = self.add(&[s, tan]);
                self.call(fs.ln, sum)
            }
        } else if Some(inner_op) == csc {
            if squared {
                // -cot u
                let c = self.cx.graph.node(cot?, &[u]);
                self.neg(c)
            } else {
                // ln(tan(u/2))
                let half = self.frac(1, 2);
                let hu = self.mul(&[half, u]);
                let t = self.call(fs.tan, hu);
                self.call(fs.ln, t)
            }
        } else if squared {
            // ∫ cot² = -cot u - u
            let c = self.cx.graph.node(cot?, &[u]);
            let sum = self.add(&[c, u]);
            self.neg(sum)
        } else {
            let s = self.call(fs.sin, u);
            self.call(fs.ln, s)
        };
        Some(self.div(antiderivative, a))
    }

    /// The single node of `f` with operator in `ops` that depends on `x`,
    /// if every dependence on `x` goes through such nodes with one common
    /// argument.
    fn single_argument(
        &self,
        f: NodeId,
        ops: &[crate::graph::OpId],
    ) -> Option<NodeId> {
        let x = self.x;
        let mut argument = None;
        let mut stack = vec![f];
        let mut seen = Vec::new();
        while let Some(n) = stack.pop() {
            if seen.contains(&n) {
                continue;
            }
            seen.push(n);
            if !self.depends(n) {
                continue;
            }
            if n == x {
                return None;
            }
            let op = self.cx.graph.op(n);
            if ops.contains(&op) {
                let &[u] = self.cx.graph.children(n) else { return None };
                match argument {
                    | None => argument = Some(u),
                    | Some(a) if a == u => {},
                    | Some(_) => return None,
                }
                continue;
            }
            stack.extend_from_slice(self.cx.graph.children(n));
        }
        argument
    }

    /// Rational functions of `e^u`, `u = a x + b`.
    pub(super) fn exp_substitution(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        let fs = self.f;
        let u = self.single_argument(f, &[fs.exp])?;
        let a = self.linear(u)?;
        let w = {
            let s = self.cx.graph.interner_mut().fresh_symbol("w");
            self.cx.graph.symbol_node(s)
        };
        let e = self.call(fs.exp, u);
        let in_w = self.cx.graph.replace_subterm(f, e, w);
        if self.depends(in_w) {
            return None;
        }
        // dx = dw / (a w)
        let aw = self.mul(&[a, w]);
        let integrand = self.div(in_w, aw);
        let integrand = self.cx.simplify(integrand);
        let primitive = super::antiderivative(self.cx, self.f, integrand, w)?;
        let back = self.cx.graph.substitute(primitive, w, e);
        Some(back)
    }

    fn fresh(
        &mut self,
        name: &str,
    ) -> NodeId {
        let s = self.cx.graph.interner_mut().fresh_symbol(name);
        self.cx.graph.symbol_node(s)
    }

    /// Fractional powers `x^(p/q)` (and `sqrt(x)`) of the variable itself:
    /// `x = t^L` with `L` the least common denominator makes the integrand
    /// free of radicals in `x`; `dx = L t^(L-1) dt`.
    pub(super) fn fractional_power_substitution(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        let (x, fs) = (self.x, self.f);
        let mut powers: Vec<(NodeId, BigRational)> = Vec::new();
        let mut stack = vec![f];
        while let Some(n) = stack.pop() {
            let op = self.cx.graph.op(n);
            let children = self.cx.graph.children(n).to_vec();
            match children.as_slice() {
                | &[b] if op == fs.sqrt && b == x => powers.push((n, BigRational::new(BigInt::one(), BigInt::from(2)))),
                | &[b, e] if op == core::POW && b == x => {
                    let r = self.number(e)?.to_rational()?;
                    if !r.is_integer() {
                        powers.push((n, r));
                    }
                },
                | _ => stack.extend(children),
            }
        }
        if powers.is_empty() {
            return None;
        }
        let l = powers.iter().fold(BigInt::one(), |acc, (_, r)| num_integer::Integer::lcm(&acc, r.denom()));
        let l_i = i64::try_from(&l).ok().filter(|&v| v <= 12)?;
        let t = self.fresh("t");
        let mut g = f;
        for (node, r) in &powers {
            let k = i64::try_from(&(r * BigRational::from_integer(l.clone())).to_integer()).ok()?;
            let k = self.int(k);
            let replacement = self.pow(t, k);
            g = self.cx.graph.replace_subterm(g, *node, replacement);
        }
        let l_node = self.int(l_i);
        let t_l = self.pow(t, l_node);
        g = self.cx.graph.substitute(g, x, t_l);
        let e = self.int(l_i - 1);
        let jac = self.pow(t, e);
        let integrand = self.mul(&[l_node, jac, g]);
        let integrand = self.cx.simplify(integrand);
        let primitive = super::antiderivative(self.cx, self.f, integrand, t)?;
        let inv = self.frac(1, l_i);
        let back = self.pow(x, inv);
        Some(self.cx.graph.substitute(primitive, t, back))
    }

    /// Functions of `ln x` (times powers of `x`): `x = e^t`, `dx = e^t dt`.
    pub(super) fn log_substitution(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        let (x, fs) = (self.x, self.f);
        let log_x = self.call(fs.ln, x);
        let t = self.fresh("t");
        let g = self.cx.graph.replace_subterm(f, log_x, t);
        if g == f {
            return None;
        }
        let e_t = self.call(fs.exp, t);
        let g = self.cx.graph.substitute(g, x, e_t);
        let integrand = self.mul(&[g, e_t]);
        let integrand = self.cx.simplify(integrand);
        let primitive = super::antiderivative(self.cx, self.f, integrand, t)?;
        Some(self.cx.graph.substitute(primitive, t, log_x))
    }

    /// `sinh`, `cosh`, `tanh` rewritten through `exp`.
    pub(super) fn hyperbolic_as_exponentials(
        &mut self,
        f: NodeId,
        depth: usize,
    ) -> Option<NodeId> {
        let fs = self.f;
        let mut g = f;
        let mut changed = false;
        let mut stack = vec![f];
        let mut seen = Vec::new();
        while let Some(n) = stack.pop() {
            if seen.contains(&n) {
                continue;
            }
            seen.push(n);
            let op = self.cx.graph.op(n);
            let children = self.cx.graph.children(n).to_vec();
            if let (true, &[u]) = (op == fs.sinh || op == fs.cosh || op == fs.tanh, children.as_slice()) {
                let plus = self.call(fs.exp, u);
                let minus_u = self.neg(u);
                let minus = self.call(fs.exp, minus_u);
                let neg_minus = self.neg(minus);
                let difference = self.add(&[plus, neg_minus]);
                let sum = self.add(&[plus, minus]);
                let half = self.frac(1, 2);
                let value = if op == fs.sinh {
                    self.mul(&[half, difference])
                } else if op == fs.cosh {
                    self.mul(&[half, sum])
                } else {
                    self.div(difference, sum)
                };
                g = self.cx.graph.replace_subterm(g, n, value);
                changed = true;
            } else {
                stack.extend(children);
            }
        }
        if !changed {
            return None;
        }
        let g = self.cx.simplify(g);
        let g = crate::rules::poly::expand_form(self.cx.graph, g).unwrap_or(g);
        let g = self.cx.simplify(g);
        self.integrate(g, depth + 1)
    }

    /// `∫ sec^n u` and `∫ csc^n u` (`n >= 2`, also written `cos u^(-n)`):
    /// `sec^n = sec^(n-2) tan/(n-1) + (n-2)/(n-1) ∫ sec^(n-2)`, and
    /// `csc^n = -csc^(n-2) cot/(n-1) + (n-2)/(n-1) ∫ csc^(n-2)`, over `a`.
    pub(super) fn secant_power(
        &mut self,
        f: NodeId,
        depth: usize,
    ) -> Option<NodeId> {
        let fs = self.f;
        let graph = &*self.cx.graph;
        let (sec, csc) = (graph.ops().lookup("sec"), graph.ops().lookup("csc"));
        let &[base, e] = graph.children(f) else {
            return None;
        };
        if graph.op(f) != core::POW {
            return None;
        }
        let k = graph.number_of(e)?.to_i64()?;
        let (is_cos, n) = match graph.op(base) {
            | op if op == fs.cos && k <= -2 => (true, -k),
            | op if op == fs.sin && k <= -2 => (false, -k),
            | op if Some(op) == sec && k >= 2 => (true, k),
            | op if Some(op) == csc && k >= 2 => (false, k),
            | _ => return None,
        };
        let &[u] = graph.children(base) else {
            return None;
        };
        let a = self.linear(u)?;
        let (c, s) = (self.call(fs.cos, u), self.call(fs.sin, u));
        // sec^(n-2) tan = sin / cos^(n-1);  -csc^(n-2) cot = -cos / sin^(n-1)
        let (num, den) = if is_cos { (s, c) } else { (c, s) };
        let e = self.int(1 - n);
        let raised = self.pow(den, e);
        let sign = self.int(if is_cos { 1 } else { -1 });
        let scale = self.frac(1, n - 1);
        let head = self.mul(&[sign, scale, num, raised]);
        let head = self.div(head, a);
        let lower = if n == 2 {
            None
        } else {
            let e = self.int(2 - n);
            let rest = self.pow(den, e);
            let rest = self.cx.simplify(rest);
            let inner = self.integrate(rest, depth + 1)?;
            let weight = self.frac(n - 2, n - 1);
            Some(self.mul(&[weight, inner]))
        };
        Some(match lower {
            | Some(l) => self.add(&[head, l]),
            | None => head,
        })
    }

    /// `∫ (B x + C) Q^(-(n + 1/2)) dx`, `n >= 1`, for a quadratic `Q`:
    /// `B/(2a) Q'` integrates directly; `J_n = ∫ Q^(-(n+1/2))` satisfies
    /// `J_n = 2(2a x + b)/((2n-1) Δ Q^(n-1/2)) + 8a(n-1)/((2n-1)Δ) J_(n-1)`
    /// with `Δ = 4ac - b²`, down to `J_1 = 2(2ax + b)/(Δ √Q)`.
    pub(super) fn quadratic_radical_power(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        let factors = if self.cx.graph.op(f) == core::MUL { self.cx.graph.children(f).to_vec() } else { vec![f] };
        let mut constants = Vec::new();
        let mut linear: Option<QPoly> = None;
        let mut radical: Option<(NodeId, i64)> = None;
        for factor in factors {
            if !self.depends(factor) {
                constants.push(factor);
                continue;
            }
            if let (true, &[q, e]) = (self.cx.graph.op(factor) == core::POW, self.cx.graph.children(factor)) {
                if let Some(r) = self.number(e).and_then(|v| v.to_rational()) {
                    if *r.denom() == BigInt::from(2) && r.is_negative() && radical.is_none() {
                        // r = -(n + 1/2)
                        let n = i64::try_from(&(-(r + BigRational::new(BigInt::one(), BigInt::from(2)))).to_integer()).ok()?;
                        if n >= 1 {
                            radical = Some((q, n));
                            continue;
                        }
                    }
                }
            }
            if linear.is_none() {
                let p = self.q_poly_of(factor)?;
                if p.len() <= 2 {
                    linear = Some(p);
                    continue;
                }
            }
            return None;
        }
        let (q_node, n) = radical?;
        let q = self.q_poly_of(q_node)?;
        let [c, b, a] = q.as_slice() else {
            return None;
        };
        let lin = linear.unwrap_or_else(|| vec![BigRational::one()]);
        let big_c = lin.first().cloned().unwrap_or_else(BigRational::zero);
        let big_b = lin.get(1).cloned().unwrap_or_else(BigRational::zero);
        let two = BigRational::from_integer(BigInt::from(2));
        let delta = BigRational::from_integer(BigInt::from(4)) * a * c - b * b;
        if delta.is_zero() {
            return None;
        }
        let mut pieces = Vec::new();
        let derivative_scale = &big_b / (&two * a);
        if !derivative_scale.is_zero() {
            // ∫ Q' Q^(-(n+1/2)) = Q^(1/2 - n)/(1/2 - n)
            let exponent = BigRational::new(BigInt::from(1 - 2 * n), BigInt::from(2));
            let e = self.cx.graph.num(Number::rat(exponent.clone()));
            let raised = self.pow(q_node, e);
            let scale = self.cx.graph.num(Number::rat(&derivative_scale / exponent));
            pieces.push(self.mul(&[scale, raised]));
        }
        let mut weight = &big_c - &derivative_scale * b;
        if !weight.is_zero() {
            let lin_node = {
                let coefficients = [b.clone(), &two * a];
                self.polynomial(&coefficients)
            };
            let mut k = n;
            while k >= 1 {
                let denominator = BigRational::from_integer(BigInt::from(2 * k - 1)) * &delta;
                let exponent = BigRational::new(BigInt::from(1 - 2 * k), BigInt::from(2));
                let e = self.cx.graph.num(Number::rat(exponent));
                let raised = self.pow(q_node, e);
                let scale = self.cx.graph.num(Number::rat(&weight * &two / &denominator));
                pieces.push(self.mul(&[scale, lin_node, raised]));
                weight = weight * BigRational::from_integer(BigInt::from(8 * (k - 1))) * a / &denominator;
                if weight.is_zero() {
                    break;
                }
                k -= 1;
            }
        }
        let body = self.add(&pieces);
        constants.push(body);
        Some(self.mul(&constants))
    }

    /// Rational functions of `sin u`, `cos u`, `tan u`, `u = a x + b`, by
    /// `t = tan(u/2)`.
    pub(super) fn weierstrass(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        let fs = self.f;
        let u = self.single_argument(f, &[fs.sin, fs.cos, fs.tan])?;
        let a = self.linear(u)?;
        let t = {
            let s = self.cx.graph.interner_mut().fresh_symbol("t");
            self.cx.graph.symbol_node(s)
        };
        let one = self.int(1);
        let two = self.int(2);
        let t2 = self.pow(t, two);
        let plus = self.add(&[one, t2]);
        let neg_t2 = self.neg(t2);
        let minus = self.add(&[one, neg_t2]);
        let two_t = self.mul(&[two, t]);
        let s_val = self.div(two_t, plus);
        let c_val = self.div(minus, plus);
        let t_val = self.div(two_t, minus);
        let mut g = f;
        for (op, value) in [(fs.sin, s_val), (fs.cos, c_val), (fs.tan, t_val)] {
            let node = self.call(op, u);
            g = self.cx.graph.replace_subterm(g, node, value);
        }
        if self.depends(g) {
            return None;
        }
        // dx = 2 dt / (a (1 + t²))
        let a_plus = self.mul(&[a, plus]);
        let jacobian = self.div(two, a_plus);
        let integrand = self.mul(&[g, jacobian]);
        let integrand = self.cx.simplify(integrand);
        let primitive = super::antiderivative(self.cx, self.f, integrand, t)?;
        let half = self.frac(1, 2);
        let hu = self.mul(&[half, u]);
        let tan_half = self.call(fs.tan, hu);
        Some(self.cx.graph.substitute(primitive, t, tan_half))
    }

    /// `∫ R(x, √Q) dx` for a quadratic `Q` over `Q`.
    pub(super) fn quadratic_radical(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        let x = self.x;
        let fs = self.f;
        // The radicand: Q^(k/2) or sqrt(Q).
        let mut radicand = None;
        let mut stack = vec![f];
        while let Some(n) = stack.pop() {
            let op = self.cx.graph.op(n);
            let children = self.cx.graph.children(n).to_vec();
            let candidate = match children.as_slice() {
                | &[q] if op == fs.sqrt => Some(q),
                | &[q, e] if op == core::POW && self.number(e).and_then(|v| v.to_rational()).is_some_and(|r| *r.denom() == BigInt::from(2)) => Some(q),
                | _ => None,
            };
            if let Some(q) = candidate.filter(|&q| self.depends(q)) {
                match radicand {
                    | None => radicand = Some(q),
                    | Some(r) if r == q => {},
                    | Some(_) => return None,
                }
            }
            stack.extend(children);
        }
        let q_node = radicand?;
        let q = self.q_poly_of(q_node)?;
        if q.len() != 3 {
            return None;
        }
        // r stands for √Q.
        let r = {
            let s = self.cx.graph.interner_mut().fresh_symbol("r");
            self.cx.graph.symbol_node(s)
        };
        let mut g = f;
        let sqrt_node = self.call(fs.sqrt, q_node);
        g = self.cx.graph.replace_subterm(g, sqrt_node, r);
        for k in [-5_i64, -3, -1, 1, 3, 5] {
            let e = self.frac(k, 2);
            let node = self.pow(q_node, e);
            let k_node = self.int(k);
            let replacement = self.pow(r, k_node);
            g = self.cx.graph.replace_subterm(g, node, replacement);
        }
        let mut gens = Gens::default();
        let gx = gens.index(self.cx.graph, x);
        let gr = gens.index(self.cx.graph, r);
        let fraction = ratio(self.cx.graph, &mut gens, g, Limits::default())?;
        if gens.len() != 2 {
            return None;
        }
        // N = N0 + N1 r with r² = Q.
        let q_poly = Poly::from_univariate(gx, &q.iter().cloned().map(Number::rat).collect::<Vec<_>>());
        let split = |p: &Poly| -> Option<(Poly, Poly)> {
            let mut even = Poly::zero();
            let mut odd = Poly::zero();
            for (k, c) in p.coefficients_in(gr).iter().enumerate() {
                let k = u32::try_from(k).ok()?;
                let factor = q_poly.pow(k / 2, 4096)?;
                let term = c.mul(&factor, 4096)?;
                if k % 2 == 0 {
                    even = even.add(&term);
                } else {
                    odd = odd.add(&term);
                }
            }
            Some((even, odd))
        };
        let (n0, n1) = split(&fraction.numer)?;
        let (d0, d1) = split(&fraction.denom)?;
        let den = d0.mul(&d0, 4096)?.sub(&d1.mul(&d1, 4096)?.mul(&q_poly, 4096)?);
        let a_num = n0.mul(&d0, 4096)?.sub(&n1.mul(&d1, 4096)?.mul(&q_poly, 4096)?);
        let b_num = n1.mul(&d0, 4096)?.sub(&n0.mul(&d1, 4096)?);
        let as_q = |p: &Poly| -> Option<QPoly> {
            let mut v: QPoly = p.univariate_in(gx)?.iter().map(Number::to_rational).collect::<Option<_>>()?;
            while v.last().is_some_and(Zero::is_zero) {
                v.pop();
            }
            Some(v)
        };
        let (a_q, b_q, den_q) = (as_q(&a_num)?, as_q(&b_num)?, as_q(&den)?);
        let mut pieces = Vec::new();
        if !a_q.is_empty() {
            let num = self.q_term(&a_q);
            let d = self.q_term(&den_q);
            let rational = self.div(num, d);
            pieces.push(super::antiderivative(self.cx, self.f, rational, x)?);
        }
        if !b_q.is_empty() {
            // ∫ B √Q = ∫ (B Q) / √Q
            let t_num = univariate::mul(&b_q, &q);
            let parts = crate::rules::poly::apart::apart(&t_num, &den_q)?;
            if !parts.quotient.is_empty() {
                pieces.push(self.polynomial_over_root(&parts.quotient, &q)?);
            }
            for piece in &parts.pieces {
                let ([alpha_den, alpha_num], [c]) = (piece.factor.as_slice(), piece.numerator.as_slice()) else {
                    return None;
                };
                // factor = alpha_num x + alpha_den = alpha_num (x - α)
                let alpha = -(alpha_den / alpha_num);
                let mut scale = c.clone();
                for _ in 0..piece.power {
                    scale /= alpha_num;
                }
                let integral = self.pole_over_root(&alpha, piece.power, &q)?;
                let scale = self.cx.graph.num(Number::rat(scale));
                pieces.push(self.mul(&[scale, integral]));
            }
        }
        Some(self.add(&pieces))
    }

    /// `∫ dx / √Q` for `Q = a x² + b x + c`.
    fn base_root_integral(
        &mut self,
        q: &[BigRational],
    ) -> Option<NodeId> {
        let fs = self.f;
        let x = self.x;
        let (c, b, a) = (&q[0], &q[1], &q[2]);
        let q_node = self.q_term(q);
        let half = self.frac(1, 2);
        let root_q = self.pow(q_node, half);
        let two = BigRational::from_integer(BigInt::from(2));
        let linear = self.q_term(&[b.clone(), &two * a]);
        if a.is_positive() {
            // ln(2 √a √Q + 2 a x + b) / √a
            let a_node = self.cx.graph.num(Number::rat(a.clone()));
            let sqrt_a = self.pow(a_node, half);
            let two_n = self.int(2);
            let t = self.mul(&[two_n, sqrt_a, root_q]);
            let sum = self.add(&[t, linear]);
            let log = self.call(fs.ln, sum);
            Some(self.div(log, sqrt_a))
        } else {
            // -asin((2 a x + b)/√(b² - 4 a c)) / √(-a)
            let disc = b * b - BigRational::from_integer(BigInt::from(4)) * a * c;
            if !disc.is_positive() {
                return None;
            }
            let disc = self.cx.graph.num(Number::rat(disc));
            let sqrt_d = self.pow(disc, half);
            let arg = self.div(linear, sqrt_d);
            let asin = self.call(fs.asin, arg);
            let minus_a = self.cx.graph.num(Number::rat(-a.clone()));
            let sqrt_ma = self.pow(minus_a, half);
            let quotient = self.div(asin, sqrt_ma);
            let _ = x;
            Some(self.neg(quotient))
        }
    }

    /// `∫ P / √Q`: `S √Q + λ ∫ dx/√Q` with `S' Q + S Q'/2 + λ = P`.
    fn polynomial_over_root(
        &mut self,
        p: &[BigRational],
        q: &[BigRational],
    ) -> Option<NodeId> {
        let n = p.len(); // degree n - 1; S has degree n - 2 (n - 1 unknowns) plus λ
        let s_len = n.saturating_sub(1);
        let unknowns = s_len + 1;
        let rows_count = n.max(1) + 1;
        let mut rows = vec![vec![BigRational::zero(); unknowns]; rows_count];
        let dq = univariate::derivative(q);
        let half = BigRational::new(BigInt::one(), BigInt::from(2));
        for i in 0..s_len {
            // S = x^i: (x^i)' Q + x^i Q'/2
            let mut mono = vec![BigRational::zero(); i + 1];
            mono[i] = BigRational::one();
            let d = univariate::derivative(&mono);
            let term1 = univariate::mul(&d, q);
            let term2: QPoly = univariate::mul(&mono, &dq).iter().map(|c| c * &half).collect();
            let total = univariate::add(&term1, &term2);
            for (k, v) in total.iter().enumerate() {
                if k < rows_count {
                    rows[k][i] += v;
                }
            }
        }
        rows[0][s_len] = BigRational::one();
        let mut rhs = vec![BigRational::zero(); rows_count];
        for (k, v) in p.iter().enumerate() {
            rhs[k] = v.clone();
        }
        let sol = solve_q(rows, rhs, unknowns)?;
        let s_poly: QPoly = sol[..s_len].to_vec();
        let lambda = sol[s_len].clone();
        let mut pieces = Vec::new();
        if s_poly.iter().any(|c| !c.is_zero()) {
            let s_term = self.q_term(&s_poly);
            let q_node = self.q_term(q);
            let half_n = self.frac(1, 2);
            let root = self.pow(q_node, half_n);
            pieces.push(self.mul(&[s_term, root]));
        }
        if !lambda.is_zero() {
            let base = self.base_root_integral(q)?;
            let l = self.cx.graph.num(Number::rat(lambda));
            pieces.push(self.mul(&[l, base]));
        }
        Some(self.add(&pieces))
    }

    /// `∫ dx / ((x - α)^m √Q)` by `t = 1/(x - α)`.
    fn pole_over_root(
        &mut self,
        alpha: &BigRational,
        m: u32,
        q: &[BigRational],
    ) -> Option<NodeId> {
        let x = self.x;
        // Q(α + 1/t) t² = Q(α) t² + Q'(α) t + a
        let q_alpha = univariate::eval(q, alpha);
        let dq_alpha = univariate::eval(&univariate::derivative(q), alpha);
        let q_star: QPoly = vec![q[2].clone(), dq_alpha, q_alpha];
        let t = {
            let s = self.cx.graph.interner_mut().fresh_symbol("t");
            self.cx.graph.symbol_node(s)
        };
        // -∫ t^(m-1) / √Q*(t) dt (for t > 0)
        let mut gens = Gens::default();
        let gt = gens.index(self.cx.graph, t);
        let star = Poly::from_univariate(gt, &q_star.iter().cloned().map(Number::rat).collect::<Vec<_>>());
        let star = to_term(self.cx.graph, &gens, &star);
        let half = self.frac(-1, 2);
        let root_inv = self.pow(star, half);
        let e = self.int(i64::from(m) - 1);
        let tp = self.pow(t, e);
        let minus = self.int(-1);
        let integrand = self.mul(&[minus, tp, root_inv]);
        let integrand = self.cx.simplify(integrand);
        let primitive = super::antiderivative(self.cx, self.f, integrand, t)?;
        let alpha_node = self.cx.graph.num(Number::rat(alpha.clone()));
        let shifted = {
            let neg = self.neg(alpha_node);
            self.add(&[x, neg])
        };
        let one = self.int(1);
        let t_of_x = self.div(one, shifted);
        Some(self.cx.graph.substitute(primitive, t, t_of_x))
    }

    /// Risch–Norman for one extension `θ = e^{g}` or `θ = ln x`.
    pub(super) fn risch_norman(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        let x = self.x;
        let fs = self.f;
        // Find θ.
        let mut theta = None;
        let mut stack = vec![f];
        while let Some(n) = stack.pop() {
            let op = self.cx.graph.op(n);
            if (op == fs.exp || op == fs.ln) && self.depends(n) {
                match theta {
                    | None => theta = Some(n),
                    | Some(t) if t == n => {},
                    | Some(_) => return None,
                }
                continue;
            }
            stack.extend_from_slice(self.cx.graph.children(n));
        }
        let theta = theta?;
        let is_exp = self.cx.graph.op(theta) == fs.exp;
        let inner = *self.cx.graph.children(theta).first()?;
        if is_exp {
            self.q_poly_of(inner)?;
        } else if inner != x {
            return None;
        }
        // f = N(x, θ) / D(x)
        let mut gens = Gens::default();
        let gx = gens.index(self.cx.graph, x);
        let gt = gens.index(self.cx.graph, theta);
        let fraction = ratio(self.cx.graph, &mut gens, f, Limits::default())?;
        if gens.len() != 2 || fraction.denom.degree_in(gt) != 0 {
            return None;
        }
        let d_q: QPoly = {
            let mut v: QPoly = fraction.denom.univariate_in(gx)?.iter().map(Number::to_rational).collect::<Option<_>>()?;
            while v.last().is_some_and(Zero::is_zero) {
                v.pop();
            }
            v
        };
        // H = D / squarefree(D): the denominator part that survives.
        let squarefree = univariate::square_free(&d_q).into_iter().fold(vec![BigRational::one()], |acc, (p, _)| univariate::mul(&acc, &p));
        let (h, _) = univariate::divrem(&d_q, &squarefree)?;
        let deg_x = fraction.numer.degree_in(gx) as usize + 2;
        let theta_degrees: Vec<u32> = if is_exp {
            let mut v: Vec<u32> = fraction.numer.coefficients_in(gt).iter().enumerate().filter(|(_, c)| !c.is_zero()).map(|(k, _)| u32::try_from(k).unwrap_or(0)).collect();
            v.sort_unstable();
            v
        } else {
            (0..=fraction.numer.degree_in(gt) + 1).collect()
        };
        let mut unknowns = Vec::new();
        let mut terms = Vec::new();
        for &j in &theta_degrees {
            for i in 0..=(deg_x + h.len()) {
                let s = self.cx.graph.interner_mut().fresh_symbol("c");
                let c = self.cx.graph.symbol_node(s);
                unknowns.push(c);
                let ei = self.int(i64::try_from(i).ok()?);
                let ej = self.int(i64::from(j));
                let xi = self.pow(x, ei);
                let tj = self.pow(theta, ej);
                terms.push(self.mul(&[c, xi, tj]));
            }
        }
        let numerator = self.add(&terms);
        let h_term = self.q_term(&h);
        let mut candidate = self.div(numerator, h_term);
        // Logarithms of the squarefree factors of D.
        let (_, factors) = univariate::factor(&squarefree);
        for (factor, _) in factors {
            let fq: QPoly = factor.into_iter().map(BigRational::from_integer).collect();
            if fq.len() < 2 {
                continue;
            }
            let s = self.cx.graph.interner_mut().fresh_symbol("c");
            let c = self.cx.graph.symbol_node(s);
            unknowns.push(c);
            let p = self.q_term(&fq);
            let log = self.call(fs.ln, p);
            let term = self.mul(&[c, log]);
            candidate = self.add(&[candidate, term]);
        }
        // F' - κ f = 0 identically, κ = 1 afterwards.
        let kappa = {
            let s = self.cx.graph.interner_mut().fresh_symbol("k");
            self.cx.graph.symbol_node(s)
        };
        unknowns.push(kappa);
        let derivative = self.derivative(candidate);
        let scaled = self.mul(&[kappa, f]);
        let negated = self.neg(scaled);
        let residual = self.add(&[derivative, negated]);
        let basis = crate::rules::ode::solve_linear_identity(self.cx.graph, residual, &unknowns)?;
        let solution = basis.into_iter().find(|v| v.last().is_some_and(|k| !k.is_zero()))?;
        let k = solution.last()?.clone();
        let mut result = candidate;
        for (u, value) in unknowns.iter().zip(&solution) {
            let v = self.cx.graph.num(Number::rat(value / &k));
            result = self.cx.graph.substitute(result, *u, v);
        }
        Some(self.cx.simplify(result))
    }

    /// `P / (x⁴ + c x² + e)` with an irreducible biquadratic denominator.
    pub(super) fn biquadratic(
        &mut self,
        f: NodeId,
    ) -> Option<NodeId> {
        let x = self.x;
        let mut gens = Gens::default();
        let gx = gens.index(self.cx.graph, x);
        let fraction = ratio(self.cx.graph, &mut gens, f, Limits::default())?;
        if gens.len() != 1 {
            return None;
        }
        let as_q = |p: &Poly| -> Option<QPoly> {
            let mut v: QPoly = p.univariate_in(gx)?.iter().map(Number::to_rational).collect::<Option<_>>()?;
            while v.last().is_some_and(Zero::is_zero) {
                v.pop();
            }
            Some(v)
        };
        let (numer, denom) = (as_q(&fraction.numer)?, as_q(&fraction.denom)?);
        if denom.len() != 5 || !denom[1].is_zero() || !denom[3].is_zero() || numer.len() > 4 {
            return None;
        }
        let lead = denom[4].clone();
        let (c, e) = (&denom[2] / &lead, &denom[0] / &lead);
        if !e.is_positive() {
            return None;
        }
        let half = self.frac(1, 2);
        let e_node = self.cx.graph.num(Number::rat(e));
        let t = self.pow(e_node, half);
        let t = self.cx.simplify(t);
        let two = self.int(2);
        let two_t = self.mul(&[two, t]);
        let c_node = self.cx.graph.num(Number::rat(c));
        let neg_c = self.neg(c_node);
        let s2 = self.add(&[two_t, neg_c]);
        let s2 = self.cx.simplify(s2);
        if self.cx.graph.eval(s2, &Env::numeric(0.0)).is_none_or(|v| v <= 0.0) {
            return None;
        }
        let s = self.pow(s2, half);
        let s = self.cx.simplify(s);
        let n: Vec<NodeId> = (0..4)
            .map(|k| {
                let v = numer.get(k).cloned().unwrap_or_else(BigRational::zero) / &lead;
                self.cx.graph.num(Number::rat(v))
            })
            .collect();
        // A = (n3 - (n2 - n0/t)/s)/2, C = (n3 + (n2 - n0/t)/s)/2,
        // B = (n0/t - (n1 - t n3)/s)/2, D = (n0/t + (n1 - t n3)/s)/2
        let n0_t = self.div(n[0], t);
        let neg_n0_t = self.neg(n0_t);
        let u = self.add(&[n[2], neg_n0_t]);
        let u_s = self.div(u, s);
        let t_n3 = self.mul(&[t, n[3]]);
        let neg_tn3 = self.neg(t_n3);
        let v = self.add(&[n[1], neg_tn3]);
        let v_s = self.div(v, s);
        let neg_us = self.neg(u_s);
        let neg_vs = self.neg(v_s);
        let a_c = self.add(&[n[3], neg_us]);
        let c_c = self.add(&[n[3], u_s]);
        let b_c = self.add(&[n0_t, neg_vs]);
        let d_c = self.add(&[n0_t, v_s]);
        let (a_c, b_c, c_c, d_c) = (self.mul(&[half, a_c]), self.mul(&[half, b_c]), self.mul(&[half, c_c]), self.mul(&[half, d_c]));
        let x2 = self.pow(x, two);
        let sx = self.mul(&[s, x]);
        let neg_sx = self.neg(sx);
        let q1 = self.add(&[x2, sx, t]);
        let q2 = self.add(&[x2, neg_sx, t]);
        let ax = self.mul(&[a_c, x]);
        let cx_ = self.mul(&[c_c, x]);
        let top2 = self.add(&[ax, b_c]);
        let top1 = self.add(&[cx_, d_c]);
        // (A x + B) q₂ + (C x + D) q₁ = N: (A x + B) sits over q₁.
        let f2 = self.div(top2, q1);
        let f1 = self.div(top1, q2);
        let f1 = self.cx.simplify(f1);
        let f2 = self.cx.simplify(f2);
        let i1 = self.rational_parametric(f1)?;
        let i2 = self.rational_parametric(f2)?;
        Some(self.add(&[i1, i2]))
    }
}

/// `∫_{-∞}^{∞} R(x) w(x) dx` with `w ∈ {1, cos(a x), sin(a x)}` by residues
/// in the upper half-plane, and Dirichlet's integral.
pub(super) fn by_residues(
    cx: &mut Cx<'_>,
    integrand: NodeId,
    x: NodeId,
    lower: NodeId,
    upper: NodeId,
) -> Option<NodeId> {
    let value = |cx: &Cx<'_>, n: NodeId| cx.graph.eval(n, &Env::numeric(0.0));
    let (lo, hi) = (value(cx, lower)?, value(cx, upper)?);
    let (sin, cos, pi, i_op, exp) = (
        cx.graph.ops().lookup("sin")?,
        cx.graph.ops().lookup("cos")?,
        cx.graph.ops().lookup("pi")?,
        cx.graph.ops().lookup("I")?,
        cx.graph.ops().lookup("exp")?,
    );
    let symbol = cx.graph.symbol_of(x)?;
    let pi = cx.graph.node(pi, &[]);
    let i_unit = cx.graph.node(i_op, &[]);
    // Split off one sin/cos(a x) factor.
    let factors = if cx.graph.op(integrand) == core::MUL { cx.graph.children(integrand).to_vec() } else { vec![integrand] };
    let mut wave = None;
    let mut rest = Vec::new();
    for fac in factors {
        let op = cx.graph.op(fac);
        if (op == sin || op == cos) && wave.is_none() {
            let &[u] = cx.graph.children(fac) else { return None };
            wave = Some((op == cos, u));
        } else {
            rest.push(fac);
        }
    }
    let r = match rest.as_slice() {
        | [] => cx.graph.int(1),
        | [one] => *one,
        | _ => cx.graph.node(core::MUL, &rest),
    };
    // a from u = a x.
    let a = match wave {
        | Some((_, u)) => {
            let parts = crate::rules::ode::coefficients_of(cx.graph, u, x)?;
            match parts.as_slice() {
                | [c0, a] if cx.is_zero(*c0) => Some(*a),
                | _ => return None,
            }
        },
        | None => None,
    };
    let (numer, denom) = crate::rules::poly::rational_function_in(cx.graph, r, x)?;
    // Dirichlet: sin(a x)/x on (0, ∞) or (-∞, ∞).
    if let (Some((false, _)), Some(a)) = (wave, a) {
        if numer.len() == 1 && denom.len() == 2 && denom[0].is_zero() && hi.is_infinite() && hi > 0.0 {
            let c = &numer[0] / &denom[1];
            let sign = cx.graph.ops().lookup("sign")?;
            let s = cx.graph.node(sign, &[a]);
            let c = cx.graph.num(Number::rat(c));
            let factor = if lo == 0.0 {
                let half = cx.graph.num(Number::fraction(1, 2)?);
                cx.graph.node(core::MUL, &[half, pi])
            } else if lo.is_infinite() && lo < 0.0 {
                pi
            } else {
                return None;
            };
            return Some(cx.graph.node(core::MUL, &[c, factor, s]));
        }
    }
    if !(lo.is_infinite() && lo < 0.0 && hi.is_infinite() && hi > 0.0) {
        return None;
    }
    // Convergence: deg D ≥ deg N + 2, or + 1 with an oscillating factor.
    let need = if wave.is_some() { 1 } else { 2 };
    if denom.len() < numer.len() + need {
        return None;
    }
    // Poles: simple roots of quadratic factors over Q, none real.
    let (_, factors) = univariate::factor(&denom);
    let mut residues = Vec::new();
    let a_value = match a {
        | Some(a) => {
            let v = value(cx, a)?;
            if v == 0.0 {
                return None;
            }
            Some(v)
        },
        | None => None,
    };
    let flip = a_value.is_some_and(|v| v < 0.0);
    // Real quadratic factors (a, b, c) as terms: quadratics over Q, and
    // biquadratics x⁴ + c x² + e split as (x² ± s x + t).
    let mut quadratics: Vec<[NodeId; 3]> = Vec::new();
    for (factor, multiplicity) in factors {
        if multiplicity != 1 {
            return None;
        }
        match factor.as_slice() {
            | [c0, b0, a0] => {
                if !(b0 * b0 - BigInt::from(4) * a0 * c0).is_negative() {
                    return None;
                }
                let terms = [a0, b0, c0].map(|v| cx.graph.num(Number::Int(v.clone())));
                quadratics.push(terms);
            },
            | [e0, z1, c2, z3, lead] if z1.is_zero() && z3.is_zero() => {
                let (c, e) = (BigRational::new(c2.clone(), lead.clone()), BigRational::new(e0.clone(), lead.clone()));
                if !e.is_positive() {
                    return None;
                }
                let half = cx.graph.num(Number::fraction(1, 2)?);
                let e_node = cx.graph.num(Number::rat(e));
                let t = cx.graph.node(core::POW, &[e_node, half]);
                let t = cx.simplify(t);
                let two = cx.graph.int(2);
                let two_t = cx.graph.node(core::MUL, &[two, t]);
                let minus_c = cx.graph.num(Number::rat(-c));
                let s2 = cx.graph.node(core::ADD, &[two_t, minus_c]);
                let s2 = cx.simplify(s2);
                if value(cx, s2).is_none_or(|v| v <= 0.0) {
                    return None;
                }
                let s_node = cx.graph.node(core::POW, &[s2, half]);
                let s_node = cx.simplify(s_node);
                let one = cx.graph.int(1);
                let minus_one = cx.graph.int(-1);
                let neg_s = cx.graph.node(core::MUL, &[minus_one, s_node]);
                quadratics.push([one, s_node, t]);
                quadratics.push([one, neg_s, t]);
            },
            | _ => return None,
        }
    }
    let n_term = poly_term(cx, &numer, x);
    let d_term = poly_term(cx, &univariate::derivative(&denom), x);
    for [a0, b0, c0] in quadratics {
        // α = (-b + i √(4ac - b²))/(2a), in the half-plane where e^{i a x}
        // decays.
        let two = cx.graph.int(2);
        let four = cx.graph.int(4);
        let minus_one = cx.graph.int(-1);
        let two_a = cx.graph.node(core::MUL, &[two, a0]);
        let inv_two_a = cx.graph.node(core::POW, &[two_a, minus_one]);
        let neg_b = cx.graph.node(core::MUL, &[minus_one, b0]);
        let re = cx.graph.node(core::MUL, &[neg_b, inv_two_a]);
        let four_ac = cx.graph.node(core::MUL, &[four, a0, c0]);
        let b2 = cx.graph.node(core::POW, &[b0, two]);
        let neg_b2 = cx.graph.node(core::MUL, &[minus_one, b2]);
        let disc = cx.graph.node(core::ADD, &[four_ac, neg_b2]);
        let half = cx.graph.num(Number::fraction(1, 2)?);
        let root = cx.graph.node(core::POW, &[disc, half]);
        let im = cx.graph.node(core::MUL, &[root, inv_two_a]);
        let sign = cx.graph.int(if flip { -1 } else { 1 });
        let i_im = cx.graph.node(core::MUL, &[sign, i_unit, im]);
        let alpha = cx.graph.node(core::ADD, &[re, i_im]);
        let alpha = cx.simplify(alpha);
        // Res = N(α) e^{i a α} / D'(α)
        let n_at = cx.graph.substitute(n_term, x, alpha);
        let d_at = cx.graph.substitute(d_term, x, alpha);
        let inv = cx.graph.node(core::POW, &[d_at, minus_one]);
        let mut res = cx.graph.node(core::MUL, &[n_at, inv]);
        if let Some(a) = a {
            let phase = cx.graph.node(core::MUL, &[i_unit, a, alpha]);
            let e = cx.graph.node(exp, &[phase]);
            res = cx.graph.node(core::MUL, &[res, e]);
        }
        residues.push(res);
    }
    let sum = cx.graph.node(core::ADD, &residues);
    let two = cx.graph.int(2);
    let sign = cx.graph.int(if flip { -1 } else { 1 });
    let total = cx.graph.node(core::MUL, &[sign, two, pi, i_unit, sum]);
    let (re, im) = (cx.graph.ops().lookup("re")?, cx.graph.ops().lookup("im")?);
    let result = match wave {
        | Some((true, _)) | None => cx.graph.node(re, &[total]),
        | Some((false, _)) => cx.graph.node(im, &[total]),
    };
    let result = cx.simplify(result);
    let _ = symbol;
    // Numeric sanity: the result must be a finite real number when it can
    // be evaluated.
    match cx.graph.eval(result, &Env::numeric(0.0)) {
        | Some(v) if !v.is_finite() => None,
        | _ => Some(result),
    }
}

fn poly_term(
    cx: &mut Cx<'_>,
    p: &[BigRational],
    x: NodeId,
) -> NodeId {
    let mut gens = Gens::default();
    let gx = gens.index(cx.graph, x);
    let numbers: Vec<Number> = p.iter().cloned().map(Number::rat).collect();
    to_term(cx.graph, &gens, &Poly::from_univariate(gx, &numbers))
}
