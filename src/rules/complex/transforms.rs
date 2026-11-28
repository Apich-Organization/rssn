//! Symbolic complex integral transforms (Mellin, inverse Mellin, Bromwich Laplace inversion,
//! Meijer G, Mellin–Barnes reduction, Hilbert transform) and Risch–Bronstein complex integration.

use super::Ops;
use super::build;
use super::eval_complex;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::Tier;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::rules::calculus::derivative;
use crate::rules::poly::best;
use num_complex::Complex64;
use std::collections::HashMap;

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let mellin = i.op(OpDescriptor::new("mellin_transform", Arity::Fixed(3))
        .flags(OpFlags::HEAVY)
        .cost(70)
        .doc(
            "`mellin_transform(f, x, s)`: symbolic Mellin transform `int_0^oo x^(s-1) f(x) dx`.",
        ))?;
    let inv_mellin = i.op(
        OpDescriptor::new("inverse_mellin_transform", Arity::Fixed(3))
            .flags(OpFlags::HEAVY)
            .cost(70)
            .doc("`inverse_mellin_transform(F, s, x)`: symbolic inverse Mellin transform."),
    )?;
    let bromwich = i.op(
        OpDescriptor::new("bromwich_inverse_laplace", Arity::Fixed(3))
            .flags(OpFlags::HEAVY)
            .cost(80)
            .doc("`bromwich_inverse_laplace(F, s, t)`: Bromwich inverse Laplace transform via sum of residues of `F(s)*exp(s*t)`."),
    )?;
    let meijer = i.op(
        OpDescriptor::new("meijer_g", Arity::Fixed(5))
            .flags(OpFlags::HEAVY)
            .cost(75)
            .doc("`meijer_g(a_list, b_list, n, m, z)`: Meijer G-function `G_{p,q}^{m,n}(z | a; b)` with elementary reductions."),
    )?;
    let mb_reduce = i.op(
        OpDescriptor::new("mellin_barnes_to_meijer_g", Arity::Fixed(5))
            .flags(OpFlags::HEAVY)
            .cost(70)
            .doc("`mellin_barnes_to_meijer_g(num_plus, num_minus, den_plus, den_minus, z)`: reduces a Mellin–Barnes Gamma-product integral to `meijer_g`."),
    )?;
    let hilbert = i.op(
        OpDescriptor::new("hilbert_transform", Arity::Fixed(3))
            .flags(OpFlags::HEAVY)
            .cost(70)
            .doc("`hilbert_transform(f, x, y)`: symbolic Hilbert transform `(1/pi) PV int_{-oo}^{oo} f(x)/(y - x) dx`."),
    )?;
    let risch_c = i.op(
        OpDescriptor::new("complex_risch_integrate", Arity::Fixed(2))
            .flags(OpFlags::HEAVY)
            .cost(90)
            .doc("`complex_risch_integrate(f, z)`: Risch–Bronstein complex antiderivative via Euler complexification and complex logarithmic residues."),
    )?;
    let ops = Ops::of(i.graph()).ok_or_else(|| {
        RuleError::Invalid {
            rule: "complex/transforms".into(),
            reason: "needs elementary",
        }
    })?;
    i.kernel(
        "complex/transforms",
        Tier::Reduce,
        ComplexTransforms {
            mellin,
            inv_mellin,
            bromwich,
            meijer,
            mb_reduce,
            hilbert,
            risch_c,
            ops,
        },
    );
    Ok(())
}

struct ComplexTransforms {
    mellin: OpId,
    inv_mellin: OpId,
    bromwich: OpId,
    meijer: OpId,
    mb_reduce: OpId,
    hilbert: OpId,
    risch_c: OpId,
    ops: Ops,
}

impl Kernel for ComplexTransforms {
    fn ops(&self) -> Vec<OpId> {
        vec![
            self.mellin,
            self.inv_mellin,
            self.bromwich,
            self.meijer,
            self.mb_reduce,
            self.hilbert,
            self.risch_c,
        ]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        self.reduce_opt(cx, node)
            .map_or(Outcome::Pass, Outcome::Equal)
    }
}

impl ComplexTransforms {
    fn reduce_opt(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Option<NodeId> {
        let op = cx.graph.op(node);
        let ch = cx.graph.children(node).to_vec();
        if op == self.mellin {
            let (f, x, s) = (*ch.first()?, *ch.get(1)?, *ch.get(2)?);
            self.mellin_reduce(cx, f, x, s)
        } else if op == self.inv_mellin {
            let (big_f, s, x) = (*ch.first()?, *ch.get(1)?, *ch.get(2)?);
            self.inv_mellin_reduce(cx, big_f, s, x)
        } else if op == self.bromwich {
            let (big_f, s, t) = (*ch.first()?, *ch.get(1)?, *ch.get(2)?);
            self.bromwich_reduce(cx, big_f, s, t)
        } else if op == self.meijer {
            let (a_list, b_list, n_node, m_node, z) = (
                *ch.first()?,
                *ch.get(1)?,
                *ch.get(2)?,
                *ch.get(3)?,
                *ch.get(4)?,
            );
            self.meijer_reduce(cx, a_list, b_list, n_node, m_node, z)
        } else if op == self.mb_reduce {
            let (np, nm, dp, dm, z) = (
                *ch.first()?,
                *ch.get(1)?,
                *ch.get(2)?,
                *ch.get(3)?,
                *ch.get(4)?,
            );
            self.mb_to_meijer(cx, np, nm, dp, dm, z)
        } else if op == self.hilbert {
            let (f, x, y) = (*ch.first()?, *ch.get(1)?, *ch.get(2)?);
            self.hilbert_reduce(cx, f, x, y)
        } else if op == self.risch_c {
            let (f, z) = (*ch.first()?, *ch.get(1)?);
            self.complex_risch(cx, f, z)
        } else {
            None
        }
    }

    fn mellin_reduce(
        &self,
        cx: &mut Cx<'_>,
        f: NodeId,
        x: NodeId,
        s: NodeId,
    ) -> Option<NodeId> {
        let gamma_op = cx.graph.ops().lookup("gamma")?;
        let f_best = best(cx.graph, f)?;
        let neg_x = build::neg(cx.graph, x);
        let exp_neg_x = build::call(cx.graph, self.ops.exp, &[neg_x]);
        if cx.graph.find(f_best) == cx.graph.find(exp_neg_x) {
            let g = build::call(cx.graph, gamma_op, &[s]);
            return Some(cx.simplify(g));
        }
        // Check 1 / (1 + x) -> pi / sin(pi * s)
        let one = cx.graph.int(1);
        let one_plus_x = build::add(cx.graph, &[one, x]);
        let neg_one = cx.graph.int(-1);
        let inv_1px = build::pow(cx.graph, one_plus_x, neg_one);
        let inv_1px_s = cx.simplify(inv_1px);
        if cx.graph.find(f_best) == cx.graph.find(inv_1px_s) {
            let pi = build::call(cx.graph, self.ops.pi, &[]);
            let pi_s = build::mul(cx.graph, &[pi, s]);
            let sin_pi_s = build::call(cx.graph, self.ops.sin, &[pi_s]);
            let inv_sin = build::pow(cx.graph, sin_pi_s, neg_one);
            let res = build::mul(cx.graph, &[pi, inv_sin]);
            return Some(cx.simplify(res));
        }
        // Check exp(-a * x) where a is independent of x -> gamma(s) * a^(-s)
        let f_op = cx.graph.op(f_best);
        if f_op == self.ops.exp {
            let arg = *cx.graph.children(f_best).first()?;
            let x_sym = cx.graph.symbol_of(x)?;
            let neg_arg = build::neg(cx.graph, arg);
            let inv_x = build::pow(cx.graph, x, neg_one);
            let a_raw = build::mul(cx.graph, &[neg_arg, inv_x]);
            let a = cx.simplify(a_raw);
            if !cx.graph.depends_on(cx.graph.find(a), x_sym) {
                let g = build::call(cx.graph, gamma_op, &[s]);
                let neg_s = build::neg(cx.graph, s);
                let a_pow = build::pow(cx.graph, a, neg_s);
                let res = build::mul(cx.graph, &[g, a_pow]);
                return Some(cx.simplify(res));
            }
        }
        None
    }

    fn inv_mellin_reduce(
        &self,
        cx: &mut Cx<'_>,
        big_f: NodeId,
        s: NodeId,
        x: NodeId,
    ) -> Option<NodeId> {
        let gamma_op = cx.graph.ops().lookup("gamma")?;
        let f_best = best(cx.graph, big_f)?;
        let gamma_s = build::call(cx.graph, gamma_op, &[s]);
        if cx.graph.find(f_best) == cx.graph.find(gamma_s) {
            let neg_x = build::neg(cx.graph, x);
            let res = build::call(cx.graph, self.ops.exp, &[neg_x]);
            return Some(cx.simplify(res));
        }
        let pi = build::call(cx.graph, self.ops.pi, &[]);
        let pi_s = build::mul(cx.graph, &[pi, s]);
        let sin_pi_s = build::call(cx.graph, self.ops.sin, &[pi_s]);
        let neg_one = cx.graph.int(-1);
        let inv_sin = build::pow(cx.graph, sin_pi_s, neg_one);
        let pi_csc_raw = build::mul(cx.graph, &[pi, inv_sin]);
        let pi_csc = cx.simplify(pi_csc_raw);
        if cx.graph.find(f_best) == cx.graph.find(pi_csc) {
            let one = cx.graph.int(1);
            let one_p_x = build::add(cx.graph, &[one, x]);
            let res = build::pow(cx.graph, one_p_x, neg_one);
            return Some(cx.simplify(res));
        }
        None
    }

    fn bromwich_reduce(
        &self,
        cx: &mut Cx<'_>,
        big_f: NodeId,
        s: NodeId,
        t: NodeId,
    ) -> Option<NodeId> {
        // Check for 1/sqrt(s) -> 1/sqrt(pi * t)
        let half = cx.graph.num(Number::fraction(1, 2)?);
        let neg_half = cx.graph.num(Number::fraction(-1, 2)?);
        let s_mhalf = build::pow(cx.graph, s, neg_half);
        if cx.graph.find(big_f) == cx.graph.find(s_mhalf) {
            let pi = build::call(cx.graph, self.ops.pi, &[]);
            let pi_t = build::mul(cx.graph, &[pi, t]);
            let inv_sqrt = build::pow(cx.graph, pi_t, neg_half);
            return Some(cx.simplify(inv_sqrt));
        }
        let _ = half;
        let poles_op = cx.graph.ops().lookup("poles")?;
        let residue_op = cx.graph.ops().lookup("residue")?;
        let poles_req = build::call(cx.graph, poles_op, &[big_f, s]);
        let poles_list = cx.simplify(poles_req);
        if cx.graph.op(poles_list) != core::LIST {
            return None;
        }
        let entries = cx.graph.children(poles_list).to_vec();
        let st = build::mul(cx.graph, &[s, t]);
        let exp_st = build::call(cx.graph, self.ops.exp, &[st]);
        let integrand = build::mul(cx.graph, &[big_f, exp_st]);
        let mut terms = Vec::with_capacity(entries.len());
        for entry in entries {
            if cx.graph.op(entry) != core::LIST {
                return None;
            }
            let p_node = *cx.graph.children(entry).first()?;
            let res_req = build::call(cx.graph, residue_op, &[integrand, s, p_node]);
            let res_val = cx.simplify(res_req);
            if cx.graph.op(res_val) == residue_op {
                return None;
            }
            terms.push(res_val);
        }
        let sum = build::add(cx.graph, &terms);
        Some(cx.simplify(sum))
    }

    fn meijer_reduce(
        &self,
        cx: &mut Cx<'_>,
        a_list: NodeId,
        b_list: NodeId,
        n_node: NodeId,
        m_node: NodeId,
        z: NodeId,
    ) -> Option<NodeId> {
        if cx.graph.op(a_list) != core::LIST || cx.graph.op(b_list) != core::LIST {
            return None;
        }
        let a_ch = cx.graph.children(a_list).to_vec();
        let b_ch = cx.graph.children(b_list).to_vec();
        let n = cx.graph.number_of(n_node).and_then(Number::to_i64)?;
        let m = cx.graph.number_of(m_node).and_then(Number::to_i64)?;

        // G_{0,1}^{1,0}(z | b0) = z^{b0} * exp(-z)
        if a_ch.is_empty() && b_ch.len() == 1 && n == 0 && m == 1 {
            let b0 = *b_ch.first()?;
            let neg_z = build::neg(cx.graph, z);
            let exp_nz = build::call(cx.graph, self.ops.exp, &[neg_z]);
            let res = if cx.graph.number_of(b0).is_some_and(Number::is_zero) {
                exp_nz
            } else {
                let z_b0 = build::pow(cx.graph, z, b0);
                build::mul(cx.graph, &[z_b0, exp_nz])
            };
            return Some(cx.simplify(res));
        }

        // G_{1,1}^{1,1}(z | 0; 0) = 1 / (1 + z)
        if a_ch.len() == 1 && b_ch.len() == 1 && n == 1 && m == 1 {
            let a0 = *a_ch.first()?;
            let b0 = *b_ch.first()?;
            if cx.graph.number_of(a0).is_some_and(Number::is_zero)
                && cx.graph.number_of(b0).is_some_and(Number::is_zero)
            {
                let one = cx.graph.int(1);
                let one_pz = build::add(cx.graph, &[one, z]);
                let neg_one = cx.graph.int(-1);
                let res = build::pow(cx.graph, one_pz, neg_one);
                return Some(cx.simplify(res));
            }
        }
        None
    }

    fn mb_to_meijer(
        &self,
        cx: &mut Cx<'_>,
        np: NodeId,
        nm: NodeId,
        dp: NodeId,
        dm: NodeId,
        z: NodeId,
    ) -> Option<NodeId> {
        if cx.graph.op(np) != core::LIST
            || cx.graph.op(nm) != core::LIST
            || cx.graph.op(dp) != core::LIST
            || cx.graph.op(dm) != core::LIST
        {
            return None;
        }
        let np_ch = cx.graph.children(np).to_vec();
        let nm_ch = cx.graph.children(nm).to_vec();
        let dp_ch = cx.graph.children(dp).to_vec();
        let dm_ch = cx.graph.children(dm).to_vec();
        let m = cx.graph.int(np_ch.len() as i64);
        let n = cx.graph.int(nm_ch.len() as i64);
        let one = cx.graph.int(1);

        let mut a_vec = Vec::with_capacity(nm_ch.len() + dm_ch.len());
        for c in nm_ch {
            let sub = build::sub(cx.graph, one, c);
            a_vec.push(cx.simplify(sub));
        }
        a_vec.extend(dm_ch);

        let mut b_vec = np_ch;
        for c in dp_ch {
            let sub = build::sub(cx.graph, one, c);
            b_vec.push(cx.simplify(sub));
        }

        let a_list = cx.graph.node(core::LIST, &a_vec);
        let b_list = cx.graph.node(core::LIST, &b_vec);
        let req = build::call(cx.graph, self.meijer, &[a_list, b_list, n, m, z]);
        Some(cx.simplify(req))
    }

    fn hilbert_reduce(
        &self,
        cx: &mut Cx<'_>,
        f: NodeId,
        x: NodeId,
        y: NodeId,
    ) -> Option<NodeId> {
        let f_best = best(cx.graph, f)?;
        let cos_x = build::call(cx.graph, self.ops.cos, &[x]);
        if cx.graph.find(f_best) == cx.graph.find(cos_x) {
            return Some(build::call(cx.graph, self.ops.sin, &[y]));
        }
        let sin_x = build::call(cx.graph, self.ops.sin, &[x]);
        if cx.graph.find(f_best) == cx.graph.find(sin_x) {
            let cos_y = build::call(cx.graph, self.ops.cos, &[y]);
            return Some(build::neg(cx.graph, cos_y));
        }
        // 1 / (x^2 + 1) -> y / (y^2 + 1)
        let one = cx.graph.int(1);
        let neg_one = cx.graph.int(-1);
        let x2 = build::powi(cx.graph, x, 2);
        let x2_p1 = build::add(cx.graph, &[x2, one]);
        let lorentz_raw = build::pow(cx.graph, x2_p1, neg_one);
        let lorentz = cx.simplify(lorentz_raw);
        if cx.graph.find(f_best) == cx.graph.find(lorentz) {
            let y2 = build::powi(cx.graph, y, 2);
            let y2_p1 = build::add(cx.graph, &[y2, one]);
            let inv_y2_p1 = build::pow(cx.graph, y2_p1, neg_one);
            let res = build::mul(cx.graph, &[y, inv_y2_p1]);
            return Some(cx.simplify(res));
        }
        None
    }

    fn complex_risch(
        &self,
        cx: &mut Cx<'_>,
        f: NodeId,
        z: NodeId,
    ) -> Option<NodeId> {
        let integral_op = cx.graph.ops().lookup("integral")?;
        let poles_op = cx.graph.ops().lookup("poles")?;
        let residue_op = cx.graph.ops().lookup("residue")?;
        let z_sym = cx.graph.symbol_of(z)?;

        // 1. Try complex partial-fraction + Laurent integration for rational functions
        //    (handles both simple and repeated poles, plus polynomial quotient part).
        let poles_req = build::call(cx.graph, poles_op, &[f, z]);
        let poles_list = cx.simplify(poles_req);
        if cx.graph.op(poles_list) == core::LIST {
            let entries = cx.graph.children(poles_list).to_vec();
            if !entries.is_empty() {
                let mut int_terms = Vec::new();
                let mut ok = true;
                for entry in entries {
                    if cx.graph.op(entry) != core::LIST {
                        ok = false;
                        break;
                    }
                    let pair = cx.graph.children(entry).to_vec();
                    let p_node = *pair.first()?;
                    let mult = pair
                        .get(1)
                        .and_then(|&m| cx.graph.number_of(m))
                        .and_then(Number::to_i64)
                        .unwrap_or(1)
                        .max(1);
                    let z_minus_p = build::sub(cx.graph, z, p_node);
                    // Residue (c_{-1}) -> c_{-1} * ln(z - p)
                    let res_req = build::call(cx.graph, residue_op, &[f, z, p_node]);
                    let r_k = cx.simplify(res_req);
                    if cx.graph.op(r_k) == residue_op {
                        ok = false;
                        break;
                    }
                    if !cx.graph.number_of(r_k).is_some_and(Number::is_zero) {
                        let ln_term = build::call(cx.graph, self.ops.ln, &[z_minus_p]);
                        int_terms.push(build::mul(cx.graph, &[r_k, ln_term]));
                    }
                    // Higher-order principal terms (j = 2..=mult): c_{-j} / (z - p)^j -> c_{-j} / ((1 - j)*(z - p)^{j-1})
                    if mult > 1 {
                        if let Some((start_pow, coeffs)) =
                            crate::rules::calculus::laurent_expansion(cx, f, z, p_node, mult + 1)
                        {
                            for (idx, c_val) in coeffs.into_iter().enumerate() {
                                let pow = start_pow + (idx as i64);
                                if pow <= -2
                                    && !cx.graph.number_of(c_val).is_some_and(Number::is_zero)
                                {
                                    let j = -pow;
                                    let denom_int = cx.graph.int(1 - j);
                                    let zp = build::powi(cx.graph, z_minus_p, 1 - j);
                                    let inv_d = build::powi(cx.graph, denom_int, -1);
                                    int_terms.push(build::mul(cx.graph, &[c_val, inv_d, zp]));
                                }
                            }
                        } else {
                            ok = false;
                            break;
                        }
                    }
                }
                if ok {
                    let cand_raw = build::add(cx.graph, &int_terms);
                    let candidate = cx.simplify(cand_raw);
                    if let Some(d_cand) = derivative(cx.graph, candidate, z) {
                        let diff_raw = build::sub(cx.graph, d_cand, f);
                        let diff = cx.simplify(diff_raw);
                        if cx.graph.number_of(diff).is_some_and(Number::is_zero) {
                            return Some(candidate);
                        }
                        // Check at multiple non-singular points for complex-valued residues
                        let pts = [
                            Complex64::new(0.37, 0.51),
                            Complex64::new(-0.73, 0.29),
                            Complex64::new(1.13, -0.41),
                        ];
                        if pts.iter().all(|&pt| {
                            let m = HashMap::from([(z_sym, pt)]);
                            eval_complex(cx.graph, diff, &m).is_some_and(|dv| dv.norm() < 1e-10)
                        }) {
                            return Some(candidate);
                        }
                    }
                }
            }
        }

        // 2. Euler-complexify trig/hyperbolic and delegate to Risch / integral
        let f_euler = complexify_trig(cx.graph, f, &self.ops);
        let int_req = build::call(cx.graph, integral_op, &[f_euler, z]);
        let ans = cx.simplify(int_req);
        if cx.graph.op(ans) != integral_op {
            return Some(ans);
        }
        None
    }
}

fn complexify_trig(
    graph: &mut Graph,
    node: NodeId,
    ops: &Ops,
) -> NodeId {
    let op = graph.op(node);
    let ch = graph.children(node).to_vec();
    if ch.is_empty() {
        return node;
    }
    let new_ch: Vec<NodeId> = ch.iter().map(|&c| complexify_trig(graph, c, ops)).collect();
    if new_ch.len() == 1
        && let Some(&u) = new_ch.first()
        && (op == ops.cos || op == ops.sin || op == ops.cosh || op == ops.sinh)
    {
        let i_unit = build::call(graph, ops.unit, &[]);
        let half = graph.num(Number::fraction(1, 2).unwrap_or_else(|| Number::from(0)));
        if op == ops.cos {
            let iu = build::mul(graph, &[i_unit, u]);
            let neg_iu = build::neg(graph, iu);
            let ep = build::call(graph, ops.exp, &[iu]);
            let em = build::call(graph, ops.exp, &[neg_iu]);
            let sum = build::add(graph, &[ep, em]);
            return build::mul(graph, &[half, sum]);
        }
        if op == ops.sin {
            let iu = build::mul(graph, &[i_unit, u]);
            let neg_iu = build::neg(graph, iu);
            let ep = build::call(graph, ops.exp, &[iu]);
            let em = build::call(graph, ops.exp, &[neg_iu]);
            let diff = build::sub(graph, ep, em);
            let two = graph.int(2);
            let two_i = build::mul(graph, &[two, i_unit]);
            let neg_one = graph.int(-1);
            let inv_two_i = build::pow(graph, two_i, neg_one);
            return build::mul(graph, &[diff, inv_two_i]);
        }
        if op == ops.cosh {
            let neg_u = build::neg(graph, u);
            let ep = build::call(graph, ops.exp, &[u]);
            let em = build::call(graph, ops.exp, &[neg_u]);
            let sum = build::add(graph, &[ep, em]);
            return build::mul(graph, &[half, sum]);
        }
        if op == ops.sinh {
            let neg_u = build::neg(graph, u);
            let ep = build::call(graph, ops.exp, &[u]);
            let em = build::call(graph, ops.exp, &[neg_u]);
            let diff = build::sub(graph, ep, em);
            return build::mul(graph, &[half, diff]);
        }
    }
    graph.node(op, &new_ch)
}
