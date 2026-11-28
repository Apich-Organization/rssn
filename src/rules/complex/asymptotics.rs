//! Symbolic Laurent principal/analytic parts, Mittag-Leffler partial-fraction
//! expansion, and saddle-point asymptotic approximation.

use super::Ops;
use super::build;
use crate::graph::Arity;
use crate::graph::Cx;
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
use crate::rules::calculus::laurent_expansion;
use crate::rules::poly::best;
use num_traits::Zero;

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let prin = i.op(
        OpDescriptor::new("laurent_principal_part", Arity::Fixed(4))
            .flags(OpFlags::HEAVY)
            .cost(80)
            .doc("`laurent_principal_part(f, z, a, order)`: negative-power terms of the Laurent series of `f(z)` at `z = a`."),
    )?;
    let anal = i.op(
        OpDescriptor::new("laurent_analytic_part", Arity::Fixed(4))
            .flags(OpFlags::HEAVY)
            .cost(80)
            .doc("`laurent_analytic_part(f, z, a, order)`: non-negative-power terms of the Laurent series of `f(z)` at `z = a`."),
    )?;
    let saddle = i.op(
        OpDescriptor::new("saddle_point_asymptotic", Arity::Fixed(5))
            .flags(OpFlags::HEAVY)
            .cost(90)
            .doc("`saddle_point_asymptotic(g, phi, z, z0, lambda)`: leading saddle-point asymptotic approximation of `int g(z) exp(lambda * phi(z)) dz` near `z0`."),
    )?;
    let ml = i.op(
        OpDescriptor::new("mittag_leffler_expansion", Arity::Fixed(2))
            .flags(OpFlags::HEAVY)
            .cost(90)
            .doc("`mittag_leffler_expansion(f, z)`: Mittag-Leffler principal-part expansion over the poles of rational `f(z)`."),
    )?;
    let ops = Ops::of(i.graph()).ok_or_else(|| {
        RuleError::Invalid {
            rule: "complex/asymptotics".into(),
            reason: "needs elementary",
        }
    })?;
    i.kernel(
        "complex/asymptotics",
        Tier::Reduce,
        Asymptotics {
            prin,
            anal,
            saddle,
            ml,
            ops,
        },
    );
    Ok(())
}

struct Asymptotics {
    prin: OpId,
    anal: OpId,
    saddle: OpId,
    ml: OpId,
    ops: Ops,
}

impl Kernel for Asymptotics {
    fn ops(&self) -> Vec<OpId> {
        vec![self.prin, self.anal, self.saddle, self.ml]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let op = cx.graph.op(node);
        let children = cx.graph.children(node).to_vec();
        let res = if op == self.prin
            && let [f, z, a, order] = children.as_slice()
        {
            Self::split_laurent(cx, *f, *z, *a, *order, true)
        } else if op == self.anal
            && let [f, z, a, order] = children.as_slice()
        {
            Self::split_laurent(cx, *f, *z, *a, *order, false)
        } else if op == self.saddle
            && let [g, phi, z, z0, lambda] = children.as_slice()
        {
            self.saddle_point(cx, *g, *phi, *z, *z0, *lambda)
        } else if op == self.ml
            && let [f, z] = children.as_slice()
        {
            Self::mittag_leffler(cx, *f, *z)
        } else {
            None
        };
        res.map_or(Outcome::Pass, Outcome::Equal)
    }
}

impl Asymptotics {
    fn split_laurent(
        cx: &mut Cx<'_>,
        f: NodeId,
        z: NodeId,
        a: NodeId,
        order_node: NodeId,
        principal: bool,
    ) -> Option<NodeId> {
        let order = cx.graph.number_of(order_node).and_then(Number::to_i64)?;
        if !(0..=32).contains(&order) {
            return None;
        }
        let (start_pow, coeffs) = laurent_expansion(cx, f, z, a, order)?;
        let z_minus_a = build::sub(cx.graph, z, a);
        let mut kept = Vec::new();
        for (idx, coeff) in coeffs.into_iter().enumerate() {
            if cx.graph.number_of(coeff).is_some_and(Number::is_zero) {
                continue;
            }
            let power = start_pow + (idx as i64);
            let is_neg = power < 0;
            if is_neg != principal {
                continue;
            }
            let term = if power == 0 {
                coeff
            } else {
                let zp = build::powi(cx.graph, z_minus_a, power);
                build::mul(cx.graph, &[coeff, zp])
            };
            kept.push(term);
        }
        let sum = build::add(cx.graph, &kept);
        Some(cx.simplify(sum))
    }

    fn saddle_point(
        &self,
        cx: &mut Cx<'_>,
        g: NodeId,
        phi: NodeId,
        z: NodeId,
        z0: NodeId,
        lambda: NodeId,
    ) -> Option<NodeId> {
        let phi_best = best(cx.graph, phi)?;
        let d1_raw = derivative(cx.graph, phi_best, z)?;
        let d1 = cx.simplify(d1_raw);
        let d1_at_raw = cx.graph.substitute(d1, z, z0);
        let d1_at = cx.simplify(d1_at_raw);
        if !cx.graph.number_of(d1_at).is_some_and(Number::is_zero) {
            return None;
        }
        let d2_raw = derivative(cx.graph, d1, z)?;
        let d2 = cx.simplify(d2_raw);
        let d2_at_raw = cx.graph.substitute(d2, z, z0);
        let d2_at = cx.simplify(d2_at_raw);
        let g_at_raw = cx.graph.substitute(g, z, z0);
        let g_at = cx.simplify(g_at_raw);
        let phi_at_raw = cx.graph.substitute(phi_best, z, z0);
        let phi_at = cx.simplify(phi_at_raw);

        let neg_two = cx.graph.int(-2);
        let pi = build::call(cx.graph, self.ops.pi, &[]);
        let num = build::mul(cx.graph, &[neg_two, pi]);
        let den = build::mul(cx.graph, &[lambda, d2_at]);
        let den_inv = build::powi(cx.graph, den, -1);
        let ratio = build::mul(cx.graph, &[num, den_inv]);
        let sqrt_factor = build::call(cx.graph, self.ops.sqrt, &[ratio]);
        let l_phi = build::mul(cx.graph, &[lambda, phi_at]);
        let exp_factor = build::call(cx.graph, self.ops.exp, &[l_phi]);
        let prod = build::mul(cx.graph, &[sqrt_factor, g_at, exp_factor]);
        Some(cx.simplify(prod))
    }

    fn mittag_leffler(
        cx: &mut Cx<'_>,
        f: NodeId,
        z: NodeId,
    ) -> Option<NodeId> {
        let poles_op = cx.graph.ops().lookup("poles")?;
        let poles_req = build::call(cx.graph, poles_op, &[f, z]);
        let poles_list = cx.simplify(poles_req);
        if cx.graph.op(poles_list) != core::LIST {
            return None;
        }
        let entries = cx.graph.children(poles_list).to_vec();
        let mut parts = Vec::new();
        for entry in entries {
            if cx.graph.op(entry) != core::LIST {
                return None;
            }
            let pair = cx.graph.children(entry).to_vec();
            let p_node = *pair.first()?;
            let order = cx
                .graph
                .number_of(*pair.get(1)?)
                .and_then(Number::to_i64)
                .unwrap_or(1)
                .max(1);
            let order_node = cx.graph.int(order + 1);
            if let Some(pp) = Self::split_laurent(cx, f, z, p_node, order_node, true) {
                parts.push(pp);
            }
        }
        // Include polynomial quotient part if deg(numer) >= deg(denom)
        let mut gens = crate::rules::poly::repr::Gens::default();
        let gz = gens.index(cx.graph, z);
        if let Some(r) = crate::rules::poly::ratio(
            cx.graph,
            &mut gens,
            f,
            crate::rules::poly::repr::Limits::default(),
        ) && let (Some(num_u), Some(den_u)) =
            (r.numer.univariate_in(gz), r.denom.univariate_in(gz))
            && !den_u.is_empty()
            && num_u.len() >= den_u.len()
        {
            // Polynomial long division over Q when possible
            if let (Some(nq), Some(dq)) = (
                num_u
                    .iter()
                    .map(Number::to_rational)
                    .collect::<Option<Vec<_>>>(),
                den_u
                    .iter()
                    .map(Number::to_rational)
                    .collect::<Option<Vec<_>>>(),
            ) {
                let mut rem = nq;
                let d_deg = dq.len() - 1;
                let lead_d = dq[d_deg].clone();
                if !lead_d.is_zero() {
                    let q_len = rem.len() - d_deg;
                    let mut q_coeffs =
                        vec![num_rational::BigRational::from_integer(0.into()); q_len];
                    for k in (0..q_len).rev() {
                        let c = &rem[k + d_deg] / &lead_d;
                        q_coeffs[k] = c.clone();
                        for j in 0..=d_deg {
                            let sub = &c * &dq[j];
                            rem[k + j] = &rem[k + j] - sub;
                        }
                    }
                    for (k, qc) in q_coeffs.into_iter().enumerate() {
                        if !qc.is_zero() {
                            let cn = cx.graph.num(Number::rat(qc));
                            if k == 0 {
                                parts.push(cn);
                            } else {
                                let zk = build::powi(cx.graph, z, k as i64);
                                parts.push(build::mul(cx.graph, &[cn, zk]));
                            }
                        }
                    }
                }
            }
        }
        let sum = build::add(cx.graph, &parts);
        Some(cx.simplify(sum))
    }
}
