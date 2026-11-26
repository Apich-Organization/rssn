//! Symbolic complex dynamics: orbit iteration, periodic point multipliers,
//! attracting fixed point checks, and formal Böttcher Laurent series.

use super::Ops;
use super::build;
use super::eval_complex;
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
use crate::graph::SymbolId;
use crate::graph::Tier;
use crate::graph::rule::Installer;
use crate::rules::calculus::derivative;
use crate::rules::poly::best;
use num_complex::Complex64;
use std::collections::HashMap;

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let orbit = i.op(OpDescriptor::new("orbit_iterate", Arity::Fixed(4))
        .flags(OpFlags::HEAVY)
        .cost(50)
        .doc("`orbit_iterate(f, z, z0, n)`: symbolic `n`-fold iterate `f^{o n}(z0)`."))?;
    let mult = i.op(
        OpDescriptor::new("multiplier", Arity::Fixed(4))
            .flags(OpFlags::HEAVY)
            .cost(60)
            .doc("`multiplier(f, z, z0, period)`: symbolic periodic-point multiplier `d/dz (f^{o p})(z) |_{z=z0}`."),
    )?;
    let attract = i.op(
        OpDescriptor::new("is_attracting_fixed_point", Arity::Fixed(3))
            .flags(OpFlags::HEAVY)
            .cost(60)
            .doc("`is_attracting_fixed_point(f, z, z0)`: returns `1` if `f(z0) = z0` and `|f'(z0)| < 1`, `0` otherwise."),
    )?;
    let bottcher = i.op(
        OpDescriptor::new("bottcher_series", Arity::Fixed(3))
            .flags(OpFlags::HEAVY)
            .cost(70)
            .doc("`bottcher_series(c, z, order)`: formal Laurent series of the Böttcher coordinate `phi_c(z)` for `z^2 + c` up to `z^(-order)`."),
    )?;
    let ops = Ops::of(i.graph()).ok_or_else(|| {
        RuleError::Invalid {
            rule: "complex/dynamics".into(),
            reason: "needs elementary".into(),
        }
    })?;
    i.kernel(
        "complex/dynamics",
        Tier::Reduce,
        ComplexDynamics {
            orbit,
            mult,
            attract,
            bottcher,
            ops,
        },
    );
    Ok(())
}

struct ComplexDynamics {
    orbit: OpId,
    mult: OpId,
    attract: OpId,
    bottcher: OpId,
    ops: Ops,
}

impl Kernel for ComplexDynamics {
    fn ops(&self) -> Vec<OpId> {
        vec![self.orbit, self.mult, self.attract, self.bottcher]
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

impl ComplexDynamics {
    fn reduce_opt(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Option<NodeId> {
        let op = cx.graph.op(node);
        let ch = cx.graph.children(node).to_vec();
        if op == self.orbit {
            let (f, z, z0, n_node) = (*ch.first()?, *ch.get(1)?, *ch.get(2)?, *ch.get(3)?);
            let n = cx.graph.number_of(n_node).and_then(Number::to_i64)?;
            if !(0..=16).contains(&n) {
                return None;
            }
            let f_best = best(cx.graph, f)?;
            let mut cur = z0;
            for _ in 0..n {
                cur = cx.graph.substitute(f_best, z, cur);
                cur = cx.simplify(cur);
            }
            Some(cur)
        } else if op == self.mult {
            let (f, z, z0, p_node) = (*ch.first()?, *ch.get(1)?, *ch.get(2)?, *ch.get(3)?);
            let p = cx.graph.number_of(p_node).and_then(Number::to_i64)?;
            if !(1..=32).contains(&p) {
                return None;
            }
            let f_best = best(cx.graph, f)?;
            let df = derivative(cx.graph, f_best, z)?;
            let df_s = cx.simplify(df);
            // Chain rule product: prod_{k=0}^{p-1} f'(z_k) where z_0 = z0, z_{k+1} = f(z_k)
            let mut factors = Vec::with_capacity(p as usize);
            let mut cur_z = z0;
            for k in 0..p {
                let df_at = cx.graph.substitute(df_s, z, cur_z);
                factors.push(cx.simplify(df_at));
                if k + 1 < p {
                    let next_z = cx.graph.substitute(f_best, z, cur_z);
                    cur_z = cx.simplify(next_z);
                }
            }
            let prod = build::mul(cx.graph, &factors);
            Some(cx.simplify(prod))
        } else if op == self.attract {
            let (f, z, z0) = (*ch.first()?, *ch.get(1)?, *ch.get(2)?);
            let z_sym = cx.graph.symbol_of(z)?;
            let empty = HashMap::new();
            let z0_val = eval_complex(cx.graph, z0, &empty)?;
            let mut bindings: HashMap<SymbolId, Complex64> = HashMap::new();
            bindings.insert(z_sym, z0_val);
            let f_z0 = eval_complex(cx.graph, f, &bindings)?;
            if (f_z0 - z0_val).norm() > 1e-9 {
                return Some(cx.graph.int(0));
            }
            let f_best = best(cx.graph, f)?;
            let df = derivative(cx.graph, f_best, z)?;
            let df_z0 = eval_complex(cx.graph, df, &bindings)?;
            if df_z0.norm() < 1.0 - 1e-12 {
                Some(cx.graph.int(1))
            } else {
                Some(cx.graph.int(0))
            }
        } else if op == self.bottcher {
            let (c, z, order_node) = (*ch.first()?, *ch.get(1)?, *ch.get(2)?);
            self.bottcher_series_reduce(cx, c, z, order_node)
        } else {
            None
        }
    }

    fn bottcher_series_reduce(
        &self,
        cx: &mut Cx<'_>,
        c: NodeId,
        z: NodeId,
        order_node: NodeId,
    ) -> Option<NodeId> {
        let _ = self.ops;
        let order = cx.graph.number_of(order_node).and_then(Number::to_i64)? as usize;
        if order > 15 {
            return None;
        }
        // phi_c(z) = z + sum_{k=1}^order a_k * z^(-k)
        // Only odd k = 2p - 1 are non-zero: a_1 = c/2,
        // a_{2p+1} = 1/2 * (sum_{r=1}^p C(-r, p-r) c^{p-r} a_r - sum_{j=1}^{2p-1} a_j a_{2p-j})
        let zero = cx.graph.int(0);
        let mut a: Vec<NodeId> = vec![zero; order.saturating_add(1)];
        if order >= 1 {
            let half = cx.graph.num(Number::fraction(1, 2)?);
            let a1_raw = build::mul(cx.graph, &[half, c]);
            let a1 = cx.simplify(a1_raw);
            if let Some(slot) = a.get_mut(1) {
                *slot = a1;
            }
        }
        let mut p = 1_usize;
        while 2 * p + 1 <= order {
            let mut rhs_terms = Vec::new();
            for r in 1..=p {
                let s = p - r;
                // C(-r, s) = (-1)^s * C(r + s - 1, s)
                let binom = signed_neg_binom(r as i64, s as i64);
                let ar = *a.get(r)?;
                if binom != 0 && !cx.graph.number_of(ar).is_some_and(Number::is_zero) {
                    let b_node = cx.graph.int(binom);
                    let c_pow = if s == 0 {
                        cx.graph.int(1)
                    } else {
                        build::powi(cx.graph, c, s as i64)
                    };
                    rhs_terms.push(build::mul(cx.graph, &[b_node, c_pow, ar]));
                }
            }
            let rhs = build::add(cx.graph, &rhs_terms);
            let mut conv_terms = Vec::new();
            for j in 1..=(2 * p - 1) {
                let aj = *a.get(j)?;
                let al = *a.get(2 * p - j)?;
                if !cx.graph.number_of(aj).is_some_and(Number::is_zero)
                    && !cx.graph.number_of(al).is_some_and(Number::is_zero)
                {
                    conv_terms.push(build::mul(cx.graph, &[aj, al]));
                }
            }
            let conv = build::add(cx.graph, &conv_terms);
            let diff = build::sub(cx.graph, rhs, conv);
            let half = cx.graph.num(Number::fraction(1, 2)?);
            let coeff_raw = build::mul(cx.graph, &[half, diff]);
            let coeff = cx.simplify(coeff_raw);
            if let Some(slot) = a.get_mut(2 * p + 1) {
                *slot = coeff;
            }
            p += 1;
        }

        let mut terms = vec![z];
        for (k, &ak) in a.iter().enumerate().skip(1) {
            if !cx.graph.number_of(ak).is_some_and(Number::is_zero) {
                let z_neg_k = build::powi(cx.graph, z, -(k as i64));
                terms.push(build::mul(cx.graph, &[ak, z_neg_k]));
            }
        }
        let sum = build::add(cx.graph, &terms);
        Some(cx.simplify(sum))
    }
}

fn signed_neg_binom(
    r: i64,
    s: usize_or_i64,
) -> i64 {
    let mut num = 1_i64;
    let mut den = 1_i64;
    for j in 0..s {
        num *= -r - j;
        den *= j + 1;
    }
    num / den
}

#[allow(non_camel_case_types)]
type usize_or_i64 = i64;
