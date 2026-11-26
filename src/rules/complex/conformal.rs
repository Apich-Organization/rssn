//! Symbolic conformal mappings, Schwarzian derivative, conformality check,
//! Rouché comparison, and Cauchy root bound.

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
use std::f64::consts::PI;

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let joukowsky = i.op(OpDescriptor::new("joukowsky", Arity::Fixed(1))
        .flags(OpFlags::HEAVY)
        .cost(10)
        .doc("`joukowsky(z)`: Joukowsky conformal map `(z + 1/z) / 2`."))?;
    let koebe = i.op(OpDescriptor::new("koebe", Arity::Fixed(1))
        .flags(OpFlags::HEAVY)
        .cost(10)
        .doc("`koebe(z)`: Koebe univalent extremal map `z / (1 - z)^2`."))?;
    let cayley = i.op(
        OpDescriptor::new("cayley_transform", Arity::Fixed(1))
            .flags(OpFlags::HEAVY)
            .cost(10)
            .doc("`cayley_transform(z)`: Cayley transform `(z - I) / (z + I)` mapping the upper half-plane to the unit disk."),
    )?;
    let cross_ratio = i.op(
        OpDescriptor::new("cross_ratio", Arity::Fixed(4))
            .flags(OpFlags::HEAVY)
            .cost(15)
            .doc("`cross_ratio(z1, z2, z3, z4)`: Möbius-invariant cross ratio `((z1 - z3)*(z2 - z4)) / ((z1 - z4)*(z2 - z3))`."),
    )?;
    let schwarzian = i.op(OpDescriptor::new("schwarzian_derivative", Arity::Fixed(2))
        .flags(OpFlags::HEAVY)
        .cost(60)
        .doc(
            "`schwarzian_derivative(f, z)`: Schwarzian derivative `(f'''/f') - (3/2)*(f''/f')^2`.",
        ))?;
    let is_conf = i.op(
        OpDescriptor::new("is_conformal_at", Arity::Fixed(3))
            .flags(OpFlags::HEAVY)
            .cost(50)
            .doc("`is_conformal_at(f, z, z0)`: returns `1` if holomorphic `f` has `f'(z0) != 0`, `0` if `f'(z0) == 0`."),
    )?;
    let rouche = i.op(
        OpDescriptor::new("rouche_compare", Arity::Fixed(5))
            .flags(OpFlags::HEAVY)
            .cost(80)
            .doc("`rouche_compare(f, g, z, center, radius)`: returns `1` if `|f(z)| > |g(z)|` on `|z - center| = radius`, `0` if `|f(z)| < |g(z)|`."),
    )?;
    let cauchy_bnd = i.op(
        OpDescriptor::new("cauchy_root_bound", Arity::Fixed(2))
            .flags(OpFlags::HEAVY)
            .cost(60)
            .doc("`cauchy_root_bound(poly, z)`: Cauchy upper bound `1 + max_{k<n} |c_k / c_n|` on the roots of polynomial `poly(z)`."),
    )?;
    let ops = Ops::of(i.graph()).ok_or_else(|| {
        RuleError::Invalid {
            rule: "complex/conformal".into(),
            reason: "needs elementary".into(),
        }
    })?;
    i.kernel(
        "complex/conformal",
        Tier::Reduce,
        ConformalKernel {
            joukowsky,
            koebe,
            cayley,
            cross_ratio,
            schwarzian,
            is_conf,
            rouche,
            cauchy_bnd,
            ops,
        },
    );
    Ok(())
}

struct ConformalKernel {
    joukowsky: OpId,
    koebe: OpId,
    cayley: OpId,
    cross_ratio: OpId,
    schwarzian: OpId,
    is_conf: OpId,
    rouche: OpId,
    cauchy_bnd: OpId,
    ops: Ops,
}

impl Kernel for ConformalKernel {
    fn ops(&self) -> Vec<OpId> {
        vec![
            self.joukowsky,
            self.koebe,
            self.cayley,
            self.cross_ratio,
            self.schwarzian,
            self.is_conf,
            self.rouche,
            self.cauchy_bnd,
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

impl ConformalKernel {
    fn reduce_opt(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Option<NodeId> {
        let op = cx.graph.op(node);
        let ch = cx.graph.children(node).to_vec();
        if op == self.joukowsky {
            let z = *ch.first()?;
            let neg_one = cx.graph.int(-1);
            let inv_z = build::pow(cx.graph, z, neg_one);
            let sum = build::add(cx.graph, &[z, inv_z]);
            let half = cx.graph.num(Number::fraction(1, 2)?);
            let res = build::mul(cx.graph, &[half, sum]);
            Some(cx.simplify(res))
        } else if op == self.koebe {
            let z = *ch.first()?;
            let one = cx.graph.int(1);
            let one_m_z = build::sub(cx.graph, one, z);
            let neg_two = cx.graph.int(-2);
            let inv_sq = build::pow(cx.graph, one_m_z, neg_two);
            let res = build::mul(cx.graph, &[z, inv_sq]);
            Some(cx.simplify(res))
        } else if op == self.cayley {
            let z = *ch.first()?;
            let i_unit = build::call(cx.graph, self.ops.unit, &[]);
            let num = build::sub(cx.graph, z, i_unit);
            let den = build::add(cx.graph, &[z, i_unit]);
            let neg_one = cx.graph.int(-1);
            let inv_den = build::pow(cx.graph, den, neg_one);
            let res = build::mul(cx.graph, &[num, inv_den]);
            Some(cx.simplify(res))
        } else if op == self.cross_ratio {
            let (z1, z2, z3, z4) = (*ch.first()?, *ch.get(1)?, *ch.get(2)?, *ch.get(3)?);
            let a = build::sub(cx.graph, z1, z3);
            let b = build::sub(cx.graph, z2, z4);
            let c = build::sub(cx.graph, z1, z4);
            let d = build::sub(cx.graph, z2, z3);
            let num = build::mul(cx.graph, &[a, b]);
            let den = build::mul(cx.graph, &[c, d]);
            let neg_one = cx.graph.int(-1);
            let inv_den = build::pow(cx.graph, den, neg_one);
            let res = build::mul(cx.graph, &[num, inv_den]);
            Some(cx.simplify(res))
        } else if op == self.schwarzian {
            let (f, z) = (*ch.first()?, *ch.get(1)?);
            self.schwarzian_reduce(cx, f, z)
        } else if op == self.is_conf {
            let (f, z, z0) = (*ch.first()?, *ch.get(1)?, *ch.get(2)?);
            self.is_conformal_reduce(cx, f, z, z0)
        } else if op == self.rouche {
            let (f, g, z, center, radius) = (
                *ch.first()?,
                *ch.get(1)?,
                *ch.get(2)?,
                *ch.get(3)?,
                *ch.get(4)?,
            );
            self.rouche_reduce(cx, f, g, z, center, radius)
        } else if op == self.cauchy_bnd {
            let (poly, z) = (*ch.first()?, *ch.get(1)?);
            self.cauchy_bound_reduce(cx, poly, z)
        } else {
            None
        }
    }

    fn schwarzian_reduce(
        &self,
        cx: &mut Cx<'_>,
        f: NodeId,
        z: NodeId,
    ) -> Option<NodeId> {
        let f_best = best(cx.graph, f)?;
        let d1_raw = derivative(cx.graph, f_best, z)?;
        let d1 = cx.simplify(d1_raw);
        if cx.graph.number_of(d1).is_some_and(Number::is_zero) {
            return None;
        }
        let d2_raw = derivative(cx.graph, d1, z)?;
        let d2 = cx.simplify(d2_raw);
        let d3_raw = derivative(cx.graph, d2, z)?;
        let d3 = cx.simplify(d3_raw);

        let neg_one = cx.graph.int(-1);
        let inv_d1 = build::pow(cx.graph, d1, neg_one);
        let term1 = build::mul(cx.graph, &[d3, inv_d1]);
        let ratio = build::mul(cx.graph, &[d2, inv_d1]);
        let ratio_sq = build::powi(cx.graph, ratio, 2);
        let three_halves = cx.graph.num(Number::fraction(3, 2)?);
        let term2 = build::mul(cx.graph, &[three_halves, ratio_sq]);
        let diff = build::sub(cx.graph, term1, term2);
        Some(cx.simplify(diff))
    }

    fn is_conformal_reduce(
        &self,
        cx: &mut Cx<'_>,
        f: NodeId,
        z: NodeId,
        z0: NodeId,
    ) -> Option<NodeId> {
        let f_best = best(cx.graph, f)?;
        let d1_raw = derivative(cx.graph, f_best, z)?;
        let d1 = cx.simplify(d1_raw);
        let at_raw = cx.graph.substitute(d1, z, z0);
        let at_val = cx.simplify(at_raw);
        let empty = HashMap::new();
        let val = eval_complex(cx.graph, at_val, &empty)?;
        if !val.re.is_finite() || !val.im.is_finite() {
            return None;
        }
        if val.norm() > 1e-12 {
            Some(cx.graph.int(1))
        } else {
            Some(cx.graph.int(0))
        }
    }

    fn rouche_reduce(
        &self,
        cx: &mut Cx<'_>,
        f: NodeId,
        g: NodeId,
        z: NodeId,
        center: NodeId,
        radius: NodeId,
    ) -> Option<NodeId> {
        let z_sym = cx.graph.symbol_of(z)?;
        let empty = HashMap::new();
        let c_val = eval_complex(cx.graph, center, &empty)?;
        let r_val = eval_complex(cx.graph, radius, &empty)?;
        if r_val.im.abs() > 1e-12 || r_val.re <= 0.0 {
            return None;
        }
        let r = r_val.re;
        let n_samples = 64_usize;
        let mut min_diff = f64::INFINITY;
        let mut max_diff = f64::NEG_INFINITY;
        for j in 0..n_samples {
            let theta = 2.0 * PI * (j as f64) / (n_samples as f64);
            let pt = c_val + Complex64::from_polar(r, theta);
            let mut bindings: HashMap<SymbolId, Complex64> = HashMap::new();
            bindings.insert(z_sym, pt);
            let fv = eval_complex(cx.graph, f, &bindings)?;
            let gv = eval_complex(cx.graph, g, &bindings)?;
            if !fv.re.is_finite() || !fv.im.is_finite() || !gv.re.is_finite() || !gv.im.is_finite()
            {
                return None;
            }
            let diff = fv.norm() - gv.norm();
            min_diff = min_diff.min(diff);
            max_diff = max_diff.max(diff);
        }
        if min_diff > 1e-10 {
            Some(cx.graph.int(1))
        } else if max_diff < -1e-10 {
            Some(cx.graph.int(0))
        } else {
            None
        }
    }

    fn cauchy_bound_reduce(
        &self,
        cx: &mut Cx<'_>,
        poly: NodeId,
        z: NodeId,
    ) -> Option<NodeId> {
        let z_sym = cx.graph.symbol_of(z)?;
        let p_best = best(cx.graph, poly)?;
        let zero = cx.graph.int(0);
        let empty = HashMap::new();
        let mut cur = p_best;
        let mut coeffs = Vec::new();
        let mut fact = 1.0_f64;
        for k in 0..=32_usize {
            let at0 = cx.graph.substitute(cur, z, zero);
            let at0_s = cx.simplify(at0);
            let val = eval_complex(cx.graph, at0_s, &empty)?;
            coeffs.push(val / fact);
            let next_raw = derivative(cx.graph, cur, z)?;
            let next = cx.simplify(next_raw);
            if cx.graph.number_of(next).is_some_and(Number::is_zero)
                && !cx.graph.depends_on(cx.graph.find(next), z_sym)
            {
                break;
            }
            cur = next;
            fact *= (k.saturating_add(1)) as f64;
            if k == 32 {
                return None;
            }
        }
        while coeffs.last().is_some_and(|c| c.norm() < 1e-14) {
            coeffs.pop();
        }
        let bnd = crate::kernels::complex_roots::cauchy_root_bound(&coeffs)?;
        let num = Number::from(bnd);
        Some(cx.graph.num(num))
    }
}
