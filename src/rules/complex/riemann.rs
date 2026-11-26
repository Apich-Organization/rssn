//! Riemann surface invariants and symbolic branch-cut jump evaluation.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::ToPrimitive;

use super::Ops;
use super::build;
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
use crate::rules::poly::best;

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let genus_op = i.op(
        OpDescriptor::new("riemann_genus", Arity::Fixed(2))
            .flags(OpFlags::HEAVY)
            .cost(50)
            .doc("`riemann_genus(n, d)`: topological genus of the superelliptic Riemann surface `w^n = prod_{k=1}^d (z - r_k)` with distinct roots by the Riemann-Hurwitz formula."),
    )?;
    let jump_op = i.op(
        OpDescriptor::new("branch_cut_jump", Arity::Fixed(3))
            .flags(OpFlags::HEAVY)
            .cost(50)
            .doc("`branch_cut_jump(f, z, p)`: symbolic discontinuity `lim_{eps -> 0+} (f(p + I*eps) - f(p - I*eps))` across a principal branch cut at `z = p`."),
    )?;
    let ops = Ops::of(i.graph()).ok_or_else(|| {
        RuleError::Invalid {
            rule: "complex/riemann".into(),
            reason: "needs elementary",
        }
    })?;
    i.kernel(
        "complex/riemann-genus",
        Tier::Reduce,
        RiemannGenus { op: genus_op },
    );
    i.kernel(
        "complex/branch-cut-jump",
        Tier::Reduce,
        BranchCutJump { op: jump_op, ops },
    );
    Ok(())
}

struct RiemannGenus {
    op: OpId,
}

impl Kernel for RiemannGenus {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let &[n_node, d_node] = cx.graph.children(node) else {
            return Outcome::Pass;
        };
        let (Some(n_int), Some(d_int)) = (
            cx.graph.number_of(n_node).and_then(Number::to_i64),
            cx.graph.number_of(d_node).and_then(Number::to_i64),
        ) else {
            return Outcome::Pass;
        };
        if n_int < 2 || d_int < 1 {
            return Outcome::Pass;
        }
        let n = BigInt::from(n_int);
        let d = BigInt::from(d_int);
        let g = n.gcd(&d);
        let total_ram = &d * (&n - 1) + (&n - g);
        let genus: BigInt = BigInt::from(1) - &n + total_ram / 2;
        if let Some(g_i64) = genus.to_i64() {
            Outcome::Equal(cx.graph.int(g_i64))
        } else {
            Outcome::Equal(cx.graph.num(Number::Int(genus)))
        }
    }
}

struct BranchCutJump {
    op: OpId,
    ops: Ops,
}

/// Extracts `a, b` if `arg` is `a * z + b` with `a` free of `z`.
fn linear_in(
    graph: &mut Graph,
    arg: NodeId,
    z: NodeId,
) -> Option<(NodeId, NodeId)> {
    let sym = graph.symbol_of(z)?;
    if graph.depends_on(graph.find(arg), sym) {
        if graph.same(arg, z) {
            return Some((graph.int(1), graph.int(0)));
        }
        let op = graph.op(arg);
        let children = graph.children(arg).to_vec();
        match (op, children.as_slice()) {
            | (core::ADD, terms) => {
                let mut a_parts = Vec::new();
                let mut b_parts = Vec::new();
                for &t in terms {
                    if !graph.depends_on(graph.find(t), sym) {
                        b_parts.push(t);
                    } else if let Some((ai, bi)) = linear_in(graph, t, z) {
                        a_parts.push(ai);
                        b_parts.push(bi);
                    } else {
                        return None;
                    }
                }
                Some((build::add(graph, &a_parts), build::add(graph, &b_parts)))
            },
            | (core::MUL, factors) => {
                let mut z_factor = None;
                let mut const_factors = Vec::new();
                for &f in factors {
                    if !graph.depends_on(graph.find(f), sym) {
                        const_factors.push(f);
                    } else if z_factor.is_none() {
                        z_factor = Some(f);
                    } else {
                        return None;
                    }
                }
                let (za, zb) = linear_in(graph, z_factor?, z)?;
                let c = build::mul(graph, &const_factors);
                Some((build::mul(graph, &[c, za]), build::mul(graph, &[c, zb])))
            },
            | _ => None,
        }
    } else {
        Some((graph.int(0), arg))
    }
}

impl Kernel for BranchCutJump {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let &[f, z, p] = cx.graph.children(node) else {
            return Outcome::Pass;
        };
        let Some(f_best) = best(cx.graph, f) else {
            return Outcome::Pass;
        };
        let Some(p_best) = best(cx.graph, p) else {
            return Outcome::Pass;
        };
        let op = cx.graph.op(f_best);
        let children = cx.graph.children(f_best).to_vec();
        let pi = build::call(cx.graph, self.ops.pi, &[]);
        let i_unit = build::call(cx.graph, self.ops.unit, &[]);
        let two = cx.graph.int(2);

        if op == self.ops.ln
            && let [arg] = children.as_slice()
        {
            let Some((a, b)) = linear_in(cx.graph, *arg, z) else {
                return Outcome::Pass;
            };
            let ap = build::mul(cx.graph, &[a, p_best]);
            let val = build::add(cx.graph, &[ap, b]);
            let val_s = cx.simplify(val);
            if cx.graph.number_of(val_s).is_some_and(|n| n.to_f64() < 0.0) {
                let jump = build::mul(cx.graph, &[two, i_unit, pi]);
                return Outcome::Equal(jump);
            }
            if cx.graph.number_of(val_s).is_some_and(|n| n.to_f64() > 0.0) {
                return Outcome::Equal(cx.graph.int(0));
            }
            return Outcome::Pass;
        }

        if (op == self.ops.sqrt || op == core::POW)
            && let Some((arg, alpha)) = match (op, children.as_slice()) {
                | (op_id, &[arg]) if op_id == self.ops.sqrt => {
                    let half = cx
                        .graph
                        .num(Number::fraction(1, 2).unwrap_or_else(|| Number::from(0)));
                    Some((arg, half))
                },
                | (core::POW, &[arg, alpha]) => Some((arg, alpha)),
                | _ => None,
            }
        {
            let Some((a, b)) = linear_in(cx.graph, arg, z) else {
                return Outcome::Pass;
            };
            let ap = build::mul(cx.graph, &[a, p_best]);
            let val = build::add(cx.graph, &[ap, b]);
            let val_s = cx.simplify(val);
            if cx.graph.number_of(val_s).is_some_and(|n| n.to_f64() < 0.0) {
                // (-r)^alpha_above - (-r)^alpha_below = 2*I*sin(pi*alpha) * |val|^alpha
                let neg_val = build::neg(cx.graph, val_s);
                let mag_pow = build::pow(cx.graph, neg_val, alpha);
                let pi_alpha = build::mul(cx.graph, &[pi, alpha]);
                let sin_pa = build::call(cx.graph, self.ops.sin, &[pi_alpha]);
                let jump = build::mul(cx.graph, &[two, i_unit, sin_pa, mag_pow]);
                return Outcome::Equal(jump);
            }
            if cx.graph.number_of(val_s).is_some_and(|n| n.to_f64() > 0.0) {
                return Outcome::Equal(cx.graph.int(0));
            }
        }
        Outcome::Pass
    }
}
