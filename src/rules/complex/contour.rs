//! Symbolic contour descriptors, residue at infinity, Cauchy principal value,
//! and residue-theorem real & trigonometric integral evaluators.

use num_complex::Complex64;

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
use crate::graph::Tier;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::rules::poly::best;

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    i.op(OpDescriptor::new("arc", Arity::Fixed(4)).doc(
        "`arc(c, r, t0, t1)`: circular arc contour `c + r * exp(I*t)` for `t` from `t0` to `t1`.",
    ))?;
    let res_inf = i.op(
        OpDescriptor::new("residue_at_infinity", Arity::Fixed(2))
            .flags(OpFlags::HEAVY)
            .cost(100)
            .doc("`residue_at_infinity(f, z)`: residue of `f(z)` at `z = infinity`, defined as `residue(-f(1/w)/w^2, w, 0)`."),
    )?;
    let cauchy_pv = i.op(
        OpDescriptor::new("cauchy_pv", Arity::Fixed(5))
            .flags(OpFlags::HEAVY)
            .cost(100)
            .doc("`cauchy_pv(f, x, a, b, pole)`: Cauchy principal value `P.V. int_a^b f(x) dx` across a simple pole in `(a, b)`."),
    )?;
    let res_real = i.op(
        OpDescriptor::new("residue_real_integral", Arity::Fixed(2))
            .flags(OpFlags::HEAVY)
            .cost(100)
            .doc("`residue_real_integral(f, x)`: evaluates `int_{-oo}^{oo} f(x) dx` by summing upper half-plane residues."),
    )?;
    let res_trig = i.op(
        OpDescriptor::new("residue_trig_integral", Arity::Fixed(2))
            .flags(OpFlags::HEAVY)
            .cost(100)
            .doc("`residue_trig_integral(f, theta)`: evaluates `int_0^{2*pi} f(theta) dtheta` by unit-circle residue substitution."),
    )?;
    let ops = Ops::of(i.graph()).ok_or_else(|| {
        RuleError::Invalid {
            rule: "complex/contour".into(),
            reason: "needs elementary",
        }
    })?;
    i.kernel(
        "complex/contour-ext",
        Tier::Reduce,
        ContourExt {
            res_inf,
            cauchy_pv,
            res_real,
            res_trig,
            ops,
        },
    );
    Ok(())
}

struct ContourExt {
    res_inf: OpId,
    cauchy_pv: OpId,
    res_real: OpId,
    res_trig: OpId,
    ops: Ops,
}

impl Kernel for ContourExt {
    fn ops(&self) -> Vec<OpId> {
        vec![self.res_inf, self.cauchy_pv, self.res_real, self.res_trig]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let op = cx.graph.op(node);
        let children = cx.graph.children(node).to_vec();
        let res = if op == self.res_inf
            && let [f, z] = children.as_slice()
        {
            self.residue_at_infinity(cx, *f, *z)
        } else if op == self.cauchy_pv
            && let [f, x, a, b, p] = children.as_slice()
        {
            self.cauchy_principal_value(cx, *f, *x, *a, *b, *p)
        } else if op == self.res_real
            && let [f, x] = children.as_slice()
        {
            self.real_integral(cx, *f, *x)
        } else if op == self.res_trig
            && let [f, theta] = children.as_slice()
        {
            self.trig_integral(cx, *f, *theta)
        } else {
            None
        };
        res.map_or(Outcome::Pass, Outcome::Equal)
    }
}

impl ContourExt {
    fn residue_at_infinity(
        &self,
        cx: &mut Cx<'_>,
        f: NodeId,
        z: NodeId,
    ) -> Option<NodeId> {
        let residue_op = cx.graph.ops().lookup("residue")?;
        let w_sym = cx.graph.interner_mut().fresh_symbol("w");
        let w = cx.graph.symbol_node(w_sym);
        let inv_w = build::powi(cx.graph, w, -1);
        let inv_w2 = build::powi(cx.graph, w, -2);
        let fw = cx.graph.substitute(f, z, inv_w);
        let neg_fw = build::neg(cx.graph, fw);
        let integrand = build::mul(cx.graph, &[neg_fw, inv_w2]);
        let zero = cx.graph.int(0);
        let req = build::call(cx.graph, residue_op, &[integrand, w, zero]);
        let ans = cx.simplify(req);
        (cx.graph.op(ans) != residue_op).then_some(ans)
    }

    fn cauchy_principal_value(
        &self,
        cx: &mut Cx<'_>,
        f: NodeId,
        x: NodeId,
        a: NodeId,
        b: NodeId,
        pole: NodeId,
    ) -> Option<NodeId> {
        // Handle simple rational functions of the form c / (x - pole) + smooth part
        // by subtracting r / (x - pole) where r = residue(f, x, pole).
        let residue_op = cx.graph.ops().lookup("residue")?;
        let int_op = cx.graph.ops().lookup("defint")?;
        let res_req = build::call(cx.graph, residue_op, &[f, x, pole]);
        let r = cx.simplify(res_req);
        if cx.graph.op(r) == residue_op {
            return None;
        }
        let x_m_p = build::sub(cx.graph, x, pole);
        let x_m_p_inv = build::powi(cx.graph, x_m_p, -1);
        let sing = build::mul(cx.graph, &[r, x_m_p_inv]);
        let smooth_raw = build::sub(cx.graph, f, sing);
        let smooth = cx.simplify(smooth_raw);
        let pv_sing = if cx.graph.number_of(r).is_some_and(Number::is_zero) {
            cx.graph.int(0)
        } else {
            let b_m_p = build::sub(cx.graph, b, pole);
            let a_m_p = build::sub(cx.graph, a, pole);
            let a_m_p_inv = build::powi(cx.graph, a_m_p, -1);
            let ratio = build::mul(cx.graph, &[b_m_p, a_m_p_inv]);
            let abs_ratio = build::call(cx.graph, self.ops.abs, &[ratio]);
            let ln_ratio = build::call(cx.graph, self.ops.ln, &[abs_ratio]);
            build::mul(cx.graph, &[r, ln_ratio])
        };
        let smooth_int = if cx.graph.number_of(smooth).is_some_and(Number::is_zero) {
            cx.graph.int(0)
        } else {
            let req = build::call(cx.graph, int_op, &[smooth, x, a, b]);
            let out = cx.simplify(req);
            if cx.graph.op(out) == int_op {
                return None;
            }
            out
        };
        let total = build::add(cx.graph, &[pv_sing, smooth_int]);
        Some(cx.simplify(total))
    }

    fn real_integral(
        &self,
        cx: &mut Cx<'_>,
        f: NodeId,
        x: NodeId,
    ) -> Option<NodeId> {
        let poles_op = cx.graph.ops().lookup("poles")?;
        let residue_op = cx.graph.ops().lookup("residue")?;
        let poles_req = build::call(cx.graph, poles_op, &[f, x]);
        let poles_list = cx.simplify(poles_req);
        if cx.graph.op(poles_list) != core::LIST {
            return None;
        }
        let entries = cx.graph.children(poles_list).to_vec();
        let mut upper_residues = Vec::new();
        let none = std::collections::HashMap::new();
        for entry in entries {
            if cx.graph.op(entry) != core::LIST {
                return None;
            }
            let pair = cx.graph.children(entry).to_vec();
            let p_node = *pair.first()?;
            let p_val = eval_complex(cx.graph, p_node, &none)?;
            if p_val.im.abs() < 1e-12 {
                // Real pole on axis: not a standard improper integral without PV.
                return None;
            }
            if p_val.im > 0.0 {
                let res_req = build::call(cx.graph, residue_op, &[f, x, p_node]);
                let r = cx.simplify(res_req);
                if cx.graph.op(r) == residue_op {
                    return None;
                }
                upper_residues.push(r);
            }
        }
        let sum = build::add(cx.graph, &upper_residues);
        let two = cx.graph.int(2);
        let pi = build::call(cx.graph, self.ops.pi, &[]);
        let i_unit = build::call(cx.graph, self.ops.unit, &[]);
        let ans_raw = build::mul(cx.graph, &[two, pi, i_unit, sum]);
        let ans = cx.simplify(ans_raw);

        // Numerical sanity check against R = 1000 real integral when evaluable.
        let xs = cx.graph.symbol_of(x)?;
        if let Some(ans_val) = eval_complex(cx.graph, ans, &none) {
            let num_val = crate::kernels::integrate::gauss_kronrod_any(
                |t| {
                    let mut m = std::collections::HashMap::new();
                    m.insert(xs, Complex64::new(t, 0.0));
                    eval_complex(cx.graph, f, &m).map_or(0.0, |v| v.re)
                },
                -500.0,
                500.0,
                1e-6,
                24,
            )
            .value;
            if (ans_val.re - num_val).abs() > 1e-2 * (1.0 + ans_val.norm())
                || ans_val.im.abs() > 1e-6
            {
                return None;
            }
        }
        Some(ans)
    }

    fn trig_integral(
        &self,
        cx: &mut Cx<'_>,
        f: NodeId,
        theta: NodeId,
    ) -> Option<NodeId> {
        let contour_op = cx.graph.ops().lookup("contour_integral")?;
        let circle_op = cx.graph.ops().lookup("circle")?;
        let f_best = best(cx.graph, f)?;
        let z_sym = cx.graph.interner_mut().fresh_symbol("zc");
        let z = cx.graph.symbol_node(z_sym);
        let i_unit = build::call(cx.graph, self.ops.unit, &[]);
        let inv_z = build::powi(cx.graph, z, -1);
        let two = cx.graph.int(2);
        let half = cx.graph.num(Number::fraction(1, 2)?);
        let two_i = build::mul(cx.graph, &[two, i_unit]);
        let two_i_inv = build::powi(cx.graph, two_i, -1);

        let z_plus_inv = build::add(cx.graph, &[z, inv_z]);
        let cos_z = build::mul(cx.graph, &[half, z_plus_inv]);
        let z_minus_inv = build::sub(cx.graph, z, inv_z);
        let sin_z = build::mul(cx.graph, &[two_i_inv, z_minus_inv]);
        let sub2 = replace_trig_calls(
            cx.graph,
            f_best,
            self.ops.cos,
            self.ops.sin,
            theta,
            cos_z,
            sin_z,
        );
        let theta_sym = cx.graph.symbol_of(theta)?;
        if cx.graph.depends_on(cx.graph.find(sub2), theta_sym) {
            return None;
        }
        let integrand_raw = build::mul(cx.graph, &[sub2, inv_z]);
        let mut gens = crate::rules::poly::repr::Gens::default();
        let integrand = if let Some(r) = crate::rules::poly::ratio(
            cx.graph,
            &mut gens,
            integrand_raw,
            crate::rules::poly::repr::Limits::default(),
        ) {
            let numer = crate::rules::poly::repr::to_term(cx.graph, &gens, &r.numer);
            let denom = crate::rules::poly::repr::to_term(cx.graph, &gens, &r.denom);
            let neg_one = cx.graph.int(-1);
            let inv_den = build::pow(cx.graph, denom, neg_one);
            let rat = build::mul(cx.graph, &[numer, inv_den]);
            cx.simplify(rat)
        } else {
            cx.simplify(integrand_raw)
        };
        let zero = cx.graph.int(0);
        let one = cx.graph.int(1);
        let unit_circle = build::call(cx.graph, circle_op, &[zero, one]);
        let req = build::call(cx.graph, contour_op, &[integrand, z, unit_circle]);
        let ans = cx.simplify(req);
        if cx.graph.op(ans) == contour_op {
            return None;
        }
        let neg_i = build::neg(cx.graph, i_unit);
        let scaled = build::mul(cx.graph, &[neg_i, ans]);
        Some(cx.simplify(scaled))
    }
}

fn replace_trig_calls(
    graph: &mut crate::graph::Graph,
    node: NodeId,
    cos_op: OpId,
    sin_op: OpId,
    theta: NodeId,
    cos_z: NodeId,
    sin_z: NodeId,
) -> NodeId {
    let op = graph.op(node);
    let ch = graph.children(node).to_vec();
    if ch.len() == 1
        && let Some(&arg) = ch.first()
        && graph.find(arg) == graph.find(theta)
    {
        if op == cos_op {
            return cos_z;
        }
        if op == sin_op {
            return sin_z;
        }
    }
    if ch.is_empty() {
        return node;
    }
    let mut changed = false;
    let mut new_ch = Vec::with_capacity(ch.len());
    for c in ch {
        let nc = replace_trig_calls(graph, c, cos_op, sin_op, theta, cos_z, sin_z);
        if nc != c {
            changed = true;
        }
        new_ch.push(nc);
    }
    if changed {
        graph.node(op, &new_ch)
    } else {
        node
    }
}
