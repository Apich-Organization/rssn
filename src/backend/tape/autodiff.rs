//! Reverse accumulation with a directional derivative of every adjoint.
//! Local analytic rules share the dual/hyper-dual interpreter, so new tape
//! instructions cannot silently acquire a different derivative convention.
use super::Dual;
use super::HyperDual;
use super::Inst;
use super::Reg;
use super::Tape;
use super::apply;
use super::int_pow;
use super::operands;

impl Tape {
    /// Vector-Jacobian product `weightsᵀ J`, with respect to inputs.
    ///
    /// One value pass and one reverse sweep, with storage linear in tape size.
    /// Parameters are held constant. Repeated outputs and gathered inputs add
    /// their adjoints. Only selected branches receive adjoints; comparison and
    /// gather-index derivatives are zero. Opaque calls use the same central
    /// differences as forward-mode AD. A wrong weight count, invalid selected
    /// output, or NaN active output returns a vector of NaNs.
    #[must_use]
    pub fn vjp(
        &self,
        inputs: &[f64],
        params: &[f64],
        weights: &[f64],
    ) -> Vec<f64> {
        self.reverse(inputs, params, weights, None).0
    }

    /// Gradient of output `k`, in one value pass and one reverse sweep.
    #[must_use]
    pub fn gradient_reverse(
        &self,
        k: usize,
        inputs: &[f64],
        params: &[f64],
    ) -> Vec<f64> {
        let Some(weights) = self.output_seed(k) else {
            return vec![f64::NAN; inputs.len()];
        };
        self.vjp(inputs, params, &weights)
    }

    /// `H_k(inputs) direction` without assembling a Hessian.
    ///
    /// Forward-over-reverse AD uses one directional value pass and one reverse
    /// sweep of dual adjoints. Cost and storage are linear in tape size for
    /// bounded-arity analytic instructions, independent of input count.
    /// Opaque calls have no second-order rule: affected entries are NaN.
    /// Non-smooth functions use the documented hyper-dual conventions.
    /// An invalid output or direction length returns a vector of NaNs.
    #[must_use]
    pub fn hessian_vector_product(
        &self,
        k: usize,
        inputs: &[f64],
        params: &[f64],
        direction: &[f64],
    ) -> Vec<f64> {
        let Some(weights) = self.output_seed(k) else {
            return vec![f64::NAN; inputs.len()];
        };
        self.reverse(inputs, params, &weights, Some(direction)).1
    }

    fn output_seed(
        &self,
        k: usize,
    ) -> Option<Vec<f64>> {
        let mut weights = vec![0.0; self.outputs.len()];
        *weights.get_mut(k)? = 1.0;
        Some(weights)
    }

    fn reverse(
        &self,
        inputs: &[f64],
        params: &[f64],
        weights: &[f64],
        direction: Option<&[f64]>,
    ) -> (Vec<f64>, Vec<f64>) {
        let invalid = || (vec![f64::NAN; inputs.len()], vec![f64::NAN; inputs.len()]);
        if weights.len() != self.outputs.len() || direction.is_some_and(|v| v.len() != inputs.len())
        {
            return invalid();
        }
        let mut values = Vec::with_capacity(self.insts.len());
        let mut out = vec![Dual(0.0, 0.0); self.outputs.len()];
        if let Some(direction) = direction {
            let seeded: Vec<_> = inputs
                .iter()
                .zip(direction)
                .map(|(&x, &d)| Dual(x, d))
                .collect();
            let params: Vec<_> = params.iter().map(|&x| Dual(x, 0.0)).collect();
            self.run(&seeded, &params, &mut values, &mut out);
        } else {
            let mut plain = Vec::with_capacity(self.insts.len());
            let mut plain_out = vec![0.0; self.outputs.len()];
            self.run(inputs, params, &mut plain, &mut plain_out);
            values.extend(plain.into_iter().map(|x| Dual(x, 0.0)));
            for (dst, x) in out.iter_mut().zip(plain_out) {
                dst.0 = x;
            }
        }
        let mut adjoints = vec![Dual(0.0, 0.0); self.insts.len()];
        for ((&reg, &weight), value) in self.outputs.iter().zip(weights).zip(&out) {
            if weight == 0.0 {
                continue;
            }
            if value.0.is_nan() || weight.is_nan() {
                return invalid();
            }
            let Some(adjoint) = adjoints.get_mut(reg as usize) else {
                return invalid();
            };
            adjoint.0 += weight;
        }
        let mut gradient = vec![0.0; inputs.len()];
        let mut product = vec![0.0; inputs.len()];
        let mut local = Vec::new();
        for (index, inst) in self.insts.iter().enumerate().rev() {
            let bar = adjoints[index];
            if bar.0 == 0.0 && bar.1 == 0.0 {
                continue;
            }
            match *inst {
                | Inst::Input(i) => accumulate_input(&mut gradient, &mut product, i as usize, bar),
                | Inst::Gather(base, reg) => {
                    let offset = values.get(reg as usize).map_or(f64::NAN, |x| x.0).trunc();
                    if offset.is_finite() && (0.0..=1e18).contains(&offset) {
                        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                        if let Some(i) = (base as usize).checked_add(offset as usize)
                            && i < self.inputs
                        {
                            accumulate_input(&mut gradient, &mut product, i, bar);
                        }
                    }
                },
                | Inst::Const(_) | Inst::Param(_) | Inst::Cmp(..) => {},
                | Inst::Select(c, a, b) => {
                    let condition = values.get(c as usize).map_or(f64::NAN, |x| x.0);
                    let selected = if condition == 0.0 { b } else { a };
                    accumulate_register(&mut adjoints, selected, bar);
                },
                | _ => {
                    let registers = operands(inst, &self.args);
                    local.clear();
                    local.extend(registers.iter().map(|&r| {
                        values
                            .get(r as usize)
                            .copied()
                            .unwrap_or(Dual(f64::NAN, f64::NAN))
                    }));
                    for (operand, &reg) in registers.iter().enumerate() {
                        let partial = local_partial(inst, &local, operand, direction.is_some());
                        let contribution = Dual(
                            product_unless_zero(bar.0, partial.0),
                            product_unless_zero(bar.1, partial.0)
                                + product_unless_zero(bar.0, partial.1),
                        );
                        accumulate_register(&mut adjoints, reg, contribution);
                    }
                },
            }
        }
        (gradient, product)
    }
}

pub(super) const fn product_unless_zero(
    weight: f64,
    value: f64,
) -> f64 {
    if weight == 0.0 {
        0.0
    } else {
        weight * value
    }
}

fn accumulate_input(
    gradient: &mut [f64],
    product: &mut [f64],
    i: usize,
    bar: Dual,
) {
    if let (Some(g), Some(h)) = (gradient.get_mut(i), product.get_mut(i)) {
        *g += bar.0;
        *h += bar.1;
    }
}

fn accumulate_register(
    adjoints: &mut [Dual],
    reg: Reg,
    contribution: Dual,
) {
    if let Some(value) = adjoints.get_mut(reg as usize) {
        value.0 += contribution.0;
        value.1 += contribution.1;
    }
}

#[allow(clippy::single_call_fn)]
fn local_partial(
    inst: &Inst,
    values: &[Dual],
    operand: usize,
    second_order: bool,
) -> Dual {
    // Evaluate powers directly at zero as well as negative bases with a
    // constant integer exponent; logarithmic differentiation creates 0/0.
    match *inst {
        | Inst::PowI(_, n) => {
            let x = values[0];
            let (_, slope, curve) = integer_power_derivatives(x.0, n);
            return Dual(
                slope,
                if second_order {
                    product_unless_zero(x.1, curve)
                } else {
                    0.0
                },
            );
        },
        | Inst::PowF(..) => {
            let (x, b) = (values[0], values[1]);
            let [dx, db, dxx, dxb, dbb] = power_partials(x.0, b.0);
            let (slope, curve) = if operand == 0 {
                (
                    dx,
                    product_unless_zero(x.1, dxx) + product_unless_zero(b.1, dxb),
                )
            } else {
                (
                    db,
                    product_unless_zero(x.1, dxb) + product_unless_zero(b.1, dbb),
                )
            };
            return Dual(slope, if second_order { curve } else { 0.0 });
        },
        | _ => {},
    }
    if second_order {
        let seeded: Vec<_> = values
            .iter()
            .enumerate()
            .map(|(i, x)| HyperDual::new(x.0, if i == operand { 1.0 } else { 0.0 }, x.1, 0.0))
            .collect();
        let result = apply(inst, &seeded);
        Dual(result.first(), result.mixed())
    } else {
        let seeded: Vec<_> = values
            .iter()
            .enumerate()
            .map(|(i, x)| Dual(x.0, if i == operand { 1.0 } else { 0.0 }))
            .collect();
        Dual(apply(inst, &seeded).1, 0.0)
    }
}

// Shared with forward and hyper-dual power rules, including i32::MIN and
// polynomial derivatives at zero (where logarithmic differentiation fails).
pub(super) fn integer_power_derivatives(
    x: f64,
    n: i32,
) -> (f64, f64, f64) {
    let power = |shift: i32| {
        n.checked_sub(shift).map_or_else(
            || x.powf(f64::from(n) - f64::from(shift)),
            |p| int_pow(x, p),
        )
    };
    let slope = if n == 0 {
        0.0
    } else {
        f64::from(n) * power(1)
    };
    let curve = if n == 0 || n == 1 {
        0.0
    } else {
        f64::from(n) * (f64::from(n) - 1.0) * power(2)
    };
    (int_pow(x, n), slope, curve)
}

// These exact exponent identities prevent 0/0; nearby exponents must not use them.
#[allow(clippy::float_cmp)]
pub(super) fn power_partials(
    x: f64,
    b: f64,
) -> [f64; 5] {
    let dx = if b == 0.0 {
        0.0
    } else if b == 1.0 {
        1.0
    } else {
        b * x.powf(b - 1.0)
    };
    let dxx = if b == 0.0 || b == 1.0 {
        0.0
    } else {
        b * (b - 1.0) * x.powf(b - 2.0)
    };
    let value = x.powf(b);
    let log = x.ln();
    let db = if x == 0.0 && b > 0.0 {
        0.0
    } else {
        value * log
    };
    let dbb = if x == 0.0 && b > 0.0 {
        0.0
    } else {
        value * log * log
    };
    let dxb = if x == 0.0 && b > 1.0 {
        0.0
    } else {
        x.powf(b - 1.0) * (1.0 + b * log)
    };
    [dx, db, dxx, dxb, dbb]
}
