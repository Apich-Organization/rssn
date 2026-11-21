//! Reverse accumulation with a directional derivative of every adjoint.
//! Local analytic rules share the dual/hyper-dual interpreter, so new tape
//! instructions cannot silently acquire a different derivative convention.
use super::Dual;
use super::HyperDual;
use super::Inst;
use super::Reg;
use super::Scalar;
use super::Tape;
use super::apply;
use super::int_pow;
use super::with_operands;

struct ReverseWork {
    adjoints: Vec<Dual>,
    product: Vec<f64>,
    local: Vec<Dual>,
}

impl ReverseWork {
    fn new(
        registers: usize,
        inputs: usize,
    ) -> Self {
        Self {
            adjoints: vec![Dual(0.0, 0.0); registers],
            product: vec![0.0; inputs],
            local: Vec::new(),
        }
    }
}

struct OutputDependencies {
    seen: Vec<usize>,
    pending: Vec<usize>,
    descending: Vec<usize>,
}

impl OutputDependencies {
    fn collect(
        &mut self,
        tape: &Tape,
        output: usize,
        row: usize,
    ) -> bool {
        let limit = tape.insts.len() / 4;
        if limit < 16 {
            return false;
        }
        self.pending.clear();
        self.descending.clear();
        self.seen[output] = row;
        self.pending.push(output);
        while let Some(index) = self.pending.pop() {
            if self.descending.len() == limit {
                return false;
            }
            self.descending.push(index);
            with_operands(&tape.insts[index], &tape.args, |operands| {
                for &operand in operands {
                    let operand = operand as usize;
                    if self.seen[operand] != row {
                        self.seen[operand] = row;
                        self.pending.push(operand);
                    }
                }
            });
        }
        self.descending.sort_unstable_by(|a, b| b.cmp(a));
        true
    }
}

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
        let mut work = ReverseWork::new(self.insts.len(), inputs.len());
        for ((&reg, &weight), value) in self.outputs.iter().zip(weights).zip(&out) {
            if weight == 0.0 {
                continue;
            }
            if value.0.is_nan() || weight.is_nan() {
                return invalid();
            }
            let Some(adjoint) = work.adjoints.get_mut(reg as usize) else {
                return invalid();
            };
            adjoint.0 += weight;
        }
        let mut gradient = vec![0.0; inputs.len()];
        self.reverse_sweep(
            &values,
            &mut work,
            &mut gradient,
            direction.is_some(),
            (0..self.insts.len()).rev(),
        );
        (gradient, work.product)
    }

    /// Reuses the checked elementary tape's finite-point value pass for every Jacobian row.
    pub(in crate::backend) fn jacobian_reverse_checked(
        &self,
        registers: Vec<f64>,
        columns: usize,
        cells: usize,
    ) -> Vec<f64> {
        const MIN_ROWS_PER_BLOCK: usize = 8;
        let values: Vec<_> = registers
            .into_iter()
            .map(|value| Dual(value, 0.0))
            .collect();
        let mut data = vec![0.0; cells];
        let fill = |start: usize, rows: &mut [f64]| {
            let mut work = ReverseWork::new(self.insts.len(), columns);
            let mut dependencies = OutputDependencies {
                seen: if self.insts.len() >= 64 {
                    vec![0; self.insts.len()]
                } else {
                    Vec::new()
                },
                pending: Vec::new(),
                descending: Vec::new(),
            };
            for (offset, (row, &output)) in rows
                .chunks_exact_mut(columns)
                .zip(&self.outputs[start..])
                .enumerate()
            {
                let register = output as usize;
                if values.get(register).is_none_or(|value| value.0.is_nan()) {
                    row.fill(f64::NAN);
                    continue;
                }
                work.product.fill(0.0);
                work.adjoints[register].0 = 1.0;
                if dependencies.collect(self, register, offset + 1) {
                    self.reverse_sweep(
                        &values,
                        &mut work,
                        row,
                        false,
                        dependencies.descending.iter().copied(),
                    );
                } else {
                    self.reverse_sweep(&values, &mut work, row, false, (0..self.insts.len()).rev());
                }
            }
        };
        let work = self.outputs.len().saturating_mul(self.insts.len());
        let workers = if work >= 16_384 {
            crate::execution::cpu_threads()
        } else {
            1
        };
        if workers > 1 && self.outputs.len() > MIN_ROWS_PER_BLOCK {
            let rows_per_block = self
                .outputs
                .len()
                .div_ceil(workers.saturating_mul(2))
                .max(MIN_ROWS_PER_BLOCK);
            crate::kernels::par::for_chunks_mut(
                &mut data,
                columns * rows_per_block,
                2,
                |block, rows| {
                    fill(block * rows_per_block, rows);
                },
            );
        } else {
            fill(0, &mut data);
        }
        data
    }

    fn reverse_sweep(
        &self,
        values: &[Dual],
        work: &mut ReverseWork,
        gradient: &mut [f64],
        second_order: bool,
        indices: impl Iterator<Item = usize>,
    ) {
        let ReverseWork {
            adjoints,
            product,
            local,
        } = work;
        for index in indices {
            let inst = &self.insts[index];
            let bar = adjoints[index];
            if bar.0 == 0.0 && bar.1 == 0.0 {
                continue;
            }
            adjoints[index] = Dual(0.0, 0.0);
            match *inst {
                | Inst::Input(i) => accumulate_input(gradient, product, i as usize, bar),
                | Inst::Gather(base, reg) => {
                    let offset = values.get(reg as usize).map_or(f64::NAN, |x| x.0).trunc();
                    if offset.is_finite() && (0.0..=1e18).contains(&offset) {
                        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                        if let Some(i) = (base as usize).checked_add(offset as usize)
                            && i < self.inputs
                        {
                            accumulate_input(gradient, product, i, bar);
                        }
                    }
                },
                | Inst::Const(_) | Inst::Param(_) | Inst::Cmp(..) => {},
                | Inst::Select(c, a, b) => {
                    let condition = values.get(c as usize).map_or(f64::NAN, |x| x.0);
                    let selected = if condition == 0.0 { b } else { a };
                    accumulate_register(adjoints, selected, bar);
                },
                | _ => {
                    with_operands(inst, &self.args, |registers| {
                        local.clear();
                        local.extend(
                            registers.iter().map(|&r| {
                                values.get(r as usize).copied().unwrap_or_else(Dual::nan)
                            }),
                        );
                        for (operand, &reg) in registers.iter().enumerate() {
                            let partial = local_partial(inst, local, operand, second_order);
                            let contribution = Dual(
                                product_unless_zero(bar.0, partial.0),
                                product_unless_zero(bar.1, partial.0)
                                    + product_unless_zero(bar.0, partial.1),
                            );
                            accumulate_register(adjoints, reg, contribution);
                        }
                    });
                },
            }
        }
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
        let result = apply_seeded(inst, values, |i, x| {
            HyperDual::new(x.0, if i == operand { 1.0 } else { 0.0 }, x.1, 0.0)
        });
        Dual(result.first(), result.mixed())
    } else {
        let result = apply_seeded(inst, values, |i, x| {
            Dual(x.0, if i == operand { 1.0 } else { 0.0 })
        });
        Dual(result.1, 0.0)
    }
}

fn apply_seeded<T: Scalar>(
    inst: &Inst,
    values: &[Dual],
    mut seed: impl FnMut(usize, Dual) -> T,
) -> T {
    let mut inline = [T::nan(); 3];
    let seeded = values.iter().copied().enumerate();
    if let Some(inline) = inline.get_mut(..values.len()) {
        for (slot, (index, value)) in inline.iter_mut().zip(seeded) {
            *slot = seed(index, value);
        }
        apply(inst, inline)
    } else {
        let heap: Vec<_> = seeded.map(|(index, value)| seed(index, value)).collect();
        apply(inst, &heap)
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
