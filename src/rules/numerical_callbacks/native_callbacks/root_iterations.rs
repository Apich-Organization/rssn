//! Dimensioned numerical root iterations in a normalized coordinate frame.
#![allow(
    clippy::arithmetic_side_effects,
    reason = "operation layouts have fixed arity and the physical frame checks finite floating-point values"
)]

use super::CallCtx;
use super::Callbacks;
use super::Complex64;
use super::DimExpr;
use super::EvalDyn;
use super::Installer;
use super::OpError;
use super::RuleError;
use super::Signature;
use super::TypeExpr;
use super::Value;
use super::arity;
use super::complex_type;
use super::count;
use super::decode;
use super::dimension;
use super::domain;
use super::equal;
use super::finite;
use super::function;
use super::native;
use super::ratio;
use super::real;
use super::record;
use super::register;
use super::scale;
use super::typed;
use super::value;
use super::var;
use crate::kernel::signature::NatExpr;
use crate::kernel::value::TensorValue;

const MAX_ITERATIONS: usize = 4096;

#[derive(Clone, Copy)]
enum Method {
    Muller,
    Halley,
    Homotopy,
    DelvesLyness,
}

impl Method {
    const fn seeds(self) -> usize {
        if matches!(self, Self::Muller) {
            3
        } else {
            1
        }
    }

    const fn start(self) -> usize {
        match self {
            | Self::Muller => 8,
            | Self::Halley => 6,
            | Self::Homotopy => 7,
            | Self::DelvesLyness => 5,
        }
    }

    const fn derivative_orders(self) -> &'static [usize] {
        match self {
            | Self::Muller => &[0],
            | Self::Halley => &[0, 1, 2],
            | Self::Homotopy => &[0, 1, 0],
            | Self::DelvesLyness => &[0, 1],
        }
    }

    fn run(
        self,
        frame: &Frame<'_, '_>,
        seeds: &[Complex64],
        steps: usize,
        corrector: usize,
        tolerance: f64,
    ) -> Result<Vec<Complex64>, String> {
        let nan = Complex64::new(f64::NAN, f64::NAN);
        let first = seeds.first().copied().unwrap_or(nan);
        let f = |z| frame.sample(0, z, None);
        let df = |z| frame.sample(1, z, None);
        let ddf = |z| frame.sample(2, z, None);
        match self {
            | Self::Muller => {
                let [z0, z1, z2] = seeds else {
                    return Err("Muller needs exactly three guesses".into());
                };
                native::muller_root(f, *z0, *z1, *z2, steps, tolerance).map(|z| vec![z])
            },
            | Self::Halley => {
                native::halley_root(f, df, ddf, first, steps, tolerance).map(|z| vec![z])
            },
            | Self::Homotopy => {
                native::homotopy_track_1d(
                    |z, t| frame.sample(0, z, Some(t)),
                    |z, t| frame.sample(1, z, Some(t)),
                    |z, t| frame.sample(2, z, Some(t)),
                    first,
                    steps,
                    corrector,
                    tolerance,
                )
                .map(|z| vec![z])
            },
            | Self::DelvesLyness => native::delves_lyness_roots(f, df, first, 1.0, steps),
        }
    }
}

struct Frame<'a, 'ctx> {
    ctx: &'a CallCtx<'ctx>,
    calls: Callbacks<'a, 'ctx>,
    axis: DimExpr,
    origin: Complex64,
    step: f64,
    outputs: Vec<(DimExpr, f64)>,
}

impl<'a, 'ctx> Frame<'a, 'ctx> {
    #[allow(
        clippy::single_call_fn,
        reason = "the frame constructor authenticates the callback dimensions before evaluation"
    )]
    fn new(
        method: Method,
        args: &[Value],
        ctx: &'a CallCtx<'ctx>,
        signature: &Signature,
        budget: usize,
    ) -> Result<Self, OpError> {
        let calls = Callbacks::new(ctx, args, signature, method.start(), budget)?;
        let base = calls
            .signatures
            .first()
            .ok_or_else(|| domain("complex root callback missing"))?;
        let axis = dimension(
            base.params
                .first()
                .ok_or_else(|| domain("complex root axis missing"))?,
        )?;
        let amplitude = dimension(&base.ret)?;
        let seeds = method.seeds();
        let origin = decode(
            args.get(seeds - 1)
                .ok_or_else(|| domain("complex root seed missing"))?,
            &axis,
            seeds - 1,
            ctx,
        )?;
        let step = scale(
            args.get(seeds)
                .ok_or_else(|| domain("coordinate scale missing"))?,
            &axis,
            seeds,
            ctx,
        )?;
        let base_scale = scale(
            args.get(seeds + 1)
                .ok_or_else(|| domain("value scale missing"))?,
            &amplitude,
            seeds + 1,
            ctx,
        )?;

        let mut outputs = Vec::with_capacity(calls.signatures.len());
        for (signature, order) in calls.signatures.iter().zip(method.derivative_orders()) {
            let mut output = amplitude.clone();
            let mut output_scale = base_scale;
            for _ in 0..*order {
                output = ratio(&output, &axis)?;
                output_scale = finite(output_scale / step)?;
            }
            if output_scale <= 0.0 {
                return Err(domain("native callback reference scale underflowed"));
            }
            equal(&dimension(&signature.ret)?, &output)?;
            outputs.push((output, output_scale));
        }
        Ok(Self {
            ctx,
            calls,
            axis,
            origin,
            step,
            outputs,
        })
    }

    fn physical(
        &self,
        point: Complex64,
    ) -> Complex64 {
        self.origin + self.step * point
    }

    fn sample(
        &self,
        index: usize,
        point: Complex64,
        time: Option<f64>,
    ) -> Complex64 {
        let nan = Complex64::new(f64::NAN, f64::NAN);
        let z = self.physical(point);
        if !z.is_finite() {
            return nan;
        }
        let Some((dimension, scale)) = self.outputs.get(index) else {
            return nan;
        };
        let z = value(z, self.axis.clone());
        // At most two arguments, so the sample path needs no heap allocation.
        let arguments = [z, Value::Float(time.unwrap_or_default())];
        let length = if time.is_some() { 2 } else { 1 };
        let Some(arguments) = arguments.get(..length) else {
            return nan;
        };
        self.calls
            .call(index, arguments, |v| {
                let z = decode(&v, dimension, index, self.ctx)? / *scale;
                Ok(Complex64::new(finite(z.re)?, finite(z.im)?))
            })
            .unwrap_or(nan)
    }
}

impl EvalDyn for Method {
    fn signature(&self) -> Signature {
        let axis = var("D");
        let amplitude = var("E");
        let point = complex_type(axis.clone());
        let mut params = vec![point.clone(); self.seeds()];
        params.extend([real(axis.clone()), real(amplitude.clone())]);
        if !matches!(self, Self::DelvesLyness) {
            params.push(TypeExpr::Float); // normalized residual tolerance
        }
        params.push(TypeExpr::Int); // iterations, nominal homotopy steps, or circle samples
        if matches!(self, Self::Homotopy) {
            params.push(TypeExpr::Int); // corrector iteration bound
        }
        params.push(TypeExpr::Int); // joint callback evaluation budget
        for order in self.derivative_orders() {
            let output = if *order == 0 {
                amplitude.clone()
            } else {
                let mut factors = vec![amplitude.clone()];
                factors.extend((0..*order).map(|_| DimExpr::Inverse(Box::new(axis.clone()))));
                DimExpr::Mul(factors)
            };
            let mut callback = vec![point.clone()];
            if matches!(self, Self::Homotopy) {
                callback.push(TypeExpr::Float);
            }
            params.push(function(callback, complex_type(output)));
        }
        let result = if matches!(self, Self::DelvesLyness) {
            (
                "roots",
                TypeExpr::Tensor {
                    rank: NatExpr::Lit(1),
                    shape: None,
                    elem: Box::new(point),
                },
            )
        } else {
            ("root", point)
        };
        typed(
            params,
            TypeExpr::Record(vec![
                (result.0.into(), result.1),
                ("evaluations".into(), TypeExpr::Int),
            ]),
            &["D", "E"],
        )
    }

    #[allow(
        clippy::indexing_slicing,
        reason = "method arity is checked first, so the fixed position offsets are valid"
    )]
    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let signature = self.signature();
        arity(args, signature.params.len())?;
        let budget = count(&args[self.start() - 1], self.start() - 1)?;
        let steps_at = self.seeds()
            + if matches!(self, Self::DelvesLyness) {
                2
            } else {
                3
            };
        let steps = count(&args[steps_at], steps_at)?;
        let corrector = if matches!(self, Self::Homotopy) {
            count(&args[steps_at + 1], steps_at + 1)?
        } else {
            0
        };
        if steps > MAX_ITERATIONS
            || corrector > MAX_ITERATIONS
            || (matches!(self, Self::DelvesLyness) && steps.saturating_mul(2) > budget)
        {
            return Err(OpError::BudgetExceeded);
        }
        let tolerance = if matches!(self, Self::DelvesLyness) {
            0.0
        } else {
            let index = self.seeds() + 2;
            super::super::positive(&args[index], index, ctx)?
        };
        let frame = Frame::new(*self, args, ctx, &signature, budget)?;
        let seeds = args[..self.seeds()]
            .iter()
            .enumerate()
            .map(|(index, seed)| {
                let z = (decode(seed, &frame.axis, index, ctx)? - frame.origin) / frame.step;
                Ok(Complex64::new(finite(z.re)?, finite(z.im)?))
            })
            .collect::<Result<Vec<_>, OpError>>()?;
        let result = self.run(&frame, &seeds, steps, corrector, tolerance);
        let evaluations = frame.calls.used.get();
        let roots: Result<_, OpError> = result.map_err(OpError::NotReduced).and_then(|roots| {
            roots
                .into_iter()
                .map(|point| {
                    let z = frame.physical(point);
                    finite(z.re)?;
                    finite(z.im)?;
                    Ok(value(z, frame.axis.clone()))
                })
                .collect::<Result<Vec<_>, OpError>>()
        });
        frame.calls.finish()?;
        crate::cancellation::checkpoint()?;
        let mut roots = roots?;
        let (label, payload) = if matches!(self, Self::DelvesLyness) {
            (
                "roots",
                Value::Tensor(TensorValue::of(
                    vec![roots.len()],
                    roots,
                    complex_type(frame.axis),
                )),
            )
        } else {
            (
                "root",
                roots
                    .pop()
                    .ok_or_else(|| domain("native iteration returned no root"))?,
            )
        };
        Ok(record(vec![
            (label, payload),
            ("evaluations", Value::int(evaluations)),
        ]))
    }
}

#[allow(
    clippy::single_call_fn,
    reason = "registration entrypoint for measured root iteration methods"
)]
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, path, docs, method) in [
        (
            "num_complex_root_muller",
            "crate::kernels::complex_roots::muller_root",
            "Muller root using three measured guesses, coordinate scale, residual scale, normalized residual tolerance, iteration limit, callback budget and one analytic callback.",
            Method::Muller,
        ),
        (
            "num_complex_root_halley",
            "crate::kernels::complex_roots::halley_root",
            "Halley root using a measured seed, coordinate and value scales, normalized residual tolerance, iteration limit, callback budget, f, f' and f''. The derivative units follow the axis dimension.",
            Method::Halley,
        ),
        (
            "num_complex_root_homotopy",
            "crate::kernels::complex_roots::homotopy_track_1d",
            "Follows a measured complex zero from t=0 to t=1. Takes seed, coordinate and value scales, normalized residual tolerance, nominal step count, correction limit, callback budget, H, dH/dz and dH/dt.",
            Method::Homotopy,
        ),
        (
            "num_complex_root_delves_lyness",
            "crate::kernels::complex_roots::delves_lyness_roots",
            "Approximates zeros of a pole-free function inside a sampled circle. Takes its measured center and radius, residual reference scale, sample count, callback budget, f and f'.",
            Method::DelvesLyness,
        ),
    ] {
        register(
            i,
            name,
            path,
            docs,
            (method.start() + method.derivative_orders().len()) as u8,
            method,
        )?;
    }
    Ok(())
}
