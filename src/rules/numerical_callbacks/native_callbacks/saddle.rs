//! Checked saddle-point callbacks with physical coordinate and amplitude units.
#![allow(
    clippy::arithmetic_side_effects,
    reason = "the fixed arities bound index arithmetic and numeric scale operations are checked"
)]

use super::CallCtx;
use super::Callbacks;
use super::DimExpr;
use super::EvalDyn;
use super::Installer;
use super::OpError;
use super::RuleError;
use super::Signature;
use super::TypeExpr;
use super::Value;
use super::arity;
use super::complex::complex_type;
use super::complex::decode;
use super::complex::dimension;
use super::complex::value;
use super::count;
use super::domain;
use super::finite;
use super::function;
use super::positive;
use super::real;
use super::record;
use super::register;
use super::scalar;
use super::scale;
use super::typed;
use super::var;
use crate::kernel::value::VariantValue;
use crate::kernels::steepest_descent as native;
use crate::rules::numerical_abi::codecs::Arg;
use num_complex::Complex64 as C;

const NAN: C = C::new(f64::NAN, f64::NAN);

#[derive(Clone, Copy)]
enum Method {
    Find,
    Leading,
    Order2,
    Integrate,
}

impl Method {
    const fn callback_start(self) -> usize {
        if matches!(self, Self::Integrate) {
            6
        } else {
            4
        }
    }

    const fn callback_count(self) -> usize {
        match self {
            | Self::Find => 3,
            | Self::Leading => 1,
            | Self::Order2 | Self::Integrate => 2,
        }
    }
}

fn phase_derivative(
    axis: &DimExpr,
    order: usize,
) -> DimExpr {
    DimExpr::Mul(
        (0..order)
            .map(|_| DimExpr::Inverse(Box::new(axis.clone())))
            .collect(),
    )
}

#[allow(
    clippy::single_call_fn,
    reason = "the exact structural record is shared by both saddle input and output contracts"
)]
fn saddle_schema(axis: &DimExpr) -> TypeExpr {
    TypeExpr::Record(vec![
        ("z0".into(), complex_type(axis.clone())),
        ("phase".into(), TypeExpr::Complex),
        ("curvature".into(), complex_type(phase_derivative(axis, 2))),
        ("angle".into(), TypeExpr::Float),
        ("degeneracy".into(), TypeExpr::Int),
    ])
}

struct Frame<'a, 'ctx> {
    ctx: &'a CallCtx<'ctx>,
    callbacks: Callbacks<'a, 'ctx>,
    axis: DimExpr,
    origin: C,
    width: f64,
}

impl Frame<'_, '_> {
    fn sample(
        &self,
        index: usize,
        coordinate: C,
        result_dimension: &DimExpr,
    ) -> C {
        let point = self.origin + self.width * coordinate;
        if !point.is_finite() {
            return NAN;
        }
        self.callbacks
            .call(index, &[value(point, self.axis.clone())], |result| {
                let sampled: C = decode(&result, result_dimension, index, self.ctx)?;
                Ok(C::new(finite(sampled.re)?, finite(sampled.im)?))
            })
            .unwrap_or(NAN)
    }

    #[allow(
        clippy::single_call_fn,
        reason = "returns a typed option of the complete physical saddle record"
    )]
    fn encode_saddle(
        &self,
        saddle: Option<native::SaddlePoint>,
    ) -> Result<Value, OpError> {
        let Some(saddle) = saddle else {
            return Ok(Value::Variant(VariantValue {
                tag: "None".into(),
                payload: Vec::new(),
            }));
        };
        let point = self.origin + self.width * saddle.z0;
        let curvature = saddle.phi_d2 / self.width / self.width;
        for component in [point.re, point.im, curvature.re, curvature.im] {
            finite(component)?;
        }
        let result = record(vec![
            ("z0", value(point, self.axis.clone())),
            ("phase", value(saddle.phi_val, DimExpr::Dimless)),
            (
                "curvature",
                value(curvature, phase_derivative(&self.axis, 2)),
            ),
            ("angle", Value::Float(saddle.steepest_angle)),
            ("degeneracy", Value::int(saddle.degeneracy)),
        ]);
        Ok(Value::Variant(VariantValue {
            tag: "Some".into(),
            payload: vec![result],
        }))
    }
}

fn saddle_field<'a>(
    source: &'a Value,
    field: &str,
) -> Result<&'a Value, OpError> {
    let Value::Record(fields) = source else {
        return Err(domain("expected a saddle point record"));
    };
    fields
        .get(field)
        .ok_or_else(|| domain(format!("saddle record lacks {field}")))
}

impl EvalDyn for Method {
    fn signature(&self) -> Signature {
        let d = var("D");
        let a = var("A");
        let product = DimExpr::Mul(vec![a.clone(), d.clone()]);
        let point = complex_type(d.clone());
        let saddle = saddle_schema(&d);
        let mut params = vec![
            if matches!(self, Self::Find) {
                point.clone()
            } else {
                saddle.clone()
            },
            real(d.clone()),
            TypeExpr::Float,
        ];
        if matches!(self, Self::Integrate) {
            params.extend([real(d.clone()), real(product.clone())]);
        }
        params.push(TypeExpr::Int);
        if !matches!(self, Self::Leading) {
            params.push(function(vec![point.clone()], TypeExpr::Complex));
        }
        if matches!(self, Self::Find) {
            params.extend([
                function(vec![point.clone()], complex_type(phase_derivative(&d, 1))),
                function(vec![point], complex_type(phase_derivative(&d, 2))),
            ]);
            typed(
                params,
                TypeExpr::Variant(vec![("None".into(), vec![]), ("Some".into(), vec![saddle])]),
                &["D"],
            )
        } else {
            params.push(function(vec![point], complex_type(a)));
            typed(params, complex_type(product), &["D", "A"])
        }
    }

    #[allow(
        clippy::indexing_slicing,
        clippy::too_many_lines,
        reason = "the four fixed method arities share one typed and budgeted coordinate frame"
    )]
    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let signature = self.signature();
        arity(args, signature.params.len())?;
        let start = self.callback_start();
        let budget = count(&args[start - 1], start - 1)?;
        let callbacks = Callbacks::new(ctx, args, &signature, start, budget)?;
        let axis = dimension(&callbacks.signatures[0].params[0])?;
        let width = scale(&args[1], &axis, 1, ctx)?;
        let control = positive(&args[2], 2, ctx)?;
        let source = &args[0];
        let origin = decode(
            if matches!(self, Self::Find) {
                source
            } else {
                saddle_field(source, "z0")?
            },
            &axis,
            0,
            ctx,
        )?;
        let frame = Frame {
            ctx,
            callbacks,
            axis,
            origin,
            width,
        };
        let phase = |z| frame.sample(0, z, &DimExpr::Dimless);
        if matches!(self, Self::Find) {
            let derivative = phase_derivative(&frame.axis, 1);
            let curvature = phase_derivative(&frame.axis, 2);
            let saddle = native::find_saddle_point(
                phase,
                |z| frame.sample(1, z, &derivative) * width,
                |z| frame.sample(2, z, &curvature) * width * width,
                C::new(0.0, 0.0),
                control,
            );
            let value = frame.encode_saddle(saddle);
            frame.callbacks.finish()?;
            crate::cancellation::checkpoint()?;
            return value;
        }

        let phi_val = decode(saddle_field(source, "phase")?, &DimExpr::Dimless, 0, ctx)?;
        let phi_d2: C = decode(
            saddle_field(source, "curvature")?,
            &phase_derivative(&frame.axis, 2),
            0,
            ctx,
        )?;
        let saddle = native::SaddlePoint {
            z0: C::new(0.0, 0.0),
            phi_val,
            phi_d2: phi_d2 * width * width,
            steepest_angle: scalar(saddle_field(source, "angle")?, &DimExpr::Dimless, 0, ctx)?,
            degeneracy: usize::decode(saddle_field(source, "degeneracy")?, 0)?,
        };
        let amplitude_index = self.callback_count() - 1;
        let amplitude_dim = dimension(&frame.callbacks.signatures[amplitude_index].ret)?;
        let amplitude = |z| frame.sample(amplitude_index, z, &amplitude_dim);
        let integral_dim = DimExpr::Mul(vec![amplitude_dim.clone(), frame.axis.clone()]);
        let approximation = match self {
            | Self::Leading => native::saddle_point_leading(amplitude, &saddle, control),
            | Self::Order2 => {
                native::saddle_point_asymptotic_order2(phase, amplitude, &saddle, control)
            },
            | Self::Integrate => {
                let half = finite(scale(&args[3], &frame.axis, 3, ctx)? / width)?;
                let tolerance = finite(scale(&args[4], &integral_dim, 4, ctx)? / width)?;
                native::steepest_descent_integrate(
                    phase, amplitude, &saddle, control, half, tolerance,
                )
            },
            | Self::Find => return Err(domain("saddle search does not integrate")),
        };
        frame.callbacks.finish()?;
        crate::cancellation::checkpoint()?;
        let physical = approximation * width;
        if !physical.is_finite() {
            return Err(OpError::NotReduced(
                "saddle approximation or contour integration did not converge".into(),
            ));
        }
        Ok(value(physical, integral_dim))
    }
}

#[allow(
    clippy::single_call_fn,
    reason = "registration entrypoint for measured saddle-point operators"
)]
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, rust_path, docs, method) in [
        (
            "num_saddle_find",
            "crate::kernels::steepest_descent::find_saddle_point",
            "Finds a slope zero from a measured complex guess, coordinate reference, normalized slope tolerance, callback budget, phase and its first two derivatives. Phase is dimensionless.",
            Method::Find,
        ),
        (
            "num_saddle_leading",
            "crate::kernels::steepest_descent::saddle_point_leading",
            "Leading oriented Gaussian saddle term. Takes the saddle, coordinate reference, positive lambda, callback budget and analytic amplitude.",
            Method::Leading,
        ),
        (
            "num_saddle_order2",
            "crate::kernels::steepest_descent::saddle_point_asymptotic_order2",
            "Oriented Gaussian approximation with the first inverse-lambda correction. The phase and amplitude must be analytic on the sampling disk of radius one tenth of the coordinate reference.",
            Method::Order2,
        ),
        (
            "num_saddle_integrate",
            "crate::kernels::steepest_descent::steepest_descent_integrate",
            "Numeric finite path integral from the negative to positive branch of an oriented simple saddle. Takes saddle, coordinate reference, positive lambda, measured half length, absolute tolerance in the result units, callback budget, phase and amplitude.",
            Method::Integrate,
        ),
    ] {
        register(
            i,
            name,
            rust_path,
            docs,
            (method.callback_start() + method.callback_count()) as u8,
            method,
        )?;
    }
    Ok(())
}
