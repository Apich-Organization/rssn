//! Serial native matrix visitors: handles stay outside the callback value schema.
use super::super::CallCtx;
use super::super::EvalDyn;
use super::super::Installer;
use super::super::OpError;
use super::super::RuleError;
use super::super::Signature;
use super::super::TypeExpr;
use super::super::Value;
use super::super::common::register;
use super::super::native_callbacks::Callbacks;
use super::super::native_callbacks::count;
use super::super::native_callbacks::function;
use super::super::native_callbacks::typed;
use crate::rules::numerical_abi::native_analysis_jacobian::MatrixObject;
use crate::rules::numerical_abi::native_object_api::Handle;
use crate::rules::numerical_abi::native_object_api::NativeArg;
use crate::rules::numerical_abi::native_object_api::NativeType;
#[derive(Clone, Copy)]
struct Visit {
    mass: bool,
}
impl EvalDyn for Visit {
    fn signature(&self) -> Signature {
        typed(
            vec![
                TypeExpr::Opaque(MatrixObject::KIND.into()),
                TypeExpr::Int,
                function(
                    vec![TypeExpr::Int, TypeExpr::Int, TypeExpr::Float],
                    TypeExpr::Bool,
                ),
            ],
            TypeExpr::Bool,
            &[],
        )
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let run = Callbacks::new(ctx, args, &self.signature(), 2, count(&args[1], 1)?)?;
        let matrix = <Handle<MatrixObject> as NativeArg>::decode(&args[0], 0, ctx)?;
        let visit = |i: usize, j: usize, value: f64| {
            let _ = run.call(
                0,
                &[Value::int(i), Value::int(j), Value::Float(value)],
                |v| {
                    let Value::Bool(ack) = v else {
                        return Err(OpError::Domain(
                            "native visitor must return Bool acknowledgement".into(),
                        ));
                    };
                    Ok(ack)
                },
            );
        };
        matrix.with(ctx, |matrix| {
            if self.mass {
                matrix.mass().for_each(visit);
            } else {
                matrix.jac().for_each(visit);
            }
            Ok(())
        })?;
        run.finish()?;
        Ok(Value::Bool(true))
    }
}
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, path, mass) in [
        ("mass", "crate::kernels::jacobian::Mass<'_>::for_each", true),
        ("jacobian", "crate::kernels::jacobian::Jac::for_each", false),
    ] {
        register(
            i,
            &format!("num_analysis_{name}_for_each"),
            path,
            "Genuine serial native matrix visitor on normalized coefficients. Ordinary native handle is authenticated by its session bank and never enters the callback. Visitor receives exact row/column Int and Float value,returns an ignored Bool acknowledgement. Shared evaluation budget and original source errors stop subsequent source invocations; native iteration order is preserved.",
            3,
            Visit { mass },
        )?;
    }
    Ok(())
}
