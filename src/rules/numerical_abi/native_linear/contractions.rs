//! Genuine scalar callback monomorphizations for general contraction and simplex integration.
use super::super::native_backend::compiled_snapshot;
use super::Array;
use super::ArrayD;
use super::DimExpr;
use super::Dynamic;
use super::Handle;
use super::Input;
use super::Installer;
use super::IxDyn;
use super::LIMIT;
use super::NativeType;
use super::NewObject;
use super::Normalized;
use super::OpError;
use super::Output;
use super::RuleError;
use super::api;
use super::array;
use super::err;
use super::object_binding;
use super::product;
use crate::kernel::CallCtx;
use crate::kernel::signature::Exp;
use crate::kernel::signature::Rational;
use crate::kernels::simplex_quadrature::MAX_DIMENSION;
use crate::kernels::simplex_quadrature::SimplexMap;
object_binding!(rectangular,"num_linear_contract_rectangular","crate::kernels::contractions::contract_rectangular","Native generic contraction specialized to f64. The compiled scalar callback receives row, inner index, column, has_previous, previous in that order; index/presence inputs are dimensionless and previous/output use the explicit reference dimension. Indices are zero based and coherent physical values are normalized only at the callback boundary.",false,(rows:usize,inner:usize,cols:usize,zero:Normalized<f64>,program:Input,reference:Input)->NewObject<Array>,|ctx,rows,inner,cols,zero,program,reference|{
 if rows>512||inner>512||cols>512||rows.checked_mul(cols).is_none_or(|n|n>LIMIT)||rows.checked_mul(inner).and_then(|n|n.checked_mul(cols)).is_none_or(|n|n>1_000_000){return Err(err("rectangular contraction work/allocation budget exceeded"));}
 let (scale,d)=api::reference(&reference.0)?;let budget=rows*inner*cols;let callback=compiled_snapshot(ctx,&program.0,budget.max(1))?;if callback.input_dimensions.len()!=5{return Err(err("contraction callback needs five inputs"));}for unit in &callback.input_dimensions[..4]{api::same_dimension(unit,&DimExpr::Dimless)?;}api::same_dimension(&callback.input_dimensions[4],&d)?;api::same_dimension(&callback.output_dimension,&d)?;
 let values=crate::kernels::contractions::contract_rectangular(rows,inner,cols,||zero.0,|row,index,column,previous:Option<f64>|{callback.evaluate_coherent(&[row as f64,index as f64,column as f64,if previous.is_some(){1.0}else{0.0},previous.unwrap_or(0.0)*scale]).map(|y|y/scale)}).map_err(|e|err(format!("{e:?}")))?;
 Ok(array(ArrayD::from_shape_vec(IxDyn(&[rows,cols]),values).map_err(err)?,scale,d))
});
struct PreparedSimplex {
    map: SimplexMap,
    coordinate_reference: f64,
    coordinate_dimension: DimExpr,
}

impl NativeType for PreparedSimplex {
    const KIND: &'static str = "rssn/native/simplex-map/v1";
    const MAX_BYTES: usize = 8192;

    fn validate(&self) -> Result<(), OpError> {
        if self.coordinate_reference.is_finite() && self.coordinate_reference > 0.0 {
            Ok(())
        } else {
            Err(err("invalid simplex coordinate reference"))
        }
    }
}

impl PreparedSimplex {
    fn new(
        vertices: &[Vec<f64>],
        coordinate_reference: &Input,
    ) -> Result<Self, OpError> {
        let (coordinate_reference, coordinate_dimension) = api::reference(&coordinate_reference.0)?;
        let map = SimplexMap::new(vertices).map_err(|error| err(format!("{error:?}")))?;
        Ok(Self {
            map,
            coordinate_reference,
            coordinate_dimension,
        })
    }

    fn integral(
        &self,
        ctx: &CallCtx<'_>,
        value_reference: &Input,
        program: &Input,
        tolerances: [f64; 2],
        max_evaluations: usize,
    ) -> Result<Dynamic, OpError> {
        let n = self.map.dimension();
        if max_evaluations == 0 || max_evaluations > 10_000_000 {
            return Err(err("simplex evaluation budget invalid"));
        }
        let (y, dy) = api::reference(&value_reference.0)?;
        let callback = compiled_snapshot(ctx, &program.0, max_evaluations)?;
        if callback.input_dimensions.len() != n {
            return Err(err(
                "simplex callback arity differs from coordinate dimension",
            ));
        }
        for input in &callback.input_dimensions {
            api::same_dimension(input, &self.coordinate_dimension)?;
        }
        api::same_dimension(&callback.output_dimension, &dy)?;
        let function = |points: &[f64]| {
            let mut coherent = [0.0; MAX_DIMENSION];
            for (value, point) in coherent.iter_mut().zip(points) {
                *value = point * self.coordinate_reference;
            }
            callback.sample(&coherent[..n]) / y
        };
        let result = self
            .map
            .integrate(function, tolerances[0], tolerances[1], max_evaluations);
        callback.finish()?;
        let result = result.map_err(|error| err(format!("{error:?}")))?;
        let dimension = self.measure_dimension(dy)?;
        let mut scales = [self.coordinate_reference; MAX_DIMENSION + 1];
        scales[n] = y;
        let measure = |value| {
            api::quantity(
                api::scale_product(value, &scales[..=n], &[])?,
                dimension.clone(),
            )
        };
        Ok(Dynamic(api::record(vec![
            ("value", measure(result.value)?),
            ("error", measure(result.error)?),
            ("evaluations", result.evaluations.encode()?),
            ("converged", result.converged.encode()?),
        ])))
    }

    fn measure_dimension(
        &self,
        value: DimExpr,
    ) -> Result<DimExpr, OpError> {
        product(
            value,
            self.coordinate_dimension
                .clone()
                .pow(Exp::Rat(Rational::int(self.map.dimension() as i64))),
        )
    }
}

object_binding!(
    simplex,
    "num_linear_integrate_simplex",
    "crate::kernels::simplex_quadrature::integrate_simplex",
    "Native adaptive simplex integration of a retained compiled scalar field, with explicit independent coordinate and ordinate references. Each vertex is normalized by the coordinate reference; value/error dimensions are ordinate times coordinate^dimension. Native evaluation and dimension limits are preserved with a checked finite callback/work budget.",
    false,
    (
        vertices: Normalized<Vec<Vec<f64>>>,
        coordinate_reference: Input,
        value_reference: Input,
        program: Input,
        absolute_tolerance: Normalized<f64>,
        relative_tolerance: Normalized<f64>,
        max_evaluations: usize,
    ) -> Dynamic,
    |ctx, vertices, coordinate_reference, value_reference, program, absolute_tolerance, relative_tolerance, max_evaluations| {
        PreparedSimplex::new(&vertices.0, &coordinate_reference)?.integral(
            ctx, &value_reference, &program, [absolute_tolerance.0, relative_tolerance.0], max_evaluations,
        )
    }
);
object_binding!(
    simplex_map_new,
    "num_linear_simplex_map_new",
    "crate::kernels::simplex_quadrature::SimplexMap::new",
    "Prepares the affine geometry of an n-dimensional simplex from n+1 rows of normalized coordinates. The positive coordinate reference records the physical unit and scale for volume and repeated integrations.",
    false,
    (vertices: Normalized<Vec<Vec<f64>>>, coordinate_reference: Input) -> NewObject<PreparedSimplex>,
    |_, vertices, coordinate_reference| {
        let plan = PreparedSimplex::new(&vertices.0, &coordinate_reference)?;
        let dimensions = vec![plan.coordinate_dimension.clone()];
        Ok(NewObject::new(plan, dimensions))
    }
);
object_binding!(
    simplex_map_volume,
    "num_linear_simplex_map_volume",
    "crate::kernels::simplex_quadrature::SimplexMap::volume",
    "Returns the unsigned volume of the prepared simplex in coordinate units raised to its spatial dimension.",
    false,
    (simplex: Handle<PreparedSimplex>) -> Dynamic,
    |ctx, simplex| simplex.with(ctx, |simplex| {
        let scales = [simplex.coordinate_reference; MAX_DIMENSION];
        let physical = api::scale_product(simplex.map.volume(), &scales[..simplex.map.dimension()], &[])?;
        Ok(Dynamic(api::quantity(physical, simplex.measure_dimension(DimExpr::Dimless)?)?))
    })
);
object_binding!(
    simplex_map_integrate,
    "num_linear_integrate_simplex_map",
    "crate::kernels::simplex_quadrature::SimplexMap::integrate",
    "Integrates a retained compiled scalar field over a prepared simplex using adaptive cubature. Geometry is reused. Give a positive value reference and normalized tolerances. Value and error have field dimension times coordinate dimension raised to the spatial rank.",
    false,
    (
        simplex: Handle<PreparedSimplex>,
        value_reference: Input,
        program: Input,
        absolute_tolerance: Normalized<f64>,
        relative_tolerance: Normalized<f64>,
        max_evaluations: usize,
    ) -> Dynamic,
    |ctx, simplex, value_reference, program, absolute_tolerance, relative_tolerance, max_evaluations| {
        simplex.with(ctx, |simplex| simplex.integral(
            ctx, &value_reference, &program, [absolute_tolerance.0, relative_tolerance.0], max_evaluations,
        ))
    }
);
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &rectangular::BINDING,
            &simplex::BINDING,
            &simplex_map_new::BINDING,
            &simplex_map_volume::BINDING,
            &simplex_map_integrate::BINDING,
        ],
        "numerical_matrix",
    )
}
