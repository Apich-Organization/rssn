// Source-derived adapters for native coding, topology and normalized chaos.
// These codecs only validate representations and encode native kernel output.
const NATIVE_VALUE_LIMIT: usize = 1_048_576;

integer_codec!(u16);
integer_output!(u16);

fn native_fields<'a>(
    value: &'a Value,
    names: &[&str],
) -> Result<&'a RecordValue, OpError> {
    let Value::Record(record) = value else {
        return Err(OpError::Domain("expected a native numerical record".into()));
    };
    if record.fields.len() != names.len()
        || names.iter().any(|name| {
            record
                .fields
                .iter()
                .filter(|(field, _)| field == name)
                .count()
                != 1
        })
    {
        return Err(OpError::Domain(
            "native numerical record has missing, duplicate or extra fields".into(),
        ));
    }
    Ok(record)
}
fn native_field<'a>(
    record: &'a RecordValue,
    name: &str,
) -> Result<&'a Value, OpError> {
    record
        .fields
        .iter()
        .find(|(field, _)| field == name)
        .map(|(_, value)| value)
        .ok_or_else(|| OpError::Domain(format!("native numerical record lacks {name}")))
}
macro_rules! native_vector_codec {
    ($($element:ty),+) => {$(
        impl Arg for Vec<$element> {
            fn ty() -> TypeExpr {TypeExpr::Tensor{rank:NatExpr::Lit(1),shape:None,elem:Box::new(<$element as Arg>::ty())}}
            #[allow(clippy::indexing_slicing,clippy::arithmetic_side_effects)] // Rank/storage checks precede every row/index operation.
            fn decode(value:&Value,index:usize)->Result<Self,OpError>{
                let Value::Tensor(t)=value else{return Err(OpError::Domain("expected a native rank-one tensor".into()));};
                if t.shape.len()!=1 || t.shape[0]!=t.data.len() || t.data.len()>NATIVE_VALUE_LIMIT {return Err(OpError::Domain("native vector shape or allocation budget is invalid".into()));}
                t.data.iter().map(|v|<$element as Arg>::decode(v,index)).collect()
            }
        }
        impl Output for Vec<$element> {
            fn ty()->TypeExpr{<Self as Arg>::ty()}
            fn encode(self)->Result<Value,OpError>{
                if self.len()>NATIVE_VALUE_LIMIT{return Err(OpError::Domain("native output exceeds allocation budget".into()));}
                Ok(Value::Tensor(TensorValue{shape:vec![self.len()],data:self.into_iter().map(<$element as Output>::encode).collect::<Result<_,_>>()?,elem:<$element as Output>::ty()}))
            }
        }
    )+};
}
native_vector_codec!(u8, u32, bool);
macro_rules! native_integer_grid {
    ($($element:ty),+) => {$(
        impl Arg for Vec<Vec<$element>> {
            fn ty()->TypeExpr{TypeExpr::Tensor{rank:NatExpr::Lit(2),shape:None,elem:Box::new(TypeExpr::Int)}}
            #[allow(clippy::indexing_slicing,clippy::arithmetic_side_effects)] // Rank/storage checks precede every row/index operation.
            fn decode(value:&Value,index:usize)->Result<Self,OpError>{
                let Value::Tensor(t)=value else{return Err(OpError::Domain("expected a native integer matrix".into()));};
                if t.shape.len()!=2 || t.shape.iter().any(|n|*n>NATIVE_VALUE_LIMIT) || t.shape[0].checked_mul(t.shape[1])!=Some(t.data.len()) || t.data.len()>NATIVE_VALUE_LIMIT{return Err(OpError::Domain("integer matrix shape or allocation budget is invalid".into()));}
                (0..t.shape[0]).map(|row|t.data[row*t.shape[1]..(row+1)*t.shape[1]].iter().map(|v|<$element as Arg>::decode(v,index)).collect()).collect()
            }
        }
        impl Output for Vec<Vec<$element>> {
            fn ty()->TypeExpr{<Self as Arg>::ty()}
            fn encode(self)->Result<Value,OpError>{
                let rows=self.len();let columns=self.first().map_or(0,Vec::len);
                if rows.checked_mul(columns).is_none_or(|n|n>NATIVE_VALUE_LIMIT)||self.iter().any(|r|r.len()!=columns){return Err(OpError::Domain("native integer matrix output is ragged or over budget".into()));}
                Ok(Value::Tensor(TensorValue{shape:vec![rows,columns],data:self.into_iter().flatten().map(Value::int).collect(),elem:TypeExpr::Int}))
            }
        }
    )+};
}
native_integer_grid!(u8, u32, i64);

macro_rules! native_tuple_arg {
    ($($ty:ident:$field:ident),+) => {impl<$($ty:Arg),+> Arg for ($($ty,)+){
        fn ty()->TypeExpr{TypeExpr::Record(vec![$((stringify!($field).into(),$ty::ty())),+])}
        fn decode(value:&Value,index:usize)->Result<Self,OpError>{
            let record=native_fields(value,&[$(stringify!($field)),+])?;
            Ok(($($ty::decode(native_field(record,stringify!($field))?,index)?,)+))
        }
    }};
}
native_tuple_arg!(A:item0,B:item1);
native_tuple_arg!(A:item0,B:item1,C:item2);
macro_rules! native_point_vector {
    ($tuple:ty,$width:expr,$from:expr,$into:expr) => {
        impl Arg for Vec<$tuple> {
            fn ty() -> TypeExpr {
                <Vec<[f64; $width]> as Arg>::ty()
            }

            fn decode(
                value: &Value,
                index: usize,
            ) -> Result<Self, OpError> {
                Ok(<Vec<[f64; $width]> as Arg>::decode(value, index)?
                    .into_iter()
                    .map($from)
                    .collect())
            }
        }
        impl Output for Vec<$tuple> {
            fn ty() -> TypeExpr {
                <Vec<[f64; $width]> as Output>::ty()
            }

            fn encode(self) -> Result<Value, OpError> {
                self.into_iter()
                    .map($into)
                    .collect::<Vec<[f64; $width]>>()
                    .encode()
            }
        }
    };
}
native_point_vector!((f64, f64, f64), 3, |[a, b, c]| (a, b, c), |(a, b, c)| {
    [a, b, c]
});

// Ragged integer rows use CSR, preserving simplices and graph components
// without rectangular padding that changes their mathematical meaning.
impl Arg for Vec<Vec<usize>> {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("offsets".into(), <Vec<usize> as Arg>::ty()),
            ("indices".into(), <Vec<usize> as Arg>::ty()),
        ])
    }

    #[allow(clippy::indexing_slicing)] // Monotone offsets are validated against the full storage length below.
    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let record = native_fields(value, &["offsets", "indices"])?;
        let offsets = Vec::<usize>::decode(native_field(record, "offsets")?, index)?;
        let indices = Vec::<usize>::decode(native_field(record, "indices")?, index)?;
        if offsets.is_empty()
            || offsets.len() > NATIVE_VALUE_LIMIT
            || indices.len() > NATIVE_VALUE_LIMIT
            || offsets.first() != Some(&0)
            || offsets.last() != Some(&indices.len())
            || offsets.windows(2).any(|w| w[0] > w[1])
        {
            return Err(OpError::Domain("CSR offsets must start at zero, be monotone and end at the index count within budget".into()));
        }
        Ok(offsets
            .windows(2)
            .map(|w| indices[w[0]..w[1]].to_vec())
            .collect())
    }
}
impl Output for Vec<Vec<usize>> {
    fn ty() -> TypeExpr {
        <Self as Arg>::ty()
    }

    fn encode(self) -> Result<Value, OpError> {
        let mut offsets = vec![0usize];
        let mut indices = Vec::new();
        for row in self {
            if indices
                .len()
                .checked_add(row.len())
                .is_none_or(|n| n > NATIVE_VALUE_LIMIT)
                || offsets.len() >= NATIVE_VALUE_LIMIT
            {
                return Err(OpError::Domain(
                    "CSR output exceeds its allocation budget".into(),
                ));
            }
            indices.extend(row);
            offsets.push(indices.len());
        }
        Ok(Value::Record(RecordValue {
            fields: vec![
                ("offsets".into(), offsets.encode()?),
                ("indices".into(), indices.encode()?),
            ],
        }))
    }
}

impl Arg for crate::kernels::error_correction::PolyGF256 {
    fn ty() -> TypeExpr {
        <Vec<u8> as Arg>::ty()
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        Ok(Self(Vec::<u8>::decode(value, index)?))
    }
}
impl Arg for crate::kernels::robust::Kernel {
    fn ty() -> TypeExpr {
        TypeExpr::Int
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        match u8::decode(value, index)? {
            | 0 => Ok(Self::Gaussian),
            | 1 => Ok(Self::Epanechnikov),
            | _ => {
                Err(OpError::Domain(
                    "KDE kernel must be 0 (Gaussian) or 1 (Epanechnikov)".into(),
                ))
            },
        }
    }
}
impl Arg for crate::kernels::stats::VarianceType {
    fn ty() -> TypeExpr {
        TypeExpr::Bool
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        Ok(if bool::decode(value, index)? {
            Self::Population
        } else {
            Self::Sample
        })
    }
}

#[allow(clippy::arithmetic_side_effects)] // Checked unit scales are exact rationals and cannot overflow.
impl Arg for num_rational::BigRational {
    fn ty() -> TypeExpr {
        TypeExpr::Real
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        match value {
            | Value::Int(n) => Ok(Self::from_integer(n.clone())),
            | Value::Real(crate::graph::number::Number::Int(n)) => {
                Ok(Self::from_integer(n.clone()))
            },
            | Value::Real(crate::graph::number::Number::Rat(n)) => Ok(n.clone()),
            | Value::Quantity(q) => {
                if q.unit.affine
                    || q.unit.offset.is_some()
                    || q.unit.scale.pi != 0
                    || num_traits::Zero::is_zero(&q.unit.scale.denom)
                    || num_traits::Zero::is_zero(&q.unit.scale.numer)
                {
                    return Err(OpError::Domain("exact rational kernels require a finite nonzero rational linear unit scale".into()));
                }
                Ok(Self::decode(&q.value, index)? * q.unit.scale.rational())
            },
            | Value::Dyn(d) => Self::decode(&d.value, index),
            | _ => {
                Err(OpError::Type {
                    expected: Box::new(TypeExpr::Real),
                    got: Box::new(value.type_expr()),
                    arg: index,
                })
            },
        }
    }
}
impl Output for num_rational::BigRational {
    fn ty() -> TypeExpr {
        TypeExpr::Real
    }

    fn encode(self) -> Result<Value, OpError> {
        Ok(Value::rational(self))
    }
}
native_vector_codec!(num_rational::BigRational);

impl Arg for Vec<crate::kernels::fractal_geometry_and_chaos::AffineTransform2D> {
    fn ty() -> TypeExpr {
        <Vec<[f64; 6]> as Arg>::ty()
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        Ok(<Vec<[f64; 6]> as Arg>::decode(value, index)?
            .into_iter()
            .map(|[a, b, c, d, e, f]| {
                crate::kernels::fractal_geometry_and_chaos::AffineTransform2D::new(a, b, c, d, e, f)
            })
            .collect())
    }
}
impl Output for Vec<crate::kernels::fractal_geometry_and_chaos::AffineTransform2D> {
    fn ty() -> TypeExpr {
        <Self as Arg>::ty()
    }

    fn encode(self) -> Result<Value, OpError> {
        self.into_iter()
            .map(|t| [t.a, t.b, t.c, t.d, t.e, t.f])
            .collect::<Vec<_>>()
            .encode()
    }
}
// Simplicial boundaries preserve the signs separately from ragged vertex rows.
impl Output for Vec<(i64, Vec<usize>)> {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("coefficients".into(), <Vec<i64> as Output>::ty()),
            ("simplices".into(), <Vec<Vec<usize>> as Output>::ty()),
        ])
    }

    fn encode(self) -> Result<Value, OpError> {
        let (coefficients, simplices): (Vec<_>, Vec<_>) = self.into_iter().unzip();
        Ok(Value::Record(RecordValue {
            fields: vec![
                ("coefficients".into(), coefficients.encode()?),
                ("simplices".into(), simplices.encode()?),
            ],
        }))
    }
}
// A filtration is one radius and one simplex per row. Neither padding nor
// sentinel vertices are introduced when the simplex dimensions differ.
impl Arg for Vec<(f64, Vec<usize>)> {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("radii".into(), <Vec<f64> as Arg>::ty()),
            ("simplices".into(), <Vec<Vec<usize>> as Arg>::ty()),
        ])
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let record = native_fields(value, &["radii", "simplices"])?;
        let radii = Vec::<f64>::decode(native_field(record, "radii")?, index)?;
        let simplices = Vec::<Vec<usize>>::decode(native_field(record, "simplices")?, index)?;
        if radii.len() != simplices.len() {
            return Err(OpError::Domain(
                "filtration radii and simplices must have equal row counts".into(),
            ));
        }
        Ok(radii.into_iter().zip(simplices).collect())
    }
}
impl Output for Vec<(f64, Vec<Vec<usize>>)> {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("radii".into(), <Vec<f64> as Output>::ty()),
            ("complex_offsets".into(), <Vec<usize> as Output>::ty()),
            ("simplices".into(), <Vec<Vec<usize>> as Output>::ty()),
        ])
    }

    fn encode(self) -> Result<Value, OpError> {
        let mut radii = Vec::new();
        let mut complex_offsets = vec![0usize];
        let mut simplices = Vec::new();
        for (radius, complex) in self {
            radii.push(radius);
            if simplices
                .len()
                .checked_add(complex.len())
                .is_none_or(|n| n > NATIVE_VALUE_LIMIT)
            {
                return Err(OpError::Domain(
                    "filtration output exceeds allocation budget".into(),
                ));
            }
            simplices.extend(complex);
            complex_offsets.push(simplices.len());
        }
        Ok(Value::Record(RecordValue {
            fields: vec![
                ("radii".into(), radii.encode()?),
                ("complex_offsets".into(), complex_offsets.encode()?),
                ("simplices".into(), simplices.encode()?),
            ],
        }))
    }
}
impl Output for Vec<crate::kernels::topology::PersistenceDiagram> {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("dimensions".into(), <Vec<usize> as Output>::ty()),
            ("offsets".into(), <Vec<usize> as Output>::ty()),
            ("births".into(), <Vec<f64> as Output>::ty()),
            ("deaths".into(), <Vec<f64> as Output>::ty()),
            ("essential".into(), <Vec<bool> as Output>::ty()),
        ])
    }

    fn encode(self) -> Result<Value, OpError> {
        let mut dimensions = Vec::new();
        let mut offsets = vec![0usize];
        let mut births = Vec::new();
        let mut deaths = Vec::new();
        let mut essential = Vec::new();
        for diagram in self {
            dimensions.push(diagram.dimension);
            for interval in diagram.intervals {
                if births.len() >= NATIVE_VALUE_LIMIT {
                    return Err(OpError::Domain(
                        "persistence output exceeds allocation budget".into(),
                    ));
                }
                let infinite = crate::kernels::exact_eq(interval.death, f64::INFINITY);
                births.push(interval.birth);
                deaths.push(if infinite {
                    0.0
                } else {
                    interval.death
                });
                essential.push(infinite);
            }
            offsets.push(births.len());
        }
        Ok(Value::Record(RecordValue {
            fields: vec![
                ("dimensions".into(), dimensions.encode()?),
                ("offsets".into(), offsets.encode()?),
                ("births".into(), births.encode()?),
                ("deaths".into(), deaths.encode()?),
                ("essential".into(), essential.encode()?),
            ],
        }))
    }
}

fn native_budget(
    factors: &[usize],
    limit: usize,
) -> Result<usize, OpError> {
    if factors.iter().any(|n| *n > limit) {
        return Err(OpError::Domain(
            "native numerical extent exceeds its finite budget".into(),
        ));
    }
    factors.iter().try_fold(1usize, |total, n| {
        total
            .checked_mul(*n)
            .filter(|v| *v <= limit)
            .ok_or_else(|| {
                OpError::Domain("native numerical work/allocation budget exceeded".into())
            })
    })
}
#[allow(clippy::indexing_slicing)] // The CodecExport arity check runs before family preflight.
fn native_count(
    args: &[Value],
    index: usize,
) -> Result<usize, OpError> {
    usize::decode(&args[index], index)
}
fn native_tensor_shape(value: &Value) -> Result<&[usize], OpError> {
    let Value::Tensor(t) = value else {
        return Err(OpError::Domain("expected a native tensor argument".into()));
    };
    if t.shape.iter().any(|extent| *extent > NATIVE_VALUE_LIMIT)
        || t.shape.iter().try_fold(1usize, |a, b| a.checked_mul(*b)) != Some(t.data.len())
        || t.data.len() > NATIVE_VALUE_LIMIT
    {
        return Err(OpError::Domain(
            "native tensor shape or allocation budget is invalid".into(),
        ));
    }
    Ok(&t.shape)
}
#[allow(clippy::arithmetic_side_effects)] // Simplex length is bounded before shifts.
fn native_complex_budget(complex: &[Vec<usize>]) -> Result<(), OpError> {
    let mut faces = 0usize;
    for simplex in complex {
        if simplex.len() > 16 {
            return Err(OpError::Domain(
                "native topology simplex dimension exceeds finite face budget".into(),
            ));
        }
        faces = faces
            .checked_add(1usize << simplex.len())
            .filter(|n| *n <= 65_536)
            .ok_or_else(|| {
                OpError::Domain("native topology closure exceeds finite face budget".into())
            })?;
    }
    Ok(())
}
#[allow(clippy::arithmetic_side_effects)] // Extents/dimensions and every count product are bounded before incrementing.
fn native_rips_budget(
    points: usize,
    dimension: usize,
) -> Result<(), OpError> {
    if dimension > 8 {
        return Err(OpError::Domain(
            "native topology dimension must not exceed eight".into(),
        ));
    }
    let mut term = 1usize;
    let mut sum = 0usize;
    for size in 1..=(dimension + 2).min(points) {
        term = term
            .checked_mul(points + 1 - size)
            .and_then(|n| n.checked_div(size))
            .ok_or_else(|| OpError::Domain("native topology simplex budget exceeded".into()))?;
        sum = sum
            .checked_add(term)
            .filter(|n| *n <= 4096)
            .ok_or_else(|| OpError::Domain("native topology simplex budget exceeded".into()))?;
    }
    native_budget(&[points, points], NATIVE_VALUE_LIMIT)?;
    Ok(())
}
#[allow(
    clippy::too_many_lines,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects
)] // One preflight dispatch audits family work/shape invariants after exact arity validation.
fn native_scalar_preflight(
    name: &str,
    args: &[Value],
) -> Result<(), OpError> {
    native_discrete_preflight(name, args)?;
    if let Some(operation) = name.strip_prefix("num_fractal_geometry_and_chaos_") {
        if operation.starts_with("generate_")
            && (operation.ends_with("set")
                || matches!(
                    operation,
                    "generate_burning_ship" | "generate_multibrot" | "generate_newton_fractal"
                ))
        {
            let pixels = native_budget(
                &[native_count(args, 0)?, native_count(args, 1)?],
                NATIVE_VALUE_LIMIT,
            )?;
            let iterations = if matches!(operation, "generate_julia_set" | "generate_multibrot") {
                5
            } else {
                4
            };
            native_budget(&[pixels, native_count(args, iterations)?], 50_000_000)?;
        } else if matches!(operation, "mandelbrot_escape_time" | "julia_escape_time") {
            native_budget(&[native_count(args, args.len() - 1)?], 50_000_000)?;
        } else if matches!(
            operation,
            "generate_lorenz_attractor"
                | "generate_lorenz_attractor_custom"
                | "generate_rossler_attractor"
                | "lyapunov_exponent_lorenz"
        ) {
            native_budget(&[native_count(args, 2)?, 3], NATIVE_VALUE_LIMIT)?;
        } else if matches!(operation, "generate_henon_map" | "generate_tinkerbell_map") {
            native_budget(&[native_count(args, 1)?, 2], NATIVE_VALUE_LIMIT)?;
        } else if operation == "logistic_map_iterate" {
            native_budget(&[native_count(args, 2)?], NATIVE_VALUE_LIMIT)?;
        } else if operation == "logistic_bifurcation" {
            let kept = native_count(args, 3)?;
            let discarded = native_count(args, 2)?;
            let steps = kept
                .checked_add(discarded)
                .ok_or_else(|| OpError::Domain("logistic work count overflows".into()))?;
            native_budget(&[native_count(args, 1)?, steps], 50_000_000)?;
            native_budget(&[native_count(args, 1)?, kept, 2], NATIVE_VALUE_LIMIT)?;
        } else if operation == "lyapunov_exponent_logistic" {
            let steps = native_count(args, 2)?
                .checked_add(native_count(args, 3)?)
                .ok_or_else(|| OpError::Domain("logistic work count overflows".into()))?;
            native_budget(&[steps], 50_000_000)?;
        } else if matches!(
            operation,
            "box_counting_dimension" | "correlation_dimension"
        ) {
            let shape = native_tensor_shape(&args[0])?;
            let points = shape.first().copied().unwrap_or(0);
            native_budget(&[points, points], NATIVE_VALUE_LIMIT)?;
            native_budget(&[points, native_count(args, 1)?], 50_000_000)?;
        } else if operation == "orbit_density" {
            native_budget(
                &[native_count(args, 1)?, native_count(args, 2)?],
                NATIVE_VALUE_LIMIT,
            )?;
        } else if operation == "generate_ifs_fractal" {
            let rows = native_tensor_shape(&args[0])?.first().copied().unwrap_or(0);
            let probabilities = Vec::<f64>::decode(&args[1], 1)?;
            let total: f64 = probabilities.iter().sum();
            if rows == 0
                || rows != probabilities.len()
                || probabilities.iter().any(|p| *p < 0.0)
                || !total.is_finite()
                || total <= 0.0
            {
                return Err(OpError::Domain(
                    "IFS requires matching nonnegative probabilities with a positive finite sum"
                        .into(),
                ));
            }
            let work = native_count(args, 3)?
                .checked_add(native_count(args, 4)?)
                .ok_or_else(|| OpError::Domain("IFS work count overflows".into()))?;
            native_budget(&[work], 50_000_000)?;
            native_budget(&[native_count(args, 3)?, 2], NATIVE_VALUE_LIMIT)?;
        }
    }
    if let Some(operation) = name.strip_prefix("num_error_correction_") {
        if let Some(first) = args.first()
            && let Value::Tensor(_) = first
        {
            native_tensor_shape(first)?;
        }
        if matches!(
            operation,
            "reed_solomon_encode"
                | "reed_solomon_decode"
                | "reed_solomon_check"
                | "calculate_syndromes"
        ) {
            let length = native_tensor_shape(&args[0])?.first().copied().unwrap_or(0);
            let parity = native_count(args, 1)?;
            if parity > 255
                || length > 255
                || (operation == "reed_solomon_encode" && length + parity > 255)
            {
                return Err(OpError::Domain(
                    "GF256 Reed-Solomon block/parity counts must fit 255 symbols".into(),
                ));
            }
        }
        if operation == "minimum_distance" {
            let shape = native_tensor_shape(&args[0])?;
            if shape.len() != 2 {
                return Err(OpError::Domain(
                    "codewords must be a rectangular byte matrix".into(),
                ));
            }
            native_budget(&[shape[0], shape[0], shape[1]], 50_000_000)?;
        }
        if operation == "convolutional_encode" {
            native_budget(
                &[
                    native_tensor_shape(&args[0])?.first().copied().unwrap_or(0),
                    2,
                ],
                NATIVE_VALUE_LIMIT,
            )?;
        }
        if matches!(operation, "bch_encode" | "bch_decode") {
            let length = native_tensor_shape(&args[0])?.first().copied().unwrap_or(0);
            let t = native_count(args, 1)?;
            if length > 65_535 || t == 0 || t > 32 {
                return Err(OpError::Domain("BCH requires at most 65535 bits and correction capability 1..=32 within the native work budget".into()));
            }
            let bits = Vec::<u8>::decode(&args[0], 0)?;
            if bits.iter().any(|bit| *bit > 1) {
                return Err(OpError::Domain(
                    "BCH input must contain only binary digits".into(),
                ));
            }
        }
        if matches!(operation, "interleave" | "deinterleave")
            && (native_count(args, 1)? == 0 || native_count(args, 1)? > NATIVE_VALUE_LIMIT)
        {
            return Err(OpError::Domain(
                "interleaving depth must be positive and within the finite work budget".into(),
            ));
        }
        if matches!(operation, "chien_search" | "forney_algorithm") {
            for (index, value) in args.iter().enumerate() {
                if native_tensor_shape(value)?.first().copied().unwrap_or(0) > 255 {
                    return Err(OpError::Domain(format!(
                        "GF256 polynomial argument {index} exceeds 255 coefficients"
                    )));
                }
            }
        }
    }
    if name == "num_fractal_geometry_and_chaos_orbit_entropy" {
        let rows = Vec::<Vec<usize>>::decode(&args[0], 0)?;
        rows.iter().flatten().try_fold(0usize, |sum, count| {
            sum.checked_add(*count)
                .ok_or_else(|| OpError::Domain("orbit histogram total exceeds usize range".into()))
        })?;
    }
    if name.starts_with("num_special_exact_") {
        native_budget(&[native_count(args, 0)?], 4096)?;
    }
    if let Some(operation) = name.strip_prefix("num_topology_") {
        if operation == "euclidean_distance"
            && native_tensor_shape(&args[0])? != native_tensor_shape(&args[1])?
        {
            return Err(OpError::Domain(
                "Euclidean points must have equal extents".into(),
            ));
        }
        if matches!(
            operation,
            "vietoris_rips_complex"
                | "betti_numbers_at_radius"
                | "persistent_homology_rips"
                | "compute_persistence"
                | "vietoris_rips_filtration"
        ) {
            let shape = native_tensor_shape(&args[0])?;
            let dim = if operation == "vietoris_rips_filtration" {
                1
            } else if operation == "compute_persistence" {
                native_count(args, 3)?
            } else {
                native_count(args, 2)?
            };
            native_rips_budget(shape.first().copied().unwrap_or(0), dim)?;
            if numeric(&args[1], 1)? < 0.0 {
                return Err(OpError::Domain("topology radii must be nonnegative".into()));
            }
            if operation == "vietoris_rips_filtration" {
                native_budget(
                    &[native_count(args, 2)?.saturating_add(1), 4096],
                    NATIVE_VALUE_LIMIT,
                )?;
            }
            if operation == "persistent_homology_rips" && native_count(args, 3)? > 4096 {
                return Err(OpError::Domain(
                    "persistence simplex limit exceeds adapter budget".into(),
                ));
            }
        } else if matches!(operation, "grid_complex" | "torus_complex") {
            native_budget(
                &[native_count(args, 0)?, native_count(args, 1)?, 16],
                65_536,
            )?;
        } else if operation == "persistent_homology" {
            let filtration = Vec::<(f64, Vec<usize>)>::decode(&args[0], 0)?;
            if filtration.len() > 4096 || native_count(args, 1)? > 8 {
                return Err(OpError::Domain(
                    "explicit persistence filtration exceeds adapter budget".into(),
                ));
            }
            native_complex_budget(
                &filtration
                    .into_iter()
                    .map(|(_, simplex)| simplex)
                    .collect::<Vec<_>>(),
            )?;
        } else if matches!(operation, "boundary_matrix" | "int_matmul") {
            if operation == "boundary_matrix" {
                let a = Vec::<Vec<usize>>::decode(&args[0], 0)?;
                let b = Vec::<Vec<usize>>::decode(&args[1], 1)?;
                native_complex_budget(&a)?;
                native_complex_budget(&b)?;
                native_budget(&[a.len(), b.len()], NATIVE_VALUE_LIMIT)?;
            } else {
                let a = Vec::<Vec<i64>>::decode(&args[0], 0)?;
                let b = Vec::<Vec<i64>>::decode(&args[1], 1)?;
                let inner = a.first().map_or(0, Vec::len);
                if inner != b.len() {
                    return Err(OpError::Domain(
                        "integer matrix product inner extents differ".into(),
                    ));
                }
                native_budget(&[a.len(), b.first().map_or(0, Vec::len), inner], 50_000_000)?;
                let max_a = a
                    .iter()
                    .flatten()
                    .map(|x| x.unsigned_abs())
                    .max()
                    .unwrap_or(0);
                let max_b = b
                    .iter()
                    .flatten()
                    .map(|x| x.unsigned_abs())
                    .max()
                    .unwrap_or(0);
                if u128::from(max_a)
                    .checked_mul(u128::from(max_b))
                    .and_then(|v| v.checked_mul(inner as u128))
                    .is_none_or(|n| n > i64::MAX as u128)
                {
                    return Err(OpError::Domain(
                        "integer matrix product may exceed i64 accumulation range".into(),
                    ));
                }
            }
        } else if operation == "transpose" {
            let shape = native_tensor_shape(&args[0])?;
            if shape.get(1).copied() != Some(native_count(args, 1)?) {
                return Err(OpError::Domain(
                    "integer matrix transpose column count disagrees with storage".into(),
                ));
            }
        } else if operation == "integer_rank" {
            let shape = native_tensor_shape(&args[0])?;
            if shape.len() != 2 {
                return Err(OpError::Domain("integer rank requires a matrix".into()));
            }
            native_budget(&[shape[0], shape[1], shape[0].min(shape[1])], 50_000_000)?;
        } else if operation == "simplex_boundary" {
            if native_tensor_shape(&args[0])?.first().copied().unwrap_or(0) > 16 {
                return Err(OpError::Domain("simplex exceeds face-work budget".into()));
            }
        } else {
            let complex = Vec::<Vec<usize>>::decode(&args[0], 0)?;
            if operation == "find_connected_components" {
                return Ok(());
            }
            native_complex_budget(&complex)?;
            native_budget(&[complex.len(), complex.len()], 16_777_216)?;
        }
    }
    if name == "num_robust_kde"
        && (numeric(&args[2], 2)? <= 0.0 || native_tensor_shape(&args[0])?.first() == Some(&0))
    {
        return Err(OpError::Domain(
            "KDE requires a positive bandwidth and nonempty sample".into(),
        ));
    }
    if name == "num_descriptive_histogram" {
        let bins = native_count(args, 1)?;
        if bins == 0 || bins > NATIVE_VALUE_LIMIT {
            return Err(OpError::Domain(
                "histogram bin count must fit the adapter allocation budget".into(),
            ));
        }
    }
    Ok(())
}

binding!(native_rs_decode,"num_error_correction_reed_solomon_decode","numerical_vector","Decode/correct a finite GF256 Reed-Solomon block, returning a new corrected block. The source's mutable buffer is local to this adapter; inputs are unchanged.","crate::kernels::error_correction::reed_solomon_decode",DimensionPolicy::Dimensionless,
(data:Vec<u8>,parity:usize)->Result<Vec<u8>,OpError>,|mut data:Vec<u8>,parity:usize|{crate::kernels::error_correction::reed_solomon_decode(&mut data,parity).map_err(OpError::Domain)?;Ok(data)});
binding!(native_median,"num_descriptive_median","numerical_stats","Median with type-7 interpolation; sorting occurs in adapter-owned storage, preserving the caller's vector.","crate::kernels::descriptive::median",DimensionPolicy::Preserve,
(data:Vec<f64>)->Result<f64,OpError>,|mut data:Vec<f64>|crate::kernels::descriptive::median(&mut data).map_err(|e|OpError::Domain(e.to_string())));
binding!(native_median_exact,"num_descriptive_median_exact","numerical_stats","Exact rational median; sorting occurs in local storage and no float conversion is performed.","crate::kernels::descriptive::median_exact",DimensionPolicy::Preserve,
(data:Vec<num_rational::BigRational>)->Result<num_rational::BigRational,OpError>,|mut data:Vec<num_rational::BigRational>|crate::kernels::descriptive::median_exact(&mut data).map_err(|e|OpError::Domain(e.to_string())));
struct NativeSampleGroups(Vec<Vec<f64>>);
impl Arg for NativeSampleGroups {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("offsets".into(), <Vec<usize> as Arg>::ty()),
            ("values".into(), <Vec<f64> as Arg>::ty()),
        ])
    }

    #[allow(clippy::indexing_slicing)] // Monotone offsets are validated against value storage.
    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let record = native_fields(value, &["offsets", "values"])?;
        let offsets = Vec::<usize>::decode(native_field(record, "offsets")?, index)?;
        let values = Vec::<f64>::decode(native_field(record, "values")?, index)?;
        if offsets.len() < 3
            || offsets.len() > NATIVE_VALUE_LIMIT
            || values.len() > NATIVE_VALUE_LIMIT
            || offsets.first() != Some(&0)
            || offsets.last() != Some(&values.len())
            || offsets.windows(2).any(|w| w[0] >= w[1])
        {
            return Err(OpError::Domain(
                "ANOVA requires at least two nonempty groups with valid CSR offsets".into(),
            ));
        }
        Ok(Self(
            offsets
                .windows(2)
                .map(|w| values[w[0]..w[1]].to_vec())
                .collect(),
        ))
    }
}
binding!(native_anova,"num_stats_one_way_anova","numerical_stats","One-way ANOVA on independently sized groups encoded as CSR offsets/values. Input groups have one shared measurement unit; F statistic and p-value are dimensionless.","crate::kernels::stats::one_way_anova",DimensionPolicy::NativeGroupsAnova,
(data:crate::rules::numerical_abi::codecs::NativeSampleGroups)->(f64,f64),|mut data:super::NativeSampleGroups|crate::kernels::stats::one_way_anova(&mut data.0.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>()));
binding!(native_line_fit,"num_stats_simple_linear_regression","numerical_stats","Simple native linear regression on separately measured equal-length predictor and response vectors; slope is response/predictor and intercept is response.","crate::kernels::stats::simple_linear_regression",DimensionPolicy::NativeLineFit,
(x:Vec<f64>,y:Vec<f64>)->Result<(f64,f64),OpError>,|x:Vec<f64>,y:Vec<f64>|{if x.len()!=y.len()||x.len()<2{return Err(OpError::Domain("regression needs equal-length samples with at least two points".into()));}Ok(crate::kernels::stats::simple_linear_regression(&x.into_iter().zip(y).collect::<Vec<_>>()))});

// Point clouds deliberately use homogeneous N x 2 tensors, distinct from the
// sampled-function Vec<(x,y)> codec that requires monotonically increasing x.
binding!(native_fractal_box,"num_fractal_geometry_and_chaos_box_counting_dimension","numerical_vector","Box-counting dimension of a finite N x 2 point tensor; both coordinate axes share units. Scale/work counts are bounded before native execution.","crate::kernels::fractal_geometry_and_chaos::box_counting_dimension",DimensionPolicy::UnitlessResult,
(points:Vec<[f64;2]>,scales:usize)->f64,|points:Vec<[f64;2]>,scales:usize|crate::kernels::fractal_geometry_and_chaos::box_counting_dimension(&points.into_iter().map(|[x,y]|(x,y)).collect::<Vec<_>>(),scales));
binding!(native_fractal_correlation,"num_fractal_geometry_and_chaos_correlation_dimension","numerical_vector","Correlation dimension of a finite N x 2 coordinate tensor. Pair allocation and radius-work counts are bounded.","crate::kernels::fractal_geometry_and_chaos::correlation_dimension",DimensionPolicy::UnitlessResult,
(points:Vec<[f64;2]>,radii:usize)->f64,|points:Vec<[f64;2]>,radii:usize|crate::kernels::fractal_geometry_and_chaos::correlation_dimension(&points.into_iter().map(|[x,y]|(x,y)).collect::<Vec<_>>(),radii));
binding!(native_fractal_henon,"num_fractal_geometry_and_chaos_generate_henon_map","numerical_vector","Native Henon map in normalized dimensionless coordinates. Returns an N x 2 point tensor with a checked finite allocation budget.","crate::kernels::fractal_geometry_and_chaos::generate_henon_map",DimensionPolicy::Dimensionless,
(start:(f64,f64),steps:usize,a:f64,b:f64)->Vec<[f64;2]>,|start:(f64,f64),steps:usize,a:f64,b:f64|crate::kernels::fractal_geometry_and_chaos::generate_henon_map(start,steps,a,b).into_iter().map(|(x,y)|[x,y]).collect());
binding!(native_fractal_tinkerbell,"num_fractal_geometry_and_chaos_generate_tinkerbell_map","numerical_vector","Native Tinkerbell map in normalized dimensionless coordinates. Returns an N x 2 tensor under a finite path budget.","crate::kernels::fractal_geometry_and_chaos::generate_tinkerbell_map",DimensionPolicy::Dimensionless,
(start:(f64,f64),steps:usize,a:f64,b:f64,c:f64,d:f64)->Vec<[f64;2]>,|start:(f64,f64),steps:usize,a:f64,b:f64,c:f64,d:f64|crate::kernels::fractal_geometry_and_chaos::generate_tinkerbell_map(start,steps,a,b,c,d).into_iter().map(|(x,y)|[x,y]).collect());
binding!(native_fractal_bifurcation,"num_fractal_geometry_and_chaos_logistic_bifurcation","numerical_vector","Logistic bifurcation samples in normalized coordinates; transient and retained iterations share a finite work budget. Returns columns r and x.","crate::kernels::fractal_geometry_and_chaos::logistic_bifurcation",DimensionPolicy::Dimensionless,
(range:(f64,f64),parameters:usize,transient:usize,points:usize,x0:f64)->Vec<[f64;2]>,|range:(f64,f64),parameters:usize,transient:usize,points:usize,x0:f64|crate::kernels::fractal_geometry_and_chaos::logistic_bifurcation(range,parameters,transient,points,x0).into_iter().map(|(r,x)|[r,x]).collect());
binding!(native_fractal_ifs,"num_fractal_geometry_and_chaos_generate_ifs_fractal","numerical_vector","Native deterministic seeded affine chaos game in normalized coordinates. Affine transforms are N x 6 tensors (a,b,c,d,e,f), probabilities must match and have positive finite sum. Returns an N x 2 point tensor.","crate::kernels::fractal_geometry_and_chaos::generate_ifs_fractal",DimensionPolicy::Dimensionless,
(transforms:Vec<crate::kernels::fractal_geometry_and_chaos::AffineTransform2D>,probabilities:Vec<f64>,start:(f64,f64),points:usize,skip:usize)->Vec<[f64;2]>,|transforms:Vec<crate::kernels::fractal_geometry_and_chaos::AffineTransform2D>,probabilities:Vec<f64>,start:(f64,f64),points:usize,skip:usize|crate::kernels::fractal_geometry_and_chaos::generate_ifs_fractal(&transforms,&probabilities,start,points,skip).into_iter().map(|(x,y)|[x,y]).collect());
binding!(native_fractal_density,"num_fractal_geometry_and_chaos_orbit_density","numerical_vector","Normalized phase-space point histogram. N x 2 input points and bounded bin counts produce CSR integer rows (offsets/indices); entries are counts, preserving empty rows without padding.","crate::kernels::fractal_geometry_and_chaos::orbit_density",DimensionPolicy::Dimensionless,
(points:Vec<[f64;2]>,x_bins:usize,y_bins:usize,x_range:(f64,f64),y_range:(f64,f64))->Vec<Vec<usize>>,|points:Vec<[f64;2]>,x_bins:usize,y_bins:usize,x_range:(f64,f64),y_range:(f64,f64)|crate::kernels::fractal_geometry_and_chaos::orbit_density(&points.into_iter().map(|[x,y]|(x,y)).collect::<Vec<_>>(),x_bins,y_bins,x_range,y_range));

/// Native source functions requiring local mutable storage or separated units.
pub const NATIVE_SCALAR_EXPORTS: &[&CodecExport] = &[
    &native_rs_decode::EXPORT,
    &native_median::EXPORT,
    &native_median_exact::EXPORT,
    &native_anova::EXPORT,
    &native_line_fit::EXPORT,
    &native_fractal_box::EXPORT,
    &native_fractal_correlation::EXPORT,
    &native_fractal_henon::EXPORT,
    &native_fractal_tinkerbell::EXPORT,
    &native_fractal_bifurcation::EXPORT,
    &native_fractal_ifs::EXPORT,
    &native_fractal_density::EXPORT,
];

include!("native_discrete_codecs.rs");
