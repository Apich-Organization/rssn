// A sampled function separates axis coordinates and ordinates at the ABI.
// Unlike a homogeneous two-column tensor, this record preserves both units.
fn sampled_type(rank: u64) -> TypeExpr {
    TypeExpr::Record(vec![
        ("item0".into(), tensor_type(rank)),
        ("item1".into(), tensor_type(rank)),
    ])
}
fn sampled_fields(value: &Value) -> Result<(&Value, &Value), OpError> {
    let Value::Record(record) = value else {
        return Err(OpError::Domain(
            "sampled function requires coordinate/value fields item0 and item1".into(),
        ));
    };
    if record.fields.len() != 2 {
        return Err(OpError::Domain(
            "sampled function requires exactly two fields".into(),
        ));
    }
    let get = |name: &str| {
        let mut fields = record.fields.iter().filter(|(n, _)| n == name);
        let result = fields
            .next()
            .ok_or_else(|| OpError::Domain("missing sampled function field".into()))?;
        if fields.next().is_some() {
            return Err(OpError::Domain("duplicate sampled function field".into()));
        }
        Ok(&result.1)
    };
    Ok((get("item0")?, get("item1")?))
}
fn samples(
    xs: Vec<f64>,
    ys: Vec<f64>,
) -> Result<Vec<(f64, f64)>, OpError> {
    if xs.len() != ys.len()
        || xs.iter().chain(&ys).any(|v| !v.is_finite())
        || xs.windows(2).any(|w| w[0] >= w[1])
    {
        return Err(OpError::Domain("sample coordinates and values must be finite, equally sized, with strictly increasing coordinates".into()));
    }
    Ok(xs.into_iter().zip(ys).collect())
}
impl Arg for Vec<(f64, f64)> {
    fn ty() -> TypeExpr {
        sampled_type(1)
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let (x, y) = sampled_fields(value)?;
        samples(
            <Vec<f64> as Arg>::decode(x, index)?,
            <Vec<f64> as Arg>::decode(y, index)?,
        )
    }
}
impl Output for Vec<(f64, f64)> {
    fn ty() -> TypeExpr {
        <Self as Arg>::ty()
    }

    fn encode(self) -> Result<Value, OpError> {
        let (x, y): (Vec<_>, Vec<_>) = self.into_iter().unzip();
        (x, y).encode()
    }
}
impl Arg for Vec<Vec<(f64, f64)>> {
    fn ty() -> TypeExpr {
        sampled_type(2)
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let (x, y) = sampled_fields(value)?;
        let x = <Vec<Vec<f64>> as Arg>::decode(x, index)?;
        let y = <Vec<Vec<f64>> as Arg>::decode(y, index)?;
        if x.len() != y.len() {
            return Err(OpError::Domain(
                "sample basis coordinates and ordinates must have equal row counts".into(),
            ));
        }
        x.into_iter().zip(y).map(|(x, y)| samples(x, y)).collect()
    }
}
impl Output for Vec<Vec<(f64, f64)>> {
    fn ty() -> TypeExpr {
        <Self as Arg>::ty()
    }

    fn encode(self) -> Result<Value, OpError> {
        let (x, y): (Vec<Vec<f64>>, Vec<Vec<f64>>) =
            self.into_iter().map(|row| row.into_iter().unzip()).unzip();
        (x, y).encode()
    }
}
impl Arg for Vec<u64> {
    fn ty() -> TypeExpr {
        <Vec<usize> as Arg>::ty()
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let Value::Tensor(t) = value else {
            return Err(OpError::Domain("expected integer vector".into()));
        };
        if t.shape != [t.data.len()] {
            return Err(OpError::Domain("malformed integer vector".into()));
        }
        t.data
            .iter()
            .map(|v| <u64 as Arg>::decode(v, index))
            .collect()
    }
}
impl Output for Vec<u64> {
    fn ty() -> TypeExpr {
        <Self as Arg>::ty()
    }

    fn encode(self) -> Result<Value, OpError> {
        Ok(Value::Tensor(TensorValue {
            shape: vec![self.len()],
            data: self.into_iter().map(Value::int).collect(),
            elem: TypeExpr::Int,
        }))
    }
}
