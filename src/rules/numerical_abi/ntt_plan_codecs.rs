// A portable description selects an immutable RSSN plan from the shared cache.
#[derive(Clone, Copy)]
struct PrimeNttPlan {
    modulus: u64,
    length: usize,
}

record!(PrimeNttPlan {
    modulus: u64,
    length: usize,
});

impl Arg for PrimeNttPlan {
    fn ty() -> TypeExpr {
        <Self as Output>::ty()
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let record = native_fields(value, &["modulus", "length"])?;
        Ok(Self {
            modulus: u64::decode(native_field(record, "modulus")?, index)?,
            length: usize::decode(native_field(record, "length")?, index)?,
        })
    }
}

binding!(
    finite_field_ntt_plan, "num_finite_field_ntt_plan", "numerical_vector",
    "Validates an odd prime and a power-of-two transform length dividing p-1. Returns a portable plan description backed by the shared NTT cache.",
    "crate::kernels::finite_field::ntt_plan", DimensionPolicy::Dimensionless,
    (modulus: u64, length: usize) -> Result<super::PrimeNttPlan, OpError>,
    |modulus: u64, length: usize| {
        super::checked_ntt_plan(modulus, length)?;
        Ok(super::PrimeNttPlan { modulus, length })
    }
);

binding!(
    finite_field_ntt_transform, "num_finite_field_ntt_transform", "numerical_vector",
    "Applies the cached forward or inverse number-theoretic transform to ordinary unsigned coefficients. The input must match the plan length. Coefficients are reduced modulo the prime and the inverse is normalized.",
    "crate::kernels::finite_field::ntt_in_place_mont", DimensionPolicy::Dimensionless,
    (values: Vec<u64>, plan: super::PrimeNttPlan, inverse: bool) -> Result<Vec<u64>, OpError>,
    |values: Vec<u64>, plan: super::PrimeNttPlan, inverse: bool| {
        super::apply_ntt_plan(values, plan, inverse)
    }
);

fn checked_ntt_plan(
    modulus: u64,
    length: usize,
) -> Result<std::sync::Arc<finite_field::NttPrimePlan>, OpError> {
    if length >= super::native_object_api::VALUE_LIMIT || !length.is_power_of_two() {
        return Err(OpError::Domain(
            "NTT plan length must be a positive power of two below the native tensor size limit"
                .into(),
        ));
    }
    crate::cancellation::checkpoint()?;
    finite_field::ntt_plan(modulus, length).ok_or_else(|| {
        OpError::Domain("NTT requires an odd prime modulus and a length dividing p-1".into())
    })
}

fn apply_ntt_plan(
    mut values: Vec<u64>,
    plan: PrimeNttPlan,
    inverse: bool,
) -> Result<Vec<u64>, OpError> {
    if values.len() != plan.length {
        return Err(OpError::Domain(
            "NTT input length differs from the plan length".into(),
        ));
    }
    let native = checked_ntt_plan(plan.modulus, plan.length)?;
    let mont = native.mont;
    for value in &mut values {
        *value = mont.to_mont(*value);
    }
    finite_field::ntt_in_place_mont(&mut values, &native, inverse);
    crate::cancellation::checkpoint()?;
    for value in &mut values {
        *value = mont.from_mont(*value);
    }
    Ok(values)
}
