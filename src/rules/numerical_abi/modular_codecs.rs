// Exact modular adapters share the native kernels' bounded integer arithmetic.
use crate::kernels::finite_field;
use crate::kernels::number_theory;

fn require_residue_modulus(modulus: u64) -> Result<(), OpError> {
    if modulus < 2 {
        return Err(OpError::Domain("residue modulus must be at least 2".into()));
    }
    Ok(())
}

fn require_modular_vector_pair(
    left: &[u64],
    right: &[u64],
    modulus: u64,
) -> Result<(), OpError> {
    require_residue_modulus(modulus)?;
    if left.len() != right.len() {
        return Err(OpError::Domain(format!(
            "modular vector lengths differ: {} and {}",
            left.len(),
            right.len()
        )));
    }
    Ok(())
}

fn require_nonzero_modulus(modulus: u128) -> Result<(), OpError> {
    if modulus == 0 {
        return Err(OpError::Domain("modulus must be positive".into()));
    }
    Ok(())
}

binding!(
    finite_field_batch_inverse,
    "num_finite_field_batch_inverse_mod",
    "numerical_vector",
    "`num_finite_field_batch_inverse_mod(values, modulus)` returns every modular inverse using a single inversion and pooled storage. The modulus must be at least 2 and every element must be invertible.",
    "crate::kernels::finite_field::batch_inverse_mod",
    DimensionPolicy::Dimensionless,
    (values: Vec<u64>, modulus: u64) -> Result<Vec<u64>, OpError>,
    |values: Vec<u64>, modulus: u64| {
        super::require_residue_modulus(modulus)?;
        super::finite_field::batch_inverse_mod(&values, modulus).ok_or_else(|| {
            OpError::Domain("every element must be coprime to the modulus".into())
        })
    }
);
binding!(
    finite_field_batch_inverse_skip_zeros,
    "num_finite_field_batch_inverse_skip_zeros",
    "numerical_vector",
    "`num_finite_field_batch_inverse_skip_zeros(values, modulus)` returns each inverse, with 0 in the position of any noninvertible element. The modulus must be at least 2.",
    "crate::kernels::finite_field::batch_inverse_skip_zeros",
    DimensionPolicy::Dimensionless,
    (values: Vec<u64>, modulus: u64) -> Result<Vec<u64>, OpError>,
    |values: Vec<u64>, modulus: u64| {
        super::require_residue_modulus(modulus)?;
        Ok(super::finite_field::batch_inverse_skip_zeros(&values, modulus))
    }
);
binding!(
    finite_field_slice_add,
    "num_finite_field_slice_add_mod",
    "numerical_vector",
    "`num_finite_field_slice_add_mod(a, b, modulus)` adds two equal-length vectors of exact u64 integers, returning residues in [0, modulus). The modulus must be at least 2.",
    "crate::kernels::finite_field::slice_add_mod",
    DimensionPolicy::Dimensionless,
    (a: Vec<u64>, b: Vec<u64>, modulus: u64) -> Result<Vec<u64>, OpError>,
    |a: Vec<u64>, b: Vec<u64>, modulus: u64| {
        super::require_modular_vector_pair(&a, &b, modulus)?;
        let mut out = vec![0; a.len()];
        super::finite_field::slice_add_mod(&a, &b, modulus, &mut out);
        Ok(out)
    }
);
binding!(
    finite_field_slice_sub,
    "num_finite_field_slice_sub_mod",
    "numerical_vector",
    "`num_finite_field_slice_sub_mod(a, b, modulus)` subtracts two equal-length vectors of exact u64 integers, returning residues in [0, modulus). The modulus must be at least 2.",
    "crate::kernels::finite_field::slice_sub_mod",
    DimensionPolicy::Dimensionless,
    (a: Vec<u64>, b: Vec<u64>, modulus: u64) -> Result<Vec<u64>, OpError>,
    |a: Vec<u64>, b: Vec<u64>, modulus: u64| {
        super::require_modular_vector_pair(&a, &b, modulus)?;
        let mut out = vec![0; a.len()];
        super::finite_field::slice_sub_mod(&a, &b, modulus, &mut out);
        Ok(out)
    }
);
binding!(
    finite_field_slice_scale,
    "num_finite_field_slice_scale_mod",
    "numerical_vector",
    "`num_finite_field_slice_scale_mod(values, scalar, modulus)` multiplies an exact u64 vector by a scalar modulo an integer of at least 2.",
    "crate::kernels::finite_field::slice_scale_mod",
    DimensionPolicy::Dimensionless,
    (values: Vec<u64>, scalar: u64, modulus: u64) -> Result<Vec<u64>, OpError>,
    |values: Vec<u64>, scalar: u64, modulus: u64| {
        super::require_residue_modulus(modulus)?;
        let mut out = vec![0; values.len()];
        super::finite_field::slice_scale_mod(&values, scalar, modulus, &mut out);
        Ok(out)
    }
);
binding!(
    finite_field_slice_axpy,
    "num_finite_field_slice_axpy_mod",
    "numerical_vector",
    "`num_finite_field_slice_axpy_mod(alpha, x, y, modulus)` returns (alpha*x + y) reduced elementwise. The vectors must have equal lengths and the modulus must be at least 2.",
    "crate::kernels::finite_field::slice_axpy_mod",
    DimensionPolicy::Dimensionless,
    (alpha: u64, x: Vec<u64>, y: Vec<u64>, modulus: u64) -> Result<Vec<u64>, OpError>,
    |alpha: u64, x: Vec<u64>, mut y: Vec<u64>, modulus: u64| {
        super::require_modular_vector_pair(&x, &y, modulus)?;
        super::finite_field::slice_axpy_mod(alpha, &x, &mut y, modulus);
        Ok(y)
    }
);
binding!(
    finite_field_dot,
    "num_finite_field_dot_mod",
    "numerical_vector",
    "`num_finite_field_dot_mod(a, b, modulus)` returns the exact modular dot product of two equal-length vectors. The modulus must be at least 2.",
    "crate::kernels::finite_field::dot_mod",
    DimensionPolicy::Dimensionless,
    (a: Vec<u64>, b: Vec<u64>, modulus: u64) -> Result<u64, OpError>,
    |a: Vec<u64>, b: Vec<u64>, modulus: u64| {
        super::require_modular_vector_pair(&a, &b, modulus)?;
        Ok(super::finite_field::dot_mod(&a, &b, modulus))
    }
);
binding!(
    number_theory_mul_mod,
    "num_number_theory_mul_mod",
    "numerical_special",
    "`num_number_theory_mul_mod(a, b, modulus)` multiplies two u64 values modulo a positive u64 integer without overflow.",
    "crate::kernels::number_theory::mul_mod",
    DimensionPolicy::Dimensionless,
    (a: u64, b: u64, modulus: u64) -> Result<u64, OpError>,
    |a: u64, b: u64, modulus: u64| {
        super::require_nonzero_modulus(modulus.into())?;
        Ok(super::number_theory::mul_mod(a, b, modulus))
    }
);
binding!(
    number_theory_mul_mod_u128,
    "num_number_theory_mul_mod_u128",
    "numerical_special",
    "`num_number_theory_mul_mod_u128(a, b, modulus)` multiplies two u128 values modulo a positive u128 integer without overflow.",
    "crate::kernels::number_theory::mul_mod_u128",
    DimensionPolicy::Dimensionless,
    (a: u128, b: u128, modulus: u128) -> Result<u128, OpError>,
    |a: u128, b: u128, modulus: u128| {
        super::require_nonzero_modulus(modulus)?;
        Ok(super::number_theory::mul_mod_u128(a, b, modulus))
    }
);
binding!(
    number_theory_mod_pow_u128,
    "num_number_theory_mod_pow_u128",
    "numerical_special",
    "`num_number_theory_mod_pow_u128(base, exponent, modulus)` raises a u128 base to a u128 exponent modulo a positive u128 integer without overflow.",
    "crate::kernels::number_theory::mod_pow_u128",
    DimensionPolicy::Dimensionless,
    (base: u128, exponent: u128, modulus: u128) -> Result<u128, OpError>,
    |base: u128, exponent: u128, modulus: u128| {
        super::require_nonzero_modulus(modulus)?;
        Ok(super::number_theory::mod_pow_u128(base, exponent, modulus))
    }
);
