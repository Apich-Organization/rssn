// A portable Euclidean disk is the homogeneous measured tensor [real, imag, radius].
pub(in crate::rules) struct CheckedComplexDisk(
    pub(in crate::rules) crate::kernels::complex_ball::ComplexBall,
);

impl Arg for CheckedComplexDisk {
    fn ty() -> TypeExpr {
        <[f64; 3] as Arg>::ty()
    }

    fn decode(
        value: &Value,
        at: usize,
    ) -> Result<Self, OpError> {
        decode_disk_cell(<[f64; 3] as Arg>::decode(value, at)?).map(Self)
    }
}

impl Output for CheckedComplexDisk {
    fn ty() -> TypeExpr {
        <[f64; 3] as Output>::ty()
    }

    fn encode(self) -> Result<Value, OpError> {
        disk_coordinates(self.0).encode()
    }
}

impl CheckedComplexDisk {
    fn new(
        center: num_complex::Complex64,
        radius: f64,
    ) -> Result<Self, OpError> {
        crate::kernels::complex_ball::ComplexBall::new(center, radius)
            .map(Self)
            .ok_or_else(|| OpError::Domain("complex disk radius must be nonnegative".into()))
    }
}

binding!(
    complex_disk_new,
    "num_complex_disk_new",
    "numerical_special",
    "`num_complex_disk_new(center, radius)` constructs a disk `[real, imag, radius]` from a finite complex center and a finite nonnegative radius with the same physical dimension.",
    "crate::kernels::complex_ball::ComplexBall::new",
    DimensionPolicy::EqualPreserve,
    (center: num_complex::Complex64, radius: f64) -> Result<super::CheckedComplexDisk, OpError>,
    |center: num_complex::Complex64, radius: f64| super::CheckedComplexDisk::new(center, radius)
);

binding!(
    complex_disk_of,
    "num_complex_disk_of",
    "numerical_special",
    "`num_complex_disk_of(center)` returns the exact point disk of a finite complex center, preserving its dimension.",
    "crate::kernels::complex_ball::ComplexBall::of",
    DimensionPolicy::Preserve,
    (center: num_complex::Complex64) -> super::CheckedComplexDisk,
    |center: num_complex::Complex64| super::CheckedComplexDisk(crate::kernels::complex_ball::ComplexBall::of(center))
);

binding!(
    complex_disk_point,
    "num_complex_disk_point",
    "numerical_special",
    "`num_complex_disk_point(real, imag)` returns an exact point disk from finite real and imaginary coordinates with the same dimension.",
    "crate::kernels::complex_ball::ComplexBall::point",
    DimensionPolicy::EqualPreserve,
    (real: f64, imag: f64) -> super::CheckedComplexDisk,
    |real: f64, imag: f64| super::CheckedComplexDisk(crate::kernels::complex_ball::ComplexBall::point(real, imag))
);

binding!(
    complex_disk_zero,
    "num_complex_disk_zero",
    "numerical_special",
    "`num_complex_disk_zero()` returns a dimensionless point disk at 0.",
    "crate::kernels::complex_ball::ComplexBall::zero",
    DimensionPolicy::Dimensionless,
    () -> super::CheckedComplexDisk,
    || super::CheckedComplexDisk(crate::kernels::complex_ball::ComplexBall::zero())
);

binding!(
    complex_disk_one,
    "num_complex_disk_one",
    "numerical_special",
    "`num_complex_disk_one()` returns a dimensionless point disk at 1.",
    "crate::kernels::complex_ball::ComplexBall::one",
    DimensionPolicy::Dimensionless,
    () -> super::CheckedComplexDisk,
    || super::CheckedComplexDisk(crate::kernels::complex_ball::ComplexBall::one())
);

binding!(
    complex_disk_i,
    "num_complex_disk_i",
    "numerical_special",
    "`num_complex_disk_i()` returns a dimensionless point disk at the imaginary unit.",
    "crate::kernels::complex_ball::ComplexBall::i",
    DimensionPolicy::Dimensionless,
    () -> super::CheckedComplexDisk,
    || super::CheckedComplexDisk(crate::kernels::complex_ball::ComplexBall::i())
);

binding!(
    complex_disk_center,
    "num_complex_disk_center",
    "numerical_special",
    "`num_complex_disk_center(disk)` returns the disk center as a complex quantity.",
    "crate::kernels::complex_ball::ComplexBall::center",
    DimensionPolicy::Preserve,
    (disk: super::CheckedComplexDisk) -> num_complex::Complex64,
    |disk: super::CheckedComplexDisk| disk.0.center()
);

binding!(
    complex_disk_radius,
    "num_complex_disk_radius",
    "numerical_special",
    "`num_complex_disk_radius(disk)` returns the disk radius with the coordinates’ physical dimension.",
    "crate::kernels::complex_ball::ComplexBall::radius",
    DimensionPolicy::Preserve,
    (disk: super::CheckedComplexDisk) -> f64,
    |disk: super::CheckedComplexDisk| disk.0.radius()
);

binding!(
    complex_disk_contains,
    "num_complex_disk_contains",
    "numerical_special",
    "`num_complex_disk_contains(disk, point)` tests whether the point may lie in the closed disk. Uncertainty due to rounding at its boundary counts as possible membership.",
    "crate::kernels::complex_ball::ComplexBall::contains",
    DimensionPolicy::EqualInputs,
    (disk: super::CheckedComplexDisk, point: num_complex::Complex64) -> bool,
    |disk: super::CheckedComplexDisk, point: num_complex::Complex64| disk.0.contains(point)
);

binding!(
    complex_disk_contains_zero,
    "num_complex_disk_contains_zero",
    "numerical_special",
    "`num_complex_disk_contains_zero(disk)` tests whether zero may lie in the closed disk.",
    "crate::kernels::complex_ball::ComplexBall::contains_zero",
    DimensionPolicy::UnitlessResult,
    (disk: super::CheckedComplexDisk) -> bool,
    |disk: super::CheckedComplexDisk| disk.0.contains_zero()
);

binding!(
    complex_disk_inside_strictly,
    "num_complex_disk_inside_strictly",
    "numerical_special",
    "`num_complex_disk_inside_strictly(inner, outer)` certifies that the inner disk lies strictly in the interior of the outer disk.",
    "crate::kernels::complex_ball::ComplexBall::inside_strictly",
    DimensionPolicy::EqualInputs,
    (inner: super::CheckedComplexDisk, outer: super::CheckedComplexDisk) -> bool,
    |inner: super::CheckedComplexDisk, outer: super::CheckedComplexDisk| inner.0.inside_strictly(&outer.0)
);

binding!(
    complex_disk_overlaps,
    "num_complex_disk_overlaps",
    "numerical_special",
    "`num_complex_disk_overlaps(a, b)` tests whether two disks with the same coordinate dimension may intersect.",
    "crate::kernels::complex_ball::ComplexBall::overlaps",
    DimensionPolicy::EqualInputs,
    (a: super::CheckedComplexDisk, b: super::CheckedComplexDisk) -> bool,
    |a: super::CheckedComplexDisk, b: super::CheckedComplexDisk| a.0.overlaps(&b.0)
);

binding!(
    complex_disk_modulus_lower,
    "num_complex_disk_modulus_lower",
    "numerical_special",
    "`num_complex_disk_modulus_lower(disk)` returns an outward lower bound on the modulus of every point in the disk.",
    "crate::kernels::complex_ball::ComplexBall::modulus_lower",
    DimensionPolicy::Preserve,
    (disk: super::CheckedComplexDisk) -> f64,
    |disk: super::CheckedComplexDisk| disk.0.modulus_lower()
);

binding!(
    complex_disk_modulus_upper,
    "num_complex_disk_modulus_upper",
    "numerical_special",
    "`num_complex_disk_modulus_upper(disk)` returns an outward upper bound on the modulus of every point in the disk.",
    "crate::kernels::complex_ball::ComplexBall::modulus_upper",
    DimensionPolicy::Preserve,
    (disk: super::CheckedComplexDisk) -> f64,
    |disk: super::CheckedComplexDisk| disk.0.modulus_upper()
);

binding!(
    complex_disk_neg,
    "num_complex_disk_neg",
    "numerical_special",
    "`num_complex_disk_neg(disk)` negates a disk and preserves its coordinate dimension.",
    "crate::kernels::complex_ball::ComplexBall::neg_ball",
    DimensionPolicy::Preserve,
    (disk: super::CheckedComplexDisk) -> super::CheckedComplexDisk,
    |disk: super::CheckedComplexDisk| super::CheckedComplexDisk(disk.0.neg_ball())
);

binding!(
    complex_disk_conj,
    "num_complex_disk_conj",
    "numerical_special",
    "`num_complex_disk_conj(disk)` reflects the disk in the real axis and preserves its coordinate dimension.",
    "crate::kernels::complex_ball::ComplexBall::conj",
    DimensionPolicy::Preserve,
    (disk: super::CheckedComplexDisk) -> super::CheckedComplexDisk,
    |disk: super::CheckedComplexDisk| super::CheckedComplexDisk(disk.0.conj())
);

binding!(
    complex_disk_add,
    "num_complex_disk_add",
    "numerical_special",
    "`num_complex_disk_add(a, b)` encloses every sum of points in two disks of the same dimension.",
    "crate::kernels::complex_ball::ComplexBall::add_ball",
    DimensionPolicy::EqualPreserve,
    (a: super::CheckedComplexDisk, b: super::CheckedComplexDisk) -> super::CheckedComplexDisk,
    |a: super::CheckedComplexDisk, b: super::CheckedComplexDisk| super::CheckedComplexDisk(a.0.add_ball(&b.0))
);

binding!(
    complex_disk_sub,
    "num_complex_disk_sub",
    "numerical_special",
    "`num_complex_disk_sub(a, b)` encloses every difference of points in two disks of the same dimension.",
    "crate::kernels::complex_ball::ComplexBall::sub_ball",
    DimensionPolicy::EqualPreserve,
    (a: super::CheckedComplexDisk, b: super::CheckedComplexDisk) -> super::CheckedComplexDisk,
    |a: super::CheckedComplexDisk, b: super::CheckedComplexDisk| super::CheckedComplexDisk(a.0.sub_ball(&b.0))
);

binding!(
    complex_disk_mul,
    "num_complex_disk_mul",
    "numerical_special",
    "`num_complex_disk_mul(a, b)` encloses products of points in two disks. The coordinate dimensions multiply.",
    "crate::kernels::complex_ball::ComplexBall::mul_ball",
    DimensionPolicy::Product,
    (a: super::CheckedComplexDisk, b: super::CheckedComplexDisk) -> super::CheckedComplexDisk,
    |a: super::CheckedComplexDisk, b: super::CheckedComplexDisk| super::CheckedComplexDisk(a.0.mul_ball(&b.0))
);

binding!(
    complex_disk_square,
    "num_complex_disk_square",
    "numerical_special",
    "`num_complex_disk_square(disk)` encloses the squares of points in a disk. The result has squared coordinate dimension.",
    "crate::kernels::complex_ball::ComplexBall::sqr",
    DimensionPolicy::Square,
    (disk: super::CheckedComplexDisk) -> super::CheckedComplexDisk,
    |disk: super::CheckedComplexDisk| super::CheckedComplexDisk(disk.0.sqr())
);

binding!(
    complex_disk_scale,
    "num_complex_disk_scale",
    "numerical_special",
    "`num_complex_disk_scale(disk, scale)` multiplies the disk by a finite real scalar. The physical dimensions multiply.",
    "crate::kernels::complex_ball::ComplexBall::scale_real",
    DimensionPolicy::Product,
    (disk: super::CheckedComplexDisk, scale: f64) -> super::CheckedComplexDisk,
    |disk: super::CheckedComplexDisk, scale: f64| super::CheckedComplexDisk(disk.0.scale_real(scale))
);

binding!(
    complex_disk_recip,
    "num_complex_disk_recip",
    "numerical_special",
    "`num_complex_disk_recip(disk)` encloses reciprocals with inverse coordinate dimension. Returns None if the disk cannot be separated from zero.",
    "crate::kernels::complex_ball::ComplexBall::recip",
    DimensionPolicy::Inverse,
    (disk: super::CheckedComplexDisk) -> Option<super::CheckedComplexDisk>,
    |disk: super::CheckedComplexDisk| disk.0.recip().map(super::CheckedComplexDisk)
);

binding!(
    complex_disk_div,
    "num_complex_disk_div",
    "numerical_special",
    "`num_complex_disk_div(numerator, denominator)` encloses quotients with numerator over denominator dimension. Returns None if the denominator disk cannot be separated from zero.",
    "crate::kernels::complex_ball::ComplexBall::div_ball",
    DimensionPolicy::Quotient,
    (numerator: super::CheckedComplexDisk, denominator: super::CheckedComplexDisk) -> Option<super::CheckedComplexDisk>,
    |numerator: super::CheckedComplexDisk, denominator: super::CheckedComplexDisk| numerator.0.div_ball(&denominator.0).map(super::CheckedComplexDisk)
);

// First-order jets of functions of a dimensionless complex variable.
pub(in crate::rules) struct CheckedComplexDiskJet(
    pub(in crate::rules) crate::kernels::complex_ball::ComplexBallJet,
);

impl Arg for CheckedComplexDiskJet {
    fn ty() -> TypeExpr {
        <Self as Output>::ty()
    }

    fn decode(
        value: &Value,
        at: usize,
    ) -> Result<Self, OpError> {
        use super::native_object_api::field;
        let Value::Record(record) = value else {
            return Err(OpError::Type {
                expected: Box::new(<Self as Arg>::ty()),
                got: Box::new(value.type_expr()),
                arg: at,
            });
        };
        if record.fields.len() != 2 {
            return Err(OpError::Domain(
                "complex disk jet requires value and derivative".into(),
            ));
        }
        Ok(Self(crate::kernels::complex_ball::ComplexBallJet {
            val: CheckedComplexDisk::decode(field(value, "value")?, at)?.0,
            deriv: CheckedComplexDisk::decode(field(value, "derivative")?, at)?.0,
        }))
    }
}

impl Output for CheckedComplexDiskJet {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("value".into(), <CheckedComplexDisk as Output>::ty()),
            ("derivative".into(), <CheckedComplexDisk as Output>::ty()),
        ])
    }

    fn encode(self) -> Result<Value, OpError> {
        Ok(Value::Record(RecordValue {
            fields: vec![
                ("value".into(), CheckedComplexDisk(self.0.val).encode()?),
                (
                    "derivative".into(),
                    CheckedComplexDisk(self.0.deriv).encode()?,
                ),
            ],
        }))
    }
}

pub(in crate::rules) struct CheckedComplexDiskVector(
    pub(in crate::rules) Vec<crate::kernels::complex_ball::ComplexBall>,
);

impl Arg for CheckedComplexDiskVector {
    fn ty() -> TypeExpr {
        <Vec<[f64; 3]> as Arg>::ty()
    }

    fn decode(
        value: &Value,
        at: usize,
    ) -> Result<Self, OpError> {
        let disks = <Vec<[f64; 3]> as Arg>::decode(value, at)?;
        disks
            .into_iter()
            .map(decode_disk_cell)
            .collect::<Result<_, _>>()
            .map(Self)
    }
}

impl Output for CheckedComplexDiskVector {
    fn ty() -> TypeExpr {
        <Vec<[f64; 3]> as Output>::ty()
    }

    fn encode(self) -> Result<Value, OpError> {
        self.0
            .into_iter()
            .map(disk_coordinates)
            .collect::<Vec<_>>()
            .encode()
    }
}

pub(in crate::rules) struct CheckedComplexDiskMatrix(
    pub(in crate::rules) Vec<Vec<crate::kernels::complex_ball::ComplexBall>>,
);

impl Arg for CheckedComplexDiskMatrix {
    fn ty() -> TypeExpr {
        tensor_type(3)
    }

    fn decode(
        value: &Value,
        at: usize,
    ) -> Result<Self, OpError> {
        let (shape, scalars) = tensor(value, 3, at)?;
        let [rows, cols, 3] = shape.as_slice() else {
            return Err(OpError::Domain(
                "complex disk matrices need shape [rows, columns, 3]".into(),
            ));
        };
        if *rows == 0 || *cols == 0 {
            return Err(OpError::Domain(
                "complex disk matrices need nonempty row and column axes".into(),
            ));
        }
        let width = cols.checked_mul(3).ok_or(OpError::BudgetExceeded)?;
        let matrix = scalars
            .chunks_exact(width)
            .map(|row| {
                let (disks, _) = row.as_chunks::<3>();
                disks.iter().copied().map(decode_disk_cell).collect()
            })
            .collect::<Result<_, _>>()?;
        Ok(Self(matrix))
    }
}

impl Output for CheckedComplexDiskMatrix {
    fn ty() -> TypeExpr {
        tensor_type(3)
    }

    fn encode(self) -> Result<Value, OpError> {
        let (rows, cols) = self.shape();
        output_tensor(
            vec![rows, cols, 3],
            self.0
                .into_iter()
                .flatten()
                .flat_map(disk_coordinates)
                .collect(),
        )
    }
}

impl CheckedComplexDiskMatrix {
    fn shape(&self) -> (usize, usize) {
        (self.0.len(), self.0.first().map_or(0, Vec::len))
    }
}

fn decode_disk_cell(
    [real, imag, radius]: [f64; 3]
) -> Result<crate::kernels::complex_ball::ComplexBall, OpError> {
    Ok(CheckedComplexDisk::new(num_complex::Complex64::new(real, imag), radius)?.0)
}

const fn disk_coordinates(disk: crate::kernels::complex_ball::ComplexBall) -> [f64; 3] {
    [disk.mid.re, disk.mid.im, disk.rad]
}

fn check_disk_matrix_work(
    rows: usize,
    inner: usize,
    cols: usize,
) -> Result<(), OpError> {
    const WORK_LIMIT: usize = 1_048_576;
    // Three measured scalars can carry three additional unit witnesses.
    let output_values = rows
        .saturating_mul(cols)
        .saturating_mul(6)
        .saturating_add(1);
    if output_values > super::native_object_api::VALUE_LIMIT
        || rows.saturating_mul(inner).saturating_mul(cols) > WORK_LIMIT
    {
        return Err(OpError::Domain(
            "complex disk matrix exceeds the native value or work limit".into(),
        ));
    }
    crate::cancellation::checkpoint()?;
    Ok(())
}

fn disk_matrix_result<T>(value: Option<T>) -> Result<T, OpError> {
    crate::cancellation::checkpoint()?;
    value.ok_or_else(|| OpError::Domain("complex disk matrix dimensions do not agree".into()))
}

binding!(
    complex_disk_jet_variable,
    "num_complex_disk_jet_variable",
    "numerical_special",
    "`num_complex_disk_jet_variable(disk)` initializes a dimensionless jet with a unit derivative. The record stores value and derivative disks.",
    "crate::kernels::complex_ball::ComplexBallJet::variable",
    DimensionPolicy::Dimensionless,
    (disk: super::CheckedComplexDisk) -> super::CheckedComplexDiskJet,
    |disk: super::CheckedComplexDisk| super::CheckedComplexDiskJet(crate::kernels::complex_ball::ComplexBallJet::variable(disk.0))
);

binding!(
    complex_disk_jet_constant,
    "num_complex_disk_jet_constant",
    "numerical_special",
    "`num_complex_disk_jet_constant(disk)` initializes a dimensionless jet with a zero derivative. The record stores value and derivative disks.",
    "crate::kernels::complex_ball::ComplexBallJet::constant",
    DimensionPolicy::Dimensionless,
    (disk: super::CheckedComplexDisk) -> super::CheckedComplexDiskJet,
    |disk: super::CheckedComplexDisk| super::CheckedComplexDiskJet(crate::kernels::complex_ball::ComplexBallJet::constant(disk.0))
);

binding!(
    complex_disk_jet_add,
    "num_complex_disk_jet_add",
    "numerical_special",
    "`num_complex_disk_jet_add(a, b)` adds two dimensionless jets, enclosing both value and first derivative.",
    "crate::kernels::complex_ball::ComplexBallJet::add",
    DimensionPolicy::Dimensionless,
    (a: super::CheckedComplexDiskJet, b: super::CheckedComplexDiskJet) -> super::CheckedComplexDiskJet,
    |a: super::CheckedComplexDiskJet, b: super::CheckedComplexDiskJet| super::CheckedComplexDiskJet(a.0.add(&b.0))
);

binding!(
    complex_disk_jet_sub,
    "num_complex_disk_jet_sub",
    "numerical_special",
    "`num_complex_disk_jet_sub(a, b)` subtracts two dimensionless jets, enclosing both value and first derivative.",
    "crate::kernels::complex_ball::ComplexBallJet::sub",
    DimensionPolicy::Dimensionless,
    (a: super::CheckedComplexDiskJet, b: super::CheckedComplexDiskJet) -> super::CheckedComplexDiskJet,
    |a: super::CheckedComplexDiskJet, b: super::CheckedComplexDiskJet| super::CheckedComplexDiskJet(a.0.sub(&b.0))
);

binding!(
    complex_disk_jet_mul,
    "num_complex_disk_jet_mul",
    "numerical_special",
    "`num_complex_disk_jet_mul(a, b)` multiplies two dimensionless jets, enclosing both value and first derivative.",
    "crate::kernels::complex_ball::ComplexBallJet::mul",
    DimensionPolicy::Dimensionless,
    (a: super::CheckedComplexDiskJet, b: super::CheckedComplexDiskJet) -> super::CheckedComplexDiskJet,
    |a: super::CheckedComplexDiskJet, b: super::CheckedComplexDiskJet| super::CheckedComplexDiskJet(a.0.mul(&b.0))
);

binding!(
    complex_disk_jet_div,
    "num_complex_disk_jet_div",
    "numerical_special",
    "`num_complex_disk_jet_div(numerator, denominator)` divides dimensionless jets with the quotient rule. Returns None when a denominator cannot be separated from zero.",
    "crate::kernels::complex_ball::ComplexBallJet::div",
    DimensionPolicy::Dimensionless,
    (numerator: super::CheckedComplexDiskJet, denominator: super::CheckedComplexDiskJet) -> Option<super::CheckedComplexDiskJet>,
    |numerator: super::CheckedComplexDiskJet, denominator: super::CheckedComplexDiskJet| numerator.0.div(&denominator.0).map(super::CheckedComplexDiskJet)
);

binding!(
    complex_disk_jet_powi,
    "num_complex_disk_jet_powi",
    "numerical_special",
    "`num_complex_disk_jet_powi(value, exponent)` raises a dimensionless jet to a signed integer power and encloses its first derivative. The exponent and exponent minus one must fit i32.",
    "crate::kernels::complex_ball::ComplexBallJet::powi",
    DimensionPolicy::Dimensionless,
    (value: super::CheckedComplexDiskJet, exponent: i32) -> Result<Option<super::CheckedComplexDiskJet>, OpError>,
    |value: super::CheckedComplexDiskJet, exponent: i32| {
        if exponent == i32::MIN {
            return Err(OpError::Domain("jet exponent minus one must fit i32".into()));
        }
        Ok(value.0.powi(exponent).map(super::CheckedComplexDiskJet))
    }
);

binding!(
    complex_disk_matvec,
    "num_complex_disk_matvec",
    "numerical_special",
    "`num_complex_disk_matvec(matrix, vector)` encloses a nonempty disk matrix of shape `[rows, columns, 3]` times a disk vector of shape `[columns, 3]`. Physical dimensions multiply.",
    "crate::kernels::complex_ball::mat_vec_ball",
    DimensionPolicy::Product,
    (matrix: super::CheckedComplexDiskMatrix, vector: super::CheckedComplexDiskVector) -> Result<super::CheckedComplexDiskVector, OpError>,
    |matrix: super::CheckedComplexDiskMatrix, vector: super::CheckedComplexDiskVector| {
        let (rows, cols) = matrix.shape();
        super::check_disk_matrix_work(rows, cols, 1)?;
        super::disk_matrix_result(crate::kernels::complex_ball::mat_vec_ball(&matrix.0, &vector.0)).map(super::CheckedComplexDiskVector)
    }
);

binding!(
    complex_disk_matmul,
    "num_complex_disk_matmul",
    "numerical_special",
    "`num_complex_disk_matmul(a, b)` encloses a product of disk matrices. Each input has nonempty row and column axes, and every disk uses a trailing coordinate axis of length 3; physical dimensions multiply.",
    "crate::kernels::complex_ball::mat_mul_ball",
    DimensionPolicy::Product,
    (a: super::CheckedComplexDiskMatrix, b: super::CheckedComplexDiskMatrix) -> Result<super::CheckedComplexDiskMatrix, OpError>,
    |a: super::CheckedComplexDiskMatrix, b: super::CheckedComplexDiskMatrix| {
        let (rows, inner) = a.shape();
        let (_, cols) = b.shape();
        super::check_disk_matrix_work(rows, inner, cols)?;
        super::disk_matrix_result(crate::kernels::complex_ball::mat_mul_ball(&a.0, &b.0)).map(super::CheckedComplexDiskMatrix)
    }
);

binding!(
    complex_disk_solve,
    "num_complex_disk_solve",
    "numerical_special",
    "`num_complex_disk_solve(matrix, rhs)` encloses the solution to a square disk-valued system when residual certification succeeds. The solution dimension is right-hand side over matrix; failed certification gives None.",
    "crate::kernels::complex_ball::solve_linear_ball",
    DimensionPolicy::LeastSquares,
    (matrix: super::CheckedComplexDiskMatrix, rhs: super::CheckedComplexDiskVector) -> Result<Option<super::CheckedComplexDiskVector>, OpError>,
    |matrix: super::CheckedComplexDiskMatrix, rhs: super::CheckedComplexDiskVector| {
        let (rows, cols) = matrix.shape();
        if rows != cols || rows != rhs.0.len() {
            return Err(OpError::Domain("complex disk solve needs a square matrix and matching right-hand side".into()));
        }
        super::check_disk_matrix_work(rows, rows.saturating_mul(4), cols)?;
        let solution = crate::kernels::complex_ball::solve_linear_ball(&matrix.0, &rhs.0);
        crate::cancellation::checkpoint()?;
        Ok(solution.map(super::CheckedComplexDiskVector))
    }
);

binding!(
    complex_disk_jet_exp,
    "num_complex_disk_jet_exp",
    "numerical_special",
    "`num_complex_disk_jet_exp(value)` encloses the dimensionless exponential and its first derivative. A non-finite enclosure reports an error.",
    "crate::kernels::complex_ball::ComplexBallJet::exp",
    DimensionPolicy::Dimensionless,
    (value: super::CheckedComplexDiskJet) -> super::CheckedComplexDiskJet,
    |value: super::CheckedComplexDiskJet| super::CheckedComplexDiskJet(value.0.exp())
);

binding!(
    complex_disk_jet_sin,
    "num_complex_disk_jet_sin",
    "numerical_special",
    "`num_complex_disk_jet_sin(value)` encloses the dimensionless sine and its first derivative. A non-finite enclosure reports an error.",
    "crate::kernels::complex_ball::ComplexBallJet::sin",
    DimensionPolicy::Dimensionless,
    (value: super::CheckedComplexDiskJet) -> super::CheckedComplexDiskJet,
    |value: super::CheckedComplexDiskJet| super::CheckedComplexDiskJet(value.0.sin())
);

binding!(
    complex_disk_jet_cos,
    "num_complex_disk_jet_cos",
    "numerical_special",
    "`num_complex_disk_jet_cos(value)` encloses the dimensionless cosine and its first derivative. A non-finite enclosure reports an error.",
    "crate::kernels::complex_ball::ComplexBallJet::cos",
    DimensionPolicy::Dimensionless,
    (value: super::CheckedComplexDiskJet) -> super::CheckedComplexDiskJet,
    |value: super::CheckedComplexDiskJet| super::CheckedComplexDiskJet(value.0.cos())
);

binding!(
    complex_disk_jet_sinh,
    "num_complex_disk_jet_sinh",
    "numerical_special",
    "`num_complex_disk_jet_sinh(value)` encloses the dimensionless hyperbolic sine and its first derivative. A non-finite enclosure reports an error.",
    "crate::kernels::complex_ball::ComplexBallJet::sinh",
    DimensionPolicy::Dimensionless,
    (value: super::CheckedComplexDiskJet) -> super::CheckedComplexDiskJet,
    |value: super::CheckedComplexDiskJet| super::CheckedComplexDiskJet(value.0.sinh())
);

binding!(
    complex_disk_jet_cosh,
    "num_complex_disk_jet_cosh",
    "numerical_special",
    "`num_complex_disk_jet_cosh(value)` encloses the dimensionless hyperbolic cosine and its first derivative. A non-finite enclosure reports an error.",
    "crate::kernels::complex_ball::ComplexBallJet::cosh",
    DimensionPolicy::Dimensionless,
    (value: super::CheckedComplexDiskJet) -> super::CheckedComplexDiskJet,
    |value: super::CheckedComplexDiskJet| super::CheckedComplexDiskJet(value.0.cosh())
);

binding!(
    complex_disk_jet_ln,
    "num_complex_disk_jet_ln",
    "numerical_special",
    "`num_complex_disk_jet_ln(value)` encloses the dimensionless principal logarithm and its first derivative. Returns None when a finite analytic enclosure cannot be certified.",
    "crate::kernels::complex_ball::ComplexBallJet::ln",
    DimensionPolicy::Dimensionless,
    (value: super::CheckedComplexDiskJet) -> Option<super::CheckedComplexDiskJet>,
    |value: super::CheckedComplexDiskJet| value.0.ln().map(super::CheckedComplexDiskJet)
);

binding!(
    complex_disk_jet_sqrt,
    "num_complex_disk_jet_sqrt",
    "numerical_special",
    "`num_complex_disk_jet_sqrt(value)` encloses the dimensionless principal square root and its first derivative. Returns None when a finite analytic enclosure cannot be certified.",
    "crate::kernels::complex_ball::ComplexBallJet::sqrt",
    DimensionPolicy::Dimensionless,
    (value: super::CheckedComplexDiskJet) -> Option<super::CheckedComplexDiskJet>,
    |value: super::CheckedComplexDiskJet| value.0.sqrt().map(super::CheckedComplexDiskJet)
);

binding!(
    complex_disk_jet_tan,
    "num_complex_disk_jet_tan",
    "numerical_special",
    "`num_complex_disk_jet_tan(value)` encloses the dimensionless tangent and its first derivative. Returns None when a finite analytic enclosure cannot be certified.",
    "crate::kernels::complex_ball::ComplexBallJet::tan",
    DimensionPolicy::Dimensionless,
    (value: super::CheckedComplexDiskJet) -> Option<super::CheckedComplexDiskJet>,
    |value: super::CheckedComplexDiskJet| value.0.tan().map(super::CheckedComplexDiskJet)
);

binding!(
    complex_disk_jet_tanh,
    "num_complex_disk_jet_tanh",
    "numerical_special",
    "`num_complex_disk_jet_tanh(value)` encloses the dimensionless hyperbolic tangent and its first derivative. Returns None when a finite analytic enclosure cannot be certified.",
    "crate::kernels::complex_ball::ComplexBallJet::tanh",
    DimensionPolicy::Dimensionless,
    (value: super::CheckedComplexDiskJet) -> Option<super::CheckedComplexDiskJet>,
    |value: super::CheckedComplexDiskJet| value.0.tanh().map(super::CheckedComplexDiskJet)
);

binding!(
    complex_disk_jet_asin,
    "num_complex_disk_jet_asin",
    "numerical_special",
    "`num_complex_disk_jet_asin(value)` encloses the dimensionless principal inverse sine and its first derivative. Returns None when a finite analytic enclosure cannot be certified.",
    "crate::kernels::complex_ball::ComplexBallJet::asin",
    DimensionPolicy::Dimensionless,
    (value: super::CheckedComplexDiskJet) -> Option<super::CheckedComplexDiskJet>,
    |value: super::CheckedComplexDiskJet| value.0.asin().map(super::CheckedComplexDiskJet)
);

binding!(
    complex_disk_jet_acos,
    "num_complex_disk_jet_acos",
    "numerical_special",
    "`num_complex_disk_jet_acos(value)` encloses the dimensionless principal inverse cosine and its first derivative. Returns None when a finite analytic enclosure cannot be certified.",
    "crate::kernels::complex_ball::ComplexBallJet::acos",
    DimensionPolicy::Dimensionless,
    (value: super::CheckedComplexDiskJet) -> Option<super::CheckedComplexDiskJet>,
    |value: super::CheckedComplexDiskJet| value.0.acos().map(super::CheckedComplexDiskJet)
);

binding!(
    complex_disk_jet_atan,
    "num_complex_disk_jet_atan",
    "numerical_special",
    "`num_complex_disk_jet_atan(value)` encloses the dimensionless principal inverse tangent and its first derivative. Returns None when a finite analytic enclosure cannot be certified.",
    "crate::kernels::complex_ball::ComplexBallJet::atan",
    DimensionPolicy::Dimensionless,
    (value: super::CheckedComplexDiskJet) -> Option<super::CheckedComplexDiskJet>,
    |value: super::CheckedComplexDiskJet| value.0.atan().map(super::CheckedComplexDiskJet)
);

binding!(
    complex_disk_jet_asinh,
    "num_complex_disk_jet_asinh",
    "numerical_special",
    "`num_complex_disk_jet_asinh(value)` encloses the dimensionless principal inverse hyperbolic sine and its first derivative. Returns None when a finite analytic enclosure cannot be certified.",
    "crate::kernels::complex_ball::ComplexBallJet::asinh",
    DimensionPolicy::Dimensionless,
    (value: super::CheckedComplexDiskJet) -> Option<super::CheckedComplexDiskJet>,
    |value: super::CheckedComplexDiskJet| value.0.asinh().map(super::CheckedComplexDiskJet)
);

binding!(
    complex_disk_jet_acosh,
    "num_complex_disk_jet_acosh",
    "numerical_special",
    "`num_complex_disk_jet_acosh(value)` encloses the dimensionless principal inverse hyperbolic cosine and its first derivative. Returns None when a finite analytic enclosure cannot be certified.",
    "crate::kernels::complex_ball::ComplexBallJet::acosh",
    DimensionPolicy::Dimensionless,
    (value: super::CheckedComplexDiskJet) -> Option<super::CheckedComplexDiskJet>,
    |value: super::CheckedComplexDiskJet| value.0.acosh().map(super::CheckedComplexDiskJet)
);

binding!(
    complex_disk_jet_atanh,
    "num_complex_disk_jet_atanh",
    "numerical_special",
    "`num_complex_disk_jet_atanh(value)` encloses the dimensionless principal inverse hyperbolic tangent and its first derivative. Returns None when a finite analytic enclosure cannot be certified.",
    "crate::kernels::complex_ball::ComplexBallJet::atanh",
    DimensionPolicy::Dimensionless,
    (value: super::CheckedComplexDiskJet) -> Option<super::CheckedComplexDiskJet>,
    |value: super::CheckedComplexDiskJet| value.0.atanh().map(super::CheckedComplexDiskJet)
);

binding!(
    complex_disk_jet_pow,
    "num_complex_disk_jet_pow",
    "numerical_special",
    "`num_complex_disk_jet_pow(base, exponent)` encloses the principal complex power of two dimensionless jets and its first derivative. Returns None when a principal logarithm cannot be certified.",
    "crate::kernels::complex_ball::ComplexBallJet::powc",
    DimensionPolicy::Dimensionless,
    (base: super::CheckedComplexDiskJet, exponent: super::CheckedComplexDiskJet) -> Option<super::CheckedComplexDiskJet>,
    |base: super::CheckedComplexDiskJet, exponent: super::CheckedComplexDiskJet| base.0.powc(&exponent.0).map(super::CheckedComplexDiskJet)
);

binding!(
    complex_disk_sec,
    "num_complex_disk_sec",
    "numerical_special",
    "`num_complex_disk_sec(disk)` encloses the dimensionless secant. Returns None if the denominator cannot be certified nonzero.",
    "crate::kernels::complex_ball::ComplexBall::sec",
    DimensionPolicy::Dimensionless,
    (disk: super::CheckedComplexDisk) -> Option<super::CheckedComplexDisk>,
    |disk: super::CheckedComplexDisk| disk.0.sec().map(super::CheckedComplexDisk)
);

binding!(
    complex_disk_csc,
    "num_complex_disk_csc",
    "numerical_special",
    "`num_complex_disk_csc(disk)` encloses the dimensionless cosecant. Returns None if the denominator cannot be certified nonzero.",
    "crate::kernels::complex_ball::ComplexBall::csc",
    DimensionPolicy::Dimensionless,
    (disk: super::CheckedComplexDisk) -> Option<super::CheckedComplexDisk>,
    |disk: super::CheckedComplexDisk| disk.0.csc().map(super::CheckedComplexDisk)
);

binding!(
    complex_disk_cot,
    "num_complex_disk_cot",
    "numerical_special",
    "`num_complex_disk_cot(disk)` encloses the dimensionless cotangent. Returns None if the denominator cannot be certified nonzero.",
    "crate::kernels::complex_ball::ComplexBall::cot",
    DimensionPolicy::Dimensionless,
    (disk: super::CheckedComplexDisk) -> Option<super::CheckedComplexDisk>,
    |disk: super::CheckedComplexDisk| disk.0.cot().map(super::CheckedComplexDisk)
);

binding!(
    complex_disk_sech,
    "num_complex_disk_sech",
    "numerical_special",
    "`num_complex_disk_sech(disk)` encloses the dimensionless hyperbolic secant. Returns None if the denominator cannot be certified nonzero.",
    "crate::kernels::complex_ball::ComplexBall::sech",
    DimensionPolicy::Dimensionless,
    (disk: super::CheckedComplexDisk) -> Option<super::CheckedComplexDisk>,
    |disk: super::CheckedComplexDisk| disk.0.sech().map(super::CheckedComplexDisk)
);

binding!(
    complex_disk_csch,
    "num_complex_disk_csch",
    "numerical_special",
    "`num_complex_disk_csch(disk)` encloses the dimensionless hyperbolic cosecant. Returns None if the denominator cannot be certified nonzero.",
    "crate::kernels::complex_ball::ComplexBall::csch",
    DimensionPolicy::Dimensionless,
    (disk: super::CheckedComplexDisk) -> Option<super::CheckedComplexDisk>,
    |disk: super::CheckedComplexDisk| disk.0.csch().map(super::CheckedComplexDisk)
);

binding!(
    complex_disk_coth,
    "num_complex_disk_coth",
    "numerical_special",
    "`num_complex_disk_coth(disk)` encloses the dimensionless hyperbolic cotangent. Returns None if the denominator cannot be certified nonzero.",
    "crate::kernels::complex_ball::ComplexBall::coth",
    DimensionPolicy::Dimensionless,
    (disk: super::CheckedComplexDisk) -> Option<super::CheckedComplexDisk>,
    |disk: super::CheckedComplexDisk| disk.0.coth().map(super::CheckedComplexDisk)
);

struct DiskFiniteness(crate::kernels::complex_ball::ComplexBall);

impl Arg for DiskFiniteness {
    fn ty() -> TypeExpr {
        <CheckedComplexDisk as Arg>::ty()
    }

    fn decode(
        value: &Value,
        at: usize,
    ) -> Result<Self, OpError> {
        let Value::Tensor(tensor) = value else {
            return Err(OpError::Domain("expected complex disk coordinates".into()));
        };
        let [real, imag, radius] = tensor.data.as_slice() else {
            return Err(OpError::Domain(
                "expected three complex disk coordinates".into(),
            ));
        };
        if tensor.shape != [3] {
            return Err(OpError::Domain("expected a rank-one complex disk".into()));
        }
        let scalar = |v: &Value| {
            match v {
                | Value::Float(x) => Ok(*x),
                | v => <f64 as Arg>::decode(v, at),
            }
        };
        Ok(Self(decode_disk_cell([
            scalar(real)?,
            scalar(imag)?,
            scalar(radius)?,
        ])?))
    }
}

binding!(
    complex_disk_is_finite,
    "num_complex_disk_is_finite",
    "numerical_special",
    "Tests whether well-formed disk coordinates and radius are all finite. An infinite radius gives false; malformed disks give an error.",
    "crate::kernels::complex_ball::ComplexBall::is_finite",
    DimensionPolicy::UnitlessResult,
    (disk: super::DiskFiniteness) -> bool,
    |disk: super::DiskFiniteness| disk.0.is_finite()
);

binding!(
    complex_disk_erf,
    "num_complex_disk_erf",
    "numerical_special",
    "Encloses the complex error function for a dimensionless disk, using an explicit series remainder or a whole-disk integral bound.",
    "crate::kernels::complex_ball::ComplexBall::erf",
    DimensionPolicy::Dimensionless,
    (disk: super::CheckedComplexDisk) -> Option<super::CheckedComplexDisk>,
    |disk: super::CheckedComplexDisk| disk.0.erf().map(super::CheckedComplexDisk)
);

binding!(
    complex_disk_log_gamma,
    "num_complex_disk_log_gamma",
    "numerical_special",
    "Encloses principal log Gamma on a dimensionless disk off the nonpositive real cut. Uses outward recurrence and a bounded Stirling remainder.",
    "crate::kernels::complex_ball::ComplexBall::ln_gamma",
    DimensionPolicy::Dimensionless,
    (disk: super::CheckedComplexDisk) -> Option<super::CheckedComplexDisk>,
    |disk: super::CheckedComplexDisk| disk.0.ln_gamma().map(super::CheckedComplexDisk)
);

binding!(
    complex_disk_gamma,
    "num_complex_disk_gamma",
    "numerical_special",
    "Encloses complex Gamma for a dimensionless disk on the supported log-Gamma domain. None indicates an unavailable finite enclosure.",
    "crate::kernels::complex_ball::ComplexBall::gamma",
    DimensionPolicy::Dimensionless,
    (disk: super::CheckedComplexDisk) -> Option<super::CheckedComplexDisk>,
    |disk: super::CheckedComplexDisk| disk.0.gamma().map(super::CheckedComplexDisk)
);

binding!(
    complex_disk_bessel_j,
    "num_complex_disk_bessel_j",
    "numerical_special",
    "Encloses the principal Bessel function J of nonnegative real or negative integral order. Fractional orders use the principal branch on dimensionless disks.",
    "crate::kernels::complex_ball::ComplexBall::bessel_j",
    DimensionPolicy::Dimensionless,
    (disk: super::CheckedComplexDisk, order: f64) -> Option<super::CheckedComplexDisk>,
    |disk: super::CheckedComplexDisk, order: f64| disk.0.bessel_j(order).map(super::CheckedComplexDisk)
);
