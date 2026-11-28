// Portable contour geometry. The zero-based arc indices refer to segment positions.
struct ContourPath {
    starts: Vec<num_complex::Complex64>,
    ends: Vec<num_complex::Complex64>,
    arc_segments: Vec<usize>,
    arc_centers: Vec<num_complex::Complex64>,
    arc_radii: Vec<f64>,
    theta_start: Vec<f64>,
    theta_end: Vec<f64>,
    closed: bool,
}
record!(ContourPath {
    starts: Vec<num_complex::Complex64>,
    ends: Vec<num_complex::Complex64>,
    arc_segments: Vec<usize>,
    arc_centers: Vec<num_complex::Complex64>,
    arc_radii: Vec<f64>,
    theta_start: Vec<f64>,
    theta_end: Vec<f64>,
    closed: bool
});

fn contour_points_join(
    a: num_complex::Complex64,
    b: num_complex::Complex64,
    extent: f64,
) -> bool {
    let scale = extent
        .max(a.re.abs())
        .max(a.im.abs())
        .max(b.re.abs())
        .max(b.im.abs());
    scale <= 0.0 || (a / scale - b / scale).norm() <= 64.0 * f64::EPSILON
}

fn portable_contour(
    segments: &[crate::kernels::contour_deform::ContourSegment],
    closed: bool,
    scale: f64,
) -> Result<ContourPath, OpError> {
    use crate::kernels::contour_deform::ContourSegment;
    use num_complex::Complex64;
    let arc_count = segments
        .iter()
        .filter(|segment| matches!(segment, ContourSegment::Arc { .. }))
        .count();
    let mut path = ContourPath {
        starts: Vec::with_capacity(segments.len()),
        ends: Vec::with_capacity(segments.len()),
        arc_segments: Vec::with_capacity(arc_count),
        arc_centers: Vec::with_capacity(arc_count),
        arc_radii: Vec::with_capacity(arc_count),
        theta_start: Vec::with_capacity(arc_count),
        theta_end: Vec::with_capacity(arc_count),
        closed,
    };
    let mut first_extent = 0.0;
    let mut previous_extent = 0.0_f64;
    for (index, segment) in segments.iter().copied().enumerate() {
        let (start, end, extent) = match segment {
            | ContourSegment::Line { start, end } => {
                let (start, end) = (start * scale, end * scale);
                let delta = end - start;
                (start, end, delta.re.abs().max(delta.im.abs()))
            },
            | ContourSegment::Arc {
                center,
                radius,
                theta0,
                theta1,
            } => {
                let center = center * scale;
                let radius = radius * scale;
                if radius <= 0.0
                    || !radius.is_finite()
                    || (theta1 - theta0).abs() <= 0.0
                    || !(theta1 - theta0).is_finite()
                    || ![center.re, center.im].into_iter().all(|axis| {
                        let lower = axis - radius;
                        let upper = axis + radius;
                        lower.is_finite() && upper.is_finite() && lower < axis && upper > axis
                    })
                {
                    return Err(OpError::Domain(
                        "contour arc is not representable in Float coordinates".into(),
                    ));
                }
                path.arc_segments.push(index);
                path.arc_centers.push(center);
                path.arc_radii.push(radius);
                path.theta_start.push(theta0);
                path.theta_end.push(theta1);
                (
                    center + Complex64::from_polar(radius, theta0),
                    center + Complex64::from_polar(radius, theta1),
                    radius,
                )
            },
        };
        if !(start.is_finite() && end.is_finite() && (end - start).is_finite()) {
            return Err(OpError::Domain(
                "contour segment is not representable in Float coordinates".into(),
            ));
        }
        if path.ends.last().is_some_and(|&previous| {
            !contour_points_join(previous, start, extent.max(previous_extent))
        }) {
            return Err(OpError::Domain(
                "contour segment endpoints do not join in Float coordinates".into(),
            ));
        }
        if index == 0 {
            first_extent = extent;
        }
        previous_extent = extent;
        path.starts.push(start);
        path.ends.push(end);
    }
    if closed
        && !path
            .starts
            .first()
            .zip(path.ends.last())
            .is_some_and(|(&a, &b)| contour_points_join(a, b, first_extent.max(previous_extent)))
    {
        return Err(OpError::Domain(
            "closed contour endpoints do not coincide in Float coordinates".into(),
        ));
    }
    crate::cancellation::checkpoint()?;
    Ok(path)
}

fn contour_cut_angle(
    inner: f64,
    outer: f64,
    angle: f64,
    gap: f64,
) -> Result<f64, OpError> {
    use std::f64::consts::PI;
    use std::f64::consts::TAU;
    if inner <= 0.0 || outer <= inner || gap <= 0.0 || gap >= PI || angle.abs() > f64::from(1 << 20)
    {
        return Err(OpError::Domain("contour needs 0 < inner radius < outer radius, a gap in (0, pi), and a cut angle in +/-2^20 radians".into()));
    }
    let normalized = angle.rem_euclid(TAU);
    let turn = normalized + TAU;
    let above = normalized + gap;
    let below = turn - gap;
    if !(above > normalized && below < turn && above < below) {
        return Err(OpError::Domain(
            "contour gap is indistinguishable in Float angles".into(),
        ));
    }
    Ok(normalized)
}

fn checked_indented_contour(
    radius: f64,
    poles: Vec<f64>,
    indent: f64,
    upper: bool,
) -> Result<ContourPath, OpError> {
    use crate::kernels::contour_deform::indented_semicircle_contour;
    if radius <= 0.0 || indent <= 0.0 || indent >= radius || poles.len() > 4096 {
        return Err(OpError::Domain(
            "indented contour needs 0 < indent < radius and at most 4096 poles".into(),
        ));
    }
    let indent = indent / radius;
    let poles: Vec<f64> = poles.into_iter().map(|pole| pole / radius).collect();
    if poles
        .iter()
        .any(|&pole| pole <= -1.0 + indent || pole >= 1.0 - indent)
    {
        return Err(OpError::Domain(
            "indented contour poles need a full indent of clearance from the diameter endpoints"
                .into(),
        ));
    }
    // The kernel's absolute geometric tolerance now works in a unit-radius
    // frame. Require a returned arc for every pole rather than losing an indent.
    let segments = indented_semicircle_contour(1.0, &poles, indent, upper);
    let path = portable_contour(&segments, true, radius)?;
    if path.arc_segments.len() != poles.len() + 1 {
        return Err(OpError::Domain(
            "contour poles or indents cannot be separated in Float coordinates".into(),
        ));
    }
    Ok(path)
}

/// Sizes all obstacle tensors before conversion and allocation.
pub(in crate::rules) fn contour_deformation_preflight(args: &[Value]) -> Result<(), OpError> {
    use crate::kernels::contour_deform::MAX_OBSTACLES;
    let count = |index| -> Result<usize, OpError> {
        match args.get(index) {
            | Some(Value::Tensor(tensor)) if tensor.shape == [tensor.data.len()] => {
                Ok(tensor.data.len())
            },
            | _ => {
                Err(OpError::Domain(
                    "contour obstacles need rank-one tensors".into(),
                ))
            },
        }
    };
    let (poles, origins, angles) = (count(2)?, count(3)?, count(5)?);
    if origins != angles {
        return Err(OpError::Domain(
            "branch-ray origins and angles have different lengths".into(),
        ));
    }
    if poles.saturating_add(origins) > MAX_OBSTACLES {
        return Err(OpError::BudgetExceeded);
    }
    crate::cancellation::checkpoint()?;
    Ok(())
}

pub(in crate::rules) fn contour_complex_positions(
    value: &Value,
    dim: &crate::kernel::signature::DimExpr,
    index: usize,
    scale: f64,
) -> Result<Vec<Complex64>, OpError> {
    use super::runtime_dimensions::canonical;
    use super::runtime_dimensions::dimension;
    if dimension(value)? != canonical(dim)? {
        return Err(OpError::Domain(
            "contour obstacle coordinates have the wrong physical dimension".into(),
        ));
    }
    Ok(Vec::<Complex64>::decode(value, index)?
        .into_iter()
        .map(|z| z / scale)
        .collect())
}

pub(in crate::rules) fn checked_deformation(
    start: Complex64,
    end: Complex64,
    poles: &[Complex64],
    origins: &[Complex64],
    clearance: f64,
    angles: &[f64],
) -> Result<Vec<crate::kernels::contour_deform::ContourSegment>, OpError> {
    use crate::kernels::contour_deform::BranchRay;
    use crate::kernels::contour_deform::DeformationError;
    use crate::kernels::contour_deform::deformation_checked;
    let rays: Vec<_> = origins
        .iter()
        .zip(angles)
        .map(|(&origin, &angle)| BranchRay { origin, angle })
        .collect();
    deformation_checked(start, end, poles, &rays, clearance).map_err(|error| {
        match error {
            | DeformationError::WorkLimit => OpError::BudgetExceeded,
            | DeformationError::Cancelled => crate::cancellation::Cancelled.into(),
            | error => OpError::Domain(error.to_string()),
        }
    })
}

binding!(
    contour_deform,
    "num_contour_deform_checked",
    "numerical_special",
    "Constructs an open contour from start to end with positive pole clearance and no branch-ray crossing. Endpoint, pole, ray origin, and clearance dimensions agree. Angles are bare radians. Up to 64 poles and rays are allowed in total. Returns an error when the bounded construction cannot verify a path.",
    "crate::kernels::contour_deform::deform_path_around_obstacles",
    DimensionPolicy::ContourFiveCoordinates,
    (start: num_complex::Complex64, end: num_complex::Complex64, poles: Vec<num_complex::Complex64>, origins: Vec<num_complex::Complex64>, clearance: f64, angles: Vec<f64>) -> Result<super::ContourPath, OpError>,
    |start: num_complex::Complex64, end: num_complex::Complex64, poles: Vec<num_complex::Complex64>, origins: Vec<num_complex::Complex64>, clearance: f64, angles: Vec<f64>| {
        let segments = super::checked_deformation(start, end, &poles, &origins, clearance, &angles)?;
        super::portable_contour(&segments, false, 1.0)
    }
);

binding!(
    contour_keyhole,
    "num_contour_keyhole",
    "numerical_special",
    "`num_contour_keyhole(branch_point, inner, outer, cut_angle, gap)` returns a closed, positively oriented keyhole path around a cut. The branch point and both radii share units; angles are bare radians.",
    "crate::kernels::contour_deform::keyhole_contour",
    DimensionPolicy::ContourThreeCoordinates,
    (branch_point: num_complex::Complex64, inner: f64, outer: f64, cut_angle: f64, gap: f64) -> Result<super::ContourPath, OpError>,
    |branch_point: num_complex::Complex64, inner: f64, outer: f64, cut_angle: f64, gap: f64| {
        let angle = super::contour_cut_angle(inner, outer, cut_angle, gap)?;
        super::portable_contour(&crate::kernels::contour_deform::keyhole_contour(branch_point, inner, outer, angle, gap), true, 1.0)
    }
);

binding!(
    contour_hankel,
    "num_contour_hankel",
    "numerical_special",
    "`num_contour_hankel(inner, outer, cut_angle, gap)` returns an open Hankel path around a cut through zero. The radii share units; angles are bare radians.",
    "crate::kernels::contour_deform::hankel_contour",
    DimensionPolicy::ContourTwoCoordinates,
    (inner: f64, outer: f64, cut_angle: f64, gap: f64) -> Result<super::ContourPath, OpError>,
    |inner: f64, outer: f64, cut_angle: f64, gap: f64| {
        let angle = super::contour_cut_angle(inner, outer, cut_angle, gap)?;
        super::portable_contour(&crate::kernels::contour_deform::hankel_contour(inner, outer, angle, gap), false, 1.0)
    }
);

binding!(
    contour_dogbone,
    "num_contour_dogbone",
    "numerical_special",
    "`num_contour_dogbone(a, b, radius)` returns the closed contour enclosing the finite cut from a to b at a positive offset. Coincident endpoints give a circle.",
    "crate::kernels::contour_deform::dogbone_contour",
    DimensionPolicy::ContourThreeCoordinates,
    (a: num_complex::Complex64, b: num_complex::Complex64, radius: f64) -> Result<super::ContourPath, OpError>,
    |a: num_complex::Complex64, b: num_complex::Complex64, radius: f64| {
        if radius <= 0.0 || !(b - a).is_finite() {
            return Err(OpError::Domain("dogbone contour needs a positive radius and representable endpoint difference".into()));
        }
        super::portable_contour(&crate::kernels::contour_deform::dogbone_contour(a, b, radius), true, 1.0)
    }
);

binding!(
    contour_indented_semicircle,
    "num_contour_indented_semicircle",
    "numerical_special",
    "`num_contour_indented_semicircle(radius, poles, indent, upper)` returns a closed semicircle with a bypass for each real pole. Internal geometry uses a unit-radius frame and returns the original coordinate units.",
    "crate::kernels::contour_deform::indented_semicircle_contour",
    DimensionPolicy::ContourThreeCoordinates,
    (radius: f64, poles: Vec<f64>, indent: f64, upper: bool) -> Result<super::ContourPath, OpError>,
    |radius: f64, poles: Vec<f64>, indent: f64, upper: bool| super::checked_indented_contour(radius, poles, indent, upper)
);

binding!(
    contour_bromwich,
    "num_contour_bromwich",
    "numerical_special",
    "`num_contour_bromwich(c, y_max)` returns the open upward interval from c-i*y_max to c+i*y_max, with y_max positive.",
    "crate::kernels::contour_deform::bromwich_contour",
    DimensionPolicy::ContourTwoCoordinates,
    (c: f64, y_max: f64) -> Result<super::ContourPath, OpError>,
    |c: f64, y_max: f64| {
        if y_max <= 0.0 {
            return Err(OpError::Domain("Bromwich contour needs a positive imaginary extent".into()));
        }
        super::portable_contour(&crate::kernels::contour_deform::bromwich_contour(c, y_max), false, 1.0)
    }
);

binding!(
    contour_rectangle,
    "num_contour_rectangle",
    "numerical_special",
    "`num_contour_rectangle(x0, x1, y0, y1)` returns a closed counterclockwise rectangle. Both real bounds and both imaginary bounds must strictly increase and share units.",
    "crate::kernels::contour_deform::rectangular_contour",
    DimensionPolicy::ContourFourCoordinates,
    (x0: f64, x1: f64, y0: f64, y1: f64) -> Result<super::ContourPath, OpError>,
    |x0: f64, x1: f64, y0: f64, y1: f64| {
        if x0 >= x1 || y0 >= y1 {
            return Err(OpError::Domain("rectangular contour needs strictly increasing bounds on both axes".into()));
        }
        super::portable_contour(&crate::kernels::contour_deform::rectangular_contour(x0, x1, y0, y1), true, 1.0)
    }
);

binding!(
    contour_sector,
    "num_contour_sector",
    "numerical_special",
    "`num_contour_sector(radius, angle)` returns a closed sector with positive radius and a nonzero signed sweep of at most one turn in radians.",
    "crate::kernels::contour_deform::sector_contour",
    DimensionPolicy::ContourOneCoordinate,
    (radius: f64, angle: f64) -> Result<super::ContourPath, OpError>,
    |radius: f64, angle: f64| {
        if radius <= 0.0 || angle.abs() <= 0.0 || angle.abs() > std::f64::consts::TAU {
            return Err(OpError::Domain("sector contour needs a positive radius and a nonzero sweep of at most one turn".into()));
        }
        super::portable_contour(&crate::kernels::contour_deform::sector_contour(radius, angle), true, 1.0)
    }
);

/// The measured contour record used by the checked quadrature callbacks.
pub(in crate::rules) fn contour_path_type() -> TypeExpr {
    super::dimensions::apply(
        Signature::new(vec![TypeExpr::Float], ContourPath::ty()),
        DimensionPolicy::ContourOneCoordinate,
    )
    .ret
}

fn checked_contour_field<'a>(
    path: &'a Value,
    name: &str,
    expected: &crate::kernel::signature::DimExpr,
) -> Result<&'a Value, OpError> {
    use super::native_object_api::field;
    use super::runtime_dimensions::canonical;
    use super::runtime_dimensions::dimension;
    let field = field(path, name)?;
    let Value::Tensor(tensor) = field else {
        return Err(OpError::Domain(format!(
            "contour field {name} must be a tensor"
        )));
    };
    if tensor.data.len() > 4096 {
        return Err(OpError::BudgetExceeded);
    }
    if dimension(field)? != canonical(expected)? {
        return Err(OpError::Domain(format!(
            "contour field {name} has the wrong physical dimension"
        )));
    }
    Ok(field)
}

/// Reconstructs the piecewise path, checking measured coordinates and arc metadata.
pub(in crate::rules) fn decode_contour_path(
    value: &Value,
    dim: &crate::kernel::signature::DimExpr,
    scale: f64,
) -> Result<Vec<crate::kernels::contour_deform::ContourSegment>, OpError> {
    use super::native_object_api::field;
    use crate::kernel::signature::DimExpr;
    use crate::kernels::contour_deform::ContourSegment;
    use num_complex::Complex64;

    let Value::Record(record) = value else {
        return Err(OpError::Domain(
            "contour must be a portable path record".into(),
        ));
    };
    if record.fields.len() != 8 {
        return Err(OpError::Domain("contour path needs eight fields".into()));
    }
    let starts = Vec::<Complex64>::decode(checked_contour_field(value, "starts", dim)?, 0)?;
    let ends = Vec::<Complex64>::decode(checked_contour_field(value, "ends", dim)?, 0)?;
    let arc_segments = Vec::<usize>::decode(
        checked_contour_field(value, "arc_segments", &DimExpr::Dimless)?,
        0,
    )?;
    let arc_centers =
        Vec::<Complex64>::decode(checked_contour_field(value, "arc_centers", dim)?, 0)?;
    let arc_radii = Vec::<f64>::decode(checked_contour_field(value, "arc_radii", dim)?, 0)?;
    let theta_start = Vec::<f64>::decode(
        checked_contour_field(value, "theta_start", &DimExpr::Dimless)?,
        0,
    )?;
    let theta_end = Vec::<f64>::decode(
        checked_contour_field(value, "theta_end", &DimExpr::Dimless)?,
        0,
    )?;
    let closed = bool::decode(field(value, "closed")?, 0)?;
    let arcs = arc_segments.len();
    if starts.len() != ends.len()
        || starts.len() > 4096
        || [
            arc_centers.len(),
            arc_radii.len(),
            theta_start.len(),
            theta_end.len(),
        ]
        .into_iter()
        .any(|size| size != arcs)
        || arc_segments
            .last()
            .is_some_and(|&index| index >= starts.len())
        || arc_segments.windows(2).any(|pair| pair[0] >= pair[1])
    {
        return Err(OpError::Domain(
            "contour endpoint and arc metadata lengths or positions disagree".into(),
        ));
    }
    let mut arc_info = arc_segments
        .into_iter()
        .zip(arc_centers)
        .zip(arc_radii)
        .zip(theta_start)
        .zip(theta_end)
        .peekable();
    let mut segments = Vec::with_capacity(starts.len());
    for (index, (&start, &end)) in starts.iter().zip(&ends).enumerate() {
        let start = start / scale;
        let end = end / scale;
        let segment = if let Some(((((arc_index, center), radius), theta0), theta1)) =
            arc_info.peek().copied()
            && arc_index == index
        {
            arc_info.next();
            let center = center / scale;
            let radius = radius / scale;
            if !contour_points_join(
                start,
                center + Complex64::from_polar(radius, theta0),
                radius,
            ) || !contour_points_join(
                end,
                center + Complex64::from_polar(radius, theta1),
                radius,
            ) {
                return Err(OpError::Domain(
                    "contour arc metadata does not meet its declared endpoints".into(),
                ));
            }
            ContourSegment::Arc {
                center,
                radius,
                theta0,
                theta1,
            }
        } else {
            ContourSegment::Line { start, end }
        };
        segments.push(segment);
    }
    // Reuse the constructors' continuity checks on the canonical coordinates.
    portable_contour(&segments, closed, 1.0)?;
    Ok(segments)
}
