// Checked charts and geodesic distances over the existing complex-coordinate kernels.
use crate::kernels::complex_coords as hyperbolic;

binding!(
    uhp_to_disk,
    "num_uhp_to_disk",
    "numerical_special",
    "Maps a dimensionless complex point in the open upper half-plane to the open Poincare unit disk. A nonpositive input imaginary part or a result that rounds outside the open disk produces a domain error.",
    "crate::kernels::complex_coords::uhp_to_disk",
    DimensionPolicy::Dimensionless,
    (z: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |z: num_complex::Complex64| {
        if z.im <= 0.0 {
            return Err(OpError::Domain("the upper half-plane requires a positive imaginary part".into()));
        }
        super::hyperbolic::uhp_to_disk(z)
            .filter(|w| super::hyperbolic::poincare_disk_gap(*w) > 0.0)
            .ok_or_else(|| OpError::Domain("the Cayley image is outside the representable open unit disk".into()))
    }
);
binding!(
    disk_to_uhp,
    "num_disk_to_uhp",
    "numerical_special",
    "Maps a dimensionless complex point in the open Poincare unit disk to the open upper half-plane. An input on or outside the unit circle or a nonfinite result produces a domain error.",
    "crate::kernels::complex_coords::disk_to_uhp",
    DimensionPolicy::Dimensionless,
    (w: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |w: num_complex::Complex64| {
        super::hyperbolic::disk_to_uhp(w)
            .filter(|z| z.im > 0.0)
            .ok_or_else(|| OpError::Domain("the inverse Cayley transform requires an interior disk point and a representable upper-half-plane image".into()))
    }
);
binding!(
    poincare_disk_distance,
    "num_poincare_disk_distance",
    "numerical_special",
    "Returns the dimensionless hyperbolic geodesic distance between two complex points inside the open Poincare unit disk.",
    "crate::kernels::complex_coords::poincare_disk_distance",
    DimensionPolicy::Dimensionless,
    (z1: num_complex::Complex64, z2: num_complex::Complex64) -> Result<f64, OpError>,
    |z1: num_complex::Complex64, z2: num_complex::Complex64| {
        super::hyperbolic::poincare_disk_distance(z1, z2)
            .ok_or_else(|| OpError::Domain("both points must lie inside the open Poincare disk".into()))
    }
);
binding!(
    uhp_distance,
    "num_uhp_distance",
    "numerical_special",
    "Returns the dimensionless hyperbolic geodesic distance between two complex points with positive imaginary parts.",
    "crate::kernels::complex_coords::uhp_distance",
    DimensionPolicy::Dimensionless,
    (z1: num_complex::Complex64, z2: num_complex::Complex64) -> Result<f64, OpError>,
    |z1: num_complex::Complex64, z2: num_complex::Complex64| {
        super::hyperbolic::uhp_distance(z1, z2)
            .ok_or_else(|| OpError::Domain("both points must lie inside the open upper half-plane".into()))
    }
);
