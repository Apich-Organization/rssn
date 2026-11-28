use crate::kernels::complex::Mobius;
use crate::kernels::conformal;

struct MobiusFixedPoints {
    points: Vec<Complex64>,
    all: bool,
}
record!(MobiusFixedPoints { points: Vec<Complex64>, all: bool });

binding!(
    conformal_disk_automorphism,
    "num_conformal_disk_automorphism",
    "numerical_special",
    "Builds a Möbius automorphism of the unit disk with |alpha| < 1 and angle theta in radians.",
    "crate::kernels::conformal::disk_automorphism",
    DimensionPolicy::Dimensionless,
    (alpha: num_complex::Complex64, theta: f64) -> Result<super::Mobius, OpError>,
    |alpha: num_complex::Complex64, theta: f64| { super::conformal::disk_automorphism(alpha, theta).map_err(OpError::Domain) }
);

binding!(
    conformal_mobius_class,
    "num_conformal_classify_mobius",
    "numerical_special",
    "Classifies finite, projectively equivalent Möbius coefficients as identity, elliptic, parabolic, hyperbolic or loxodromic.",
    "crate::kernels::conformal::classify_mobius",
    DimensionPolicy::Dimensionless,
    (map: super::Mobius) -> Result<String, OpError>,
    |map: super::Mobius| {
        use super::conformal::MobiusClass;
        let class = super::conformal::classify_mobius(&map).map_err(OpError::Domain)?;
        Ok(match class {
            MobiusClass::Identity => "identity",
            MobiusClass::Elliptic => "elliptic",
            MobiusClass::Parabolic => "parabolic",
            MobiusClass::Hyperbolic => "hyperbolic",
            MobiusClass::Loxodromic => "loxodromic",
        }.into())
    }
);

binding!(
    conformal_fixed_points,
    "num_conformal_mobius_fixed_points",
    "numerical_special",
    "Returns finite complex fixed points in the points field. For the identity all is true and points is empty.",
    "crate::kernels::conformal::mobius_fixed_points",
    DimensionPolicy::Dimensionless,
    (map: super::Mobius) -> Result<super::MobiusFixedPoints, OpError>,
    |map: super::Mobius| {
        let class = super::conformal::classify_mobius(&map).map_err(OpError::Domain)?;
        Ok(super::MobiusFixedPoints {
            points: super::conformal::mobius_fixed_points(&map),
            all: class == super::conformal::MobiusClass::Identity,
        })
    }
);

binding!(
    conformal_bergman_disk,
    "num_conformal_bergman_kernel_disk",
    "numerical_special",
    "Bergman kernel of the unit disk. Both dimensionless complex inputs must lie inside the open disk.",
    "crate::kernels::conformal::bergman_kernel_disk",
    DimensionPolicy::Dimensionless,
    (z: num_complex::Complex64, w: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |z: num_complex::Complex64, w: num_complex::Complex64| {
        if z.norm() >= 1.0 || w.norm() >= 1.0 {
            return Err(OpError::Domain("Bergman kernel requires two points in the open unit disk".into()));
        }
        Ok(super::conformal::bergman_kernel_disk(z, w))
    }
);
