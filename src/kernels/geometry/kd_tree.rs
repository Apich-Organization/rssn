//! Euclidean neighbor search in a k-d tree.

use super::GeometryError;
use std::cmp::Ordering;

// Keep ordinary squared distances on the fast path. Multiplying coordinate
// differences by 2^±600 before squaring covers the remaining Float range.
const SCALE_EXPONENT: u64 = 600;
const GROW: f64 = f64::from_bits((1023 + SCALE_EXPONENT) << 52);
const SHRINK: f64 = f64::from_bits((1023 - SCALE_EXPONENT) << 52);
const NORMAL_FLOOR_GROWN: f64 = f64::MIN_POSITIVE * GROW * GROW;
const NORMAL_CEILING_SHRUNK: f64 = f64::MAX * SHRINK * SHRINK;

// In this interval even the smallest nonzero coordinate difference has
// magnitude at least 2^-502. The square of a difference is at most 2^902,
// which leaves ample room to sum even usize::MAX squared components.
const ORDINARY_MIN: f64 = f64::from_bits((1023 - 450) << 52);
const ORDINARY_MAX: f64 = f64::from_bits((1023 + 450) << 52);

pub(super) fn ordinary_coordinate(value: f64) -> bool {
    let absolute = value.abs();
    absolute == 0.0 || (ORDINARY_MIN..=ORDINARY_MAX).contains(&absolute)
}

trait DistanceKey: Copy {
    const UNBOUNDED: Self;
    fn between(
        left: &[f64],
        right: &[f64],
    ) -> Self;
    fn compare(
        self,
        other: Self,
    ) -> Ordering;
    fn less_than(
        self,
        other: Self,
    ) -> bool;
    fn distance(self) -> f64;

    fn radius(radius: f64) -> Self {
        if radius.is_infinite() {
            Self::UNBOUNDED
        } else {
            Self::between(&[radius], &[0.0])
        }
    }
}

impl DistanceKey for f64 {
    const UNBOUNDED: Self = Self::INFINITY;

    #[inline]
    fn between(
        left: &[f64],
        right: &[f64],
    ) -> Self {
        left.iter().zip(right).map(|(a, b)| (a - b) * (a - b)).sum()
    }

    fn compare(
        self,
        other: Self,
    ) -> Ordering {
        self.total_cmp(&other)
    }

    fn less_than(
        self,
        other: Self,
    ) -> bool {
        self < other
    }

    fn distance(self) -> f64 {
        self.sqrt()
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
enum Range {
    Small,
    Ordinary,
    Large,
    Unbounded,
}

#[derive(Clone, Copy, Debug)]
struct SquaredDistance {
    range: Range,
    squared: f64,
}

impl DistanceKey for SquaredDistance {
    const UNBOUNDED: Self = Self {
        range: Range::Unbounded,
        squared: f64::INFINITY,
    };

    #[inline]
    fn between(
        left: &[f64],
        right: &[f64],
    ) -> Self {
        let squared = f64::between(left, right);
        if squared.is_normal() {
            return Self {
                range: Range::Ordinary,
                squared,
            };
        }
        Self::rescaled(left, right, squared.is_finite())
    }

    fn compare(
        self,
        other: Self,
    ) -> Ordering {
        self.range
            .cmp(&other.range)
            .then_with(|| self.squared.total_cmp(&other.squared))
    }

    fn less_than(
        self,
        other: Self,
    ) -> bool {
        self.compare(other).is_lt()
    }

    fn distance(self) -> f64 {
        let root = self.squared.sqrt();
        match self.range {
            | Range::Small => root * SHRINK,
            | Range::Ordinary | Range::Unbounded => root,
            | Range::Large => root * GROW,
        }
    }
}

impl SquaredDistance {
    #[cold]
    fn rescaled(
        left: &[f64],
        right: &[f64],
        small: bool,
    ) -> Self {
        if small {
            let squared: f64 = left
                .iter()
                .zip(right)
                .map(|(a, b)| {
                    let delta = (a - b) * GROW;
                    delta * delta
                })
                .sum();
            // Rounded subnormal terms may have put the raw sum just below
            // the ordinary range even when the recomputed sum belongs to it.
            if squared >= NORMAL_FLOOR_GROWN {
                Self {
                    range: Range::Ordinary,
                    squared: squared * SHRINK * SHRINK,
                }
            } else {
                Self {
                    range: Range::Small,
                    squared,
                }
            }
        } else {
            let squared: f64 = left
                .iter()
                .zip(right)
                .map(|(a, b)| {
                    let delta = a - b;
                    let scaled = if delta.is_finite() {
                        delta * SHRINK
                    } else {
                        a * SHRINK - b * SHRINK
                    };
                    scaled * scaled
                })
                .sum();
            if squared <= NORMAL_CEILING_SHRUNK {
                Self {
                    range: Range::Ordinary,
                    squared: squared * GROW * GROW,
                }
            } else {
                Self {
                    range: Range::Large,
                    squared,
                }
            }
        }
    }
}

#[derive(Clone, Copy)]
struct Neighbor<D> {
    index: usize,
    distance: D,
}

impl<D: DistanceKey> Neighbor<D> {
    fn into_result(self) -> (usize, f64) {
        (self.index, self.distance.distance())
    }
}

trait Neighbors<D> {
    fn consider(
        &mut self,
        neighbor: Neighbor<D>,
    );
    fn limit(&self) -> D;
}

struct Nearest<D>(Neighbor<D>);

impl<D: DistanceKey> Neighbors<D> for Nearest<D> {
    fn consider(
        &mut self,
        neighbor: Neighbor<D>,
    ) {
        if neighbor.distance.less_than(self.0.distance) {
            self.0 = neighbor;
        }
    }

    fn limit(&self) -> D {
        self.0.distance
    }
}

struct KNearest<D> {
    k: usize,
    found: Vec<Neighbor<D>>,
}

impl<D: DistanceKey> Neighbors<D> for KNearest<D> {
    fn consider(
        &mut self,
        neighbor: Neighbor<D>,
    ) {
        if self.found.len() < self.k || neighbor.distance.less_than(self.limit()) {
            let index = self
                .found
                .binary_search_by(|item| item.distance.compare(neighbor.distance))
                .unwrap_or_else(|insertion| insertion);
            self.found.insert(index, neighbor);
            if self.found.len() > self.k {
                self.found.pop();
            }
        }
    }

    fn limit(&self) -> D {
        if self.found.len() < self.k {
            D::UNBOUNDED
        } else {
            self.found.last().map_or(D::UNBOUNDED, |item| item.distance)
        }
    }
}

struct WithinRadius<D> {
    radius: D,
    found: Vec<Neighbor<D>>,
}

impl<D: DistanceKey> Neighbors<D> for WithinRadius<D> {
    fn consider(
        &mut self,
        neighbor: Neighbor<D>,
    ) {
        if !self.radius.less_than(neighbor.distance) {
            self.found.push(neighbor);
        }
    }

    fn limit(&self) -> D {
        self.radius
    }
}

#[derive(Clone, Debug)]
struct KdNode {
    point_index: usize,
    axis: usize,
    left: Option<usize>,
    right: Option<usize>,
}

/// A k-d tree for exact nearest-neighbor, k-nearest, and radius queries in `R^d`.
#[derive(Clone, Debug)]
pub struct KdTree {
    points: Vec<Vec<f64>>,
    dim: usize,
    root: Option<usize>,
    nodes: Vec<KdNode>,
    ordinary_coordinates: bool,
}

impl KdTree {
    /// Builds a k-d tree from `points`.
    ///
    /// # Errors
    ///
    /// * [`GeometryError::Empty`] if `points` is empty or 0-dimensional.
    /// * [`GeometryError::ShapeMismatch`] if points have differing dimensions.
    /// * [`GeometryError::Domain`] if a coordinate is non-finite.
    pub fn build(points: Vec<Vec<f64>>) -> Result<Self, GeometryError> {
        let (dim, ordinary_coordinates) = Self::validate_defining_points(&points)?;
        let mut indices: Vec<usize> = (0..points.len()).collect();
        let mut nodes = Vec::with_capacity(points.len());
        let root = Self::build_recursive(&points, &mut indices, 0, dim, &mut nodes);
        Ok(Self {
            points,
            dim,
            root,
            nodes,
            ordinary_coordinates,
        })
    }

    /// Validates points stored by the portable model or consumed by the tree.
    pub(crate) fn validate_defining_points(
        points: &[Vec<f64>]
    ) -> Result<(usize, bool), GeometryError> {
        if points.is_empty() {
            return Err(GeometryError::Empty);
        }
        let dim = points[0].len();
        if dim == 0 {
            return Err(GeometryError::Empty);
        }
        for pt in points {
            if pt.len() != dim {
                return Err(GeometryError::ShapeMismatch(format!(
                    "expected dimension {dim}, found {}",
                    pt.len()
                )));
            }
        }
        let mut ordinary_coordinates = true;
        for &coordinate in points.iter().flatten() {
            if !coordinate.is_finite() {
                return Err(GeometryError::Domain(
                    "coordinates must be finite".to_string(),
                ));
            }
            ordinary_coordinates = ordinary_coordinates && ordinary_coordinate(coordinate);
        }
        Ok((dim, ordinary_coordinates))
    }

    fn build_recursive(
        points: &[Vec<f64>],
        indices: &mut [usize],
        depth: usize,
        dim: usize,
        nodes: &mut Vec<KdNode>,
    ) -> Option<usize> {
        if indices.is_empty() {
            return None;
        }
        let axis = depth % dim;
        let mid = indices.len() / 2;
        indices.select_nth_unstable_by(mid, |&a, &b| points[a][axis].total_cmp(&points[b][axis]));
        let point_index = indices[mid];
        let (left_slice, right_with_mid) = indices.split_at_mut(mid);
        let right_slice = &mut right_with_mid[1..];
        let left = Self::build_recursive(points, left_slice, depth + 1, dim, nodes);
        let right = Self::build_recursive(points, right_slice, depth + 1, dim, nodes);
        let node_idx = nodes.len();
        nodes.push(KdNode {
            point_index,
            axis,
            left,
            right,
        });
        Some(node_idx)
    }

    /// Finds the nearest point to `query`, returning `(point_index, euclidean_distance)`.
    ///
    /// Distance is infinite if the mathematical distance exceeds the Float range.
    ///
    /// # Errors
    ///
    /// [`GeometryError::ShapeMismatch`] if `query.len()` does not match the tree dimension.
    /// [`GeometryError::Domain`] if a query coordinate is non-finite.
    pub fn nearest(
        &self,
        query: &[f64],
    ) -> Result<Option<(usize, f64)>, GeometryError> {
        if self.validate_query(query)? {
            Ok(self.nearest_validated::<f64>(query))
        } else {
            Ok(self.nearest_validated::<SquaredDistance>(query))
        }
    }

    /// Finds the nearest point for every query in input order, using the same
    /// selection rule as [`Self::nearest`]. Large batches use the execution
    /// policy for parallel mapping. A query does not allocate its own buffer.
    ///
    /// # Errors
    ///
    /// The error from the first invalid query, as in [`Self::nearest`].
    pub fn nearest_batch(
        &self,
        queries: &[Vec<f64>],
    ) -> Result<Vec<Option<(usize, f64)>>, GeometryError> {
        let mut ordinary = true;
        for query in queries {
            ordinary &= self.validate_query(query)?;
        }
        let nearest = if ordinary {
            Self::nearest_validated::<f64>
        } else {
            Self::nearest_validated::<SquaredDistance>
        };
        Ok(crate::kernels::par::map_slice(queries, 256, |query| {
            nearest(self, query)
        }))
    }

    fn nearest_validated<D: DistanceKey>(
        &self,
        query: &[f64],
    ) -> Option<(usize, f64)> {
        let mut nearest = Nearest(Neighbor {
            index: 0,
            distance: D::UNBOUNDED,
        });
        self.search(self.root?, query, &mut nearest);
        Some(nearest.0.into_result())
    }

    /// Finds up to `k` nearest points to `query`, sorted by ascending distance.
    ///
    /// # Errors
    ///
    /// As in [`Self::nearest`].
    pub fn k_nearest(
        &self,
        query: &[f64],
        k: usize,
    ) -> Result<Vec<(usize, f64)>, GeometryError> {
        let ordinary = self.validate_query(query)?;
        let k = k.min(self.points.len());
        if k == 0 {
            return Ok(Vec::new());
        }
        if ordinary {
            Ok(self.k_nearest_validated::<f64>(query, k))
        } else {
            Ok(self.k_nearest_validated::<SquaredDistance>(query, k))
        }
    }

    fn k_nearest_validated<D: DistanceKey>(
        &self,
        query: &[f64],
        k: usize,
    ) -> Vec<(usize, f64)> {
        let mut nearest = KNearest {
            k,
            found: Vec::<Neighbor<D>>::with_capacity(k + 1),
        };
        if let Some(root) = self.root {
            self.search(root, query, &mut nearest);
        }
        nearest
            .found
            .into_iter()
            .map(Neighbor::into_result)
            .collect()
    }

    /// Finds all points within Euclidean `radius` of `query`, sorted by ascending distance.
    ///
    /// A negative radius returns an empty result. Positive infinity includes every point.
    ///
    /// # Errors
    ///
    /// As in [`Self::nearest`]. A NaN radius returns [`GeometryError::Domain`].
    pub fn radius_search(
        &self,
        query: &[f64],
        radius: f64,
    ) -> Result<Vec<(usize, f64)>, GeometryError> {
        let ordinary = self.validate_query(query)?;
        if radius.is_nan() {
            return Err(GeometryError::Domain("radius must not be NaN".to_string()));
        }
        if radius < 0.0 {
            return Ok(Vec::new());
        }
        if ordinary {
            Ok(self.radius_validated::<f64>(query, radius))
        } else {
            Ok(self.radius_validated::<SquaredDistance>(query, radius))
        }
    }

    fn radius_validated<D: DistanceKey>(
        &self,
        query: &[f64],
        radius: f64,
    ) -> Vec<(usize, f64)> {
        let mut nearest = WithinRadius {
            radius: D::radius(radius),
            found: Vec::new(),
        };
        if let Some(root) = self.root {
            self.search(root, query, &mut nearest);
        }
        nearest.found.sort_by(|a, b| a.distance.compare(b.distance));
        nearest
            .found
            .into_iter()
            .map(Neighbor::into_result)
            .collect()
    }

    fn search<D: DistanceKey>(
        &self,
        node_index: usize,
        query: &[f64],
        neighbors: &mut impl Neighbors<D>,
    ) {
        let node = &self.nodes[node_index];
        let point = &self.points[node.point_index];
        let distance = D::between(query, point);
        neighbors.consider(Neighbor {
            index: node.point_index,
            distance,
        });

        let axis = node.axis;
        let target = query[axis];
        let split = point[axis];
        let (near, far) = if target <= split {
            (node.left, node.right)
        } else {
            (node.right, node.left)
        };
        if let Some(near_index) = near {
            self.search(near_index, query, neighbors);
        }
        if let Some(far_index) = far
            && !neighbors.limit().less_than(D::between(&[target], &[split]))
        {
            self.search(far_index, query, neighbors);
        }
    }

    fn validate_query(
        &self,
        query: &[f64],
    ) -> Result<bool, GeometryError> {
        let mut ordinary = self.ordinary_coordinates;
        for &coordinate in query {
            if !coordinate.is_finite() {
                return Err(GeometryError::Domain(
                    "query coordinates must be finite".to_string(),
                ));
            }
            ordinary = ordinary && ordinary_coordinate(coordinate);
        }
        if query.len() != self.dim {
            return Err(GeometryError::ShapeMismatch(format!(
                "query dimension {} does not match tree dimension {}",
                query.len(),
                self.dim,
            )));
        }
        Ok(ordinary)
    }
}
