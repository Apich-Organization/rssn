//! Separable Q1 Poisson solves with homogeneous Dirichlet boundaries.

use std::f64::consts::PI;
use std::sync::Arc;

use num_complex::Complex64;
use rustfft::Fft;
use rustfft::FftDirection;

use crate::kernels::par::for_chunks_mut;
use crate::kernels::par::map_indexed;
use crate::kernels::reuse::f64_scratch;
use crate::kernels::reuse::fft_plan;
use crate::kernels::reuse::fft_scratch_for;

use super::GaussQuadrature;

const CELL_GRAIN: usize = 256;
const AXIS_GRAIN: usize = 8192;
const LINES_PER_JOB: usize = 8;
const CORNERS_2D: [usize; 4] = [0, 1, 3, 2];
const CORNERS_3D: [usize; 8] = [0, 1, 2, 3, 4, 5, 6, 7];

struct SineAxis {
    extent: usize,
    stride: usize,
    fft: Arc<dyn Fft<f64>>,
    scale: f64,
}

impl SineAxis {
    fn line(
        &self,
        row: &mut [f64],
        stride: usize,
        work: &mut [Complex64],
    ) {
        let (input, scratch) = work.split_at_mut(self.fft.len());
        let half = self.extent + 1;
        input[0] = Complex64::ZERO;
        input[half] = Complex64::ZERO;
        for i in 1..half {
            let v = row[(i - 1) * stride];
            input[i] = Complex64::new(v, 0.0);
            input[2 * half - i] = Complex64::new(-v, 0.0);
        }
        self.fft.process_with_scratch(input, scratch);
        for i in 1..half {
            row[(i - 1) * stride] = -self.scale * input[i].im;
        }
    }

    fn work_len(&self) -> usize {
        self.fft.len() + self.fft.get_inplace_scratch_len()
    }

    fn contiguous(
        &self,
        data: &mut [f64],
    ) {
        let n = self.extent;
        for_chunks_mut(data, n * LINES_PER_JOB, 2, |_, block| {
            let mut work = fft_scratch_for(self.work_len());
            work.resize(self.work_len(), Complex64::ZERO);
            for line in block.chunks_exact_mut(n) {
                self.line(line, 1, &mut work);
            }
        });
    }

    fn apply(
        &self,
        data: &mut [f64],
    ) {
        let n = self.extent;
        let s = self.stride;
        let group = n * s;
        if !crate::execution::should_split(data.len(), AXIS_GRAIN) {
            let mut work = fft_scratch_for(self.work_len());
            work.resize(self.work_len(), Complex64::ZERO);
            for block in data.chunks_exact_mut(group) {
                for offset in 0..s {
                    self.line(&mut block[offset..], s, &mut work);
                }
            }
        } else if s == 1 {
            self.contiguous(data);
        } else {
            // Pack strided lines so parallel FFT tasks have independent mutable slices.
            let mut packed = f64_scratch(data.len());
            packed.resize(data.len(), 0.0);
            for_chunks_mut(&mut packed, n * LINES_PER_JOB, 2, |job, block| {
                let mut work = fft_scratch_for(self.work_len());
                work.resize(self.work_len(), Complex64::ZERO);
                for (i, line) in block.chunks_exact_mut(n).enumerate() {
                    let line_index = job * LINES_PER_JOB + i;
                    let first = (line_index / s) * group + line_index % s;
                    for (k, value) in line.iter_mut().enumerate() {
                        *value = data[first + s * k];
                    }
                    self.line(line, 1, &mut work);
                }
            });
            for_chunks_mut(data, 1024, 2, |offset, block| {
                for (i, value) in block.iter_mut().enumerate() {
                    let flat = offset * 1024 + i;
                    let line = (flat / group) * s + flat % s;
                    let mode = (flat % group) / s;
                    *value = packed[line * n + mode];
                }
            });
        }
    }
}

struct Sample {
    offset: [f64; 3],
    weights: [f64; 8],
}

/// Prepared Q1 finite-element Poisson solver on a uniform Cartesian unit square or cube.
///
/// The one-dimensional mass and stiffness matrices share a discrete sine
/// basis. The plan retains the tensor-product stiffness spectrum and Fourier
/// plans. Repeated solves only assemble the load and apply the sine transforms.
pub struct PreparedCartesianPoisson {
    dim: usize,
    cells: [usize; 3],
    nodes: [usize; 3],
    interior: [usize; 3],
    spacing: [f64; 3],
    samples: Vec<Sample>,
    axes: Vec<SineAxis>,
    spectrum: Vec<f64>,
}

impl PreparedCartesianPoisson {
    /// Prepares a grid with two or three positive element counts.
    /// Coordinates and output storage use x-fast order on the unit domain.
    ///
    /// # Errors
    /// Rejects invalid extents and grids too large for addressable storage.
    pub fn new(elements: &[usize]) -> Result<Self, String> {
        let dim = elements.len();
        if !(2..=3).contains(&dim) || elements.contains(&0) {
            return Err("Cartesian Poisson needs two or three positive element counts".into());
        }
        let mut cells = [1; 3];
        cells[..dim].copy_from_slice(elements);
        let mut nodes = [1; 3];
        let mut interior = [1; 3];
        let mut spacing = [1.0; 3];
        for (axis, &n) in elements.iter().enumerate() {
            nodes[axis] = n
                .checked_add(1)
                .ok_or("Cartesian Poisson extent overflows")?;
            interior[axis] = n - 1;
            spacing[axis] = 1.0 / n as f64;
        }
        let node_count = nodes
            .iter()
            .try_fold(1_usize, |size, &n| size.checked_mul(n));
        if node_count
            .is_none_or(|n| n > isize::MAX as usize / size_of::<Result<[f64; 8], String>>())
        {
            return Err("Cartesian Poisson grid exceeds addressable storage".into());
        }
        let mut axes = Vec::new();
        let mut spectrum = Vec::new();
        if !interior.contains(&0) {
            let mut stride = 1;
            for axis in 0..dim {
                let extent = interior[axis];
                let n = cells[axis];
                if extent > 1 {
                    let fft = fft_plan(2 * n, FftDirection::Forward);
                    axes.push(SineAxis {
                        extent,
                        stride,
                        fft,
                        scale: (0.5 / n as f64).sqrt(),
                    });
                }
                stride *= extent;
            }
            let eigenvalues: [Vec<_>; 3] = std::array::from_fn(|axis| {
                if axis >= dim {
                    return vec![(0.0, 1.0)];
                }
                (1..cells[axis])
                    .map(|mode| {
                        let angle = PI * (mode as f64 / cells[axis] as f64);
                        let half_sine = (0.5 * angle).sin();
                        let stiffness = (4.0 * half_sine * half_sine) / spacing[axis];
                        let mass = spacing[axis] * ((2.0 + angle.cos()) / 3.0);
                        (stiffness, mass)
                    })
                    .collect()
            });
            spectrum.reserve_exact(stride);
            for &(kz, mz) in &eigenvalues[2] {
                for &(ky, my) in &eigenvalues[1] {
                    for &(kx, mx) in &eigenvalues[0] {
                        spectrum.push(kx * my * mz + mx * ky * mz + mx * my * kz);
                    }
                }
            }
        }
        let points = GaussQuadrature::new().points;
        let count = 1 << dim;
        let mut samples = Vec::with_capacity(count);
        for sample in 0..count {
            let coordinate: [f64; 3] = std::array::from_fn(|a| points[(sample >> a) & 1]);
            let mut weights = [0.0; 8];
            for (corner, weight) in weights[..count].iter_mut().enumerate() {
                *weight = (0..dim).fold(1.0 / count as f64, |weight, axis| {
                    let side = if (corner >> axis) & 1 == 0 {
                        -1.0
                    } else {
                        1.0
                    };
                    weight * (1.0 + side * coordinate[axis])
                });
            }
            samples.push(Sample {
                offset: coordinate.map(|x| f64::midpoint(1.0, x)),
                weights,
            });
        }
        Ok(Self {
            dim,
            cells,
            nodes,
            interior,
            spacing,
            samples,
            axes,
            spectrum,
        })
    }

    /// Nodal shape in x, y and optional z order, including the boundary.
    #[must_use]
    pub fn nodal_shape(&self) -> &[usize] {
        &self.nodes[..self.dim]
    }

    pub(crate) const fn heap_bytes(&self) -> usize {
        self.samples.capacity() * size_of::<Sample>()
            + self.axes.capacity() * size_of::<SineAxis>()
            + self.spectrum.capacity() * size_of::<f64>()
    }

    /// Solves for a new source. The callback runs once per element quadrature point.
    ///
    /// # Errors
    /// Reports non-finite source values, cancellation and unrepresentable results.
    pub fn solve<F>(
        &self,
        source: F,
    ) -> Result<Vec<f64>, String>
    where
        F: Fn([f64; 3]) -> f64 + Send + Sync,
    {
        self.solve_source(&source, false)
    }

    pub(super) fn solve_legacy<F>(
        &self,
        source: F,
    ) -> Result<Vec<f64>, String>
    where
        F: Fn([f64; 3]) -> f64 + Send + Sync,
    {
        self.solve_source(&source, true)
    }

    const fn cell(
        &self,
        flat: usize,
    ) -> [usize; 3] {
        let [nx, ny, _] = self.cells;
        [flat % nx, (flat / nx) % ny, flat / (nx * ny)]
    }

    fn assemble_plane_loads(
        &self,
        element_loads: &[Result<[f64; 8], String>],
        rhs: &mut [f64],
    ) -> Result<(), String> {
        // Validate in element order before starting independent nodal sums.
        if let Some(error) = element_loads.iter().find_map(|load| load.as_ref().err()) {
            return Err(error.clone());
        }
        let width = self.interior[0].max(1);
        let columns = self.cells[0];
        let coefficient = |cell: usize, corner| {
            element_loads[cell]
                .as_ref()
                .expect("validated element load")[corner]
        };
        for_chunks_mut(rhs, width, CELL_GRAIN.div_ceil(width), |row, values| {
            let first = row * columns;
            for (x, value) in values.iter_mut().enumerate() {
                let base = first + x;
                let mut sum = 0.0;
                // Each node retains the ascending cell order of the scatter.
                sum += coefficient(base, 3);
                sum += coefficient(base + 1, 2);
                sum += coefficient(base + columns, 1);
                sum += coefficient(base + columns + 1, 0);
                *value = sum;
            }
        });
        Ok(())
    }

    fn solve_source<F>(
        &self,
        source: &F,
        repeated: bool,
    ) -> Result<Vec<f64>, String>
    where
        F: Fn([f64; 3]) -> f64 + Sync,
    {
        let corners = if self.dim == 2 {
            &CORNERS_2D[..]
        } else {
            &CORNERS_3D[..]
        };
        let cell_count = self.cells.iter().product();
        let determinant = self.spacing.iter().product::<f64>() / corners.len() as f64;
        let element_loads: Vec<Result<[f64; 8], String>> =
            map_indexed(cell_count, CELL_GRAIN, |i| {
                crate::cancellation::checkpoint().map_err(|e| e.to_string())?;
                let cell = self.cell(i);
                let mut load = [0.0; 8];
                for sample in &self.samples {
                    let point = std::array::from_fn(|axis| {
                        if axis < self.dim {
                            (cell[axis] as f64 + sample.offset[axis]) * self.spacing[axis]
                        } else {
                            0.0
                        }
                    });
                    let shared = (!repeated).then(|| source(point));
                    for &corner in corners {
                        let f = shared.unwrap_or_else(|| source(point));
                        load[corner] += sample.weights[corner] * f * determinant;
                    }
                }
                if corners.iter().any(|&i| !load[i].is_finite()) {
                    return Err("Cartesian Poisson source produced a non-finite load".into());
                }
                Ok(load)
            });

        let mut rhs = vec![0.0; self.spectrum.len()];
        let [nx, ny, _] = self.interior;
        if self.dim == 2 {
            self.assemble_plane_loads(&element_loads, &mut rhs)?;
        } else {
            for (i, load) in element_loads.into_iter().enumerate() {
                let load = load?;
                let cell = self.cell(i);
                for &corner in corners {
                    let node =
                        std::array::from_fn::<_, 3, _>(|axis| cell[axis] + ((corner >> axis) & 1));
                    if (0..self.dim).all(|a| node[a] > 0 && node[a] < self.cells[a]) {
                        let z = node[2].saturating_sub(1);
                        rhs[(node[0] - 1) + nx * ((node[1] - 1) + ny * z)] += load[corner];
                    }
                }
            }
        }
        let mut output = vec![0.0; self.nodes.iter().product()];
        if rhs.is_empty() {
            return Ok(output);
        }
        let mut peak = 0.0_f64;
        for &r in &rhs {
            if !r.is_finite() {
                return Err("Cartesian Poisson assembled load exceeds Float range".into());
            }
            peak = peak.max(r.abs());
        }
        if peak == 0.0 {
            return Ok(output);
        }
        let scale = if (1e-140..=1e140).contains(&peak) {
            1.0
        } else {
            peak
        };
        if !crate::kernels::exact_eq(scale, 1.0) {
            for value in &mut rhs {
                *value /= scale;
            }
        }
        for axis in &self.axes {
            crate::cancellation::checkpoint().map_err(|e| e.to_string())?;
            axis.apply(&mut rhs);
        }
        for (r, &k) in rhs.iter_mut().zip(&self.spectrum) {
            *r /= k;
        }
        for axis in self.axes.iter().rev() {
            crate::cancellation::checkpoint().map_err(|e| e.to_string())?;
            axis.apply(&mut rhs);
        }
        let [width, height, _] = self.nodes;
        for (row, line) in rhs.chunks_exact(nx).enumerate() {
            let y = row % ny + 1;
            let z = if self.dim == 2 {
                0
            } else {
                row / ny + 1
            };
            let start = 1 + width * (y + height * z);
            for (out, &value) in output[start..start + nx].iter_mut().zip(line) {
                *out = value * scale;
                if !out.is_finite() {
                    return Err("Cartesian Poisson solution exceeds Float range".into());
                }
            }
        }
        Ok(output)
    }
}

#[cfg(test)]
mod tests {
    use super::PreparedCartesianPoisson;

    #[test]
    fn nodal_assembly_retains_incident_element_order() {
        let plan = PreparedCartesianPoisson::new(&[2, 2]).unwrap();
        let mut loads = [[0.0; 8]; 4];
        loads[0][3] = 1.0e30;
        loads[1][2] = -1.0e30;
        loads[2][1] = 1.0;
        loads[3][0] = -0.0;
        let mut rhs = [0.0];
        plan.assemble_plane_loads(&loads.map(Ok), &mut rhs).unwrap();
        assert_eq!(rhs, [1.0]);
    }

    #[test]
    fn nodal_assembly_keeps_the_first_failed_element_on_boundary_only_grids() {
        let plan = PreparedCartesianPoisson::new(&[1, 2]).unwrap();
        let result =
            plan.assemble_plane_loads(&[Err("first".into()), Err("second".into())], &mut []);
        assert_eq!(result, Err("first".into()));
        assert_eq!(
            plan.assemble_plane_loads(&[Ok([0.0; 8]), Ok([0.0; 8])], &mut []),
            Ok(())
        );
    }
}
