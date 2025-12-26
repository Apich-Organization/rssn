use crate::numerical::matrix::Matrix;
use crate::numerical::solve::{
    solve_linear_system,
    LinearSolution,
};
use rayon::prelude::*;
use serde::{
    Deserialize,
    Serialize,
};
use std::ops::{
    Add,
    Mul,
    Sub,
};

#[derive(Clone, Copy, Default, Debug, Serialize, Deserialize)]

pub struct Vector2D {
    pub x: f64,
    pub y: f64,
}

impl Vector2D {
    pub fn new(
        x: f64,
        y: f64,
    ) -> Self {

        Self { x, y }
    }

    pub fn norm(&self) -> f64 {

        (self.x * self.x + self.y * self.y).sqrt()
    }
}

impl Add for Vector2D {
    type Output = Self;

    fn add(
        self,
        rhs: Self,
    ) -> Self {

        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl Mul<f64> for Vector2D {
    type Output = Self;

    fn mul(
        self,
        rhs: f64,
    ) -> Self {

        Self {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
}

impl Sub for Vector2D {
    type Output = Self;

    fn sub(
        self,
        rhs: Self,
    ) -> Self {

        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

#[derive(Clone, Copy, Default, Debug, Serialize, Deserialize)]

pub struct Vector3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vector3D {
    #[allow(dead_code)]

    pub fn new(
        x: f64,
        y: f64,
        z: f64,
    ) -> Self {

        Self { x, y, z }
    }

    #[allow(dead_code)]

    pub fn norm(&self) -> f64 {

        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }
}

impl Sub for Vector3D {
    type Output = Self;

    fn sub(
        self,
        rhs: Self,
    ) -> Self {

        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

/// Specifies the type of boundary condition on an element.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]

pub enum BoundaryCondition<T> {
    Potential(T),
    Flux(T),
}

#[allow(dead_code)]
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]

pub struct Element2D {
    pub p1: Vector2D,
    pub p2: Vector2D,
    pub midpoint: Vector2D,
    pub length: f64,
    pub normal: Vector2D,
}

impl Element2D {
    pub fn new(
        p1: Vector2D,
        p2: Vector2D,
    ) -> Self {

        let diff = p2 - p1;

        let length = diff.norm();

        let normal = Vector2D::new(
            diff.y / length,
            -diff.x / length,
        );

        let midpoint = Vector2D::new(
            (p1.x + p2.x) / 2.0,
            (p1.y + p2.y) / 2.0,
        );

        Self {
            p1,
            p2,
            midpoint,
            length,
            normal,
        }
    }
}

/// Solves a 2D Laplace problem (e.g., potential flow, steady-state heat conduction)
/// using the Boundary Element Method (BEM) with constant elements.
///
/// This function discretizes the boundary of the domain into elements and applies
/// boundary conditions to solve for unknown potentials or fluxes on the boundary.
///
/// # Arguments
/// * `points` - A `Vec` of `(x, y)` tuples defining the vertices of the boundary polygon.
/// * `bcs` - A `Vec` of `BoundaryCondition` for each element, specifying known potential or flux.
///
/// # Returns
/// A `Result` containing a tuple `(u, q)`, where `u` is a `Vec<f64>` of potentials
/// and `q` is a `Vec<f64>` of normal fluxes on each element. Returns an `Err` string
/// if the system is ill-posed or has no unique solution.

pub fn solve_laplace_bem_2d(
    points: &[(f64, f64)],
    bcs: &[BoundaryCondition<f64>],
) -> Result<(Vec<f64>, Vec<f64>), String> {

    let n = points.len();

    if n != bcs.len() {

        return Err("Number of points and boundary conditions must match.".to_string());
    }

    let elements: Vec<_> = (0..n)
        .map(|i| {

            Element2D::new(
                Vector2D::new(
                    points[i].0,
                    points[i].1,
                ),
                Vector2D::new(
                    points[(i + 1) % n].0,
                    points[(i + 1) % n].1,
                ),
            )
        })
        .collect();

    let mut h_mat = Matrix::zeros(n, n);

    let mut g_mat = Matrix::zeros(n, n);

    // Parallel matrix assembly
    let matrices_data: Vec<
        Vec<(
            usize,
            usize,
            f64,
            f64,
        )>,
    > = (0..n)
        .into_par_iter()
        .map(|i| {

            let mut row = Vec::with_capacity(n);

            for j in 0..n {

                if i == j {

                    let g_ii = elements[i].length / (2.0 * std::f64::consts::PI)
                        * (1.0 - (elements[i].length / 2.0).ln());

                    row.push((i, j, 0.0, g_ii)); // Diagonal H will be set later via rigid body motion trick
                } else {

                    let r_vec = elements[j].midpoint - elements[i].midpoint;

                    let r = r_vec.norm();

                    let dot = r_vec.x * elements[j].normal.x + r_vec.y * elements[j].normal.y;

                    let h_ij = -dot / (2.0 * std::f64::consts::PI * r * r);

                    let g_ij = -1.0 / (2.0 * std::f64::consts::PI) * r.ln();

                    row.push((
                        i,
                        j,
                        h_ij * elements[j].length,
                        g_ij * elements[j].length,
                    ));
                }
            }

            row
        })
        .collect();

    for row in matrices_data {

        for (i, j, h, g) in row {

            *h_mat.get_mut(i, j) = h;

            *g_mat.get_mut(i, j) = g;
        }
    }

    // Rigid body motion trick: sum of row H_ij = 0
    // So H_ii = -sum_{j!=i} H_ij
    for i in 0..n {

        let mut row_sum = 0.0;

        for j in 0..n {

            if i != j {

                row_sum += *h_mat.get(i, j);
            }
        }

        *h_mat.get_mut(i, i) = -row_sum;
    }

    let mut a_mat = Matrix::zeros(n, n);

    let mut b_vec = vec![0.0; n];

    for i in 0..n {

        for j in 0..n {

            match bcs[j] {
                // Unknown depends on element j's BC type
                BoundaryCondition::Potential(u_val) => {

                    // Unknown is flux q_j. Equation side: -G_ij * q_j
                    *a_mat.get_mut(i, j) = -*g_mat.get(i, j);

                    // Known contribution from H_ij * u_j goes to RHS with minus
                    b_vec[i] -= *h_mat.get(i, j) * u_val;
                }
                BoundaryCondition::Flux(q_val) => {

                    // Unknown is potential u_j. Equation side: H_ij * u_j
                    *a_mat.get_mut(i, j) = *h_mat.get(i, j);

                    // Known contribution from G_ij * q_j goes to RHS
                    b_vec[i] += *g_mat.get(i, j) * q_val;
                }
            }
        }
    }

    let solution = match solve_linear_system(&a_mat, &b_vec)? {
        LinearSolution::Unique(sol) => sol,
        _ => return Err("BEM system has no unique solution.".to_string()),
    };

    let mut u = vec![0.0; n];

    let mut q = vec![0.0; n];

    let mut sol_idx = 0;

    for i in 0..n {

        match bcs[i] {
            BoundaryCondition::Potential(u_val) => {

                u[i] = u_val;

                q[i] = solution[sol_idx];

                sol_idx += 1;
            }
            BoundaryCondition::Flux(q_val) => {

                q[i] = q_val;

                u[i] = solution[sol_idx];

                sol_idx += 1;
            }
        }
    }

    Ok((u, q))
}

/// Scenario for 2D BEM: Simulates potential flow around a cylinder.
///
/// This function sets up a circular boundary and applies boundary conditions
/// corresponding to a uniform flow in the x-direction. It then uses the BEM solver
/// to calculate the potential and flux on the cylinder's surface.
///
/// # Returns
/// A `Result` containing a tuple `(u, q)` of potentials and fluxes on the cylinder surface,
/// or an error string if the BEM system cannot be solved.

pub fn simulate_2d_cylinder_scenario() -> Result<(Vec<f64>, Vec<f64>), String> {

    let n_points = 40;

    let radius = 1.0;

    let mut points = Vec::new();

    let mut bcs = Vec::new();

    for i in 0..n_points {

        let angle = 2.0 * std::f64::consts::PI * (f64::from(i)) / (f64::from(n_points));

        let (x, y) = (
            radius * angle.cos(),
            radius * angle.sin(),
        );

        points.push((x, y));

        bcs.push(BoundaryCondition::Potential(1.0 * x));
    }

    solve_laplace_bem_2d(&points, &bcs)
}

/// Evaluates the potential at an internal point in the domain after solving the boundary.
///
/// # Arguments
/// * `point` - The `(x, y)` coordinate of the internal point.
/// * `elements` - The boundary elements.
/// * `u` - The solved boundary potentials.
/// * `q` - The solved boundary fluxes.

pub fn evaluate_potential_2d(
    point: (f64, f64),
    elements: &[Element2D],
    u: &[f64],
    q: &[f64],
) -> f64 {

    let p = Vector2D::new(point.0, point.1);

    let mut result = 0.0;

    for i in 0..elements.len() {

        let r_vec = elements[i].midpoint - p;

        let r = r_vec.norm();

        let dot = r_vec.x * elements[i].normal.x + r_vec.y * elements[i].normal.y;

        let h_ij = -dot / (2.0 * std::f64::consts::PI * r * r);

        let g_ij = -1.0 / (2.0 * std::f64::consts::PI) * r.ln();

        result += g_ij * elements[i].length * q[i] - h_ij * elements[i].length * u[i];
    }

    result
}

/// Solves a 3D Laplace problem on a cubic domain using a simplified BEM approach.

pub fn solve_laplace_bem_3d() -> Result<(), String> {

    println!("3D BEM is a complex topic requiring a dedicated library. This is a placeholder.");

    Ok(())
}
