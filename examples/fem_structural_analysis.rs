use std::time::Instant;

#[cfg(feature = "output")]
use ndarray::Array2;
#[cfg(feature = "output")]
use plotters::prelude::*;
use rssn::numerical::matrix::Matrix;
use rssn::numerical::physics_fea::Material;
use rssn::numerical::physics_fea::TriangleElement2D;
use rssn::numerical::physics_fea::assemble_2d_stiffness_matrix;
use rssn::numerical::physics_fea::compute_element_strain;
use rssn::numerical::physics_fea::create_rectangular_mesh;
use rssn::numerical::physics_fea::solve_static_structural;
#[cfg(feature = "output")]
use rssn::output::plotting::PlotConfig;
#[cfg(feature = "output")]
use rssn::output::plotting::plot_heatmap_2d;
#[cfg(feature = "output")]
use rssn::output::plotting::plot_surface_2d;

fn main() {

    println!(
        "Starting FEM Structural \
         Analysis Example (Bridge \
         Structure)..."
    );

    // 1. Define Bridge Constants
    const LENGTH: f64 = 10.0;

    const HEIGHT: f64 = 2.0;

    const NX: usize = 40; // Elements in X
    const NY: usize = 8; // Elements in Y
    const LOAD_FORCE: f64 = -1e5; // Downward force N

    // 2. Material Properties (Steel)
    let material = Material::steel();

    println!(
        "Material: Steel (E={:.2e}, \
         v={:.2})",
        material.youngs_modulus,
        material.poissons_ratio
    );

    // 3. Generate Mesh
    println!(
        "Generating mesh ({}x{} \
         elements)...",
        NX, NY
    );

    let (nodes, raw_elements) =
        create_rectangular_mesh(
            LENGTH,
            HEIGHT,
            NX,
            NY,
        );

    // Convert raw elements to TriangleElement2D structs
    let elements: Vec<
        TriangleElement2D,
    > = raw_elements
        .iter()
        .map(|node_indices| {

            let n1 = node_indices[0];

            let n2 = node_indices[1];

            let n3 = node_indices[2];

            let coords = [
                (
                    nodes[n1].x,
                    nodes[n1].y,
                ),
                (
                    nodes[n2].x,
                    nodes[n2].y,
                ),
                (
                    nodes[n3].x,
                    nodes[n3].y,
                ),
            ];

            TriangleElement2D::new(
                *node_indices,
                coords,
                0.1, // Thickness
                material,
                true, // Plane stress
            )
        })
        .collect();

    // 4. Assemble Global Stiffness Matrix
    println!(
        "Assembling stiffness matrix \
         (Nodes: {}, DOFs: {})...",
        nodes.len(),
        nodes.len() * 2
    );

    let start_assemble = Instant::now();

    // Pre-calculate local K matrices for assembly
    let element_matrices: Vec<(
        Matrix<f64>,
        [usize; 6],
    )> = elements
        .iter()
        .map(|el| {

            let k_local = el
                .local_stiffness_matrix(
                );

            let dof_map = [
                el.nodes[0] * 2,
                el.nodes[0] * 2 + 1,
                el.nodes[1] * 2,
                el.nodes[1] * 2 + 1,
                el.nodes[2] * 2,
                el.nodes[2] * 2 + 1,
            ];

            (k_local, dof_map)
        })
        .collect();

    let global_k =
        assemble_2d_stiffness_matrix(
            nodes.len() * 2,
            &element_matrices,
        );

    println!(
        "Assembly completed in {:?}",
        start_assemble.elapsed()
    );

    // 5. Apply Boundary Conditions and Loads
    // Fix left and right bottom corners (Pinned support)
    let mut fixed_dofs = Vec::new();

    let num_nodes_x = NX + 1;

    // Bottom-left corner nodes (fix a small region for stability)
    for i in 0 ..= 2 {

        let node_idx = i;

        fixed_dofs
            .push((node_idx * 2, 0.0)); // Fix X
        fixed_dofs.push((
            node_idx * 2 + 1,
            0.0,
        )); // Fix Y
    }

    // Bottom-right corner nodes
    for i in
        (num_nodes_x - 3) .. num_nodes_x
    {

        let node_idx = i;

        fixed_dofs
            .push((node_idx * 2, 0.0)); // Fix X
        fixed_dofs.push((
            node_idx * 2 + 1,
            0.0,
        )); // Fix Y
    }

    // Apply Uniform Distributed Load on Top Surface
    let mut forces =
        vec![0.0; nodes.len() * 2];

    for i in 0 ..= NX {

        // Top row nodes
        let node_idx =
            NY * num_nodes_x + i;

        forces[node_idx * 2 + 1] =
            LOAD_FORCE
                / (NX as f64 + 1.0);
    }

    // 6. Solve System
    println!(
        "Solving system of linear \
         equations..."
    );

    let start_solve = Instant::now();

    let displacements =
        solve_static_structural(
            global_k,
            forces,
            &fixed_dofs,
        )
        .unwrap();

    println!(
        "Solved in {:?}",
        start_solve.elapsed()
    );

    // 7. Post-Processing & Visualization
    #[cfg(feature = "output")]
    {

        println!("Generatin plots...");

        // Reconstruct displacement field grid for plotting (approximated)
        // We'll create a 2D array of displacement magnitude

        let mut disp_magnitude =
            Array2::zeros((
                NY + 1,
                NX + 1,
            ));

        let mut von_mises =
            Array2::zeros((NY, NX)); // Element-based

        let mut max_disp = 0.0;

        // Node displacements
        for j in 0 ..= NY {

            for i in 0 ..= NX {

                let node_id =
                    j * num_nodes_x + i;

                let ux = displacements
                    [node_id * 2];

                let uy = displacements
                    [node_id * 2 + 1];

                let mag = (ux * ux
                    + uy * uy)
                    .sqrt();

                disp_magnitude
                    [[NY - j, i]] = mag; // Note: Y-axis flip for image coord system

                if mag > max_disp {

                    max_disp = mag;
                }
            }
        }

        // Element stresses (Von Mises) - taking simple average of elements in a grid cell (approx)
        // We iterate raw elements pair-wise as they form rectangles
        for j in 0 .. NY {

            for i in 0 .. NX {

                // Get the two triangles for this grid cell
                let el_idx1 =
                    (j * NX + i) * 2;

                let el_idx2 =
                    el_idx1 + 1;

                let s1 = compute_element_stress(&elements[el_idx1], &displacements);

                let vm1 = TriangleElement2D::von_mises_stress(&s1);

                let s2 = compute_element_stress(&elements[el_idx2], &displacements);

                let vm2 = TriangleElement2D::von_mises_stress(&s2);

                von_mises
                    [[NY - 1 - j, i]] =
                    (vm1 + vm2) / 2.0;
            }
        }

        // Plot Displacement Magnitude
        let mut plot_config =
            PlotConfig::default();

        plot_config.caption =
            "Structural Displacement \
             Magnitude"
                .to_string();

        plot_config.width = 3840;

        plot_config.height = 2160; // Aspect ratio of bridge

        // Heatmap of displacement
        let path_disp =
            "fem_bridge_displacement.\
             png";

        if let Err(e) = plot_heatmap_2d(
            &disp_magnitude,
            path_disp,
            Some(plot_config.clone()),
        ) {

            eprintln!(
                "Error plotting \
                 displacement: {}",
                e
            );
        } else {

            println!(
                "Saved {}",
                path_disp
            );
        }

        // Plot Stress
        plot_config.caption =
            "Von Mises Stress \
             Distribution"
                .to_string();

        let path_stress =
            "fem_bridge_stress.png";

        if let Err(e) = plot_heatmap_2d(
            &von_mises,
            path_stress,
            Some(plot_config.clone()),
        ) {

            eprintln!(
                "Error plotting \
                 stress: {}",
                e
            );
        } else {

            println!(
                "Saved {}",
                path_stress
            );
        }

        // 3D Visual of Displacement (Fun exaggeration)
        plot_config.caption =
            "3D Displacement \
             Visualization \
             (Exaggerated)"
                .to_string();

        plot_config.height = 800;

        let path_3d =
            "fem_bridge_3d_disp.png";

        // Scale up values for visibility in 3D plot
        let disp_viz = disp_magnitude
            .mapv(|v| {
                v * 10.0 / max_disp
            });

        if let Err(e) = plot_surface_2d(
            &disp_viz,
            path_3d,
            Some(plot_config),
        ) {

            eprintln!(
                "Error plotting 3D \
                 surface: {}",
                e
            );
        } else {

            println!(
                "Saved {}",
                path_3d
            );
        }
    }

    #[cfg(not(feature = "output"))]

    println!(
        "Output disabled. Max node \
         displacement: {:?}",
        displacements
            .iter()
            .fold(0.0f64, |a, &b| {
                a.max(b.abs())
            })
    );
}

fn compute_element_stress(
    element: &TriangleElement2D,
    global_displacements: &[f64],
) -> Vec<f64> {

    let mut el_disp =
        Vec::with_capacity(6);

    for &node_idx in &element.nodes {

        el_disp.push(
            global_displacements
                [node_idx * 2],
        );

        el_disp.push(
            global_displacements
                [node_idx * 2 + 1],
        );
    }

    let b_mat = element.b_matrix();

    let strain = compute_element_strain(
        &b_mat,
        &el_disp,
    );

    let d_mat =
        element.constitutive_matrix();

    // Stress = D * Strain
    // Custom mult for 3x3 * 3x1
    let mut stress = vec![0.0; 3];

    for i in 0 .. 3 {

        for j in 0 .. 3 {

            stress[i] += d_mat
                .get(i, j)
                * strain[j];
        }
    }

    stress
}
