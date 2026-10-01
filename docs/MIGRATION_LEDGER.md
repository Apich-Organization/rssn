# Migration ledger

Every public function of the pre-rewrite tree (commit `84f3ee89`), and where its
functionality lives now. A row is closed only when the new home is named and a test
exercises it, or when it is dropped with a reason. Status values: `pending`, `done`,
`partial` (say what is missing), `dropped` (say why).

Read a legacy implementation with `git show 84f3ee89:<file>`.

## `src/compute/config.rs` (19)

| legacy function | new home | status |
|---|---|---|
| `new` | replaced by `Config::new`/`Config::with` (rule-set selection) in src/api.rs | dropped |
| `with_calculus` | replaced by `Config::new`/`Config::with` (rule-set selection) in src/api.rs | dropped |
| `with_ode` | replaced by `Config::new`/`Config::with` (rule-set selection) in src/api.rs | dropped |
| `with_pde` | replaced by `Config::new`/`Config::with` (rule-set selection) in src/api.rs | dropped |
| `with_algebra` | replaced by `Config::new`/`Config::with` (rule-set selection) in src/api.rs | dropped |
| `with_matrix` | replaced by `Config::new`/`Config::with` (rule-set selection) in src/api.rs | dropped |
| `with_transforms` | replaced by `Config::new`/`Config::with` (rule-set selection) in src/api.rs | dropped |
| `with_special_functions` | replaced by `Config::new`/`Config::with` (rule-set selection) in src/api.rs | dropped |
| `with_optimization` | replaced by `Config::new`/`Config::with` (rule-set selection) in src/api.rs | dropped |
| `with_stats` | replaced by `Config::new`/`Config::with` (rule-set selection) in src/api.rs | dropped |
| `with_physics` | replaced by `Config::new`/`Config::with` (rule-set selection) in src/api.rs | dropped |
| `with_all` | replaced by `Config::new`/`Config::with` (rule-set selection) in src/api.rs | dropped |
| `target_symbolic` | replaced by `api::Target` in src/api.rs | dropped |
| `target_numerical` | replaced by `api::Target` in src/api.rs | dropped |
| `bind` | replaced by `Config::bind` in src/api.rs | dropped |
| `with_budget` | replaced by `Config::budget` in src/api.rs | dropped |
| `with_initial_condition` | replaced by the conditions argument of `dsolve` (rules::ode) | dropped |
| `with_ode_range` | replaced by the arguments of `odeint` (rules::ode) | dropped |
| `with_ode_steps` | replaced by the arguments of `odeint` (rules::ode) | dropped |

## `src/compute/operators.rs` (14)

| legacy function | new home | status |
|---|---|---|
| `d` | `api::diff` / `rules::calculus` `diff` (test derivatives_symbolic) | done |
| `integral` | Expr request builders; requests are built by `Session::call` / `Term::apply` and registered by the rule sets | dropped |
| `definite_integral` | Expr request builders; requests are built by `Session::call` / `Term::apply` and registered by the rule sets | dropped |
| `indefinite_integral` | Expr request builders; requests are built by `Session::call` / `Term::apply` and registered by the rule sets | dropped |
| `limit` | Expr request builders; requests are built by `Session::call` / `Term::apply` and registered by the rule sets | dropped |
| `gradient` | Expr request builders; requests are built by `Session::call` / `Term::apply` and registered by the rule sets | dropped |
| `var` | Expr request builders; requests are built by `Session::call` / `Term::apply` and registered by the rule sets | dropped |
| `ode` | Expr request builders; requests are built by `Session::call` / `Term::apply` and registered by the rule sets | dropped |
| `solve` | Expr request builders; requests are built by `Session::call` / `Term::apply` and registered by the rule sets | dropped |
| `eq` | Expr request builders; requests are built by `Session::call` / `Term::apply` and registered by the rule sets | dropped |
| `matrix_inv` | Expr request builders; requests are built by `Session::call` / `Term::apply` and registered by the rule sets | dropped |
| `det` | Expr request builders; requests are built by `Session::call` / `Term::apply` and registered by the rule sets | dropped |
| `eigenvalues` | Expr request builders; requests are built by `Session::call` / `Term::apply` and registered by the rule sets | dropped |
| `charpoly` | Expr request builders; requests are built by `Session::call` / `Term::apply` and registered by the rule sets | dropped |

## `src/input/parser.rs` (1)

| legacy function | new home | status |
|---|---|---|
| `parse_expr` | `Session::parse` (api.rs; used by every rule-set test) | done |

## `src/jit/engine.rs` (7)

| legacy function | new home | status |
|---|---|---|
| `register_memory_region` | old JIT engine; replaced by `backend::Backend`/`Interpreter` (`Term::compile`) | dropped |
| `clear_memory_regions` | old JIT engine; replaced by `backend::Backend`/`Interpreter` (`Term::compile`) | dropped |
| `allow_call_target` | old JIT engine; replaced by `backend::Backend`/`Interpreter` (`Term::compile`) | dropped |
| `build_sandbox_context` | old JIT engine; replaced by `backend::Backend`/`Interpreter` (`Term::compile`) | dropped |
| `new` | old JIT engine; replaced by `backend::Backend`/`Interpreter` (`Term::compile`) | dropped |
| `register_custom_op` | old JIT engine; custom operators are `OpDescriptor::eval` registrations | dropped |
| `compile` | old JIT engine; replaced by `Term::compile` / `backend::Interpreter` (test reduce_then_compile) | dropped |

## `src/nightly/matrix.rs` (32)

| legacy function | new home | status |
|---|---|---|
| `new` | `kernels::matrix::Matrix::new` (tests/kernels/matrix.rs `construction_and_access`) | done |
| `with_backend` | `kernels::matrix::Matrix::with_backend` (no direct test; used by the faer paths) | partial |
| `set_backend` | `kernels::matrix::Matrix::set_backend` (tests/kernels/matrix.rs `set_backend_switches_the_backend_in_place`) | done |
| `zeros` | `kernels::matrix::Matrix::zeros` (tests/kernels/matrix.rs `construction_and_access`) | done |
| `get` | `kernels::matrix::Matrix::get` (tests/kernels/matrix.rs `construction_and_access`) | done |
| `get_mut` | `kernels::matrix::Matrix::get_mut` (tests/kernels/matrix.rs `construction_and_access`) | done |
| `rows` | `kernels::matrix::Matrix::rows` (tests/kernels/matrix.rs `construction_and_access`) | done |
| `cols` | `kernels::matrix::Matrix::cols` (tests/kernels/matrix.rs `construction_and_access`) | done |
| `decompose` | `kernels::matrix::Matrix::decompose` (tests/kernels/matrix.rs `faer_lu_and_qr_reconstruct_the_matrix`) | done |
| `data` | `kernels::matrix::Matrix::data` (tests/kernels/matrix.rs `construction_and_access`) | done |
| `into_data` | `kernels::matrix::Matrix::into_data` (tests/kernels/matrix.rs `construction_and_access`) | done |
| `get_cols` | `kernels::matrix::Matrix::get_cols` (tests/kernels/matrix.rs `construction_and_access`) | done |
| `rref` | `kernels::matrix::Matrix::rref` (tests/kernels/matrix.rs `rref_and_null_space`) | done |
| `transpose` | `kernels::matrix::Matrix::transpose` (tests/kernels/matrix.rs `transpose`) | done |
| `mul_strassen` | `kernels::matrix::Matrix::mul_strassen` (tests/kernels/matrix.rs `strassen_matches_naive_product`) | done |
| `determinant` | `kernels::matrix::Matrix::determinant` (tests/kernels/matrix.rs `determinant`) | done |
| `lu_decomposition` | `kernels::matrix::Matrix::lu_decomposition` (tests/kernels/matrix.rs `lu_decomposition_of_a_2x2_matrix`) | done |
| `determinant_block` | `kernels::matrix::Matrix::determinant_block` (tests/kernels/matrix.rs `block_determinant_matches_lu`) | done |
| `determinant_lu` | `kernels::matrix::Matrix::determinant_lu` (tests/kernels/matrix.rs `block_determinant_matches_lu`) | done |
| `inverse` | `kernels::matrix::Matrix::inverse` (tests/kernels/matrix.rs `inverse_of_3x3_known`) | done |
| `null_space` | `kernels::matrix::Matrix::null_space` (tests/kernels/matrix.rs `rref_and_null_space`) | done |
| `rank` | `kernels::matrix::Matrix::rank` (tests/kernels/matrix.rs `trace_and_rank`) | done |
| `trace` | `kernels::matrix::Matrix::trace` (tests/kernels/matrix.rs `trace_and_rank`) | done |
| `is_symmetric` | `kernels::matrix::Matrix::is_symmetric` (tests/kernels/matrix.rs `identity_orthogonal_symmetric_diagonal`) | done |
| `is_diagonal` | `kernels::matrix::Matrix::is_diagonal` (tests/kernels/matrix.rs `identity_orthogonal_symmetric_diagonal`) | done |
| `frobenius_norm` | `kernels::matrix::Matrix::frobenius_norm` (tests/kernels/matrix.rs `norms`) | done |
| `l1_norm` | `kernels::matrix::Matrix::l1_norm` (tests/kernels/matrix.rs `norms`) | done |
| `linf_norm` | `kernels::matrix::Matrix::linf_norm` (tests/kernels/matrix.rs `norms`) | done |
| `identity` | `kernels::matrix::Matrix::identity` (tests/kernels/matrix.rs `identity_orthogonal_symmetric_diagonal`) | done |
| `is_identity` | `kernels::matrix::Matrix::is_identity` (tests/kernels/matrix.rs `identity_orthogonal_symmetric_diagonal`) | done |
| `is_orthogonal` | `kernels::matrix::Matrix::is_orthogonal` (tests/kernels/matrix.rs `identity_orthogonal_symmetric_diagonal`) | done |
| `jacobi_eigen_decomposition` | `kernels::matrix::Matrix::jacobi_eigen_decomposition` (tests/kernels/matrix.rs `jacobi_eigen_decomposition_of_symmetric_matrix`) | done |

## `src/numerical/calculus.rs` (4)

| legacy function | new home | status |
|---|---|---|
| `partial_derivative` | `kernels::calculus::partial_derivative` | done |
| `gradient` | `kernels::calculus::gradient` | done |
| `jacobian` | `kernels::calculus::jacobian` | done |
| `hessian` | `kernels::calculus::hessian` | done |

## `src/numerical/calculus_of_variations.rs` (2)

| legacy function | new home | status |
|---|---|---|
| `evaluate_action` |  | pending |
| `euler_lagrange` |  | pending |

## `src/numerical/combinatorics.rs` (9)

| legacy function | new home | status |
|---|---|---|
| `factorial` | `kernels::combinatorics::factorial` | done |
| `permutations` | `kernels::combinatorics::permutations` | done |
| `combinations` | `kernels::combinatorics::combinations` | done |
| `solve_recurrence_numerical` | `kernels::combinatorics::solve_recurrence_numerical` | done |
| `stirling_second` | `kernels::combinatorics::stirling_second` | done |
| `bell` | `kernels::combinatorics::bell` | done |
| `catalan` | `kernels::combinatorics::catalan` | done |
| `rising_factorial` | `kernels::combinatorics::rising_factorial` | done |
| `falling_factorial` | `kernels::combinatorics::falling_factorial` | done |

## `src/numerical/complex_analysis.rs` (11)

| legacy function | new home | status |
|---|---|---|
| `contour_integral` | `kernels::complex::contour_integral` (unit test residues_and_contours) | done |
| `residue` | `kernels::complex::residue` (test residues_and_contours) | done |
| `count_zeros_poles` | `kernels::complex::count_zeros_poles` (test argument_principle) | done |
| `complex_derivative` | `kernels::complex::complex_derivative` (test derivatives) | done |
| `new` | `kernels::complex::Mobius` (test mobius) | done |
| `apply` | `kernels::complex::Mobius` (test mobius) | done |
| `compose` | `kernels::complex::Mobius` (test mobius) | done |
| `inverse` | `kernels::complex::Mobius` (test mobius) | done |
| `contour_integral_expr` | `rules::complex::analysis` `contour_integral` (path or circle; numerically checked) | done |
| `residue_expr` | `rules::complex::analysis` `residue` | done |
| `eval_complex_expr` | `Graph::eval_complex` (ComplexEval attributes) | done |

## `src/numerical/computer_graphics.rs` (51)

| legacy function | new home | status |
|---|---|---|
| `new` | `kernels::computer_graphics::Point2D::new` | done |
| `distance_to` | `kernels::computer_graphics::Point2D::distance_to` | done |
| `to_vector` | `kernels::computer_graphics::Point3D::to_vector` | done |
| `magnitude` | `kernels::computer_graphics::Vector2D::magnitude` | done |
| `normalize` | `kernels::computer_graphics::Vector2D::normalize` | done |
| `rotate` | `kernels::computer_graphics::Vector2D::rotate` | done |
| `perpendicular` | `kernels::computer_graphics::Vector2D::perpendicular` | done |
| `magnitude_squared` | `kernels::computer_graphics::Vector3D::magnitude_squared` | done |
| `rgb` | `kernels::computer_graphics::Color::rgb` | done |
| `clamp` | `kernels::computer_graphics::Color::clamp` | done |
| `lerp` | `kernels::computer_graphics::Color::lerp` | done |
| `dot_product_2d` | `kernels::computer_graphics::dot_product_2d` | done |
| `dot_product` | `kernels::computer_graphics::dot_product` | done |
| `cross_product` | `kernels::computer_graphics::cross_product` | done |
| `reflect` | `kernels::computer_graphics::reflect` | done |
| `refract` | `kernels::computer_graphics::refract` | done |
| `slerp` | `kernels::computer_graphics::slerp` | done |
| `angle_between` | `kernels::computer_graphics::angle_between` | done |
| `project` | `kernels::computer_graphics::project` | done |
| `translation_matrix` | `kernels::computer_graphics::translation_matrix` | done |
| `scaling_matrix` | `kernels::computer_graphics::scaling_matrix` | done |
| `uniform_scaling_matrix` | `kernels::computer_graphics::uniform_scaling_matrix` | done |
| `rotation_matrix_x` | `kernels::computer_graphics::rotation_matrix_x` | done |
| `rotation_matrix_y` | `kernels::computer_graphics::rotation_matrix_y` | done |
| `rotation_matrix_z` | `kernels::computer_graphics::rotation_matrix_z` | done |
| `rotation_matrix_axis` | `kernels::computer_graphics::rotation_matrix_axis` | done |
| `shearing_matrix` | `kernels::computer_graphics::shearing_matrix` | done |
| `perspective_matrix` | `kernels::computer_graphics::perspective_matrix` | done |
| `orthographic_matrix` | `kernels::computer_graphics::orthographic_matrix` | done |
| `look_at_matrix` | `kernels::computer_graphics::look_at_matrix` | done |
| `identity_matrix` | `kernels::computer_graphics::identity_matrix` | done |
| `identity` | `kernels::computer_graphics::Quaternion::identity` | done |
| `from_axis_angle` | `kernels::computer_graphics::Quaternion::from_axis_angle` | done |
| `from_euler` | `kernels::computer_graphics::Quaternion::from_euler` | done |
| `conjugate` | `kernels::computer_graphics::Quaternion::conjugate` | done |
| `inverse` | `kernels::computer_graphics::Quaternion::inverse` | done |
| `multiply` | `kernels::computer_graphics::Quaternion::multiply` | done |
| `rotate_vector` | `kernels::computer_graphics::Quaternion::rotate_vector` | done |
| `to_matrix` | `kernels::computer_graphics::Quaternion::to_matrix` | done |
| `at` | `kernels::computer_graphics::Ray::at` | done |
| `ray_sphere_intersection` | `kernels::computer_graphics::ray_sphere_intersection` | done |
| `ray_plane_intersection` | `kernels::computer_graphics::ray_plane_intersection` | done |
| `ray_triangle_intersection` | `kernels::computer_graphics::ray_triangle_intersection` | done |
| `bezier_quadratic` | `kernels::computer_graphics::bezier_quadratic` | done |
| `bezier_cubic` | `kernels::computer_graphics::bezier_cubic` | done |
| `catmull_rom` | `kernels::computer_graphics::catmull_rom` | done |
| `degrees_to_radians` | `kernels::computer_graphics::degrees_to_radians` | done |
| `radians_to_degrees` | `kernels::computer_graphics::radians_to_degrees` | done |
| `transform_point` | `kernels::computer_graphics::transform_point` | done |
| `transform_vector` | `kernels::computer_graphics::transform_vector` | done |
| `barycentric_coordinates` | `kernels::computer_graphics::barycentric_coordinates` | done |

## `src/numerical/convergence.rs` (5)

| legacy function | new home | status |
|---|---|---|
| `sum_series_numerical` | `kernels::series::sum_to_infinity` (unit tests in kernels/series.rs) | done |
| `aitken_acceleration` |  | pending |
| `find_sequence_limit` |  | pending |
| `richardson_extrapolation` | private `richardson` inside `kernels::series::sum_to_infinity` (sums only, doubling checkpoints) | partial |
| `wynn_epsilon` | private `wynn_epsilon` inside `kernels::series::sum_to_infinity`; not callable on its own | partial |

## `src/numerical/coordinates.rs` (3)

| legacy function | new home | status |
|---|---|---|
| `transform_point` |  | pending |
| `numerical_jacobian` |  | pending |
| `transform_point_pure` |  | pending |

## `src/numerical/differential_geometry.rs` (5)

| legacy function | new home | status |
|---|---|---|
| `metric_tensor_at_point` |  | pending |
| `christoffel_symbols` |  | pending |
| `riemann_tensor` |  | pending |
| `ricci_tensor` |  | pending |
| `ricci_scalar` |  | pending |

## `src/numerical/elementary.rs` (25)

| legacy function | new home | status |
|---|---|---|
| `eval_expr` | `Term::eval` (api.rs; test derivatives_agree_across_phases) | done |
| `eval_expr_single` | `Term::eval` (api.rs; test derivatives_agree_across_phases) | done |
| `sin` | `rules::elementary` operator with float semantics (test numeric_phase) | done |
| `cos` | `rules::elementary` operator with float semantics (test numeric_phase) | done |
| `tan` | `rules::elementary` operator with float semantics (test numeric_phase) | done |
| `asin` | `rules::elementary` operator with float semantics (test numeric_phase) | done |
| `acos` | `rules::elementary` operator with float semantics (test numeric_phase) | done |
| `atan` | `rules::elementary` operator with float semantics (test numeric_phase) | done |
| `atan2` | `rules::elementary` operator with float semantics (test numeric_phase) | done |
| `sinh` | `rules::elementary` operator with float semantics (test numeric_phase) | done |
| `cosh` | `rules::elementary` operator with float semantics (test numeric_phase) | done |
| `tanh` | `rules::elementary` operator with float semantics (test numeric_phase) | done |
| `asinh` | `rules::elementary` operator with float semantics (test numeric_phase) | done |
| `acosh` | `rules::elementary` operator with float semantics (test numeric_phase) | done |
| `atanh` | `rules::elementary` operator with float semantics (test numeric_phase) | done |
| `abs` | `rules::elementary` operator with float semantics (test numeric_phase) | done |
| `sqrt` | `rules::elementary` operator with float semantics (test numeric_phase) | done |
| `ln` | `rules::elementary` operator with float semantics (test numeric_phase) | done |
| `log` | `rules::elementary` operator with float semantics (test numeric_phase) | done |
| `exp` | `rules::elementary` operator with float semantics (test numeric_phase) | done |
| `pow` | core `pow` operator (`graph::op`) with float semantics; folded in `rules::arith` (test numbers_fold_exactly) | done |
| `floor` |  | pending |
| `ceil` |  | pending |
| `round` |  | pending |
| `signum` | `rules::special` `sign` (test step_like_values) | done |

## `src/numerical/error_correction.rs` (36)

| legacy function | new home | status |
|---|---|---|
| `new` | `kernels::error_correction::PolyGF256::new` | done |
| `degree` | `kernels::error_correction::PolyGF256::degree` | done |
| `eval` | `kernels::error_correction::PolyGF256::eval` | done |
| `poly_add` | `kernels::error_correction::PolyGF256::poly_add` | done |
| `poly_sub` | `kernels::error_correction::PolyGF256::poly_sub` | done |
| `poly_mul` | `kernels::error_correction::PolyGF256::poly_mul` | done |
| `poly_div` | `kernels::error_correction::PolyGF256::poly_div` | done |
| `derivative` | `kernels::error_correction::PolyGF256::derivative` | done |
| `scale` | `kernels::error_correction::PolyGF256::scale` | done |
| `normalize` | `kernels::error_correction::PolyGF256::normalize` | done |
| `reed_solomon_encode` | `kernels::error_correction::reed_solomon_encode` | done |
| `reed_solomon_decode` | `kernels::error_correction::reed_solomon_decode` | done |
| `reed_solomon_check` | `kernels::error_correction::reed_solomon_check` | done |
| `calculate_syndromes` | `kernels::error_correction::calculate_syndromes` | done |
| `chien_search` | `kernels::error_correction::chien_search` | done |
| `forney_algorithm` | `kernels::error_correction::forney_algorithm` | done |
| `hamming_distance_numerical` | `kernels::error_correction::hamming_distance_numerical` | done |
| `hamming_weight_numerical` | `kernels::error_correction::hamming_weight_numerical` | done |
| `hamming_encode_numerical` | `kernels::error_correction::hamming_encode_numerical` | done |
| `hamming_decode_numerical` | `kernels::error_correction::hamming_decode_numerical` | done |
| `hamming_check_numerical` | `kernels::error_correction::hamming_check_numerical` | done |
| `bch_encode` | `kernels::error_correction::bch_encode` | done |
| `bch_decode` | `kernels::error_correction::bch_decode` | done |
| `crc32_compute_numerical` | `kernels::error_correction::crc32_compute_numerical` | done |
| `crc32_verify_numerical` | `kernels::error_correction::crc32_verify_numerical` | done |
| `crc32_update_numerical` | `kernels::error_correction::crc32_update_numerical` | done |
| `crc32_finalize_numerical` | `kernels::error_correction::crc32_finalize_numerical` | done |
| `crc16_compute` | `kernels::error_correction::crc16_compute` | done |
| `crc8_compute` | `kernels::error_correction::crc8_compute` | done |
| `interleave` | `kernels::error_correction::interleave` | done |
| `deinterleave` | `kernels::error_correction::deinterleave` | done |
| `convolutional_encode` | `kernels::error_correction::convolutional_encode` | done |
| `minimum_distance` | `kernels::error_correction::minimum_distance` | done |
| `code_rate` | `kernels::error_correction::code_rate` | done |
| `error_correction_capability` | `kernels::error_correction::error_correction_capability` | done |
| `error_detection_capability` | `kernels::error_correction::error_detection_capability` | done |

## `src/numerical/finite_field.rs` (8)

| legacy function | new home | status |
|---|---|---|
| `new` | `kernels::finite_field::PrimeFieldElement::new` | done |
| `inverse` | `kernels::finite_field::PrimeFieldElement::inverse` | done |
| `pow` | `kernels::finite_field::PrimeFieldElement::pow` | done |
| `gf256_add` | `kernels::finite_field::gf256_add` | done |
| `gf256_mul` | `kernels::finite_field::gf256_mul` | done |
| `gf256_inv` | `kernels::finite_field::gf256_inv` | done |
| `gf256_div` | `kernels::finite_field::gf256_div` | done |
| `gf256_pow` | `kernels::finite_field::gf256_pow` | done |

## `src/numerical/fractal_geometry_and_chaos.rs` (27)

| legacy function | new home | status |
|---|---|---|
| `new` | `kernels::fractal_geometry_and_chaos::FractalData::new` | done |
| `get` | `kernels::fractal_geometry_and_chaos::FractalData::get` | done |
| `set` | `kernels::fractal_geometry_and_chaos::FractalData::set` | done |
| `generate_mandelbrot_set` | `kernels::fractal_geometry_and_chaos::generate_mandelbrot_set` | done |
| `mandelbrot_escape_time` | `kernels::fractal_geometry_and_chaos::mandelbrot_escape_time` | done |
| `generate_julia_set` | `kernels::fractal_geometry_and_chaos::generate_julia_set` | done |
| `julia_escape_time` | `kernels::fractal_geometry_and_chaos::julia_escape_time` | done |
| `generate_burning_ship` | `kernels::fractal_geometry_and_chaos::generate_burning_ship` | done |
| `generate_multibrot` | `kernels::fractal_geometry_and_chaos::generate_multibrot` | done |
| `generate_newton_fractal` | `kernels::fractal_geometry_and_chaos::generate_newton_fractal` | done |
| `generate_lorenz_attractor` | `kernels::fractal_geometry_and_chaos::generate_lorenz_attractor` | done |
| `generate_lorenz_attractor_custom` | `kernels::fractal_geometry_and_chaos::generate_lorenz_attractor_custom` | done |
| `generate_rossler_attractor` | `kernels::fractal_geometry_and_chaos::generate_rossler_attractor` | done |
| `generate_henon_map` | `kernels::fractal_geometry_and_chaos::generate_henon_map` | done |
| `generate_tinkerbell_map` | `kernels::fractal_geometry_and_chaos::generate_tinkerbell_map` | done |
| `logistic_map_iterate` | `kernels::fractal_geometry_and_chaos::logistic_map_iterate` | done |
| `logistic_bifurcation` | `kernels::fractal_geometry_and_chaos::logistic_bifurcation` | done |
| `lyapunov_exponent_logistic` | `kernels::fractal_geometry_and_chaos::lyapunov_exponent_logistic` | done |
| `lyapunov_exponent_lorenz` | `kernels::fractal_geometry_and_chaos::lyapunov_exponent_lorenz` | done |
| `box_counting_dimension` | `kernels::fractal_geometry_and_chaos::box_counting_dimension` | done |
| `correlation_dimension` | `kernels::fractal_geometry_and_chaos::correlation_dimension` | done |
| `orbit_density` | `kernels::fractal_geometry_and_chaos::orbit_density` | done |
| `orbit_entropy` | `kernels::fractal_geometry_and_chaos::orbit_entropy` | done |
| `apply` | `kernels::fractal_geometry_and_chaos::AffineTransform2D::apply` | done |
| `generate_ifs_fractal` | `kernels::fractal_geometry_and_chaos::generate_ifs_fractal` | done |
| `sierpinski_triangle_ifs` | `kernels::fractal_geometry_and_chaos::sierpinski_triangle_ifs` | done |
| `barnsley_fern_ifs` | `kernels::fractal_geometry_and_chaos::barnsley_fern_ifs` | done |

## `src/numerical/functional_analysis.rs` (8)

| legacy function | new home | status |
|---|---|---|
| `l1_norm` | `kernels::functional_analysis::l1_norm` | done |
| `l2_norm` | `kernels::functional_analysis::l2_norm` | done |
| `infinity_norm` | `kernels::functional_analysis::infinity_norm` | done |
| `inner_product` | `kernels::functional_analysis::inner_product` | done |
| `project` | `kernels::functional_analysis::project` | done |
| `normalize` | `kernels::functional_analysis::normalize` | done |
| `gram_schmidt` | `kernels::functional_analysis::gram_schmidt` | done |
| `gram_schmidt_orthonormal` | `kernels::functional_analysis::gram_schmidt_orthonormal` | done |

## `src/numerical/geometric_algebra.rs` (8)

| legacy function | new home | status |
|---|---|---|
| `new` | `kernels::geometric_algebra::Multivector3D::new` | done |
| `reverse` | `kernels::geometric_algebra::Multivector3D::reverse` | done |
| `conjugate` | `kernels::geometric_algebra::Multivector3D::conjugate` | done |
| `norm_sq` | `kernels::geometric_algebra::Multivector3D::norm_sq` | done |
| `norm` | `kernels::geometric_algebra::Multivector3D::norm` | done |
| `inv` | `kernels::geometric_algebra::Multivector3D::inv` | done |
| `wedge` | `kernels::geometric_algebra::Multivector3D::wedge` | done |
| `dot` | `kernels::geometric_algebra::Multivector3D::dot` | done |

## `src/numerical/graph.rs` (10)

| legacy function | new home | status |
|---|---|---|
| `new` | `kernels::graph::Graph::new` | done |
| `add_edge` | `kernels::graph::Graph::add_edge` | done |
| `num_nodes` | `kernels::graph::Graph::num_nodes` | done |
| `adj` | `kernels::graph::Graph::adj` | done |
| `dijkstra` | `kernels::graph::dijkstra` | done |
| `bfs` | `kernels::graph::bfs` | done |
| `page_rank` | `kernels::graph::page_rank` | done |
| `floyd_warshall` | `kernels::graph::floyd_warshall` | done |
| `connected_components` | `kernels::graph::connected_components` | done |
| `minimum_spanning_tree` | `kernels::graph::minimum_spanning_tree` | done |

## `src/numerical/indefinite_sum.rs` (11)

| legacy function | new home | status |
|---|---|---|
| `try_closed_form_sum` |  | pending |
| `eval_antidiff` |  | pending |
| `eval_normalized` |  | pending |
| `new` |  | pending |
| `eval` |  | pending |
| `eval_indefinite_product_numerical` |  | pending |
| `series_antidiff` |  | pending |
| `compute_taylor_coeffs_numerical` |  | pending |
| `eval_indefinite_sum_numerical` |  | pending |
| `expr_contains_var` |  | pending |
| `extract_linear_coeff` |  | pending |

## `src/numerical/integrate.rs` (6)

| legacy function | new home | status |
|---|---|---|
| `trapezoidal_rule` | `kernels::integrate::trapezoidal_rule` | done |
| `simpson_rule` | `kernels::integrate::simpson_rule` | done |
| `adaptive_quadrature` | `kernels::integrate::adaptive_quadrature` | done |
| `romberg_integration` | `kernels::integrate::romberg_integration` | done |
| `gauss_legendre_quadrature` | `kernels::integrate::gauss_legendre_quadrature` | done |
| `quadrature` | `rules::calculus` `Quadrature` kernel (`defint`, test numeric_quadrature_without_a_closed_form) and `kernels::integrate::gauss_kronrod` | done |

## `src/numerical/interpolate.rs` (4)

| legacy function | new home | status |
|---|---|---|
| `lagrange_interpolation` | `kernels::interpolate::lagrange_interpolation` | done |
| `cubic_spline_interpolation` | `kernels::interpolate::cubic_spline_interpolation` | done |
| `bezier_curve` | `kernels::interpolate::bezier_curve` | done |
| `b_spline` | `kernels::interpolate::b_spline` | done |

## `src/numerical/matrix.rs` (32)

| legacy function | new home | status |
|---|---|---|
| `new` | `kernels::matrix::Matrix::new` | done |
| `with_backend` | `kernels::matrix::Matrix::with_backend` | done |
| `set_backend` | `kernels::matrix::Matrix::set_backend` | done |
| `zeros` | `kernels::matrix::Matrix::zeros` | done |
| `get` | `kernels::matrix::Matrix::get` | done |
| `get_mut` | `kernels::matrix::Matrix::get_mut` | done |
| `rows` | `kernels::matrix::Matrix::rows` | done |
| `cols` | `kernels::matrix::Matrix::cols` | done |
| `decompose` | `kernels::matrix::Matrix::decompose` | done |
| `data` | `kernels::matrix::Matrix::data` | done |
| `into_data` | `kernels::matrix::Matrix::into_data` | done |
| `get_cols` | `kernels::matrix::Matrix::get_cols` | done |
| `rref` | `kernels::matrix::Matrix::rref` | done |
| `transpose` | `kernels::matrix::Matrix::transpose` | done |
| `mul_strassen` | `kernels::matrix::Matrix::mul_strassen` | done |
| `determinant` | `kernels::matrix::Matrix::determinant` | done |
| `lu_decomposition` | `kernels::matrix::Matrix::lu_decomposition` | done |
| `determinant_block` | `kernels::matrix::Matrix::determinant_block` | done |
| `determinant_lu` | `kernels::matrix::Matrix::determinant_lu` | done |
| `inverse` | `kernels::matrix::Matrix::inverse` | done |
| `null_space` | `kernels::matrix::Matrix::null_space` | done |
| `rank` | `kernels::matrix::Matrix::rank` | done |
| `trace` | `kernels::matrix::Matrix::trace` | done |
| `is_symmetric` | `kernels::matrix::Matrix::is_symmetric` | done |
| `is_diagonal` | `kernels::matrix::Matrix::is_diagonal` | done |
| `frobenius_norm` | `kernels::matrix::Matrix::frobenius_norm` | done |
| `l1_norm` | `kernels::matrix::Matrix::l1_norm` | done |
| `linf_norm` | `kernels::matrix::Matrix::linf_norm` | done |
| `identity` | `kernels::matrix::Matrix::identity` | done |
| `is_identity` | `kernels::matrix::Matrix::is_identity` | done |
| `is_orthogonal` | `kernels::matrix::Matrix::is_orthogonal` | done |
| `jacobi_eigen_decomposition` | `kernels::matrix::Matrix::jacobi_eigen_decomposition` | done |

## `src/numerical/multi_valued.rs` (8)

| legacy function | new home | status |
|---|---|---|
| `newton_method_complex` | `kernels::complex::newton` | done |
| `complex_log_k` | `rules::complex::branches` `log_branch` (ComplexEval) | done |
| `complex_sqrt_k` | `rules::complex::branches` `sqrt_branch` (ComplexEval) | done |
| `complex_pow_k` | `rules::complex::branches` `power_branch` (ComplexEval) | done |
| `complex_nth_root_k` | `rules::complex::branches` `root_branch` (ComplexEval) | done |
| `complex_arcsin_k` | `rules::complex::branches` `asin_branch` (ComplexEval) | done |
| `complex_arccos_k` | `rules::complex::branches` `acos_branch` (ComplexEval) | done |
| `complex_arctan_k` | `rules::complex::branches` `atan_branch` (ComplexEval) | done |

## `src/numerical/number_theory.rs` (8)

| legacy function | new home | status |
|---|---|---|
| `gcd` | `kernels::number_theory::gcd` | done |
| `mod_pow` | `kernels::number_theory::mod_pow` | done |
| `mod_inverse` | `kernels::number_theory::mod_inverse` | done |
| `is_prime_miller_rabin` | `kernels::number_theory::is_prime_miller_rabin` | done |
| `lcm` | `kernels::number_theory::lcm` | done |
| `phi` | `kernels::number_theory::phi` | done |
| `factorize` | `kernels::number_theory::factorize` | done |
| `primes_sieve` | `kernels::number_theory::primes_sieve` | done |

## `src/numerical/ode.rs` (5)

| legacy function | new home | status |
|---|---|---|
| `solve_ode_system` | `rules::ode` `odeint` with list arguments (test numeric_integration_agrees_with_closed_forms) / `kernels::ode::solve_adaptive` | done |
| `solve_ode_euler` | `kernels::ode::solve_fixed` with `OdeSolverMethod::Euler` (test euler_exponential) | done |
| `solve_ode_heun` | `kernels::ode::solve_fixed` with `OdeSolverMethod::Heun` (test heun_exponential) | done |
| `solve_ode_system_rk4` | `kernels::ode::solve_fixed` with `OdeSolverMethod::RungeKutta4` (test rk4_harmonic_oscillator_half_period) | done |
| `solve_ode_system_rk4_named` | `rules::ode` `odeint` (symbols are named in the request) / `kernels::ode::solve_fixed` | done |

## `src/numerical/optimize.rs` (8)

| legacy function | new home | status |
|---|---|---|
| `new` | `kernels::optimize::LinearRegression::new` | done |
| `solve_with_gradient_descent` | `kernels::optimize::EquationOptimizer::solve_with_gradient_descent` | done |
| `auto_solve_conjugate_gradient` | `kernels::optimize::EquationOptimizer::auto_solve_conjugate_gradient` | partial (no test: needs a custom argmin `Operator` impl, more than a few lines) |
| `solve_with_bfgs` | `kernels::optimize::EquationOptimizer::solve_with_bfgs` | done |
| `solve_with_pso` | `kernels::optimize::EquationOptimizer::solve_with_pso` | done |
| `auto_solve` | `kernels::optimize::EquationOptimizer::auto_solve` | partial (no test: `problem: P` is the `Array1<f64>` alias, which is never a `CostFunction`, so it cannot be called) |
| `print_optimization_result` | `kernels::optimize::ResultAnalyzer::print_optimization_result` | done |
| `analyze_convergence` | `kernels::optimize::ResultAnalyzer::analyze_convergence` | done |

## `src/numerical/pde.rs` (1)

| legacy function | new home | status |
|---|---|---|
| `pde_solver` | `kernels::pde::pde_solver` | done |

## `src/numerical/physics.rs` (45)

| legacy function | new home | status |
|---|---|---|
| `new` | `sim::classical::Particle3D::new` | done |
| `kinetic_energy` | `sim::classical::Particle3D::kinetic_energy` | done |
| `momentum` | `sim::classical::Particle3D::momentum` | done |
| `simulate_particle_motion` | `sim::classical::simulate_particle_motion` | done |
| `simulate_ising_model` | `sim::classical::simulate_ising_model` | done |
| `solve_1d_schrodinger` | `sim::classical::solve_1d_schrodinger` | done |
| `solve_2d_schrodinger` | `sim::classical::solve_2d_schrodinger` | done |
| `solve_3d_schrodinger` | `sim::classical::solve_3d_schrodinger` | done |
| `solve_heat_equation_1d_crank_nicolson` | `sim::classical::solve_heat_equation_1d_crank_nicolson` | done |
| `solve_wave_equation_1d` | `sim::classical::solve_wave_equation_1d` | done |
| `projectile_motion_with_drag` | `sim::classical::projectile_motion_with_drag` | done |
| `simple_harmonic_oscillator` | `sim::classical::simple_harmonic_oscillator` | done |
| `damped_harmonic_oscillator` | `sim::classical::damped_harmonic_oscillator` | done |
| `simulate_n_body` | `sim::classical::simulate_n_body` | done |
| `gravitational_potential_energy` | `sim::classical::gravitational_potential_energy` | done |
| `total_kinetic_energy` | `sim::classical::total_kinetic_energy` | done |
| `coulomb_force` | `sim::classical::coulomb_force` | done |
| `electric_field_point_charge` | `sim::classical::electric_field_point_charge` | done |
| `electric_potential_point_charge` | `sim::classical::electric_potential_point_charge` | done |
| `magnetic_field_infinite_wire` | `sim::classical::magnetic_field_infinite_wire` | done |
| `lorentz_force` | `sim::classical::lorentz_force` | done |
| `cyclotron_radius` | `sim::classical::cyclotron_radius` | done |
| `ideal_gas_pressure` | `sim::classical::ideal_gas_pressure` | done |
| `ideal_gas_volume` | `sim::classical::ideal_gas_volume` | done |
| `ideal_gas_temperature` | `sim::classical::ideal_gas_temperature` | done |
| `maxwell_boltzmann_speed_distribution` | `sim::classical::maxwell_boltzmann_speed_distribution` | done |
| `maxwell_boltzmann_mean_speed` | `sim::classical::maxwell_boltzmann_mean_speed` | done |
| `maxwell_boltzmann_rms_speed` | `sim::classical::maxwell_boltzmann_rms_speed` | done |
| `blackbody_power` | `sim::classical::blackbody_power` | done |
| `wien_displacement_wavelength` | `sim::classical::wien_displacement_wavelength` | done |
| `lorentz_factor` | `sim::classical::lorentz_factor` | done |
| `time_dilation` | `sim::classical::time_dilation` | done |
| `length_contraction` | `sim::classical::length_contraction` | done |
| `relativistic_momentum` | `sim::classical::relativistic_momentum` | done |
| `relativistic_kinetic_energy` | `sim::classical::relativistic_kinetic_energy` | done |
| `relativistic_total_energy` | `sim::classical::relativistic_total_energy` | done |
| `mass_energy` | `sim::classical::mass_energy` | done |
| `relativistic_velocity_addition` | `sim::classical::relativistic_velocity_addition` | done |
| `quantum_harmonic_oscillator_energy` | `sim::classical::quantum_harmonic_oscillator_energy` | done |
| `hydrogen_energy_level` | `sim::classical::hydrogen_energy_level` | done |
| `de_broglie_wavelength` | `sim::classical::de_broglie_wavelength` | done |
| `heisenberg_position_uncertainty` | `sim::classical::heisenberg_position_uncertainty` | done |
| `photon_energy` | `sim::classical::photon_energy` | done |
| `photon_wavelength` | `sim::classical::photon_wavelength` | done |
| `compton_wavelength` | `sim::classical::compton_wavelength` | done |

## `src/numerical/physics_cfd.rs` (31)

| legacy function | new home | status |
|---|---|---|
| `new` | `sim::physics_cfd::FluidProperties::new` | done |
| `air` | `sim::physics_cfd::FluidProperties::air` | done |
| `water` | `sim::physics_cfd::FluidProperties::water` | done |
| `kinematic_viscosity` | `sim::physics_cfd::FluidProperties::kinematic_viscosity` | done |
| `thermal_diffusivity` | `sim::physics_cfd::FluidProperties::thermal_diffusivity` | done |
| `prandtl_number` | `sim::physics_cfd::FluidProperties::prandtl_number` | done |
| `reynolds_number` | `sim::physics_cfd::reynolds_number` | done |
| `mach_number` | `sim::physics_cfd::mach_number` | done |
| `froude_number` | `sim::physics_cfd::froude_number` | done |
| `cfl_number` | `sim::physics_cfd::cfl_number` | done |
| `check_cfl_stability` | `sim::physics_cfd::check_cfl_stability` | done |
| `diffusion_number` | `sim::physics_cfd::diffusion_number` | done |
| `solve_advection_1d` | `sim::physics_cfd::solve_advection_1d` | done |
| `solve_diffusion_1d` | `sim::physics_cfd::solve_diffusion_1d` | done |
| `solve_poisson_2d_jacobi` | `sim::physics_cfd::solve_poisson_2d_jacobi` | done |
| `solve_poisson_2d_gauss_seidel` | `sim::physics_cfd::solve_poisson_2d_gauss_seidel` | done |
| `solve_poisson_2d_sor` | `sim::physics_cfd::solve_poisson_2d_sor` | done |
| `solve_advection_diffusion_1d` | `sim::physics_cfd::solve_advection_diffusion_1d` | done |
| `solve_burgers_1d` | `sim::physics_cfd::solve_burgers_1d` | done |
| `compute_vorticity` | `sim::physics_cfd::compute_vorticity` | done |
| `compute_stream_function` | `sim::physics_cfd::compute_stream_function` | done |
| `velocity_from_stream_function` | `sim::physics_cfd::velocity_from_stream_function` | done |
| `compute_divergence` | `sim::physics_cfd::compute_divergence` | done |
| `compute_gradient` | `sim::physics_cfd::compute_gradient` | done |
| `compute_laplacian` | `sim::physics_cfd::compute_laplacian` | done |
| `lid_driven_cavity_simple` | `sim::physics_cfd::lid_driven_cavity_simple` | done |
| `apply_dirichlet_bc` | `sim::physics_cfd::apply_dirichlet_bc` | done |
| `apply_neumann_bc` | `sim::physics_cfd::apply_neumann_bc` | done |
| `max_velocity_magnitude` | `sim::physics_cfd::max_velocity_magnitude` | done |
| `l2_norm` | `sim::physics_cfd::l2_norm` | done |
| `max_abs` | `sim::physics_cfd::max_abs` | done |

## `src/numerical/physics_fea.rs` (27)

| legacy function | new home | status |
|---|---|---|
| `new` | `sim::physics_fea::Material::new` | done |
| `steel` | `sim::physics_fea::Material::steel` | done |
| `aluminum` | `sim::physics_fea::Material::aluminum` | done |
| `copper` | `sim::physics_fea::Material::copper` | done |
| `shear_modulus` | `sim::physics_fea::Material::shear_modulus` | done |
| `bulk_modulus` | `sim::physics_fea::Material::bulk_modulus` | done |
| `distance_to` | `sim::physics_fea::Node2D::distance_to` | done |
| `local_stiffness_matrix` | `sim::physics_fea::LinearElement1D::local_stiffness_matrix` | done |
| `assemble_global_stiffness_matrix` | `sim::physics_fea::assemble_global_stiffness_matrix` | done |
| `solve_static_structural` | `sim::physics_fea::solve_static_structural` | done |
| `area` | `sim::physics_fea::TriangleElement2D::area` | done |
| `constitutive_matrix` | `sim::physics_fea::TriangleElement2D::constitutive_matrix` | done |
| `b_matrix` | `sim::physics_fea::TriangleElement2D::b_matrix` | done |
| `compute_stress` | `sim::physics_fea::TriangleElement2D::compute_stress` | done |
| `von_mises_stress` | `sim::physics_fea::TriangleElement2D::von_mises_stress` | done |
| `transformation_matrix` | `sim::physics_fea::BeamElement2D::transformation_matrix` | done |
| `global_stiffness_matrix` | `sim::physics_fea::BeamElement2D::global_stiffness_matrix` | done |
| `mass_matrix` | `sim::physics_fea::BeamElement2D::mass_matrix` | done |
| `conductivity_matrix` | `sim::physics_fea::ThermalElement1D::conductivity_matrix` | done |
| `apply_boundary_conditions_penalty` | `sim::physics_fea::apply_boundary_conditions_penalty` | done |
| `assemble_2d_stiffness_matrix` | `sim::physics_fea::assemble_2d_stiffness_matrix` | done |
| `compute_element_strain` | `sim::physics_fea::compute_element_strain` | done |
| `principal_stresses` | `sim::physics_fea::principal_stresses` | done |
| `max_shear_stress` | `sim::physics_fea::max_shear_stress` | done |
| `safety_factor_von_mises` | `sim::physics_fea::safety_factor_von_mises` | done |
| `create_rectangular_mesh` | `sim::physics_fea::create_rectangular_mesh` | done |
| `refine_mesh` | `sim::physics_fea::refine_mesh` | done |

## `src/numerical/physics_md.rs` (27)

| legacy function | new home | status |
|---|---|---|
| `new` | `sim::physics_md::Particle::new` | done |
| `with_charge` | `sim::physics_md::Particle::with_charge` | done |
| `kinetic_energy` | `sim::physics_md::Particle::kinetic_energy` | done |
| `momentum` | `sim::physics_md::Particle::momentum` | done |
| `speed` | `sim::physics_md::Particle::speed` | done |
| `distance_to` | `sim::physics_md::Particle::distance_to` | done |
| `lennard_jones_interaction` | `sim::physics_md::lennard_jones_interaction` | done |
| `integrate_velocity_verlet` | `sim::physics_md::integrate_velocity_verlet` | done |
| `morse_interaction` | `sim::physics_md::morse_interaction` | done |
| `harmonic_interaction` | `sim::physics_md::harmonic_interaction` | done |
| `coulomb_interaction` | `sim::physics_md::coulomb_interaction` | done |
| `soft_sphere_interaction` | `sim::physics_md::soft_sphere_interaction` | done |
| `total_kinetic_energy` | `sim::physics_md::total_kinetic_energy` | done |
| `total_momentum` | `sim::physics_md::total_momentum` | done |
| `center_of_mass` | `sim::physics_md::center_of_mass` | done |
| `temperature` | `sim::physics_md::temperature` | done |
| `pressure` | `sim::physics_md::pressure` | done |
| `remove_com_velocity` | `sim::physics_md::remove_com_velocity` | done |
| `velocity_rescale` | `sim::physics_md::velocity_rescale` | done |
| `berendsen_thermostat` | `sim::physics_md::berendsen_thermostat` | done |
| `apply_pbc` | `sim::physics_md::apply_pbc` | done |
| `minimum_image_distance` | `sim::physics_md::minimum_image_distance` | done |
| `radial_distribution_function` | `sim::physics_md::radial_distribution_function` | done |
| `mean_square_displacement` | `sim::physics_md::mean_square_displacement` | done |
| `initialize_velocities_maxwell_boltzmann` | `sim::physics_md::initialize_velocities_maxwell_boltzmann` | done |
| `create_cubic_lattice` | `sim::physics_md::create_cubic_lattice` | done |
| `create_fcc_lattice` | `sim::physics_md::create_fcc_lattice` | done |

## `src/numerical/polynomial.rs` (9)

| legacy function | new home | status |
|---|---|---|
| `new` | `kernels::polynomial::Polynomial::new` | done |
| `eval` | `kernels::polynomial::Polynomial::eval` | done |
| `find_roots` | `kernels::polynomial::Polynomial::find_roots` | done |
| `derivative` | `kernels::polynomial::Polynomial::derivative` | done |
| `long_division` | `kernels::polynomial::Polynomial::long_division` | done |
| `degree` | `kernels::polynomial::Polynomial::degree` | done |
| `is_zero` | `kernels::polynomial::Polynomial::is_zero` | done |
| `integral` | `kernels::polynomial::Polynomial::integral` | done |
| `div_scalar` | `kernels::polynomial::Polynomial::div_scalar` | done |

## `src/numerical/real_roots.rs` (4)

| legacy function | new home | status |
|---|---|---|
| `sturm_sequence` | `kernels::real_roots::sturm_sequence` | done |
| `isolate_real_roots` | `kernels::real_roots::isolate_real_roots` | done |
| `refine_root_bisection` | `kernels::real_roots::refine_root_bisection` | done |
| `find_roots` | `kernels::real_roots::find_roots` | done |

## `src/numerical/series.rs` (3)

| legacy function | new home | status |
|---|---|---|
| `taylor_coefficients` | `rules::calculus` `taylor` (exact coefficients; test taylor_and_laurent) | done |
| `evaluate_power_series` | `kernels::polynomial::Polynomial::eval` for finite coefficient lists; no centre argument | partial |
| `sum_series` | `kernels::series::sum_range` (unit test in kernels/series.rs) | done |

## `src/numerical/signal.rs` (5)

| legacy function | new home | status |
|---|---|---|
| `fft` | `kernels::signal::fft` | done |
| `convolve` | `kernels::signal::convolve` | done |
| `cross_correlation` | `kernels::signal::cross_correlation` | done |
| `hann_window` | `kernels::signal::hann_window` | done |
| `hamming_window` | `kernels::signal::hamming_window` | done |

## `src/numerical/solve.rs` (5)

| legacy function | new home | status |
|---|---|---|
| `solve_linear_system` | `kernels::solve::solve_linear_system` | done |
| `solve_nonlinear_system` | `kernels::solve::solve_nonlinear_system` | done |
| `solve_root_newton` | `kernels::solve::solve_root_newton` | done |
| `solve_root_bisection` | `kernels::solve::solve_root_bisection` | done |
| `solve_root` | `kernels::solve::solve_root` | done |

## `src/numerical/sparse.rs` (14)

| legacy function | new home | status |
|---|---|---|
| `csr_from_triplets` | `kernels::sparse::csr_from_triplets` | done |
| `sp_mat_vec_mul` | `kernels::sparse::sp_mat_vec_mul` | done |
| `to_csr` | `kernels::sparse::to_csr` | done |
| `to_dense` | `kernels::sparse::to_dense` | done |
| `rank` | `kernels::sparse::rank` | done |
| `transpose` | `kernels::sparse::transpose` | done |
| `trace` | `kernels::sparse::trace` | done |
| `is_symmetric` | `kernels::sparse::is_symmetric` | done |
| `is_diagonal` | `kernels::sparse::is_diagonal` | done |
| `frobenius_norm` | `kernels::sparse::frobenius_norm` | done |
| `l1_norm` | `kernels::sparse::l1_norm` | done |
| `linf_norm` | `kernels::sparse::linf_norm` | done |
| `to_csmat` | `kernels::sparse::SparseMatrixData::to_csmat` | done |
| `solve_conjugate_gradient` | `kernels::sparse::solve_conjugate_gradient` | done |

## `src/numerical/special.rs` (37)

| legacy function | new home | status |
|---|---|---|
| `gamma_numerical` | `kernels::special::gamma_numerical` | done |
| `ln_gamma_numerical` | `kernels::special::ln_gamma_numerical` | done |
| `digamma_numerical` | `kernels::special::digamma_numerical` | done |
| `lower_incomplete_gamma` | `kernels::special::lower_incomplete_gamma` | done |
| `upper_incomplete_gamma` | `kernels::special::upper_incomplete_gamma` | done |
| `regularized_lower_gamma` | `kernels::special::regularized_lower_gamma` | done |
| `regularized_upper_gamma` | `kernels::special::regularized_upper_gamma` | done |
| `beta_numerical` | `kernels::special::beta_numerical` | done |
| `ln_beta_numerical` | `kernels::special::ln_beta_numerical` | done |
| `incomplete_beta` | `kernels::special::incomplete_beta` | done |
| `regularized_beta` | `kernels::special::regularized_beta` | done |
| `erf_numerical` | `kernels::special::erf_numerical` | done |
| `erfc_numerical` | `kernels::special::erfc_numerical` | done |
| `inverse_erf_numerical` | `kernels::special::inverse_erf_numerical` | done |
| `bessel_j0` | `kernels::special::bessel_j0` | done |
| `bessel_j1` | `kernels::special::bessel_j1` | done |
| `bessel_y0` | `kernels::special::bessel_y0` | done |
| `bessel_y1` | `kernels::special::bessel_y1` | done |
| `bessel_i0` | `kernels::special::bessel_i0` | done |
| `bessel_i1` | `kernels::special::bessel_i1` | done |
| `legendre_p` | `kernels::special::legendre_p` | done |
| `chebyshev_t` | `kernels::special::chebyshev_t` | done |
| `chebyshev_u` | `kernels::special::chebyshev_u` | done |
| `hermite_h` | `kernels::special::hermite_h` | done |
| `laguerre_l` | `kernels::special::laguerre_l` | done |
| `factorial` | `kernels::special::factorial` | done |
| `double_factorial` | `kernels::special::double_factorial` | done |
| `binomial` | `kernels::special::binomial` | done |
| `riemann_zeta` | `kernels::special::riemann_zeta` | done |
| `sinc` | `kernels::special::sinc` | done |
| `logit` | `kernels::special::logit` | done |
| `sigmoid` | `kernels::special::sigmoid` | done |
| `softplus` | `kernels::special::softplus` | done |
| `bernoulli_number` | `kernels::special::bernoulli_number` | done |
| `bernoulli_poly` | `kernels::special::bernoulli_poly` | done |
| `hurwitz_zeta` | `kernels::special::hurwitz_zeta` | done |
| `polygamma_numerical` | `kernels::special::polygamma_numerical` | done |

## `src/numerical/stats.rs` (30)

| legacy function | new home | status |
|---|---|---|
| `mean` | `kernels::stats::mean` | done |
| `variance_with_type` | `kernels::stats::variance_with_type` | done |
| `variance` | `kernels::stats::variance` | done |
| `std_dev` | `kernels::stats::std_dev` | done |
| `median` | `kernels::stats::median` | done |
| `percentile` | `kernels::stats::percentile` | done |
| `covariance` | `kernels::stats::covariance` | done |
| `correlation` | `kernels::stats::correlation` | done |
| `pdf` | `kernels::stats::NormalDist::pdf` | done |
| `cdf` | `kernels::stats::NormalDist::cdf` | done |
| `new` | `kernels::stats::UniformDist::new` | done |
| `pmf` | `kernels::stats::BinomialDist::pmf` | done |
| `simple_linear_regression` | `kernels::stats::simple_linear_regression` | done |
| `min` | `kernels::stats::min` | done |
| `max` | `kernels::stats::max` | done |
| `skewness` | `kernels::stats::skewness` | partial (only covered by an #[ignore]d test that exposes a library bug) |
| `kurtosis` | `kernels::stats::kurtosis` | done |
| `one_way_anova` | `kernels::stats::one_way_anova` | done |
| `two_sample_t_test` | `kernels::stats::two_sample_t_test` | done |
| `shannon_entropy` | `kernels::stats::shannon_entropy` | done |
| `geometric_mean` | `kernels::stats::geometric_mean` | done |
| `harmonic_mean` | `kernels::stats::harmonic_mean` | done |
| `range` | `kernels::stats::range` | done |
| `iqr` | `kernels::stats::iqr` | done |
| `z_scores` | `kernels::stats::z_scores` | done |
| `mode` | `kernels::stats::mode` | done |
| `welch_t_test` | `kernels::stats::welch_t_test` | done |
| `chi_squared_test` | `kernels::stats::chi_squared_test` | done |
| `coefficient_of_variation` | `kernels::stats::coefficient_of_variation` | done |
| `standard_error` | `kernels::stats::standard_error` | done |

## `src/numerical/tensor.rs` (7)

| legacy function | new home | status |
|---|---|---|
| `tensordot` | `kernels::tensor::tensordot` | done |
| `outer_product` | `kernels::tensor::outer_product` | done |
| `tensor_vec_mul` | `kernels::tensor::tensor_vec_mul` | done |
| `inner_product` | `kernels::tensor::inner_product` | done |
| `contract` | `kernels::tensor::contract` | done |
| `norm` | `kernels::tensor::norm` | done |
| `to_arrayd` | `kernels::tensor::TensorData::to_arrayd` | done |

## `src/numerical/testing.rs` (8)

| legacy function | new home | status |
|---|---|---|
| `solve` | `rules::solve` `solve` (tests linear_and_quadratic, transcendental_by_inversion) | done |
| `solve_polynomial` | `rules::solve` `solve` (test higher_degree_by_factorisation) | done |
| `extract_polynomial_coeffs` | `rules::poly` `coeff` / `poly::repr::Poly::coefficients_in` (test degree_and_coefficients) | done |
| `solve_transcendental_numerical` | `rules::solve` `nsolve` (test numeric_phase) / `kernels::solve::solve_root` | done |
| `solve_linear_system_numerical` | `kernels::solve::solve_linear_system` (tests/kernels/solve.rs) | done |
| `solve_linear_system_symbolic` | `rules::linalg` `linsolve` (test linear_systems) | done |
| `solve_system` | `rules::solve` `solve` with lists (tests linear_systems, polynomial_systems) | done |
| `solve_nonlinear_system_numerical` | `kernels::solve::solve_nonlinear_system` (tests/kernels/solve.rs) | done |

## `src/numerical/topology.rs` (5)

| legacy function | new home | status |
|---|---|---|
| `find_connected_components` |  | pending |
| `vietoris_rips_complex` |  | pending |
| `betti_numbers_at_radius` |  | pending |
| `compute_persistence` |  | pending |
| `euclidean_distance` |  | pending |

## `src/numerical/transforms.rs` (4)

| legacy function | new home | status |
|---|---|---|
| `fft` | `kernels::transforms::fft` | done |
| `ifft` | `kernels::transforms::ifft` | done |
| `fft_slice` | `kernels::transforms::fft_slice` | done |
| `ifft_slice` | `kernels::transforms::ifft_slice` | done |

## `src/numerical/vector.rs` (18)

| legacy function | new home | status |
|---|---|---|
| `vec_add` | `kernels::vector::vec_add` | done |
| `vec_sub` | `kernels::vector::vec_sub` | done |
| `scalar_mul` | `kernels::vector::scalar_mul` | done |
| `dot_product` | `kernels::vector::dot_product` | done |
| `norm` | `kernels::vector::norm` | done |
| `l1_norm` | `kernels::vector::l1_norm` | done |
| `linf_norm` | `kernels::vector::linf_norm` | done |
| `lp_norm` | `kernels::vector::lp_norm` | done |
| `normalize` | `kernels::vector::normalize` | done |
| `cross_product` | `kernels::vector::cross_product` | done |
| `distance` | `kernels::vector::distance` | done |
| `angle` | `kernels::vector::angle` | done |
| `project` | `kernels::vector::project` | done |
| `reflect` | `kernels::vector::reflect` | done |
| `lerp` | `kernels::vector::lerp` | done |
| `is_orthogonal` | `kernels::vector::is_orthogonal` | done |
| `is_parallel` | `kernels::vector::is_parallel` | done |
| `cosine_similarity` | `kernels::vector::cosine_similarity` | done |

## `src/numerical/vector_calculus.rs` (7)

| legacy function | new home | status |
|---|---|---|
| `gradient` | `rules::linalg` `grad` (test vector_calculus) | done |
| `divergence` | `rules::linalg` `div` (test vector_calculus) | done |
| `divergence_expr` | `rules::linalg` `div` (test vector_calculus) | done |
| `curl` | `rules::linalg` `curl` (test vector_calculus) | done |
| `curl_expr` | `rules::linalg` `curl` (test vector_calculus) | done |
| `laplacian` | `rules::linalg` `laplacian` (test vector_calculus) | done |
| `directional_derivative` | `rules::linalg` `directional` (test vector_calculus) | done |

## `src/output/io.rs` (14)

| legacy function | new home | status |
|---|---|---|
| `write_npy_file` | `io::write_npy_file` | done |
| `read_npy_file` | `io::read_npy_file` | done |
| `write_csv_file` | `io::write_csv_file` | done |
| `read_csv_file` | `io::read_csv_file` | done |
| `write_json_file` | `io::write_json_file` | done |
| `read_json_file` | `io::read_json_file` | done |
| `save_expr_as_npy` | Expr serialisation; numeric arrays use `io::{write,read}_{npy,csv,json}_file` | dropped |
| `save_expr_as_csv` | Expr serialisation; numeric arrays use `io::{write,read}_{npy,csv,json}_file` | dropped |
| `load_csv_as_expr` | Expr serialisation; numeric arrays use `io::{write,read}_{npy,csv,json}_file` | dropped |
| `save_expr_as_json` | Expr serialisation; numeric arrays use `io::{write,read}_{npy,csv,json}_file` | dropped |
| `load_json_as_expr` | Expr serialisation; numeric arrays use `io::{write,read}_{npy,csv,json}_file` | dropped |
| `load_npy_as_expr` | Expr serialisation; numeric arrays use `io::{write,read}_{npy,csv,json}_file` | dropped |
| `save_expr` | Expr serialisation; numeric arrays use `io::{write,read}_{npy,csv,json}_file` | dropped |
| `load_expr` | Expr serialisation; numeric arrays use `io::{write,read}_{npy,csv,json}_file` | dropped |

## `src/output/latex.rs` (3)

| legacy function | new home | status |
|---|---|---|
| `to_latex` |  | pending |
| `to_latex_prec_with_parens` |  | pending |
| `to_greek` |  | pending |

## `src/output/plotting.rs` (9)

| legacy function | new home | status |
|---|---|---|
| `plot_function_2d` |  | pending |
| `plot_series_2d` |  | pending |
| `plot_vector_field_2d` |  | pending |
| `plot_surface_3d` |  | pending |
| `plot_surface_2d` |  | pending |
| `plot_parametric_curve_3d` |  | pending |
| `plot_vector_field_3d` |  | pending |
| `plot_3d_path_from_points` |  | pending |
| `plot_heatmap_2d` |  | pending |

## `src/output/pretty_print.rs` (1)

| legacy function | new home | status |
|---|---|---|
| `pretty_print` | `Graph::display` (graph/print.rs tests; every rule-set test compares printed answers) | done |

## `src/output/typst.rs` (1)

| legacy function | new home | status |
|---|---|---|
| `to_typst` |  | pending |

## `src/physics/physics_bem.rs` (6)

| legacy function | new home | status |
|---|---|---|
| `new` | `sim::physics_bem::Vector2D::new` | done |
| `norm` | `sim::physics_bem::Vector2D::norm` | done |
| `solve_laplace_bem_2d` | `sim::physics_bem::solve_laplace_bem_2d` | done |
| `simulate_2d_cylinder_scenario` | `sim::physics_bem::simulate_2d_cylinder_scenario` | done |
| `evaluate_potential_2d` | `sim::physics_bem::evaluate_potential_2d` | done |
| `solve_laplace_bem_3d` | `sim::physics_bem::solve_laplace_bem_3d` | done |

## `src/physics/physics_cnm.rs` (5)

| legacy function | new home | status |
|---|---|---|
| `solve_schrodinger_1d_cn` | `sim::physics_cnm::solve_schrodinger_1d_cn` | done |
| `solve_heat_equation_1d_cn` | `sim::physics_cnm::solve_heat_equation_1d_cn` | done |
| `simulate_1d_heat_conduction_cn_scenario` | `sim::physics_cnm::simulate_1d_heat_conduction_cn_scenario` | done |
| `solve_heat_equation_2d_cn_adi` | `sim::physics_cnm::solve_heat_equation_2d_cn_adi` | done |
| `simulate_2d_heat_conduction_cn_adi_scenario` | `sim::physics_cnm::simulate_2d_heat_conduction_cn_adi_scenario` | done |

## `src/physics/physics_em.rs` (8)

| legacy function | new home | status |
|---|---|---|
| `solve_forward_euler` | `sim::physics_em::solve_forward_euler` | done |
| `solve_midpoint_euler` | `sim::physics_em::solve_midpoint_euler` | done |
| `solve_heun_euler` | `sim::physics_em::solve_heun_euler` | done |
| `solve_semi_implicit_euler` | `sim::physics_em::solve_semi_implicit_euler` | done |
| `simulate_oscillator_forward_euler_scenario` | `sim::physics_em::simulate_oscillator_forward_euler_scenario` | done |
| `simulate_gravity_semi_implicit_euler_scenario` | `sim::physics_em::simulate_gravity_semi_implicit_euler_scenario` | done |
| `solve_backward_euler_linear` | `sim::physics_em::solve_backward_euler_linear` | done |
| `simulate_stiff_decay_scenario` | `sim::physics_em::simulate_stiff_decay_scenario` | done |

## `src/physics/physics_fdm.rs` (16)

| legacy function | new home | status |
|---|---|---|
| `from_data` | `sim::physics_fdm::FdmGrid::from_data` | done |
| `new` | `sim::physics_fdm::FdmGrid::new` | done |
| `with_value` | `sim::physics_fdm::FdmGrid::with_value` | done |
| `dimensions` | `sim::physics_fdm::FdmGrid::dimensions` | done |
| `as_slice` | `sim::physics_fdm::FdmGrid::as_slice` | done |
| `as_mut_slice` | `sim::physics_fdm::FdmGrid::as_mut_slice` | done |
| `len` | `sim::physics_fdm::FdmGrid::len` | done |
| `is_empty` | `sim::physics_fdm::FdmGrid::is_empty` | done |
| `solve_heat_equation_2d` | `sim::physics_fdm::solve_heat_equation_2d` | done |
| `solve_wave_equation_2d` | `sim::physics_fdm::solve_wave_equation_2d` | done |
| `solve_wave_equation_3d` | `sim::physics_fdm::solve_wave_equation_3d` | done |
| `solve_poisson_2d` | `sim::physics_fdm::solve_poisson_2d` | done |
| `solve_burgers_1d` | `sim::physics_fdm::solve_burgers_1d` | done |
| `solve_advection_diffusion_1d` | `sim::physics_fdm::solve_advection_diffusion_1d` | done |
| `simulate_2d_heat_conduction_scenario` | `sim::physics_fdm::simulate_2d_heat_conduction_scenario` | done |
| `simulate_2d_wave_propagation_scenario` | `sim::physics_fdm::simulate_2d_wave_propagation_scenario` | done |

## `src/physics/physics_fem.rs` (6)

| legacy function | new home | status |
|---|---|---|
| `solve_poisson_1d` | `sim::physics_fem::solve_poisson_1d` | done |
| `simulate_1d_poisson_scenario` | `sim::physics_fem::simulate_1d_poisson_scenario` | done |
| `solve_poisson_2d` | `sim::physics_fem::solve_poisson_2d` | done |
| `simulate_2d_poisson_scenario` | `sim::physics_fem::simulate_2d_poisson_scenario` | done |
| `solve_poisson_3d` | `sim::physics_fem::solve_poisson_3d` | done |
| `simulate_3d_poisson_scenario` | `sim::physics_fem::simulate_3d_poisson_scenario` | done |

## `src/physics/physics_fvm.rs` (13)

| legacy function | new home | status |
|---|---|---|
| `new` | `sim::physics_fvm::Mesh::new` | done |
| `num_cells` | `sim::physics_fvm::Mesh::num_cells` | done |
| `lax_friedrichs_flux` | `sim::physics_fvm::lax_friedrichs_flux` | done |
| `minmod` | `sim::physics_fvm::minmod` | done |
| `van_leer` | `sim::physics_fvm::van_leer` | done |
| `solve_advection_1d` | `sim::physics_fvm::solve_advection_1d` | done |
| `simulate_1d_advection_scenario` | `sim::physics_fvm::simulate_1d_advection_scenario` | done |
| `solve_burgers_1d` | `sim::physics_fvm::solve_burgers_1d` | done |
| `solve_shallow_water_1d` | `sim::physics_fvm::solve_shallow_water_1d` | done |
| `solve_advection_2d` | `sim::physics_fvm::solve_advection_2d` | done |
| `simulate_2d_advection_scenario` | `sim::physics_fvm::simulate_2d_advection_scenario` | done |
| `solve_advection_3d` | `sim::physics_fvm::solve_advection_3d` | done |
| `simulate_3d_advection_scenario` | `sim::physics_fvm::simulate_3d_advection_scenario` | done |

## `src/physics/physics_mm.rs` (6)

| legacy function | new home | status |
|---|---|---|
| `new` | `sim::physics_mm::Vector2D::new` | done |
| `compute_density_pressure` | `sim::physics_mm::SPHSystem::compute_density_pressure` | done |
| `compute_forces` | `sim::physics_mm::SPHSystem::compute_forces` | done |
| `integrate` | `sim::physics_mm::SPHSystem::integrate` | done |
| `update` | `sim::physics_mm::SPHSystem::update` | done |
| `simulate_dam_break_2d_scenario` | `sim::physics_mm::simulate_dam_break_2d_scenario` | done |

## `src/physics/physics_mtm.rs` (4)

| legacy function | new home | status |
|---|---|---|
| `solve_poisson_1d_multigrid` | `sim::physics_mtm::solve_poisson_1d_multigrid` | done |
| `simulate_1d_poisson_multigrid_scenario` | `sim::physics_mtm::simulate_1d_poisson_multigrid_scenario` | done |
| `solve_poisson_2d_multigrid` | `sim::physics_mtm::solve_poisson_2d_multigrid` | done |
| `simulate_2d_poisson_multigrid_scenario` | `sim::physics_mtm::simulate_2d_poisson_multigrid_scenario` | done |

## `src/physics/physics_rkm.rs` (7)

| legacy function | new home | status |
|---|---|---|
| `solve_rk4` | `sim::physics_rkm::solve_rk4` | done |
| `new` | `sim::physics_rkm::DormandPrince54::new` | done |
| `solve` | `sim::physics_rkm::DormandPrince54::solve` | done |
| `simulate_lorenz_attractor_scenario` | `sim::physics_rkm::simulate_lorenz_attractor_scenario` | done |
| `simulate_damped_oscillator_scenario` | `sim::physics_rkm::simulate_damped_oscillator_scenario` | done |
| `simulate_vanderpol_scenario` | `sim::physics_rkm::simulate_vanderpol_scenario` | done |
| `simulate_lotka_volterra_scenario` | `sim::physics_rkm::simulate_lotka_volterra_scenario` | done |

## `src/physics/physics_sim/fdtd_electrodynamics.rs` (2)

| legacy function | new home | status |
|---|---|---|
| `run_fdtd_simulation` | `sim::models::fdtd_electrodynamics::run_fdtd_simulation` | done |
| `simulate_and_save_final_state` | `sim::models::fdtd_electrodynamics::simulate_and_save_final_state` | done |

## `src/physics/physics_sim/geodesic_relativity.rs` (3)

| legacy function | new home | status |
|---|---|---|
| `effective_potential` | `sim::models::geodesic_relativity::GeodesicParameters::effective_potential` | done |
| `run_geodesic_simulation` | `sim::models::geodesic_relativity::run_geodesic_simulation` | done |
| `simulate_black_hole_orbits_scenario` | `sim::models::geodesic_relativity::simulate_black_hole_orbits_scenario` | partial (no test: scenario writes .csv/.npy files into the working directory; underlying run_* function is tested) |

## `src/physics/physics_sim/gpe_superfluidity.rs` (2)

| legacy function | new home | status |
|---|---|---|
| `run_gpe_ground_state_finder` | `sim::models::gpe_superfluidity::run_gpe_ground_state_finder` | done |
| `simulate_bose_einstein_vortex_scenario` | `sim::models::gpe_superfluidity::simulate_bose_einstein_vortex_scenario` | partial (no test: scenario writes .csv/.npy files into the working directory; underlying run_* function is tested) |

## `src/physics/physics_sim/ising_statistical.rs` (2)

| legacy function | new home | status |
|---|---|---|
| `run_ising_simulation` | `sim::models::ising_statistical::run_ising_simulation` | done |
| `simulate_ising_phase_transition_scenario` | `sim::models::ising_statistical::simulate_ising_phase_transition_scenario` | partial (no test: scenario writes .csv/.npy files into the working directory; underlying run_* function is tested) |

## `src/physics/physics_sim/linear_elasticity.rs` (3)

| legacy function | new home | status |
|---|---|---|
| `element_stiffness_matrix` | `sim::models::linear_elasticity::element_stiffness_matrix` | done |
| `run_elasticity_simulation` | `sim::models::linear_elasticity::run_elasticity_simulation` | done |
| `simulate_cantilever_beam_scenario` | `sim::models::linear_elasticity::simulate_cantilever_beam_scenario` | partial (no test: scenario writes .csv/.npy files into the working directory; underlying run_* function is tested) |

## `src/physics/physics_sim/navier_stokes_fluid.rs` (3)

| legacy function | new home | status |
|---|---|---|
| `run_channel_flow` | `sim::models::navier_stokes_fluid::run_channel_flow` | done |
| `run_lid_driven_cavity` | `sim::models::navier_stokes_fluid::run_lid_driven_cavity` | done |
| `simulate_lid_driven_cavity_scenario` | `sim::models::navier_stokes_fluid::simulate_lid_driven_cavity_scenario` | done |

## `src/physics/physics_sim/schrodinger_quantum.rs` (2)

| legacy function | new home | status |
|---|---|---|
| `run_schrodinger_simulation` | `sim::models::schrodinger_quantum::run_schrodinger_simulation` | done |
| `simulate_double_slit_scenario` | `sim::models::schrodinger_quantum::simulate_double_slit_scenario` | partial (no test: scenario writes .csv/.npy files into the working directory; underlying run_* function is tested) |

## `src/physics/physics_sm.rs` (10)

| legacy function | new home | status |
|---|---|---|
| `fft2d` | `sim::physics_sm::fft2d` | done |
| `ifft2d` | `sim::physics_sm::ifft2d` | done |
| `solve_advection_diffusion_1d` | `sim::physics_sm::solve_advection_diffusion_1d` | done |
| `simulate_1d_advection_diffusion_scenario` | `sim::physics_sm::simulate_1d_advection_diffusion_scenario` | done |
| `solve_advection_diffusion_2d` | `sim::physics_sm::solve_advection_diffusion_2d` | done |
| `simulate_2d_advection_diffusion_scenario` | `sim::physics_sm::simulate_2d_advection_diffusion_scenario` | done |
| `fft3d` | `sim::physics_sm::fft3d` | done |
| `ifft3d` | `sim::physics_sm::ifft3d` | partial (only covered by an #[ignore]d test that exposes a library bug) |
| `solve_advection_diffusion_3d` | `sim::physics_sm::solve_advection_diffusion_3d` | partial (only covered by an #[ignore]d test that exposes a library bug) |
| `simulate_3d_advection_diffusion_scenario` | `sim::physics_sm::simulate_3d_advection_diffusion_scenario` | done |

## `src/plugins/manager.rs` (7)

| legacy function | new home | status |
|---|---|---|
| `empty` | plugin manager (FFI plumbing); extension is by `RuleSet`/`OpDescriptor` | dropped |
| `new` | plugin manager (FFI plumbing); extension is by `RuleSet`/`OpDescriptor` | dropped |
| `execute_plugin` | plugin manager (FFI plumbing); extension is by `RuleSet`/`OpDescriptor` | dropped |
| `register_plugin` | plugin manager (FFI plumbing); extension is by `RuleSet`/`OpDescriptor` | dropped |
| `get_loaded_plugin_names` | plugin manager (FFI plumbing); extension is by `RuleSet`/`OpDescriptor` | dropped |
| `unload_plugin` | plugin manager (FFI plumbing); extension is by `RuleSet`/`OpDescriptor` | dropped |
| `get_plugin_metadata` | plugin manager (FFI plumbing); extension is by `RuleSet`/`OpDescriptor` | dropped |

## `src/plugins/plugin_c.rs` (1)

| legacy function | new home | status |
|---|---|---|
| `new` | plugin manager (FFI plumbing) | dropped |

## `src/symbolic/cad.rs` (1)

| legacy function | new home | status |
|---|---|---|
| `cad` |  | pending |

## `src/symbolic/calculus.rs` (16)

| legacy function | new home | status |
|---|---|---|
| `substitute` | `graph::Graph::substitute` (graph/subst.rs test replaces_free_occurrences) | done |
| `differentiate` | `rules::calculus` `diff` (test derivatives_symbolic) | done |
| `integrate` | `rules::calculus` `integral` (tests table_integrals, substitution, by_parts) | done |
| `integrate_internal` | `rules::calculus` `integral` (tests table_integrals, substitution, by_parts) | done |
| `substitute_expr` | `graph::Graph::substitute` (graph/subst.rs test replaces_free_occurrences) | done |
| `evaluate_at_point` | `Term::eval` with bindings (api.rs) | done |
| `definite_integrate` | `rules::calculus` `defint` (test definite_integrals_agree_across_phases) | done |
| `check_analytic` |  | pending |
| `find_poles` |  | pending |
| `calculate_residue` | `kernels::complex::residue` numeric closure kernel only | partial |
| `is_inside_contour` |  | pending |
| `path_integrate` | `rules::linalg` `line_integral` for real parametrised curves; `kernels::complex::contour_integral` numerically | partial |
| `factorial` | `rules::combinatorics` `factorial` (test factorials) | done |
| `improper_integral` | `defint` over `oo` limits, numeric only via `Quadrature` (test numeric_quadrature_without_a_closed_form); no residue-theorem symbolic result | partial |
| `limit` | `rules::calculus` `limit` (tests limits_by_continuity_and_cancellation, limits_at_infinity) | done |
| `limit_internal` | `rules::calculus` `limit` (tests limits_by_continuity_and_cancellation, limits_at_infinity) | done |

## `src/symbolic/calculus_of_variations.rs` (5)

| legacy function | new home | status |
|---|---|---|
| `euler_lagrange` |  | pending |
| `euler_lagrange_internal` |  | pending |
| `solve_euler_lagrange` |  | pending |
| `solve_euler_lagrange_internal` |  | pending |
| `hamiltons_principle` |  | pending |

## `src/symbolic/cas_foundations.rs` (12)

| legacy function | new home | status |
|---|---|---|
| `get_term_factors` | Expr helper; terms are normalised by `rules::arith` Collect | dropped |
| `build_expr_from_factors` | Expr helper; terms are normalised by `rules::arith` Collect | dropped |
| `normalize` | `rules::arith` Collect window pass (test like_terms_are_collected) | done |
| `expand` | `rules::poly` `expand` (test expand_is_pinned_against_cheaper_spellings) | done |
| `factorize` | `rules::poly` `factor` (tests factor_univariate, factor_multivariate_pulls_out_content) | done |
| `factorize_internal` | `rules::poly` `factor` (tests factor_univariate, factor_multivariate_pulls_out_content) | done |
| `risch_integrate` | `rules::calculus` `integral` heuristic stages (no Risch algorithm) | partial |
| `grobner_basis` | `rules::poly` `groebner` (test groebner_bases) | done |
| `cylindrical_algebraic_decomposition` |  | pending |
| `simplify_with_relations` |  | pending |
| `simplify_with_relations_internal` |  | pending |
| `normalize_with_relations` |  | pending |

## `src/symbolic/classical_mechanics.rs` (19)

| legacy function | new home | status |
|---|---|---|
| `new` |  | pending |
| `newtons_second_law` |  | pending |
| `momentum` | `sim::classical::momentum` (f64 only; tests/sim/classical.rs test_particle3d_momentum) | partial |
| `kinetic_energy` | `sim::classical::kinetic_energy` (f64 only; test_particle3d_kinetic_energy) | partial |
| `potential_energy_gravity_uniform` |  | pending |
| `potential_energy_gravity_universal` |  | pending |
| `potential_energy_spring` |  | pending |
| `work_constant_force` |  | pending |
| `work_line_integral` |  | pending |
| `power` |  | pending |
| `torque` |  | pending |
| `angular_momentum` |  | pending |
| `centripetal_acceleration` |  | pending |
| `moment_of_inertia_point_mass` |  | pending |
| `rotational_kinetic_energy` |  | pending |
| `lagrangian` |  | pending |
| `hamiltonian` |  | pending |
| `euler_lagrange_equation` |  | pending |
| `poisson_bracket` |  | pending |

## `src/symbolic/combinatorics.rs` (10)

| legacy function | new home | status |
|---|---|---|
| `expand_binomial` | `expand((a+b)^n)` of the poly rules, literal n | partial: symbolic n needs a sum operator |
| `permutations` | `permutations` | done |
| `combinations` | `binomial` | done |
| `solve_recurrence` | `rsolve(eq, a(n), list(init))`: rational/quadratic-irrational roots, repeated roots, polynomial*exponential forcing | done |
| `get_sequence_from_gf` | `gf_coeffs(f, x, n)` | done |
| `apply_inclusion_exclusion` | `inclusion_exclusion` | done |
| `find_period` | `period` (length when no shorter period; legacy gave None) | done |
| `catalan_number` | `catalan` | done |
| `stirling_number_second_kind` | `stirling2` | done |
| `bell_number` | `bell` | done |

## `src/symbolic/complex_analysis.rs` (23)

| legacy function | new home | status |
|---|---|---|
| `new` | `rules::complex::analysis` Möbius maps as matrices / `kernels::complex::Mobius::new` | done |
| `continue_along_path` | `rules::complex::analysis` `continue_along` (disc-by-disc check, Taylor at the end point) | done |
| `get_final_expression` | `rules::complex::analysis` `continue_along` result | done |
| `estimate_radius_of_convergence` | `rules::complex::analysis` `radius_of_convergence` (exact distance to nearest pole) | done |
| `complex_distance` | `rules::complex::analysis` `distance` | done |
| `classify_singularity` | `rules::complex::analysis` `singularity` (regular/removable/pole order/essential) | done |
| `laurent_series` | `rules::calculus` `laurent` (complex expansion points via `rules::complex`) | done |
| `calculate_residue` | `rules::complex::analysis` `residue` (derivative formula, numerically checked) | done |
| `calculate_residue_internal` | `rules::complex::analysis` `residue` | done |
| `contour_integral_residue_theorem` | `rules::complex::analysis` `contour_integral(f, z, circle(c, r))` | done |
| `contour_integral_residue_theorem_internal` | `rules::complex::analysis` `contour_integral` | done |
| `identity` | `identity(2)` / `kernels::complex::Mobius::identity` | done |
| `apply` | `rules::complex::analysis` `mobius_apply` | done |
| `compose` | `rules::complex::analysis` `mobius_compose` | done |
| `inverse` | `rules::complex::analysis` `mobius_inverse` | done |
| `cauchy_integral_formula` | `rules::complex::analysis` `cauchy_integral` | done |
| `cauchy_integral_formula_internal` | `rules::complex::analysis` `cauchy_integral` | done |
| `cauchy_derivative_formula` | `rules::complex::analysis` `cauchy_derivative` | done |
| `cauchy_derivative_formula_internal` | `rules::complex::analysis` `cauchy_derivative` | done |
| `complex_exp` | `rules::complex` parts kernel: exp(a + I b) = e^a (cos b + I sin b) | done |
| `complex_log` | `rules::complex` parts kernel: ln z = ln|z| + I arg z | done |
| `complex_arg` | `rules::complex` `arg` | done |
| `complex_modulus` | `rules::complex` `abs` of complex terms | done |

## `src/symbolic/computer_graphics.rs` (22)

| legacy function | new home | status |
|---|---|---|
| `translation_2d` |  | pending |
| `translation_3d` |  | pending |
| `rotation_2d` |  | pending |
| `rotation_3d_x` |  | pending |
| `rotation_3d_y` |  | pending |
| `rotation_3d_z` |  | pending |
| `scaling_2d` |  | pending |
| `scaling_3d` |  | pending |
| `perspective_projection` |  | pending |
| `orthographic_projection` |  | pending |
| `look_at` |  | pending |
| `evaluate` |  | pending |
| `derivative` |  | pending |
| `split` |  | pending |
| `new` |  | pending |
| `apply_transformation` |  | pending |
| `compute_normals` |  | pending |
| `triangulate` |  | pending |
| `shear_2d` |  | pending |
| `reflection_2d` |  | pending |
| `reflection_3d` |  | pending |
| `rotation_axis_angle` |  | pending |

## `src/symbolic/convergence.rs` (1)

| legacy function | new home | status |
|---|---|---|
| `analyze_convergence` | `rules::calculus` `converges` (ratio test only) | partial |

## `src/symbolic/coordinates.rs` (12)

| legacy function | new home | status |
|---|---|---|
| `transform_point` |  | pending |
| `transform_expression` |  | pending |
| `get_transform_rules` |  | pending |
| `get_to_cartesian_rules` |  | pending |
| `transform_contravariant_vector` |  | pending |
| `transform_covariant_vector` |  | pending |
| `transform_tensor2` |  | pending |
| `symbolic_mat_mat_mul` |  | pending |
| `get_metric_tensor` |  | pending |
| `transform_divergence` |  | pending |
| `transform_curl` |  | pending |
| `transform_gradient` |  | pending |

## `src/symbolic/core/api.rs` (54)

| legacy function | new home | status |
|---|---|---|
| `new_constant` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `new_variable` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `new_bigint` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `new_rational` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `new_pi` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `new_e` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `new_infinity` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `new_negative_infinity` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `new_matrix` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `new_predicate` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `new_forall` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `new_exists` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `new_interval` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `new_derivative` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `new_derivativen` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `new_indefinite_sum` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `new_indefinite_product` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `new_sparse_polynomial` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `new_custom_zero` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `new_custom_string` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `new_custom_arc_three` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `new_custom_arc_four` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `new_custom_arc_five` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `is_dag` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `to_dag` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `to_dag_form` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `to_ast` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `register_dynamic_op` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `get_dynamic_op_properties` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `sin` | `api::{sin,cos,tan,exp,ln,abs,sqrt,asin,acos,atan,sinh,cosh,tanh}` (api.rs `functions!`; test operators_build_canonical_terms) | done |
| `cos` | `api::{sin,cos,tan,exp,ln,abs,sqrt,asin,acos,atan,sinh,cosh,tanh}` (api.rs `functions!`; test operators_build_canonical_terms) | done |
| `tan` | `api::{sin,cos,tan,exp,ln,abs,sqrt,asin,acos,atan,sinh,cosh,tanh}` (api.rs `functions!`; test operators_build_canonical_terms) | done |
| `exp` | `api::{sin,cos,tan,exp,ln,abs,sqrt,asin,acos,atan,sinh,cosh,tanh}` (api.rs `functions!`; test operators_build_canonical_terms) | done |
| `ln` | `api::{sin,cos,tan,exp,ln,abs,sqrt,asin,acos,atan,sinh,cosh,tanh}` (api.rs `functions!`; test operators_build_canonical_terms) | done |
| `log` | `Term::apply("log"/"pow")` / `Term::pow` (api.rs) | done |
| `abs` | `api::{sin,cos,tan,exp,ln,abs,sqrt,asin,acos,atan,sinh,cosh,tanh}` (api.rs `functions!`; test operators_build_canonical_terms) | done |
| `sqrt` | `api::{sin,cos,tan,exp,ln,abs,sqrt,asin,acos,atan,sinh,cosh,tanh}` (api.rs `functions!`; test operators_build_canonical_terms) | done |
| `pow` | `Term::apply("log"/"pow")` / `Term::pow` (api.rs) | done |
| `asin` | `api::{sin,cos,tan,exp,ln,abs,sqrt,asin,acos,atan,sinh,cosh,tanh}` (api.rs `functions!`; test operators_build_canonical_terms) | done |
| `acos` | `api::{sin,cos,tan,exp,ln,abs,sqrt,asin,acos,atan,sinh,cosh,tanh}` (api.rs `functions!`; test operators_build_canonical_terms) | done |
| `atan` | `api::{sin,cos,tan,exp,ln,abs,sqrt,asin,acos,atan,sinh,cosh,tanh}` (api.rs `functions!`; test operators_build_canonical_terms) | done |
| `sinh` | `api::{sin,cos,tan,exp,ln,abs,sqrt,asin,acos,atan,sinh,cosh,tanh}` (api.rs `functions!`; test operators_build_canonical_terms) | done |
| `cosh` | `api::{sin,cos,tan,exp,ln,abs,sqrt,asin,acos,atan,sinh,cosh,tanh}` (api.rs `functions!`; test operators_build_canonical_terms) | done |
| `tanh` | `api::{sin,cos,tan,exp,ln,abs,sqrt,asin,acos,atan,sinh,cosh,tanh}` (api.rs `functions!`; test operators_build_canonical_terms) | done |
| `to_f64` | `Term::as_f64` (api.rs) | done |
| `is_positive` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `is_negative` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `is_fractional` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `is_integer` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `is_float` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `is_zero` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `is_one` | Expr/DAG/FFI plumbing (terms are built through `Session`/`Term`) | dropped |
| `as_f64` | `Term::as_f64` (api.rs) | done |
| `simplify` | `Session::simplify` (test simplification) | done |

## `src/symbolic/core/ast_impl.rs` (5)

| legacy function | new home | status |
|---|---|---|
| `re` | Expr plumbing (term inspection through `Graph`) | dropped |
| `im` | Expr plumbing (term inspection through `Graph`) | dropped |
| `to_f64` | Expr plumbing (term inspection through `Graph`) | dropped |
| `op` | Expr plumbing (term inspection through `Graph`) | dropped |
| `children` | Expr plumbing (term inspection through `Graph`) | dropped |

## `src/symbolic/core/dag_mgr.rs` (3)

| legacy function | new home | status |
|---|---|---|
| `new` | DAG manager; hash-consing is `graph::store` | dropped |
| `get_or_create_normalized` | DAG manager; hash-consing is `graph::store` | dropped |
| `get_or_create` | DAG manager; hash-consing is `graph::store` | dropped |

## `src/symbolic/core/expr_impl.rs` (4)

| legacy function | new home | status |
|---|---|---|
| `pre_order_walk` | Expr plumbing; traversal is `graph` window/extract, normal form is `rules::arith` | dropped |
| `post_order_walk` | Expr plumbing; traversal is `graph` window/extract, normal form is `rules::arith` | dropped |
| `in_order_walk` | Expr plumbing; traversal is `graph` window/extract, normal form is `rules::arith` | dropped |
| `normalize` | Expr plumbing; traversal is `graph` window/extract, normal form is `rules::arith` | dropped |

## `src/symbolic/core/to_expr.rs` (4)

| legacy function | new home | status |
|---|---|---|
| `to_expr` | Expr plumbing | dropped |
| `new` | Expr plumbing | dropped |
| `clone_box_dist` | Expr plumbing | dropped |
| `clone_box_quant` | Expr plumbing | dropped |

## `src/symbolic/cryptography.rs` (15)

| legacy function | new home | status |
|---|---|---|
| `is_infinity` |  | pending |
| `x` |  | pending |
| `y` |  | pending |
| `new` |  | pending |
| `is_on_curve` |  | pending |
| `negate` |  | pending |
| `double` |  | pending |
| `add` |  | pending |
| `scalar_mult` |  | pending |
| `generate_keypair` |  | pending |
| `generate_shared_secret` |  | pending |
| `point_compress` |  | pending |
| `point_decompress` |  | pending |
| `ecdsa_sign` |  | pending |
| `ecdsa_verify` |  | pending |

## `src/symbolic/differential_geometry.rs` (7)

| legacy function | new home | status |
|---|---|---|
| `exterior_derivative` |  | pending |
| `wedge_product` |  | pending |
| `boundary` |  | pending |
| `generalized_stokes_theorem` |  | pending |
| `gauss_theorem` |  | pending |
| `stokes_theorem` |  | pending |
| `greens_theorem` |  | pending |

## `src/symbolic/discrete_groups.rs` (4)

| legacy function | new home | status |
|---|---|---|
| `cyclic_group` |  | pending |
| `dihedral_group` |  | pending |
| `symmetric_group` |  | pending |
| `klein_four_group` |  | pending |

## `src/symbolic/electromagnetism.rs` (8)

| legacy function | new home | status |
|---|---|---|
| `new` |  | pending |
| `lorentz_force` | `sim::classical::lorentz_force` (f64 only; test_lorentz_force) | partial |
| `electric_field_from_potentials` |  | pending |
| `electric_field_from_potential` |  | pending |
| `magnetic_field_from_vector_potential` |  | pending |
| `poynting_vector` |  | pending |
| `energy_density` |  | pending |
| `coulombs_law` | `sim::classical::electric_field_point_charge` / `coulomb_force` (f64 only) | partial |

## `src/symbolic/elementary.rs` (34)

| legacy function | new home | status |
|---|---|---|
| `sin` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `cos` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `tan` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `sinh` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `cosh` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `tanh` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `ln` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `exp` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `sqrt` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `pow` | core `pow` operator (`Term::pow`, api.rs) | done |
| `infinity` | `oo` operator of `rules::calculus` (test limits_at_infinity) | done |
| `negative_infinity` | `oo` operator of `rules::calculus` (test limits_at_infinity) | done |
| `log_base` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `cot` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `sec` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `csc` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `acot` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `asec` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `acsc` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `coth` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `sech` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `csch` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `asinh` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `acosh` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `atanh` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `acoth` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `asech` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `acsch` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `atan2` | `rules::elementary` operator (tests exact_special_values, reciprocal_and_inverse_families, numeric_phase) | done |
| `pi` | `rules::elementary` `pi` and `E` operators (test exact_special_values) | done |
| `e` | `rules::elementary` `pi` and `E` operators (test exact_special_values) | done |
| `expand` | `rules::poly` `expand` (polynomial expansion); sum-angle expansion only in the Explore tier of `rules::elementary` | partial |
| `expand_internal` | `rules::poly` `expand` (polynomial expansion); sum-angle expansion only in the Explore tier of `rules::elementary` | partial |
| `binomial_coefficient` | `rules::combinatorics` `binomial` (test binomials) | done |

## `src/symbolic/error_correction.rs` (13)

| legacy function | new home | status |
|---|---|---|
| `hamming_distance` |  | pending |
| `hamming_weight` |  | pending |
| `hamming_encode` |  | pending |
| `hamming_check` |  | pending |
| `hamming_decode` |  | pending |
| `rs_encode` |  | pending |
| `rs_check` |  | pending |
| `rs_error_count` |  | pending |
| `rs_decode` |  | pending |
| `crc32_compute` |  | pending |
| `crc32_verify` |  | pending |
| `crc32_update` |  | pending |
| `crc32_finalize` |  | pending |

## `src/symbolic/error_correction_helper.rs` (23)

| legacy function | new home | status |
|---|---|---|
| `new` |  | pending |
| `from_bigint` |  | pending |
| `is_zero` |  | pending |
| `is_one` |  | pending |
| `inverse` |  | pending |
| `pow` |  | pending |
| `gf256_exp` |  | pending |
| `gf256_log` |  | pending |
| `gf256_add` |  | pending |
| `gf256_mul` |  | pending |
| `gf256_inv` |  | pending |
| `gf256_div` |  | pending |
| `gf256_pow` |  | pending |
| `poly_eval_gf256` |  | pending |
| `poly_add_gf256` |  | pending |
| `poly_mul_gf256` |  | pending |
| `poly_scale_gf256` |  | pending |
| `poly_derivative_gf256` |  | pending |
| `poly_gcd_gf256` |  | pending |
| `poly_div_gf256` |  | pending |
| `poly_add_gf` |  | pending |
| `poly_mul_gf` |  | pending |
| `poly_div_gf` |  | pending |

## `src/symbolic/finite_field.rs` (10)

| legacy function | new home | status |
|---|---|---|
| `serialize` |  | pending |
| `deserialize` |  | pending |
| `new` |  | pending |
| `inverse` |  | pending |
| `degree` |  | pending |
| `long_division` |  | pending |
| `add` |  | pending |
| `sub` |  | pending |
| `mul` |  | pending |
| `div` |  | pending |

## `src/symbolic/fractal_geometry_and_chaos.rs` (12)

| legacy function | new home | status |
|---|---|---|
| `new` |  | pending |
| `apply` |  | pending |
| `similarity_dimension` |  | pending |
| `new_mandelbrot_family` |  | pending |
| `iterate` |  | pending |
| `orbit` |  | pending |
| `fixed_points` |  | pending |
| `stability_index` |  | pending |
| `find_fixed_points` |  | pending |
| `analyze_stability` |  | pending |
| `lyapunov_exponent` |  | pending |
| `lorenz_system` |  | pending |

## `src/symbolic/functional_analysis.rs` (13)

| legacy function | new home | status |
|---|---|---|
| `new` |  | pending |
| `apply` |  | pending |
| `inner_product` | `kernels::functional_analysis::inner_product` (discrete samples only, no symbolic integrals) | partial |
| `inner_product_internal` | `kernels::functional_analysis::inner_product` (discrete samples only, no symbolic integrals) | partial |
| `norm` | `kernels::functional_analysis::{l1_norm,l2_norm,infinity_norm}` (discrete samples; no general L^p) | partial |
| `norm_internal` | `kernels::functional_analysis::{l1_norm,l2_norm,infinity_norm}` (discrete samples; no general L^p) | partial |
| `banach_norm` | `kernels::functional_analysis::{l1_norm,l2_norm,infinity_norm}` (discrete samples; no general L^p) | partial |
| `banach_norm_internal` | `kernels::functional_analysis::{l1_norm,l2_norm,infinity_norm}` (discrete samples; no general L^p) | partial |
| `are_orthogonal` |  | pending |
| `project` | `kernels::functional_analysis::project` (discrete samples) | partial |
| `project_internal` | `kernels::functional_analysis::project` (discrete samples) | partial |
| `gram_schmidt` | `kernels::functional_analysis::gram_schmidt` (discrete samples) | partial |
| `gram_schmidt_orthonormal` | `kernels::functional_analysis::gram_schmidt_orthonormal` (discrete samples) | partial |

## `src/symbolic/geometric_algebra.rs` (11)

| legacy function | new home | status |
|---|---|---|
| `new` | `Multivector3D` fields (numeric G3 only) | partial |
| `scalar` | `Multivector3D` fields (numeric G3 only) | partial |
| `vector` | `Multivector3D` fields (numeric G3 only) | partial |
| `geometric_product` | `kernels::geometric_algebra::Multivector3D` `Mul` (numeric G3 only) | partial |
| `grade_projection` |  | pending |
| `outer_product` | `Multivector3D::wedge` (numeric G3 only) | partial |
| `inner_product` | `Multivector3D::dot` (numeric G3 only) | partial |
| `reverse` | `Multivector3D::reverse` (numeric G3 only) | partial |
| `magnitude` | `Multivector3D::norm` (numeric G3 only) | partial |
| `dual` |  | pending |
| `normalize` |  | pending |

## `src/symbolic/graph.rs` (15)

| legacy function | new home | status |
|---|---|---|
| `new` |  | pending |
| `nodes` |  | pending |
| `node_count` |  | pending |
| `is_directed` |  | pending |
| `add_node` |  | pending |
| `add_edge` |  | pending |
| `get_node_id` |  | pending |
| `neighbors` |  | pending |
| `out_degree` |  | pending |
| `in_degree` |  | pending |
| `get_edges` |  | pending |
| `add_hyperedge` |  | pending |
| `to_adjacency_matrix` |  | pending |
| `to_incidence_matrix` |  | pending |
| `to_laplacian_matrix` |  | pending |

## `src/symbolic/graph_algorithms.rs` (26)

| legacy function | new home | status |
|---|---|---|
| `dfs` |  | pending |
| `bfs` |  | pending |
| `connected_components` |  | pending |
| `is_connected` |  | pending |
| `strongly_connected_components` |  | pending |
| `has_cycle` |  | pending |
| `find_bridges_and_articulation_points` |  | pending |
| `kruskal_mst` |  | pending |
| `edmonds_karp_max_flow` |  | pending |
| `dinic_max_flow` |  | pending |
| `bellman_ford` |  | pending |
| `min_cost_max_flow` |  | pending |
| `is_bipartite` |  | pending |
| `bipartite_maximum_matching` |  | pending |
| `prim_mst` |  | pending |
| `topological_sort_kahn` |  | pending |
| `topological_sort_dfs` |  | pending |
| `topological_sort` |  | pending |
| `bipartite_minimum_vertex_cover` |  | pending |
| `hopcroft_karp_bipartite_matching` |  | pending |
| `blossom_algorithm` |  | pending |
| `shortest_path_unweighted` |  | pending |
| `dijkstra` |  | pending |
| `floyd_warshall` |  | pending |
| `spectral_analysis` |  | pending |
| `algebraic_connectivity` |  | pending |

## `src/symbolic/graph_isomorphism_and_coloring.rs` (3)

| legacy function | new home | status |
|---|---|---|
| `are_isomorphic_heuristic` |  | pending |
| `greedy_coloring` |  | pending |
| `chromatic_number_exact` |  | pending |

## `src/symbolic/graph_operations.rs` (8)

| legacy function | new home | status |
|---|---|---|
| `induced_subgraph` |  | pending |
| `union` |  | pending |
| `intersection` |  | pending |
| `cartesian_product` |  | pending |
| `tensor_product` |  | pending |
| `complement` |  | pending |
| `disjoint_union` |  | pending |
| `join` |  | pending |

## `src/symbolic/grobner.rs` (4)

| legacy function | new home | status |
|---|---|---|
| `poly_division_multivariate` | `poly::groebner::reduce` (exercised by groebner tests) | done |
| `subtract_poly` | `poly::repr::Poly::sub` (test arithmetic_laws) | done |
| `buchberger` | `poly::groebner::groebner` (tests circle_and_line_eliminate_to_a_univariate_polynomial, linear_system) | done |
| `reduced_basis` | `poly::groebner::groebner` returns the reduced basis (test the_reduced_basis_does_not_depend_on_the_generators) | done |

## `src/symbolic/group_theory.rs` (9)

| legacy function | new home | status |
|---|---|---|
| `new` |  | pending |
| `multiply` |  | pending |
| `inverse` |  | pending |
| `is_abelian` |  | pending |
| `element_order` |  | pending |
| `conjugacy_classes` |  | pending |
| `center` |  | pending |
| `is_valid` |  | pending |
| `character` |  | pending |

## `src/symbolic/handles.rs` (8)

| legacy function | new home | status |
|---|---|---|
| `insert` | FFI handle table | dropped |
| `get` | FFI handle table | dropped |
| `clone_expr` | FFI handle table | dropped |
| `free` | FFI handle table | dropped |
| `exists` | FFI handle table | dropped |
| `count` | FFI handle table | dropped |
| `clear` | FFI handle table | dropped |
| `get_all_handles` | FFI handle table | dropped |

## `src/symbolic/integral_equations.rs` (7)

| legacy function | new home | status |
|---|---|---|
| `new` |  | pending |
| `solve_neumann_series` |  | pending |
| `solve_separable_kernel` |  | pending |
| `solve_successive_approximations` |  | pending |
| `solve_by_differentiation` |  | pending |
| `solve_airfoil_equation` |  | pending |
| `solve_airfoil_equation_internal` |  | pending |

## `src/symbolic/integration.rs` (8)

| legacy function | new home | status |
|---|---|---|
| `integrate_rational_function` | `rules::calculus` `integral` partial-fraction stage, exact over Q (test rational_functions) | done |
| `risch_norman_integrate` | `rules::calculus` `integral` staged heuristics; no Risch-Norman undetermined-coefficients ansatz | partial |
| `integrate_poly_exp` | `integral` by_parts stage handles polynomial*exp(ax) (test by_parts); no general exponential extension | partial |
| `poly_from_coeffs` | `poly::repr::Poly::from_univariate` | done |
| `partial_fraction_integrate` | `rules::calculus` `integral` partial-fraction stage, exact over Q (test rational_functions) | done |
| `hermite_integrate_rational` | `rules::calculus` `integral` partial-fraction stage, exact over Q (test rational_functions) | done |
| `integrate_rational_function_expr` | `rules::calculus` `integral` partial-fraction stage, exact over Q (test rational_functions) | done |
| `poly_derivative_symbolic` | `poly::repr::Poly::derivative` / `diff` (test derivatives_symbolic) | done |

## `src/symbolic/lie_groups_and_algebras.rs` (11)

| legacy function | new home | status |
|---|---|---|
| `lie_bracket` |  | pending |
| `lie_bracket_internal` |  | pending |
| `exponential_map` |  | pending |
| `adjoint_representation_group` |  | pending |
| `adjoint_representation_algebra` |  | pending |
| `commutator_table` |  | pending |
| `check_jacobi_identity` |  | pending |
| `so3_generators` |  | pending |
| `so3` |  | pending |
| `su2_generators` |  | pending |
| `su2` |  | pending |

## `src/symbolic/logic.rs` (6)

| legacy function | new home | status |
|---|---|---|
| `simplify_logic` | `rules::logic` `simplify_logic` (test minimal_sum_of_products) | done |
| `to_cnf` | `rules::logic` `cnf` (test conjunctive_normal_form) | done |
| `to_cnf_internal` | `rules::logic` `cnf` (test conjunctive_normal_form) | done |
| `to_dnf` | `rules::logic` `dnf` (test disjunctive_normal_form) | done |
| `to_dnf_internal` | `rules::logic` `dnf` (test disjunctive_normal_form) | done |
| `is_satisfiable` | `rules::logic` `satisfiable` (test satisfiability) | done |

## `src/symbolic/matrix.rs` (31)

| legacy function | new home | status |
|---|---|---|
| `get_matrix_dims` | `rules::linalg` `dims` | done |
| `create_empty_matrix` | `rules::linalg` `zeros` | done |
| `identity_matrix` | `rules::linalg` `identity` | done |
| `add_matrices` | `rules::linalg` `madd` | done |
| `sub_matrices` | `rules::linalg` `madd(A, smul(-1, B))` | done |
| `mul_matrices` | `rules::linalg` `matmul` | done |
| `mul_matrices_internal` | `rules::linalg` `matmul` | done |
| `scalar_mul_matrix` | `rules::linalg` `smul` | done |
| `transpose_matrix` | `rules::linalg` `transpose` | done |
| `transpose_internal` | `rules::linalg` `transpose` | done |
| `determinant` | `rules::linalg` `det` (exact rref / memoised Laplace) | done |
| `determinant_internal` | `rules::linalg` `det` | done |
| `inverse_matrix` | `rules::linalg` `inverse` (exact Gauss–Jordan / adjugate) | done |
| `inverse_internal` | `rules::linalg` `inverse` | done |
| `solve_linear_system` | `rules::linalg` `linsolve` (parametric when underdetermined) | done |
| `solve_linear_system_internal` | `rules::linalg` `linsolve` | done |
| `trace` | `rules::linalg` `trace` | done |
| `trace_internal` | `rules::linalg` `trace` | done |
| `characteristic_polynomial` | `rules::linalg` `charpoly` (Faddeev–LeVerrier) | done |
| `characteristic_polynomial_internal` | `rules::linalg` `charpoly` | done |
| `lu_decomposition` | `rules::linalg` `lu` → list(P, L, U) | done |
| `qr_decomposition` | `rules::linalg` `qr` (Gram–Schmidt) | done |
| `rref` | `rules::linalg` `rref` | done |
| `rref_internal` | `rules::linalg` `rref` | done |
| `null_space` | `rules::linalg` `nullspace` | done |
| `null_space_internal` | `rules::linalg` `nullspace` | done |
| `eigen_decomposition` | `rules::linalg` `eigenvals`, `eigenvects` (closed-form roots only) | done |
| `svd_decomposition` | `rules::linalg` `svd` (numeric, faer) | partial |
| `rank` | `rules::linalg` `rank` | done |
| `gaussian_elimination` | `rules::linalg` `rref` | done |
| `is_zero_matrix` | `rules::linalg` `rank(A) = 0` | done |

## `src/symbolic/multi_valued.rs` (12)

| legacy function | new home | status |
|---|---|---|
| `arg` | `rules::complex` `arg` | done |
| `abs` | `rules::complex` `abs` | done |
| `general_log` | `rules::complex::branches` `log_branch` | done |
| `general_sqrt` | `rules::complex::branches` `sqrt_branch` | done |
| `general_power` | `rules::complex::branches` `power_branch` | done |
| `general_nth_root` | `rules::complex::branches` `root_branch` | done |
| `general_arcsin` | `rules::complex::branches` `asin_branch` | done |
| `general_arccos` | `rules::complex::branches` `acos_branch` | done |
| `general_arctan` | `rules::complex::branches` `atan_branch` | done |
| `general_arcsinh` | `rules::complex::branches` `asinh_branch` | done |
| `general_arccosh` | `rules::complex::branches` `acosh_branch` | done |
| `general_arctanh` | `rules::complex::branches` `atanh_branch` | done |

## `src/symbolic/number_theory.rs` (15)

| legacy function | new home | status |
|---|---|---|
| `expr_to_sparse_poly` | poly::repr::from_term | dropped: helper for the legacy Expr; polynomials are handled by poly::repr |
| `solve_diophantine` | `diophantine(eq, list(vars))`: linear (parametric t), Pell +-1, Pythagorean | done |
| `solve_diophantine_internal` | `diophantine` | done |
| `solve_pell_from_poly` | `diophantine` / `pell(n)` | done |
| `is_neg_one` |  | dropped: Expr predicate |
| `is_two` |  | dropped: Expr predicate |
| `extended_gcd_inner` | `egcd` | done |
| `chinese_remainder` | `crt` | done |
| `chinese_remainder_internal` | `crt` | done |
| `is_prime` | `isprime` | done |
| `is_prime_internal` | `isprime` | done |
| `sqrt_continued_fraction` | `cfrac(sqrt(n))` | done |
| `get_convergent` | `convergents(x, k)` | done |
| `extended_gcd` | `egcd` | done |
| `extended_gcd_internal` | `egcd` | done |

## `src/symbolic/numeric.rs` (2)

| legacy function | new home | status |
|---|---|---|
| `evaluate_numerical` | `Term::eval` / `Answer::as_f64` (api.rs; test derivatives_agree_across_phases) | done |
| `evaluate_complex` |  | pending |

## `src/symbolic/ode.rs` (22)

| legacy function | new home | status |
|---|---|---|
| `solve_ode` | `rules::ode` `dsolve` (tests first_order_linear_and_separable ... cauchy_euler) | done |
| `solve_ode_internal` | `rules::ode` `dsolve` (tests first_order_linear_and_separable ... cauchy_euler) | done |
| `solve_ode_system` | `rules::ode` `odeint` (numeric systems only; no symbolic system solver) | partial |
| `solve_ode_system_internal` | `rules::ode` `odeint` (numeric systems only; no symbolic system solver) | partial |
| `solve_separable_ode` | `rules::ode` `dsolve` (test first_order_linear_and_separable) | done |
| `solve_separable_ode_internal` | `rules::ode` `dsolve` (test first_order_linear_and_separable) | done |
| `solve_first_order_linear_ode` | `rules::ode` `dsolve` (test first_order_linear_and_separable) | done |
| `solve_first_order_linear_ode_internal` | `rules::ode` `dsolve` (test first_order_linear_and_separable) | done |
| `solve_bernoulli_ode` | `rules::ode` `dsolve` (test bernoulli_riccati_homogeneous_exact) | done |
| `solve_bernoulli_ode_internal` | `rules::ode` `dsolve` (test bernoulli_riccati_homogeneous_exact) | done |
| `solve_riccati_ode` | `rules::ode` `dsolve` (test bernoulli_riccati_homogeneous_exact) | done |
| `solve_riccati_ode_internal` | `rules::ode` `dsolve` (test bernoulli_riccati_homogeneous_exact) | done |
| `solve_cauchy_euler_ode` | `rules::ode` `dsolve` (test cauchy_euler) | done |
| `solve_cauchy_euler_ode_internal` | `rules::ode` `dsolve` (test cauchy_euler) | done |
| `solve_by_reduction_of_order` |  | pending |
| `solve_by_reduction_of_order_internal` |  | pending |
| `solve_exact_ode` | `rules::ode` `dsolve` (test bernoulli_riccati_homogeneous_exact) | done |
| `solve_exact_ode_internal` | `rules::ode` `dsolve` (test bernoulli_riccati_homogeneous_exact) | done |
| `solve_ode_by_series` |  | pending |
| `solve_ode_by_series_internal` |  | pending |
| `solve_ode_by_fourier` |  | pending |
| `solve_ode_by_fourier_internal` |  | pending |

## `src/symbolic/optimize.rs` (4)

| legacy function | new home | status |
|---|---|---|
| `find_extrema` |  | pending |
| `hessian_matrix` | `rules::linalg` `hessian` (test vector_calculus) | done |
| `hessian_matrix_internal` | `rules::linalg` `hessian` (test vector_calculus) | done |
| `find_constrained_extrema` |  | pending |

## `src/symbolic/pde.rs` (37)

| legacy function | new home | status |
|---|---|---|
| `solve_pde` |  | pending |
| `solve_pde_internal` |  | pending |
| `solve_pde_by_separation_of_variables` |  | pending |
| `solve_pde_by_separation_of_variables_internal` |  | pending |
| `classify_pde_heuristic` |  | pending |
| `solve_pde_by_characteristics` |  | pending |
| `solve_pde_by_characteristics_internal` |  | pending |
| `solve_pde_by_greens_function` |  | pending |
| `solve_pde_by_greens_function_internal` |  | pending |
| `solve_second_order_pde` |  | pending |
| `solve_second_order_pde_internal` |  | pending |
| `solve_wave_equation_1d_dalembert` |  | pending |
| `solve_wave_equation_1d_dalembert_internal` |  | pending |
| `solve_heat_equation_1d` |  | pending |
| `solve_heat_equation_1d_internal` |  | pending |
| `solve_laplace_equation_2d` |  | pending |
| `solve_laplace_equation_2d_internal` |  | pending |
| `solve_wave_equation_3d` |  | pending |
| `solve_wave_equation_3d_internal` |  | pending |
| `solve_heat_equation_3d` |  | pending |
| `solve_heat_equation_3d_internal` |  | pending |
| `solve_laplace_equation_3d` |  | pending |
| `solve_laplace_equation_3d_internal` |  | pending |
| `solve_poisson_equation_2d` |  | pending |
| `solve_poisson_equation_2d_internal` |  | pending |
| `solve_poisson_equation_3d` |  | pending |
| `solve_poisson_equation_3d_internal` |  | pending |
| `solve_helmholtz_equation` |  | pending |
| `solve_helmholtz_equation_internal` |  | pending |
| `solve_schrodinger_equation` |  | pending |
| `solve_schrodinger_equation_internal` |  | pending |
| `solve_klein_gordon_equation` |  | pending |
| `solve_klein_gordon_equation_internal` |  | pending |
| `solve_burgers_equation` |  | pending |
| `solve_burgers_equation_internal` |  | pending |
| `solve_with_fourier_transform` |  | pending |
| `solve_with_fourier_transform_internal` |  | pending |

## `src/symbolic/poly_factorization.rs` (11)

| legacy function | new home | status |
|---|---|---|
| `factor_gf` | private GF(p) helper inside `rules::poly::univariate`, reachable only through `factor` over Q; no public GF(p) API | partial |
| `poly_derivative_gf` | private GF(p) helper inside `rules::poly::univariate`, reachable only through `factor` over Q; no public GF(p) API | partial |
| `square_free_factorization_gf` |  | pending |
| `berlekamp_factorization` |  | pending |
| `berlekamp_zassenhaus` | `rules::poly::univariate::factor` (Berlekamp-Zassenhaus over Q; tests known_factorisations, factorisation_reproduces_the_input) | done |
| `cantor_zassenhaus` | private GF(p) helper inside `rules::poly::univariate`, reachable only through `factor` over Q; no public GF(p) API | partial |
| `distinct_degree_factorization` | private GF(p) helper inside `rules::poly::univariate`, reachable only through `factor` over Q; no public GF(p) API | partial |
| `poly_gcd_gf` | private GF(p) helper inside `rules::poly::univariate`, reachable only through `factor` over Q; no public GF(p) API | partial |
| `poly_pow_mod` | private GF(p) helper inside `rules::poly::univariate`, reachable only through `factor` over Q; no public GF(p) API | partial |
| `poly_mul_scalar` | private GF(p) helper inside `rules::poly::univariate`, reachable only through `factor` over Q; no public GF(p) API | partial |
| `poly_extended_gcd` | private GF(p) helper inside `rules::poly::univariate`, reachable only through `factor` over Q; no public GF(p) API | partial |

## `src/symbolic/polynomial.rs` (24)

| legacy function | new home | status |
|---|---|---|
| `add_poly` | `poly::repr::Poly::{add,mul,scale}` (test arithmetic_laws) | done |
| `mul_poly` | `poly::repr::Poly::{add,mul,scale}` (test arithmetic_laws) | done |
| `differentiate_poly` | `poly::repr::Poly::derivative` | done |
| `contains_var` | Expr helper; `free_of` guard / `poly::repr::from_term` / Poly keeps no zero terms | dropped |
| `is_polynomial` | `rules::poly` `degree` (test degree_and_coefficients; stays unreduced for non-polynomials) | done |
| `polynomial_degree` | `rules::poly` `degree` (test degree_and_coefficients; stays unreduced for non-polynomials) | done |
| `leading_coefficient` | `rules::poly` `coeff` (test degree_and_coefficients) | done |
| `polynomial_long_division` | `rules::poly` `quo`/`rem` and `univariate::divrem` (test division_and_gcd) | done |
| `polynomial_long_division_internal` | `rules::poly` `quo`/`rem` and `univariate::divrem` (test division_and_gcd) | done |
| `to_polynomial_coeffs_vec` | `poly::repr::Poly::univariate_in` (test univariate_views) | done |
| `from_coeffs_to_expr` | `poly::repr::Poly::from_univariate` / `to_term` (test expansion) | done |
| `polynomial_long_division_coeffs` | `rules::poly` `quo`/`rem` and `univariate::divrem` (test division_and_gcd) | done |
| `expr_to_sparse_poly` | Expr helper; `free_of` guard / `poly::repr::from_term` / Poly keeps no zero terms | dropped |
| `eval` | `Term::eval` / `univariate::eval` | done |
| `poly_mul_scalar_expr` | `poly::repr::Poly::{add,mul,scale}` (test arithmetic_laws) | done |
| `gcd` | `rules::poly` `pgcd` (test division_and_gcd) | done |
| `degree` | `rules::poly` `degree` (test degree_and_coefficients; stays unreduced for non-polynomials) | done |
| `leading_term` | `poly::repr::Poly::terms` ordered by monomial (no dedicated accessor) | partial |
| `long_division` | `rules::poly` `quo`/`rem` and `univariate::divrem` (test division_and_gcd) | done |
| `get_coeffs_as_vec` | `poly::repr::Poly::univariate_in` (test univariate_views) | done |
| `get_coeff_for_power` | `rules::poly` `coeff` (test degree_and_coefficients) | done |
| `prune_zeros` | Expr helper; `free_of` guard / `poly::repr::from_term` / Poly keeps no zero terms | dropped |
| `poly_from_coeffs` | `poly::repr::Poly::from_univariate` / `to_term` (test expansion) | done |
| `sparse_poly_to_expr` | `poly::repr::to_term` (test expansion) | done |

## `src/symbolic/proof.rs` (7)

| legacy function | new home | status |
|---|---|---|
| `verify_equation_solution` |  | pending |
| `verify_indefinite_integral` |  | pending |
| `verify_definite_integral` |  | pending |
| `verify_ode_solution` |  | pending |
| `verify_matrix_inverse` |  | pending |
| `verify_derivative` |  | pending |
| `verify_limit` |  | pending |

## `src/symbolic/quantum_field_theory.rs` (8)

| legacy function | new home | status |
|---|---|---|
| `dirac_adjoint` |  | pending |
| `feynman_slash` |  | pending |
| `scalar_field_lagrangian` |  | pending |
| `qed_lagrangian` |  | pending |
| `qcd_lagrangian` |  | pending |
| `propagator` |  | pending |
| `scattering_cross_section` |  | pending |
| `feynman_propagator_position_space` |  | pending |

## `src/symbolic/quantum_mechanics.rs` (21)

| legacy function | new home | status |
|---|---|---|
| `bra_ket` |  | pending |
| `bra_ket_internal` |  | pending |
| `new` |  | pending |
| `apply` |  | pending |
| `commutator` |  | pending |
| `commutator_internal` |  | pending |
| `expectation_value` |  | pending |
| `expectation_value_internal` |  | pending |
| `uncertainty` |  | pending |
| `probability_density` |  | pending |
| `hamiltonian_free_particle` |  | pending |
| `hamiltonian_harmonic_oscillator` |  | pending |
| `angular_momentum_z` |  | pending |
| `pauli_matrices` |  | pending |
| `spin_operator` |  | pending |
| `solve_time_independent_schrodinger` |  | pending |
| `time_dependent_schrodinger_equation` |  | pending |
| `dirac_equation` |  | pending |
| `klein_gordon_equation` |  | pending |
| `first_order_energy_correction` |  | pending |
| `scattering_amplitude` |  | pending |

## `src/symbolic/radicals.rs` (2)

| legacy function | new home | status |
|---|---|---|
| `simplify_radicals` |  | pending |
| `denest_sqrt` |  | pending |

## `src/symbolic/real_roots.rs` (4)

| legacy function | new home | status |
|---|---|---|
| `sturm_sequence` | `kernels::real_roots::sturm_sequence` (tests/kernels/real_roots.rs) | done |
| `count_real_roots_in_interval` | `kernels::real_roots::isolate_real_roots` returns intervals (count = length); no interval-count function | partial |
| `isolate_real_roots` | `kernels::real_roots::isolate_real_roots` (tests/kernels/real_roots.rs) | done |
| `eval_expr` | Expr helper; `Term::eval` | dropped |

## `src/symbolic/relativity.rs` (13)

| legacy function | new home | status |
|---|---|---|
| `lorentz_factor` | `sim::classical::lorentz_factor` (f64 only; test_lorentz_factor_high_speed) | partial |
| `lorentz_transformation_x` |  | pending |
| `velocity_addition` | `sim::classical::relativistic_velocity_addition` (f64 only) | partial |
| `mass_energy_equivalence` | `sim::classical::mass_energy` (f64 only) | partial |
| `relativistic_momentum` | `sim::classical::relativistic_momentum` (f64 only) | partial |
| `doppler_effect` |  | pending |
| `schwarzschild_radius` |  | pending |
| `gravitational_time_dilation` |  | pending |
| `einstein_tensor` |  | pending |
| `geodesic_acceleration` |  | pending |
| `lorentz_transformation` | legacy alias of lorentz_transformation_x | dropped |
| `einstein_field_equations` | legacy placeholder, superseded by einstein_tensor | dropped |
| `geodesic_equation` | legacy placeholder, superseded by geodesic_acceleration | dropped |

## `src/symbolic/rewriting.rs` (2)

| legacy function | new home | status |
|---|---|---|
| `apply_rules_to_normal_form` | `Session::with_rules` + `Engine` saturation (rules fixed at session build, not per call) | partial |
| `knuth_bendix` |  | pending |

## `src/symbolic/series.rs` (14)

| legacy function | new home | status |
|---|---|---|
| `taylor_series` | `rules::calculus` `taylor` (test taylor_and_laurent) | done |
| `taylor_series_internal` | `rules::calculus` `taylor` (test taylor_and_laurent) | done |
| `calculate_taylor_coefficients` | `rules::calculus` `taylor` (test taylor_and_laurent) | done |
| `laurent_series` | `rules::calculus` `laurent` (test taylor_and_laurent) | done |
| `laurent_series_internal` | `rules::calculus` `laurent` (test taylor_and_laurent) | done |
| `fourier_series` | `rules::calculus` `fourier_series` (test fourier_series_of_simple_functions) | done |
| `fourier_series_internal` | `rules::calculus` `fourier_series` (test fourier_series_of_simple_functions) | done |
| `summation` | `rules::calculus` `sum`/`product` (tests sums_and_products, numeric_sums) | done |
| `summation_internal` | `rules::calculus` `sum`/`product` (tests sums_and_products, numeric_sums) | done |
| `product` | `rules::calculus` `sum`/`product` (tests sums_and_products, numeric_sums) | done |
| `product_internal` | `rules::calculus` `sum`/`product` (tests sums_and_products, numeric_sums) | done |
| `analyze_convergence` | `rules::calculus` `converges` (ratio test only) | partial |
| `asymptotic_expansion` |  | pending |
| `analytic_continuation` |  | pending |

## `src/symbolic/simplify.rs` (9)

| legacy function | new home | status |
|---|---|---|
| `simplify` | `Session::simplify` (test simplification) | done |
| `heuristic_simplify` | `Session::simplify` (test simplification) | done |
| `is_zero` | Expr predicate / rule naming | dropped |
| `is_infinite` | Expr predicate / rule naming | dropped |
| `is_one` | Expr predicate / rule naming | dropped |
| `as_f64` | `Term::as_f64` (api.rs) | done |
| `is_numeric` | Expr predicate / rule naming | dropped |
| `get_name` | Expr predicate / rule naming | dropped |
| `collect_and_order_terms` | `rules::arith` Collect (test like_terms_are_collected) | done |

## `src/symbolic/simplify_dag.rs` (3)

| legacy function | new home | status |
|---|---|---|
| `simplify` | `Session::simplify` (test simplification) | done |
| `pattern_match` | `graph::pattern` e-matching (graph/pattern.rs tests) | done |
| `substitute_patterns` | `graph::pattern` instantiation / `graph::subst` (graph/pattern.rs, subst.rs tests) | done |

## `src/symbolic/solid_state_physics.rs` (13)

| legacy function | new home | status |
|---|---|---|
| `new` |  | pending |
| `volume` |  | pending |
| `reciprocal_lattice_vectors` |  | pending |
| `bloch_theorem` |  | pending |
| `energy_band` |  | pending |
| `density_of_states_3d` |  | pending |
| `fermi_energy_3d` |  | pending |
| `drude_conductivity` |  | pending |
| `hall_coefficient` |  | pending |
| `debye_frequency` |  | pending |
| `einstein_heat_capacity` |  | pending |
| `plasma_frequency` |  | pending |
| `london_penetration_depth` |  | pending |

## `src/symbolic/solve.rs` (9)

| legacy function | new home | status |
|---|---|---|
| `solve` | `rules::solve` `solve` (tests linear_and_quadratic, higher_degree_by_factorisation, transcendental_by_inversion) | done |
| `solve_internal` | `rules::solve` `solve` (tests linear_and_quadratic, higher_degree_by_factorisation, transcendental_by_inversion) | done |
| `solve_system` | `rules::solve` `solve` with lists (tests linear_systems, polynomial_systems) | done |
| `solve_system_internal` | `rules::solve` `solve` with lists (tests linear_systems, polynomial_systems) | done |
| `solve_system_parcial` | `rules::solve` `solve` with lists (tests linear_systems, polynomial_systems) | done |
| `solve_linear_system_mat` | `rules::linalg` `linsolve` (test linear_systems) | done |
| `solve_linear_system` | `rules::linalg` `linsolve` (test linear_systems) | done |
| `solve_linear_system_gauss` | `rules::linalg` `linsolve` (test linear_systems) | done |
| `extract_polynomial_coeffs` | `rules::poly` `coeff` / `Poly::coefficients_in` (test degree_and_coefficients) | done |

## `src/symbolic/special.rs` (28)

| legacy function | new home | status |
|---|---|---|
| `gamma_numerical` | `kernels::special::gamma_numerical (test_gamma)` | done |
| `ln_gamma_numerical` | `kernels::special::ln_gamma_numerical (test_ln_gamma)` | done |
| `digamma_numerical` | `kernels::special::digamma_numerical (test_digamma)` | done |
| `beta_numerical` | `kernels::special::beta_numerical (test_beta)` | done |
| `ln_beta_numerical` | `kernels::special::ln_beta_numerical` (test in tests/kernels/special.rs) | done |
| `regularized_incomplete_beta` | `kernels::special::regularized_beta (test_regularized_beta)` | done |
| `erf_numerical` | `kernels::special::erf_numerical (test_erf)` | done |
| `erfc_numerical` | `kernels::special::erfc_numerical (test_erfc)` | done |
| `inverse_erf` | `kernels::special::inverse_erf_numerical (test_inverse_erf)` | done |
| `inverse_erfc` |  | pending |
| `factorial` | `kernels::special::factorial (test_factorial)` | done |
| `double_factorial` | `kernels::special::double_factorial (test_double_factorial)` | done |
| `binomial` | `kernels::special::binomial (test_binomial)` | done |
| `rising_factorial` | `rules::combinatorics` `rising`/`falling` (test falling_and_rising) | done |
| `falling_factorial` | `rules::combinatorics` `rising`/`falling` (test falling_and_rising) | done |
| `bessel_j0` | `kernels::special::bessel_j0 (test_bessel_j0)` | done |
| `bessel_j1` | `kernels::special::bessel_j1 (test_bessel_j1)` | done |
| `bessel_y0` | `kernels::special::bessel_y0 (test_bessel_y0)` | done |
| `bessel_y1` | `kernels::special::bessel_y1 (bessel_y1_small_argument_reference_values)` | done |
| `bessel_i0` | `kernels::special::bessel_i0 (test_bessel_i0)` | done |
| `bessel_i1` | `kernels::special::bessel_i1 (modified_bessel_recurrence_and_edges)` | done |
| `bessel_k0` |  | pending |
| `bessel_k1` |  | pending |
| `sinc` | `kernels::special::sinc (test_sinc)` | done |
| `zeta` | `kernels::special::riemann_zeta (test_riemann_zeta)` | done |
| `ln_factorial` |  | pending |
| `regularized_gamma_p` | `kernels::special::regularized_lower_gamma (incomplete_gamma_closed_forms)` | done |
| `regularized_gamma_q` | `kernels::special::regularized_upper_gamma (incomplete_gamma_closed_forms)` | done |

## `src/symbolic/special_functions.rs` (26)

| legacy function | new home | status |
|---|---|---|
| `gamma` | `rules::special` `gamma` (test gamma_family_values) | done |
| `ln_gamma` | `rules::special` `lgamma` (test gamma_family_values) | done |
| `beta` | `rules::special` `beta` (test gamma_family_values) | done |
| `digamma` | `rules::special` `digamma` (test gamma_family_values) | done |
| `polygamma` | `kernels::special::polygamma_numerical` only; no `polygamma` operator | partial |
| `erf` | `rules::special` `erf`/`erfc` (test error_function_values) | done |
| `erfc` | `rules::special` `erf`/`erfc` (test error_function_values) | done |
| `erfi` |  | pending |
| `zeta` | `rules::special` `zeta` (test zeta_values) | done |
| `bessel_j` | `rules::special` `besselj`/`bessely`/`besseli` (test bessel_values; integer orders 0 and 1 only) | partial |
| `bessel_y` | `rules::special` `besselj`/`bessely`/`besseli` (test bessel_values; integer orders 0 and 1 only) | partial |
| `bessel_i` | `rules::special` `besselj`/`bessely`/`besseli` (test bessel_values; integer orders 0 and 1 only) | partial |
| `bessel_k` |  | pending |
| `legendre_p` | `rules::special` `legendre`/`laguerre`/`hermite`/`chebyshevt`/`chebyshevu` (tests polynomials_of_any_degree_match_the_recurrence) | done |
| `laguerre_l` | `rules::special` `legendre`/`laguerre`/`hermite`/`chebyshevt`/`chebyshevu` (tests polynomials_of_any_degree_match_the_recurrence) | done |
| `generalized_laguerre` |  | pending |
| `hermite_h` | `rules::special` `legendre`/`laguerre`/`hermite`/`chebyshevt`/`chebyshevu` (tests polynomials_of_any_degree_match_the_recurrence) | done |
| `chebyshev_t` | `rules::special` `legendre`/`laguerre`/`hermite`/`chebyshevt`/`chebyshevu` (tests polynomials_of_any_degree_match_the_recurrence) | done |
| `chebyshev_u` | `rules::special` `legendre`/`laguerre`/`hermite`/`chebyshevt`/`chebyshevu` (tests polynomials_of_any_degree_match_the_recurrence) | done |
| `bessel_differential_equation` |  | pending |
| `legendre_differential_equation` |  | pending |
| `legendre_rodrigues_formula` |  | pending |
| `laguerre_differential_equation` |  | pending |
| `hermite_differential_equation` |  | pending |
| `hermite_rodrigues_formula` |  | pending |
| `chebyshev_differential_equation` |  | pending |

## `src/symbolic/stats.rs` (10)

| legacy function | new home | status |
|---|---|---|
| `mean` | `mean` | done |
| `mean_internal` | `mean` | done |
| `variance` | `variance` | done |
| `variance_internal` | `variance` | done |
| `std_dev` | `std` | done |
| `std_dev_internal` | `std` | done |
| `covariance` | `covariance` | done |
| `covariance_internal` | `covariance` | done |
| `correlation` | `correlation` | done |
| `correlation_internal` | `correlation` | done |

## `src/symbolic/stats_inference.rs` (3)

| legacy function | new home | status |
|---|---|---|
| `one_sample_t_test_symbolic` | `t_test` | done |
| `two_sample_t_test_symbolic` | `welch_test` | done |
| `z_test_symbolic` | `z_test` | done |

## `src/symbolic/stats_information_theory.rs` (9)

| legacy function | new home | status |
|---|---|---|
| `shannon_entropy` | `entropy` (nats), `entropy_bits` | done |
| `shannon_entropy_internal` | `entropy_bits` | done |
| `kl_divergence` | `kl_divergence` (nats) | done |
| `cross_entropy` | `cross_entropy` (nats) | done |
| `joint_entropy` | `joint_entropy` | done |
| `conditional_entropy` | `conditional_entropy` | done |
| `mutual_information` | `mutual_information` | done |
| `gini_impurity` | `gini` | done |
| `gini_impurity_internal` | `gini` | done |

## `src/symbolic/stats_regression.rs` (3)

| legacy function | new home | status |
|---|---|---|
| `simple_linear_regression_symbolic` | `linear_regression` | done |
| `nonlinear_regression_symbolic` | `nonlinear_regression` (builds a `solve` request; needs the solve rules) | done |
| `polynomial_regression_symbolic` | `polynomial_regression` | partial: exact/float data only, symbolic data not supported |

## `src/symbolic/tensor.rs` (15)

| legacy function | new home | status |
|---|---|---|
| `new` | `kernels::tensor::TensorData` / `ndarray` arrays (f64 only; no symbolic components) | partial |
| `rank` | `kernels::tensor::TensorData` / `ndarray` arrays (f64 only; no symbolic components) | partial |
| `get` | `kernels::tensor::TensorData` / `ndarray` arrays (f64 only; no symbolic components) | partial |
| `add` | `kernels::tensor::TensorData` / `ndarray` arrays (f64 only; no symbolic components) | partial |
| `sub` | `kernels::tensor::TensorData` / `ndarray` arrays (f64 only; no symbolic components) | partial |
| `scalar_mul` | `kernels::tensor::TensorData` / `ndarray` arrays (f64 only; no symbolic components) | partial |
| `outer_product` | `kernels::tensor::outer_product` (f64 ndarray only) | partial |
| `contract` | `kernels::tensor::contract` (f64 ndarray only) | partial |
| `to_matrix_expr` | `kernels::tensor::TensorData` / `ndarray` arrays (f64 only; no symbolic components) | partial |
| `raise_index` |  | pending |
| `lower_index` |  | pending |
| `christoffel_symbols_first_kind` |  | pending |
| `christoffel_symbols_second_kind` |  | pending |
| `riemann_curvature_tensor` |  | pending |
| `covariant_derivative_vector` |  | pending |

## `src/symbolic/thermodynamics.rs` (13)

| legacy function | new home | status |
|---|---|---|
| `first_law_thermodynamics` |  | pending |
| `ideal_gas_law` | `sim::classical::{ideal_gas_pressure,ideal_gas_volume,ideal_gas_temperature}` (f64 only) | partial |
| `enthalpy` |  | pending |
| `helmholtz_free_energy` |  | pending |
| `gibbs_free_energy` |  | pending |
| `boltzmann_entropy` |  | pending |
| `carnot_efficiency` |  | pending |
| `boltzmann_distribution` |  | pending |
| `partition_function` |  | pending |
| `fermi_dirac_distribution` |  | pending |
| `bose_einstein_distribution` |  | pending |
| `work_isothermal_expansion` |  | pending |
| `verify_maxwell_relation_helmholtz` |  | pending |

## `src/symbolic/topology.rs` (19)

| legacy function | new home | status |
|---|---|---|
| `new` |  | pending |
| `dimension` |  | pending |
| `boundary` |  | pending |
| `symbolic_boundary` |  | pending |
| `add_term` |  | pending |
| `add_simplex` |  | pending |
| `get_simplices_by_dim` |  | pending |
| `get_boundary_matrix` |  | pending |
| `get_symbolic_boundary_matrix` |  | pending |
| `apply_boundary_operator` |  | pending |
| `apply_symbolic_boundary_operator` |  | pending |
| `compute_euler_characteristic` |  | pending |
| `verify_boundary_property` |  | pending |
| `verify_coboundary_property` |  | pending |
| `compute_homology_betti_number` |  | pending |
| `compute_cohomology_betti_number` |  | pending |
| `create_grid_complex` |  | pending |
| `create_torus_complex` |  | pending |
| `vietoris_rips_filtration` |  | pending |

## `src/symbolic/transforms.rs` (28)

| legacy function | new home | status |
|---|---|---|
| `fourier_time_shift` | `rules::transforms` `fourier` (shift/modulation theorems applied inside the kernel) | done |
| `fourier_frequency_shift` | `rules::transforms` `fourier` modulation by exp(I a t), cos, sin | done |
| `fourier_scaling` | `rules::transforms` `fourier` (linear arguments in the table) | partial |
| `fourier_differentiation` | `rules::transforms` `fourier` multiplication-by-t theorem; derivative theorem not applied | partial |
| `laplace_time_shift` | `rules::transforms` `laplace` heaviside(t - c) shift theorem | done |
| `laplace_differentiation` | `rules::transforms` `laplace` of `diff(y(t), t)` = s L[y] - y(0) (first order) | done |
| `laplace_frequency_shift` | `rules::transforms` `laplace` exp(a t) shift theorem | done |
| `laplace_scaling` | `rules::transforms` `laplace` (linear arguments in the table) | done |
| `laplace_integration` | `rules::transforms` `laplace` of `defint(g, u, 0, t)` = G/s | done |
| `z_time_shift` | `rules::transforms` `ztransform` of `kronecker(n - k)`; general shifted sequences not detected | partial |
| `z_scaling` | `rules::transforms` `ztransform` a^n g(n) = G(z/a) | done |
| `z_differentiation` | `rules::transforms` `ztransform` n g(n) = -z G'(z) | done |
| `fourier_transform` | `rules::transforms` `fourier` (table + theorems, numerically checked; test fourier_transforms) | done |
| `fourier_transform_internal` | `rules::transforms` `fourier` | done |
| `inverse_fourier_transform` | `rules::transforms` `inverse_fourier` (duality) | done |
| `inverse_fourier_transform_internal` | `rules::transforms` `inverse_fourier` | done |
| `laplace_transform` | `rules::transforms` `laplace` (tests laplace_table_and_theorems, round trips) | done |
| `laplace_transform_internal` | `rules::transforms` `laplace` | done |
| `inverse_laplace_transform` | `rules::transforms` `inverse_laplace` (exact partial fractions over Q; symbolic single linear/quadratic factors) | done |
| `inverse_laplace_transform_internal` | `rules::transforms` `inverse_laplace` | done |
| `z_transform` | `rules::transforms` `ztransform` (unilateral; legacy summed over all integers) | done |
| `z_transform_internal` | `rules::transforms` `ztransform` | done |
| `inverse_z_transform` | `rules::transforms` `inverse_ztransform` (partial fractions of F/z; complex-pole quadratics) | done |
| `inverse_z_transform_internal` | `rules::transforms` `inverse_ztransform` | done |
| `partial_fraction_decomposition` | `rules::poly` `apart(f, x)` over Q | done |
| `partial_fraction_decomposition_internal` | `rules::poly::apart::apart` | done |
| `convolution_fourier` | `rules::transforms` `convolve` builds the convolution integral; the transform-product theorem is not applied | partial |
| `convolution_laplace` | `rules::transforms` `convolve(f, g, t)` = ∫_0^t f(u) g(t-u) du (test convolution) | done |

## `src/symbolic/unit_unification.rs` (1)

| legacy function | new home | status |
|---|---|---|
| `unify_expression` |  | pending |

## `src/symbolic/vector.rs` (15)

| legacy function | new home | status |
|---|---|---|
| `new` | `rules::linalg` `list(...)` | dropped |
| `magnitude` | `rules::linalg` `norm` | done |
| `dot` | `rules::linalg` `dot` | done |
| `cross` | `rules::linalg` `cross` | done |
| `normalize` | `rules::linalg` `normalize` | done |
| `scalar_mul` | `rules::linalg` `smul` / product with list entries | done |
| `angle` | `rules::linalg` `angle` | done |
| `project_onto` | `rules::linalg` `project` | done |
| `to_expr` | `rules::linalg` (vectors are terms) | dropped |
| `gradient` | `rules::linalg` `grad` | done |
| `divergence` | `rules::linalg` `div` | done |
| `curl` | `rules::linalg` `curl` | done |
| `laplacian` | `rules::linalg` `laplacian` | done |
| `directional_derivative` | `rules::linalg` `directional` | done |
| `partial_derivative_vector` | `rules::linalg` `jacobian` column / `diff` per entry | done |

## `src/symbolic/vector_calculus.rs` (8)

| legacy function | new home | status |
|---|---|---|
| `line_integral_scalar` | `rules::linalg` `line_integral` | done |
| `line_integral_scalar_internal` | `rules::linalg` `line_integral` | done |
| `line_integral_vector` | `rules::linalg` `line_integral_vec` | done |
| `line_integral_vector_internal` | `rules::linalg` `line_integral_vec` | done |
| `surface_integral` | `rules::linalg` `surface_integral` (nested quadrature numerically) | done |
| `surface_integral_internal` | `rules::linalg` `surface_integral` | done |
| `volume_integral` | `rules::linalg` `volume_integral` | done |
| `volume_integral_internal` | `rules::linalg` `volume_integral` | done |
