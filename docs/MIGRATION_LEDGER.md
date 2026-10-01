# Migration ledger

Every public function of the pre-rewrite tree (commit `84f3ee89`), and where its
functionality lives now. A row is closed only when the new home is named and a test
exercises it, or when it is dropped with a reason. Status values: `pending`, `done`,
`partial` (say what is missing), `dropped` (say why).

Read a legacy implementation with `git show 84f3ee89:<file>`.

A row with status `partial (no test: ...)` has code but no test naming it, so
it is not closed. A status followed by `(in progress: ...)` is being implemented on
the transforms/complex/finite-field/units branch and is not to be started here.

## Summary (updated 2026-10-01)

| status | rows |
|---|---|
| done | 1453 |
| dropped | 140 |
| partial | 63 (14 missing functionality, 49 missing a test; 4 in progress) |
| pending | 19 (9 in progress) |
| total | 1675 |

The remaining `pending` and `partial` rows are listed by domain in
`docs/MIGRATION_GAPS.md`, which is regenerated from this file.

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
| `with_backend` | `kernels::matrix::Matrix::with_backend` (test tests/kernels/matrix.rs, Faer backend) | done |
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
| `evaluate_action` | `rules::variational` operator `action(L, y(x), path, x, a, b)` (tests rules::variational::tests) | done |
| `euler_lagrange` | `rules::variational` operator `euler_lagrange` (tests rules::variational::tests) (numeric evaluation of the resulting expression through the session) | done |

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
| `aitken_acceleration` | `kernels::convergence::aitken_acceleration` (test tests/kernels/convergence.rs) | done |
| `find_sequence_limit` | `kernels::convergence::find_sequence_limit` (test tests/kernels/convergence.rs) | done |
| `richardson_extrapolation` | `kernels::convergence::richardson_extrapolation` (test tests/kernels/convergence.rs) | done |
| `wynn_epsilon` | `kernels::convergence::wynn_epsilon` (test tests/kernels/convergence.rs) | done |

## `src/numerical/coordinates.rs` (3)

| legacy function | new home | status |
|---|---|---|
| `transform_point` | `kernels::coordinates::transform_point` (test tests/kernels/coordinates.rs) | done |
| `numerical_jacobian` | `kernels::coordinates::numerical_jacobian` (test tests/kernels/coordinates.rs) | done |
| `transform_point_pure` | `kernels::coordinates::transform_point_pure` (test tests/kernels/coordinates.rs) | done |

## `src/numerical/differential_geometry.rs` (5)

| legacy function | new home | status |
|---|---|---|
| `metric_tensor_at_point` | `kernels::differential_geometry::metric_tensor_at_point` (test tests/kernels/differential_geometry.rs) | done |
| `christoffel_symbols` | `kernels::differential_geometry::christoffel_symbols` (test tests/kernels/differential_geometry.rs) | done |
| `riemann_tensor` | `kernels::differential_geometry::riemann_tensor` (test tests/kernels/differential_geometry.rs) | done |
| `ricci_tensor` | `kernels::differential_geometry::ricci_tensor` (test tests/kernels/differential_geometry.rs) | done |
| `ricci_scalar` | `kernels::differential_geometry::ricci_scalar` (test tests/kernels/differential_geometry.rs) | done |

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
| `floor` | `rules::special` operator `floor` (tests rounding_functions, floor_ceil_round_rewrites) | done |
| `ceil` | `rules::special` operator `ceil` (tests rounding_functions, floor_ceil_round_rewrites) | done |
| `round` | `rules::special` operator `round` (tests rounding_functions, floor_ceil_round_rewrites) | done |
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
| `try_closed_form_sum` | `rules::calculus` `antidifference(t, k)` (Gosper: polynomials, c^k, factorials, binomials, rational terms) and `harmonic(n)` for sums of 1/k; missing: closed forms for sin(a k)/cos(a k), ln k (lgamma) and k^a (Hurwitz zeta) | partial |
| `eval_antidiff` | `kernels::indefinite_sum::indefinite_sum` (test tests/kernels/indefinite_sum.rs) (Euler-Maclaurin and Taylor/Bernoulli strategies; the legacy Abel-Plana engine is not reproduced, the normalisation `F(h) = 0` is built in) | done |
| `eval_normalized` | `kernels::indefinite_sum::indefinite_sum` (test tests/kernels/indefinite_sum.rs) (Euler-Maclaurin and Taylor/Bernoulli strategies; the legacy Abel-Plana engine is not reproduced, the normalisation `F(h) = 0` is built in) | done |
| `new` | `kernels::indefinite_sum::indefinite_sum` (test tests/kernels/indefinite_sum.rs) (Euler-Maclaurin and Taylor/Bernoulli strategies; the legacy Abel-Plana engine is not reproduced, the normalisation `F(h) = 0` is built in) | done |
| `eval` | `kernels::indefinite_sum::indefinite_sum` (test tests/kernels/indefinite_sum.rs) (Euler-Maclaurin and Taylor/Bernoulli strategies; the legacy Abel-Plana engine is not reproduced, the normalisation `F(h) = 0` is built in) | done |
| `eval_indefinite_product_numerical` | `kernels::indefinite_sum::indefinite_product` (test tests/kernels/indefinite_sum.rs) | done |
| `series_antidiff` | `kernels::indefinite_sum::series_antidiff` (test tests/kernels/indefinite_sum.rs) | done |
| `compute_taylor_coeffs_numerical` | `kernels::indefinite_sum::numeric_taylor_coefficients` (test tests/kernels/indefinite_sum.rs) | done |
| `eval_indefinite_sum_numerical` | `kernels::indefinite_sum::indefinite_sum` (test tests/kernels/indefinite_sum.rs) | done |
| `expr_contains_var` | helper on the legacy `Expr`; the new kernels take closures and the rule sets use `Graph::free_symbols`/polynomial coefficient extraction | dropped |
| `extract_linear_coeff` | helper on the legacy `Expr`; the new kernels take closures and the rule sets use `Graph::free_symbols`/polynomial coefficient extraction | dropped |

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
| `auto_solve_conjugate_gradient` | `kernels::optimize::EquationOptimizer::auto_solve_conjugate_gradient` (tests auto_solve_conjugate_gradient_minimises_sphere, ..._improves_rosenbrock) | done |
| `solve_with_bfgs` | `kernels::optimize::EquationOptimizer::solve_with_bfgs` | done |
| `solve_with_pso` | `kernels::optimize::EquationOptimizer::solve_with_pso` | done |
| `auto_solve` | `kernels::optimize::EquationOptimizer::auto_solve` (tests auto_solve_runs_a_population_solver, auto_solve_runs_steepest_descent_with_line_search) | done |
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
| `evaluate_power_series` | `kernels::series::evaluate_power_series` (with centre; test tests/kernels/series.rs) | done |
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
| `skewness` | `kernels::stats::skewness` (test skewness_sign) | done |
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
| `find_connected_components` | `kernels::topology::find_connected_components` (test tests/kernels/topology.rs) | done |
| `vietoris_rips_complex` | `kernels::topology::vietoris_rips_complex` (test tests/kernels/topology.rs) | done |
| `betti_numbers_at_radius` | `kernels::topology::betti_numbers_at_radius` (test tests/kernels/topology.rs) | done |
| `compute_persistence` | `kernels::topology::compute_persistence` (test tests/kernels/topology.rs) | done |
| `euclidean_distance` | `kernels::topology::euclidean_distance` (test tests/kernels/topology.rs) | done |

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
| `to_latex` | `io::latex::to_latex` (test src/io/latex.rs constructs) | done |
| `to_latex_prec_with_parens` | `io::latex::to_latex_prec_with_parens` (exercised by the term walker; test constructs) | done |
| `to_greek` | `io::latex::to_greek` (test greek_names) | done |

## `src/output/plotting.rs` (9)

| legacy function | new home | status |
|---|---|---|
| `plot_function_2d` | `io::plot::plot_function_2d` (SVG output; tests in src/io/plot.rs) | done |
| `plot_series_2d` | `io::plot::plot_series_2d` (SVG output; tests in src/io/plot.rs) | done |
| `plot_vector_field_2d` | `io::plot::plot_vector_field_2d` (SVG output; tests in src/io/plot.rs) | done |
| `plot_surface_3d` | `io::plot::plot_surface_3d` (SVG output; tests in src/io/plot.rs) | done |
| `plot_surface_2d` | `io::plot::plot_surface_2d` (SVG output; tests in src/io/plot.rs) | done |
| `plot_parametric_curve_3d` | `io::plot::plot_parametric_curve_3d` (SVG output; tests in src/io/plot.rs) | done |
| `plot_vector_field_3d` | `io::plot::plot_vector_field_3d` (SVG output; tests in src/io/plot.rs) | done |
| `plot_3d_path_from_points` | `io::plot::plot_3d_path_from_points` (SVG output; tests in src/io/plot.rs) | done |
| `plot_heatmap_2d` | `io::plot::plot_heatmap_2d` (SVG output; tests in src/io/plot.rs) | done |

## `src/output/pretty_print.rs` (1)

| legacy function | new home | status |
|---|---|---|
| `pretty_print` | `Graph::display` (graph/print.rs tests; every rule-set test compares printed answers) | done |

## `src/output/typst.rs` (1)

| legacy function | new home | status |
|---|---|---|
| `to_typst` | `io::typst::to_typst` (test src/io/typst.rs) | done |

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
| `simulate_black_hole_orbits_scenario` | `sim::models::geodesic_relativity::simulate_black_hole_orbits_scenario` | done |

## `src/physics/physics_sim/gpe_superfluidity.rs` (2)

| legacy function | new home | status |
|---|---|---|
| `run_gpe_ground_state_finder` | `sim::models::gpe_superfluidity::run_gpe_ground_state_finder` | done |
| `simulate_bose_einstein_vortex_scenario` | `sim::models::gpe_superfluidity::simulate_bose_einstein_vortex_scenario` | done |

## `src/physics/physics_sim/ising_statistical.rs` (2)

| legacy function | new home | status |
|---|---|---|
| `run_ising_simulation` | `sim::models::ising_statistical::run_ising_simulation` | done |
| `simulate_ising_phase_transition_scenario` | `sim::models::ising_statistical::simulate_ising_phase_transition_scenario` | done |

## `src/physics/physics_sim/linear_elasticity.rs` (3)

| legacy function | new home | status |
|---|---|---|
| `element_stiffness_matrix` | `sim::models::linear_elasticity::element_stiffness_matrix` | done |
| `run_elasticity_simulation` | `sim::models::linear_elasticity::run_elasticity_simulation` | done |
| `simulate_cantilever_beam_scenario` | `sim::models::linear_elasticity::simulate_cantilever_beam_scenario` | done |

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
| `simulate_double_slit_scenario` | `sim::models::schrodinger_quantum::simulate_double_slit_scenario` | done |

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
| `ifft3d` | `sim::physics_sm::ifft3d` | done |
| `solve_advection_diffusion_3d` | `sim::physics_sm::solve_advection_diffusion_3d` | done |
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
| `cad` | `rules::poly::algebra` operator `cad(polys, vars)`: sample points of the cells for R^1 and R^2 only; no 3+ variables (legacy projected and lifted to any dimension), no cell adjacency/structure output | partial |

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
| `check_analytic` |  | pending (in progress: transforms/complex/finite-field/units branch) |
| `find_poles` | `rules::complex::analysis` operator `poles(f, z)` (test rules::complex::analysis::tests) | done |
| `calculate_residue` | `rules::complex::analysis` operator `residue(f, z, a)` (test rules::complex::analysis::tests) | done |
| `is_inside_contour` |  | pending (in progress: transforms/complex/finite-field/units branch) |
| `path_integrate` | `rules::complex::analysis` operator `contour_integral(f, z, path(g(t), t, t0, t1))` (test rules::complex::analysis::tests) | done |
| `factorial` | `rules::combinatorics` `factorial` (test factorials) | done |
| `improper_integral` | `defint` over `oo` limits (antiderivative or `Quadrature`, test numeric_quadrature_without_a_closed_form); `poles`/`residue` exist, but no operator evaluates `∫_{-oo}^{oo}` of a rational function as `2 pi I * sum of upper-half-plane residues` | partial |
| `limit` | `rules::calculus` `limit` (tests limits_by_continuity_and_cancellation, limits_at_infinity) | done |
| `limit_internal` | `rules::calculus` `limit` (tests limits_by_continuity_and_cancellation, limits_at_infinity) | done |

## `src/symbolic/calculus_of_variations.rs` (5)

| legacy function | new home | status |
|---|---|---|
| `euler_lagrange` | `rules::variational` operator `euler_lagrange` (tests rules::variational::tests) | done |
| `euler_lagrange_internal` | `rules::variational` operator `euler_lagrange` (tests rules::variational::tests) | done |
| `solve_euler_lagrange` | `rules::variational` operator `solve_euler_lagrange` (tests rules::variational::tests) | done |
| `solve_euler_lagrange_internal` | `rules::variational` operator `solve_euler_lagrange` (tests rules::variational::tests) | done |
| `hamiltons_principle` | `rules::variational` operator `hamiltons_principle` | partial (no test: `hamiltons_principle` is defined but no test exercises it) |

## `src/symbolic/cas_foundations.rs` (12)

| legacy function | new home | status |
|---|---|---|
| `get_term_factors` | Expr helper; terms are normalised by `rules::arith` Collect | dropped |
| `build_expr_from_factors` | Expr helper; terms are normalised by `rules::arith` Collect | dropped |
| `normalize` | `rules::arith` Collect window pass (test like_terms_are_collected) | done |
| `expand` | `rules::poly` `expand` (test expand_is_pinned_against_cheaper_spellings) | done |
| `factorize` | `rules::poly` `factor` (tests factor_univariate, factor_multivariate_pulls_out_content) | done |
| `factorize_internal` | `rules::poly` `factor` (tests factor_univariate, factor_multivariate_pulls_out_content) | done |
| `risch_integrate` | legacy function was a deprecated placeholder returning an unevaluated `RischIntegrate(..)` variable; integration is `rules::calculus` `integral` | dropped |
| `grobner_basis` | `rules::poly` `groebner` (test groebner_bases) | done |
| `cylindrical_algebraic_decomposition` | legacy function was a deprecated placeholder; the real implementation is `cad` (rules::poly::algebra) | dropped |
| `simplify_with_relations` | `rules::poly::algebra` operator `simplify_with_relations(e, relations, vars)` (tests rules::poly::algebra::tests) | done |
| `simplify_with_relations_internal` | `rules::poly::algebra` operator `simplify_with_relations(e, relations, vars)` (tests rules::poly::algebra::tests) | done |
| `normalize_with_relations` | `rules::poly::algebra` operator `normal_form(p, polys, vars)` (tests rules::poly::algebra::tests) and `simplify_with_relations` | done |

## `src/symbolic/classical_mechanics.rs` (19)

| legacy function | new home | status |
|---|---|---|
| `new` | `rules::physics` operator `kinematics` (test rules::physics::tests) | done |
| `newtons_second_law` | `rules::physics` operator `newtons_second_law` | partial (no test: `newtons_second_law` is defined but no test exercises it) |
| `momentum` | `rules::physics` operator `momentum` (test rules::physics::tests) | done |
| `kinetic_energy` | `rules::physics` operator `kinetic_energy` (test rules::physics::tests) | done |
| `potential_energy_gravity_uniform` | `rules::physics` operator `potential_energy_gravity_uniform` | partial (no test: `potential_energy_gravity_uniform` is defined but no test exercises it) |
| `potential_energy_gravity_universal` | `rules::physics` operator `potential_energy_gravity_universal` (test rules::physics::tests) | done |
| `potential_energy_spring` | `rules::physics` operator `potential_energy_spring` (test rules::physics::tests) | done |
| `work_constant_force` | `rules::physics` operator `work_constant_force` (test rules::physics::tests) | done |
| `work_line_integral` | `rules::physics` operator `work_line_integral` | partial (no test: `work_line_integral` is defined but no test exercises it) |
| `power` | `rules::physics` operator `power` (test rules::physics::tests) | done |
| `torque` | `rules::physics` operator `torque` (test rules::physics::tests) | done |
| `angular_momentum` | `rules::physics` operator `angular_momentum` | partial (no test: `angular_momentum` is defined but no test exercises it) |
| `centripetal_acceleration` | `rules::physics` operator `centripetal_acceleration` | partial (no test: `centripetal_acceleration` is defined but no test exercises it) |
| `moment_of_inertia_point_mass` | `rules::physics` operator `moment_of_inertia_point_mass` | partial (no test: `moment_of_inertia_point_mass` is defined but no test exercises it) |
| `rotational_kinetic_energy` | `rules::physics` operator `rotational_kinetic_energy` | partial (no test: `rotational_kinetic_energy` is defined but no test exercises it) |
| `lagrangian` | `rules::physics` operator `lagrangian` (test rules::physics::tests) | done |
| `hamiltonian` | `rules::physics` operator `hamiltonian` | partial (no test: `hamiltonian` is defined but no test exercises it) |
| `euler_lagrange_equation` | `rules::physics` operator `euler_lagrange_equation` (test rules::physics::tests) | done |
| `poisson_bracket` | `rules::physics` operator `poisson_bracket` (test rules::physics::tests) | done |

## `src/symbolic/combinatorics.rs` (10)

| legacy function | new home | status |
|---|---|---|
| `expand_binomial` | `expand((a+b)^n)` for literal n; symbolic n has no `expand_binomial(a, b, n)` operator returning `sum(binomial(n, k) a^(n-k) b^k, k, 0, n)` (`sum` and `binomial` exist; binomial sums are closed by rules::calculus recurrence guessing) | partial |
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
| `translation_2d` | `rules::discrete::graphics` operator `translation_2d` (tests rules::discrete::graphics::tests) | done |
| `translation_3d` | `rules::discrete::graphics` operator `translation_3d` (tests rules::discrete::graphics::tests) | done |
| `rotation_2d` | `rules::discrete::graphics` operator `rotation_2d` (tests rules::discrete::graphics::tests) | done |
| `rotation_3d_x` | `rules::discrete::graphics` operator `rotation_3d_x` (tests rules::discrete::graphics::tests) | done |
| `rotation_3d_y` | `rules::discrete::graphics` operator `rotation_3d_y` (tests rules::discrete::graphics::tests) | done |
| `rotation_3d_z` | `rules::discrete::graphics` operator `rotation_3d_z` (tests rules::discrete::graphics::tests) | done |
| `scaling_2d` | `rules::discrete::graphics` operator `scaling_2d` (tests rules::discrete::graphics::tests) | done |
| `scaling_3d` | `rules::discrete::graphics` operator `scaling_3d` (tests rules::discrete::graphics::tests) | done |
| `perspective_projection` | `rules::discrete::graphics` operator `perspective` (tests rules::discrete::graphics::tests) | done |
| `orthographic_projection` | `rules::discrete::graphics` operator `orthographic` (tests rules::discrete::graphics::tests) | done |
| `look_at` | `rules::discrete::graphics` operator `look_at` (tests rules::discrete::graphics::tests) | done |
| `evaluate` | `rules::discrete::graphics` operator `bezier` (and `bspline` for B-spline curves) (tests rules::discrete::graphics::tests) | done |
| `derivative` | `rules::discrete::graphics` operator `bezier_derivative` (tests rules::discrete::graphics::tests) | done |
| `split` | `rules::discrete::graphics` operator `bezier_split` (tests rules::discrete::graphics::tests) | done |
| `new` | `rules::discrete::graphics` operator `mesh_transform` (a mesh is the pair of a vertex list and a face list; see also `mesh_normals`, `mesh_triangulate`) (tests rules::discrete::graphics::tests) | done |
| `apply_transformation` | `rules::discrete::graphics` operator `mesh_transform` (and `apply_transform` for single points) (tests rules::discrete::graphics::tests) | done |
| `compute_normals` | `rules::discrete::graphics` operator `mesh_normals` (tests rules::discrete::graphics::tests) | done |
| `triangulate` | `rules::discrete::graphics` operator `mesh_triangulate` (tests rules::discrete::graphics::tests) | done |
| `shear_2d` | `rules::discrete::graphics` operator `shear_2d` (tests rules::discrete::graphics::tests) | done |
| `reflection_2d` | `rules::discrete::graphics` operator `reflection_2d` (tests rules::discrete::graphics::tests) | done |
| `reflection_3d` | `rules::discrete::graphics` operator `reflection_3d` (tests rules::discrete::graphics::tests) | done |
| `rotation_axis_angle` | `rules::discrete::graphics` operator `rotation_axis_angle` (tests rules::discrete::graphics::tests) | done |

## `src/symbolic/convergence.rs` (1)

| legacy function | new home | status |
|---|---|---|
| `analyze_convergence` | `rules::calculus` operator `converges(term, k)` (divergence, ratio, root, alternating, limit-comparison and Cauchy-condensation tests; tests in rules::calculus::tests) | done |

## `src/symbolic/coordinates.rs` (12)

| legacy function | new home | status |
|---|---|---|
| `transform_point` | `rules::geometry` operator `transform_point` (tests rules::geometry::tests) | done |
| `transform_expression` | `rules::geometry` operator `transform_expression` (tests rules::geometry::tests) | done |
| `get_transform_rules` | `rules::geometry` operator `to_cartesian` and `from_cartesian` (tests rules::geometry::tests) | done |
| `get_to_cartesian_rules` | `rules::geometry` operator `to_cartesian` (tests rules::geometry::tests) | done |
| `transform_contravariant_vector` | `rules::geometry` operator `transform_vector` (tests rules::geometry::tests) | done |
| `transform_covariant_vector` | `rules::geometry` operator `transform_covector` (tests rules::geometry::tests) | done |
| `transform_tensor2` | `rules::geometry` operator `transform_tensor2` (tests rules::geometry::tests) | done |
| `symbolic_mat_mat_mul` | `rules::linalg` operator `matmul` | done |
| `get_metric_tensor` | `rules::geometry` operator `coordinate_metric` (tests rules::geometry::tests) | done |
| `transform_divergence` | `rules::geometry` operator `div_in` (tests rules::geometry::tests) | done |
| `transform_curl` | `rules::geometry` operator `curl_in` (tests rules::geometry::tests) | done |
| `transform_gradient` | `rules::geometry` operator `grad_in` (tests rules::geometry::tests) | done |

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
| `is_infinity` | `rules::discrete::crypto` operator `ec_is_infinity` (tests rules::discrete::crypto::tests) | done |
| `x` | `rules::discrete::crypto` operator `ec_x` (tests rules::discrete::crypto::tests) | done |
| `y` | `rules::discrete::crypto` operator `ec_y` (tests rules::discrete::crypto::tests) | done |
| `new` | `rules::discrete::crypto` operator `ec_curve` (curve term; points are `list(x, y)`) (tests rules::discrete::crypto::tests) | done |
| `is_on_curve` | `rules::discrete::crypto` operator `ec_on_curve` (tests rules::discrete::crypto::tests) | done |
| `negate` | `rules::discrete::crypto` operator `ec_neg` (tests rules::discrete::crypto::tests) | done |
| `double` | `rules::discrete::crypto` operator `ec_double` (tests rules::discrete::crypto::tests) | done |
| `add` | `rules::discrete::crypto` operator `ec_add` (tests rules::discrete::crypto::tests) | done |
| `scalar_mult` | `rules::discrete::crypto` operator `ec_mul` (tests rules::discrete::crypto::tests) | done |
| `generate_keypair` | `rules::discrete::crypto` operator `ecdh_public` (the private key `d` is an argument, not drawn at random) (tests rules::discrete::crypto::tests) | done |
| `generate_shared_secret` | `rules::discrete::crypto` operator `ecdh_shared` (tests rules::discrete::crypto::tests) | done |
| `point_compress` | `rules::discrete::crypto` operator `ec_compress` (tests rules::discrete::crypto::tests) | done |
| `point_decompress` | `rules::discrete::crypto` operator `ec_decompress` (tests rules::discrete::crypto::tests) | done |
| `ecdsa_sign` | `rules::discrete::crypto` operator `ecdsa_sign` (nonce `k` supplied) (tests rules::discrete::crypto::tests) | done |
| `ecdsa_verify` | `rules::discrete::crypto` operator `ecdsa_verify` (tests rules::discrete::crypto::tests) | done |

## `src/symbolic/differential_geometry.rs` (7)

| legacy function | new home | status |
|---|---|---|
| `exterior_derivative` | `rules::geometry` operator `exterior_d` (tests rules::geometry::tests) | done |
| `wedge_product` | `rules::geometry` operator `wedge` (tests rules::geometry::tests) | done |
| `boundary` | no symbolic region/manifold boundary term `∂M`; the concrete theorems take explicit bounds (rectangle, box, parametrised surface) | pending |
| `generalized_stokes_theorem` | no operator stating `∫_M dω = ∫_∂M ω` for a symbolic manifold; `greens_theorem`, `gauss_theorem`, `stokes_theorem` cover the concrete cases | pending |
| `gauss_theorem` | `rules::geometry` operator `gauss_theorem` (tests rules::geometry::tests) | done |
| `stokes_theorem` | `rules::geometry` operator `stokes_theorem` (tests rules::geometry::tests) | done |
| `greens_theorem` | `rules::geometry` operator `greens_theorem` (tests rules::geometry::tests) | done |

## `src/symbolic/discrete_groups.rs` (4)

| legacy function | new home | status |
|---|---|---|
| `cyclic_group` | `rules::discrete::groups` operator `cyclic_group` (tests rules::discrete::groups::tests) | done |
| `dihedral_group` | `rules::discrete::groups` operator `dihedral_group` (tests rules::discrete::groups::tests) | done |
| `symmetric_group` | `rules::discrete::groups` operator `symmetric_group` (tests rules::discrete::groups::tests) | done |
| `klein_four_group` | `rules::discrete::groups` operator `klein_four_group` (tests rules::discrete::groups::tests) | done |

## `src/symbolic/electromagnetism.rs` (8)

| legacy function | new home | status |
|---|---|---|
| `new` | `rules::physics` operator `maxwell_equations` (test rules::physics::tests) | done |
| `lorentz_force` | `rules::physics` operator `lorentz_force` (test rules::physics::tests) | done |
| `electric_field_from_potentials` | `rules::physics` operator `electric_field_from_potentials` | partial (no test: `electric_field_from_potentials` is defined but no test exercises it) |
| `electric_field_from_potential` | `rules::physics` operator `electric_field_from_potential` (test rules::physics::tests) | done |
| `magnetic_field_from_vector_potential` | `rules::physics` operator `magnetic_field_from_vector_potential` (test rules::physics::tests) | done |
| `poynting_vector` | `rules::physics` operator `poynting_vector` (test rules::physics::tests) | done |
| `energy_density` | `rules::physics` operator `em_energy_density` | partial (no test: `em_energy_density` is defined but no test exercises it) |
| `coulombs_law` | `rules::physics` operator `coulombs_law` | partial (no test: `coulombs_law` is defined but no test exercises it) |

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
| `expand` | `rules::poly` `expand` (polynomials) and `expand_trig` (sum and multiple angles, rules::poly::algebra tests) | done |
| `expand_internal` | `rules::poly` `expand` (polynomials) and `expand_trig` (sum and multiple angles, rules::poly::algebra tests) | done |
| `binomial_coefficient` | `rules::combinatorics` `binomial` (test binomials) | done |

## `src/symbolic/error_correction.rs` (13)

| legacy function | new home | status |
|---|---|---|
| `hamming_distance` | `rules::discrete::coding` operator `hamming_distance` (tests rules::discrete::coding::tests) | done |
| `hamming_weight` | `rules::discrete::coding` operator `hamming_weight` (tests rules::discrete::coding::tests) | done |
| `hamming_encode` | `rules::discrete::coding` operator `hamming_encode` (tests rules::discrete::coding::tests) | done |
| `hamming_check` | `rules::discrete::coding` operator `hamming_check` (tests rules::discrete::coding::tests) | done |
| `hamming_decode` | `rules::discrete::coding` operator `hamming_decode` (tests rules::discrete::coding::tests) | done |
| `rs_encode` | `rules::discrete::coding` operator `rs_encode` (tests rules::discrete::coding::tests) | done |
| `rs_check` | `rules::discrete::coding` operator `rs_check` (tests rules::discrete::coding::tests) | done |
| `rs_error_count` | `rules::discrete::coding` operator `rs_error_count` (tests rules::discrete::coding::tests) | done |
| `rs_decode` | `rules::discrete::coding` operator `rs_decode` (tests rules::discrete::coding::tests) | done |
| `crc32_compute` | `rules::discrete::coding` operator `crc32` (tests rules::discrete::coding::tests) | done |
| `crc32_verify` | `rules::discrete::coding` operator `crc32_verify` (tests rules::discrete::coding::tests) | done |
| `crc32_update` | `rules::discrete::coding` operator `crc32_update` (tests rules::discrete::coding::tests) | done |
| `crc32_finalize` | `rules::discrete::coding` operator `crc32_finalize` (tests rules::discrete::coding::tests) | done |

## `src/symbolic/error_correction_helper.rs` (23)

| legacy function | new home | status |
|---|---|---|
| `new` | `rules::discrete::finite_field` operator `gf_add` (the field GF(p) is the modulus argument `p`; there is no field object) (tests rules::discrete::finite_field::tests) | done |
| `from_bigint` | `rules::discrete::finite_field` operator `gf_add` (arbitrary-precision modulus `p`) (tests rules::discrete::finite_field::tests) | done |
| `is_zero` | `rules::discrete::finite_field` operator `gf_is_zero` (tests rules::discrete::finite_field::tests) | done |
| `is_one` | `rules::discrete::finite_field` operator `gf_is_one` (tests rules::discrete::finite_field::tests) | done |
| `inverse` | `rules::discrete::finite_field` operator `gf_inv` (tests rules::discrete::finite_field::tests) | done |
| `pow` | `rules::discrete::finite_field` operator `gf_pow` (tests rules::discrete::finite_field::tests) | done |
| `gf256_exp` | `rules::discrete::coding` operator `gf256_exp` (tests rules::discrete::coding::tests) | done |
| `gf256_log` | `rules::discrete::coding` operator `gf256_log` (tests rules::discrete::coding::tests) | done |
| `gf256_add` | `rules::discrete::coding` operator `gf256_add` (tests rules::discrete::coding::tests) | done |
| `gf256_mul` | `rules::discrete::coding` operator `gf256_mul` (tests rules::discrete::coding::tests) | done |
| `gf256_inv` | `rules::discrete::coding` operator `gf256_inv` (tests rules::discrete::coding::tests) | done |
| `gf256_div` | `rules::discrete::coding` operator `gf256_div` (tests rules::discrete::coding::tests) | done |
| `gf256_pow` | `rules::discrete::coding` operator `gf256_pow` (tests rules::discrete::coding::tests) | done |
| `poly_eval_gf256` | `rules::discrete::coding` operator `gf256_poly_eval` (tests rules::discrete::coding::tests) | done |
| `poly_add_gf256` | `rules::discrete::coding` operator `gf256_poly_add` (tests rules::discrete::coding::tests) | done |
| `poly_mul_gf256` | `rules::discrete::coding` operator `gf256_poly_mul` (tests rules::discrete::coding::tests) | done |
| `poly_scale_gf256` | `rules::discrete::coding` operator `gf256_poly_scale` (tests rules::discrete::coding::tests) | done |
| `poly_derivative_gf256` | `rules::discrete::coding` operator `gf256_poly_derivative` (tests rules::discrete::coding::tests) | done |
| `poly_gcd_gf256` | `rules::discrete::coding` operator `gf256_poly_gcd` (tests rules::discrete::coding::tests) | done |
| `poly_div_gf256` | `rules::discrete::coding` operator `gf256_poly_divmod` (and `gf256_poly_mod`) (tests rules::discrete::coding::tests) | done |
| `poly_add_gf` | `rules::discrete::finite_field` operator `gfp_add` (tests rules::discrete::finite_field::tests) | done |
| `poly_mul_gf` | `rules::discrete::finite_field` operator `gfp_mul` (tests rules::discrete::finite_field::tests) | done |
| `poly_div_gf` | `rules::discrete::finite_field` operator `gfp_divmod` (tests rules::discrete::finite_field::tests) | done |

## `src/symbolic/finite_field.rs` (10)

| legacy function | new home | status |
|---|---|---|
| `serialize` | serde helpers for `Arc<PrimeField>` fields of the legacy element structs; field elements are plain terms (integers, coefficient lists) now | dropped |
| `deserialize` | serde helpers for `Arc<PrimeField>` fields of the legacy element structs; field elements are plain terms (integers, coefficient lists) now | dropped |
| `new` | `rules::discrete::finite_field` operator `gfx_reduce` (elements of GF(p)[x]/(m) reduced by `gfx_reduce`; `gfp_norm` for polynomials) (tests rules::discrete::finite_field::tests) | done |
| `inverse` | `rules::discrete::finite_field` operator `gfx_inv` (and `gf_inv`) (tests rules::discrete::finite_field::tests) | done |
| `degree` | `rules::discrete::finite_field` operator `gfp_degree` (tests rules::discrete::finite_field::tests) | done |
| `long_division` | `rules::discrete::finite_field` operator `gfp_divmod` (tests rules::discrete::finite_field::tests) | done |
| `add` | `rules::discrete::finite_field` operator `gfx_add` (and `gf_add`, `gfp_add`) (tests rules::discrete::finite_field::tests) | done |
| `sub` | `rules::discrete::finite_field` operator `gfx_sub` (and `gf_sub`, `gfp_sub`) (tests rules::discrete::finite_field::tests) | done |
| `mul` | `rules::discrete::finite_field` operator `gfx_mul` (and `gf_mul`, `gfp_mul`) (tests rules::discrete::finite_field::tests) | done |
| `div` | `rules::discrete::finite_field` operator `gfx_div` (and `gf_div`) (tests rules::discrete::finite_field::tests) | done |

## `src/symbolic/fractal_geometry_and_chaos.rs` (12)

| legacy function | new home | status |
|---|---|---|
| `new` | `rules::discrete::fractal` operator `ifs_apply` (an IFS is a list of coordinate formulas; `ifs_generate` runs the chaos game) (tests rules::discrete::fractal::tests) | done |
| `apply` | `rules::discrete::fractal` operator `ifs_apply` (tests rules::discrete::fractal::tests) | done |
| `similarity_dimension` | `rules::discrete::fractal` operator `similarity_dimension` (and `moran_dimension`) (tests rules::discrete::fractal::tests) | done |
| `new_mandelbrot_family` | `rules::discrete::fractal` operator `mandelbrot_iterate` (the family `z^2 + c` with symbolic `c`; see also `mandelbrot_orbit`) (tests rules::discrete::fractal::tests) | done |
| `iterate` | `rules::discrete::fractal` operator `mandelbrot_iterate` (tests rules::discrete::fractal::tests) | done |
| `orbit` | `rules::discrete::fractal` operator `mandelbrot_orbit` (tests rules::discrete::fractal::tests) | done |
| `fixed_points` | `rules::discrete::fractal` operator `mandelbrot_fixed_points` (tests rules::discrete::fractal::tests) | done |
| `stability_index` | `rules::discrete::fractal` operator `mandelbrot_stability` (tests rules::discrete::fractal::tests) | done |
| `find_fixed_points` | `rules::discrete::fractal` operator `map_fixed_points` (and `complex_map_fixed_points`) (tests rules::discrete::fractal::tests) | done |
| `analyze_stability` | `rules::discrete::fractal` operator `map_stability` (and `complex_map_stability`) (tests rules::discrete::fractal::tests) | done |
| `lyapunov_exponent` | `rules::discrete::fractal` operator `lyapunov_exponent` (tests rules::discrete::fractal::tests) | done |
| `lorenz_system` | `rules::discrete::fractal` operator `lorenz` (tests rules::discrete::fractal::tests) | done |

## `src/symbolic/functional_analysis.rs` (13)

| legacy function | new home | status |
|---|---|---|
| `new` | function spaces are the trailing `x, a, b` arguments of `rules::functional` `inner_product`/`l2_norm`/`lp_norm` (Hilbert and Banach space structs have no term form); operators are built from `op_mul`, `op_d`, `op_int`, ... | done |
| `apply` | `rules::functional` operator `op_apply` (tests rules::functional::tests) | done |
| `inner_product` | `rules::functional` operator `inner_product(f, g, x, a, b)` (tests rules::functional::tests) (discrete samples: kernels::functional_analysis::inner_product) | done |
| `inner_product_internal` | `rules::functional` operator `inner_product(f, g, x, a, b)` (tests rules::functional::tests) (discrete samples: kernels::functional_analysis::inner_product) | done |
| `norm` | `rules::functional` operator `l2_norm` (tests rules::functional::tests) and `lp_norm` (discrete samples: kernels::functional_analysis::{l1_norm,l2_norm,infinity_norm}) | done |
| `norm_internal` | `rules::functional` operator `l2_norm` (tests rules::functional::tests) and `lp_norm` (discrete samples: kernels::functional_analysis::{l1_norm,l2_norm,infinity_norm}) | done |
| `banach_norm` | `rules::functional` operator `lp_norm(f, p, x, a, b)` (tests rules::functional::tests) | done |
| `banach_norm_internal` | `rules::functional` operator `lp_norm(f, p, x, a, b)` (tests rules::functional::tests) | done |
| `are_orthogonal` | `rules::functional` operator `are_orthogonal` (tests rules::functional::tests) | done |
| `project` | `rules::functional` operator `project_onto` (tests rules::functional::tests) (discrete samples: kernels::functional_analysis::project) | done |
| `project_internal` | `rules::functional` operator `project_onto` (tests rules::functional::tests) (discrete samples: kernels::functional_analysis::project) | done |
| `gram_schmidt` | `rules::functional` operator `gram_schmidt` (tests rules::functional::tests) (discrete samples: kernels::functional_analysis::gram_schmidt) | done |
| `gram_schmidt_orthonormal` | `rules::functional` operator `gram_schmidt_orthonormal` (tests rules::functional::tests) (discrete samples: kernels::functional_analysis::gram_schmidt_orthonormal) | done |

## `src/symbolic/geometric_algebra.rs` (11)

| legacy function | new home | status |
|---|---|---|
| `new` | `rules::geometric_algebra` operator `mv` (tests rules::geometric_algebra::tests) (any signature, symbolic coefficients) | done |
| `scalar` | `rules::geometric_algebra` operator `mv_scalar` (tests rules::geometric_algebra::tests) | done |
| `vector` | `rules::geometric_algebra` operator `mv_vector` (tests rules::geometric_algebra::tests) | done |
| `geometric_product` | `rules::geometric_algebra` operator `ga_gp` (tests rules::geometric_algebra::tests) | done |
| `grade_projection` | `rules::geometric_algebra` operator `ga_grade` (tests rules::geometric_algebra::tests) | done |
| `outer_product` | `rules::geometric_algebra` operator `ga_wedge` (tests rules::geometric_algebra::tests) | done |
| `inner_product` | `rules::geometric_algebra` operator `ga_inner` (tests rules::geometric_algebra::tests) | done |
| `reverse` | `rules::geometric_algebra` operator `ga_reverse` (tests rules::geometric_algebra::tests) | done |
| `magnitude` | `rules::geometric_algebra` operator `ga_norm` (tests rules::geometric_algebra::tests) | done |
| `dual` | `rules::geometric_algebra` operator `ga_dual` (tests rules::geometric_algebra::tests) | done |
| `normalize` | `rules::geometric_algebra` operator `ga_normalize` (tests rules::geometric_algebra::tests) | done |

## `src/symbolic/graph.rs` (15)

| legacy function | new home | status |
|---|---|---|
| `new` | `rules::discrete::graphs` operator `graph` / `digraph` (terms `graph(n, edges)`) (tests rules::discrete::graphs::tests) | done |
| `nodes` | `rules::discrete::graphs` operator `graph_nodes` (tests rules::discrete::graphs::tests) | done |
| `node_count` | `rules::discrete::graphs` operator `graph_node_count` (tests rules::discrete::graphs::tests) | done |
| `is_directed` | `rules::discrete::graphs` operator `graph_is_directed` (tests rules::discrete::graphs::tests) | done |
| `add_node` | `rules::discrete::graphs` operator `graph_add_node` (tests rules::discrete::graphs::tests) | done |
| `add_edge` | `rules::discrete::graphs` operator `graph_add_edge` (tests rules::discrete::graphs::tests) | done |
| `get_node_id` | `rules::discrete::graphs` operator `graph_node_id` (tests rules::discrete::graphs::tests) | done |
| `neighbors` | `rules::discrete::graphs` operator `graph_neighbors` (tests rules::discrete::graphs::tests) | done |
| `out_degree` | `rules::discrete::graphs` operator `graph_out_degree` (tests rules::discrete::graphs::tests) | done |
| `in_degree` | `rules::discrete::graphs` operator `graph_in_degree` (tests rules::discrete::graphs::tests) | done |
| `get_edges` | `rules::discrete::graphs` operator `graph_edges` (tests rules::discrete::graphs::tests) | done |
| `add_hyperedge` | `rules::discrete::graphs` operator `graph_add_hyperedge` (tests rules::discrete::graphs::tests) | done |
| `to_adjacency_matrix` | `rules::discrete::graphs` operator `graph_adjacency` (tests rules::discrete::graphs::tests) | done |
| `to_incidence_matrix` | `rules::discrete::graphs` operator `graph_incidence` (tests rules::discrete::graphs::tests) | done |
| `to_laplacian_matrix` | `rules::discrete::graphs` operator `graph_laplacian` (tests rules::discrete::graphs::tests) | done |

## `src/symbolic/graph_algorithms.rs` (26)

| legacy function | new home | status |
|---|---|---|
| `dfs` | `rules::discrete::graphs` operator `graph_dfs` (tests rules::discrete::graphs::tests) | done |
| `bfs` | `rules::discrete::graphs` operator `graph_bfs` (tests rules::discrete::graphs::tests) | done |
| `connected_components` | `rules::discrete::graphs` operator `graph_components` (tests rules::discrete::graphs::tests) | done |
| `is_connected` | `rules::discrete::graphs` operator `graph_is_connected` (tests rules::discrete::graphs::tests) | done |
| `strongly_connected_components` | `rules::discrete::graphs` operator `graph_scc` (tests rules::discrete::graphs::tests) | done |
| `has_cycle` | `rules::discrete::graphs` operator `graph_has_cycle` (tests rules::discrete::graphs::tests) | done |
| `find_bridges_and_articulation_points` | `rules::discrete::graphs` operator `graph_bridges` (tests rules::discrete::graphs::tests) | done |
| `kruskal_mst` | `rules::discrete::graphs` operator `graph_kruskal` (tests rules::discrete::graphs::tests) | done |
| `edmonds_karp_max_flow` | `rules::discrete::graphs` operator `graph_edmonds_karp` (tests rules::discrete::graphs::tests) | done |
| `dinic_max_flow` | `rules::discrete::graphs` operator `graph_dinic` (tests rules::discrete::graphs::tests) | done |
| `bellman_ford` | `rules::discrete::graphs` operator `graph_bellman_ford` (tests rules::discrete::graphs::tests) | done |
| `min_cost_max_flow` | `rules::discrete::graphs` operator `graph_min_cost_flow` (tests rules::discrete::graphs::tests) | done |
| `is_bipartite` | `rules::discrete::graphs` operator `graph_is_bipartite` (tests rules::discrete::graphs::tests) | done |
| `bipartite_maximum_matching` | `rules::discrete::graphs` operator `graph_bipartite_matching` (tests rules::discrete::graphs::tests) | done |
| `prim_mst` | `rules::discrete::graphs` operator `graph_prim` (tests rules::discrete::graphs::tests) | done |
| `topological_sort_kahn` | `rules::discrete::graphs` operator `graph_toposort_kahn` (tests rules::discrete::graphs::tests) | done |
| `topological_sort_dfs` | `rules::discrete::graphs` operator `graph_toposort_dfs` (tests rules::discrete::graphs::tests) | done |
| `topological_sort` | `rules::discrete::graphs` operator `graph_toposort` (tests rules::discrete::graphs::tests) | done |
| `bipartite_minimum_vertex_cover` | `rules::discrete::graphs` operator `graph_vertex_cover` (tests rules::discrete::graphs::tests) | done |
| `hopcroft_karp_bipartite_matching` | `rules::discrete::graphs` operator `graph_hopcroft_karp` (tests rules::discrete::graphs::tests) | done |
| `blossom_algorithm` | `rules::discrete::graphs` operator `graph_blossom` (tests rules::discrete::graphs::tests) | done |
| `shortest_path_unweighted` | `rules::discrete::graphs` operator `graph_shortest_path_unweighted` (tests rules::discrete::graphs::tests) | done |
| `dijkstra` | `rules::discrete::graphs` operator `graph_dijkstra` (tests rules::discrete::graphs::tests) | done |
| `floyd_warshall` | `rules::discrete::graphs` operator `graph_floyd_warshall` (tests rules::discrete::graphs::tests) | done |
| `spectral_analysis` | `rules::discrete::graphs` operator `spectral_analysis` (tests rules::discrete::graphs::tests) | done |
| `algebraic_connectivity` | `rules::discrete::graphs` operator `algebraic_connectivity` (tests rules::discrete::graphs::tests) | done |

## `src/symbolic/graph_isomorphism_and_coloring.rs` (3)

| legacy function | new home | status |
|---|---|---|
| `are_isomorphic_heuristic` | `rules::discrete::graphs` operator `graph_isomorphic_heuristic` (tests rules::discrete::graphs::tests) | done |
| `greedy_coloring` | `rules::discrete::graphs` operator `graph_greedy_coloring` (tests rules::discrete::graphs::tests) | done |
| `chromatic_number_exact` | `rules::discrete::graphs` operator `graph_chromatic_number` (tests rules::discrete::graphs::tests) | done |

## `src/symbolic/graph_operations.rs` (8)

| legacy function | new home | status |
|---|---|---|
| `induced_subgraph` | `rules::discrete::graphs` operator `induced_subgraph` (tests rules::discrete::graphs::tests) | done |
| `union` | `rules::discrete::graphs` operator `graph_union` (tests rules::discrete::graphs::tests) | done |
| `intersection` | `rules::discrete::graphs` operator `graph_intersection` (tests rules::discrete::graphs::tests) | done |
| `cartesian_product` | `rules::discrete::graphs` operator `graph_cartesian` (tests rules::discrete::graphs::tests) | done |
| `tensor_product` | `rules::discrete::graphs` operator `graph_tensor` (tests rules::discrete::graphs::tests) | done |
| `complement` | `rules::discrete::graphs` operator `graph_complement` (tests rules::discrete::graphs::tests) | done |
| `disjoint_union` | `rules::discrete::graphs` operator `graph_disjoint_union` (tests rules::discrete::graphs::tests) | done |
| `join` | `rules::discrete::graphs` operator `graph_join` (tests rules::discrete::graphs::tests) | done |

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
| `new` | `rules::discrete::groups` operator `group_from_table` (tests rules::discrete::groups::tests) | done |
| `multiply` | `rules::discrete::groups` operator `group_mul` (tests rules::discrete::groups::tests) | done |
| `inverse` | `rules::discrete::groups` operator `group_inverse` (tests rules::discrete::groups::tests) | done |
| `is_abelian` | `rules::discrete::groups` operator `group_is_abelian` (tests rules::discrete::groups::tests) | done |
| `element_order` | `rules::discrete::groups` operator `group_element_order` (tests rules::discrete::groups::tests) | done |
| `conjugacy_classes` | `rules::discrete::groups` operator `group_conjugacy_classes` (tests rules::discrete::groups::tests) | done |
| `center` | `rules::discrete::groups` operator `group_center` (tests rules::discrete::groups::tests) | done |
| `is_valid` | `rules::discrete::groups` operator `group_is_valid` (tests rules::discrete::groups::tests) | done |
| `character` | `rules::discrete::groups` operator `group_character` (tests rules::discrete::groups::tests) | done |

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
| `new` | `rules::functional` operator `fredholm_solve` (tests rules::functional::tests) (equation parameters are the arguments `f, lambda, K, x, t, a, b`) | done |
| `solve_neumann_series` | `rules::functional` operator `fredholm_neumann` (tests rules::functional::tests) | done |
| `solve_separable_kernel` | `rules::functional` operator `fredholm_separable` (tests rules::functional::tests) | done |
| `solve_successive_approximations` | `rules::functional` operator `volterra_successive` (tests rules::functional::tests) | done |
| `solve_by_differentiation` | `rules::functional` operator `volterra_to_ode` (tests rules::functional::tests) (with `volterra_solve`) | done |
| `solve_airfoil_equation` | `rules::functional` operator `airfoil_equation` | partial (no test: `airfoil_equation` is defined but no test exercises it) |
| `solve_airfoil_equation_internal` | `rules::functional` operator `airfoil_equation` | partial (no test: `airfoil_equation` is defined but no test exercises it) |

## `src/symbolic/integration.rs` (8)

| legacy function | new home | status |
|---|---|---|
| `integrate_rational_function` | `rules::calculus` `integral` partial-fraction stage, exact over Q (test rational_functions) | done |
| `risch_norman_integrate` | `rules::calculus` `integral` staged heuristics (table, partial fractions, substitution, by parts, verified by differentiation); no Risch-Norman undetermined-coefficients ansatz | partial |
| `integrate_poly_exp` | `integral` by_parts stage integrates polynomial*exp(a x) (test by_parts) and `erfi`/`erf` for Gaussians; no general exp-extension (towers `exp(g(x))` with polynomial-in-t coefficients, `g` non-linear) integrator | partial |
| `poly_from_coeffs` | `poly::repr::Poly::from_univariate` | done |
| `partial_fraction_integrate` | `rules::calculus` `integral` partial-fraction stage, exact over Q (test rational_functions) | done |
| `hermite_integrate_rational` | `rules::calculus` `integral` partial-fraction stage, exact over Q (test rational_functions) | done |
| `integrate_rational_function_expr` | `rules::calculus` `integral` partial-fraction stage, exact over Q (test rational_functions) | done |
| `poly_derivative_symbolic` | `poly::repr::Poly::derivative` / `diff` (test derivatives_symbolic) | done |

## `src/symbolic/lie_groups_and_algebras.rs` (11)

| legacy function | new home | status |
|---|---|---|
| `lie_bracket` | `rules::lie` operator `lie_bracket` (tests rules::lie::tests) | done |
| `lie_bracket_internal` | `rules::lie` operator `lie_bracket` (tests rules::lie::tests) | done |
| `exponential_map` | `rules::lie` operator `exp_map` (tests rules::lie::tests) | done |
| `adjoint_representation_group` | `rules::lie` operator `adjoint_group` (tests rules::lie::tests) | done |
| `adjoint_representation_algebra` | `rules::lie` operator `adjoint_rep` (tests rules::lie::tests) | done |
| `commutator_table` | `rules::lie` operator `commutator_table` (tests rules::lie::tests) | done |
| `check_jacobi_identity` | `rules::lie` operator `check_jacobi` (tests rules::lie::tests) | done |
| `so3_generators` | `rules::lie` operator `so3_basis` (tests rules::lie::tests) | done |
| `so3` | `rules::lie` operator `so3_basis` (tests rules::lie::tests) (the algebra is `so3_basis()` with `structure_constants`/`killing_form`) | done |
| `su2_generators` | `rules::lie` operator `su2_basis` (tests rules::lie::tests) with `-i sigma/2` convention (legacy `+i sigma/2` flipped the sign of the structure constants) | done |
| `su2` | `rules::lie` operator `su2_basis` (tests rules::lie::tests) (the algebra is `su2_basis()` with `structure_constants`/`killing_form`) | done |

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
| `svd_decomposition` | `rules::linalg` `svd` (numeric, faer); no exact symbolic SVD | partial |
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
| `evaluate_complex` | `rules::complex::eval_complex` / `Graph::eval_complex` (complex evaluation of terms; test rules::complex::tests) | done |

## `src/symbolic/ode.rs` (22)

| legacy function | new home | status |
|---|---|---|
| `solve_ode` | `rules::ode` `dsolve` (tests first_order_linear_and_separable ... cauchy_euler) | done |
| `solve_ode_internal` | `rules::ode` `dsolve` (tests first_order_linear_and_separable ... cauchy_euler) | done |
| `solve_ode_system` | `rules::ode` `odeint` (numeric systems only); `dsolve` takes a single equation: no symbolic solver for linear/first-order systems | partial |
| `solve_ode_system_internal` | `rules::ode` `odeint` (numeric systems only); `dsolve` takes a single equation: no symbolic solver for linear/first-order systems | partial |
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
| `solve_ode_by_fourier` | needs the Fourier derivative theorem (transforms branch) to turn the ODE into an algebraic equation | pending |
| `solve_ode_by_fourier_internal` | needs the Fourier derivative theorem (transforms branch) to turn the ODE into an algebraic equation | pending |

## `src/symbolic/optimize.rs` (4)

| legacy function | new home | status |
|---|---|---|
| `find_extrema` | `rules::optimize` operator `find_extrema(f, vars)` (test rules::optimize::tests) | done |
| `hessian_matrix` | `rules::linalg` `hessian` (test vector_calculus) | done |
| `hessian_matrix_internal` | `rules::linalg` `hessian` (test vector_calculus) | done |
| `find_constrained_extrema` | `rules::optimize` operator `find_constrained_extrema(f, constraints, vars)` (Lagrange multipliers; test rules::optimize::tests) | done |

## `src/symbolic/pde.rs` (37)

| legacy function | new home | status |
|---|---|---|
| `solve_pde` | `rules::pde` operator `solve_pde` (exercised through the method dispatch of `pdsolve`, test rules::pde::tests) | done |
| `solve_pde_internal` | `rules::pde` operator `solve_pde` (exercised through the method dispatch of `pdsolve`, test rules::pde::tests) | done |
| `solve_pde_by_separation_of_variables` | `rules::pde` operator `solve_pde_by_separation_of_variables` (exercised through the method dispatch of `pdsolve`, test rules::pde::tests) | done |
| `solve_pde_by_separation_of_variables_internal` | `rules::pde` operator `solve_pde_by_separation_of_variables` (exercised through the method dispatch of `pdsolve`, test rules::pde::tests) | done |
| `classify_pde_heuristic` | `rules::pde` operator `pde_classify` (test rules::pde::tests) | done |
| `solve_pde_by_characteristics` | `rules::pde` operator `solve_pde_by_characteristics` (test rules::pde::tests) | done |
| `solve_pde_by_characteristics_internal` | `rules::pde` operator `solve_pde_by_characteristics` (test rules::pde::tests) | done |
| `solve_pde_by_greens_function` | `rules::pde` operator `solve_pde_by_greens_function` (exercised through the method dispatch of `pdsolve`, test rules::pde::tests) | done |
| `solve_pde_by_greens_function_internal` | `rules::pde` operator `solve_pde_by_greens_function` (exercised through the method dispatch of `pdsolve`, test rules::pde::tests) | done |
| `solve_second_order_pde` | `rules::pde` operator `solve_second_order_pde` (exercised through the method dispatch of `pdsolve`, test rules::pde::tests) | done |
| `solve_second_order_pde_internal` | `rules::pde` operator `solve_second_order_pde` (exercised through the method dispatch of `pdsolve`, test rules::pde::tests) | done |
| `solve_wave_equation_1d_dalembert` | `rules::pde` operator `solve_wave_equation_1d_dalembert` (test rules::pde::tests) | done |
| `solve_wave_equation_1d_dalembert_internal` | `rules::pde` operator `solve_wave_equation_1d_dalembert` (test rules::pde::tests) | done |
| `solve_heat_equation_1d` | `rules::pde` operator `solve_heat_equation_1d` (test rules::pde::tests) | done |
| `solve_heat_equation_1d_internal` | `rules::pde` operator `solve_heat_equation_1d` (test rules::pde::tests) | done |
| `solve_laplace_equation_2d` | `rules::pde` operator `solve_laplace_equation_2d` (exercised through the method dispatch of `pdsolve`, test rules::pde::tests) | done |
| `solve_laplace_equation_2d_internal` | `rules::pde` operator `solve_laplace_equation_2d` (exercised through the method dispatch of `pdsolve`, test rules::pde::tests) | done |
| `solve_wave_equation_3d` | `rules::pde` operator `solve_wave_equation_3d` (test rules::pde::tests) | done |
| `solve_wave_equation_3d_internal` | `rules::pde` operator `solve_wave_equation_3d` (test rules::pde::tests) | done |
| `solve_heat_equation_3d` | `rules::pde` operator `solve_heat_equation_3d` (test rules::pde::tests) | done |
| `solve_heat_equation_3d_internal` | `rules::pde` operator `solve_heat_equation_3d` (test rules::pde::tests) | done |
| `solve_laplace_equation_3d` | `rules::pde` operator `solve_laplace_equation_3d` (exercised through the method dispatch of `pdsolve`, test rules::pde::tests) | done |
| `solve_laplace_equation_3d_internal` | `rules::pde` operator `solve_laplace_equation_3d` (exercised through the method dispatch of `pdsolve`, test rules::pde::tests) | done |
| `solve_poisson_equation_2d` | `rules::pde` operator `solve_poisson_equation_2d` (test rules::pde::tests) | done |
| `solve_poisson_equation_2d_internal` | `rules::pde` operator `solve_poisson_equation_2d` (test rules::pde::tests) | done |
| `solve_poisson_equation_3d` | `rules::pde` operator `solve_poisson_equation_3d` (test rules::pde::tests) | done |
| `solve_poisson_equation_3d_internal` | `rules::pde` operator `solve_poisson_equation_3d` (test rules::pde::tests) | done |
| `solve_helmholtz_equation` | `rules::pde` operator `solve_helmholtz_equation` (test rules::pde::tests) | done |
| `solve_helmholtz_equation_internal` | `rules::pde` operator `solve_helmholtz_equation` (test rules::pde::tests) | done |
| `solve_schrodinger_equation` | `rules::pde` operator `solve_schrodinger_equation` (test rules::pde::tests) | done |
| `solve_schrodinger_equation_internal` | `rules::pde` operator `solve_schrodinger_equation` (test rules::pde::tests) | done |
| `solve_klein_gordon_equation` | `rules::pde` operator `solve_klein_gordon_equation` (test rules::pde::tests) | done |
| `solve_klein_gordon_equation_internal` | `rules::pde` operator `solve_klein_gordon_equation` (test rules::pde::tests) | done |
| `solve_burgers_equation` | `rules::pde` operator `solve_burgers_equation` (test rules::pde::tests) | done |
| `solve_burgers_equation_internal` | `rules::pde` operator `solve_burgers_equation` (test rules::pde::tests) | done |
| `solve_with_fourier_transform` | `rules::pde` operator `solve_with_fourier_transform` (exercised through the method dispatch of `pdsolve`, test rules::pde::tests) | done |
| `solve_with_fourier_transform_internal` | `rules::pde` operator `solve_with_fourier_transform` (exercised through the method dispatch of `pdsolve`, test rules::pde::tests) | done |

## `src/symbolic/poly_factorization.rs` (11)

| legacy function | new home | status |
|---|---|---|
| `factor_gf` | no public GF(p) operator: only the private helper inside `rules::poly::univariate` used by `factor` over Q | pending (in progress: transforms/complex/finite-field/units branch) |
| `poly_derivative_gf` | no public GF(p) operator: only the private helper inside `rules::poly::univariate` used by `factor` over Q | pending (in progress: transforms/complex/finite-field/units branch) |
| `square_free_factorization_gf` | no public GF(p) operator: only the private helper inside `rules::poly::univariate` used by `factor` over Q | pending (in progress: transforms/complex/finite-field/units branch) |
| `berlekamp_factorization` | no public GF(p) operator: only the private helper inside `rules::poly::univariate` used by `factor` over Q | pending (in progress: transforms/complex/finite-field/units branch) |
| `berlekamp_zassenhaus` | `rules::poly::univariate::factor` (Berlekamp-Zassenhaus over Q; tests known_factorisations, factorisation_reproduces_the_input) | done |
| `cantor_zassenhaus` | no public GF(p) operator: only the private helper inside `rules::poly::univariate` used by `factor` over Q | pending (in progress: transforms/complex/finite-field/units branch) |
| `distinct_degree_factorization` | no public GF(p) operator: only the private helper inside `rules::poly::univariate` used by `factor` over Q | pending (in progress: transforms/complex/finite-field/units branch) |
| `poly_gcd_gf` | `rules::discrete::finite_field` operator `gfp_gcd(f, g, p)` (tests rules::discrete::finite_field::tests) | done |
| `poly_pow_mod` | `rules::discrete::finite_field` operator `gfx_pow(a, e, m, p)` (tests rules::discrete::finite_field::tests) (power modulo the polynomial `m` over GF(p)) | done |
| `poly_mul_scalar` | `rules::discrete::finite_field` operator `gfp_mul(f, list(c), p)` (tests rules::discrete::finite_field::tests) (scalar as a constant polynomial) | done |
| `poly_extended_gcd` | `rules::discrete::finite_field` operator `gfp_egcd(f, g, p)` (tests rules::discrete::finite_field::tests) | done |

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
| `leading_term` | `rules::poly::algebra` operator `leading_term(p, vars[, order])` (tests rules::poly::algebra::tests) | done |
| `long_division` | `rules::poly` `quo`/`rem` and `univariate::divrem` (test division_and_gcd) | done |
| `get_coeffs_as_vec` | `poly::repr::Poly::univariate_in` (test univariate_views) | done |
| `get_coeff_for_power` | `rules::poly` `coeff` (test degree_and_coefficients) | done |
| `prune_zeros` | Expr helper; `free_of` guard / `poly::repr::from_term` / Poly keeps no zero terms | dropped |
| `poly_from_coeffs` | `poly::repr::Poly::from_univariate` / `to_term` (test expansion) | done |
| `sparse_poly_to_expr` | `poly::repr::to_term` (test expansion) | done |

## `src/symbolic/proof.rs` (7)

| legacy function | new home | status |
|---|---|---|
| `verify_equation_solution` | `rules::verify` operator `verify_solution` (tests rules::verify::tests) | done |
| `verify_indefinite_integral` | `rules::verify` operator `verify_integral` (tests rules::verify::tests) | done |
| `verify_definite_integral` | `rules::verify` operator `verify_definite_integral` (tests rules::verify::tests) | done |
| `verify_ode_solution` | `rules::verify` operator `verify_ode_solution` (tests rules::verify::tests) | done |
| `verify_matrix_inverse` | `rules::verify` operator `verify_inverse` (tests rules::verify::tests) | done |
| `verify_derivative` | `rules::verify` operator `verify_derivative` (tests rules::verify::tests) | done |
| `verify_limit` | `rules::verify` operator `verify_limit` (tests rules::verify::tests) | done |

## `src/symbolic/quantum_field_theory.rs` (8)

| legacy function | new home | status |
|---|---|---|
| `dirac_adjoint` | `rules::physics` operator `dirac_adjoint` (test rules::physics::tests) | done |
| `feynman_slash` | `rules::physics` operator `feynman_slash` (test rules::physics::tests) | done |
| `scalar_field_lagrangian` | `rules::physics` operator `scalar_field_lagrangian` (test rules::physics::tests) | done |
| `qed_lagrangian` | `rules::physics` operator `qed_lagrangian` | partial (no test: `qed_lagrangian` is defined but no test exercises it) |
| `qcd_lagrangian` | `rules::physics` operator `qcd_lagrangian` | partial (no test: `qcd_lagrangian` is defined but no test exercises it) |
| `propagator` | `rules::physics` operator `propagator` | partial (no test: `propagator` is defined but no test exercises it) |
| `scattering_cross_section` | `rules::physics` operator `scattering_cross_section` | partial (no test: `scattering_cross_section` is defined but no test exercises it) |
| `feynman_propagator_position_space` | `rules::physics` operator `feynman_propagator_position_space` | partial (no test: `feynman_propagator_position_space` is defined but no test exercises it) |

## `src/symbolic/quantum_mechanics.rs` (21)

| legacy function | new home | status |
|---|---|---|
| `bra_ket` | `rules::physics` operator `braket` (test rules::physics::tests) | done |
| `bra_ket_internal` | `rules::physics` operator `braket_on` | partial (no test: `braket_on` is defined but no test exercises it) |
| `new` | `rules::physics` operator `op_mul` (with `op_d`, `op_add`, `op_compose`, ... from rules::functional) (test rules::physics::tests) | done |
| `apply` | `rules::physics` operator `qm_apply` (test rules::physics::tests) | done |
| `commutator` | `rules::physics` operator `commutator` (test rules::physics::tests) | done |
| `commutator_internal` | `rules::physics` operator `commutator` (test rules::physics::tests) | done |
| `expectation_value` | `rules::physics` operator `expectation_value` (test rules::physics::tests) | done |
| `expectation_value_internal` | `rules::physics` operator `expectation_value` (test rules::physics::tests) | done |
| `uncertainty` | `rules::physics` operator `uncertainty` (test rules::physics::tests) | done |
| `probability_density` | `rules::physics` operator `probability_density` (test rules::physics::tests) | done |
| `hamiltonian_free_particle` | `rules::physics` operator `hamiltonian_free_particle` | partial (no test: `hamiltonian_free_particle` is defined but no test exercises it) |
| `hamiltonian_harmonic_oscillator` | `rules::physics` operator `hamiltonian_harmonic_oscillator` (test rules::physics::tests) | done |
| `angular_momentum_z` | `rules::physics` operator `angular_momentum_z` | partial (no test: `angular_momentum_z` is defined but no test exercises it) |
| `pauli_matrices` | `rules::physics` operator `pauli_matrices` | partial (no test: `pauli_matrices` is defined but no test exercises it) |
| `spin_operator` | `rules::physics` operator `spin_operator` (test rules::physics::tests) | done |
| `solve_time_independent_schrodinger` | `rules::physics` operator `solve_time_independent_schrodinger` | partial (no test: `solve_time_independent_schrodinger` is defined but no test exercises it) |
| `time_dependent_schrodinger_equation` | `rules::physics` operator `time_dependent_schrodinger_equation` | partial (no test: `time_dependent_schrodinger_equation` is defined but no test exercises it) |
| `dirac_equation` | `rules::physics` operator `dirac_equation` (test rules::physics::tests) | done |
| `klein_gordon_equation` | `rules::physics` operator `klein_gordon_equation` (test rules::physics::tests) | done |
| `first_order_energy_correction` | `rules::physics` operator `first_order_energy_correction` | partial (no test: `first_order_energy_correction` is defined but no test exercises it) |
| `scattering_amplitude` | `rules::physics` operator `scattering_amplitude` | partial (no test: `scattering_amplitude` is defined but no test exercises it) |

## `src/symbolic/radicals.rs` (2)

| legacy function | new home | status |
|---|---|---|
| `simplify_radicals` |  | pending |
| `denest_sqrt` |  | pending |

## `src/symbolic/real_roots.rs` (4)

| legacy function | new home | status |
|---|---|---|
| `sturm_sequence` | `kernels::real_roots::sturm_sequence` (tests/kernels/real_roots.rs) | done |
| `count_real_roots_in_interval` | `rules::poly::algebra` operator `count_real_roots(p, x, a, b)` (tests rules::poly::algebra::tests) (Sturm sequences) | done |
| `isolate_real_roots` | `kernels::real_roots::isolate_real_roots` (tests/kernels/real_roots.rs) | done |
| `eval_expr` | Expr helper; `Term::eval` | dropped |

## `src/symbolic/relativity.rs` (13)

| legacy function | new home | status |
|---|---|---|
| `lorentz_factor` | `rules::physics` operator `lorentz_factor` (test rules::physics::tests) | done |
| `lorentz_transformation_x` | `rules::physics` operator `lorentz_transformation_x` (test rules::physics::tests) | done |
| `velocity_addition` | `rules::physics` operator `velocity_addition` (test rules::physics::tests) | done |
| `mass_energy_equivalence` | `rules::physics` operator `mass_energy_equivalence` (test rules::physics::tests) | done |
| `relativistic_momentum` | `rules::physics` operator `relativistic_momentum` (test rules::physics::tests) | done |
| `doppler_effect` | `rules::physics` operator `doppler_effect` | partial (no test: `doppler_effect` is defined but no test exercises it) |
| `schwarzschild_radius` | `rules::physics` operator `schwarzschild_radius` (test rules::physics::tests) | done |
| `gravitational_time_dilation` | `rules::physics` operator `gravitational_time_dilation` | partial (no test: `gravitational_time_dilation` is defined but no test exercises it) |
| `einstein_tensor` | `rules::physics` operator `einstein_tensor_from` | partial (no test: `einstein_tensor_from` is defined but no test exercises it) |
| `geodesic_acceleration` | `rules::physics` operator `geodesic_acceleration` (test rules::physics::tests) | done |
| `lorentz_transformation` | legacy alias of lorentz_transformation_x | dropped |
| `einstein_field_equations` | legacy placeholder, superseded by einstein_tensor | dropped |
| `geodesic_equation` | legacy placeholder, superseded by geodesic_acceleration | dropped |

## `src/symbolic/rewriting.rs` (2)

| legacy function | new home | status |
|---|---|---|
| `apply_rules_to_normal_form` | `rules::rewriting` operator `rewrite_with(t, rules, vars)` (tests rules::rewriting::tests) | done |
| `knuth_bendix` | `rules::rewriting` operator `knuth_bendix(equations, variables[, precedence])` (tests rules::rewriting::tests) | done |

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
| `analyze_convergence` | `rules::calculus` operator `converges(term, k)` (divergence, ratio, root, alternating, limit-comparison and Cauchy-condensation tests; tests in rules::calculus::tests) | done |
| `asymptotic_expansion` | `rules::calculus` definition `asymptotic(f, x, n) := laurent(f, x, oo, n)` (expansions at infinity; tests in rules::calculus::tests) | done |
| `analytic_continuation` | `rules::complex::analysis` operator `continue_along(f, z, points, order)` (test rules::complex::analysis::tests) | done |

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
| `new` | `rules::physics` operator `lattice_volume` (test rules::physics::tests) | done |
| `volume` | `rules::physics` operator `lattice_volume` (test rules::physics::tests) | done |
| `reciprocal_lattice_vectors` | `rules::physics` operator `reciprocal_lattice_vectors` (test rules::physics::tests) | done |
| `bloch_theorem` | `rules::physics` operator `bloch_wave` | partial (no test: `bloch_wave` is defined but no test exercises it) |
| `energy_band` | `rules::physics` operator `energy_band` (test rules::physics::tests) | done |
| `density_of_states_3d` | `rules::physics` operator `density_of_states_3d` | partial (no test: `density_of_states_3d` is defined but no test exercises it) |
| `fermi_energy_3d` | `rules::physics` operator `fermi_energy_3d` | partial (no test: `fermi_energy_3d` is defined but no test exercises it) |
| `drude_conductivity` | `rules::physics` operator `drude_conductivity` | partial (no test: `drude_conductivity` is defined but no test exercises it) |
| `hall_coefficient` | `rules::physics` operator `hall_coefficient` (test rules::physics::tests) | done |
| `debye_frequency` | `rules::physics` operator `debye_frequency` | partial (no test: `debye_frequency` is defined but no test exercises it) |
| `einstein_heat_capacity` | `rules::physics` operator `einstein_heat_capacity` (test rules::physics::tests) | done |
| `plasma_frequency` | `rules::physics` operator `plasma_frequency` | partial (no test: `plasma_frequency` is defined but no test exercises it) |
| `london_penetration_depth` | `rules::physics` operator `london_penetration_depth` | partial (no test: `london_penetration_depth` is defined but no test exercises it) |

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
| `inverse_erfc` | `rules::special` operator `erfcinv`, definition `inverse_erfc(x)` (tests rules::special::tests) | done |
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
| `bessel_k0` | `rules::special` definition `bessel_k0(x) := besselk(0, x)` (tests rules::special::tests) | done |
| `bessel_k1` | `rules::special` definition `bessel_k1(x) := besselk(1, x)` (tests rules::special::tests) | done |
| `sinc` | `kernels::special::sinc (test_sinc)` | done |
| `zeta` | `kernels::special::riemann_zeta (test_riemann_zeta)` | done |
| `ln_factorial` | `rules::special` definition `ln_factorial(n) := lgamma(n + 1)` (tests rules::special::tests) | done |
| `regularized_gamma_p` | `kernels::special::regularized_lower_gamma (incomplete_gamma_closed_forms)` | done |
| `regularized_gamma_q` | `kernels::special::regularized_upper_gamma (incomplete_gamma_closed_forms)` | done |

## `src/symbolic/special_functions.rs` (26)

| legacy function | new home | status |
|---|---|---|
| `gamma` | `rules::special` `gamma` (test gamma_family_values) | done |
| `ln_gamma` | `rules::special` `lgamma` (test gamma_family_values) | done |
| `beta` | `rules::special` `beta` (test gamma_family_values) | done |
| `digamma` | `rules::special` `digamma` (test gamma_family_values) | done |
| `polygamma` | `rules::special` operator `polygamma(n, x)` (exact at integers, `digamma` for n = 0, derivative rules) (tests rules::special::tests) | done |
| `erf` | `rules::special` `erf`/`erfc` (test error_function_values) | done |
| `erfc` | `rules::special` `erf`/`erfc` (test error_function_values) | done |
| `erfi` | `rules::special` operator `erfi` (integral of exp(x^2) reduces to it) (tests rules::special::tests) | done |
| `zeta` | `rules::special` `zeta` (test zeta_values) | done |
| `bessel_j` | `rules::special` operator `besselj(n, x)` (any order via kernels::special::bessel_j; half-integer closed forms) (tests rules::special::tests) | done |
| `bessel_y` | `rules::special` operator `bessely(n, x)` (any order via kernels::special::bessel_y; half-integer closed forms) (tests rules::special::tests) | done |
| `bessel_i` | `rules::special` operator `besseli(n, x)` (any order via kernels::special::bessel_i; half-integer closed forms) (tests rules::special::tests) | done |
| `bessel_k` | `rules::special` operator `besselk(n, x)` (any order via kernels::special::bessel_k; half-integer closed forms) (tests rules::special::tests) | done |
| `legendre_p` | `rules::special` `legendre`/`laguerre`/`hermite`/`chebyshevt`/`chebyshevu` (tests polynomials_of_any_degree_match_the_recurrence) | done |
| `laguerre_l` | `rules::special` `legendre`/`laguerre`/`hermite`/`chebyshevt`/`chebyshevu` (tests polynomials_of_any_degree_match_the_recurrence) | done |
| `generalized_laguerre` | `rules::special` operator `laguerre_gen(n, alpha, x)` (tests rules::special::tests) | done |
| `hermite_h` | `rules::special` `legendre`/`laguerre`/`hermite`/`chebyshevt`/`chebyshevu` (tests polynomials_of_any_degree_match_the_recurrence) | done |
| `chebyshev_t` | `rules::special` `legendre`/`laguerre`/`hermite`/`chebyshevt`/`chebyshevu` (tests polynomials_of_any_degree_match_the_recurrence) | done |
| `chebyshev_u` | `rules::special` `legendre`/`laguerre`/`hermite`/`chebyshevt`/`chebyshevu` (tests polynomials_of_any_degree_match_the_recurrence) | done |
| `bessel_differential_equation` | `rules::special` definition `bessel_differential_equation(y, x, n)` (tests rules::special::tests) | done |
| `legendre_differential_equation` | `rules::special` definition `legendre_differential_equation(y, x, n)` | partial (no test: `legendre_differential_equation` is defined but no test exercises it) |
| `legendre_rodrigues_formula` | `rules::special` definition `legendre_rodrigues(n, x)` (tests rules::special::tests) | done |
| `laguerre_differential_equation` | `rules::special` definition `laguerre_differential_equation(y, x, n)` | partial (no test: `laguerre_differential_equation` is defined but no test exercises it) |
| `hermite_differential_equation` | `rules::special` definition `hermite_differential_equation(y, x, n)` | partial (no test: `hermite_differential_equation` is defined but no test exercises it) |
| `hermite_rodrigues_formula` | `rules::special` definition `hermite_rodrigues(n, x)` (tests rules::special::tests) | done |
| `chebyshev_differential_equation` | `rules::special` definition `chebyshev_differential_equation(y, x, n)` | partial (no test: `chebyshev_differential_equation` is defined but no test exercises it) |

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
| `polynomial_regression_symbolic` | `polynomial_regression` (rules::stats): exact and float data only; symbolic (non-numeric) data points are not supported | partial |

## `src/symbolic/tensor.rs` (15)

| legacy function | new home | status |
|---|---|---|
| `new` | a tensor is a nested `list` term (see `rules::geometry` `tensor_rank`, `component`) | done |
| `rank` | `rules::geometry` operator `tensor_rank` (tests rules::geometry::tests) | done |
| `get` | `rules::geometry` operator `component` (tests rules::geometry::tests) | done |
| `add` | `rules::geometry` operator `tensor_add` (tests rules::geometry::tests) | done |
| `sub` | `rules::geometry` `tensor_add(A, tensor_smul(-1, B))` (no dedicated operator) | done |
| `scalar_mul` | `rules::geometry` operator `tensor_smul` (tests rules::geometry::tests) | done |
| `outer_product` | `rules::geometry` operator `tensor_outer` (tests rules::geometry::tests) | done |
| `contract` | `rules::geometry` operator `contract` (tests rules::geometry::tests) | done |
| `to_matrix_expr` | tensors are nested lists, so a rank-2 tensor already is a matrix term for `rules::linalg` (`matmul`, `det`, ...) | done |
| `raise_index` | `rules::geometry` operator `raise_index` (tests rules::geometry::tests) | done |
| `lower_index` | `rules::geometry` operator `lower_index` (tests rules::geometry::tests) | done |
| `christoffel_symbols_first_kind` | `rules::geometry` operator `christoffel1` | partial (no test: `christoffel1` is defined but no test exercises it) |
| `christoffel_symbols_second_kind` | `rules::geometry` operator `christoffel` (tests rules::geometry::tests) | done |
| `riemann_curvature_tensor` | `rules::geometry` operator `riemann` (tests rules::geometry::tests) | done |
| `covariant_derivative_vector` | `rules::geometry` operator `covariant_derivative` (tests rules::geometry::tests) | done |

## `src/symbolic/thermodynamics.rs` (13)

| legacy function | new home | status |
|---|---|---|
| `first_law_thermodynamics` | `rules::physics` operator `first_law_thermodynamics` | partial (no test: `first_law_thermodynamics` is defined but no test exercises it) |
| `ideal_gas_law` | `rules::physics` operator `ideal_gas_law` | partial (no test: `ideal_gas_law` is defined but no test exercises it) |
| `enthalpy` | `rules::physics` operator `enthalpy` (test rules::physics::tests) | done |
| `helmholtz_free_energy` | `rules::physics` operator `helmholtz_free_energy` | partial (no test: `helmholtz_free_energy` is defined but no test exercises it) |
| `gibbs_free_energy` | `rules::physics` operator `gibbs_free_energy` (test rules::physics::tests) | done |
| `boltzmann_entropy` | `rules::physics` operator `boltzmann_entropy` | partial (no test: `boltzmann_entropy` is defined but no test exercises it) |
| `carnot_efficiency` | `rules::physics` operator `carnot_efficiency` (test rules::physics::tests) | done |
| `boltzmann_distribution` | `rules::physics` operator `boltzmann_distribution` | partial (no test: `boltzmann_distribution` is defined but no test exercises it) |
| `partition_function` | `rules::physics` operator `partition_function` (test rules::physics::tests) | done |
| `fermi_dirac_distribution` | `rules::physics` operator `fermi_dirac_distribution` (test rules::physics::tests) | done |
| `bose_einstein_distribution` | `rules::physics` operator `bose_einstein_distribution` | partial (no test: `bose_einstein_distribution` is defined but no test exercises it) |
| `work_isothermal_expansion` | `rules::physics` operator `work_isothermal_expansion` | partial (no test: `work_isothermal_expansion` is defined but no test exercises it) |
| `verify_maxwell_relation_helmholtz` | `rules::physics` operator `verify_maxwell_relation_helmholtz` (test rules::physics::tests) | done |

## `src/symbolic/topology.rs` (19)

| legacy function | new home | status |
|---|---|---|
| `new` | `rules::discrete::topology` operator `sc_complex` (tests rules::discrete::topology::tests) | done |
| `dimension` | `rules::discrete::topology` operator `sc_dimension` (tests rules::discrete::topology::tests) | done |
| `boundary` | `rules::discrete::topology` operator `sc_boundary` (tests rules::discrete::topology::tests) | done |
| `symbolic_boundary` | `rules::discrete::topology` operator `sc_boundary`; chain coefficients are arbitrary terms (tests rules::discrete::topology::tests) | done |
| `add_term` | `rules::discrete::topology` operator `sc_chain_boundary`; a chain is `list(list(coeff, simplex), ...)` and equal simplices are combined (tests rules::discrete::topology::tests) | done |
| `add_simplex` | `rules::discrete::topology` operator `sc_complex`; `sc_complex(list(...))` closes under faces (tests rules::discrete::topology::tests) | done |
| `get_simplices_by_dim` | `rules::discrete::topology` operator `sc_simplices` (tests rules::discrete::topology::tests) | done |
| `get_boundary_matrix` | `rules::discrete::topology` operator `sc_boundary_matrix` (tests rules::discrete::topology::tests) | done |
| `get_symbolic_boundary_matrix` | `rules::discrete::topology` operator `sc_boundary_matrix` (tests rules::discrete::topology::tests) | done |
| `apply_boundary_operator` | `rules::discrete::topology` operator `sc_chain_boundary` (tests rules::discrete::topology::tests) | done |
| `apply_symbolic_boundary_operator` | `rules::discrete::topology` operator `sc_chain_boundary`; coefficients may be symbolic (tests rules::discrete::topology::tests) | done |
| `compute_euler_characteristic` | `rules::discrete::topology` operator `sc_euler_characteristic` (tests rules::discrete::topology::tests) | done |
| `verify_boundary_property` | `rules::discrete::topology` operator `sc_verify_boundary` (tests rules::discrete::topology::tests) | done |
| `verify_coboundary_property` | `rules::discrete::topology` operator `sc_verify_coboundary` (tests rules::discrete::topology::tests) | done |
| `compute_homology_betti_number` | `rules::discrete::topology` operator `sc_betti` (tests rules::discrete::topology::tests) | done |
| `compute_cohomology_betti_number` | `rules::discrete::topology` operator `sc_cohomology_betti` (tests rules::discrete::topology::tests) | done |
| `create_grid_complex` | `rules::discrete::topology` operator `sc_grid` (tests rules::discrete::topology::tests) | done |
| `create_torus_complex` | `rules::discrete::topology` operator `sc_torus` (tests rules::discrete::topology::tests) | done |
| `vietoris_rips_filtration` | `rules::discrete::topology` operator `vietoris_rips_filtration` (tests rules::discrete::topology::tests) | done |

## `src/symbolic/transforms.rs` (28)

| legacy function | new home | status |
|---|---|---|
| `fourier_time_shift` | `rules::transforms` `fourier` (shift/modulation theorems applied inside the kernel) | done |
| `fourier_frequency_shift` | `rules::transforms` `fourier` modulation by exp(I a t), cos, sin | done |
| `fourier_scaling` | `rules::transforms` `fourier` (linear arguments through the table); no general scaling theorem `F(f(a t)) = F(w/a)/|a|` | partial (in progress: transforms/complex/finite-field/units branch) |
| `fourier_differentiation` | `rules::transforms` `fourier` (multiplication-by-t theorem); derivative theorem `F(f-prime) = I w F(f)` not applied | partial (in progress: transforms/complex/finite-field/units branch) |
| `laplace_time_shift` | `rules::transforms` `laplace` heaviside(t - c) shift theorem | done |
| `laplace_differentiation` | `rules::transforms` `laplace` of `diff(y(t), t)` = s L[y] - y(0) (first order) | done |
| `laplace_frequency_shift` | `rules::transforms` `laplace` exp(a t) shift theorem | done |
| `laplace_scaling` | `rules::transforms` `laplace` (linear arguments in the table) | done |
| `laplace_integration` | `rules::transforms` `laplace` of `defint(g, u, 0, t)` = G/s | done |
| `z_time_shift` | `rules::transforms` `ztransform` of `kronecker(n - k)`; general shifted sequences `f(n - k)` not detected | partial (in progress: transforms/complex/finite-field/units branch) |
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
| `convolution_fourier` | `rules::transforms` `convolve` builds the convolution integral; the transform-product theorem is not applied | partial (in progress: transforms/complex/finite-field/units branch) |
| `convolution_laplace` | `rules::transforms` `convolve(f, g, t)` = ∫_0^t f(u) g(t-u) du (test convolution) | done |

## `src/symbolic/unit_unification.rs` (1)

| legacy function | new home | status |
|---|---|---|
| `unify_expression` |  | pending (in progress: transforms/complex/finite-field/units branch) |

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
