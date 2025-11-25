# PDE Solver Implementation Progress

## Completed Tasks

### 1. Fixed Expr::PartialEq Bug ✅

**Problem**: The `Expr::PartialEq` implementation was not properly comparing `Derivative` expressions. It only compared the inner expression but ignored the variable name, causing `d/dx(d/dx(u))` to incorrectly compare equal to `d/dt(d/dt(u))`.

**Root Cause**: The `get_children_internal()` method for `Derivative(a, var)` only returned `[a]`, omitting the variable name. The generic fallback in `PartialEq` then compared only the children, missing the variable distinction.

**Solution**: Added explicit pattern matching in `Expr::PartialEq` for:

- `Derivative(e, v)` - compares both expression and variable name
- `Solve(e, v)`, `ConvergenceAnalysis(e, v)`, `ForAll(v, e)`, `Exists(v, e)` - also fixed for consistency

**Impact**: PDE solvers now correctly extract coefficients, producing accurate solutions like `u = F(x - 2t) + G(x + 2t)` instead of `u = F(x) + G(x)`.

### 2. Implemented PDE Solver Infrastructure ✅

**Manual Parsing Approach**:

- Created `collect_terms(expr)` - recursively collects terms from Add/Sub/Neg expressions
- Created `extract_coefficient(term, var)` - extracts coefficient of a variable from a term, handling DAG unwrapping
- Replaced pattern matching with manual term collection and coefficient extraction

**Implemented Solvers**:

1. **Method of Characteristics** (`solve_pde_by_characteristics`):
   - Solves first-order linear PDEs: `a*u_x + b*u_y = c`
   - Returns: `u = u_p(x) + F(y - (b/a)*x)` where `u_p` is a particular solution

2. **D'Alembert's Formula** (`solve_wave_equation_1d_dalembert`):
   - Solves 1D wave equation: `u_tt = c^2*u_xx`
   - Returns: `u = F(x - ct) + G(x + ct)`

### 3. Created Comprehensive Tests ✅

**Test File**: `tests/symbolic_pde_test.rs`

**Test Cases**:

1. `test_method_of_characteristics`: Tests `u_x + u_y = 1`
   - Expected: `u = x + F(y - x)`
   - Status: ✅ Passing

2. `test_wave_equation_dalembert`: Tests `u_tt - 4*u_xx = 0` (c = 2)
   - Expected: `u = F(x - 2t) + G(x + 2t)`
   - Status: ✅ Passing

## Next Steps

### 1. Add More PDE Types 🔄

- Heat equation (parabolic): `u_t = α*u_xx`
- Laplace equation (elliptic): `u_xx + u_yy = 0`
- Poisson equation: `u_xx + u_yy = f(x,y)`
- 2D wave equation: `u_tt = c^2*(u_xx + u_yy)`
- Schrödinger equation: `i*u_t = -u_xx + V(x)*u`

### 2. Expand Test Coverage 📝

- Add tests for each new PDE type
- Test with various boundary conditions
- Test with non-constant coefficients
- Property-based tests for solution verification

### 3. Implement FFI APIs 🔌

- Create `src/ffi_apis/symbolic_pde_ffi/` module
- Implement handle, JSON, and Bincode APIs
- Follow the pattern from `symbolic_ode_ffi` and `symbolic_calculus_ffi`
- Register in `src/ffi_apis/mod.rs`

## Technical Notes

### DAG Handling

All PDE solvers now properly handle DAG-wrapped expressions by:

1. Unwrapping DAG nodes in `solve_pde` before dispatching
2. Unwrapping DAG nodes in `collect_terms` during recursion
3. Unwrapping DAG nodes in `extract_coefficient` before comparison

### Coefficient Extraction

The `extract_coefficient` function handles:

- Direct matches: `u_x` → coefficient `1`
- Multiplication: `4*u_x` or `u_x*4` → coefficient `4`
- Negation: `-u_x` → coefficient `-1`
- Negated multiplication: `-(4*u_x)` → coefficient `-4`
- DAG-wrapped expressions at any level

### Known Limitations

1. Only handles linear PDEs with constant coefficients
2. Separation of variables not yet implemented
3. Boundary conditions not yet integrated into solutions
4. No numerical PDE solvers yet
