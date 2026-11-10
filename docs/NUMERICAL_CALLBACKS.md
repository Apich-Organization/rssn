# Numerical host callbacks

The `numerical_callback` rule set registers serial, effectful numerical
operations through the same manifest and dynamic-call protocol as other RSSN
operators. Each trailing function argument is a `CallbackRef` registered in the
call's `CallCtx`; the bank snapshots its concrete signature on registration.
The name-based and pinned `OpKey` entry points use the same evaluator.

| Operator | Ordinary arguments | Callback | Result fields |
| --- | --- | --- | --- |
| `num_integrate` | `a`, `b`, absolute tolerance, relative tolerance, panel budget | `Float[D] -> Float[E]` | `value: Float[D*E]`, `error: Float[D*E]`, `evaluations`, `converged` |
| `num_brent` | `a`, `b`, positive argument tolerance, positive iteration budget | `Float[D] -> Float[E]` | `root: Float[D]`, `evaluations` |
| `num_least_squares` | initial parameter vector, residual tolerance, gradient tolerance, step tolerance, initial damping, iteration budget, evaluation budget | `Tensor<1, Float[D]> -> Tensor<1, Float[E]>` | `x: Tensor<1, Float[D]>`, `residual_norm: Float[E]`, `iterations`, `evaluations`, `jacobian_evaluations`, `accepted_steps`, `rejected_steps`, `converged`, `status` |
| `num_ode_dopri5` | initial time, final time, state vector, absolute state tolerance, relative tolerance, maximal time step, attempt budget | `(Float[D], Tensor<1, Float[E]>) -> Tensor<1, Float[E/D]>` | `times: Tensor<1, Float[D]>`, `states: Tensor<2, Float[E]>`, `evaluations`, `rejected` |
| `num_ode_rosenbrock23` | same IVP arguments | same vector RHS | same IVP result |
| `num_ode_radau5` | same IVP arguments | same vector RHS | same IVP result |
| `num_ode_bdf` | same IVP arguments | same vector RHS | same IVP result |

The table's numeric coefficients are finite floating-point values. Linear
quantities are converted to coherent coefficients before entering a solver;
affine points are rejected. Domain, state and derivative dimensions are checked
against the actual callback signature, independently of the result type a host
later requests. An exact zero absolute tolerance contributes no absolute error
and is dimension-polymorphic. Relative and absolute tolerances are nonnegative
and cannot both vanish. ODE maximal steps are positive and finite at this ABI;
its initial step is chosen automatically within the requested interval.

`num_least_squares` is damped nonlinear least squares (Levenberg–Marquardt
style) over a scoped residual callback. The parameter and residual vectors have
independent dimensions `D` and `E`; the residual tolerance has dimension `E`, the
other three tolerances and the damping are dimensionless, and the parameter
vector must be nonempty. The Jacobian is formed by finite differences of the
callback (an analytic Jacobian is available only through the Rust kernel
`kernels::rootfind::least_squares`). `status` is one of `residual_converged`,
`gradient_converged`, `iteration_limit`, `evaluation_limit` or `stagnated`, and
`converged` is true only for the first two.

The state vector of an ODE operator must be nonempty. Every derivative has the same width. Concrete
callback extents, when present in its signature, are checked before executing
the host, and every returned tensor is checked for rank, shape/storage agreement,
width and physical dimension. The state result stores one row per accepted time
point. All four methods support reverse time. Native event, mass-matrix, DAE and
dense-output facilities remain available through the Rust solver APIs.

The first host or callback-result conversion error propagates unchanged, including
`OpError::HostDiagnostic`'s code. Once it fails, the operation never invokes the
host again. ODE adapters also stop adaptive retries through an internal cancellation
check; this does not unwind through host code. `evaluations` counts actual host
callback invocations. A root result does not evaluate its residual once more just
to populate a report. Brent's budget counts iterations after endpoint samples;
ODE budgets count attempted steps, including rejected attempts. Exhaustion is an
explicit `NotReduced` error. Quadrature retains its explicit `converged` result.

The quadrature and ODE estimators are heuristic floating-point error estimates.
Brent requires a bracketed sign change (or an endpoint root); its tolerance is in
the domain's units and makes no promise about residual magnitude.

## Shared contracts

`numerical_callbacks::common::CallbackRun` owns first-error preservation and
observable invocation counts for all adapters. Registration supplies the same
stateful/diverging effects, session-bound threading and interpreter backend
metadata to each operator. There is no symbolic rewrite, scalar JIT, GPU or C
function-pointer interpretation for these host handles.

`tools/rssn-api-audit` generates record codecs from public Rust structs. An
exhaustive typed destructuring pattern makes a source field change a compilation
error until its policy is reviewed. The quadrature `Integral` record and SVD
record each derive their manifest field types and runtime dimensional decoration
from one generated field policy. Generated codec availability and callable export
coverage are reported separately: a record layout by itself does not establish
its solver's dimensional contract.
