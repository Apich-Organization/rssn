\name{sim_run}
\alias{sim_run}
\alias{sim_names}
\alias{sim_params}
\alias{sim_defaults}
\alias{SCENARIOS}
\alias{SCENARIO_PARAMS}
\title{Simulations}
\usage{
sim_run(scenario, ...)

sim_names()

sim_params(name)

sim_defaults(name)
}
\arguments{
  \item{scenario}{The scenario name; see \code{sim_names()}.}
  \item{...}{The parameters, as named arguments; a single named list is
    accepted too.}
  \item{name}{The scenario name.}
}
\value{
  \code{sim_run()} returns the scenario's result as a named list.
  \code{sim_names()} a character vector; \code{sim_params()} a named list of
  \code{list(type = , default = )}; \code{sim_defaults()} the parameters with
  the defaults the engine declares.
}
\description{
  The engine's simulation catalog is reached through one entry point: a scenario
  name and a JSON object of parameters.  \code{SCENARIOS} is the catalog -- name,
  group and the parameters each request type reads -- generated from the
  engine's own source.
}
\details{
  Only literal defaults are recovered from the engine's source, so a parameter
  without one is left out of \code{sim_defaults()}: it is a starting point, not
  a complete schema.
}
\examples{
sim_run("ising", width = 8, height = 8, temperature = 2.0, mc_steps = 10, seed = 3)

sim_names()                       # 29 scenarios
sim_params("poisson_multigrid")
sim_defaults("poisson_multigrid")

sim_run("poisson_multigrid", dim = 2, n = 33, tol = 1e-10)$converged
}
