# A tour of the simulation catalog, from R.
#
#     RSSN_LIB=/path/to/rssn.dll Rscript inst/examples/simulation.R

library(rssn)

section <- function(title) cat("\n==", title, "==\n")

section("the catalog")
cat(length(sim_names()), "scenarios:\n")
cat(strwrap(paste(sim_names(), collapse = ", "), width = 78, prefix = "  "), sep = "\n")

section("what a scenario reads")
print(sim_params("poisson_multigrid"))
print(sim_defaults("poisson_multigrid"))

section("an Ising lattice")
result <- sim_run("ising", width = 8, height = 8, temperature = 2.0,
                  mc_steps = 10, seed = 3)
cat("magnetization:", result$magnetization, "\n")
cat("spins:\n")
spins <- matrix(unlist(result$spins), nrow = 8, byrow = TRUE)
print(spins)

section("a multigrid solve")
result <- sim_run("poisson_multigrid", dim = 2, n = 33, tol = 1e-10)
cat("converged:", result$converged, " levels:", result$levels,
    " max error:", result$max_error, "\n")

section("a stochastic model")
result <- sim_run("spin_model", dim = 2, n = 16, beta = 0.5, sweeps = 100,
                  model = "ising", update = "wolff", seed = 1)
cat("energy per site:", result$energy_per_site,
    " |magnetisation| per site:", result$abs_magnetisation_per_site, "\n")

section("a failure")
condition <- tryCatch({
    sim_run("no_such_scenario")
    NULL
}, rssn_error = function(condition) condition)
cat("status:", condition$status_name, "\n")
