# Simulations: a scenario by name, with JSON parameters.
#
# The engine's simulation catalog is reached through one entry point: a scenario
# name and a JSON object of parameters. `{}` selects the defaults of the
# scenarios that have them, and `SCENARIOS` / `scenario_params()` describe what
# each one reads.
#
#     result <- sim_run("ising", width = 8, height = 8, temperature = 2.0,
#                       mc_steps = 10, seed = 3)
#     names(result)                       # "magnetization", "spins", ...
#     sim_run("poisson_multigrid", dim = 2, n = 33, tol = 1e-10)$converged

#' Run a simulation scenario
#'
#' @param scenario The scenario name; see [sim_names()].
#' @param ... The parameters, as named arguments (see [sim_params()] for the
#'   ones a scenario reads). A single named list is accepted too.
#' @return The scenario's JSON result, decoded: a named list.
#' @examples
#' sim_run("ising", width = 4, height = 4, temperature = 2, mc_steps = 5, seed = 7)
#' @export
sim_run <- function(scenario, ...) {
    .rssn_require()
    parameters <- list(...)
    if (length(parameters) == 1L && is.null(names(parameters)) && is.list(parameters[[1L]])) {
        parameters <- parameters[[1L]]
    }
    text <- .Call(C_rssn_sim_run, as.character(scenario), json_encode(parameters))
    json_decode(text)
}

#' Every scenario name
#'
#' @return A character vector, in catalog order.
#' @export
sim_names <- function() {
    vapply(SCENARIOS, function(scenario) scenario$name, character(1))
}

#' The parameters a scenario's request type reads
#'
#' @param name The scenario name.
#' @return A named list: parameter name -> `list(type = , default = )`. Empty
#'   when the catalog does not describe the scenario.
#' @export
sim_params <- function(name) {
    SCENARIO_PARAMS[[as.character(name)]] %||% list()
}

#' The parameters of a scenario with their declared defaults
#'
#' Only literals are recovered from the engine's source, so a field without a
#' literal default is left out: this is a starting point for [sim_run()], not a
#' complete schema.
#'
#' @param name The scenario name.
#' @return A named list, ready to pass to [sim_run()].
#' @examples
#' sim_defaults("poisson_multigrid")
#' @export
sim_defaults <- function(name) {
    parameters <- sim_params(name)
    out <- list()
    for (field in names(parameters)) {
        if (!is.null(parameters[[field]]$default)) {
            out[[field]] <- parameters[[field]]$default
        }
    }
    out
}
