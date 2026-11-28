# The computation configuration: target, bindings, assumptions, dimensions.
#
# A configuration is a handle with reference semantics, and every setter
# returns it, so a target reads as a pipe:
#
#     cfg <- rssn_config() |> cfg_numeric(1e-12) |> cfg_bind("x", 0.5)
#     compute(rssn_diff(rssn_sin(x), x), cfg)$value     # 0.877582561890373
#
# `cfg_bind()` binds a symbol for numeric evaluation, `cfg_assume()` declares a
# fact, `cfg_dimension()` declares a dimension and `cfg_dimension_mode()`
# decides what the engine does with an inconsistency.

#' The fact bits `cfg_assume()` accepts
#'
#' A named integer vector: `real` (1), `nonzero` (2), `nonnegative` (4),
#' `nonpositive` (8), `integer` (16), and the two composites `positive` (7) and
#' `negative` (11), which is how the engine defines them.
#'
#' @examples
#' assumption[["integer"]]   # 16L
#' @export
assumption <- c(real = 1L, nonzero = 2L, nonnegative = 4L, nonpositive = 8L,
                integer = 16L, positive = 7L, negative = 11L)

#' The dimension checking modes
#'
#' A named integer vector: `off` (0), `advisory` (1) and `strict` (2).
#' @export
dimension_mode <- c(off = 0L, advisory = 1L, strict = 2L)

#' A computation configuration
#'
#' An owned `RssnConfig`: the target of a computation, the symbol bindings, the
#' assumptions and the dimension declarations. It is single-threaded, like the
#' session it is used with, and is released by `close()` or by the garbage
#' collector.
#'
#' @return An object of class `rssn_config`.
#' @examples
#' cfg <- rssn_config() |> cfg_numeric(1e-12) |> cfg_bind("x", 0)
#' @export
rssn_config <- function() {
    .rssn_require()
    structure(list(ptr = .Call(C_rssn_config_new)), class = "rssn_config")
}

#' The raw external pointer of a configuration
#' @param config A `rssn_config`.
#' @return The external pointer the C bridge holds.
#' @export
config_ptr <- function(config) {
    if (!inherits(config, "rssn_config")) stop("expected a rssn_config")
    config$ptr
}

# `NULL` stays `NULL`: the C bridge reads it as "the default configuration".
.rssn_config_ptr <- function(config) {
    if (is.null(config)) return(NULL)
    config_ptr(config)
}

#' @export
print.rssn_config <- function(x, ...) {
    cat("<rssn config>\n")
    invisible(x)
}

#' @export
close.rssn_config <- function(con, ...) {
    if (inherits(con, "rssn_config") && !is.null(con$ptr)) {
        .Call(C_rssn_config_free, con$ptr)
        con$ptr <- NULL
    }
    invisible(NULL)
}

#' Ask for a closed form
#'
#' The default target: reduce the request symbolically.
#'
#' @param config A `rssn_config`.
#' @return `config`, invisibly.
#' @export
cfg_symbolic <- function(config) {
    .Call(C_rssn_config_symbolic, config_ptr(config))
    invisible(config)
}

#' Ask for a numeric answer
#'
#' @param config A `rssn_config`.
#' @param tolerance The tolerance the answer must meet.
#' @return `config`, invisibly.
#' @examples
#' cfg <- rssn_config() |> cfg_numeric(1e-12)
#' @export
cfg_numeric <- function(config, tolerance = 1e-12) {
    .Call(C_rssn_config_numeric, config_ptr(config), as.numeric(tolerance))
    invisible(config)
}

#' Bind a symbol to a number for numeric evaluation
#'
#' @param config A `rssn_config`.
#' @param symbol The symbol.
#' @param value A number.
#' @return `config`, invisibly.
#' @examples
#' cfg <- rssn_config() |> cfg_numeric(1e-12) |> cfg_bind("x", 0.5)
#' @export
cfg_bind <- function(config, symbol, value) {
    .Call(C_rssn_config_bind, config_ptr(config), as.character(symbol), as.numeric(value))
    invisible(config)
}

#' Declare a fact about a symbol
#'
#' @param config A `rssn_config`.
#' @param symbol The symbol.
#' @param facts An [assumption] bit set, a name such as `"positive"`, or a
#'   character vector of names.
#' @return `config`, invisibly.
#' @examples
#' cfg <- rssn_config() |> cfg_assume("n", "integer")
#' @export
cfg_assume <- function(config, symbol, facts) {
    .Call(C_rssn_config_assume, config_ptr(config), as.character(symbol), .rssn_facts(facts))
    invisible(config)
}

.rssn_facts <- function(facts) {
    if (is.numeric(facts)) return(as.integer(sum(facts)))
    if (is.character(facts)) {
        unknown <- setdiff(tolower(facts), names(assumption))
        if (length(unknown)) {
            stop(sprintf("unknown assumption '%s'; use one of %s", unknown[[1L]],
                         paste(names(assumption), collapse = ", ")))
        }
        return(as.integer(sum(assumption[tolower(facts)])))
    }
    stop("cfg_assume() takes a name, a set of names, or an assumption bit set")
}

#' Limit the engine's search
#'
#' @param config A `rssn_config`.
#' @param max_iterations Outer iterations; `0` keeps the current limit.
#' @param max_nodes Nodes; `0` keeps the current limit.
#' @return `config`, invisibly.
#' @export
cfg_budget <- function(config, max_iterations = 0, max_nodes = 0) {
    .Call(C_rssn_config_budget, config_ptr(config), as.numeric(max_iterations),
          as.numeric(max_nodes))
    invisible(config)
}

#' Set the dimension checking mode
#'
#' @param config A `rssn_config`.
#' @param mode A [dimension_mode] value or name: `"off"`, `"advisory"` or
#'   `"strict"`.
#' @return `config`, invisibly.
#' @export
cfg_dimension_mode <- function(config, mode = "advisory") {
    if (is.character(mode)) {
        if (!mode[[1L]] %in% names(dimension_mode)) {
            stop(sprintf("unknown dimension mode '%s'; use one of %s", mode[[1L]],
                         paste(names(dimension_mode), collapse = ", ")))
        }
        mode <- dimension_mode[[mode[[1L]]]]
    }
    .Call(C_rssn_config_dimension_mode, config_ptr(config), as.integer(mode))
    invisible(config)
}

#' Disable hard dimensional errors even in strict mode
#'
#' @param config A `rssn_config`.
#' @return `config`, invisibly.
#' @export
cfg_bypass_dimensions <- function(config) {
    .Call(C_rssn_config_bypass_dimensions, config_ptr(config))
    invisible(config)
}

#' Declare the dimension of a symbol
#'
#' @param config A `rssn_config`.
#' @param symbol The symbol.
#' @param text The dimension as text: `"m/s"`, `"Length/Time"`, `"any"`.
#' @return `config`, invisibly.
#' @examples
#' cfg <- rssn_config() |> cfg_dimension("x", "m") |> cfg_dimension("t", "s")
#' @export
cfg_dimension <- function(config, symbol, text) {
    .Call(C_rssn_config_dimension, config_ptr(config), as.character(symbol),
          as.character(text))
    invisible(config)
}

#' Register a base dimension the engine does not have
#'
#' @param config A `rssn_config`.
#' @param name The base dimension, for example `"Currency"`.
#' @return `config`, invisibly.
#' @export
cfg_define_base_dimension <- function(config, name) {
    .Call(C_rssn_config_define_base_dimension, config_ptr(config), as.character(name))
    invisible(config)
}
