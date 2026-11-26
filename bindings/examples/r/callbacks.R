# Higher-order operators: R functions the engine calls.
#
#     RSSN_LIB=/path/to/rssn.dll Rscript bindings/examples/r/callbacks.R
#
# Some operators take a function and call it -- root finders, quadrature,
# optimisers, bootstrap.  In R any function works: the package turns it into the
# C function pointer the engine knows how to call, and converts the arguments
# and the result for you.

library(rssn)

root_of_two <- function() {
    # `sqrt(2)` as the root of `x^2 - 2`, by Brent's method.
    s <- rssn_session()
    answer <- call_value("num_brent",
                         1, 2,        # a, b: a bracketing interval
                         1e-12,       # tolerance
                         64,          # iteration budget
                         function(x) x * x - 2,
                         session = s)
    cat(sprintf("num_brent:        sqrt(2) = %.15g\n", value_float(answer[[1]])))
    cat(sprintf("                  evaluations = %d\n", as.integer(value_float(answer[[2]]))))
}

through_the_facade <- function() {
    # The generated facade routes a higher-order operator the same way.
    s <- rssn_session()
    answer <- rssn_num_brent(1, 2, 1e-12, 64, function(x) x * x - 2, session = s)
    cat(sprintf("rssn_num_brent:   %.15g\n", value_float(answer[[1]])))
}

minimise <- function() {
    # A two-dimensional objective and its gradient, as plain R functions.
    #
    # `num_bfgs(x0, tolerance, max_iterations, max_evaluations, objective,
    # gradient)`: the objective receives the point as a numeric vector and
    # answers a number; the gradient answers a vector of the same length.
    s <- rssn_session()
    objective <- function(point) sum((point - 1)^2)
    gradient <- function(point) 2 * (point - 1)
    answer <- call_value("num_bfgs", c(3, -2), 1e-10, 100, 1000,
                         objective, gradient, session = s)
    point <- as.numeric(c(value_float(answer[[1]][[1]]), value_float(answer[[1]][[2]])))
    cat(sprintf("num_bfgs:         minimum at (%.6f, %.6f)\n", point[[1]], point[[2]]))
    cat(sprintf("                  value = %.3g\n", value_float(answer[[2]])))
}

explicit_callback <- function() {
    # A callback that outlives one call, with its signature spelled out.
    #
    # `callback_signatures()` says what the engine expects at an argument
    # position; `rssn_callback()` registers one explicitly, which is useful when
    # the same function goes to several calls.
    s <- rssn_session()
    declared <- callback_signatures("num_brent", 5, session = s)
    cat(sprintf("declared at num_brent/5: %s\n", json_encode(declared)))

    residual <- rssn_callback(function(x) x * x - 2, s, signature = "(float) -> float")
    answer <- call_value("num_brent", 1, 2, 1e-12, 64, residual, session = s)
    cat(sprintf("explicit callback: %.15g\n", value_float(answer[[1]])))
    close(residual)
}

failures_cross <- function() {
    # An error in the callback aborts the operator with its message.
    s <- rssn_session()
    condition <- tryCatch({
        call_value("num_brent", 1, 2, 1e-12, 64,
                   function(x) stop("no thanks"), session = s)
        NULL
    }, rssn_error = function(condition) condition)
    if (!is.null(condition)) {
        cat(sprintf("callback failure: %s: %s\n", condition$status_name, condition$message))
    }
}

root_of_two()
through_the_facade()
minimise()
explicit_callback()
failures_cross()
