# Sessions, terms, answers and computation.
#
# The model is the engine's: an operator applied to arguments builds an
# *unevaluated request*, and `compute()` reduces it. `rssn_diff(sin(x) * x, x)`
# is a term you can print, differentiate again or substitute into; only
# `compute()` runs the rules.
#
#     s <- rssn_session()
#     x <- sym("x", s)
#     t <- rssn_diff(sin(x) * cos(x), x)   # a request; prints "diff(sin(x)*cos(x), x)"
#     compute(t)$text                      # "cos(2*x)"
#
# A session is single-threaded: its store sits behind a `RefCell` in Rust, so
# two threads using one session at once is a borrow conflict. Use one session
# per thread, or `compute_all()` for parallel requests inside one session.

# ---------------------------------------------------------------------------
# sessions
# ---------------------------------------------------------------------------

#' A rssn session
#'
#' A session owns every term made through it. It is single-threaded, and it is
#' released by `close()`, by the garbage collector, or with the session itself.
#'
#' @return An object of class `rssn_session`.
#' @examples
#' s <- rssn_session()
#' x <- sym("x", s)
#' compute(rssn_diff(x^3, x))$text   # "3*x^2"
#' close(s)
#' @export
rssn_session <- function() {
    .rssn_require()
    structure(list(ptr = .Call(C_rssn_session_new)), class = "rssn_session")
}

#' The raw external pointer of a session
#' @param session A `rssn_session`.
#' @return The external pointer the C bridge holds.
#' @export
session_ptr <- function(session) {
    if (!inherits(session, "rssn_session")) {
        stop("expected a rssn_session")
    }
    session$ptr
}

#' @export
print.rssn_session <- function(x, ...) {
    cat("<rssn session>\n")
    invisible(x)
}

#' @export
format.rssn_session <- function(x, ...) "<rssn session>"

#' Release a session and every term made through it
#'
#' Safe to call twice. This extends the `close()` generic, so `close(s)` needs
#' no namespace qualification.
#'
#' @param con A `rssn_session`.
#' @param ... Ignored.
#' @return `NULL`, invisibly.
#' @export
close.rssn_session <- function(con, ...) {
    if (inherits(con, "rssn_session") && !is.null(con$ptr)) {
        .Call(C_rssn_session_free, con$ptr)
        con$ptr <- NULL
    }
    invisible(NULL)
}

#' @rdname close.rssn_session
#' @export
rssn_close <- function(session) {
    close(session)
}

#: The calling thread's session, created on first use.
.rssn_default <- local({
    value <- NULL
    list(
        get = function() {
            if (is.null(value) || is.null(value$ptr)) {
                value <<- rssn_session()
            }
            value
        },
        set = function(session) {
            previous <- value
            value <<- session
            previous
        }
    )
})

#' The current default session
#'
#' The session the operator functions use when no argument carries one. It is
#' created on first use and can be replaced by [with_session()].
#'
#' @return A `rssn_session`.
#' @export
default_session <- function() {
    .rssn_default$get()
}

#' Run code with a given default session
#'
#' @param session A `rssn_session`.
#' @param code An expression, evaluated in the caller's environment.
#' @return The value of `code`.
#' @examples
#' s <- rssn_session()
#' x <- with_session(s, sym("x"))
#' with_session(s, compute(rssn_diff(x^2, x))$text)   # "2*x"
#' @export
with_session <- function(session, code) {
    if (!inherits(session, "rssn_session")) stop("expected a rssn_session")
    expression <- substitute(code)
    previous <- .rssn_default$set(session)
    on.exit(.rssn_default$set(previous), add = TRUE)
    eval(expression, parent.frame())
}

# ---------------------------------------------------------------------------
# terms
# ---------------------------------------------------------------------------

.rssn_term <- function(session, handle) {
    structure(as.integer(handle), class = "rssn_term", session = session)
}

#' The raw handle of a term
#' @param term A `rssn_term`.
#' @return The `uint32_t` handle as an integer.
#' @export
handle <- function(term) {
    as.integer(unclass(term))
}

.rssn_session_of <- function(term) attr(term, "session")

#' The symbol `name`
#'
#' @param symbol The symbol's name.
#' @param session The session to build the term in; the default session when
#'   omitted.
#' @return A `rssn_term`.
#' @examples
#' x <- sym("x")
#' compute(x^2)$text
#' @export
sym <- function(symbol, session = NULL) {
    session <- session %||% default_session()
    .rssn_term(session, .Call(C_rssn_term_sym, session$ptr, as.character(symbol)))
}

#' Parse infix text into a term
#'
#' @param source Infix text such as `"diff(sin(x)^2, x)"`.
#' @param session The session to build the term in.
#' @return A `rssn_term`.
#' @examples
#' parse_term("x^2 + 2*x + 1")
#' @export
parse_term <- function(source, session = NULL) {
    session <- session %||% default_session()
    .rssn_term(session, .Call(C_rssn_term_parse, session$ptr, as.character(source)))
}

#' An exact integer term
#' @param value A number with no fractional part.
#' @param session The session to build the term in.
#' @return A `rssn_term`.
#' @export
int_term <- function(value, session = NULL) {
    session <- session %||% default_session()
    .rssn_term(session, .Call(C_rssn_term_int, session$ptr, as.numeric(value)))
}

#' A floating point term
#' @param value A number.
#' @param session The session to build the term in.
#' @return A `rssn_term`.
#' @export
real_term <- function(value, session = NULL) {
    session <- session %||% default_session()
    .rssn_term(session, .Call(C_rssn_term_float, session$ptr, as.numeric(value)))
}

#' An exact rational term
#'
#' @param numerator An integer.
#' @param denominator An integer; zero is an error.
#' @param session The session to build the term in.
#' @return A `rssn_term`.
#' @examples
#' rational(1, 3)
#' @export
rational <- function(numerator, denominator = 1, session = NULL) {
    session <- session %||% default_session()
    .rssn_term(session,
               .Call(C_rssn_term_rational, session$ptr, as.numeric(numerator),
                     as.numeric(denominator)))
}

#' Coerce a value to a term
#'
#' A `rssn_term` is returned as is (it must belong to `session`); a number
#' becomes an exact integer or a float; a string is parsed as infix text.
#'
#' @param x A term, a number, or infix text.
#' @param session The session to build the term in.
#' @return A `rssn_term`.
#' @export
term <- function(x, session = NULL) {
    if (inherits(x, "rssn_term")) {
        if (!is.null(session) && !identical(.rssn_session_of(x)$ptr, session$ptr)) {
            stop("the term belongs to another session")
        }
        return(x)
    }
    session <- session %||% default_session()
    if (is.character(x)) return(parse_term(x, session))
    if (is.logical(x)) return(int_term(as.integer(x), session))
    if (is.numeric(x)) {
        if (length(x) != 1L) stop("a term is a scalar; use a list for a tensor")
        if (!is.na(x) && x == round(x) && abs(x) < 2^53) return(int_term(x, session))
        return(real_term(x, session))
    }
    stop(sprintf("cannot use an object of class '%s' as a term",
                 paste(class(x), collapse = "/")))
}

# -- renderings -------------------------------------------------------------

#' The term as infix text
#'
#' The rendering the engine parses back: `term_text()` is what `as.character()`
#' returns, so a term round-trips through [parse_term()].
#'
#' @param term A `rssn_term`.
#' @return A string.
#' @examples
#' term_text(sym("x")^2)
#' @export
term_text <- function(term) {
    .Call(C_rssn_term_to_string, .rssn_session_of(term)$ptr, handle(term))
}

#' The term as LaTeX math
#' @param term A `rssn_term`.
#' @return A string.
#' @examples
#' term_latex(sym("x")^2)    # "x^{2}"
#' @export
term_latex <- function(term) {
    .Call(C_rssn_term_to_latex, .rssn_session_of(term)$ptr, handle(term))
}

#' The term as Typst math
#' @param term A `rssn_term`.
#' @return A string.
#' @export
term_typst <- function(term) {
    .Call(C_rssn_term_to_typst, .rssn_session_of(term)$ptr, handle(term))
}

#' The term as a multi-line Unicode drawing
#' @param term A `rssn_term`.
#' @return A string.
#' @export
term_pretty <- function(term) {
    .Call(C_rssn_term_to_pretty, .rssn_session_of(term)$ptr, handle(term))
}

#' @export
as.character.rssn_term <- function(x, ...) term_text(x)

#' @export
format.rssn_term <- function(x, ...) term_text(x)

#' @export
print.rssn_term <- function(x, ...) {
    cat("<rssn term> ", text(x), "\n", sep = "")
    invisible(x)
}

#' The value of a literal term
#'
#' @param term A `rssn_term`.
#' @return A number, or `NA` when the term is symbolic.
#' @export
term_value <- function(term) {
    .Call(C_rssn_term_to_f64, .rssn_session_of(term)$ptr, handle(term))
}

#' The numeric value of a term, computing it when it is not a literal
#'
#' @param term A `rssn_term`.
#' @param tolerance The tolerance of the numeric target; `1e-12` by default.
#' @return A number.
#' @export
to_float <- function(term, tolerance = 1e-12) {
    literal <- term_value(term)
    if (!is.na(literal)) return(literal)
    answer <- compute(term, rssn_config() |> cfg_numeric(tolerance))
    if (is.na(answer$value)) {
        stop("the term has no numeric value")
    }
    answer$value
}

#' @export
as.numeric.rssn_term <- function(x, ...) to_float(x)

# -- arithmetic, as requests ------------------------------------------------

#: `Ops` group members the engine names itself.
.OPS_NAMES <- c(
    "+" = "add", "-" = "sub", "*" = "mul", "/" = "divide", "^" = "pow",
    "==" = "eq", "!=" = "ne", "<" = "lt", "<=" = "le", ">" = "gt", ">=" = "ge",
    "&" = "and", "|" = "or"
)

#' @export
Ops.rssn_term <- function(e1, e2) {
    generic <- .Generic
    if (missing(e2)) {
        switch(generic,
            "-" = return(.rssn_apply("neg", list(e1), NULL)),
            "+" = return(e1),
            "!" = return(.rssn_apply("not", list(e1), NULL)),
            stop(sprintf("unary '%s' is not defined for a term", generic))
        )
    }
    operator <- .OPS_NAMES[[generic]]
    if (is.null(operator)) {
        stop(sprintf("'%s' has no rssn operator; call op_apply() by name", generic))
    }
    .rssn_apply(operator, list(e1, e2), NULL)
}

#: `Math` group members the engine has, by engine name. `log()` is the one with
#: two arities: R's `log(x)` is the natural logarithm, which the engine spells
#: `ln`, and R's `log(x, base)` is the engine's `log(base, x)` -- base first.
.MATH_NAMES <- c(
    abs = "abs", sign = "sign", sqrt = "sqrt", floor = "floor", ceiling = "ceil",
    round = "round", exp = "exp", cos = "cos", sin = "sin",
    tan = "tan", acos = "acos", asin = "asin", atan = "atan", cosh = "cosh",
    sinh = "sinh", tanh = "tanh", acosh = "acosh", asinh = "asinh",
    atanh = "atanh", gamma = "gamma", lgamma = "lgamma", digamma = "digamma"
)

#' @export
Math.rssn_term <- function(x, ...) {
    extra <- list(...)
    if (.Generic == "log") {
        if (length(extra) == 0L) return(.rssn_apply("ln", list(x), NULL))
        # `log(x, base)` in R, `log(base, x)` in the engine.
        return(.rssn_apply("log", list(extra[[1L]], x), NULL))
    }
    operator <- .MATH_NAMES[[.Generic]]
    if (is.null(operator)) {
        stop(sprintf("'%s()' has no rssn operator; call op_apply() by name", .Generic))
    }
    .rssn_apply(operator, c(list(x), extra), NULL)
}

#' @export
Summary.rssn_term <- function(..., na.rm = FALSE) {
    arguments <- list(...)
    if (.Generic != "sum") {
        stop(sprintf("'%s()' has no rssn operator", .Generic))
    }
    .rssn_apply("sum", arguments, NULL)
}

# ---------------------------------------------------------------------------
# answers
# ---------------------------------------------------------------------------

#' The answer of a computation
#'
#' A list with the reduced term and everything the engine reported about it:
#' `term`, `text`, `value` (`NA` when unknown), `error`, `reduced`,
#' `iterations`, `dimension` (`NULL` when none) and `diagnostic_codes`.
#'
#' @param term A `rssn_term`.
#' @param text The answer as infix text.
#' @param value The numeric value, or `NA`.
#' @param error The estimated absolute error, or `NA`.
#' @param reduced Whether the request was fully reduced.
#' @param iterations Outer iterations the engine performed.
#' @param dimension The inferred dimension as text, or `NULL`.
#' @param diagnostic_codes Stable codes of the advisory dimensional findings.
#' @return An object of class `rssn_answer`.
#' @export
answer <- function(term, text, value = NA_real_, error = NA_real_, reduced = TRUE,
                   iterations = 0, dimension = NULL, diagnostic_codes = character()) {
    structure(
        list(term = term, text = text, value = value, error = error, reduced = reduced,
             iterations = iterations, dimension = dimension,
             diagnostic_codes = diagnostic_codes),
        class = "rssn_answer"
    )
}

#' @export
is_answer <- function(x) inherits(x, "rssn_answer")

#' @export
print.rssn_answer <- function(x, ...) {
    cat("<rssn answer> ", x$text, "\n", sep = "")
    if (!is.na(x$value)) {
        cat("  value      : ", format(x$value, digits = 15), "\n", sep = "")
        if (!is.na(x$error) && x$error > 0) {
            cat("  error      : ", format(x$error, digits = 3), "\n", sep = "")
        }
    }
    if (!isTRUE(x$reduced)) cat("  reduced    : FALSE\n")
    if (!is.null(x$dimension)) cat("  dimension  : ", x$dimension, "\n", sep = "")
    if (length(x$diagnostic_codes)) {
        cat("  findings   : ", paste(x$diagnostic_codes, collapse = ", "), "\n", sep = "")
    }
    invisible(x)
}

#' @export
format.rssn_answer <- function(x, ...) x$text

# Reads the thread's dimension report, then the text of the answer term.
.rssn_answer_from <- function(session, term, has_value, value, error, reduced, iterations,
                              config = NULL) {
    findings <- .Call(C_rssn_dimension_count)
    codes <- if (findings > 0L) {
        vapply(seq_len(findings) - 1L, function(index) .Call(C_rssn_dimension_code, index),
               character(1))
    } else {
        character()
    }
    inferred <- .Call(C_rssn_dimension, session$ptr, term, .rssn_config_ptr(config))
    answer(
        term = .rssn_term(session, term),
        text = .Call(C_rssn_term_to_string, session$ptr, term),
        value = if (has_value) value else NA_real_,
        error = if (has_value) error else NA_real_,
        reduced = reduced,
        iterations = iterations,
        dimension = if (nzchar(inferred)) inferred else NULL,
        diagnostic_codes = codes
    )
}

#' Reduce a request to the target of a configuration
#'
#' @param x A term, or infix text.
#' @param config A `rssn_config`; the symbolic default when omitted.
#' @param session The session to parse `x` in, when `x` is text.
#' @return An `rssn_answer`.
#' @examples
#' compute("diff(sin(x)*cos(x), x)")$text                       # "cos(2*x)"
#' compute("diff(sin(x)*cos(x), x)",
#'         rssn_config() |> cfg_numeric(1e-12) |> cfg_bind("x", 0))$value
#' @export
compute <- function(x, config = NULL, session = NULL) {
    term <- if (inherits(x, "rssn_term")) {
        x
    } else {
        term(x, session %||% default_session())
    }
    session <- .rssn_session_of(term)
    raw <- .Call(C_rssn_compute, session$ptr, handle(term), .rssn_config_ptr(config))
    .rssn_answer_from(session, raw$term, raw$has_value, raw$value, raw$error,
                      raw$reduced, raw$iterations)
}

#' Parse and reduce text with the default symbolic configuration
#'
#' @param source Infix text.
#' @param session The session to parse it in.
#' @return An `rssn_answer`.
#' @examples
#' simplify("(x^2 - 1)/(x - 1)")$text    # "x + 1"
#' @export
simplify <- function(source, session = NULL) {
    session <- session %||% default_session()
    term <- .rssn_term(session, .Call(C_rssn_simplify, session$ptr, as.character(source)))
    compute(term)
}

#' Compute several requests in parallel inside the library
#'
#' @param requests A list of terms or infix strings.
#' @param config A `rssn_config`.
#' @param session The session to parse the requests in.
#' @param errors What to do with a request that failed: `"raise"` (the default)
#'   raises the first failure, `"collect"` puts the condition in its place and
#'   `"ignore"` puts `NULL` there.
#' @return A list of `rssn_answer` objects, one per request.
#' @examples
#' compute_all(list("diff(x^3, x)", "integral(x^2, x)", "gcd(4, 6)"))
#' @export
compute_all <- function(requests, config = NULL, session = NULL,
                        errors = c("raise", "collect", "ignore")) {
    errors <- match.arg(errors)
    requests <- as.list(requests)
    if (length(requests) == 0L) return(list())
    terms <- lapply(requests, function(request) term(request, session))
    session <- .rssn_session_of(terms[[1L]])
    handles <- vapply(terms, handle, integer(1))
    raw <- .Call(C_rssn_compute_all, session$ptr, handles, .rssn_config_ptr(config))
    first_error <- NULL
    out <- vector("list", length(terms))
    for (index in seq_along(terms)) {
        status <- raw$status[[index]]
        if (status == 0L) {
            out[[index]] <- .rssn_answer_from(
                session, raw$term[[index]], raw$has_value[[index]], raw$value[[index]],
                raw$error[[index]], raw$reduced[[index]], raw$iterations[[index]])
        } else {
            condition <- .rssn_condition(status)
            if (is.null(first_error)) first_error <- condition
            out[[index]] <- switch(errors, raise = NULL, collect = condition, ignore = NULL)
        }
    }
    if (errors == "raise" && !is.null(first_error)) stop(first_error)
    out
}

#' The compute record as a list
#'
#' The engine's `RssnComputeRecord`, decoded: `term`, `value`, `error`,
#' `reduced`, `iterations`, `dimension` and `diagnostics` (the prose of each
#' advisory finding).
#'
#' @param x A term or infix text.
#' @param config A `rssn_config`.
#' @param session The session to parse `x` in.
#' @return A named list.
#' @export
compute_json <- function(x, config = NULL, session = NULL) {
    term <- if (inherits(x, "rssn_term")) x else term(x, session %||% default_session())
    session <- .rssn_session_of(term)
    json_decode(.Call(C_rssn_compute_json, session$ptr, handle(term),
                      .rssn_config_ptr(config)))
}

#' The compute record in the deterministic `bincode-next` format
#'
#' @param x A term or infix text.
#' @param config A `rssn_config`.
#' @param session The session to parse `x` in.
#' @return A raw vector.
#' @export
compute_bincode <- function(x, config = NULL, session = NULL) {
    term <- if (inherits(x, "rssn_term")) x else term(x, session %||% default_session())
    session <- .rssn_session_of(term)
    .Call(C_rssn_compute_bincode, session$ptr, handle(term), .rssn_config_ptr(config))
}

#' The batch record in the `bincode-next` format
#'
#' @param requests A list of terms or infix strings.
#' @param config A `rssn_config`.
#' @param session The session to parse the requests in.
#' @return A raw vector.
#' @export
compute_all_bincode <- function(requests, config = NULL, session = NULL) {
    terms <- lapply(as.list(requests), function(request) term(request, session))
    session <- .rssn_session_of(terms[[1L]])
    handles <- vapply(terms, handle, integer(1))
    .Call(C_rssn_compute_all_bincode, session$ptr, handles, .rssn_config_ptr(config))
}

#' Evaluate a term numerically without running any rule
#'
#' @param term A `rssn_term`.
#' @param ... Named numbers binding the free symbols, such as `x = 0.5`.
#' @return A number.
#' @examples
#' eval_term(parse_term("sin(x) + 1"), x = 0)
#' @export
eval_term <- function(term, ...) {
    bindings <- list(...)
    if (length(bindings) == 0L) {
        return(.Call(C_rssn_term_eval, .rssn_session_of(term)$ptr, handle(term),
                     character(), numeric()))
    }
    if (is.null(names(bindings)) || any(!nzchar(names(bindings)))) {
        stop("eval_term() takes named bindings, such as x = 0.5")
    }
    .Call(C_rssn_term_eval, .rssn_session_of(term)$ptr, handle(term), names(bindings),
          as.numeric(bindings))
}

#' The inferred dimension of a request
#'
#' @param x A term or infix text.
#' @param config A `rssn_config`, usually with `cfg_dimension()` declarations.
#' @param session The session to parse `x` in.
#' @return A string; `""` when no dimension was inferred.
#' @examples
#' cfg <- rssn_config() |> cfg_dimension("x", "m")
#' x <- sym("x")
#' dimension(rssn_diff(x^2, x), cfg)   # "Length"
#' @export
dimension <- function(x, config = NULL, session = NULL) {
    term <- if (inherits(x, "rssn_term")) x else term(x, session %||% default_session())
    .Call(C_rssn_dimension, .rssn_session_of(term)$ptr, handle(term),
          .rssn_config_ptr(config))
}

#' The dimensional findings of the last computation on this thread
#'
#' @return A character vector of stable codes (`"mismatch"`, `"unknown_unit"`,
#'   `"unannotated"`, ...); the prose of each finding is in the `"diagnostics"`
#'   field of [compute_json()].
#' @export
dimension_report <- function() {
    .rssn_require()
    findings <- .Call(C_rssn_dimension_count)
    if (findings == 0L) return(character())
    vapply(seq_len(findings) - 1L, function(index) .Call(C_rssn_dimension_code, index),
           character(1))
}

# ---------------------------------------------------------------------------
# the operator facade's entry points
# ---------------------------------------------------------------------------

# Coerces one argument of an operator to a term handle in `session`.
.rssn_handle <- function(x, session) {
    if (inherits(x, "rssn_term")) {
        if (!identical(.rssn_session_of(x)$ptr, session$ptr)) {
            stop("the term belongs to another session")
        }
        return(handle(x))
    }
    handle(term(x, session))
}

# The one function every generated operator calls.
.rssn_apply <- function(op, args, session) {
    if (is.null(session)) {
        for (argument in args) {
            if (inherits(argument, "rssn_term")) {
                session <- .rssn_session_of(argument)
                break
            }
        }
        if (is.null(session)) session <- default_session()
    }
    handles <- vapply(args, .rssn_handle, integer(1), session = session)
    .rssn_term(session, .Call(C_rssn_term_apply, session$ptr, op, handles))
}

#' Apply an operator by name
#'
#' The escape hatch behind the generated `rssn_<name>()` functions: it reaches
#' any operator at run time, including one whose name is computed.
#'
#' @param op The operator name, as the manifest spells it (`"diff"`).
#' @param ... The arguments.
#' @param session The session to build the request in.
#' @return A `rssn_term`.
#' @examples
#' op_apply("diff", parse_term("x^3"), parse_term("x"))
#' @export
op_apply <- function(op, ..., session = NULL) {
    .rssn_apply(op, list(...), session)
}

#: `%||%`, the null-coalescing operator base R only gained in 4.4.
`%||%` <- function(left, right) if (is.null(left)) right else left
