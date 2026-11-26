# Dynamic values: the cell ABI behind operators that take native objects.
#
# The term API carries only values that are nodes of the expression graph. An
# operator that returns a native object -- an AD tape, a compiled backend, a
# solver handle -- returns a *value* instead, and an operator that takes one can
# only be called with a value. Those are the operators `OPS` marks with
# `ffi = FALSE`; they are reached with `call_value()`.
#
#     s <- rssn_session()
#     backend <- call_value("num_backend_current", session = s)
#     value_kind(backend)                                                 # "opaque"
#     value_text(call_value("num_backend_name", backend, session = s))     # "interpreter"
#
#     tape <- call_value("num_ad_tape_new", 16, session = s)
#     x <- call_value("num_ad_tape_var", tape, 1.5, session = s)
#     z <- call_value("num_ad_var_mul", call_value("num_ad_var_sin", x, session = s), x,
#                     session = s)
#     value_float(call_value("num_ad_var_value", z, session = s))   # 1.4962424799060816
#
# A value is an environment with a finalizer, not a bare integer: the cell id
# the engine hands out can be reused after it is freed, so an object that
# remembers whether it still owns its cell is the only safe representation.
# `close()` releases the cell early; otherwise the garbage collector does it,
# and the session releases whatever is left when it is closed.

#: The kind codes `rssn_value_kind` reports, in declaration order.
VALUE_KINDS <- c(
    "int", "real", "float", "complex", "bool", "string", "symbol", "unit",
    "quantity", "tensor", "sparse", "record", "variant", "term", "function",
    "dyn", "opaque"
)

.rssn_value <- function(session, cell, owned = TRUE) {
    value <- new.env(parent = emptyenv())
    value$cell <- as.integer(cell)
    value$session <- session
    value$owned <- owned
    class(value) <- "rssn_value"
    if (owned) {
        reg.finalizer(value, .rssn_finalize_value, onexit = TRUE)
    }
    value
}

.rssn_finalize_value <- function(value) {
    if (isTRUE(value$owned) && value$cell != 0L) {
        session <- value$session
        if (!is.null(session) && !is.null(session$ptr)) {
            try(.Call(C_rssn_value_free, session$ptr, value$cell), silent = TRUE)
        }
        value$cell <- 0L
        value$owned <- FALSE
    }
    invisible(NULL)
}

# The cell id of a value, or a clear error once it has been released.
value_cell <- function(value) {
    if (value$cell == 0L) stop("this value has been released")
    value$cell
}

.rssn_value_session <- function(value) value$session

#' @export
print.rssn_value <- function(x, ...) {
    if (x$cell == 0L) {
        cat("<rssn value> released\n")
    } else {
        cat("<rssn value> ", value_kind(x), ": ", value_text(x), "\n", sep = "")
    }
    invisible(x)
}

#' @export
format.rssn_value <- function(x, ...) {
    if (x$cell == 0L) "released" else value_text(x)
}

#' @export
as.character.rssn_value <- function(x, ...) format(x)

#' @export
as.numeric.rssn_value <- function(x, ...) value_float(x)

#' @export
length.rssn_value <- function(x) value_len(x)

#' @export
`[[.rssn_value` <- function(x, index) {
    index <- as.integer(index)
    length <- value_len(x)
    if (index < 1L || index > length) stop("subscript out of bounds")
    .rssn_value(x$session, .Call(C_rssn_value_at, x$session$ptr, value_cell(x), index - 1L))
}

#' Release a value cell early
#'
#' The cell is released by the garbage collector anyway, and the session
#' releases whatever is left when it is closed; `close()` is for a long loop
#' that would otherwise hold many cells at once. Safe to call twice.
#'
#' @param con A `rssn_value`.
#' @param ... Ignored.
#' @return `NULL`, invisibly.
#' @export
close.rssn_value <- function(con, ...) {
    .rssn_finalize_value(con)
    invisible(NULL)
}

#' The kind of a value
#'
#' @param value A `rssn_value`.
#' @return One of [VALUE_KINDS]: `"float"`, `"tensor"`, `"opaque"`, ...
#' @export
value_kind <- function(value) {
    code <- .Call(C_rssn_value_kind, value$session$ptr, value_cell(value))
    if (code < length(VALUE_KINDS)) VALUE_KINDS[[code + 1L]] else paste0("kind", code)
}

#' The numeric value of a value
#' @param value A `rssn_value`.
#' @return A number.
#' @export
value_float <- function(value) {
    .Call(C_rssn_value_to_f64, value$session$ptr, value_cell(value))
}

#' The text of a string or symbol value, or any value's display form
#' @param value A `rssn_value`.
#' @return A string.
#' @export
value_text <- function(value) {
    .Call(C_rssn_value_to_text, value$session$ptr, value_cell(value))
}

#' The number of elements of a list, tensor or record value
#' @param value A `rssn_value`.
#' @return An integer.
#' @export
value_len <- function(value) {
    .Call(C_rssn_value_len, value$session$ptr, value_cell(value))
}

#' The value in its `serde` JSON form, decoded
#' @param value A `rssn_value`.
#' @return A list, vector or scalar.
#' @export
value_json <- function(value) {
    json_decode(.Call(C_rssn_value_json, value$session$ptr, value_cell(value)))
}

#' The value as a term, when it has a graph form
#' @param value A `rssn_value`.
#' @return A `rssn_term`.
#' @export
value_term <- function(value) {
    .rssn_term(value$session, .Call(C_rssn_value_to_term, value$session$ptr, value_cell(value)))
}

#' A dense extraction of a tensor value
#'
#' The engine stores tensors row-major and R stores arrays column-major, so the
#' extents are permuted on the way out: the answer reads like the matrix it
#' describes.
#'
#' @param value A `rssn_value`.
#' @return A numeric array.
#' @export
value_tensor <- function(value) {
    raw <- .Call(C_rssn_value_to_tensor, value$session$ptr, value_cell(value))
    .rssn_array(raw$data, raw$dim)
}

#' A dense extraction of a term's list value
#'
#' The shortcut for `value_tensor(value_new(x))`: a matrix or a nested list
#' term becomes the array it describes, in reading order.
#'
#' @param x A `rssn_term` or an `rssn_answer`.
#' @return A numeric array.
#' @examples
#' tensor(compute(parse_term("matmul(list(list(1,2),list(3,4)), list(list(1,0),list(0,1)))")))
#' @export
tensor <- function(x) {
    if (inherits(x, "rssn_answer")) x <- x$term
    value_tensor(value_new(x, .rssn_session_of(x)))
}

# Row-major data plus extents as an R array (column-major, permuted back).
.rssn_array <- function(data, dim) {
    dim <- as.integer(dim)
    if (length(dim) == 0L) return(if (length(data)) data[[1L]] else NA_real_)
    if (length(dim) == 1L) return(array(data, dim = dim))
    array(aperm(array(data, dim = rev(dim)), rev(seq_along(dim))), dim = dim)
}

#' The number of live value cells of a session
#'
#' A diagnostic for leak hunting: it counts the cells the session still owns.
#'
#' @param session A `rssn_session`.
#' @return An integer.
#' @export
value_count <- function(session = NULL) {
    session <- session %||% default_session()
    .Call(C_rssn_value_count, session$ptr)
}

# ---------------------------------------------------------------------------
# R objects to cells
# ---------------------------------------------------------------------------

#' A value cell for an R object
#'
#' `x` as a cell of `session`: a `rssn_value` (returned as is), a `rssn_term`
#' (copied as a term value), a number, a string, a logical, a complex number, a
#' list (a record when it is named, a list value otherwise), a vector (a list of
#' its elements, which keeps each element's kind) or an array (a dense tensor).
#'
#' @param x The object.
#' @param session The session that owns the cell.
#' @return A `rssn_value`.
#' @examples
#' s <- rssn_session()
#' value_kind(value_new(c(1, 2.5, 3), s))   # "tensor"
#' @export
value_new <- function(x, session = NULL) {
    session <- session %||% default_session()
    .rssn_value(session, .rssn_cell(x, session))
}

#' @rdname value_new
#' @export
value_of <- value_new

.rssn_cell <- function(x, session, numeric = "auto") {
    if (inherits(x, "rssn_value")) {
        if (!identical(x$session$ptr, session$ptr)) {
            stop("the value belongs to another session")
        }
        return(value_cell(x))
    }
    if (inherits(x, "rssn_term")) {
        if (!identical(.rssn_session_of(x)$ptr, session$ptr)) {
            stop("the term belongs to another session")
        }
        return(.Call(C_rssn_value_new_term, session$ptr, handle(x)))
    }
    if (inherits(x, "rssn_callback")) {
        if (!identical(x$session$ptr, session$ptr)) {
            stop("the callback belongs to another session")
        }
        return(callback_cell(x))
    }
    if (is.function(x)) {
        stop("a function is a value only as an operator argument: pass it to call_value()")
    }
    if (is.logical(x)) {
        if (length(x) == 1L) {
            if (is.na(x)) stop("NA has no value form")
            return(.Call(C_rssn_value_new_bool, session$ptr, x))
        }
        return(.rssn_cell_list(x, session))
    }
    if (is.numeric(x)) {
        if (length(x) == 1L) return(.rssn_cell_number(x, session, numeric))
        if (length(dim(x)) >= 2L) return(.rssn_cell_tensor(x, session))
        return(.rssn_cell_list(x, session))
    }
    if (is.complex(x)) {
        if (length(x) == 1L) {
            return(.Call(C_rssn_value_new_complex, session$ptr, Re(x), Im(x)))
        }
        return(.rssn_cell_list(x, session))
    }
    if (is.character(x)) {
        if (length(x) == 1L) {
            if (is.na(x)) stop("NA has no value form")
            return(.Call(C_rssn_value_new_text, session$ptr, x))
        }
        return(.rssn_cell_list(x, session))
    }
    if (is.list(x)) {
        names <- names(x)
        if (!is.null(names) && any(nzchar(names))) {
            if (any(!nzchar(names))) stop("a record needs a name for every field")
            cells <- vapply(x, .rssn_cell, integer(1), session = session)
            return(.Call(C_rssn_value_new_record, session$ptr, names, cells))
        }
        return(.rssn_cell_list(x, session))
    }
    stop(sprintf("cannot use an object of class '%s' as a value",
                 paste(class(x), collapse = "/")))
}

# R has no integer literal: `1` is a double and `1L` is an integer, and the
# engine's value ABI is strict about the difference (an operator that wants a
# `Float` refuses an `Int`, and the other way round). The R type decides by
# default; `numeric` overrides it for one argument, which is what the retry in
# `.rssn_value_call()` uses when an operator's parameter turns out to be the
# other kind.
.rssn_cell_number <- function(x, session, numeric = "auto") {
    if (is.na(x)) stop("NA has no value form")
    if (identical(numeric, "int")) {
        if (x != round(x)) stop("the operator wants an integer and this is not one")
        return(.Call(C_rssn_value_new_int, session$ptr, as.numeric(x)))
    }
    if (identical(numeric, "float")) {
        return(.Call(C_rssn_value_new_float, session$ptr, as.numeric(x)))
    }
    if (is.integer(x)) {
        return(.Call(C_rssn_value_new_int, session$ptr, as.numeric(x)))
    }
    .Call(C_rssn_value_new_float, session$ptr, as.numeric(x))
}

.rssn_cell_list <- function(x, session) {
    cells <- vapply(seq_along(x), function(index) .rssn_cell(x[[index]], session), integer(1))
    .Call(C_rssn_value_new_list, session$ptr, cells)
}

.rssn_cell_tensor <- function(x, session) {
    extents <- dim(x)
    rank <- length(extents)
    data <- if (rank == 1L) as.numeric(x) else as.numeric(aperm(x, rev(seq_len(rank))))
    .Call(C_rssn_value_new_tensor, session$ptr, data, as.integer(extents))
}

# ---------------------------------------------------------------------------
# cells to R objects (what a callback receives)
# ---------------------------------------------------------------------------

.rssn_to_r <- function(session, cell, owned = FALSE) {
    value <- .rssn_value(session, cell, owned = owned)
    kind <- value_kind(value)
    if (kind == "int") return(as.integer(round(value_float(value))))
    if (kind == "real" || kind == "float") return(value_float(value))
    if (kind == "bool") return(isTRUE(value_json(value)[["bool"]]))
    if (kind == "complex") {
        body <- value_json(value)[["complex"]]
        if (!is.null(body[["Exact"]])) {
            return(complex(real = .rssn_number(body[["Exact"]][["re"]]),
                           imaginary = .rssn_number(body[["Exact"]][["im"]])))
        }
        floating <- body[["Float"]]
        return(complex(real = as.numeric(floating[["re"]]),
                       imaginary = as.numeric(floating[["im"]])))
    }
    if (kind == "string" || kind == "symbol") return(value_text(value))
    if (kind == "tensor") {
        return(tryCatch(value_tensor(value), error = function(condition) {
            lapply(seq_len(value_len(value)), function(index) {
                .rssn_to_r(session, value[[index]]$cell, owned = TRUE)
            })
        }))
    }
    if (kind == "record") {
        fields <- list()
        for (index in seq_len(value_len(value))) {
            name <- .Call(C_rssn_value_field_name, session$ptr, cell, index - 1L)
            fields[[name]] <- .rssn_to_r(session, value[[index]]$cell, owned = TRUE)
        }
        return(fields)
    }
    if (kind == "variant") {
        tag <- .Call(C_rssn_value_tag, session$ptr, cell)
        payload <- lapply(seq_len(value_len(value)), function(index) {
            .rssn_to_r(session, value[[index]]$cell, owned = TRUE)
        })
        return(list(tag = tag, value = payload))
    }
    value_json(value)
}

# A `Number` JSON (`list(Float = 1)`, `list(Int = 3)`, `list(Rat = c(1, 3))`).
.rssn_number <- function(payload) {
    if (is.numeric(payload)) return(as.numeric(payload))
    if (is.list(payload) && length(payload) == 1L) {
        tag <- names(payload)[[1L]]
        body <- payload[[1L]]
        if (tag %in% c("Float", "float", "Int", "int")) return(as.numeric(body))
        if (tag %in% c("Rat", "rat")) {
            if (is.list(body)) return(body$numer / body$denom)
            return(body[[1L]] / body[[2L]])
        }
    }
    stop("not a number")
}

# ---------------------------------------------------------------------------
# calling operators through the value ABI
# ---------------------------------------------------------------------------

#' Call an operator through the value ABI
#'
#' This is the path to the operators whose arguments have no term form: native
#' object handles, records, host callbacks. An ordinary R function in the
#' argument list becomes a host callback, with the signature the operator
#' declares for that position.
#'
#' @param op The operator name.
#' @param ... The arguments: values, terms, numbers, strings, logicals, lists,
#'   arrays, and functions.
#' @param session The session to call it in; taken from the first value or term
#'   argument when omitted.
#' @return A `rssn_value`.
#' @examples
#' s <- rssn_session()
#' answer <- call_value("num_brent", 1, 2, 1e-12, 64, function(x) x * x - 2,
#'                      session = s)
#' value_float(answer[[1]])   # 1.4142135623730951
#' @export
call_value <- function(op, ..., session = NULL) {
    .rssn_call(op, list(...), session)
}

# The entry point the generated facade uses for a higher-order operator: the
# same call, with the arguments already collected in a list.
.rssn_call <- function(op, arguments, session) {
    if (is.null(session)) {
        for (argument in arguments) {
            if (inherits(argument, "rssn_value")) {
                session <- argument$session
                break
            }
            if (inherits(argument, "rssn_term")) {
                session <- .rssn_session_of(argument)
                break
            }
        }
        if (is.null(session)) session <- default_session()
    }
    .rssn_value(session, .rssn_value_call(op, arguments, session))
}

# Calls the operator, registering every function argument as a callback and
# retrying when the engine reports that an argument's numeric kind is the other
# one.
#
# R has no integer literal, so `64` reaches the engine as a `Float` while the
# parameter may be an `Int` (and `1L` the other way round). Rather than make
# every caller learn which parameters are integers, the first failure that names
# an argument and a kind flips that one argument and tries again -- at most a
# handful of times, and only for a mismatch the engine itself described. The
# message is read for this one recovery, never to decide what failed: the status
# code is what `rssn_error` carries.
.rssn_value_call <- function(op, arguments, session) {
    kinds <- rep("auto", length(arguments))
    outcome <- NULL
    for (attempt in seq_len(6L)) {
        outcome <- tryCatch(
            list(value = .rssn_value_call_once(op, arguments, session, kinds)),
            error = function(condition) list(error = condition)
        )
        if (!is.null(outcome$value)) return(outcome$value)
        flipped <- .rssn_flip_kind(conditionMessage(outcome$error), kinds)
        if (is.null(flipped)) break
        kinds <- flipped
    }
    stop(outcome$error)
}

# The per-argument numeric kind the engine asked for, or `NULL` when the message
# is not about one (`argument 3 has type Float, expected Int`).
.rssn_flip_kind <- function(message, kinds) {
    match <- regmatches(
        message,
        regexec("argument ([0-9]+) has type [A-Za-z0-9_<>]+, expected (Int|Float)", message)
    )[[1L]]
    if (length(match) < 3L) return(NULL)
    index <- as.integer(match[[2L]]) + 1L
    if (index < 1L || index > length(kinds)) return(NULL)
    wanted <- if (match[[3L]] == "Int") "int" else "float"
    if (identical(kinds[[index]], wanted)) return(NULL)
    kinds[[index]] <- wanted
    kinds
}

# One attempt: the callbacks registered, the other arguments converted, one call.
.rssn_value_call_once <- function(op, arguments, session, kinds) {
    positions <- which(vapply(arguments, is.function, logical(1)))
    if (length(positions) == 0L) {
        cells <- integer(length(arguments))
        for (index in seq_along(arguments)) {
            cells[[index]] <- .rssn_cell(arguments[[index]], session, kinds[[index]])
        }
        return(.Call(C_rssn_value_call, session$ptr, op, cells))
    }
    options <- lapply(positions, function(position) {
        signatures <- callback_signatures(op, position, session = session)
        if (length(signatures) == 0L) list(NULL) else signatures
    })
    combinations <- expand.grid(lapply(options, seq_along), KEEP.OUT.ATTRS = FALSE,
                                stringsAsFactors = FALSE)
    first_error <- NULL
    for (row in seq_len(min(nrow(combinations), 8L))) {
        callbacks <- list()
        result <- tryCatch({
            cells <- integer(length(arguments))
            for (index in seq_along(arguments)) {
                slot <- match(index, positions)
                if (is.na(slot)) {
                    cells[[index]] <- .rssn_cell(arguments[[index]], session, kinds[[index]])
                } else {
                    declared <- options[[slot]][[combinations[row, slot]]]
                    callback <- if (is.null(declared)) {
                        rssn_callback(arguments[[index]], session = session, op = op,
                                      index = index)
                    } else {
                        rssn_callback(arguments[[index]], session = session,
                                      signature = declared)
                    }
                    callbacks[[length(callbacks) + 1L]] <- callback
                    cells[[index]] <- callback_cell(callback)
                }
            }
            .Call(C_rssn_value_call, session$ptr, op, cells)
        }, error = function(condition) {
            if (is.null(first_error)) first_error <<- condition
            NULL
        })
        for (callback in callbacks) close(callback)
        if (!is.null(result)) return(result)
    }
    if (is.null(first_error)) stop("no callback signature to try for this operator")
    stop(first_error)
}
