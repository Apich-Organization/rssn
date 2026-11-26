# Host callbacks: an R function the engine calls.
#
# 189 operators are higher-order: they take a function and call it once per
# quadrature node, optimisation step or root-finder probe (`num_bfgs`,
# `num_brent`, `num_analysis_integrate_adaptive_scaled`, `num_bootstrap`, ...).
# In R any function works, because this file turns one into the C function
# pointer the engine knows how to call.
#
#     s <- rssn_session()
#     answer <- call_value("num_brent", 1, 2, 1e-12, 64,
#                          function(x) x * x - 2, session = s)
#     value_float(answer[[1]])        # 1.4142135623730951
#
# The function receives native R objects -- numbers, integers, logicals,
# complex numbers, strings, arrays for tensors and named lists for records --
# and whatever it returns is converted back into a value. An explicit
# `rssn_callback()` is for a function that goes to several calls, or one whose
# type has to be spelled out.
#
# The engine calls the function on the thread that made the call, so the R API
# is used from the thread it belongs to. An error inside the function is caught,
# reported to the engine as `InvalidArgument`, and the enclosing call raises an
# `rssn_error` condition whose message is the R condition's.

#: The signatures the compact grammar accepts, by spelling.
.SIGNATURE_SCALARS <- c(
    "int" = "int", "integer" = "int",
    "real" = "real", "rational" = "real",
    "float" = "float", "double" = "float",
    "complex" = "complex",
    "bool" = "bool", "boolean" = "bool",
    "string" = "string", "str" = "string", "text" = "string",
    "symbol" = "symbol"
)

#' Convert a compact callback signature to the JSON the engine reads
#'
#' `"(float, float) -> float"` becomes `'{"params":["float","float"],"ret":"float"}'`.
#' The arrow may be `->` or `=>`; `[]` appends a tensor rank (`"float[]"` is a
#' vector, `"float[][]"` a matrix). A string that is already JSON, or a list
#' that already holds the object, is passed through.
#'
#' @param text A compact signature, a JSON string, or a list.
#' @return A length-one JSON string.
#' @examples
#' parse_signature("(float, float) -> float")
#' @export
parse_signature <- function(text) {
    if (is.list(text)) return(json_encode(text))
    text <- trimws(as.character(text))
    if (startsWith(text, "{") || startsWith(text, "[")) return(text)
    arrow <- if (grepl("->", text, fixed = TRUE)) "->" else "=>"
    if (!grepl(arrow, text, fixed = TRUE)) {
        stop(sprintf("a callback signature needs '%s': %s", arrow, text))
    }
    parts <- strsplit(text, arrow, fixed = TRUE)[[1L]]
    left <- trimws(parts[[1L]])
    right <- trimws(paste(parts[-1L], collapse = arrow))
    if (startsWith(left, "(") && endsWith(left, ")")) {
        left <- substr(left, 2L, nchar(left) - 1L)
    }
    parameters <- if (nzchar(trimws(left))) {
        lapply(.rssn_split_top(left), .rssn_type_expression)
    } else {
        list()
    }
    json_encode(list(params = parameters, ret = .rssn_type_expression(right)))
}

# Splits on commas that are not inside brackets.
.rssn_split_top <- function(text) {
    parts <- character()
    depth <- 0L
    start <- 1L
    characters <- strsplit(text, "")[[1L]]
    for (index in seq_along(characters)) {
        character <- characters[[index]]
        if (character %in% c("(", "[", "{")) depth <- depth + 1L
        if (character %in% c(")", "]", "}")) depth <- depth - 1L
        if (character == "," && depth == 0L) {
            parts <- c(parts, trimws(substr(text, start, index - 1L)))
            start <- index + 1L
        }
    }
    tail <- trimws(substr(text, start, nchar(text)))
    if (nzchar(tail)) parts <- c(parts, tail)
    parts
}

# One compact type (`real`, `float[]`) as the JSON `TypeExpr` shape.
.rssn_type_expression <- function(token) {
    name <- trimws(token)
    rank <- 0L
    while (endsWith(name, "[]")) {
        rank <- rank + 1L
        name <- trimws(substr(name, 1L, nchar(name) - 2L))
    }
    base <- .SIGNATURE_SCALARS[[tolower(name)]]
    if (is.null(base)) {
        stop(sprintf("unknown type '%s' in a callback signature; use one of %s, or append []",
                     name, paste(unique(.SIGNATURE_SCALARS), collapse = ", ")))
    }
    if (rank == 0L) return(base)
    list(tensor = list(rank = list(lit = rank), shape = NULL, elem = base))
}

#' A host callback: an R function registered with a session
#'
#' Create one when the same function goes to several calls, or when its type has
#' to be spelled out. Give either `signature` -- a compact `"(float) -> float"`,
#' a JSON string, or the list form -- or the operator and the argument position
#' the callback fills, in which case the engine supplies the type from the
#' operator's own declaration.
#'
#' @param function_ The function. It receives one R object per argument.
#' @param session The session the callback belongs to.
#' @param signature The function's type, when it is not taken from an operator.
#' @param op The operator the callback is for.
#' @param index The argument position of `op` the callback fills, counting from
#'   one.
#' @return An object of class `rssn_callback`.
#' @examples
#' s <- rssn_session()
#' f <- rssn_callback(function(x) x * x - 2, s, signature = "(float) -> float")
#' answer <- call_value("num_brent", 1, 2, 1e-12, 64, f, session = s)
#' close(f)
#' @export
rssn_callback <- function(function_, session = NULL, signature = NULL, op = NULL,
                          index = 1L) {
    if (!is.function(function_)) stop("a callback must be a function")
    session <- session %||% default_session()
    registered <- .rssn_register_callback(session, function_, signature, op, index)
    callback <- new.env(parent = emptyenv())
    callback$session <- session
    callback$cell <- as.integer(registered[[1L]])
    callback$registry <- as.integer(registered[[2L]])
    callback$function_ <- function_
    callback$signature <- signature
    callback$op <- op
    callback$index <- index
    class(callback) <- "rssn_callback"
    reg.finalizer(callback, .rssn_finalize_callback, onexit = TRUE)
    callback
}

.rssn_register_callback <- function(session, function_, signature, op, index) {
    # The closure the C trampoline calls: it receives the argument cell ids and
    # answers one cell id, or a negative status after reporting the message.
    closure <- function(cells) {
        tryCatch({
            arguments <- lapply(as.integer(cells), function(cell) {
                .rssn_to_r(session, cell, owned = FALSE)
            })
            result <- do.call(function_, arguments)
            as.integer(.rssn_cell(result, session))
        }, error = function(condition) {
            message <- paste0(class(condition)[[1L]], ": ", conditionMessage(condition))
            try(.Call(C_rssn_callback_fail, message), silent = TRUE)
            -1L
        })
    }
    if (!is.null(signature)) {
        .Call(C_rssn_callback_new, session$ptr, parse_signature(signature), closure)
    } else if (!is.null(op)) {
        .Call(C_rssn_callback_new_for, session$ptr, as.character(op), as.integer(index) - 1L,
              closure)
    } else {
        stop("a callback needs a signature, or an operator and an argument position")
    }
}

.rssn_finalize_callback <- function(callback) {
    if (callback$cell != 0L) {
        session <- callback$session
        if (!is.null(session) && !is.null(session$ptr)) {
            try(.Call(C_rssn_callback_free, session$ptr, callback$cell), silent = TRUE)
        }
        callback$cell <- 0L
    }
    if (callback$registry != 0L) {
        try(.Call(C_rssn_registry_release, as.numeric(callback$registry)), silent = TRUE)
        callback$registry <- 0L
    }
    invisible(NULL)
}

#' The cell id of a callback
#' @param callback A `rssn_callback`.
#' @return An integer cell id.
#' @export
callback_cell <- function(callback) {
    if (callback$cell == 0L) stop("this callback has been released")
    callback$cell
}

#' @export
print.rssn_callback <- function(x, ...) {
    cat("<rssn callback> ", if (x$cell == 0L) "released" else paste0("cell ", x$cell), "\n",
        sep = "")
    invisible(x)
}

#' Retire a callback
#'
#' Releases the callback's cell and the registry slot that keeps the R function
#' alive. Safe to call twice.
#'
#' @param con A `rssn_callback`.
#' @param ... Ignored.
#' @return `NULL`, invisibly.
#' @export
close.rssn_callback <- function(con, ...) {
    .rssn_finalize_callback(con)
    invisible(NULL)
}

#' The signatures an operator declares for a callback at a position
#'
#' Usually one; an operator with an objective *and* a gradient declares two
#' different function types at different positions, and a few declare more than
#' one at the same position. Positions count from one, as in R.
#'
#' @param op The operator name.
#' @param index The argument position, counting from one.
#' @param session The session to ask.
#' @return A list of signature objects, each of which `rssn_callback()` accepts.
#' @examples
#' callback_signatures("num_brent", 5)
#' @export
callback_signatures <- function(op, index, session = NULL) {
    session <- session %||% default_session()
    if (index < 1L) stop("argument positions count from one")
    text <- .Call(C_rssn_callback_signatures, session$ptr, as.character(op),
                  as.integer(index) - 1L)
    decoded <- json_decode(text)
    if (is.list(decoded) && length(decoded) == 0L) return(list())
    if (is.list(decoded) && is.null(names(decoded))) return(decoded)
    list(decoded)
}

#' The number of live host callbacks of a session
#' @param session A `rssn_session`.
#' @return An integer; a diagnostic for leak hunting.
#' @export
callback_count <- function(session = NULL) {
    session <- session %||% default_session()
    .Call(C_rssn_callback_count, session$ptr)
}
