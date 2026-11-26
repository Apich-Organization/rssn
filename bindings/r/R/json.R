# A small JSON codec.
#
# The engine hands JSON to the host in four places: the compute record
# (`compute_json`), a value's `serde` form (`value_json`), the manifest
# (`manifest`) and a simulation's result (`sim_run`). R has no JSON parser in
# base, and this package deliberately has no dependencies, so the codec lives
# here -- about the same size as the one the Julia binding carries.
#
# `json_decode` answers the natural R object: an object becomes a named list, an
# array an unnamed list, a number a double (an integer when it has no fraction),
# `null` `NULL`. `json_encode` is the inverse and is what `sim_run` uses to send
# parameters.
#
# The decoder is written for the manifest, which is about a megabyte: a string is
# scanned with one anchored regular expression rather than character by
# character, and the containers grow by doubling rather than by appending.

#' Decode a JSON string
#'
#' @param text A length-one character vector.
#' @return A named list (object), an unnamed list (array), a double, an integer,
#'   a string, a logical, or `NULL`.
#' @examples
#' json_decode("{\"a\": [1, 2.5, null]}")
#' @export
json_decode <- function(text) {
    if (length(text) != 1L || is.na(text)) {
        stop("json_decode() takes one string")
    }
    state <- new.env(parent = emptyenv())
    state$text <- text
    state$position <- 1L
    state$length <- nchar(text, type = "chars")
    value <- .json_value(state)
    .json_skip_space(state)
    if (state$position <= state$length) {
        stop(sprintf("trailing text in JSON at character %d", state$position))
    }
    value
}

#' Encode an R object as a JSON string
#'
#' @param value A list, vector, or scalar. Named lists become objects, unnamed
#'   lists and vectors become arrays, `NULL` becomes `null`.
#' @param auto_unbox Whether a length-one vector becomes a scalar (`TRUE`, the
#'   default, which is what the engine's parameter objects want).
#' @return A length-one character vector.
#' @examples
#' json_encode(list(a = c(1, 2), b = TRUE))
#' @export
json_encode <- function(value, auto_unbox = TRUE) {
    paste0(.json_write(value, auto_unbox), collapse = "")
}

# -- decoding ---------------------------------------------------------------

.json_peek <- function(state) {
    if (state$position > state$length) return("")
    substr(state$text, state$position, state$position)
}

.json_take <- function(state) {
    character <- .json_peek(state)
    state$position <- state$position + 1L
    character
}

.json_skip_space <- function(state) {
    while (state$position <= state$length) {
        character <- .json_peek(state)
        if (character != " " && character != "\t" && character != "\n" && character != "\r") {
            break
        }
        state$position <- state$position + 1L
    }
    invisible(NULL)
}

.json_expect <- function(state, expected) {
    .json_skip_space(state)
    character <- .json_take(state)
    if (character != expected) {
        stop(sprintf("expected '%s' in JSON at character %d", expected, state$position))
    }
    invisible(NULL)
}

.json_value <- function(state) {
    .json_skip_space(state)
    character <- .json_peek(state)
    if (character == "") stop("unexpected end of JSON")
    if (character == "{") return(.json_object(state))
    if (character == "[") return(.json_array(state))
    if (character == "\"") return(.json_string(state))
    if (character == "t" || character == "f" || character == "n") return(.json_word(state))
    .json_number(state)
}

# A growable list, so an array or an object of n members costs n appends rather
# than n copies.
.json_buffer <- function() {
    buffer <- new.env(parent = emptyenv())
    buffer$items <- vector("list", 32L)
    buffer$count <- 0L
    buffer
}

.json_push <- function(buffer, value) {
    buffer$count <- buffer$count + 1L
    if (buffer$count > length(buffer$items)) {
        length(buffer$items) <- length(buffer$items) * 2L
    }
    # Single-bracket assignment with a one-element list, so a `NULL` value is
    # stored rather than deleting the slot (`[[<-` with `NULL` would delete it).
    buffer$items[buffer$count] <- list(value)
    invisible(NULL)
}

.json_collect <- function(buffer) {
    if (buffer$count == 0L) return(list())
    buffer$items[seq_len(buffer$count)]
}

.json_object <- function(state) {
    .json_expect(state, "{")
    values <- .json_buffer()
    keys <- .json_buffer()
    .json_skip_space(state)
    if (.json_peek(state) == "}") {
        state$position <- state$position + 1L
        return(structure(list(), names = character()))
    }
    repeat {
        .json_skip_space(state)
        key <- .json_string(state)
        .json_expect(state, ":")
        .json_push(keys, key)
        .json_push(values, .json_value(state))
        .json_skip_space(state)
        separator <- .json_take(state)
        if (separator == "}") break
        if (separator != ",") {
            stop(sprintf("expected ',' or '}' in JSON at character %d", state$position))
        }
    }
    out <- .json_collect(values)
    names(out) <- unlist(.json_collect(keys), use.names = FALSE)
    out
}

.json_array <- function(state) {
    .json_expect(state, "[")
    values <- .json_buffer()
    .json_skip_space(state)
    if (.json_peek(state) == "]") {
        state$position <- state$position + 1L
        return(list())
    }
    repeat {
        .json_push(values, .json_value(state))
        .json_skip_space(state)
        separator <- .json_take(state)
        if (separator == "]") break
        if (separator != ",") {
            stop(sprintf("expected ',' or ']' in JSON at character %d", state$position))
        }
    }
    .json_collect(values)
}

# A JSON string, scanned with one anchored regular expression: a string is a
# quote, any run of characters that are neither a quote nor a backslash, and an
# escape (a backslash and one more character) allowed anywhere.
.JSON_STRING_PATTERN <- '^"(?:[^"\\\\]|\\\\.)*"'

# Strings are matched inside a window that grows when a token runs past its end,
# so a megabyte of document costs one pass rather than one copy of the remaining
# text per token.
.JSON_WINDOW <- 4096L

.json_string <- function(state) {
    .json_skip_space(state)
    start <- state$position
    size <- .JSON_WINDOW
    repeat {
        end <- min(start + size - 1L, state$length)
        window <- substr(state$text, start, end)
        match <- regexpr(.JSON_STRING_PATTERN, window, perl = TRUE)
        if (match > 0L) {
            width <- attr(match, "match.length")
            if (start + width - 1L < end || end >= state$length) {
                token <- substr(window, 2L, width - 1L)
                state$position <- start + width
                return(if (grepl("\\", token, fixed = TRUE)) .json_unescape(token) else token)
            }
        }
        if (end >= state$length) {
            stop(sprintf("unterminated string in JSON at character %d", start))
        }
        size <- size * 2L
    }
}

.JSON_ESCAPES <- c("\\" = "\\", "\"" = "\"", "/" = "/", "b" = "\b", "f" = "\f",
                   "n" = "\n", "r" = "\r", "t" = "\t")

.json_unescape <- function(token) {
    # The simple escapes first, then `\uXXXX`, which is rare and needs a
    # conversion per occurrence.
    for (escape in names(.JSON_ESCAPES)) {
        token <- gsub(paste0("\\", escape), .JSON_ESCAPES[[escape]], token, fixed = TRUE)
    }
    repeat {
        position <- regexpr("[\\\\]u[0-9a-fA-F]{4}", token)
        if (position < 0L) break
        code <- strtoi(substr(token, position + 2L, position + 5L), 16L)
        token <- paste0(substr(token, 1L, position - 1L), intToUtf8(code),
                        substr(token, position + 6L, nchar(token)))
    }
    token
}

.json_word <- function(state) {
    for (word in list(list("true", TRUE), list("false", FALSE), list("null", NULL))) {
        width <- nchar(word[[1]])
        if (substr(state$text, state$position, state$position + width - 1L) == word[[1]]) {
            state$position <- state$position + width
            return(word[[2]])
        }
    }
    stop(sprintf("unexpected token in JSON at character %d", state$position))
}

.JSON_NUMBER_PATTERN <- "^-?(?:0|[1-9][0-9]*)(?:\\.[0-9]+)?(?:[eE][-+]?[0-9]+)?"

.json_number <- function(state) {
    .json_skip_space(state)
    start <- state$position
    end <- min(start + 63L, state$length)
    window <- substr(state$text, start, end)
    match <- regexpr(.JSON_NUMBER_PATTERN, window, perl = TRUE)
    if (match < 0L) stop(sprintf("expected a value in JSON at character %d", start))
    width <- attr(match, "match.length")
    text <- substr(window, 1L, width)
    state$position <- start + width
    if (!grepl("[.eE]", text)) {
        number <- suppressWarnings(as.integer(text))
        if (!is.na(number)) return(number)
    }
    as.numeric(text)
}

# -- encoding ---------------------------------------------------------------

.json_write <- function(value, auto_unbox) {
    if (is.null(value)) return("null")
    if (is.function(value)) stop("a function cannot be encoded as JSON")
    if (inherits(value, "rssn_term") || inherits(value, "rssn_value")) {
        return(.json_write(as.character(value), auto_unbox))
    }
    if (is.logical(value)) {
        if (length(value) == 1L && auto_unbox) {
            return(if (is.na(value)) "null" else if (value) "true" else "false")
        }
        return(.json_array_of(value, auto_unbox))
    }
    if (is.numeric(value) || is.integer(value)) {
        if (length(value) == 1L && auto_unbox) return(.json_number_text(value))
        return(.json_array_of(value, auto_unbox))
    }
    if (is.complex(value)) {
        stop("a complex number has no JSON form; encode it as a list(re = , im = )")
    }
    if (is.character(value)) {
        if (length(value) == 1L && auto_unbox) {
            return(if (is.na(value)) "null" else .json_quote(value))
        }
        return(.json_array_of(value, auto_unbox))
    }
    if (is.list(value)) {
        names <- names(value)
        if (!is.null(names) && any(nzchar(names))) {
            fields <- vapply(seq_along(value), function(index) {
                paste0(.json_quote(names[[index]]), ":", .json_write(value[[index]], auto_unbox))
            }, character(1))
            return(paste0("{", paste0(fields, collapse = ","), "}"))
        }
        return(.json_array_of(value, auto_unbox))
    }
    stop(sprintf("cannot encode an object of class '%s' as JSON",
                 paste(class(value), collapse = "/")))
}

.json_array_of <- function(value, auto_unbox) {
    items <- vapply(seq_along(value), function(index) {
        paste0(.json_write(value[[index]], auto_unbox), collapse = "")
    }, character(1))
    paste0("[", paste0(items, collapse = ","), "]")
}

.json_number_text <- function(value) {
    if (is.na(value)) return("null")
    if (is.integer(value)) return(as.character(value))
    if (!is.finite(value)) return("null")
    if (value == round(value) && abs(value) < 1e15) {
        return(format(value, scientific = FALSE, trim = TRUE))
    }
    sprintf("%.17g", value)
}

.json_quote <- function(text) {
    text <- gsub("\\", "\\\\", text, fixed = TRUE)
    text <- gsub("\"", "\\\"", text, fixed = TRUE)
    text <- gsub("\n", "\\n", text, fixed = TRUE)
    text <- gsub("\r", "\\r", text, fixed = TRUE)
    text <- gsub("\t", "\\t", text, fixed = TRUE)
    text <- gsub("\b", "\\b", text, fixed = TRUE)
    text <- gsub("\f", "\\f", text, fixed = TRUE)
    paste0("\"", text, "\"")
}
