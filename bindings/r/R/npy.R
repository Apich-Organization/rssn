#' Encode or decode NumPy array bytes
#'
#' The native codec reads Float32 and Float64 arrays and writes little-endian
#' Float64. File access is handled by the caller.
#'
#' @param x A numeric array or a flat row-major vector of samples.
#' @param shape Array extents. Defaults to `dim(x)` for an array and
#'   `length(x)` for a vector. Use `integer()` for a scalar.
#' @param session The RSSN session. Defaults to the package's default session.
#' @return `rssn_encode_npy()` returns a raw vector. `rssn_decode_npy()`
#'   returns an owned binary64 `rssn_value` tensor.
#' @export
rssn_encode_npy <- function(x, shape = dim(x) %||% length(x), session = NULL) {
    if (!is.numeric(x)) stop("supply numeric NumPy samples")
    if (!is.numeric(shape) || any(!is.finite(shape) | shape < 0 |
                                  shape != trunc(shape) | shape > .Machine$integer.max)) {
        stop("NumPy shape must contain nonnegative integers representable by R")
    }
    if (!is.null(dim(x)) && !identical(as.integer(shape), dim(x))) {
        stop("array extents must match the NumPy shape")
    }
    session <- session %||% default_session()
    extents <- as.integer(shape)
    dimensions <- value_tensor_cells(extents, length(extents), session,
                                     exemplar = if (!length(extents)) 0L)
    on.exit(close(dimensions))
    storage <- if (length(dim(x)) < 2L) x else aperm(x, rev(seq_along(dim(x))))
    samples <- as.numeric(storage)
    data <- .rssn_value(session, .Call(C_rssn_value_new_tensor, session$ptr,
                                      samples, length(samples)))
    on.exit(close(data), add = TRUE)
    encoded <- call_value("num_npy_encode", dimensions, data, session = session)
    on.exit(close(encoded), add = TRUE)
    .rssn_npy_unhex(value_text(encoded))
}

#' @rdname rssn_encode_npy
#' @param payload A raw vector containing NumPy array bytes.
#' @export
rssn_decode_npy <- function(payload, session = NULL) {
    if (!is.raw(payload)) stop("supply NumPy bytes as a raw vector")
    session <- session %||% default_session()
    hex <- paste(sprintf("%02x", as.integer(payload)), collapse = "")
    decoded <- call_value("num_npy_decode", hex, session = session)
    on.exit(close(decoded))
    dimensions <- decoded[["shape"]]
    on.exit(close(dimensions), add = TRUE)
    data <- decoded[["data"]]
    on.exit(close(data), add = TRUE)
    extents <- vapply(seq_len(length(dimensions)), function(index) {
        extent <- dimensions[[index]]
        on.exit(close(extent))
        result <- suppressWarnings(as.integer(value_text(extent)))
        if (is.na(result) || result < 0L) stop("NumPy extent exceeds R's array index range")
        result
    }, integer(1L))
    .rssn_value(session, .Call(C_rssn_value_new_tensor, session$ptr,
                               as.numeric(value_tensor(data)), extents))
}

.rssn_npy_unhex <- function(hex) {
    digits <- match(charToRaw(hex), charToRaw("0123456789abcdef")) - 1L
    if (!length(digits) || length(digits) %% 2L || anyNA(digits)) {
        stop("invalid RSSN NumPy byte payload")
    }
    high <- seq.int(1L, length(digits), by = 2L)
    as.raw(bitwShiftL(digits[high], 4L) + digits[high + 1L])
}
