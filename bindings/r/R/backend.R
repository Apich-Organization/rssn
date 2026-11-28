#' Select an RSSN execution backend
#'
#' The active backend is process-wide. Use an explicit backend handle when
#' operating on a program through the native value API.
#'
#' @param session The RSSN session. Defaults to the package's default session.
#' @return `rssn_available_backends()` returns a character vector of usable
#'   selectors. `rssn_current_backend()` and `rssn_set_default_backend()` return
#'   the effective backend name. `rssn_open_backend()` returns an `rssn_value`.
#'   `rssn_release_backend()` returns `NULL` invisibly.
#' @export
rssn_available_backends <- function(session = NULL) {
    choices <- call_value("num_backend_available", session = session)
    on.exit(close(choices))
    vapply(seq_len(length(choices)), function(index) {
        choice <- choices[[index]]
        on.exit(close(choice))
        value_text(choice)
    }, character(1L))
}

#' @rdname rssn_available_backends
#' @export
rssn_current_backend <- function(session = NULL) {
    selected <- tryCatch(
        call_value("num_backend_current_name", session = session),
        rssn_error = function(error) {
            if (!identical(error$status_name, "unknown_operator")) stop(error)
            NULL
        }
    )
    if (!is.null(selected)) {
        on.exit(close(selected))
        return(value_text(selected))
    }

    backend <- call_value("num_backend_current", session = session)
    on.exit(tryCatch(rssn_release_backend(backend), finally = close(backend)))
    name <- call_value("num_backend_name", backend)
    on.exit(close(name), add = TRUE, after = FALSE)
    value_text(name)
}

#' @rdname rssn_available_backends
#' @param name A selector from `rssn_available_backends()`. Use `"auto"` to
#'   restore the native library's automatic policy.
#' @export
rssn_open_backend <- function(name, session = NULL) {
    .rssn_check_backend_name(name)
    call_value("num_backend_named", name, session = session)
}

#' @rdname rssn_available_backends
#' @param backend A native handle from `rssn_open_backend()`. Release makes
#'   every copy of the handle invalid.
#' @export
rssn_release_backend <- function(backend) {
    result <- call_value("num_backend_release", backend)
    close(result)
    invisible(NULL)
}

#' @rdname rssn_available_backends
#' @export
rssn_set_default_backend <- function(name, session = NULL) {
    .rssn_check_backend_name(name)
    selected <- call_value("num_backend_select_default", name, session = session)
    on.exit(close(selected))
    value_text(selected)
}

.rssn_check_backend_name <- function(name) {
    if (!is.character(name) || length(name) != 1L || is.na(name)) {
        stop("supply one backend selector string")
    }
}
