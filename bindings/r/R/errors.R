# Conditions: the status vocabulary of `rssn.h`.
#
# Every failure carries the `RssnStatus` code the engine reported, as a
# condition of class `c("rssn_error", "error", "condition")` with a `status`
# element. The code is the contract; the message is for humans and the engine
# does not promise its wording.
#
#     tryCatch(parse_term("sin("), rssn_error = function(e) e$status)   # 3L
#     tryCatch(parse_term("sin("), rssn_error = function(e) e$status_name)  # "parse"

#: The status codes of `rssn.h`, in declaration order.
.STATUS_NAMES <- c(
    "ok", "null_argument", "invalid_utf8", "parse", "unknown_operator",
    "not_numeric", "inconsistent", "rule", "no_such_term", "no_value",
    "simulation", "panic", "dimension", "invalid_argument", "cancelled"
)

#' Build the condition for a status code without raising it
#'
#' @param status An integer status code.
#' @param message The message; the thread's last error when omitted.
#' @return A condition of class `rssn_error`.
#' @export
rssn_condition <- function(status, message = NULL) {
    status <- as.integer(status)
    if (is.null(message)) {
        message <- if (.rssn_state$loaded) .Call(C_rssn_last_error) else ""
    }
    if (!nzchar(message)) message <- "the rssn call failed"
    structure(
        class = c("rssn_error", "error", "condition"),
        list(
            message = message,
            call = NULL,
            status = status,
            status_name = if (status >= 0L && status < length(.STATUS_NAMES)) {
                .STATUS_NAMES[[status + 1L]]
            } else {
                "unknown"
            }
        )
    )
}

.rssn_condition <- rssn_condition
