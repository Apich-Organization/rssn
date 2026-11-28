# Build an exact rational matrix from row-major cells.
library(rssn)

session <- rssn_session()
foreign <- rssn_session()
verify <- function() {
    on.exit(close(session))
    on.exit(close(foreign), add = TRUE)
    wide <- value_integer_decimal("18446744073709551617", session)
    fraction <- call_value("divide", 2L, 3L, session = session)
    matrix <- value_tensor_cells(list(wide, fraction, 3L, 4L), c(2L, 2L), session)
    close(wide)
    close(fraction)
    a <- call_value("num_native_qmatrix_new", matrix, 1L, session = session)
    determinant <- call_value("num_native_qlinalg_det", a, session = session)
    stopifnot(value_text(determinant) == "73786976294838206466")
    cat("Exact determinant:", value_text(determinant), "\n")

    scalar <- value_tensor_cells(list(7L), integer(), session)
    empty <- value_tensor_cells(list(), c(2L, 0L, 3L), session, exemplar = 0L)
    stopifnot(value_text(scalar[[1L]]) == "7", value_len(empty) == 0L)
    stopifnot(value_json(empty)$tensor$elem == "int")

    fails <- function(expr, text) {
        tryCatch({ force(expr); FALSE },
                 error = function(e) grepl(text, conditionMessage(e), fixed = TRUE))
    }
    stopifnot(fails(value_tensor_cells(list(1L), c(2L, 2L), session), "does not match"))
    stopifnot(fails(value_tensor_cells(list(list(1L), list(2L, 3L)), 2L, session),
                    "flat row-major"))
    alien <- value_integer_decimal("7", foreign)
    stopifnot(fails(value_tensor_cells(list(alien), 1L, session), "another session"))
}

verify()
