library(rssn)

check <- function(label, condition) {
    if (!isTRUE(condition)) stop("FAILED: ", label, call. = FALSE)
}

check("empty object", identical(json_decode("{}"), structure(list(), names = character())))
check("empty array", identical(json_decode("[]"), list()))
check("integer", identical(json_decode("1"), 1L))
check("double", identical(json_decode("1.5"), 1.5))
check("exponent", identical(json_decode("-2e-3"), -0.002))
check("true", isTRUE(json_decode("true")))
check("null", is.null(json_decode("null")))
check("newline escape", identical(json_decode("\"a\\nb\""), "a\nb"))
check("quote escape", identical(json_decode("\"q\\\"q\""), "q\"q"))
check("backslash escape", identical(json_decode("\"back\\\\slash\""), "back\\slash"))
check("unicode escape", identical(json_decode("\"\\u00e9\""), "\u00e9"))
check("nesting", identical(json_decode("{\"a\": [1, 2.5, \"x\"]}")$a[[3]], "x"))
check("null member", is.null(json_decode("{\"a\": null}")$a))
check("encoding", identical(json_encode(list(a = c(1, 2))), "{\"a\":[1,2]}"))
check("round trip", identical(json_decode(json_encode(list(a = 1, b = "x", c = TRUE))),
                              list(a = 1L, b = "x", c = TRUE)))
check("nested arrays", identical(json_decode("[[1,2],[3]]"), list(list(1L, 2L), list(3L))))

# A large document: 5000 records, which exercises the growable buffers.
record <- "{\"n\": 1, \"v\": [1.5, null, \"x\"], \"ok\": true}"
document <- paste0("[", paste(rep(record, 5000L), collapse = ","), "]")
start <- Sys.time()
decoded <- json_decode(document)
elapsed <- as.numeric(Sys.time() - start, units = "secs")
cat(sprintf("5000 records: %.2fs
", elapsed))
check("large array length", length(decoded) == 5000L)
check("large array member", identical(decoded[[5000L]]$v[[1L]], 1.5))
check("large array null", is.null(decoded[[1L]]$v[[2L]]))
check("large array string", identical(decoded[[1L]]$v[[3L]], "x"))
check("round trip of the document", identical(json_decode(json_encode(decoded)), decoded))
cat("all json checks pass
")
