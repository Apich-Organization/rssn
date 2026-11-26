#!/usr/bin/env Rscript
# Run the cross-language conformance table against the R binding.
#
#     R CMD INSTALL bindings/r
#     RSSN_LIB=/path/to/rssn.dll Rscript bindings/tests/run_conformance.R
#
# `cases.tsv` is the same file the Python, Julia and C++ runners use, so a
# wrapper that disagrees with the engine -- or with another wrapper -- fails
# here. The expectations are recorded by the Python runner
# (`run_conformance.py --record`).

suppressPackageStartupMessages(library(rssn))

arguments <- commandArgs(trailingOnly = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                       value = TRUE)[1])))
cases_path <- file.path(here, "cases.tsv")
if (length(arguments) && !startsWith(arguments[[1]], "--")) {
    cases_path <- arguments[[1]]
}

read_cases <- function(path) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    lines <- lines[!grepl("^\\s*#", lines) & nzchar(trimws(lines))]
    header <- strsplit(lines[[1]], "\t", fixed = TRUE)[[1]]
    rows <- lapply(lines[-1], function(line) {
        fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
        length(fields) <- length(header)
        fields[is.na(fields)] <- ""
        names(fields) <- header
        as.list(fields)
    })
    rows
}

literal <- function(text) {
    text <- trimws(text)
    if (grepl("[.eE]", text)) return(as.numeric(text))
    value <- suppressWarnings(as.integer(text))
    if (is.na(value)) as.numeric(text) else value
}

config_for <- function(row) {
    target <- trimws(row$target)
    bindings <- trimws(row$bindings)
    if (target != "numeric" && !nzchar(bindings)) return(NULL)
    config <- rssn_config()
    if (target == "numeric") {
        tolerance <- if (nzchar(trimws(row$tolerance))) as.numeric(row$tolerance) else 1e-12
        cfg_numeric(config, tolerance)
    }
    for (pair in strsplit(bindings, ";", fixed = TRUE)[[1]]) {
        if (!nzchar(trimws(pair))) next
        parts <- strsplit(pair, "=", fixed = TRUE)[[1]]
        if (length(parts) < 2) next
        cfg_bind(config, trimws(parts[[1]]), as.numeric(trimws(parts[[2]])))
    }
    config
}

status_of <- function(condition) {
    if (inherits(condition, "rssn_error")) condition$status_name else "exception"
}

# Runs one case; answers list(text=, value=, status=).
run_case <- function(row) {
    kind <- trimws(row$kind)
    payload <- row$payload
    if (kind == "manifest") {
        digest <- manifest_hash()
        stopifnot(nzchar(digest), digest != "0")
        return(list(text = "", value = digest, status = "ok"))
    }
    if (kind == "units") {
        return(list(text = "", value = as.character(length(UNITS)), status = "ok"))
    }
    if (kind == "constants") {
        return(list(text = "", value = as.character(length(CONSTANTS)), status = "ok"))
    }
    if (kind == "value") {
        # The second ABI: an operator called with native objects rather than
        # terms. The display form is the engine's own, so it is stable across
        # languages.
        parts <- strsplit(payload, " ", fixed = TRUE)[[1]]
        arguments <- if (length(parts) > 1) {
            lapply(strsplit(paste(parts[-1], collapse = " "), ",", fixed = TRUE)[[1]], literal)
        } else {
            list()
        }
        answer <- do.call(call_value, c(list(parts[[1]]), arguments))
        return(list(text = value_text(answer), value = "", status = "ok"))
    }
    if (kind == "sim") {
        parts <- strsplit(payload, " ", fixed = TRUE)[[1]]
        parameters <- if (length(parts) > 1) json_decode(paste(parts[-1], collapse = " "))
                      else list()
        output <- sim_run(parts[[1]], parameters)
        stopifnot(is.list(output), length(output) > 0)
        return(list(text = paste(sort(names(output)), collapse = ","), value = "", status = "ok"))
    }

    session <- rssn_session()
    on.exit(close(session), add = TRUE)

    if (kind == "apply") {
        values <- lapply(strsplit(row$bindings, ",", fixed = TRUE)[[1]], function(entry) {
            if (!nzchar(trimws(entry))) return(NULL)
            as.numeric(trimws(entry))
        })
        values <- Filter(Negate(is.null), values)
        built <- do.call(op_apply, c(list(payload), values, list(session = session)))
        return(list(text = term_text(built), value = "", status = "ok"))
    }
    if (kind == "callback2") {
        # Two callbacks in one call: the objective and its gradient `num_bfgs`
        # declares, from the same start point every runner uses.
        parts <- strsplit(payload, " ", fixed = TRUE)[[1]]
        arguments <- if (length(parts) > 1) {
            lapply(strsplit(paste(parts[-1], collapse = " "), ",", fixed = TRUE)[[1]], literal)
        } else {
            list()
        }
        objective <- function(point) sum((point - 1)^2)
        gradient <- function(point) 2 * (point - 1)
        answer <- do.call(call_value,
                          c(list(parts[[1]]), list(c(3, -2)), arguments,
                            list(objective, gradient, session = session)))
        return(list(text = sprintf("%.10g", value_float(answer[[1]][[1]])), value = "", status = "ok"))
    }
    if (kind == "callback") {
        parts <- strsplit(payload, " ", fixed = TRUE)[[1]]
        arguments <- if (length(parts) > 1) {
            lapply(strsplit(paste(parts[-1], collapse = " "), ",", fixed = TRUE)[[1]], literal)
        } else {
            list()
        }
        answer <- do.call(call_value,
                          c(list(parts[[1]]), arguments,
                            list(function(x) x * x - 2, session = session)))
        root <- value_float(answer[[1]])
        return(list(text = sprintf("%.10g", root), value = "", status = "ok"))
    }
    answer <- if (kind == "simplify") simplify(payload, session = session)
              else compute(parse_term(payload, session), config_for(row))
    value <- ""
    if (!is.na(answer$value) && trimws(row$target) == "numeric") {
        value <- sprintf("%.17g", answer$value)
    }
    list(text = answer$text, value = value, status = "ok")
}

cases <- read_cases(cases_path)
failures <- 0L
for (row in cases) {
    wanted <- if (nzchar(trimws(row$expect_status))) trimws(row$expect_status) else "ok"
    observed <- tryCatch(run_case(row),
                         error = function(condition) {
                             list(text = "", value = "", status = status_of(condition))
                         })
    name <- row$name
    if (!identical(observed$status, wanted)) {
        cat(sprintf("FAIL %s: status %s != %s\n", name, observed$status, wanted))
        failures <- failures + 1L
        next
    }
    if (wanted != "ok") {
        cat(sprintf("ok   %s (%s)\n", name, wanted))
        next
    }
    if (!identical(observed$text, row$expect_text)) {
        cat(sprintf("FAIL %s: text '%s' != '%s'\n", name, observed$text, row$expect_text))
        failures <- failures + 1L
        next
    }
    if (nzchar(row$expect_value) && nzchar(observed$value)) {
        want <- as.numeric(row$expect_value)
        got <- as.numeric(observed$value)
        if (!is.na(want) && !is.na(got) && abs(want - got) > 1e-9 * max(1, abs(want))) {
            cat(sprintf("FAIL %s: value %s != %s\n", name, observed$value, row$expect_value))
            failures <- failures + 1L
            next
        }
    }
    cat(sprintf("ok   %s\n", name))
}
cat(sprintf("%d/%d cases pass\n", length(cases) - failures, length(cases)))
if (failures > 0L) quit(status = 1L)
