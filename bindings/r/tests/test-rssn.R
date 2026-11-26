# Tests for the R binding.
#
#     R CMD INSTALL bindings/r
#     RSSN_LIB=/path/to/rssn.dll Rscript bindings/r/tests/test-rssn.R
#
# They need a built rssn shared library (see `bindings/README.md`). The
# cross-language half of the suite is `bindings/tests/cases.tsv`, which
# `bindings/tests/run_conformance.R` checks; this file covers what the table
# cannot see: the R coercions, the `Ops` and `Math` groups, the JSON codec, the
# callback trampoline and the value objects.

library(rssn)

failures <- 0L
checks <- 0L

expect <- function(condition, label) {
    checks <<- checks + 1L
    if (!isTRUE(condition)) {
        failures <<- failures + 1L
        cat("FAIL ", label, "\n", sep = "")
    }
}

expect_equal <- function(actual, expected, label) {
    checks <<- checks + 1L
    if (!identical(actual, expected)) {
        failures <<- failures + 1L
        cat("FAIL ", label, ": ", paste(format(actual), collapse = " "), " != ",
            paste(format(expected), collapse = " "), "\n", sep = "")
    }
}

expect_error <- function(expression, status, label) {
    checks <<- checks + 1L
    condition <- tryCatch({
        expression
        NULL
    }, rssn_error = function(condition) condition)
    if (is.null(condition)) {
        failures <<- failures + 1L
        cat("FAIL ", label, ": no error\n", sep = "")
        return(invisible(NULL))
    }
    if (!is.null(status) && !identical(condition$status_name, status)) {
        failures <<- failures + 1L
        cat("FAIL ", label, ": status ", condition$status_name, " != ", status, "\n", sep = "")
    }
    invisible(NULL)
}

# -- sessions and terms -----------------------------------------------------

session <- rssn_session()
x <- sym("x", session)

expect_equal(term_text(x), "x", "a symbol renders as its name")
expect_equal(term_text(x * x), "x*x", "a product is a request")
expect_equal(compute(x * x)$text, "x^2", "compute reduces the request")
expect_equal(term_text(sym("y", session) + 1), "y + 1", "addition")
expect_equal(term_text(2 * x), "2*x", "scalar on the left")
expect_equal(term_text(x - 1), "x - 1", "subtraction")
expect_equal(term_text(1 - x), "1 - x", "subtraction on the left")
expect_equal(term_text(-x), "-x", "negation")
expect_equal(term_text(x / 2), "x/2", "quotient")
expect_equal(term_text(2 / x), "2/x", "quotient on the left")
expect_equal(term_text(x^2), "x^2", "power")
expect_equal(term_text(2^x), "2^x", "power on the left")
expect_equal(compute((x^2 - 1) / (x - 1))$text, "x + 1", "cancellation")

# The `Math` group generics reach the engine's own operators.
expect_equal(term_text(sin(x)), "sin(x)", "sin")
expect_equal(compute(rssn_diff(sin(x), x))$text, "cos(x)", "diff")
expect_equal(term_text(sqrt(x)), "sqrt(x)", "sqrt")
expect_equal(term_text(exp(x)), "exp(x)", "exp")
expect_equal(term_text(abs(x)), "abs(x)", "abs")
expect_equal(term_text(floor(x)), "floor(x)", "floor")
expect_equal(term_text(ceiling(x)), "ceil(x)", "ceiling maps to ceil")

# Coercion of numbers, rationals and text.
expect_equal(term_text(term(2, session)), "2", "an integer term")
expect_equal(term_text(term(0.5, session)), "0.5", "a float")
expect_equal(term_text(rational(1, 3, session)), "1/3", "an exact fraction")
expect_equal(term_text(term("1/4", session)), "1/4", "text is parsed")
expect_equal(term_text(rational(2, 3, session)), "2/3", "an exact rational")
expect_equal(term_value(term(7, session)), 7, "the value of a literal")
expect(is.na(term_value(x)), "a symbol has no literal value")

# Renderings.
answer <- compute(rssn_diff(sin(x) * cos(x), x))
expect_equal(answer$text, "cos(2*x)", "the derivative")
expect(grepl("cos", term_latex(answer$term)), "LaTeX")
expect(nzchar(term_typst(answer$term)), "Typst")
expect(nzchar(term_pretty(answer$term)), "the Unicode drawing")
expect(isTRUE(answer$reduced), "the answer is reduced")

# A tensor comes back in reading order.
matrix <- compute(parse_term("matmul(list(list(1,2),list(3,4)), list(list(0,1),list(1,0)))",
                             session))$term
entries <- value_tensor(value_new(matrix, session))
expect_equal(dim(entries), c(2L, 2L), "the tensor's extents")
expect_equal(entries[1, 2], 1, "row-major reading order (1,2)")
expect_equal(entries[2, 1], 4, "row-major reading order (2,1)")

# Batch computation, and the failure policy.
answers <- compute_all(list("diff(x^3, x)", "integral(x^2, x)", "gcd(4, 6)"))
expect_equal(length(answers), 3L, "three answers")
expect(all(vapply(answers, is_answer, logical(1))), "all answers")
numeric_config <- rssn_config() |> cfg_numeric(1e-12)
collected <- compute_all(list("1 + 1", "diff(x^2, x)"), numeric_config, errors = "collect")
expect(is_answer(collected[[1L]]), "the first request succeeds")
expect(inherits(collected[[2L]], "rssn_error"), "the second is collected")
expect_error(compute_all(list("diff(x^2, x)"), numeric_config), "not_numeric", "raise by default")

# The records and the wire format.
record <- compute_json("diff(sin(x), x)", session = session)
expect_equal(record$term, "cos(x)", "the compute record")
expect(inherits(compute_bincode("diff(sin(x), x)", session = session), "raw"), "bincode bytes")
expect(length(compute_all_bincode(list("1 + 1", "gcd(2, 4)"))) > 4, "the batch record")

expect_equal(eval_term(parse_term("sin(x)", session), x = 0), 0, "eval without rules")

# Dimensions.
dimension_config <- rssn_config() |> cfg_dimension("x", "m")
expect_equal(dimension(rssn_diff(x^2, x), dimension_config), "Length", "the inferred dimension")
expect_equal(dimension(parse_term("2 + 2", session)), "1", "a dimensionless answer")
strict <- rssn_config() |> cfg_dimension("x", "m") |> cfg_dimension("y", "s") |>
    cfg_dimension_mode("strict")
expect_error(compute(parse_term("x + y", session), strict), "dimension", "a strict conflict")

# -- errors -----------------------------------------------------------------

expect_error(parse_term("sin(", session), "parse", "a parse failure")
expect_error(op_apply("definitely_not_an_operator", 1, session = session), "unknown_operator",
             "an unknown operator")
expect_error(compute(parse_term("diff(x^2, x)", session), numeric_config), "not_numeric",
             "a symbolic answer under a numeric target")
expect_error(rational(1, 0, session), NULL, "a zero denominator")

# -- values -----------------------------------------------------------------

backend <- call_value("num_backend_current", session = session)
expect_equal(value_kind(backend), "opaque", "the backend is opaque")
expect_equal(value_text(call_value("num_backend_name", backend, session = session)),
             "interpreter", "the backend's name")
expect_error(value_term(backend), NULL, "an opaque handle has no term form")

tape <- call_value("num_ad_tape_new", 16, session = session)
variable <- call_value("num_ad_tape_var", tape, 1.5, session = session)
product <- call_value("num_ad_var_mul",
                      call_value("num_ad_var_sin", variable, session = session),
                      variable, session = session)
expect(abs(value_float(call_value("num_ad_var_value", product, session = session)) -
           sin(1.5) * 1.5) < 1e-12, "the AD tape's value")
gradient <- call_value("num_ad_var_backward", product, session = session)
expect(abs(value_float(call_value("num_ad_var_grad_in", variable, gradient,
                                  session = session)) - (cos(1.5) * 1.5 + sin(1.5))) < 1e-12,
       "the AD tape's gradient")
expect(is.list(value_json(backend)), "a value's JSON form")

# A list argument keeps its element kinds.
list_value <- value_new(list(1L, 2.5, 3L), session)
expect_equal(value_kind(list_value), "tensor", "a list becomes a tensor")
expect_equal(value_len(list_value), 3L, "three elements")
expect_equal(value_float(list_value[[2L]]), 2.5, "the second element")

# Cells are released with the session, and `close()` releases one early.
scratch <- rssn_session()
expect_equal(value_count(scratch), 0L, "no cells yet")
value <- call_value("num_backend_current", session = scratch)
expect_equal(value_count(scratch), 1L, "one cell")
close(value)
expect_equal(value_count(scratch), 0L, "released")
close(scratch)

# -- host callbacks ---------------------------------------------------------

brent <- call_value("num_brent", 1, 2, 1e-12, 64, function(t) t * t - 2, session = session)
expect(abs(value_float(brent[[1]]) - sqrt(2)) < 1e-10, "the root of x^2 - 2")

explicit <- rssn_callback(function(t) t * t - 2, session, signature = "(float) -> float")
first <- call_value("num_brent", 1, 2, 1e-12, 64, explicit, session = session)
expect(abs(value_float(first[[1]]) - sqrt(2)) < 1e-10, "an explicit callback")
close(explicit)
expect_equal(callback_count(session), 0L, "the callback is retired")

signatures <- callback_signatures("num_brent", 5, session = session)
expect_equal(length(signatures), 1L, "one declared signature")
expect_equal(signatures[[1L]]$params[[1L]], "float", "the declared parameter")
expect_equal(length(callback_signatures("num_brent", 1, session = session)), 0L,
             "a non-function position declares nothing")

expect_equal(parse_signature("(float, float) -> float"),
             '{"params":["float","float"],"ret":"float"}', "the compact signature grammar")
expect_equal(parse_signature("() -> bool"), '{"params":[],"ret":"bool"}', "a nullary signature")

# An error inside the callback aborts the operator with its message.
refusal <- tryCatch({
    call_value("num_brent", 1, 2, 1e-12, 64, function(t) stop("no thanks"), session = session)
    NULL
}, rssn_error = function(condition) condition)
expect(inherits(refusal, "rssn_error"), "a failing callback raises")
expect(grepl("no thanks", refusal$message), "the R message crosses")
expect_equal(callback_count(session), 0L, "the implicit callback is retired too")

# Two callbacks in one call: an objective and its gradient.
point <- call_value("num_bfgs", c(3, -2), 1e-10, 100, 1000,
                    function(p) sum((p - 1)^2),
                    function(p) 2 * (p - 1),
                    session = session)
expect(abs(value_float(point[[1]][[1]]) - 1) < 1e-6, "the minimiser's first coordinate")
expect(abs(value_float(point[[1]][[2]]) - 1) < 1e-6, "the minimiser's second coordinate")

# -- metadata ---------------------------------------------------------------

expect(length(op_names()) == 5411L, "every operator is registered")
expect_equal(op_help("diff")$rule_set, "calculus", "the operator reference")
expect(nrow(op_table(pattern = "brent")) >= 3L, "the search")
expect_equal(OPS[["diff"]], c(2L, 2L), "the arity table")
expect_equal(OPS[["add"]], c(0L, NA_integer_), "a variadic operator")
expect(nzchar(manifest_hash()), "the manifest hash")
expect(nrow(manifest()$ops) == 5411L, "the manifest")
expect(nzchar(manifest()$rssn_version), "the manifest version")
expect(length(manifest()$rule_sets) > 30L, "the rule sets")

# -- units, constants, scenarios --------------------------------------------

expect_equal(length(UNITS), 82L, "the unit table")
expect_equal(unit_get("m")$dimension, "Length", "a unit's dimension")
expect_equal(unit_get("metre")$name, "m", "an alias")
expect(abs(unit_scale(unit_get("ft")) - 0.3048) < 1e-12, "the foot's exact scale")
expect_equal(length(CONSTANTS), 30L, "the constant table")
expect_equal(constant_get("c")$symbol, "c", "a constant by symbol")
expect(abs(constant_value(constant_get("c")) - 299792458) < 1e-6, "the speed of light")
expect("ising" %in% sim_names(), "the scenario catalog")
expect_equal(sim_defaults("poisson_multigrid")$dim, 2L, "a scenario's defaults")
result <- sim_run("ising", width = 4, height = 4, temperature = 2, mc_steps = 5, seed = 7)
expect(all(c("spins", "magnetization") %in% names(result)), "the Ising result")

# -- the JSON codec ---------------------------------------------------------

expect_equal(json_decode("{}"), structure(list(), names = character()), "an empty object")
expect_equal(json_decode("[]"), list(), "an empty array")
expect_equal(json_decode("1"), 1L, "an integer")
expect_equal(json_decode("1.5"), 1.5, "a double")
expect_equal(json_decode("true"), TRUE, "a boolean")
expect(is.null(json_decode("null")), "null")
expect_equal(json_decode("\"a\\nb\""), "a\nb", "escapes")
expect_equal(json_decode("{\"a\": [1, 2.5, \"x\"]}")$a[[3L]], "x", "nesting")
expect_equal(json_encode(list(a = c(1, 2))), '{"a":[1,2]}', "encoding")
expect(is.null(json_decode("{\"a\": null}")$a), "a null field survives")
expect_equal(json_decode(json_encode(list(a = 1, b = "x", c = TRUE))),
             list(a = 1L, b = "x", c = TRUE), "a round trip")

# -- the default session ----------------------------------------------------

expect(identical(default_session(), default_session()), "the default session is stable")
with_session(session, expect(identical(default_session(), session), "with_session installs it"))
expect(!identical(default_session(), session), "and restores the previous one")

cat(sprintf("%d/%d checks pass\n", checks - failures, checks))
if (failures > 0L) quit(status = 1L)
