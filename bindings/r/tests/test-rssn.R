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

# The installed bridge has its own namespace and cannot serve as the kernel.
bridge <- normalizePath(getLoadedDLLs()[["rssn"]][["path"]], winslash = "/")
installed <- system.file("libs", package = "rssn")
candidates <- rssn_library_candidates()
expect(!installed %in% candidates, "automatic discovery inspects installed libraries by file")
expect(!bridge %in% normalizePath(candidates[file.exists(candidates)], winslash = "/"),
       "automatic discovery skips the R bridge")

# Check a nested source working directory without relying on Cargo artifacts.
check_bazel_candidates <- function() {
    previous <- getwd()
    root <- tempfile("rssn-checkout-")
    child <- file.path(root, "examples")
    dir.create(child, recursive = TRUE)
    file.create(file.path(root, "rssn.h"), file.path(root, "MODULE.bazel"))
    on.exit({
        setwd(previous)
        unlink(root, recursive = TRUE)
    })
    setwd(child)
    candidates <- rssn_library_candidates()
    targets <- c("ffi_distribution", "ffi", "ffi_distribution_portable", "ffi_dev")
    bazel <- match(file.path("..", "bazel-bin", targets, "ffi"), candidates)
    cargo <- match(file.path(".", "target", "debug"), candidates)
    expect(!anyNA(bazel) && !is.na(cargo) &&
           all(diff(bazel) > 0) && max(bazel) < cargo,
           "automatic discovery prefers Bazel library directories in order")
}
check_bazel_candidates()

# A directory resolves to the same library. A subsequent failed rebind leaves
# the function table and any existing session usable.
bound_library <- rssn_library()
if (file.exists(bound_library)) {
    expect_equal(normalizePath(rssn_load(dirname(bound_library))), normalizePath(bound_library),
                 "load from a directory")
}

# -- sessions and terms -----------------------------------------------------

session <- rssn_session()
x <- sym("x", session)

expect_equal(term_text(x), "x", "a symbol renders as its name")
bad_library <- tempfile(fileext = .Platform$dynlib.ext)
writeLines("invalid rssn library", bad_library)
failed_load <- tryCatch(rssn_load(bad_library), error = identity)
unlink(bad_library)
expect(inherits(failed_load, "error"), "a bad library reports a loading error")
bridge_path <- getLoadedDLLs()[["rssn"]][["path"]]
incompatible <- tryCatch(rssn_load(bridge_path), error = identity)
expect(inherits(incompatible, "error") && grepl("no symbol", conditionMessage(incompatible)),
       "a shared library without RSSN's ABI is rejected")
expect_equal(term_text(x), "x", "a failed rebind leaves the session usable")
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
expect(is.na(answer$value), "a symbolic answer has no numeric value")
expect_equal(compute("1 + 1")$value, 2, "a numeric answer has a value")

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
quantity_requests <- list("quantity(3, m) + quantity(4, s)", "quantity(2, m)")
metadata <- function(ans) list(ans$dimension, ans$diagnostic_codes)
individual <- lapply(quantity_requests, compute)
expect("mismatch" %in% individual[[1L]]$diagnostic_codes, "advisory mismatch is reported")
expect(identical(lapply(compute_all(quantity_requests), metadata), lapply(individual, metadata)),
       "batch metadata matches each individual request")
numeric_config <- rssn_config() |> cfg_numeric(1e-12)
collected <- compute_all(list("1 + 1", "diff(x^2, x)"), numeric_config, errors = "collect")
expect(is_answer(collected[[1L]]), "the first request succeeds")
expect(inherits(collected[[2L]], "rssn_error"), "the second is collected")
expect_error(compute_all(list("diff(x^2, x)"), numeric_config), "not_numeric", "raise by default")
ignored <- compute_all(list("1 + 1", "diff(x^2, x)", "3 + 4"), numeric_config, errors = "ignore")
expect_equal(length(ignored), 3L, "ignore keeps every batch position")
expect(is.null(ignored[[2L]]), "ignore leaves NULL at the failed position")
expect_equal(ignored[[3L]]$value, 7, "a later success keeps its position")
mixed <- compute_all(list("1 + 1", x))
expect_equal(mixed[[2L]]$text, "x", "mixed input uses the supplied term's session")
other <- rssn_session()
foreign <- sym("foreign", other)
for (batch in list(compute_all, compute_all_bincode)) {
    condition <- tryCatch(batch(list(x, foreign)), error = identity)
    expect(inherits(condition, "error") && grepl("another session", conditionMessage(condition)),
           "a batch rejects terms from two sessions")
}
expect(is.raw(compute_all_bincode(list(), session = session)), "empty batch record")

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
             rssn_current_backend(session), "the backend's name")
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

moments <- call_value("num_timeseries_rolling_moments", c(1, 2, 3), 2L, 1L, session = session)
expect_equal(value_tensor(moments[["variance"]]), array(c(0.5, 0.5)), "a record field by name")
expect_equal(value_tensor(moments[[2L]]), array(c(0.5, 0.5)), "a record field by position")
missing_field <- tryCatch(moments[["var"]], error = identity)
expect(inherits(missing_field, "error") &&
       startsWith(conditionMessage(missing_field), "no record field named "),
       "record field names must match exactly")

for (original in list(array(9, dim = c(1L, 1L, 1L)), array(numeric(), dim = c(0L, 2L, 3L)),
                      array(as.double(seq_len(24L)), dim = c(2L, 3L, 4L)))) {
    stored <- value_new(original, session)
    expect_equal(value_tensor(stored), original, "array round trips with its extents")
    close(stored)
}

digits <- "-1267650600228229401496703205383"
exact <- value_integer_decimal(digits, session)
expect_equal(value_kind(exact), "int", "decimal integer kind")
expect_equal(value_text(exact), digits, "decimal integer precision")
close(exact)

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

center <- 0.3 + 0.4i
sampled <- complex()
derivative <- rssn_num_complex_derivative_scaled(center, 1.0, 1.0, 64L, function(z) {
    sampled <<- c(sampled, z)
    z * z
}, session = session)
parts <- value_json(derivative)$complex$float
expect(Mod(complex(real = parts$re, imaginary = parts$im) - 2 * center) < 1e-8,
       "a callback receives complex arguments across the ABI")
expect(length(sampled) > 0L && all(Mod(sampled - center) < 1e-4),
       "the callback receives the operator's complex samples")
close(derivative)

# Two callbacks in one call: an objective and its gradient.
point <- call_value("num_bfgs", c(3, -2), 1e-10, 100, 1000,
                    function(p) sum((p - 1)^2),
                    function(p) 2 * (p - 1),
                    session = session)
expect(abs(value_float(point[[1]][[1]]) - 1) < 1e-6, "the minimiser's first coordinate")
expect(abs(value_float(point[[1]][[2]]) - 1) < 1e-6, "the minimiser's second coordinate")

# -- metadata ---------------------------------------------------------------

reference <- Sys.getenv("RSSN_CONFORMANCE_OP_NAMES")
if (nzchar(reference)) {
    live_names <- readLines(reference, encoding = "UTF-8")
    expect_equal(sort(op_names()), sort(live_names), "every live operator has an R binding")
}
expect_equal(op_help("diff")$rule_set, "calculus", "the operator reference")
expect(nrow(op_table(pattern = "brent")) >= 3L, "the search")
expect_equal(OPS[["diff"]], c(2L, 2L), "the arity table")
expect_equal(OPS[["add"]], c(0L, NA_integer_), "a variadic operator")
expect(nzchar(manifest_hash()), "the manifest hash")
expect_equal(sort(manifest()$ops$name), sort(op_names()), "the packaged manifest names")
expect(nzchar(manifest()$rssn_version), "the manifest version")
expect(length(manifest()$rule_sets) > 30L, "the rule sets")

(function() {
    s <- rssn_session()
    on.exit(close(s))
    retained <- value_new(0, s)
    independent <- value_new(retained)
    close(independent)
    expect_equal(value_float(retained), 0, "value_new keeps its source alive")
    callback <- rssn_callback(function(x) retained, s, signature = "(float) -> float")
    answer <- rssn_num_brent(1, 2, 1e-12, 64L, callback)
    expect(identical(answer$session, s), "callback's session inferred")
    close(answer)
    close(callback)
    before <- value_count(s)
    record <- value_new(list(a = list(retained, 2L), b = TRUE), s)
    close(record)
    expect_equal(value_count(s), before, "temporary record and list members released")
    expect(inherits(try(value_new(list(1, NA_character_), s), silent = TRUE), "try-error"),
           "invalid nested member is rejected")
    expect_equal(value_count(s), before, "partial construction releases earlier members")
    expect_equal(value_float(retained), 0, "callback return keeps its source alive")
    close(retained)
})()

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

# -- checked dynamic facade --------------------------------------------------
request <- rssn_num_special_ext_lerch(-1234.0, 2.0, 1.0, session = session)
expect(inherits(request, "rssn_term"), "the dynamic facade still builds a term by default")
lerch <- rssn_num_special_ext_lerch(-1234.0, 2.0, 1.0, session = session, checked = TRUE)
polylog <- rssn_num_special_ext_polylog(2.0, -1234.0, session = session, checked = TRUE)
expect(abs(value_float(lerch) - 0.021861589304464496) < 1e-15, "a checked real Lerch call")
expect(abs(value_float(polylog) / -1234.0 - value_float(lerch)) < 1e-15, "checked polylog agrees")
point <- rssn_num_uhp_to_disk(1i, session = session, checked = TRUE)
expect_equal(value_json(point)$complex$float, list(re = 0.0, im = 0.0), "a checked complex call")
recovered <- rssn_num_disk_to_uhp(point, checked = TRUE)
expect_equal(value_json(recovered)$complex$float, list(re = 0.0, im = 1.0), "infer the session from a value")
expect_error(rssn_num_disk_to_uhp(1 + 0i, session = session, checked = TRUE),
             "invalid_argument", "checked facade preserves domain errors")

# -- backend selection and NumPy arrays -------------------------------------
(function() {
    s <- rssn_session()
    other <- rssn_session()
    on.exit({
        rssn_set_default_backend("auto", s)
        close(other)
        close(s)
    })
    before <- value_text(call_value("num_native_usage", session = s)[["objects"]])
    expect(nzchar(rssn_current_backend(s)), "query the active backend")
    expect_equal(value_text(call_value("num_native_usage", session = s)[["objects"]]),
                 before, "the effective query leaves no native objects")
    expect("interpreter" %in% rssn_available_backends(s), "the portable backend is available")
    handle <- rssn_open_backend("interpreter", s)
    expect_equal(value_text(call_value("num_backend_name", handle)),
                 "interpreter", "open an explicit interpreter backend")
    rssn_release_backend(handle)
    expect_error(call_value("num_backend_name", handle), NULL, "released backend is invalid")
    close(handle)

    expect_equal(rssn_set_default_backend("interpreter", s), "interpreter", "select an explicit default")
    expect_equal(rssn_current_backend(other), "interpreter", "sessions share the process default")
    expect_error(rssn_set_default_backend("not-an-rssn-backend", other), NULL, "unknown backend fails")
    expect_equal(rssn_current_backend(other), "interpreter", "failed selection leaves the default intact")

    if (isTRUE(rssn_features()$npy)) {
        samples <- matrix(c(-0.0, 1.0, 2.0, 3.0, Inf, NaN), nrow = 2L, byrow = TRUE)
        payload <- rssn_encode_npy(samples, session = s)
        expect_equal(as.integer(payload[seq_len(6L)]), c(147L, 78L, 85L, 77L, 80L, 89L), "NumPy magic")
        restored <- rssn_decode_npy(payload, s)
        array <- value_tensor(restored)
        expect_equal(array, samples, "NumPy shape, row order and IEEE samples")
        expect(1 / array[1L, 1L] == -Inf, "the NumPy codec preserves negative zero")
        close(restored)

        scalar <- rssn_decode_npy(rssn_encode_npy(-0.0, integer(), s), s)
        expect(is.null(dim(value_tensor(scalar))) && 1 / value_tensor(scalar) == -Inf,
               "the NumPy codec supports a scalar")
        close(scalar)

        empty <- rssn_decode_npy(rssn_encode_npy(numeric(), c(2L, 0L, 3L), s), s)
        expect_equal(dim(value_tensor(empty)), c(2L, 0L, 3L), "the NumPy codec preserves empty extents")
        close(empty)
        expect_error(rssn_decode_npy(as.raw(c(147, 78)), s), NULL, "truncated NumPy bytes fail")
    }
})()

# -- the default session ----------------------------------------------------

expect(identical(default_session(), default_session()), "the default session is stable")
with_session(session, expect(identical(default_session(), session), "with_session installs it"))
expect(!identical(default_session(), session), "and restores the previous one")

cat(sprintf("%d/%d checks pass\n", checks - failures, checks))
if (failures > 0L) quit(status = 1L)
