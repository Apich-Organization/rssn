# A tour of the numeric side of rssn, from R: targets, tensors and dimensions.
#
#     RSSN_LIB=/path/to/rssn.dll Rscript inst/examples/numeric.R

library(rssn)

section <- function(title) cat("\n==", title, "==\n")

session <- rssn_session()
x <- sym("x", session)

section("a numeric target")
cfg <- rssn_config() |> cfg_numeric(1e-12) |> cfg_bind("x", 0.5)
answer <- compute(rssn_diff(sin(x), x), cfg)
cat(sprintf("d/dx sin(x) at 0.5 = %.15g (error %.3g)\n", answer$value, answer$error))

section("evaluation without any rule")
cat("sin(0) + 1 =", eval_term(parse_term("sin(x) + 1", session), x = 0), "\n")

section("linear algebra and tensors")
matrix_expression <- parse_term(
    "matmul(list(list(1,2),list(3,4)), list(list(0,1),list(1,0)))", session)
entries <- tensor(compute(matrix_expression))
print(entries)

for (source in c("det(list(list(2,0),list(0,3)))",
                 "trace(list(list(1,2),list(3,4)))",
                 "norm(list(3,4))",
                 "inverse(list(list(1,2),list(3,4)))")) {
    cat(sprintf("%-34s -> %s\n", source, compute(parse_term(source, session))$text))
}

section("a term's value as an array, through the value ABI")
value <- value_new(compute(matrix_expression)$term, session)
cat("kind:", value_kind(value), " extents:", paste(dim(value_tensor(value)), collapse = "x"), "\n")

section("dimensions")
dimension_config <- rssn_config() |>
    cfg_dimension("x", "m") |>
    cfg_dimension("t", "s")
cat("d/dx x^2   ->", dimension(rssn_diff(x^2, x), dimension_config), "\n")
cat("2 + 2      ->", dimension(parse_term("2 + 2", session)), "\n")

strict <- rssn_config() |>
    cfg_dimension("x", "m") |>
    cfg_dimension("y", "s") |>
    cfg_dimension_mode("strict")
condition <- tryCatch({
    compute(parse_term("x + y", session), strict)
    NULL
}, rssn_error = function(condition) condition)
cat("x + y (strict) ->", condition$status_name, ":", condition$message, "\n")

section("units and constants")
cat("ft ->", unit_scale(unit_get("ft")), "m\n")
cat("c  ->", constant_value(constant_get("c")), "m/s\n")
print(head(units_table(), 3))

close(session)
