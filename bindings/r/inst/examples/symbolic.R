# A tour of the symbolic side of rssn, from R.
#
#     RSSN_LIB=/path/to/rssn.dll Rscript inst/examples/symbolic.R
#
# The model is the engine's: an operator applied to arguments builds an
# *unevaluated request*, and `compute()` reduces it.

library(rssn)

section <- function(title) cat("\n==", title, "==\n")

session <- rssn_session()
x <- sym("x", session)

section("a request is not a result")
t <- rssn_diff(sin(x) * cos(x), x)
cat("request :", term_text(t), "\n")
cat("answer  :", compute(t)$text, "\n")

section("calculus")
for (source in c("diff(x^3, x)",
                 "diff(diff(x^3, x), x)",
                 "integral(x^2, x)",
                 "defint(x^2, x, 0, 1)",
                 "limit(sin(x)/x, x, 0)",
                 "taylor(sin(x), x, 0, 5)",
                 "sum(1/(k*(k+1)), k, 1, n)")) {
    cat(sprintf("%-32s -> %s\n", source, compute(parse_term(source, session))$text))
}

section("algebra")
for (source in c("expand((x+1)^3)", "factor(x^2 - 1)", "solve(x^2 - 4, x)",
                 "(x^2 - 1)/(x - 1)")) {
    cat(sprintf("%-32s -> %s\n", source, compute(parse_term(source, session))$text))
}

section("arithmetic and the Math group")
cat("2*x + 1                ->", term_text(2 * x + 1), "\n")
cat("x/2                    ->", term_text(x / 2), "\n")
cat("sin(x)^2 + cos(x)^2    ->", compute(sin(x)^2 + cos(x)^2)$text, "\n")
cat("sqrt(x)                ->", term_text(sqrt(x)), "\n")

section("renderings")
answer <- compute(rssn_integral(x^2, x))
cat("text   :", term_text(answer$term), "\n")
cat("latex  :", term_latex(answer$term), "\n")
cat("typst  :", term_typst(answer$term), "\n")
cat("pretty :\n")
cat(term_pretty(parse_term("integral(x^2, x)", session)), "\n")

section("a batch, in parallel inside the library")
answers <- compute_all(list("diff(x^3, x)", "integral(x^2, x)", "gcd(4, 6)"))
for (answer in answers) cat(" ", answer$text, "\n")

close(session)
