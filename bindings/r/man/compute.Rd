\name{compute}
\alias{compute}
\alias{simplify}
\alias{compute_all}
\alias{compute_json}
\alias{compute_bincode}
\alias{compute_all_bincode}
\alias{eval_term}
\alias{dimension}
\alias{dimension_report}
\alias{answer}
\alias{is_answer}
\alias{print.rssn_answer}
\alias{format.rssn_answer}
\title{Reduce a request}
\usage{
compute(x, config = NULL, session = NULL)

simplify(source, session = NULL)

compute_all(requests, config = NULL, session = NULL,
            errors = c("raise", "collect", "ignore"))

compute_json(x, config = NULL, session = NULL)

compute_bincode(x, config = NULL, session = NULL)

compute_all_bincode(requests, config = NULL, session = NULL)

eval_term(term, ...)

dimension(x, config = NULL, session = NULL)

dimension_report()

answer(term, text, value = NA_real_, error = NA_real_, reduced = TRUE,
       iterations = 0, dimension = NULL, diagnostic_codes = character())

is_answer(x)
}
\arguments{
  \item{x}{A term, or infix text.}
  \item{source}{Infix text.}
  \item{config}{A \code{rssn_config}; the symbolic default when omitted.}
  \item{session}{The session to parse text in.}
  \item{requests}{A list of terms or infix strings.}
  \item{errors}{What to do with a request that failed: \code{"raise"} raises
    the first failure, \code{"collect"} puts the condition in its place and
    \code{"ignore"} puts \code{NULL} there.}
  \item{term}{A \code{rssn_term}.}
  \item{...}{Named numbers binding the free symbols, such as \code{x = 0.5}.}
  \item{text, value, error, reduced, iterations, dimension, diagnostic_codes}{
    The fields of an answer.}
}
\value{
  \code{compute()}, \code{simplify()} and \code{compute_all()} return an
  \code{rssn_answer} (a list of them, for \code{compute_all()}).  The record
  functions return a list or a raw vector; \code{eval_term()} a number;
  \code{dimension()} a string.
}
\description{
  \code{compute()} reduces a request to the target of a configuration: a closed
  form by default, a number with a numeric target.  An answer carries the
  reduced term, its text, the numeric value and error when the engine knows
  them, whether the request was fully reduced, the iterations it took, the
  inferred dimension and the stable codes of the advisory dimensional findings.
}
\details{
  \code{compute_all()} runs the requests in parallel inside the library, which
  is the supported way to use several cores with one session.

  \code{eval_term()} evaluates a term numerically without running any rule.

  \code{dimension()} infers a dimension from the declarations a configuration
  carries (\code{\link{cfg_dimension}()}); \code{dimension_report()} lists the
  codes of the findings of the last computation on this thread, whose prose is
  in the \code{"diagnostics"} field of \code{compute_json()}.
}
\examples{
compute("diff(sin(x)*cos(x), x)")$text                  # "cos(2*x)"
compute("diff(sin(x)*cos(x), x)",
        rssn_config() |> cfg_numeric(1e-12) |> cfg_bind("x", 0))$value

simplify("(x^2 - 1)/(x - 1)")$text                      # "x + 1"
length(compute_all(list("diff(x^3, x)", "gcd(4, 6)")))

cfg <- rssn_config() |> cfg_dimension("x", "m")
x <- sym("x")
dimension(rssn_diff(x^2, x), cfg)                       # "Length"

eval_term(parse_term("sin(x) + 1"), x = 0)              # 1
}
\seealso{\code{\link{rssn_config}}, \code{\link{term}}}
