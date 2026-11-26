\name{term}
\alias{term}
\alias{sym}
\alias{parse_term}
\alias{int_term}
\alias{real_term}
\alias{rational}
\alias{handle}
\alias{term_text}
\alias{term_latex}
\alias{term_typst}
\alias{term_pretty}
\alias{term_value}
\alias{to_float}
\alias{as.character.rssn_term}
\alias{format.rssn_term}
\alias{print.rssn_term}
\alias{as.numeric.rssn_term}
\alias{Ops.rssn_term}
\alias{Math.rssn_term}
\alias{Summary.rssn_term}
\title{Terms: the nodes of a session's expression graph}
\usage{
sym(symbol, session = NULL)

parse_term(source, session = NULL)

int_term(value, session = NULL)

real_term(value, session = NULL)

rational(numerator, denominator = 1, session = NULL)

term(x, session = NULL)

handle(term)

term_text(term)

term_latex(term)

term_typst(term)

term_pretty(term)

term_value(term)

to_float(term, tolerance = 1e-12)

\method{as.character}{rssn_term}(x, ...)

\method{as.numeric}{rssn_term}(x, ...)

\method{Ops}{rssn_term}(e1, e2)

\method{Math}{rssn_term}(x, ...)

\method{Summary}{rssn_term}(..., na.rm = FALSE)
}
\arguments{
  \item{symbol}{The symbol's name.}
  \item{source}{Infix text such as \code{"diff(sin(x)^2, x)"}.}
  \item{value}{A number.}
  \item{numerator, denominator}{Integers; a zero denominator is an error.}
  \item{x, term}{A \code{rssn_term}.}
  \item{e1, e2}{The operands of a binary operator.}
  \item{session}{The session to build the term in; the default session when
    omitted.}
  \item{tolerance}{The tolerance of the numeric target \code{to_float()} uses.}
  \item{...}{Further arguments: the rest of a \code{Math} call, or the
    arguments of \code{sum()}.}
  \item{na.rm}{Ignored.}
}
\value{
  A \code{rssn_term} (an integer handle with the session attached), or for the
  accessors a string, a number, or \code{NA}.
}
\description{
  A term is a \code{uint32_t} handle into one session's graph, carried as an
  integer vector with the session as an attribute, so a term is copyable and
  cheap.

  Operators are \emph{requests}, not results: \code{rssn_diff(sin(x) * x, x)}
  builds the unevaluated derivative and \code{\link{compute}()} reduces it.
  \code{term_text()} shows the request.
}
\details{
  Arithmetic and the \code{Math} group work on a term directly, through the
  engine's own \code{add}, \code{sub}, \code{mul}, \code{divide}, \code{pow} and
  \code{neg} operators:

  \preformatted{
  x^2 + 2*x + 1     # a request; compute() gives "x^2 + 2*x + 1"
  sin(x)            # Math.rssn_term; the engine's sin
  sqrt(x), exp(x), log(x), abs(x), floor(x), ceiling(x), gamma(x), ...
  }

  The engine has no arithmetic quotient named \code{div} (that is the
  divergence of a vector field), so \code{/} is \code{divide}.

  A comparison (\code{==}, \code{<}, ...) builds a term too: the engine's
  \code{eq}, \code{lt}, \code{gt} and their siblings.

  \code{term()} coerces: a number becomes an exact integer when it has no
  fractional part, a string is parsed as infix text, and a \code{rssn_term}
  belonging to another session is refused.
}
\examples{
x <- sym("x")

term_text(x^2 + 1)                        # "x^2 + 1"
compute(x^2 + 1)$text                     # "x^2 + 1"
term_text(sin(x))                         # "sin(x)"
compute(rssn_diff(sin(x), x))$text        # "cos(x)"
term_value(term(7))                       # 7
to_float(parse_term("2 + 3"))             # 5

term_latex(x^2)                           # "x^{2}"
term_typst(x^2)
term_pretty(parse_term("integral(x^2, x)"))
}
\seealso{\code{\link{compute}}, \code{\link{rssn_session}}, \code{\link{op_apply}}}
