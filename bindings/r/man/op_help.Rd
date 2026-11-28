\name{op_help}
\alias{op_help}
\alias{op_table}
\alias{op_names}
\alias{op_apply}
\alias{OPERATORS}
\alias{OPS}
\alias{print.rssn_operator}
\title{The operator reference}
\usage{
op_help(name)

op_table(pattern = NULL, category = NULL, rule_set = NULL, limit = 200L)

op_names()

op_apply(op, ..., session = NULL)

OPERATORS

OPS
}
\arguments{
  \item{name}{An operator name, as the engine spells it (\code{"diff"}).}
  \item{pattern}{A regular expression matched against the name and the
    documentation; \code{NULL} keeps every operator.}
  \item{category}{One manifest category (\code{"calculus"},
    \code{"numerical_matrix"}, ...); \code{NULL} keeps all of them.}
  \item{rule_set}{One rule set; \code{NULL} keeps all of them.}
  \item{limit}{The maximum number of rows to return.}
  \item{op}{The operator name to apply.}
  \item{...}{The arguments of the operator.}
  \item{session}{The session to build the request in.}
}
\value{
  \code{op_help()} returns a one-row data frame of class \code{rssn_operator};
  \code{op_table()} a data frame of the matching operators; \code{op_names()} a
  sorted character vector; \code{op_apply()} a \code{rssn_term} (or a value, for
  a higher-order operator).
}
\description{
  Every registered operator has a function here, named \code{rssn_<name>()}:
  \code{rssn_diff()}, \code{rssn_integral()}, \code{rssn_factorint()}, and 5408
  more.  The prefix keeps the generated names out of base R's way, so nothing
  masks \code{base::diff}.

  \code{OPERATORS} is the reference the generated functions came from: one row
  per operator with its signature, arity, category, rule set, what it returns
  and its documentation.  It is generated from the engine's kernel manifest, so
  it cannot drift from the engine.  \code{OPS} is the arity table the facade
  uses.
}
\details{
  \code{op_apply()} is the escape hatch: it reaches any operator by name at run
  time, including one whose name is computed.

  A higher-order operator -- one that takes a function -- returns a value rather
  than a term, and is called through the value ABI; \code{OPERATORS$returns} is
  \code{"value"} for exactly those.

  Arithmetic and the \code{Math} group generics cover the mathematical
  operators directly on a term, so \code{sin(x)}, \code{x^2} and \code{sqrt(x)}
  need no generated function.
}
\examples{
op_help("diff")
op_help("num_brent")

head(op_table(pattern = "brent"))
nrow(op_table(category = "calculus"))

length(op_names())                       # count for this binding release
op_apply("diff", sym("x")^3, sym("x"))
}
\seealso{\code{\link{term}}, \code{\link{compute}}}
