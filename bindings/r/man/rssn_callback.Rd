\name{rssn_callback}
\alias{rssn_callback}
\alias{callback_cell}
\alias{callback_signatures}
\alias{callback_count}
\alias{parse_signature}
\alias{close.rssn_callback}
\alias{print.rssn_callback}
\title{Host callbacks: an R function the engine calls}
\usage{
rssn_callback(function_, session = NULL, signature = NULL, op = NULL,
              index = 1L)

callback_cell(callback)

callback_signatures(op, index, session = NULL)

callback_count(session = NULL)

parse_signature(text)

\method{close}{rssn_callback}(con, ...)
}
\arguments{
  \item{function_}{The function.  It receives one R object per argument.}
  \item{session}{The session the callback belongs to.}
  \item{signature}{The function's type: a compact \code{"(float) -> float"}, a
    JSON string, or the list form.}
  \item{op}{The operator the callback is for.}
  \item{index}{The argument position of \code{op} the callback fills, counting
    from one.}
  \item{callback}{A \code{rssn_callback}.}
  \item{con}{A \code{rssn_callback}.}
  \item{...}{Ignored.}
  \item{text}{A compact signature.}
}
\value{
  \code{rssn_callback()} returns an object of class \code{rssn_callback};
  \code{callback_signatures()} a list of signature objects;
  \code{parse_signature()} a JSON string.
}
\description{
  189 operators are higher-order: they take a function and call it once per
  quadrature node, optimisation step or root-finder probe.  Any R function
  works, because the package turns one into the C function pointer the engine
  calls.  The function receives native R objects -- numbers, integers, logicals,
  complex numbers, strings, arrays for tensors and named lists for records -- and
  whatever it returns is converted back into a value.
}
\details{
  \code{\link{call_value}()} registers a function argument automatically, for
  the duration of one call.  An explicit callback is for a function that goes to
  several calls, or one whose type has to be spelled out.

  An error inside the function is caught, reported to the engine as
  \code{InvalidArgument}, and the enclosing call raises an \code{rssn_error}
  condition whose message is the R condition's.

  The engine calls the function on the thread that made the call, so the R API
  is used from the thread it belongs to.
}
\examples{
s <- rssn_session()

f <- rssn_callback(function(x) x * x - 2, s, signature = "(float) -> float")
answer <- call_value("num_brent", 1, 2, 1e-12, 64, f, session = s)
value_float(answer[[1]])          # 1.4142135623730951
close(f)

callback_signatures("num_brent", 5, session = s)
parse_signature("(float, float) -> float")
}
\seealso{\code{\link{call_value}}}
