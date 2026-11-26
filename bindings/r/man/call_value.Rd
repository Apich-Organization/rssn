\name{call_value}
\alias{call_value}
\alias{value_new}
\alias{value_of}
\alias{value_kind}
\alias{value_float}
\alias{value_text}
\alias{value_json}
\alias{value_term}
\alias{value_tensor}
\alias{tensor}
\alias{value_len}
\alias{value_count}
\alias{value_cell}
\alias{VALUE_KINDS}
\alias{close.rssn_value}
\alias{print.rssn_value}
\alias{format.rssn_value}
\alias{as.character.rssn_value}
\alias{as.numeric.rssn_value}
\alias{length.rssn_value}
\alias{[[.rssn_value}
\title{Dynamic values: the cell ABI}
\usage{
call_value(op, ..., session = NULL)

value_new(x, session = NULL)

value_of(x, session = NULL)

value_kind(value)

value_float(value)

value_text(value)

value_json(value)

value_term(value)

value_tensor(value)

tensor(x)

value_len(value)

value_cell(value)

value_count(session = NULL)

VALUE_KINDS

\method{close}{rssn_value}(con, ...)

\method{as.numeric}{rssn_value}(x, ...)

\method{as.character}{rssn_value}(x, ...)

\method{length}{rssn_value}(x)

\method{[[}{rssn_value}(x, index)
}
\arguments{
  \item{op}{The operator name.}
  \item{...}{The arguments: values, terms, numbers, strings, logicals, lists,
    arrays, and functions.}
  \item{session}{The session to call it in; taken from the first value or term
    argument when omitted.}
  \item{x}{The object to store as a cell, or a \code{rssn_value}.}
  \item{value}{A \code{rssn_value}.}
  \item{con}{A \code{rssn_value}.}
  \item{index}{The position of an element, counting from one.}
  \item{term}{A \code{rssn_term} or an \code{rssn_answer}, for \code{tensor()}.}
}
\value{
  \code{call_value()} and \code{value_new()} return an \code{rssn_value};
  the readers return a number, a string, a list, a term, an array or an
  integer.
}
\description{
  The term API carries only values that are nodes of the expression graph.  An
  operator that returns a native object -- an AD tape, a compiled backend, a
  solver handle -- returns a \emph{value} instead, and an operator that takes
  one can only be called with a value.  \code{call_value()} is that path, and
  it is also the one that carries a host callback.

  A value is an environment with a finalizer holding a cell id: the engine can
  reuse a freed id, so an object that remembers whether it still owns its cell
  is the only safe representation.  \code{close()} releases it early; otherwise
  the garbage collector does, and the session releases whatever is left.
}
\details{
  A function argument becomes a host callback, registered with the signature the
  operator declares for that position:

  \preformatted{
  call_value("num_brent", 1, 2, 1e-12, 64, function(x) x * x - 2)
  }

  R has no integer literal -- \code{1} is a double and \code{1L} is an integer
  -- and the value ABI is strict about the difference.  The first failure that
  names an argument and the kind it wanted flips that one argument and retries,
  so \code{64} reaches an \code{Int} parameter as an integer without the caller
  having to write \code{64L}.

  \code{VALUE_KINDS} lists the kind names \code{value_kind()} reports.
}
\examples{
s <- rssn_session()

backend <- call_value("num_backend_current", session = s)
value_kind(backend)                                        # "opaque"
value_text(call_value("num_backend_name", backend, session = s))

tape <- call_value("num_ad_tape_new", 16, session = s)
x <- call_value("num_ad_tape_var", tape, 1.5, session = s)
z <- call_value("num_ad_var_mul", call_value("num_ad_var_sin", x, session = s), x,
                session = s)
value_float(call_value("num_ad_var_value", z, session = s))  # sin(1.5) * 1.5

answer <- call_value("num_brent", 1, 2, 1e-12, 64, function(t) t * t - 2,
                     session = s)
value_float(answer[[1]])                                   # 1.4142135623730951

value_kind(value_new(c(1, 2.5, 3), s))                     # "tensor"

tensor(compute(parse_term("matmul(list(list(1,2),list(3,4)), list(list(1,0),list(0,1)))")))
}
\seealso{\code{\link{rssn_callback}}, \code{\link{term}}}
