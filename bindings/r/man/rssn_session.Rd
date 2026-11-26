\name{rssn_session}
\alias{rssn_session}
\alias{close.rssn_session}
\alias{rssn_close}
\alias{default_session}
\alias{with_session}
\alias{session_ptr}
\alias{print.rssn_session}
\alias{format.rssn_session}
\title{A rssn session}
\usage{
rssn_session()

session_ptr(session)

rssn_close(session)

\method{close}{rssn_session}(con, ...)

default_session()

with_session(session, code)
}
\arguments{
  \item{session}{A \code{rssn_session}.}
  \item{con}{A \code{rssn_session}.}
  \item{code}{An expression, evaluated in the caller's environment.}
  \item{...}{Ignored.}
}
\value{
  \code{rssn_session()} returns an object of class \code{rssn_session}.
  \code{default_session()} returns the calling thread's session, creating it on
  first use.  \code{with_session()} returns the value of \code{code}.
  \code{session_ptr()} returns the external pointer the C bridge holds.
}
\description{
  A session owns every term made through it.  It is single-threaded: its store
  sits behind a \code{RefCell} in Rust, so two threads using one session at once
  is a borrow conflict.  Use one session per thread, or \code{\link{compute_all}()}
  for parallel requests inside one session.
}
\details{
  The operator functions take the session from their first term argument, and
  fall back to \code{default_session()}.  \code{with_session()} scopes that
  fallback, which is what makes the unqualified form read naturally:

  \preformatted{
  s <- rssn_session()
  x <- sym("x", s)
  with_session(s, compute(rssn_diff(x^2, x))$text)   # "2*x"
  }

  A session is released by \code{close()}, which extends the \code{close()}
  generic, or by the garbage collector.
}
\examples{
s <- rssn_session()
x <- sym("x", s)
compute(rssn_diff(x^3, x))$text   # "3*x^2"
close(s)

x <- with_session(rssn_session(), sym("x"))
}
\seealso{\code{\link{term}}, \code{\link{compute}}}
