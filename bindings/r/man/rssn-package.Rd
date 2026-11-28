\name{rssn-package}
\alias{rssn-package}
\alias{rssn}
\docType{package}
\title{Scientific Computing with the rssn Engine}
\description{
  A computer algebra system, a set of numerical kernels and a physics
  simulation engine behind one small API.  This package binds the rssn C ABI:
  a \emph{session} owns an expression graph, \emph{operators} are requests that
  \code{\link{compute}()} reduces, and every one of the engine's registered
  registered operators has a function here.

  The bindings are generated from the engine's own kernel manifest, so the
  operator table, the unit table, the physical constants and the simulation
  catalog cannot drift from the engine.
}
\details{
  The workflow is the engine's own, in R's idiom:

  \enumerate{
    \item A session (\code{\link{rssn_session}()}) owns every term.
    \item Build a request: \code{sym("x")}, arithmetic (\code{x^2 + 1}), the
      \code{Math} group (\code{sin(x)}, \code{sqrt(x)}) or a generated
      operator (\code{rssn_diff(f, x)}, \code{rssn_integral(f, x)}).
    \item Reduce it with \code{\link{compute}()}, to a closed form or, with a
      numeric target, to a number.
  }

  Operators whose arguments have no term form -- AD tapes, compiled backends,
  solver handles -- and the higher-order operators that take a function are
  reached through the value ABI with \code{\link{call_value}()}.  An ordinary R
  function passed there becomes a host callback.

  \code{\link{op_help}()} prints one operator's documentation and
  \code{\link{op_table}()} searches the whole reference.
}
\section{The library}{
  The package does not link against rssn: it loads the shared library at run
  time from \code{RSSN_LIB}, \code{options(rssn.library = )}, the package's own
  \code{libs/} directory, or the build directories of a source checkout.
  \code{\link{rssn_library}()} reports what it bound and
  \code{\link{rssn_load}()} rebinds.
}
\seealso{
  \code{\link{rssn_session}}, \code{\link{term}}, \code{\link{compute}},
  \code{\link{rssn_config}}, \code{\link{call_value}}, \code{\link{sim_run}},
  \code{\link{op_help}}
}
\keyword{package}
