\name{rssn_available_backends}
\alias{rssn_available_backends}
\alias{rssn_current_backend}
\alias{rssn_open_backend}
\alias{rssn_release_backend}
\alias{rssn_set_default_backend}
\title{Select an RSSN execution backend}
\usage{
rssn_available_backends(session = NULL)
rssn_current_backend(session = NULL)
rssn_open_backend(name, session = NULL)
rssn_release_backend(backend)
rssn_set_default_backend(name, session = NULL)
}
\arguments{
  \item{session}{The RSSN session. Defaults to the package's default session.}
  \item{name}{A selector returned by \code{rssn_available_backends()}.
    Use \code{"auto"} to restore the native library's automatic policy.}
  \item{backend}{A native handle from \code{rssn_open_backend()}.
    Release makes every copy of the native handle invalid.}
}
\value{
  \code{rssn_available_backends()} returns a character vector of usable selectors.
  \code{rssn_current_backend()} and \code{rssn_set_default_backend()} return the
  effective backend name. \code{rssn_open_backend()} returns an
  \code{rssn_value}. \code{rssn_release_backend()} returns \code{NULL} invisibly.
}
\description{
  Query and select backends in the loaded RSSN library. The default is
  process-wide. An explicit backend handle can be passed to a native
  compilation call.
}
