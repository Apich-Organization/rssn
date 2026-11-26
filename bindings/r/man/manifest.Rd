\name{manifest}
\alias{manifest}
\alias{manifest_json}
\alias{manifest_hash}
\alias{MANIFEST_META}
\alias{rssn_version}
\alias{rssn_library}
\alias{rssn_library_candidates}
\alias{rssn_load}
\alias{json_decode}
\alias{json_encode}
\alias{rssn_condition}
\title{The library, the manifest, and the JSON codec}
\usage{
rssn_version()

rssn_library()

rssn_library_candidates()

rssn_load(path = NULL)

manifest()

manifest_json()

manifest_hash()

json_decode(text)

json_encode(value, auto_unbox = TRUE)

rssn_condition(status, message = NULL)
}
\arguments{
  \item{path}{A file, a directory holding the library, or \code{NULL} to
    restore the automatic search.}
  \item{text}{A length-one JSON string.}
  \item{value}{A list, vector or scalar.}
  \item{auto_unbox}{Whether a length-one vector becomes a scalar.}
  \item{status}{An integer status code.}
  \item{message}{The message; the thread's last error when omitted.}
}
\value{
  \code{rssn_version()} and \code{rssn_library()} return a string;
  \code{rssn_library_candidates()} the search order; \code{manifest()} a named
  list; \code{manifest_hash()} a decimal string; the codec functions a decoded
  object or a string; \code{rssn_condition()} an \code{rssn_error} condition.
}
\description{
  The package does not link against rssn: it loads the shared library at run
  time.  The search order is \code{options(rssn.library = )}, \code{RSSN_LIB},
  the package's own \code{libs/} directory, the build directories of a source
  checkout, and finally the platform loader's own path.

  The manifest is the engine's own description of itself: the operator table,
  the units, the constants and the feature flags.  The bindings are generated
  from it, so \code{OPERATORS} and \code{OPS} cannot drift from the engine.
  \code{manifest()} assembles it from those tables -- a data frame of operators
  rather than a list of lists -- and \code{manifest_json()} hands over the
  engine's own serialisation, which is about 18 MB and slow to parse in R.

  R has no JSON parser in base, and this package deliberately has no
  dependencies, so the codec is included.  \code{json_decode()} answers a named
  list for an object, an unnamed list for an array, a number, a string, a
  logical or \code{NULL}; \code{json_encode()} is the inverse.
}
\details{
  \code{manifest_hash()} is a \code{uint64_t}, which R's double cannot hold
  exactly, so it is returned as decimal text.

  \code{rssn_condition()} builds the condition the C bridge raises for a status
  code, which is what makes \code{tryCatch(expr, rssn_error = function(e)
  e$status)} work.
}
\examples{
rssn_version()
rssn_library()

nrow(manifest()$ops)             # 5411
manifest()$rssn_version
manifest_hash()

json_decode("{\"a\": [1, 2.5, null]}")
json_encode(list(a = c(1, 2), b = TRUE))

tryCatch(parse_term("sin("), rssn_error = function(e) e$status_name)   # "parse"
}
