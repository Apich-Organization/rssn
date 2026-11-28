\name{op_info}
\alias{op_info}
\title{Read one live operator contract}
\usage{
op_info(name)
}
\arguments{
\item{name}{The operator's name or alias.}
}
\value{
A named list with the live manifest entry, including overload signatures,
signature sources, dimension rules, and a \code{keys} list. Each key can be
passed to \code{\link{call_key}}.
}
\description{
Queries the loaded RSSN library for an operator without decoding the entire
manifest.
}
