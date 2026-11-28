\name{call_key}
\alias{call_key}
\title{Call an exact operator overload}
\usage{
call_key(key, ..., session = NULL)
}
\arguments{
\item{key}{An element of \code{op_info(name)$keys}.}
\item{...}{Values accepted by the overload. Use R integer literals such as
\code{2L} for integer parameters and registered callbacks for functions.}
\item{session}{The RSSN session, inferred from the arguments when omitted.}
}
\value{
An owned \code{rssn_value}.
}
\description{
Executes the overload named by the key's rule set, operator name and signature
digest. RSSN checks the supplied values against that overload's contract.
}
