\name{value_integer_decimal}
\alias{value_integer_decimal}
\title{Construct an exact integer from signed decimal digits}
\usage{
value_integer_decimal(digits, session = NULL)
}
\arguments{
\item{digits}{One character string with decimal digits and an optional sign.}
\item{session}{The RSSN session.}
}
\value{An owned \code{rssn_value} containing an exact integer.}
\description{
Constructs an integer cell directly from its decimal spelling. Read the same
integer back with \code{value_text()}.
}
