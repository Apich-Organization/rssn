\name{value_with_coherent_dimension}
\alias{value_with_coherent_dimension}
\title{Attach a concrete dimension to a coherent Float value}
\usage{value_with_coherent_dimension(value, dimension)}
\arguments{
\item{value}{An RSSN Float or dense Float tensor in coherent units.}
\item{dimension}{A dimension object such as \code{list(base = "Length")}, or \code{"dimless"}.}
}
\value{A separate RSSN value cell with the same tensor shape.}
\description{Attaches a concrete dimension to an existing coherent Float value.}
