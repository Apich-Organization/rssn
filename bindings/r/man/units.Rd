\name{unit_get}
\alias{unit_get}
\alias{unit_names}
\alias{units_table}
\alias{unit_scale}
\alias{UNITS}
\alias{UNIT_BY_ALIAS}
\alias{constant_get}
\alias{constants_table}
\alias{constant_value}
\alias{CONSTANTS}
\alias{CONSTANT_BY_SYMBOL}
\title{The unit and constant tables}
\usage{
unit_get(name)

unit_names()

units_table()

unit_scale(unit)

constant_get(name)

constants_table()

constant_value(constant)

UNITS

UNIT_BY_ALIAS

CONSTANTS

CONSTANT_BY_SYMBOL
}
\arguments{
  \item{name}{A unit name or alias (\code{"m"}, \code{"metre"}, \code{"ft"}),
    or a constant display name or symbol.}
  \item{unit}{A unit record from \code{unit_get()}.}
  \item{constant}{A constant record from \code{constant_get()}.}
}
\value{
  \code{unit_get()} and \code{constant_get()} return a record (a named list);
  \code{units_table()} and \code{constants_table()} a data frame;
  \code{unit_scale()} and \code{constant_value()} a number.
}
\description{
  The engine's unit and constant tables travel with its manifest, so the
  package carries them exactly: 82 units with their exact rational scale
  (\code{numer / denom * pi^pi}), affine offset, aliases and category, and 30
  physical constants with their exact value, dimension, unit and standard
  uncertainty.
}
\details{
  A scale is exact: the rational part is kept as decimal text in
  \code{numer} / \code{denom} and the \code{pi} power separately, so a
  conversion never loses precision on the way in.  \code{unit_scale()} is the
  numeric approximation.
}
\examples{
unit_get("m")$dimension          # "Length"
unit_get("metre")$name           # "m"
unit_scale(unit_get("ft"))       # 0.3048

length(UNITS)                    # 82
head(units_table())

constant_get("c")$symbol         # "c"
constant_value(constant_get("c"))  # 299792458
constant_get("G")$uncertainty
}
