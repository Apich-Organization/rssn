\name{value_tensor_cells}
\alias{value_tensor_cells}
\title{Construct a shaped tensor from dynamic cells}
\usage{value_tensor_cells(values, shape, session = NULL, exemplar = NULL)}
\arguments{
  \item{values}{A flat vector or list of values in row-major order.}
  \item{shape}{Nonnegative integer extents, outermost axis first.}
  \item{session}{The RSSN session, inferred from the elements when omitted.}
  \item{exemplar}{A value providing the element type of an empty tensor.}
}
\value{A separate RSSN tensor cell.}
\description{Preserves the kinds of existing cells and exact integers. Rank
zero takes one element. A shape containing a zero extent takes no elements and
accepts an exemplar for its element type.}
\examples{
s <- rssn_session()
large <- value_integer_decimal("18446744073709551617", s)
m <- value_tensor_cells(list(large, 2L), c(1L, 2L), s)
value_text(m[[1L]])
close(m)
close(large)
close(s)
}
