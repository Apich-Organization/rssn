\name{rssn_encode_npy}
\alias{rssn_encode_npy}
\alias{rssn_decode_npy}
\title{Encode or decode NumPy array bytes}
\usage{
rssn_encode_npy(x, shape = dim(x) \%||\% length(x), session = NULL)
rssn_decode_npy(payload, session = NULL)
}
\arguments{
  \item{x}{A numeric array or a flat row-major vector of samples.}
  \item{shape}{Array extents. Defaults to \code{dim(x)} for an array and
    \code{length(x)} for a vector. Use \code{integer()} for a scalar.}
  \item{payload}{A raw vector containing NumPy array bytes.}
  \item{session}{The RSSN session. Defaults to the package's default session.}
}
\value{
  \code{rssn_encode_npy()} returns a raw vector. \code{rssn_decode_npy()}
  returns an owned binary64 \code{rssn_value} tensor.
}
\description{
  The native codec reads Float32 and Float64 arrays and writes little-endian
  Float64. File access is handled by the caller.
}
