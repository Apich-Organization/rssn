\name{rssn_config}
\alias{rssn_config}
\alias{config_ptr}
\alias{cfg_symbolic}
\alias{cfg_numeric}
\alias{cfg_bind}
\alias{cfg_assume}
\alias{cfg_budget}
\alias{cfg_dimension_mode}
\alias{cfg_bypass_dimensions}
\alias{cfg_dimension}
\alias{cfg_define_base_dimension}
\alias{assumption}
\alias{dimension_mode}
\alias{close.rssn_config}
\alias{print.rssn_config}
\title{The computation configuration}
\usage{
rssn_config()

config_ptr(config)

cfg_symbolic(config)

cfg_numeric(config, tolerance = 1e-12)

cfg_bind(config, symbol, value)

cfg_assume(config, symbol, facts)

cfg_budget(config, max_iterations = 0, max_nodes = 0)

cfg_dimension_mode(config, mode = "advisory")

cfg_bypass_dimensions(config)

cfg_dimension(config, symbol, text)

cfg_define_base_dimension(config, name)

assumption

dimension_mode

\method{close}{rssn_config}(con, ...)
}
\arguments{
  \item{config}{A \code{rssn_config}.}
  \item{con}{A \code{rssn_config}.}
  \item{tolerance}{The tolerance a numeric answer must meet.}
  \item{symbol}{The symbol a binding, an assumption or a dimension is about.}
  \item{value}{The number a symbol is bound to.}
  \item{facts}{An \code{\link{assumption}} bit set, a name such as
    \code{"positive"}, or a character vector of names.}
  \item{max_iterations, max_nodes}{Search limits; zero keeps the current one.}
  \item{mode}{A \code{\link{dimension_mode}} value or name: \code{"off"},
    \code{"advisory"} or \code{"strict"}.}
  \item{text}{The dimension as text: \code{"m/s"}, \code{"Length/Time"},
    \code{"any"}.}
  \item{name}{A base dimension the engine does not have, such as
    \code{"Currency"}.}
  \item{...}{Ignored.}
}
\value{
  Every setter returns \code{config}, invisibly, so a target reads as a pipe.
  \code{rssn_config()} returns an object of class \code{rssn_config};
  \code{config_ptr()} the external pointer the C bridge holds.
}
\description{
  A configuration is the target of a computation, the symbol bindings, the
  assumptions and the dimension declarations.  It is single-threaded, like the
  session it is used with, and is released by \code{close()} or by the garbage
  collector.
}
\details{
  \code{assumption} is a named integer vector of the fact bits the engine
  accepts: \code{real}, \code{nonzero}, \code{nonnegative}, \code{nonpositive},
  \code{integer}, and the composites \code{positive} (7) and \code{negative}
  (11).  \code{dimension_mode} holds \code{off}, \code{advisory} and
  \code{strict}.
}
\examples{
cfg <- rssn_config() |> cfg_numeric(1e-12) |> cfg_bind("x", 0.5)
x <- sym("x")
compute(rssn_diff(sin(x), x), cfg)$value     # 0.877582561890373

cfg_assume(rssn_config(), "n", "integer")
assumption[["integer"]]                      # 16

strict <- rssn_config() |>
    cfg_dimension("x", "m") |>
    cfg_dimension("y", "s") |>
    cfg_dimension_mode("strict")
try(compute(parse_term("x + y"), strict))
}
\seealso{\code{\link{compute}}, \code{\link{dimension}}}
