\name{HX.cost.function}
\alias{HX.cost.function}
\title{HX.cost.function}
\description{Computes cost function to be optimized by HX.optimize.field. It depends on the homogeneity of fluence and penalized heavily particle numbers below the minimum limit.}
\usage{HX.cost.function( d.mm, focus.FWHM.mm, fluence.cm2, field.size.mm, N.min)}
\arguments{
  \item{d.mm}{Value for distance d (in mm).}
  \item{focus.FWHM.mm}{ FWHM size of double-Gaussian beam spot (in mm).}
  \item{field.size.mm}{ Side length (in mm) of the square field in which homgenous fluence shall be obtained. Actual irradiation field will be bigger.}
  \item{fluence.cm2}{Particle fluence (in 1/cm2) within the homogenous field area.}
  \item{N.min}{Minimal number of particle per raster spot. Given by machine constraints.}
}
\value{
  Single number.
}

\examples{
#
}
