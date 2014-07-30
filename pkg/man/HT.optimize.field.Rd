\name{HX.optimize.field}
\alias{HX.optimize.field}
\title{HX.optimize.field}
\description{Optimizes distance between raster spot positions to yield a square-shaped inner field area of optimal homogeneity. The beam is assumed to be double-gaussien, the distance between the spot positions is same in x and y.}
\usage{HX.optimize.field( focus.FWHM.mm, field.size.mm, fluence.cm2, N.min, d.start.mm, plot = FALSE)}
\arguments{
  \item{focus.FWHM.mm}{ FWHM size of double-Gaussian beam spot (in mm).}
  \item{field.size.mm}{ Side length (in mm) of the square field in which homgenous fluence shall be obtained. Actual irradiation field will be bigger.}
  \item{fluence.cm2}{Particle fluence (in 1/cm2) within the homogenous field area.}
  \item{N.min}{Minimal number of particle per raster spot. Given by machine constraints.}
  \item{d.start.mm}{Start value for distance d (in mm), a reasonable value is focus / 3.}
}
\value{
  It is also adviced to follow the optimization routine output to check if it works ok.
}

\examples{
#
}
