\name{HX.compute.field}
\alias{HX.compute.field}
\title{HX.compute.field}
\description{Computes fluence in voxels of given resolution in a field part for given beam spot grid.}
\usage{HX.compute.field( beam.spot.grid, field.mm = NULL, resolution.mm)}
\arguments{
  \item{beam.spot.grid}{ Data frame holding beam spot positions, foci, and particle numbers. See \link{beam.spot.grid}}.
  \item{field.mm}{Vector of four elements, holding coordinates of field corners in mm (x.min, y.min, x.max, y.max). If NULL, the field size will be automatically chosen such that the entire irradiation field is covered.}
  \item{resolution.mm}{Side length of voxels in mm.}
}
\value{
  Matrix (for the sake of speed).
}

\examples{
#
}
