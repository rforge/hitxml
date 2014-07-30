\name{beam.spot.grid}
\alias{beam.spot.grid}
\title{beam.spot.grid}
\description{Data frame that holds beam spot positions, foci and particle numbers per spot. It anticipates the structure in the input/output HIT XML file.}

\details{
  Data frame with five columns and as many rows as beam spot positions:
  x.mm --> x position of beam spots in mm
  y.mm --> y position of beam spots in mm
  focus.X.FWHM.mm --> focus (FWHM) in x in mm
  focus.Y.FWHM.mm --> focus (FWHM) in y in mm
  N.particles --> number of particles per beam spots
}
