% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/dataSpectrum.R
\name{Mass.Stopping.Power.MeV.cm2.g}
\alias{Mass.Stopping.Power.MeV.cm2.g}
\title{Mass stopping power for all energy bins in a spectrum}
\usage{
Mass.Stopping.Power.MeV.cm2.g(x, stopping.power.source, target.material)
}
\arguments{
\item{x}{Object of class \code{\link{dataSpectrum}} or list of objects of this class}

\item{stopping.power.source}{Descriptor for source of stopping power data (\code{\link[libamtrack]{stopping.power.source}})}

\item{target.material}{Descriptor for target material (\code{\link[libamtrack]{material.no}})}
}
\description{
(List of) vector(s) of mass stopping power values (in the order of particle number and energies)
of the given spectrum/a.
}
\details{
Uses \code{\link[libamtrack]{libamtrack}} stopping power function
}

