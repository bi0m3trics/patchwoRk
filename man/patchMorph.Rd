% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/patchMorph.R
\name{patchMorph}
\alias{patchMorph}
\title{patchMorph}
\usage{
patchMorph(data_in, res = -1, suitThresh = -1, gapThresh = -1,
  spurThresh = -1, suitVals = -1, gapVals = -1, spurVals = -1,
  proj4 = -1, ...)
}
\arguments{
\item{data_in}{A RasterLayer or a SpatialPolygonsDataFrame. Map of suitable/non-suitable habitat}

\item{suitThresh}{A interger. A threshold value over which some organism may perceive the area as
suitable habitat (resulting in a binary map of suitable and non-suitable pixels)}

\item{gapThresh}{A interger. The gap diameter of non-suitable land cover within a habitat patch
that should be considered part of the patch if small enough}

\item{spurThresh}{A interger. The width of a section of narrow, unsuitable edge habitat extending
out from a larger, wider patch that is too thin to be considered part of suitable habitat}

\item{suitVals}{Integer vector. A vector of size = 3 specifying the lower suitability threshold,
the upper suitability threshold, and the total number of values to be evaluated.}

\item{gapVals}{Integer vector. A vector of size = 3 specifying the lower gap threshold, the upper
gap threshold, and the total number of values to be evaluated.}

\item{spurVals}{Integer vector. A vector of size = 3 specifying the lower spur threshold, the upper
spur threshold, and the total number of values to be evaluated.}
}
\value{
A RasterLayer or a list of RasterLayers of the same dimensions as data_in where 1's are
suitbale habitat and 0's are unsutiable habitat. In the case of PM_Hierarchy, patchMorph returns
a list of RasterLayers (one per suitability-gap-spur combination) outcomes, otherwise it returns
a single RasterLayer of the single resulting suitability-gap-spur outcome.
}
\description{
patchMorph
}
\examples{
myFile <- system.file("extdata", "mixedconifer.tif", package="patchwoRk")
myRas <- raster(myFile)

pm.result.single <- patchMorph(data_in = myRas, suitThresh = 1, gapThresh = 2, spurThresh = 2)
plot(pm.result.single, main="PatchMorph Results (Gap-2 & Spur-2)")

pm.layered.result.map <- patchMorph(data_in = myRas, suitVals = c(0, 1, 2),
gapVals = c(2, 6, 3), spurVals = c(2, 6, 3))
names(pm.layered.result)
plot(pm.layered.result.[1], main=names(pm.layered.result.[1]))

}
\references{
Girvetz EH, and Greco SE. 2007. How to define a patch: a spatial model for hierarchically
delineating organism-specific habitat patches. Landscape Ecology 22: 1131-1142.
}
