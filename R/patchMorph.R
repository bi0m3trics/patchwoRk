#' patchMorph
#'
#' @param data_in A RasterLayer or a SpatialPolygonsDataFrame. Map of suitable/non-suitable habitat
#' @param suitThresh A interger. A threshold value over which some organism may perceive the area as
#' suitable habitat (resulting in a binary map of suitable and non-suitable pixels)
#' @param gapThresh A interger. The gap diameter of non-suitable land cover within a habitat patch
#' that should be considered part of the patch if small enough
#' @param spurThresh A interger. The width of a section of narrow, unsuitable edge habitat extending
#' out from a larger, wider patch that is too thin to be considered part of suitable habitat
#' @param suitVals Integer vector. A vector of size = 3 specifying the lower suitability threshold,
#' the upper suitability threshold, and the total number of values to be evaluated.
#' @param gapVals Integer vector. A vector of size = 3 specifying the lower gap threshold, the upper
#' gap threshold, and the total number of values to be evaluated.
#' @param spurVals Integer vector. A vector of size = 3 specifying the lower spur threshold, the upper
#' spur threshold, and the total number of values to be evaluated.
#' @return A RasterLayer or a list of RasterLayers of the same dimensions as data_in where 1's are
#' suitbale habitat and 0's are unsutiable habitat. In the case of PM_Hierarchy, patchMorph returns
#' a list of RasterLayers (one per suitability-gap-spur combination) outcomes, otherwise it returns
#' a single RasterLayer of the single resulting suitability-gap-spur outcome.
#' @references
#' Girvetz EH, and Greco SE. 2007. How to define a patch: a spatial model for hierarchically
#' delineating organism-specific habitat patches. Landscape Ecology 22: 1131-1142.
#' @examples
#' myFile <- system.file("extdata", "mixedconifer.tif", package="patchwoRk")
#' myRas <- raster(myFile)
#'
#' pm.result.single <- patchMorph(data_in = myRas, suitThresh = 1, gapThresh = 2, spurThresh = 2)
#' plot(pm.result.single, main="PatchMorph Results (Gap-2 & Spur-2)")
#'
#' pm.layered.result.map <- patchMorph(data_in = myRas, suitVals = c(0, 1, 2),
#' gapVals = c(2, 6, 3), spurVals = c(2, 6, 3))
#' names(pm.layered.result)
#' plot(pm.layered.result.[1], main=names(pm.layered.result.[1]))
#'
#' @export
patchMorph <- function(data_in, res=-1, suitThresh=-1, gapThresh=-1, spurThresh=-1, suitVals=-1, gapVals=-1, spurVals=-1, proj4=-1,...)
{
  if(length(suitThresh) == 1)
    class(data_in) <- "RasterLayer"
  if(length(suitVals) > 1)
    class(data_in) <- "pmMulti"
  UseMethod("patchMorph", data_in)
}

#' @describeIn patchMorph.RasterLayer Input is a RasterLayer, and only a single suitability, gap, and spur
#' values is specified, for which the only that outcomes is returned
#' @method patchMorph RasterLayer
#' @export
patchMorph.RasterLayer <- function(data_in, suitThresh = 1, gapThresh = 2, spurThresh = 2)
{
  if(!is.numeric(c(suitThresh, gapThresh, spurThresh)))
    stop("suitThresh, gapThresh, and spurThresh must be numeric.")
  if(gapThresh < 2 | spurThresh < 2)
    stop("Gap/Spur threshold is too small! Must be at least twice the raster resolution.")

  ## Set up the crs, the extent, and a NA mask for the original raster
  r.crs <- crs(data_in)
  r.e<-extent(data_in)
  e.mask<-mask(data_in, data_in, maskvalue=0, updatevalue=1)

  ## Get the associated kernels
  gapKernel  <- getCircleKernel(as.integer(gapThresh / 2))
  spurKernel <- getCircleKernel(as.integer(spurThresh / 2))

  ## Set up the suitability threshold, with 0 as non-suitable, and 1 as suitable
  data_in[data_in < suitThresh] <- 0
  data_in[data_in >= suitThresh] <- 1

  ## Get the euclidean distances to suitable habitat, and ensure the extent is the same as original
  xy.suit <- rasterToPoints(data_in, function(x) x == 1) # data in RANN::nn2
  xy.nonsuit <- rasterToPoints(data_in, function(x) x == 0) # query in RANN::nn2
  if(nrow(xy.suit) == 0 | nrow(xy.nonsuit) == 0) {
    stop("Gap parameter resulted there being no euclidean distances to suitable habitat")
  } else {
    euc_dists<-RANN::nn2(data=xy.suit, query=xy.nonsuit, k=1, treetype="kd",
                         searchtype="standard")$nn.dists
    data_in<-raster::rasterFromXYZ(cbind(xy.nonsuit[,1:2], euc_dists/res(data_in)[1]))
    crs(data_in)<-r.crs
    data_in[is.na(values(data_in))] = 0
  }
  data_in<-extend(data_in, r.e, value=NA)
  data_in<-mask(data_in, e.mask)

  ## Apply a focal maximum
  data_in <- focal(data_in, gapKernel, fun=max, na.rm=TRUE, pad=FALSE)
  data_in<-mask(data_in, e.mask)

  cat("Processing gap threshold diameter:", ncol(gapKernel)-1,"pixels\n")
  ## Reclassify based on the gap threshold
  data_in[data_in <= (ncol(gapKernel)+1)/2] <- 1
  data_in[data_in > (ncol(gapKernel)+1)/2] <- 0

  ## Check to make sure there's still non-suitable pixels in the raster, othwewise return data_in
  if( (sum(data_in[values(data_in)==1]) + sum(is.na(values(data_in))) ) == ( nrow(data_in)*ncol(data_in)) ) return(data_in)

  ## Get the euclidean distances to non-suitable habitat, and ensure the extent is the same as
  ## original
  xy.nonsuit <- rasterToPoints(data_in, function(x) x == 0)
  xy.suit <- rasterToPoints(data_in, function(x) x == 1)
  if(nrow(xy.nonsuit) == 0 | nrow(xy.suit) == 0) {
    stop("Gap parameter resulted there being no euclidean distances to non-suitable habitat")
  } else {
    euc_dists<-RANN::nn2(data=xy.nonsuit, query=xy.suit, k=1, treetype="kd",
                         searchtype="standard")$nn.dists
    data_in<-raster::rasterFromXYZ(cbind(xy.suit[,1:2], euc_dists/res(data_in)[1]))
    crs(data_in)<-r.crs
    data_in[is.na(values(data_in))] = 0
  }
  data_in<-extend(data_in, r.e, value=NA)
  data_in<-mask(data_in, e.mask)

  ## Apply a focal maximum
  data_in <- focal(data_in, spurKernel, fun=max, na.rm=TRUE, pad=FALSE)
  data_in <- mask(data_in, e.mask)

  cat("Processing spur threshold diameter:",ncol(spurKernel)-1, "pixels\n")
  ## Reclassify based on the spur threshold
  data_in[data_in <= (ncol(spurKernel)+1)/2] <- 0
  data_in[data_in > (ncol(spurKernel)+1)/2] <- 1

  return(data_in)
}

#' @describeIn patchMorph.pmMulti Input is a pm_Multi, and numeric vecotrs of suitability, gap, and
#' spur values are specified and all possible outcomes are returned.
#' @method patchMorph pmMulti
#' @export
patchMorph.pmMulti <- function(data_in, suitVals=c(0,1,2), gapVals=c(2,8,4), spurVals=c(2,8,4), ...)
{
  if(class(data_in) == "SpatialPolygonsDataFrame")
    data_in <- polygonToRaster(data_in, ...)
  else if(class(data_in) == "matrix" || class(data_in) == "data.frame")
    data_in <- matrixToRaster(data_in, ...)
  else if(class(data_in) != "RasterLayer")
    stop("data_in must be of class RasterLayer, SpatialPolygonsDataFrame, matrix, or data.frame")

  ## As in the paper, we won't consider the lowest suitability class
  suitSeq <- seq(suitVals[1], suitVals[2], (suitVals[2] - suitVals[1]) / (suitVals[3] - 1))[-1]
  gapSeq  <- seq(gapVals[1], gapVals[2], (gapVals[2] - gapVals[1]) / (gapVals[3] - 1))
  spurSeq <- seq(spurVals[1], spurVals[2], (spurVals[2] - spurVals[1]) / (spurVals[3] - 1))

  outList <- list()
  cl <- parallel::makeCluster(pmax(1, parallel::detectCores()-2))
  doParallel::registerDoParallel(cl)

  outList<-
    foreach (i = 1:length(suitSeq)) %:% # nesting operator
    foreach (j = 1:length(gapSeq)) %:% # nesting operator
    foreach (k = 1:length(spurSeq), .packages=c('RANN', 'raster', 'foreach', 'doParallel'), .combine=c,
             .verbose=FALSE,
             .export = c("getCircleKernel", "matrixToRaster","matrixToRaster.default",
                         "matrixToRaster.RasterLayer","patchMorph","patchMorph.pmMulti",
                         "patchMorph.RasterLayer","patchMorph.SpatialPolygonsDataFrame",
                         "patchMorphSummary","polygonToRaster")) %dopar% {
                           patchMorph(data_in, suitThresh = suitSeq[i], gapThresh = gapSeq[j], spurThresh = spurSeq[k])
                         }
  stopCluster(cl)
  outList <- unlist(outList)
  labels <- expand.grid(q=suitSeq,r=gapSeq, s=spurSeq)
  labels <- apply(labels,1, function(x) paste("suit", x[1], "gap", x[2], "spur", x[3], sep=""))
  names(outList)<-labels

  return(outList)
}
