#' patchMorphSummary
#'
#' @param pm_list A list of patchMorph RasterLayers (pm_list).
#' @param suit A integer. The values of suitability over which to provide summaries. If left blank, all
#' values in the pm_laist object will be used.
#' @param gap A integer. The gap values over which to provide summaries. If left blank, all
#' values in the pm_laist object will be used.
#' @param spur A integer. The spur valuesover which to provide summaries. If left blank, all
#' values in the pm_laist object will be used.
#' @references
#' Girvetz EH, and Greco SE. 2007. How to define a patch: a spatial model for hierarchically
#' delineating organism-specific habitat patches. Landscape Ecology 22: 1131-1142.
#' @return A rasterLayer. Sums the members of a patchMorph list (pm_list) from pm_hierarchy
#  case of PatchMorph() where the suitability, gap, and spur thresholds are provided
#' @examples
#' myFile <- system.file("extdata", "mixedconifer.tif", package="patchwoRk")
#' myRas <- raster(myFile)
#'
#' pm.layered.result.map <- patchMorph(data_in = myRas, suitVals = c(0, 1, 2),
#' gapVals = c(2, 6, 3), spurVals = c(2, 6, 3))
#' pm.layered.sum. <- patchMorphSummary(pm.layered.result)
#'
#' plot(myRas, main="Original Raster")
#' plot(pm.layered.sum, main="PatchMorph Multi Results (Gap-2:10 & Spur-2:10)")
#' @export
patchMorphSummary <- function(pm_list, suit=-1, gap=-1, spur=-1)
{
  if(!is.list(pm_list))
    stop("pm_list must be a list")
  if(class(pm_list[[1]]) != "RasterLayer")
    stop("pm_list must be a list of RasterLayers")

  pm_mat_list <- pm_list
  for(i in 1:length(pm_list))
    pm_mat_list[[i]] <- matrix(pm_list[[i]]@data@values, ncol=ncol(pm_list[[i]]), byrow=T)

  if(!is.numeric(c(suit, gap, spur)))
    stop("suit, gap, and spur must be numeric values")

  if(suit > -1)
    pm_mat_list <- pm_mat_list[grep(paste("suit", suit, sep=""), names(pm_mat_list))]
  if(gap > -1)
    pm_mat_list <- pm_mat_list[grep(paste("gap", gap, sep=""), names(pm_mat_list))]
  if(spur > -1)
    pm_mat_list <- pm_mat_list[grep(paste("spur", spur, sep=""), names(pm_mat_list))]

  return(matrixToRaster(Reduce("+", pm_mat_list), pm_list[[1]]))
}
