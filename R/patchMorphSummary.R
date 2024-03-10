#' patchMorphSummary
#'
#' @param pm_list A list of patchMorph SpatRasters (pm_list) all in the same
#' units and with the same resolution.
#' @param suit A integer. The values of suitability over which to provide
#' summaries. If left blank, all values in the suit object will be used.
#' Assumed to be increasing and in the same units as data_in
#' @param gap A integer. The gap values over which to provide summaries. If left
#' blank, all values in the gap object will be used. Assumed to be increasing and
#' in the same units as data_in
#' @param spur A integer. The spur valuesover which to provide summaries. If
#' left blank, all values in the spur object will be used. Assumed to be
#' increasing and in the same units as data_in
#' @references
#' Girvetz EH, and Greco SE. 2007. How to define a patch: a spatial model for hierarchically
#' delineating organism-specific habitat patches. Landscape Ecology 22: 1131-1142.
#' @return A rasterLayer. Sums the members of a patchMorph list (pm_list) from pm_hierarchy
#  case of PatchMorph() where the suitability, gap, and spur thresholds are provided
#' @examples
#' myFile <- system.file("extdata", "mixedconifer.tif", package="patchwoRk")
#' myRas <- rast(myFile)
#'
#' pm.layered.result <- patchMorph(data_in = myRas, buffer = max(c(2, 6, 3)),
#' suitVals = c(0, 1, 2), gapVals = c(2, 6, 3),
#' spurVals = c(2, 6, 3), verbose = TRUE)
#' pm.layered.sum <- patchMorphSummary(pm.layered.result)
#'
#' plot(myRas, main="Original Raster")
#' plot(pm.layered.sum, main="PatchMorph Multi Results (Gap-2:10 & Spur-2:10)")
#' @export
patchMorphSummary <- function(pm_list, suit=-1, gap=-1, spur=-1)
{
  if (!is.list(pm_list))
    stop("pm_list must be a list")

  if (!all(sapply(pm_list, inherits, "SpatRaster")))
    stop("pm_list must be a list of SpatRasters")

  pm_mat_list <- vector("list", length(pm_list))
  for (i in seq_along(pm_list)) {
    # Convert SpatRaster to matrix
    pm_mat_list[[i]] <- terra::as.matrix(pm_list[[i]])
  }

  if (!is.numeric(c(suit, gap, spur)))
    stop("suit, gap, and spur must be numeric values")

  # Filtering based on suit, gap, and spur values if provided
  if (suit > -1)
    pm_mat_list <- pm_mat_list[grep(paste("suit", suit, sep=""), names(pm_mat_list))]
  if (gap > -1)
    pm_mat_list <- pm_mat_list[grep(paste("gap", gap, sep=""), names(pm_mat_list))]
  if (spur > -1)
    pm_mat_list <- pm_mat_list[grep(paste("spur", spur, sep=""), names(pm_mat_list))]

  # Sum the matrices and convert back to SpatRaster
  summed_matrix <- Reduce(`+`, pm_mat_list)
  summed_raster <- terra::rast(pm_list[[1]])
  terra::values(summed_raster) <- summed_matrix

  return(summed_raster)
}
