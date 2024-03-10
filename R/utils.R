getCircleKernel <- function(radius)
{
  kernel_side <- 2 * as.integer(radius) + 1
  kernel_y <- matrix(rep(radius:-radius, kernel_side), ncol=kernel_side)
  kernel_x <- -t(kernel_y)
  kernel   <- matrix(as.matrix(dist(cbind(as.vector(kernel_x), as.vector(kernel_y))))[as.integer((kernel_side^2) / 2) + 1,], ncol=kernel_side)
  kernel[kernel <= radius] <- 0
  kernel[kernel > 0]  <- 1
  kernel <- 1 - kernel
  return(kernel)
}

matrixToRaster.SpatRaster <- function(data_in, ras_match) {
  # Create an empty SpatRaster with the same properties as ras_match
  result_raster <- rast(ras_match)

  # Check if the dimensions of the matrix and the SpatRaster match
  if (nrow(data_in) != nrow(result_raster) || ncol(data_in) != ncol(result_raster)) {
    stop("Dimensions of matrix do not match the dimensions of the SpatRaster template")
  }

  # Assign the matrix values to the SpatRaster
  values(result_raster) <- as.vector(data_in)

  return(result_raster)
}

matrixToRaster <- function(data_in, res, crs = "+proj=utm +zone=12 +ellps=GRS80") {
  nrows <- nrow(data_in)
  ncols <- ncol(data_in)
  xmn <- 0
  xmx <- ncols * res
  ymn <- 0
  ymx <- nrows * res
  r <- rast(nrows = nrows, ncols = ncols, xmin = xmn, xmax = xmx, ymin = ymn, ymax = ymx, crs = crs)
  values(r) <- as.vector(data_in)
  return(r)
}

matrixToRaster <- function(data_in, ras_match, ...)
{
  UseMethod("matrixToRaster", ras_match)
}

is_strictly_increasing <- function(x) {
  all(diff(x) > 0)
}
