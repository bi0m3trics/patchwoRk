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

matrixToRaster.RasterLayer <- function(data_in, ras_match)
{
  return(raster(data_in, template=ras_match))
}

matrixToRaster.default <- function(data_in, res, proj4=CRS("+proj=utm +zone=12 +ellps=GRS80"))
{
  ras_match <- raster(nrows=nrow(data_in), ncols=ncol(data_in),
                      xmn=0, xmx=ncol(data_in) * res,
                      ymn=0, ymx=nrow(data_in) * res, crs=proj4)
  return(matrixToRaster(data_in, ras_match))
}

matrixToRaster <- function(data_in, ras_match, ...)
{
  UseMethod("matrixToRaster", ras_match)
}
