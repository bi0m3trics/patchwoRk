# R
library(terra)

getCircularKernel <- function(radius)
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

# Set parameters
input_raster <- "J:/fragPatch/FtValleyClipCHM.tif"
outDir <- "J:/fragPatch/patchMorph"

# Create lists
suitList <- c(0,2,50)
gapList <- seq(2, 16, by = 2)
spurList <- seq(4, 24, by = 4)

# Load raster
r <- rast(input_raster)

# Loop over suitabilities
for (cntSuit in suitList) {
  suit <- cntSuit
  outDir1 <- file.path(outDir, paste0("suit", suit))
  dir.create(outDir1)
  
  # Loop over gaps
  for (cntGap in gapList) {
    gapDist <- cntGap
    outDir2 <- file.path(outDir1, paste0("gap", gapDist))
    dir.create(outDir2)
    
    # Loop over spurs
    for (cntSpur in spurList) {
      spurDist <- cntSpur
      outDir3 <- file.path(outDir2, paste0("spur", spurDist))
      dir.create(outDir3)
      
      # Reclassify raster based on suitability
      rSuit <- classify(r, matrix(c(-Inf, suit, NA, suit, Inf, 1), ncol=3, byrow=TRUE))
      
      # Calculate Euclidean distance
      rDist <- terra::distance(rSuit)
      
      spurKernal <- getCircularKernel((spurDist/2))
      
      # Apply focal statistics
      rFocal <- focal(rDist, w = spurKernal, fun = max, na.rm = TRUE)
      
      # Reclassify raster based on spur distance
      rSpur <- classify(rFocal, matrix(c(-Inf, spurDist/2, NA, spurDist/2, Inf, 1), ncol=3, byrow=TRUE), right = FALSE)
      
      # Calculate Euclidean distance
      rDist2 <- terra::distance(rSpur)
      
      gapKernal <- getCircularKernel((gapDist/2))
      
      # Apply focal statistics
      rFocal2 <- focal(rDist2, w = gapKernal, fun = max, na.rm = TRUE)
      
      # Reclassify raster based on gap distance
      rGap <- classify(rFocal2, matrix(c(-Inf, gapDist/2, 0, gapDist/2, Inf, 1), ncol=3, byrow=TRUE), right = FALSE)
      
      # Save output raster
      writeRaster(rGap, filename = file.path(outDir3, paste0("su", cntSuit, "sp", cntSpur, "gp", cntGap, ".tif")), overwrite = TRUE)
    }
  }
}

merge_rasters <- function(suitList, gapList, spurList, outDir) {
  # Loop over suitabilities
  for (suit in suitList) {
    outDir1 <- file.path(outDir, paste0("suit", suit))
    
    # Initialize a list to store the rasters
    rasters <- list()
    
    # Loop over gaps
    for (gapDist in gapList) {
      outDir2 <- file.path(outDir1, paste0("gap", gapDist))
      
      # Loop over spurs
      for (spurDist in spurList) {
        outDir3 <- file.path(outDir2, paste0("spur", spurDist))
        
        # Construct the file path
        filePath <- file.path(outDir3, paste0("su", suit, "sp", spurDist, "gp", gapDist, ".tif"))
        
        # Check if the file exists
        if (!file.exists(filePath)) {
          print(paste("File does not exist:", filePath))
          next
        }
        
        # Load the raster
        r <- rast(filePath)
        
        # Add the raster to the list
        rasters[[length(rasters) + 1]] <- r
      }
    }
    
    # Merge the rasters
    rasters_collection <- rast(rasters)
    
    # Sum the rasters
    rSum <- sum(rasters_collection)
    
    # Save the merged raster
    writeRaster(rSum, filename = file.path(outDir1, paste0("psum_su", suit, ".tif")), overwrite = TRUE)
  }
}

merge_rasters(suitList = suitList, gapList = gapList, spurList = spurList, outDir = outDir)
plot(rast("J:/fragPatch/patchMorph/suit2/psum_su2.tif"))
