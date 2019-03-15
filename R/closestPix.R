#' @title Find the pixel closest to the border of a polygon
#' 
#' @description Given a set of locations (x and y coordinates), find the location of the closest pixel in a raster, which will be inside the border of polygon. This function is designed to account for locations that are outside the extent of the raster.
#' 
#' @param sp A \code{SpatialPolygons} or a \code{SpatialPolygonsDataFrame}
#' @param mesh An \code{inla.mesh} object
#' @param X A \code{\link{raster}} that includes the explanatory variables to consider for the analysis. This could also be a \code{\link{stack}} or a \code{\link{brick}}.
#' 
#' @details 
#' 
#' The time it takes to run this function is directly related to the number of vertex (points) in the mesh; more specifically the ones outside the region of interest.
#' 
#' @importFrom sp SpatialPoints
#' @importFrom raster crs
#' @importFrom raster raster
#' @importFrom raster res
#' @importFrom raster extend
#' @importFrom raster extent
#' @importFrom raster rasterize
#' @importFrom raster distance
#' @importFrom raster brick
#' @importFrom raster crop
#' @importFrom raster mask
#' @importFrom raster values
#'
#' @export
closestPix <- function(sp, mesh, X){
  ### Construct SpatialPoints from mesh edges
  loc <- SpatialPoints(coords = mesh$loc[,1:2],
                       proj4string = crs(sp))
  
  ### Extract values from X at loc
  locVal <- extract(X, loc)
  
  #==================================================
  # If there are some elements in locVal that are NAs
  #==================================================
  # Find locVal with NAs
  locValNA <- which(is.na(locVal[,1]))
  nlocValNA <- length(locValNA)
  
  if(nlocValNA > 0){
    loc <- SpatialPoints(coords = coordinates(loc)[locValNA,1:2],
                         proj4string = crs(sp))
    
    for(i in 1:nlocValNA){
      buffer <- 20000
      try(locExtract <- extract(X,
                                loc[i,],
                                buffer = buffer,
                                fun = mean), silent = TRUE)
      
      while(any(is.na(locExtract))){
        locExtract <- extract(X, loc[i,], 
                              buffer = buffer,
                              fun = mean)
        buffer <- buffer + buffer
      }
      
      if(!is.null(locExtract)){
        locVal[locValNA[i],] <- locExtract
      }
      
      locExtract <- NULL
    }
    
    #========================================================
    # If there are still some elements in locVal that are NAs
    #========================================================
    # Find locVal with NAs
    locValNA <- which(is.na(locVal[,1]))
    nlocValNA <- length(locValNA)
    
    if(nlocValNA > 0){
      ### Calculate distance for the whole area for a single point
      distList <- vector("list", length = nlocValNA)
      
      distRaster <- raster(xmn = min(mesh$loc[,1]) - 1, 
                           xmx = max(mesh$loc[,1]) + 1, 
                           ymn = min(mesh$loc[,2]) - 1, 
                           ymx = max(mesh$loc[,2]) + 1, 
                           resolution = res(X))
      
      extentMesh <- extent(min(mesh$loc[,1]) - 1, 
                           max(mesh$loc[,1]) + 1,
                           min(mesh$loc[,2]) - 1,
                           max(mesh$loc[,2]) + 1)
      
      distRaster <- raster(X)
      distRaster <- extend(distRaster, extentMesh)
      distRaster <- crop(distRaster, extentMesh)
      
      for(i in 1:nlocValNA){
        distRasterize <- rasterize(matrix(mesh$loc[locValNA[i],1:2],ncol=2),distRaster)
        distList[[i]] <- distance(distRasterize)
      }
      
      ### Build brick
      distBrick <- brick(distList)
      
      ### Crop to X size
      rasterMask <- brick(X, nl = nlocValNA)
      
      distExtend <- extend(distBrick, extent(X))
      distCrop <- crop(distExtend, extent(X))
      distMask <- mask(distCrop, sp)
      
      pixels <- values(distMask)
      minPixel <- apply(pixels,2, which.min)
      loc <- SpatialPoints(coords = coordinates(X)[minPixel,],
                           proj4string = crs(sp))
      
      for(i in 1:nlocValNA){
        buffer <- 5000
        locExtract <- extract(X,
                              loc[i,],
                              buffer = buffer,
                              fun = mean)
        
        while(any(is.na(locExtract))){
          locExtract <- extract(X, loc[i,], 
                                buffer = buffer,
                                fun = mean)
          buffer <- buffer + buffer
        }
        locVal[locValNA[i],] <- locExtract
      }
    }
  }
  
  ### return
  results <- list(mesh = mesh, Xmesh = locVal)
  class(results) <- "closestPix"
  
  return(results)
}
