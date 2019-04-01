#' @title Prepare data for model construction
#' @name explanaMesh
#' 
#' @description Prepare and organise the data for model prediction. This function is meant to be used as a starting point for all spatial and spatiotemporal model in this package. This function also gathers explanatory variables values associated to the edges within and outside the area of interest.
#' 
#' @param sp A \code{SpatialPolygons} or a \code{SpatialPolygonsDataFrame}
#' @param mesh An \code{inla.mesh} object
#' @param X A \code{\link{raster}} that includes the explanatory variables to consider for the analysis. This could also be a \code{\link{stack}} or a \code{\link{brick}}.
#' @param verbose Logical. Whether or not to print on the screen (five times) the number of edges that have been assigned values.
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
explanaMesh <- function(sp, mesh, X, verbose = TRUE){
  #============================================================
  # If the edges are not within the raster boundaries associate
  # them with a close by pixel
  #
  # Some of this code was borrows from : 
  # https://stackoverflow.com/questions/26652629/extracting-a-value-from-a-raster-for-a-specific-point-based-on-the-closest-cell
  #============================================================
  # Extract the points
  pts <- mesh$loc[,1:2]
  
  # Find the quadrant around the raster in which the pixel is located
  quad <- colSums(sapply(pts[, 1], '>', bbox(X)[1, ])) * 3 + 
         colSums(sapply(pts[, 2], '>', bbox(X)[2, ]))
  
  Quad <- split(as.data.frame(pts), quad)
  
  # Find new coordinates (by quadrants)
  new.pts <- lapply(names(Quad), function(x) {
    switch(x, 
           '0' = xyFromCell(X, ncell(X) - ncol(X) + 1)[rep(1, nrow(Quad[[x]])), ],
           '1' = xyFromCell(X, cellFromXY(X, cbind(xmin(X), Quad[[x]][, 2]))),
           '2' = xyFromCell(X, 1)[rep(1, nrow(Quad[[x]])), ],
           '3' = xyFromCell(X, cellFromXY(X, cbind(Quad[[x]][, 1], ymin(X)))),
           '4' = {
                  xy <- as.matrix(Quad[[x]])
                  dimnames(xy) <- list(NULL, c('x', 'y'))
                  xy
                 },
           '5' = xyFromCell(X, cellFromXY(X, cbind(Quad[[x]][, 1], ymax(X)))),
           '6' = xyFromCell(X, ncell(X))[rep(1, nrow(Quad[[x]])), ],
           '7' = xyFromCell(X, cellFromXY(X, cbind(xmax(X), Quad[[x]][, 2]))),
           '8' = xyFromCell(X, ncol(X))[rep(1, nrow(Quad[[x]])), ]
           )
        }
      )
  
  # Reorganize to a data.frame 
  new.pts <- unsplit(mapply(function(x, y) {
    row.names(x) <- row.names(y)
    as.data.frame(x)
  }, new.pts, Quad, SIMPLIFY=FALSE), quad)
  
  # Construct SpatialPoints from mesh edges
  loc <- SpatialPoints(coords = new.pts,
                       proj4string = crs(sp))
  
  # Extract values from X at loc
  locVal <- extract(X, loc)
  
  #==================================================
  # If there are some elements in locVal that are NAs
  #==================================================
  # Find locVal with NAs
  locValNA <- sort(unique(which(is.na(locVal), arr.ind = TRUE)[,1]))
  nlocValNA <- length(locValNA)
  
  if(nlocValNA > 0){
    loc <- SpatialPoints(coords = coordinates(loc)[locValNA,1:2],
                         proj4string = crs(sp))
    
    for(i in 1:nlocValNA){
      extX <- extent(X)
      extDist <- dist(matrix(c(extX@xmin,extX@xmax,
                               extX@ymin,extX@ymax),nrow = 2))
      buffer <- extDist/10
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
      
      if(!is.null(locExtract)){
        locVal[locValNA[i],] <- locExtract
      }
      
      locExtract <- NULL
    
      if(verbose){
        if(any(i == round(seq(0, nlocValNA, length = 7)[-c(1,7)]))){
          print(paste(i," out of ",nlocValNA,"elements considered",
                      sep = ""))
        }
      }
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
        distRasterize <- rasterize(matrix(mesh$loc[locValNA[i],
                                                   1:2],
                                          ncol=2),distRaster)
        
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
  
  # return
  results <- list(mesh = mesh, Xmesh = locVal, X = X)
  class(results) <- "explanaMesh"
  
  return(results)
}
