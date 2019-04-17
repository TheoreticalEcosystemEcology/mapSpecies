#' @title Prepare data for model construction
#' @name explanaMesh
#' 
#' @description Prepare and organise the data for model prediction. This function is meant to be used as a starting point for all spatial and spatiotemporal model in this package. This function also gathers explanatory variables values associated to the edges within and outside the area of interest.
#' 
#' @param sPoly A \code{\link[sp]{SpatialPolygons}} or a \code{\link[sp]{SpatialPolygonsDataFrame}} defining the sampling region.
#' @param mesh An \code{inla.mesh} object
#' @param X A \code{\link{raster}} that includes the explanatory variables to consider for the analysis. This could also be a \code{\link{stack}} or a \code{\link{brick}}.
#' @param verbose Logical. Whether or not to print on the screen (five times) the number of edges that have been assigned values.
#' 
#' @details 
#' 
#' The time it takes to run this function is directly related to the number of vertex (points) in the mesh; more specifically the ones outside the region of interest.
#' 
#' Also, this function checks the projections of \code{sPoly}, \code{mesh} and \code{X} to make sure they match. Note that if all projections are NAs, the function will run but warnings messages will be printed on screen related to this. 
#' 
#' @return 
#' 
#' An object of class \code{explanaMesh} that includes a list of all the objects used as arguments (except for the verbose call) and \code{Xmesh} the values of the explanatory variables for all edges of the mesh.
#' 
#' @importFrom sp SpatialPoints
#' @importFrom sp bbox
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
#' @importFrom raster xyFromCell
#' @importFrom raster ncell
#' @importFrom raster cellFromXY
#' @importFrom raster xmin
#' @importFrom raster ymin
#' @importFrom raster xmax
#' @importFrom raster ymax
#' @importFrom stats dist
#'
#' @keywords manip
#'
#' @export
explanaMesh <- function(sPoly, mesh, X, verbose = TRUE){
  #=================
  # Check projection
  #=================
  # Extract projection of objects
  projSp <- crs(sPoly)@projargs
  projX <- crs(X)@projargs
  
  if(is.null(mesh$crs)){
    projMesh <- NA
  }else{
    projMesh <- mesh$crs@projargs
  }
  
  proj <- c(projSp, projMesh, projX)
  
  # Check if all projections are NAs
  if(all(is.na(proj))){
    warning("All objects have NA as projection, make sure this is right")
  }else{
    if(any(is.na(proj))){
      stop("At least one object has NA as projection")
    }else{
      if(!((proj[1] == proj[2]) & (proj[1] == proj[3]) & (proj[2] == proj[3]))){
        stop("The projection for at least one of the spatial object is different from the others")
      }
    }
  }

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
           '0' = matrix(xyFromCell(X, ncell(X) - ncol(X) + 1)[rep(1, nrow(Quad[[x]])), ], ncol = 2),
           '1' = matrix(xyFromCell(X, cellFromXY(X, cbind(xmin(X), Quad[[x]][, 2]))), ncol = 2),
           '2' = matrix(xyFromCell(X, 1)[rep(1, nrow(Quad[[x]])), ], ncol = 2),
           '3' = matrix(xyFromCell(X, cellFromXY(X, cbind(Quad[[x]][, 1], ymin(X)))), ncol = 2),
           '4' = matrix(unlist(Quad[[x]]), ncol = 2),
           '5' = matrix(xyFromCell(X, cellFromXY(X, cbind(Quad[[x]][, 1], ymax(X)))), ncol = 2),
           '6' = matrix(xyFromCell(X, ncell(X))[rep(1, nrow(Quad[[x]])), ], ncol = 2),
           '7' = matrix(xyFromCell(X, cellFromXY(X, cbind(xmax(X), Quad[[x]][, 2]))), ncol = 2),
           '8' = matrix(xyFromCell(X, ncol(X))[rep(1, nrow(Quad[[x]])), ], ncol = 2)
           )
        }
      )
  
  # Reorganize to a data.frame 
  new.pts <- unsplit(mapply(function(x, y) {
                            row.names(x) <- row.names(y)
                            as.data.frame(x)
                           }, 
                     new.pts, Quad, SIMPLIFY=FALSE), quad)
  
  # Construct SpatialPoints from mesh edges
  loc <- SpatialPoints(coords = new.pts,
                       proj4string = crs(sPoly))
  
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
                         proj4string = crs(sPoly))
    
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
          print(paste(i," out of ",nlocValNA," edges with NAs considered",
                      sep = ""))
        }
      }
    }
  }
  
  # Convert to data.frame
  locVal <- as.data.frame(locVal)
  
  # Check for factors in X and convert locVal accordingly
  if(class(X) =="RasterBrick"){
    X <- stack(X)
  }
  
  Xfactor <- unlist(lapply(X@layers,function(x) x@data@isfactor))

  if(any(Xfactor)){
    for(i in which(Xfactor)){
      locVal[,i] <- round(locVal[,i])
      locVal[,i] <- as.factor(locVal[,i])
      levels(locVal[,i]) <- X@layers[[i]]@data@attributes[[1]]$ID
    }
  }

  # return
  results <- list(sPoly = sPoly, mesh = mesh, X = X, Xmesh = locVal)
  class(results) <- "explanaMesh"
  
  return(results)
}

