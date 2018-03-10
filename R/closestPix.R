#' @title Find the pixel closest to the border of a polygon
#' 
#' @description Given a set of locations (x and y coordinates), find the location of the closest pixel in a raster, which will be inside the border of polygon. This function is designed to account for locations that are outside the extent of the raster.
#' 
#' @param sp A \code{SpatialPolygons} or a \code{SpatialPolygonsDataFrame}
#' @param mesh An \code{inla.mesh} object
#' @param raster A \code{\link{raster}} that includes the explanatory variables to consider for the analysis. This could also be a \code{\link{stack}} or a \code{\link{brick}}.
#' 
#' @importFrom sp SpatialPoints
#' @importFrom raster crs
#' @importFrom rgeos gWithin
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
closestPix <- function(sp, mesh, raster){
  ### Construct SpatialPoints from mesh edges
  loc <- SpatialPoints(coords = mesh$loc[,1:2],
                       proj4string = crs(sp))
  
  ### Select the edges outside polygon
  sel <- (!gWithin(loc, sp, byid = TRUE))
  
  ### Reconstruct SpatialPoints from selected mesh edges
  loc <- SpatialPoints(coords = mesh$loc[sel,1:2],
                       proj4string = crs(sp))

  ### build raster to calculate distance
  distRaster <- raster(xmn = min(mesh$loc[,1])-1, 
                       xmx = max(mesh$loc[,1])+1, 
                       ymn = min(mesh$loc[,2])-1, 
                       ymx = max(mesh$loc[,2])+1,
                       resolution = res(raster))
  
  ### Extend raster to the size of distRaster
  rasterBroad <- extend(raster,extent(distRaster))
  
  ### Calculate distance for the whole area for a single point
  distRasterList <- vector("list", length = sum(sel))
  
  for(i in 1:sum(sel)){
    distRasterize <- rasterize(matrix(mesh$loc[i,1:2],ncol=2),distRaster)
    distRasterList[[i]] <- distance(distRasterize)
  }
  
  ### Build brick
  distRasterBrick <- brick(distRasterList)
  
  ### Crop to raster size
  distRasterBrickCrop <- crop(distRasterBrick, extent(raster))

  ### Mask distRaster with sp
  distMask <- mask(distRasterBrickCrop,sp)
  pixels <- values(distMask)
  
  minPixel <- apply(pixels,2, which.min)
  names(minPixel)<-NULL
  
  ### return
  res <- list(mesh = mesh, meshLoc = loc, meshSel=sel, minPixel = minPixel)
  class(res) <- "closestPix"
  
  return(res)
}
