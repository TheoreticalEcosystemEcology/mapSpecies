#' @title Distance to geometry
#' 
#' @description Calculates for all pixels of a map the distance to a particular geometry. This function is essential a wrapper over the \code{\link{distance}} function of the raster package
#'  
#' @param sp A \code{SpatialPolygons}, \code{SpatialPolygonsDataFrame}, \code{SpatialLines}, \code{SpatialLinesDataFrame}, \code{SpatialPoints} or \code{SpatialPointsDataFrame} object
#' @param raster A \code{\link{raster}} that includes the explanatory variables to consider for the analysis. This could also be a \code{\link{stack}} or a \code{\link{brick}}. 
#' 
#' @importFrom 
#' @importFrom 
#' @importFrom 
#' 
#' @export
distRaster <- function(sp, raster){
  ### Check to make sure projection match
  if(!identical(proj4string(sp),proj4string(raster))){
    stop("The projection for 'sp' is not the same as the one for 'raster'")
  }
  
  ### Rasterize sp object using raster as a reference
  spRaster <- rasterize(sp, raster)
  
  ### Calculate distance to geometry (sp) for all pixels of spRaster
  spRasterDist <- distance(spRaster)
  
  ### Return
  return(spRasterDist)
}