#' @title Find the closest pixel in a raster given a set of coordinates
#' 
#' @description Given a set of locations (x and y coordinates), find the location of the closest pixel in a raster. This function is designed to account for locations that are outside the extent of the raster. This function is meant to be used internally.
#' 
#' @param x A vector of X coordinates
#' @param y A vector of Y coordinates
#' @param raster A \code{\link{raster}}. This could also be a \code{\link{stack}} or a \code{\link{brick}}.
#' 
#' @importFrom raster xmin
#' @importFrom raster ymin
#' @importFrom raster xres
#' @importFrom raster yres
#'
closestPix <- function(x, y, raster){
  ### Find dimension of raster
  rasterRowCol <- dim(raster)[1:2]
  
  ### Extract location in conceptual raster (assuming all values were available)
  xPointCon <- round(1 + (x - xmin(raster))/xres(X))
  yPointCon <- round(1 + (y - ymin(raster))/yres(X))
  
  ### Extract location in real raster
  xLoc <- pmax.int(1, pmin.int(xPointCon, rasterRowCol[1]))
  yLoc <- pmax.int(1, pmin.int(yPointCon, rasterRowCol[2]))
  
  ### Organise to fit in the range of the data
  xy <- cbind(xLoc/rasterRowCol[1], yLoc/rasterRowCol[2])
  
  ### Name columns
  colnames(xy) <- c("x","y")
  
  ### return
  return(xy)
}
