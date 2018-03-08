#' @title Calculate weight for the integration of the point process
#' 
#' @param sp A \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame} object.
#' @param mesh An \code{inla.mesh} object
#' 
#' @importFrom rgeos gIntersects
#' @importFrom rgeos gArea
#' @export
#' 
#' @keywords models
#' 
weightPP <- function(sp, mesh){
  #------------------------------------------------
  ### Construct a voronoid from the triangular mesh
  #------------------------------------------------
  voro <- inla.mesh.dual(mesh)
  
  #--------------------------------------------------------------
  ### Find the intersection between the polygons in the dual mesh
  ### and the location domain
  #--------------------------------------------------------------

  ### Calculate weight
  weight <- numeric()
  
  for(i in 1:length(voro)){
    if(gIntersects(voro[i,], sp)){
      weight[i] <- gArea(gIntersection(voro[i,], sp))
    }else{
      weight[i] <- 0
    }
  }
  
  ### Check to make sure there are integration points with 0 weight
  if(all(weight > 0)){
    stop("There needs to be some weights that are of 0")
  }
  
  ### Return mesh
  res <- list(mesh = mesh, weight = weight)
  class(res) <- "weightPP"
  return(res)
}