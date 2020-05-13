#' @title Calculate weight for the integration of the point process
#' 
#' @description This function calculates the weight that should be given to the 
#' 
#' @param sPoly A \code{SpatialPolygons*} object.
#' @param mesh An \code{inla.mesh} object
#' 
#' @return 
#' 
#' An object of class \code{ppWeight} which returns a vector of weights. 
#' 
#' @importFrom rgeos gIntersects
#' @importFrom rgeos gArea
#' @importFrom rgeos gIntersection
#'
#' @export
#' 
#' @keywords models
#' 
ppWeight <- function(sPoly, mesh){
  #-------------------------------------------------
  ### Construct a dual mesh from the triangular mesh
  #-------------------------------------------------
  dmesh <- inla.mesh.dual(mesh)
  crs(dmesh) <- mesh$crs
  
  #--------------------------------------------------------------
  ### Find the intersection between the polygons in the dual mesh
  ### and the location domain
  #--------------------------------------------------------------

  ### Calculate weight
  weight <- numeric()
  
  for(i in 1:length(dmesh)){
    if(gIntersects(dmesh[i,], sPoly)){
      weight[i] <- gArea(gIntersection(dmesh[i,], sPoly))
      print(i)
    }else{
      weight[i] <- 0
    }
  }
  
  ### Check to make sure there are integration points with 0 weight
  if(all(weight > 0)){
    stop("There needs to be some weights that are of 0")
  }
  
  ### Return mesh
  attributes(weight) <- mesh
  
  class(weight) <- "ppWeight"
  
  return(weight)
}