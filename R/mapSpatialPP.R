#' @title Construct maps for spatial point process models
#' 
#' @description Constructs mean and standard deviation maps for spatial point process models
#' 
#' @param model An object of class \code{inla} obtained from the \code{\link{spatialPP}} function.
#' @param mesh An \code{inla.mesh} object.
#' @param resolution A vector of length 2 defining the number of pixels in rows and columns of the map.
#' @param type Either "mean", "sd", "0.025quant", "0.5quant", "0.975quant" or "mode". Defines the map to be drawn.
#' @param sp A spatial polygon to isolate the region of interest. If none is given, a map is drawn for the entire region covered by the mesh. 
#' 
#' @importFrom INLA inla.mesh.projector
#' @importFrom INLA inla.mesh.project
#' @importFrom INLA inla.stack.index
#' @importFrom raster raster
#' @importFrom raster mask
#' 
#' @export
#' 
mapSpatialPP <- function(model, mesh, resolution, 
                         type = c("mean", "sd", "0.025quant", 
                                  "0.5quant", "0.975quant",
                                  "mode"), sp = NULL){
  ### General check
  if(length(type) > 1){
    stop("Only one type should be defined")
  }
  
  ### Define map basis
  mapBasis <- inla.mesh.projector(mesh, dims = resolution)
  
  ### Find the mesh edges on which predictions should be made
  ID <- inla.stack.index(stk.c, tag="pred")$data
  
  ### Calculate prediction
  mapPred <- inla.mesh.project(mapBasis, 
                               model$summary.fitted.values[[type]][ID])
  
  ### Transform map into a raster
  mapRaster <- raster(t(mapPred[,ncol(mapPred):1]),
                      xmn = min(mapBasis$x), xmx = max(mapBasis$x), 
                      ymn = min(mapBasis$y), ymx = max(mapBasis$y))
  
  ### Isolate region of interest
  if(!is.null(sp)){
    mapMask <- mask(log(mapRaster), RefSP)
  }
  
  ### Return
  return(mapMask)
}