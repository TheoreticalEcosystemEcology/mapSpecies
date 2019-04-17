#' @title Construct maps for spatial point process models
#' 
#' @description Constructs mean and standard deviation maps for spatial point process models
#' 
#' @param spaceTimeModel An object of class \code{spatialTimePP}.
#' @param dims A vector of length 2 defining the number of pixels in rows and columns of the map.
#' @param type Either "mean", "sd", "0.025quant", "0.5quant", "0.975quant" or "mode". Defines the map to be drawn.
#' @param sPoly A spatial polygon to isolate the region of interest. If none is given, a map is drawn for the entire region covered by the mesh. 
#' 
#' @importFrom INLA inla.mesh.projector
#' @importFrom INLA inla.mesh.project
#' @importFrom INLA inla.stack.index
#' @importFrom raster raster
#' @importFrom raster mask
#' @importFrom raster stack
#'
#' @keywords hplot
#' 
#' @export
#' 
mapSpaceTime <- function(spaceTimeModel, dims, 
                         type = c("mean", "sd", "0.025quant", 
                                "0.5quant", "0.975quant",
                                "mode"), sPoly = NULL){
  ### General check
  if(length(type) > 1){
    stop("Only one type should be defined")
  }
  
  ### Number of spatial edges
  nSpaceEdges <- spaceTimeModel$meshSpace$n
  nTimeEdges <- spaceTimeModel$meshTime$n
  
  ### Define map basis
  mapBasis <- inla.mesh.projector(spaceTimeModel$meshSpace, dims = dims)
  
  ### Find the mesh edges on which predictions should be made
  ID <- inla.stack.index(spaceTimeModel$Stack, tag="pred")$data
  
  ### Calculate prediction
  mapPred <- vector("list", length = nTimeEdges)
  
  for(i in 1:nTimeEdges){
    mapPredBase <- inla.mesh.project(mapBasis, 
                                 spaceTimeModel$model$summary.fitted.values[[type]][ID][1:nSpaceEdges+(i-1)*nSpaceEdges])
    
    mapPred[[i]] <- raster(t(mapPredBase[,ncol(mapPredBase):1]),
                           xmn = min(mapBasis$x), xmx = max(mapBasis$x), 
                           ymn = min(mapBasis$y), ymx = max(mapBasis$y))
  }

  ### Transform map into a raster
  mapRaster <- stack(mapPred)
  
  ### Isolate region of interest
  if(!is.null(sPoly)){
    mapRaster <- mask(mapRaster, sPoly)
  }
  
  ### Return
  return(mapRaster)
}