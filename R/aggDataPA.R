#' @importFrom deldir deldir
#' @importFrom deldir tile.list
#' @importFrom sp Polygon
#' @importFrom sp Polygons
#' @importFrom sp SpatialPolygons
#' @importFrom sp over
#' 
aggDataPA <- function(xy, meshSpace){
  
  #================
  ### Basic objects
  #================
  nSpaceEdges <- meshSpace$n

  #=========================================
  ### Find Dirichlet tesselation of the mesh
  #=========================================
  ### Build Delaunay and Dirichlet
  DelauDirich <- deldir(meshSpace$loc[,1], meshSpace$loc[,2])
  
  ### Dirichlet tiles
  tiles <- tile.list(DelauDirich)
  nTiles <- length(tiles)
  
  #=========================================================
  ### Construct SpatialPolygons where each polygon is a tile
  #=========================================================
  poly <- vector("list", length = nTiles)
  for(i in 1:nTiles){
    ### Tiles coordinates
    xyTiles <- cbind(tiles[[i]]$x, tiles[[i]]$y)
    nxyTiles <- nrow(xyTiles)
    poly[[i]] <- Polygons(list(Polygon(xyTiles[c(1:nxyTiles, 1),])), i)
  }

  spatialPolys <- SpatialPolygons(poly)
  
  #============================================
  ### Aggregate data points over space and time
  #============================================
  ### Organise spatial data for aggregation
  spaceAgg <- factor(over(SpatialPoints(cbind(xy$x, xy$y)),
                          spatialPolys), levels=1:nTiles)

  ### Aggregate over space
  spaceAggDF <- as.data.frame(table(spaceAgg))
  spaceAggDF <- apply(spaceAggDF, 2, 
                        function(x) as.integer(as.character(x)))  
  
  spaceAggDF[,2] <- ifelse(spaceAggDF[,2] > 0, 1, 0)
  
  ### Return
  res <- as.data.frame(spaceAggDF)
  colnames(res)[1] <- c("space")
    
  return(res)
}