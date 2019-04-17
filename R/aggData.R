#' @importFrom deldir deldir
#' @importFrom deldir tile.list
#' @importFrom sp Polygon
#' @importFrom sp Polygons
#' @importFrom sp SpatialPolygons
#' @importFrom sp over
#' 
aggData <- function(xyt, meshSpace, meshTime=NULL){
  
  #================
  ### Basic objects
  #================
  nSpaceEdges <- meshSpace$n
  
  if(!is.null(meshTime)){
    nTimeEdges <- meshTime$n
  }
  
  #=========================================
  ### Find Dirichlet tesselation of the mesh
  #=========================================
  ### Build Delaunay and Dirichlet
  DelauDirich <- deldir(meshSpace$loc[,1], meshSpace$loc[,2],
                        suppressMsge = TRUE)
  
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
  spaceAgg <- factor(over(SpatialPoints(cbind(xyt$x, xyt$y)),
                          spatialPolys), levels=1:nTiles)

  ### Organise temporal data for aggregation
  if(is.null(meshTime)){
    ### Aggregate over space
    spaceAggDF <- as.data.frame(table(spaceAgg))
    spaceAggDF <- apply(spaceAggDF, 2, 
                          function(x) as.integer(as.character(x)))  
    
    ### Return
    res <- as.data.frame(spaceAggDF)
    colnames(res)[1] <- c("space")
  }else{
    timeEdges <- sort(c(meshTime$loc[c(1,nTimeEdges)],
                        meshTime$loc[2:nTimeEdges-1]/2 + 
                          meshTime$loc[2:nTimeEdges]/2))
    
    timeAgg <- factor(findInterval(xyt$t, timeEdges), 
                      levels = 1:nTimeEdges)
    
    ### Aggregate over space and time
    spaceTimeAgg <- as.data.frame(table(spaceAgg, timeAgg))
    spaceTimeAgg <- apply(spaceTimeAgg, 2, 
                          function(x) as.integer(as.character(x)))  
  
    ### Return
    res <- as.data.frame(spaceTimeAgg)
    colnames(res)[1:2] <- c("space", "time")
  }
  
  return(res)
}