#' @title Spatial model for presence-absence data
#'
#' @description Spatial model for presence absence data using INLA. This function is essentially a sophisticated wrapper over \code{inla}
#'
#' @param formula A formula that only relates the response \code{y} and some (or all) of the explanatory variables \code{X}. A paricularity of the is formula is that the response has to be defined as \code{y}.
#' @param spdf A SpatialPointsDataFrame with a vector of 0 and 1 in the data.frame part \code{stack} (or \code{brick}) combining all of the explanatory variables to consider.
#' @param X A raster \code{stack} (or \code{brick}) combining all of the explanatory variables to consider.
#' @param sp A \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame} 
#' @param closestPix An object of class \code{closestPix}
#' @param smooth A single value ranging between 0 and 2 passed to \code{inla.spde2.pcmatern}. It defines the smoothness of the Matern SPDE model. Default is set at 2.
#' @param prior.range A vector of length 2, with (range0, Prange) specifying that P(ρ < ρ_0) = p_ρ, where ρ is the spatial range of the random field. If Prange is NA, then range0 is used as a fixed range value. Default is c(0.05, 0.01).
#' @param prior.sigma A vector of length 2, with (sigma0, Psigma) specifying that P(σ > σ_0) = p_σ, where σ is the marginal standard deviation of the field. If Psigma is NA, then sigma0 is used as a fixed range value.  Default is c(1, 0.01).
#' @param \dots Arguments passed to \code{inla}
#'
#' @export
#' 
#' @keywords models
spatialPA <- function(formula, spdf, X,  closestPix, 
                      smooth = 2,
                      prior.range = c(0.05, 0.01),
                      prior.sigma = c(1, 0.01), ...){
  
  #================
  ### Basic objects
  #================
  ny <- nrow(y)
  nEdges <- closestPix$mesh$n
  
  #==============
  ### Define SPDE
  #==============
  SPDE <- inla.spde2.pcmatern(mesh=closestPix$mesh, 
                              alpha=smooth,
                              prior.range=prior.range,
                              prior.sigma=prior.range)

  #==========================
  ### Define response objects
  #==========================
  ### Coordinates use for the estimation
  xy <- coordinates(spdf)
  
  #===========================
  ### Define covariate objects
  #===========================
  ### Organize data into a data.frame
  ### Note that y here is bogus, it will be removed later
  refData <- data.frame(y = 1,X@data@values)
  colnames(refData)[1] <- names(spdf)
  
  ### Organize X so that it follows the formula
  Xorg <- model.matrix(formula,model.frame(formula, 
                                           data = refData, 
                                           na.action = NULL))[,-1]
  
  ### Construct a brick out of Xorg
  xyXorg <- cbind(coordinates(X),Xorg)
  Xbrick <- rasterFromXYZ(xyXorg)
  
  ### Extract covariate values for model estimation
  meshLoc <- closestPix$mesh$loc[,1:2]
  meshLoc[closestPix$meshSel,] <- coordinates(Xbrick)[closestPix$minPixel,]
  
  locEst <- SpatialPoints(coords = rbind(meshLoc,xy))
  XEst <- extract(Xbrick, locEst)
  
  ### Extract covariate values for model prediction
  locPred <- SpatialPoints(coords = meshLoc)
  XPred <- extract(Xbrick, locPred)
 
  #=====================
  ### Construct A matrix
  #=====================
  ### For estimation
  AEst <- inla.spde.make.A(closestPix$mesh, loc = xy)

  ### For prediction
  APred<-inla.spde.make.A(Mesh)
  
  #====================================
  ### Build stack object for estimation
  #====================================
  ### For estimation
  resp <- unlist(spdf@data)
  names(resp) <- NULL
  resp <- list(y = resp)
  names(resp) <- names(spdf)
  
  StackEst <- inla.stack(data = resp, A = list(1, AEst), 
                         effects = list(list(Intercept = 1, 
                                             X = XEst), 
                                        list(i = 1:nEdges)),
                         tag = "est")

  ### For prediction
  respNA <- list(y = NA)
  names(respNA) <- names(spdf)
  
  StackPred <- inla.stack(data = respNA, A = list(1, APred), 
                          effects = list(list(Intercept = 1, 
                                              X = XPred),
                                         list(i = 1:nEdges)), 
                          tag = "pred")
  
  ### Combine both stack objects
  Stack <- inla.stack(StackEst, StackPred)
  
  #===============
  ### Build models
  #===============
  formule <- formula(y ~ 0 + Intercept + X + f(i, model=SPDE))
  
  model <- inla(formule, family = "binomial", 
                data = inla.stack.data(Stack),
                control.family = list(link = "logit"),
                control.predictor = list(A = inla.stack.A(Stack), 
                                         link = 1),
                E = inla.stack.data(Stack)$e, ...)
  

  ### Return model
  res <- list(model = model, Stack = Stack, mesh = closestPix$mesh)

  class(res) <- "spatialPA"
  return(res)
}