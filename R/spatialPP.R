#' @title Spatial point process model
#'
#' @description Spatial point process model using INLA. This function is essentially a sophisticated wrapper over \code{inla}
#'
#' @param formula A formula that only relates the response \code{y} and some (or all) of the explanatory variables \code{X}. A paricularity of the is formula is that the response has to be defined as \code{y}.
#' @param y A 2-columns matrix defining the location of where a species was found.
#' @param X A raster \code{stack} (or \code{brick}) combining all of the explanatory variables to consider.
#' @param weightPP An object of class \code{weightPP}
#' @param closestPix An object of class \code{closestPix}
#' @param smooth A single value ranging between 0 and 2 passed to \code{inla.spde2.pcmatern}. It defines the smoothness of the Matern SPDE model. Default is set at 2.
#' @param prior.range A vector of length 2, with (range0, Prange) specifying that P(ρ < ρ_0) = p_ρ, where ρ is the spatial range of the random field. If Prange is NA, then range0 is used as a fixed range value. Default is c(0.05, 0.01).
#' @param prior.sigma A vector of length 2, with (sigma0, Psigma) specifying that P(σ > σ_0) = p_σ, where σ is the marginal standard deviation of the field. If Psigma is NA, then sigma0 is used as a fixed range value.  Default is c(1, 0.01).
#' @param \dots Arguments passed to \code{inla}
#'
#' @return
#' 
#' A list including an \code{\link{inla}} object, an \code{\link{inla.mesh}} object.
#'
#' @importFrom INLA inla.spde2.pcmatern
#' @importFrom sp coordinates
#' @importFrom raster rasterFromXYZ
#' @importFrom sp SpatialPoints
#' @importFrom raster extract
#' @importFrom INLA inla.spde.make.A
#' @importFrom Matrix Diagonal
#' @importFrom Matrix rBind
#' @importFrom INLA inla.stack
#' @importFrom INLA inla
#' 
#' @export
#' 
#' @keywords models
spatialPP <- function(formula, y, X, weightPP, closestPix, 
                      smooth = 2,
                      prior.range = c(0.05, 0.01),
                      prior.sigma = c(1, 0.01), ...){

  #==============
  ### Basic check
  #==============
  if(as.character(formula[[2]]) != "y"){
    stop("'y' should be used to define the response variable")
  }
  
  ### Check if the mesh in weightPP and closestPix are the same
  if(!identical(weightPP$mesh$graph, closestPix$mesh$graph)){
    stop("'weightPP' and 'closestPix' were constructed using different mesh")
  }
  
  #================
  ### Basic objects
  #================
  ny <- nrow(y)
  nEdges <- weightPP$mesh$n

  #==============
  ### Define SPDE
  #==============
  SPDE <- inla.spde2.pcmatern(mesh=weightPP$mesh, alpha=smooth,
                              prior.range=prior.range,
                              prior.sigma=prior.range)

  #==========================
  ### Define response objects
  #==========================
  ### Pseudo-absences are the number of edges on the mesh
  ### Occurences are the number of points
  yPP <- rep(0:1, c(nEdges, ny))
  
  ### weight associated to pseudo-absences (w) and occurrences (0)
  ePP <- c(weightPP$weight, rep(0, ny))
  
  #===========================
  ### Define covariate objects
  #===========================
  ### Organize data into a data.frame
  ### Note that y here is bogus, it will be removed later
  refData <-  data.frame(y = 1,values(X))
  
  ### Organize X so that it follows the formula
  Xorg <- model.matrix(formula,model.frame(formula, 
                                           data = refData, 
                                           na.action = NULL))[,-1]
  
  ### Construct a brick out of Xorg
  xyXorg <- cbind(coordinates(X),Xorg)
  Xbrick <- rasterFromXYZ(xyXorg)

  ### Extract covariate values for model estimation
  meshLoc <- weightPP$mesh$loc[,1:2]
  meshLoc[closestPix$meshSel,] <- coordinates(Xbrick)[closestPix$minPixel,]

  locEst <- SpatialPoints(coords = rbind(meshLoc,y))
  XEst <- extract(Xbrick, locEst)
  
  ### Extract covariate values for model prediction
  locPred <- SpatialPoints(coords = meshLoc)
  XPred <- extract(Xbrick, locPred)
  
  #===========================
  ### Define projection matrix
  #===========================
  ### For inference
  ProjInfer <- inla.spde.make.A(weightPP$mesh, y)
  
  ### For integration
  ProjInter <- Diagonal(nEdges, rep(1, nEdges))
  
  ### Combine both projection matrix
  A <- rBind(ProjInter, ProjInfer)
  
  #=====================
  ### Build stack object
  #=====================
  StackEst <- inla.stack(data = list(y = yPP, e = ePP), 
                         A = list(1, A), 
                         effects = list(list(Intercept = 1, 
                                             X = XEst), 
                                        list(i = 1:nEdges)), 
                         tag = "est")
  
  StackPred <- inla.stack(data = list(y = NA, e = NA),
                              A = list(1, ProjInter), 
                              effects = list(list(Intercept = 1, 
                                                  X = XPred), 
                                             list(i = 1:nEdges)),
                              tag = "pred")
  
  Stack <- inla.stack(StackEst, StackPred)
  
  #===============
  ### Build models
  #===============
  formule <- formula(y ~ 0 + Intercept + X + f(i, model=SPDE))
  
  model <- inla(formule, family = "poisson", 
                data = inla.stack.data(Stack),
                control.predictor = list(A = inla.stack.A(Stack), 
                                         link = 1),
                E = inla.stack.data(Stack)$e, ...)
  
  ### Return model
  res <- list(model = model, Stack = Stack, mesh = closestPix$mesh)
  
  class(res) <- "spatialPP"
  return(res)
}
