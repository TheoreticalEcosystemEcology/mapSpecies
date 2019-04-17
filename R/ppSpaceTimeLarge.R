#' @title Spatio-temporal point process model for large data
#' @name ppSpaceTimeLarge
#'
#' @description Spatio-temporal point process model using INLA designed for data that are composed of many occurrences. This function is essentially a sophisticated wrapper over \code{inla}
#'
#' @param formula A formula that only relates the response \code{y} and some (or all) of the explanatory variables \code{X}. A paricularity of the is formula is that the response has to be defined as \code{y}.
#' @param y A 3-columns matrix defining the x and y coordinates of where a species was found and the last column represents the time (generally in days of the year) when a species was found.
#' @param sPoly A \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame} 
#' @param ppWeight An object of class \code{\link{ppWeight}}
#' @param explanaMesh An object of class \code{explanaMesh}
#' @param meshTime An object of class \code{\link{inla.mesh.1d}}
#' @param smooth A single value ranging between 0 and 2 passed to \code{inla.spde2.pcmatern}. It defines the smoothness of the Matern SPDE model. Default is set at 2.
#' @param prior.range A vector of length 2, with (range0, Prange) specifying that P(ρ < ρ_0) = p_ρ, where ρ is the spatial range of the random field. If Prange is NA, then range0 is used as a fixed range value. Default is c(0.05, 0.01).
#' @param prior.sigma A vector of length 2, with (sigma0, Psigma) specifying that P(σ > σ_0) = p_σ, where σ is the marginal standard deviation of the field. If Psigma is NA, then sigma0 is used as a fixed range value.  Default is c(1, 0.01).
#' @param prior.pccor A vector of length 2, with (cor, Pcor) specifying that P(cor > cor_0) = p_cor, where cor is the temporal autocorrelation. Default is c(0.7, 0.7).
#' @param \dots Arguments passed to \code{inla}
#'
#' @return
#' 
#' A list including an \code{\link{inla}} object, an\code{\link{inla.stack}} object, an \code{\link{inla.mesh.2d}} object for the spatial component of the model and an \code{\link{inla.mesh.1d}} object for the temporal component of the model.
#' 
#' @importFrom INLA inla.mesh.fem
#' @importFrom INLA inla.stack.A
#' @importFrom INLA inla.stack.data
#' @importFrom stats model.matrix
#' @importFrom stats model.frame
#' @importFrom Matrix rBind
#'
#' 
#' @keywords models
ppSpaceTimeLarge <- function(formula, y, ppWeight, explanaMesh, 
                             meshTime, smooth = 2,
                             prior.range = c(0.05, 0.01),
                             prior.sigma = c(1, 0.01),
                             prior.pccor = c(0.7, 0.7), ...){
  
  #==============
  ### Basic check
  #==============
  if(as.character(formula[[2]]) != "y"){
    stop("'y' should be used to define the response variable")
  }
  
  ### Check if the mesh in ppWeight and explanaMesh are the same
  if(!identical(ppWeight$mesh$graph, explanaMesh$mesh$graph)){
    stop("'ppWeight' and 'explanaMesh' were constructed using different mesh")
  }
  
  if(colnames(y)[1] != "x" && colnames(y)[2] != "y" && colnames(y)[3] != "t"){
    stop("'y' needs to be formated so that the column names are, in order, 'x', 'y' and 't'")
    ### This will have to be changed to a STF object in the future ###
    ### This will have to be changed to a STF object in the future ###
    ### This will have to be changed to a STF object in the future ###
    ### This will have to be changed to a STF object in the future ###
    ### This will have to be changed to a STF object in the future ###
    ### This will have to be changed to a STF object in the future ###
  }
  y <- as.data.frame(y)
  
  #================
  ### Basic objects
  #================
  ny <- nrow(y)
  nSpaceEdges <- ppWeight$mesh$n
  nTimeEdges <- meshTime$n
  
  #======================================
  ### Aggregate spatial and temporal data
  #======================================
  spaceTimeAgg <- aggData(y, ppWeight$mesh, meshTime)
  
  #==============
  ### Define SPDE
  #==============
  SPDE <- inla.spde2.pcmatern(mesh=ppWeight$mesh, alpha=smooth,
                              prior.range=prior.range,
                              prior.sigma=prior.range)
  
  #==========================
  ### Define response objects
  #==========================
  ### Pseudo-absences are the number of edges on the mesh
  ### Occurences are the number of points
  yPP <- spaceTimeAgg$Freq
  #---------------------------------------------------------------
  ### weight associated to pseudo-absences (w) and occurrences (0)
  #---------------------------------------------------------------
  #________________________________________________________________
  ### Define spatiotemporal volume to distribute weight across time 
  #________________________________________________________________
  timeWeight <- diag(inla.mesh.fem(meshTime)$c0)
  
  ePP <- ppWeight$weight[spaceTimeAgg$space] * timeWeight[spaceTimeAgg$time]
  
  #===========================
  ### Define covariate objects
  #===========================
  ### Organize data into a data.frame
  ### Note that y here is bogus, it will be removed later
  refData <-  data.frame(y = 1,values(explanaMesh$X))
  
  ### Organize X so that it follows the formula
  Xorg <- model.matrix(formula,model.frame(formula, 
                                           data = refData, 
                                           na.action = NULL))[,-1]
  
  ### Construct a brick out of Xorg
  xyXorg <- cbind(coordinates(explanaMesh$X),Xorg)
  Xbrick <- rasterFromXYZ(xyXorg)
  
  ### Extract covariate values for model estimation
  meshLoc <- ppWeight$mesh$loc[,1:2]
  meshLocBase <- meshLoc
  
  for(i in 1:(nTimeEdges-1)){
    meshLoc <- rbind(meshLoc,meshLocBase)
  }
  
  locEst <- SpatialPoints(coords = meshLoc)
  XEst <- extract(Xbrick, locEst)

  #===========================
  ### Define projection matrix
  #===========================
  ### For inference
  ProjInfer <- inla.spde.make.A(ppWeight$mesh, 
                                loc = ppWeight$mesh$loc[spaceTimeAgg$space,],
                                group = spaceTimeAgg$time,
                                group.mesh = meshTime)
  
  IDSpaceTime <- inla.spde.make.index('i', nSpaceEdges, n.group=nTimeEdges)
  
  #=====================
  ### Build stack object
  #=====================
  StackEst <- inla.stack(data = list(y = yPP, e = ePP), 
                         A = list(1, ProjInfer), 
                         effects = list(list(Intercept = 1, 
                                             X = XEst),
                                        IDSpaceTime), 
                         tag = "est")
  
  StackPred <- inla.stack(data = list(y = NA, e = NA),
                          A = list(1, ProjInfer), 
                          effects = list(list(Intercept = 1, 
                                              X = XEst), 
                                         IDSpaceTime),
                          tag = "pred")
  
  Stack <- inla.stack(StackEst, StackPred)

  #===============
  ### Build models
  #===============
  pcTimeAutoCor <- list(prior='pccor1', param=prior.pccor)
  formule <- formula(y ~ 0 + Intercept + X + f(i, model=SPDE, 
                                               group = i.group, 
                                               control.group = list(model = "ar1",
                                                                    hyper = list(theta = pcTimeAutoCor))))
  
  
  model <- inla(formule, family = "poisson", 
                data = inla.stack.data(Stack),
                control.predictor = list(A = inla.stack.A(Stack), 
                                         link = 1),
                E = inla.stack.data(Stack)$e,
                control.inla=list(strategy='gaussian'), ...)
  
  ### Return model
  res <- list(model = model, Stack = Stack, 
              meshSpace = explanaMesh$mesh, meshTime = meshTime)
  
  class(res) <- "ppSpaceTime"
  
  return(res)
}