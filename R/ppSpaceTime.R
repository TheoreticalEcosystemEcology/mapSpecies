#' @title Spatio-temporal point process model
#' @name ppSpaceTime
#'
#' @description Spatio-temporal point process model using INLA. This function is essentially a sophisticated wrapper over \code{inla}
#'
#' @param formula A formula that only relates the response \code{y} and some (or all) of the explanatory variables \code{X}. A paricularity of the is formula is that the response has to be defined as \code{y}.
#' @param ST An object of class \code{\link[spacetime]{STF}}*, \code{\link[spacetime]{STI}}* or \code{\link[spacetime]{STS}}*.
#' @param ppWeight An object of class \code{\link{ppWeight}}
#' @param explanaMesh An object of class \code{dataPrep}
#' @param meshTime An object of class \code{\link{inla.mesh.1d}}
#' @param timeRes A character string defining the temporal resolution to use to perform the analysis. The character string choices are given in the Details section of \code{\link{DateTimeClasses}} help.
#' @param offset A character string defining the explanatory variable in \code{explanaMesh$X} to use as offset.
#' @param smooth A single value ranging between 0 and 2 passed to \code{inla.spde2.pcmatern}. It defines the smoothness of the Matern SPDE model. Default is set at 2.
#' @param prior.range A vector of length 2, with (range0, Prange) specifying that P(ρ < ρ_0) = p_ρ, where ρ is the spatial range of the random field. If Prange is NA, then range0 is used as a fixed range value. Default is c(0.05, 0.01).
#' @param prior.sigma A vector of length 2, with (sigma0, Psigma) specifying that P(σ > σ_0) = p_σ, where σ is the marginal standard deviation of the field. If Psigma is NA, then sigma0 is used as a fixed range value.  Default is c(1, 0.01).
#' @param prior.pccor A vector of length 2, with (cor, Pcor) specifying that P(cor > cor_0) = p_cor, where cor is the temporal autocorrelation. Default is c(0.7, 0.7).
#' @param many Logical. Whether the data in \code{sPoints} is large or not. See details. Default is \code{FALSE}.
#' @param \dots Arguments passed to \code{inla}
#'
#' @details 
#'
#'  If the argument \code{many = TRUE}, the estimation and the prediction will be carried out solely at the mesh edges, whereas when \code{many = FALSE} the estimation will be carried out at the mesh edges and at the sampled location. When the number of samples is very large (e.g. tens of thousands of samples or more) using \code{many = TRUE} can be much more computationally efficient. However, there is a precision trade-off. When \code{many = TRUE}, each sample is associated to an edge and the model is constructed using the number of samples associated to an edge as an importance value. In doing so, some precision is lost at the expense of speed. 
#'
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
#' @importFrom zoo index 
#' 
#' @export
#' 
#' @keywords models
ppSpaceTime <- function(formula, 
                        ST, 
                        ppWeight, 
                        explanaMesh, 
                        meshTime, 
                        timeRes, 
                        offset = NULL, 
                        smooth = 2,
                        prior.range = c(0.05, 0.01),
                        prior.sigma = c(1, 0.01),
                        prior.pccor = c(0.7, 0.7),
                        many = FALSE, 
                        ...){
  
  #==============
  ### Basic check
  #==============
  if(as.character(formula[[2]]) != "y"){
    stop("'y' should be used to define the response variable")
  }
  
  ### Check if the mesh in ppWeight and explanaMesh are the same
  if(!identical(attributes(ppWeight)$graph, explanaMesh$mesh$graph)){
    stop("'ppWeight' and 'explanaMesh' were constructed using different mesh")
  }
  
  # Check class of ST
  if(!any(class(ST) == c("STS", "STI", "STF", "STSDF", "STIDF", "STFDF"))){
    stop("'ST' needs to be an object of class ST*")
  }
  
  # Make sure timeRes has the proper structure
  if(!any(timeRes == c("all", "sec", "min", "hour",
                       "mday", "mon", "year", "wday",
                       "yday"))){
    stop("The first entry of 'timeRes' needs to be either one of
         'all', 'sec', 'min', 'hour','mday', 'mon', 'year', 'wday','yday'")
  }
  
  if(!is.null(offset)){
    explanForm <- as.character(formula)[3]
    if(length(grep(offset,explanForm)) > 0){
      stop("'offset' is also in 'formula'")
    }
  }
  
  #================
  ### Basic objects
  #================
  ny <- length(ST)
  nSpaceEdges <- explanaMesh$mesh$n
  nTimeEdges <- meshTime$n
  
  #==============
  ### Define SPDE
  #==============
  SPDE <- inla.spde2.pcmatern(mesh=explanaMesh$mesh, 
                              alpha=smooth,
                              prior.range=prior.range,
                              prior.sigma=prior.range)
  
  #==========================
  ### Define response objects
  #==========================
  ### Pseudo-absences are the number of edges on the mesh
  ### Occurences are the number of points
  yPP <- rep(0:1, c(nSpaceEdges * nTimeEdges, ny))
  
  #---------------------------------------------------------------
  ### weight associated to pseudo-absences (w) and occurrences (0)
  #---------------------------------------------------------------
  #________________________________________________________________
  ### Define spatiotemporal volume to distribute weight across time 
  #________________________________________________________________
  spaceTimeWeight <- rep(ppWeight, nTimeEdges) * 
                     rep(diag(inla.mesh.fem(meshTime)$c0), 
                         nSpaceEdges)
  
  ePP <- c(spaceTimeWeight, rep(0, ny))
  
  #===========================
  ### Define covariate objects
  #===========================
  ### Organize data into a data.frame
  refData <- as.data.frame(values(explanaMesh$X))
  
  Xfactor <- unlist(lapply(explanaMesh$Xmesh, is.factor))
  if(any(Xfactor)){
    for(i in which(Xfactor)){
      refData[,i] <- as.factor(refData[,i])
      levels(refData[,i]) <- levels(explanaMesh$Xmesh[,i])
    }
  }
  
  refData <- data.frame(y = 1, refData)

  ### Organize X so that it follows the formula
  Xorg <- model.matrix(formula,model.frame(formula, 
                                           data = refData, 
                                           na.action = NULL))[,-1]
  
  ### Construct a brick out of Xorg
  xyXorg <- cbind(coordinates(explanaMesh$X),Xorg)
  Xbrick <- rasterFromXYZ(xyXorg)
  names(Xbrick) <- colnames(Xorg)

  if(!is.null(offset)){
    Xbrick <- addLayer(Xbrick, explanaMesh$X[[offset]])
  }
  
  ### Extract covariate values for model estimation
  meshLoc <- explanaMesh$mesh$loc[,1:2]
  meshLocBase <- meshLoc
  
  for(i in 1:(nTimeEdges-1)){
    meshLoc <- rbind(meshLoc,meshLocBase)
  }
  
  xyMesh <- rbind(meshLoc,
                  coordinates(ST))
  rownames(xyMesh) <- 1:nrow(xyMesh)
  
  locEst <- SpatialPoints(coords = xyMesh)
  XEst <- extract(Xbrick, locEst)
  XEst <- as.data.frame(cbind(Intercept = 1, XEst))
  
  ### Extract covariate values for model prediction
  locPred <- SpatialPoints(coords = explanaMesh$mesh$loc[,1:2])
  
  prefData <- data.frame(y = 1, explanaMesh$Xmesh)

  XPred <- model.matrix(formula,model.frame(formula,
                                            data = prefData,
                                            na.action = NULL))[,-1]
  
  if(!is.null(offset)){
    XPred <- cbind(XPred, explanaMesh$Xmesh[,offset])
    colnames(XPred)[ncol(XPred)] <- offset
  }
  
  XPred <- as.data.frame(XPred)
  
  # Define temporal object
  time <- as.POSIXlt(index(ST@time))
  timeSel <- time[,timeRes]
  
  #===========================
  ### Define projection matrix
  #===========================
  ### For inference
  ProjInfer <- inla.spde.make.A(explanaMesh$mesh, 
                                loc = coordinates(ST), 
                                n.group = length(meshTime$n),
                                group = timeSel, 
                                group.mesh = meshTime)

  IDSpaceTime <- inla.spde.make.index('i', nSpaceEdges, n.group=nTimeEdges)
  
  ### For integration
  ProjInter <- Diagonal(nSpaceEdges * nTimeEdges)
  
  ### Combine both projection matrix
  A <- rbind(ProjInter, ProjInfer)
  
  #=====================
  ### Build stack object
  #=====================
  Stack <- inla.stack(data = list(y = yPP, e = ePP), 
#  StackEst <- inla.stack(data = list(y = yPP, e = ePP), 
                         A = list(A, 1), 
                         effects = list(IDSpaceTime,
                                        list(Intercept = 1, 
                                             XEst)), 
                         tag = "est")
  
#  StackPred <- inla.stack(data = list(y = NA, e = NA),
#                          A = list(ProjInter, 1), 
#                          effects = list(IDSpaceTime,
#                                         list(Intercept = 1, 
#                                              X = XPredTime)),
#                          tag = "pred")
  
#  Stack <- inla.stack(StackEst, StackPred)
  
  #===============
  ### Build models
  #===============
  pcTimeAutoCor <- list(prior='pccor1', param=prior.pccor)
  
  if(is.null(offset)){
    forEffect <- paste(colnames(XPred),collapse = "+")

    formule <- formula(paste("y ~ 0 + Intercept +",
                             forEffect,
                             "+ f(i, model=SPDE, 
                                  group = IDSpaceTime$i.group, 
                                  control.group = list(model = 'ar1',
                                                       hyper = list(theta = pcTimeAutoCor)))", collapse = "+"))
  }else{
    forEffect <- paste(colnames(XPred)[colnames(XPred) != offset],
                       collapse = "+")
    
    formule <- formula(paste("y ~ 0 + Intercept + ",
                             forEffect," + offset(",offset,") + f(i, model=SPDE, 
                                               group = IDSpaceTime$i.group, 
                                               control.group = list(model = 'ar1',
                                                                    hyper = list(theta = pcTimeAutoCor)))", collapse = "+"))
  }  
  
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