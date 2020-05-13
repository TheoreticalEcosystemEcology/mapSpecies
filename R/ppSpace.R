#' @title Spatial point process model
#' @name ppSpace
#'
#' @description Spatial point process model using INLA. This function is essentially a specialized wrapper over \code{inla}
#'
#' @param formula A formula that only relates the response \code{y} and some (or all) of the explanatory variables \code{X}. A paricularity of the is formula is that the response has to be defined as \code{y}.
#' @param sPoints A \code{SpatialPoint*} object that includes the sample location of the modelled species.
#' @param ppWeight An object of class \code{ppWeight}
#' @param explanaMesh An object of class \code{explanaMesh}
#' @param smooth A single value ranging between 0 and 2 passed to \code{inla.spde2.pcmatern}. It defines the smoothness of the Matern SPDE model. Default is set at 2.
#' @param prior.range A vector of length 2, with (range0, Prange) specifying that P(ρ < ρ_0) = p_ρ, where ρ is the spatial range of the random field. If Prange is NA, then range0 is used as a fixed range value. Default is c(0.05, 0.01).
#' @param prior.sigma A vector of length 2, with (sigma0, Psigma) specifying that P(σ > σ_0) = p_σ, where σ is the marginal standard deviation of the field. If Psigma is NA, then sigma0 is used as a fixed range value.  Default is c(1, 0.01).
#' @param many Logical. Whether the data in \code{sPoints} is large or not. See details. Default is \code{FALSE}.
#' @param \dots Arguments passed to \code{inla}
#'
#' @details 
#' 
#' If the argument \code{many = TRUE}, the estimation and the prediction will be carried out solely at the mesh edges, whereas when \code{many = FALSE} the estimation will be carried out at the mesh edges and at the sampled location. When the number of samples is very large (e.g. tens of thousands of samples or more) using \code{many = TRUE} can be much more computationally efficient. However, there is a precision trade-off. When \code{many = TRUE}, each sample is associated to an edge and the model is constructed using the number of samples associated to an edge as an importance value. In doing so, some spatial precision is lost at the expense of speed. 
#'
#' @return
#' 
#' An object of class \code{ppSpace} that includes a model output, which is the model output of INLA.
#' 
#' In addition, it includes a series of attributes:
#' 
#'	  \item{\code{formula}}{The formula used to construct the model}
#'	  \item{\code{sPointsDF}}{A \code{SpatialPointDataFrame} object that includes the sample location and associated data of the modelled species.}
#'	  \item{\code{XEst}}{A matrix with all the explanatory variables used to construct the model. If there were factors in the original set of explanatory variables \code{X}, in \code{XEst}, they were decomposed into dummy variables. The values in \code{XEst} are the one from the sampled location.}
#'	  \item{\code{XPred}}{A matrix with all the explanatory variables used to construct the model. If there were factors in the original set of explanatory variables \code{X}, in \code{XPred}, they were decomposed into dummy variables. The values in \code{XPred} were gathered at the mesh edges. When \code{many = TRUE}, the values in \code{XPred} are exactly the same as the values in \code{XEst}}
#'	  \item{\code{mesh}}{An object of class \code{inla.mesh}. It is the mesh used to construct the model.}
#'	  \item{\code{Stack}}{An object of class \code{inla.data.stack}. It is a stack object constructed internally.}
#'
#' @importFrom sp coordinates
#' @importFrom raster rasterFromXYZ
#' @importFrom sp SpatialPoints
#' @importFrom raster extract
#' @importFrom INLA inla
#' @importFrom INLA inla.spde2.pcmatern
#' @importFrom INLA inla.spde.make.A
#' @importFrom INLA inla.spde.make.index
#' @importFrom INLA inla.stack
#' @importFrom INLA inla.stack.data
#' @importFrom INLA inla.stack.A
#' @importFrom Matrix Diagonal
#' @importFrom stats model.matrix
#' @importFrom stats model.frame
#' 
#' @export
#' 
#' @keywords models
ppSpace <- function(formula,
                    sPoints, 
                    ppWeight, 
                    explanaMesh, 
                    smooth = 2,
                    prior.range = c(0.05, 0.01),
                    prior.sigma = c(1, 0.01), 
                    many = FALSE,
                    ...){

  #==============
  ### Basic check
  #==============
  if(as.character(formula[[2]]) != "y"){
    stop("'y' should be used to define the response variable")
  }
  
  ### Check if the mesh in ppWeight and dataPred are the same
  if(!identical(attributes(ppWeight)$graph, explanaMesh$mesh$graph)){
    stop("'ppWeight' and 'explanaMesh' were constructed using different mesh")
  }
  
  #================
  ### Basic objects
  #================
  nsPoints <- length(sPoints)
  nEdges <- explanaMesh$mesh$n
  xy <- coordinates(sPoints)

  #==============
  ### Define SPDE
  #==============
  SPDE <- inla.spde2.pcmatern(mesh=explanaMesh$mesh, alpha=smooth,
                              prior.range=prior.range,
                              prior.sigma=prior.range)

  #==========================
  ### Define response objects
  #==========================
  if(many){
    ### Aggregate spatial data
    xyDF <- as.data.frame(xy)
    spaceAgg <- aggData(xyDF, explanaMesh$mesh)
    
    ### Pseudo-absences are the number of edges on the mesh
    ### Occurences are the number of points
    yPP <- spaceAgg$Freq
    
    ### weight associated to pseudo-absences (w) and occurrences (0)
    ePP <- ppWeight[spaceAgg$space]
  }else{
    ### Pseudo-absences are the number of edges on the mesh
    ### Occurences are the number of points
    yPP <- rep(0:1, c(nEdges, nsPoints))
    
    ### weight associated to pseudo-absences (w) and occurrences (0)
    ePP <- c(ppWeight, rep(0, nsPoints))
  }
  
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
  
  if(many){
    ### Extract covariate values for model estimation
    locEst <- SpatialPoints(coords = explanaMesh$mesh$loc[,1:2])
    XEst <- extract(Xbrick, locEst)
  }else{
    ### Extract covariate values for model estimation
    meshxy <- rbind(explanaMesh$mesh$loc[,1:2],xy)
    rownames(meshxy) <- 1:nrow(meshxy)
    
    locEst <- SpatialPoints(coords = meshxy)
    XEst <- extract(Xbrick, locEst)
    
    ### Extract covariate values for model prediction
    locPred <- SpatialPoints(coords = explanaMesh$mesh$loc[,1:2])
    XPred <- explanaMesh$Xmesh
  }
  
  #==================================================
  ### Define projection matrix and build stack object
  #==================================================
  if(many){
    ### Projection matrix
    ProjInfer <- inla.spde.make.A(explanaMesh$mesh,
                                  loc = explanaMesh$mesh$loc[spaceAgg$space,])
    
    IDSpace <- inla.spde.make.index('i', nEdges)
  
    ### Build stack objects
    StackEst <- inla.stack(data = list(y = yPP, e = ePP), 
                           A = list(1, ProjInfer), 
                           effects = list(list(Intercept = 1, 
                                               X = XEst),
                                          IDSpace), 
                           tag = "est")
    
    StackPred <- inla.stack(data = list(y = NA, e = NA),
                            A = list(1, ProjInfer), 
                            effects = list(list(Intercept = 1, 
                                                X = XEst), 
                                           IDSpace),
                            tag = "pred")
  }else{
    #--------------------
    ### Projection matrix
    #--------------------
    ### For inference
    ProjInfer <- inla.spde.make.A(explanaMesh$mesh, xy)
    
    ### For integration
    ProjInter <- Diagonal(nEdges, rep(1, nEdges))
    
    ### Combine both projection matrix
    A <- rbind(ProjInter, ProjInfer)

    #----------------------
    ### Build stack objects
    #----------------------
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
  }
  
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
  
  nameRes <- names(model)
  
  #===============
  ### Return model
  #===============
  ### add a series of attributes to the result object
  if(many){
    attributes(model) <- list(formula = formula,
                              sPoints = sPoints,
                              XEst = XEst,
                              XPred = XEst,
                              mesh = explanaMesh$mesh,
                              Stack = Stack)
  }else{
    attributes(model) <- list(formula = formula,
                              sPoints = sPoints,
                              XEst = XEst,
                              XPred = XPred,
                              mesh = explanaMesh$mesh,
                              Stack = Stack)
  }
  
  names(model) <- nameRes
  
  class(model) <- "ppSpace"
  
  return(model)
}
