#' @title Spatial point process model for large data
#' @name ppSpaceLarge
#'
#' @description Spatial point process model using INLA designed for data that are composed of many occurrences. This function is essentially a sophisticated wrapper over \code{inla}
#'
#' @param formula A formula that only relates the response \code{y} and some (or all) of the explanatory variables \code{X}. A paricularity of the is formula is that the response has to be defined as \code{y}.
#' @param y A 2-columns matrix defining the location of where a species was found.
#' @param weightPP An object of class \code{weightPP}
#' @param explanaMesh An object of class \code{dataPrep}
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
ppSpaceLarge <- function(formula, y, weightPP, explanaMesh, 
                      smooth = 2,
                      prior.range = c(0.05, 0.01),
                      prior.sigma = c(1, 0.01), ...){
  
  #==============
  ### Basic check
  #==============
  if(as.character(formula[[2]]) != "y"){
    stop("'y' should be used to define the response variable")
  }
  
  ### Check if the mesh in weightPP and explanaMesh are the same
  if(!identical(weightPP$mesh$graph, explanaMesh$mesh$graph)){
    stop("'weightPP' and 'explanaMesh' were constructed using different mesh")
  }
  
  #================
  ### Basic objects
  #================
  ny <- nrow(y)
  nEdges <- weightPP$mesh$n
  
  #=========================
  ### Aggregate spatial data
  #=========================
  colnames(y) <- c("x", "y")
  y <- as.data.frame(y)
  spaceAgg <- aggData(y, weightPP$mesh)
  
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
  yPP <- spaceAgg$Freq
  
  ### weight associated to pseudo-absences (w) and occurrences (0)
  ePP <- weightPP$weight[spaceAgg$space]
  
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
  locEst <- SpatialPoints(coords = weightPP$mesh$loc[,1:2])
  XEst <- extract(Xbrick, locEst)

  #===========================
  ### Define projection matrix
  #===========================
  ### For inference
  ProjInfer <- inla.spde.make.A(weightPP$mesh, 
                                loc = weightPP$mesh$loc[spaceAgg$space,])
  
  
  IDSpace <- inla.spde.make.index('i', nEdges)
  
  #=====================
  ### Build stack object
  #=====================
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
  res <- list(model = model, Stack = Stack, mesh = dataPrep$mesh)
  
  class(res) <- "ppSpace"
  return(res)
}