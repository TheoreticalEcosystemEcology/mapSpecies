#' @title Build a raster of expert knowledge prior information 
#' 
#' @description Rasterize an expert maps and weight its importance using the samples and their distance to the expert map with a five parameters logistic curve.
#' 
#' @param expert A \code{\link[sp]{SpatialPolygons}} or a \code{\link[sp]{SpatialPolygonsDataFrame}} defining the expert map.
#' @param sPointDF A \code{\link[sp]{SpatialPointsDataFrame}} defining the sampled location with a single variable.
#' @param raster An object of class \code{RasterLayer}.
#' @param family This argument defines which reference distribution should be used to estimate the parameter of the logistic curve.
#' @param link The link function to use if \code{family = "binomial"}. This argument is ignored if \code{family} is not \code{"binomial"}.
#' @param iniParam A named vector defining the initial parameters of the logistc curve used to weight the importance of the expert map. The name of the parameters need to be \code{upper}, \code{lower}, \code{rate}, \code{shift} and \code{skew}. As default, the values of 1 are given to each parameters.    
#' @param \dots arguments pass to \code{link[gnlm]{bnlr}} if \code{family = "binomial"} or \code{link[gnlm]{gnlr}} otherwise.
#'
#' @details 
#' 
#' This function uses the five parameters logistic curve proposed by Richards (1959) and suggested by Merow et al. (2017) as a way to weight the expert knowledge. The equation associated to this logistic curve is 
#' 
#' \deqn{W(x) = u - \frac{u - l}{\left(1 + e^{-r(x-k)}\right)^{1/s}}}{W(x) = u - (u - l)/((1 + exp(-r(x-k)))^(1/s))}.
#' 
#' where eqn{u} and eqn{l} are the upper and lower asymptotes of the logistic curve, eqn{r} is a rate that gives flexibility to the curve from a sharpe step to a flat surface and eqn{s} is a measure of skewness that adjust the symmetry of the decay on the edge of the expert map. As for eqn{k}, it shifts the curve inside or outside the expert map. Finally, eqn{x} is the distance to the expert map.
#' 
#' The detailed mathematics of the logistics curve is presented in the \code{uniSpace} vignette.
#'
#' The five parameters logistic equation are estimated using nonlinear modelling with the help of the \code{gnlm} R package.
#' 
#' This function is designed to be with presence-absence, abundance and continuous data, but not presence-only. To estimate the parameter of the logistic curve using presence-only data use the \code{bossMaps} R package.
#' 
#' The function tends to send warnings message that stems from the calls to either \code{\link[gnlm]{bnlr}} or \code{\link[gnlm]{gnlr}}. They essentially inform users of different choices automatically made in the code.
#' 
#' @return 
#' 
#' An object of class \code{RasterLayer} that includes expert maps in its object.
#'
#' @references 
#'
#' Merow, C., A. M. Wilson, and W. Jetz. 2017. Integrating occurrence data and expert maps for improved species range predictions, \emph{Global Ecology and Biogeography} \strong{26}:243–258.
#'
#' Richards, F. J. 1959. A flexible growth function for empirical use. \emph{Journal of Experimental Botany} \strong{10}:290–301
#'
#' @importFrom raster rasterize
#' @importFrom raster distance
#' @importFrom gnlm bnlr
#' @importFrom gnlm gnlr
#' @importFrom bossMaps logistic
#' @importFrom raster raster
#' 
#' @export
#' 
#' @keywords models
offsetExpert <- function(expert, sPointDF, raster,
                         family, link, iniParam = c(upper = 1, 
                                                    lower = 1, 
                                                    rate = 1, 
                                                    shift = 1, 
                                                    skew = 1), ...){
 
  #-----------------------------------------------------------
  # Find the distance of each pixel to the expert map boundary
  #-----------------------------------------------------------
  # Rasterize the expert polygon
  expertRast <- rasterize(expert, raster)

  # Calculate the distance
  expertDist <- distance(expertRast)
  
  #--------------------------------
  # Parameterize the logistic curve
  #--------------------------------
  resp <- unlist(sPointDF@data)
  x <- extract(expertDist, sPointDF)
  
  # Estimate the parameter of the logistic curve
  if(family == "binomial"){
    # If data is 0-1
    y <- cbind(resp, 1 - resp)
    
    # Force parameter variable and parameters in the function environment
    holdEnv <- list(y=y, x=x, iniParam=iniParam)
    attach(holdEnv)
    on.exit(detach("holdEnv"))
    
    logisticParam <- gnlm::bnlr(y = y, link = link,
                          mu = ~ (upper - ((upper - lower)/(1 + exp(-rate * (x - shift)))^(1/skew))), 
                          pmu = iniParam, ...)

  }else{
    # Otherwise
    y <- unlist(sPointDF@data)
    
    # Force parameter variable and parameters in the function environment
    holdEnv <- list(y=y, x=x, iniParam=iniParam)
    attach(holdEnv)
    on.exit(detach("holdEnv"))
    
    logisticParam <- gnlr(y = y, distribution = family,
                          mu = ~ (upper - ((upper - lower)/(1 + exp(-rate * (x - shift)))^(1/skew))), 
                          pmu = iniParam, ...)
  }

  # Construct a raster map using the logistic curve
  expertLogistic <- logistic(x = expertDist,
                             upper = logisticParam$coefficients[1],
                             lower = logisticParam$coefficients[2],
                             rate = logisticParam$coefficients[3],
                             shift = logisticParam$coefficients[4],
                             skew = logisticParam$coefficients[5])
  
  coeff <- logisticParam$coefficients
  
  names(coeff) <- c("upper", "lower", "rate", "shift", "skew")
  
  res <- list(expert = expertLogistic,
              coefficients = coeff)
  
  class(res) <- "offsetExpert"

  return(res)
}
