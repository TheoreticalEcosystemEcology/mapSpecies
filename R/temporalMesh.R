#' @title Build a 1 dimensional mesh for a temporal model
#'
#' @description Builds a 1 dimensional mesh for a temporal model. This is essential a generic wrapper for inla.mesh.1d
#'
#' @param knots A vector defining the temporal location of the knots. Default is NULL
#' @param \dots Arguments passed to \code{\link{inla.mesh.1d}}
#'
#' @details 
#' 
#' If \code{knots} is \code{NULL}, the mesh is constructed for one year with a daily resolution where the knots are located at the beginning of each month and a last knot is at the end of the year.
#'
#' @return
#'
#' An \code{inla.mesh} object
#'
#' @seealso
#'
#' \code{\link{inla.mesh.1d}}
#'
#' @author F. Guillaume Blanchet
#'
#' @importFrom INLA inla.mesh.1d
#' @export
#'
#' @keywords graphs
#'
temporalMesh <- function(knots = NULL, ...){
  ### Define knots if argument is NULL
  if(is.null(knots)){
    knots <- c(as.numeric(seq(as.Date("2000/1/1"),
                              as.Date("2000/12/31"), "month")) - 
                 10956, as.numeric(as.Date("2000/12/31"))-10956)
  }
  
  ### Construct mesh
  mesh <- inla.mesh.1d(knots)
  
  ### Return
  return(mesh)
}
