#' @title Build triangulation for a spatial polygon
#'
#' @description Builds a triangulation for a spatial polygon where a set of coordinates (usually fictive) is included to define the size of the triangles. In short, where there are less points the triangles will be larger insurance a more reliable model estimation in areas where there are less sampling efforts.
#'
#' @param sp A \code{SpatialPolygons}, \code{SpatialPolygonsDataFrame}, \code{SpatialLines}, \code{SpatialLinesDataFrame}, \code{SpatialPoints} or \code{SpatialPointsDataFrame} object
#'
#' @param coords A 2-columns matrix of point locations
#' @param generic Logical. Defines a generic mesh (see details).
#' @param \dots Arguments passed to \code{\link{inla.mesh.2d}}
#'
#' @details
#'
#' Although the \code{generic} argument can build a mesh designed for any spatial domain, the parameters used to build the mesh were chosen for North America.
#'
#' NOT YET FULLY TESTED
#' NOT YET FULLY TESTED
#' NOT YET FULLY TESTED
#' NOT YET FULLY TESTED
#' NOT YET FULLY TESTED
#'
#'
#' @return
#'
#' An \code{inla.mesh} object
#'
#' @seealso
#'
#' \code{\link{inla.mesh.2d}}
#'
#' @author F. Guillaume Blanchet
#'
#' @importFrom INLA inla.mesh.2d
#' @export
#'
#' @keywords graphs
#'
spatialMesh <- function(sp, coords=NULL, generic = TRUE,...){
  ### General checks
  if(ncol(coords) != 2){
    stop("'coords' needs to have two columns")
  }

  ### Generic construction of the mesh
  if(generic){
    ### Excluding coordinates
    if(is.null(coords)){
      mesh <- inla.mesh.2d(boundary=inla.sp2segment(sp),
                           max.edge=c(1,3), cutoff=c(0.5, 2))

    ### Including coordinates
    }else{
      mesh <- inla.mesh.2d(boundary=inla.sp2segment(sp),
                           loc.domain = coords,
                           max.edge=1, cutoff=0.3)
    }
  ### Detailled construction of the mesh
  }else{
    ### Excluding coordinates
    if(is.null(coords)){
      mesh <- inla.mesh.2d(boundary=inla.sp2segment(sp), ...)

    ### Including coordinates
    }else{
      mesh <- inla.mesh.2d(boundary=inla.sp2segment(sp),
                           loc.domain = coords, ...)
    }
  }

  ### Return mesh
  return(mesh)
}
