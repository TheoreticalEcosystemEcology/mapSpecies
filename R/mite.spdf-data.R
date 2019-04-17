#' Mite data
#'
#' Oribatid mite data with explanatory variables
#'
#' @docType data
#'
#' @usage data(mite.spdf)
#'
#' @details These data were gathered in a transect of 2.6 m by 10 m located on the border of a small Laurentide lake, Lake Geai, located in the Station de Biologie des Laurentides. The short side of the transect is bordering the lake. 
#' 
#' The \code{mite.spdf} data is a \code{SpatialPointDataFrame} contaning 35 pseudospecies sampled at 70 locations in the 2.6 m by 10 m transect. 
#'
#' The \code{mite.envRaster} data is a \code{rasterStack} object contaning 5 environmental variables gather across the 2.6 m by 10 m transect. Originally, substrate density (SubsDens) and soil water content (WatrCont) were only gathered at the sampled locations; they were interpolated using kriging for values to be available across the full area of the transect. The other thress variables (Substrate, Shrub and Topo) were rasterized using Figure 1 of Borcard and Legendre (1994) as reference.
#'
#' The projection used for these data preserves the fine (meter) scale, which is useful for some of the illustration presented in the vignettes.
#' 
#'
#' @references
#' 
#' Borcard, D., and P. Legendre. 1994. Environmental control and spatial structure in ecological communities: an example using oribatid mites (Acari, Oribatei). \emph{Environmental and Ecological Statistics} \strong{1}:37–61.
#'
#' @keywords datasets
"mite.spdf"
