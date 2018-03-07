context("test-mapspecies.R")

test_that("PP model works", {
  require(spatstat)
  
  #================
  ### Simulate data
  #================
  ### Define window size
  win <- owin(c(0,3), c(0,3))
  
  ### Define the number of pixel in the window (300 X 300)
  spatstat.options(npixel=300)
  
  ### Define the mean of the LGRF
  beta0 <- 3
  
  ### Expected (average) number of point in the LGRF
  areaPoly <- diff(range(win$x)) * diff(range(win$y))
  exp(beta0) * areaPoly
  
  ### Variance
  sigma2x <- 0.2
  
  ### 1/Scale
  kappa <- 2
  
  ### Build covariate
  y0 <- x0 <- seq(win$xrange[1], win$xrange[2],
                  length=spatstat.options()$npixel)
  gridcov <- outer(x0, y0, function(x,y) cos(x) - sin(y-2))
  
  gridcovRaster <- raster(gridcov[nrow(gridcov):1,], xmn = 0, 
                          xmx = 3, ymn = 0, ymx = 3)
  
  ### Define parameter associated to this covariate
  beta1 <- -0.5
  
  ### Expected (average) number of point in the LGRF
  sum(exp(beta0 + beta1*gridcov) * diff(x0[1:2])*diff(y0[1:2]))
  
  ### Simulate the data
  spatstat.options(npixel=300)
  set.seed(1)
  lg.s.c <- rLGCP('matern', im(beta0 + beta1*gridcov,
                               xcol=x0, yrow=y0), 
                  var = sigma2x, scale = 1/kappa, nu = 1, win = win)
  
  ### Point pattern
  xy.c <- cbind(lg.s.c$x, lg.s.c$y)[,2:1]
  n.c <- nrow(xy.c)
  
  ### Draw point patterns on the full simulated LGRF
  par(mfrow=c(1,2), mar=c(2,2,1,1), mgp=c(1,0.5,0))
  image.plot(list(x=x0, y=y0, z=gridcov), main='Covariate', asp=1)
  image.plot(list(x=x0, y=y0, z=log(attr(lg.s.c, 'Lambda')$v)),
             main='log-Lambda', asp=1) 
  points(xy.c, pch=19)
  
  #=============
  ### Build mesh
  #=============
  loc.d <- 3*t(matrix(c(0,0,1,0,1,1,0,1,0,0), 2))
  domainSP <- SpatialPolygons(list(Polygons(list(Polygon(loc.d)),
                                            '0')))
  
  mesh <- spatialMesh(coords = loc.d,generic=FALSE, offset=c(0.3, 1),
                      max.edge=c(0.3, 0.7), cutoff=0.05)
  

  ### Build model
  ppModel <- spatialPP(y ~ layer, y = xy.c, X = gridcovRaster, 
                       sp = domainSP, mesh = mesh, smooth = 2,
                       prior.range = c(0.05, 0.1),
                       prior.sigma = c(1, 0.01))
  
  ppPred <- mapSpatialPP(ppModel, mesh = mesh, 
                       resolution = c(500,500), type = "mean",sp = domainSP)
  
  image.plot(list(x=x0, y=y0, z=attr(lg.s.c, 'Lambda')$v),
             main='Lambda', asp=1)
  contour(ppPred,add=TRUE,lwd=2,asp=1)

})


