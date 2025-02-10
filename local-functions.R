library(spatstat.linnet)

#'  Compute the inverse of distance function of a point pattern on a linear network.
#'
#'  @description
#'  Adaptation of the method distfunlpp.R in spatstat.linnet
#'  See: https://github.com/spatstat/spatstat.linnet/blob/master/R/distfunlpp.R
#'  Date: 2025-02-10
#' 
#'  @param X A point pattern on a linear network (object of class "lpp").
#'  @param k An integer. The distance to the kth nearest point will be computed.
#'  @return A `function` with arguments `x,y` and optional arguments `seg, tp`. 
#'  It also belongs to the class `"linfun"` which has methods for plot, print etc.
#' 
#'  @seealso [spatstat.linnet::distfun.lpp()]
#' 

distfun.inverse.lpp <- function(X, ..., k=1) {
  stopifnot(inherits(X, "lpp"))
  force(X)
  force(k)
  stopifnot(length(k) == 1)
  L <- as.linnet(X)
  f <- function(x, y=NULL, seg=NULL, tp=NULL, ...) {
    # L is part of the environment
    Y <- as.lpp(x=x, y=y, seg=seg, tp=tp, L=L)
    d <- nncross.lpp(Y, X, what="dist", k=k)
    d <- 1/d # inverse of distance
    return(d)
  }
  f <- linfun(f, L)
  assign("k", k, envir=environment(f))
  assign("X", X, envir=environment(f))
  attr(f, "extrargs") <- list(k=k)
  class(f) <- c("distfunlpp", class(f))    
  return(f)
}