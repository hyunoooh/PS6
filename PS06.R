#========================================
# Statistical Programming
# Problem Set
# 
# Hyunjoo Oh
#========================================

rm(list=ls())

# Goal 1: Increase Dimensionality

sg.int <- function(g,..., lower, upper){
  # load SparseGrid library
  require("SparseGrid")
  # set lower bound of integration:
  # a numeric vector containing the largest integers of lower values
  lower <- floor(lower)
  # set upper bound of integration:
  # a numeric vector containing the smallest integers of upper values
  upper <- ceiling(upper)
  
  # if any lower bound is greater than upper bound, show the error message
  if (any(lower > upper)) stop("lower must be smaller than upper")
  # create a matrix of all combinations of variables
  gridss <- as.matrix(expand.grid(seq(lower[1], upper[1]-1, by=1),
                                seq(lower[2], upper[2]-1, by=1)))
  # create integration grid with the least number of nodes
  sp.grid <- createIntegrationGrid(type = 'KPU', dimension = 2, k = 5 )
                                  # KPU is a nested rule for unweighted integral over [0,1]
                                  # k is accuracy level. 
  # create nodes for lower bounds
  nodes <- gridss[1,]+sp.grid$nodes
  # create weights for lower bounds using weights in sp.grid
  weights <- sp.grid$weights
  
  # create nodes for combination of all points in integration
  for (i in 2:nrow(gridss)){
    # create node for i by binding by rows
    nodes <- rbind(nodes, gridss[i,]+sp.grid$nodes)  
    # create weights for combination i
    weights <- c(weights, sp.grid$weights)
  }
  # apply function over each set of nodes
  gx.sp <- apply(X = nodes, MARGIN = 1, FUN = g, ...)
  # multiply weights 
  val.sp <- gx.sp %*%weights
  # return final values
  val.sp
}

# Goal 2: Parallel

sg.int.parallel <- function(g,..., lower, upper, dim){
  # load SparseGrid library; load parallel library
  require("SparseGrid")
  require("parallel")
  # calculate the number of cpu cores
  no_cores <- detectCores()
  # Initiate cluster
  cl <- makeCluster(spec = no_cores)

  # set lower bound of integration:
  # a numeric vector containing the largest integers of lower values
  lower <- floor(lower)
  # set upper bound of integration:
  # a numeric vector containing the smallest integers of upper values
  upper <- ceiling(upper)
    # if any lower bound is greater than upper bound, show the error message
  if (any(lower > upper)) stop("lower must be smaller than upper")

  # create all lower/upper bounds of integration
  all.bounds <- parLapply(cl, 1:dim, function(x){
                           seq(lower[x], upper[x]-1, by=1)
                           })
  # create a matrix of all combinations of variables
  gridss <- as.matrix(expand.grid(all.bounds))
  
  # create integration grid with the least number of nodes
  sp.grid <- createIntegrationGrid(type = 'KPU', dimension = dim, k = 5 )
                                  # KPU is a nested rule for unweighted integral over [0,1]
                                  # dimension of integration is user input
                                  # k is accuracy level. 
  # create nodes for lower bounds
  nodes <- parLapply(cl, 1:nrow(gridss), function(x){
                       gridss[x,]+sp.grid$nodes
                     })
  # create weights for lower bounds using weights in sp.grid
  # we need to create weight values as many as number of rows of gridss
  weights <- rep(sp.grid$weights, nrow(gridss))
  
  ######## 
  # create nodes for combination of all points in integration
  for (i in 2:nrow(gridss)){
    # create node for i by binding by rows
    nodes <- rbind(nodes, gridss[i,]+sp.grid$nodes)  
    # create weights for combination i
    weights <- c(weights, sp.grid$weights)
  }
  # apply function over each set of nodes
  gx.sp <- apply(X = nodes, MARGIN = 1, FUN = g, ...)
  # multiply weights 
  val.sp <- gx.sp %*%weights
  # return final values
  val.sp
}

# Goal 3: Unit Testing
# Goal 4: Measure Speed
# Goal 5: Package cubature

# 
library(testthat)

