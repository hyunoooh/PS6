#========================================
# Statistical Programming
# Problem Set
# 
# Hyunjoo Oh
#========================================

rm(list=ls())

# We want to find the area under the curve (default = one dimension)
# We want to find the volumn under a multi-dimensional hyperplane 
# (increase dimensionality >> multiple dinemsions)


### Goal 0: Add comments to the integration function
# default = one dimension

sg.int <- function(g,..., lower, upper){
  # load SparseGrid library
  require("SparseGrid")
  
  # validity test:
  # if the length of lower and upper bounds are not equal, show the error message
  if (length(lower) != length(upper)) stop("the length of lower and upper must be the same")
  # if the length of lower and upper bounds are greater than 2, show the error message
  if (length(lower) > 2) stop("Use sg.int.multidim() function, instead")
  if (length(upper) > 2) stop("Use sg.int.multidim() function, instead")
  
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
  val.sp <- gx.sp %*% weights
  # return final numeric values
  as.numeric(val.sp)
}


### Goal 1: Increase Dimensionality

sg.int.multidim <- function(g, ..., lower, upper, dim){
  # load SparseGrid library
  require("SparseGrid")
  
  # validity checking:
    # if the length of lower and upper bounds are not equal, show the error message
    if (length(lower) != length(upper)) stop("the length of lower and upper must be the same")
    # if lengths of lower and upper bounds are not equal to dimesion, show the error message
    if (length(lower) != dim) stop("the length of lower must be equal to dim")
    if (length(upper) != dim) stop("the length of upper must be equal to dim")
    # if input dim is not "integer, greater than 0," show the error message
    if (dim%%1 != 0 | dim <= 0) stop("dim must be integer, greater than 0")
  
  # set lower bound of integration:
  # a numeric vector containing the largest integers of lower values
  lower <- floor(lower)
  # set upper bound of integration:
  # a numeric vector containing the smallest integers of upper values
  upper <- ceiling(upper)
    # if any lower bound is greater than upper bound, show the error message
    if (any(lower > upper)) stop("lower must be smaller than upper")
  
  # create all lower/upper bounds of integration
  # because it's a list, we can use lapply function
  all.bounds <- lapply(1:dim, function(x){seq(lower[x], upper[x]-1, by=1)})
  # create a matrix of all combinations of variables
  gridss <- as.matrix(expand.grid(all.bounds))
  
  # create integration grid with the least number of nodes
  sp.grid <- createIntegrationGrid(type = 'KPU', dimension = dim, k = 5 )
                                # KPU is a nested rule for unweighted integral over [0,1]
                                # dimension of integration is user input
                                # k is accuracy level. 
  
  # create nodes for lower bounds
  nodes <- gridss[1,] + sp.grid$nodes
  # create weights for lower bounds using weights in sp.grid
  # we need to create weight values as many as number of rows of gridss
  weights <- rep(sp.grid$weights, nrow(gridss))
  
  # create nodes for combination of all points in integration
  # because it's a list, we can use lapply function
  nodes <- lapply(1:nrow(gridss), function(x){gridss[x,] + sp.grid$nodes})
  # create node by binding by rows
  nodes <- do.call(rbind, nodes)  

  # apply function g over each set of nodes
  gx.sp <- apply(X = nodes, MARGIN = 1, FUN = g, ...)
  # multiply weights 
  val.sp <- gx.sp %*% weights
  # return final numeric values
  as.numeric(val.sp)
}

### Goal 2: Parallel

sg.int.parallel <- function(g, ..., lower, upper, dim){
  # load SparseGrid library; load parallel library
  require("SparseGrid")
  require("parallel")
  
  # calculate the number of cpu cores on the current host 
  no_cores <- detectCores()
  # initiate cluster
  cl <- makeCluster(no_cores)
  
  # validity checking:
    # if the length of lower and upper bounds are not equal, show the error message 
    if (length(lower) != length(upper)) stop("the length of lower and upper must be the same")
    # if lengths of lower and upper bounds are not equal to dimesion, show the error message
    if (length(lower) != dim) stop("the length of lower must be equal to dim")
    if (length(upper) != dim) stop("the length of upper must be equal to dim")
    # if input dim is not "integer, greater than 0," show the error message
    if (dim%%1 != 0 | dim <= 0) stop("dim must be integer, greater than 0")
  
  # set lower bound of integration:
  # a numeric vector containing the largest integers of lower values
  lower <- floor(lower)
  # set upper bound of integration:
  # a numeric vector containing the smallest integers of upper values
  upper <- ceiling(upper)
    
  # validity checking:
    # if any lower bound is greater than upper bound, show the error message
    if (any(lower > upper)) stop("lower must be smaller than upper")
  
  # create all lower/upper bounds of integration
                # parLappy applies operations of list parallelization using clusters
  all.bounds <- parLapply(cl, 1:dim, 
                          function(x){seq(lower[x], upper[x]-1, by=1)})
  # create a matrix of all combinations of variables
  gridss <- as.matrix(expand.grid(all.bounds))
  
  # create integration grid with the least number of nodes
  sp.grid <- createIntegrationGrid(type = 'KPU', dimension = dim, k = 5 )
                                  # KPU is a nested rule for unweighted integral over [0,1]
                                  # dimension of integration is user input
                                  # k is accuracy level. 
  
  # create weights for lower bounds using weights in sp.grid
  # we need to create weight values as many as number of rows of gridss
  weights <- rep(sp.grid$weights, nrow(gridss))
  
  # create nodes for combination of all points in integration
  # because it's a list, we want to use lapply function;
  # parLappy applies operations of list parallelization using clusters
  nodes <- parLapply(cl = cl, 
                     X = 1:nrow(gridss), 
                     fun = function(x){gridss[x,] + sp.grid$nodes})
  # create nodes by binding by rows
  require("plyr")
  nodes <- do.call("rbind", nodes)
  
  # apply function g over each set of nodes
  # foy applying operations using clusters, we can use parApply
  gx.sp <- parApply(cl = cl, 
                    X = nodes, 
                    MARGIN = 1, 
                    FUN = g, ...)
  # multiply weights 
  val.sp <- gx.sp %*% weights
  # Once we are done close the cluster so that resources 
  # such as memory are returned to the operating system.
  stopCluster(cl)
  # return final value as numeric
  as.numeric(val.sp)
}



### Goal 3: Unit Testing
library(testthat)

# sample for testing
# dimension = 2
test.f1 <- function(x) x[1] + x[2]^2
# dimension = 3
test.f2 <- function(x) x[1] + 2*x[2]^2 + 3*x[3]^3
# dimension = 4
test.f3 <- function(x) 2*x[1] + 3*x[2]^2 + 4*x[3]^3 + x[4]^4

# test basic function: sg.int()
sg.int(test.f1, lower = c(-1, -1), upper = c(1, 1))
# sg.int(test.f2, lower = c(-2, -2, -2), upper = c(2, 2, 2)) # wrong dimension >> error testing
# sg.int(test.f3, lower = c(-1, -1, -1, -1), upper = c(1, 1, 1, 1))  # wrong dimension >> error testing

# test integration function for multidimension: sg.int.multidim()
sg.int.multidim(test.f1, lower = c(-1, -1), upper = c(1, 1), dim = 2)
sg.int.multidim(test.f2, lower = c(-2, -2, -2), upper = c(2, 2, 2), dim = 3)
sg.int.multidim(test.f3, lower = c(-1, -1, -1, -1), upper = c(1, 1, 1, 1), dim = 4)  

# test integration function for parallelization: sg.int.parallel()
sg.int.parallel(test.f1, lower = c(-1, -1), upper = c(1, 1), dim = 2)
sg.int.parallel(test.f2, lower = c(-2, -2, -2), upper = c(2, 2, 2), dim = 3)
sg.int.parallel(test.f3, lower = c(-1, -1, -1, -1), upper = c(1, 1, 1, 1), dim = 4) 

    ### Goal 5: Package cubature
    # Load the package "cubature" to use the function adaptIntegrate()
    library(cubature)
    ?adaptIntegrate
    # we can use this function in unit-testing:

# test.f1 with sg.int(), sg.int.multidim(), sg.int.parallel
expect_equal(sg.int(test.f1, lower = c(-1, -1), upper = c(1, 1)),
             adaptIntegrate(f=test.f1, lowerLimit = c(-1, -1), upperLimit = c(1, 1))$integral,
             tolerance = 0.0001)
expect_equal(sg.int.multidim(test.f1, lower = c(-1, -1), upper = c(1, 1), dim = 2),
             adaptIntegrate(f=test.f1, lowerLimit = c(-1, -1), upperLimit = c(1, 1))$integral,
             tolerance = 0.0001)
expect_equal(sg.int.parallel(test.f1, lower = c(-1, -1), upper = c(1, 1), dim = 2),
             adaptIntegrate(f=test.f1, lowerLimit = c(-1, -1), upperLimit = c(1, 1))$integral,
             tolerance = 0.0001)

# test.f2 and test.f3 with sg.int.multidim(), sg.int.parallel()
expect_equal(sg.int.multidim(test.f2, lower = c(-2, -2, -2), upper = c(2, 2, 2), dim = 3),
             adaptIntegrate(f=test.f2, lowerLimit = c(-2, -2, -2), upperLimit = c(2, 2, 2))$integral,
             tolerance = 0.0001)
expect_equal(sg.int.multidim(test.f3, lower = c(-1, -1, -1, -1), upper = c(1, 1, 1, 1), dim = 4),
             adaptIntegrate(f=test.f3, lowerLimit = c(-1, -1, -1, -1), upperLimit = c(1, 1, 1, 1))$integral,
             tolerance = 0.0001)

expect_equal(sg.int.parallel(test.f2, lower = c(-2, -2, -2), upper = c(2, 2, 2), dim = 3),
             adaptIntegrate(f=test.f2, lowerLimit = c(-2, -2, -2), upperLimit = c(2, 2, 2))$integral,
             tolerance = 0.0001)
expect_equal(sg.int.parallel(test.f3, lower = c(-1, -1, -1, -1), upper = c(1, 1, 1, 1), dim = 4),
             adaptIntegrate(f=test.f3, lowerLimit = c(-1, -1, -1, -1), upperLimit = c(1, 1, 1, 1))$integral,
             tolerance = 0.0001)

# error testing
expect_error(sg.int(test.f2, lower = c(-2, -2, -2), upper = c(2, 2, 2))) # lengths do not match
expect_error(sg.int(test.f3, lower = c(-1, -1, -1, -1), upper = c(1, 1, 1, 1))) # wrong dimension
expect_error(sg.int.parallel(test.f2, lower = c(-2, -2), upper = c(2, 2), dim = 3)) # lengths do not match
expect_error(sg.int.parallel(test.f2, lower = c(-2, -2), upper = c(2, 2), dim = 4.5)) # wrong value for dim
expect_error(sg.int.parallel(test.f3, lower = c(-1, -1, -1, -1), upper = c(1, 1), dim = 4)) # length and dim do not match
expect_error(sg.int.parallel(test.f3, lower = c(-1, -1, -1, -1), upper = c(1, 1, 1, 1), dim = "4")) # dim is not numeric




# Goal 4: Measure Speed
library(microbenchmark)

# using 2D example
basic.vs.multidim <- microbenchmark(
  "sg.int" = sg.int(g = test.f1, lower = c(-1, -1), upper = c(1, 1)),
  "sg.int.multidim" = sg.int.multidim(g = test.f1, lower = c(-1, -1), upper = c(1, 1), dim=2),
  times = 100)
basic.vs.multidim
plot(basic.vs.multidim)

basic.vs.parallel <- microbenchmark(
  "sg.int" = sg.int(g = test.f1, lower = c(-1, -1), upper = c(1, 1)),
  "sg.int.parallel" = sg.int.parallel(g = test.f1, lower = c(-1, -1), upper = c(1, 1), dim=2),
  times = 100)
basic.vs.parallel
plot(basic.vs.parallel)

multidim.vs.parallel <- microbenchmark(
  "sg.int.multidim" = sg.int.multidim(g = test.f1, lower = c(-1, -1), upper = c(1, 1), dim=2),
  "sg.int.parallel" = sg.int.parallel(g = test.f1, lower = c(-1, -1), upper = c(1, 1), dim=2),
  times = 100)
multidim.vs.parallel
plot(multidim.vs.parallel)

# using 4D example 
multidim.vs.parallel.4d <- microbenchmark(
  "sg.int.multidim" = sg.int.multidim(test.f3, lower = c(-1, -1, -1, -1), upper = c(1, 1, 1, 1), dim = 4),
  "sg.int.parallel" = sg.int.parallel(test.f3, lower = c(-1, -1, -1, -1), upper = c(1, 1, 1, 1), dim = 4),
  times = 100)
multidim.vs.parallel.4d
plot(multidim.vs.parallel.4d)


# Goal 5: Package cubature
# cf. see above (Goal 3: Unit testing)
# integrate the same functions and see if the result is faster:

# Speed testing: 
# sg.int.multidim vs. adaptIntegrate: our function takes a lot more time (slower...)
multidim.vs.adaptIntegrate <- microbenchmark(
  "sg.int.multidim" = sg.int.multidim(test.f3, lower = c(-1, -1, -1, -1), upper = c(1, 1, 1, 1), dim = 4),
  "adaptIntegrate" = adaptIntegrate(f=test.f3, lowerLimit = c(-1, -1, -1, -1), upperLimit = c(1, 1, 1, 1)),
  times = 100)
multidim.vs.adaptIntegrate
plot(multidim.vs.adaptIntegrate)

# sg.int.parallel vs. adaptIntegrate: our function takes a lot more time (much slower...)
parallel.vs.adaptIntegrate <- microbenchmark(
  "sg.int.parallel" = sg.int.parallel(test.f3, lower = c(-1, -1, -1, -1), upper = c(1, 1, 1, 1), dim = 4),
  "adaptIntegrate" = adaptIntegrate(f=test.f3, lowerLimit = c(-1, -1, -1, -1), upperLimit = c(1, 1, 1, 1)),
  times = 100)
parallel.vs.adaptIntegrate
plot(parallel.vs.adaptIntegrate)



### Additional Goal: Finding maximum

# input given function
# task 1: with x
task.fun.1 <- function(x) {
  sin((x[1]^2)/2 - ((x[2]^2)/4)*cos(2*x[1] - exp(x[2])))
}

# optimization based on "L-BFGS-B" method
# The code for method "L-BFGS-B" is based on Fortran code 
# by Zhu, Byrd, Lu-Chen and Nocedal obtained from Netlib 
# (file ‘opt/lbfgs_bcm.shar’: another version is in ‘toms/778’).
max.fun.1 <- optim(par = c(2, 1), fn = task.fun,
                 lower = c(-1, -1), upper = c(3, 3),
                 method = "L-BFGS-B",
                 control = list(fnscale = -1))
# check what's in the function
max.fun.1$par
max.fun.1$value
max.fun.1$counts
max.fun.1$convergence
max.fun.1$message

# input given function
# task 2: with x, y
task.fun.2 <- function(x, y) {
  sin((x^2)/2 - ((y^2)/4)*cos(2*x - exp(y)))
}

# optimize over x for y within the interval [-1, 3]
task.slice.y <- lapply(seq(-1, 3, by=0.0001),
                      function(z){
                        optimize(f = task.fun.2, y = z,
                                 lower = -1, upper = 3,
                                 maximum = TRUE)$maximum
                      })
# find the maximum
task.maximum <- unlist(task.slice.y)
head(task.maximum, 100)
which.max(task.maximum)
task.maximum[which.max(task.maximum)]
