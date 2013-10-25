source('/home/darkil73/Documents/Dropbox/Cours-ENSAE/Simulation et Monte Carlo/Projet/gamma.R')

##########################################################################
# PROCESS SIMULATION                                                     #
##########################################################################

varianceGammaVector <- function(t,sigma,theta,kappa){
  process <- c(0)
  sum <- 0
  for (i in 1:t){
    increment <- gammaVariable(1,1/kappa)/kappa
    norm <- rnorm(1)
    increment <- sigma*norm*sqrt(increment) + theta*increment
    sum <- sum + increment
    process <- c(process,sum)    
  }
  return (process)
}

# plot(1:10001,varianceGammaVector(10000,0.1,0.1,2), type = "s", xlab = "temps t", ylab = "X(t)")

##########################################################################
# SPEED TESTS                                                            #
##########################################################################

speedVarianceGamma <- function(param1,param2,t,fixed="kappa"){
  if (fixed=="sigma"){
    sigma = 1
    theta = param1
    kappa = param2
  }
  else{
    if (fixed=="theta"){
      theta = 1
      sigma = param1
      kappa = param2
    }
    else{
      kappa = 1
      sigma = param1
      theta = param2
    }
  }
  return (system.time(varianceGammaVector(t,sigma,theta,kappa))[3])
}

speedVarianceGammaMatrix <- function(sigma=c(1),theta=c(1),kappa=c(1),t=100,fixed=c("kappa")){

  if (fixed=="sigma"){
    param1 = theta
    param2 = kappa
  }
  else{
    if (fixed=="theta"){
      param1 = sigma
      param2 = kappa
    }
    else{
      param1 = sigma
      param2 = theta
    }
  }
  
  nparam2 = length(param2)
  nparam1 = length(param1)
  
  param1Rows <- replicate(nparam2,param1)
  param2Cols <- replicate(nparam1,param2,simplify=FALSE)
  
  param1Matrix <- matrix(param1Rows,nparam1)
  param2Matrix <- do.call(rbind,param2Cols)
  
  nt <- length(t)
  result <- vector("list",nt)
  
  for (i in 1:nt){ 
    speedVarianceGammaFixedt <- function(parameter1,parameter2){
      return(speedVarianceGamma(parameter1,parameter2,t[i],fixed))
    }
    speedM <- matrix(mapply(speedVarianceGammaFixedt, param1Matrix, param2Matrix),ncol=nparam2)
    result[[i]] <- speedM
  }

  return(result)
}

# speedVarianceGammaMatrix(sigma=seq(0.1,10,by=0.1),theta=seq(0.1,10,by=0.1),t=100,fixed="kappa")
speedVK <- speedVarianceGammaMatrix(sigma=seq(0.1,5.1,by=0.5),theta=seq(0.5,5.1,by=0.5),t=100,fixed="kappa")
speedVS <- speedVarianceGammaMatrix(sigma=seq(0.1,10,by=0.1),theta=seq(0.1,10,by=0.1),t=100,fixed="sigma")
speedVT <- speedVarianceGammaMatrix(sigma=c(1),kappa=seq(0.1,10.1,by=0.5),t=100,fixed="theta")
##########################################################################
# EFFICIENCY TESTS                                                       #
##########################################################################

library(cubature)

cdfVarianceGammaIncrements <- function(z,param1,param2){
  
  kappa <- param1
  theta <- param2
  
  integrand <- function(M){
    x <- M[1]
    y <- M[2]
    z <- M[3]
    return (dnorm(x)*dgamma(y^2/x^2,shape=1/kappa,rate=kappa)*dgamma((z-y)/theta,shape=1/kappa,rate=kappa))
  }
  
  cdf <- adaptIntegrate(integrand,
                             c(-100,-100,-100),
                             c(100,100,z),
                             maxEval=10000,
                             absError = 1e-05)$integral  
  return(cdf)
}

varianceGammaEfficiencyKS <- function(sigma,theta,kappa,t,N=100){
  
  simulatedProcess <- varianceGammaVector(t,sigma,theta,kappa)
  simulatedIncrements <- diff(simulatedProcess,1)
  
  xMin = floor(min(simulatedIncrements))
  xMax = floor(max(simulatedIncrements))+1
  sub <- sapply(0:N, FUN=function(i) xMin+i*(xMax-xMin)/N)
  
  cdfEmpirical <- sapply(sub,
                         FUN=function(x) 
                           length(simulatedIncrements[simulatedIncrements<x])/(t-1)
                         )
  
  difference <- abs(sapply(sub,
                           FUN=function(x) 
                             abs(cdfVarianceGammaIncrements(x,kappa,theta))
                           )
                    -cdfEmpirical
                    )
  
  statistic <- (sqrt(t-1))*max(difference)
  
  if (statistic > 1.36){
    result = "rejet"
  }
  
  else{
    result = "non-rejet"
  }
  
  return(result)
}

#varianceGammaEfficiency(1,1,1,100)