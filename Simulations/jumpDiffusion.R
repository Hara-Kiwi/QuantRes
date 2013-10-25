source('/home/darkil73/Documents/Dropbox/Cours-ENSAE/Simulation et Monte Carlo/Projet/poisson.R')

##########################################################################
# PROCESS SIMULATION                                                     #
##########################################################################

jumpDiffusionVector <- function(lambda,Tlim,tlim){
  N <- poissonProcessPointwise(lambda,Tlim) # Tlim = 2-10 * tlim
  U <- runif(n=N,min=0,max=Tlim)
  Y <- rnorm(n=N,mean=0,sd=1)
  sortedU <- sort(x=U,index.return = TRUE)
  process <- c()
  sum <- 0
  i <- 1
  for (k in 1:tlim){
    while (sortedU$x[i] < k){
      sum <- sum + Y[sortedU$ix[i]]
      i <- i +1
    }
    G <- rnorm(n=1,mean=0,sd=1)
    sum <- sum + G[1]
    process <- c(process,sum)
  }
  return (process)
}

# lambda = 1
# tlim = 1000
# par(mfrow=c(2,1))
# plot(1:tlim,jumpDiffusionVector(lambda,2*tlim,tlim),type='s',xlab="temps t", ylab="X(t)")

##########################################################################
# SPEED TESTS                                                            #
##########################################################################

speedJumpDiffusion <- function(lambda,t){
  return (system.time(jumpDiffusionVector(lambda,2*t,t))[3])
}

speedJumpDiffusionMatrix <- function(lambda,t){
  
  nt = length(t)
  nlambda = length(lambda)
  
  lambdaRows <- replicate(nt,lambda)
  tCols <- replicate(nlambda,t,simplify=FALSE)
  
  lambdaMatrix <- matrix(lambdaRows,nlambda)
  tMatrix <- do.call(rbind,tCols)
  
  speedM <- matrix(mapply(speedJumpDiffusion, lambdaMatrix, tMatrix),ncol=nt)
  
  return(speedM)
}

# speedJumpDiffusionMatrix(c(0.01,0.1,0.5,0.9,1,1.1,2,5,10,20,50),c(10,100,1000))
# par(mfrow=c(1,1))
# speed1 <- speedJumpDiffusionMatrix(seq(1,100,by=1),c(100))[,1]
# plot(speed1,xlab="lambda",ylab="temps de calcul",type='l',col='purple')

##########################################################################
# EFFICIENCY TESTS                                                       #
##########################################################################

cdfJumpDiffusionIncrements <- function(y,lambda){
  jumpDiffusionIncrementsTruncatedSeries <- function(x){
    jumpDiffusionIncrementsSeriesTerms <- function(k){
      return (dnorm(x,mean=0,sd=k+1)*dpois(k,lambda))
    }
    return(sum(sapply(0:(10*lambda),FUN="jumpDiffusionIncrementsSeriesTerms")))
  }
  return(integrate(Vectorize(jumpDiffusionIncrementsTruncatedSeries), -Inf, y)$value)
}

jumpDiffusionEfficiency <- function(lambda,t){
  simulatedProcess <- jumpDiffusionVector(lambda,20*t,t)
  simulatedIncrements <- diff(simulatedProcess,1)
  ks.test(simulatedIncrements, "cdfJumpDiffusionIncrements", lambda)
}

jumpDiffusionEfficiency(1,1000)