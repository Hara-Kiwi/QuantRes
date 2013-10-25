source('/home/darkil73/Documents/Dropbox/Cours-ENSAE/Simulation et Monte Carlo/Projet/poisson.R')

##########################################################################
# PROCESS SIMULATION                                                     #
##########################################################################

# First algorithm (algorithm 6.1)

coumpoundPoissonPointwise1 <- function(lambda, t){
  sum <- 0
  Y <- c()
  while (sum<t){
    newExp <- rexp(n = 1, rate = lambda)
    sum <- sum + newExp
    Y <- c(Y,rnorm(n = 1, mean = 0, sd = 1))
  }
  return (sum(Y))
}

compoundPoissonVector1 <- function(lambda, Nlim, tlim){
  c <- countingList(lambda,Nlim)
  jumps <- rnorm(n = Nlim, mean = 0, sd = 1)
  process <- c()
  n <- 0
  m <- 0
  sum <- 0
  for (t in 1:tlim){
    while (c[n+1]<=t){
      n <- n+1
    }
    if (n != m){
      sum <- sum + sum(jumps[(m+1):n])
      m <- n
    }
    process <- c(process,sum)
  }
  return(process)
}

# tlim <- 100
# Nlim <- Nlim = computeNlim(tlim)
# plot(x = 1:tlim, y=compoundPoissonVector1(Nlim,tlim), type = "s", xlab = "temps t", ylab = "X(t)")  

# Second algorithm (algorithm 6.2) utilisant des va uniformes

compoundPoissonPointwise2 <- function(lambda,t,tlim){
  N <- rpois(1,lambda*tlim)
  U <- runif(n=N,min=0,max=tlim)
  Y <- rnorm(n=N,mean=0,sd=1)
  sum <- 0
  for (k in 1:N){
    if (U[k]<t){
      sum <- sum + Y[k]
    }
  }
  return (sum)  
}

compoundPoissonVector2 <- function(lambda,Tlim,tlim){
  N <- poissonProcessPointwise(lambda,Tlim)
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
    process <- c(process,sum)
  }
  return (process)
}

#par(mfrow=c(2,2))
#lambda=10
#tlim <- 1000
#plot(x = 1:tlim, y=compoundPoissonVector2(lambda,2*tlim,tlim), type = "s", xlab = "temps t", ylab = "X(t)")

##########################################################################
# SPEED TESTS                                                            #
##########################################################################

speedCompoundPoisson <- function(lambda,t,method=c("algo1","algo2")){
  switch(
    method,
    "algo1" = system.time(compoundPoissonVector1(lambda,computeNlim(lambda,t),t))[3],
    "algo2" = system.time(compoundPoissonVector2(lambda,2*t,t))[3]
  )
}

speedCompoundPoissonMatrix <- function(lambda,t,method=c("1","2")){
  
  nt = length(t)
  nlambda = length(lambda)
  
  lambdaRows <- replicate(nt,lambda)
  tCols <- replicate(nlambda,t,simplify=FALSE)
  
  lambdaMatrix <- matrix(lambdaRows,nlambda)
  tMatrix <- do.call(rbind,tCols)
  
  speedM <- matrix(mapply(speedCompoundPoisson, lambdaMatrix, tMatrix, method),ncol=nt)
  
  return(speedM)
}

# speedCompoundPoissonMatrix(c(0.01,0.1,0.5,0.9,1,1.1,2,5,10,20,50),c(10,100,1000,10000),"algo2")
# par(mfrow=c(2,1))
# speed1 <- speedCompoundPoissonMatrix(seq(1,100,by=1),c(100),"algo1")[,1]
# plot(speed1,xlab="lambda",ylab="temps de calcul",type='l',col='purple')
# speed2 <- speedCompoundPoissonMatrix(seq(1,100,by=1),c(100),"algo2")[,1]
# plot(speed2,xlab="lambda",ylab="temps de calcul",type='l',col='purple')

##########################################################################
# EFFICIENCY TESTS                                                       #
##########################################################################

cdfCompoundPoissonIncrements <- function(y,lambda){
  coumpoundPoissonIncrementsTruncatedSeries <- function(x){
    compoundPoissonIncrementsSeriesTerms <- function(k){
      return (dnorm(x,mean=0,sd=k)*dpois(k,lambda))
    }
    return(sum(sapply(0:(10*lambda),FUN="compoundPoissonIncrementsSeriesTerms")))
  }
  return(integrate(Vectorize(coumpoundPoissonIncrementsTruncatedSeries), -Inf, y)$value)
}

compoundPoissonEfficiency <- function(lambda,t,method=c("algo1","algo2")){
  simulatedProcess <- switch(
                        method,
                        "algo1" = compoundPoissonVector1(lambda,computeNlim(lambda,t),t),
                        "algo2" = compoundPoissonVector2(lambda,10*t,t)
                      )
  simulatedIncrements <- diff(simulatedProcess,1)
  ks.test(simulatedIncrements, "cdfCompoundPoissonIncrements", lambda)
}

compoundPoissonEfficiency(1,1000,method=c("algo1"))