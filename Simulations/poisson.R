##########################################################################
# PROCESS SIMULATION                                                     #
##########################################################################

# Simulate n exponential random values
countingList <- function(lambda,n){
  c <- rexp(n, rate = lambda)
  c <- cumsum(c)
  return(c)
}

# The computational time of simulation grows with lambda, this is due to the simulation
# of exponential random values

# l'espérance du dernier élément de la liste simulée est N/lambda, et l'écart-type
# de sqrt(N)/lambda), et pour construire les vecteurs on veut que tlim soit plus petit
# que ce dernier élément de la liste, pour économiser la mémoire, on pourra donc
# choisir tlim = Nlim/lambda - 10sqrt(Nlim)/lambda, ce qui devrait être suffisant,
# soit Nlim = (10 + sqrt(100 + 4*lambda*tlim))^2 / 4 

computeNlim <- function(lambda,tlim){
  return ((10 + sqrt(100 + 4*lambda*tlim))^2)/4
}

# Now simulate Poisson process from definition

poissonProcessPointwise <- function(lambda,t){
  N <- computeNlim(lambda,t)
  c <- countingList(lambda,N)
  cpt = 0
  for (k in 1:N){
    if (c[k]<=t) cpt = cpt + 1
  }
  return(cpt)
}

# Simulate a vector of Poisson process values N(t) for plotting

poissonProcessVector <- function(lambda,tlim){
  Nlim <- computeNlim(lambda,tlim)
  c = countingList(lambda,Nlim)  
  poissonVector <- c()
  n <- 0
  for (t in 1:tlim){
    while (c[n+1]<=t){
      n <- n+1
    }
    poissonVector <- c(poissonVector,n)
  }
  return(poissonVector) # This vector contains values [N(1),...,N(tlim)]
}

jumpTimes <- function(lambda,tSub){
  Nlim <- computeNlim(lambda,tSub[length(tSub)])
  c = countingList(lambda,Nlim)  
  poissonVector <- c()
  jumpTimes <- c()
  n <- 0
  for (t in tSub){
    while (c[n+1]<=t){
      n <- n+1
    }
    poissonVector <- c(poissonVector,n)
    if (length(poissonVector) > 2){
      if (poissonVector[length(poissonVector)-1]<poissonVector[length(poissonVector)]){
        jumpTimes <- c(jumpTimes,t)
      }   
    }
  }
  return(jumpTimes)
}

# Plot the Poisson process as a function of t

#par(mfrow=c(2,2))
#lambda = 10
#tlim <- 100
# plot(x = 1:tlim, y=poissonProcessVector(lambda,tlim), type = "s", xlab = "temps t", ylab = "N(t)")  

##########################################################################
# SPEED TESTS                                                            #
##########################################################################

speedPoisson <- function(lambda,t){
  return(system.time(poissonProcessVector(lambda,t))[3])
}

speedPoissonMatrix <- function(lambda,t){
  
  nt = length(t)
  nlambda = length(lambda)
  
  lambdaRows <- replicate(nt,lambda)
  tCols <- replicate(nlambda,t,simplify=FALSE)
  
  lambdaMatrix <- matrix(lambdaRows,nlambda)
  tMatrix <- do.call(rbind,tCols)
  
  speedM <- matrix(mapply(speedPoisson, lambdaMatrix, tMatrix),ncol=nt)
  
  return(speedM)
}

# speedPoissonMatrix(c(0.01,0.1,0.5,0.9,1,1.1,2,5,10,20,50),c(10,100,1000,10000))
#speed <- speedPoissonMatrix(c(1,10),seq(1000,10000,by=100))
#par(mfrow=c(2,1))
#fast <- speed[1,]
#slow <- speed[2,]
#plot(x=seq(1000,10000,by=100),fast,type='l',col='blue',xlab="t",ylab="temps de calcul")
#plot(x=seq(1000,10000,by=100),slow,type='l',col='red',xlab="t",ylab="temps de calcul")
#speed <- speedPoissonMatrix(seq(1,100,by=1),c(100))[,1]
#par(mfrow=c(1,1))
#plot(speed,xlab="lambda",ylab="temps de calcul",type='l',col='purple')

##########################################################################
# EFFICIENCY TESTS                                                       #
##########################################################################

poissonEfficiency <- function(lambda,t){
  simulatedProcess <- poissonProcessVector(lambda,t)  
  simulatedIncrements <- diff(simulatedProcess,1)
  ks.test(simulatedIncrements, "ppois",lambda)
}

#pvalue <- c()
#for (i in 1:100){
#  pvalue[i] <- poissonEfficiency(4,1000)$p.value
#}
#length(pvalue[pvalue>0.05]) normal car la p-value suit une loi uniforme
#ks.test(pvalue,'punif')