##########################################################################
# PROCESS SIMULATION                                                     #
##########################################################################

# Gamma(alpha,beta)(x) = beta^alpha/Gamma(alpha) x^{alpha-1} e^{-beta*x}
# Un processus Gamma est un processus de Lévy à incréments iid de loi Gamma
# Simuler le processus revient à simuler les incréments
# Pour tout t, X_t a une loi Gamma(ct,1/lambda) = Gamma(alpha,beta)
# Propriété d'échelle : si X suit Gamma(alpha,beta) alors beta*X suit Gamma(alpha,1)
# On simule donc des Gamma(alpha,1), de densité x^{a-1}/Gamma(a) * e^{-x}

# Algorithme 1 : à utiliser si a <= 1

gammaVariableInf1 <- function(n,a){
  res <- c()
  for (i in 1:n){
    x <- y <- 1
    while (x+y>1){
      u <- runif(1,0,1)
      v <- runif(1,0,1)
      x <- u^{1/a}
      y <- v^{1/(1-a)}
    }
    e <- rexp(1,1)
    res <- c(res,x*e/(x+y))
  }
  return(res)
}

# Algorithme 2 : si a > 1 

gammaVariableSup1 <- function(n,a){
  res <- c()
  b <- a - 1
  c <- 3*a - 0.75
  for (i in 1:n){
    repeat{
      u <- runif(1,0,1)
      v <- runif(1,0,1)
      w <- u*(1-u)
      y <- sqrt(c/w)*(u-0.5)
      x <- b + y
      if (x<0) next
      z <- 64*(w^3)*(v^3)
      if (log(z)<2*(b*log(x/b)-y)) break 
    }
    res <- c(res,x)
  }
  return(res)
}

# Cette fonction rassemble les deux précédentes

gammaVariable <- function(n,a){
  case = a>1 
  if (case){
    res <- gammaVariableSup1(n,a)
  }
  else{
    res <- gammaVariableInf1(n,a)
  }
  return(res)
}
# par(mfrow=c(2,1))
# x <- gammaVariable(10000,3)
# xmax = floor(max(x))+1
# hist(x,breaks=seq(0,xmax,by=0.05),freq=F,xlim = c(0,xmax),las=1)
# lines(seq(0,xmax,by=0.01),dgamma(seq(0,xmax,by=0.01),shape = 3, rate = 1),col = 'red')
# ks.test(x,"pgamma",shape=3)

# Rappelons maintenant que si (X_t) suit un processus Gamma(t;alpha,beta) alors
# ses incréments suivent la loi Gamma(alpha^2/beta,beta/alpha)

gammaProcess <- function(t,alpha,beta){
  a <- alpha^2/beta
  b <- beta/alpha
  process <- c(0)
  sum <- 0
  for (i in 1:t){
    increment <- gammaVariable(1,a)/b
    sum <- sum + increment
    process <- c(process,sum)    
  }
  return (process)
}

# par(mfrow=c(2,2))
# plot(1:101,gammaProcess(100,10,10), type = "s", xlab = "temps t", ylab = "X(t)")

##########################################################################
# SPEED TESTS                                                            #
##########################################################################

speedGammaProcess <- function(alpha,beta,t){
  return (system.time(gammaProcess(t,alpha,beta))[3])
}

speedGammaProcessMatrix <- function(alpha,beta,t){
  
  nbeta = length(beta)
  nalpha = length(alpha)
  
  alphaRows <- replicate(nbeta,alpha)
  betaCols <- replicate(nalpha,beta,simplify=FALSE)
  
  alphaMatrix <- matrix(alphaRows,nalpha)
  betaMatrix <- do.call(rbind,betaCols)
  
  nt <- length(t)
  result <- vector("list",nt)
  
  for (i in 1:nt){ 
    speedGammaProcessFixedt <- function(alpha,beta){
      return(speedGammaProcess(alpha,beta,t[i]))
    }
    speedM <- matrix(mapply(speedGammaProcessFixedt, alphaMatrix, betaMatrix),ncol=nbeta)
    result[[i]] <- speedM
  }
  
  return(result)
}

#speedGammaProcessMatrix(c(0.01,0.1,0.5,0.9,1,1.1,2,5,10,20,50),
#                        c(0.01,0.1,0.5,0.9,1,1.1,2,5,10,20,50),
#                        c(10000)
#                        )
# par(mfrow=c(1,1))
# speed1 <- speedGammaProcessMatrix(seq(0.1,5.1,by=0.5),seq(0.1,5.1,by=0.1),c(10))[,1]
# plot(speed1,xlab="alpha",ylab="temps de calcul",type='l',col='purple')

##########################################################################
# EFFICIENCY TESTS                                                       #
##########################################################################

gammaEfficiency <- function(alpha,beta,t){
  simulatedProcess <- gammaProcess(t,alpha,beta)
  simulatedIncrements <- diff(simulatedProcess,1)
  ks.test(simulatedIncrements, "pgamma", shape=alpha^2/beta, rate=beta/alpha)
}
