
#################################################
# Application Poisson composÃ©                  #
#################################################


#Générateur d'un actif

generateur<- function(lambda,mu,Tlim,tlim,init){
  N <- poissonProcessPointwise(lambda,Tlim)
  U <- runif(n=N,min=0,max=Tlim)
  Y <- rnorm(n=N,mean=0,sd=mu)
  sortedU <- sort(x=U,index.return = TRUE)
  process <- c(init)
  sum <- 0
  i <- 1
  for (k in 1:tlim){
    while (sortedU$x[i] < k){
      sum <- sum + Y[sortedU$ix[i]]
      i <- i +1
    }
    process <- c(process,init*exp(sum))
  }
  process<-process[1:(length(process)-1)]
  return (process)
}
###############################################################

#Simulation action de la forme S(t)=S(0)*exp(L(t)) 
# on essaye de simuler un index type cac 40 (on commence à 2990) 

tlim=365*100
Tlim=80000

x<-1:tlim
y<-generateur(1,0.004,Tlim,tlim,3880)  #Paramètre estimé auparavant
plot(x,y,type='l')

##################################################################

#Pricer d'options asiatique(call) avec Poisson composé 


#fonction pricing

Pricer<-function(T,init){   #init=S(0) T€[0,mat] (matu résid) K=strike
prix<-max(sum(generateur(1,0.005,2*T,T,init))/T-K,0)
for (i in 1:(iter-1)){
prix<-prix+max(sum(generateur(1,0.005,2*T,T,init))/T-K,0)
}
prix<-prix/iter
prix<-exp(-taux*T/100)*prix
return(prix)
}
##############################################################
#fonction pricing avec reduction de variance

generateurv<- function(lambda,mu,Tlim,tlim,init){
  N <- poissonProcessPointwise(lambda,Tlim)
  U <- runif(n=N,min=0,max=Tlim)
  Y <- rnorm(n=N,mean=0,sd=mu)
  sortedU <- sort(x=U,index.return = TRUE)
  process1 <- c(init)
  process2 <- c(init)
  sum1 <- 0
  sum2 <-0
  i <- 1
  for (k in 1:tlim){
    while (sortedU$x[i] < k){
      sum1 <- sum1 + Y[sortedU$ix[i]]
      sum2<- sum2 -Y[sortedU$ix[i]]
      i <- i +1
    }
    process1 <- c(process1,init*exp(sum1))
    process2 <- c(process2,init*exp(sum2))

  }
  process1<-process1[1:(length(process1)-1)]
  process2<-process2[1:(length(process2)-1)]
  process<-cbind(process1,process2)
  return (process)
}

Pricerv<-function(T,init){   #init=S(0) T€[0,mat] (matu résid) K=strike
prix<-(max(sum(generateurv(1,0.005,2*T,T,init)[,1])/T-K,0)+max(sum(generateurv(1,0.001,2*T,T,init)[,2])/T-K,0))/2
for (i in 1:(iterv-1)){
prix<-prix+(max(sum(generateurv(1,0.005,2*T,T,init)[,1])/T-K,0)+max(sum(generateurv(1,0.001,2*T,T,init)[,2])/T-K,0))/2
}
prix<-prix/iterv
prix<-exp(-taux*T/100)*prix
return(prix)
}


# test , on veut pricer call de maturité 30 strike =50 
K=50
#la maturité est en jour
mat<-10
taux<-0.05
iterv=10  #calcul de l'espérance avec N=100
iter=10

x<-1:1000
y<-1:100
z<- matrix(data = NA, nrow = 100, ncol = 1000)
for (i in 1:100){
for (j in 1:1000){
z[i,j]=Pricerv(1000*i,j)
}
}

x<-1:1000
y<-1:100
f<- matrix(data = NA, nrow = 100, ncol = 1000)
for (i in 1:100){
for (j in 1:1000){
f[i,j]=Pricer(1000*i,j)
}
}




op<-par(mfrow=c(1,2))
persp(x = seq(0, 1, length.out = nrow(z)),
      y = seq(0, 1, length.out = ncol(z)),
       z,xlab ="Maturité" , ylab = "Sous-jacent", zlab ="Prix de l'option",
      main = "Surface de pricing avec réduction variance", theta = 30, phi=10,
, col="blue",border="red" ,shade=0.5)
persp(x = seq(0, 1, length.out = nrow(f)),
      y = seq(0, 1, length.out = ncol(f)),
       f,xlab ="Temps" , ylab = "Sous-jacent", zlab ="Prix de l'option",
      main = "Surface de pricing d'un call asiatique",sub="sans reduction de variance
(n=100)", theta = 30, phi=10,
, col="blue",border="red" ,shade=0.5)
par(op)




#calcul du gain de la variance

var1<-matrix(data = NA, nrow = 100, ncol = 1000)
var2<-matrix(data = NA, nrow = 100, ncol = 1000)



for (i in 1:100){
for (j in 1:1000){
c<-max(sum(generateur(1,0.005,2*100*i,100*i,j))/T-K,0)
var1[i,j]<-(c-f[i,j])^2 
var2[i,j]<-(c-z[i,j])^2
for (h in 1:99){
c<-max(sum(generateur(1,0.005,2*100*i,100*i,j))/T-K,0)
var1[i,j]=var1[i,j]+(c-f[i,j])^2   #v1 doit >v2
var2[i,j]=var2[i,j]+(c-z[i,j])^2
}
var1<-var1/100
var2<-var2/100
}
}

gain<-sum(var1-var2)/sum(var1)