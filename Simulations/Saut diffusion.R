#################################################
# Application Saut diffusion                    # 
#################################################


#Générateur d'un actif (à partir d'un seuil)

generateursd<- function(lambda,sigma1,sigma2,Tlim,tlim,init){
  N <- poissonProcessPointwise(lambda,Tlim) # Tlim = 2-10 * tlim
  U <- runif(n=N,min=0,max=Tlim)
  Y <- rnorm(n=N,mean=0,sd=sigma1)
  sortedU <- sort(x=U,index.return = TRUE)
  process <- c(init)
  sum <- 0
  i <- 1
  for (k in 1:tlim){
    while (sortedU$x[i] < k){
      sum <- sum + Y[sortedU$ix[i]]
      i <- i +1
    }
    G <- rnorm(n=1,mean=0,sd=sigma2)
    sum <- sum + G[1]
    process <- c(process,init*exp(sum))
  }
  process <- process[1:(length(process)-1)]

  return (process)
}



###########################################################################

#Pricing d'une option (call) lookback (Strike fixe)

Pricersd<-function(T,init){   #init=S(0) T=maturité (matu résid) K=strike
prix<-max(max(generateursd(1,sigma1,sigma2,2*T,T,init))-K,0)
for (i in 0:(iter-1)){
prix<-prix+max(max(generateursd(1,sigma1,sigma2,2*T,T,init))-K,0)
}
prix<-prix/iter
prix<-exp(-taux*T/100)*prix
return(prix)
}

# test , on veut pricer call de maturité 30 strike =50 
K=50
#la maturité est en jour
mat<-5
taux<-0.05
sigma1=0.01
sigma2=0.01
iterv=10  #calcul de l'espérance avec N=100
iter=10

x<-1:1000
y<-1:100
yt<-10*y
z<- matrix(data = NA, nrow = 10, ncol = 100)
for (i in 1:100){
for (j in 1:100){
z[i,j]=Pricersd(i*100,j)
}
}

persp(x = seq(0, 1000, length.out = nrow(z)),
      y = seq(0, 100, length.out = ncol(z)),
       z,xlab ="Maturité" , ylab = "Sous-jacent", zlab ="Prix de l'option",
      main = "Surface de pricing d'un Call Lookback",
sub="sans reduction de variance
(n=100)", theta = 30, phi=10,
, col="green",shade=0.75)

########################################################################

#Reduction de la variance :

generateursdv<- function(lambda,sigma1,sigma2,Tlim,tlim,init){
  N <- poissonProcessPointwise(lambda,Tlim) # Tlim = 2-10 * tlim
  U <- runif(n=N,min=0,max=Tlim)
  Y <- rnorm(n=N,mean=0,sd=sigma1)
  sortedU <- sort(x=U,index.return = TRUE)
  process11 <- c(init)
  process12 <- c(init)
  process21 <- c(init)
  process22 <- c(init)    
  sum11 <- 0
  sum12 <- 0
  sum21 <- 0
  sum22 <- 0
  i <- 1
  for (k in 1:tlim){
    while (sortedU$x[i] < k){
      sum11 <- sum11 + Y[sortedU$ix[i]]
      sum12 <- sum12 + Y[sortedU$ix[i]]
      sum21 <- sum21 - Y[sortedU$ix[i]]
      sum22 <- sum22 - Y[sortedU$ix[i]]
      i <- i +1
    }
    G <- rnorm(n=1,mean=0,sd=sigma2)
    sum11 <- sum11 + G[1]
    sum12<-sum12 - G[1]
    sum21 <- sum21 + G[1]
    sum22<-sum22 - G[1]
    process11 <- c(process11,init*exp(sum11))
    process12 <- c(process12,init*exp(sum12))
    process21 <- c(process21,init*exp(sum21))
    process22 <- c(process22,init*exp(sum22))
  }
  process11 <- process11[1:(length(process11)-1)]
  process12 <- process12[1:(length(process12)-1)]
  process21 <- process21[1:(length(process21)-1)]
  process22 <- process22[1:(length(process22)-1)]
  process<-cbind(process11,process12,process21,process22)
  return (process)
}

Pricersdv<-function(T,init){   #init=S(0) T=maturité (matu résid) K=strike
prix<-(max(max(generateursdv(1,sigma1,sigma2,2*T,T,init)[,1])-K,0)+
max(max(generateursdv(1,sigma1,sigma2,2*T,T,init)[,2])-K,0)+
max(max(generateursdv(1,sigma1,sigma2,2*T,T,init)[,3])-K,0)+
max(max(generateursdv(1,sigma1,sigma2,2*T,T,init)[,4])-K,0))/4
for (i in 0:(iterv-1)){
prix<-prix+(max(max(generateursdv(1,sigma1,sigma2,2*T,T,init)[,1])-K,0)+
max(max(generateursdv(1,sigma1,sigma2,2*T,T,init)[,2])-K,0)+
max(max(generateursdv(1,sigma1,sigma2,2*T,T,init)[,3])-K,0)+
max(max(generateursdv(1,sigma1,sigma2,2*T,T,init)[,4])-K,0))/4
}
prix<-prix/iter
prix<-exp(-taux*T/100)*prix
return(prix)
}

#la maturité est en jour
mat<-5
taux<-0.05
sigma1=0.01
sigma2=0.01
iterv=10000 #calcul de l'espérance avec N=100
iter=10000

x<-1:1000
y<-1:100
yt<-10*y
f<- matrix(data = NA, nrow = 10, ncol = 100)
for (i in 1:100){
for (j in 1:1000){
f[i,j]=Pricersdv(i*100,j)
}
}

persp(x = seq(0, 1000, length.out = nrow(z)),
      y = seq(0, 100, length.out = ncol(z)),
       f,xlab ="Maturité" , ylab = "Sous-jacent", zlab ="Prix de l'option",
      main = "Surface de pricing Call lookback",
sub="avec reduction de variance", theta = 30, phi=10,
, col="white",border="black" ,shade=0.5)


#######################################################

 #Calcul du gain de variance

#comparaison graphique :

op<-par(mfrow=c(1,2))
persp(x = seq(0, 1, length.out = nrow(z)),
      y = seq(0, 1, length.out = ncol(z)),
       z,xlab ="Maturité" , ylab = "Sous-jacent", zlab ="Prix de l'option",
      main = "Surface de pricing d'un Call Lookback",sub="Sans reduction de variance
(n=20)", theta = 60, phi=10,
, col="green" ,shade=0.75)
persp(x = seq(0, 1, length.out = nrow(f)),
      y = seq(0, 1, length.out = ncol(f)),
       f,xlab ="Maturité" , ylab = "Sous-jacent", zlab ="Prix de l'option",
      main = "Surface de pricing d'un Call Lookback"
,sub="Avec reduction de variance
(n=20)", theta = 60, phi=10,
, col="green",shade=0.75)
par(op)


#calcul du gain :


var1<-matrix(data = NA, nrow = 100, ncol = 1000)
var2<-matrix(data = NA, nrow = 100, ncol = 1000)



for (i in 1:100){
for (j in 1:1000){
c<-max(sum(generateur(1,0.005,2*100*i,100*i,j))/T-K,0)
var1[i,j]<-(c-f[i,j])^2 
var2[i,j]<-(c-z[i,j])^2
for (h in 1:9){
c<-max(sum(generateur(1,0.005,2*100*i,100*i,j))/T-K,0)
var1[i,j]=var1[i,j]+(c-f[i,j])^2   #v1 doit >v2
var2[i,j]=var2[i,j]+(c-z[i,j])^2
}
var1<-var1/10
var2<-var2/10
}
}

gain<-sum(var1-var2)/sum(var1)

