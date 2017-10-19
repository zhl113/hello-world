
library(MASS)
##########################generating data############################
omega.f<-function(p, M, alpha, d){
  A=matrix(0,p,p)
  for (i in 2:p){
    
    for (j in 2:(i-1)){
      A[i,j]= M * abs(i-j)^{- alpha} - M * (abs(i-j)+1)^{- alpha}
      A[i,j]= A[i,j] + M * (abs(i-1))^{- alpha} / (i-1) /2
    }
    A[i,1]=M * (abs(i-1))^{- alpha} / (i-1) * i /2
  }
  
  A[row(A)==col(A)]=1
  D=diag(rep(d^{-1},p))
  omega=t(A)%*%D%*%A
  return (omega)}

#########################zhao method###################################
index <- function (n, p){
  if (n < 1) {return (1)}
  if (n > p) {return (p)}
  return (n)
}
restrict <- function (m, eta){
  for (i in 1: length(m))
  {
    if (m[i] < eta^-1) m[i]=eta^-1
    if (m[i] > eta) m[i]=eta
  }
  return(m)
}
proj <- function (m, eta){
  s <- eigen(m)
  D <- diag(restrict(s$values, eta))
  return (s$vectors %*% D %*% solve(s$vectors))
}


tapering.k.est <-function(k, data, eta){
  if (k < 2) { k = 1 }
  p=dim(data)[2]
  n=dim(data)[1]
  est=matrix(0,p,p)
  for (i in (2-k):p){
    subindex=c(index((i-k),p):index((i+2*k-1),p))
    subset=data[,subindex]
    hat.cov= 1/n* t(subset)%*%subset
    tilde.cov=proj(hat.cov, eta)
    #tilde.cov = hat.cov
    center=c((-index(i-k,p)+index(i,p)+1):(index(i+k-1,p)-index(i-k,p)+1))
    tilde.pre=matrix(0,p,p)
    tilde.pre[c(index(i,p):(index(i,p)+length(center)-1)),c(index(i,p):(index(i,p)+length(center)-1))]=solve(tilde.cov)[center,center]
    est=est+tilde.pre
  }
  return(est)
}

est.z<-function(k, data, eta){
  estimator=1/(k-floor(k/2))*(tapering.k.est(k, data,eta)-tapering.k.est(floor(k/2), data,eta))
  return(estimator)
}
###############################bickle method#####################################
index.f<-function(index,p){
  if (index > p) return(p)
  if (index < 1) return(1)
  return (index)
}

mse.f<-function(lm,n,k){
  return (sum(lm$residuals^2)/(n-k))
}

est.b <- function(k, data){
  n=dim(data)[1]
  p=dim(data)[2]
  A <- matrix(0,p,p)
  A[row(A)==col(A)] <- 1
  D=rep(0,p)
  D[1]=var(data[,1])/(n-1)*n
  for (i in 2:p){
    y <- data[,i]
    x <- data[,index.f(i-k, p):index.f(i-1, p)]
    lm <- lm(y ~ x - 1)
    D[i] <- mse.f(lm, n, (index.f(i-1, p)-index.f(i-k, p)))
    coeff <- lm$coefficients
    temp <- matrix(0, p, p)
    temp[i,index.f(i-k, p):index.f(i-1, p)]=-coeff
    A <- A + temp
  }
  D=diag(D^(-1))
  return (t(A)%*%D%*%A)}
###############################bandwidth calculator##########################
bd.z.f<-function (n,alpha) {return (floor(n^{1/(2*alpha)}))}
bd.b.f<-function (n,alpha,p) {return (floor((n/log(p))^{1/(2*alpha+2)}))}
bd.e.f<-function (n,alpha) {return (floor(n^{1/(2*alpha+1)}))}

##############################lets play#######################################
dist.f <- function(matrix){ return (max(abs(eigen(matrix)$values)))}

eta=10000
p.list=c(250)
alpha.list=c(0.8)
M.list=c(3)
d=1
n.list=c(1000)
total=1
begin = 5
interval = 10

# eta=10000
# p=200
# alpha=0.6
# M=2
# d=1
# n=100
# total=5

# eta=10000
# p.list=c(250,500,1000,2000)
# alpha.list=c(0.6,0.7,0.8,0.9,1)
# M.list=c(0.8,1,1.2)
# d=1
# n.list=c(250,500,1000,2000)
# total=100

result=rep(0,12)
for (p in p.list){
  for (alpha in alpha.list){
    for (M in M.list) {
      for (n in n.list){
          start.step.time = Sys.time()
          cat('alpha', alpha, 'M', M, 'n', n, 'p', p, '\n')
        omega=omega.f(p, M, alpha, d)
        sigma=solve(omega)
        bd.z=bd.z.f(n,alpha)
        bd.b=bd.b.f(n,alpha,p)
        bd.e=bd.e.f(n,alpha)
        omega.eigen=dist.f(omega)
        
        bd.z.m.z1=rep(0,total)
        bd.b.m.z1=rep(0,total)
        bd.b.m.b1=rep(0,total)
        bd.e.m.z1=rep(0,total)
        bd.p1=rep(0,total)
        bd.inv1=rep(0,total)
        sequ = seq(begin, 200, interval)
        # sequ = seq(begin, p + begin*2, interval)
        bd.m1 = matrix(0,total,length(sequ))
        for (k in 1:total){
          data=mvrnorm(n,rep(0,p),sigma)
          bd.z.m.z1[k]=dist.f(est.z(bd.z, data, eta)-omega)
          # bd.b.m.z1[k]=dist.f(est.z(bd.b, data, eta)-omega)
          # bd.b.m.b1[k]=dist.f(est.b(bd.b, data)-omega)
          bd.e.m.z1[k]=dist.f(est.z(bd.e, data, eta)-omega)
          # bd.p1[k]=dist.f(est.z(p, data, eta)-omega)
          # bd.inv1[k]=dist.f(est.z(2*p, data, eta)-omega)
          for (j in sequ){
            u = (j - begin) / interval + 1
            bd.m1[k,u]=dist.f(est.z(j, data, eta)-omega)
          }
          cat(omega.eigen, bd.z, bd.e, bd.b, bd.z.m.z1[k], bd.e.m.z1[k], bd.p1[k], bd.inv1[k],bd.b.m.z1[k], bd.b.m.b1[k], '\n')
        }
        
        bd.z.m.z=mean(bd.z.m.z1)
        bd.b.m.z=mean(bd.b.m.z1)
        bd.b.m.b=mean(bd.b.m.b1)
        bd.e.m.z=mean(bd.e.m.z1)
        bd.p=mean(bd.p1)
        bd.inv=mean(bd.inv1)
        bd.m = colMeans(bd.m1)
        
        current=c(p, n, alpha, M, omega.eigen, bd.z, bd.e, bd.b, bd.z.m.z, bd.e.m.z, bd.b.m.z, bd.b.m.b)
        result=rbind(result,current)
        
          end.step.time = Sys.time()
          cat(end.step.time - start.step.time, '\n')
          
      }}}}
colnames(result)=c("p", "n", "alpha", "M", "omega.eigen" ,  "bd.z", "bd.e", "bd.b",   "bd.z.m.z", "bd.e.m.z", "bd.b.m.z", "bd.b.m.b")
result=result[-1,]
display=round(result,4)
display

whatineed=display[display[,9] < display[,10], 2:4]
whatineed
#write.csv(display,'result.csv')
#################################test test###############################
