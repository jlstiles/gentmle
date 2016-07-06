library(tmle)
library(foreach)
library(doSNOW)
library(ggplot2)
library(SuperLearner)
library(parallel)
library(ggplot2)
library(reshape)
library(plyr)
library(origami)
library(geepack)

#
gendata.blip=function(n){
  U1 = runif(n,0,1)
  W1= -1*as.numeric(U1<.5)+1*as.numeric(U1>=.5)
  # W1=rnorm(n)
  W2=rnorm(n)
  W3=rbinom(n,1,.3)
  W4=rnorm(n)
  A=rbinom(n,1,g0(W1,W2,W3,W4))
  Y=rbinom(n,1,Q0(A,W1,W2,W3,W4))
  blip = Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4)
  return(list(df=data.frame(A,W1,W2,W3,W4,Y),blip=blip))
}


integral = function(blip,y){
  x=seq(-1,1,.0005)
  min=which(abs(y-x)==min(abs(y-x)))
  y = x[min]
  int.vals = vapply(x[1:min],FUN = function(a) (mean(blip>a)
                                               +mean(blip>(a+.0005))),
                    FUN.VALUE=1)
  sum(int.vals)*.0005/2
  }

# ttt = simPrim(0,10000)
# mean(ttt[[6]])

# n=1000
simPrim = function(y,n) {
  df = gendata.blip(n)
  data = df$df

  data1 = data
  data1$A = 1
  data0 = data
  data0$A = 0
  newdata = rbind(data,data0,data1)

  QAWfit = glm(Y~A+W1+W2+W3+W4+A*W4,data=data,family='binomial')
  QAW = predict(QAWfit,type='response')
  Q0W = predict(QAWfit, newdata = newdata[(n+1):(2*n),], type='response')
  Q1W = predict(QAWfit, newdata = newdata[(2*n+1):(3*n),], type='response')
  Q = cbind(QAW,Q0W,Q1W)

  Binit = Q[,"Q1W"]-Q[,"Q0W"]
  int.init = ifelse(Binit>y,y+1,Binit+1)
  psi.init = mean(int.init)

  gfit = glm(A~.,data=data[,1:5],family='binomial')
  g1W = predict(gfit,type='response')
  depsilon = .001
  gbounds <- Qbounds <- c(10^-9, 1-10^-9)

  onestepPrim <- oneStepPrim(Y = data$Y, A = data$A, Q = Q, g1W = g1W,
                             depsilon = depsilon, max_iter = max(1000, 2/depsilon),
                             gbounds = gbounds, Qbounds = Qbounds,y=y)
  psi = onestepPrim[[1]]
  SE = sqrt(onestepPrim[[2]])
  onestepPrimCI = c(psi-1.96*SE,psi+1.96*SE)
  cover = (truth>=onestepPrimCI[1])*(truth<=onestepPrimCI[2])
  # c(cover,psi.init, psi,onestepPrimCI)
  return(list(cover,psi.init,psi,onestepPrimCI,onestepPrim$epsilon,onestepPrim$IC))
}




#-----------------------------------One-Step TMLE for ATT parameter  ----------------------------------------
oneStepPrim <- function(Y, A, Q, g1W, y,depsilon, max_iter, gbounds, Qbounds){
  # max_iter=500
  # y=.03
  # Y=data$Y
  # A=data$A
  n <- length(Y)
  B = Q[,"Q1W"]-Q[,"Q0W"]
  int = ifelse(B>y,y+1,B+1)
  # int
  # B
  calcLoss <- function(Q) -mean(Y * log(Q[,"QAW"]) + (1-Y) * log(1 - Q[,"QAW"]))
  psi.prev <- psi  <- mean(int)
  H1.AW =  (B<=y)*(A/g1W-(1-A)/(1-g1W))
  IC.prev <- IC.cur <- H1.AW*(Y-Q[, "QAW"])+int-psi
  deriv <-  mean(H1.AW* (Y-Q[, "QAW"]))
  if (deriv > 0) { depsilon <- -depsilon}
  loss.prev <- Inf
  loss.cur <-  calcLoss(Q)
  if(is.nan(loss.cur) | is.na(loss.cur) | is.infinite(loss.cur)) {
    loss.cur <- Inf
    loss.prev <- 0
  }
  iter <-  0
  # max_iter=10
  # psi
  while (loss.prev > loss.cur & iter < max_iter){
    IC.prev <-  IC.cur
    Q.prev <- Q
    g1W <- .bound(g1W, gbounds)
    H1 <- cbind(HAW = (B<=y)*(A/g1W-(1-A)/(1-g1W)),
                H0W = -(B<=y)*(1-A)/(1-g1W),
                H1W = (B<=y)*A/g1W)
    Q  <- .bound(plogis(qlogis(Q.prev) - depsilon * H1), Qbounds)
    psi.prev <- psi
    # psi
    B = Q[,"Q1W"]-Q[,"Q0W"]
    int = ifelse(B>y,y+1,B+1)
    psi <- mean(int)
    loss.prev <- loss.cur
    loss.cur <- calcLoss(Q)
    IC.cur <- (B<=y)*(A/g1W-(1-A)/(1-g1W))*(Y-Q[, "QAW"])+
      int-psi
    if(is.nan(loss.cur) | is.infinite(loss.cur) | is.na(loss.cur)) {loss.cur <- Inf}
    iter <- iter + 1
    # print(psi.prev)
  }
  return(list(psi = psi.prev, var.psi = var(IC.prev)/n,  epsilon = (iter-1)*depsilon,IC=IC.prev))
}

#------------- bound function --------------
.bound <- function(x, bds){
  x[x > max(bds)] <- max(bds)
  x[x < min(bds)] <- min(bds)
  x
}

g0=function(W1,W2,W3,W4) {plogis(-.8*W1+.39*W2+.08*W3-.12*W4-.15)}
g0=function(W1,W2,W3,W4) {plogis(-.28*W1+1*W2+.08*W3-.12*W4-1)}
Q0=function(A,W1,W2,W3,W4) {plogis(3*A-1*W1*A+1.2*W2-1.9*W3+.4*W4*A*W2-cos(W4))}
Q0=function(A,W1,W2,W3,W4) {plogis(.5*A-1*W1+1.2*W2-1.9*W3+.4*W4)}
Q0=function(A,W1,W2,W3,W4) {plogis(.5*A-1*W1+1.2*W2-1.9*W3+.4*W4-A*W4)}

n=1000000
data.blip=gendata.blip(n)
hist(data.blip$blip,breaks=200)

y=0.1
n = 1000
truth = integral(data.blip$blip,y)
truth


detectCores()
cl = makeCluster(8, type = "SOCK")
registerDoSNOW(cl)

ptm = proc.time()
allresults=foreach(i=1:1000,
                   .packages=c("tmle","plyr","SuperLearner","geepack","BlipVariance"))%dopar%{simPrim(y,n)}
proc.time()-ptm

#calculate performance
results=do.call(rbind,allresults)
head(results,20)
mean(results[,1])

# results
# results = allresults
perf=function(ests,truth){
  n=length(ests)
  var=((n-1)/n)*var(ests)
  bias=mean(ests)-truth
  mse=mean((ests-truth)^2)
  c(var=var,bias=bias,mse=mse)
}

