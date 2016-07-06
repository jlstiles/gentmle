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
library(mvtnorm)

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

truth = vapply(1:length(truth),FUN = function(x) truth[[x]],FUN.VALUE=1)

n=1000
cc = 1:51
yy = y[cc]
test = simPrimMD(yy,n)
test[[2]]$r
test[[2]]$epsilon

apply(test[[2]]$IC,2,mean)

mcmc = rmvnorm(100000,rep(0,length(yy)),sigma  = test[[2]]$r,method="chol")
maxes = apply(mcmc,1,FUN = function(x) max(abs(x)))

zscore = quantile(maxes,.95)
zscore

SEs = apply(test[[2]]$IC,2,sd)/sqrt(n)
halfwidth = zscore*SEs

truths = truth[cc]
testpsi = vapply(1:length(truths),FUN = function(x) test[[2]]$psi[[x]],FUN.VALUE=3)
testpsi

sqrt(sum(apply(test[[2]]$IC,2,mean)^2))
lower95CI = testpsi-halfwidth
upper95CI = testpsi+halfwidth

plotter = data.frame(x=rep(yy,4),y=c(truths,testpsi,lower95CI,upper95CI),
                     type=c(rep("truth",length(truths)),rep("est",length(truths)),
                            rep("lower95CI",length(truths)),rep("upper95CI",length(truths))))

peter = ggplot(plotter,aes(x=x,y=y,colour=type))+geom_line()
peter = peter + labs(x="blip",y="primitive",title="primitive curve")
peter

# n=1000
simPrimMD = function(y,n) {
  # n=1000
  # y = c(-.1,0,.1,.2,.3)
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
  
  int.init = vapply(y,FUN = function(y) ifelse(Binit>y,y+1,Binit+1), FUN.VALUE=rep(1,n))
  psi.init = apply(int.init,2,mean)

  gfit = glm(A~.,data=data[,1:5],family='binomial')
  g1W = predict(gfit,type='response')
  depsilon = .001
  gbounds <- Qbounds <- c(10^-9, 1-10^-9)

  onestepPrim <- oneStepPrimMD(Y = data$Y, A = data$A, Q = Q, g1W = g1W,
                             depsilon = depsilon, max_iter = max(1000, 10/depsilon),
                             gbounds = gbounds, Qbounds = Qbounds,y=y)
  # psi = onestepPrim[[1]]
  # SE = sqrt(onestepPrim[[2]])
  # onestepPrimCI = c(psi-1.96*SE,psi+1.96*SE)
  # cover = (truth>=onestepPrimCI[1])*(truth<=onestepPrimCI[2])
  # c(cover,psi.init, psi,onestepPrimCI)
  return(list(psi.init,onestepPrim))
}
# 
# 
# 
# 
# #-----------------------------------One-Step TMLE for ATT parameter  ----------------------------------------
oneStepPrimMD <- function(Y, A, Q, g1W, y,depsilon, max_iter, gbounds, Qbounds){
  # max_iter=500
  # Y = data$Y
  # A = data$A
  # Q = Q
  # g1W = g1W
  # depsilon = .001
  # max_iter = max(1000, 2/depsilon)
  # gbounds = gbounds
  # Qbounds = Qbounds
  
  n <- length(Y)
  B = Q[,"Q1W"]-Q[,"Q0W"]
  
  int = vapply(y, FUN = function(y) ifelse(B>y,y+1,B+1),FUN.VALUE = rep(1,n))
  psi = apply(int,2,mean)
  psi.prev <- psi
  
  # B
  calcLoss <- function(Q) -mean(Y * log(Q[,"QAW"]) + (1-Y) * log(1 - Q[,"QAW"]))
  
  H1 =  vapply(y,FUN = function(y) (B<=y)*(A/g1W-(1-A)/(1-g1W)),FUN.VALUE=rep(1,n))
  # head(H1)
  
  # matrix gotten from previous step times Y-Q--not matrix mult
  deriv_delta <-  colMeans(H1* (Y-Q[, "QAW"]))
  norm = sqrt(sum(deriv_delta^2))
  # norm
  
  # multidimensional IC over cols of H1.AW (as "x")
  IC.prev <- IC.cur <- apply(t(H1),1,FUN = function(x) x*(Y-Q[, "QAW"]))+
    t(apply(int,1,FUN=function(x) x-psi))
  # head(IC.cur)
  # dim(IC.cur)
  
  # Happens after deriv, a scalar.  (Some each params deriv)
  # if (deriv > 0) { depsilon <- -depsilon}
  loss.prev <- Inf
  loss.cur <-  calcLoss(Q)
  if(is.nan(loss.cur) | is.na(loss.cur) | is.infinite(loss.cur)) {
    loss.cur <- Inf
    loss.prev <- 0
  }
  iter <-  0
  loss.cur
  # max_iter=10
  # psi
  while (loss.prev > loss.cur & iter < max_iter){
    # Have multidimension IC here, written as matrix.
    IC.prev <-  IC.cur
    
    # Do this via vectorized code
    Q.prev <- Q
    g1W <- .bound(g1W, gbounds)
    
    HAW_delta <- vapply(y,FUN = function(y) (B<=y)*(A/g1W-(1-A)/(1-g1W)), FUN.VALUE = rep(1,n))
    H0W_delta = vapply(y,FUN = function(y) -(B<=y)*(1-A)/(1-g1W), FUN.VALUE = rep(1,n))
    H1W_delta = vapply(y,FUN = function(y) (B<=y)*A/g1W, FUN.VALUE = rep(1,n))
 
    deriv_delta = colMeans(IC.cur)
    norm = sqrt(sum(deriv_delta^2))
    
    HAW = apply(HAW_delta,1,FUN=function(x) sum(x*deriv_delta))/norm
    H0W = apply(H0W_delta,1,FUN=function(x) sum(x*deriv_delta))/norm
    H1W = apply(H1W_delta,1,FUN=function(x) sum(x*deriv_delta))/norm
    
    H1 = cbind(HAW=HAW,H0W=H0W,H1W=H1W)
    
    # Need all the ICs for this H1 = t(mean(Dstar))Dstar/normj
    Q  <- .bound(plogis(qlogis(Q.prev) + depsilon * H1), Qbounds) 
    psi.prev <- psi
    B = Q[,"Q1W"]-Q[,"Q0W"]
    
    psi.prev <- psi
    int = vapply(y, FUN = function(y) ifelse(B>y,y+1,B+1),FUN.VALUE = rep(1,n))
    psi = apply(int,2,mean)
    
    # and here, vectorize
    H_delta =  vapply(y,FUN = function(y) (B<=y)*(A/g1W-(1-A)/(1-g1W)),FUN.VALUE=rep(1,n))
    
    # matrix gotten from previous step times Y-Q--not matrix mult
    deriv_delta <-  colMeans(H_delta* (Y-Q[, "QAW"]))
    norm = sqrt(sum(deriv_delta^2))
    
    # multidimensional IC over cols of H1.AW (as "x")
    IC.prev <- IC.cur 
    IC.cur <- apply(t(H_delta),1,FUN = function(x) x*(Y-Q[, "QAW"]))+
      t(apply(int,1,FUN=function(x) x-psi))
    # head(IC.cur)
    # dim(IC.cur)
    
    # Same as ever
    loss.prev <- loss.cur
    loss.cur <- calcLoss(Q)
    loss.cur
    if(is.nan(loss.cur) | is.infinite(loss.cur) | is.na(loss.cur)) {loss.cur <- Inf}
    iter <- iter + 1
    # print(psi.prev)
  }
  return(list(psi = psi.prev, 
              var.psi = apply(IC.prev,2,FUN = function(x) var(x)/n),  
              epsilon = (iter-1)*depsilon, r = cor(IC.prev),IC = IC.prev))
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

y = seq(-.1,.4,by=.01)
y
truth = mclapply(y,FUN = function(y) integral(data.blip$blip,y))


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

save.image(file="primitiveMD.Rdata")
