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
library(BlipVariance)

##
detectCores()
cl = makeCluster(8, type = "SOCK")
registerDoSNOW(cl)

ptm = proc.time()
allresults=foreach(i=1:1000,
                   .packages=c("tmle","plyr","SuperLearner","geepack","BlipVariance"))%dopar%{sim_onestep(1000)}
proc.time()-ptm

#calculate performance
results=do.call(rbind,allresults)
# results
# results = allresults
perf=function(ests,truth){
  n=length(ests)
  var=((n-1)/n)*var(ests)
  bias=mean(ests)-truth
  mse=mean((ests-truth)^2)
  c(var=var,bias=bias,mse=mse)
}

results = as.data.frame(apply(results,2, as.numeric))

t(apply(results[,c(1,5)], 2, perf,var0))
head(results)
mean(results$cover1)
mean(results$cover)

hist(results[,1],breaks=50)
hist(results[,5],breaks=50)

g0=function(W1,W2,W3,W4) {plogis(-.8*W1+.39*W2+.08*W3-.12*W4-.15)}
g0=function(W1,W2,W3,W4) {plogis(-.28*W1+1*W2+.08*W3-.12*W4-1)}
Q0=function(A,W1,W2,W3,W4) {plogis(3*A-.1*W1+1.2*W2-1.9*W3+2*W4)}

gendata=function(n){
  U1 = runif(n,0,1)
  W1= -1*as.numeric(U1<.5)+1*as.numeric(U1>=.5)
  # W1=rnorm(n)
  W2=rnorm(n)
  W3=rbinom(n,1,.3)
  W4=rnorm(n)
  A=rbinom(n,1,g0(W1,W2,W3,W4))
  Y=rbinom(n,1,Q0(A,W1,W2,W3,W4))
  data.frame(A,W1,W2,W3,W4,Y)
}

testdata=gendata(1000000)
v1=testdata$W1
v2=testdata$W2
v3=testdata$W3
v4=testdata$W4
ATE0=mean(Q0(1,v1,v2,v3,v4)-Q0(0,v1,v2,v3,v4))
var0=mean((Q0(1,v1,v2,v3,v4)-Q0(0,v1,v2,v3,v4)-ATE0)^2)
var0
ATE0

sim_onestep = function(n){
  n = 1000
  data <- gendata(n)
  data1=data
  data1$A=1
  data0=data
  data0$A=0
  newdata = rbind(data,data1,data0)
  X = data
  X$Y=NULL
  newdata$Y=NULL
  SL.library = "SL.glm"
  Qfit=SuperLearner(data$Y,X,newX=newdata, family = binomial(),
                    SL.library=SL.library,  method = "method.NNLS",
                    id = NULL, verbose = FALSE, control = list(),
                    cvControl = list(V=10), obsWeights = NULL)
  gfit=glm(A~W1+W2+W3+W4,data=data,family=binomial(link="logit"))
  gk=predict(gfit,type="response")
  hist(gk,breaks=50)
  Qk = Qfit$SL.predict[1:n]
  Q1k = Qfit$SL.predict[(n+1):(2*n)]
  Q0k = Qfit$SL.predict[(2*n+1):(3*n)]

  initdata = data.frame(Qk=Qk,Q1k=Q1k,Q0k=Q0k,gk=gk,A=data$A,Y=data$Y)
  sigmaATE.info = gentmle(initdata,sigmaATE_estimate,sigmaATE_update,max_iter=100)
  ci_regular = ci_gentmle(sigmaATE.info)[2:5]
  Q = cbind(QAW=Qk,Q0W=Q0k,Q1W=Q1k)
  depsilon = .001
  gbounds <- Qbounds <- c(10^-9, 1-10^-9)

  res.onestep <- oneStepblipvar(Y = data$Y, A = data$A, Q = Q, g1W = gk,
                            depsilon = depsilon, max_iter = max(1000, 2/depsilon),
                            gbounds = gbounds, Qbounds = Qbounds)

  sd = sqrt(res.onestep$var.psi)
  est = res.onestep$psi
  ci_onestep = c(est = est, sd = sd,lower = est - 1.96*sd, upper = est + 1.96*sd)
  cover = as.numeric((ci_regular[3]<=var0)*(ci_regular[4]>=var0))
  cover1 = as.numeric((ci_onestep[3]<=var0)*(ci_onestep[4]>=var0))
  ci_regular
  ci_onestep
  return(c(ci_regular,ci_onestep,cover=cover, cover1=cover1))
}


#-----------------------------------One-Step TMLE for ATT parameter  ----------------------------------------
oneStepblipvar <- function(Y, A, Q, g1W, depsilon, max_iter, gbounds, Qbounds){
	n <- length(Y)
	calcLoss <- function(Q) -mean(Y * log(Q[,"QAW"]) + (1-Y) * log(1 - Q[,"QAW"]))
	psi.prev <- psi  <- var((Q[,"Q1W"] - Q[, "Q0W"]))
	H1.AW =  2*(Q[,"Q1W"]-Q[,"Q0W"]-mean(Q[,"Q1W"]-Q[,"Q0W"]))*(A/g1W-(1-A)/(1-g1W))
	IC.prev <- IC.cur <- H1.AW*(Y-Q[, "QAW"])+
	  (Q[,"Q1W"]-Q[,"Q0W"]-mean(Q[,"Q1W"]-Q[,"Q0W"]))^2-psi

	deriv <-  mean(H1.AW* (Y-Q[, "QAW"]))

	if (deriv > 0) { depsilon <- -depsilon}
	loss.prev <- Inf
 	loss.cur <-  calcLoss(Q)
 	if(is.nan(loss.cur) | is.na(loss.cur) | is.infinite(loss.cur)) {
 		loss.cur <- Inf
 		loss.prev <- 0
 	}
 	iter <-  0
 	while (loss.prev > loss.cur & iter < max_iter){
		IC.prev <-  IC.cur
		Q.prev <- Q
		g1W <- .bound(g1W, gbounds)
 		H1 <- cbind(HAW = 2*(Q.prev[,"Q1W"]-Q.prev[,"Q0W"]-mean(Q.prev[,"Q1W"]-Q.prev[,"Q0W"]))*(A/g1W-(1-A)/(1-g1W)),
					  H0W =   - 2*(Q.prev[,"Q1W"]-Q.prev[,"Q0W"]-mean(Q.prev[,"Q1W"]-Q.prev[,"Q0W"]))*(1-A)/(1-g1W),
					  H1W =    2*(Q.prev[,"Q1W"]-Q.prev[,"Q0W"]-mean(Q.prev[,"Q1W"]-Q.prev[,"Q0W"]))*A/g1W)

 		Q  <- .bound(plogis(qlogis(Q.prev) - depsilon * H1), Qbounds)

 		psi.prev <- psi
 		psi <- var((Q[,"Q1W"] - Q[, "Q0W"]))
 		loss.prev <- loss.cur
 		loss.cur <- calcLoss(Q)
 		IC.cur <- 2*(Q[,"Q1W"]-Q[,"Q0W"]-mean(Q[,"Q1W"]-Q[,"Q0W"]))*(A/g1W-(1-A)/(1-g1W))*(Y-Q[, "QAW"])+
 		  (Q[,"Q1W"]-Q[,"Q0W"]-mean(Q[,"Q1W"]-Q[,"Q0W"]))^2-psi
 		if(is.nan(loss.cur) | is.infinite(loss.cur) | is.na(loss.cur)) {loss.cur <- Inf}
 		iter <- iter + 1
 		# print(psi.prev)
 	}
 	 	return(list(psi = psi.prev, var.psi = var(IC.prev)/n,  epsilon = (iter-1)*depsilon))
 }

#------------- bound function --------------
.bound <- function(x, bds){
	x[x > max(bds)] <- max(bds)
	x[x < min(bds)] <- min(bds)
	x
}

