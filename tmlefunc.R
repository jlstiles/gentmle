# tmle(Y, A, W, Z=NULL, Delta = rep(1,length(Y)), Q = NULL, Q.Z1 = NULL, Qform = NULL,
#      Qbounds = NULL, Q.SL.library = c("SL.glm", "SL.step", "SL.glm.interaction"),
#      cvQinit = FALSE, g1W = NULL, gform = NULL, gbound = 0.025,  pZ1=NULL,
#      g.Zform = NULL, pDelta1 = NULL, g.Deltaform = NULL,
#      g.SL.library = c("SL.glm", "SL.step", "SL.glm.interaction"),
#      family = "gaussian", fluctuation = "logistic", alpha = 0.995, id=1:length(Y),
#      verbose = FALSE)
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

n = 1000
data = gendata(n)
Qkfit= glm(Y~.,data=data,family='binomial')
gkfit = glm(A~W1+W2+W3+W4,data=data,family='binomial')

data1=data
data1$A=1
data0=data
data0$A=0
newdata = rbind(data,data1,data0)

Qk = predict(Qkfit,type='response')
Q0k = predict(Qkfit,newdata = newdata[(2*n+1):(3*n),],type='response')
Q1k = predict(Qkfit,newdata = newdata[(n+1):(2*n),],type='response')
gk=predict(gkfit,type='response')
initdata = data.frame(Qk=Qk,Q1k=Q1k,Q0k=Q0k,gk=gk,A=data$A,Y=data$Y)

tmle = function(Y,A,W,Q.SL.library,g.SL.library,parameter) {
   parameter = "eyATE"
   esti = paste0(parameter,"_estimate")
   upda = paste0(parameter,"_update")
   funas.name(esti)
   fit = gentmle(initdata, estimate_fun=eyATE_estimate, update_fun=eyATE_update, max_iter = 100)
   return(fit)
}
tmle(parameter="eyATE")
head(initdata)

