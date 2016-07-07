library(Hmisc)
library(reshape2)
library(plyr)
library(dplyr)
library(sas7bdat)
library(memoise)
library(ltmle)
library(foreach)
library(doSNOW)
library(parallel)
library(mnlogit)

set.seed(252)
n = 1000000
W0 = rbinom(n,1,0.5)
W1 = runif(n,0,1)

p1 = exp(-0.6+0.4*W0-0.3*W1)
p2 = exp(W0-W1)
p3 = exp(-.4*W0+W1)
p4 = exp(.3-W0+.4*W1)
p=p1+p2+p3+p4

# These are the probs of getting each level
p1=p1/p
p2=p2/p
p3=p3/p
p4=p4/p

# Initializing a matrix so each row can be filled with a treatment for each subject
Arows = matrix(rep(NA,(4*n)),ncol=4)

# Fills in the treatment for each subject
for (a in 1:n){
  Arows[a,]=rmultinom(1,1,c(p1[a],p2[a],p3[a],p4[a]))
}

df = data.frame(W0,W1,Arows)
head(df)
#
# df$id = 1:1000
# df.melt = melt(df,id=c("id","W0","W1"))
#
# df.melt = df.melt[order(df.melt$id),]
# df.melt[1:10,]


f1 = glm(Arows[,1]~W0+W1,family='binomial')

ind2 = which(Arows[,1]==0)
ind2.not = which(Arows[,1]==1)
df2 = data.frame(W0,W1,A = Arows[,2])
df2 = df2[ind2,]
f2 = glm(A~W0+W1,data=df2,family='binomial')

ind3 = which(Arows[,1]==0&Arows[,2]==0)
ind3.not = which(!(Arows[,1]==0&Arows[,2]==0))
df3 = data.frame(W0,W1,A = Arows[,3])
df3 = df3[ind3,]
f3 = glm(A~W0+W1,data=df3,family = 'binomial')

newdata = data.frame(W0=W0, W1=W1)
cat4= (1-predict(f3,type='response',newdata=newdata))*
  (1-predict(f2,type='response',newdata=newdata))*
  (1-predict(f1,type='response',newdata=newdata))
cat3 = predict(f3,type='response',newdata=newdata)*
  (1-predict(f2,type='response',newdata=newdata))*
  (1-predict(f1,type='response',newdata=newdata))
cat2 = predict(f2,type='response',newdata=newdata)*
  (1-predict(f1,type='response',newdata=newdata))
cat1 = predict(f1,type='response',newdata=newdata)

preds = data.frame(cat1,cat2,cat3,cat4)

fmn = multinom(Arows~W0+W1)
pred <- predict(fmn, newdata = data.frame(W0=W0,W1=W1), type = "prob")

preds[1:20,]-
pred[1:20,]

Qbar0 <- function(A, W) {

  W1 <- W[, 1]
  W2 <- W[, 2]
  W3 <- W[, 3]
  W4 <- W[, 4]
  Qbar <- plogis(ifelse(W4 > 0, (A == 1) + (A == 1) * (5 * W1^2 - 4.45), (A == 2) + (A == 3) +
                          (A == 2) * (4 * W2) + (A == 3) * (5 * W3)))
  return(Qbar)
}

g0 <- function(W) {
  W1 <- W[, 1]
  W2 <- W[, 2]
  W3 <- W[, 3]
  W4 <- W[, 4]

  # rep(0.5, nrow(W))
  A1 <- plogis(W1)
  A2 <- plogis(W2)
  A3 <- plogis(W3)
  A <- cbind(A1, A2, A3)

  # make sure A sums to 1
  A <- normalize_rows(A)
}

gen_data <- function(n = 1000, p = 4) {
  # n=1000
  # p=4
  W <- matrix(rnorm(n * p), nrow = n)
  colnames(W) <- paste("W", seq_len(p), sep = "")
  pA <- g0(W)
  # head(pA)
  A <- factor(apply(pA, 1, function(pAi) which(rmultinom(1, 1, pAi) == 1)))
  A_vals <- vals_from_factor(A)
  u <- runif(n)
  Y <- as.numeric(u < Qbar0(A, W))
  Q0aW <- sapply(A_vals, Qbar0, W)
  # head(Q0aW)
  d0 <- apply(Q0aW, 1, which.max)
  d0[1:10]
  Yd0 <- as.numeric(u < Qbar0(d0, W))
  data.frame(W, A, Y, Q0aW, d0, Yd0)
}

data <- gen_data(1e+05, 5)

Anode <- "A"
Wnodes <- grep("^W", names(data), value = T)

Q_fit <- glm(data$Y ~ ., data[, c("A", Wnodes)], family = binomial(link = "logit"))
g_fit <- multinomial_SuperLearner(data$A, data[, Wnodes])

A_vals <- vals_from_factor(data$A)

Q_a <- sapply(A_vals, function(A_val) {
  newdata <- data[, c(Anode, Wnodes)]
  newdata[, Anode] <- A_val
  predict(Q_fit, newdata, type = "response")
})

pA <- predict(g_fit, newdata = data[, Wnodes])$pred

# A sample gstar--treat with A=1, if the patient received A=1 or 2, otherwise leave alone
gstar <- function(A, gk) {
  ifelse(A == 1, gk[, 1] + gk[, 2], ifelse(A == 3, gk[, 3], 0))
}
