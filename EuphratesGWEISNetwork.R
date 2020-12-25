#-----FunMap-----
library(mvtnorm)
salttable <- read.csv("./path/saltstable.csv",header = T,row.names = 1)
cktable <- read.csv("./path/ckstable.csv",header = T,row.names = 1)
genotable <- read.csv("./path/genostable.csv",header = T,row.names = 1)

get_miu = function(miu_par,t){
  miu = miu_par[1]/(1 + miu_par[2] * exp(-miu_par[3] * t)) - (miu_par[4] * exp(-miu_par[5] * t))
  miu
}

get_SAD1 <- function(x,t){
  n=length(t)
  tmp1 <- (1-x[1]^2) 
  tmp2 <- (1-x[3]^2) 
  sig1 <- array(1, dim=c(n,n))   
  for(i in 1:n)   
  {     
    sig1[i,i:n] <- x[1]^( c(i:n) - i ) * (1-x[1]^(2*i))/tmp1 
    sig1[i:n,i] <- sig1[i,i:n]   
  } 
  sig1 <- sig1*x[2]^2
  sig2 <- array(1, dim=c(n,n))   
  for(i in 1:n){     
    sig2[i,i:n] <- x[3]^( c(i:n) - i ) * (1-x[3]^(2*i))/tmp2    
    sig2[i:n,i] <- sig2[i,i:n]   
  }  
  sig2 <- sig2*x[4]^2
  sig12 <- array(0, dim=c(n,n))
  
  sigma1 <- cbind(sig1,sig12)
  sigma2 <- cbind(sig12,sig2)
  sigma <- rbind(sigma1,sigma2)
  return(sigma)
}

merge_miu <- function(x,t){
  merge_miu <- c(get_miu(x[1:5],t),get_miu(x[6:10],t))
  merge_miu
}

L0 = function(par,t,y){
  miu = merge_miu(par[1:10],t)
  SAD1 = get_SAD1(par[11:14],t)
  L0 = -sum(dmvnorm(y,miu,SAD1,log = T))
  L0
}
L1 = function(par,t,marker){
  ckgeno1 <- cktable[which(marker==0),]
  saltgeno1 <- salttable[which(marker==0),]
  ckgeno2 <- cktable[which(marker==1),]
  saltgeno2 <- salttable[which(marker==1),]
  ckgeno3 <- cktable[which(marker==2),]
  saltgeno3 <- salttable[which(marker==2),]
  SAD1 = get_SAD1(par[31:34],t)
  miu1 = merge_miu(par[1:10],t)
  miu2 = merge_miu(par[11:20],t)
  miu3 = merge_miu(par[21:30],t)
  l1_1 <- sum(dmvnorm(cbind(ckgeno1,saltgeno1),miu1,SAD1,log = T))
  l1_2 <- sum(dmvnorm(cbind(ckgeno2,saltgeno2),miu2,SAD1,log = T))
  l1_3 <- sum(dmvnorm(cbind(ckgeno3,saltgeno3),miu3,SAD1,log = T))
  L1 <- -(l1_1 + l1_2 + l1_3)
  L1
}
L2 = function(par,t,marker){
  ckgeno1 <- cktable[which(marker==0),]
  saltgeno1 <- salttable[which(marker==0),]
  ckgeno2 <- cktable[which(marker==1),]
  saltgeno2 <- salttable[which(marker==1),]
  SAD1 = get_SAD1(par[21:24],t)
  miu1 = merge_miu(par[1:10],t)
  miu2 = merge_miu(par[11:20],t)
  l1_1 <- sum(dmvnorm(cbind(ckgeno1,saltgeno1),miu1,SAD1,log = T))
  l1_2 <- sum(dmvnorm(cbind(ckgeno2,saltgeno2),miu2,SAD1,log = T))
  L2 <- -(l1_1 + l1_2)
  L2
}
t = seq(13,78,5)

LR = function(marker){
  t = t
  ckgeno1 <- cktable[which(marker==0),]
  saltgeno1 <- salttable[which(marker==0),]
  ckgeno2 <- cktable[which(marker==1),]
  saltgeno2 <- salttable[which(marker==1),]
  ckgeno3 <- cktable[which(marker==2),]
  saltgeno3 <- salttable[which(marker==2),]
  if(length(which(marker==9))==0){
    y_all = cbind(cktable,salttable)
  }else{
    y_all <- cbind(cktable[-which(marker == 9),],salttable[-which(marker == 9),])
  }
  if(2 %in% marker){
    NH0 = optim(c(71,21,0.062,7.8,0.0165,71.8,16.3,0.0336,12,0.05,0.001,50,0.001,50),L0,t=t,y=y_all,method = "BFGS",control=list(maxit=50000))
    pars <- c(NH0$par[1:10],NH0$par[1:10],NH0$par[1:10],NH0$par[11:12],NH0$par[11:12])
    NH1 <- optim(pars,L1,t=t,marker=marker,method="BFGS",control=list(maxit=50000))
    LR<- 2*(NH0$value - NH1$value)
    cat(i,NH0$value,NH1$value,LR,"\n")
    LR
  }else{
    NH0 = optim(c(71,21,0.062,7.8,0.0165,71.8,16.3,0.0336,12,0.05,0.001,50,0.001,50),L0,t=t,y=y_all,method = "BFGS",control=list(maxit=50000))
    pars <- c(NH0$par[1:10],NH0$par[1:10],NH0$par[11:12],NH0$par[11:12])
    NH1 <- optim(pars,L2,t=t,marker=marker,method="BFGS",control=list(maxit=50000))
    LR<- 2*(NH0$value - NH1$value)
    cat(i,NH0$value,NH1$value,LR,"\n")
    LR
  }
}

lr = rep(0,(dim(genotable)[1]))
for (i in 1:dim(genotable)[1]) {
  lr[i] = LR(as.numeric(genotable[i,]))
}
save(lr,file = "filename.Rdata")
#-----FunCluster-----
library(MASS)
setwd("./path")

SsquareDiff <<- 1.0e-5
IncLimit <<- 3
REPEAT_LIMIT <<- 200

LGD_J <<- 5
START_J <<- 8
END_J <<- 15

LIKELIHOOD_DIFF <<- 0.5
rhoIncreasement <<- 0.002
rhoStart <<-0.8

datafile.L <<- "ck_gene_effect.csv.csv"
datafile.H <<- "salt_gene_effect.csv.csv"

set.seed(Sys.time())

# Legendre Polynominals
LgdP <- expression( tt,
                    ( 3* tt^2 - 1 )/2 , 
                    ( 5 * tt^3 - 3* tt )/2, 
                    ( 35 * tt^4 - 30 * tt^2 + 3)/8,
                    ( 63 * tt^5 - 70 * tt^3 + 15 * tt )/8,
                    ( 231 * tt^6 - 315 * tt^4 + 105 * tt^2 - 5)/16  )

GetMR <- function(rho,times)
{
  MR <- matrix(1,length(times),length(times))
  for ( i in 1:length(times)){
    for(j in 1:length(times)){
      MR[i,j]= rho^(abs(times[j] - times[i]))
    }
  }
  return (MR)
}
GetMX <- function(times,r)
{
  tnum = length(times) 
  X <- matrix(1,tnum,r+1)
  
  for(t in 1:tnum ){
    tt <- -1 + 2*(times[t] - times[1])/(times[tnum] - times[1])
    for(i in 1:r){
      X[t,i+1] <- eval(LgdP[i])
    }
  }
  return (X)
}
GetInitPij <- function(N,J)
{
  P <- matrix(1/J,N,J)
  for (i in 1:N){
    P[i,] <- rnorm(J, mean=1/J, sd= 0.5 * 1/J )
    P[i,] <- P[i,]/sum(P[i,])
  }
  
  return (P)
}
GetMeanMatrix <- function(J,times,P,X,Asdata,InvMSigema)
{
  m <- matrix(NA,J,length(times))
  N <- length(Asdata[,1])
  r <- length(X[1,])
  
  xInvSigema <- t(X) %*% InvMSigema
  xInvSigemax <- xInvSigema%*% X
  
  for( j in 1:J){
    ud <- matrix(0, r, r)
    for( i in 1: N){
      ud <- ud + P[i,j]*xInvSigemax
    }
    ubd <- matrix(0, r, 1)
    for( i in 1: N){
      ubd <- ubd + P[i,j]*( xInvSigema %*% (Asdata[i,]) )
    }
    uj <- ginv(ud) %*% ubd
    m[j,] <- X %*% uj
  }
  return(m)
}
GetNewSsquare <- function(Asdata,m,MR,times,P,J)
{
  N <- length(Asdata[,1])
  
  InvMR <- ginv(MR)
  newSsquare <- 0
  for(i in 1:N){
    SumJ <- 0
    for(j in 1:J){
      yi_mj <- Asdata[i,]-m[j,]
      SumJ <- SumJ + P[i,j] * ((yi_mj) %*% InvMR %*% (yi_mj) )
    }
    newSsquare <- newSsquare + SumJ
  }
  newSsquare <- as.numeric(newSsquare/(length(times)*N))
  
  return(newSsquare)
}
GetNewRho.b <- function(rho,rhoDir)
{
  
  newrho <- as.numeric(rho + rhoIncreasement*rhoDir)
  if (newrho > 1) newrho <- 1
  if (newrho < 0) newrho <- 0
  
  return (newrho)
}
GetNewRho<- function(Asdata,m,MR,times,P,J,rho,Ssquare)
{
  N <- length(Asdata[,1])
  newrho <- 0
  for(i in 1:N){
    SumJ <- 0
    for(j in 1:J){
      yi_mj <- Asdata[i,]-m[j,]
      Item1 <- (1/(1 - rho*rho))*((yi_mj) %*% MR %*% (yi_mj) )
      Item2 <- 0
      for(k in 2:(length(times)-1) )
        Item2 <- Item2 + (yi_mj[k]^2)
      Item2 <- Item2 * rho
      Item3 <- 0
      for(k in 1:(length(times)-1) )
        Item2 <- Item3 + yi_mj[k] * yi_mj[k+1]
      SumJ <- SumJ + P[i,j] * (Item1 + Item2 - Item3)
    }
    newrho <- newrho + SumJ
  }
  newrho <- as.numeric(newrho/( (length(times)-1)* N * Ssquare))
  
  if(abs(newrho) >= 1) return( sign(newrho)*.5)
  else return(newrho)
}
GetLikelihood <- function(Asdata.H, m.H, MSigema.H, Asdata.L, m.L, MSigema.L, omiga,P,J,times)
{
  N <- length(Asdata.H[,1])
  
  InvMSigema.H <- ginv(MSigema.H)
  InvMSigema.L <- ginv(MSigema.L)
  
  DetMSigema.H <- det(MSigema.H)
  DetMSigema.L <- det(MSigema.L)
  
  LogDetMSigema.H <- log(DetMSigema.H)/2
  LogDetMSigema.L <- log(DetMSigema.L)/2
  
  LogM2Pi <- length(times)*log(2*pi)
  
  oneterm <- function(i, j) {
    f <- function(i,j)
      P[i,j]*(log(omiga[j]) - LogM2Pi - LogDetMSigema.H - LogDetMSigema.L
              - ( ((Asdata.H[i,]-m.H[j,])) %*% InvMSigema.H %*% (Asdata.H[i,]-m.H[j,])) /2
              - ( ((Asdata.L[i,]-m.L[j,])) %*% InvMSigema.L %*% (Asdata.L[i,]-m.L[j,])) /2)
    mapply(f, i, j)
  }
  tmp <- outer(1:N, 1:J, oneterm)  
  tmp[!is.finite(tmp)] <- min(tmp[is.finite(tmp)])
  return(sum(tmp))
}
StepE <- function(Asdata.H, m.H, MSigema.H, Asdata.L, m.L, MSigema.L, omiga,P,J,times)
{
  InvMSigema.H <- ginv(MSigema.H)
  InvMSigema.L <- ginv(MSigema.L)
  
  N <- length(Asdata.H[,1])
  
  for( i in 1:N){
    Fi <- rep(0,J)
    for( j in 1:J){
      yi_mj.H <- Asdata.H[i,]-m.H[j,]
      Fi.H = exp( ( (yi_mj.H) %*% InvMSigema.H %*% (yi_mj.H) ) / -2)
      
      yi_mj.L <- Asdata.L[i,]-m.L[j,]
      Fi.L = exp( ( (yi_mj.L) %*% InvMSigema.L %*% (yi_mj.L) ) / -2)
      
      Fi[j] = Fi.H * Fi.L
    }
    OmigaF <- omiga %*% Fi
    P[i,] <- (omiga * Fi) / c(OmigaF)
  }
  if (all(is.nan(P))!=TRUE) {
    for (q in 1:length(P[1,])) {
      sit <- which(is.nan(P[,q]) == TRUE)
      P[sit,q] <- P[which.min(P[,q]),q]
    }
  }
  return(P)
}

StepM <- function(Asdata,m,MR,times,Ssquare,P,rho,rhoDir,J,rpt)
{
  newSsquare <- GetNewSsquare(Asdata,m,MR,times,P,J)
  if (rpt > 0)
    newrho <- GetNewRho(Asdata,m,MR,times,P,J,rho,Ssquare)
  else
    newrho <- rho
  
  return( c(newSsquare, newrho))
}
StepM.b <- function(Asdata.H,Asdata.L,m.H,m.L,MR,times,P,rho,rhoDir,J,rpt)
{
  newSsquare.H <- GetNewSsquare(Asdata.H,m.H,MR,times,P,J)
  newSsquare.L <- GetNewSsquare(Asdata.L,m.L,MR,times,P,J)
  if (rpt > 0){
    newrho <- GetNewRho.b(rho,rhoDir)
  }else{
    newrho <- rho
  }
  return( c(newSsquare.H, newSsquare.L, newrho))
}
RunEM.Joint <- function(Asdata.H,Asdata.L,times,X,P,MR,rho,MSigema.H,MSigema.L,omiga,m.H,m.L,J,r){
  rpt <- 1
  Likelihood <- -Inf
  
  rhoDir <- 1
  rhoIncCount <- 0
  
  while(TRUE){
    OldLikelihood <- Likelihood
    
    P <- StepE(Asdata.H, m.H, MSigema.H, Asdata.L, m.L, MSigema.L, omiga,P,J,times)
    
    newpars <- StepM.b(Asdata.H, Asdata.L, m.H, m.L, MR, times,P,rho,rhoDir,J,rpt)
    Ssquare.H <- newpars[1]
    Ssquare.L <- newpars[2]
    rho <- newpars[3]
    
    MR <- GetMR(rho,times)
    
    MSigema.H <- Ssquare.H * MR
    InvMSigema.H <- ginv(MSigema.H)
    
    MSigema.L <- Ssquare.L * MR
    InvMSigema.L <- ginv(MSigema.L)
    
    N <- length(Asdata.H[,1])
    omiga <- colSums(P)/ N 
    
    m.H <- GetMeanMatrix(J,times,P,X,Asdata.H,InvMSigema.H)
    m.L <- GetMeanMatrix(J,times,P,X,Asdata.L,InvMSigema.L)
    
    Likelihood <- GetLikelihood(Asdata.H, m.H, MSigema.H, Asdata.L, m.L, MSigema.L, omiga,P,J,times)
    
    if ( Likelihood >= OldLikelihood){
      rhoIncCount <- 0
    }else{
      rhoIncCount <- rhoIncCount + 1
      if (rhoIncCount >= IncLimit){
        rhoIncCount <- 0
        rhoDir <- rhoDir * -1
      }
    }
    
    cat("J:",J,"	rpt:", rpt, "\n")
    cat("Ssquare.H:", Ssquare.H, "    Ssquare.L:", Ssquare.L, "\n")
    cat("rho:", rho, "\n")
    cat("omiga:",omiga,"\n")
    cat("Likelihood:",Likelihood,"\n\n")

    if( (abs(abs(OldLikelihood - Likelihood) - Likelihood) < LIKELIHOOD_DIFF) ) { 
      cat("quit due to likelihood\n")
      cat("LIKELIHOOD_DIFF:",LIKELIHOOD_DIFF,"\n")
      break
    }
    if( rpt >= REPEAT_LIMIT ){
      cat("quit due to rpt\n")
      break
    }
    rpt <- rpt + 1
  }
  OpiPfileName <- sprintf("OpiP%02d.LGD%d.csv",J,r)
  write.csv(P,OpiPfileName,row.names = FALSE)
  
  OpiMfileName <- sprintf("OpiM%02d.LGD%d.H.csv",J,r)
  write.csv(m.H,OpiMfileName,row.names = FALSE)
  
  OpiMfileName <- sprintf("OpiM%02d.LGD%d.L.csv",J,r)
  write.csv(m.L,OpiMfileName,row.names = FALSE)
  
  return(c(rho,Ssquare.H,Ssquare.L,Likelihood))
}

InitAndRunEM.Joint <- function(J,r) 
{
  cat("r:",r,"\n")
  Asdata.H <- read.csv(datafile.H,row.names = 1)
  Asdata.H <- Asdata.H[,-1]
  Asdata.H <- as.matrix(Asdata.H)
  colnames(Asdata.H) <- NULL
  
  Asdata.L <- read.csv(datafile.L,row.names = 1)
  Asdata.L <- Asdata.L[,-1]
  Asdata.L <- as.matrix(Asdata.L)
  colnames(Asdata.L) <- NULL
  
  times <- seq(13,78,5)
  
  rho <- rhoStart
  Ssquare <- 20
  
  X <- GetMX(times,r)
  N <- length(Asdata.H[,1])
  
  P <- GetInitPij(N,J)
  
  MR <- GetMR(rho,times)
  MSigema <- Ssquare * MR
  InvMSigema <- ginv(MSigema)
  DetMSigema <- (det(MSigema))^0.5
  
  omiga <- colSums(P)/ N
  
  m.H <- GetMeanMatrix(J,times,P,X,Asdata.H,InvMSigema)
  m.L <- GetMeanMatrix(J,times,P,X,Asdata.L,InvMSigema)
  
  EMResults <- RunEM.Joint(Asdata.H, Asdata.L,times,X,P,MR,rho,MSigema,MSigema,omiga,m.H,m.L,J,r)
  
  return(EMResults)
  
}
EMResults <- matrix(0, END_J - START_J + 1, 5)
for(j in START_J:END_J){
  EMResults[j - START_J + 1,1] <- j
  EMResults[j - START_J + 1,2:5] <- InitAndRunEM.Joint(J=j,r=LGD_J)
}
resultFileName <- sprintf("EMResults%02d~%02d.LGD%d.csv",START_J,END_J,LGD_J)
write.csv(EMResults,resultFileName,row.names = FALSE)
#-----NetRestructure-----
library(glmnet)

data_effect <- read.csv("./path/ck_gene_effect.csv")
score <- matrix(NA,length(data_effect[,1]),1)
clusterscore <- read.csv(OpiPfileName,row.names = 1)
for (i in 1:length(clusterscore[,1])) {
  score[i] <- which.max(clusterscore[i,])
}
mod = 1

means <- matrix(NA,length(which(score == mod)),14)
for (i in 1:length(table(score))) {
  sit <- which(score == i)
  means[i,] = apply(data_effect[sit,],2,mean)
}
means <- t(means)
colnames(means) <- c(1:length(table(score)))
name <- c(1:length(means[1,]))
marker_list <- list()

for (col in mod) {
  ridge1_cv <- cv.glmnet(x = means[,-col], y = means[,col],type.measure = "mse",nfold = 10,alpha = 0,grouped=FALSE)
  best_ridge_coef <- as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1]
  fit_res <- cv.glmnet(x = means[,-col], y = means[,col],type.measure = "mse",nfold = 10,alpha = 1,penalty.factor = 1 / abs(best_ridge_coef),keep = TRUE,grouped=FALSE)
  best_alasso_coef1 <- coef(fit_res, s = fit_res$lambda.min)
  marker_list_one <- list()
  marker_list_one[[1]] <- name[col]#第一个列表是直接qtl的名字
  marker_list_one[[2]] <- as.numeric(best_alasso_coef1@Dimnames[[1]][best_alasso_coef1@i[-1]+1])#第二个列表是间接qtl的名字
  marker_list_one[[3]] <- best_alasso_coef1@x[-1]#第三个列表是变量选择系数
  marker_list[[col]] <- marker_list_one
  #proc.time() - tim
}

load("./path/Effect.Rdata")

get_LOPm <- function(X){
  len = length(X)
  LOP <- function(r){
    t <- seq(-1,1,2/(len-1))
    temp <- rep(0,len)
    for (m in 0:as.integer(r/2)) {
      temp <- temp  + (-1)^m*gamma(2*r - 2*m + 1)/(2^r*gamma(m+1)*gamma(r-m+1)*gamma(r-2*m + 1)) * t^(r-2*m)
    }
    return(temp)
  }
  LOPm <- cbind(LOP(0),LOP(1),LOP(2),LOP(3),LOP(4),LOP(5),LOP(6) )
  return(LOPm[,1:ORDER])
}
ORDER <- 6
library(mvtnorm)
f1 <- function(x,t){
  y = t[1]/(1 + t[2] * exp(-t[3] * x)) - (t[4] * exp(-t[5] * x))
  return(y)
}
fy <- function(t,X){
  if(t[1] == 2){
    e1 = f1(X,t[2:6])
    e2 = f1(X,t[8:12])
    additive = 0.5 * (e1 - e2)
    dominant = 0
    y <- (2 * t[7] * t[13] * (additive + (t[7] - t[13]) * dominant)^2 + 4 * t[7]^2 * t[13]^2 * dominant^2) ^ 0.5
  }else{
    e1 = f1(X,t[2:6])
    e2 = f1(X,t[8:12])
    e3 = f1(X,t[14:18])
    additive = 0.5 * (e1 - e3)
    dominant = e2 - 0.5 * (e1 + e3)
    y <- (2 * t[7] * t[13] * (additive + (t[7] - t[13]) * dominant)^2 + 4 * t[7]^2 * t[13]^2 * dominant^2) ^ 0.5
  }
}
get_origin <- function(dy,X,y0){
  y0 <- c(y0)
  for (i in 2:(length(X)-1)) {
    slope <- dy[i-1]
    y_before <- y0[length(y0)]
    add <- y_before + slope*(X[2]-X[1])
    y0 <- c(y0,add)
  }
  return(y0)
}
fl_new <- function(t,X,dep,ind,dep_per,ind_per,LOPm){
  ydep <- matrix(NA,length(dep_per[,1]),length(X))
  for (i in 1:length(dep_per[,1])) {
    ydep[i,] <- fy(c(dep_per[i,]),X) 
  }
  ydep <- apply(ydep, 2, mean)
  d <- 2*(ydep[-1] - ydep[-length(ydep)])/(X[2]-X[1])
  tm <- matrix(t,ncol=ORDER,byrow = T)
  temp1 <- LOPm[-1,]%*%t(tm) # mp * m个线
  temp1 <- temp1*matrix(rep(ydep[-1],length(ind)+1),ncol = length(ind)+1,byrow = F)
  num = 0 
  fy1 <- matrix(NA,length(ind),length(X))
  for (i in 1:length(ind)) {
    sit = which(score == ind[i])
    leng = length(sit)
    raw = matrix(NA,leng,length(X))
    if (i == 1) {
      parameter = ind_per[1:leng,]
      for (o in 1:leng) {
        raw[o,] = fy(parameter[o,] ,X)
      }
      fy1[i,] <- apply(raw,2,mean)
      num = num + leng
    }else{
      parameter = ind_per[(num+1):(num +leng),]
      for (o in 1:leng) {
        raw[o,] = fy(parameter[o,] ,X)
      }
      fy1[i,] <- apply(raw,2,mean)
      num = num + leng
    }
  }
  for (i in 1:length(ind)) {
    yind <- fy1[i,]
    temp1[,i+1] <- temp1[,i+1]*yind[-1]
  }
  x0 <- X[-length(X)]
  temp0 <- LOPm[-length(X),]%*%t(tm)
  temp0 <- temp0*matrix(rep(ydep[-length(X)],length(ind)+1),ncol = length(ind)+1,byrow = F)
  for (i in 1:length(ind)) {
    yind <- fy1[i,]
    temp0[,i+1] <- temp0[,i+1]*yind[-length(X)]
  }
  #----------------
  d_mat <- LOPm%*%t(tm)
  for (i in 1:length(tm[,1])) {
    if (i == 1) {
      d_mat[,i] <- d_mat[,i]*ydep
    }else{
      d_mat[,i] <- d_mat[,i]*fy1[i-1,]
    }
  }
  # 将d_mat转化为o_mat
  o_mat <- c()
  for (i in 1:length(d_mat[1,]) ) {
    if(i == 1){
      o_mat <- cbind(o_mat,get_origin(d_mat[,i],X,ydep[1]))
    }else{
      o_mat <- cbind(o_mat,get_origin(d_mat[,i],X,0))
    }
  }
  o_mat <- as.data.frame(o_mat)
  y <- ydep[-length(ydep)]
  e <- colSums(t(o_mat))
  ssr <- sum((y-e)^2)
  sst <- sum((y-mean(y))^2)
  r <- 1-(ssr/sst)
  #----------------
  return(sum((d - colSums(t(temp1 + temp0)) )^2) + abs(1-r))
}

for (col in length(marker_list) ) {
  cat(col,'\n')
  dep <- marker_list[[col]][[1]]
  dep <- which(score == dep)
  
  ind <- marker_list[[col]][[2]]
  sit = list()
  for (i in 1:length(ind)) {
    sit[[i]] = which(score == ind[i])
  }
  sit = unlist(sit)
  
  if( length(sit) == 0 ){
    next
  }
  dep_per <- matrix(NA, ncol = 18,nrow = length(dep))
  for (i in 1:length(dep)) {
    dep_per[i,1:length(picture[[dep[i]]])] <- picture[[dep[i]]]
  }
  
  ind_per <- matrix(NA, ncol = 18,nrow = length(sit))
  for (i in 1:length(sit)) {
    ind_per[i,1:length(picture[[sit[i]]])] <- picture[[sit[i]]]
  }
  
  X <- seq(0,78,78/((length(ind)+1)*ORDER*4))
  
  # 等式的右边有两部分,1~n/0~n-1
  t0 <- rep(0.001,(length(ind)+1)*ORDER)
  
  itimes <- 1
  repeat{
    s1 <- optim(t0,fl_new,method = 'Nelder-Mead',X = X,dep = dep,ind = ind, 
                dep_per = dep_per, ind_per = ind_per, LOPm = get_LOPm(X))
    r1 <- s1$par
    s2 <- optim(r1,fl_new,method = 'Nelder-Mead',X = X,dep = dep,ind = ind, 
                dep_per = dep_per, ind_per = ind_per, LOPm = get_LOPm(X))
    cat(col,'-',itimes,s2$value,'\n')
    itimes <- itimes + 1
    if(all( abs(r1-s2$par) == 0 )||itimes == 10){ #*** itimes越高精度越高,计算速度越慢,有条件部署在集群时,应该尽可能大与1000 ***#
      break
    }else{
      t0 <- s2$par
    }
  }
  
  marker_list[[col]][[4]] <- matrix(s2$par,ncol=ORDER,byrow=TRUE)
  tm <-  matrix(s2$par,ncol=ORDER,byrow=TRUE)
  
  ydep <- matrix(NA,length(dep_per[,1]),length(X))
  for (i in 1:length(dep_per[,1])) {
    ydep[i,] <- fy(c(dep_per[i,]),X) 
  }
  ydep <- apply(ydep, 2, mean)
  d_mat <- get_LOPm(X)%*%t(tm)
  num = 0 
  fy1 <- matrix(NA,length(ind),length(X))
  for (i in 1:length(ind)) {
    sit = which(score == ind[i])
    leng = length(sit)
    raw = matrix(NA,leng,length(X))
    if (i == 1) {
      parameter = ind_per[1:leng,]
      for (o in 1:leng) {
        raw[o,] = fy(parameter[o,] ,X)
      }
      fy1[i,] <- apply(raw,2,mean)
      num = num + leng
    }else{
      parameter = ind_per[(num+1):(num +leng),]
      for (o in 1:leng) {
        raw[o,] = fy(parameter[o,] ,X)
      }
      fy1[i,] <- apply(raw,2,mean)
      num = num + leng
    }
  }
  for (i in 1:length(tm[,1])) {
    if (i == 1) {
      d_mat[,i] <- d_mat[,i]*ydep
    }else{
      d_mat[,i] <- d_mat[,i]*fy1[i-1,]
    }
  }
  
  # 将d_mat转化为o_mat
  o_mat <- c()
  for (i in 1:length(d_mat[1,]) ) {
    if(i == 1){
      o_mat <- cbind(o_mat,get_origin(d_mat[,i],X,ydep[1]))
    }else{
      o_mat <- cbind(o_mat,get_origin(d_mat[,i],X,0))
    }
  }
  
  o_mat <- as.data.frame(o_mat)
  
  if( dim(o_mat)[2] <= 2 ){
    marker_list[[col]][[5]] <- sum(o_mat[,length(o_mat[1,])]*(X[2]-X[1]))
  }else{
    marker_list[[col]][[5]] <- colSums(o_mat[,2:(length(ind)+1)]*(X[2]-X[1]))
  }
}
filename <- paste0("salt",mod,".Rdata")
save(marker_list,file = filename)