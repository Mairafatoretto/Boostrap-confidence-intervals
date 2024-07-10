
##################################################################
# Fatoretto et al.  - Method 1 -#
##################################################################

#packages

set.seed(23)
require(lme4)
require(HDInterval)
require(pbmcapply)
require(mvtnorm)
require(glmmTMB)
require(pbmcapply)


#parameters used to simulation

I <- 4#nº Blocks
J <- 14 #nº Isolates
K <- 5 #nº Dose

m <- 100 #total of conidias


#X matrix
BLOCK <- rep(c(1:I),each = J*K)
FUNG <- rep(c(1:J),each = 5,I)
DOSE <- rep(c(0,2,4,6,8), I*J)
IND <- c(1:(I*J*K))
sig_ind <-0.9
cov.matrix <- matrix(c(2,-0.4,-0.4,0.1), ncol=2)


#Simulated parameters
Bl1 <- 8
Bl2 <- 7
Bl3 <- 8
Bl4 <- 8  
#FUNG1296 <- 0
FUNG1741 <- -1.5
FUNG1998 <- -1.5
FUNG2778 <- -1.5
FUNG3300 <- -2.5
FUNG3302 <- -2.5
FUNG3307 <- -3
FUNGA152.I <- -2.5
FUNGC23.I <- -2
FUNGE149.I <- -2
FUNGE71.I <- -2.5
FUNGIT1.I2 <- -1.17
FUNGIT2.I10 <- -1.4
FUNGSB4.I7 <- -1.5
SLOPE <- -0.9



#BETA
beta <- c(Bl1,Bl2,Bl3,Bl4,FUNG1741,FUNG1998,FUNG2778,FUNG3300,FUNG3302,FUNG3307,FUNGA152.I, FUNGC23.I,FUNGE149.I,FUNGE71.I,
          FUNGIT1.I2,FUNGIT2.I10,FUNGSB4.I7,SLOPE)



#X matrix
Cl <- kronecker(diag(1, I), array(1, dim=c(J*K,1)))
Bl1 <- matrix(rep(t(kronecker(diag(1, J), array(1, dim=c(K,1)))),I) , ncol =  ncol(kronecker(diag(1, J), array(1, dim=c(K,1)))) , byrow = TRUE )[,-1]
X <- cbind(Cl,Bl1,DOSE)


#Z matrix
Z <- diag(1,nrow=I*J*K)
Z2 <- cbind(kronecker(diag(1, I*J), array(1, dim=c(K,1))))


nboot <- 1000
ncores <- 4

#simulation
coverage <-  pbmclapply(1:100, function(i){  
  
  b <- rnorm(I*J*K, 0, sig_ind)
  b2 <- rmvnorm(I*J, mean=c(0,0), sigma=cov.matrix,method="chol")
  
  intercept <- b2[,1]
  slope <- rep(b2[,2],each=K)
  
  eta <- X%*%beta+Z2%*%intercept+DOSE*slope+Z%*%b
  pi <- plogis(eta)
  y <- rbinom(I * J * K, m, pi)
  
  dataset <- data.frame(IND,BLOCK,FUNG,DOSE,y,m)
  dataset$FUNG <- as.factor(dataset$FUNG)
  dataset$BLOCK <- as.factor(dataset$BLOCK)
  dataset$IND <- as.factor(dataset$IND)
  
  resp <- with(dataset,cbind(y,m-y))
  
  fit_model <- glmmTMB(cbind(y,m-y)~ -1+BLOCK + FUNG+DOSE+(DOSE|BLOCK:FUNG) +(1|IND),family=binomial,data=dataset)
  
  
  #function to bootrstrap resampling
  sampler <- function(dat, clustervar, replace = TRUE, reps = 1) {
    cid <- unique(dat[, clustervar[1]]) 
    ncid <- length(cid)
    recid <- sample(cid, size = ncid * reps, replace = TRUE) 
    if (replace) {
      rid <- lapply(seq_along(recid), function(i) {
        cbind(NewID = i, RowID = sample(which(dat[, clustervar] == recid[i]),
                                        size = length(which(dat[, clustervar] == recid[i])), replace = TRUE))
      })
    } else {
      rid <- lapply(seq_along(recid), function(i) {
        cbind(NewID = i, RowID = which(dat[, clustervar] == recid[i]))
      })
    }
    dat <- as.data.frame(do.call(rbind, rid))
    dat$Replicate <- factor(cut(dat$NewID, breaks = c(1, ncid * 1:reps), include.lowest = TRUE,
                                labels = FALSE))
    dat$NewID <- factor(dat$NewID)
    return(dat)
  }
  
  #M1
  tmp <- sampler(dataset, "IND", replace=FALSE,reps = nboot)
  
  #For M2 and M3 change for this.
  #tmp <- sampler(dataset, "BLOCK",reps = nboot)
  
  
  bigdata <- cbind(tmp, dataset[tmp$RowID, ])
  
  
  #newdata
  resp = with(bigdata,cbind(y,m-y))
  
  
  #fiiting new datas for M1 and M2
  fit_boot_model <- function(i) {
    object <- try(glmmTMB(cbind(y,m-y)~ -1+BLOCK+FUNG+DOSE+(DOSE|BLOCK:FUNG) +(1|IND),   
                          family=binomial,data=droplevels(subset                           (bigdata,bigdata$Replicate==i))),
                  silent = TRUE)
    if (class(object) == "try-error")
      return(object)
    c(fixef(object), VarCorr(object))
  }
  
  #for M3 change for this
  #fit_boot_model <- function(i) {
  #  object <- try(glmmTMB(cbind(y,m-y)~ -1+NewID+FUNG+DOSE+(DOSE|BLOCK:FUNG) +(1|IND),   
  #                        family=binomial,data=droplevels(subset                           (bigdata,bigdata$Replicate==i))),
  #                silent = TRUE)
  #  if (class(object) == "try-error")
  #    return(object)
  #  c(fixef(object), VarCorr(object))
  #}
  
  
  
  start <- proc.time()
  res <-pbmclapply(1:nboot, FUN = fit_boot_model,mc.cores = 8L)   
  end <- proc.time()
  
  
  res <- res[lapply(res,function(x) length(grep("Error",x,value=FALSE))) == 0]
  
  
  
  NewID1 <- c(1:length(res))
  NewID2 <- c(1:length(res))
  NewID3 <- c(1:length(res))
  NewID4 <- c(1:length(res))
  FUNG2 <- c(1:length(res))
  FUNG3 <- c(1:length(res))
  FUNG4 <- c(1:length(res))
  FUNG5 <- c(1:length(res))
  FUNG6 <- c(1:length(res))
  FUNG7 <- c(1:length(res))
  FUNG8 <- c(1:length(res))
  FUNG9 <- c(1:length(res))
  FUNG10 <- c(1:length(res))
  FUNG11 <- c(1:length(res))
  FUNG12 <- c(1:length(res))
  FUNG13 <- c(1:length(res))
  FUNG14 <- c(1:length(res))
  DOSE <- c(1:length(res))
  ind <- c(1:length(res))
  sigma_I <- c(1:length(res))
  sigma_S <- c(1:length(res))
  sigma_IS <- c(1:length(res))
  
  for(i in 1:length(res)) {
    NewID1[i] <- res[[i]][[1]][1]
    NewID2[i] <- res[[i]][[1]][2]
    NewID3[i] <- res[[i]][[1]][3]
    NewID4[i] <- res[[i]][[1]][4]
    FUNG2[i] <- res[[i]][[1]][5]
    FUNG3[i] <- res[[i]][[1]][6]
    FUNG4[i] <- res[[i]][[1]][7]
    FUNG5[i] <- res[[i]][[1]][8]
    FUNG6[i] <- res[[i]][[1]][9]
    FUNG7[i] <- res[[i]][[1]][10]
    FUNG8[i] <- res[[i]][[1]][11]
    FUNG9[i] <- res[[i]][[1]][12]
    FUNG10[i] <- res[[i]][[1]][13]
    FUNG11[i] <- res[[i]][[1]][14]
    FUNG12[i] <- res[[i]][[1]][15]
    FUNG13[i] <- res[[i]][[1]][16]
    FUNG14[i] <- res[[i]][[1]][17]
    DOSE[i] <- res[[i]][[1]][18]
    
  }
  
  for(i in 1:length(res)) {
    ind[i] <- attr(res[[i]][[4]][['IND']],'stddev')
    sigma_I[i] <- attr(res[[i]][[4]][['BLOCK:FUNG']],'stddev')["(Intercept)"]
    sigma_S[i] <- attr(res[[i]][[4]][['BLOCK:FUNG']],'stddev')["DOSE"]
    sigma_IS[i] <- attr(res[[i]][[4]][['BLOCK:FUNG']],'correlation')[2,1]
  }
  
  
  
  
  fixed_boot <- data.frame(NewID1,NewID2,NewID3,NewID4,FUNG2,FUNG3,FUNG4,FUNG5,FUNG6,FUNG7,FUNG8,FUNG9,FUNG10,FUNG11,FUNG12,FUNG13,FUNG14,DOSE)
  
  
  random_boot <- data.frame(ind,sigma_I,sigma_S,sigma_IS)
  
  
  #confidence intervals
  list(
    "random" <- VarCorr(fit_model),
    "fixed" <- fixef(fit_model),
    "coefs" <- data.frame(fixed_boot,random_boot),
    "hpd" = t(apply(data.frame(fixed_boot,random_boot), 2, hdi)),
    "qtl" = t(apply(data.frame(fixed_boot,random_boot), 2, quantile, probs = c(0.025, 0.975)))
  )
},mc.cores=ncores)

  


#######################################################################
##############Calculating coverage rate################################
#######################################################################

total <- coverage[1:100]


for(i in 1:22) {                   
  assign(paste0("coverage.", i), c(1:100))
}

#update random intervals
theta_random <- lapply(1:100, function(i){  
  list(
    "hpd" = t(apply(total[[i]][[3]][22]*total[[i]][[3]][21]*total[[i]][[3]][20], 2, hdi)),
    "qtl" = t(apply(total[[i]][[3]][22]*total[[i]][[3]][21]*total[[i]][[3]][20], 2, quantile, probs = c(0.025, 0.975)))    
  )
}) 


orig_values <- list(8,7,8,8,-1.5,-1.5,-1.5,-2.5,-2.5,-3,-2.5,
-2,-2,-2.5,-1.17,-1.4,-1.5,-0.9,0.9,1.4,0.33,-0.4)

coverage_total. <- c(1:21)

for (j in 1:21){
  Object = get(paste0("coverage.", j))
  for(i in 1:100){
    if(total[[i]][['qtl']][j, 1]<=orig_values[[j]] & total[[i]][['qtl']][j, 2]>=orig_values[[j]]){
      Object[i] <- 1
    }
    else{ 
      Object[i] <- 0
    }
  }
  assign(paste0("Coverage.", j), Object)
  coverage_total.[j] <- sum(get(paste0("Coverage.", j)))/100 
}

Coverage.22 <- c(1:100)
for(i in 1:100){ 
  if (theta_random[[i]][['qtl']][1,1]<=-0.4 & theta_random[[i]][['qtl']][1,2]>=-0.4)
  {
    Coverage.22[i]<- 1
  }
  else{ 
    Coverage.22[i] <-  0
  }
  coverage.22<- sum(Coverage.22)/100
}

#getting coverage rate
coverage_total.
coverage.22

