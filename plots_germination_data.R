#################################################################
##########Figure 6, 7 and 8####################################
################################################################
rm(list = ls())
#packages
require(HDInterval)
require(glmmTMB)
require(pbmcapply)
require(ggplot2)
require(reshape2)


cex = 1.7
cores = c('tomato','gold','orange','brown','red','black','darkgrey','purple','green','blue','hotpink','magenta','darkblue','darkgreen')
cores2 = c(1,gray(0.1),gray(0.3),gray(0.5),gray(0.7),1,gray(0.1),gray(0.3),gray(0.5),gray(0.7),1,gray(0.1),gray(0.3),gray(0.5))
shapes = c(1,2,3,4,5,6,7,8,9,14,15,16,17,18)
tracos = c('solid','dashed','dotted','dotdash','longdash','dashed','dotted','dotdash','longdash','solid','dotted','dotdash','longdash','dashed')
tam_leg = 18 
tam_val = 15 
tam_tit=14
fonte = 'Calibri'


#reading dataset 
dataset = read.table('germination_data.txt',header=TRUE,dec=',',sep='\t')

dataset$BLOCK <- as.factor(dataset$BLOCK)
dataset$FUNG <- as.factor(dataset$FUNG)
dataset$resp = with(dataset,cbind(GERM,TOTAL-GERM));
dataset$IND = gl(nrow(dataset),1)

#fitting the model with the original dataset
a <- glmmTMB(resp~-1+BLOCK+FUNG+DOSE+(DOSE|BLOCK:FUNG) +(1|IND),REML=TRUE,family=binomial,data=dataset)
  
  
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


tmp <- sampler(dataset, "BLOCK",reps = 1000)


bigdata <- cbind(tmp, dataset[tmp$RowID, ])


#newdata
resp = with(bigdata,cbind(GERM,TOTAL-GERM))
  
  
#fiiting the model with the resampling dataset
fit_boot_model <- function(i) {
  object <- try(glmmTMB(cbind(GERM,TOTAL-GERM)~ -1+NewID+FUNG+DOSE+(DOSE|BLOCK:FUNG) +(1|IND),  REML=TRUE, 
                        family=binomial,data=droplevels(subset                           (bigdata,bigdata$Replicate==i))),
                silent = TRUE)
  if (class(object) == "try-error")
    return(object)
  c(fixef(object), VarCorr(object))
}

  
  
start <- proc.time()
res <-pbmclapply(1:1000, FUN = fit_boot_model,mc.cores = 8L)    
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
  FUNG2[i] <- res[[i]][[1]]['FUNGFUNG2']
  FUNG3[i] <- res[[i]][[1]]['FUNGFUNG3']
  FUNG4[i] <- res[[i]][[1]]['FUNGFUNG4']
  FUNG5[i] <- res[[i]][[1]]['FUNGFUNG5']
  FUNG6[i] <- res[[i]][[1]]['FUNGFUNG6']
  FUNG7[i] <- res[[i]][[1]]['FUNGFUNG7']
  FUNG8[i] <- res[[i]][[1]]['FUNGFUNG8']
  FUNG9[i] <- res[[i]][[1]]['FUNGFUNG9']
  FUNG10[i] <- res[[i]][[1]]['FUNGFUNG10']
  FUNG11[i] <- res[[i]][[1]]['FUNGFUNG11']
  FUNG12[i] <- res[[i]][[1]]['FUNGFUNG12']
  FUNG13[i] <- res[[i]][[1]]['FUNGFUNG13']
  FUNG14[i] <- res[[i]][[1]]['FUNGFUNG14']
  DOSE[i] <- res[[i]][[1]][['DOSE']]
  
}  
  
fixed <- data.frame(NewID1,NewID2,NewID3,NewID4,FUNG2,FUNG3,FUNG4,FUNG5,FUNG6,FUNG7,FUNG8,FUNG9,FUNG10,FUNG11,FUNG12,FUNG13,FUNG14,DOSE)

  
  
for(i in 1:length(res)) {
  ind[i] <- attr(res[[i]][[4]][['IND']],'stddev')
  sigma_I[i] <- attr(res[[i]][[4]][['BLOCK:FUNG']],'stddev')["(Intercept)"]
  sigma_S[i] <- attr(res[[i]][[4]][['BLOCK:FUNG']],'stddev')["DOSE"]
  sigma_IS[i] <- attr(res[[i]][[4]][['BLOCK:FUNG']],'correlation')[2,1]
}


random <- data.frame(ind,sigma_I,sigma_S,sigma_IS)
bigres <- data.frame(fixed,random)
  
###################################################################################################figure 6#########################################
#########################################################################
fungus = bigres[,5:17]


x = seq(0,8,by=0.01)
reta = (1/4)*(mean(bigres$NewID1) + mean(bigres$NewID2)  +mean(bigres$NewID3)  +mean(bigres$NewID4)) + mean(bigres$DOSE)*x
sigmoid1 = exp(reta)/(1+exp(reta))
pred = data.frame(reta,Fungus="FUNG1",dose=x)


for (i in 1:13) {
  reta = (1/4)*(mean(bigres$NewID1) + mean(bigres$NewID2)  +mean(bigres$NewID3)  +mean(bigres$NewID4)) +  apply(fungus[i],2,FUN=mean) +  mean(bigres$DOSE)*x
  pred = rbind(pred,data.frame(reta,Fungus=names(fungus[i]),dose=x))
}

ord <- c("FUNG1","FUNG2","FUNG3","FUNG4","FUNG5","FUNG6","FUNG7",
         "FUNG8","FUNG9"
         ,"FUNG10","FUNG11","FUNG12","FUNG13","FUNG14")


pred$Fungus <- factor(pred$Fungus,levels=ord)

#Figure 6a  
ggplot(pred, aes(x=dose, y=exp(reta)/(1+exp(reta)), colour=Fungus,fill=Fungus))	+geom_line(size=0.7)+ 
  theme_bw()+theme(plot.title = element_text(hjust = 0.5,size=8,vjust=3), text = element_text(family = fonte),
                   panel.grid = element_blank(), 
                   legend.background = element_rect(color=1,size=0.3),
                   legend.title = element_text(size=tam_leg), 
                   legend.text = element_text(size=tam_val), 
                   axis.title.x=element_text( size=tam_leg), 
                   axis.text.x=element_text(size=tam_val),
                   axis.title.y=element_text( size=tam_leg), 
                   axis.text.y=element_text(size=tam_val),
                   legend.key.height=unit(0.8,"line")) +
  xlab('\nExposure time (hours)\n') +
  ylab('Predicted proportion\n') + scale_color_manual(values=cores,labels = c("1 ","2","3","4","5","6",                                               "7","8","9","10","11","12","13","14"),name = "Isolates")    

fitted = data.frame(valor=c(dataset$GERM/dataset$TOTAL,exp(pred$reta)/(1+exp(pred$reta))), tipo=c(rep('Observado',280),rep('Estimado',length(pred$dose))), fungo=c(as.character(dataset$FUNG),as.character(pred$Fungus)), dose=c(dataset$DOSE,pred$dose))

fitted$fungo <- factor(fitted$fungo,levels=ord)

labels <- c("FUNG1"="1 ","FUNG2"="2","FUNG3"="3","FUNG4"="4","FUNG5"="5","FUNG6"="6","FUNG7"="7","FUNG8"="8","FUNG9"="9","FUNG10"="10","FUNG11"="11","FUNG12"="12","FUNG13"="13","FUNG14"="14")

#Figure 6b
ggplot(subset(fitted, tipo=='Observado'), 
       aes(x=dose, y=valor, colour=fungo)) + 
  geom_point(size=3) +
  geom_line(data=subset(fitted, tipo=='Estimado'), 
            aes(x=dose,y=valor,color=fungo),size=1) + 
  facet_wrap(~fungo,ncol=4,labeller=labeller(fungo = labels)) +
  xlab('\nExposure time (hours)') +
  ylab('Proportion\n') + 
  ggtitle('') + theme_bw()+
  theme(strip.text = element_text(size=tam_leg,),legend.position="",plot.title = element_text(hjust = 0.5,size=15,vjust=3), text = element_text(family = fonte),
        panel.grid = element_blank(), 
        legend.background = element_rect(color=1,size=0.3),
        legend.title = element_text(size=tam_leg), 
        legend.text = element_text(size=tam_val), 
        axis.title.x=element_text( size=tam_leg), 
        axis.text.x=element_text(size=tam_val),
        axis.title.y=element_text( size=tam_leg), 
        axis.text.y=element_text(size=tam_val),
        legend.key.height=unit(0.8,"line"))+
  scale_color_manual(values=cores,labels = c("1 ","2","3","4","5","6","7","8","9","10","11","12","13","14"),name = "Isolates")                                                                                                                                                            

################################################################
##################Figure 7####################################
################################################################

#isolate 1
bigres$boot <-c(gl(nrow(bigres),1))
reta= (1/4)*(bigres$NewID1[1] + bigres$NewID2[1] + bigres$NewID3[1] +bigres$NewID4[1]) + bigres$DOSE[1]*x
pred = data.frame(reta,boot=as.factor(bigres$boot[1]),dose=x)

for (i in 2:1000) {
  reta= (1/4)*(bigres$NewID1[i] + bigres$NewID2[i] + bigres$NewID3[i] +bigres$NewID4[i])+ bigres$DOSE[i]*x
  pred=rbind(pred,data.frame(reta,boot=as.factor(bigres$boot[i]),dose=x))
} 
    
Pred_mat <-  t(mapply(dcast(pred,boot~dose, value.var = "reta")[,-1], FUN=as.numeric))  
   
    
# Pack the estimates for plotting
Estims_plot <- cbind(
  x = x, 
  as.data.frame(t(apply(Pred_mat, 1, function(t) c(
    ci_lower_est = quantile(t, probs = 0.025, names = FALSE), 
    ci_upper_est = quantile(t, probs = 0.975, names = FALSE)
  ))))
)  
   
    
ggplot() + 
  geom_line(alpha=0.5,data=pred, aes(x=dose, y=exp(reta)/(1+exp(reta)), fill=boot)) +
  geom_line(data=Estims_plot, aes(x = x, y = exp(ci_lower_est)/(1+exp(ci_lower_est))),size=1, color = "blue") +
  geom_line(data=Estims_plot, aes(x = x, y = exp(ci_upper_est)/(1+exp(ci_upper_est))), size=1,color = "blue") +
  ggtitle("ISOLATE 1")+ 
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,size=tam_tit,vjust=3), text = element_text(family = fonte), panel.grid = element_blank(), legend.background = element_rect(color=1,size=0.3), legend.title = element_text(size=tam_leg), legend.text = element_text(size=tam_val), axis.title.x=element_text( size=tam_leg), axis.text.x=element_text(size=tam_leg), axis.title.y=element_text( size=tam_leg), axis.text.y=element_text(size=tam_val),legend.key.height=unit(0.8,"line")) + xlab('\nExposure time (hours)\n') + ylab('Predicted proportion\n') 

  
#isolate 2
reta= (1/4)*(bigres$NewID1[1] + bigres$NewID2[1] + bigres$NewID3[1] +bigres$NewID4[1])+ bigres$FUNG2[1]+ bigres$DOSE[1]*x#For others isolates change FUNG2
pred = data.frame(reta,boot=as.factor(bigres$boot[1]),dose=x)

for (i in 2:1000) {
  reta= (1/4)*(bigres$NewID1[i] + bigres$NewID2[i] + bigres$NewID3[i] +bigres$NewID4[i])+  bigres$FUNG2[i]+bigres$DOSE[i]*x#For others isolates change FUNG2
  pred=rbind(pred,data.frame(reta,boot=as.factor(bigres$boot[i]),dose=x))
} 
  
  
pred_mat <-  t(mapply(dcast(pred,boot~dose, value.var = "reta")[,-1], FUN=as.numeric))
  
  
# Pack the estimates for plotting
Estims_plot <- cbind(
  x = x, 
  as.data.frame(t(apply(pred_mat, 1, function(t) c(
    ci_lower_est = quantile(t, probs = 0.025, names = FALSE), 
    ci_upper_est = quantile(t, probs = 0.975, names = FALSE)
  ))))
)
  
ggplot() + 
  geom_line(alpha=0.5,data=pred, aes(x=dose, y=exp(reta)/(1+exp(reta)), fill=boot)) +
    geom_line(data=Estims_plot, aes(x = x, y = exp(ci_lower_est)/(1+exp(ci_lower_est))),size=1, color = "blue") +
    geom_line(data=Estims_plot, aes(x = x, y = exp(ci_upper_est)/(1+exp(ci_upper_est))), size=1,color = "blue") +
      ggtitle("ISOLATE 2") + #For others isolates change the plot title here 
    theme_bw()+ 
    theme(plot.title = element_text(hjust = 0.5,size=tam_tit,vjust=3), text = element_text(family = fonte), panel.grid = element_blank(), legend.background = element_rect(color=1,size=0.3), legend.title = element_text(size=tam_leg), legend.text = element_text(size=tam_val), axis.title.x=element_text( size=tam_leg), axis.text.x=element_text(size=tam_leg), axis.title.y=element_text( size=tam_leg), axis.text.y=element_text(size=tam_val),legend.key.height=unit(0.8,"line")) + xlab('\nExposure time (hours)\n') + ylab('Predicted proportion\n') 
  
  
  
  

############################################################
##################Figure 8 #################################
############################################################

  
  
#isolate 1
bigres$boot <-c(gl(nrow(bigres),1))
reta= (1/4)*(bigres$NewID1[1] + bigres$NewID2[1] + bigres$NewID3[1] +bigres$NewID4[1]) + bigres$DOSE[1]*x
pred = data.frame(reta,boot=as.factor(bigres$boot[1]),dose=x)

for (i in 2:1000) {
  reta= (1/4)*(bigres$NewID1[i] + bigres$NewID2[i] + bigres$NewID3[i] +bigres$NewID4[i])+ bigres$DOSE[i]*x
  pred=rbind(pred,data.frame(reta,boot=as.factor(bigres$boot[i]),dose=x))
} 

require(reshape2)
pred_mat <-  t(mapply(dcast(pred,boot~dose, value.var = "reta")[,-1], FUN=as.numeric))


# Pack the estimates for plotting
Estims_plot <- cbind(
  x = x, 
  as.data.frame(t(apply(pred_mat, 1, function(t) c(
    ci_lower_est = quantile(t, probs = 0.025, names = FALSE), 
    ci_upper_est = quantile(t, probs = 0.975, names = FALSE)
  ))))
)  
  
  
#isolate 2
bigres$boot <-c(gl(nrow(bigres),1))
reta2= (1/4)*(bigres$NewID1[1] + bigres$NewID2[1] + bigres$NewID3[1] +bigres$NewID4[1])+ bigres$FUNG2[1]+ bigres$DOSE[1]*x
pred2 = data.frame(reta2,boot=as.factor(bigres$boot[1]),dose=x)
  
for (i in 2:1000) {
  reta2= (1/4)*(bigres$NewID1[i] + bigres$NewID2[i] + bigres$NewID3[i] +bigres$NewID4[i])+  bigres$FUNG2[i]+bigres$DOSE[i]*x
  pred2=rbind(pred2,data.frame(reta2,boot=as.factor(bigres$boot[i]),dose=x))
} 

  
pred_mat2 <-  t(mapply(dcast(pred2,boot~dose, value.var = "reta2")[,-1], FUN=as.numeric))
  
  
# Pack the estimates for plotting
Estims_plot2 <- cbind(
  x = x, 
  as.data.frame(t(apply(pred_mat2, 1, function(t) c(
    ci_lower_est = quantile(t, probs = 0.025, names = FALSE), 
    ci_upper_est = quantile(t, probs = 0.975, names = FALSE)
  ))))
)


#comparing isolate 1 and isolte2
ggplot()+
  geom_ribbon(data=Estims_plot2, aes(color="cyan4",x=x, ymin = exp(ci_lower_est)/(1+exp(ci_lower_est)), ymax = exp(ci_upper_est)/(1+exp(ci_upper_est))),linetype = 2,alpha=0.2) +
  geom_ribbon(data=Estims_plot, aes(color="brown1",x=x, ymin = exp(ci_lower_est)/(1+exp(ci_lower_est)), ymax = exp(ci_upper_est)/(1+exp(ci_upper_est))),linetype = 2,alpha=0.2) +
  scale_color_identity(
    labels = c(" ISOLATE 1", " ISOLATE 2"),
    guide = "legend")+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        legend.position=c(0.3,0.2),legend.title = element_blank(), 
        legend.text =element_text(size=10), axis.title.x=element_text(size=tam_leg), 
        axis.text.x=element_text(size=tam_val), axis.title.y=element_text(size=tam_leg),
        axis.text.y=element_text(size=tam_val),legend.key.height=unit(0.5,"line")) +
  xlab('\nExposure time (hours)')+ylab('Confidence intervals\n')




  
  
  
  
  