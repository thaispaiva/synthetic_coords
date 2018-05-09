###################################
## looking at the MCMC results

require(fields)

## loading the results
load("/net/nfs1/s/grad/tvp/Dropbox/PHD/Project - Jerry/Death data/results.RData")
dim(avg.lambda)

top4

## ploting the posterior means for the top4 combinations
image.plot(xvec, yvec, matrix(avg.lambda[,top4[1]],size-1,size-1), col=rev(heat.colors(30)), xlab="x", ylab="y", main=paste("Posterior mean of intensity",top4[1]),useRaster=T)
points(jitter(sp[comb==top4[1],],amount=0.1),pch=19,cex=0.8)

image.plot(xvec, yvec, matrix(avg.lambda[,top4[2]],size-1,size-1), col=rev(heat.colors(30)), xlab="x", ylab="y", main=paste("Posterior mean of intensity",top4[2]),useRaster=T)
points(jitter(sp[comb==top4[2],],amount=0.1),pch=19,cex=0.8)

image.plot(xvec, yvec, matrix(avg.lambda[,top4[3]],size-1,size-1), col=rev(heat.colors(30)), xlab="x", ylab="y", main=paste("Posterior mean of intensity",top4[3]),useRaster=T)
points(jitter(sp[comb==top4[3],],amount=0.1),pch=20,cex=0.8)

image.plot(xvec, yvec, matrix(avg.lambda[,top4[4]],size-1,size-1), col=rev(heat.colors(30)), xlab="x", ylab="y", main=paste("Posterior mean of intensity",top4[4]),useRaster=T)
points(jitter(sp[comb==top4[4],],amount=0.1),pch=19,cex=0.8)

## looking at the covariates values for the top4
colnames(X) = c("cancer","sex","race","age","edu")
colnames(b.matrix) = c("comb","cancer","sex","race","age","edu")
b.matrix[top4,]

aux.age = c("<60","60-75","75-85",">85")
aux.edu = c("<high school","high school","some college",">college")

## saving the posterior mean plots

setwd("~/Dropbox/PHD/Project - Jerry/Death data/figures mcmc")

for(i in top4){
  pdf(paste("mean_lambda",i,sep="_",".pdf"), width=15, height=15, pointsize=24)
  image.plot(xvec, yvec, matrix(avg.lambda[,i],size-1,size-1), col=rev(heat.colors(30)), xlab="x", ylab="y", main=paste("Posterior mean of intensity",i),useRaster=T)
  points(jitter(sp[comb==i,],amount=0.5),pch=20,cex=0.8)
  text(7.6,0.2,labels=paste(" not cancer \n female \n white \n",aux.age[b.matrix[i,5]],"\n",aux.edu[b.matrix[i,6]]),adj = c(0,0))
  dev.off()
}


####################################################################################################

## load the synthetic datasets

setwd("~/Project - Jerry/Death data")

synt = list(1); synt = c(synt,2:5)
for(i in 1:5){
  synt[[i]] = read.table(paste("syndata",i,".txt"))
}

####
#### SYNTHETIC LOCATIONS HAVE SWITCHED X AND Y COLUMNS !!!!!!!!!!!!!!!!!!!!
####

par(mfrow=c(3,2))
for(i in 1:5){
  plot(synt[[i]][synt[[i]][,1]==top4[1],3:2]) ## synthetic locations
}
plot(sp[comb==top4[1],1:2]) ## original locations

## saving the synthetic locations

setwd("~/Dropbox/PHD/Project - Jerry/Death data/figures mcmc")

for(i in top4){
  pdf(paste("synt_lambda",i,sep="_",".pdf"), width=30, height=30, pointsize=24)
  par(mfrow=c(2,2))
  for(j in 1:4){
    image.plot(xvec, yvec, matrix(avg.lambda[,i],size-1,size-1), col=rev(heat.colors(30)), xlab="x", ylab="y", useRaster=T)
    points(synt[[j]][synt[[j]][,1]==i,3:2],pch=20,cex=0.8)
  }
  mtext(paste("Synthetic locations of intensity",i), side=3, outer=TRUE, line=-3, cex=2)
  dev.off()  
}

## save the original dataset in the same format as the synthetic datasets

setwd("~/Project - Jerry/Death data")
orig = matrix(c(comb,sp),n,3,byrow=F)
orig = orig[order(comb),]
write(t(orig),file="origdata.txt",ncolumns=3)

###################################################################################################

orig = matrix(c(comb,sp),n,3,byrow=F)
orig = orig[order(comb),]

## plot original vs synthetic locations
plot(orig[,2],syn.data[[1]][,3]) #,xlim=c(2,8),ylim=c(2,8))
plot(orig[,3],syn.data[[1]][,2])

## EDA of original vs synthetic points

setwd("~/Dropbox/PHD/Project - Jerry/Death data/figures synt")

## cancer and not cancer

CANCER_reord = CANCER[order(comb)]

pdf("cancer_orig.pdf", width=15, height=15, pointsize=24, paper="special")
plot(jitter(orig[,2],amount=0.1),jitter(orig[,3],amount=0.1),pch=20,col=CANCER_reord+1,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
# plot(orig[,2],orig[,3],pch=20,col=CANCER_reord+1,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
title("Cancer and immune system - Original data")
legend("topright",col=c(2,1),pch=20,leg=c(paste("cancer (",sum(CANCER_reord),")"),paste("others (",sum(!CANCER_reord),")")))
dev.off()

pdf("cancer_synt.pdf", width=30, height=30, pointsize=24, paper="special")
par(mfrow = c(2,2))
for(i in 1:4){
  ## THE SYNTHETIC DATASET ROWS AND COLS ARE STILL SWITCHED!!!!!!
  plot(syn.data[[i]][,3],syn.data[[i]][,2],pch=20,col=CANCER_reord+1,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
  mtext(paste("Cancer and immune system - Synthetic data"), side=3, outer=TRUE, line=-3, cex=2)
}
dev.off()

## race

RACE = dg.durham2$RACE[order(comb)] - 1  

pdf("race_orig.pdf", width=15, height=15, pointsize=24, paper="special")
plot(jitter(orig[,2],amount=0.1),jitter(orig[,3],amount=0.1),pch=20,col=RACE+1,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
title("Race - Original data")
legend("topright",col=c(1,2),pch=20,leg=c(paste("white (",sum(!RACE),")"),paste("black (",sum(RACE),")")))
dev.off()

pdf("race_synt.pdf", width=30, height=30, pointsize=24, paper="special")
par(mfrow = c(2,2))
for(i in 1:4){
  ## THE SYNTHETIC DATASET ROWS AND COLS ARE STILL SWITCHED!!!!!!
  plot(syn.data[[i]][,3],syn.data[[i]][,2],pch=20,col=RACE+1,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
  mtext(paste("Race - Synthetic data"), side=3, outer=TRUE, line=-3, cex=2)
}
dev.off()

####################################################################################################


