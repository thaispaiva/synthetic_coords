## Looking at the results of the posterior risk measures

require(fields)
require(xtable)

setwd("/home/grad/tvp/Dropbox/PHD/Project - Jerry/Death data/y tilda/size20")

## load the workspace
# load("/net/nfs1/s/grad/tvp/Dropbox/PHD/Project - Jerry/Death data/results2.RData")
load("/net/nfs1/s/grad/tvp/Project - Jerry/Death data/y tilda/size20/results_size20.RData")

colnames(X) = c("cancer","sex","race","age","edu")
colnames(b.matrix) = c("comb","cancer","sex","race","age","edu")

## calculate tau and tau.bar

tau.r = numeric(n)
for(i in 1:n){
  tau.r[i] = ifelse(which(post[i,]==max(post[i,]))==di[i],1,0)
}
tau.bar = mean(tau.r)
cat("\n\n tau.bar = ",tau.bar)

## saving tau & di
setwd("/net/nfs1/s/grad/tvp/Dropbox/PHD/Project - Jerry/Death data/y tilda/size20/risk")
write(t(tau.r),file="tau.txt",ncol=1)
write(t(di),file="di.txt",ncol=1)

## histogram of the prob at the true grid cell
pi.di = numeric(n)
for(i in 1:n){
  pi.di[i] = post[i,di[i]]
}
pdf("hist_pidi.pdf", width=12, height=12, pointsize=24)
hist(pi.di,main="Histogram of pi^{di}")
dev.off()


## select some individuals to look at their values

# ## first, individuals that have tau=1
# tau1 = which(tau.r>0)
# length(tau1)
# tau1
# # i = sample(tau1,1)  # select different observations
# 
# ## check how many points had tau.r=1 in each combination
# table(comb[tau1])[table(comb[tau1])>0]

####################
## theta.far = proportion of similar points that are within the neighborhood

W.aux = W
diag(W.aux) = 1  ## to include the true gridcell when calculating the proportions

theta.far = numeric(n)
for(i in 1:n){
  theta.far[i] = sum( counts[ comb[i] , which(W.aux[di[i],]>0) ] )/sum( counts[comb[i],] )
}

## boxplot of theta.far based on tau.r

pdf("theta_far.pdf", width=12, height=12, pointsize=24)
boxplot(theta.far~tau.r,ylab="theta.far",xlab="tau.r",main="Proportion of similar points that 
 are within the neighborhood, based on tau.r")
dev.off()

# ## plot the points with colors based on tau.r
# rnd.comb = sample(b,4) # select a few combinations
# 
# par(mfrow=c(2,2))
# for(i in 1:4){
#   plot(sp[comb==rnd.comb[i],],pch=20,col=tau.r[comb==rnd.comb[i]]+1)
# }

# dev.off()

## look at the shortest distance from the points that have tau.r=1 to the nearest neighbor

ne.dist = numeric(n)
for(i in 1:n){
  aux = sp[-i,]
  ne.dist[i] = min(c(nearest.dist(matrix(sp[i,],1,2),matrix(aux[comb[-i]==comb[i],],ncol=2),delta=15)))
  # return the nearest distance from each point to the other points in the same combination
}
# hist(ne.dist)

pdf("ne_dist.pdf", width=12, height=12, pointsize=24)
boxplot(ne.dist~tau.r,ylab="ne.dist",xlab="tau.r",main="Distance to nearest similar neighbor,
 based on tau.r")
dev.off()

#################

## look at the unique observations

# sort(table(comb))[2:6]
# comb.uni = c(95,92,96,97)
# 
# ## plotting the posterior mean of the combinations with few obs
# pdf("post_uni.pdf", width=24, height=24, pointsize=24)
# par(mfrow=c(2,2))
# i = 1
#   image.plot(xvec, yvec, matrix(avg.lambda[,comb.uni[i]],size-1,size-1), col=rev(heat.colors(30)), xlab="s1", ylab="s2", main=paste("Posterior mean of intensity",comb.uni[i]),useRaster=T)
#   points(jitter(sp[comb==comb.uni[i],1:2],amount=0.1),pch=19,cex=0.8)
#   text(7.2,0.2,labels=paste(" cancer \n male \n black \n >85 \n college "),adj = c(0,0))
# i = 2
#   image.plot(xvec, yvec, matrix(avg.lambda[,comb.uni[i]],size-1,size-1), col=rev(heat.colors(30)), xlab="s1", ylab="s2", main=paste("Posterior mean of intensity",comb.uni[i]),useRaster=T)
#   points(jitter(sp[comb==comb.uni[i],1:2],amount=0.1),pch=19,cex=0.8)
#   text(7.2,0.2,labels=paste(" cancer \n male \n black \n 75-85 \n >college "),adj = c(0,0))
# i = 3
#   image.plot(xvec, yvec, matrix(avg.lambda[,comb.uni[i]],size-1,size-1), col=rev(heat.colors(30)), xlab="s1", ylab="s2", main=paste("Posterior mean of intensity",comb.uni[i]),useRaster=T)
#   points(jitter(sp[comb==comb.uni[i],1:2],amount=0.1),pch=19,cex=0.8)
#   text(7.2,0.2,labels=paste(" cancer \n male \n black \n >85 \n >college "),adj = c(0,0))
# i = 4
#   image.plot(xvec, yvec, matrix(avg.lambda[,comb.uni[i]],size-1,size-1), col=rev(heat.colors(30)), xlab="s1", ylab="s2", main=paste("Posterior mean of intensity",comb.uni[i]),useRaster=T)
#   points(jitter(sp[comb==comb.uni[i],1:2],amount=0.1),pch=19,cex=0.8)
#   text(7.2,0.2,labels=paste(" cancer \n female \n white \n <60 \n <high school"),adj = c(0,0))
# dev.off()


## plotting the original and synthetic locations

# i = 1
# pdf(paste("uni_orig_synt",i,".pdf",sep=""), width=24, height=24, pointsize=24)
# image.plot(xvec, yvec, matrix(avg.lambda[,comb.uni[i]],size-1,size-1), col=rev(heat.colors(30)), xlab="s1", ylab="s2", main=paste("Posterior mean of intensity",comb.uni[i]),useRaster=T)
# points(jitter(sp[comb==comb.uni[i],1:2],amount=0.1),pch=19,cex=0.8)
# text(7.2,0.2,labels=paste(" cancer \n male \n black \n >85 \n college "),adj = c(0,0))
# for(j in 1:5){
# #   points(jitter(syn.data[[j]][syn.data[[j]][,1]==comb.uni[i],3:2],amount=0.1),cex=0.8)  
#   points(jitter(syn.data[[j]][syn.data[[j]][,1]==comb.uni[i],2:3],amount=0.1),cex=0.8)
# }
# dev.off()
# 
# i = 2
# pdf(paste("uni_orig_synt",i,".pdf",sep=""), width=24, height=24, pointsize=24)
# image.plot(xvec, yvec, matrix(avg.lambda[,comb.uni[i]],size-1,size-1), col=rev(heat.colors(30)), xlab="s1", ylab="s2", main=paste("Posterior mean of intensity",comb.uni[i]),useRaster=T)
# points(jitter(sp[comb==comb.uni[i],1:2],amount=0.1),pch=19,cex=0.8)
# text(7.2,0.2,labels=paste(" cancer \n male \n black \n >85 \n college "),adj = c(0,0))
# for(j in 1:5){
#   points(jitter(syn.data[[j]][syn.data[[j]][,1]==comb.uni[i],2:3],amount=0.1),cex=0.8)  
# }
# dev.off()
# 
# i = 3
# pdf(paste("uni_orig_synt",i,".pdf",sep=""), width=24, height=24, pointsize=24)
# image.plot(xvec, yvec, matrix(avg.lambda[,comb.uni[i]],size-1,size-1), col=rev(heat.colors(30)), xlab="s1", ylab="s2", main=paste("Posterior mean of intensity",comb.uni[i]),useRaster=T)
# points(jitter(sp[comb==comb.uni[i],1:2],amount=0.1),pch=19,cex=0.8)
# text(7.2,0.2,labels=paste(" cancer \n male \n black \n >85 \n college "),adj = c(0,0))
# for(j in 1:5){
#   points(jitter(syn.data[[j]][syn.data[[j]][,1]==comb.uni[i],2:3],amount=0.1),cex=0.8)  
# }
# dev.off()
# 
# i = 4
# pdf(paste("uni_orig_synt",i,".pdf",sep=""), width=24, height=24, pointsize=24)
# image.plot(xvec, yvec, matrix(avg.lambda[,comb.uni[i]],size-1,size-1), col=rev(heat.colors(30)), xlab="s1", ylab="s2", main=paste("Posterior mean of intensity",comb.uni[i]),useRaster=T)
# points(jitter(sp[comb==comb.uni[i],1:2],amount=0.1),pch=19,cex=0.8)
# text(7.2,0.2,labels=paste(" cancer \n male \n black \n >85 \n college "),adj = c(0,0))
# for(j in 1:5){
#   points(jitter(syn.data[[j]][syn.data[[j]][,1]==comb.uni[i],2:3],amount=0.1),cex=0.8)  
# }
# dev.off()

## plot the risk and the points

## first comb, with two unique obs
# j = 1
# i = which(comb==comb.uni[j])
# 
# pdf(paste("uni_obs_comb",comb.uni[j],".pdf",sep=""), width=36, height=12, pointsize=24)
# par(mfrow=c(1,3))
# for(k in 1:length(i)){
#   ## plot the risk and the points
#   image.plot(xvec,yvec,matrix(post[i[k],],size-1,size-1),col=rev(heat.colors(30)),main=paste("Obs ",i[k]))
#   points(sp[i[k],1],sp[i[k],2],pch=20)
#   points(sp[comb==comb[i[k]],1:2])
# }
# ## plot the average lambda
# image.plot(xvec, yvec, matrix(avg.lambda[,comb.uni[j]],size-1,size-1), col=rev(heat.colors(30)), xlab="s1", ylab="s2", main=paste("Posterior mean of intensity",comb.uni[j]),useRaster=T)
# points(sp[comb==comb.uni[j],1:2],pch=19,cex=0.8)
# dev.off()
# 
# ## second and third comb, with three unique obs
# j = 2
# i = which(comb==comb.uni[j])
# 
# pdf(paste("uni_obs_comb",comb.uni[j],".pdf",sep=""), width=24, height=24, pointsize=24)
# par(mfrow=c(2,2))
# for(k in 1:length(i)){
#   ## plot the risk and the points
#   image.plot(xvec,yvec,matrix(post[i[k],],size-1,size-1),col=rev(heat.colors(30)),main=paste("Obs ",i[k]))
#   points(sp[i[k],1],sp[i[k],2],pch=20)
#   points(sp[comb==comb[i[k]],1:2])
# }
# ## plot the average lambda
# image.plot(xvec, yvec, matrix(avg.lambda[,comb.uni[j]],size-1,size-1), col=rev(heat.colors(30)), xlab="s1", ylab="s2", main=paste("Posterior mean of intensity",comb.uni[j]),useRaster=T)
# points(sp[comb==comb.uni[j],1:2],pch=19,cex=0.8)
# dev.off()
# 
# j = 3
# i = which(comb==comb.uni[j])
# 
# pdf(paste("uni_obs_comb",comb.uni[j],".pdf",sep=""), width=24, height=24, pointsize=24)
# par(mfrow=c(2,2))
# for(k in 1:length(i)){
#   ## plot the risk and the points
#   image.plot(xvec,yvec,matrix(post[i[k],],size-1,size-1),col=rev(heat.colors(30)),main=paste("Obs ",i[k]))
#   points(sp[i[k],1],sp[i[k],2],pch=20)
#   points(sp[comb==comb[i[k]],1:2])
# }
# ## plot the average lambda
# image.plot(xvec, yvec, matrix(avg.lambda[,comb.uni[j]],size-1,size-1), col=rev(heat.colors(30)), xlab="s1", ylab="s2", main=paste("Posterior mean of intensity",comb.uni[j]),useRaster=T)
# points(sp[comb==comb.uni[j],1:2],pch=19,cex=0.8)
# dev.off()
# 
# ## fourth comb, with four unique obs
# j = 4
# i = which(comb==comb.uni[j])
# 
# pdf(paste("uni_obs_comb",comb.uni[j],".pdf",sep=""), width=36, height=24, pointsize=24)
# par(mfrow=c(2,3))
# for(k in 1:length(i)){
#   ## plot the risk and the points
#   image.plot(xvec,yvec,matrix(post[i[k],],size-1,size-1),col=rev(heat.colors(30)),main=paste("Obs ",i[k]))
#   points(sp[i[k],1],sp[i[k],2],pch=20)
#   points(sp[comb==comb[i[k]],1:2])
# }
# ## plot the average lambda
# image.plot(xvec, yvec, matrix(avg.lambda[,comb.uni[j]],size-1,size-1), col=rev(heat.colors(30)), xlab="s1", ylab="s2", main=paste("Posterior mean of intensity",comb.uni[j]),useRaster=T)
# points(sp[comb==comb.uni[j],1:2],pch=19,cex=0.8)
# dev.off()

###############
## look at some of the combinations with a lot of tau.r=1

# table(comb[tau1])[table(comb[tau1])>0]
# 
# i = 45
# image.plot(xvec,yvec,matrix(avg.lambda[,i],size-1,size-1),col=rev(heat.colors(30)),main=paste("Comb ",i))
# points(sp[comb==i,1:2],pch=21)
# points(sp[comb==i & tau.r==1,1:2],pch=20)
# 
# i = 46
# image.plot(xvec,yvec,matrix(avg.lambda[,i],size-1,size-1),col=rev(heat.colors(30)),main=paste("Comb ",i))
# points(sp[comb==i,1:2],pch=21)
# points(sp[comb==i & tau.r==1,1:2],pch=20)

#################################################################################################################

## theta - measure of risk when tau.r=1

theta.max = numeric(n)

for(i in 1:n){
  if(tau.r[i]>0){
    theta.max[i] = 1/sum( counts[ comb[i] , which(post[i,]==max(post[i,])) ] )
  }
}
hist(theta.max)
hist(theta.max[theta.max>0])

table(theta.max)

xtable(t(table(theta.max)))

# boxplot(theta.far~theta.max)
# boxplot(theta.far[theta.max>0]~theta.max[theta.max>0])
# boxplot(ne.dist~theta.max)

################################################################################################################

## delta measure - distance between true location and the maximum

delta = numeric(n)

for(i in 1:n){
  max.cell = vec2mat(which(post[i,]==max(post[i,])),size)
  max.coord = c((xvec[max.cell[1,2]+1]+xvec[max.cell[1,2]])/2,(yvec[max.cell[1,1]+1]+yvec[max.cell[1,1]])/2)
  delta[i] = dist(rbind(sp[i,],max.coord))
}

pdf("hist_delta.pdf", width=12, height=12, pointsize=24)
hist(delta)
dev.off()

summary(delta)
xtable(t(as.matrix(summary(delta))))

# distribution of delta_i when tau.r==0 (the max isn't the true gridcell, but how far is it?)

summary(delta[tau.r==0])
xtable(t(as.matrix(summary(delta[tau.r==0]))))

# proportion of delta_i < 1
prop_less1 = sum(delta[tau.r==0] < 1)/length(delta[tau.r==0])
prop_less1

pdf("hist_delta_taur0.pdf", width=12, height=12, pointsize=24)
hist(delta[tau.r==0],main=expression(paste("Histogram of ",delta,"|",tau[r]==0)),freq=F,col=c(1,rep("white",11)),
     labels=c(round(prop_less1,3),rep(" ",11)),
     xlab=expression(paste(delta,"|",tau[r]==0)))
dev.off()

# distribution of delta_i

boxplot(delta~tau.r)

#########################

## Pr(true match|alone)

# first, find the cases that are alone in their original grid cell
numb.neigh = numeric(n)
for(i in 1:n){
  numb.neigh[i] = sum(di[comb==comb[i]]==di[i]) - 1
}
hist(numb.neigh)
table(numb.neigh)
sum(numb.neigh==0)/n  # 40% of the points are alone

# points that are alone
alone = which(numb.neigh==0) 

# proportion of tau.r=1 for the alone points
sum(tau.r[alone])/length(alone)
round(sum(tau.r[alone])/length(alone),4)

1-sum(tau.r[alone])/length(alone)
round(1-sum(tau.r[alone])/length(alone),4)

## WHAT IF WE LOOK AT ALL THE NEIGHBORS, NOT ONLY THE SIMILAR ONES?

# numb.allneigh = numeric(n)
# for(i in 1:n){
#   numb.allneigh[i] = sum(di==di[i]) - 1
# }
# hist(numb.allneigh)
# table(numb.allneigh)
# sum(numb.allneigh==0)/n  # 0.4% of the points are alone
# # points that are alone
# alone = which(numb.allneigh==0) 
# # proportion of tau.r=1 for the alone points
# sum(tau.r[alone])/length(alone)
# 1-sum(tau.r[alone])/length(alone)



