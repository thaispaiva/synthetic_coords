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

# post is the n x m matrix with the posterior probability distributions for each observation over the entire grid
# di is the true grid cells

tau.r = numeric(n)
for(i in 1:n){
  tau.r[i] = ifelse(which(post[i,]==max(post[i,]))==di[i],1,0)
}
tau.bar = mean(tau.r)
cat("\n\n tau.bar = ",tau.bar)

write(t(tau.r),file="/home/grad/tvp/Dropbox/PHD/Project - Jerry/Death data/y tilda/summary all sizes/risk/tau_size20.txt",ncol=1)

## saving tau & di
# setwd("/net/nfs1/s/grad/tvp/Dropbox/PHD/Project - Jerry/Death data/y tilda/size20/risk")
# write(t(tau.r),file="tau.txt",ncol=1)
# write(t(di),file="di.txt",ncol=1)

## histogram of the prob at the true grid cell
# pi.di = numeric(n)
# for(i in 1:n){
#   pi.di[i] = post[i,di[i]]
# }
# pdf("hist_pidi.pdf", width=12, height=12, pointsize=24)
# hist(pi.di,main="Histogram of pi^{di}")
# dev.off()

####################
## theta.far = proportion of similar points that are within the neighborhood

W.aux = W
diag(W.aux) = 1  ## to include the true gridcell when calculating the proportions

theta.far = numeric(n)
for(i in 1:n){
  theta.far[i] = sum( counts[ comb[i] , which(W.aux[di[i],]>0) ] )/sum( counts[comb[i],] )
}

## boxplot of theta.far based on tau.r

# pdf("theta_far.pdf", width=12, height=12, pointsize=24)
# boxplot(theta.far~tau.r,ylab="theta.far",xlab="tau.r",main="Proportion of similar points that 
#  are within the neighborhood, based on tau.r")
# dev.off()

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

# pdf("ne_dist.pdf", width=12, height=12, pointsize=24)
# boxplot(ne.dist~tau.r,ylab="ne.dist",xlab="tau.r",main="Distance to nearest similar neighbor,
#  based on tau.r")
# dev.off()

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

write(t(theta.max),file="/home/grad/tvp/Dropbox/PHD/Project - Jerry/Death data/y tilda/summary all sizes/risk/theta_size20.txt",ncol=1)

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

write(t(delta),file="/home/grad/tvp/Dropbox/PHD/Project - Jerry/Death data/y tilda/summary all sizes/risk/delta_size20.txt",ncol=1)

# pdf("hist_delta.pdf", width=12, height=12, pointsize=24)
# hist(delta)
# dev.off()

summary(delta)
xtable(t(as.matrix(summary(delta))))

# distribution of delta_i when tau.r==0 (the max isn't the true gridcell, but how far is it?)

summary(delta[tau.r==0])
xtable(t(as.matrix(summary(delta[tau.r==0]))))

# proportion of delta_i < 1
prop_less1 = sum(delta[tau.r==0] < 1)/length(delta[tau.r==0])
round(prop_less1,4)

# pdf("hist_delta_taur0.pdf", width=12, height=12, pointsize=24)
# hist(delta[tau.r==0],main=expression(paste("Histogram of ",delta,"|",tau[r]==0)),freq=F,col=c(1,rep("white",11)),
#      labels=c(round(prop_less1,3),rep(" ",11)),
#      xlab=expression(paste(delta,"|",tau[r]==0)))
# dev.off()

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

write(t(numb.neigh),file="/home/grad/tvp/Dropbox/PHD/Project - Jerry/Death data/y tilda/summary all sizes/risk/nneigh_size20.txt",ncol=1)

# points that are alone
alone = which(numb.neigh==0) 

# proportion of tau.r=1 for the alone points
sum(tau.r[alone])/length(alone)
round(sum(tau.r[alone])/length(alone),4)

1-sum(tau.r[alone])/length(alone)
round(1-sum(tau.r[alone])/length(alone),4)




