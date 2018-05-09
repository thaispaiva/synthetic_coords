######################
## UTILITY MEASURES ##
######################

require(spBayes)
require(coda)
require(xtable)

########################
## SPATIAL REGRESSION ##
########################

## Summary of results of the spatial regression

# ## load the workspaces and the glm results
load("/net/nfs1/s/grad/tvp/Project - Jerry/Death data/y tilda/size20/utility/utility_size20_orig_nonadapt.RData")
glm.orig = glm.sp

glm.synt2 = list(1); glm.synt2 = c(glm.synt2,2:10)
load("/net/nfs1/s/grad/tvp/Project - Jerry/Death data/y tilda/size20/utility/utility_size20_syntp1_nonadapt.RData")
glm.synt2[[1]] = glm.synt[[1]]; glm.synt2[[2]] = glm.synt[[2]];
load("/net/nfs1/s/grad/tvp/Project - Jerry/Death data/y tilda/size20/utility/utility_size20_syntp2_nonadapt.RData")
glm.synt2[[3]] = glm.synt[[3]]; glm.synt2[[4]] = glm.synt[[4]]
load("/net/nfs1/s/grad/tvp/Project - Jerry/Death data/y tilda/size20/utility/utility_size20_syntp3_nonadapt.RData")
glm.synt2[[5]] = glm.synt[[5]]; glm.synt2[[6]] = glm.synt[[6]]
load("/net/nfs1/s/grad/tvp/Project - Jerry/Death data/y tilda/size20/utility/utility_size20_syntp4_nonadapt.RData")
glm.synt2[[7]] = glm.synt[[7]]; glm.synt2[[8]] = glm.synt[[8]]
load("/net/nfs1/s/grad/tvp/Project - Jerry/Death data/y tilda/size20/utility/utility_size20_syntp5_nonadapt.RData")
glm.synt2[[9]] = glm.synt[[9]]; glm.synt2[[10]] = glm.synt[[10]]
# 
# save.image("/net/nfs1/s/grad/tvp/Dropbox/PHD/Project - Jerry/Death data/y tilda/size20/utility/summary_utility.RData")

# load("/home/grad/tvp/Project - Jerry/Death data/y tilda/size20/utility/summary_utility.RData")
# rm(glm.synt)

## posterior means and HPD intervals for original data

# print(summary(mcmc(glm.orig$p.samples)))
# plot(mcmc(glm.orig$p.samples))
# # burn-in = 50000
# print(summary(mcmc(glm.orig$p.samples,start=50001,end=100000)))

mean.orig = apply(glm.orig$p.samples[50001:100000,],MAR=2,FUN=mean)
n.par = length(mean.orig)

ci.orig = matrix(0,2,n.par)
for(i in 1:n.par){
  ci.orig[,i] = c(HPDinterval(mcmc(glm.orig$p.samples[50001:100000,i])))
}

res.orig = rbind(mean.orig,ci.orig)
row.names(res.orig) = c("mean","LB","UB")
xtable(res.orig)

## posterior means and HPD intervals for synthetic data

# plot(glm.synt2[[1]]$p.samples) # burn-in = 20000
# plot(glm.synt2[[2]]$p.samples) # burn-in = 40000
# plot(glm.synt2[[3]]$p.samples) # burn-in = 20000
# plot(glm.synt2[[4]]$p.samples) # burn-in = 40000
# plot(glm.synt2[[5]]$p.samples) # burn-in = 20000
# plot(glm.synt2[[6]]$p.samples) # burn-in = 20000
# plot(glm.synt2[[7]]$p.samples) # burn-in = 40000
# plot(glm.synt2[[8]]$p.samples) # burn-in = 40000
# plot(glm.synt2[[9]]$p.samples) # burn-in = 30000
# plot(glm.synt2[[10]]$p.samples) # burn-in = 40000

# combine the samples - overall burn-in = 40000
synt.samples = matrix(0,500000,11)
synt.samples = rbind(glm.synt2[[1]]$p.samples[50001:100000,],glm.synt2[[2]]$p.samples[50001:100000,],
      glm.synt2[[3]]$p.samples[50001:100000,],glm.synt2[[4]]$p.samples[50001:100000,],
      glm.synt2[[5]]$p.samples[50001:100000,],glm.synt2[[6]]$p.samples[50001:100000,],
      glm.synt2[[7]]$p.samples[50001:100000,],glm.synt2[[8]]$p.samples[50001:100000,],
      glm.synt2[[9]]$p.samples[50001:100000,],glm.synt2[[10]]$p.samples[50001:100000,])

mean.synt = apply(synt.samples,MAR=2,FUN=mean)

ci.synt = matrix(0,2,n.par)
for(i in 1:n.par){
  ci.synt[,i] = c(HPDinterval(mcmc(synt.samples[,i])))
}
res.synt = rbind(mean.synt,ci.synt)
row.names(res.synt) = c("mean","LB","UB")
xtable(res.synt)

## plot to compare the estimations from original vs synthetic

setwd("/net/nfs1/s/grad/tvp/Dropbox/PHD/Project - Jerry/Death data/y tilda/size20/utility")

pdf("spreg_plot.pdf", width=12, height=12, pointsize=20)

aux = range(res.orig,res.synt)
plot(1:11-0.1,mean.orig,ylim=aux,xlim=c(0.5,11.5),pch=20,xaxt="n",xlab="parameters",ylab=" ")
segments(1:11-0.1,ci.orig[1,],1:11-0.1,ci.orig[2,],lwd=3)
points(1:11+0.1,mean.synt,pch=20,col="grey")
segments(1:11+0.1,ci.synt[1,],1:11+0.1,ci.synt[2,],lwd=3,col="grey")
axis(1,1:11,c("Intercept","sex","race","age2","age3","age4","edu2","edu3","edu4","sigma.sq","phi"),cex.axis=0.7)
legend("topleft",legend=c("original","synthetic"),pch=20,lwd=2,col=c("black","grey"))
title("Posterior means and 95% HPD intervals
obtained with original and synthetic datasets")

dev.off()

######################################################################################################################
