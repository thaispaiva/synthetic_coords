
#####################
## POINT ESTIMATES ##
#####################

require(xtable)

setwd("/net/nfs1/s/grad/tvp/Dropbox/PHD/Project - Jerry/Death data/y tilda/size20/utility")

## save the data that we need (dont need to load the entire workspace just for this part)
## then, we just need to load the data as follows...

# write(t(cbind(death,coords.zip)),"death_coords_zip.txt",ncolumns=12)
death.orig = as.matrix(read.table("death_coords_zip.txt"))
colnames(death.orig) = c("ytil","sex","race","age2","age3","age4","edu2","edu3","edu4","lat","lon","zip")
head(death.orig)

n = length(death.orig[,1])
n.syn = 10
death.synt=list(1); death.synt=c(death.synt,2:n.syn)
for(i in 1:n.syn){
  #   write(t(cbind(death.synt[[i]],sub.synt[[i]][,2:4])),paste("synt",i,"_coords_zip.txt",sep=""),ncolumns=12)
  death.synt[[i]] = as.matrix(read.table(paste("synt",i,"_coords_zip.txt",sep="")))
  colnames(death.synt[[i]]) = c("ytil","sex","race","age2","age3","age4","edu2","edu3","edu4","lat","lon","zip")
}

## redo the synthetic zip code classification

require(class)

for(i in 1:n.syn){
  aux.zip = knn(death.orig[,10:11],death.synt[[i]][,10:11],death.orig[,12])
  aux.lev = levels(aux.zip)
  death.synt[[i]][,12] = as.numeric(aux.lev[as.numeric(aux.zip)])
}

###############################

## PROPORTION OF BLACK PEOPLE

list.zip = c(27701,27703,27704,27705,27707,27712,27713,27517,27522,27613)
n.zip = length(list.zip)

## proportion of black people ORIGINAL DATA

p.orig = sd.orig = numeric(n.zip); ci.orig = matrix(0,n.zip,2)

for(i in 1:n.zip){
  # mean, sd of the mean and confidence interval
  p.orig[i] = mean(death.orig[death.orig[,12]==list.zip[i],3])
  sd.orig[i] = sd(death.orig[death.orig[,12]==list.zip[i],3])/sqrt(length(death.orig[death.orig[,12]==list.zip[i],3]))
  ci.orig[i,] = c(p.orig[i] - 1.96*sd.orig[i],p.orig[i] + 1.96*sd.orig[i])
}
cbind(ci.orig,p.orig)

# table(death.orig[,12])
# table(death.synt[[1]][,12])
# table(death.synt[[10]][,12])

## proportion of black people SYNTHETIC DATA

p.synt = v.synt = matrix(0,n.zip,n.syn)
q.synt = sd.synt = numeric(n.zip); ci.synt = matrix(0,n.zip,2)

for(i in 1:n.zip){
  
  # obtain the estimates for each synthetic data set (mean and var of the mean)
  for(j in 1:n.syn){
    p.synt[i,j] = mean(death.synt[[j]][death.synt[[j]][,12]==list.zip[i],3])
    v.synt[i,j] = var(death.synt[[j]][death.synt[[j]][,12]==list.zip[i],3])/length(death.synt[[j]][death.synt[[j]][,12]==list.zip[i],3])
  }
  
  # combine the estimates
  q.synt[i] = mean(p.synt[i,])
  # components of the variance
  bm = var(p.synt[i,])
  vm = mean(v.synt[i,])
  Tp = bm/n.syn + vm
  sd.synt[i] = sqrt(Tp)
  # degrees of freedom for the t distribution of q
  rm = bm/(n.syn*vm)
  vp = (n.syn-1)*(1+1/rm)^2
  # confidence intervals
  ci.synt[i,] = c(q.synt[i] + qt(0.025,df=vp)*sqrt(Tp), q.synt[i] + qt(0.975,df=vp)*sqrt(Tp))
}

cbind(ci.synt,q.synt)

pblack.summary = cbind(p.orig,ci.orig,q.synt,ci.synt)
pblack.summary = pblack.summary[1:7,]
pblack.summary

## plotting ORIGINAL AND SYNTHETIC

plot(1:7-0.2,pblack.summary[,1],xlim=c(0.5,7.5),ylim=c(0,1),pch=19,xaxt="n",xlab="zip code",ylab="proportion of black people")
segments(x0=1:7-0.2,y0=pblack.summary[,2],x1=1:7-0.2,y1=pblack.summary[,3],lwd=3)
points(1:7+0.2,pblack.summary[,4],pch=19,col="blue")
segments(x0=1:7+0.2,y0=pblack.summary[,5],x1=1:7+0.2,y1=pblack.summary[,6],lwd=3,col="blue")
axis(side=1,at=1:7,labels=list.zip[1:7])

## plotting zipcodes together

# get the colors for each zip code
aux = col2rgb(c("black","red","green","coral","turquoise","blue","hotpink","grey","grey","grey"))
cores = rgb(aux[1,],aux[2,],aux[3,],max=255)

aux.cor = numeric(n)
for(i in 1:n){
  aux.cor[i] = cores[which(death.orig[i,12]==list.zip)]
}

# plot the map with everyone

pdf("zip_pblack.pdf", width=30, height=15, pointsize=28)

par(mfrow=c(1,2))

plot(death.orig[,10:11],col=aux.cor,pch=20,xlab="x",ylab="y",xlim=c(0,10),ylim=c(0,10),
     main=expression(paste("Observations based on zipcode")))
legend("topright",col=cores,leg=list.zip[1:7],pch=20)

plot(1:7-0.2,pblack.summary[,1],xlim=c(0.5,7.5),ylim=c(0,1),pch=20,xaxt="n",xlab="zip code",ylab="proportion",cex=2,
     col=cores,main=expression(paste("Estimated proportion and 95% CI of black people")))
segments(x0=1:7-0.2,y0=pblack.summary[,2],x1=1:7-0.2,y1=pblack.summary[,3],lwd=5,col=cores)
points(1:7+0.2,pblack.summary[,4],pch=8,col=cores,cex=1.5)
segments(x0=1:7+0.2,y0=pblack.summary[,5],x1=1:7+0.2,y1=pblack.summary[,6],lwd=3,col=cores)
axis(side=1,at=1:7,labels=list.zip[1:7])
legend("topright",col="black",pch=c(20,8),legend=c("original","synthetic"))

dev.off()

# plot with race=1 (black) highlighted

pdf("zip2_pblack.pdf", width=30, height=15, pointsize=28)

par(mfrow=c(1,2))

plot(death.orig[,10:11],col=aux.cor,pch=20,xlab="x",ylab="y",xlim=c(0,10),ylim=c(0,10),
     main=expression(paste("Observations based on zipcode")))
legend("topright",col=cores,leg=list.zip[1:7],pch=20)

# if race==1
points(death.orig[death.orig[,3]==1,10:11],col=rgb(100,100,100,alpha=150,max=255),pch=20)

plot(1:7-0.2,pblack.summary[,1],xlim=c(0.5,7.5),ylim=c(0,1),pch=20,xaxt="n",xlab="zip code",ylab="proportion",cex=2,
     col=cores,main=expression(paste("Estimated proportion and 95% CI of black people")))
segments(x0=1:7-0.2,y0=pblack.summary[,2],x1=1:7-0.2,y1=pblack.summary[,3],lwd=5,col=cores)
points(1:7+0.2,pblack.summary[,4],pch=8,col=cores,cex=1.5)
segments(x0=1:7+0.2,y0=pblack.summary[,5],x1=1:7+0.2,y1=pblack.summary[,6],lwd=3,col=cores)
axis(side=1,at=1:7,labels=list.zip[1:7])
legend("topright",col="black",pch=c(20,8),legend=c("original","synthetic"))

dev.off()

###### FOR THE PAPER

## plot the proportion of black by zipcode

pdf("pblack_BW.pdf", width=15, height=15, pointsize=28)

plot(1:7-0.2,pblack.summary[,1],xlim=c(0.5,7.5),ylim=c(0,1),pch=20,xaxt="n",xlab="zip code",ylab="proportion",cex=1.5,
     col="black",main=expression(paste("Estimated proportion and 95% CI of black people")))
segments(x0=1:7-0.2,y0=pblack.summary[,2],x1=1:7-0.2,y1=pblack.summary[,3],lwd=4,col="black")
points(1:7+0.2,pblack.summary[,4],pch=18,col="gray50",cex=1.5)
segments(x0=1:7+0.2,y0=pblack.summary[,5],x1=1:7+0.2,y1=pblack.summary[,6],lwd=4,col="gray50")
axis(side=1,at=1:7,labels=list.zip[1:7])
legend("topright",col=c("black","gray50"),pch=c(20,18),legend=c("original","synthetic"))

dev.off()

## table of the estimates and sd

tab.black = cbind(p.orig,sd.orig,q.synt,sd.synt)[1:7,]

xtable(tab.black,digits=3)


###############################

## PROPORTION OF YTIL

## proportion of ytil ORIGINAL DATA

p.orig = sd.orig = numeric(n.zip); ci.orig = matrix(0,n.zip,2)

for(i in 1:n.zip){
  # mean, sd of the mean and confidence interval
  p.orig[i] = mean(death.orig[death.orig[,12]==list.zip[i],1])
  sd.orig[i] = sd(death.orig[death.orig[,12]==list.zip[i],1])/sqrt(length(death.orig[death.orig[,12]==list.zip[i],1]))
  ci.orig[i,] = c(p.orig[i] - 1.96*sd.orig[i],p.orig[i] + 1.96*sd.orig[i])
}
cbind(ci.orig,p.orig)

## proportion of ytil SYNTHETIC DATA

p.synt = v.synt = matrix(0,n.zip,n.syn)
q.synt = sd.synt = numeric(n.zip); ci.synt = matrix(0,n.zip,2)

for(i in 1:n.zip){
  
  # obtain the estimates for each synthetic data set (mean and var of the mean)
  for(j in 1:n.syn){
    p.synt[i,j] = mean(death.synt[[j]][death.synt[[j]][,12]==list.zip[i],1])
    v.synt[i,j] = var(death.synt[[j]][death.synt[[j]][,12]==list.zip[i],1])/length(death.synt[[j]][death.synt[[j]][,12]==list.zip[i],1])
  }
  
  # combine the estimates
  q.synt[i] = mean(p.synt[i,])
  # components of the variance
  bm = var(p.synt[i,])
  vm = mean(v.synt[i,])
  Tp = bm/n.syn + vm
  sd.synt[i] = sqrt(Tp)
  # degrees of freedom for the t distribution of q
  rm = bm/(n.syn*vm)
  vp = (n.syn-1)*(1+1/rm)^2
  # confidence intervals
  ci.synt[i,] = c(q.synt[i] + qt(0.025,df=vp)*sqrt(Tp), q.synt[i] + qt(0.975,df=vp)*sqrt(Tp))
}

cbind(ci.synt,q.synt)
cbind(q.synt,sd.synt)

pblack.summary = cbind(p.orig,ci.orig,q.synt,ci.synt)
pblack.summary = pblack.summary[1:7,]
pblack.summary

## plotting ORIGINAL AND SYNTHETIC

plot(1:7-0.2,pblack.summary[,1],xlim=c(0.5,7.5),ylim=c(0,1),pch=19,xaxt="n",xlab="zip code",ylab="proportion of ytil")
segments(x0=1:7-0.2,y0=pblack.summary[,2],x1=1:7-0.2,y1=pblack.summary[,3],lwd=3)
points(1:7+0.2,pblack.summary[,4],pch=19,col="blue")
segments(x0=1:7+0.2,y0=pblack.summary[,5],x1=1:7+0.2,y1=pblack.summary[,6],lwd=3,col="blue")
axis(side=1,at=1:7,labels=list.zip[1:7])

###################################
## with the zip code map and colors

# get the colors for each zip code
aux = col2rgb(c("black","red","green","coral","turquoise","blue","hotpink","grey","grey","grey"))
cores = rgb(aux[1,],aux[2,],aux[3,],max=255)

aux.cor = numeric(n)
for(i in 1:n){
  aux.cor[i] = cores[which(death.orig[i,12]==list.zip)]
}

# plot the map with everyone

pdf("zip_pytil.pdf", width=30, height=15, pointsize=28)

par(mfrow=c(1,2))

plot(death.orig[,10:11],col=aux.cor,pch=20,xlab="x",ylab="y",xlim=c(0,10),ylim=c(0,10),
     main=expression(paste("Observations based on zipcode")))
legend("topright",col=cores,leg=list.zip[1:7],pch=20)

plot(1:7-0.2,pblack.summary[,1],xlim=c(0.5,7.5),ylim=c(0,1),pch=20,xaxt="n",xlab="zip code",ylab="proportion",cex=2,
     col=cores,main=expression(paste("Estimated proportion and 95% CI of ",tilde(y)==1)))
segments(x0=1:7-0.2,y0=pblack.summary[,2],x1=1:7-0.2,y1=pblack.summary[,3],lwd=5,col=cores)
points(1:7+0.2,pblack.summary[,4],pch=8,col=cores,cex=1.5)
segments(x0=1:7+0.2,y0=pblack.summary[,5],x1=1:7+0.2,y1=pblack.summary[,6],lwd=3,col=cores)
axis(side=1,at=1:7,labels=list.zip[1:7])
legend("topright",col="black",pch=c(20,8),legend=c("original","synthetic"))

dev.off()

# plot with ytil=1 highlighted

pdf("zip2_pytil.pdf", width=30, height=15, pointsize=28)

par(mfrow=c(1,2))

plot(death.orig[,10:11],col=aux.cor,pch=20,xlab="x",ylab="y",xlim=c(0,10),ylim=c(0,10),
     main=expression(paste("Observations based on zipcode")))
legend("topright",col=cores,leg=list.zip[1:7],pch=20)

# if ytil==1
points(death.orig[death.orig[,1]==1,10:11],col=rgb(100,100,100,alpha=150,max=255),pch=20)

plot(1:7-0.2,pblack.summary[,1],xlim=c(0.5,7.5),ylim=c(0,1),pch=20,xaxt="n",xlab="zip code",ylab="proportion",cex=2,
     col=cores,main=expression(paste("Estimated proportion and 95% CI of ",tilde(y)==1)))
segments(x0=1:7-0.2,y0=pblack.summary[,2],x1=1:7-0.2,y1=pblack.summary[,3],lwd=5,col=cores)
points(1:7+0.2,pblack.summary[,4],pch=8,col=cores,cex=1.5)
segments(x0=1:7+0.2,y0=pblack.summary[,5],x1=1:7+0.2,y1=pblack.summary[,6],lwd=3,col=cores)
axis(side=1,at=1:7,labels=list.zip[1:7])
legend("topright",col="black",pch=c(20,8),legend=c("original","synthetic"))

dev.off()

###### FOR THE PAPER

## plot the proportion of ytil by zipcode

pdf("pytil_BW.pdf", width=15, height=15, pointsize=28)

plot(1:7-0.2,pblack.summary[,1],xlim=c(0.5,7.5),ylim=c(0,1),pch=20,xaxt="n",xlab="zip code",ylab="proportion",cex=1.5,
     col="black",main=expression(paste("Estimated proportion and 95% CI of ",tilde(y)==1)))
segments(x0=1:7-0.2,y0=pblack.summary[,2],x1=1:7-0.2,y1=pblack.summary[,3],lwd=4,col="black")
points(1:7+0.2,pblack.summary[,4],pch=18,col="gray50",cex=1.5)
segments(x0=1:7+0.2,y0=pblack.summary[,5],x1=1:7+0.2,y1=pblack.summary[,6],lwd=4,col="gray50")
axis(side=1,at=1:7,labels=list.zip[1:7])
legend("topright",col=c("black","gray50"),pch=c(20,18),legend=c("original","synthetic"))

dev.off()

## table of the estimates and sd

tab.ytil = cbind(p.orig,sd.orig,q.synt,sd.synt)[1:7,]

tab.aux = cbind(tab.black,tab.ytil)
row.names(tab.aux) = list.zip[1:7]

xtable(tab.aux,digits=3)


list.zip


