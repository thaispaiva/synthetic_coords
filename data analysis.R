###################################################
###         Data Analysis - deaths in NC        ###
###################################################

setwd("~/Dropbox/PHD/Project - Jerry/Death data")

require(maptools)
require(spdep)
require(ars)

source("~/Dropbox/PHD/Project - Jerry/Gibbs sampler/mat2vec.R")
source("~/Dropbox/PHD/Project - Jerry/Gibbs sampler/neigh_list.R")
source("~/Dropbox/PHD/Project - Jerry/Gibbs sampler/gensyn.R")

## loading the workspace
load("/net/nfs1/s/grad/tvp/Dropbox/PHD/Project - Jerry/Death data/data analysis.RData")
# load("~/My Dropbox/PHD/Project - Jerry/Death data/dgDurham.Rdata")
attach(dg.durham2)

####################################################################################################

## loading the data
load("~/Dropbox/PHD/Project - Jerry/Death data/dgDurham.Rdata")
# load("~/My Dropbox/PHD/Project - Jerry/Death data/dgDurham.Rdata")
attach(dg.durham)

## looking at the variables
names(dg.durham)

table(SEX) # 1:male; 2:female
table(RACE) 
summary(AGE_UNITS)
hist(AGE_UNITS)

table(CO_OCC) # country of occurrence
table(CITY_OCC) # city of occurrence
table(CO_RES) # country of residence
table(CITY_RES) # city of residence

table(AUTOPSY) # 1:autopsy performed; 2:autopsy not performed; 9: unknown
table(FINDINGS) # 1:autopsy findings considered to determine cause of death; or 2:not considered; 9:unknown/no autopsy 

table(MARITAL_ST) # 1:never married; 2:married; 3:widowed; 4:divorced
table(ATTENDANT) # 1:physician; 2:medical examiner; 3:coroner

table(HISPANIC) # hispanic origin

summary(EDUCATION)
hist(EDUCATION) # 0-12:elementary/secondary ed; 13-16: 1-4 years college; 17: 5+ years college

table(ZIP) # zipcode

table(HOSP_STATU) # hospital status

class(RES_ADD) # residence address
length(levels(RES_ADD))
length(RES_ADD)

class(RES_CITY) # residence city
length(levels(RES_CITY))
table(RES_CITY)[table(RES_CITY)!=0] # all Durham

table(MANNER_DEA) # manner of death
# 1:accident; 2:suicide; 3:homicide; 4:pending; 5:could not determine; 6:self-inflicted; 7:natural

class(ICD10) # underlying cause of death code
table(ICD10)[table(ICD10)!=0]

class(X1_MCOD) # first mentioned cause of death
table(X1_MCOD)[table(X1_MCOD)!=0]

## ploting the coordinates
plot(dg.durham$POINT_X,dg.durham$POINT_Y)

plot(POINT_X,POINT_Y,pch=19,cex=0.5)

table(dg.durham$RES_CITY)
class(dg.durham$RES_CITY)


###################################################################################################

## plotting the different categories

sp = matrix(c(POINT_X,POINT_Y),length(POINT_X),2,byrow=F)
prop = (max(POINT_X) - min(POINT_X))/(max(POINT_Y) - min(POINT_Y))
## make them into the [0,10]x[0,10] region
sp[,1] = sp[,1]-(min(sp[,1])-0.1)
sp[,1] = 10*prop*sp[,1]/(max(sp[,1])+0.1)
sp[,2] = sp[,2]-(min(sp[,2])-0.1)
sp[,2] = 10*sp[,2]/(max(sp[,2])+0.1)
range(sp)
plot(sp,xlim=c(0,10),ylim=c(0,10))
# check new proportion
(max(sp[,1]) - min(sp[,1]))/(max(sp[,2]) - min(sp[,2]))
prop

# pdf("sex1.pdf", width=24, height=24, pointsize=24, paper="special")
# par(mfrow=c(2,2))
# for(i in 1:4){
# cond = SEX==1 & MARITAL_ST==i
# cor.ref = numeric(10); cor.ref[1] = rainbow(n=8,start=0.1)[1]; cor.ref[2:3] = c("black","red");
# cor.ref[4:10] = rainbow(n=8,start=0.1)[2:8]
# cor = cor.ref[RACE[cond]+1]
# pch.ref = c(15,20,17,rep(15,7))
# ch = pch.ref[RACE[cond]+1]
# plot(sp[cond,],col=cor,pch=ch,xlim=c(0,13),ylim=c(0,10),xlab="x",ylab="y")
# title(paste("Sex=",1,", Mar=",i))
# legend("topright",legend=0:9,pch=pch.ref,col=cor.ref,title="Race")
# }
# dev.off()
# 
# pdf("sex2.pdf", width=24, height=24, pointsize=24, paper="special")
# par(mfrow=c(2,2))
# for(i in 1:4){
# cond = SEX==2 & MARITAL_ST==i
# cor.ref = numeric(10); cor.ref[1] = rainbow(n=8,start=0.1)[1]; cor.ref[2:3] = c("black","red");
# cor.ref[4:10] = rainbow(n=8,start=0.1)[2:8]
# cor = cor.ref[RACE[cond]+1]
# pch.ref = c(15,20,17,rep(15,7))
# ch = pch.ref[RACE[cond]+1]
# plot(sp[cond,],col=cor,pch=ch,xlim=c(0,13),ylim=c(0,10),xlab="x",ylab="y")
# title(paste("Sex=",1,", Mar=",i))
# legend("topright",legend=0:9,pch=pch.ref,col=cor.ref,title="Race")
# }
# dev.off()

####################################################################################################

## selecting some observations

# removing the unknown age == 99
sum(AGE_UNITS == 99)/n  # 3% of observations

sum(RACE!=1 & RACE!=2)/n  # 8% of observations

sum(AGE_UNITS!=99 & (RACE==1 | RACE==2))/n  # include 98% of the original observations

dg.durham2 = subset(dg.durham,subset=(AGE_UNITS!=99 & (RACE==1 | RACE==2)),select=c(1,2,3,10,13,19,21,22))

detach(dg.durham)
attach(dg.durham2)
names(dg.durham2)

## creating the categorical variable for education
EDU = ifelse(EDUCATION<12,1,ifelse(EDUCATION==12,2,ifelse(EDUCATION>16,4,3)))
table(EDU)

## creating the categorical variable for age
AGE = ifelse(AGE_UNITS<60,1,ifelse(AGE_UNITS<75,2,ifelse(AGE_UNITS<85,3,4)))
table(AGE)

## looking at the cause of death
fir.ICD = substr(ICD10,start=1,stop=1)  # getting the first character from the death classification

## cancer and immune system
sum(fir.ICD == "C" | fir.ICD == "D")/length(fir.ICD)  # 24%
CANCER = fir.ICD == "C" | fir.ICD == "D"

plot(POINT_X,POINT_Y,pch=20,col=CANCER+1)

## circulatory
sum(fir.ICD == "I")/length(fir.ICD)  # 30%
CIRCULATORY = fir.ICD == "I"

plot(POINT_X,POINT_Y,pch=20,col=CIRCULATORY+1)

## respiratory
sum(fir.ICD == "J")/length(fir.ICD)  # 9.7%
RESP = fir.ICD == "J"

plot(POINT_X,POINT_Y,pch=20,col=RESP+1) ## no spatial relation

## external causes
sum(fir.ICD == "V" | fir.ICD == "W" | fir.ICD == "Y")/length(fir.ICD)  # 3%
EXTERNAL = fir.ICD == "V" | fir.ICD == "W" | fir.ICD == "Y"

plot(POINT_X,POINT_Y,pch=20,col=EXTERNAL+1) ## no spatial relation

## spatial coordinates
sp = matrix(c(POINT_X,POINT_Y),length(POINT_X),2,byrow=F)
prop = (max(POINT_X) - min(POINT_X))/(max(POINT_Y) - min(POINT_Y))
## make them into the [0,10]x[0,10] region
sp[,1] = sp[,1]-(min(sp[,1])-0.1)
sp[,1] = 10*prop*sp[,1]/(max(sp[,1])+0.1)
sp[,2] = sp[,2]-(min(sp[,2])-0.1)
sp[,2] = 10*sp[,2]/(max(sp[,2])+0.1)
range(sp)
plot(sp,xlim=c(0,10),ylim=c(0,10))
# check new proportion
(max(sp[,1]) - min(sp[,1]))/(max(sp[,2]) - min(sp[,2]))
prop

## plot the potential response variables again

setwd("~/Dropbox/PHD/Project - Jerry/Death data/figures eda")

pdf("cancer.pdf", width=12, height=12, pointsize=24, paper="special")
plot(sp[,1],sp[,2],pch=20,col=CANCER+1,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
title("Cancer and immune system")
legend("topright",col=c(2,1),pch=20,leg=c(paste("cancer (",sum(CANCER),")"),paste("others (",sum(!CANCER),")")))
dev.off()

pdf("sub_cancer.pdf", width=24, height=12, pointsize=24, paper="special")
par(mfrow=c(1,2))
plot(sp[CANCER,1],sp[CANCER,2],pch=20,col=2,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
title("Cancer and immune system - cancer")
legend("bottomright",leg=paste(sum(CANCER)),pch=20,col=2)
plot(sp[!CANCER,1],sp[!CANCER,2],pch=20,col=1,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
title("Cancer and immune system - others")
legend("bottomright",leg=paste(sum(!CANCER)),pch=20,col=1)
dev.off()

pdf("circulatory.pdf", width=12, height=12, pointsize=24, paper="special")
plot(sp[,1],sp[,2],pch=20,col=CIRCULATORY+1,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
title("Circulatory system")
legend("topright",col=c(2,1),pch=20,leg=c(paste("circulatory (",sum(CIRCULATORY),")"),paste("others (",sum(!CIRCULATORY),")")))
dev.off()

pdf("sub_circulatory.pdf", width=24, height=12, pointsize=24, paper="special")
par(mfrow=c(1,2))
plot(sp[CIRCULATORY,1],sp[CIRCULATORY,2],pch=20,col=2,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
title("Circulatory system - circulatory")
legend("bottomright",leg=paste(sum(CIRCULATORY)),pch=20,col=2)
plot(sp[!CIRCULATORY,1],sp[!CIRCULATORY,2],pch=20,col=1,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
title("Circulatory system - others")
legend("bottomright",leg=paste(sum(!CIRCULATORY)),pch=20,col=1)
dev.off()


## plot the other variables

pdf("sex.pdf", width=12, height=12, pointsize=24, paper="special")
plot(sp[,1],sp[,2],pch=20,col=SEX,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
title("Points plotted by sex")
legend("topright",col=c(1,2),pch=20,leg=c(paste("male (",sum(SEX==1),")"),paste("female (",sum(SEX==2),")")))
dev.off()

pdf("sub_sex.pdf", width=24, height=12, pointsize=24, paper="special")
par(mfrow=c(1,2))
plot(sp[SEX==1,1],sp[SEX==1,2],pch=20,col=1,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
title("Points plotted by sex - male")
legend("bottomright",col=1,pch=20,leg=paste(sum(SEX==1)))
plot(sp[SEX==2,1],sp[SEX==2,2],pch=20,col=2,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
title("Points plotted by sex - female")
legend("bottomright",col=2,pch=20,leg=paste(sum(SEX==2)))
dev.off()

pdf("race.pdf", width=12, height=12, pointsize=24, paper="special")
plot(sp[,1],sp[,2],pch=20,col=RACE,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
title("Points plotted by race")
legend("topright",col=c(1,2),pch=20,leg=c(paste("white (",sum(RACE==1),")"),paste("black (",sum(RACE==2),")")))
dev.off()

pdf("sub_race.pdf", width=24, height=12, pointsize=24, paper="special")
par(mfrow=c(1,2))
plot(sp[RACE==1,1],sp[RACE==1,2],pch=20,col=1,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
title("Points plotted by race - white")
legend("bottomright",col=1,pch=20,leg=paste(sum(RACE==1)))
plot(sp[RACE==2,1],sp[RACE==2,2],pch=20,col=2,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
title("Points plotted by race - black")
legend("bottomright",col=2,pch=20,leg=paste(sum(RACE==2)))
dev.off()

pdf("age.pdf", width=12, height=12, pointsize=24, paper="special")
plot(sp[,1],sp[,2],pch=20,col=AGE,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
title("Points plotted by age")
legend("topright",col=c(1,2,3,4),pch=20,leg=c(paste("<60 (",sum(AGE==1),")"),
                                              paste("60-75 (",sum(AGE==2),")"),
                                              paste("75-85 (",sum(AGE==3),")"),
                                              paste(">85 (",sum(AGE==4),")")))
dev.off()

pdf("sub_age.pdf", width=24, height=24, pointsize=24, paper="special")
par(mfrow=c(2,2))
for(i in 1:4){
  plot(sp[AGE==i,1],sp[AGE==i,2],pch=20,col=i,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
  title(paste("Age",i))
  legend("bottomright",col=i,pch=20,leg=paste(sum(AGE==i)))
}
dev.off()

pdf("edu.pdf", width=12, height=12, pointsize=24, paper="special")
plot(sp[,1],sp[,2],pch=20,col=EDU,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
title("Points plotted by education")
legend("topright",col=c(1,2,3,4),pch=20,leg=c(paste("<12 (",sum(EDU==1),")"),
                                              paste("12 (",sum(EDU==2),")"),
                                              paste("13-16 (",sum(EDU==3),")"),
                                              paste(">16 (",sum(EDU==4),")"))) 
dev.off()

pdf("sub_edu.pdf", width=24, height=24, pointsize=24, paper="special")
par(mfrow=c(2,2))
for(i in 1:4){
  plot(sp[EDU==i,1],sp[EDU==i,2],pch=20,col=i,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
  title(paste("Edu",i))
  legend("bottomright",col=i,pch=20,leg=paste(sum(EDU==i)))
}
dev.off()

pdf("marital.pdf", width=12, height=12, pointsize=24, paper="special")
plot(sp[,1],sp[,2],pch=20,col=MARITAL_ST,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
title("Points plotted by marital status")
legend("topright",col=c(1,2,3,4),pch=20,leg=c(paste("never married (",sum(MARITAL_ST==1),")"),
                                              paste("married (",sum(MARITAL_ST==2),")"),
                                              paste("widowed (",sum(MARITAL_ST==3),")"),
                                              paste("divorced (",sum(MARITAL_ST==4),")")))
dev.off()

pdf("sub_marital.pdf", width=24, height=24, pointsize=24, paper="special")
par(mfrow=c(2,2))
for(i in 1:4){
  plot(sp[MARITAL_ST==i,1],sp[MARITAL_ST==i,2],pch=20,col=i,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
  title(paste("Marital status",i))
  legend("bottomright",col=i,pch=20,leg=paste(sum(MARITAL_ST==i)))
}
dev.off()

## investigating some interactions

## cancer/circulatory vs age
pdf("cancer_age.pdf", width=24, height=24, pointsize=24, paper="special")
par(mfrow=c(2,2))
for(i in 1:4){
  plot(sp[AGE==i,1],sp[AGE==i,2],pch=20,col=CANCER[AGE==i]+1,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
  title(paste("Age ",i))
  legend("topright",col=c(2,1),pch=20,leg=c("cancer","others"))
}
dev.off()

pdf("circulatory_age.pdf", width=24, height=24, pointsize=24, paper="special")
par(mfrow=c(2,2))
for(i in 1:4){
  plot(sp[AGE==i,1],sp[AGE==i,2],pch=20,col=CIRCULATORY[AGE==i]+1,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
  title(paste("Age ",i))
  legend("topright",col=c(2,1),pch=20,leg=c("circulatory","others"))
}
dev.off()

## cancer/circulatory vs education
pdf("cancer_edu.pdf", width=24, height=24, pointsize=24, paper="special")
par(mfrow=c(2,2))
for(i in 1:4){
  plot(sp[EDU==i,1],sp[EDU==i,2],pch=20,col=CANCER[EDU==i]+1,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
  title(paste("Education ",i))
  legend("topright",col=c(2,1),pch=20,leg=c("cancer","others"))
}
dev.off()

pdf("circulatory_edu.pdf", width=24, height=24, pointsize=24, paper="special")
par(mfrow=c(2,2))
for(i in 1:4){
  plot(sp[EDU==i,1],sp[EDU==i,2],pch=20,col=CIRCULATORY[EDU==i]+1,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
  title(paste("Education ",i))
  legend("topright",col=c(2,1),pch=20,leg=c("circulatory","others"))
}
dev.off()

## cancer/circulatory vs marital status
pdf("cancer_marital.pdf", width=24, height=24, pointsize=24, paper="special")
par(mfrow=c(2,2))
for(i in 1:4){
  plot(sp[MARITAL_ST==i,1],sp[MARITAL_ST==i,2],pch=20,col=CANCER[MARITAL_ST==i]+1,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
  title(paste("Marital Status ",i))
  legend("topright",col=c(2,1),pch=20,leg=c("cancer","others"))
}
dev.off()

pdf("circulatory_marital.pdf", width=24, height=24, pointsize=24, paper="special")
par(mfrow=c(2,2))
for(i in 1:4){
  plot(sp[MARITAL_ST==i,1],sp[MARITAL_ST==i,2],pch=20,col=CIRCULATORY[MARITAL_ST==i]+1,xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
  title(paste("Marital Status ",i))
  legend("topright",col=c(2,1),pch=20,leg=c("circulatory","others"))
}
dev.off()

## age vs marital status
pdf("age_marital.pdf", width=24, height=24, pointsize=24, paper="special")
par(mfrow=c(2,2))
for(i in 1:4){
  plot(sp[AGE==i,1],sp[AGE==i,2],pch=20,col=MARITAL_ST[AGE==i],xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
  title(paste("Age ",i))
  legend("topright",col=1:4,pch=20,leg=c("never married","married","widowed","divorced"))
}
dev.off()

## age vs education
pdf("age_education.pdf", width=24, height=24, pointsize=24, paper="special")
par(mfrow=c(2,2))
for(i in 1:4){
  plot(sp[AGE==i,1],sp[AGE==i,2],pch=20,col=EDU[AGE==i],xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
  title(paste("Age ",i))
  legend("topright",col=1:4,pch=20,leg=c("<12","12","13-16",">16"))
}
dev.off()

## sex vs education
pdf("sex_education.pdf", width=24, height=12, pointsize=24, paper="special")
par(mfrow=c(1,2))
for(i in 1:2){
  plot(sp[SEX==i,1],sp[SEX==i,2],pch=20,col=EDU[SEX==i],xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
  title(paste("Sex ",i))
  legend("topright",col=1:4,pch=20,leg=c("<12","12","13-16",">16"))
}
dev.off()

## race vs education
pdf("race_education.pdf", width=24, height=12, pointsize=24, paper="special")
par(mfrow=c(1,2))
for(i in 1:2){
  plot(sp[RACE==i,1],sp[RACE==i,2],pch=20,col=EDU[RACE==i],xlim=c(0,10),ylim=c(0,10),xlab="x",ylab="y")
  title(paste("Race ",i))
  legend("topright",col=1:4,pch=20,leg=c("<12","12","13-16",">16"))
}
dev.off()

####################################################################################################

setwd("~/Dropbox/PHD/Project - Jerry/Death data")

## first, let us select some variables

names(dg.durham2)
X = matrix(c(CANCER,SEX,RACE,AGE,EDU,MARITAL_ST),length(SEX),6,byrow=F)
# X = matrix(c(CIRCULATORY,SEX,RACE,AGE,EDU,MARITAL_ST),length(SEX),6,byrow=F)

############## other X
X = X[,1:5]

## create the variables with the dimensions
n = dim(X)[1]; p = dim(X)[2]

vx = list(1); vx = c(vx,2:p)  # list with the unique values of each x
for(i in 1:p){
  vx[[i]] = sort(unique(X[,i]))
}

nx = numeric(p)  # number of values each x can assume
for(i in 1:p){
  nx[i] = length(vx[[i]])
}

nb = prod(nx)
b = 1:nb

size = 21
nm = (size-1)*(size-1)
m = 1:nm

## create the grid points
xvec = yvec = seq(0,10,length.out=size)

## creating matrix with the combinations
b.matrix = matrix(0,nb,p+1)
b.matrix[,1] = b
b.matrix[,2] = rep(vx[[1]],times=1,each=prod(nx[2:p]))
for(i in 2:p){
  b.matrix[,i+1] = rep(vx[[i]],times=prod(nx[1:(i-1)]),each=prod(nx[(i+1):p]))
}

## identifying the combination index for all observations
comb = numeric(n)
for(i in 1:n){
  for(j in 1:nb){
    if(all(X[i,]==b.matrix[j,2:(p+1)])){
      comb[i] = j
    }
  }
}

comb = factor(comb,levels=b)
hist(as.numeric(table(comb)))
summary(as.numeric(table(comb)))

n.comb = table(comb)
length(n.comb)
sum(n.comb==0)/nb
sum(n.comb<5)/nb
sum(n.comb<2)/nb

## matrix with the counts of observed points at each gridcell
nib = matrix(0,nb,nm)
for(j in 1:nm){
  aux = vec2mat(j,size)
  row = aux[1,1]
  col = aux[1,2]
  
  for(i in 1:nb){
    nib[i,j] = sum(comb==i & (xvec[col]<sp[,1] & sp[,1]<xvec[col+1]) &
      (yvec[row]<sp[,2] & sp[,2]<yvec[row+1]) )
  }
}
sum(nib)
apply(nib,MAR=1,FUN="sum")
counts = nib

## plot some of the combinations

# top 4
sort(n.comb,dec=T)[1:4]
# top4 = c(179,183,187,243)  # combinations with the highest frequencies
top4 = c(45,46,47,42)

par(mfrow=c(2,2))
for(i in 1:4){
  image(xvec,yvec,matrix(counts[top4[i],],size-1,size-1,byrow=F),col=rev(heat.colors(30)))
  points(jitter(sp[comb==top4[i],1]),jitter(sp[comb==top4[i],2]),pch=19,cex=0.5)
}
dev.off()

# ## looking at something weird... 12 repeated (?) measurements -> no! they're nursing homes
# aux = which(sp[,2]>8.5 & sp[,1]>5 & comb==top4[i])
# sp[aux,]
# dg.durham2[aux,]
# 
# aux = which(sp[,1]>4 & sp[,1]<4.5 & sp[,2]>6.5 & sp[,2]<7 & comb==top4[i])
# sp[aux,]
# dg.durham2[aux,]

####################################################################################################

## MCMC
## include only the main effects

## create neighborhood list
#neigh = neigh.list(nm,size)
neigh = cell2nb((size-1),(size-1),type="queen")
ni = card(neigh)
W = nb2mat(neigh,style="B")
sum(apply(W,MAR=2,FUN=sum) == ni) == (size-1)^2

## creating the index aux's for the parameters alpha, beta, phi

# list of the elements in alpha used in each combination
ind.a = list(1); ind.a = c(ind.a,2:nb)
for(i in b){
  ind.a[[i]] = numeric(0)
  j=1
  if(any(b.matrix[i,j+1]==vx[[j]][-1])){
    ind.a[[i]] = c(ind.a[[i]],which(b.matrix[i,j+1]==vx[[j]][-1]))
  }
  for(j in 2:p){
    if(any(b.matrix[i,j+1]==vx[[j]][-1])){
      ind.a[[i]] = c(ind.a[[i]],sum(nx[1:(j-1)]-1)+which(b.matrix[i,j+1]==vx[[j]][-1]))
    }
  }
}

# list with the combinations that have each alpha
sub.a = list(1); sub.a = c(sub.a,2:sum(nx-1))
for(i in 1:sum(nx-1)){
  sub.a[[i]] = numeric(0)
  for(j in b){
    if(i %in% ind.a[[j]])
      sub.a[[i]] = c(sub.a[[i]],j)
  }
}

## save workspace to use in the MCMC
save.image("/net/nfs1/s/grad/tvp/Dropbox/PHD/Project - Jerry/Death data/data analysis.RData")


###################################################################################################################

## no interaction terms --> no beta's

# # list of the elements in beta used in each combination
# ind.b = list(1); ind.b = c(ind.b,2:nb)
# for(i in b){
#   ind.b[[i]] = numeric(0)
#   if(sum(b.matrix[i,2]==vy[-1])){  
#     if(sum(b.matrix[i,3]==vx1[-1])){ 
#       ind.b[[i]] = c(ind.b[[i]],(which(b.matrix[i,2]==vy[-1])-1)*nx1+which(b.matrix[i,3]==vx1[-1])) 
#     }
#     if(sum(b.matrix[i,4]==vx2[-1])){  
#       ind.b[[i]] = c(ind.b[[i]],(ny-1)*(nx1-1)+which(b.matrix[i,4]==vx2[-1]))
#     }
#   }
#   if(sum(b.matrix[i,3]==vx1[-1])){
#     if(sum(b.matrix[i,4]==vx2[-1])){
#       ind.b[[i]] = c(ind.b[[i]],(ny-1)*(nx1-1+nx2-1)+(which(b.matrix[i,3]==vx1[-1])-1)*nx2+which(b.matrix[i,4]==vx2[-1]))
#     }
#   }
# }
# 
# # list with the combinations that have each beta
# sub.b = list(1); sub.b = c(sub.b,2:((ny-1)*(nx1-1)+(ny-1)*(nx2-1)+(nx1-1)*(nx2-1)))
# for(i in 1:length(sub.b)){
#   sub.b[[i]] = numeric(0)
#   for(j in 1:length(ind.b)){
#     if(sum(i==ind.b[[j]]))
#       sub.b[[i]] = c(sub.b[[i]],j);
#   }
# }

#############################
## generate initial values ##
#############################

# sim = 20000
# burn = 10000
# n.syn = 10; int.syn = 1000
# list.syn = seq(burn+1,burn+(n.syn-1)*int.syn+1,by=int.syn) # indexes of the simulations to save
# syn.data = list(1)
# syn.data = c(syn.data,2:n.syn)    # list of the synthetic datasets
# 
# v = 5  # prior variance for the diffuse normals
# a.tau = b.tau = 0.1  # prior parameters for tau
# m.bar = mean(ni)
# b.s2 = 0.1; a.s2 = m.bar*(0.7^2)*b.s2  # equation (5.48) from the spatial book
# 
# # mu - overall intercept - prior: normal(0,v)
# mu = numeric(sim)
# # mu[1] = rnorm(1,0,sd=sqrt(v))
# mu[1] = 1
# 
# # alpha - main effects - diffuse normal prior
# alpha = matrix(0,sum(nx-1),sim)
# # alpha[,1] = rnorm(ny-1+nx1-1+nx2-1,0,sd=v)
# alpha[,1] = 0
# 
# # # beta - interactions - diffuse normal prior
# # beta = matrix(0,(ny-1)*(nx1-1)+(ny-1)*(nx2-1)+(nx1-1)*(nx2-1),sim)
# # # beta[,1] = rnorm((ny-1)*(nx1-1)+(ny-1)*(nx2-1)+(nx1-1)*(nx2-1),0,sd=v)
# # beta[,1] = 0
# 
# # theta_i - overall car model 
# theta = matrix(0,nm,sim)
# # make it sum to zero
# theta[,1] = theta[,1] - sum(theta[,1])/nm
# cur.theta = theta[,1]
# 
# # phi - main effects' car model
# phi = array(data=0,dim=c(nm,sim,sum(nx-1)))
# # make it sum to zero
# for(i in 1:dim(phi)[3]){
#   phi[,1,i] = phi[,1,i] - mean(phi[,1,i])
# }
# #apply(phi,MAR=3,FUN=sum)
# cur.phi = phi[,1,]
# 
# # tau - precision of the car model - gamma prior
# tau = numeric(sim)
# # tau[1] = rgamma(1,a.tau,rate=b.tau)
# tau[1] = 1
# 
# # tau_phi - precision of the other car models - gamma prior
# tau.phi = matrix(0,dim(phi)[3],sim)
# tau.phi[,1] = 1
# 
# # inv.s2 precision of the random errors - gamma prior
# inv.s2 = numeric(sim)
# # inv.s2[1] = rgamma(1,a.s2,rate=b.s2)
# inv.s2[1] = 1
# 
# # epsilon_i - normal random errors
# epsilon = array(data=0,dim=c(nm,sim,nb))
# 
# # gamma
# gamma = array(data=0,dim=c(nm,sim,nb))
# for(i in b){
#   if(length(ind.a[[i]])==0){
# #     if(length(ind.b[[i]])==0){
#       gamma[,1,i] = log(n) + mu[1] + theta[,1] + epsilon[,1,i]
# #     }
# #     else{
# #       gamma[,1,i] = log(n) + mu[1] + sum(beta[ind.b[[i]],1]) + theta[,1] + epsilon[,1,i]
# #     }
#   }
#   else{
# #     if(length(ind.b[[i]])==0){
#       gamma[,1,i] = log(n) + mu[1] + sum(alpha[ind.a[[i]],1]) + theta[,1] + ifelse(length(ind.a[[i]])>1,apply(phi[,1,ind.a[[i]]],MAR=1,FUN=sum),phi[,1,ind.a[[i]]]) + epsilon[,1,i]
# #     }
# #     else{
# #       gamma[,1,i] = log(n) + mu[1] + sum(alpha[ind.a[[i]],1]) + sum(beta[ind.b[[i]],1]) + theta[,1] + ifelse(length(ind.a[[i]])>1,apply(phi[,1,ind.a[[i]]],MAR=1,FUN=sum),phi[,1,ind.a[[i]]]) + epsilon[,1,i]
# #     }
#   }
# }
# cur.gamma = gamma[,1,]
# 
# ###########################
# ## start the simulations ##
# ###########################
# 
# ## functions for the adaptive rejection sampling steps
# 
# # distribution of mu
# f.mu = function(x,n.f,v.f,sumeta.f){
#   x*n.f - (1/(2*v.f))*x^2 - exp(x)*sumeta.f
# }
# 
# fprima.mu = function(x,n.f,v.f,sumeta.f){
#   n.f - x/v.f - exp(x)*sumeta.f
# }
# 
# # distribution of theta_i
# f.the = function(x,c.f,sumeta.f,ni.f,tau.f,bar.f){
#   x*c.f - exp(x)*sumeta.f - (1/2)*ni.f*tau.f*(x-bar.f)^2
# }
# 
# fprima.the = function(x,c.f,sumeta.f,ni.f,tau.f,bar.f){
#   c.f - exp(x)*sumeta.f - ni.f*tau.f*(x-bar.f)
# }
# 
# # distribution of epsilon_i
# f.eps = function(x,c.f,sumeta.f,tau.f){
#   x*c.f - exp(x)*sumeta.f - (1/2)*tau.f*x^2
# }
# 
# fprima.eps = function(x,c.f,sumeta.f,tau.f){
#   c.f - exp(x)*sumeta.f - tau.f*x
# }
# 
# 
# for(s in 2:sim){
#   
#   cat("Iteration ",s,"\n")
#   
#   ## update mu
#   eta = cur.gamma - mu[s-1]
#   sumeta = sum(exp(eta))
#   mu[s] = ars(1,f.mu,fprima.mu,lb=T,xlb=-100,ub=T,xub=100,n.f=n,v.f=v,sumeta.f=sumeta)
#   #   mu[s] = ars(1,f.mu,fprima.mu,lb=T,xlb=mu[s-1]-10,ub=T,xub=mu[s-1]+10,n.f=n,v.f=v,sumeta.f=sumeta)
#   #   mu[s] = ars(1,f.mu,fprima.mu,lb=T,xlb=-70,ub=T,xub=70,n.f=n,v.f=v,sumeta.f=sumeta)
#   cur.gamma = eta + mu[s]
#   
#   ## update alpha
#   for(i in 1:dim(alpha)[1]){
#     n.aux = sum(nib[sub.a[[i]],])
#     eta = cur.gamma[,sub.a[[i]]] - alpha[i,(s-1)]
#     sum.aux = sum(exp(eta))
#     alpha[i,s] = ars(1,f.mu,fprima.mu,lb=T,xlb=-100,ub=T,xub=100,n.f=n.aux,v.f=v,sumeta.f=sum.aux)
#     #     alpha[i,s] = ars(1,f.mu,fprima.mu,lb=T,xlb=alpha[i,(s-1)]-10,ub=T,xub=alpha[i,(s-1)]+10,n.f=n.aux,v.f=v,sumeta.f=sum.aux)
#     #     alpha[i,s] = ars(1,f.mu,fprima.mu,lb=T,xlb=-70,ub=T,xub=70,n.f=n.aux,v.f=v,sumeta.f=sum.aux)
#     cur.gamma[,sub.a[[i]]] = eta + alpha[i,s]
#   }
#   
# #   ## update beta
# #   for(i in 1:dim(beta)[1]){
# #     n.aux = sum(nib[sub.b[[i]],])
# #     eta = cur.gamma[,sub.b[[i]]] - beta[i,(s-1)]
# #     sum.aux = sum(exp(eta))
# #     beta[i,s] = ars(1,f.mu,fprima.mu,lb=T,xlb=-100,ub=T,xub=100,n.f=n.aux,v.f=v,sumeta.f=sum.aux)
# #     #     beta[i,s] = ars(1,f.mu,fprima.mu,lb=T,xlb=beta[i,(s-1)]-10,ub=T,xub=beta[i,(s-1)]+10,n.f=n.aux,v.f=v,sumeta.f=sum.aux)
# #     #     beta[i,s] = ars(1,f.mu,fprima.mu,lb=T,xlb=-70,ub=T,xub=70,n.f=n.aux,v.f=v,sumeta.f=sum.aux)
# #     cur.gamma[,sub.b[[i]]] = eta + beta[i,s]
# #   }
#   
#   ## update theta
#   c.aux = apply(counts,MAR=2,FUN=sum)
#   eta = cur.gamma - theta[,(s-1)] # since the sumeta.f does not depend on the updated theta's
#   for(i in m){
#     sum.aux = sum(exp(eta[i,]))
#     bar = W[i,]%*%cur.theta/ni[i]
#     cur.theta[i] = ars(1,f.the,fprima.the,lb=T,xlb=-100,ub=T,xub=100,c.f=c.aux[i],sumeta.f=sum.aux,ni.f=ni[i],tau.f=tau[s-1],bar.f=bar)
#     #     cur.theta[i] = ars(1,f.the,fprima.the,lb=T,xlb=theta[i,(s-1)]-10,ub=T,xub=theta[i,(s-1)]+10,c.f=c.aux[i],sumeta.f=sum.aux,ni.f=ni[i],tau.f=tau[s-1],bar.f=bar)
#     #     cur.theta[i] = ars(1,f.the,fprima.the,lb=T,xlb=-70,ub=T,xub=70,c.f=c.aux[i],sumeta.f=sum.aux,ni.f=ni[i],tau.f=tau[s-1],bar.f=bar)
#   }
#   ## make it sum to zero
#   theta[,s] = cur.theta - sum(cur.theta)/nm
#   cur.gamma = eta + theta[,s]
#   
#   ## update phi
#   for(k in 1:dim(phi)[3]){
#     c.aux = apply(counts[sub.a[[k]],],MAR=2,FUN=sum)
#     eta = cur.gamma[,sub.a[[k]]] - phi[,(s-1),k]
#     for(i in m){
#       sum.aux = sum(exp(eta[i,]))
#       bar = W[i,]%*%cur.phi[,k]/ni[i]
#       cur.phi[i,k] = ars(1,f.the,fprima.the,lb=T,xlb=-100,ub=T,xub=100,c.f=c.aux[i],sumeta.f=sum.aux,ni.f=ni[i],tau.f=tau.phi[k,(s-1)],bar.f=bar)
#       #       cur.phi[i,k] = ars(1,f.the,fprima.the,lb=T,xlb=phi[i,(s-1),k]-10,ub=T,xub=phi[i,(s-1),k]+10,c.f=c.aux[i],sumeta.f=sum.aux,ni.f=ni[i],tau.f=tau.phi[k,(s-1)],bar.f=bar)
#       #       cur.phi[i,k] = ars(1,f.the,fprima.the,lb=T,xlb=-70,ub=T,xub=70,c.f=c.aux[i],sumeta.f=sum.aux,ni.f=ni[i],tau.f=tau.phi[k,(s-1)],bar.f=bar)
#     }
#     ## make it sum to zero
#     phi[,s,k] = cur.phi[,k] - sum(cur.phi[,k])/nm
#     cur.gamma[,sub.a[[k]]] = eta + phi[,s,k]
#   }
#   
#   ## update epsilon
#   for(i in m){
#     for(j in b){
#       eta = cur.gamma[i,j] - epsilon[i,(s-1),j]
#       sum.aux = exp(eta)
#       epsilon[i,s,j] = ars(1,f.eps,fprima.eps,lb=T,xlb=-100,ub=T,xub=100,c.f=counts[j,i],sumeta.f=sum.aux,tau.f=inv.s2[s-1])
#       #       epsilon[i,s,j] = ars(1,f.eps,fprima.eps,lb=T,xlb=epsilon[i,(s-1),j]-10,ub=T,xub=epsilon[i,(s-1),j]+10,c.f=counts[j,i],sumeta.f=sum.aux,tau.f=inv.s2[s-1])
#       #       epsilon[i,s,j] = ars(1,f.eps,fprima.eps,lb=T,xlb=-70,ub=T,xub=70,c.f=counts[j,i],sumeta.f=sum.aux,tau.f=inv.s2[s-1])
#       cur.gamma[i,j] = eta + epsilon[i,s,j]
#     }
#   }
#   
#   ## update gamma
#   gamma[,s,] = cur.gamma
#   
#   ## update tau
#   sum.theta = t(theta[,s])%*%(diag(ni)-W)%*%(theta[,s])
#   tau[s] = rgamma(1,a.tau+nm/2,rate=b.tau+sum.theta/2)  
#   
#   ## update tau.phi
#   for(k in 1:dim(tau.phi)[1]){
#     sum.phi = t(phi[,s,k])%*%(diag(ni)-W)%*%(phi[,s,k])
#     tau.phi[k,s] = rgamma(1,a.tau+nm/2,rate=b.tau+sum.phi/2)
#   }
#   
#   ## update inv.s2
#   sum.s2 = sum(epsilon[,s,]^2)
#   inv.s2[s] = rgamma(1,a.s2+(nm*nb)/2,rate=b.s2+sum.s2/2)
#   
#   #   par(mfrow=c(2,2))
#   #   aux = exp(gamma[,s])
#   #   image.plot(xvec,yvec,matrix(aux,size-1,size-1),col=rev(heat.colors(10)),useRaster=T)
#   #   points(sp,pch=19,cex=0.5)
#   #   ts.plot(mu[1:s],main="mu")
#   #   ts.plot(inv.s2[1:s],main="inv.s2")
#   #   ts.plot(tau[1:s],main="tau")
#   
#   ## generate synthetic dataset
#   
#   if(sum(s==list.syn)==1){
#     i.syn = which(s==list.syn)
#     gamma.s = matrix(gamma[,s,],nm,nb)
#     
#     ## create the synthetic dataset and save it
#     syn.data[[i.syn]] = gen.syn(gamma.s,nb,size,n.comb,xvec,yvec)
#     write(t(syn.data[[i.syn]]),file=paste("syndata",i.syn,".txt"),ncolumns=3)
#     
#     ## plot the current lambda and the synthetic dataset
#     postscript(paste("syn",i.syn,".eps"), horizontal=F, width=12, height=12, pointsize=24, paper="special")
#     par(mfrow=c(4,3))
#     for(i in b){
#       image.plot(xvec,yvec,matrix(exp(gamma.s[,i]),size-1,size-1),col=rev(heat.colors(10)),useRaster=T)
#       points(syn.data[[i.syn]][which(syn.data[[i.syn]][,1]==i),2:3],pch=19,cex=0.5)
#       title(paste("Combination ",i))
#     }
#     dev.off()
#   }
#   
# }

