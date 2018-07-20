#### Thais Paiva - jul/2018

##################################################################################

## load packages
require(maptools)
require(spdep)
require(ars)
require(fields)

## load data
load("data analysis.RData")

## change some of the variables from the previous workspace
size = 21
nm = (size-1)*(size-1)
m = 1:nm

## create the grid points
xvec = yvec = seq(0,10,length.out=size)

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

## create neighborhood list
#neigh = neigh.list(nm,size)
neigh = cell2nb((size-1),(size-1),type="queen")
ni = card(neigh)
W = nb2mat(neigh,style="B")
sum(apply(W,MAR=2,FUN=sum) == ni) == (size-1)^2


#############################
## generate initial values ##
#############################

sim = 5001
burn = 1000
n.syn = 5; int.syn = 1000
list.syn = seq(burn+1,burn+(n.syn-1)*int.syn+1,by=int.syn) # indexes of the simulations to save
syn.data = list(1)
syn.data = c(syn.data,2:n.syn)    # list of the synthetic datasets

v = 5  # prior variance for the diffuse normals
a.tau = b.tau = 0.1  # prior parameters for tau
m.bar = mean(ni)
b.s2 = 0.1; a.s2 = m.bar*(0.7^2)*b.s2  # equation (5.48) from the spatial book

# mu - overall intercept - prior: normal(0,v)
mu = numeric(sim)
# mu[1] = rnorm(1,0,sd=sqrt(v))
mu[1] = 1

# alpha - main effects - diffuse normal prior
alpha = matrix(0,sum(nx-1),sim)
# alpha[,1] = rnorm(ny-1+nx1-1+nx2-1,0,sd=v)
alpha[,1] = 0

# theta_i - overall car model 
theta = matrix(0,nm,sim)
# make it sum to zero
theta[,1] = theta[,1] - sum(theta[,1])/nm
cur.theta = theta[,1]

# phi - main effects' car model
phi = array(data=0,dim=c(nm,sim,sum(nx-1)))
# make it sum to zero
for(i in 1:dim(phi)[3]){
  phi[,1,i] = phi[,1,i] - mean(phi[,1,i])
}
#apply(phi,MAR=3,FUN=sum)
cur.phi = phi[,1,]

# tau - precision of the car model - gamma prior
tau = numeric(sim)
# tau[1] = rgamma(1,a.tau,rate=b.tau)
tau[1] = 1

# tau_phi - precision of the other car models - gamma prior
tau.phi = matrix(0,dim(phi)[3],sim)
tau.phi[,1] = 1

# inv.s2 precision of the random errors - gamma prior
inv.s2 = numeric(sim)
# inv.s2[1] = rgamma(1,a.s2,rate=b.s2)
inv.s2[1] = 1

# epsilon_i - normal random errors
epsilon = array(data=0,dim=c(nm,1,nb))

# gamma
gamma = array(data=0,dim=c(nm,1,nb))
for(i in b){
  if(length(ind.a[[i]])==0){
    gamma[,1,i] = log(n) + mu[1] + theta[,1] + epsilon[,1,i]
  }
  else{
    gamma[,1,i] = log(n) + mu[1] + sum(alpha[ind.a[[i]],1]) + theta[,1] + ifelse(length(ind.a[[i]])>1,apply(phi[,1,ind.a[[i]]],MAR=1,FUN=sum),phi[,1,ind.a[[i]]]) + epsilon[,1,i]
  }
}
cur.gamma = gamma[,1,]

# the offset log(n) is added here, when we create the gammas

## saving the lambda's used to generate the synthetic datasets
lambda = array(data=0,dim=c(nm,sim-burn,nb))
# lambda = list(1); lambda = c(lambda,2:nb)
# for(i in 1:nb){
#   lambda[[i]] = matrix(0,nm,sim-burn)
# }

## running average for lambda = exp(gamma)
avg.lambda = matrix(0,nm,nb)

## synthetic counts
c.til = array(0,dim=c(nb,nm,n.syn))



###########################
## start the simulations ##
###########################

## functions for the adaptive rejection sampling steps

# distribution of mu
f.mu = function(x,n.f,v.f,sumeta.f){
  x*n.f - (1/(2*v.f))*x^2 - exp(x)*sumeta.f
}

fprima.mu = function(x,n.f,v.f,sumeta.f){
  n.f - x/v.f - exp(x)*sumeta.f
}

# distribution of theta_i
f.the = function(x,c.f,sumeta.f,ni.f,tau.f,bar.f){
  x*c.f - exp(x)*sumeta.f - (1/2)*ni.f*tau.f*(x-bar.f)^2
}

fprima.the = function(x,c.f,sumeta.f,ni.f,tau.f,bar.f){
  c.f - exp(x)*sumeta.f - ni.f*tau.f*(x-bar.f)
}

# distribution of epsilon_i
f.eps = function(x,c.f,sumeta.f,tau.f){
  x*c.f - exp(x)*sumeta.f - (1/2)*tau.f*x^2
}

fprima.eps = function(x,c.f,sumeta.f,tau.f){
  c.f - exp(x)*sumeta.f - tau.f*x
}


for(s in 2:sim){
  
  cat("\n Iteration ",s,"\n")
  
  ## update mu
  eta = cur.gamma - mu[s-1]
  sumeta = sum(exp(eta))
  mu[s] = ars(1,f.mu,fprima.mu,lb=T,xlb=-100,ub=T,xub=100,n.f=n,v.f=v,sumeta.f=sumeta)
  cur.gamma = eta + mu[s]
  
  ## update alpha
  for(i in 1:dim(alpha)[1]){
    n.aux = sum(nib[sub.a[[i]],])
    eta = cur.gamma[,sub.a[[i]]] - alpha[i,(s-1)]
    sum.aux = sum(exp(eta))
    alpha[i,s] = ars(1,f.mu,fprima.mu,lb=T,xlb=-100,ub=T,xub=100,n.f=n.aux,v.f=v,sumeta.f=sum.aux)
    cur.gamma[,sub.a[[i]]] = eta + alpha[i,s]
  }
  
  ## update theta
  c.aux = apply(counts,MAR=2,FUN=sum)
  eta = cur.gamma - theta[,(s-1)] # since the sumeta.f does not depend on the updated theta's
  for(i in m){
    sum.aux = sum(exp(eta[i,]))
    bar = W[i,]%*%cur.theta/ni[i]
    cur.theta[i] = ars(1,f.the,fprima.the,lb=T,ns=1000,xlb=-100,ub=T,xub=100,c.f=c.aux[i],sumeta.f=sum.aux,ni.f=ni[i],tau.f=tau[s-1],bar.f=bar)
  }
  ## make it sum to zero
  theta[,s] = cur.theta - sum(cur.theta)/nm
  cur.gamma = eta + theta[,s]
  
  ## update phi
  for(k in 1:dim(phi)[3]){
    c.aux = apply(counts[sub.a[[k]],],MAR=2,FUN=sum)
    eta = cur.gamma[,sub.a[[k]]] - phi[,(s-1),k]
    for(i in m){
      sum.aux = sum(exp(eta[i,]))
      bar = W[i,]%*%cur.phi[,k]/ni[i]
      cur.phi[i,k] = ars(1,f.the,fprima.the,lb=T,ns=1000,xlb=-100,ub=T,xub=100,c.f=c.aux[i],sumeta.f=sum.aux,ni.f=ni[i],tau.f=tau.phi[k,(s-1)],bar.f=bar)
    }
    ## make it sum to zero
    phi[,s,k] = cur.phi[,k] - sum(cur.phi[,k])/nm
    cur.gamma[,sub.a[[k]]] = eta + phi[,s,k]
  }
  
  ## update epsilon
  for(i in m){
    for(j in b){
      eta = cur.gamma[i,j] - epsilon[i,1,j]
      sum.aux = exp(eta)
      epsilon[i,1,j] = ars(1,f.eps,fprima.eps,lb=T,xlb=-100,ub=T,xub=100,c.f=counts[j,i],sumeta.f=sum.aux,tau.f=inv.s2[s-1])
      cur.gamma[i,j] = eta + epsilon[i,1,j]
    }
  }
  
  ## update gamma
  gamma[,1,] = cur.gamma
  
  ## update tau
  sum.theta = t(theta[,s])%*%(diag(ni)-W)%*%(theta[,s])
  tau[s] = rgamma(1,a.tau+nm/2,rate=b.tau+sum.theta/2)  
  
  ## update tau.phi
  for(k in 1:dim(tau.phi)[1]){
    sum.phi = t(phi[,s,k])%*%(diag(ni)-W)%*%(phi[,s,k])
    tau.phi[k,s] = rgamma(1,a.tau+nm/2,rate=b.tau+sum.phi/2)
  }
  
  ## update inv.s2
  sum.s2 = sum(epsilon[,1,]^2)
  inv.s2[s] = rgamma(1,a.s2+(nm*nb)/2,rate=b.s2+sum.s2/2)
  
  ## after burn-in, we can generate synthetic data sets & keep the lambdas for the risk measures
  
  if(s>burn){
    
    ## update running average of lambda  
    avg.lambda = avg.lambda*(s-burn-1)/(s-burn) + exp(gamma[,1,])/(s-burn)
    
    ## keep lambdas
    lambda[,s-burn,] = exp(gamma)
    
    if(any(list.syn==s)){
      
      ## generate synthetic dataset
      i.syn = which(list.syn==s)  # index of the synthetic data set
      gamma.s = matrix(gamma[,1,],nm,nb)
      
      syn.data[[i.syn]] = gen.syn(gamma.s,nb,size,n.comb,xvec,yvec)
      write(t(syn.data[[i.syn]]),file=paste("syndata",i.syn,".txt",sep=""),ncolumns=3)
      
      ## calculate the synthetic counts on each grid cell
      for(j in 1:nm){
        aux = vec2mat(j,size); row = aux[1,1]; col = aux[1,2]
        
        for(i in 1:nb){
          c.til[i,j,i.syn] = sum( syn.data[[i.syn]][,1]==i &
                                    (xvec[col]<syn.data[[i.syn]][,2] & syn.data[[i.syn]][,2]<xvec[col+1]) &
                                    (yvec[row]<syn.data[[i.syn]][,3] & syn.data[[i.syn]][,3]<yvec[row+1]) )
        }
      }
    }
  }
  
  #   cat("\n Simulations done \n")
  
}

#####################################################################################

## calculating the risk

cat("\n Calculating the risk \n")

## true gridcells
di = numeric(n)
for(i in 1:n){
  di[i] = mat2vec(r=min(which(sp[i,2]<yvec))-1,c=min(which(sp[i,1]<xvec))-1,size) # the x coord is going to give the col, and the y coord the row
}

## prior probabilities - low risk prior (uniform over all grid cells)
log.prior = matrix(log(1/nm),n,nm)
apply(exp(log.prior),MAR=1,FUN=sum)

S = sim-burn

## calculate the posterior probabilities based on the risk measure formula we derived
log.post = post = matrix(0,n,nm)
for(i in 1:n){  # for each person
  cat("\n Ind ",i)
  comb.i = comb[i]
  for(j in 1:n.syn){  # for each synthetic dataset
    aux.ij = apply(lambda[,,comb.i]^c.til[comb.i,,j],MAR=2,FUN=prod)
    for(k in m){
      log.post[i,k] = log.post[i,k]+log((1/S)/sum(lambda[k,,comb.i]/lambda[di[i],,comb.i])*
                                          sum( aux.ij*lambda[k,,comb.i]/lambda[di[i],,comb.i] ))  # product of the synthetic datasets
    }
    #     cat("...",100*j/n.syn,"%")
  }
  # now we have to multiply by the prior
  for(k in m){
    if(exp(log.prior[i,k])>0){
      log.post[i,k] = log.post[i,k]+log.prior[i,k]
    }
  }
  # finally, normalizing the posterior probabilities
  Const = min(log.post[i,])
  post[i,] = exp(log.post[i,]-Const)/sum(exp(log.post[i,]-Const))
  log.post[i,] = log(post[i,])
}

setwd("/net/nfs1/s/grad/tvp/Project - Jerry/Death data/size20")

# saving the posterior
write(t(post),file="post.txt",ncolumns=nm)

## saving the workspace
save.image("/net/nfs1/s/grad/tvp/Project - Jerry/Death data/size20/results_size20.RData")

## posterior analysis

## write results

cat("\n Saving the files and plots")

setwd("/net/nfs1/s/grad/tvp/Dropbox/PHD/Project - Jerry/Death data/size20")

# scatter plots
pdf("ts_mu.pdf", width=12, height=12, pointsize=24)
ts.plot(mu[(burn+1):sim])
dev.off()

pdf("ts_mu_all.pdf", width=12, height=12, pointsize=24)
ts.plot(mu)
dev.off()

pdf("ts_alpha.pdf", width=48, height=36, pointsize=24)
par(mfrow=c(4,3))
for(i in 1:dim(alpha)[1]){
  ts.plot(alpha[i,(burn+1):sim],main=paste("alpha ",i))
}
dev.off()

pdf("ts_alpha_all.pdf", width=48, height=36, pointsize=24)
par(mfrow=c(4,3))
for(i in 1:dim(alpha)[1]){
  ts.plot(alpha[i,],main=paste("alpha ",i))
}
dev.off()

# postscript("ts_beta.eps", horizontal=F, width=36, height=24, pointsize=24, paper="special")
# par(mfrow=c(2,3))
# for(i in 1:dim(beta)[1]){
#   ts.plot(beta[i,(burn+1):sim],main=paste("beta ",i))
# }
# dev.off()
# 
# postscript("ts_beta_all.eps", horizontal=F, width=36, height=24, pointsize=24, paper="special")
# par(mfrow=c(2,3))
# for(i in 1:dim(beta)[1]){
#   ts.plot(beta[i,],main=paste("beta ",i))
# }
# dev.off()

pdf("ts_invs2.pdf", width=12, height=12, pointsize=24)
ts.plot(inv.s2[(burn+1):sim])
dev.off()

pdf("ts_invs2_all.pdf", width=12, height=12, pointsize=24)
ts.plot(inv.s2)
dev.off()

pdf("ts_tau.pdf", width=12, height=12, pointsize=24)
ts.plot(tau[(burn+1):sim],main="tau")
dev.off()

pdf("ts_tau_all.pdf", width=12, height=12, pointsize=24)
ts.plot(tau,main="tau")
dev.off()

pdf("ts_tauphi.pdf", width=48, height=36, pointsize=24)
par(mfrow=c(4,3))
for(i in 1:dim(tau.phi)[1]){
  ts.plot(tau.phi[i,(burn+1):sim],main=paste("tau.phi ",i))
}
dev.off()

pdf("ts_tauphi_all.pdf", width=48, height=36, pointsize=24)
par(mfrow=c(4,3))
for(i in 1:dim(tau.phi)[1]){
  ts.plot(tau.phi[i,],main=paste("tau.phi ",i))
}
dev.off()

for(i in top4){
  #   post.g = matrix(apply(gamma[,(burn+1):sim,i],MAR=1,FUN=mean),(size-1),(size-1))
  #   postscript(paste("mean_gamma",i,".eps"), horizontal=F, width=12, height=12, pointsize=24, paper="special")
  #   image.plot(xvec, yvec, post.g, col=rev(heat.colors(30)), xlab="s1", ylab="s2", main=paste("Posterior mean of intensity",i),useRaster=T)
  #   dev.off()
  
  pdf(paste("mean_lambda",i,".pdf",sep=""), width=12, height=12, pointsize=24)
  image.plot(xvec, yvec, matrix(avg.lambda[,i],size-1,size-1), col=rev(heat.colors(30)), xlab="s1", ylab="s2", main=paste("Posterior mean of intensity",i),useRaster=T)
  points(sp[comb==i,],pch=19,cex=0.8)
  dev.off()
  
  pdf(paste("synt",i,".pdf",sep=""), width=24, height=24, pointsize=24)
  par(mfrow=c(2,2))
  # first, plot the average lambda w/ original data
  image.plot(xvec, yvec, matrix(avg.lambda[,i],size-1,size-1), col=rev(heat.colors(30)), xlab="s1", ylab="s2", main=paste("Posterior mean of intensity",i),useRaster=T)
  points(sp[comb==i,],pch=19,cex=0.8)
  # then, plot three synthetic datasets
  for(j in 1:3){
    image.plot(xvec, yvec, matrix(avg.lambda[,i],size-1,size-1), col=rev(heat.colors(30)), xlab="s1", ylab="s2", main=paste("Posterior mean of intensity",i),useRaster=T)
    points(syn.data[[j]][syn.data[[j]][,1]==i,2:3],pch=19,cex=0.8)
  }
  dev.off()
}

# files

# write(t(mu),ncol=sim,file="mu.txt")
# write(t(alpha),ncol=sim,file="alpha.txt")
# write(t(beta),ncol=sim,file="beta.txt")

# for(i in b){
#   write(exp(t(gamma[,(burn+1):sim,i])),ncol=sim-burn,file=paste("lambda",i,".txt"))
# }

# write(t(tau),ncol=sim,file="tau.txt")
# for(i in 1:dim(tau.phi)[1]){
#   write(t(tau.phi[i,]),ncol=sim,file=paste("tauphi",i,".txt"))
# }
# write(t(inv.s2),ncol=sim,file="invs2.txt")

cat("\n\n Done! :) \n")
