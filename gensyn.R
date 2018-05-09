## function to generate synthetic datasets 

source("~/Dropbox/PHD/Project - Jerry/Gibbs sampler/mat2vec.R")

## supposing gamma = array(0,dim=c(nm,sim,nb))
## gamma.s = gamma[,s,]

gen.syn = function(gamma.s,nb,size,n.comb,xvec,yvec){
  
  synt.data = numeric(0)
  
  for(j in 1:nb){
    lambda = exp(gamma.s[,j])    # we need the exponential of the gamma 
    lambda = lambda/sum(lambda)    # normalizing  
    
    # select a gridcell with prob lambda
    slct = sample(1:((size-1)^2),n.comb[j],replace=T,prob=lambda)
    row.col = vec2mat(slct,size)
    row = row.col[,1]; col = row.col[,2]
    
    # then, generate the new location uniformly inside the selected gridcell
#     aux = matrix(c(rep(j,n.comb[j]),runif(n.comb[j],xvec[row],xvec[row+1]),runif(n.comb[j],yvec[col],yvec[col+1])),n.comb[j],3,byrow=F)
    ## CHANGING THE ORDER TO SEE IF IT WILL WORK
        aux = matrix(c(rep(j,n.comb[j]),runif(n.comb[j],xvec[col],xvec[col+1]),runif(n.comb[j],yvec[row],yvec[row+1])),n.comb[j],3,byrow=F)
#     synt.data = rbind(synt.data,aux[,c(1,3:2)])  
        synt.data = rbind(synt.data,aux)
    ## THE COLUMNS AND ROWS ARE SWITCHED, THAT'S WHY IT'S SAVING THE COLUMNS IN A DIFFERENT ORDER
    ## maybe this will fix it
  }
    
  return(synt.data)
}