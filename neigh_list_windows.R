#### Function to create neighbors list
#### Thais Paiva, Sep 2011

## create neigthbors list from a vector of grid indexes

## input: nm = length of the grid (as a vector, i=1:m)
##        size = size of the grid (nm = (size-1)*(size-1))

## output: neigh = list of size nm
##                 neigh[[i]] contains the vector index of the grid cells
##                 that are neighbors of the i-th grid cell
##                 neighbors are considered as if they share vertex


neigh.list = function(nm,size){
  source("C:/Users/Thais/Documents/My Dropbox/PHD/Project - Jerry/Gibbs sampler/mat2vec.R")
  #source("~/Dropbox/PHD/Project - Jerry/Gibbs sampler/mat2vec.R")
  neigh = list(1)
  neigh = c(neigh,2:nm)

  for(i in 1:nm){
    aux = vec2mat(i,size)
    r = aux[1,1]
    c = aux[1,2]

    if(r==1){
      if(c==1){
        neigh[[i]] = mat2vec(c(r,r+1,r+1),c(c+1,c,c+1),size)
      }
      else{
        if(c==(size-1)){
          neigh[[i]] = mat2vec(c(r,r+1,r+1),c(c-1,c-1,c),size)
        }
        else{
          neigh[[i]] = mat2vec(c(r,r,r+1,r+1,r+1),c(c-1,c+1,c-1,c,c+1),size)
        }
      }
    }
    else{
      if(r==(size-1)){
        if(c==1){
          neigh[[i]] = mat2vec(c(r-1,r-1,r),c(c,c+1,c+1),size)
        }
        else{
          if(c==(size-1)){
            neigh[[i]] = mat2vec(c(r-1,r-1,r),c(c-1,c,c-1),size)
          }
          else{
            neigh[[i]] = mat2vec(c(r-1,r-1,r-1,r,r),c(c-1,c,c+1,c-1,c+1),size)
          }
        }
      }
      else{
        if(c==1){
          neigh[[i]] = mat2vec(c(r-1,r-1,r,r+1,r+1),c(c,c+1,c+1,c,c+1),size)
        }
        else{
          if(c==(size-1)){
            neigh[[i]] = mat2vec(c(r-1,r-1,r,r+1,r+1),c(c-1,c,c-1,c-1,c),size)
          }
          else{
            neigh[[i]] = mat2vec(c(r-1,r-1,r-1,r,r,r+1,r+1,r+1),
                   c(c-1,c,c+1,c-1,c+1,c-1,c,c+1),size)
          }
        }
      }
    }
  }
  return(neigh)
}



