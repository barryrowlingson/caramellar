create.W <- function(delta, nbs, ids, inv=FALSE) {
    p = match(ids, nbs$ids) # as needed p is the same order as inverse.perm in carmc!

    nbs$Adj = nbs$Adj[p,p]
    nbs$Distances = nbs$Distances[p,p]
    
    n = length(nbs$ids)
    W <- Matrix(0,n, n,sparse=TRUE)
    for(i in 1:n) { 
      if(inv==FALSE){
        W[i,][nbs$Adj[i,]] <- 1/(1+(nbs$Distances[i,][nbs$Adj[i,]]/delta)^2)
      }else{
        W[i,][nbs$Adj[i,]] =    (1+(nbs$Distances[i,][nbs$Adj[i,]]/delta)^2);
      }
    }
    W <- forceSymmetric(W)
    W = as.spam(as.matrix(W))
    return(W)	
}
