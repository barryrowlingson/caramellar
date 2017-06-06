IDspacetime <- function(sts,premise_ids){
    #premise_ids = unique(sts[,1]) # use spatial_ids which have been reordered for Cholesky decomposition
    n.x = length(premise_ids)
    ndata = nrow(sts)
    ID.space.time <- rep(NA,ndata)
    #time <- sort(unique(rr[,c("time")]))
    time = sort(unique(sts[,2]))
    coords.time <- expand.grid(1:n.x,time)
    n.x.t <- nrow(coords.time)
    
    for(i in 1:n.x.t) {
	ind.space <- which(sts[,1]==premise_ids[coords.time[i,1]] & 
	                   sts[,2]==coords.time[i,2])
        if(length(ind.space)>0) {
            ID.space.time[ind.space] <- i	
        }    	                               
    }
    return(ID.space.time)

}

