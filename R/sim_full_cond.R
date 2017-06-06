create.sim.full.cond <- function(D, n.x.t, mcmcbeta, ID.space.time){
    ndata = nrow(D)
    p <- ncol(D)
    ind.S <- sapply(1:n.x.t,function(i) which(ID.space.time==i))

    C.S <- Matrix(0,ndata,n.x.t,sparse=TRUE)
    for(i in 1:n.x.t) {	
	if(length(ind.S[[i]])>0) C.S[ind.S[[i]],i] <- 1
    }
    
    sel.beta <- 1:p
    sel.S <- (p+1):(n.x.t+p)
    n.S <- sapply(1:n.x.t,function(i) length(ind.S[[i]]))
    
    mean.beta <- mcmcbeta$mean(p)
    cov.beta.inv <- mcmcbeta$cov.inv(p)
    cov.beta.S <- Matrix(t(D)%*%C.S,sparse=TRUE)
    
    cov.beta.beta <- cov.beta.inv+t(D)%*%D
    
    cond.inv.cov <- Matrix(0,n.x.t+p,n.x.t+p,sparse=TRUE)
    cond.inv.cov[sel.beta,sel.beta] <- cov.beta.beta
    cond.inv.cov[sel.beta,sel.S] <- cov.beta.S
    cond.inv.cov[sel.S,sel.beta] <- t(cov.beta.S)
    cond.inv.cov <- forceSymmetric(cond.inv.cov)
    
    sim.full.cond <- function(V,par.sim) {
        cond.inv.cov[sel.S,sel.S] <- par.sim$Sigma.inv
        diag(cond.inv.cov[sel.S,sel.S]) <- diag(cond.inv.cov[sel.S,sel.S])+n.S           
        beta.V <- t(D)%*%V
        S.V <- sapply(1:n.x.t,function(i) {
            if(length(ind.S[[i]]>0)) {
                return(sum(V[ind.S[[i]]]))
            } else {
                return(0)
            }
        }) 
         
        b <- c(cov.beta.inv%*%mean.beta+t(D)%*%V,S.V)
        
        mean.vec <- as.numeric(solve(cond.inv.cov,b,sparse=TRUE))
        
        mu.beta <- mean.vec[sel.beta]
        mu.S <- mean.vec[sel.S]
   
        cond.sim <- backsolve(chol(cond.inv.cov),rnorm(p+n.x.t))+
            c(mu.beta,mu.S)
        sim.res <- list()
        sim.res$beta <- cond.sim[sel.beta]
        sim.res$S <- cond.sim[sel.S] 
        return(sim.res) 
    }
    return(sim.full.cond)
}
