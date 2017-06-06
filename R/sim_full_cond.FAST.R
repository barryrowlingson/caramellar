
create.sim.full.cond.FAST <- function(D, n.x.t, mcmcbeta, ID.space.time){
  
    p <- ncol(D)
    ind.S <- sapply(1:n.x.t,function(i) which(ID.space.time==i))

    C.S <- Matrix(0,nrow(D),n.x.t,sparse=TRUE)
    for(i in 1:n.x.t) {
	    if(length(ind.S[[i]])>0) C.S[ind.S[[i]],i] <- 1
    }

    sel.beta <- 1:p                                       
    sel.S <- (p+1):(n.x.t+p)                              
    n.S <- sapply(1:n.x.t,function(i) length(ind.S[[i]])) 

    mean.beta <- mcmcbeta$mean(p)         
    cov.beta.inv <- mcmcbeta$cov.inv(p)
    cov.beta.S <- Matrix(t(D)%*%C.S,sparse=TRUE)

    cond.inv.cov <- Matrix(0,n.x.t+p,n.x.t+p,sparse=TRUE)
    cond.inv.cov[sel.beta,sel.beta] <- cov.beta.inv+t(D)%*%D
    cond.inv.cov[sel.beta,sel.S] <- cov.beta.S
    cond.inv.cov[sel.S,sel.beta] <- t(cov.beta.S)
    cond.inv.cov <- forceSymmetric(cond.inv.cov)

    ### for fast computation change 'cond.inv.cov.FAST' into a spam object
    cond.inv.cov.FAST = as(cond.inv.cov, "dgCMatrix")
    cond.inv.cov.FAST = as.spam.dgCMatrix(cond.inv.cov.FAST)
    padrow = spam(0,p,n.x.t)
    padcol = spam(0,n.x.t+p,p)
    RR = (n.x.t+p):1
    
    chnp = NULL    
    chwp = NULL
    ch_test = TRUE
    
    rm(C.S,cov.beta.S,cond.inv.cov)

    sim.full.cond.FAST <- function(V,Q.t,Q,sigma2) {

        SigInv = kronecker.spam(Q.t,Q)/sigma2 # N.b. sparse structure: display(chol.spam(sInv[n.x.t:1,n.x.t:1],pivot=FALSE))
        diag.spam(SigInv) = diag.spam(SigInv)+n.S
        SigInvBig = cbind.spam( padcol, rbind.spam(padrow,SigInv) )
        cond.inv.cov.ALL = cond.inv.cov.FAST + SigInvBig

        S.V = sapply( ind.S, function(x) if(length(x)>0){sum(V[x])}else{0} )
      
        if(ch_test){
          chwp <<- chol.spam(cond.inv.cov.ALL       ,pivot=TRUE)       # only do this first time loop is called
          chnp <<- chol.spam(cond.inv.cov.ALL[RR,RR],pivot=FALSE)
          chol_wp = chwp
          chol_np = chnp
          ch_test <<- FALSE
        }else{
          chol_wp = update.spam.chol.NgPeyton(chwp, cond.inv.cov.ALL)
          chol_np = update.spam.chol.NgPeyton(chnp, cond.inv.cov.ALL[RR,RR])
        }
        
        b = c( as.numeric(cov.beta.inv%*%as.matrix(mean.beta)+t(D)%*%as.matrix(V)), S.V )
        mean.vec = as.numeric( solve.spam(chol_wp,as.vector(b)) ) # the following is slower since spam does pivot each time: mean.vec = as.numeric(solve.spam(cond.inv.cov.ALL,as.vector(b)))
        cond.sim = (backsolve.spam( chol_np, as.vector(rnorm(p+n.x.t)) ) )[RR] + mean.vec
        
        return( list(beta=cond.sim[sel.beta], S=cond.sim[sel.S]) )
    }
    return(sim.full.cond.FAST)
}
