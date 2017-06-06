create.Q.t <- function(phi, nt, nb_time) { ### speed up gained by using array instead of matrix!
  Q.t = array(0,dim=c(nt,nt))             
  Q.t[1,1] <- Q.t[nt,nt] <- 1+{phi^2}/{1-phi^2}
  diag(Q.t)[-c(1,nt)] <- 1+2*{phi^2}/{1-phi^2}
  Q.t[nb_time] <- -phi/{1-phi^2}
  Q.t <- Q.t*{1-phi^2}
  Q.t = as.spam(Q.t)
  return(Q.t)
}

create.Q <- function(rho,W){
  Q <- -rho*W
  diag.spam(Q) = rowSums.spam(W)
  return(Q)
}

log.posterior <- function(S,Q,Q.t,l.sigma2,t.phi,t.rho,logdet,log.prior.sigma2,log.prior.phi,log.prior.rho,n.x,n.t,n.x.t) { # added by ACH - faster kronecker product method using Kron(B',A)Vec(X)=Vec(AXB)
  phi = 1/{1+exp(-t.phi)}
  rho = 1/{1+exp(-t.rho)}
  sigma2 = exp(l.sigma2)
  r=log.prior.sigma2(sigma2) + l.sigma2 +
    log.prior.phi(phi)+t.phi-2*log(1+exp(t.phi)) +
    log.prior.rho(rho)+t.rho-2*log(1+exp(t.rho)) +
    -{-n.x*log(1-phi^2)-n.t*logdet+n.x.t*l.sigma2 +
        as.numeric(S%*%as.vector(as.matrix(Q%*%matrix(S,ncol=n.t)%*%t.spam(Q.t))))/sigma2}/2
  return(r)
}

make_lim <- function(model){
    y <- model$y
    ind0 <- which(y==0)
    ind1 <- which(y==1)
    lim = matrix(NA,nrow=length(y),ncol=2)
    lim[ind0,1] <- -Inf
    lim[ind0,2] <- 0
    lim[ind1,1] <- 0
    lim[ind1,2] <- Inf
    lim
}

## compute the indexes of the 1-off-diagonal elements of an n.t x n.t matrix.
## Z=Matrix(0,5,5)
## Z[make_nb_time(5)]=1
## Z
## 5 x 5 sparse Matrix of class "dgCMatrix"
##
## [1,] . 1 . . .
## [2,] 1 . 1 . .
## [3,] . 1 . 1 .
## [4,] . . 1 . 1
## [5,] . . . 1 .

make_nb_time <- function(n.t){
    nb_time <- NULL
    for(i in 1:(n.t-1)) {
        nb_time <- c(i*n.t+i,(i-1)*n.t+i+1,nb_time)
    }
    nb_time
}



##' Space-time Conditional Autoregressive MCMC
##'
##' Fit the space-time conditional AR model via MCMC
##' @title Conditional Autoregressive Space-time MCMC
##' @param cases data frame of case data
##' @param space_time formula of location_id~date
##' @param linear_model formula of caseflag~covariates
##' @param neighbours neighbourhood structure
##' @param mc_control control mcmc loops, thinning, burn-in
##' @param mcmc priors, proposals, etc
##' @param note a text note to add to the document
##' @param progress passed to txtProgressBar to create a progress bar
##' @return an object of class "carmc"
##' @author Emanuele Giorgi, Barry Rowlingson, Alison Hale
##' @export
carmc <- function(cases, space_time, linear_model, neighbours,
                  mc_control,
                  mcmc = mcmc.defaults(),
                  note = NA,
                  progress){
  
    op = spam.options()        # save inital options
    spam.options(cholsymmetrycheck=FALSE,cholpivotcheck=FALSE,safemode=c(FALSE,FALSE,FALSE))
  
    if(!igraph::is.connected(igraph::graph.adjacency(neighbours$Adjacencies))){
        stop("Non-connected neighbourhood structure")
    }

    out <- list()
    out$note = note

    out$timing = list(
        clock=list(start=Sys.time()),
        cpu = list(start=proc.time())
        )

    n.sim <- mc_control$n.sim
    burnin <- mc_control$burnin
    thin <- mc_control$thin
    ## sample acceptances etc at this resolution
    hthin <- mc_control$hthin
    nhsamples <- n.sim %/% hthin

    if(!missing(progress)){
        progress$min = 0
        progress$max=n.sim
        pb = do.call(txtProgressBar, progress)
    }

    ## get the location IDs and case dates
    spaces_times <- model.frame(space_time, data=cases)
    realDates <- spaces_times[,2] # Date objects
    spaces_times[,2] <- as.numeric(spaces_times[,2]) # convert Date to number
    spatial_ids <- unique(spaces_times[,1])

    AdjMat       = neighbours$Adjacencies/(length(spatial_ids)^2) + diag(1,length(spatial_ids))
    inverse.perm = rev(chol.spam(as.spam(AdjMat),pivot=TRUE)@pivot) # this oringially was inverse.perm = rev((chol.spam(as.spam(AdjMat[length(spatial_ids):1,length(spatial_ids):1]),pivot=TRUE))@pivot) # which is wrong!
    spatial_ids  = neighbours$ids[inverse.perm]  # n.b. chol should be sparse: display(chol.spam(as.spam(Q[length(spatial_ids):1,length(spatial_ids):1]),pivot=FALSE))
    rm(AdjMat,inverse.perm)                      # Permutations work as follows: A=chol(as.spam(AdjMat),pivot=TRUE); ord=A@pivot; B=chol(as.spam(AdjMat[ord,ord]),pivot=FALSE); sum(A@entries-Z@entries);
    
    if(!all(sort(spatial_ids) == sort(neighbours$ids))){
      stop("spatial id and neighbour ids mismatch")
    }    
    
    times <- spaces_times[,2]
    n.x <- length(spatial_ids)

    ## Construct all the W matrices
    delta.set <- mcmc$delta$values
    n.delta <- length(delta.set)
    ind.delta.curr <- mcmc$delta$start_at

    W.ind.ref = 1
    W.inv.ref = create.W(delta.set[W.ind.ref], neighbours, spatial_ids, inv=TRUE)
    delta.ref = delta.set[W.ind.ref]
    W.prop = W.inv.ref # initalise
    W.curr = W.inv.ref

    time <- sort(unique(times))
    n.t <- length(time)
    nb_time <- make_nb_time(n.t)

    ndata <- nrow(cases)

    time <- sort(unique(times))

    ID.space.time = IDspacetime(spaces_times,spatial_ids)

    n.x.t <- n.x * n.t

    ## get initial beta coeffecients from simple GLM
    mod.glm <- glm(formula=linear_model,family=binomial(link="probit"),data=cases,x=TRUE)

    D <- mod.glm$x
    p <- ncol(D)

    beta.curr <- coef(mod.glm)
    beta.names=names(beta.curr)
    mu.curr <- as.numeric(D%*%beta.curr)

    lim <- make_lim(mod.glm)

    V.curr <- truncnorm::rtruncnorm(ndata,a=lim[,1],b=lim[,2],
                         mean=mu.curr,
                         sd=1)

    rho.curr <- mcmc$rho$init # initial rho

    delta.factor = (delta.set[ind.delta.curr]/delta.ref)^2
    W.curr@entries = delta.factor/(delta.factor-1+W.inv.ref@entries)

    Q.curr <- create.Q(rho.curr,W.curr)
    Q.det.curr = determinant.spam(Q.curr)$modulus

    phi.curr <- mcmc$phi$init # initial phi
    Q.t.curr <- create.Q.t(phi.curr, n.t, nb_time)

    sigma2.curr <- mcmc$sigma2$init
    l.sigma2.curr <- log(sigma2.curr)

    sim.full.cond.FAST = create.sim.full.cond.FAST(D, n.x.t, mcmc$beta, ID.space.time)
    cond.sim.curr = sim.full.cond.FAST(V.curr,Q.t.curr,Q.curr,sigma2.curr)    

    V.curr <- rtruncnorm(ndata,a=lim[,1],b=lim[,2],
                         mean=D%*%cond.sim.curr$beta+
                         cond.sim.curr$S[ID.space.time],
                         sd=1)

    ## prior for rho:
    log.prior.rho <- mcmc$rho$log.prior

    t.rho.curr <- log(rho.curr/(1-rho.curr))

    ## prior for phi:
    log.prior.phi <- mcmc$phi$log.prior
    t.phi.curr <- log(phi.curr/(1-phi.curr))

    ## prior for sigma2
    log.prior.sigma2 = mcmc$sigma2$log.prior

    n.samples <- (n.sim-burnin)/thin
    
    out$S <- matrix(NA,nrow=n.samples,ncol=n.x.t)
    out$beta <- matrix(NA,nrow=n.samples,ncol=ncol(D))
    out$delta <- rep(NA,n.samples)
    out$sigma2 <- rep(NA,n.samples)
    out$rho <- rep(NA,n.samples)
    out$phi <- rep(NA,n.samples)
    
    ## track proposal sd and acceptance here
    out$prop = list(
        sigma2 = list(
            h.l = rep(NA,nhsamples),
            a = rep(NA,nhsamples)
            ),
        rho = list(
            h.t = rep(NA,nhsamples),
            a = rep(NA,nhsamples)
            ),
        phi = list(
            h.t = rep(NA,nhsamples),
            a = rep(NA,nhsamples)
            )
        )

    ## starting value of sd of proposals for sigma2, phi, rho

    h.l.sigma2 <- mcmc$sigma2$h.l
    acc.l.sigma2 <- 0

    h.t.phi <- mcmc$phi$h.t
    acc.t.phi <- 0

    h.t.rho <- mcmc$rho$h.t
    acc.t.rho <- 0

    lp.curr = log.posterior(cond.sim.curr$S,
                            Q.curr,Q.t.curr,
                            l.sigma2.curr,t.phi.curr,t.rho.curr,
                            Q.det.curr,
                            log.prior.sigma2,log.prior.phi,log.prior.rho,n.x,n.t,n.x.t)

    dstep = mcmc$delta$step_size

### MCMC simulations ###
    
    for(i in 1:n.sim) {
        hsave = ((i-1)%%hthin)==0
        hpos = 1 + ((i-1)%/%hthin)
        
    ### Update delta
        ## dstep is the proposal step max size
        set.delta.curr <- max(ind.delta.curr-dstep,1):min(ind.delta.curr+dstep,n.delta) # suppose we what the step between successive delta's to be something other than one?
        set.delta.curr = set.delta.curr[set.delta.curr!=ind.delta.curr]
        ind.delta.prop <- sample(set.delta.curr,1)

        set.delta.prop <- max(ind.delta.prop-dstep,1):min(ind.delta.prop+dstep,n.delta)
        set.delta.prop = set.delta.prop[set.delta.prop!=ind.delta.prop]

        dp.curr <- -log(length(set.delta.curr))
        dp.prop <- -log(length(set.delta.prop))

        delta.factor = {delta.set[ind.delta.prop]/delta.ref}^2
        W.prop@entries = delta.factor/{delta.factor-1+W.inv.ref@entries}

        Q.prop <- create.Q(rho.curr, W.prop)
        Q.det.prop = determinant.spam(Q.prop)$modulus

        lp.prop = log.posterior(cond.sim.curr$S,
                                Q.prop,Q.t.curr,
                                l.sigma2.curr,t.phi.curr,t.rho.curr,
                                Q.det.prop,
                                log.prior.sigma2,log.prior.phi,log.prior.rho,n.x,n.t,n.x.t)

# cppdelta=log_posterior_cpp(cond.sim.curr$S,
#                         par.sim.prop$Sigma.inv,
#                         sigma2.curr,
#                         l.sigma2.curr,
#                         t.phi.curr,
#                         t.rho.curr,
#                         par.sim.prop$log_det.Q,
#                         n.x,n.t,n.x.t)   
# cat(cppdelta-lp.prop,"  ,  ")

        if(log(runif(1)) < lp.prop+dp.prop-lp.curr-dp.curr) {
            ind.delta.curr <- ind.delta.prop
            W.curr  <- W.prop
            lp.curr =  lp.prop # this orginially read lp.curr <- lp.curr #which is surely wrong!
            Q.curr  <- Q.prop
            Q.det.curr = Q.det.prop
        }
        rm(ind.delta.prop,Q.prop,Q.det.prop,dp.curr,dp.prop,set.delta.curr,set.delta.prop)

    ### Update rho
	      t.rho.prop <- t.rho.curr+h.t.rho*rnorm(1)
	      rho.prop = 1/{1+exp(-t.rho.prop)}
	      Q.prop <- create.Q(rho.prop, W.curr)
	      Q.det.prop = determinant.spam(Q.prop)$modulus

        lp.prop = log.posterior(cond.sim.curr$S,
                                Q.prop,Q.t.curr,
                                l.sigma2.curr,t.phi.curr,t.rho.prop,
                                Q.det.prop,
                                log.prior.sigma2,log.prior.phi,log.prior.rho,n.x,n.t,n.x.t)
        
        log.ratio.rho <- lp.prop-lp.curr

        if(log(runif(1)) < log.ratio.rho) {
            t.rho.curr <- t.rho.prop
            rho.curr   <- rho.prop
            lp.curr    <- lp.prop
            Q.curr     <- Q.prop
            Q.det.curr = Q.det.prop
            acc.t.rho <- acc.t.rho+1
        }
        rm(lp.prop,t.rho.prop,rho.prop,Q.prop,Q.det.prop)
        ##  proposal adjustment
        h.t.rho <- mcmc$rho$adjust(h.t.rho, acc.t.rho, i)
        if(hsave){
            out$prop$rho$h.t[hpos] = h.t.rho
            out$prop$rho$a[hpos] = acc.t.rho/i
        }

    ### Update phi
	      t.phi.prop <- t.phi.curr+h.t.phi*rnorm(1)
	      phi.prop = 1/{1+exp(-t.phi.prop)}
	      Q.t.prop <- create.Q.t(phi.prop, n.t, nb_time)

        #par.sim.prop <- create.par.sim(sigma2.curr,Q.curr,Q.t.prop,log_det.Q=FALSE)
        #par.sim.prop$log_det.Q <- par.sim.curr$log_det.Q ### I've assumed this line is correct!
        
	      lp.prop = log.posterior(cond.sim.curr$S,
                                Q.curr,Q.t.prop,
                                l.sigma2.curr,t.phi.prop,t.rho.curr,
	                              Q.det.curr,
                                log.prior.sigma2,log.prior.phi,log.prior.rho,n.x,n.t,n.x.t)

        log.ratio.phi <- lp.prop-lp.curr

        if(log(runif(1)) < log.ratio.phi) {
            t.phi.curr <- t.phi.prop
            phi.curr   <- phi.prop
            lp.curr    <- lp.prop
            Q.t.curr   <- Q.t.prop
            acc.t.phi  <- acc.t.phi+1
        }
        rm(lp.prop,t.phi.prop,phi.prop,Q.t.prop)
        ## proposal adjustment
        h.t.phi <- mcmc$phi$adjust(h.t.phi, acc.t.phi, i)
        if(hsave){
            out$prop$phi$h.t[hpos] = h.t.phi
            out$prop$phi$a[hpos] = acc.t.phi/i
        }

    ### Update sigma2
	      l.sigma2.prop <- l.sigma2.curr+h.l.sigma2*rnorm(1)
	      sigma2.prop <- exp(l.sigma2.prop)

	      lp.prop = log.posterior(cond.sim.curr$S,
                                Q.curr,Q.t.curr,
                                l.sigma2.prop,t.phi.curr,t.rho.curr,
	                              Q.det.curr,
                                log.prior.sigma2,log.prior.phi,log.prior.rho,n.x,n.t,n.x.t)
	      
        log.ratio.sigma2 <- lp.prop-lp.curr

        if(log(runif(1)) < log.ratio.sigma2) {
            l.sigma2.curr <- l.sigma2.prop
            sigma2.curr   <- sigma2.prop
            lp.curr       <- lp.prop
            acc.l.sigma2 <- acc.l.sigma2+1
        }
        rm(lp.prop,l.sigma2.prop,sigma2.prop)
        ## adjust proposal
        h.l.sigma2 <- mcmc$phi$adjust(h.l.sigma2, acc.l.sigma2, i)
        if(hsave){
            out$prop$sigma2$h.l[hpos] = h.l.sigma2
            out$prop$sigma2$a[hpos] = acc.l.sigma2/i
        }

    ### Update V, U, beta
        cond.sim.curr = sim.full.cond.FAST(V.curr,Q.t.curr,Q.curr,sigma2.curr)
        mu.curr <- as.numeric(D%*%cond.sim.curr$beta)
	      V.curr <- rtruncnorm(ndata,a=lim[,1],b=lim[,2],
                             mean=mu.curr+cond.sim.curr$S[ID.space.time],
                             sd=1)

	      lp.curr = log.posterior(cond.sim.curr$S,
                                Q.curr,Q.t.curr,
                                l.sigma2.curr,t.phi.curr,t.rho.curr,
	                              Q.det.curr,
                                log.prior.sigma2,log.prior.phi,log.prior.rho,n.x,n.t,n.x.t)

        if(i > burnin & (i-burnin)%%thin==0) {
            Pos = {i-burnin}/thin
            out$S[Pos,] <- cond.sim.curr$S
            out$beta[Pos,]  <- cond.sim.curr$beta
            out$sigma2[Pos] <- sigma2.curr
            out$rho[Pos] <- rho.curr
            out$phi[Pos] <- phi.curr
            out$delta[Pos] <- delta.set[ind.delta.curr]

            if(as.duration(Sys.time() - out$timing$clock$start) > mc_control$tmax){
                print(Sys.time())
                print(out$timing$clock$start)
                print(mc_control$tmax)
                completed = 1:((i-burnin)/thin)
                warning("mcmc run timed out, only ",length(completed)," samples saved")
                out$completed = completed
                out$S=out$S[completed,,drop=FALSE]
                out$beta = out$beta[completed,,drop=FALSE]
                out$sigma2 = out$sigma2[completed,drop=FALSE]
                out$rho = out$rho[completed,drop=FALSE]
                out$phi = out$phi[completed,drop=FALSE]
                out$delta = out$delta[completed,drop=FALSE]
                break;
            }
        }

        if(!missing(progress)){
            setTxtProgressBar(pb,i)
        }
    }

    out$completed = i

    colnames(out$beta)=beta.names
    out$space = spatial_ids
    out$nSpace = length(spatial_ids)
    out$nTimes = n.t
    out$startDate = min(realDates)
    out$S = mc2array(out$S, out$nSpace, out$nTimes)
    attr(out$S, "ids") = spatial_ids
    attr(out$S, "start_date") = out$startDate

    out$mc_control=mc_control
    out$iteration <- function(index){
        (index * thin) + burnin
    }
    out$timing$clock$end=Sys.time()
    out$timing$cpu$end=proc.time()
    
    on.exit(spam.options(op))   # restore default spam options on exit
    
    attr(out, "ids") = spatial_ids
    class(out) <- "carmc"
    return(out)

}

carmcdates <- function(x){
    seq(as.Date(x$startDate), by=1, len=x$nTimes)
}

print.carmc <- function(x,...){
    cat("\n** carmc output **\n\n")
    if(!is.null(x$note)){
        if(!is.na(x$note)){
            cat("------------\n",x$note,"\n------------\n\n")
        }
    }
    cat("Dates: ",as.character(x$startDate)," to ",as.character(max(carmcdates(x))),"\n",sep="")
    cat(length(x$space)," locations\n\n",sep="")
    cat("MCMC run with ",x$mc_control$n.sim," iterations, burn: ",x$mc_control$burnin," thin: ", x$mc_control$thin,"\n",sep="")
    cat("\nTotal samples: ",length(x$delta),"\n",sep="")
    cat("\nCPU Time:\n")
    print(x$timing$cpu$end - x$timing$cpu$start)
    cat("\nClock Time:\n")
    cat("   Start: ",as.character(x$timing$clock$start),"\n",sep="")
    cat(" Run for: ",format(x$timing$clock$end-x$timing$clock$start),"\n",sep="")
    cat("\n")
    cat("Covariate betas:\n")
    print(summary(x$beta))
    cat("\nVariables:\n")
    print(summary(data.frame(delta = x$delta,
                             sigma2=x$sigma2,
                             phi=x$phi,
                             rho=x$rho)))
}

mc2array <- function(samples, nsites, ntimes){
    stopifnot(nsites*ntimes == ncol(samples))
    array(exp(samples), c(nrow(samples), nsites, ntimes))
}

location <- function(thing, ids){
    indexes = match(ids, location_ids(thing))
    if(all(is.na(indexes))){
        stop("No valid location ids in ",paste(ids,collapse=","))
    }
    if(any(is.na(indexes))){
        bad <- ids[is.na(indexes)]
        warning("Invalid location id(s): ",paste(bad, collapse=","))
    }

    indexes
}

location_ids <- function(thing){
    attr(thing,"ids")
}
