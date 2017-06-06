
adjust_proposal <- function(h, acc, i){
    max(10e-07,h + (0.001*i^(-0.0001))*(acc/i-0.45))
}

use_dlnorm <- function(mean.l, sd.l){
    force(mean.l)
    force(sd.l)
    D = function(s){
        dlnorm(s, mean.l, sd.l, log=TRUE)
    }
    D
}

use_dbeta <- function(a, b){
    force(a)
    force(b)
    D = function(s){
        dbeta(s,a,b,log=TRUE)
    }
    D
}

##' Specification for carmc MCMC
##'
##' MCMC model parameters
##' @title default carmc MCMC model parameters
##' @param sigma2 control for sigma-2
##' @param rho control for rho
##' @param phi control for phi
##' @param beta control for beta
##' @param delta control for delta
##' @return a list of control parameters
##' @author Barry Rowlingson
##' @export
mcmc.defaults <- function(
###
### sigma-squared
###
    
    sigma2=list(
        log.prior=use_dlnorm(-5, 3),
        init=exp(-5),
        h.l = 0.5,
        adjust = adjust_proposal
        ),

###
### rho
###

    rho = list(
        init = 0.98,
        log.prior = use_dbeta(1,1),
        h.t = 0.5,
        adjust = adjust_proposal
        ),

###
### phi
###

    phi = list(
        init = 0.7,
        log.prior=use_dbeta(1,1),
        h.t = 0.5,
        adjust = adjust_proposal
        ),

###
### beta
###
    beta = list(
        mean  = function(p){
            rep(0,p)
        },
        cov.inv = function(p){
            diag(1/1000,p)
        }
        ),
###
### delta
###
    delta = list(
        values=seq(1,100,1),
        start_at=10,
        step_size = 2
        )
    ){

    force(sigma2)
    force(rho)
    force(phi)
    force(beta)
    force(delta)
    
    list(sigma2 = sigma2,
         rho = rho,
         phi = phi,
         beta = beta,
         delta = delta)
    
}
