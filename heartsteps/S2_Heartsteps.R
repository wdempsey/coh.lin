## This code is similar to Peng's S2 but with Hearsteps changes. I changed
## the dimensions of the baseline coefficients
## the baseline generation distribution
## The coefficients for the baselines to get w_i (alpha0, beta0)

## Note: Peng's code previously didn't work with the baseline dimension being
## larger than 1, I had to change this.

rm(list = ls())
require(doParallel)
#require(igraph)
require(MASS)
registerDoParallel(6);
setwd("/Users/tnecamp/Dropbox/icml paper/new_stuff_NIPS/copying pengs nips/nips/sim/code")
source("functions_v2_heartsteps_states.R")  ## import functions for running coh.lin

nrep <- 1
### Simulation S2. Construct the Graph using Baseline

# distance measure in baseline
distb = function(b1, b2){
  
  eucl(b1-b2)
  
}
n <- 20;
dimB <- 4;
dimS <- 2;
p <- 2+2*dimS
sigma <- 1


nT <- 200*n
delta <- (1/nT)
BF <- sqrt(p);

#wid.all <- c(0.1, 0.25, 0.5, 1.0, 1.5, 2)
REG.all <- NULL
#for(wid in wid.all){
 wid = 2
  
  cat("********************eta", wid, " *******************\n")
  
  #dat <- foreach(icount(nrep), .inorder = FALSE) %dopar% {
   
    #alpha0 <- runif(dimB, min = -1, max = 1)  ## Intercept here
    #beta0 <- matrix(runif((p-1)*dimB, min = -1, max = 1), p-1, dimB)  ## here changed 5 to p-1, ## this should actually be a matrix that is p-1 by dimB ## changed it to a matrix
    
    coeff_vec = c(1.462, 0.362, 0.3, 0.003, 0.231, -0.038, 4.206, 0.197, -1.126, 0.007, 0.368, 0.079, 0.002, -0.033, -0.02, -0.004, -0.354, 0.273, -0.006, -0.034, 0.001, 0.011, -0.091, -0.205) ## take from linear MAB data gen
    
    alpha0 = c(coeff_vec[1], coeff_vec[4], coeff_vec[5], coeff_vec[6])
    beta0 = c(coeff_vec[7], coeff_vec[2], coeff_vec[3], coeff_vec[8], coeff_vec[9],
    			coeff_vec[10], coeff_vec[13], coeff_vec[16], coeff_vec[19], coeff_vec[22],
    			coeff_vec[11], coeff_vec[14], coeff_vec[17], coeff_vec[20], coeff_vec[23],
    			coeff_vec[12], coeff_vec[15], coeff_vec[18], coeff_vec[21], coeff_vec[24])
    beta0 = matrix(beta0, p-1, dimB)
    
    
    #LipC <- eucl(c(alpha0, beta0))
    #Bw <- (1+wid) * eucl(c(alpha0, beta0))
    # The eucl norm here works only when it is a vector. In the case where dimB > 1, we use have a matrix
    # we must use a new norm now, the norm must maintain the property in the paper, right now I'm just 
    # using the norm from a previous sim, I hope this works
    # it does work, a little too well!
    
    LipC <- norm(rbind(alpha0, (beta0)), type = "2")
    Bw <- (1+wid)*norm(rbind(alpha0, (beta0)), type = "2")
    
    # generate baseline and coef. w
    #b.all <- matrix(runif(n*dimB, 1-wid, 1+wid), n, dimB)   ## These are the baseline variables
    
   #  ## make people exactly the same
    #  age = rep(rpois(1,35), n)		## participant's age, male = 1
	  # gender = rep(rbinom(1,1,.351), n)		## participant's gender
	  # activity_base = rep(rnorm(1, mean = 0, sd = 2.3), n)			## baseline acitvity level, mean centered, might need to change to zero inflated thing
	
	## make people much closer
  age = sample(c(32, 33, 34, 35, 36, 37, 38), n, replace = TRUE )		## participant's age, male = 1
	gender = rbinom(n,1,.1)	## participant's gender
	activity_base = rnorm(n, mean = 0, sd = .5)		## baseline acitvity level, mean centered, might need to change to zero inflated thing
	
	# # ## original distn for baseline - this leads to a lot of variance in baseline distributions
  # ## leading to a very sparse graph
   # age = rpois(n,35)		## participant's age, male = 1
	 # gender = rbinom(n,1,.351)		## participant's gender
	 # activity_base = rnorm(n, mean = 0, sd = 2.3) 	## baseline acitvity level, mean centered, might need to change to zero inflated thing
	
	b.all = cbind(rep(1,n), age, gender, activity_base)
    
    
    alpha.all <- apply(b.all, 1, function(b) alphab(b,alpha0)) ## This is to multiply the baseline with the intercept
    beta.all <- t(apply(b.all, 1, function(b) betab(b, beta0))); ## This is to multiply the baseline with the non-intercept coefficients
    w.mat <- cbind(alpha.all, beta.all); # row i = coef. for user i 
    w.vec <- c(t(w.mat))
    

    
    reg.all <- cor.all <- c()
    
    # separate UCB
    temp <- ucb(w.mat, nT, lambda=0, gamma=1, delta=delta, L=matrix(0, n, n), Bdw = 0, Bw);
    temp1 = temp
    reg.all <- c(reg.all, temp$regret[nT])
    cor.all <- c(cor.all, temp$corr)
    
    # optimal graph and run coh.lin on it
    epsopt <- opteps (b.all, LipC, eps.range=c(0, 1), lam.range=c(0,10), gam.range = c(1e-2, 10), nT, Bw, BF, delta)
    L <- GraphConstr(b.all, eps.nw = epsopt, LipConst = LipC, F)
    
    Bdw0 = LipC * epsopt;
    gam.lam.opt <- optlamgam(nT, lam.range= c(0, 10), L, gam.range=c(1e-2, 10), Bdw0, Bw, BF, delta);
    lamopt <- gam.lam.opt[1];
    gammaopt <- gam.lam.opt[2];
    
    temp <- ucb(w.mat, nT, lambda=lamopt, gamma=gammaopt, delta=delta, L, Bdw = Bdw0, Bw);
    reg.all <- c(reg.all, temp$regret[nT])
    cor.all <- c(cor.all, temp$corr)


    ##  Plot the regret of sep ucbs vs coh.lin, note the regret is in log step count
	  plot(temp1$reg, cex = .1)
	  points(temp$reg, col = 2, cex = .1)

    ##  Plot the percentage of corrrect decisions made of sep ucbs vs coh.lin
	  plot(temp1$corr, cex = .1, ylim = c(0,1))
	  points(temp$corr, col = 2, cex = .1)
