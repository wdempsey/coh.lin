## This is similar to Peng's functions_v2, however here I 
## had to change the state distribution  to mimic the 
## Heartsteps data set


logdet = function(M) log(det(M))

eucl = function(x) sqrt(sum(x^2))

phi = function(i, x){
  
  temp <- rep(0, n*p);
  temp[((i-1)*p+1): (i*p)] <- x
  
  return(temp)
}


# Constrcut the graph (based on baselie) 
GraphConstr = function(b.all, eps.nw, LipConst, plot=F){
  
  Adj.mat <- matrix(0, n, n)
  
  for(i in 1:n){
    for(j in 1:n){
      
      if(i == j){
        
        next
        
      }else{
        
        Adj.mat[i, j] <- ifelse( distb(b.all[i, ], b.all[j, ]) < LipConst * eps.nw, 1, 0)
        
        
      }
      
    }
  }
  
  
  
  if(plot){
    plot(graph_from_adjacency_matrix(Adj.mat, "undirected"))
  }
  
  
  Deg.mat <- diag(apply(Adj.mat, 1, sum))
  L <- Deg.mat - Adj.mat;
  
  return(L)
  
}

# Regret Bound 
# Regret.Bnd = function(nT, lambda, gamma, L, Bdw, Bw, BF, delta){
#   
#   
#   Lo <- L %x% diag(p);
#   dbar <- mean(diag(L))
#   
#   LL <- BF # bound for features
#   
#   #Smooth.para <- eucl(Lo %*% w.vec);
#   Smooth.para <- Bdw * sqrt(sum(diag(L)^2))
#   
#   #Bound.para <- eucl(w.vec)
#   Bound.para <- sqrt(n) * Bw
#   
#   B <- sqrt(lambda * (sum(diag(L))/2) * Bdw^2 + gamma * n * Bw^2)
#   #B = (lambda/sqrt(gamma)) * Smooth.para + sqrt(gamma) * Bound.para
#   
#   V0 <- lambda * Lo + gamma * diag(n*p);
#   
#   fT <- ((n*p)* log(gamma + lambda*dbar + (nT*LL^2)/(n*p))) - logdet(V0);
#   #fT <- n * p * log(1 + nT*LL^2/(n*p*gamma));
#   #fT <- n * p * log(1 + max(eigen(solve(V0))$values)*nT*LL^2/(n*p));
#   
#   
#   temp1 <- sigma * sqrt(2*log(1/delta) + fT) + B
#   temp2 <- sqrt(fT);
#   
#   regret.bnd.t <- sqrt(nT) * temp1 * temp2;
#   
#   return(regret.bnd.t)
# }


Regret.Bnd = function(nT, lambda, gamma, L, Bdw, Bw, BF, delta){
  
  
  Lo <- L %x% diag(p);
  dbar <- mean(diag(L))
  
  LL <- BF # bound for features
  
  #Smooth.para <- eucl(Lo %*% w.vec);
  Smooth.para <- Bdw * sqrt(sum(diag(L)^2))
  
  #Bound.para <- eucl(w.vec)
  Bound.para <- sqrt(n) * Bw
  
  B <- sqrt(lambda * (sum(diag(L))/2) * Bdw^2 + gamma * n * Bw^2)
  #B = (lambda/sqrt(gamma)) * Smooth.para + sqrt(gamma) * Bound.para
  
  V0 <- lambda * Lo + gamma * diag(n*p);
  
  
#   ft = function(t){
#     #((n*p)* log(gamma + lambda*dbar + (t*LL^2)/(n*p))) - logdet(V0);
#     n * p * log(1 + t*LL^2/(n*p*gamma))
#   }
#   
#   t_th_term = function(t){
#     
#     tmp <- sigma * sqrt(2*log(1/delta) + ft(t)) + B;
#     return(tmp^2)
#   }
#   
#   regret.bnd.t <- sqrt(sum(sapply(1:nT, t_th_term)) * ft(nT))
#   
   
  fT <- ((n*p)* log(gamma + lambda*dbar + (nT*LL^2)/(n*p))) - logdet(V0);
#   #fT <- n * p * log(1 + nT*LL^2/(n*p*gamma));
#   #fT <- n * p * log(1 + max(eigen(solve(V0))$values)*nT*LL^2/(n*p));
#   
  temp1 <- sigma * sqrt(2*log(1/delta) + fT) + B
  temp2 <- sqrt(fT);
  regret.bnd.t <- sqrt(nT) * temp1 * temp2;
  
  return(regret.bnd.t)
}


# select the optimal lambda based on graph
optlambda = function(nT, lam.range, L, gamma, Bdw, Bw, BF, delta){
  
  fn = function(x) Regret.Bnd(nT, lambda = x, gamma = gamma, L, Bdw, Bw, BF, delta)
  opt <- optim(0, fn,  method = "L-BFGS-B", lower = 0, upper = lam.range[2])$par
  
  k <- 0 
  while(abs(opt-lam.range[2]) < 0.01 & k <= 3){
    
    # close to the upper bound 
    
    opt0 <- opt
    opt <- optim(opt, fn,  method = "L-BFGS-B", lower = 0, upper = 2*lam.range[2])$par
    k <- k + 1
  }
  if(k == 3) warning("SHIT in selecting lambda")
  return(opt)
  
}

optlamgam = function(nT, lam.range, L, gam.range, Bdw, Bw, BF, delta){
  
  fn = function(x) Regret.Bnd(nT, lambda = x[1], gamma = x[2], L, Bdw, Bw, BF, delta)
  opt <- optim(c(1, 1), fn,  method = "L-BFGS-B", lower = c(lam.range[1], gam.range[1]), 
               upper = c(lam.range[2], gam.range[2]))$par
  
#   k <- 0 
#   while(abs(opt-lam.range[2]) < 0.01 & k <= 3){
#     
#     # close to the upper bound 
#     
#     opt0 <- opt
#     opt <- optim(opt, fn,  method = "L-BFGS-B", lower = 0, upper = 2*lam.range[2])$par
#     k <- k + 1
#   }
#   if(k == 3) warning("SHIT in selecting lambda")
  return(opt)
  
}

# selct the optimal graph using basline
opteps = function(b.all, LC, eps.range, lam.range, gam.range = gam.range, nT, Bw, BF, delta, plot = F){
  
  
  obj = function(eps){
    
    L <- GraphConstr(b.all, eps.nw = eps, LipConst = LC, plot=F) # no edge
    
    lamgamopt <- optlamgam(nT, lam.range= lam.range, L, gam.range=gam.range, LC * eps, Bw, BF, delta);
    
    opt.reg <- Regret.Bnd(nT, lambda = lamgamopt[1], gamma = lamgamopt[2], L, Bdw = LC * eps, Bw, BF, delta)
    
    return(opt.reg)
  }
  
  fn = function(x){
    
    # input the range. Cut into 101 points and select the one that minmizes the obj
    
    x.seq <- seq(x[1], x[2], length.out = 21)
    temp <- sapply(x.seq, obj);
    opt <- x.seq[which.min(temp)]
    
    if(plot) plot(x = x.seq, y = temp, xlab = "Thres", ylab = "regret", type = "o")
    
    return(opt)
  }
  
  opt <- fn(eps.range)
  
  
  if(opt==eps.range[1]){
    
    # the range is too large. 
    stop("Something is wrong--picking 0 as the optimal threshold")
    
    
  }
  
  
  j <- 0
  while(opt == eps.range[2] && j <= 10){
    
    eps.range[2] <- 2 * eps.range[2]
    opt <- fn(eps.range)
    j <- j+1
    
  }
  if(j == 10){
    
    stop("The upper limit is not large enough and cannot fixed by mulitplied by 1000")
    
  }
  
  # in the middle
  stopifnot(all(opt<eps.range[2] & opt > eps.range[1]))
  
  k <- 0
  dif <- 1e4
  while(k <= 10 & dif > 1e-3){
    
    eps.range <- c(opt - diff(eps.range)/10, opt + diff(eps.range)/10)
    
    opt.nxt <- fn(eps.range)
    dif <- abs(opt.nxt - opt)
    opt <- opt.nxt
    k <- k + 1
  }
  if(k == 10) warning("SHIT in selecting epsilon")
  
  
  return(opt)
}

# evaluate the epsilon (thres)
evaleps = function(b.all, LC, eps, lam.range, nT, gamma, Bw, BF, delta, plot=F){
  
  
  obj = function(eps){
    
    L <- GraphConstr(b.all, eps.nw = eps, LipConst = LC, plot=F) # no edge
    
    lamgamopt <- optlamgam(nT, lam.range= c(0, 10), L, gam.range=c(1e-2, 10), LC * eps, Bw, BF, delta);

    opt.reg <- Regret.Bnd(nT, lambda = lamgamopt[1], gamma = lamgamopt[2], L, Bdw = LC * eps, Bw, BF, delta)
    
    return(opt.reg)
  }
  
  
  if(length(eps) == 1){
    
    return(obj(eps))
    
  }else{
    
    temp <- sapply(eps, obj)
    if(plot) plot(x = eps, y = temp, xlab = "Thres", ylab = "regret", type = "o", pch = 16, lty = 2)
    
    return(temp)
    
  }
  
}

# LIN-UCB with known network
ucb = function(w.mat, nT, lambda, gamma, delta, L, Bdw, Bw, print=F){
  
  x.feat = function(s, a) c(1, a, s, s*a)
  
  ## this is where I changed it to reflect hearsteps, there are two spots to change
  genStates = function() {#runif(dimS, -1, 1)
  	pre_steps = rnorm(1, mean = 2.7, sd = 3.03)		## previous step count
	engaged_level = runif(1)			## measure of engagement
	state_vec = c(pre_steps, engaged_level)
  	return(state_vec)
  	}
  
  ExpRwrd = function(s, a, i)  sum(x.feat(s, a) * w.mat[i, ])
  
  genRewards = function(s, a, i) ExpRwrd(s, a, i) + rnorm(1, sd=sigma)
  

  
  # initialization
  Lo <- L %x% diag(p);
  
#   #Smooth.para <- eucl(Lo %*% w.vec);
#   Smooth.para <- Bdw * sqrt(sum(diag(L)^2))
#   
#   #Bound.para <- eucl(w.vec)
#   Bound.para <- sqrt(n) * Bw
  
  B <- sqrt(lambda * (sum(diag(L))/2) * Bdw^2 + gamma * n * Bw^2)
  #B = (lambda/sqrt(gamma)) * Smooth.para + sqrt(gamma) * Bound.para;

  
  
  V <- V0 <- lambda * Lo + gamma * diag(n*p)
  b <- rep(0, n*p);
  
  
  #   c0 <- max(apply(w.mat, 1, function(x) sum(abs(x))))
  #   init.beta <- sigma * sqrt(2*log(1/delta)) + B;
  #   print(c(c0, init.beta))
  #   stopifnot(init.beta > c0) 
  
  
  
  # save 
  Regret <- c()
  Correct <- c()
  
  for(t in 1:nT){
    
    
    # user index
    it <- sample(n, 1);
    
    # context for user it
    St <- genStates()
    
    # UCB to select At
    V.inv <- solve(V)
    w <- V.inv %*% b;
    beta.delta <- sigma * sqrt(2*log(1/delta) + logdet(V)-logdet(V0)) + B;
    phit.a <- cbind(phi(it, x.feat(s = St, a = 0)), phi(it, x.feat(s = St, a = 1)))
    At <- which.max(t(phit.a) %*% w + apply(phit.a, 2, function(x) sqrt(t(x) %*% V.inv %*% x)) * beta.delta)-1
    phit <- phit.a[, At + 1]
    
#     cat("Time =",t, "\n")
#     cat("Point estimates =", c(t(phit.a) %*% w), "\n")
#     cat("UCB", c(apply(phit.a, 2, function(x) sqrt(t(x) %*% V.inv %*% x)) * beta.delta), "\n")
    # reward
    Rt = genRewards(s = St, a = At, i = it)
    #Rt = sum(x.feat(s = St, a = At) * w.mat[it, ]) + rnorm(1, sd = sigma)
    
    # updates for UCB
    V <- V + phit %*% t(phit);
    b <- b + Rt *  phit
    
    # optmimal action
    rwrd <- c(ExpRwrd(s = St, a = 0, i = it), ExpRwrd(s = St, a = 1, i = it))
    #rwrd <- c(sum(x.feat(s = St, a = 0) * w.mat[it, ]), sum(x.feat(s = St, a = 1) * w.mat[it, ]))
    At.opt <- ifelse(rwrd[1] <= rwrd[2], 1, 0);
    
    # pseudo regret 
    rt =  max(rwrd) -  ExpRwrd(s = St, a = At, i = it)
    
    
    
    # save 
    Regret <- c(Regret, rt);
    Correct <- c(Correct, (At == At.opt))
    
    
    #
    if(t%% (nT/10) == 0){
      if(print) {cat((t/nT)*100, "%", "\n")}
      
    }
  }
  
  list(corr = cumsum(Correct)/(1:nT),  regret = cumsum(Regret), lambda = lambda)
  
}

# Coef. function of baseline (linear)
alphab = function(b, alpha)  sum(b*alpha)

betab = function(b, beta) c( beta %*% b)


# GoB
GoB = function(w.mat, nT, delta, L, Bdw, Bw, print=F){
  
  x.feat = function(s, a) c(1, a, s, s*a)
  
  genStates = function() {#runif(dimS, -1, 1)
  	pre_steps = rnorm(1, mean = 2.7, sd = 3.03)		## previous step count
	engaged_level = runif(1)			## measure of engagement
	state_vec = c(pre_steps, engaged_level)
  	return(state_vec)
  	}
  
  ExpRwrd = function(s, a, i)  sum(x.feat(s, a) * w.mat[i, ])
  
  genRewards = function(s, a, i) ExpRwrd(s, a, i) + rnorm(1, sd=sigma)
  
  

  Lo <- L %x% diag(p);
  Ao <- Lo + diag(n*p)
  Ao.invsqrt <- mat.sqrt(solve(Ao));
  
  norm.Util =  sqrt(n*Bw^2 + Bdw^2 * (sum(diag(L)/2)));
  
  # initialization
  M <- M0 <- diag(n*p);
  b <- rep(0, n*p);
  

  # save 
  Regret <- c()
  Correct <- c()
  
  for(t in 1:nT){
    
    
    # user index
    it <- sample(n, 1);
    
    # context for user it
    St <- genStates()
    
    # UCB to select At
    M.inv <- solve(M)
    w <- M.inv %*% b;
    beta.delta <- sigma * sqrt(2*log(1/delta) + logdet(M)) + norm.Util;
    phit.a <- Ao.invsqrt %*% cbind(phi(it, x.feat(s = St, a = 0)), phi(it, x.feat(s = St, a = 1)))
    At <- which.max(t(phit.a) %*% w + apply(phit.a, 2, function(x) sqrt(t(x) %*% M.inv %*% x)) * beta.delta)-1
    phit <- phit.a[, At + 1]
    
    
#     cat("Time =",t, "\n")
#     cat("Point estimates =", c(t(phit.a) %*% w), "\n")
#     cat("UCB", c(apply(phit.a, 2, function(x) sqrt(t(x) %*% M.inv %*% x)) * beta.delta), "\n")

    
    # reward
    Rt = genRewards(s = St, a = At, i = it)
    #Rt = sum(x.feat(s = St, a = At) * w.mat[it, ]) + rnorm(1, sd = sigma)
    
    # updates for UCB
    M <- M + phit %*% t(phit);
    b <- b + Rt *  phit
    
    # optmimal action
    rwrd <- c(ExpRwrd(s = St, a = 0, i = it), ExpRwrd(s = St, a = 1, i = it))
    #rwrd <- c(sum(x.feat(s = St, a = 0) * w.mat[it, ]), sum(x.feat(s = St, a = 1) * w.mat[it, ]))
    At.opt <- ifelse(rwrd[1] <= rwrd[2], 1, 0);
    
    # pseudo regret 
    rt =  max(rwrd) -  ExpRwrd(s = St, a = At, i = it)
    
    
    
    # save 
    Regret <- c(Regret, rt);
    Correct <- c(Correct, (At == At.opt))
    
    
    #
    if(t%% (nT/10) == 0){
      if(print) {cat((t/nT)*100, "%", "\n")}
      
    }
  }
  
  list(corr = sum(Correct)/nT,  regret = cumsum(Regret))
  
}

mat.sqrt = function(M){
  
  eigen.decm <- eigen(M);
  M.sqrt <- eigen.decm$vectors %*% diag(sqrt(eigen.decm$values)) %*% t(eigen.decm$vectors)
  
}



