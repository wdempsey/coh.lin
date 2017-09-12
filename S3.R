rm(list = ls())
require(doParallel)
#require(igraph)
require(MASS)
registerDoParallel(20);
source("functions.R")

nrep <- 200
### Simulation S3. Sensitivity Analysis Guessing the Lip Constant 

# distance measure in baseline
distb = function(b1, b2){
  
  eucl(b1-b2)
  
}

n <- 20;
dimB <- 1;
dimS <- 2;
p <- 2+2*dimS
sigma <- 1

nT <- 100*n
delta <- (1/nT)
lam.range = c(0, 10)
BF <- sqrt(p);


LLL.all <- c(0.1, 0.5, 1, 5, 10)
wid.all <- c(0.1, 0.25, 0.5, 1.0, 1.5, 2)



REG.all <- NULL
for(wid in wid.all){
  
  cat("********************eta", wid, " *******************\n")
  
  dat <- foreach(icount(nrep), .inorder = FALSE) %dopar% {
    
    
    
    alpha0 <- runif(1, min = -1, max = 1)
    beta0 <- matrix(runif(5, min = -1, max = 1), p-1, 1)
    
    # generate baseline and coef. w
    b.all <- matrix(runif(n*dimB, 1-wid, 1+wid), n, dimB)
    alpha.all <- apply(b.all, 1, function(b) alphab(b,alpha0))
    beta.all <- t(apply(b.all, 1, function(b) betab(b, beta0)));
    w.mat <- cbind(alpha.all, beta.all); # row i = coef. for user i 
    w.vec <- c(t(w.mat))
    
    
    reg.all <- NULL
    
    for(LLL in LLL.all){
      
      # wrong Lip Const
      LipC <- eucl(c(alpha0, beta0))*LLL
      Bw <- (1+wid) * eucl(c(alpha0, beta0))
      
      epsopt <- opteps (b.all, LipC, eps.range=c(0, 1), lam.range=c(0,10), gam.range = c(1e-2, 10), nT, Bw, BF, delta)
      L <- GraphConstr(b.all, eps.nw = epsopt, LipConst = LipC, F)
      
      Bdw0 = LipC * epsopt;
      gam.lam.opt <- optlamgam(nT, lam.range= c(0, 10), L, gam.range=c(1e-2, 10), Bdw0, Bw, BF, delta);
      lamopt <- gam.lam.opt[1];
      gammaopt <- gam.lam.opt[2];
      
      temp <- ucb(w.mat, nT, lambda=lamopt, gamma=gammaopt, delta=delta, L, Bdw = Bdw0, Bw);
      
      
      ####
      reg.all <- rbind(reg.all, c(wid, LLL, temp$regret[nT], epsopt, lamopt, gammaopt))
      
      
    }
    
    reg.all
  }
  
  
  dat <- Reduce('+', dat)/length(dat)
  
  
  print(dat)
  REG.all <- rbind(REG.all, dat[, c(1, 2, 3)])
  save(dat, file = paste("S3_wid", wid, ".RData" ,sep=""))
  
}

save(REG.all, file = "S3.Rdata");
print(REG.all);



load("/Users/Peng/Dropbox/Tim:Walter:Peng/network MAB/nips/sim/result/untitled folder/S3.Rdata")
dat <- data.frame(REG.all);
colnames(dat) <- c("eta", "L", "regret")

dat[, 1] <- round(dat[, 1], 2)
dat[, 2] <- round(dat[, 2], 1)

opt.r <- NULL
for(LLL in unique(dat[, 2])){
  print(subset(dat, L == LLL))
  opt.r <- rbind(opt.r, subset(dat, L == LLL)[, 3])
}

load("/Users/Peng/Dropbox/Tim:Walter:Peng/network MAB/nips/sim/result/untitled folder/S2.Rdata")
sep.r <- REG.all[, 2]
opt.r[3, ] <-  REG.all[, 3]

plot(opt.r[3, ],  col = 5, lty = 1, pch = 16, type = "o", ylim = c(100, 600),
     xaxt='n', xlab = "eta in generating baseline", ylab = "Regret")
points(opt.r[1, ], col = 2, type = "o", lty = 1, pch = 16,  cex = 1.2)
points(opt.r[2, ], col = 4, type = "o", lty = 1, pch = 16,  cex = 1.2)
points(opt.r[4, ], col = 3, type = "o", lty = 1, pch = 16,  cex = 1.2)
points(opt.r[5, ], col = 6, type = "o", lty = 1, pch = 16,  cex = 1.5)



points(sep.r, col = 1, type = "o", lty = 1, pch = 16,  cex = 1.2)

axis(1, at = 1:nrow(REG.all), labels = REG.all[,1])

legend("bottomright", 
       col=c(2, 4, 5, 3, 6, 1), lty = 1, pch = 16, cex = 1, 
       legend =c("0.1Lc", "0.5Lc", "Lc", "5Lc", "10Lc", "Sep"), 
       text.width = 0.40, pt.cex = 1.0, pt.lwd = 1.0, ncol=2)



