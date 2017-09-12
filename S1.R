rm(list = ls())
require(doParallel)
#require(igraph)
require(MASS)
registerDoParallel(20);
# source("/Users/Peng/Dropbox/Tim:Walter:Peng/network MAB/code/peng/functions.R")
source("functions.R")

nrep <- 200
### Simulation S1. Fixed the graph


n <- 20;
dimB <- 1;
dimS <- 2;
p <- 2+2*dimS
sigma <- 1

L <- matrix(-1, n, n)
diag(L) <- n-1

nT <- 100*n
delta <- (1/nT)
BF <- sqrt(p) # bound for features




wid.all <- c(0.01, 0.05, 0.1, 0.25, 0.5, 1)

REG.all <- NULL
for(wid in wid.all){
  
  cat("********************eta", wid, " *******************\n")
  
  # separate 
  
  dat <- foreach(icount(nrep), .inorder = FALSE) %dopar% {
    
    alpha0 <- runif(1, min = -1, max = 1)
    beta0 <- matrix(runif(5, min = -1, max = 1), p-1, 1)
    
    
    Bw <- (1+wid) * eucl(c(alpha0, beta0))
    Bdw0 <- 2*wid * eucl(c(alpha0, beta0))
    
    gam.lam.opt <- optlamgam(nT, lam.range= c(0, 10), L, gam.range=c(1e-2, 10), Bdw0, Bw, BF, delta);
    lamopt <- gam.lam.opt[1];
    gammaopt <- gam.lam.opt[2];
    
    
    
    b.all <- matrix(runif(n*dimB, 1-wid, 1+wid), n, dimB)
    alpha.all <- apply(b.all, 1, function(b) alphab(b,alpha0))
    beta.all <- t(apply(b.all, 1, function(b) betab(b, beta0)));
    w.mat <- cbind(alpha.all, beta.all); # row i = coef. for user i 
    w.vec <- c(t(w.mat))
    
    reg.all <- cor.all <- c()
    temp <- ucb(w.mat, nT, lambda=0, gamma=1, delta=delta, L, Bdw = Bdw0, Bw);
    reg.all <- c(reg.all, temp$regret[nT])
    cor.all <- c(cor.all, temp$corr)
    
    temp <- ucb(w.mat, nT, lambda=lamopt, gamma=gammaopt, delta=delta, L, Bdw = Bdw0, Bw);
    reg.all <- c(reg.all, temp$regret[nT])
    cor.all <- c(cor.all, temp$corr)
    
    temp <- GoB(w.mat, nT, delta=delta, L, Bdw = Bdw0, Bw)
    reg.all <- c(reg.all, temp$regret[nT])
    cor.all <- c(cor.all, temp$corr)
    
    list(regret = reg.all, corr = cor.all, lamopt = lamopt, gammaopt = gammaopt)
  }
  
  avg.corr <- Reduce('+', lapply(dat, "[[", "corr"))/length(dat);
  avg.reg <- Reduce('+', lapply(dat, "[[", "regret"))/length(dat);
  
  dat <- list(corr = avg.corr, regret = avg.reg, 
              gamma = lapply(dat, "[[", "gammaopt"),
              lambda = lapply(dat, "[[", "lamopt"))
  
  
  save(dat, file = paste("S1_wid", wid,".RData" ,sep=""))
  
  
  # separate 
  
  cat("separate LinUCB:",
      ", Accuracy =", round(dat$corr[1], 3), 
      ", Regret =", round(dat$regret[1], 3), "\n")
  
  
  
  # optimal 
  
  cat("Optimal LinUCB:", 
      ", Accuracy =", round(dat$corr[2], 3), 
      ", Regret =", round(dat$regret[2], 3), "\n")
  
  # GoB
  
  cat("GoB",
      ", Accuracy =", round(dat$corr[3], 3), 
      ", Regret =", round(dat$regret[3], 3), "\n")
  
  
  REG.all <- rbind(REG.all, c(wid, dat$regret))
  
}

save(REG.all, file = "S1.Rdata");

load("/Users/Peng/Dropbox/Tim:Walter:Peng/network MAB/nips/sim/result/untitled folder/S1.Rdata")
temp <- REG.all[3:6, 2] # sep
pdf("S1.pdf",width=8,height=6) 
plot(REG.all[, 2], ylim = c(80, 520), type = "o", pch = 16, xlab = "eta in generating baseline", xaxt="n", ylab = "Regret")
axis(1, at = 1:nrow(REG.all), labels = REG.all[,1])
points(REG.all[, 3], col = 2, type = "o", pch = 16)
points(REG.all[, 4], col = 4, type = "o", pch = 16)
legend("bottomright", col = c(1, 2, 4), pch = 16, legend = c("Sep.Lin", "Coh.Lin", "GoB.Lin"), text.width = 0.8)


dev.off()
