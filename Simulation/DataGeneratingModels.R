#########################################################################################################
# A Within-Group Approach to Ensemble Machine Learning Methods for Causal Inference in Multilevel Studies
# by Youmi Suk 

# ::: Data Generating Models :::
twolevel.uniformU <- function(ids=1, Smpl.size = "150.30", m.val=2) {
  
  # ::::: 1) generate the number of observations in each cluster :::::
  
  size <- strsplit(Smpl.size, "\\.")[[1]]  # the number of cluster and cluster sizes
  
  J <- as.numeric(size[1]) # level-2 unit, the number of cluster
  n.clus <- as.numeric(size[2]) # level-1 unit, cluster sizes
  
  if (n.clus < 20) {  # variation of cluster sizes 
    sd.clus <- 0.5
  } else if (n.clus >= 20 & n.clus < 30) { 
    sd.clus <- 1 
  } else if (n.clus >= 30) {
    sd.clus <- 2
  }
  
  N <- round(rnorm(J, n.clus, sd.clus))                   # normally distributed cluster sizes 
  id <- as.factor(rep(ids:(J + ids - 1), N))              # cluster id
  
  # ::::: 2) generate level-1 covariates, Xs with cluster-specific means :::::
  
  totalN <- length(id)
  X1 <- runif(totalN, -1, 1)
  X2 <- runif(totalN, 0, 1)
  
  # ::::: 3) generate level-2 covariates, W :::::
  
  W1 <- runif(J, -1, 1)
  W2 <- runif(J, -2, 2) # W2 is an unmeasured cluster-level confounder U_j that follows a uniform distribution.	  
  names(W1) <- names(W2)  <- levels(id) 
  
  pop <- data.frame(id, X1, X2, W1=W1[id], W2=W2[id])
  
  # ::::: 4) generate selection probabilities and potential outcome ::::::::
  
  E <- rnorm(sum(N), 0, 1)   # error terms for pot.  
  
  pop$lps <- -0.6 + 0.3* pop$X1 + 0.3*(pop$X2 +pop$W1 + pop$W2) +  0.4 * I(pop$X2 < 0.3)   # ps logit
  pop$Y0 <- 70  + 2*(pop$X1 + pop$X2 + pop$W1 + I(pop$X2  < 0.3) + pop$W2) + E  
  pop$Y1 <- pop$Y0 + 2 + 2*pop$X2 + 2*pop$W1 + m.val*(pop$W2)^3
  
  
  pop$ps <- 1 / (1 + exp(-pop$lps))   # propensity scores
  pop$Z <- rbinom(nrow(pop), 1, pop$ps) # treatment indicator  
  pop$Y <- with(pop, ifelse(Z == 1, Y1, Y0))
  
  pop
} 
