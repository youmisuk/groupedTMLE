#########################################################################################################
# A Within-Group Approach to Ensemble Machine Learning Methods for Causal Inference in Multilevel Studies
# Youmi Suk 

# :: load packages
library(dplyr)
library(lme4)
library(tmle)

# :: load or write functions
source("DataGeneratingModels.R") # load a function for data generating processes (DGP) 

# write a function for the clustered IPW estimator from Li et al (2013). Propensity score weighting with multilevel data. Statistics in Medicine, 32(19), 3373–3387. doi: 10.1002/sim.5786
li.clusterATE <- function(ps, Y, Z, id) {
  
  # ps ... propensity scores
  # Y ... outcome
  # Z ... treatment
  # id ... cluster id
  
  dataf <- data.frame(Y, Z, ps, id)
  dataf$temp.wt <- with(dataf, Z/ps + (1-Z)/(1-ps))
  temp.dat <- dataf %>% group_by(id) %>% summarize(clusterATE=IPW_ate(temp.wt, Y, Z), clusterSE=IPW_se(temp.wt, Y, Z), wt=sum(temp.wt), .groups = 'drop') 
  temp.ate = sum(temp.dat$clusterATE * temp.dat$wt, na.rm = TRUE) / sum(temp.dat$wt[!is.na(temp.dat$clusterATE)])
  temp.se = sqrt(sum(temp.dat$wt^2*temp.dat$clusterSE^2, na.rm = TRUE)) / sum(temp.dat$wt[is.finite(temp.dat$clusterSE)], na.rm = TRUE)
  
  return(list(ATE=c(temp.ate, temp.se), clusterATE=temp.dat$clusterATE, clusterSE=temp.dat$clusterSE, sumIPW = temp.dat$wt)) # Li's method
}

# write a function to get the weighted mean of groups
GroupWeightedMean <- function(x, w) {
  weighted.mean(x=x, w=w)
}

# write a function to estimate the average treatment effect (ATE) from an inverse propensity weighting (IPW) estimator
IPW_ate = function(IPW, Y, Z) {
  ATE = sum((IPW*Y)[Z==1])/sum(IPW[Z==1]) - sum((IPW*Y)[Z==0])/sum(IPW[Z==0])
  return(ATE)
}

# write a function to estimate the standard error (SE) of the ATE from an IPW estimator
IPW_se = function(IPW, Y, Z) {
  SE = sqrt(var(Y)*(sum((IPW[Z==1])^2)/(sum(IPW[Z==1]))^2 + sum((IPW[Z==0])^2)/(sum(IPW[Z==0]))^2))
  return(SE)
}


# :::: Simulation codes ::::

args = commandArgs(trailingOnly=TRUE) # arguments in a command line for HTcondor
a = as.numeric(args[1]) # the first argument, the number of groups
rep.c = as.numeric(args[2]) # the second argument, the replicate number

n.group = as.numeric(a) # the number of groups
n.rep = 2 # the number of replicates 

# create matrices that will save simulation results
marginalest_reps <- clusterest_reps <- groupest_reps <- matrix(NA, nrow=n.rep, ncol=4)
sampleATE <- rep(NA, n.rep)
adjTMLEest_reps <- matrix(NA, nrow=n.rep, ncol=5)

# conduct simulations
for (i in 1:n.rep) {
  set.seed(i+2000+rep.c)
  
  # change the DGPs: Design 1 (m.val=0) versus Design 2 & 3 (m.val=2)
  dat <- twolevel.uniformU(m.val=0, Smpl.size = "170.15") # create data from Design 1
#  dat <- twolevel.uniformU(m.val=2, Smpl.size = "170.15") # create data from Design 2 & 3  
  
  sampleATE[i] <-  mean(dat$Y1-dat$Y0) # the sample mean
  
  trt.prop = tapply(dat$Z, dat$id, mean) # the treatment prevalence
  
  km.out = kmeans(scale(trt.prop[dat$id]), n.group, nstart=50) # conduct k-mean clustering

  dat$group = km.out$cluster # group membership
  
  groupTMLE.est.group <- matrix(NA, nrow=n.group, ncol=5)   
  
  if (sum(table(dat$Z, dat$group) == 0) == 0) {
    
    for (g in 1:n.group) {
      
      temp.group.re.PS <- NULL
      temp.group.fe.PS <- NULL
      temp.group.re.mis.PS <-  NULL
      temp.group.fe.mis.PS <- NULL
      
      temp.dat <- dat[dat$group == g,]
      temp.dat$id <- factor(temp.dat$id)
      
      if (length(levels(temp.dat$id)) > 1) {
        
        temp.group.re <- tryCatch({
          glmer(Z~ X1 + X2 + W1  + I(X2 < 0.3) +  (1|id), data = temp.dat, family=binomial(link="logit"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5), calc.derivs = FALSE))
        }, error = function(e) {
          print(paste(e))
          glm(Z ~  X1 + X2 + W1  + I(X2 < 0.3), data = temp.dat, family=binomial)
        })
        
        temp.group.re.mis <- tryCatch({
          glmer(Z ~  X1 + X2 + W1 + (1|id), data = temp.dat, family=binomial(link="logit"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5), calc.derivs = FALSE))
        }, error = function(e) {
          print(paste(e))
          glm(Z ~  X1 + X2 + W1, data = temp.dat, family=binomial)
        })
        
        temp.group.fe <- glm(Z ~ X1 + X2  + I(X2 < 0.3) + factor(id), data = temp.dat, family = binomial)
        temp.group.fe.mis <- glm(Z ~ X1 + X2 + as.factor(id), data = temp.dat, family = binomial) 
        
        temp.group.re.PS <- predict(temp.group.re, type='response')
        temp.group.fe.PS <- predict(temp.group.fe, type='response')
        temp.group.re.mis.PS <-  predict(temp.group.re.mis, type='response')
        temp.group.fe.mis.PS <- predict(temp.group.fe.mis, type='response')
        
        temp.out.TMLE.g <- tmle(Y=temp.dat$Y, A=temp.dat$Z, W=temp.dat[,c("X1", "X2", "W1")])
        temp.out.TMLE.g.ps.re <- tmle(Y=temp.dat$Y, A=temp.dat$Z, W=temp.dat[,c("X1", "X2", "W1")], g1W=temp.group.re.PS)
        temp.out.TMLE.g.ps.fe <- tmle(Y=temp.dat$Y, A=temp.dat$Z, W=temp.dat[,c("X1", "X2", "W1")], g1W=temp.group.fe.PS)
        temp.out.TMLE.g.ps.re.mis <- tmle(Y=temp.dat$Y, A=temp.dat$Z, W=temp.dat[,c("X1", "X2", "W1")], g1W=temp.group.re.mis.PS)
        temp.out.TMLE.g.ps.fe.mis <- tmle(Y=temp.dat$Y, A=temp.dat$Z, W=temp.dat[,c("X1", "X2", "W1")], g1W=temp.group.fe.mis.PS)
        
      } else {
        temp.group.re <- glm(Z ~  X1 + X2  + I(X2 < 0.3), data = temp.dat, family=binomial)
        temp.group.fe <- glm(Z ~  X1 + X2  + I(X2 < 0.3), data = temp.dat, family=binomial)   
        temp.group.re.mis <- glm(Z ~  X1 + X2 , data = temp.dat, family=binomial)
        temp.group.fe.mis <- glm(Z ~  X1 + X2 , data = temp.dat, family=binomial)          
        
        temp.group.re.PS <- predict(temp.group.re, type='response')
        temp.group.fe.PS <- predict(temp.group.fe, type='response')
        temp.group.re.mis.PS <-  predict(temp.group.re.mis, type='response')
        temp.group.fe.mis.PS <- predict(temp.group.fe.mis, type='response')
        
        temp.out.TMLE.g <- tmle(Y=temp.dat$Y, A=temp.dat$Z, W=temp.dat[,c("X1", "X2")])
        temp.out.TMLE.g.ps.re <- tmle(Y=temp.dat$Y, A=temp.dat$Z, W=temp.dat[,c("X1", "X2")], g1W=temp.group.re.PS)
        temp.out.TMLE.g.ps.fe <- tmle(Y=temp.dat$Y, A=temp.dat$Z, W=temp.dat[,c("X1", "X2")], g1W=temp.group.fe.PS)
        temp.out.TMLE.g.ps.re.mis <- tmle(Y=temp.dat$Y, A=temp.dat$Z, W=temp.dat[,c("X1", "X2")], g1W=temp.group.re.mis.PS)
        temp.out.TMLE.g.ps.fe.mis <- tmle(Y=temp.dat$Y, A=temp.dat$Z, W=temp.dat[,c("X1", "X2")], g1W=temp.group.fe.mis.PS)       
      }
      
      dat$group.re.PS[dat$group == g] <- temp.group.re.PS
      dat$group.fe.PS[dat$group == g] <- temp.group.fe.PS
      dat$group.re.mis.PS[dat$group == g] <- temp.group.re.mis.PS
      dat$group.fe.mis.PS[dat$group == g] <- temp.group.fe.mis.PS          
      
      temp.TMLE_ate <- summary(temp.out.TMLE.g)$estimates$ATE$psi
      temp.TMLE_ate.ps.re <- summary(temp.out.TMLE.g.ps.re)$estimates$ATE$psi
      temp.TMLE_ate.ps.fe <- summary(temp.out.TMLE.g.ps.fe)$estimates$ATE$psi
      temp.TMLE_ate.ps.re.mis <- summary(temp.out.TMLE.g.ps.re.mis)$estimates$ATE$psi
      temp.TMLE_ate.ps.fe.mis <- summary(temp.out.TMLE.g.ps.fe.mis)$estimates$ATE$psi     
      
      groupTMLE.est.group[g,] <- c(temp.TMLE_ate, temp.TMLE_ate.ps.re, temp.TMLE_ate.ps.fe, temp.TMLE_ate.ps.re.mis, temp.TMLE_ate.ps.fe.mis)
    }
    
    # marginal IPW estimator
    marginalest_reps[i, ] <- c(li.clusterATE(dat$group.re.PS,  dat$Y, dat$Z, 1)$ATE[1],
                               li.clusterATE(dat$group.fe.PS,  dat$Y, dat$Z, 1)$ATE[1],
                               li.clusterATE(dat$group.re.mis.PS,  dat$Y, dat$Z, 1)$ATE[1],
                               li.clusterATE(dat$group.fe.mis.PS,  dat$Y, dat$Z, 1)$ATE[1])
    
    # clustered IPW estimator
    clusterest_reps[i, ] <- c(li.clusterATE(dat$group.re.PS,  dat$Y, dat$Z, dat$id)$ATE[1],
                              li.clusterATE(dat$group.fe.PS,  dat$Y, dat$Z, dat$id)$ATE[1],
                              li.clusterATE(dat$group.re.mis.PS,  dat$Y, dat$Z, dat$id)$ATE[1],
                              li.clusterATE(dat$group.fe.mis.PS,  dat$Y, dat$Z, dat$id)$ATE[1])
    
    # grouped IPW estimator
    groupest_reps[i, ] <-  c(li.clusterATE(dat$group.re.PS,  dat$Y, dat$Z, dat$group)$ATE[1],
                             li.clusterATE(dat$group.fe.PS,  dat$Y, dat$Z, dat$group)$ATE[1],
                             li.clusterATE(dat$group.re.mis.PS,  dat$Y, dat$Z, dat$group)$ATE[1],
                             li.clusterATE(dat$group.fe.mis.PS,  dat$Y, dat$Z, dat$group)$ATE[1])
    
    # TMLE estimator
    adjTMLEest_reps[i, ] <- 
      
      c(apply(groupTMLE.est.group, 2, function(x) GroupWeightedMean(x, w=table(dat$group))))
    
  } 
  
}

allrlst <- cbind(as.numeric(sampleATE), marginalest_reps, clusterest_reps, groupest_reps, adjTMLEest_reps)

allrlst
