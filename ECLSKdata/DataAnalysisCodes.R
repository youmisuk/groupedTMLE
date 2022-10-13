#########################################################################################################
# A Within-Group Approach to Ensemble Machine Learning Methods for Causal Inference in Multilevel Studies
# : The ECLS-K data analysis 
# by Youmi Suk 

# This is our complete data. 
# The original ECLS-K 1998-99 dataset is available at https://nces.ed.gov/ecls/dataproducts.asp. 
# For more information on ECLS-K data, see Tourangeau et al. (2009) from https://nces.ed.gov/pubs2009/2009003.pdf.

# The variables in the dataset include:

# :: ID
# S7_ID    : school ID

# :: treatment
# C7DESMTH : whether students took the algebra or a higher-level math course (= 1) or not (= 0)

# :: outcome
# C7R4MSCL_s : math IRT scale score

# :: student-level covariates
# C6R4MSCL : prior math IRT scale score
# P6EXPECT : expected degree of child
# GENDER : male or female
# RACE : race/ethnicity
# WKSESL : socio-economic status measure
# WKPOV_R : poverty level
# WKMOMED : mother's education level
# P7HFAMIL : family type

# :: school-level covariates
# S7PUPRI : public or private school
# R7REGION : census region - WEST, MIDWEST, NORTHEAST, SOUTH
# R7URBAN : school location - urban, suburb, small town
####################################################################################################

# :: load packages/sources
library(dplyr)
library(lme4)
library(tmle)

# :: write useful functions
# write a function that creates inverse-propensity weights (IPW).
IPWfun <- function(ps, Z) Z/ps + (1-Z)/(1-ps) 

# write a generic function that implements the marginal IPW, clustered IPW, and grouped GIPW estimators.
groupedestimator <- function(ps, Y, Z, id) {
  
  dataf <- data.frame(Y, Z, ps, id)
  dataf$temp.wt <- temp.wt <- with(dataf, Z/ps + (1-Z)/(1-ps))
  
  IPW_ate = function(IPW, Y, Z) {
    ATE = sum((IPW*Y)[Z==1])/sum(IPW[Z==1]) - sum((IPW*Y)[Z==0])/sum(IPW[Z==0])
    return(ATE)
  }
  
  IPW_se = function(IPW, Y, Z, homo=FALSE) {
    
    if (homo== FALSE) {
      sigsq1 = var(Y[Z==1])
      sigsq0 = var(Y[Z==0])
    } else {
      sigsq1 = var(Y)
      sigsq0 = var(Y) 
    }
    
    SE = sqrt(sigsq1*(sum((IPW[Z==1])^2)/(sum(IPW[Z==1]))^2) + sigsq0*(sum((IPW[Z==0])^2)/(sum(IPW[Z==0]))^2))    
    return(SE)
  }
  
  if (length(unique(id)) == 1) {
    temp.ate = IPW_ate(temp.wt, Y, Z)
    temp.se = IPW_se(temp.wt, Y, Z, homo=FALSE)
    temp.dat = NULL
    
  } else {
    
    temp.dat <- dataf %>% group_by(id) %>% summarize(groupATE=IPW_ate(temp.wt, Y, Z), groupSE=IPW_se(temp.wt, Y, Z, homo=TRUE), wt=sum(temp.wt), .groups = 'drop') 
    temp.ate = sum(temp.dat$groupATE * temp.dat$wt, na.rm = TRUE) / sum(temp.dat$wt[!is.na(temp.dat$groupATE)])
    temp.se = sqrt(sum(temp.dat$wt^2*temp.dat$groupSE^2, na.rm = TRUE))  / sum(temp.dat$wt[is.finite(temp.dat$groupSE)], na.rm = TRUE)
    
  }
  
  return(list(ATE=c(temp.ate, temp.se), groupATE=temp.dat$groupATE, groupSE=temp.dat$groupSE, sumIPW = temp.dat$wt)) # Li's method
}

# write a function to get the weighted mean of groups
GroupWeightedMean <- function(x, w) {
  weighted.mean(x, w)
}

# :: load data
dat <- read.csv("ECLSK_Algebra_complete.csv") # load complete data

dat$WKSESL_s <- as.numeric(scale(dat$WKSESL)) # create the standardized variable of the SES
dat$C6R4MSCL_s <- as.numeric(scale(dat$C6R4MSCL)) # create the standardized variable of the prior achievement scores
dat$C7R4MSCL_s <- as.numeric(scale(dat$C7R4MSCL)) # create the standardized variable of the post-treatment achievement scores

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ::: Within-Group: When the number of groups is equal to 8 ::: ####
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# :: create 8 groups of clusters
set.seed(1)
n.group = 8 # the number of groups = 8
trt.prop = tapply(dat$C7DESMTH, dat$S7_ID, mean)
km.out = kmeans(scale(trt.prop[dat$S7_ID]), n.group, nstart=50)  # run k-means clustering

groupTMLE.est.group <- matrix(NA, nrow=n.group, ncol=3)   

dat$group = km.out$cluster # save the group membership in the data.

# :: implement different estimators: marginal IPW, clustered IPW, and grouped IPW, grouped TMLE estimators with group-specific propensity scores.

# create model matrices
covs.dum_c <- model.matrix(~ C6R4MSCL_s + WKSESL_s + P6EXPECT + GENDER +  RACE + WKPOV_R + WKMOMED + P7HFAMIL+ S7PUPRI + R7REGION + R7URBAN + 0,  data =dat) # create a model matrix with all confounders for tmle
ps.matrix <- model.matrix(~ C6R4MSCL_s +  I(C6R4MSCL_s^2) + WKSESL_s + P6EXPECT + GENDER + RACE + WKPOV_R  + C6R4MSCL_s*WKMOMED +  P7HFAMIL + S7PUPRI + R7REGION + R7URBAN, data = dat) # a model matrix with all confounders 
ps.matrix2 <- model.matrix(~ C6R4MSCL_s +  I(C6R4MSCL_s^2) + WKSESL_s + P6EXPECT + GENDER + RACE + WKPOV_R  + C6R4MSCL_s*WKMOMED +  P7HFAMIL + R7REGION, data = dat) # a model matrix with level-1 confounders only 

groupTMLE.est.group <- matrix(NA, nrow=n.group, ncol=3)   

for (g in 1:n.group) {
  
  # create subsets for group g.
  temp.dat <- dat[dat$group == g,]
  temp.covs.dum_c <- covs.dum_c[dat$group == g,]
  temp.dat$S7_ID <- factor(temp.dat$S7_ID)
  temp.ps.matrix = ps.matrix[dat$group == g,]
  temp.ps.matrix2 = ps.matrix2[dat$group == g,]  
  
  # fit group-specific random effects propensity scores
  temp.group.re <- tryCatch({
    glmer(C7DESMTH ~ temp.ps.matrix +  (1 |S7_ID), data = temp.dat, family=binomial(link="logit"), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5), calc.derivs = FALSE))
  }, error = function(e) {
    print(paste(e))
    glm(C7DESMTH ~ temp.ps.matrix, data = temp.dat, family=binomial)
  })
  
  # fit group-specific fixed effects propensity scores
  temp.group.fe <- glm(C7DESMTH ~ temp.ps.matrix2 + S7_ID, data = temp.dat, family = binomial)
  
  # get predicted values and save them in the data    
  temp.group.re.PS <- predict(temp.group.re, type='response')
  temp.group.fe.PS <- predict(temp.group.fe, type='response')
  dat$group.re.PS[dat$group == g] <- temp.group.re.PS
  dat$group.fe.PS[dat$group == g] <- temp.group.fe.PS
  
  # fit group-specific tmle
  temp.out.TMLE.g <- tmle(Y=temp.dat$C7R4MSCL_s, A=temp.dat$C7DESMTH, W=temp.covs.dum_c) # with default propensity scores
  temp.out.TMLE.g.ps.re <- tmle(Y=temp.dat$C7R4MSCL_s, A=temp.dat$C7DESMTH, W=temp.covs.dum_c, g1W=temp.group.re.PS) # with random effects propensity scores
  temp.out.TMLE.g.ps.fe <- tmle(Y=temp.dat$C7R4MSCL_s, A=temp.dat$C7DESMTH, W=temp.covs.dum_c, g1W=temp.group.fe.PS) # with fixed effects propensity scores
  
  temp.TMLE_ate <- summary(temp.out.TMLE.g)$estimates$ATE$psi
  temp.TMLE_ate.ps.re <- summary(temp.out.TMLE.g.ps.re)$estimates$ATE$psi
  temp.TMLE_ate.ps.fe <- summary(temp.out.TMLE.g.ps.fe)$estimates$ATE$psi
  
  groupTMLE.est.group[g,] <- c(temp.TMLE_ate, temp.TMLE_ate.ps.re, temp.TMLE_ate.ps.fe)
  
  # save info for sensitivity analysis
  dat$groupTMLE.default.PS[dat$group == g] <- temp.out.TMLE.g$g$g1W
  dat$groupTMLE.default.out[dat$group == g] <- ifelse(dat$C7DESMTH[dat$group == g] == 0, temp.out.TMLE.g$Qstar[,1], temp.out.TMLE.g$Qstar[,2])
  dat$groupTMLE.re.out[dat$group == g] <- ifelse(dat$C7DESMTH[dat$group == g] == 0, temp.out.TMLE.g.ps.re$Qstar[,1], temp.out.TMLE.g.ps.re$Qstar[,2])
  dat$groupTMLE.fe.out[dat$group == g] <- ifelse(dat$C7DESMTH[dat$group == g] == 0, temp.out.TMLE.g.ps.fe$Qstar[,1], temp.out.TMLE.g.ps.fe$Qstar[,2])
  
}

# marginal IPW estimator, clustered IPW estimator, grouped IPW estimator
IPW_ATE_SE <- data.frame(rbind(groupedestimator(dat$group.re.PS,  dat$C7R4MSCL_s, dat$C7DESMTH, 1)$ATE, # marginal IPW with random effects propensity scores
                        groupedestimator(dat$group.fe.PS,  dat$C7R4MSCL_s, dat$C7DESMTH, 1)$ATE, # marginal IPW with fixed effects propensity scores
                        groupedestimator(dat$group.re.PS,  dat$C7R4MSCL_s, dat$C7DESMTH, dat$S7_ID)$ATE, # clustered IPW with random effects propensity scores
                        groupedestimator(dat$group.fe.PS,  dat$C7R4MSCL_s, dat$C7DESMTH, dat$S7_ID)$ATE, # clustered IPW with fixed effects propensity scores
                        groupedestimator(dat$group.re.PS,  dat$C7R4MSCL_s, dat$C7DESMTH, dat$group)$ATE, # grouped IPW with random effects propensity scores
                        groupedestimator(dat$group.fe.PS,  dat$C7R4MSCL_s, dat$C7DESMTH, dat$group)$ATE)) # grouped IPW with fixed effects propensity scores

colnames(IPW_ATE_SE) <- c("Estimate", "SE")

# grouped TMLE estimator
TMLE_ATE <- apply(groupTMLE.est.group, 2, function(x) GroupWeightedMean(x, w=table(dat$group)))

# ::::::::::::::::::::::
# ::: Summary of ATE ::: ####
# ::::::::::::::::::::::

# :: summary of results for IPW estimators
IPW_ATE_SE$LCI <- IPW_ATE_SE$Estimate - 1.96*IPW_ATE_SE$SE # lower CI (LCI)
IPW_ATE_SE$UCI <- IPW_ATE_SE$Estimate + 1.96*IPW_ATE_SE$SE # upper CI (UCI)

rownames(IPW_ATE_SE) <- c("MarignalIPW+RePS", "MarignalIPW+FePS", "ClusteredIPW+RePS", "ClusteredIPW+FePS", "GroupedIPW+RePS", "GroupedIPW+FePS")
IPW_ATE_SE  # print results from IPW estimators when K = 8 (i.e., within-cluster approach)

  
# :: summary of results for TMLE estimators
# We note that cluster bootstrap SE and confidence intervals are used for TMLE estimators.
TMLE_ATE_SE <- data.frame(Estimate=TMLE_ATE)
TMLE_ATE_SE$SE <- c(0.03209, 0.03617,	0.03784) # cluster bootstrap SE
TMLE_ATE_SE$LCI <- c(0.08964,	0.05011,	0.05205) # cluster bootstrap LCI
TMLE_ATE_SE$UCI <- c(0.21415,	0.19185,	0.19841) # cluster bootstrap UCI

rownames(TMLE_ATE_SE) <- c("GroupedTMLE+Default", "GroupedTMLE+RePS", "GroupedTMLE+FePS")
TMLE_ATE_SE # print results from TMLE estimators when K = 8 (i.e., within-cluster approach) 

# combine IPW results and TMLE results
ATErlst <- rbind(IPW_ATE_SE, TMLE_ATE_SE)
ATErlst

# ::::::::::::::::::::::::::::
# ::: Sensitivity Analysis ::: ####
# ::::::::::::::::::::::::::::

# extract predictions for the propensity score and the outcome
pred_rlst <- data.frame(groupTMLE.default_trt = dat$group.re.PS, 
                        groupTMLE.RePS_trt= dat$group.re.PS, 
                        groupTMLE.FePS_trt= dat$group.fe.PS,
                        groupTMLE.default_out = dat$groupTMLE.default.out,
                        groupTMLE.RePS_out= dat$groupTMLE.re.out,
                        groupTMLE.FePS_out=  dat$groupTMLE.fe.out)

# compute S^2
S_sq <- apply((dat$C7R4MSCL_s - pred_rlst[, 4:6])^2, 2, mean)/apply((dat$C7DESMTH - pred_rlst[, 1:3])^2, 2, mean) 
S_sq

# compute bias, B
eta_1 <- 0.05
eta_2 <- 0.05
B <- sqrt(S_sq*eta_1*eta_2/(1-eta_2)) 
B

# compute adjusted estimates when positive bias is present.
TMLE_ATE_SE[c("Estimate", "LCI", "UCI")] + B

# compute adjusted estimates when negative bias is present.
TMLE_ATE_SE[c("Estimate", "LCI", "UCI")] - B

