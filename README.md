# A Within-Group Approach to Ensemble Machine Learning Methods for Causal Inference in Multilevel Studies

Youmi Suk

## Overview

Machine learning (ML) methods for causal inference have gained popularity due to their flexibility to predict the outcome model and the propensity score. In this paper, we provide a within-group approach for ML-based causal inference methods to robustly estimate average treatment effects in multilevel studies when there is cluster-level unmeasured confounding. We focus on one particular ML-based causal inference method based on the targeted maximum likelihood estimation (TMLE) with an ensemble learner called SuperLearner. Through our simulation studies, we observe that training TMLE within groups of similar clusters helps remove bias from cluster-level unmeasured confounders. Also, using within-group propensity scores estimated from fixed effects logistic regression increases the robustness of the proposed within-group TMLE method. Even if the propensity scores are partially misspecified, the within-group TMLE still produces robust ATE estimates due to double robustness with flexible modeling, unlike parametric-based inverse propensity weighting methods. We demonstrate our proposed methods and conduct sensitivity analyses against the number of groups and individual-level unmeasured confounding to evaluate the effect of taking an eighth-grade algebra course on math achievement in the Early Childhood Longitudinal Study.

For more details of our proposed methods, see [our paper](https://psyarxiv.com/8s7ut/). 
Here, we provide `R` codes to reproduce our simulation study and replicate our data analysis using the kindergarten cohort of ECLS (ECLS-K) data. 

## Simulation Study

* `DataGeneratingModels.R`  

   This `R` file includes data generating codes for multilevel observational data with cluster-level unmeasured confounding.

* `SimulationCodes.R`
 
   This `R` file includes simulation codes with our proposed modifcations for TMLE where `m.val` controls whether there is a cross-level interaction between a cluster-level unmeasured confounder and a treatment variable. For more information on simulation condtions, see [our paper](https://psyarxiv.com/8s7ut/).


## ECLS-K Data Study

* `ECLSK_Algebra_complete.csv`

  This is our complete data. The original ECLSK 1998-99 dataset is available at https://nces.ed.gov/ecls/dataproducts.asp. For more information on ECLS-K data, see [Tourangeau et al. (2009)](https://nces.ed.gov/pubs2009/2009003.pdf).

* `DataAnalysisCodes.R` 
 
   This `R` file can be used to replicate our data analysis.
   