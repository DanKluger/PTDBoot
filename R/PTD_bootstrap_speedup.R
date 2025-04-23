
#' Speedup to the Predict-Then-Debias Bootstrap
#'
#' This function runs a computational speedup to the Predict-Then-Debias Bootstrap, which corresponds
#' to Algorithm 3 in Kluger et al. (2025). It takes in 3 component datasets and as well as an algorithm that produces a statistical estimate of a quantity of interest and its
#' covariance matrix and runs a quicker version of the Predict-Then-Debias bootstrap.
#' PLEASE NOTE: the 3 input component datasets, `predicted_data_completeSamp`, `predicted_data_incompleteSamp`, and `true_data_completeSamp` should all be dataframes with the exact same column names.
#'
#' @param EstAlgorithm A a function that takes a data sample (called "dfInp"), a  vector of sample weights (called "weightsInp") and a logical "calcVCOV" as inputs.
#' When calcVCOV=FALSE, this function should run the statistical estimation algorithm and return a vector of point estimates for the quantity of interest.
#' When calcVCOV=TRUE, this function should run the statistical estimation algorithm and return a list with a vector of point estimates for the quantity of interest (called "estimate") and its estimated covariance matrix (called "VCOV").
#' If sample weighting is unnecessary, "weightsInp" still needs to be officially defined as an input argument to `EstAlgorithm()`, but it can be ignored or unused in `EstAlgorithm()`.
#' @param true_data_completeSamp An n x p data frame with the gold standard data.
#' @param predicted_data_completeSamp An n x p data frame where all variables that are not widely available are replaced by their widely available proxies
#' @param predicted_data_incompleteSamp An (N-n) x p data where all variables that are not widely available are replaced by their widely available proxies. This dataset corresponds to the samples where gold standard measurement of some variables are unavailable.
#' @param prob_lab_completeSamp An n x 1 vector with the probabilities that each complete sample would have been labelled and therefore assigned to the complete sample  (default is NULL).
#' @param prob_lab_incompleteSamp An (N-n) x 1 vector with the probabilities that each incomplete sample would have been labelled and therefore assigned to the complete sample  (default is NULL).
#' @param B The number of bootstrap draws B (the default is 2000. Larger B cannot hurt, but can lead to slower runtime).
#' @param alpha Statistical significance level. The default is 0.05, but note that the corresponding paper, Kluger et al. (2025), used 0.1.
#' @param TuningScheme Character string of tuning schemes to use. "Optimal" uses an estimate of the optimal tuning matrix."Diagonal" (the default) uses an estimate of the optimal tuning matrix among diagonal tuning matrices. "None" sets the tuning matrix to the identity matrix and corresponds to the untuned Predict-Then-Debias estimator.
#'
#' @return Returns a list with
#' @return  -The Predict-Then-Debias estimate of the regression coefficients.
#' @return  -The (1-alpha) x 100 percent confidence intervals for the regression coefficients (based on the Predict-Then-Debias bootstrap).
#' @return  -The tuning matrix that was ultimately used (this will typically be an estimate of the optimal or optimal diagonal tuning matrix).
#' @importFrom stats as.formula binomial cov gaussian glm poisson quantile rbinom rnorm var
#' @importFrom sandwich sandwich
#' @export
PTD_bootstrap_speedup <- function(EstAlgorithm, true_data_completeSamp, predicted_data_completeSamp, predicted_data_incompleteSamp,B=2000,alpha=0.05,TuningScheme="Diagonal",
                                  prob_lab_completeSamp=NULL,prob_lab_incompleteSamp=NULL){

  ##############  Warning messages ########################################
  if(sum(c(prob_lab_completeSamp,prob_lab_incompleteSamp)==0)>0){print("Warning: the probabilities that each sample
  is labelled should be strictly positive for statistical guarantees")}

  ####################  Computing preliminary variables of interest############################

  #Complete and total sample sizes
  n <- nrow(true_data_completeSamp)
  N <- n + nrow(predicted_data_incompleteSamp)

  #Calculate inverse probability weights
  if( is.null(prob_lab_incompleteSamp) & is.null(prob_lab_completeSamp) ){
    w_complete_samp <- rep(N/n,n) #Default is that all samples get the same weight
    w_incomplete_samp <- rep(N/(N-n),N-n) #Default is that all samples get the same weight
  } else if( (!is.null(prob_lab_incompleteSamp)) & (!is.null(prob_lab_completeSamp)) ){
    w_complete_samp <- 1/prob_lab_completeSamp #Weight set to inverse of probability of being labelled
    w_incomplete_samp <- 1/(1-prob_lab_incompleteSamp) #Weight set to inverse of probability of not being labelled
  } else {
    print("Warning: either set both labelling probabilities to be NULL or both to be nonnull vectors (of lengths n and N-n)")
  }

  #Point estimators on original sample
  hatThetaComplete <- EstAlgorithm(dfInp=true_data_completeSamp, weightsInp=w_complete_samp,calcVCOV=FALSE)
  hatGammaComplete <- EstAlgorithm(dfInp=predicted_data_completeSamp, weightsInp=w_complete_samp,calcVCOV=FALSE)
  hatGammaIncFit <- EstAlgorithm(dfInp=predicted_data_incompleteSamp, weightsInp=w_incomplete_samp,calcVCOV=TRUE)
  hatGammaIncomplete <- hatGammaIncFit$estimate
  CovHatGammaIncEst <- hatGammaIncFit$VCOV

  ################################### Calculating Bootstrapped Component Estimators ########################################

  hatThetaComplete_Boot <- matrix(NA, nrow = B, ncol = length(hatThetaComplete))
  hatGammaComplete_Boot <- matrix(NA, nrow = B, ncol = length(hatGammaComplete))
  n_completeBoot <- rbinom(B,size = N,prob = n/N)
  for(b in 1:B){
    idx_complete_b <- sample(1:n,size = n_completeBoot[b],replace = T)
    hatThetaComplete_Boot[b,] <-  EstAlgorithm(dfInp=true_data_completeSamp[idx_complete_b,], weightsInp=w_complete_samp[idx_complete_b],calcVCOV=FALSE)
    hatGammaComplete_Boot[b,] <-  EstAlgorithm(dfInp=predicted_data_completeSamp[idx_complete_b,], weightsInp=w_complete_samp[idx_complete_b],calcVCOV=FALSE)
  }

  #Draw hatGammaIncomplete from its asymptotic approximation
  IIDGaussianNoise <- matrix(rnorm(B*length(hatGammaIncomplete)),nrow = length(hatGammaIncomplete),ncol = B)
  GaussianNoise <- t(chol(CovHatGammaIncEst)) %*% IIDGaussianNoise
  hatGammaIncompleteMat <- matrix(rep(hatGammaIncomplete,B),ncol = B,byrow = F)
  hatGammaIncomplete_ApproxBoot <- t(hatGammaIncompleteMat+GaussianNoise)

  ################################### Calculating Optimal Tuning Matrix ###################################################
  CovHat_theta_gamma_comp <- cov(hatThetaComplete_Boot,hatGammaComplete_Boot)
  VarHat_gamma_comp  <- var(hatGammaComplete_Boot)

  if(TuningScheme=="Optimal"){
    hatOmega <- CovHat_theta_gamma_comp %*% solve(VarHat_gamma_comp+CovHatGammaIncEst)
  } else if(TuningScheme=="Diagonal"){
    hatOmega <- diag(diag(CovHat_theta_gamma_comp)/diag(VarHat_gamma_comp+CovHatGammaIncEst))
  }else if(TuningScheme=="None"){
    hatOmega <- diag(length(hatGammaComplete)) #Set to identity matrix
  }



  ################################### Calculating PTD point estimate and CIs #################################################################
  #PTD point estimate
  PTD_estimate <- as.vector(hatThetaComplete + hatOmega %*% (hatGammaIncomplete-hatGammaComplete))
  #Bootstrap draws of PTD estimator
  PTD_estimator_boot <- hatThetaComplete_Boot + t(hatOmega %*% (t(hatGammaIncomplete_ApproxBoot-hatGammaComplete_Boot)))
  #Percentile Bootstrap Confidence Intervals
  PTD_Boot_CIs <-  t(sapply(data.frame(PTD_estimator_boot), FUN = function(x) quantile(x,probs = c(alpha/2,1-alpha/2))))

  #Assigning names to outputs
  colnames(hatOmega) <- rownames(hatOmega)  <- names(hatThetaComplete)
  names(PTD_estimate) <- names(hatThetaComplete)
  rownames(PTD_Boot_CIs) <- names(hatThetaComplete)
  colnames(PTD_Boot_CIs) <- c("lower","upper")

  return(list(PTD_estimate=PTD_estimate,PTD_Boot_CIs=PTD_Boot_CIs,Tuning_matrix=hatOmega))
}
