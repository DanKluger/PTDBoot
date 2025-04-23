#Algorithms 6 from Kluger et al. (2025), 'Prediction-Powered Inference with Imputed Covariates and Nonuniform Sampling'
#arXiv URL: https://arxiv.org/abs/2501.18577.



##########################################################


#' Stratified Predict-Then-Debias Bootstrap
#'
#' This function runs a stratified bootstrap modification to the Predict-Then-Debias Bootstrap, which corresponds
#' to Algorithm 6 in Kluger et al. (2025). The algorithm is designed for settings where there is a small number of large strata. It takes in 3 component datasets, 2 strata ID vectors, as well as an algorithm that produces a statistical estimate of a quantity of interest runs a modification of the Predict-Then-Debias bootstrap where within each strata complete and incomplete samples are both resampled with replacement.
#' PLEASE NOTE: the 3 input component datasets, `predicted_data_completeSamp`, `predicted_data_incompleteSamp`, and `true_data_completeSamp` should all be dataframes with the exact same column names.
#'
#'
#' @param EstAlgorithm A function that takes a data sample (called "dfInp") and vector of sample weights (called "weightsInp") as inputs, runs a statistical estimation algorithm and subsequently
#'  returns a vector of estimates of the quantities of interest as outputs. If sample weighting is unnecessary, "weightsInp" still needs to be officially defined as an input argument to `EstAlgorithm()`, but it can be ignored or unused in `EstAlgorithm()`.
#' @param true_data_completeSamp An n x p data frame with the gold standard data.
#' @param predicted_data_completeSamp An n x p data frame where all variables that are not widely available are replaced by their widely available proxies
#' @param predicted_data_incompleteSamp An (N-n) x p data where all variables that are not widely available are replaced by their widely available proxies. This dataset corresponds to the samples where gold standard measurement of some variables are unavailable.
#' @param StrataID_completeSamp An n x 1 vector of integers or character strings giving the corresponding Strata IDs for each of the n complete samples.
#' @param StrataID_incompleteSamp An (N-n) x 1 vector of integers or character strings giving the corresponding Strata IDs for each of the N-n incomplete samples.
#' @param TotalStrataSizes A named vector with the number of elements in each strata from the large population from which the stratified sample was obtained. The names should correspond to the unique strata IDs and the values should correspond to the overall strata sizes.
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
PTD_stratified_bootstrap <- function(EstAlgorithm, true_data_completeSamp, predicted_data_completeSamp, predicted_data_incompleteSamp,
                                     StrataID_completeSamp,StrataID_incompleteSamp,TotalStrataSizes,
                                  B=2000,alpha=0.05,TuningScheme="Diagonal"){

  print("Warning: Please note that this method is designed for settings where there is a small number of large strata.")

  ####################  Computing preliminary variables of interest############################

  #Complete and total sample sizes
  n <- nrow(true_data_completeSamp)
  N <- n + nrow(predicted_data_incompleteSamp)



  #Formatting strata IDs and counts for stratified sampling and calculating inverse probability weights
  uniqueStrataIDs <- unique(c(StrataID_completeSamp,StrataID_incompleteSamp))
  idxCompleteInEachStrata <- list()
  idxIncompleteInEachStrata <- list()

  w_complete_samp <- rep(NA,n)
  w_incomplete_samp <- rep(NA,N-n)
  for(j in 1:length(uniqueStrataIDs)){
    #Storing indices corresponding to each strata
    idxCompleteInEachStrata[[j]] <- which(StrataID_completeSamp==uniqueStrataIDs[j])
    idxIncompleteInEachStrata[[j]] <- which(StrataID_incompleteSamp==uniqueStrataIDs[j])
    #calculate probability that an element in strata j is both sampled and complete
    n_complete_j <- length(idxCompleteInEachStrata[[j]])
    strat_size_j <- TotalStrataSizes[names(TotalStrataSizes)==uniqueStrataIDs[j]]
    prob_samp_complete_j <- n_complete_j/strat_size_j
    #calculate probability that an element in strata j is both sampled and incomplete
    n_inc_j <- length(idxIncompleteInEachStrata[[j]])
    prob_samp_incomplete_j <- n_inc_j/strat_size_j

    #setting inverse probability weights for samples in strata j
    w_complete_samp[idxCompleteInEachStrata[[j]]] <- 1/prob_samp_complete_j
    w_incomplete_samp[idxIncompleteInEachStrata[[j]]] <- 1/prob_samp_incomplete_j
  }

  #Point estimators on original sample
  hatThetaComplete <- EstAlgorithm(dfInp=true_data_completeSamp, weightsInp=w_complete_samp)
  hatGammaComplete <- EstAlgorithm(dfInp=predicted_data_completeSamp, weightsInp=w_complete_samp)
  hatGammaIncomplete <- EstAlgorithm(dfInp=predicted_data_incompleteSamp, weightsInp=w_incomplete_samp)

  ################################### Calculating Bootstrapped Component Estimators ########################################

  hatThetaComplete_Boot <- matrix(NA, nrow = B, ncol = length(hatThetaComplete))
  hatGammaComplete_Boot <- matrix(NA, nrow = B, ncol = length(hatGammaComplete))
  hatGammaIncomplete_Boot <- matrix(NA, nrow = B, ncol = length(hatGammaIncomplete))

  for(b in 1:B){

    #Resample within each strata with replacement
    idx_complete_b <- vector()
    idx_inc_b <- vector()
    for(j in 1:length(uniqueStrataIDs)){
      clust_j_complete <- sample(x = idxCompleteInEachStrata[[j]],size = length(idxCompleteInEachStrata[[j]]),replace = T)
      idx_complete_b <- c(idx_complete_b,clust_j_complete)

      clust_j_inc <- sample(x = idxIncompleteInEachStrata[[j]],size = length(idxIncompleteInEachStrata[[j]]),replace = T)
      idx_inc_b <- c(idx_inc_b,clust_j_inc)
    }

    hatThetaComplete_Boot[b,] <-  EstAlgorithm(dfInp=true_data_completeSamp[idx_complete_b,], weightsInp=w_complete_samp[idx_complete_b])
    hatGammaComplete_Boot[b,] <-  EstAlgorithm(dfInp=predicted_data_completeSamp[idx_complete_b,], weightsInp=w_complete_samp[idx_complete_b])
    hatGammaIncomplete_Boot[b,] <-  EstAlgorithm(dfInp=predicted_data_incompleteSamp[idx_inc_b,], weightsInp=w_incomplete_samp[idx_inc_b])
  }

  ################################### Calculating Optimal Tuning Matrix ###################################################
  CovHat_theta_gamma_comp <- cov(hatThetaComplete_Boot,hatGammaComplete_Boot)
  VarHat_gamma_comp  <- var(hatGammaComplete_Boot)
  VarHat_gamma_incomplete <- var(hatGammaIncomplete_Boot)

  if(TuningScheme=="Optimal"){
    hatOmega <- CovHat_theta_gamma_comp %*% solve(VarHat_gamma_comp+VarHat_gamma_incomplete)
  } else if(TuningScheme=="Diagonal"){
    hatOmega <- diag(diag(CovHat_theta_gamma_comp)/diag(VarHat_gamma_comp+VarHat_gamma_incomplete))
  }else if(TuningScheme=="None"){
    hatOmega <- diag(length(hatGammaComplete)) #Set to identity matrix
  }

  ################################### Calculating PTD point estimator and CIs #################################################################
  #PTD point estimator
  PTD_estimate <- as.vector(hatThetaComplete + hatOmega %*% (hatGammaIncomplete-hatGammaComplete))
  #Bootstrap draws of PTD estimator
  PTD_estimator_boot <- hatThetaComplete_Boot + t(hatOmega %*% (t(hatGammaIncomplete_Boot-hatGammaComplete_Boot)))
  #Percentile Bootstrap Confidence Intervals
  PTD_Boot_CIs <-  t(sapply(data.frame(PTD_estimator_boot), FUN = function(x) quantile(x,probs = c(alpha/2,1-alpha/2))))

  #Assigning names to outputs
  colnames(hatOmega) <- rownames(hatOmega)  <- names(hatThetaComplete)
  names(PTD_estimate) <- names(hatThetaComplete)
  rownames(PTD_Boot_CIs) <- names(hatThetaComplete)
  colnames(PTD_Boot_CIs) <- c("lower","upper")

  return(list(PTD_estimate=PTD_estimate,PTD_Boot_CIs=PTD_Boot_CIs,Tuning_matrix=hatOmega))
}
