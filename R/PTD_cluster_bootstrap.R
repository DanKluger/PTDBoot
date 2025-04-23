
#' Clustered Predict-Then-Debias Bootstrap
#'
#' This function runs a cluster  bootstrap modification to the Predict-Then-Debias Bootstrap, which corresponds
#' to Algorithm 5 in Kluger et al. (2025). It takes in 3 component datasets, 2 cluster ID vectors, as well as an algorithm that produces a statistical estimate of a quantity of interest runs a modification of the Predict-Then-Debias bootstrap where clusters (as opposed to individual samples) are resampled with replacement.
#' PLEASE NOTE: the 3 input component datasets, `predicted_data_completeSamp`, `predicted_data_incompleteSamp`, and `true_data_completeSamp` should all be dataframes with the exact same column names.
#'
#'
#' @param EstAlgorithm A function that takes a data sample (called "dfInp") and vector of sample weights (called "weightsInp") as inputs, runs a statistical estimation algorithm and subsequently
#'  returns a vector of estimates of the quantities of interest as outputs. If sample weighting is unnecessary, "weightsInp" still needs to be officially defined as an input argument to `EstAlgorithm()`, but it can be ignored or unused in `EstAlgorithm()`.
#' @param true_data_completeSamp An n x p data frame with the gold standard data.
#' @param predicted_data_completeSamp An n x p data frame where all variables that are not widely available are replaced by their widely available proxies
#' @param predicted_data_incompleteSamp An (N-n) x p data where all variables that are not widely available are replaced by their widely available proxies. This dataset corresponds to the samples where gold standard measurement of some variables are unavailable.
#' @param clusterID_completeSamp An n x 1 vector of integers or character strings giving the corresponding Cluster IDs for each of the n complete samples.
#' @param clusterID_incompleteSamp An (N-n) x 1 vector of integers or character strings giving the corresponding Cluster IDs for each of the N-n incomplete samples.
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

PTD_cluster_bootstrap <- function(EstAlgorithm, true_data_completeSamp, predicted_data_completeSamp, predicted_data_incompleteSamp,
                                  clusterID_completeSamp,clusterID_incompleteSamp,
                                  B=2000,alpha=0.05,TuningScheme="Diagonal",prob_lab_completeSamp=NULL,prob_lab_incompleteSamp=NULL){

  ##############  Warning messages ########################################
  if(sum(c(prob_lab_completeSamp,prob_lab_incompleteSamp)==0)>0){print("Warning: the probabilities that each sample
  is labelled should be strictly positive for statistical guarantees")}

  if(length(intersect(clusterID_completeSamp,clusterID_incompleteSamp))>0){
    print("Warning: some clusters contain a mixture of complete samples and incomplete samples. The Predict-Then-Debias cluster bootstrap will still run but it has not been tested and its validity has not been verified in such settings")
    }

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

  #Formatting clusters for cluster sampling
  uniqueClusterIDs <- unique(c(clusterID_completeSamp,clusterID_incompleteSamp))
  idxCompleteInEachCluster <- list()
  idxIncompleteInEachCluster <- list()

  for(j in 1:length(uniqueClusterIDs)){
    idxCompleteInEachCluster[[j]] <- which(clusterID_completeSamp==uniqueClusterIDs[j])
    idxIncompleteInEachCluster[[j]] <- which(clusterID_incompleteSamp==uniqueClusterIDs[j])
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

    #Resample clusters with replacement
    ClustersResamp <- sample(1:length(uniqueClusterIDs),size = length(uniqueClusterIDs),replace = T)
    idx_complete_b <- vector()
    idx_inc_b <- vector()
    for(j in 1:length(ClustersResamp)){
      idx_complete_b <- c(idx_complete_b,idxCompleteInEachCluster[[ClustersResamp[j]]])
      idx_inc_b <- c(idx_inc_b,idxIncompleteInEachCluster[[ClustersResamp[j]]])
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
