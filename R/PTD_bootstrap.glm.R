
#' Predict-Then-Debias bootstrap for generalized linear models
#'
#'   Implementation of the Predict-Then-Debias bootstrap for settings where
#'   the estimand is the regression coefficient in a generalized linear model (GLM).
#'   This function calls more general functions such as \code{PTD_bootstrap()}
#'   and does not require the user to specify a function as an input argument. PLEASE NOTE: the first 3 arguments, predicted_data_completeSamp, predicted_data_incompleteSamp, and true_data_completeSamp should all be dataframes with the exact same column names.
#'
#' @param true_data_completeSamp An n x p data frame with the gold standard data.
#' @param predicted_data_completeSamp An n x p data frame where all variables that are not widely available are replaced by their widely available proxies
#' @param predicted_data_incompleteSamp An (N-n) x p data where all variables that are not widely available are replaced by their widely available proxies. This dataset corresponds to the samples where gold standard measurement of some variables are unavailable.
#' @param regFormula.glm A character string giving the regression specification. Standard formatting that is used for the \code{lm()} and \code{glm()} functions in R should be used.
#' @param GLM_type A character string that should be set to either "linear", "logistic", or "Poisson" to indicate whether the estimands of interest are the regression coefficients of a linear regression model, logistic regression model or a Poisson regression model, respectively.
#' @param B The number of bootstrap draws B (the default is 2000. Larger B cannot hurt, but can lead to slower runtime).
#' @param alpha Statistical significance level. The default is 0.05, but note that the corresponding paper, Kluger et al. (2025), used 0.1.
#' @param TuningScheme Character string of tuning schemes to use. "Optimal" uses an estimate of the optimal tuning matrix."Diagonal" (the default) uses an estimate of the optimal tuning matrix among diagonal tuning matrices. "None" sets the tuning matrix to the identity matrix and corresponds to the untuned Predict-Then-Debias estimator.
#' @param prob_lab_completeSamp An n x 1 vector with the probabilities that each complete sample would have been labelled and therefore assigned to the complete sample  (default is NULL).
#' @param prob_lab_incompleteSamp An (N-n) x 1 vector with the probabilities that each incomplete sample would have been labelled and therefore assigned to the complete sample  (default is NULL).
#' @param GroupID_completeSamp  An n x 1 vector of integers or character strings giving the corresponding Cluster IDs (if clustered=TRUE) or the Strata IDs (if stratified=TRUE) for each of the n complete samples. When clustered=FALSE and stratified=FALSE, this variable is ignored.
#' @param GroupID_incompleteSamp An (N-n) x 1 vector of integers or character strings giving the corresponding Cluster IDs (if clustered=TRUE) or the Strata IDs (if stratified=TRUE) for each of the N-n incomplete samples. When clustered=FALSE and stratified=FALSE, this variable is ignored.
#' @param TotalStrataSizes A named vector with the number of elements in each strata from the large population from which the stratified sample was obtained. The names should correspond to the unique strata IDs and the values should correspond to the overall strata sizes. Only need to set this when stratified=TRUE.
#' @param clustered A TRUE/FALSE logical of whether to use the clustered Predict-Then-Debias Bootstrap (Algorithm 5 in  Kluger et al. (2025)) to account for data that was labelled according to cluster sampling.
#' @param stratified A TRUE/FALSE logical of whether to use the stratified Predict-Then-Debias Bootstrap (Algorithm 6 in  Kluger et al. (2025)) to account for data that was collected according to stratified sampling (with a small number of large strata).
#' @param speedup A TRUE/FALSE logical of whether to speedup the bootstrap procedure using Algorithm 3 in Kluger et al. (2025).
#'
#' @return Returns a list with
#' @return  -The Predict-Then-Debias estimate of the regression coefficients.
#' @return  -The (1-alpha) x 100 percent confidence intervals for the regression coefficients (based on the Predict-Then-Debias bootstrap).
#' @return  -The tuning matrix that was ultimately used (this will typically be an estimate of the optimal or optimal diagonal tuning matrix).
#' @importFrom stats as.formula binomial cov gaussian glm poisson quantile rbinom rnorm var
#' @importFrom sandwich sandwich
#' @export
PTD_bootstrap.glm <- function(true_data_completeSamp, predicted_data_completeSamp, predicted_data_incompleteSamp,
                              regFormula.glm,GLM_type,B=2000,alpha=0.05,TuningScheme="Diagonal",
                              clustered=FALSE,stratified=FALSE,speedup=FALSE,
                              prob_lab_completeSamp=NULL,prob_lab_incompleteSamp=NULL,
                              GroupID_completeSamp=NULL,GroupID_incompleteSamp=NULL,TotalStrataSizes=NULL){

  ########################### Defining a function that runs the specified GLM ########################


  EstAlgorithmInp <- function(dfInp,weightsInp,calcVCOV=FALSE){

    if(GLM_type %in% c("linear","Linear","Gaussian","gaussian","Normal","normal","OLS","ols")){
      familyLocal <- gaussian(link = "identity")
    } else if(GLM_type %in% c("logistic","Logistic","logit")){
      familyLocal <- binomial(link = "logit")
    } else if(GLM_type %in% c("Poisson","poisson")){
      familyLocal <- poisson(link = "log")
    }

    if(GLM_type %in% c("logistic","Logistic","logit")){
      #Weighted logistic regression software returns warnings when the weights are non-integers
      #These warnings can be ignored for our purposes
      glmFit <- suppressWarnings(
        glm(formula = as.formula(regFormula.glm),family = familyLocal,data = dfInp,weights = weightsInp)
      )
    } else{
      glmFit <- glm(formula = as.formula(regFormula.glm),family = familyLocal,data = dfInp,weights = weightsInp)
    }
    if(calcVCOV){
      VCOVSandwich <- sandwich(glmFit)
      return(list(estimate=glmFit$coefficients,VCOV=VCOVSandwich))
    }else{
      return(glmFit$coefficients)
    }
  }

  ################### Running the specified variant of the Predict-Then-Debias bootstrap ####################

  if(clustered){
    #Clustered Predict-Then-Debias bootstrap (Algorithm 5 from Kluger et al. (2025))

    if(speedup){print("Speedup has not been implemented the clustered Predict-Then-Debias bootstrap in this package, so the standard version is being run. A speedup is possible though, see the reproducibility repo or the remarks at the end of Section B.1 in Kluger et al. (2025)")}

    out <- PTD_cluster_bootstrap(EstAlgorithmInp, true_data_completeSamp, predicted_data_completeSamp, predicted_data_incompleteSamp,
                                 prob_lab_completeSamp=prob_lab_completeSamp,prob_lab_incompleteSamp=prob_lab_incompleteSamp,
                                 clusterID_completeSamp=GroupID_completeSamp,clusterID_incompleteSamp=GroupID_incompleteSamp,B=B,alpha=alpha,TuningScheme=TuningScheme)
  } else if(stratified){
    #Stratified Predict-Then-Debias bootstrap (Algorithm 6 from Kluger et al. (2025))

    if((!is.null(prob_lab_completeSamp))| (!is.null(prob_lab_completeSamp))){
      print("Warning: The user specified labelling probabilities are being ignored. The stratified Predict-Then-Debias bootstrap automatically calculates inverse probability weights based on number complete and incomplete samples in each strata.")
    }
    if(speedup){print("A speedup has not been implemented the stratified Predict-Then-Debias bootstrao. Running the standard version.")}

    out <- PTD_stratified_bootstrap(EstAlgorithmInp, true_data_completeSamp, predicted_data_completeSamp, predicted_data_incompleteSamp,
                                    StrataID_completeSamp=GroupID_completeSamp,StrataID_incompleteSamp=GroupID_incompleteSamp,
                                    B=B,alpha=alpha,TuningScheme=TuningScheme,TotalStrataSizes=TotalStrataSizes)

  }else if(speedup){

    #Speedup of Predict-Then-Debias bootstrap (Algorithm 3 from Kluger et al. (2025))
    out <- PTD_bootstrap_speedup(EstAlgorithmInp, true_data_completeSamp, predicted_data_completeSamp, predicted_data_incompleteSamp,
                                 prob_lab_completeSamp=prob_lab_completeSamp,prob_lab_incompleteSamp=prob_lab_incompleteSamp,B=B,alpha=alpha,TuningScheme=TuningScheme)
  }else{

    #Standard Predict-Then-Debias bootstrap (Algorithm 1 or 2 from Kluger et al. (2025))
    out <- PTD_bootstrap(EstAlgorithmInp, true_data_completeSamp, predicted_data_completeSamp, predicted_data_incompleteSamp,
                         prob_lab_completeSamp=prob_lab_completeSamp,prob_lab_incompleteSamp=prob_lab_incompleteSamp,B=B,alpha=alpha,TuningScheme=TuningScheme)
  }

  ## Add cluster bootstrap, and stratified bootstrap

  return(out)
}
