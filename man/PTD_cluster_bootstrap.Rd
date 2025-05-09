% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PTD_cluster_bootstrap.R
\name{PTD_cluster_bootstrap}
\alias{PTD_cluster_bootstrap}
\title{Clustered Predict-Then-Debias Bootstrap}
\usage{
PTD_cluster_bootstrap(
  EstAlgorithm,
  true_data_completeSamp,
  predicted_data_completeSamp,
  predicted_data_incompleteSamp,
  clusterID_completeSamp,
  clusterID_incompleteSamp,
  B = 2000,
  alpha = 0.05,
  TuningScheme = "Diagonal",
  prob_lab_completeSamp = NULL,
  prob_lab_incompleteSamp = NULL
)
}
\arguments{
\item{EstAlgorithm}{A function that takes a data sample (called "dfInp") and vector of sample weights (called "weightsInp") as inputs, runs a statistical estimation algorithm and subsequently
returns a vector of estimates of the quantities of interest as outputs. If sample weighting is unnecessary, "weightsInp" still needs to be officially defined as an input argument to \code{EstAlgorithm()}, but it can be ignored or unused in \code{EstAlgorithm()}.}

\item{true_data_completeSamp}{An n x p data frame with the gold standard data.}

\item{predicted_data_completeSamp}{An n x p data frame where all variables that are not widely available are replaced by their widely available proxies}

\item{predicted_data_incompleteSamp}{An (N-n) x p data where all variables that are not widely available are replaced by their widely available proxies. This dataset corresponds to the samples where gold standard measurement of some variables are unavailable.}

\item{clusterID_completeSamp}{An n x 1 vector of integers or character strings giving the corresponding Cluster IDs for each of the n complete samples.}

\item{clusterID_incompleteSamp}{An (N-n) x 1 vector of integers or character strings giving the corresponding Cluster IDs for each of the N-n incomplete samples.}

\item{B}{The number of bootstrap draws B (the default is 2000. Larger B cannot hurt, but can lead to slower runtime).}

\item{alpha}{Statistical significance level. The default is 0.05, but note that the corresponding paper, Kluger et al. (2025), used 0.1.}

\item{TuningScheme}{Character string of tuning schemes to use. "Optimal" uses an estimate of the optimal tuning matrix."Diagonal" (the default) uses an estimate of the optimal tuning matrix among diagonal tuning matrices. "None" sets the tuning matrix to the identity matrix and corresponds to the untuned Predict-Then-Debias estimator.}

\item{prob_lab_completeSamp}{An n x 1 vector with the probabilities that each complete sample would have been labelled and therefore assigned to the complete sample  (default is NULL).}

\item{prob_lab_incompleteSamp}{An (N-n) x 1 vector with the probabilities that each incomplete sample would have been labelled and therefore assigned to the complete sample  (default is NULL).}
}
\value{
Returns a list with

-The Predict-Then-Debias estimate of the regression coefficients.

-The (1-alpha) x 100 percent confidence intervals for the regression coefficients (based on the Predict-Then-Debias bootstrap).

-The tuning matrix that was ultimately used (this will typically be an estimate of the optimal or optimal diagonal tuning matrix).
}
\description{
This function runs a cluster  bootstrap modification to the Predict-Then-Debias Bootstrap, which corresponds
to Algorithm 5 in Kluger et al. (2025). It takes in 3 component datasets, 2 cluster ID vectors, as well as an algorithm that produces a statistical estimate of a quantity of interest runs a modification of the Predict-Then-Debias bootstrap where clusters (as opposed to individual samples) are resampled with replacement.
PLEASE NOTE: the 3 input component datasets, \code{predicted_data_completeSamp}, \code{predicted_data_incompleteSamp}, and \code{true_data_completeSamp} should all be dataframes with the exact same column names.
}
