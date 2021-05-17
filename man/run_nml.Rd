\name{run_nml}
\alias{run_nml}
\docType{package}
\title{
Computes NML penalty for multinomial model
}
\description{
Simple function using Rcpp
}
\arguments{
\item{fun}{Model likelihood function. Should have three arguments The first should be the vector of parameter values (which should be unbounded). Second should be the frequency data vector. The third should be a sample size vector of the same length as the data vector (each value indicating the sample size associated with the multinomial distribution of the value in data vector with the same index).}
\item{parl}{Number of free parameters in the model}
\item{ks}{Vector, with each element indicating the number of categories in a multinomial distribution. Note that the length of this vector will indicate the number of multinomial distributions.}
    \item{Ns}{Vector, with each element indicating the sample size of each multinomial distribution. Needs to have the same length as argument ks.}
    \item{fits}{Number of fit runs per dataset (only best fit is retained). Default is 2. NML requires the determination of the maximum likelihood for each dataset considered. Therefore, you should set the number of fits to a number that is very likely to result in the best possible fit.}
    \item{batchsize}{Number of samples used (per core) in each update. Default is 1000.}
    \item{thin}{Number of consecutive samples discarded in a sequence. Default is 1.  Note that batchsize is preserved, which means that by increasing the thinning, one is increasing the total number of data samples taken. For example, if batchsize = 1000 and thin = 5, then 5000 samples are taken (per core), but only 1000 are ultimately retained for fitting.}
    \item{burn}{Scalar, number of initial samples to be taken and discarded before sampling the first batch. For each multinomial distribution, the initial frequency samples are taken from a noninformative Dirichlet-multinomial distribution with parameters alpha = 1.}
\item{cores}{Scalar, number of cores to be used in the computation of the penalty. Default is NULL, which corresponds to 80\% of the cores available.}
    \item{precision}{Scalar, length of the 99.9\% confidence interval of the log-penalty. Default value is 0.20. This value is deemed reasonable given that any penalized-fit differences below 0.5 are pretty much irrelevant for purposes of model selection.}
	\item{G2}{Logical. Whether the model returns the G-squared statistic or negative log-likelihood instead. Default value is TRUE. If set to FALSE, then the function must return a negative log-likelihood.}
	\item{stopcrit}{Character. Which stopping criteria is to be used when determining convergence. The default is the variance between the log-means (betw_logmean).}
	\item{aggmethod}{Character. Which chain aggregation procedure is to be used. The default is to take the log of the average likelihoods (logmean).}
}
