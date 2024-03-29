% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pval_infl.R
\name{pval_infl}
\alias{pval_infl}
\title{Calculate inflation factor from p-values}
\usage{
pval_infl(pvals, df = 1)
}
\arguments{
\item{pvals}{The vector of association p-values to analyze.
This function assumes all p-values are provided (a mix of null and alternative tests).
\code{NA} values are allowed in input and removed.
Non-\code{NA} values outside of [0, 1] will trigger an error.}

\item{df}{The degrees of freedom of the assumed chi-squared distribution (default 1).}
}
\value{
The inflation factor
}
\description{
The inflation factor is defined as the median association test statistic divided by the expected median under the null hypothesis, which is typically assumed to have a chi-squared distribution.
This function takes a p-value distribution and maps its median back to the chi-squared value (using the quantile function) in order to compute the inflation factor in the chi-squared scale.
The full p-value distribution (a mix of null and alternative cases) is used to calculate the desired median value (the true \code{causal_loci} is not needed, unlike \code{\link[=pval_srmsd]{pval_srmsd()}}).
}
\examples{
# simulate truly null p-values, which should be uniform
pvals <- runif(10)
# calculate desired measure
pval_infl( pvals )

}
\seealso{
\code{\link[=pval_srmsd]{pval_srmsd()}}, a more robust measure of null p-value accuracy, but which requires knowing the true causal loci.

\code{\link[=pval_type_1_err]{pval_type_1_err()}} for classical type I error rate estimates.
}
