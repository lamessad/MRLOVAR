#' Inverse-Variance Weighted (IVW) Estimation
#'
#' @description This function performs Inverse-Variance Weighted (IVW) estimation of causal effect using linear regression. It estimates the causal effect of an exposure on an outcome by weighting the observations inversely proportional to their variance. The function returns the estimate of the causal effect, its standard error, and p-value.
#'
#' @param betaY Numeric vector of effect sizes for the outcome.
#' @param betaX Numeric vector of effect sizes for the exposure.
#' @param betaYse Numeric vector of standard errors for the outcome.
#' @param ld_cor Numeric matrix for the LD correlation among SNPs.
#'
#' @details The function uses linear regression to estimate the causal effect of the exposure on the outcome, with weights based on the inverse of the variance of the outcome's effect sizes. The standard error is adjusted based on the model's residual standard error. The p-value is computed from the normal distribution of the test statistic.
#'
#' @importFrom stats pnorm
#'
#' @return A named vector with the following elements:
#' \item{Estimate}{The estimate of the causal effect.}
#' \item{SE}{The standard error of the estimate.}
#' \item{P}{The p-value for the test.}
#'
IVW <- function(betaY, betaX, betaYse, ld_cor) {

  nsnps <- length(betaY)
  omega <- betaYse %o% betaYse * ld_cor
  thetaIVW <- as.numeric(solve(t(betaX) %*% solve(omega) %*% betaX) %*% (t(betaX) %*% solve(omega) %*% betaY))
  rse <- betaY - thetaIVW * betaX
  thetaIVWse <- sqrt(solve(t(betaX) %*% solve(omega) %*% betaX)) * max(sqrt(t(rse) %*% solve(omega) %*% rse / (nsnps - 1)), 1)

  pvalue <- 2 * pnorm(-abs(thetaIVW / thetaIVWse))

  return(list(
    Estimate = thetaIVW,
    SE = thetaIVWse,
    P = pvalue))
}
