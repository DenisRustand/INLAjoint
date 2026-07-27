#' Compare INLAjoint models by pairwise WAIC differences
#'
#' @description
#' Pairwise comparison of two or more fitted joint models based on the
#' widely applicable information criterion (Watanabe, 2010). For each
#' pair, the WAIC difference is computed from the pointwise WAIC
#' contributions, which provides a standard error for the difference and an
#' approximate z-test of equal predictive accuracy.
#'
#' @details
#' Since WAIC decomposes as a sum over observations, the difference between
#' two models can be computed observation by observation. For a pair of
#' models with pointwise contributions \eqn{w_i^{(1)}}{w_i(1)} and
#' \eqn{w_i^{(2)}}{w_i(2)}, the function reports
#' \eqn{\Delta = \sum_i (w_i^{(1)} - w_i^{(2)})}{dWAIC = sum_i (w_i(1) - w_i(2))}
#' with standard error
#' \eqn{\sqrt{n}\, \mathrm{sd}(w_i^{(1)} - w_i^{(2)})}{sqrt(n) * sd(w_i(1) - w_i(2))}.
#' This paired construction accounts for the correlation between the two
#' model' contributions for the same observations and follows the standard-error
#' recommendation of Vehtari, Gelman and Gabry (2017) for pointwise
#' predictive criteria. The standardized statistic
#' \eqn{z = \Delta / \mathrm{se}}{z = dWAIC / se} is interpreted using
#' a two-sided normal approximation.
#'
#' A positive \code{dWAIC} for a pair (\code{model1}, \code{model2}) means
#' \code{model2} has the lower WAIC (better estimated predictive accuracy).
#' The \code{better} column always names the preferred model, whether or not
#' the difference is significant.
#'
#' Although the compared models are Bayesian, WAIC is a point estimate of
#' out-of-sample predictive accuracy, and the z-test assesses the
#' sampling uncertainty of that estimate. It is a frequentist
#' evaluation of predictive performance, not a posterior probability that
#' one model is better.
#'
#' @section Limitations:
#' The comparison is most informative when \code{dWAIC} is large relative
#' to its standard error. When models have similar predictive performance,
#' the normal approximation can be poorly calibrated and the standard error
#' too small (Sivula et al., 2025). Small, marginally significant differences
#' should therefore not be taken as clear evidence of superiority.
#'
#' The standard error assumes independent pointwise differences. With repeated
#' measurements, these may remain correlated within subjects, making the
#' standard error optimistic. Pairwise differencing can reduce this problem
#' by cancelling variation shared by both models, but may not remove it.
#'
#' @param ... Two or more objects of class \code{INLAjoint}. Models must be
#'   fitted to the same data so that contributions are aligned. Observations
#'   with a missing contribution in a pair are dropped from the comparison.
#'
#' @return A \code{data.frame} with one row per pair of models and columns:
#' \describe{
#'   \item{model1, model2}{names of the compared models.}
#'   \item{dWAIC}{summed pointwise WAIC difference (model1 minus model2);
#'     positive values favor \code{model2}.}
#'   \item{se}{standard error of \code{dWAIC}.}
#'   \item{z}{\code{dWAIC / se}.}
#'   \item{p}{two-sided p-value from the normal approximation.}
#'   \item{better}{name of the model with the lower WAIC.}
#'   \item{sig}{significance codes: \code{***} for p < 0.001, \code{**}
#'   for p < 0.01, \code{*} for p < 0.05.}
#' }
#'
#' @references
#' Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC.
#' \emph{Statistics and Computing}, 27(5), 1413--1432.
#'
#' Sivula, T., Magnusson, M., Matamoros, A. A. and Vehtari, A. (2025).
#' Uncertainty in Bayesian leave-one-out cross-validation based model
#' comparison. \emph{Bayesian Analysis}.
#' \doi{10.1214/25-BA1569}.
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and
#' widely applicable information criterion in singular learning theory.
#' \emph{Journal of Machine Learning Research}, 11, 3571--3594.
#'
#' @examples
#' \dontrun{
#' # m1, m2, m3 fitted with joint() on the same data
#' joint.compare(m1, m2, m3)
#' }
#'
#' @export

joint.compare <- function(...) {
  models <- list(...)
  nm <- as.character(match.call()[-1])
  k <- length(models)
  if (k < 2) stop("need at least two models")
  for (i in seq_len(k)) {
    if (!inherits(models[[i]], "INLAjoint"))
      stop(paste0("argument ", i, " is not an INLAjoint object"))
    if (is.null(models[[i]]$waic_vec))
      stop(paste0("model '", nm[i], "' has no waic_vec field"))
  }
  pairs <- expand.grid(m1 = seq_len(k), m2 = seq_len(k))
  pairs <- pairs[pairs$m1 < pairs$m2, ]
  out <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(p) {
    i <- pairs$m1[p]
    j <- pairs$m2[p]
    d <- models[[i]]$waic_vec - models[[j]]$waic_vec
    dw <- sum(na.omit(d))
    se <- sqrt(length(na.omit(d))) * sd(na.omit(d))
    z <- dw / se
    p <- 2 * pnorm(-abs(z))
    sig <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
    better <- if (dw > 0) nm[j] else nm[i]
    # if (p >= 0.05) better <- ""
    data.frame(model1 = nm[i], model2 = nm[j],
               dWAIC = round(dw, 1), se = round(se, 1),
               z = round(z, 2), p = round(p, 3),
               better = better, sig = sig,
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

