#' Compare INLAjoint models
#'
#' @description Pairwise comparison of joint models based on pointwise WAIC differences.
#'   Positive dWAIC means model2 fits better (lower WAIC).
#' @param ... Two or more objects of class \code{INLAjoint}.
#' @return A data.frame with one row per pair of models.
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
    dw <- sum(d)
    se <- sqrt(length(d)) * sd(d)
    z <- dw / se
    p <- 2 * pnorm(-abs(z))
    sig <- if (p < 0.01) "**" else if (p < 0.05) "*" else ""
    better <- if (dw > 0) nm[j] else nm[i]
    if (p >= 0.05) better <- ""
    data.frame(model1 = nm[i], model2 = nm[j],
               dWAIC = round(dw, 1), se = round(se, 1),
               z = round(z, 2), p = round(p, 3),
               better = better, sig = sig,
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}
