#' Compute time-dependent predictive accuracy metrics (AUC and Brier Score)
#'
#' @description Computes time-dependent AUC and Brier score from predicted
#' risks and observed outcomes.
#'
#' The function handles censoring through the IPCW (Inverse Probability of Censoring Weighted)
#' method, which reweights contributions by the inverse of the Kaplan-Meier
#' estimate of the censoring distribution.
#'
#' @param risk a numeric matrix of dimension \eqn{n \times T}{n x T} giving the predicted risk
#'   (i.e., event probability) for each subject (row) at each evaluation time (column).
#'   For survival predictions, this should be \eqn{1 - S(t)}. For competing risks, this should
#'   be the cause-specific cumulative incidence function (CIF).
#' @param pred_times a numeric vector of length \eqn{T} giving the evaluation time points
#'   corresponding to the columns of \code{risk}.
#' @param time a numeric vector of length \eqn{n} giving the observed event or censoring time
#'   for each subject.
#' @param event a numeric or integer vector of length \eqn{n} giving the event indicator for
#'   each subject (1 = event occurred, 0 = censored). For competing risks, use the cause-specific
#'   event indicator (1 = event of interest, 0 = censored or competing event; see Details).
#' @param metrics character vector specifying which metrics to compute. Options are \code{"auc"},
#'   \code{"brier"}, or both (default is \code{c("auc", "brier")}).
#'
#' @details
#'
#' \strong{Brier Score (IPCW):}
#'
#' The time-dependent Brier score at evaluation time \eqn{t} is defined as:
#' \deqn{BS(t) = \frac{1}{n}\sum_{i=1}^{n}\left[\frac{\hat{\pi}_i(t)^2 \cdot I(T_i > t)}{\hat{G}(t)} +
#' \frac{(1 - \hat{\pi}_i(t))^2 \cdot I(T_i \le t) \cdot \delta_i}{\hat{G}(T_i)}\right]}
#' where \eqn{\hat{\pi}_i(t)} is the predicted risk for subject \eqn{i} at time \eqn{t},
#' \eqn{T_i} is the observed time, \eqn{\delta_i} is the event indicator, and \eqn{\hat{G}(\cdot)}
#' is the Kaplan-Meier estimate of the censoring survival function.
#'
#' \strong{AUC (IPCW):}
#'
#' The time-dependent AUC at evaluation time \eqn{t} measures the concordance
#' between predicted risks and observed outcomes over all case-control pairs:
#' \deqn{AUC(t) = \frac{\sum_{i}\sum_{j} I(T_i \le t)\delta_i I(T_j > t)
#' I(\hat{\pi}_i(t) > \hat{\pi}_j(t)) w_{ij}(t)}{\sum_{i}\sum_{j} I(T_i \le t)
#' \delta_i I(T_j > t) w_{ij}(t)}}
#' where the IPCW weight is \eqn{w_{ij}(t) = 1 / (\hat{G}(T_i) \cdot \hat{G}(t))}.
#' Ties in predicted risk contribute 0.5.
#'
#' \strong{Competing Risks:}
#'
#' For competing risks, pass the cause-specific CIF as the \code{risk} matrix, and set
#' the \code{event} indicator to the binary indicator for the cause of interest
#' (1 = cause of interest occurred, 0 = censored or competing event). Each cause is
#' scored separately with its own call to \code{score_predictions()}.
#'
#' @return An object of class \code{"score_predictions"}, which is a list containing:
#' \describe{
#'   \item{scores}{A data.frame with columns \code{time} and the requested metrics
#'     (\code{AUC} and/or \code{Brier}).}
#'   \item{metrics}{Character vector of metrics computed.}
#'   \item{pred_times}{The evaluation time points.}
#'   \item{n_subjects}{Number of subjects.}
#'   \item{n_events}{Total number of events observed.}
#'   \item{n_censored}{Total number of censored observations.}
#'   \item{event_rate}{Proportion of subjects with events.}
#'   \item{n_cases}{Number of cases (events before each \code{pred_times}) used for AUC.}
#'   \item{n_controls}{Number of controls (event-free at each \code{pred_times}) used for AUC.}
#' }
#'
#' @seealso \code{\link{predict.INLAjoint}} for computing predictions.
#'
#' @examples
#' \donttest{
#' if(requireNamespace("INLA")){
#'
#' # simulate longitudinal + survival data
#' set.seed(1)
#' n <- 500
#' sex <- rbinom(n, 1, 0.5)
#' b <- MASS::mvrnorm(n, mu = c(0, 0),
#'                    Sigma = matrix(c(4, 1.2, 1.2, 0.6), 2, 2))
#' long_data <- do.call(rbind, lapply(1:n, function(i) {
#'   times <- seq(0, 10, by = 0.5)
#'   Y <- 5 - 0.3 * times + 0.8 * sex[i] + b[i,1] + b[i,2] * times +
#'        rnorm(length(times), 0, 0.5)
#'   data.frame(id = i, time = times, Y = Y, sex = sex[i])
#' }))
#' surv_data <- do.call(rbind, lapply(1:n, function(i) {
#'   tg <- seq(0, 15, length.out = 500)
#'   Yt <- 5 - 0.3*tg + 0.8*sex[i] + b[i,1] + b[i,2]*tg
#'   h <- exp(-4 + 0.4 * Yt)
#'   H <- cumsum(c(0, (h[-1]+h[-length(h)])/2 * diff(tg)))
#'   et <- tg[which(H >= rexp(1))[1]]
#'   st <- min(ifelse(is.na(et), 10, et), 10)
#'   data.frame(id = i, survtime = st, event = as.integer(st < 10), sex = sex[i])
#' }))
#' long_data <- merge(long_data, surv_data[, c("id","survtime")], by="id")
#' long_data <- long_data[long_data$time <= long_data$survtime, ]
#' long_data$survtime <- NULL
#'
#' # train/test split
#' train_ids <- 1:480
#' test_ids <- 481:500
#' train_long <- long_data[long_data$id %in% train_ids, ]
#' train_surv <- surv_data[surv_data$id %in% train_ids, ]
#' test_long  <- long_data[long_data$id %in% test_ids, ]
#' test_surv  <- surv_data[surv_data$id %in% test_ids, ]
#'
#' # fit joint model on training data
#' fit <- joint(
#'   formSurv = INLA::inla.surv(survtime, event) ~ sex,
#'   dataSurv = train_surv,
#'   formLong = list(Y ~ time + sex + (1 + time | id)),
#'   dataLong = train_long,
#'   id = "id", timeVar = "time", assoc = "CV"
#' )
#'
#' # predict for test subjects from a landmark
#' landmark <- 2
#' horizon <- 9
#' eval_times <- seq(4, 8)
#' valid_ids <- test_surv$id[test_surv$survtime >= landmark]
#' newdata <- test_long[test_long$id %in% valid_ids & test_long$time <= landmark, ]
#' obs_per_id <- table(newdata$id)
#' valid_ids <- as.integer(names(obs_per_id[obs_per_id >= 2]))
#' newdata <- newdata[newdata$id %in% valid_ids, ]
#' test_surv_valid <- test_surv[test_surv$id %in% valid_ids, ]
#'
#' # fine grid that includes the integer evaluation times
#' fine_grid <- sort(unique(c(seq(landmark, horizon, length.out = 50), eval_times)))
#' pred <- predict(fit, newData = newdata, horizon = horizon,
#'                 Csurv = landmark, survival = TRUE,
#'                 timePoints = fine_grid)
#'
#' # reshape predictions into a risk matrix
#' predS <- pred$PredS
#' ids <- unique(predS$id)
#' all_times <- sort(unique(predS$time))
#' surv_mat <- matrix(predS$Surv_quant0.5, nrow = length(ids),
#'                    ncol = length(all_times), byrow = TRUE)
#' risk_mat <- 1 - surv_mat
#'
#' # select integer evaluation times for scoring
#' keep <- match(eval_times, all_times)
#' risk_mat <- risk_mat[, keep]
#' pred_times <- eval_times
#'
#' # compute scores
#' sc <- score_predictions(
#'   risk = risk_mat,
#'   pred_times = pred_times,
#'   time = test_surv_valid$survtime,
#'   event = test_surv_valid$event
#' )
#' print(sc)
#' plot(sc)
#' }
#' }
#'
#' @export
score_predictions <- function(risk, pred_times, time, event,
                              metrics = c("auc", "brier")) {
  metrics <- match.arg(metrics, c("auc", "brier"), several.ok = TRUE)

  # validate risk input
  if (is.numeric(risk) && is.null(dim(risk))) {
    risk <- matrix(risk, ncol = 1)
  }
  if (!is.matrix(risk) || !is.numeric(risk)) {
    stop("'risk' must be a numeric matrix (n subjects x T time points).")
  }

  n <- nrow(risk)
  n_times <- ncol(risk)

  if (length(pred_times) != n_times) {
    stop("Length of 'pred_times' (", length(pred_times),
         ") must match number of columns in 'risk' (", n_times, ").")
  }
  if (length(time) != n) {
    stop("Length of 'time' (", length(time), ") must match number of rows in 'risk' (",
         n, ").")
  }
  if (length(event) != n) {
    stop("Length of 'event' (", length(event), ") must match number of rows in 'risk' (",
         n, ").")
  }
  if (any(is.na(time)) || any(is.na(event))) {
    stop("'time' and 'event' must not contain NA values.")
  }
  if (!all(event %in% c(0, 1))) {
    stop("'event' must be a binary indicator (0 = censored, 1 = event). ",
         "For competing risks, provide the cause-specific binary indicator.")
  }

  # Kaplan-Meier estimate of the censoring survival distribution
  # (flip event indicator: censoring becomes the "event")
  cens_event <- 1L - as.integer(event)
  ord <- order(time)
  t_ord <- time[ord]
  d_ord <- cens_event[ord]
  unique_times <- sort(unique(t_ord))
  n_at_risk <- length(t_ord)
  surv <- 1.0
  G_time <- 0
  G_surv <- 1
  for (tk in unique_times) {
    at_tk <- (t_ord == tk)
    n_cens <- sum(at_tk & d_ord == 1)
    n_total <- sum(at_tk)
    if (n_at_risk > 0 && n_cens > 0) {
      surv <- surv * (1 - n_cens / n_at_risk)
    }
    surv <- max(surv, 1e-8)
    G_time <- c(G_time, tk)
    G_surv <- c(G_surv, surv)
    n_at_risk <- n_at_risk - n_total
  }
  G <- data.frame(time = G_time, surv = G_surv)

  # compute metrics at each evaluation time
  auc_vals <- brier_vals <- rep(NA, n_times)
  n_cases <- n_controls <- integer(n_times)

  for (j in seq_len(n_times)) {
    t_eval <- pred_times[j]
    risk_j <- risk[, j]

    cases <- which(time <= t_eval & event == 1)
    controls <- which(time > t_eval)
    n_cases[j] <- length(cases)
    n_controls[j] <- length(controls)

    # evaluate censoring distribution G at needed time points
    G_t <- G$surv[pmax(findInterval(t_eval, G$time), 1L)]

    if ("brier" %in% metrics) {
      G_Ti <- G$surv[pmax(findInterval(time, G$time), 1L)]
      survived <- (time > t_eval)
      term1 <- (risk_j^2) * survived / G_t
      event_before <- (time <= t_eval) & (event == 1)
      term2 <- ((1 - risk_j)^2) * event_before / G_Ti
      brier_vals[j] <- sum(term1 + term2, na.rm = TRUE) / n
    }

    if ("auc" %in% metrics) {
      if (length(cases) == 0 || length(controls) == 0) {
        auc_vals[j] <- NA_real_
      } else {
        G_Ti_cases <- G$surv[pmax(findInterval(time[cases], G$time), 1L)]
        risk_cases <- risk_j[cases]
        risk_controls <- risk_j[controls]
        w_cases <- 1 / (G_Ti_cases * G_t)
        conc_matrix <- outer(risk_cases, risk_controls, function(a, b) {
          ifelse(a > b, 1, ifelse(a == b, 0.5, 0))
        })
        w_matrix <- matrix(w_cases, nrow = length(cases), ncol = length(controls))
        numerator <- sum(conc_matrix * w_matrix)
        denominator <- sum(w_matrix)
        auc_vals[j] <- if (denominator == 0) NA_real_ else numerator / denominator
      }
    }
  }

  scores <- data.frame(time = pred_times)
  if ("auc" %in% metrics) scores$AUC <- auc_vals
  if ("brier" %in% metrics) scores$Brier <- brier_vals

  result <- list(
    scores = scores,
    metrics = metrics,
    pred_times = pred_times,
    n_subjects = n,
    n_events = sum(event == 1),
    n_censored = sum(event == 0),
    event_rate = mean(event == 1),
    n_cases = n_cases,
    n_controls = n_controls
  )
  class(result) <- "score_predictions"
  return(result)
}


#' Print method for score_predictions objects
#'
#' @param x an object of class \code{"score_predictions"}.
#' @param digits number of digits to display. Default is 4.
#' @param ... additional arguments (ignored).
#'
#' @export
print.score_predictions <- function(x, digits = 4, ...) {
  cat("Predictive Accuracy Scores (IPCW-corrected)\n")
  cat("=============================================\n")
  cat("Subjects:", x$n_subjects, " | Events:", x$n_events,
      " | Censored:", x$n_censored,
      " | Event rate:", round(x$event_rate * 100, 1), "%\n\n")

  df <- x$scores
  if ("AUC" %in% colnames(df)) df$AUC <- round(df$AUC, digits)
  if ("Brier" %in% colnames(df)) df$Brier <- round(df$Brier, digits)
  df$Cases <- x$n_cases
  df$Controls <- x$n_controls

  print(df, row.names = FALSE)

  cat("\nSummary:\n")
  if ("AUC" %in% colnames(x$scores)) {
    auc_valid <- x$scores$AUC[!is.na(x$scores$AUC)]
    if (length(auc_valid) > 0) {
      cat("  Mean AUC:", round(mean(auc_valid), digits),
          " | Range: [", round(min(auc_valid), digits), ",",
          round(max(auc_valid), digits), "]\n")
    }
  }
  if ("Brier" %in% colnames(x$scores)) {
    brier_valid <- x$scores$Brier[!is.na(x$scores$Brier)]
    if (length(brier_valid) > 0) {
      cat("  Mean Brier:", round(mean(brier_valid), digits),
          " | Range: [", round(min(brier_valid), digits), ",",
          round(max(brier_valid), digits), "]\n")
    }
  }
  invisible(x)
}


#' Plot method for score_predictions objects
#'
#' @description Produces time-dependent AUC and/or Brier score plots.
#'
#' @param x an object of class \code{"score_predictions"}.
#' @param which character vector specifying which plots to produce. Options are
#'   \code{"auc"}, \code{"brier"}, or both. Default uses the metrics available in \code{x}.
#' @param add logical; if \code{TRUE}, add lines to an existing plot instead of creating
#'   a new one. Useful for overlaying results from different models. Default is \code{FALSE}.
#' @param col color for the lines and points. Default is \code{"blue"}.
#' @param lwd line width. Default is 2.
#' @param pch point character. Default is 19 (filled circle).
#' @param main optional title(s) for the plot(s). If a single string is provided, it is
#'   used for all panels. If a named list with elements \code{"auc"} and/or \code{"brier"},
#'   each is used for the respective panel.
#' @param legend_label character string for the legend label when using \code{add = TRUE}.
#' @param ... additional graphical parameters passed to \code{plot()} or \code{lines()}.
#'
#' @export
plot.score_predictions <- function(x, which = NULL, add = FALSE,
                                   col = "blue", lwd = 2, pch = 19,
                                   main = NULL, legend_label = NULL, ...) {
  if (is.null(which)) which <- x$metrics
  which <- match.arg(which, c("auc", "brier"), several.ok = TRUE)

  which <- intersect(which, x$metrics)
  if (length(which) == 0) {
    stop("No matching metrics to plot.")
  }

  both <- length(which) == 2
  if (both && !add) {
    oldpar <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
    on.exit(par(oldpar))
  }

  times <- x$pred_times

  if ("auc" %in% which) {
    auc <- x$scores$AUC
    auc_main <- "Time-Dependent AUC (IPCW)"
    if (!is.null(main)) {
      if (is.list(main) && "auc" %in% names(main)) auc_main <- main$auc
      else if (is.character(main) && length(main) == 1) auc_main <- main
    }

    if (add) {
      lines(times, auc, type = "o", col = col, lwd = lwd, pch = pch, ...)
    } else {
      plot(times, auc, type = "o", col = col, lwd = lwd, pch = pch,
           ylim = c(0, 1), xlab = "Time", ylab = "AUC",
           main = auc_main, ...)
      abline(h = 0.5, lty = 2, col = "gray")
    }
  }

  if ("brier" %in% which) {
    brier <- x$scores$Brier
    brier_main <- "Time-Dependent Brier Score (IPCW)"
    if (!is.null(main)) {
      if (is.list(main) && "brier" %in% names(main)) brier_main <- main$brier
      else if (is.character(main) && length(main) == 1) brier_main <- main
    }

    brier_valid <- brier[!is.na(brier)]
    ylim_max <- if (length(brier_valid) > 0) max(brier_valid) * 1.1 else 0.25

    if (add) {
      lines(times, brier, type = "o", col = col, lwd = lwd, pch = pch, ...)
    } else {
      plot(times, brier, type = "o", col = col, lwd = lwd, pch = pch,
           ylim = c(0, ylim_max), xlab = "Time", ylab = "Brier Score",
           main = brier_main, ...)
    }
  }

  invisible(x)
}
