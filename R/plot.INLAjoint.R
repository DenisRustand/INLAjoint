#' Plot the output from a multivariate joint model for longitudinal and/or survival data
#'
#' @description This function provide plots for the output of a
#' multivariate joint model for longitudinal and/or survival data.
#' The output can be stored into an object and manipulated as
#' a list of ggplot outputs, described bellow.
#'
#' @param x an object with the output of the the \link{joint} function
#' @param ... Extra arguments including:
#' \code{sdcor}: logical indicating if the random effects
#' correlation are to be shown. If FALSE the covariance is shown.
#' \code{priors}: logical indicating if the priors are added to the posterior marginals plots.
#' @return return a named list of \code{ggplot} objects containing:
#' \describe{
#'  \item{\code{Outcomes}}{
#'  as a list of length equal the number of longitudinal
#'  outcomes plus the number of survival outcomes, each one including
#'  the plot for the posterior marginals of the associated fixed effects
#'  and residual or baseline variance (or standard error).
#'  Each element contains the plot for the posterior marginal.}
#'  \item{\code{Covariances}}{
#'  the plots for the posterior marginal distribution of the
#'  covariance parameters.}
#'  \item{\code{Associations}}{
#'  the plots for the posterior marginal distribution of the
#'  association parameters.}
#'  \item{\code{Random}}{
#'  The plot for the fitted baseline risk functions as shown
#'  as the posterior mean and credible interval.}
#'  }
#'
#' @import ggplot2
#' @importFrom grDevices dev.new
#' @export

plot.INLAjoint <- function(x, ...) {
  arguments <- list(...)
  if(is.null(arguments$sdcor)) sdcor=F else sdcor=arguments$sdcor
  if(is.null(arguments$priors)) priors=F else priors=arguments$priors
  stopifnot(is.logical(sdcor))
  stopifnot(is.logical(priors))
  y <- NULL
  group <- NULL
  type <- NULL
  kid <- NULL
  ID <- NULL
  lower <- NULL
  upper <- NULL
  out <- list(
        Outcomes=NULL, Covariances=NULL,
        Associations=NULL, Baseline=NULL, Random=NULL)
  if(priors){
    x.fixed.intercept.prior <- seq(x$priors_used$priorFixed$mean.intercept-2*sqrt(1/x$priors_used$priorFixed$prec.intercept),
                                   x$priors_used$priorFixed$mean.intercept+2*sqrt(1/x$priors_used$priorFixed$prec.intercept), len=500)
    fixed.intercept.prior <- dnorm(x.fixed.intercept.prior, x$priors_used$priorFixed$mean.intercept, sqrt(1/x$priors_used$priorFixed$prec.intercept))
    x.fixed.prior <- seq(x$priors_used$priorFixed$mean-2*sqrt(1/x$priors_used$priorFixed$prec),
                             x$priors_used$priorFixed$mean+2*sqrt(1/x$priors_used$priorFixed$prec), len=500)
    fixed.prior <- dnorm(x.fixed.prior, x$priors_used$priorFixed$mean, sqrt(1/x$priors_used$priorFixed$prec))
    x.assoc.prior <- seq(x$priors_used$priorAssoc$mean-2*sqrt(1/x$priors_used$priorAssoc$prec),
                         x$priors_used$priorAssoc$mean+2*sqrt(1/x$priors_used$priorAssoc$prec), len=500)
    assoc.prior <- dnorm(x.assoc.prior, x$priors_used$priorAssoc$mean, sqrt(1/x$priors_used$priorAssoc$prec))
    x.priorSRE_ind.prior <- seq(x$priors_used$priorSRE_ind$mean-2*sqrt(1/x$priors_used$priorSRE_ind$prec),
                                x$priors_used$priorSRE_ind$mean+2*sqrt(1/x$priors_used$priorSRE_ind$prec), len=500)
    priorSRE_ind.prior <- dnorm(x.priorSRE_ind.prior, x$priors_used$priorSRE_ind$mean, sqrt(1/x$priors_used$priorSRE_ind$prec))

    # residual error (gaussian and lognormal)
    if(sdcor){
      limitsreserr <- qgamma(c(.99, .01), shape = exp(1), rate = 1/exp(5e-05))^(-1)
      x.reserr.prior <- seq(limitsreserr[1], limitsreserr[2], length = 501)
      reserr.prior <- sqrt(exp(dgamma(1/x.reserr.prior, shape = exp(1), rate = 1/exp(5e-05), log=TRUE) - 2 * log(x.reserr.prior)))
    }else{
      limitsreserr <- qgamma(c(.99, .01), shape = exp(1), rate = 1/exp(5e-05))^(-1)
      x.reserr.prior <- seq(limitsreserr[1], limitsreserr[2], length = 501)
      reserr.prior <-  sqrt(exp(dgamma(1/x.reserr.prior, shape = exp(1), rate = 1/exp(5e-05), log=TRUE) - 2 * log(x.reserr.prior)))
    }

    # Random effects variance prior ~invGamma(alpha = r/2, beta = R/2):
    N <- 1e4
    invsample <- rWishart(N, df=x$priors_used$priorRandom$r, Sigma=solve(diag(x$priors_used$priorRandom$R, 2)))
    sample <- lapply(seq.int(N), function(x) solve(invsample[,,x]))
    # remove dependency on MCMCpack
    # invsample <- lapply(seq.int(N), function(x) MCMCpack::rwish(x$priors_used$priorRandom$r, solve(diag(x$priors_used$priorRandom$R, 2))))
    # sample <- lapply(invsample, solve)
    cov.prior <- sapply(sample, function(x) x[1,2]) # covariance
    if(sdcor){
      limits <- qgamma(c(.999, .001), shape = x$priors_used$priorRandom$r/2, rate = 1/(x$priors_used$priorRandom$R/2))^(-1)
      x.sd.prior <- seq(limits[1], limits[2], length = 501)
      sd.prior <- sqrt(exp(dgamma(1/x.sd.prior, shape = x$priors_used$priorRandom$r/2, rate = 1/(x$priors_used$priorRandom$R/2), log=TRUE) - 2 * log(x.sd.prior)))
      # correlation
      denom <- apply(sqrt(sapply(sample, diag)), 2, prod)
      stopifnot(abs(cov.prior) < denom)
      corr.prior <- cov.prior/denom
    }else{
      limits <- qgamma(c(.99, .01), shape = x$priors_used$priorRandom$r/2, rate = 1/(x$priors_used$priorRandom$R/2))^(-1)
      x.var.prior <- seq(limits[1], limits[2], length = 501)
      var.prior <- exp(dgamma(1/x.var.prior, shape = x$priors_used$priorRandom$r/2, rate = 1/(x$priors_used$priorRandom$R/2), log=TRUE) - 2 * log(x.var.prior))
    }
  }
  trimMarginal <- function(m, p=0.001) {
          m <- INLA::inla.smarginal(m)
          ab <- INLA::inla.qmarginal(c(p, 1-p), m)
          ii <- which((m$x>=ab[1])&(m$x<=ab[2]))
          return(list(x=m$x[ii], y=m$y[ii]))
  }
  joinMarginals <- function(listofmarginals, trim=TRUE, p=0.001) {
      ## join the marginals into one long data with indicator
      if(trim) {
          Reduce('rbind', lapply(1:length(listofmarginals), function(k) {
              data.frame(m=k, trimMarginal(listofmarginals[[k]], p=p))
          }))
      } else {
          Reduce('rbind', lapply(1:length(listofmarginals), function(k) {
              data.frame(m=k, as.data.frame(listofmarginals[[k]]))
          }))
      }
  }

  out.patt <- '_[LS][0-9]+$'
  hhid <- sapply(x$internal.marginals.hyperpar, attr, 'hyperid')
  if(length(hhid)>0) hid <- sapply(strsplit(hhid, '|', fixed=TRUE), tail, 1) else hid <- NULL
  if(length(grep("Intercept_S", names(x$marginals.fixed)))>0){
    JRM <- x$marginals.fixed[-grep("Intercept_S", names(x$marginals.fixed))]
  }else{
    JRM <- x$marginals.fixed
  }
  x.n <- length(x.names <- names(JRM))
  if(x.n>0) {
      x.psub <- regexpr(out.patt, x.names)
      x.group <- substring(x.names, x.psub+1)
      xMargs <- joinMarginals(JRM)
      xMargs$Effect <- factor(x.names[xMargs$m], x.names, x.names)
      xMargs$Outcome <- x.group[xMargs$m]
      lfamilies0 <- c(
          'gaussian', 'lognormal', 'gamma', 'beta'
      ) ## others to be added
      hl.jj <- which(unlist(x$famLongi) %in% lfamilies0)
      nhl <- length(hl.jj)
      hs.jj <- grep('baseline[1-9]', hid)
      nhs <- length(hs.jj)
      if((nhl+nhs)>0) {
          thMargs <- joinMarginals(lapply(
              x$internal.marginals.hyperpar[c(seq_len(nhl), hs.jj)],
              function(m) INLA::inla.tmarginal(function(x) exp(-x/(1+sdcor)), m)),
              trim=FALSE)
          thMargs$Effect <-
              paste0(c(rep(paste0(
                  'Residual ',
                   c('Var.', 'S.D.')[sdcor+1], '_L'), nhl),
                rep(paste0(
                    'Baseline',
                    c('Var.', 'S.D.')[sdcor+1], '_S'), nhs)),
                c(seq_len(nhl), seq_len(nhs)))[thMargs$m]
          thMargs$Outcome <-
              c(paste0(rep('L', nhl), seq_len(nhl)),
                paste0(rep('S', nhs), seq_len(nhs))
                )[thMargs$m]
          xMargs <- as.data.frame(rbind(
              xMargs, thMargs))
      }
      if(priors){
        xMargs$group <- "posterior"
        for(outc in unique(xMargs$Outcome)){
          DatOutc <- xMargs[xMargs$Outcome==outc,]
          for(effc in unique(DatOutc$Effect)){
            if(substr(effc, 1, 10)=="Intercept_"){
              addPrior <- data.frame("m"=1, "x"=x.fixed.intercept.prior, "y"=fixed.intercept.prior, "Effect"=effc, "Outcome"=outc, "group"="prior")
              xMargs <- rbind(xMargs, addPrior)
            }else if(substr(effc, 1, 8)=="Residual"){
              if(substr(effc, 10, 11)=="Va"){ # variance
                addPrior <- data.frame("m"=1, "x"=x.reserr.prior, "y"=reserr.prior, "Effect"=effc, "Outcome"=outc, "group"="prior")
                xMargs <- rbind(xMargs, addPrior)
              }else{
                addPrior <- data.frame("m"=1, "x"=x.reserr.prior, "y"=reserr.prior, "Effect"=effc, "Outcome"=outc, "group"="prior")
                xMargs <- rbind(xMargs, addPrior)
              }
            }else if(substr(effc, 1, 8)=="Baseline"){
              # pc prior here?
            }else{
              addPrior <- data.frame("m"=1, "x"=x.fixed.prior, "y"=fixed.prior, "Effect"=effc, "Outcome"=outc, "group"="prior")
              xMargs <- rbind(xMargs, addPrior)
            }
          }
        }
        out$Outcomes <- lapply(
          split(xMargs, xMargs$Outcome), function(d)
            ggplot(d, aes(x=x, y=y, group=group, colour=group)) +
            xlab('') +
            ylab('Density') +
            geom_line() +
            facet_wrap(~Effect, scales='free'))
      }else{
        out$Outcomes <- lapply(
          split(xMargs, xMargs$Outcome), function(d)
            ggplot(d, aes(x=x, y=y)) +
            xlab('') +
            ylab('Density') +
            geom_line() +
            facet_wrap(~Effect, scales='free'))
      }
  }
  nhk <- length(hd.idx <- grep('^Theta[0-9]+ for ', names(hid)))
  if(nhk>0) {
      k1 <- grep('Theta1 for ', names(hid))
      class(x) <- 'inla'
      out$Covariances <- vector('list', length(k1))
      names(out$Covariances) <- paste0('L', 1:length(k1))
      for (l in 1:length(k1)) {
        kdsamples <- INLA::inla.iidkd.sample(
            2e4, x, hid[k1[l]], return.cov=!sdcor)
        k <- nrow(kdsamples[[1]])
        if(k>0) {
            kdsamples <- sapply(kdsamples, as.vector)
            ii.m <- matrix(1:(k*k), k)
            kdens <- Reduce('rbind', lapply(1:nrow(kdsamples), function(j) {
                ldu <- 2-(j%in%diag(ii.m))
                if(ldu==1) {
                    dd <- density(log(kdsamples[j,]))
                    dd$x <- exp(dd$x)
                    dd$y <- dd$y/exp(dd$x)
                } else {
                    dd <- density(kdsamples[j,])
                }
                ldu <- ldu + (j%in%(ii.m[lower.tri(ii.m)]))
                return(data.frame(
                    m=j, as.data.frame(
                        trimMarginal(dd[c('x', 'y')], 0.005)),
                    ldu=ldu))
            }))
            kdnames <- apply(expand.grid(
                x$REstruc, x$REstruc), 1,
                function(x) paste(unique(x), collapse=':'))
            kdens$Effect <- factor(
                kdnames[kdens$m],
                unique(kdnames))
            kdens$type <- factor(
                kdens$ldu, 1:3,
                c(c('Var.', 'St.Dev.')[sdcor+1],
                  rbind(c('Cov.', 'Cov.'), c('Correl.', 'Correl.'))[sdcor+1,]))
            if(priors){
              kdens$group <- "posterior"
              for(effc in unique(kdens$Effect)){
                if(kdens[which(kdens$Effect==effc)[1], "type"]=="Var."){
                  addPrior <- data.frame("m"=1, "x"=x.var.prior, "y"=var.prior, "ldu"=1, "Effect"=effc, "type"=kdens[which(kdens$Effect==effc)[1], "type"], "group"="prior")
                  kdens <- rbind(kdens, addPrior)
                }else if(kdens[which(kdens$Effect==effc)[1], "type"]=="St.Dev."){
                  addPrior <- data.frame("m"=1, "x"=x.sd.prior, "y"=sd.prior, "ldu"=1, "Effect"=effc, "type"=kdens[which(kdens$Effect==effc)[1], "type"], "group"="prior")
                  kdens <- rbind(kdens, addPrior)
                }else if(kdens[which(kdens$Effect==effc)[1], "type"]=="Cov."){
                  addPrior <- data.frame("m"=1, "x"=density(cov.prior)$x, "y"=density(cov.prior)$y, "ldu"=1, "Effect"=effc, "type"=kdens[which(kdens$Effect==effc)[1], "type"], "group"="prior")
                  kdens <- rbind(kdens, addPrior)
                }else if(kdens[which(kdens$Effect==effc)[1], "type"]=="Correl."){
                  addPrior <- data.frame("m"=1, "x"=density(corr.prior)$x, "y"=density(corr.prior)$y, "ldu"=1, "Effect"=effc, "type"=kdens[which(kdens$Effect==effc)[1], "type"], "group"="prior")
                  kdens <- rbind(kdens, addPrior)
                }
              }
              out$Covariances[[l]] <- ggplot(kdens, aes(x=x,y=y,group=group)) +
                xlab('') +
                ylab('Density') +
                geom_line(aes(color=type, linetype=group)) +
                facet_wrap(~Effect, scales='free')
            }else{
              out$Covariances[[l]] <- ggplot(kdens, aes(x=x,y=y)) +
                xlab('') +
                ylab('Density') +
                geom_line(aes(color=type, linetype=type)) +
                facet_wrap(~Effect, scales='free')
            }
        } else {
            warning('Something wrong with', kid[k1], 'happened!')
        }
      }
  }
  nhc <- length(hc.idx <- grep('Beta_intern for ', names(hid)))
  if(nhc>0) {
      cMargs <- joinMarginals(
          x$internal.marginals.hyperpar[hc.idx])
      cnames <- substring(names(x$internal.marginals.hyperpar)[hc.idx],16)
      cMargs$Effect <- factor(cnames[cMargs$m], cnames, cnames)
      out$Associations <- ggplot(cMargs, aes(x=x,y=y)) +
          xlab('') +
          ylab('Density') +
          geom_line() +
          facet_wrap(~Effect, scales='free')
  }
  rnames <- names(x$summary.random)
  nbas <- length(bas.idx <- grep(
      '^baseline[0-9]+', rnames))
  nbasP <- length(c(grep("weibullsurv", unlist(x$basRisk)), # number of parametric baseline risks
                    grep("exponentialsurv", unlist(x$basRisk))))
  if(nbas>0) {
    BaselineValues <- NULL
    for(i in 1:nbas){
      BHmean <- NULL
      BHlo <- NULL
      BHup <- NULL
      for(j in 1:length(x$marginals.random[[i]])){
        # m <- inla.smarginal(x$marginals.random[[i]][[j]])
        # ab <- inla.qmarginal(c(0.001, 0.999), m)
        # ii <- which((m$x>=ab[1]) & (m$x<=ab[2]))
        # m$x <- m$x[ii]
        # m$y <- m$y[ii]
        Mm <- INLA::inla.qmarginal(c(0.025, 0.5, 0.975), x$marginals.random[[i]][[j]])
        BHmean <- c(BHmean, exp(Mm[2]))
        BHlo <- c(BHlo, exp(Mm[1]))
        BHup <- c(BHup, exp(Mm[3]))
      }
      BaselineValues <- rbind(BaselineValues,
                              cbind(time=x$summary.random[[paste0("baseline",i,".hazard")]]$ID,
                              mean=BHmean,
                              lower=BHlo,
                              upper=BHup,
                              S=i))
    }
      jsr <- Reduce('rbind', lapply(1:nbas, function(k) {
          data.frame(x$summary.random[[bas.idx[k]]],
                     S=paste0('S',k))
      }))
      colnames(jsr)  <- gsub('X0.', 'q', colnames(jsr), fixed=TRUE)
      out$Baseline <- ggplot(jsr, aes(x=ID)) +
          geom_ribbon(aes(ymin=BaselineValues[,"lower"],
                          ymax=BaselineValues[,"upper"]),
                      fill='grey70') +
          geom_line(aes(y=BaselineValues[,"mean"])) +
          xlab('Time') +
          ylab('Baseline risk') +
          facet_wrap(~S,  scales='free')
  }
  if(nbasP>0){
    HW0 <- function(t, lambda, alpha){ # risk function Weibull variant 0 (also exponential for alpha=1)
      res = lambda*alpha*t^(alpha-1)
    }
    HW1 <- function(t, lambda, alpha){ # risk function Weibull variant 1
      res = lambda*alpha*(lambda*t)^(alpha-1)
    }
    BHM <- NULL # baseline risk marginals
    nbl <- 1 # to keep track of baseline risk in case of multiple parametric survival outcomes
    nbl2 <- 1
    BaselineValues <- NULL
    if(x$.args$control.compute$config){
      NSAMPLES = 500
      SEL <- sapply(paste0(rownames(x$summary.fixed)[grep("Intercept_S", rownames(x$summary.fixed))]), function(x) x=1, simplify=F)
      SAMPLES <- INLA::inla.posterior.sample(NSAMPLES, x, selection=SEL)
    }else{
      NSAMPLES = 500
      SEL <- sapply(paste0(rownames(x$summary.fixed)[grep("Intercept_S", rownames(x$summary.fixed))]), function(x) x=1, simplify=F)
      SAMPLES <- INLA::inla.posterior.sample(NSAMPLES, x, selection=SEL)
      message("The parametric baseline risk is plotted without uncertainty, to get uncertainty, rerun the model with control=list(..., config=TRUE)")
    }
    for(i in 1:(nbas+nbasP)){
      if(nbasP==1){
        maxTime <- max(na.omit(x$.args$data$Yjoint[which(class(x$.args$data$Yjoint)=="inla.surv")][i]$time))
      }else if(nbasP>1){
        maxTime <- max(na.omit(x$.args$data$Yjoint[which(sapply(x$.args$data$Yjoint, class)=="inla.surv")][[i]]$time))
      }
      timePts <- seq(0, maxTime, len=500)
      timePts2 <- timePts#c(timePts[2], timePts[-1]) # avoid computing parametric baseline risk at time 0 since it's always 0
      Variant_i <- x$.args$control.family[[i]]$variant
      if(x$basRisk[[i]]=="exponentialsurv"){
        BHM <- append(BHM, list(INLA::inla.tmarginal(function(x) exp(x),
                                               x$marginals.fixed[grep("Intercept_S", names(x$marginals.fixed))][[i]])))
        names(BHM)[nbl] <- paste0("Exponential (rate)_S", i)
        if(x$.args$control.compute$config){ # config set to TRUE then we can compute uncertainty
          curves_SMP <- sapply(1:NSAMPLES, function(x) HW0(t=timePts2, lambda=exp(SAMPLES[[x]]$latent[nbl]), alpha=1))
          QUANT <- apply(curves_SMP, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))
          values_i <- t(QUANT)
        }else{
          values_i <- HW0(t=timePts2, lambda=exp(x$summary.fixed[grep("Intercept_S", names(x$marginals.fixed))][[i]]), alpha=1)
        }
        name_i <- paste0("Exponential baseline risk (S", nbl, ")")
      }else if(x$basRisk[[i]]=="weibullsurv"){
        BHM <- append(BHM, list(INLA::inla.tmarginal(function(x) exp(x),
                                               x$marginals.fixed[grep("Intercept_S", names(x$marginals.fixed))][[i]])))
        names(BHM)[nbl2] <- paste0("Weibull (scale)_S", i)
        BHM <- append(BHM, list(x$marginals.hyperpar[grep("weibull", names(x$marginals.hyperpar))][[nbl]]))
        names(BHM)[nbl2+1] <- paste0("Weibull (shape)_S", i)
        if(x$.args$control.compute$config){ # config set to TRUE then we can compute uncertainty
          if(Variant_i==0){
            curves_SMP <- sapply(1:NSAMPLES, function(x) HW0(t=timePts2, lambda=exp(SAMPLES[[x]]$latent[nbl]), alpha=SAMPLES[[x]]$hyperpar[nbl]))
          }else if(Variant_i==1){
            curves_SMP <- sapply(1:NSAMPLES, function(x) HW1(t=timePts2, lambda=exp(SAMPLES[[x]]$latent[nbl]), alpha=SAMPLES[[x]]$hyperpar[nbl]))
          }
          QUANT <- apply(curves_SMP, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))
          values_i <- t(QUANT)
        }else{
          if(Variant_i==0){
            values_i <- HW0(t=timePts2, lambda=exp(x$summary.fixed[grep("Intercept_S", rownames(x$summary.fixed)), "mean"][[i]]),
                            alpha=x$summary.hyperpar[grep("weibull", names(x$marginals.hyperpar)), "mean"][[nbl]])
          }else if(Variant_i==1){
            values_i <- HW1(t=timePts2, lambda=exp(x$summary.fixed[grep("Intercept_S", rownames(x$summary.fixed)), "mean"][[i]]),
                            alpha=x$summary.hyperpar[grep("weibull", names(x$marginals.hyperpar)), "mean"][[nbl]])
          }
        }
        name_i <- paste0("Weibull baseline risk (S", nbl, ")")
      }
      nbl2 <- nbl2+2
      nbl <- nbl+1
      BaselineValues <- rbind(BaselineValues, data.frame(timePts, values_i, name_i))
    }
    if(x$.args$control.compute$config){ # with uncertainty
      colnames(BaselineValues) <- c("x", "lower", "y", "upper", "Effect")
      out$Baseline <- ggplot(data=BaselineValues) +
        geom_ribbon(aes(x=x, ymin=lower,
                        ymax=upper),
                    fill='grey70') +
        geom_line(aes(x=x, y=y)) +
        xlab('Time') +
        ylab('Baseline risk') +
        facet_wrap(~Effect,  scales='free')
    }else{ # only means
      colnames(BaselineValues) <- c("x", "y", "Effect")
      out$Baseline <- ggplot(BaselineValues, aes(x=x, y=y)) +
        geom_line(aes(y=BaselineValues[,"y"])) +
        xlab('Time') +
        ylab('Baseline risk') +
        facet_wrap(~Effect,  scales='free')
    }
    sMargs <- joinMarginals(BHM)
    snames <- names(BHM)
    sMargs$Effect <- factor(snames[sMargs$m], snames, snames)
    out$BaselineParam <- ggplot(sMargs, aes(x=x,y=y)) +
      xlab('') +
      ylab('Density') +
      geom_line() +
      facet_wrap(~Effect, scales='free')
  }
  out <- out[!sapply(out, is.null)]
  class(out) <- c("plot.INLAjoint", "list")
  return(out)
}
print.plot.INLAjoint <- function(x, ...) {
  for(i in 1:length(x)) {
    if(i>1) dev.new()
    if(all(class(x[[i]]) == "list")) {
      for(j in 1:length(x[[i]])) {
        if(j>1) dev.new()
        print(x[[i]][[j]])
      }
    } else {
      print(x[[i]])
    }
  }
}
