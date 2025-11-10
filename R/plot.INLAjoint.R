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
#' @importFrom grDevices dev.new dev.interactive devAskNewPage dev.flush
#' @export

plot.INLAjoint <- function(x, ...) {
  arguments <- list(...)
  if("run" %in% names(x)) if(!x$run) stop("Please run the model (with function `joint.run()`)")
  if(is.null(arguments$sdcor)) sdcor=F else sdcor=arguments$sdcor
  if(is.null(arguments$priors)) priors=F else priors=arguments$priors
  if(is.null(arguments$NL_fun)) NL_fun=F else priors=arguments$NL_fun
  stopifnot(is.logical(sdcor))
  stopifnot(is.logical(priors))
  methodNL="sampling"
  NsampleNL=1000
  if(is.null(arguments$NLeffectonly)) NLeffectonly=FALSE else NLeffectonly=arguments$NLeffectonly
  if(is.null(arguments$hr)) hr=FALSE else hr=arguments$hr
  # NLeffectonly <- F
  y <- NULL
  group <- NULL
  type <- NULL
  kid <- NULL
  ID <- NULL
  lower <- NULL
  upper <- NULL
  group_factor <- NULL
  fitted_lower <- NULL
  fitted_upper <- NULL
  fitted_mean <- NULL
  out <- list(
        Outcomes=NULL, Covariances=NULL,
        Associations=NULL, Baseline=NULL, Random=NULL)
  if(priors){
    x.fixed.intercept.prior <- seq(x$priors_used$priorFixed$mean-2*sqrt(1/x$priors_used$priorFixed$prec),
                                   x$priors_used$priorFixed$mean+2*sqrt(1/x$priors_used$priorFixed$prec), len=500)
    fixed.intercept.prior <- dnorm(x.fixed.intercept.prior, x$priors_used$priorFixed$mean, sqrt(1/x$priors_used$priorFixed$prec))
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
                c(hl.jj, seq_len(nhs)))[thMargs$m]#seq_len(nhl)
          thMargs$Outcome <-
              c(paste0(rep('L', nhl), hl.jj),#seq_len(nhl)
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
            ggplot(d, aes(x=x, y=y, group=group, colour=group, linetype=group)) +
            xlab('') +
            ylab('Density') +
            geom_line() +
            facet_wrap(~Effect, scales='free') +
            theme_minimal())
      }else{
        out$Outcomes <- lapply(
          split(xMargs, xMargs$Outcome), function(d)
            ggplot(d, aes(x=x, y=y)) +
            xlab('') +
            ylab('Density') +
            geom_line() +
            facet_wrap(~Effect, scales='free')+
            theme_minimal())
      }
  }
  nhk <- length(hd.idx <- unique(c(grep('^Theta[0-9]+ for ', names(hid)), grep('Log precision for ID', names(hid)))))
  if(nhk>0) {
      k11 <- grep('Theta1 for ', names(hid))
      k12 <- grep('Log precision for ID', names(hid))
      k1 <- sort(c(k11, k12))
      class(x) <- 'inla'
      out$Covariances <- vector('list', length(k1))
      names(out$Covariances) <- paste0('L', 1:length(k1))
      for (l in 1:length(k1)) {
        if(k1[l] %in% k11){ # iddkd
          kdsamples <- INLA::inla.iidkd.sample(
            2e4, x, hid[k1[l]], return.cov=!sdcor)
          k <- nrow(kdsamples[[1]])
          if(k>0) {
            kdsamples <- sapply(kdsamples, as.vector)
            ii.m <- matrix(1:(k*k), k)
            kdens <- Reduce('rbind', lapply(1:nrow(kdsamples), function(j) {
              ldu <- 2-(j%in%diag(ii.m))
              # if(ldu==1) {
              #   dd <- density(log(kdsamples[j,]))
              #   dd$x <- exp(dd$x)
              #   dd$y <- dd$y/exp(dd$x)
              # } else {
                dd <- density(kdsamples[j,])
              # }
              ldu <- ldu + (j%in%(ii.m[lower.tri(ii.m)]))
              return(data.frame(
                m=j, as.data.frame(
                  trimMarginal(dd[c('x', 'y')], 0.005)),
                ldu=ldu))
            }))
            if(l>9) shiftRE <- 3 else shiftRE <- 2
            struc_k_start <- grep(substr(hid[k1[l]], 3, nchar(hid[k1[l]])), x$REstruc)
            #which(substring(x$REstruc, nchar(x$REstruc)-shiftRE, nchar(x$REstruc))==paste0("_L", l))
            struc_k <- x$REstruc[struc_k_start:(struc_k_start+(k-1))]
            kdnames <- apply(expand.grid(
              struc_k, struc_k), 1,
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
                facet_wrap(~Effect, scales='free') +
                theme_minimal()
            }else{
              out$Covariances[[l]] <- ggplot(kdens, aes(x=x,y=y)) +
                xlab('') +
                ylab('Density') +
                geom_line(aes(color=type, linetype=type)) +
                facet_wrap(~Effect, scales='free') +
                theme_minimal()
            }
          } else {
            warning('Something wrong with', kid[k1], 'happened!')
          }
        }else if(k1[l] %in% k12){ # iid
          if(k1[l] == k12[1]) kdsamples <- 1/exp(INLA::inla.hyperpar.sample(2e4, x, intern=TRUE)[, k12])
          if(sdcor) kdsamples <- sqrt(kdsamples)
          k <- nrow(kdsamples)
          if(l>9) shiftRE <- 3 else shiftRE <- 2
          struc_k <- x$REstruc[which(substring(x$REstruc, nchar(x$REstruc)-shiftRE, nchar(x$REstruc))==paste0("_L", l) |
                                       substring(x$REstruc, nchar(x$REstruc)-shiftRE, nchar(x$REstruc))==paste0("_S", l))]
          if(!is.null(dim(kdsamples))){
            k_dens <- density(kdsamples[, which(k12==k1[l])])
          }else{
            k_dens <- density(kdsamples)
          }
          if(sdcor) typeRE <- rep("St.Dev.", length(k_dens$x)) else typeRE <- rep("Var.", length(k_dens$x))
          if(all.equal(struc_k, character(0))==TRUE) struc_k <- "Random effect"
          if(all.equal(typeRE, character(0))==TRUE) typeRE <- ""
          kdens <- data.frame("x"=k_dens$x, "y"=k_dens$y, "Effect"=rep(struc_k, length(k_dens$x)), "type"=typeRE)
          if(priors){
            kdens$group <- "posterior"
            for(effc in unique(kdens$Effect)){
              if(kdens[which(kdens$Effect==effc)[1], "type"]=="Var."){
                addPrior <- data.frame("x"=x.var.prior, "y"=var.prior, "Effect"=effc, "type"=kdens[which(kdens$Effect==effc)[1], "type"], "group"="prior")
                kdens <- rbind(kdens, addPrior)
              }else if(kdens[which(kdens$Effect==effc)[1], "type"]=="St.Dev."){
                addPrior <- data.frame("x"=x.sd.prior, "y"=sd.prior, "Effect"=effc, "type"=kdens[which(kdens$Effect==effc)[1], "type"], "group"="prior")
                kdens <- rbind(kdens, addPrior)
              }else if(kdens[which(kdens$Effect==effc)[1], "type"]=="Cov."){
                addPrior <- data.frame("x"=density(cov.prior)$x, "y"=density(cov.prior)$y, "Effect"=effc, "type"=kdens[which(kdens$Effect==effc)[1], "type"], "group"="prior")
                kdens <- rbind(kdens, addPrior)
              }else if(kdens[which(kdens$Effect==effc)[1], "type"]=="Correl."){
                addPrior <- data.frame("x"=density(corr.prior)$x, "y"=density(corr.prior)$y, "Effect"=effc, "type"=kdens[which(kdens$Effect==effc)[1], "type"], "group"="prior")
                kdens <- rbind(kdens, addPrior)
              }
            }
            out$Covariances[[l]] <- ggplot(kdens, aes(x=x,y=y,group=group)) +
              xlab('') +
              ylab('Density') +
              geom_line(aes(color=type, linetype=group)) +
              facet_wrap(~Effect, scales='free') +
              theme_minimal()
          }else{
            out$Covariances[[l]] <- ggplot(kdens, aes(x=x,y=y)) +
              xlab('') +
              ylab('Density') +
              geom_line(aes(color=type, linetype=type)) +
              facet_wrap(~Effect, scales='free') +
              theme_minimal()
          }
        }

      }
  }
  nhc <- length(hc.idx <- grep('Beta_intern for ', names(hid)))
  if(nhc>0) {
    if(priors){
      cMargs <- joinMarginals(
        x$internal.marginals.hyperpar[hc.idx])
      cnames <- substring(names(x$internal.marginals.hyperpar)[hc.idx],16)
      cMargs$Effect <- factor(cnames[cMargs$m], cnames, cnames)
      cMargs$group <- "posterior"
      for(effa in unique(cMargs$Effect)){
        DatEffa <- cMargs[cMargs$Effect==effa,]
        addPrior <- data.frame("m"=1, "x"=x.assoc.prior, "y"=assoc.prior,
                               "Effect"=effa,  "group"="prior")
        cMargs <- rbind(cMargs, addPrior)
      }
      out$Associations <- ggplot(cMargs, aes(x=x,y=y, group=group, colour=group, linetype=group)) +
        xlab('') +
        ylab('Density') +
        geom_line() +
        facet_wrap(~Effect, scales='free') +
        theme_minimal()
    }else{
      cMargs <- joinMarginals(
        x$internal.marginals.hyperpar[hc.idx])
      cnames <- substring(names(x$internal.marginals.hyperpar)[hc.idx],16)
      cMargs$Effect <- factor(cnames[cMargs$m], cnames, cnames)
      out$Associations <- ggplot(cMargs, aes(x=x,y=y)) +
        xlab('') +
        ylab('Density') +
        geom_line() +
        facet_wrap(~Effect, scales='free') +
        theme_minimal()
    }
  }
  rnames <- names(x$summary.random)
  nbas <- length(bas.idx <- grep(
      '^baseline[0-9]+', rnames))
  nbasP <- length(basP.idx <- c(grep("weibullsurv", unlist(x$basRisk)), # number of parametric baseline risks
                                grep("exponentialsurv", unlist(x$basRisk)),
                                grep("dgompertzsurv", unlist(x$basRisk)),
                                grep("gompertzsurv", unlist(x$basRisk))))
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
      if(!is.null(x$strata[[i]])){
        Nstrata <- length(which(x$summary.random[[paste0("baseline",i,".hazard")]]$ID==0))
        Pstrata <- rep(1:Nstrata, each=length(x$summary.random[[paste0("baseline",i,".hazard")]]$ID)/Nstrata)
        BaselineValues <- rbind(BaselineValues,
                                cbind(time=x$summary.random[[paste0("baseline",i,".hazard")]]$ID,
                                      mean=BHmean,
                                      lower=BHlo,
                                      upper=BHup,
                                      S=as.numeric(paste0(i, ".", Pstrata))))
        jsr <- Reduce('rbind', lapply(1:nbas, function(k) {
          data.frame(x$summary.random[[bas.idx[k]]],
                     S=paste0('S',k, '.', Pstrata))
        }))
      }else{
        BaselineValues <- rbind(BaselineValues,
                                cbind(time=x$summary.random[[paste0("baseline",i,".hazard")]]$ID,
                                      mean=BHmean,
                                      lower=BHlo,
                                      upper=BHup,
                                      S=i))
        jsr <- Reduce('rbind', lapply(1:nbas, function(k) {
          data.frame(x$summary.random[[bas.idx[k]]],
                     S=paste0('S',k))
        }))
      }
    }
      colnames(jsr)  <- gsub('X0.', 'q', colnames(jsr), fixed=TRUE)
      out$Baseline <- ggplot(jsr, aes(x=ID)) +
          geom_ribbon(aes(ymin=BaselineValues[,"lower"],
                          ymax=BaselineValues[,"upper"]),
                      fill='grey70') +
          geom_line(aes(y=BaselineValues[,"mean"])) +
          xlab('Time') +
          ylab('Baseline risk') +
          facet_wrap(~S,  scales='free') +
        theme_minimal()
  }
  # if exists NLcovName then look at vector NLassoc et Lassoc et pour chaque
  # TRUE in Lassoc check NL assoc et grab les parametres correspondants
  # reconstruire d'abord les Linear puis les splines (vecotr avec mean slope
  # pspline pour linear et ajouter thera dans la second etape)
  if(length(grep('scopy', names(hid)))>0){
    NL_data <- NULL
    NL_data2 <- NULL
    NL_mean <- grep('(scopy mean)', names(hid))
    NL_slope <- grep('(scopy slope)', names(hid))
    NL_theta <- grep('(scopy theta)', names(hid))
    hc.idxNL <- grep('scopy', names(hid))
    numNL_mean <- length(substr(names(hid)[NL_mean], 1, 6))
    numNL_slope <- length(substr(names(hid)[NL_slope], 1, 6))
    numNL_theta <- length(substr(names(hid)[NL_theta], 1, 6))
    NLeffid <- sapply(names(x$summary.random), function(x) grep(x, names(hid)[hc.idxNL]))
    NLeff <- names(NLeffid)[which(sapply(NLeffid, length)>0)]

    if(methodNL=="sampling") Hnl <- INLA::inla.hyperpar.sample(NsampleNL, x)

    for(effNL in NLeff){
      k_NL <- as.integer(strsplit(strsplit(effNL, "_L")[[1]][2], "_S")[[1]][1])
      if(length(grep("CV", effNL)>0)){
        x_NLid <- grep(paste0("uv", k_NL), names(x$summary.random))
      }else if(length(grep("CS", effNL)>0)){
        x_NLid <- grep(paste0("us", k_NL), names(x$summary.random))
      }else if(length(grep("SRE", effNL)>0)){
        x_NLid <- grep(paste0("usre", k_NL), names(x$summary.random))
      } # CV_CS not done here
      # x_NLid <- grep(effNL, names(x$summary.random))
      xval <- x$summary.random[[x_NLid]]$mean# x$cov_NL[[k_NL]] #
      xval2 <- seq(min(xval), max(xval), len=1000)# seq(min(x$cov_NL[[k_NL]]), max(x$cov_NL[[k_NL]]), len=1000)
      # xval2 <- seq(range(x$summary.random[[x_NLid2]]$mean), len=1000)# seq(min(x$cov_NL[[k_NL]]), max(x$cov_NL[[k_NL]]), len=1000)
      if(methodNL=="analytical"){
        stop("WIP")
        # sf_NL <- smooth.spline(xval, x$summary.random[[effNL]]$mean)
        # sfupp_NL <- smooth.spline(xval, x$summary.random[[effNL]]$'0.025quant')
        # sflow_NL <- smooth.spline(xval, x$summary.random[[effNL]]$'0.975quant')
        # # sfupp_NL <- smooth.spline(xval, x$summary.random[[effNL]]$mean+1.96*x$summary.random[[effNL]]$sd)
        # # sflow_NL <- smooth.spline(xval, x$summary.random[[effNL]]$mean-1.96*x$summary.random[[effNL]]$sd)
        # NL_data <- rbind(NL_data, cbind("x"=sf_NL$x,
        #                                 "y"=sf_NL$y,
        #                                 "upper"=sfupp_NL$y,
        #                                 "lower"=sflow_NL$y, "Effect"=effNL))

      }else if(methodNL=="sampling"){
        nb <- length(grep(effNL, names(hid))) # number of splines parameters
        # xx.loc <- min(xval) + (max(xval)-min(xval)) * (0:(nb - 1))/(nb - 1)
        xx.loc <- x$range[[which(NLeff == effNL)]][1] + diff(x$range[[which(NLeff == effNL)]]) * seq(0, 1, len = nb)
        prop <- INLAjoint.scopy.define(nb)
        for(nsmp in 1:NsampleNL){
          # funNL <- splinefun(xx.loc, Hnl[nsmp, grep(effNL, names(hid))], method = "natural")
          funNL <- splinefun(xx.loc, prop$W %*% Hnl[nsmp, grep(effNL, names(hid))], method = "natural")
          if(NLeffectonly){
            if(!hr) NLval = funNL(xval2) else NLval = exp(funNL(xval2))
            NL_data2 <- cbind(NL_data2, NLval)
          }else{
            if(!hr) NLval = xval2*funNL(xval2) else NLval = exp(xval2*funNL(xval2))
            NL_data2 <- cbind(NL_data2, NLval)
          }
        }
        statsNL <- apply(NL_data2, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))
        NL_data <- rbind(NL_data, cbind("x"=xval2,
                                        "y"=statsNL[2,],
                                        "upper"=statsNL[1,],
                                        "lower"=statsNL[3,], "Effect"=effNL))
      }
      NL_data2 <- NULL
    }
    NL_data <- as.data.frame(NL_data)
    NL_data[, 1:4] <- apply(NL_data[, 1:4], 2, as.numeric)


    #
    #     vals <- fun(xHs[, i])
    #     lines(xHs[,i], xHs[,i]*vals, col=4, lty=2)
    #
    #     NL_data[, 1:4] <- apply(NL_data[, 1:4], 2, as.numeric)
    #     plot(NL_data$x, NL_data$y, type="o", pch=19)
    #     lines(NL_data$x, NL_data$lower, lty=2)
    #     lines(NL_data$x, NL_data$upper, lty=2)
    #
    #     plot(xval, x$summary.random[[effNL]]$mean, pch=19, cex=0.5)
    #     points(xval, x$summary.random[[effNL]]$'0.025quant', col=2, pch=19, cex=0.4)
    #     points(xval, x$summary.random[[effNL]]$'0.975quant', col=4, pch=19, cex=0.4)
    #
    #     plot(xval, x$summary.linear.predictor$mean[which(!is.na(x$.args$data$Yjoint$y1..coxph))], pch=19, cex=0.5)
    #     points(xval, x$summary.random[[effNL]]$'0.025quant', col=2, pch=19, cex=0.4)
    #     points(xval, x$summary.random[[effNL]]$'0.975quant', col=4, pch=19, cex=0.4)

    out$NL_Associations <- ggplot(NL_data, aes(x=x,y=y)) +
      xlab('') +
      ylab('Effect') +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70")+
      geom_line() +
      facet_wrap(~Effect, scales='free') +
      theme_minimal()
  }
  if(nbasP>0){
    HW0 <- function(t, lambda, alpha){ # risk function Weibull variant 0 (also exponential for alpha=1)
      res = lambda*alpha*t^(alpha-1)
    }
    HW1 <- function(t, lambda, alpha){ # risk function Weibull variant 1
      res = lambda*alpha*(lambda*t)^(alpha-1)
    }
    HG <- function(t, lambda, alpha){ # risk function Gompertz
      res = lambda*alpha*exp(alpha*t)
    }
    HDG <- function(t, lambda, alpha){ # risk function Defective Gompertz
      res = lambda*abs(alpha)*exp(abs(alpha)*t) / (exp(abs(alpha)*t) - 1 + exp(-lambda/abs(alpha)))
    }
    BHM <- NULL # baseline risk marginals
    nbl <- 1 # to keep track of baseline risk in case of multiple parametric survival outcomes
    nbl2 <- 1
    BaselineValuesP <- NULL
    if(x$.args$control.compute$config==TRUE){
      NSAMPLES = 500
      SEL <- sapply(paste0(rownames(x$summary.fixed)[grep("Intercept_S", rownames(x$summary.fixed))]), function(x) x=1, simplify=F)
      SAMPLES <- INLA::inla.posterior.sample(NSAMPLES, x, selection=SEL)
    }else if(x$.args$control.compute$config=="lite"){
      NSAMPLES = 500
      SEL <- sapply(paste0(rownames(x$summary.fixed)[grep("Intercept_S", rownames(x$summary.fixed))]), function(x) x=1, simplify=F)
      SAMPLES <- INLA::inla.rjmarginal(NSAMPLES, x)
      SAMPLESH <- INLA::inla.hyperpar.sample(NSAMPLES, x)
    }
    ctBP <- 1 # counter for baseline parametric risks
    for(i in 1:(nbas+nbasP)){
      if(i %in% basP.idx){
        # if(nbasP==1){
        if(length(x$.args$family)==1){
          maxTime <- max(na.omit(x$.args$data$Yjoint$time))
        }else{
          maxTime <- max(na.omit(x$.args$data$Yjoint[which(sapply(x$.args$data$Yjoint, class)=="inla.surv")][[ctBP]]$time))
        }
        # }else if(nbasP>1){
        #   maxTime <- max(na.omit(x$.args$data$Yjoint[which(sapply(x$.args$data$Yjoint, class)=="inla.surv")][[i]]$time))
        # }
        timePts <- seq(0, maxTime, len=500)
        timePts2 <- timePts#c(timePts[2], timePts[-1]) # avoid computing parametric baseline risk at time 0 since it's always 0
        Variant_i <- x$.args$control.family[[i]]$variant
        if(x$basRisk[[i]]=="exponentialsurv"){
          BHM <- append(BHM, list(INLA::inla.tmarginal(function(x) exp(x),
                                                 x$marginals.fixed[grep("Intercept_S", names(x$marginals.fixed))][[ctBP]])))
          names(BHM)[nbl] <- paste0("Exponential (rate)_S", i)
          if(x$.args$control.compute$config==TRUE){ # config set to TRUE then we can compute uncertainty
            curves_SMP <- sapply(1:NSAMPLES, function(x) HW0(t=timePts2, lambda=exp(SAMPLES[[x]]$latent[nbl]), alpha=1))
            QUANT <- apply(curves_SMP, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))
            values_i <- t(QUANT)
          }else{
            values_i <- HW0(t=timePts2, lambda=exp(x$summary.fixed[grep("Intercept_S", names(x$marginals.fixed)), "mean"][[ctBP]]), alpha=1)
          }
          name_i <- paste0("Exponential baseline risk (S", nbl, ")")
        }else if(x$basRisk[[i]]=="weibullsurv"){
          BHM <- append(BHM, list(INLA::inla.tmarginal(function(x) exp(x),
                                                 x$marginals.fixed[grep("Intercept_S", names(x$marginals.fixed))][[ctBP]])))
          names(BHM)[nbl2] <- paste0("Weibull (scale)_S", i)
          BHM <- append(BHM, list(x$marginals.hyperpar[grep("weibull", names(x$marginals.hyperpar))][[nbl]]))
          names(BHM)[nbl2+1] <- paste0("Weibull (shape)_S", i)
          if(x$.args$control.compute$config==TRUE){ # config set to TRUE then we can compute uncertainty
            if(Variant_i==0){
              curves_SMP <- sapply(1:NSAMPLES, function(x) HW0(t=timePts2, lambda=exp(SAMPLES[[x]]$latent[nbl]), alpha=SAMPLES[[x]]$hyperpar[grep("weibull", names(SAMPLES[[x]]$hyperpar))][nbl]))
            }else if(Variant_i==1){
              curves_SMP <- sapply(1:NSAMPLES, function(x) HW1(t=timePts2, lambda=exp(SAMPLES[[x]]$latent[nbl]), alpha=SAMPLES[[x]]$hyperpar[grep("weibull", names(SAMPLES[[x]]$hyperpar))][nbl]))
            }
            QUANT <- apply(curves_SMP, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))
            values_i <- t(QUANT)
          }else if(x$.args$control.compute$config=="lite"){
            if(Variant_i==0){
              curves_SMP <- sapply(1:NSAMPLES, function(x) HW0(t=timePts2, lambda=exp(SAMPLES$samples[,x])[grep(paste0("Intercept_S", nbl), names(SAMPLES$samples[,x]))],
                                                               alpha=SAMPLESH[x,][grep("weibullsurv", names(SAMPLESH[x,]))][nbl]))
            }else if(Variant_i==1){
              curves_SMP <- sapply(1:NSAMPLES, function(x) HW1(t=timePts2, lambda=exp(SAMPLES$samples[,x])[grep(paste0("Intercept_S", nbl), names(SAMPLES$samples[,x]))],
                                                               alpha=SAMPLESH[x,][grep("weibullsurv", names(SAMPLESH[x,]))][nbl]))
            }
            QUANT <- apply(curves_SMP, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))
            values_i <- t(QUANT)
            # if(Variant_i==0){
            #   values_i <- HW0(t=timePts2, lambda=exp(x$summary.fixed[grep("Intercept_S", rownames(x$summary.fixed)), "mean"][[i]]),
            #                   alpha=x$summary.hyperpar[grep("weibull", names(x$marginals.hyperpar)), "mean"][[nbl]])
            # }else if(Variant_i==1){
            #   values_i <- HW1(t=timePts2, lambda=exp(x$summary.fixed[grep("Intercept_S", rownames(x$summary.fixed)), "mean"][[i]]),
            #                   alpha=x$summary.hyperpar[grep("weibull", names(x$marginals.hyperpar)), "mean"][[nbl]])
            # }
          }
          name_i <- paste0("Weibull baseline risk (S", nbl, ")")
        }else if(x$basRisk[[i]]=="gompertzsurv"){
          BHM <- append(BHM, list(INLA::inla.tmarginal(function(x) exp(x),
                                                 x$marginals.fixed[grep("Intercept_S", names(x$marginals.fixed))][[ctBP]])))
          names(BHM)[nbl2] <- paste0("Gompertz (scale)_S", i)
          BHM <- append(BHM, list(x$marginals.hyperpar[grep("alpha parameter for Gompertz-surv", names(x$marginals.hyperpar))][[nbl]]))
          names(BHM)[nbl2+1] <- paste0("Gompertz (alpha)_S", i)
          if(x$.args$control.compute$config==TRUE){ # config set to TRUE then we can compute uncertainty
            curves_SMP <- sapply(1:NSAMPLES, function(x) HG(t=timePts2, lambda=exp(SAMPLES[[x]]$latent[nbl]),
                                                            alpha=SAMPLES[[x]]$hyperpar[grep("alpha parameter for Gompertz-surv", names(SAMPLES[[x]]$hyperpar))][nbl]))
            QUANT <- apply(curves_SMP, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))
            values_i <- t(QUANT)
          }else if(x$.args$control.compute$config=="lite"){
            curves_SMP <- sapply(1:NSAMPLES, function(x) HG(t=timePts2, lambda=exp(SAMPLES$samples[,x])[grep(paste0("Intercept_S", nbl), names(SAMPLES$samples[,x]))],
                                                            alpha=SAMPLESH[x,][grep("alpha parameter for Gompertz-surv", names(SAMPLESH[x,]))][nbl]))
            QUANT <- apply(curves_SMP, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))
            values_i <- t(QUANT)
          }
          name_i <- paste0("Gompertz baseline risk (S", nbl, ")")
        }else if(x$basRisk[[i]]=="dgompertzsurv"){
          BHM <- append(BHM, list(INLA::inla.tmarginal(function(x) exp(x),
                                                 x$marginals.fixed[grep("Intercept_S", names(x$marginals.fixed))][[ctBP]])))
          names(BHM)[nbl2] <- paste0("Defective Gompertz (scale)_S", i)
          BHM <- append(BHM, list(x$marginals.hyperpar[grep("alpha parameter for dGompertz-surv", names(x$marginals.hyperpar))][[nbl]]))
          names(BHM)[nbl2+1] <- paste0("Defective Gompertz (alpha)_S", i)
          if(x$.args$control.compute$config==TRUE){ # config set to TRUE then we can compute uncertainty
            curves_SMP <- sapply(1:NSAMPLES, function(x) HDG(t=timePts2, lambda=exp(SAMPLES[[x]]$latent[nbl]),
                                                             alpha=SAMPLES[[x]]$hyperpar[grep("alpha parameter for dGompertz-surv", names(SAMPLES[[x]]$hyperpar))][nbl]))
            QUANT <- apply(curves_SMP, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))
            values_i <- t(QUANT)
          }else if(x$.args$control.compute$config=="lite"){
            curves_SMP <- sapply(1:NSAMPLES, function(x) HDG(t=timePts2, lambda=exp(SAMPLES$samples[,x])[grep(paste0("Intercept_S", nbl), names(SAMPLES$samples[,x]))],
                                                             alpha=SAMPLESH[x,][grep("alpha parameter for dGompertz-surv", names(SAMPLESH[x,]))][nbl]))
            QUANT <- apply(curves_SMP, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))
            values_i <- t(QUANT)
          }
          name_i <- paste0("Defective Gompertz baseline risk (S", nbl, ")")
        }
        # Increment counters based on number of parameters
        if(x$basRisk[[i]]=="exponentialsurv"){
          nbl2 <- nbl2+1  # exponential has only 1 parameter (rate)
        }else{
          nbl2 <- nbl2+2  # weibull, gompertz, dgompertz have 2 parameters (scale + shape/alpha)
        }
        nbl <- nbl+1
        ctBP <- ctBP+1
        BaselineValuesP <- rbind(BaselineValuesP, data.frame(timePts, values_i, name_i)[-1,])
      }
    }
    if(x$.args$control.compute$config==TRUE | x$.args$control.compute$config=="lite"){ # with uncertainty #### NEEDS FIX (check JM1 example bayes surv analysis with INLA => uncertainty leads to extreme values (maybe cut tails?))
      colnames(BaselineValuesP) <- c("x", "lower", "y", "upper", "Effect")
      out$BaselineP <- ggplot(data=BaselineValuesP) +
        geom_ribbon(aes(x=x, ymin=lower,
                        ymax=upper),
                    fill='grey70') +
        geom_line(aes(x=x, y=y)) +
        xlab('Time') +
        ylab('Baseline risk') +
        facet_wrap(~Effect,  scales='free') +
        theme_minimal()
    }else{ # only means
      # colnames(BaselineValues) <- c("x", "y", "Effect")
      # out$Baseline <- ggplot(BaselineValues, aes(x=x, y=y)) +
      #   geom_line(aes(y=BaselineValues[,"y"])) +
      #   xlab('Time') +
      #   ylab('Baseline risk') +
      #   facet_wrap(~Effect,  scales='free')
    }
    sMargs <- joinMarginals(BHM, trim=FALSE)
    snames <- names(BHM)
    sMargs$Effect <- factor(snames[sMargs$m], snames, snames)
    out$ParamBaseline <- ggplot(sMargs, aes(x=x,y=y)) +
      xlab('') +
      ylab('Density') +
      geom_line() +
      facet_wrap(~Effect, scales='free') +
      theme_minimal()
  }

  # RW2 trajectory plots
  if(!is.null(x$rw2_info)) {
    out$RW2_Trajectories <- list()

    for(k in seq_along(x$rw2_info)) {
      rw2_k <- x$rw2_info[[k]]
      if(!is.null(rw2_k) && !is.null(rw2_k$group_map)) {
        n_samples <- 1000

        # Posterior samples
        original_class <- class(x)
        class(x) <- "inla"
        samples <- tryCatch({
          INLA::inla.posterior.sample(n_samples, x)
        }, error = function(e) NULL)
        class(x) <- original_class

        if(!is.null(samples)) {
          # Extract fixed effect samples
          extract_fixed_effect <- function(samples, coef_name) {
            latent_names <- rownames(samples[[1]]$latent)
            full_name <- paste0(coef_name, ":1")
            idx <- which(latent_names == full_name)
            if (length(idx) == 0) return(rep(0, length(samples)))
            sapply(samples, function(s) s$latent[idx[1], 1])
          }

          marker_suffix <- paste0("_L", k)
          time_var <- rw2_k$time_var
          n_groups <- rw2_k$n_groups
          group_map <- rw2_k$group_map

          # Get intercept samples
          intercept_samp <- extract_fixed_effect(samples, paste0("Intercept", marker_suffix))

          # Extract group-specific fixed effects from model
          group_fixed_effects <- vector("list", n_groups)
          for(g in seq_len(n_groups)) {
            group_fixed_effects[[g]] <- rep(0, n_samples)
          }

          # Parse group expression to identify grouping variables
          if(!is.null(rw2_k$group_expr)) {
            group_expr <- rw2_k$group_expr
            group_vars <- unique(unlist(strsplit(gsub("[*:()]", " ", group_expr), " ")))
            group_vars <- group_vars[nchar(group_vars) > 0]

            # Extract fixed effect samples for each grouping variable
            var_samples <- list()
            for(var in group_vars) {
              var_cols <- grep(paste0("^", var, ".*", marker_suffix, "$"),
                              names(x$marginals.fixed), value = TRUE)
              if(length(var_cols) > 0) {
                # Main effect
                main_col <- var_cols[!grepl("\\.", var_cols)]
                if(length(main_col) > 0) {
                  var_samples[[var]]$main <- extract_fixed_effect(samples, main_col[1])
                }
                # Interaction effects
                int_cols <- var_cols[grepl("\\.", var_cols)]
                for(int_col in int_cols) {
                  int_name <- gsub(marker_suffix, "", int_col)
                  var_samples[[var]][[int_name]] <- extract_fixed_effect(samples, int_col)
                }
              }
            }

            # Compute group-specific fixed effects based on group labels
            for(g in seq_len(n_groups)) {
              label <- as.character(group_map$group_label[g])
              # Parse label to identify which variables are active
              for(var in names(var_samples)) {
                if(!is.null(var_samples[[var]]$main) && grepl(var, label, fixed = TRUE)) {
                  group_fixed_effects[[g]] <- group_fixed_effects[[g]] + var_samples[[var]]$main
                }
              }
              # Add interaction terms if applicable
              for(var in names(var_samples)) {
                int_terms <- names(var_samples[[var]])[names(var_samples[[var]]) != "main"]
                for(int_term in int_terms) {
                  if(grepl(int_term, label, fixed = TRUE)) {
                    group_fixed_effects[[g]] <- group_fixed_effects[[g]] + var_samples[[var]][[int_term]]
                  }
                }
              }
            }
          }

          # Extract RW2 components
          configs <- x$misc$configs$contents
          time_l_idx <- which(configs$tag == paste0(time_var, marker_suffix))

          if(length(time_l_idx) > 0) {
            time_l_start <- configs$start[time_l_idx]
            time_l_length <- configs$length[time_l_idx]

            rw_summary <- x$summary.random[[paste0(time_var, marker_suffix)]]
            if ("ID" %in% colnames(rw_summary)) {
              n_unique_times <- length(unique(rw_summary$ID))
              if (time_l_length %% n_groups == 0) {
                n_per_group <- time_l_length / n_groups
              } else {
                n_per_group <- n_unique_times
              }
              time_grid <- unique(rw_summary$ID)[seq_len(min(n_per_group, n_unique_times))]
            } else {
              n_per_group <- round(time_l_length / n_groups)
              time_grid <- seq_len(n_per_group)
            }
            if(is.factor(time_grid)) time_grid <- as.character(time_grid)
            time_grid <- as.numeric(time_grid)

            # Compute trajectories: intercept + group fixed effects + RW2 component
            traj_list <- vector("list", n_groups)
            for(g in seq_len(n_groups)) {
              traj_list[[g]] <- matrix(0, nrow = n_samples, ncol = n_per_group)
            }

            for (i in seq_len(n_samples)) {
              rw_vals <- samples[[i]]$latent[time_l_start:(time_l_start + time_l_length - 1), 1]
              for(g in seq_len(n_groups)) {
                start_idx <- (g-1) * n_per_group + 1
                end_idx <- g * n_per_group
                rw_g <- rw_vals[start_idx:end_idx]
                traj_list[[g]][i, ] <- intercept_samp[i] + group_fixed_effects[[g]][i] + rw_g
              }
            }

            # Compute quantiles
            compute_quantiles <- function(traj_matrix) {
              list(
                median = apply(traj_matrix, 2, median),
                lower = apply(traj_matrix, 2, quantile, probs = 0.025),
                upper = apply(traj_matrix, 2, quantile, probs = 0.975)
              )
            }

            # Build plot data
            plot_data_list <- vector("list", n_groups)
            for(g in seq_len(n_groups)) {
              quant_g <- compute_quantiles(traj_list[[g]])
              plot_data_list[[g]] <- data.frame(
                time = time_grid,
                fitted_mean = quant_g$median,
                fitted_lower = quant_g$lower,
                fitted_upper = quant_g$upper,
                group_label = group_map$group_label[g]
              )
            }

            plot_data <- do.call(rbind, plot_data_list)
            plot_data$group_factor <- factor(plot_data$group_label,
                                             levels = unique(plot_data$group_label))

            # Line types per group
            n_linetypes <- length(unique(plot_data$group_factor))
            linetypes <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")[1:n_linetypes]

            p <- ggplot(plot_data, aes(x = time, color = group_factor, linetype = group_factor)) +
              geom_line(aes(y = fitted_lower), linewidth = 0.5, alpha = 0.5) +
              geom_line(aes(y = fitted_upper), linewidth = 0.5, alpha = 0.5) +
              geom_line(aes(y = fitted_mean), linewidth = 1.2) +
              scale_linetype_manual(values = linetypes) +
              labs(
                title = paste("RW2 Trajectories - Marker", k),
                x = time_var,
                y = "Fitted value",
                color = "Group",
                linetype = "Group"
              ) +
              theme_bw() +
              theme(
                legend.position = "bottom",
                plot.title = element_text(hjust = 0.5, face = "bold")
              )

            out$RW2_Trajectories[[k]] <- p
          }
        }
      }
    }

    if(length(out$RW2_Trajectories) == 0) {
      out$RW2_Trajectories <- NULL
    }
  }

  out <- out[!sapply(out, is.null)]
  class(out) <- c("plot.INLAjoint", "list")
  return(out)
}

