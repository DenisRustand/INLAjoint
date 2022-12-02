#' Plot the output from a multivariate joint model for longitudinal and/or survival data
#'
#' @description This function provide plots for the output of a
#' multivariate joint model for longitudinal and/or survival data.
#' The output can be stored into an object and manipulated as
#' a list of ggplot outputs, described bellow.
#'
#' @param jres an object with the output of the the \link{joint} function
#' @param sdcor logical indicating if the random effects
#' correlation are to be shown. If FALSE the covariance is shown.
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
#' @export

plot.INLAjoint <- function(jres, sdcor=FALSE, ...) {

    stopifnot(is.logical(sdcor))
    out <- list(
        Outcomes=NULL, Covariances=NULL,
        Associations=NULL, Baseline=NULL, Random=NULL)

    trimMarginal <- function(m, p=0.001) {
            m <- INLA:::inla.smarginal(m)
            ab <- INLA:::inla.qmarginal(c(p, 1-p), m)
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
    hhid <- sapply(jres$internal.marginals.hyperpar, attr, 'hyperid')
    hid <- sapply(strsplit(hhid, '|', fixed=TRUE), tail, 1)
    x.n <- length(x.names <- names(jres$marginals.fixed[-grep("Intercept_S", names(jres$marginals.fixed))]))
    if(x.n>0) {
        x.psub <- regexpr(out.patt, x.names)
        x.group <- substring(x.names, x.psub+1)
        xMargs <- joinMarginals(jres$marginals.fixed[-grep("Intercept_S", names(jres$marginals.fixed))])
        xMargs$Effect <- factor(x.names[xMargs$m], x.names, x.names)
        xMargs$Outcome <- x.group[xMargs$m]
        lfamilies0 <- c(
            'gaussian', 'lognormal', 'gamma', 'beta'
        ) ## others to be added
        hl.jj <- which(unlist(jres$famLongi) %in% lfamilies0)
        nhl <- length(hl.jj)
        hs.jj <- grep('baseline[1-9]', hid)
        nhs <- length(hs.jj)
        if((nhl+nhs)>0) {
            thMargs <- joinMarginals(lapply(
                jres$internal.marginals.hyperpar[c(seq_len(nhl), hs.jj)],
                function(m) inla.tmarginal(function(x) exp(-x/(1+sdcor)), m)),
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
        out$Outcomes <- lapply(
            split(xMargs, xMargs$Outcome), function(d)
                ggplot(d, aes(x=x, y=y)) +
                xlab('') +
                ylab('Density') +
                geom_line() +
                facet_wrap(~Effect, scales='free'))
    }

    nhk <- length(hd.idx <- grep('^Theta[0-9]+ for ', names(hid)))
    if(nhk>0) {
        k1 <- grep('Theta1 for ', names(hid))
        class(jres) <- 'inla'
        out$Covariances <- vector('list', length(k1))
        names(out$Covariances) <- paste0('L', 1:length(k1))
        for (l in 1:length(k1)) {
          kdsamples <- inla.iidkd.sample(
              2e4, jres, hid[k1[l]], return.cov=!sdcor)
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
                  jres$REstruc, jres$REstruc), 1,
                  function(x) paste(unique(x), collapse=':'))
              kdens$Effect <- factor(
                  kdnames[kdens$m],
                  unique(kdnames))
              kdens$type <- factor(
                  kdens$ldu, 1:3,
                  c(c('Var.', 'St.Dev.')[sdcor+1],
                    'Correl.', 'Correl.'))
              out$Covariances[[l]] <- ggplot(kdens, aes(x=x,y=y)) +
                  xlab('') +
                  ylab('Density') +
                  geom_line(aes(color=type)) +
                  facet_wrap(~Effect, scales='free')
          } else {
              warning('Something wrong with', kid[k1], 'happened!')
          }
        }
    }
    nhc <- length(hc.idx <- grep('Beta_intern for ', names(hid)))
    if(nhc>0) {
        cMargs <- joinMarginals(
            jres$internal.marginals.hyperpar[hc.idx])
        cnames <- substring(names(jres$internal.marginals.hyperpar)[hc.idx],16)
        cMargs$Effect <- factor(cnames[cMargs$m], cnames, cnames)
        out$Associations <- ggplot(cMargs, aes(x=x,y=y)) +
            xlab('') +
            ylab('Density') +
            geom_line() +
            facet_wrap(~Effect, scales='free')
    }
    rnames <- names(jres$summary.random)
    nbas <- length(bas.idx <- grep(
        '^baseline[0-9]+', rnames))
    nbasP <- length(c(grep("weibullsurv", unlist(jres$basRisk)), # number of parametric baseline risks
                      grep("exponentialsurv", unlist(jres$basRisk))))
    if(nbas>0) {
      BaselineValues <- NULL
      for(i in 1:nbas){
        BHmean <- NULL
        BHlo <- NULL
        BHup <- NULL
        for(j in 1:length(jres$marginals.random[[i]])){
          # m <- inla.smarginal(jres$marginals.random[[i]][[j]])
          # ab <- inla.qmarginal(c(0.001, 0.999), m)
          # ii <- which((m$x>=ab[1]) & (m$x<=ab[2]))
          # m$x <- m$x[ii]
          # m$y <- m$y[ii]
          Mm <- inla.qmarginal(c(0.025, 0.5, 0.975), jres$marginals.random[[i]][[j]])
          BHmean <- c(BHmean, exp(Mm[2]))
          BHlo <- c(BHlo, exp(Mm[1]))
          BHup <- c(BHup, exp(Mm[3]))
        }
        BaselineValues <- rbind(BaselineValues,
                                cbind(time=jres$summary.random[[paste0("baseline",i,".hazard")]]$ID,
                                mean=BHmean,
                                lower=BHlo,
                                upper=BHup,
                                S=i))
      }
        jsr <- Reduce('rbind', lapply(1:nbas, function(k) {
            data.frame(jres$summary.random[[bas.idx[k]]],
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
      BHM <- NULL # baseline risk marginals
      nbl <- 1 # to keep track of baseline risk in case of multiple parametric survival outcomes
      nbl2 <- 1
      for(i in 1:(nbas+nbasP)){
        if(jres$basRisk[[i]]=="exponentialsurv"){
          BHM <- append(BHM, list(inla.tmarginal(function(x) exp(x),
                                                 jres$marginals.fixed[grep("Intercept_S", names(jres$marginals.fixed))][[i]])))
          names(BHM)[nbl2] <- paste0("Exponential (rate)_S", i)
          nbl2 <- nbl2+1
        }else if(jres$basRisk[[i]]=="weibullsurv"){
          BHM <- append(BHM, list(inla.tmarginal(function(x) exp(x),
                                                 jres$marginals.fixed[grep("Intercept_S", names(jres$marginals.fixed))][[i]])))
          names(BHM)[nbl2] <- paste0("Weibull (scale)_S", i)
          BHM <- append(BHM, list(jres$marginals.hyperpar[grep("weibull", names(jres$marginals.hyperpar))][[nbl]]))
          names(BHM)[nbl2+1] <- paste0("Weibull (shape)_S", i)
          nbl2 <- nbl2+2
          nbl <- nbl+1
        }
      }
      sMargs <- joinMarginals(BHM)
      snames <- names(BHM)
      sMargs$Effect <- factor(snames[sMargs$m], snames, snames)
      out$Baseline <- ggplot(sMargs, aes(x=x,y=y)) +
        xlab('') +
        ylab('Density') +
        geom_line() +
        facet_wrap(~Effect, scales='free')
    }
    return(out[!sapply(out, is.null)])
}
