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
#'
#' @export

plot.INLAjoint <- function(jres, sdcor=FALSE, ...) {

    stopifnot(is.logical(sdcor))
    out <- list(
        Outcomes=NULL, Covariances=NULL,
        Associations=NULL, Baseline=NULL, Random=NULL)

    trimMarginal <- function(m, p=0.001) {
        ## trim the tails of a marginal
        s <- INLA:::inla.smarginal(m)
        ab <- INLA:::inla.qmarginal(c(p, 1-p), s)
        ii <- which((s$x>=ab[1])&(s$x<=ab[2]))
        return(list(x=s$x[ii], y=s$y[ii]))
    }
    joinMarginals <- function(listofmarginals) {
        ## join the marginals into one long data with indicator
        Reduce('rbind', lapply(1:length(listofmarginals), function(k) {
            data.frame(m=k, trimMarginal(listofmarginals[[k]]))
        }))
    }

    out.patt <- '_[LS][0-9]+$'
    hhid <- sapply(jres$internal.marginals.hyperpar, attr, 'hyperid')
    hid <- sapply(strsplit(hhid, '|', fixed=TRUE), tail, 1)
    x.n <- length(x.names <- names(jres$marginals.fixed))
    if(x.n>0) {
        x.psub <- regexpr(out.patt, x.names)
        x.group <- substring(x.names, x.psub+1)
        x.sname <- substr(x.names, 1, x.psub-1)
        xMargs <- joinMarginals(jres$marginals.fixed)
        xMargs$Effect <- x.sname[xMargs$m]
        xMargs$Outcome <- x.group[xMargs$m]
        lfamilies0 <- c('gaussian', 'lognomal')
        hl.jj <- which(unlist(jres$famL) %in% lfamilies0)
        nhl <- length(hl.jj)
        hs.jj <- grep('baseline[1-9]', hid)
        nhs <- length(hs.jj)
        if((nhl+nhs)>0) {
            thMargs <- joinMarginals(lapply(
                jres$internal.marginals.hyperpar[c(hl.jj, hs.jj)],
                function(m) inla.tmarginal(function(x) exp(-x/(1+sdcor)), m)))
            thMargs$Effect <-
                c(rep(paste0('Residual ', c('Var.', 'S.D.')[sdcor+1]), nhl),
                  rep(paste0('Baseline', c('Var.', 'S.D.')[sdcor+1]), nhs)
                  )[thMargs$m]
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
        kdsamples <- inla.iidkd.sample(
            2e4, jres, hid[k1], return.cov=!sdcor)
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
            out$Covariances <- ggplot(kdens, aes(x=x,y=y)) +
                xlab('') +
                ylab('Density') +
                geom_line(aes(color=type)) +
                facet_wrap(~Effect, scales='free')
        } else {
            warning('Something wrong with', kid[k1], 'happened!')
        }
    }

    nhc <- length(hc.idx <- grep('Beta_intern for ', names(hid)))
    if(nhc>0) {
        cMargs <- joinMarginals(
            jres$internal.marginals.hyperpar[hc.idx])
        cnames <- substring(names(jres$internal.marginals.hyperpar)[hc.idx],16)
        cMargs$Effect <- cnames[cMargs$m]
        out$Associations <- ggplot(cMargs, aes(x=x,y=y)) +
            xlab('') +
            ylab('Density') +
            geom_line() +
            facet_wrap(~Effect, scales='free')
    }

    rnames <- names(jres$summary.random)
    nbas <- length(bas.idx <- grep(
        '^baseline[0-9]+', rnames))
    if(nbas>0) {
        jsr <- Reduce('rbind', lapply(1:nbas, function(k) {
            data.frame(jres$summary.random[[bas.idx[k]]],
                       S=paste0('S',k))
        }))
        colnames(jsr)  <- gsub('X0.', 'q', colnames(jsr), fixed=TRUE)
        out$Random <- ggplot(jsr, aes(x=ID)) +
            geom_ribbon(aes(ymin=exp(q025quant),
                            ymax=exp(q975quant)),
                        fill='grey70') +
            geom_line(aes(y=exp(mean))) +
            xlab('Time') +
            ylab('Baseline risk') +
            facet_wrap(~S,  scales='free')
    }

    return(out[!sapply(out, is.null)])
}
