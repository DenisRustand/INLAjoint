# INLAjoint

## Joint modeling multivariate longitudinal and time-to-event outcomes with INLA

INLAjoint is a package that fits joint models for multivariate longitudinal markers (with various distributions available) and survival outcomes (possibly accounting for competing risks) with Integrated Nested Laplace Approximations (INLA). The flexible and user friendly function joint() facilitates the use of the fast and reliable inference technique implemented in INLA package for joint modeling. More details are given in the help page of the joint function (accessible via ?joint in the R console), the vignette associated to the joint() function (accessible via vignette("INLAjoint") in the R console) and the following article: https://arxiv.org/abs/2203.06256

## Install

devtools::install_github('DenisRustand/INLAjoint', build_vignettes = TRUE)

Note that INLA is required, you can install it with:
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

More details at https://www.r-inla.org/download-install

Finally, it is possible to use this package with pardiso, which provides a high performance computing environment with parallel computing support using OpenMP (not avail. on windows). See https://pardiso-project.org/r-inla/ for more informations.

