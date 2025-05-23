# INLAjoint

## Joint modeling multivariate longitudinal and time-to-event outcomes with INLA

INLAjoint is a package that fits joint models for multivariate longitudinal markers (with various distributions available) and survival outcomes (possibly accounting for competing risks and multi-state) with Integrated Nested Laplace Approximations (INLA). The flexible and user friendly function joint() facilitates the use of the fast and reliable inference technique implemented in INLA package for joint modeling. More details are given in the help page of the joint function (accessible via ?joint in the R console), the vignette associated to the joint() function (accessible via vignette("INLAjoint") in the R console).

## Install
install.packages("R.rsp") # (only for the vignette)

devtools::install_github('DenisRustand/INLAjoint', build_vignettes = TRUE)

Note that INLA is required, you can install it with:
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)

More details at https://www.r-inla.org/download-install

Contact: INLAjoint@gmail.com

References:

- Rustand, D., van Niekerk, J., Krainski, E. T., & Rue, H. (2024). Joint Modeling of Multivariate Longitudinal and Survival Outcomes with the R package INLAjoint. arXiv preprint arXiv:2402.08335.
https://arxiv.org/abs/2402.08335

- Rustand, D., van Niekerk, J., Krainski, E. T., Rue, H., & Proust-Lima, C. (2024). Fast and flexible inference for joint models of multivariate longitudinal and survival data using integrated nested Laplace approximations. Biostatistics, kxad019.
https://doi.org/10.1093/biostatistics/kxad019

- Alvares, D., Van Niekerk, J., Krainski, E.T., Rue, H. and Rustand, D., 2024. Bayesian survival analysis with INLA. Statistics in medicine, 43(20), pp.3975-4010.
https://doi.org/10.1002/sim.10160
