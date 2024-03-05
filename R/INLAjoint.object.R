#' Fitted \code{joint} object
#'
#' @description An object of class \code{INLAjoint} returned by the \code{joint}
#'   function that fits a joint model to multivariate longitudinal and
#'   time-to-event data. The following functions can apply to objects of this
#'   class: \code{plot}, \code{print}, \code{summary} and \code{priors.used}.
#'
#' @keywords joint model multivariate longitudinal survival competing risks
#' @seealso \code{\link{joint}}.
#' @return A list with the following components: \describe{
#'   \item{\code{names.fixed}}{a vector with the name of the fixed effects of
#'   the model. The corresponding submodel is indicated by the suffix including
#'   a letter and a number ("L" for longitudinal and "S" for survival).}
#'   \item{\code{summary.fixed}}{summary statistics for the fixed effects of
#'   the model. The summary statistics sorted by longitudinal and
#'   survival components are available by applying the \code{summary} function
#'   to the \code{INLAjoint} object.}
#'   \item{\code{summary.fixed}}{marginals for the fixed effects of the model.}
#'   \item{\code{mlik}}{log marginal-likelihood.}
#'   \item{\code{cpo}}{Conditional Predictive Ordinate.}
#'   \item{\code{gcpo}}{Group-Conditional Predictive Ordinate.}
#'   \item{\code{po}}{Predictive ordinate.}
#'   \item{\code{waic}}{Widely applicable Bayesian information criterion}
#'   \item{\code{model.random}}{a vector with the name of the random parameters of
#'   the model, possibly including the following components: \describe{
#'   \item{\code{RW1 model and RW2 model}}{Random walk of order 1 or 2
#'   corresponding to Bayesian smoothing splines for the baseline hazard risk}
#'   \item{\code{IID model}}{Univariate random effect.}
#'   \item{\code{IIDKD model}}{Multivariate random effects.}
#'   \item{\code{Copy}}{association parameter.}}}
#'   \item{\code{summary.random}}{summary statistics for the random parameters of the model.}
#'   \item{\code{marginals.random}}{marginals for the random parameters of the model.}
#'   \item{\code{size.random}}{size of the random parameters of the model.}
#'   \item{\code{summary.linear.predictor}}{summary statistics of the linear predictors.}
#'   \item{\code{marginals.linear.predictor}}{marginals for the linear predictors.}
#'   \item{\code{summary.fitted.values}}{summary statistics of the fitted values.}
#'   \item{\code{marginals.fitted.values}}{marginals for the fitted values.}
#'   \item{\code{size.linear.predictor}}{size of the linear predictors.}
#'   \item{\code{summary.hyperpar}}{summary statistics for the hyperparameters of
#'   the model. The summary statistics sorted by longitudinal and
#'   survival components are available by applying the \code{summary} function
#'   to the \code{INLAjoint} object. Particularly, this is the raw output of INLA
#'   and therefore the precision of the residual errors and baseline hazard functions
#'   hyperparameters are provided. Similarly, the Cholesky matrix is given for the
#'   random-effects. The summary function can easily return either variance and covariance
#'    or standard deviations and correlations for all these hyperparameters.}
#'   \item{\code{marginals.hyperpar}}{marginals for the hyperparameters of the model.}
#'   \item{\code{internal.summary.hyperpar}}{summary of the internal hyperparameters,
#'   this is similar to the summary of the hyperparameters but here they are provided as used
#'   for the computations (logarithm scale for residual error and baseline risk hyperparameters).}
#'   \item{\code{internal.marginals.hyperpar}}{marginals for the internal hyperparameters of the model.}
#'   \item{\code{misc}}{miscellaneous (as provided in the INLA output).}
#'   \item{\code{dic}}{Deviance Information Criterion.}
#'   \item{\code{mode}}{.}
#'   \item{\code{joint.hyper}}{.}
#'   \item{\code{nhyper}}{.}
#'   \item{\code{version}}{Version of INLA.}
#'   \item{\code{cpu.used}}{Computation time of INLA.}
#'   \item{\code{all.hyper}}{.}
#'   \item{\code{.args}}{.}
#'   \item{\code{call}}{INLA call.}
#'   \item{\code{selection}}{information about parameters for sampling with inla.rjmarginal.}
#'   \item{\code{cureVar}}{informations about cure fraction submodel for mixture cure survival models.}
#'   \item{\code{variant}}{information about variant for Weibull baseline hazards.}
#'   \item{\code{SurvInfo}}{some information about survival submodels (names of event
#'   indicator and event time variables as well as baseline hazard).}
#'   \item{\code{famLongi}}{list of distributions for the longitudinal markers.}
#'   \item{\code{corLong}}{boolean indicating if random effects are correlated accross markers.}
#'   \item{\code{control.link}}{informations about link function (1=default).}
#'   \item{\code{longOutcome}}{name of longitudinal outcomes.}
#'   \item{\code{survOutcome}}{name of survival outcomes.}
#'   \item{\code{assoc}}{vector with names of all association parameters (longi-surv).}
#'   \item{\code{id}}{name of the id variable.}
#'   \item{\code{timeVar}}{name of time variable.}
#'   \item{\code{range}}{information about range of X-axis values for non-linear associations.}
#'   \item{\code{REstruc}}{names of the grouped random effects for the longitudinal markers.}
#'   \item{\code{mat_k}}{contains the list of random effects covariance matrices when they are
#'   fixed as they are not part of the estimated parameters (used for displaying them in summary).}
#'   \item{\code{fixRE}}{list of the size of number of groups of random effects, each element is a
#'   boolean indicating if the random effects of the group is fixed (TRUE) or estimated (FALSE).}
#'   \item{\code{lonFacChar}}{list of factors and character covariates included in the longitudinal
#'   submodels to keep track of modalities (used internally when doing predictions to reconstruct
#'   categorical covariates).}
#'   \item{\code{survFacChar}}{same as lonFacChar but for survival submodels.}
#'   \item{\code{corRE}}{list indicating if groups of random effects are correlated within
#'   longitudinal submodels.}
#'   \item{\code{basRisk}}{list of the baseline risk used for each survival submodel.}
#'   \item{\code{priors_used}}{informations about priors used in the model, internally used
#'   to display priors in plots (with argument priors=TRUE in the call of the plot function).
#'   Note that priors can also be displayed with the function priors.used() applied to an
#'   INLAjoint object.}
#'   \item{\code{dataLong}}{name of the longitudinal dataset.}
#'   \item{\code{dataSurv}}{name of the survival dataset.}
#'   }
"INLAjoint.object" <- NULL
