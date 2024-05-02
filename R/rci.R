#' Repeated Confidence Interval
#'
#' Computes a repeated confidence interval on the hazard ratio for a group
#' sequential study with a failure-time endpoint.
#'
#' @details
#' Computes a repeated confidence interval on the hazard ratio.  The ratio is
#' of experimental/control, so values < 1 indicate that the experimental
#' treatment is better (has lower hazard).  The hazard ratio is estimated from
#' the Cox proportional hazards model.  Note that the RCI is on the hazard
#' ratio scale, not the log scale.
#'
#' The user can either specify the type of use function boundary \code{use}
#' (see \code{\link{sequse}} for definitions), the information times of the
#' previous and current interim analyses \code{inf}, and the confidence level
#' \code{conf}, or can specify the boundary critical value on the standard
#' normal scale at the current information time \code{crit}.
#'
#' @param time Vector of failure/censoring times
#' @param status Vector of failure indicators (1=failure, 0=censored)
#' @param rx Treatment variable
#' @param inf Information times of the previous and current interim analyses
#' (not needed if \code{crit} specified)
#' @param control The code of \code{rx} that corresponds to the control group
#' @param strat If given, a stratified proportional hazards model is fit, with
#' strata defined by the unique values of \code{strat}
#' @param conf The two-sided confidence coefficient (not needed if \code{crit}
#' specified)
#' @param use The type of use-function boundary (not needed if \code{crit}
#' specified)
#' @param crit Standard normal deviate to use as a multiplier in the RCI.  Not
#' needed if \code{inf} and \code{use} are specified
#'
#' @return Returns a vector of length 4 giving the lower confidence limit on
#' the hazard ratio, the estimated ratio, the upper confidence limit, and the
#' critical value on the standard normal scale at the current analysis.
#'
#' @seealso
#' \code{\link{sequse}}
#'
#' @references Jennison and Turnbull (1990). \emph{Statistical Science}
#' \strong{5}:299-317.
#'
#' @keywords survival
#'
#' @examples
#' set.seed(3)
#' ft <- c(rexp(100), rexp(100) / 0.67)
#' ct <- runif(200) * 3
#' fi <- ifelse(ct < ft, 0, 1)
#' ft <- pmin(ft, ct)
#' rx <- c(rep(0, 100), rep(1, 100))
#' rci(ft, fi, rx, inf = c(0.25, 0.54))
#' rci(ft, fi, rx, crit = 2.4)
#'
#' @export rci

rci <- function(time, status, rx, inf, control = 0, strat, conf = 0.95,
                use = 1, crit = NULL) {
  if (length(time) != length(status) | length(time) != length(rx))
    stop('lengths incompatible')
  if (is.null(crit)) {
    if (min(inf) <= 0 | max(inf) > 1)
      stop('inf out of range')
    if (conf <= 0 | conf >= 1)
      stop('conf must be in (0,1)')
    alpha <- (1 - conf) / 2
    crit <- rev(sequse(inf, alpha, use))[1L]
  }
  # rx <- as.numeric(rx == control)
  rx <- as.numeric(rx != control)
  if (missing(strat)) {
    z <- do.call('coxph', list(formula = Surv(time, status) ~ rx,
                               data = data.frame(time, status, rx)))
  } else {
    if (length(strat) != length(time))
      stop('lengths incompatible')
    z <- do.call('coxph', list(formula = Surv(time, status) ~ rx + strata(strat),
                               data = data.frame(time, status, rx, strat)))
  }
  z <- c(exp(z$coef + c(-1, 0, 1) * crit * c(sqrt(z$var))), crit)
  # names(z) <- c('lcl', 'ehr(C/E)', 'ucl', 'crit')
  names(z) <- c('lcl', 'ehr(E/C)', 'ucl', 'crit')
  z
}
