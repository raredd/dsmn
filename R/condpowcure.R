#' Perform Simulation to Calculate Conditional Power of a Logrank Test
#'
#' Performs a simulation to calculate conditional power of a logrank test in a
#' group sequential experiment with failure times generated from a cure rate
#' model.
#'
#' @details
#' Adds \code{add.acc} cases and generates additional follow-up for subjects
#' already entered but censored in the current data. Failure times are
#' generated from exponential cure rate models, with fractions
#' \code{control.cure} and \code{test.cure} that will never fail and
#' exponential distributions for the failure times of those who will eventually
#' fail.
#'
#' Additional interim analyses are performed when \code{inf.time * total.inf}
#' failures are observed. The critical values used for rejecting
#' the null at these additional interim analyses are specified in
#' \code{crit.value}. The number of additional analyses is the length of
#' \code{inf.time} and \code{crit.value}. The information of the current data
#' should not be included. Full information (\code{inf.time = 1}) must be
#' explicitly included. Information times should be in increasing order.
#'
#' After generating the additional data, the function performs the interim
#' analyses, and determines whether the upper boundary is crossed.  This is
#' repeated \code{nsamp} times, giving an estimate of the probability of
#' eventually rejecting the null under the specified distributions, conditional
#' on the current data. Low conditional power under the target alternative for
#' the study might be justification for stopping early.
#'
#' It is very important that consistent time units be used for all arguments
#' and for the current data.
#'
#' @param time Failure/censoring times in the current data
#' @param status failure indicator for current data (1=failure, 0=censored)
#' @param rx treatment variable in current data
#' @param nsamp number of samples to generate in the simulation
#' @param crit.val critical values to be used at \emph{future} analyses
#' @param control the code of \code{rx} for the control group
#' @param control.cure the cure fraction in the control group
#' @param control.rate the failure rate in the conditional exponential
#' distribution for non-cured subjects in the control group
#' @param test.cure the cure fraction in the experimental treatment group
#' @param test.rate the failure rate in the conditional exponential
#' distribution for non-cured subjects in the experimental treatment group
#' @param inf.time information times of future analyses
#' @param total.inf Number of failures at full information
#' @param add.acc Additional number of cases that will be accrued
#' @param add.acc.period Duration of time over which the \code{add.acc} cases
#' will be accrued
#' @param p.con The proportion randomized to the control group
#'
#' @return
#' Returns the proportion of samples where the upper boundary was
#' crossed, and its simulation standard error. Also prints the value of the
#' logrank statistic for the current data.
#'
#' @seealso
#' \code{\link{condpow}}; \code{\link{sequse}}; \code{\link{acondpow}}
#'
#' @references
#' Jennison and Turnbull (1990). \emph{Statistical Science} \strong{5}:299-317.
#'
#' Betensky (1997). \emph{Biometrics}, \strong{53}:794-806.
#'
#' @keywords survival design
#'
#' @import survival
#'
#' @examples
#' ## current data
#' set.seed(3)
#' ft <- c(ifelse(runif(100) < 0.6, rexp(100), 100), ifelse(runif(100) < 0.55,
#'         rexp(100) / 0.95, 100))
#' ct <- runif(200) * 4
#' fi <- ifelse(ct < ft, 0, 1)
#' ft <- pmin(ft, ct)
#' rx <- c(rep(0, 100), rep(1, 100))
#'
#' ## currently at 0.375 information -- assume no prior interim analyses
#' critv <- sequse(c(0.375, 0.7, 1))[-1]
#'
#' condpowcure(ft, fi, rx, nsamp = 10, crit.val = critv, control.cure = 0.4,
#'             control.rate = 1, test.cure = 0.6, test.rate = 0.75,
#'             inf.time = c(0.7, 1), total.inf = 200, add.acc = 200, add.acc.period = 1)
#'
#' \dontrun{
#' ## use larger nsamp in practice:
#' condpowcure(ft, fi, rx, nsamp = 1000, crit.val = critv, control.cure = 0.4,
#'             control.rate = 1, test.cure = 0.6, test.rate = 0.75,
#'             inf.time = c(0.7, 1), total.inf = 200, add.acc = 200, add.acc.period = 1)
#' ## Observed Z = 1.43
#' ## [1] 0.958
#' }
#'
#' @export

condpowcure <- function(time, status, rx, nsamp = 500, crit.val = 1.96, control = 0,
                        control.cure, control.rate, test.cure, test.rate,
                        inf.time = 1, total.inf, add.acc = 0, add.acc.period,
                        p.con = 0.5) {
  # status must be 0=censored, 1=failed
  if (length(time) != length(status) | length(time) != length(rx))
    stop('invalid data')
  if (min(control.cure, test.cure, control.rate, test.rate, total.inf) <= 0)
    stop('invalid rates')
  if (length(crit.val) != length(inf.time))
    stop('crit/inf')
  rx <- as.numeric(rx != control) + 1
  if (sum(status) > 0) {
    z <- do.call('survdiff', list(formula = Surv(time, status) ~ rx,
                                  data = data.frame(time, status, rx)))
    z <- (z$obs[1L] - z$exp[1L]) / sqrt(z$var[1L, 1L])
    cat('Observed Z = ', format(round(z, 2L)), '\n')
  } else {
    cat('No events\n')
  }
  n2 <- nn <- length(time)
  haz <- c(control.rate, test.rate)
  cure <- c(control.cure, test.cure)
  out <- 0
  if (add.acc > 0) {
    n2 <- nn + add.acc
    time <- c(time, rep(0, add.acc))
    status <- c(status, rep(0, add.acc))
    n3 <- round(add.acc * p.con, 0)
    rx <- c(rx, rep(1, n3), rep(2, add.acc - n3))
  }
  # conditional cure probability
  ccf <- cure[rx]
  ccf <- ccf / ((1 - ccf) * exp(-haz[rx] * time) + ccf)
  for (i in 1:nsamp) {
  # generate data
    enter <- rep(0, nn)
    if (add.acc > 0)
      enter <- c(enter, runif(add.acc) * add.acc.period)
    # generate increment from conditional cure rate model
    ctime <- ifelse(status == 1, 0, # already failure
                    ifelse(runif(n2) <= ccf, 1000, # cure
                           rexp(n2) / haz[rx])) + enter # new failure time
    # chronological times of future analyses
    at <- quantile(ctime, probs = c(inf.time * total.inf / n2))
    for (j in 1:length(crit.val)) {
      # increments censored at at[j] on chronological scale--then shift back to
      # analysis time scale
      newtime <- time + pmin(ctime, at[j]) - pmin(enter, at[j])
      newstat <- ifelse(ctime <= at[j], 1, 0)
      z <- do.call('survdiff', list(formula = Surv(newtime, newstat) ~ rx,
                                    data = data.frame(newtime, newstat, rx)))
      z <- (z$obs[1L] - z$exp[1L]) / sqrt(z$var[1L, 1L])
      if (z >= crit.val[j]) {
        out <- out + 1
        break
      }
    }
  }
  out / nsamp
}
