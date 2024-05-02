#' Perform Simulation to Calculate Conditional Power of a Logrank Test
#'
#' Performs a simulation to calculate conditional power of a logrank test in a
#' group sequential experiment with exponential failure times.
#'
#' @details
#' Adds \code{add.acc} cases and generates additional follow-up for subjects
#' already entered but censored in the current data. Failure times are
#' generated from exponential distributions. Additional interim analyses are
#' performed when \code{inf.time * total.inf} failures are observed. The
#' critical values used for rejecting the null at these additional interim
#' analyses are specified in \code{crit.value}. The number of additional
#' analyses is the length of \code{inf.time} and \code{crit.value}. The
#' information of the current data should not be included. Full information
#' (\code{inf.time = 1}) must be explicitly included.  Information times should
#' be in increasing order.
#'
#' After generating the additional data, the function performs the interim
#' analyses, and determines whether the upper boundary is crossed (note that
#' all analyses are performed, to calculate additional quantities, even if the
#' boundary is crossed at an early analysis). This is repeated \code{nsamp}
#' times, giving an estimate of the probability of eventually rejecting the
#' null under the specified distributions, conditional on the current data.
#' Low conditional power under the target alternative for the study might be
#' justification for stopping early.
#'
#' It is very important that consistent time units be used for all arguments
#' and for the current data.
#'
#' If \code{strat} is used and \code{add.acc > 0}, then the new observations
#' are generated in a new stratum.
#'
#' @param time Failure/censoring times in the current data
#' @param status failure indicator for current data (1=failure, 0=censored)
#' @param rx treatment variable in current data
#' @param nsamp number of samples to generate in the simulation
#' @param crit.val critical values to be used at \emph{future} analyses
#' @param control the code of \code{rx} for the control group
#' @param control.rate exponential failure rate for the control group
#' @param test.rate exponential failure rate for the experimental treatment
#' group
#' @param inf.time information times of future analyses
#' @param total.inf Number of failures at full information
#' @param add.acc Additional number of cases that will be accrued
#' @param add.acc.period Duration of time over which the \code{add.acc} cases
#' will be accrued
#' @param p.con The proportion randomized to the control group
#' @param strat Variable defining the strata for a stratified test
#'
#' @return Returns the estimated rejection probability and its simulation
#' standard error.  Also prints the value of the logrank statistic for the
#' current data, several approximations to the conditional power (based on
#' large sample normality and independent increments -- each based on slightly
#' different approximations to the distribution of the logrank statistic), the
#' average number of events in each group at full information, the average
#' increment in the logrank score from the current data to full information,
#' and the average calendar time of full information.
#'
#' @seealso
#' \code{\link{condpowcure}}; \code{\link{sequse}}; \code{\link{acondpow}}
#'
#' @references
#' Jennison and Turnbull (1990). \emph{Statistical Science} \strong{5}:299-317.
#'
#' Betensky (1997). \emph{Biometrics} \strong{53}:794-806.
#'
#' @keywords survival design
#'
#' @import survival
#'
#' @examples
#' ## current data
#' set.seed(3)
#' ft <- c(rexp(100), rexp(100) / 0.95)
#' ct <- runif(200) * 3
#' fi <- ifelse(ct < ft, 0, 1)
#' ft <- pmin(ft, ct)
#' rx <- c(rep(0, 100), rep(1, 100))
#'
#' ## currently at 0.43 information -- assume no prior interim analyses
#' critv <- sequse(c(0.43, 0.7, 1))[-1]
#'
#' condpow(ft, fi, rx, nsamp = 10, crit.val = critv, control.rate = 1, test.rate = 0.75,
#'         inf.time = c(0.7, 1), total.inf = 300, add.acc = 200, add.acc.period = 1)
#'
#' \dontrun{
#' # use larger nsamp in practice, eg
#' condpow(ft, fi, rx, nsamp = 1000, crit.val = critv, control.rate = 1, test.rate = 0.75,
#'         inf.time = c(0.7, 1), total.inf = 300, add.acc = 200, add.acc.period = 1)
#' }
#'
#' @export

condpow <- function(time, status, rx, nsamp = 500, crit.val = 1.96,
                    control = 0, control.rate, test.rate, inf.time = 1, total.inf,
                    add.acc = 0, add.acc.period, p.con = 0.5,
                    strat = rep(1, length(time))) {
  if (length(time) != length(status) | length(time) != length(rx))
    stop('invalid data')
  if (min(control.rate,test.rate,total.inf) <= 0)
    stop('invalid rates')
  if (length(crit.val) != length(inf.time))
    stop('crit/inf')
  rx <- as.numeric(rx != control) + 1
  nd <- sum(status)
  if (nd > 0) {
    z <- do.call('survdiff', list(formula = Surv(time, status) ~ rx + strata(strat),
                                  data = data.frame(time, status, rx, strat)))
    Uf <- sum(as.matrix(z$obs - z$exp)[1L, ])
    print(c(Logrank.score = Uf, Logrank.var = z$var[1L, 1L]))
    z <- Uf / sqrt(z$var[1L, 1L])
    cat('Observed Z = ', format(round(z, 3L)),'\n')
  } else {
    Uf <- 0
    cat('No events\n')
  }
  # approximations
  crv <- crit.val[length(crit.val)]
  rho <- table(rx[status == 0]) + add.acc * c(p.con, 1 - p.con)
  rho <- rho / sum(rho)
  fbar <- 1 - nd / total.inf
  tmp <- sqrt(fbar * rho[1L] * rho[2L] * total.inf)
  b1 <- crv / sqrt(fbar) - Uf / tmp + tmp * log(test.rate / control.rate)
  b2 <- crv / sqrt(fbar) - Uf / tmp + sqrt(fbar * total.inf) *
    (test.rate - control.rate) / (test.rate + control.rate)
  d1 <- rho[1L] / (rho[2L] * test.rate / control.rate + rho[1L])
  b3 <- b1 * sqrt(d1 * (1 - d1) / (rho[1L] * rho[2L]))
  b4 <- (b1 - crv / sqrt(fbar)) *
    sqrt(d1 * (1 - d1) / (rho[1L] * rho[2L])) + crv / sqrt(fbar)
  cat('approx 1: ', format(round(1 - pnorm(b1), 3L)), '\n')
  cat('Freedman: ', format(round(1 - pnorm(b2), 3L)), '\n')
  cat('Exponent: ', format(round(1 - pnorm(b3), 3L)), '\n')
  cat('Exponent2: ', format(round(1 - pnorm(b4), 3L)), '\n')
  n2 <- nn <- length(time)
  haz <- c(control.rate, test.rate)
  out3 <- out <- 0
  if (add.acc > 0) {
    strat <- unclass(as.factor(strat))
    ustrat <- sort(unique(strat))
    n2 <- nn + add.acc
    time <- c(time, rep(0, add.acc))
    status <- c(status, rep(0, add.acc))
    n3 <- round(add.acc * p.con, 0L)
    rx <- c(rx, rep(1, n3), rep(2, add.acc - n3))
    strat <- if (length(ustrat) == 1)
      c(strat, rep(ustrat, add.acc))
    else c(strat, rep(max(ustrat) + 1, add.acc))
  }
  out2 <- rep(0, 3L)
  for (i in 1:nsamp) {
    # generate data
    enter <- rep(0, nn)
    if (add.acc > 0)
      enter <- c(enter, runif(add.acc) * add.acc.period)
    # ctime=calendar time of failures
    ctime <- ifelse(status == 1, 0, rexp(n2) / haz[rx]) + enter
    at <- quantile(ctime, probs = c(inf.time * total.inf / n2))
    out3 <- out3 + max(at) # calendar time scale
    tmp <- 0
    for (j in 1:length(crit.val)) {
      newtime <- time + pmin(ctime, at[j]) - pmin(enter, at[j])
      newstat <- ifelse(ctime <= at[j], 1, 0)
      z <- do.call('survdiff',
                   list(formula = Surv(newtime, newstat) ~ rx + strata(strat),
                        data = data.frame(newtime, newstat, rx, strat)))
      zt <- sum(as.matrix(z$obs - z$exp)[1L, ]) / sqrt(z$var[1L, 1L])
      # (z$obs[1]-z$exp[1])/sqrt(z$var[1,1])
      if (zt >= crit.val[j])
        tmp <- 1
    }
    out <- out + tmp
    out2 <- out2 + c(sum(as.matrix(z$obs - z$exp)[1L, ]) - Uf,
                     apply(as.matrix(z$obs), 1L, sum))
  }
  cat('Ave # failures at full information in control group:',
      format(out2[2L] / nsamp), '\n')
  cat('Ave # failures at full information in experim. group:',
      format(out2[3L] / nsamp), '\n')
  cat('Ave increment in logrank score:', format(out2[1L] / nsamp), '\n')
  cat('Ave calendar time of full information:', format(out3/nsamp), '\n')
  c(Power = out, se = sqrt(out * (nsamp - out) / nsamp)) / nsamp
}
