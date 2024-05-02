#' Compute Approximate Conditional Power
#'
#' Computes an approximation to the conditional power of a sequential test
#' based on the asymptotic normality of the statistic.
#'
#' @details
#' Computes an approximation to the conditional probability of the test being
#' rejected at the final analysis (i.e., that the test statistic will be >
#' \code{za}), conditional on the current data, if the alternative for which
#' the sequential test has power \code{power} is true.  The approximation is
#' based on asymptotic joint normality of the test statistic at the current and
#' final analysis, and on the increment in the information being independent of
#' the current data.
#'
#' @param za standard normal critical value at final analysis
#' @param power the overall planned power of the sequential test
#' @param zs value of the test statistic on the standard normal scale at the
#' current interim analysis
#' @param infr information fraction (proportion of total information) at the
#' current interim analysis
#'
#' @return The approximate conditional power (a scalar).
#'
#' @seealso
#' \code{\link{condpow}}; \code{\link{condpowcure}}
#'
#' @keywords design
#'
#' @examples
#' acondpow(2.0, 0.85, 0.32, 0.41)
#'
#' @export

acondpow <- function(za, power, zs, infr) {
  # za=normal critical value at final analysis
  # zs=current value of logrank test stat (std normal scale)
  # infr = current information fraction
  # power = power (<1)
  1 - pnorm((za - zs * sqrt(infr)) /
              sqrt(1 - infr) - (za + qnorm(power)) * sqrt(1 - infr))
}
