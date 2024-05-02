#' Power for two sample exponential (or logrank)
#'
#' Computes power for a stratified test comparing two groups with exponential
#' failure time distributions.  The power can also be used for the logrank test
#' under proportional hazards.
#'
#' @details
#' Computes the power for a stratified test in the two sample exponential
#' problem.  The test is actually based on the optimally weighted average of
#' the log hazard ratio estimates from each stratum.  Unlike other sample size
#' functions, any of the input parameters (except alpha2) can take different
#' values for different strata.  For example, the proportion in each stratum
#' can be specified through either prop.strat or by giving difference accrual
#' rates for each stratum, and the follow-up time can be different for
#' different stratum to allow for different strata to have opened or closed at
#' different times.
#'
#' @param acc.rate Accrual rate
#' @param acc.per Accrual period
#' @param add.fu Additional follow-up after the end of accrual before data are
#' analyzed
#' @param l1 Hazard rates in group 1 (one per stratum)
#' @param r1o2 Hazard ratios for group 1 over group 2; either length 1 or
#' length = number of strata
#' @param p1 Proportion in group 1; either length 1 or length = number of
#' strata.  Can be different for different strata.
#' @param prop.strat Proportion of the total sample in each stratum
#' @param alpha2 Two-sided type I error rate
#' @param nullh Null hypothesis hazard ratio (1 over 2)
#'
#' @return A vector of length 4 giving the overall expected standard error of
#' with weighted average log hazard ratio (se), the expected value under the
#' alternative of the weighted average log hazard ratio (me), the expected
#' total number of failure events under the alternative (n.fail), and the power
#' (power).
#'
#' @seealso
#' \code{\link{seqopr}}; \code{\link{powlgrnk}}; \code{\link{correl.power}}
#'
#' @keywords design
#'
#' @examples
#' powexp2(200, 2, 3, log(2) / 4, 1.5, 0.5)
#' powexp2(200, 2, 3, log(2) / 4, 1.33, 0.5)
#' powexp2(200, 2, 3, log(2) / 4, c(1.5, 1), 0.5, c(0.7, 0.3))
#' powexp2(10, 31.5, 13.75, c(log(2) / 6, log(2) / 8), c(1.4, 1.2), 0.5, c(0.9, 0.1))
#' powexp2(c(9, 1), 31.5, 13.75, c(log(2) / 6, log(2) / 8), c(1.4, 1.2), 0.5)
#'
#' @export powexp2

powexp2 <- function(acc.rate, acc.per, add.fu, l1, r1o2, p1,
                    prop.strat = rep(1, length(l1)), alpha2 = 0.05, nullh = 1) {
  n <- acc.rate * acc.per
  zc <- qnorm(1 - alpha2 / 2)
  p.strat <- prop.strat / sum(prop.strat)
  n.strat <- if (length(n) > 1)
    n else n * p.strat / sum(p.strat)
  q <- 1 - p1
  a <- add.fu
  b <- acc.per + add.fu # censoring uniform on (a,b)
  cpx <- function(l) (exp(-a * l) - exp(-b * l)) / (l * acc.per)
  l2 <- l1 / r1o2
  fp1 <- 1 - cpx(l1)
  fp2 <- 1 - cpx(l2)
  v <- 1 / (n.strat * p1 * fp1) + 1 / (n.strat * q * fp2)
  w <- sum(1 / v)
  sa <- sqrt(1 / w)
  ma <- sum(log(r1o2) / v) / w
  # use sa under null and alt, per RGS 1981 -- also, consistent with seqopr
  # and srvpwr
  # sa <- sqrt(1/(n*p1*fp1)+1/(n*q*fp2))
  # zb <- (zc*sa-abs(log(r1o2)))/sa
  zb <- (zc * sa - abs(ma - log(nullh))) / sa
  # zb2 <- (zc*sn2-abs(log(r1o2)))/sa
  c(se = sa, me = ma, n.fail = sum(n.strat * fp1 * p1 + n.strat * fp2 * q),
    power = 1 - pnorm(zb))
}
