#' Sample Size Calculation for a One-Sample Exponential
#'
#' Computes the sample size for a given size and power for the test for a
#' specified improvement in the exponential rate parameter.
#'
#' @details
#' Accrual is assumed to be uniform over the accrual period, so censoring is
#' uniform on the interval \code{add.fu} to \code{acc.per + add.fu}. The power
#' is based on the Wald test for the log failure rate parameter unless
#' \code{nonpar=TRUE}, in which case the test is based on the nonparametric
#' estimate of the cumulative hazard at \code{tp}. Simulations suggest that the
#' distribution approximation to the latter test is anti-conservative, so the
#' test (and the corresponding power calculations should be used with extreme
#' caution). Exponential failures are assumed for the calculations, even for
#' the nonparametric test.
#'
#' Any combination of parameters determining the null and alternative event
#' rates can be specified.
#'
#' @param control.rate The hazard rate under the null hypothesis
#' @param pct.imp The percent improvement in the median failure time under the
#' alternative
#' @param acc.per The length of the time period during which patients are
#' accrued
#' @param add.fu The length of time after the end of accrual during which
#' patients continue to be followed prior to analysis
#' @param alpha The one-sided type I error rate
#' @param power The power of the test
#' @param control.med The median under the null
#' @param control.survp The survival probability at \code{tp} under the null
#' @param tp Time point for survival probability values and/or the
#' nonparametric test
#' @param pct.reduc Percent reduction in the hazard rate under the alternative
#' @param nonpar If TRUE, the test based on the nonparametric cumulative hazard
#' estimator at \code{tp} is used
#' @param alt.rate The hazard rate under the alternative
#' @param alt.med The median under the alternative
#' @param alt.survp The survival probability at \code{tp} under the alternative
#'
#' @return
#' A vector giving the sample size (n), the number of failures needed
#' (nd), the proportion of failures (fail.prob), the hazard ratio of the null
#' to alternative hazard rates (haz.ratio), the input values of \code{pct.imp}
#' and \code{control.rate}, the hazard rate under the alternative (alt.rate)
#' and the input values of \code{alpha} and \code{beta}.
#'
#' @keywords design survival
#'
#' @examples
#' surv1samp(-1.333 * log(0.96) / 5, 33.3, 4, 3, alpha = 0.1, power = 0.95)
#' surv1samp(-log(0.95) / 5, 100 * (log(0.95) / log(0.96) - 1), 4, 3,
#'           alpha = 0.1, power = 0.95)
#' surv1samp(-log(0.955) / 5, 100 * (log(0.955) / log(0.965) - 1), 4, 3,
#'           alpha = 0.1, power = 0.95)
#'
#' @export surv1samp

surv1samp <- function(control.rate = NULL, pct.imp = NULL, acc.per, add.fu,
                      alpha = 0.1, power = 0.9, control.med = NULL,
                      control.survp = NULL, tp = NULL, pct.reduc = NULL,
                      nonpar = FALSE, alt.rate = NULL, alt.med = NULL,
                      alt.survp = NULL) {
  # determine control rate
  if (!is.null(control.survp) & !is.null(tp)) {
    control.rate <- -log(control.survp) / tp
  } else if (!is.null(control.med)) {
    control.rate <- log(2) / control.med
  } else if (is.null(control.rate))
    stop('One of control.rate, control.med, or control.survp and tp must be given')
  # determine alt.rate
  if (!is.null(pct.imp)) {
    ratio <- 1 + pct.imp / 100
    alt.rate <- control.rate / ratio
  } else if (!is.null(pct.reduc)) {
    ratio <- 1 / (1 - pct.reduc / 100)
    alt.rate <- control.rate / ratio
  } else if (!is.null(alt.med)) {
    alt.rate <- log(2) / alt.med
    ratio <- control.rate / alt.rate
  } else if (!is.null(alt.survp) & !is.null(tp)) {
    alt.rate <- -log(alt.survp) / tp
    ratio <- control.rate / alt.rate
  } else if (!is.null(alt.rate)) {
    ratio <- control.rate / alt.rate
  } else
    stop('One of pct.imp, pct.reduc, alt.rate, alt.med, or alt.survp and tp must be specified')
  if (is.null(pct.imp))
    pct.imp <- 100 * (control.rate / alt.rate - 1)
  # failure prob -- compute under alternative
  pcens <- (exp(-add.fu * alt.rate) - exp(-(add.fu + acc.per) * alt.rate)) /
    (acc.per * alt.rate)
  pfail <- 1 - pcens
  if (nonpar) {
    if (is.null(tp))
      stop('tp must be specified for nonpar test')
    if (is.null(control.survp))
      control.survp <- exp(-tp * control.rate)
    if (is.null(alt.survp))
      alt.survp <- exp(-tp * alt.rate)
    # variance under the alternative
    if (tp <= add.fu) { # test hazard rate
      s1 <- exp(-tp * alt.rate)
      v <- (1 - s1) / s1
    } else {
      if (tp >= add.fu+acc.per)
        stop('tp must be < add.fu+acc.per')
      s1 <- exp(-add.fu * alt.rate)
      f1 <- function(x)
        acc.per * alt.rate * exp(alt.rate * x) / (add.fu + acc.per - x)
      v <- integrate(f1, add.fu, tp)$value + (1 - s1) / s1
    }
    n <- ((qnorm(1 - alpha) + qnorm(power)) /
            log(control.survp / alt.survp)) ^ 2 * v
    nd <- n * pfail
  } else { # exponential
    # number events
    nd <- ((qnorm(1 - alpha) + qnorm(power)) / log(ratio)) ^ 2
    n <- nd / pfail
  }

  c(n = n, nd = nd, fail.prob = pfail, haz.ratio = ratio, pct.imp = pct.imp,
    control.rate = control.rate, alt.rate = alt.rate, alpha = alpha, power = power)
}
