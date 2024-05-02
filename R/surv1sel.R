#' Compute selection probabilities in a randomized phase II study with a
#' failure time endpoint
#'
#' Uses simulations with exponential failure times and uniform accrual to
#' estimate the probability of each arm being the best in a randomized phase II
#' study, where best is defined as the lowest hazard rate in a Cox proportional
#' hazards model.
#'
#' @details
#' Assumes that the clinical trial will enroll \code{ng*npg} subjects uniformly
#' over the period \code{(0,acc.per)}, with the analysis performed when
#' \code{nevents} failures have been observed.  A Cox proportional hazards
#' model with indicator variables for the groups is fit and the group estimated
#' to have the lowest hazard rate is selected as the best.  In the function,
#' samples are repeatedly generated from exponential distributions with the
#' specified hazard rates and the proportion of the samples for which each
#' group is selected as the best is calculated.
#'
#' @param ng Number of groups
#' @param npg Number of subjects per group (the same scalar value is assumed
#' for all groups)
#' @param acc.per The planned accrual period
#' @param nevents The total planned number of events (across all arms) at the
#' time of analysis
#' @param haz A vector of length \code{ng} giving the failure hazard rates for
#' the expoential distribution in each group
#' @param nsamp The number of samples to generate
#'
#' @return
#' The vector giving the proportion of times each group is selected as
#' the best (ie has the lowest estimated hazard rate).
#'
#' @seealso
#' \code{\link{surv1samp}}
#'
#' @keywords design survival
#'
#' @examples
#' surv1sel(6, 55, 12, 40 * 6, log(2) / c(7.2, 4.8, 4.8, 4.8, 4.8, 4.8))
#'
#' @export

surv1sel <- function(ng, npg, acc.per, nevents, haz, nsamp = 1000) {
  # ng = number of groups
  # npg = number of patients per group
  # acc.per = number of months of accrual
  # nevents = total number of events at the time ov analysis
  # haz = vector of length ng giving the hazard rates for each of the ng groups
  # nsamp = number of samples generated in the simulation
  # output is a vector of length nsamp giving specifying which group has the
  #   lowest estimated failure hazard rate for each sample
  n <- ng * npg
  rx <- rep((1:ng), rep(npg, ng))
  haz2 <- haz[rx]
  rx <- as.factor(rx)
  out <- rep(0, nsamp)
  for (i in 1:nsamp) {
    entry <- runif(n) * acc.per
    failtime <- rexp(n) / haz2
    u <- entry + failtime
    uc <- quantile(u, probs = nevents / n)
    failind <- ifelse(u <= uc, 1, 0)
    failtime <- ifelse(failind == 1, failtime, uc - entry)
    z <- coxph(Surv(failtime, failind) ~ C(rx, contr.treatment))
    z <- z$coef
    if (min(z) > 0) {
      out[i] <- 1
    } else {
      out[i] <- min((1:length(z))[z==min(z)])+1
    }
  }

  table(out) / nsamp
}
