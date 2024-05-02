' Calculates the expected number of events exponential failure models
#'
#' Given a uniform accrual period, specified accrual rate and follow-up period,
#' and specified parameters for the exponential failure distribution, computes
#' the expected number of failures
#'
#' @details
#' The exponential distribution has survivor function S(t) equal to
#' \code{exp(-l*t)}, where \code{l = haz}.  \code{haz} and \code{med} are
#' related through \code{haz=log(2)/med}.  The calculations assume that
#' \code{n=acc.rate*acc.per} patients are entered uniformly over the period
#' \code{[0,acc.per]}, with follow-up continuing for an addition \code{add.fu}
#' time units, so censoring will be uniform on \code{[add.fu,add.fu+acc.per]}.
#' The failure probability under the specified failure distribution is then
#' computed, as is the expected number of failures, \code{n*fail.prob}.
#'
#' Caution: consistent time units must be used for all quantities.
#'
#' @param acc.rate The accrual rate in patients per time unit
#' @param acc.per Duration of accrual
#' @param add.fu Duration of follow-up from end of accrual to analysis
#' @param haz Exponential hazard rate for failure times.  Not needed if
#' \code{med} is specified.
#' @param med Median failure time.  Not needed if \code{haz} is specified.
#'
#' @return
#' A vector of length 3 giving the expected number of failures
#' \code{nevents}, the failure probability \code{fail.prob} and the total
#' sample size \code{n}.
#'
#' @seealso
#' \code{\link{nevents.cure}}
#'
#' @keywords design survival
#'
#' @examples
#' acc.per <- 2.2
#' add.fu <- 3.2
#' nevents(240,acc.per,add.fu,med=1.7)
#'
#' @export nevents

nevents <- function(acc.rate, acc.per, add.fu, haz = NULL, med = NULL) {
  # acc.rate=accrual per time unit
  # acc.per = accrual period
  # add.fu = amount of additional follow-up after end of accrual
  # haz = hazard rate
  # med = median failure time
  if (is.null(haz))
    haz <- log(2) / med
  n <- acc.rate * acc.per
  pcens <- (exp(-add.fu * haz) - exp(-(add.fu + acc.per) * haz)) / (acc.per * haz)
  pfail <- 1 - pcens
  c(nevents = n * pfail, fail.prob = pfail, n = n)
}

#' Calculates the expected number of events for cure rate failure models
#'
#' Given a uniform accrual period, specified accrual rate and follow-up period,
#' and specified parameters for the exponential cure rate failure model,
#' computes the expected number of failures
#'
#' @details
#' The exponential cure rate model has survivor function S(t) equal to
#' \code{p+(1-p)*exp(-l*t)}, where \code{p = cure.rate} and \code{l = haz}.
#' \code{haz} and \code{med} are related through \code{haz=log(2)/med}.
#' Optionally, a proportional hazards shift from the cure rate model with
#' survivor function given by \code{(p+(1-p)*exp(-l*t))^a}, where \code{a =
#' haz.ratio}, can be used by specifying a value for \code{haz.ratio}.  The
#' calculations assume that \code{n=acc.rate*acc.per} patients are entered
#' uniformly over the period \code{[0,acc.per]}, with follow-up continuing for
#' an addition \code{add.fu} time units, so censoring will be uniform on
#' \code{[add.fu,add.fu+acc.per]}.  The failure probability under the specified
#' failure distribution is then computed, as is the expected number of
#' failures, \code{n*fail.prob}.
#'
#' Caution: consistent time units must be used for all quantities.
#'
#' @param acc.rate The accrual rate in patients per time unit
#' @param acc.per Duration of accrual
#' @param add.fu Duration of follow-up from end of accrual to analysis
#' @param cure.rate Proportion cured
#' @param haz Exponential hazard rate for failure times for those not cured.
#' Not needed if \code{med} is specified.
#' @param med Median failure time for those not cured.  Not needed if
#' \code{haz} is specified.
#' @param haz.ratio Proportional hazards shift from the exponential cure rate
#' model
#'
#' @return
#' A vector of length 3 giving the expected number of failures
#' \code{nevents}, the failure probability \code{fail.prob} and the total
#' sample size \code{n}.
#'
#' @seealso
#' \code{\link{nevents}}
#'
#' @keywords design survival
#'
#' @examples
#' acc.per <- 2.2
#' add.fu <- 3.2
#' nevents.cure(240,acc.per,add.fu,cure.rate=.38,med=1.7)
#' nevents.cure(240,acc.per,add.fu,cure.rate=.38,med=1.7,haz.ratio=.667)
#'
#' @export nevents.cure

nevents.cure <- function(acc.rate, acc.per, add.fu, cure.rate,
                         haz = NULL, med = NULL, haz.ratio = 1) {
  # acc.rate=accrual per time unit
  # acc.per = accrual period
  # add.fu = amount of additional follow-up after end of accrual
  # cure.rate = proportion cured
  # haz = hazard rate
  # med = median failure time
  # haz.ratio=hazard ratio from null exponential cure rate model
  if (is.null(haz))
    haz <- log(2) / med
  n <- acc.rate * acc.per
  if (haz.ratio == 1) {
    pcens <- cure.rate + (1 - cure.rate) *
      (exp(-add.fu * haz) - exp(-(add.fu + acc.per) * haz)) / (acc.per * haz)
  } else {
    f <- function(x) (cure.rate + (1 - cure.rate) *
                        exp(-haz * x)) ^ haz.ratio / acc.per
    pcens <- integrate(f, add.fu, add.fu + acc.per)$value
  }
  pfail <- 1 - pcens
  c(nevents = n * pfail, fail.prob = pfail, n = n)
}
