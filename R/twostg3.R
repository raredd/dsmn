#' Compute operating characteristics of two-stage trinomial-outcome designs
#'
#' For a phase II study with three possible ordered outcomes (e.g., response,
#' stable, progression), computes the probability of not rejecting the null
#' hypothesis for a specified two-stage design.
#'
#' @details
#' Consider a phase II study with an ordered categorical outcome with three
#' categories, which are referred to here as response, stable and progression
#' (based on the best response outcome in oncology).  A treatment might improve
#' outcome by either improving the response rate or by improving the disease
#' stabilization rate (where stabilization includes best response of response
#' or stable).  In the designs considered here, the treatment is declared to be
#' active (or sufficiently active to be worth studying further) if it shows a
#' sufficient level of activity for either response or disease stabilization.
#'
#' More formally, let \code{R1} be the number of patients observed to respond
#' and \code{S1} be the number of patients observed to be stable (but not
#' responses) during the first stage of accrual, and similarly let \code{R2}
#' and \code{S2} be the number of additional cases observed to respond or
#' stabilize from the second stage.  The study is stopped after the first stage
#' if BOTH \code{R1 <= r1} and \code{R1+S1 <= s1}.  The study proceeds to the
#' second stage if either \code{R1 > r1} or \code{R1+S1 > s1}.  The treatment
#' is declared to be inactive if either the study stops after the first stage
#' or if at the end of the second stage BOTH \code{R1+R2 <= r2} and
#' \code{R1+R2+S1+S2 <= s2}.  The treatment is considered sufficiently active
#' for further investigation if either \code{R1+R2 > r2} or \code{R1+R2+S1+S2 >
#' s2}.
#'
#' Note that \code{presp} and \code{pstab} are the probabilities of the
#' corresponding disjoint cells in the trinomial model (that is, \code{pstab}
#' is the probability that a case will be stable but not a response).
#'
#' @param n1 Number of subjects enrolled in the first stage
#' @param n2 Number of additional subjects enrolled in the second stage
#' @param presp The probability of the best (response) category
#' @param pstab The probability of the intermediate (stable only) category
#' @param r1 Max number of responses that can occur at the first stage without
#' proceeding to the second stage
#' @param s1 Max number cases that can have responses or stable outcomes at
#' first stage without proceeding to the second stage
#' @param r2 Max number of responses from stages 1 and 2 combined that can be
#' observed without declaring the treatment to be effective
#' @param s2 Max number of cases from stages 1 and 2 combined that can have
#' responses or stable outcomes without declaring the treatment to be effective
#'
#' @return A vector giving the probability of stopping after the first stage
#' (\code{p.stop.1}) and the overall probability that the treatment is declared
#' inactive (\code{p.inactive}).
#'
#' @seealso
#' \code{\link{twostg}}
#'
#' @keywords design
#'
#' @examples
#' twostg3(20, 40, 0.05, 0.10, 1, 3, 5, 22)
#'
#' @export twostg3

twostg3 <- function(n1, n2, presp, pstab, r1, s1, r2, s2) {
  # n1=# at first stage
  # n2=additiona # at second stage
  # presp = Prob response
  # pstab = Prob stable (excluding response)
  # r1=max # responses at first stage and still stop
  # r2=max # total responses and still reject activity
  # s1=max # stable or response at 1st stage and still stop
  # s2=max # stable or response total and still accept H0
  # i1=# responses first stage; i2=#stable (non-response) 1st stg
  # i3=# add responses 2nd stage; i4=# add stable 2nd stg
  # rejp is the probability of concluding the drug is not active
  # (rejecting the drug), which is 1 - the probability of rejecting the
  # null hypothesis
  if (n1 < 1 | n2 < 1 | r1 < 0 | r2 < 0 | s1 < 0 | s2 < 0 | presp <= 0 |
      pstab<=0 | r1 + s1 > n1 | r2 + s2 > n2 + n1 | presp >= 1| pstab >= 1)
    stop('invalid arguments')
  p2 <- pstab / (1 - presp) # cond prob stable given not a response
  rejp1 <- rejp <- 0
  for (i1 in 0:(min(r2, n1))) {
    u1 <- dbinom(i1, n1, presp)
    for (i2 in 0:min(n1 - i1, s2 - i1)) {
      u2 <- dbinom(i2, n1 - i1, p2)
      if (i1 <= r1 & (i1 + i2) <= s1) {
        rejp1 <- rejp1 + u1 * u2
      } else {
        for (i3 in 0:min(n2, r2 - i1, s2 - i2 - i1)) {
          u3 <- dbinom(i3, n2, presp)
          for (i4 in 0:min(n2 - i3, s2 - i1 - i2 - i3)) {
            u4 <- dbinom(i4, n2 - i3, p2)
            rejp <- rejp + u1 * u2 * u3 * u4
          }
        }
      }
    }
  }
  c(p.stop.1 = rejp1, p.inactive = rejp + rejp1)
}
