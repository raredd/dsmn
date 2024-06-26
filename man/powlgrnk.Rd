% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/powlgrnknew.R
\name{powlgrnk}
\alias{powlgrnk}
\title{Computes the power of the logrank test}
\usage{
powlgrnk(
  acc.per,
  acc.rate,
  add.fu,
  alpha = 0.025,
  p.con = 0.5,
  hazcon = function(x, lc, pc, ...) {
     u <- (1 - pc) * exp(-lc * x)
     lc * u/(pc +
    u)
 },
  survcon = function(x, lc, pc, ...) pc + (1 - pc) * exp(-lc * x),
  haztst = function(x, lt, pt, ...) {
     u <- (1 - pt) * exp(-lt * x)
     lt * u/(pt +
    u)
 },
  survtst = function(x, lt, pt, ...) pt + (1 - pt) * exp(-lt * x),
  control.rate = NULL,
  test.rate = NULL,
  control.cure = 0,
  test.cure = 0,
  ...
)
}
\arguments{
\item{acc.per}{Planned duration of accrual}

\item{acc.rate}{Number of patients expected to be entered per time unit}

\item{add.fu}{Additional follow-up between the end of accrual and the time
of analysis}

\item{alpha}{The one-sided type I error of the test}

\item{p.con}{The proportion randomized to the control arm}

\item{hazcon}{A function evaluating the control group hazard function at a
vector of times}

\item{survcon}{A function evaluating the control group survivor function at
a vector of times}

\item{haztst}{A function evaluating the experimental or test group hazard
function at a vector of times}

\item{survtst}{A function evaluating the experimental or test group survivor
function at a vector of times}

\item{control.rate}{Exponential hazard rate in control group for non-cured}

\item{test.rate}{Exponential hazard rate in test group for non-cured}

\item{control.cure}{Cure fraction in control group}

\item{test.cure}{Cure fraction in test group}

\item{...}{additional arguments passed to survival and hazard functions}
}
\value{
Returns a vector giving the power of the test under the specified
conditions (\code{power}), the total sample size (\code{n}) and the
expected number of events (\code{nd}).
}
\description{
Computes the power of the two-group logrank test for arbitrary failure time
distributions in a standard clinical trials setting.
}
\details{
The calculations assume that \code{n=acc.rate*acc.per} patients are entered
uniformly over the period \code{[0,acc.per]}, with follow-up continuing for
an addition \code{add.fu} time units, so censoring will be uniform on
\code{[add.fu,add.fu+acc.per]}.  The failure probability under the specified
failure distribution is then computed, as is the expected number of
failures, \code{n*fail.prob}.

\code{hazcon} must be the hazard function corresponding to the survivor
function in \code{survcon} (and similarly for \code{haztst} and
\code{survtst}).  The program does not check for consistency of these
functions, though.  One way to check would be to compare \code{survcon(x)}
to \code{exp(-integrate(hazcon,0,x)$value)} for various values \code{x}.
The default functions are for the exponential cure rate model, which reduces
to the two-sample exponential when the cure fractions are 0.

If custom \code{hazcon}, \code{survcon}, \code{haztst}, and \code{survtst}
functions are given, the arguments should match those in the default
functions (ie \code{x,lc,pc,...} for \code{hazcon}, etc).  The values of
\code{lc}, \code{pc}, \code{lt}, and \code{pt} are specified in the
\code{control.rate}, \code{control.cure}, \code{test.rate}, and
\code{test.cure} arguments and passed through to these functions.

The calculations are performed by using the \code{integrate} function to
approximate the expectations of the logrank score and the logrank variance
estimator and the variance of the logrank score under the specified
conditions, and then using a normal approximation to the distribution of the
logrank statistic.

This function only computes power for a single analysis, and does not
consider group sequential tests.

Caution: consistent time units must be used for all quantities; eg if the
accrual rate is given in patients/month, then the hazard rates must be in
units of events/month and \code{acc.per} and \code{add.fu} must be in
months.
}
\examples{
# Exponential distributions
powlgrnk(5, 200, 3, control.rate = 0.1, test.rate = 0.075)

# Cure rate
powlgrnk(3, 200, 3, control.rate = log(2) / 3, test.rate = log(2) / 4,
         control.cure = 0.3, test.cure = 0.4)

# Exponential cure-rate with proportional hazards alternative
ht <- function(x, rat, ...) rat * hc(x, ...)
st <- function(x, rat, ...) sc(x, ...) ^ rat
sc <- function(x, pic, lic, ...) pic + (1 - pic) * exp(-lic * x)
hc <- function(x, pic, lic, ...) {
  u <- (1 - pic) * exp(-lic * x); lic * u / (pic + u)
}
powlgrnk(5, 200, 3, hazcon = hc, survcon = sc, haztst = ht,
         survtst = st, lic = log(2) / 3, pic = 0.3, rat = 0.5)

}
\seealso{
\code{\link{powlgrnk6}}
}
\keyword{design}
\keyword{survival}
