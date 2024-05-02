#' Operating Characteristics of Two-Stage Phase II Designs
#'
#' Calculates the operating characteristics of single arm, two-stage phase II
#' designs.
#'
#' @details
#' In the first stage, \code{n1} cases are accrued and the study is stopped if
#' \code{r1} or fewer responses are observed.  If more than \code{r1} responses
#' are observed then an additional \code{n2} cases are accrued.  The drug or
#' regimen will be declared ineffective if \code{r2} or fewer total responses
#' are observed, and worth further study if more than \code{r2} responses are
#' observed.  Subjects are assumed to be independent with a common response
#' probability \code{p1}.
#'
#' @param n1 Number of cases accrued in the first stage
#' @param n2 Number of additional cases accrued in the second stage
#' @param p1 Response probability
#' @param r1 max number responses in first stage where drug would still be
#' declared ineffective
#' @param r2 max number of total responses for drug to be declared ineffective
#'
#' @return \code{twostg} returns a list of class \code{twostg}, with components
#' \item{inputs }{a vector containing the input values} \item{prob.inactive}{a
#' vector giving the total probability of \code{r2} or fewer responses and the
#' probability of \code{r1} or fewer in the first stage}
#'
#' \code{print.twostg} prints a summary of the input data and the probability
#' of being declared ineffective overall and at the first stage.
#'
#' @seealso \code{\link{pickwin}}, \code{\link{rp21}}, \code{\link{simon}},
#' \code{\link{bin1samp}}
#' @keywords design
#'
#' @examples
#' twostg(14, 18, 0.1, 1, 4)
#' twostg(14, 18, 0.3, 1, 4)
#'
#' @export twostg

twostg <- function(n1, n2, p1, r1, r2) {
  # prob calculations for 2-stage phase II design
  # n1=# entered at first stage, n2=additional # entered at second stage
  # p1, response probability
  # r1=max number responses stage 1 trt 1 and still stop
  # r2=max total number responses on trt 1 to be inactive
  if (n1 < 1 | n2 < 1 | r1 < 0 | r2 < 0 | p1 <= 0 | r1 > n1 |
      r2 > n2 + n1 | p1 >= 1)
    stop('invalid arguments')
  x1 <- 0:n1
  x2 <- 0:n2
  w1 <- dbinom(x1, n1, p1)
  w3 <- dbinom(x2, n2, p1)
  # b1=marginal dist of # responses in group 1
  u1 <- c(outer(x1[-(1:(r1 + 1))], x2, '+'))
  u2 <- c(outer(w1[-(1:(r1 + 1))], w3))
  b1 <- c(w1[1:(r1 + 1)], tapply(u2, u1, sum))

  structure(
    list(inputs = c(n1 = n1, n2 = n2, p1 = p1, r1 = r1, r2 = r2),
         prob.inactive = c(total = sum(b1[1:(r2 + 1)]), sum(b1[1:(r1 + 1)]))),
    class = 'twostg'
  )
}

#' Print a summary of the output from twostg
#'
#' Prints a summary of the output from \code{\link{twostg}}
#'
#' @details
#' Prints operating characteristics and stopping rules of \code{\link{twostg}}
#'
#' @param x An object of class \code{twostg} (output of \code{twostg})
#' @param ... included for compatibility with generic.  Any additional
#' arguments are ignored.
#'
#' @return
#' No value is returned - used for side-effect of printing to console.
#'
#' @seealso \code{\link{twostg}}
#'
#' @export

print.twostg <- function(x, ...) {
  cat('RX1 is declared inactive if <=',format(x$inputs[4]),
      'responses are observed in\nthe first',format(x$inputs[1]),
      'cases or <=',format(x$inputs[5]),'responses are',
    'observed in the first',format(x$inputs[1]+x$inputs[2]),'cases\n\n')
  cat('P(RX1 declared inactive) =',format(signif(x$prob.inactive[1],3)),'\n')
  cat('P(stop at first stage) =',format(signif(x$prob.inactive[2],3)),'\n')
  invisible()
}

#' Optimal Two-Stage Single-Arm Designs
#'
#' For studies with binary endpoints, searches for two-stage sampling designs
#' that minimize the expected number of subjects under the null, subject to
#' various constraints.
#'
#' @details
#' For two-stage phase II designs for studies with binary endpoints, searches
#' over possible two-stage sampling designs to find those that minimize the
#' expected \# of subjects, subject to specified constraints. If the only
#' constraints are the type I and type II errors of the tests, then the designs
#' are the optimal designs of Simon (1989).  If a positive value of
#' \code{n1max} is specified, then only designs with \code{<= n1max} subjects
#' in the first stage are considered.  Also, only designs with \code{<= ntmax}
#' total subjects are considered.  Setting \code{ntmax} to a large value (as in
#' the default), effectively allows the search to consider all possible
#' designs.
#'
#' If \code{minimax=TRUE}, then the minimax designs of Simon (1989), which
#' minimize the maximum sample size, are considered.  As there are typically
#' multiple designs with the same minimum max sample size, the program still
#' selects among the designs in this class based on the the expected sample
#' size under the null.
#'
#' Designs which optimize one particular criterion sometimes have other
#' undesirable properties.  By specifying a value of \code{del > 0}, alternate
#' designs meeting the criteria that have expected sample sizes under the null
#' within \code{del} of the optimal design will also be retained.
#'
#' @aliases simon
#'
#' @param p0 Null hypothesis response probability
#' @param pa Alternative hypothesis response probability
#' @param n1max The maximum number of subjects entered during the first stage.
#' Ignored if \code{<= 0}.
#' @param ntmax The maximum total number of subjects.
#' @param alpha Type I error rate
#' @param beta Type II error rate
#' @param del Searches for designs where the expected number of subjects under
#' the null is within \code{del} of the minimum possible value
#' @param minimax If \code{TRUE}, only searches for designs with the total
#' sample size equal to the minimum possible value
#'
#' @return
#' Returns a list with components \item{designs }{A matrix with a row
#' giving a summary of each design meeting the criteria.  The columns are
#' \code{n1}, the number of subjects entered in the first stage; \code{r1}, the
#' cutoff for stopping at the first stage (continue if \# responses \code{>
#' r1}); \code{n2}, the additional number of subjects enrolled in the second
#' stage; \code{r2}, the cutoff for inactivity after the second stage (reject
#' null if \# responses \code{> r2}); \code{Pstop1.H0}, the probability of
#' stopping after the first stage under H0; \code{size}, the actual type I
#' error; \code{type2}, the actual type II error; \code{E.tot.n.H0}, the
#' expected \# subjects under H0.} \item{call }{The call to \code{simon()}.}
#' \item{description }{A text string giving a brief description of the columns
#' in \code{$designs}.}
#'
#' @seealso
#' \code{\link{twostg}}; \code{\link{bin1samp}}; \code{\link{pickwin}};
#' \code{\link{rp21}}
#'
#' @references
#' Simon R (1989). Optimal two-stage designs for phase II clinical
#' trials. \emph{Controlled Clinical Trials} \strong{10}:1-10.
#'
#' @keywords design
#'
#' @examples
#' simon(0.05, 0.20)
#' simon(0.05, 0.20, minimax = TRUE)
#'
#' @export simon

simon <- function(p0, pa, n1max = 0, ntmax = 1e5, alpha = 0.1, beta = 0.1,
                  del = 1, minimax = FALSE) {
  # prob calculations for 2-stage phase II design
  # optimal simon designs (minimize E(n|H0))
  # p0, pa, null & alt response probabilities
  # n1max=max # subject in first stage
  if (alpha > 0.5 | alpha <= 0 | 1 - beta <= alpha | beta <= 0 | p0 <= 0 |
      p0 >= pa | pa >= 1 | n1max > ntmax)
    stop('invalid arguments')
  # determine min sample size first stage
  n1min <- max(ceiling(log(beta) / log(1 - pa)), 2)
  if (n1min > ntmax)
    stop('no valid designs')
  # optimal one-sample design
  u1 <- bin1samp(p0, pa, alpha, beta)
  # include even if n>n1max
  z <- matrix(c(u1[1:2], 0, u1[2L], 1, u1[5:6], u1[1L]), nrow = 1L)
  if (n1min < u1[1L]) {
    n1max <- if (n1max > 0) min(n1max, u1[1L] - 1) else u1[1L] - 1
  } else if (n1min == u1[1L]) {
    n1max <- if (n1max > 0) min(n1max, u1[1L]) else u1[1L]
  } else stop('no valid designs')
  n1max <- min(n1max, ntmax)
  # determine min total sample size (use randomized decision rule on
  # boundary for single sample)
  b <- 0
  n <- u1[1L]
  while (b <= beta) {
    n <- n - 1
    u <- c(1, 1 - pbinom(0:(u1[2L] + 1), n, p0))
    index <- min((1:(u1[2L] + 3))[u <= alpha])
    pi <- (alpha - u[index]) / (u[index - 1] - u[index])
    u <- if (index > 2)
      1-pbinom((index-3):(index-2), n, pa)
    else c(1, 1 - pbinom(0, n, pa))
    b <- 1 - u[2L] - pi * (u[1L] - u[2L])
  }
  if (n > ntmax)
    stop('no valid designs')
  e0 <- u1[1L]
  for (i in n1min:n1max) { # cases 1st stage
    # stage I
    # feasible stopping rules
    x1 <- 0:i
    w1 <- dbinom(x1, i, p0)
    w2 <- dbinom(x1, i, pa)
    sub <- cumsum(w2) <= beta
    y1 <- x1[sub] # feasible 1st stage stopping rules
    for (r1 in y1) {
      j <- n - i # additional cases at 2nd stage
      if (j > 0) {
        q <- 0
        pi0 <- 1 - sum(w1[1:(r1 + 1)])
        while (q == 0) {
          x2 <- 0:j
          w3 <- dbinom(x2, j, p0)
          w4 <- dbinom(x2, j, pa)
          # b0,ba=marginal dist of # responses H0, Ha
          u3 <- c(outer(x1[-(1:(r1 + 1))], x2, '+'))
          u4 <- c(outer(w1[-(1:(r1 + 1))], w3))
          b0 <- cumsum(c(w1[1:(r1 + 1)], tapply(u4, u3, sum)))
          sub <- b0 < 1 - alpha
          r2 <- sum(sub)
          u4 <- c(outer(w2[-(1:(r1 + 1))], w4))
          ba <- cumsum(c(w2[1:(r1 + 1)], tapply(u4, u3, sum)))
          e1 <- i + pi0 * j
          if (ba[r2 + 1] <= beta) {
            q <- 1
            if (minimax) {
              if (i + j < ntmax) {
                ntmax <- i + j
                e0 <- e1 # need to reset to min E.H0 at min # cases
              } else if (i + j == ntmax)
                e0 <- min(e0, e1)
            } else {
              e0 <- min(e0, e1)
            }
            z <- rbind(z, c(i, r1, j, r2, b0[r1 + 1],
                            1 - b0[r2 + 1], ba[r2 + 1], e1))
          } else {
            j <- j + 1
            if (minimax) {
              if (i + j > ntmax)
                q <- 1
            } else {
              if (e1 > e0 + del | i + j > ntmax)
                q <- 1
            }
          }
        }
      }
    }
  }
  dimnames(z) <- list(NULL, c('n1', 'r1', 'n2', 'r2', 'Pstop1.H0',
                              'size', 'type2', 'E.tot.n.H0'))
  z <- z[z[, 1L] + z[, 3L] <= ntmax, , drop = FALSE]
  z <- z[z[, 8L] <= e0 + del, , drop = FALSE]
  list(
    designs = z[order(z[, 8L]), , drop = FALSE],
    call = match.call(),
    description = c('n1, n2 = cases 1st stage and additional # in 2nd',
                    'r1, r2 = max # responses 1st stage and total to declare trt inactive')
  )
}
