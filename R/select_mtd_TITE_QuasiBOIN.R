
#' Obtain the maximum tolerated dose (MTD) of BOIN/gBOIN/TITEBOIN/TITEgBOIN designs
#'
#'
#' @param target the target toxicity probability (example: target <- 0.30) or the target normalized ETS (example: target <- 0.47 / 1.5)
#' @param ntox   Number of patients with DLT or the sum of Normalized ETS (equivalent toxicity score)
#' @param npts   Number of patients treated at each dose level
#' @param Neli   the sample size cutoff for elimination. The default is Neli=3.
#' @param cutoff.eli the cutoff to eliminate an overly toxic dose for safety.
#'                  We recommend the default value of (\code{cutoff.eli=0.95}) for general use.
#' @param extrasafe set \code{extrasafe=TRUE} to impose a more stringent stopping rule
#' @param offset a small positive number (between 0 and 0.5) to control how strict the
#'               stopping rule is when \code{extrasafe=TRUE}. A larger value leads to a more
#'               strict stopping rule. The default value \code{offset=0.05} generally works well.
#' @param print print the additional result or not. The default value is print=FALSE
#' @param gdesign For BOIN and TITEBOIN, "FALSE" should be assigned.
#'        For gBOIN and TITEgBOIN, "TRUE" should be assigned . The default is \code{gdesign=FALSE}.

################################################################################################################################
# Isotonic transformation will transform non-monotonic points to monotomic points.                                             #
# For example, ntox=c(1.38, 1.38, 1.05, 2.72) and npts=c(9.1, 6.1, 9.1, 6.1), then p=ntox/npts=c(0.152, 0.227, 0.115, 0.445)   #
# After Isotonic transformation, it will be 0.149, 0.149, 0.149, 0.445, no decreasing, only increasing                         #
################################################################################################################################


select_mtd_TITE_QuasiBOIN <- function(target,ntox, npts, Neli=3, cutoff.eli = 0.95, extrasafe = FALSE, offset = 0.05, print = FALSE,gdesign=FALSE) {
  ## isotonic transformation using the pool adjacent violator algorithm (PAVA)
  pava <- function(x, wt = rep(1, length(x))) {
    n <- length(x)
    if (n <= 1)
      return(x)
    if (any(is.na(x)) || any(is.na(wt))) {
      stop("Missing values in 'x' or 'wt' not allowed")
    }
    lvlsets <- (1:n)
    repeat {
      viol <- (as.vector(diff(x)) < 0)
      if (!(any(viol)))
        break
      i <- min((1:(n - 1))[viol])
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i + 1]
      ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
      x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
      lvlsets[ilvl] <- lvl1
    }
    x
  }
  ## determine whether the dose has been eliminated during the trial
  y = ntox
  n = npts
  ndose = length(n)
  elimi = rep(0, ndose)
  for (i in 1:ndose) {
    if (n[i] >= Neli) {
      if (1 - pbeta(target, y[i] + 1, n[i] - y[i] + 1) > cutoff.eli) {
        elimi[i:ndose] = 1
        break
      }
    }
  }

  if (extrasafe) {
    if (n[1] >= Neli) {
      if (1 - pbeta(target, y[1] + 1, n[1] - y[1] + 1) > cutoff.eli - offset) {
        elimi[1:ndose] = 1
      }
    }
  }

  ## no dose should be selected (i.e., selectdose=99) if the first dose is already very toxic or all uneliminated doses are never used to treat patients
  if (elimi[1] == 1 || sum(n[elimi == 0]) == 0) {
    selectdose = 99
  }
  else {
    adm.set = (n != 0) & (elimi == 0)
    adm.index = which(adm.set == T)
    y.adm = y[adm.set]
    n.adm = n[adm.set]

    ## poster mean and variance of toxicity probabilities using beta(0.005, 0.005) as the prior
    phat = (y.adm + 0.005)/(n.adm + 0.01)
    phat.var = (y.adm + 0.005) * (n.adm - y.adm + 0.005)/((n.adm + 0.01)^2 * (n.adm + 0.01 + 1))

    ## perform the isotonic transformation using PAVA
    phat = pava(phat, wt = 1/phat.var)
    phat = phat + (1:length(phat)) * 1e-10  ## break ties by adding an increasingly small number
    selectd = sort(abs(phat - target), index.return = T)$ix[1]  ## select dose closest to the target as the MTD
    selectdose = adm.index[selectd]
  }
  if (print == TRUE) {
    if (selectdose == 99) {
      message("All tested doses are overly toxic. No MTD is selected! \n")
    }
    else {
      message("The MTD is dose level ", selectdose, "\n\n")
    }
    trtd = (n != 0)
    poverdose = pava(1 - pbeta(target, y[trtd] + 0.05, n[trtd] -
                                 y[trtd] + 0.05))
    phat.all = pava((y[trtd] + 0.05)/(n[trtd] + 0.1), wt = 1/((y[trtd] +
                                                                 0.05) * (n[trtd] - y[trtd] + 0.05)/((n[trtd] + 0.1)^2 *
                                                                                                       (n[trtd] + 0.1 + 1))))

    if(gdesign==TRUE){

      message("Dose    Posterior normalized ETS estimate        95%                  \n",
          sep = "")
      message("Level     Estimate                       Credible Interval   Pr(toxicity>",
          target, "|data)\n", sep = "")
      for (i in 1:ndose) {
        if (n[i] > 0) {
          message(" ", i, "        ", formatC(phat.all[i],
                                          digits = 2, format = "f"), "                            (", formatC(qbeta(0.025,
                                                                                                               y[i] + 0.05, n[i] - y[i] + 0.05), digits = 2,
                                                                                                         format = "f"), ", ", formatC(qbeta(0.975, y[i] +
                                                                                                                                              0.05, n[i] - y[i] + 0.05), digits = 2, format = "f"),
              ")            ", formatC(poverdose[i], digits = 2,
                                       format = "f"), "\n")
        }

        else {
          message(" ", i, "        ", "----", "                           (",
              "------------", ")            ", "----", "\n")
        }
      }

    }else{
      message("Dose    Posterior DLT rate        95%                  \n",
          sep = "")
      message("Level     Estimate          Credible Interval   Pr(toxicity>",
          target, "|data)\n", sep = "")
      for (i in 1:ndose) {
        if (n[i] > 0) {
          message(" ", i, "        ", formatC(phat.all[i],
                                          digits = 2, format = "f"), "                (", formatC(qbeta(0.025,
                                                                                                  y[i] + 0.05, n[i] - y[i] + 0.05), digits = 2,
                                                                                            format = "f"), ", ", formatC(qbeta(0.975, y[i] +
                                                                                                                                 0.05, n[i] - y[i] + 0.05), digits = 2, format = "f"),
              ")            ", formatC(poverdose[i], digits = 2,
                                       format = "f"), "\n")
        }

        else {
          message(" ", i, "        ", "----", "               (",
              "------------", ")            ", "----", "\n")
        }
      }
    }
    message("NOTE: no estimate is provided for the doses at which no patient was treated.")
  }
  else {
    return(selectdose)
  }
}








