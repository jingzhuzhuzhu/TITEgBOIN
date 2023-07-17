\name{next_TITE_QuasiBOIN}
\alias{next_TITE_QuasiBOIN}
\title{
Determine the dose for the next cohort of new patients for single-agent trials using Bayesian optimal interval (BOIN) (Yuan et al. 2016)/Generalized Bayesian optimal interval (gBOIN) (Mu et al. 2019)/Time-to-event bayesian optimal interval (TITEBOIN) (Lin et al. 2020)/Time-to-event generalized bayesian optimal interval (TITEgBOIN) (Takeda et al. 2022) designs.
}
\description{
Determine the dose for the next cohort of new patients for single-agent trials using Bayesian optimal interval (BOIN)/Generalized Bayesian optimal interval (gBOIN)/
Time-to-event bayesian optimal interval (TITEBOIN)/Time-to-event generalized bayesian optimal interval (TITEgBOIN) designs.
}
\usage{
next_TITE_QuasiBOIN(target,n,npend, y, ft, d, maxt=28, p.saf = 0.6 * target,
                                p.tox = 1.4 * target,elimination=NA,cutoff.eli = 0.95,
                                extrasafe = FALSE, offset=0.05,n.earlystop = 100,
                                maxpen=0.50,Neli=3,print_d = FALSE,gdesign=FALSE)

}

\arguments{
  \item{target}{
the target toxicity probability (example: target <- 0.30) or the target normalized equivalent toxicity score (ETS) (example: target <- 0.47 / 1.5).
}
  \item{n}{
number of patients treated at each dose level.
}
  \item{npend}{
for Time-to-event bayesian optimal interval (TITEBOIN)/Time-to-event generalized bayesian optimal interval (TITEgBOIN), the number of pending patients at each dose level.
for Bayesian optimal interval (BOIN)/Generalized Bayesian optimal interval (gBOIN), "NA" should be assigned.
}
  \item{y}{
number of patients with dose limiting toxicity (DLT) or the sum of Normalized equivalent toxicity score (ETS).
}
  \item{ft}{
for Time-to-event bayesian optimal interval (TITEBOIN)/Time-to-event generalized bayesian optimal interval (TITEgBOIN), Total follow-up time for pending patients for toxicity at each dose level (days). for Bayesian optimal interval (BOIN)/Generalized Bayesian optimal interval (gBOIN), "NA" should be assigned.
}
  \item{d}{
current dose level.
}
  \item{maxt}{
for Time-to-event bayesian optimal interval (TITEBOIN)/Time-to-event generalized bayesian optimal interval (TITEgBOIN), Length of assessment window for toxicity (days).
for Bayesian optimal interval (BOIN)/Generalized Bayesian optimal interval (gBOIN), "NA" should be assigned.
}
  \item{p.saf}{
the lower bound. The default value is p.saf=0.6*target.
}
  \item{p.tox}{
the upper bound. The default value is p.tox=1.4*target.
}
  \item{elimination}{
elimination of each dose (0,1 should be assigned, 0 means the dose is not eliminated,
1 means the dose is eliminated due to over toxic(elimination=NA, 0 is defauted for each dose level)).
}
  \item{cutoff.eli}{
the cutoff to eliminate an overly toxic dose for safety.
We recommend the default value of (cutoff.eli=0.95) for general use.
}
  \item{extrasafe}{
set extrasafe=TRUE to impose a more stringent stopping rule.
}
  \item{offset}{
a small positive number (between 0 and 0.5) to control how strict the
stopping rule is when extrasafe=TRUE. A larger value leads to a more
strict stopping rule. The default value offset=0.05 generally works well.
}
  \item{n.earlystop}{
the early stopping parameter. The default value is n.earlystop=100.
}
  \item{maxpen}{
for Time-to-event bayesian optimal interval (TITEBOIN)/Time-to-event generalized bayesian optimal interval (TITEgBOIN), the upper limit of the ratio of pending patients.
for Bayesian optimal interval (BOIN)/Generalized Bayesian optimal interval (gBOIN), "NA" should be assigned.
}
  \item{print_d}{
print the additional result or not. The default value is print=FALSE.
}
  \item{Neli}{
the sample size cutoff for elimination. The default is Neli=3.
}
  \item{gdesign}{
for Bayesian optimal interval (BOIN) and Time-to-event bayesian optimal interval (TITEBOIN), "FALSE" should be assigned.
for Generalized Bayesian optimal interval (gBOIN) and Time-to-event generalized bayesian optimal interval (TITEgBOIN), "TRUE" should be assigned . The default is gdesign=FALSE.
}
}

\value{
next_TITE_QuasiBOIN() returns the toxicity probability and the recommended dose level for the next cohort
including: (1) the lower Bayesian optimal boundary (lambda_e)
(2) the upper Bayesian optimal boundary (lambda_d)
(3) The number of patients or the effective sampe size (ESS) at each dose level (ESS)
(4) The dose limiting toxicity (DLT) rate or mu (the estimated quasi-Bernoulli toxicity probability) at each dose level (mu)
(5) the recommended dose level for the next cohort as a numeric value under (d)
}
\references{
1. Liu S. and Yuan, Y. (2015). Bayesian optimal interval designs for phase I clinical trials, Journal of the Royal Statistical Society: Series C , 64, 507-523.

2. Yuan, Y., Hess, K. R., Hilsenbeck, S. G., & Gilbert, M. R. (2016). Bayesian optimal interval design: a simple and well-performing design for phase I oncology trials. Clinical Cancer Research, 22(17), 4291-4301.

3. Zhou, H., Yuan, Y., & Nie, L. (2018). Accuracy, safety, and reliability of novel phase I trial designs. Clinical Cancer Research, 24(18), 4357-4364.

4. Zhou, Y., Lin, R., Kuo, Y. W., Lee, J. J., & Yuan, Y. (2021). BOIN Suite: A Software Platform to Design and Implement Novel Early-Phase Clinical Trials. JCO Clinical Cancer Informatics, 5, 91-101.

5. Takeda K, Xia Q, Liu S, Rong A. TITE-gBOIN: Time-to-event Bayesian optimal interval design to accelerate dose-finding accounting for toxicity grades. Pharm Stat. 2022 Mar;21(2):496-506. doi: 10.1002/pst.2182. Epub 2021 Dec 3. PMID: 34862715.

6. Yuan, Y., Lin, R., Li, D., Nie, L. and Warren, K.E. (2018). Time-to-event Bayesian Optimal Interval Design to Accelerate Phase I Trials. Clinical Cancer Research, 24(20): 4921-4930.

7. Rongji Mu, Ying Yuan, Jin Xu, Sumithra J. Mandrekar, Jun Yin, gBOIN: A Unified Model-Assisted Phase I Trial Design Accounting for Toxicity Grades, and Binary or Continuous End Points, Journal of the Royal Statistical Society Series C: Applied Statistics, Volume 68, Issue 2, February 2019, Pages 289–308, https://doi.org/10.1111/rssc.12263.

8. Lin R, Yuan Y. Time-to-event model-assisted designs for dose-finding trials with delayed toxicity. Biostatistics. 2020 Oct 1;21(4):807-824. doi: 10.1093/biostatistics/kxz007. PMID: 30984972; PMCID: PMC8559898.

9. Hsu C, Pan H, Mu R (2022). _UnifiedDoseFinding: Dose-Finding Methods for Non-Binary Outcomes_. R package version 0.1.9, <https://CRAN.R-project.org/package=UnifiedDoseFinding>.

}
\author{
Jing Zhu, Jun Zhang, Kentato Takeda
}

\examples{

#For Bayesian optimal interval (BOIN) design
target<-0.3
next_TITE_QuasiBOIN(target=target,n=c(3,3,4,4,4,0),npend=NA, y=c(0,0,1,1,1,0), ft=NA,
                    d=5, maxt=NA,p.saf= 0.6 * target, p.tox = 1.4 * target,elimination=NA,
                    cutoff.eli = 0.95,extrasafe = FALSE, n.earlystop = 10,
                    maxpen=NA,print_d = TRUE,gdesign=FALSE)


#For Generalized Bayesian optimal interval (gBOIN) design
target=0.47/1.5
next_TITE_QuasiBOIN(target=target,n=c(3,3,4,4,4,0),npend=NA,
                    y=c(0, 0, 0.5/1.5, 1.0/1.5, 1.5/1.5, 0),ft=NA, d=5, maxt=NA,
                    p.saf= 0.6 * target, p.tox = 1.4 * target,elimination=NA,
                    cutoff.eli = 0.95,extrasafe = FALSE, n.earlystop = 10,
                    maxpen=NA,print_d = TRUE,gdesign=TRUE)


#For Time-to-event bayesian optimal interval (TITEBOIN) design
target=0.3
next_TITE_QuasiBOIN(target=target,n=c(3,3,4,4,4,0),npend=c(0,0,0,1,2,0), y=c(0,0,1,1,1,0),
                    ft=c(0, 0, 0, 14, 28, 0),d=5, maxt=28,p.saf= 0.6 * target,
                    p.tox = 1.4 * target,elimination=NA,cutoff.eli = 0.95,
                    extrasafe = FALSE, n.earlystop = 10,maxpen=0.5,print_d = TRUE,
                    gdesign=FALSE)


#For Time-to-event generalized bayesian optimal interval (TITEgBOIN) design
target=0.47/1.5
next_TITE_QuasiBOIN(target=target,n=c(3,3,4,4,4,0),npend=c(0,0,0,1,2,0),
                    y=c(0, 0, 0.5/1.5, 1.0/1.5, 1.5/1.5, 0),ft=c(0, 0, 0, 14, 28, 0),
                    d=5, maxt=28,p.saf= 0.6 * target, p.tox = 1.4 * target,
                    elimination=NA,cutoff.eli = 0.95,extrasafe = FALSE,
                    n.earlystop = 10,maxpen=0.5,print_d = TRUE,gdesign=TRUE)

}
