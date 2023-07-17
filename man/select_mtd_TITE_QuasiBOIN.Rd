\name{select_mtd_TITE_QuasiBOIN}
\alias{select_mtd_TITE_QuasiBOIN}
\title{
Obtain the maximum tolerated dose (MTD) of Bayesian optimal interval (BOIN)/Generalized Bayesian optimal interval (gBOIN)/Time-to-event bayesian optimal interval (TITEBOIN)/Time-to-event generalized bayesian optimal interval (TITEgBOIN) designs
}
\description{
Obtain the maximum tolerated dose (MTD) of Bayesian optimal interval (BOIN) (Yuan et al. 2016)/Generalized Bayesian optimal interval (gBOIN) (Mu et al. 2019)/Time-to-event bayesian optimal interval (TITEBOIN) (Lin et al. 2020)/Time-to-event generalized bayesian optimal interval (TITEgBOIN) (Takeda et al. 2022) designs
}
\usage{
select_mtd_TITE_QuasiBOIN(target,ntox, npts, Neli=3, cutoff.eli = 0.95,
                                extrasafe = FALSE, offset = 0.05,print = FALSE,
                                gdesign=FALSE)
}
\arguments{
  \item{target}{
the target toxicity probability (example: target <- 0.30) or the target normalized equivalent toxicity score (ETS) (example: target <- 0.47 / 1.5).
}
  \item{ntox}{
number of patients with dose limiting toxicity (DLT) or the sum of normalized equivalent toxicity score (ETS).
}
  \item{npts}{
the number of patients enrolled at each dose level.
}
  \item{Neli}{
the sample size cutoff for elimination. The default is Neli=3.
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
  \item{print}{
print the additional result or not. The default value is print=FALSE.
}
  \item{gdesign}{
for Bayesian optimal interval (BOIN) and Time-to-event bayesian optimal interval (TITEBOIN), "FALSE" should be assigned.
for Generalized Bayesian optimal interval (gBOIN) and Time-to-event generalized bayesian optimal interval (TITEgBOIN), "TRUE" should be assigned . The default is gdesign=FALSE.
}
}

\value{
select_mtd_TITE_QuasiBOIN() returns the selected dose
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

#For Bayesian optimal interval (BOIN) design/Time-to-event bayesian optimal interval (TITEBOIN)
#design
target<-0.3
y<-c(0,0,1,2,3,0)
n<-c(3,3,6,9,9,0)
select_mtd_TITE_QuasiBOIN(target=target,ntox=y,npts=n,print=TRUE,gdesign=FALSE)


#For Generalized Bayesian optimal interval (gBOIN) design/Time-to-event generalized bayesian
#optimal interval (TITEgBOIN) design
target<-0.47/1.5
y<-c(0,0,2/1.5,3.5/1.5,5.5/1.5,0)
n<-c(3,3,6,9,9,0)
select_mtd_TITE_QuasiBOIN(target=target,ntox=y,npts=n,print=TRUE,gdesign=TRUE)

}



