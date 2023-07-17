\name{get_oc_TITE_QuasiBOIN}
\alias{get_oc_TITE_QuasiBOIN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate operating characteristics for finding the maxinum tolerated dose (MTD) using Bayesian optimal interval (BOIN), Generalized Bayesian optimal interval (gBOIN), Time-to-event bayesian optimal interval (TITEBOIN) and Time-to-event generalized bayesian optimal interval (TITEgBOIN) designs
}
\description{Obtain the operating characteristics of the model-assisted design for single agent trials by simulating trials
using Bayesian optimal interval (BOIN) (Yuan et al. 2016)/Generalized Bayesian optimal interval (gBOIN) (Mu et al. 2019)/Time-to-event bayesian optimal interval (TITEBOIN) (Lin et al. 2020)/Time-to-event generalized bayesian optimal interval (TITEgBOIN) designs(Takeda et al. 2022).
}
\usage{
get_oc_TITE_QuasiBOIN(target, prob, score=c(0,0.5,1.0,1.5), TITE=TRUE,
                            ncohort,cohortsize, maxt=1, accrual=3,maxpen=0.5,
                            alpha1=0.5,alpha2=0.5,n.earlystop = 100,Neli=3,
                            startdose = 1,p.saf = 0.6 * target,
                            p.tox =  1.4 * target,cutoff.eli = 0.95,
                            extrasafe = FALSE, offset = 0.05,ntrial = 1000,
                            seed=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target}{
the target toxicity probability (example: target <- 0.30) or the target normalized equivalent toxicity score (ETS) (example: target <- 0.47 / 1.5).
}
  \item{prob}{
a vector (Bayesian optimal interval (BOIN) or Time-to-event bayesian optimal interval (TITEBOIN) design)/matrix (Generalized Bayesian optimal interval (gBOIN) or Time-to-event generalized bayesian optimal interval (TITEgBOIN) design) containing the true toxicity probabilities of the investigational dose levels.
}
  \item{score}{
for Generalized Bayesian optimal interval (gBOIN)/Time-to-event generalized bayesian optimal interval (TITEgBOIN), a vector containing the relative severity of different toxicity grades
in terms of dose limiting toxicity (DLTs) in the dose-finding procedure. As default, toxicity grades of 0/1,2,3, and 4 are assigned values of 0,0.5,1,1.5.
for Bayesian optimal interval (BOIN)/Time-to-event bayesian optimal interval (TITEBOIN), "NA" should be assigned.
}
  \item{TITE}{
for Time-to-event bayesian optimal interval (TITEBOIN)/Time-to-event generalized bayesian optimal interval (TITEgBOIN), "TRUE" should be assigned.
for Bayesian optimal interval (BOIN)/Generalized Bayesian optimal interval (gBOIN), "FALSE" should be assigned.
}
  \item{ncohort}{
the total number of cohorts.
}
  \item{cohortsize}{
the cohort size.
}
  \item{maxt}{
for Time-to-event bayesian optimal interval (TITEBOIN)/Time-to-event generalized bayesian optimal interval (TITEgBOIN), the maximum follow-up time.
for Bayesian optimal interval (BOIN)/ Generalized Bayesian optimal interval (gBOIN), if you don't need to get 1,the average trial duration needed for the trial,
2, the standard deviation of average trial duration needed for the trial. Then "NA" should be assigned;
if you need to get 1,the average trial duration needed for the trial, 2, the standard deviation of
average trial duration needed for the trial. Then please specify the accrual rate and the maximum follow-up time.
}
  \item{accrual}{
for Time-to-event bayesian optimal interval (TITEBOIN)/Time-to-event generalized bayesian optimal interval (TITEgBOIN), the accrual rate, i.e., the number of patients accrued in 1 unit of time,
for Bayesian optimal interval (BOIN)/ Generalized Bayesian optimal interval (gBOIN), if you don't need to get 1,the average trial duration needed for the trial,
2, the standard deviation of average trial duration needed for the trial. Then "NA" should be assigned;
if you need to get 1,the average trial duration needed for the trial, 2, the standard deviation of
average trial duration needed for the trial, Then please specify the accrual rate and the maximum follow-up time.
}
  \item{maxpen}{
for Time-to-event bayesian optimal interval (TITEBOIN)/Time-to-event generalized bayesian optimal interval (TITEgBOIN), the upper limit of the ratio of pending patients.
for Bayesian optimal interval (BOIN)/Generalized Bayesian optimal interval (gBOIN), "NA" should be assigned.
}
  \item{alpha1}{
for Time-to-event bayesian optimal interval (TITEBOIN)/Time-to-event generalized bayesian optimal interval (TITEgBOIN), a number from (0,1) that assume toxicity outcomes occurred with probability alpha1 in the last fraction of alpha2 of the assessment window. The default is alpha1=0.5.
for Bayesian optimal interval (BOIN)/Generalized Bayesian optimal interval (gBOIN), "NA" should be assigned.
}
  \item{alpha2}{
for Time-to-event bayesian optimal interval (TITEBOIN)/Time-to-event generalized bayesian optimal interval (TITEgBOIN), a number from (0,1) that assume toxicity outcomes occurred with probability alpha1 in the last fraction of alpha2 of the assessment window. The default is alpha2=0.5.
for Bayesian optimal interval (BOIN)/Generalized Bayesian optimal interval (gBOIN), "NA" should be assigned.
}
  \item{n.earlystop}{
the early stopping parameter and the decision is to stay. If the number of patients treated at the current dose
reaches n.earlystop, stop the trial and select the maxinum tolerated dose (MTD) based on the observed data.
The default value n.earlystop=100 essentially turns off this type of early stopping.
}
  \item{Neli}{
the sample size cutoff for elimination. The default is Neli=3.
}
  \item{startdose}{
the starting dose level for the trial.
}
  \item{p.saf}{
the lower bound. The default value is p.saf=0.6*target.
}
  \item{p.tox}{
the upper bound. The default value is p.tox=1.4*target.
}
  \item{cutoff.eli}{
the cutoff to eliminate an overly toxic dose for safety. We recommend the default value of cutoff.eli=0.95 for general use.
}
  \item{extrasafe}{
set extrasafe=TRUE to impose a more stringent stopping rule.
}
  \item{offset}{
a small positive number (between 0 and 0.5) to control how strict the stopping rule is when extrasafe=TRUE. A larger value leads to a more strict stopping rule. The default value offset=0.05 generally works well.
}
  \item{ntrial}{
the total number of trials to be simulated.
}
  \item{seed}{
the seed, The default value is seed = 100
}
}
\details{
This function generates he operating characteristics of the Bayesian optimal interval (BOIN)/ Generalized Bayesian optimal interval (gBOIN)/Time-to-event bayesian optimal interval (TITEBOIN)/ Time-to-event generalized bayesian optimal interval (TITEgBOIN) designs for trials by simulating trials under the prespecified true toxicity probabilities of the investigational doses.
}
\value{
get_oc_TITE_QuasiBOIN() returns the operating characteristics of the Bayesian optimal interval (BOIN)/Generalized Bayesian optimal interval (gBOIN)/Time-to-event bayesian optimal interval (TITEBOIN)/Time-to-event generalized bayesian optimal interval (TITEgBOIN) designs as a data frame,
including: (1) the percentage of trials that the maxinum tolerated dose (MTD) is correctly selected,
(2) the percentage of patients that are correctly allocated to the maxinum tolerated dose (MTD),
(3) the percentage of overdosing selection,
(4) the percentage of overdosing allocation,
(5) selection percentage at each dose level,
(6) the number of patients treated at each dose level,
(7) the percentage of patients treated at each dose level,
(8) the number of toxicities observed at each dose level,
(9) the average number of toxicities,
(10) the average number of patients,
(11) the percentage of early stopping without selecting the maxinum tolerated dose (MTD),
(12) the average trial duration needed for the trial,
(13) the standard deviation of average trial duration needed for the trial,
(14) simulation set up data frame, include the target toxicity probability/the normalized target equivalent toxicity score (ETS); the true target toxicity probability/
     the true normalized equivalent toxicity score (ETS) at each dose level based on prob and score, and lambda_e denotes the lower Bayesian optimal boundary
     and lambda_d denotes the upper Bayesian optimal boundary.
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
\note{
We should avoid setting the values of p.saf and p.tox very close to the target.
This is because the small sample sizes of typical phase I trials prevent us from
differentiating the target toxicity rate from the rates close to it. In addition,
in most clinical applications, the target toxicity rate is often a rough guess,
and finding a dose level with a toxicity rate reasonably close to the target rate
will still be of interest to the investigator. In addition, we recommend setting
the value of priortox relatively small, for example, priortox=target/2 to accelerate
the escalation procedure.
}
\examples{

#For Bayesian optimal interval (BOIN) design and Output trial duration as an operating
#characteristics
get_oc_TITE_QuasiBOIN(target=0.3, score=NA,prob=c(0.25,0.30,0.45,0.49,0.53), TITE=FALSE,
                      ncohort=10, cohortsize=3,startdose=1,maxt=28,accrual=10,
                      maxpen=NA,alpha1=NA,alpha2=NA,cutoff.eli=0.95, ntrial=10,seed=6)


#For Bayesian optimal interval (BOIN) design and not Output trial duration as an operating
#characteristics
get_oc_TITE_QuasiBOIN(target=0.3, score=NA,prob=c(0.25,0.30,0.45,0.49,0.53), TITE=FALSE,
                      ncohort=10, cohortsize=3,startdose=1,maxt=NA,accrual=NA,
                      maxpen=NA,alpha1=NA,alpha2=NA,cutoff.eli=0.95, ntrial=10,seed=6)


#For Generalized Bayesian optimal interval (gBOIN) design and Output trial duration as an
#operating characteristics
target<-0.47/1.5
prob <- matrix(c(0.83,	0.75,	0.62,	0.51,	0.34,	0.19,
                 0.12,	0.15,	0.18,	0.19,	0.16,	0.11,
                 0.04,	0.07,	0.11,	0.14,	0.15,	0.11,
                 0.01,	0.03,	0.09,	0.16,	0.35,	0.59), ncol = 6, byrow = TRUE)
get_oc_TITE_QuasiBOIN(target=target, score=c(0,0.5,1,1.5),prob=prob, TITE=FALSE,ncohort=10,
                      cohortsize=3,startdose=1,maxt=28,accrual=10, maxpen=NA,alpha1=NA,
                      alpha2=NA,cutoff.eli=0.95, ntrial=10,seed=6)


#For Generalized Bayesian optimal interval (gBOIN) design and not Output trial duration as
#an operating characteristics
target<-0.47/1.5
prob <- matrix(c(0.83,	0.75,	0.62,	0.51,	0.34,	0.19,
                 0.12,	0.15,	0.18,	0.19,	0.16,	0.11,
                 0.04,	0.07,	0.11,	0.14,	0.15,	0.11,
                 0.01,	0.03,	0.09,	0.16,	0.35,	0.59), ncol = 6, byrow = TRUE)
get_oc_TITE_QuasiBOIN(target=target, score=c(0,0.5,1,1.5),prob=prob, TITE=FALSE,ncohort=10,
                      cohortsize=3,startdose=1,maxt=NA,accrual=NA, maxpen=NA,alpha1=NA,
                      alpha2=NA,cutoff.eli=0.95, ntrial=10,seed=6)


#For Time-to-event bayesian optimal interval (TITEBOIN) design
get_oc_TITE_QuasiBOIN(target=0.3, score=NA,prob=c(0.25,0.30,0.45,0.49,0.53), TITE=TRUE,
                      ncohort=10, cohortsize=3,startdose=1,maxt=28,accrual=10,
                      maxpen=0.5,alpha1=0.5,alpha2=0.5,cutoff.eli=0.95,
                      ntrial=10,seed=6)


#For Time-to-event generalized bayesian optimal interval (TITEgBOIN) design
target<-0.47/1.5
prob <- matrix(c(0.83,	0.75,	0.62,	0.51,	0.34,	0.19,
                 0.12,	0.15,	0.18,	0.19,	0.16,	0.11,
                 0.04,	0.07,	0.11,	0.14,	0.15,	0.11,
                 0.01,	0.03,	0.09,	0.16,	0.35,	0.59), ncol = 6, byrow = TRUE)
get_oc_TITE_QuasiBOIN(target=target, score=c(0,0.5,1,1.5),prob=prob, TITE=TRUE,ncohort=10,
                      cohortsize=3,startdose=1,maxt=28,accrual=10, maxpen=0.5,alpha1=0.5,
                      alpha2=0.5,cutoff.eli=0.95, ntrial=10,seed=6)
}
