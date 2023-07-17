
#' Determine the dose for the next cohort of new patients for single-agent trials
#'
#' @param target the target toxicity probability (example: target <- 0.30) or the target normalized ETS (example: target <- 0.47 / 1.5)
#' @param n Number of patients treated at each dose level
#' @param npend For TITEBOIN/TITEgBOIN, the number of pending patients at each dose level
#'              For BOIN/gBOIN, "NA" should be assigned
#' @param y:  Number of patients with DLT or the sum of Normalized ETS (equivalent toxicity score)
#' @param ft: For TITEBOIN/TITEgBOIN, Total follow-up time for pending patients for toxicity at each dose level (days)
#'            For BOIN/gBOIN, "NA" should be assigned
#' @param d: current dose level
#' @param maxt For TITEBOIN/TITEgBOIN, Length of assessment window for toxicity (days)
#'             For BOIN/gBOIN, "NA" should be assigned
#' @param p.saf the lower bound. The default value is p.saf=0.6*target
#' @param p.tox the upper bound. The default value is p.tox=1.4*target
#' @param elimination, elimination of each dose (0,1 should be assigned, 0 means the dose is not eliminated,
#'                     1 means the dose is eliminated due to over toxic(\code{elimination=NA}, 0 is defauted for each dose level))
#' @param cutoff.eli the cutoff to eliminate an overly toxic dose for safety.
#'                  We recommend the default value of (\code{cutoff.eli=0.95}) for general use.
#' @param extrasafe set \code{extrasafe=TRUE} to impose a more stringent stopping rule
#' @param offset a small positive number (between 0 and 0.5) to control how strict the
#'               stopping rule is when \code{extrasafe=TRUE}. A larger value leads to a more
#'               strict stopping rule. The default value \code{offset=0.05} generally works well.
#' @param n.earlystop the early stopping parameter. The default value is n.earlystop=100
#' @param maxpen  For TITEBOIN/TITEgBOIN, the upper limit of the ratio of pending patients
#'                For BOIN/gBOIN, "NA" should be assigned
#' @param Neli the sample size cutoff for elimination. The default is \code{Neli=3}.
#' @param gdesign For BOIN and TITEBOIN, "FALSE" should be assigned.
#'        For gBOIN and TITEgBOIN, "TRUE" should be assigned . The default is \code{gdesign=FALSE}.
#'
#' @return \code{next_TITE_QuasiBOIN()} returns the toxicity probability and the recommended dose level for the next cohort
#'              including: (1) the lower Bayesian optimal boundary (\code{lambda_e})
#'              (2) the upper Bayesian optimal boundary (\code{lambda_d})
#'              (3) The number of patients or ESS (the effective sampe size) at each dose level (\code{ESS})
#'              (4) The DLT rate or mu (the estimated quasi-Bernoulli toxicity probability) at each dose level (\code{mu})
#'              (5) the recommended dose level for the next cohort as a numeric value under (\code{d})


next_TITE_QuasiBOIN <- function(target, n, npend, y, ft, d, maxt=28, p.saf = 0.6 * target, p.tox =  1.4 * target,elimination=NA,
                                cutoff.eli = 0.95, extrasafe = FALSE, offset=0.05,n.earlystop = 100,maxpen=0.50,Neli=3,print_d = FALSE,gdesign=FALSE)
{


  earlystop <- 0
  stop <- 0
  elimineed <-0
  pending <- 0
  suspend<-0

  if(is.na(maxt)){maxt = 1}
  if(is.na(maxpen)){maxpen=0.5;}

  ndose <- length(y)
  if(is.na(npend[1])){npend = rep(0,ndose)}
  if(is.na(ft[1])){ft = rep(0,ndose)}

  if(is.na(elimination[1])){
    elimi = rep(0,ndose)
  }else{
    elimi=elimination
  }

  ntox.curr1<-y
  n.curr1<-n
  n.pend1<-npend
  totalt1<-ft/maxt

  #effective sampe size
  e_n<-(n.curr1-n.pend1+totalt1)
  #mu DLT rate;
  mu<-ntox.curr1/(n.curr1-n.pend1+totalt1)


  lambda1 = log((1 - p.saf)/(1 - target))/log(target *
                                                (1 - p.saf)/(p.saf * (1 - target)))
  lambda2 = log((1 - target)/(1 - p.tox))/log(p.tox * (1 -

                                                         target)/(target * (1 - p.tox)))

  ###dose elimination rule: check all the doses for elimination####
  for(dd in 1:ndose){
    #if (1-pbeta(target, ntox.curr1[dd]+1, (n.curr1[dd]-n.pend1[dd]+totalt1[dd])-ntox.curr1[dd]+1)>cutoff.eli && n.curr1[dd]>=Neli){
    if (1-pbeta(target, ntox.curr1[dd]+1, (n.curr1[dd])-ntox.curr1[dd]+1)>cutoff.eli && n.curr1[dd]>=Neli){
      elimi[dd:ndose]=1;
      break;
    }
  }

  ####current dose#####
  ntox.curr<-ntox.curr1[d]
  n.curr<-n.curr1[d]
  n.pend<-n.pend1[d]
  totalt<-totalt1[d]

  #check whether extra safey rule should be applied
  if(extrasafe)
  {
    if(d==1){
      #if(1-pbeta(target, ntox.curr+1, (n.curr-n.pend+totalt)-ntox.curr+1)>cutoff.eli-offset && n.curr>=Neli) {
      if(1-pbeta(target, ntox.curr+1, n.curr-ntox.curr+1)>cutoff.eli-offset && n.curr>=Neli) {
        earlystop = 1; elimi[1:ndose]=1;  }
    }
  }

  #############################################################################
  # Why use n instead of effective n for the safety/extra safety rule?        #
  # If there are many pending people when we check elimination                #
  # (it is possible since our logic checks elimination first and              #
  # then check whether we have enough non-pending people), using n            #
  # can reduce the chance of eliminating ??false positive?? overly-toxic dose   #
  # comparing to using effective n, which results in higher probability of    #
  # selecting MTD.                                                            #
  #############################################################################

  # check whether the current dose level should be eliminated
  if(elimi[d]==1) {
    d=which(elimi==1)[1]-1
    if(d==0){earlystop = 1;}
    elimineed=1
  }
  # else{#check whether the current dose is toxic based on observed data
  #   if((ntox.curr/n.curr)>=lambda2){
  #     if(d==1){d=d; if(n.pend>0){pending=1};} else{d=d-1};
  #     elimineed=1
  #   }
  # } # This part of code doesn't change results, delete it.

  #check whether current dose reaches the "adequate size" of n.earlystop
  if(n.curr>=n.earlystop){
    #earlystop = 1; #JZ 01/28/2023 edit: Change "earlystop=1" to "stop=1" since we will count earlystop as early termination%
    stop=1;
  }


  #####################################################################
  # Assign Dose:                                                      #
  # 1. First check if should early stop trial                         #
  # 2. Else, check if the current dose should be eliminated           #
  # 3. Else, check if should suspend enrollment                       #
  # 4. Else, make dose assignment by comparing mu with lambda_1,2     #
  #####################################################################
  stay<-0
  if(elimineed==1){

  }else{
    if(n.pend>n.curr*maxpen) {pending=1;suspend<-1}
    else
    {
      #compute the estimate of toxicity rate
      if(n.pend==0){
        phat=ntox.curr/n.curr
      } else {
        #compute the estimated toxicity rate based on the imputed data
        phat = ntox.curr/(n.curr-n.pend+totalt)
      }


      #make dose assignment decisions
      if(phat<=lambda1){
        #check whether the current dose is the highest
        if(d==ndose){d=d;stay<-1}
        else{
          if(elimi[d+1]==1){d=d; stay=1} else{ d=d+1}
        }}
      else if (phat>=lambda2 ){
        if(d==1){d=d;stay=1; if(n.pend>0){pending=1}} else{d=d-1}
      }
      else {d=d; stay=1}
    }
  }

  if(stop==1 & stay==1){
    stop=1
  }else{stop=0}

  if (print_d == FALSE) {
    ddose=list()
    ddose=list(d=d,pending=pending,n.pend=n.pend,earlystop=earlystop,stop=stop,elimineed=elimineed,elimi=elimi)
    return(ddose)
  }else{
    if (earlystop == 1){
      message("\n")
      message("Early stop because the 1st dose is too toxic\n")
    }else if(earlystop == 0){
      if(suspend==1){
        message("\n")
        message("More than 50% of patients have not finished the assessment at the current dose level, suspend the dose allocation\n")
      }else{
        if(sum(elimi)==0){
          ddose=list()
          ddose=list(lambda_e = lambda1,lambda_d = lambda2, e_n=e_n,mu=mu,d=d)
          names(ddose)[1]<-'lambda_e'
          names(ddose)[2]<-'lambda_d'
          if(gdesign=="TRUE"){
          names(ddose)[3]<-'The ESS (the effective sampe size) at each dose level'
          names(ddose)[4]<-'The estimated quasi-Bernoulli toxicity probability at each dose level'
          }else{
            names(ddose)[3]<-'The number of patients at each dose level'
            names(ddose)[4]<-'The DLT rate at each dose level'
          }
          names(ddose)[5]<-'Next dose level'
          print(ddose)
        }else{
          elimi_note<-paste0(paste0("Eliminate dose level [",min(which(elimi==1))),"] and all the doses above this dose level.")
          ddose=list()
          ddose=list(lambda_e = lambda1,lambda_d = lambda2, e_n=e_n,mu=mu,d=d,elimi=elimi_note)
          names(ddose)[1]<-'lambda_e'
          names(ddose)[2]<-'lambda_d'
          if(gdesign=="TRUE"){
          names(ddose)[3]<-'The ESS (the effective sampe size) at each dose level'
          names(ddose)[4]<-'The estimated quasi-Bernoulli toxicity probability at each dose level'
          }else{
            names(ddose)[3]<-'The number of patients at each dose level'
            names(ddose)[4]<-'The DLT rate at each dose level'
          }
          names(ddose)[5]<-'Next dose level'
          names(ddose)[6]<-'Dose elimination due to too toxic'
          print(ddose)
        }
      }
    }

  }
}








