#### FISHING FUNCTION FOR ESS ####
#LAST UPDATED DB 260429

#Messing with DK's version of fishing function

#1: dat:     # stock-specific management data (mean er and fm, HRC rules, er above/below HRCs)
#2: bm:     #starting biomass
#3: repo.loc:       # Location of the Github repo, defaults to "D:/GitHub/Multispecies_model/"

#NOTE: writing this as if we are fishing population first, then growing population (discrete-time model)
#NOTE 2: writing this so that we are applying updated exploitation rate in the year it is calculated (no time lag)



proj.catch <- function(dat, bm = NULL,repo.loc = "D:/GitHub/Multispecies_model/")
{
  #CHECK: Why fm.mn and not er.mn?
  #if not using HCR...
  if (dat$use.hcr == "NO") {
    er <- rlnorm(1, mean = log(dat$er.mn), sd = dat$er.sd)
    #way to track where stock is in relation to LRP/USR (not managing, just tracking)
    if (bm <= dat$lrp) {
      status <- "Below LRP"
      stat.num <- 1
    }
    if (dat$lrp < bm & bm < dat$usr){
      status <- "Between LRP and USR"
      stat.num <- 2
    }
    
    if (bm >= dat$usr){
      status <- "Above USR"
      stat.num <- 3
    }

  }
  #CHECK: what is "1-exp" doing? Making sure er can't go higher than 1?
  
  #if using HCR...
  #if (dat$use.hcr == "YES"){
    #if (bm <= dat$lrp) 
    #{
      #if(dat$fm.below.lrp == 0) er <- 0 
      #if(dat$fm.below.lrp > 0) er <- 1-exp(-rlnorm(1, mean = log(dat$fm.below.lrp), sd = dat$er.sd))
    #} # end the below lrp if
    #if(bm > dat$lrp  & bm < dat$usr)  er <- 1-exp(-rlnorm(1, mean = log(dat$fm.mn*(bm-dat$lrp)/dat$lrp + dat$fm.below.lrp) , sd = dat$er.sd))
    #if(bm >= dat$urp) ex.next <- 1-exp(-rlnorm(1, mean = log(dat$fm.mn), sd = dat$er.sd))
  #}
  #CHECK:not sure how to get fm at the moment - making alt version with er
  #if (dat$use.hcr == "YES"){
    #if (bm <= dat$lrp){
      #if(dat$er.below.lrp == 0) {
       # er <- 0
       # status <- "Below LRP"
       # stat.num <- 1
      #}
    #}#end of "below lrp" loop
    #if (bm > dat$lrp & bm < dat$usr) {
     # er <- 1-exp(-rlnorm(1, mean = log(dat$er.mn*(bm - dat$lrp)/dat$lrp + dat$er.below.lrp), sd = dat$er.sd))
     # status <- "Between LRP and USR"
      #stat.num <- 2
     # }
    #CHECK: go over behaviour of mean relationship to make sure I understand
    
    #if(bm >= dat$usr) {
     # er <- 1-exp(-rlnorm(1, mean = log(dat$er.mn), sd = dat$er.sd))
     # status <- "Above USR"
     # stat.num <- 3
    #}
    #CHECK: to confirm, er.mn is acting as the cap? I.E. er can't go higher than er.mn
      #is it supposed to work this way? Currently, it IS possible to get er above er.mn when 
      #between LRP and USR (longhorn sculpin)
    
    #}#end HCR loop
  
  return(list(er = er, status = status, stat.num = stat.num))
  
}#end fishing function