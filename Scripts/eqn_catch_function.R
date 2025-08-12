
# Fishing function development
#This function is called in the stock loop of trophic_model_function simulation,
#and determines the catch and net biomass of a given stock in a given projection year

# DB Notes
#1: mgmt.stock:     # stock-specific management plan details (thresholds, assessment interval, stock ID #, 
                    #and slope/intercept of equation between u and biomass)
#2: repo.loc:       # Location of the Github repo, defaults to "D:/GitHub/Multispecies_model/"

#NOTE: writing this as if we are fishing population first, then growing population (discrete-time model)
#NOTE 2: each updated exploitation rate is applied not to the current projection year, but to the next one 
#(this delay models lags in mgmt implementation)

proj.catch.eqn <- function(mgmt.stock, repo.loc = "D:/GitHub/Multispecies_model/")
{
  #isolate management criteria from stock-specific management plan
  low.threshold.stock <- mgmt.stock$low.threshold
  high.threshold.stock <- mgmt.stock$high.threshold
  num <- mgmt.stock$stock.num
  a.interval <- mgmt.stock$assessment.interval
  slope <- mgmt.stock$slope
  intercept <- mgmt.stock$intercept
  ex.curr <- mgmt.stock$ex.curr
  stock.num <- mgmt.stock$stock.num
  
  #start "stock assessment"
    #set up objects
    update <- NULL
    update.type <- NULL
    ex.next <- NULL
  
    #start "assessment year" loop
    if(t %% a.interval == 0){
      if (bm.start <= low.threshold.stock){
        ex.next <- 0
        update <- 1
        update.type <- "assessment completed; u = 0"

      }#end low.threshold loop
    if(low.threshold.stock < bm.start && bm.start < high.threshold.stock){
      #apply equation to get exploitation rate for bm.start
      ex.next <- slope*(bm.start) + intercept
      update <- 2
      update.type <- paste("assessment completed; u changed to ", ex.next)

    }#end "in between thresholds" loop
    if(bm.start >= high.threshold.stock){
      ex.next <- 0.4
      update <- 3
      update.type <- "assessment completed; u changed to 0.4"

    }#end "above threshold" loop
    }#end "assessment year" loop
    if(t %% a.interval != 0){
      update <- 0
      update.type <- "no assessment this year"
      ex.next <- ex.curr
    }#end of "regular year" loop
  #end of assessment
  
  #find removals that are equivalent to exploitation rate from stock assessment
  removals <- bm.start * ex.curr
  #apply removals to results
  tst.res <- lam.samp*(bm.start - removals)
  
  
  return(list(tst.res = tst.res, removals = removals, ex.curr = ex.curr,
              update = update, update.type = update.type, ex.next = ex.next))
  
}#end function