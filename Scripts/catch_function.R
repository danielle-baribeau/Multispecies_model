# Here we develop a multi-species model for the North Sea.
# Fishing function development
#This function will be called in the stock loop (starts on line 586 of trophic_model_function) to determine the catch of each stock in a given projection year

# DB Notes
#1: 


#NOTE: writing this as if we are fishing population first, then growing population (discrete-time model)

proj.catch <- function(mgmt, stock.num, repo.loc = "D:/GitHub/Multispecies_model/")
{
#isolate stock-specific management information  
  #get stock-specific management plan
  mgmt.stock <- mgmt.plan |> collapse::fsubset(stock == s)
  #isolate management criteria from stock-specific management plan
  ex.rate.stock <- mgmt.stock$ex.rate
  low.threshold.stock <- mgmt.stock$low.threshold
  high.threshold.stock <- mgmt.stock$high.threshold
  num <- mgmt.stock$stock.num
  interval <- mgmt.stock$assessment.interval
  step <- mgmt.stock$adjustment
  
  #set up assessment update objects
  update <- NULL
  update.type <- NULL

  #"Stock assessment" based on the biomass of the previous year
  #VERSION 1 - upper and lower threshold, with adjustments of constant magnitude that are input by user
  if(t %% interval == 0){
  if(bm.start < low.threshold.stock){
    ex.rate.stock <- ex.rate.stock - adjustment
    update <- 2
    update.type <- paste("assessment; u -",adjustment)
    #apply updated exploitation rate to mgmt plan for future years
    mgmt.plan[num, 1] <- ex.rate.stock
  }#end of "under low threshold" loop
  if(bm.start >= low.threshold.stock && bm.start < high.threshold.stock){
    update <- 1
    update.type <- "assessment; no update needed"
  }#end of "between thresholds" loop
  if(bm.start > high.threshold.stock){
    update <- 2
    update.type <- paste("assessment; u +",adjustment)
  } #end of "above high threshold" loop 
  }# end of "assessment year" loop
  if(t %% interval != 0){
    update <- 0
    update.type <- "no assessment this year"
    }#end of "regular year" loop
  #end of assessment
  
  #find removals equivalent to exploitation rate
  removals <- bm.start * ex.rate.stock
  #apply removals to results
  tst.res <- lam.samp*(bm.start - removals)
  
  
  return(list(tst.res = tst.res, removals = removals, update = update, update.type = update.type))
  
}#end function

#testing catch function
stocks = stock.lst
lambdas= eco.lambdas
n.sims=n.sims
mgmt = list(mgmt =mgmt.plan,er.mn = NULL,er.sd = NULL)
n.yrs.proj= n.yrs.proj
repo.loc=repo.loc