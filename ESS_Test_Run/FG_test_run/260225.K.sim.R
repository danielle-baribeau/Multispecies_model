## NEXT STEPS
# get time loop working (need to fix indexing so that you start with empty df
#and fill it as you go)
#don't want to always reset mother df in the loop; need to just add to it sequentially
#put 2016 stuff in mother df after you initialize it, then add on to this df in loop
#means that removals/lambdas/Ks for 2016 will be on 2018 row, etc.
#look at res.ts in the original code and example from H

# time loop
# figures
# management scenarios
# environment
# run sims
# done!



#what inputs do we need?
#community biomass time series + pacf() on this time series

#proportion of community biomass taken up by larges time series + pacf() (logit-transformed)

#for all these, need starting value, mean, start-mean and sd
com.bm <- eco.tot.bm.best$bm.eco    
com.cor <- pacf(com.bm)
com.cor.lag.1 <- com.cor$acf[1]
com.cor.lag.2 <- com.cor$acf[2]

start.com.sim <- com.bm[length(com.bm)]
mn.com.sim <- mean(com.bm)
med.com.sim <- median(com.bm)
start.com.diff <- start.com.sim - med.com.sim
sd.com.sim <- sd(log(com.bm))

L.com.prop <- unique(bm.best.b$prop.eco.bfg[bm.best.b$BFG == "L"])
L.com.prop.cor <- pacf(L.com.prop)
L.com.prop.cor.lag.1 <- L.com.prop.cor$acf[1]
L.com.prop.cor.lag.2 <- L.com.prop.cor$acf[2]

L.com.prop.logit <- logit(L.com.prop)
start.L.com.prop.logit <- L.com.prop.logit[length(L.com.prop.logit)]
mn.L.sim.logit <- mean(L.com.prop.logit)
med.L.sim.logit <- median(L.com.prop.logit)
start.L.com.prop.diff.logit <- start.L.com.prop.logit - med.L.sim.logit
sd.L.sim.logit <- sd(L.com.prop.logit)

#starting proportions for each stock
L.stocks.data <- bm.best.b |> collapse::fsubset(Year == 2016) |> collapse::fsubset(BFG == "L")
L.stocks <- as.character(unique(L.stocks.data$species))
L.bm.start <- L.stocks.data$bm.bfg[1]


M.stocks.data <- bm.best.b |> collapse::fsubset(Year == 2016) |> collapse::fsubset(BFG == "M")
M.stocks <- unique(M.stocks.data$species)


#initialize starting objects
n.yrs.proj <- 2
bm.sim.K.com <- NULL
bm.L.Ks <- NULL
arima.logit <- NULL
sim.K.stock <- NULL
bm.sim.K.stocks <- NULL
sim.stocks.K <- NULL
res.ts <- NULL
sim.bm.stock <- NULL
L.sim.bm.stock.tmp <- NULL
M.sim.bm.stock.tmp <- NULL

n.sims <- 1

lambdas <- do.call("rbind", stocks.lst)
#truncate to years we are looking at
lambdas <- lambdas |> collapse::fsubset(year %in% yrs)


#get start year biomasses (last year of input time series) into a data frame
start.res.df <- data.frame(sim = 1, year = 2017,
                           stock = c(L.stocks,M.stocks), 
                           fg = c(L.stocks.data$BFG, M.stocks.data$BFG),
                           bm.stock = c(L.stocks.data$bm.stock,M.stocks.data$bm.stock),
                           bm.fg = c(L.stocks.data$bm.bfg, M.stocks.data$bm.bfg),
                           bm.com = c(L.stocks.data$bm.eco, M.stocks.data$bm.eco),
                           removals = 0, lambda = 0,stock.K = 0, fg.K = 0, com.K = 0)


start.res.L.df <- data.frame(sim = 1, year = 2017,
                             stock = L.stocks, 
                             fg = L.stocks.data$BFG,
                             bm.stock = L.stocks.data$bm.stock,
                             bm.fg = L.stocks.data$bm.bfg,
                             bm.com = L.stocks.data$bm.eco,
                             removals = 0, lambda = 0,stock.K = 0, fg.K = 0, com.K = 0)

start.res.M.df <- data.frame(sim = 1, year = 2017,
                             stock = M.stocks, 
                             fg = M.stocks.data$BFG,
                             bm.stock = M.stocks.data$bm.stock,
                             bm.fg = M.stocks.data$bm.bfg,
                             bm.com = M.stocks.data$bm.eco,
                             removals = 0, lambda = 0,stock.K = 0, fg.K = 0, com.K = 0)

#initialize a data frame to hold results of simulation
res.sim.df <- NULL
res.yr.df <- NULL
res.yrs.df <- NULL
res.yr.df.tmp <- NULL
res.L.stock.df <- NULL
res.M.stock.df <- NULL
res.stock.df <- NULL

res.ts.yr <- NULL
#include first year biomasses into simulation results data frame
#res.df <- rbind(res.df, start.res.df)

for (i in 1:n.sims){
  #STARTING WITH YEAR 1 OF PROJECTION:
  #go through entire first year (lambdas, K-space), then continue with other years
  #browser()
  #size of community is autocorrelated
  # [[i]] should go back on this later!!
  bm.sim.K.com <- data.frame(bm = c(arima.sim(model =list(ar = com.cor.lag.1),n = n.yrs.proj,n.start=1,start.innov = start.com.diff/com.cor.lag.1,
                                              innov = c(1,rlnorm(n.yrs.proj-1,0,sd.com.sim))*med.com.sim - med.com.sim) + med.com.sim),
                             Years = 1:n.yrs.proj,sim = i)
  #everything before + med.com.sim gives you amount of change in a given year
  #then add this change to the median to get fish numbers
  
  #if no ac, arima assumes white noise
  #if ac added, arima assumes autocorrelated process with normally distributed (Gaussian) background
  #"background" = need arima to tell you how far white noise is going from the mean
  #what is gaussian distributed is the size of the peaks in the arima
  
  #* med and - med gives positive and negative deviations above/below what median started at
  #only got positive values because log scale was making deviates always positive
  #lognormal always thinks in proportions rather than multiplying
  
  #sim.com.bm[[i]] gives you community K
  
  #in terms of ecosystem change, we don't care about autocorrelation
  #if we want declining process = exponentially decaying vector, sequence of deviates between 1 and 0.8 (losing 20% of eco capacity)
  #use the same process that we used to develop the ecosystem arima to develop the deviates
  #meaning take com.sim.K output values, and see if you can get a multiplier or additive number 
  #that changes ecosystem size in a way that matches your ecological system assumption
  #multipying something by a deviate of 1 = same number
  #sequence between 1 and 0.8 that is same length as n.years
  #multiply bm.sim.K.com by that matrix of deviates
  
  #proportion of community taken up by larges is autocorrelated
  prop.sim.K.L <-inv.logit(arima.sim(model =list(ar = c(L.com.prop.cor.lag.1,L.com.prop.cor.lag.2)),n = n.yrs.proj,
                                     n.start =2, start.innov = c(start.L.com.prop.diff.logit/L.com.prop.cor.lag.1,
                                                                 start.L.com.prop.diff.logit/L.com.prop.cor.lag.1), 
                                     innov = c(0,rnorm(n.yrs.proj-1,0,sd.L.sim.logit))) + mn.L.sim.logit)
  #this is the K value for the larges
  bm.sim.K.L <- bm.sim.K.com$bm[t] * c(prop.sim.K.L)
  
  #this is the (stand-by) value for the mediums (will be updated later after L dynamics are calculated)
  bm.sim.K.M <- bm.sim.K.com$bm[t] * (1 - c(prop.sim.K.L))
  
  for (t in 1:n.yrs.proj){
    if (t == 1) res.ts.yr[[t]] <- data.frame(sim = i, 
                                year = t + 2017,
                                stock = rep(L.stocks, 2), 
                                fg = "L",
                                bm.stock = start.,
                                bm.fg = M.stocks.data$bm.bfg,
                                bm.com = M.stocks.data$bm.eco,
                                removals = 0, lambda = 0,stock.K = 0, fg.K = 0, com.K = 0)
    
    if (t > 1) res.ts.yr[[t]] <- data.frame(sim = i, year = 2017,
                                            stock = M.stocks, 
                                            fg = M.stocks.data$BFG,
                                            bm.stock = M.stocks.data$bm.stock,
                                            bm.fg = M.stocks.data$bm.bfg,
                                            bm.com = M.stocks.data$bm.eco,
                                            removals = 0, lambda = 0,stock.K = 0, fg.K = 0, com.K = 0)
    #Get stock-specific K values using a multinomial distribution:
    #get proportions of start L K that each stock uses with a multinomial distribution
    L.multinom.year <- rmultinom(1, bm.sim.K.L[t], prob = (L.stocks.data$prop.bm.stock.bfg))
    
    #L Ks are now done!
    
    #now onto L lambdas:
    count <- 0
    for(s in L.stocks)
    {
      # Reset samples
      #browser()
      count <- count + 1
      stock.lambdas <- lambdas |> collapse::fsubset(common == s)
      #bm.ts.stock <- bm.final[bm.final$Stock == s & bm.final$Year %in% first.year:last.year,]  
      #a=s
      #browser()
      
      #this needs to be indexed according to time
      if (t == 1)bm.last.yr <- start.res.L.df$bm.stock[start.res.df$stock == s]
      
      if (t > 1)bm.last.yr <- res.L.stock.df$bm.stock[res.L.stock.df$stock == s 
                                                      & res.L.stock.df$year == t+2016
                                                      & res.L.stock.df$sim == i]
      
      # Sort out which of the years are low or high bm
      # I'm using 0.5 as the cut off, other options are valid (0.4 is my fav...)
      low.vs.high <- 0.4
      low.vs.high.bm <- low.vs.high * max(stock.lambdas$year.minus.one)
      # Have to drop the final year because we don't have a lambda estimate for the final year
      low.bm.years <- which(stock.lambdas$year.minus.one[-nrow(stock.lambdas)] < low.vs.high.bm)
      if(length(low.bm.years) == 0) low.years <- F else low.years <- T
      high.bm.years <- which(stock.lambdas$year.minus.one[-nrow(stock.lambdas)] >= low.vs.high.bm)
      # ANd get what our carrying capacity is at this moment
      cur.K <- L.multinom.year[count]
      
      # So first, get a sample 
      method <- "not_sample"
      if(method == 'sample')
      {
        # Pick one of these to sample if that's how we want to roll, if we have low biomass years (as
        # defined by our cut off low.vs.high)
        if(low.years == T)
        {
          if(bm.last.yr < low.vs.high.bm) samp <- sample(low.bm.years,1)
          if(bm.last.yr >= low.vs.high.bm) samp <- sample(high.bm.years,1)
        } # end If we have low years
        if(low.years == F) samp <- sample(high.bm.years,1)
        # The simple way to do it is just to sample from the natural mortality distribution
        # Now using the right lambda, go look at trends from the stocks that are declining to see what's up there.
        lam.samp <- stock.lambdas$lambda[samp] # Get the sample years.  
      } # end the sample method.
      
      # Or do it the fun way...
      if(method != "sample")
      {
        
        # The fun way to do it is to do something multivariate! Note these are instantaneous now!!
        if(bm.last.yr < low.vs.high.bm) 
        {
          if(length(low.bm.years) >0)
          {
            lam.med <- median(stock.lambdas$lambda[low.bm.years],na.rm=T)
            lam.sd <- sd(log(stock.lambdas$lambda[low.bm.years]),na.rm=T)
            if(length(low.bm.years) == 1) lam.sd <- 0.2 # In case there is just one low biomass year
          }
          if(length(low.bm.years) ==0)
          {
            lam.med <- median(stock.lambdas$lambda,na.rm=T)
            lam.sd <- sd(log(stock.lambdas$lambda),na.rm=T)
          }
          lam.samp <- rlnorm(1,log(lam.med),lam.sd)
          #if(is.na(lam.samp)) browser()
        } # end if(bm.last.yr < low.vs.high.bm) 
        
        if(bm.last.yr >= low.vs.high.bm & bm.last.yr < cur.K) 
        {
          lam.med <- median(stock.lambdas$lambda[high.bm.years],na.rm=T)
          lam.sd <- sd(log(stock.lambdas$lambda[high.bm.years]),na.rm=T)
          lam.samp <- rlnorm(1,log(lam.med),lam.sd)
          
        } # end if(bm.last.yr < low.vs.high.bm) 
        
      } # end if(method != "sample")
      #if(is.na(lam.samp)) browser()
      while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.med),lam.sd)
      # Final one, if we are above the K, we are just doing the high biomass scenario for now.
      # Solution, sample from the lambdas at high biomass, but only take lambdas that are <= 1
      if(bm.last.yr >= cur.K) 
      {
        lam.med <- median(stock.lambdas$lambda[high.bm.years],na.rm=T)
        lam.sd <- sd(log(stock.lambdas$lambda[high.bm.years]),na.rm=T)
        #going to switch to stretched beta dist
        #beta distribution bounded by minimum value in ts and maximum value in ts
        #lam.samp <- r_(1, median lambda, min = min observed, max = max observed)
        lam.samp <- rlnorm(1,log(lam.med),lam.sd)
        #if(is.na(lam.samp)) browser()
        while(lam.samp >1) lam.samp <- rlnorm(1,log(lam.med),lam.sd)
      }
      #while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      
      ex.curr <- 0.1
      #just make this a vector for x years from now on; same for mgmt
      removals <- bm.last.yr * ex.curr
      #apply removals to results
      tst.res <- lam.samp*(bm.last.yr - removals)
      
      if (t == 1) {
        start.res.L.df$removals[start.res.L.df$stock == s] <- removals
        start.res.L.df$lambda[start.res.L.df$stock == s] <- lam.samp
        start.res.L.df$fg.K[start.res.L.df$stock == s] <- bm.sim.K.L[t]
        start.res.L.df$com.K[start.res.L.df$stock == s] <- bm.sim.K.com$bm[t]
        start.res.L.df$stock.K[start.res.L.df$stock == s] <- cur.K
        
        res.L.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "L",
                                         bm.stock = tst.res, bm.fg = 0, 
                                         bm.com = 0,
                                         removals = 0, lambda = 0, stock.K = 0,
                                         fg.K = 0, com.K = 0)
      }
      
      if (t > 1){
        res.L.stock.df$removals[res.L.stock.df$stock == s & 
                                  res.L.stock.df$year == 2016 + t] <- removals
        res.L.stock.df$lambda[res.L.stock.df$stock == s & 
                                res.L.stock.df$year == 2016 + t] <- lam.samp
        res.L.stock.df$stock.K[res.L.stock.df$stock == s & 
                                 res.L.stock.df$year == 2016 + t] <- cur.K
        res.L.stock.df$fg.K[res.L.stock.df$stock == s & 
                              res.L.stock.df$year == 2016 + t] <- bm.sim.K.L[t]
        res.L.stock.df$com.K[res.L.stock.df$stock == s & 
                               res.L.stock.df$year == 2016 + t] <- bm.sim.K.com$bm[t]
        res.L.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "L",
                                         bm.stock = tst.res, bm.fg = 0, 
                                         bm.com = 0,
                                         removals = 0, lambda = 0, stock.K = 0,
                                         fg.K = 0, com.K = 0)
      }
      
      res.L.stock.df <- rbind(res.L.stock.df, res.L.stock.df.tmp)
      #reset temp filler vector
      res.L.stock.df.tmp <- NULL
      
    }#end stock K loop
    
    #add up all larges
    #will split into LP and LB later on
    L.bm <- res.L.stock.df |> group_by(year, fg) |> 
      summarize(L.bm = sum(bm.stock, na.rm = TRUE), .groups = "keep") |> 
      arrange(year, fg) |>
      ungroup()
    
    res.L.stock.df$bm.fg <- L.bm$L.bm[t]
    
    print("done larges")
    
    
    #this is the difference between realized biomass of L and L K (theoretical size of L)
    if (t == 1) L.K.space <- bm.sim.K.L[t] - start.res.L.df$bm.fg[1]
    if (t > 1) L.K.space <- bm.sim.K.L[t] - bm.sim.K.L[t - 1]
    
    #adjust M K according to realized growth of Ls
    M.K.real <- bm.sim.K.M[t] + L.K.space
    
    #get K values for mediums for YEAR 1
    #get proportions of start M K that each stock uses with a multinomial distribution
    M.multinom.year <- rmultinom(1, M.K.real, prob = (M.stocks.data$prop.bm.stock.bfg))
    
    #now onto M lambdas for year 1:
    
    count <- 0
    for(s in M.stocks)
    {
      # Reset samples
      #browser()
      count <- count + 1
      stock.lambdas <- lambdas |> collapse::fsubset(common == s)
      #bm.ts.stock <- bm.final[bm.final$Stock == s & bm.final$Year %in% first.year:last.year,]  
      #a=s
      #browser()
      #res.ts[[s]] <- data.frame(net.bm = M.stocks.data$bm.stock[M.stocks.data$species == s],
      #removals = NA, mgmt.update = NA,
      #update.type = NA, current.u = NA, next.yrs.u = NA,
      #Stock = s,sim= i,lambda = NA,Years=M.stocks.data$Year[1],
      #fg = "Wrong",
      #K.bm = NA)
      #this needs to be indexed according to time
      if (t == 1)bm.last.yr <- start.res.df$bm.stock[start.res.df$stock == s]
      
      if (t > 1)res.M.stock.df$bm.stock[res.M.stock.df$stock == s 
                                        & res.M.stock.df$year == t+2016
                                        & res.M.stock.df$sim == i-1]
      
      # Sort out which of the years are low or high bm
      # I'm using 0.5 as the cut off, other options are valid (0.4 is my fav...)
      low.vs.high <- 0.4
      low.vs.high.bm <- low.vs.high * max(stock.lambdas$year.minus.one)
      # Have to drop the final year because we don't have a lambda estimate for the final year
      low.bm.years <- which(stock.lambdas$year.minus.one[-nrow(stock.lambdas)] < low.vs.high.bm)
      if(length(low.bm.years) == 0) low.years <- F else low.years <- T
      high.bm.years <- which(stock.lambdas$year.minus.one[-nrow(stock.lambdas)] >= low.vs.high.bm)
      # ANd get what our carrying capacity is at this moment
      cur.K <- M.multinom.year[count]
      
      # So first, get a sample 
      method <- "not_sample"
      if(method == 'sample')
      {
        # Pick one of these to sample if that's how we want to roll, if we have low biomass years (as
        # defined by our cut off low.vs.high)
        if(low.years == T)
        {
          if(bm.last.yr < low.vs.high.bm) samp <- sample(low.bm.years,1)
          if(bm.last.yr >= low.vs.high.bm) samp <- sample(high.bm.years,1)
        } # end If we have low years
        if(low.years == F) samp <- sample(high.bm.years,1)
        # The simple way to do it is just to sample from the natural mortality distribution
        # Now using the right lambda, go look at trends from the stocks that are declining to see what's up there.
        lam.samp <- stock.lambdas$lambda[samp] # Get the sample years.  
      } # end the sample method.
      
      # Or do it the fun way...
      if(method != "sample")
      {
        
        # The fun way to do it is to do something multivariate! Note these are instantaneous now!!
        if(bm.last.yr < low.vs.high.bm) 
        {
          if(length(low.bm.years) >0)
          {
            lam.med <- median(stock.lambdas$lambda[low.bm.years],na.rm=T)
            lam.sd <- sd(log(stock.lambdas$lambda[low.bm.years]),na.rm=T)
            if(length(low.bm.years) == 1) lam.sd <- 0.2 # In case there is just one low biomass year
          }
          if(length(low.bm.years) ==0)
          {
            lam.med <- median(stock.lambdas$lambda,na.rm=T)
            lam.sd <- sd(log(stock.lambdas$lambda),na.rm=T)
          }
          lam.samp <- rlnorm(1,log(lam.med),lam.sd)
          #if(is.na(lam.samp)) browser()
        } # end if(bm.last.yr < low.vs.high.bm) 
        
        if(bm.last.yr >= low.vs.high.bm & bm.last.yr < cur.K) 
        {
          lam.med <- median(stock.lambdas$lambda[high.bm.years],na.rm=T)
          lam.sd <- sd(log(stock.lambdas$lambda[high.bm.years]),na.rm=T)
          lam.samp <- rlnorm(1,log(lam.med),lam.sd)
          
        } # end if(bm.last.yr < low.vs.high.bm) 
        
      } # end if(method != "sample")
      #if(is.na(lam.samp)) browser()
      while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.med),lam.sd)
      # Final one, if we are above the K, we are just doing the high biomass scenario for now.
      # Solution, sample from the lambdas at high biomass, but only take lambdas that are <= 1
      if(bm.last.yr >= cur.K) 
      {
        lam.med <- median(stock.lambdas$lambda[high.bm.years],na.rm=T)
        lam.sd <- sd(log(stock.lambdas$lambda[high.bm.years]),na.rm=T)
        #going to switch to stretched beta dist
        #beta distribution bounded by minimum value in ts and maximum value in ts
        #lam.samp <- r_(1, median lambda, min = min observed, max = max observed)
        lam.samp <- rlnorm(1,log(lam.med),lam.sd)
        #if(is.na(lam.samp)) browser()
        while(lam.samp >1) lam.samp <- rlnorm(1,log(lam.med),lam.sd)
      }
      #while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      
      ex.curr <- 0.1
      #just make this a vector for x years from now on; same for mgmt
      removals <- bm.last.yr * ex.curr
      #apply removals to results
      tst.res <- lam.samp*(bm.last.yr - removals)
      
      if (t == 1) {
        start.res.M.df$removals[start.res.M.df$stock == s] <- removals
        start.res.M.df$lambda[start.res.M.df$stock == s] <- lam.samp
        start.res.M.df$fg.K[start.res.M.df$stock == s] <- M.K.real
        start.res.M.df$com.K[start.res.M.df$stock == s] <- bm.sim.K.com$bm[t]
        start.res.M.df$stock.K[start.res.M.df$stock == s] <- cur.K
        
        res.M.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "M",
                                         bm.stock = tst.res, bm.fg = 0, bm.com = 0,
                                         removals = 0, lambda = 0, stock.K = 0,
                                         fg.K = 0, com.K = 0)
      }
      
      if (t > 1){
        res.M.stock.df$removals[res.M.stock.df$stock == s & 
                                  res.M.stock.df$year == 2016 + t] <- removals
        res.M.stock.df$lambda[res.M.stock.df$stock == s & 
                                res.M.stock.df$year == 2016 + t] <- lam.samp
        res.M.stock.df$stock.K[res.M.stock.df$stock == s & 
                                 res.M.stock.df$year == 2016 + t] <- cur.K
        res.M.stock.df$fg.K[res.M.stock.df$stock == s & 
                              res.M.stock.df$year == 2016 + t] <- bm.sim.K.L[t]
        res.M.stock.df$com.K[res.M.stock.df$stock == s & 
                               res.M.stock.df$year == 2016 + t] <- bm.sim.K.com$bm[t]
        res.M.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "M",
                                         bm.stock = tst.res, bm.fg = 0, bm.com = 0,
                                         removals = 0, lambda = 0, stock.K = 0,
                                         fg.K = 0, com.K = 0)
      }
      
      res.M.stock.df <- rbind(res.M.stock.df, res.M.stock.df.tmp)
      #reset temp filler vector
      res.M.stock.df.tmp <- NULL
      
    }#end stock K loop
    
    #sum up M biomass
    M.bm <- res.M.stock.df |> collapse::fsubset(year == t + 2017) |> group_by(year, fg) |> 
      summarize(M.bm = sum(bm.stock, na.rm = TRUE), .groups = "keep") |> 
      arrange(year, fg) |>
      ungroup()
    
    res.M.stock.df$bm.fg <- M.bm$M.bm[t]
    print("done mediums")
    
    #get community biomass for the year
    bm.com <- res.M.stock.df[res.M.stock.df$year == t + 2017,]$bm.fg[1] + 
      res.L.stock.df[res.L.stock.df$year == t + 2017,]$bm.fg[1]
    res.L.stock.df$bm.com <- bm.com
    res.M.stock.df$bm.com <- bm.com
    
    #combine year 0 with first projection year
    if (t == 1) {
      res.L.stock.df <- rbind(start.res.L.df, res.L.stock.df)
      res.M.stock.df <- rbind(start.res.M.df, res.M.stock.df)
    }
    
    #results for year t
    res.yr.df <- rbind(res.L.stock.df, res.M.stock.df)
    
  }
  #get everything into one data frame
  res.sim.df <- rbind(res.sim.df, res.yr.df)
  
  ########################################################
  
  #large pop dynamics need to happen in separate place from mediums
  #get down to lambdas with larges, then look at mediums
  
  #next steps for larges
  #Ks - from L to individual stocks
  #population dynamics of each of the species
  #compare population dynamics back to K (sum population dynamics to get how much larges are actually using)
  #run all these steps in one go
  #can add in fishing pressure, get everything set up - but just for larges
  #ignore LP and LB for now
  #once all the large stuff is working, think about deviates on ecosystem size (environmental change)
  #THEN think about mediums
  
  ########################################################
}
res.sim.df.s <- res.sim.df[order(res.sim.df$stock),]

#IT RUNNSSSSSSSSSSSS!!!!!

#Now check if it runs properly
#it does not :(

for (j in 1:nrow(it.L.Ks)){
  #in the first year of the projection...
  if (it.L.Ks$Years[j] == 1){
    sim.K.stock <- data.frame(Year = NULL, sim = NULL, stock = NULL, prop.L.stock = NULL)
    
    
  }
  
  
  if (t == 1){
    res.L.ts[[s]] <- data.frame(sim = i, year = 2017,
                                stock = s, 
                                fg = "L",
                                bm.stock = start.res.L.df$bm.stock[start.res.df$stock == s]$bm.stock,
                                bm.fg = start.res.L.df$bm.stock[start.res.df$stock == s]$bm.fg,
                                bm.com = start.res.L.df$bm.stock[start.res.df$stock == s]$bm.com,
                                removals = 0, lambda = 0,stock.K = 0, fg.K = 0, com.K = 0)
  }
  
  if (t > 1){ res.L.ts[[s]] <- data.frame(sim = 1, year = t,
                                          stock = s, 
                                          fg = "L",
                                          bm.stock = M.stocks.data$bm.stock,
                                          bm.fg = M.stocks.data$bm.bfg,
                                          bm.com = M.stocks.data$bm.eco,
                                          removals = 0, lambda = 0,stock.K = 0, fg.K = 0, com.K = 0)
  
  }
  
  
  
  sim.stocks.K <- do.call("rbind",sim.K.stock)
  sim.L.K <- do.call("rbind",bm.L.Ks)
  sim.com.K <- do.call("rbind",bm.sim.K.com)
  