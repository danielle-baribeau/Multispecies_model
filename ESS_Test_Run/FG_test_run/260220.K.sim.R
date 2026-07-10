## NEXT STEPS
  # get time loop working (need to fix indexing so that you start with empty df
  #and fill it as you go)
    #don't want to always reset mother df in the loop; need to just add to it sequentially
    #put 2017 stuff in mother df after you initialize it, then add on to this df in loop
    #means that removals/lambdas/Ks for 2017 will be on 2018 row, etc.
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
L.stocks.data <- bm.best.b |> collapse::fsubset(Year == 2017) |> collapse::fsubset(BFG == "L")
L.stocks <- as.character(unique(L.stocks.data$species))
L.bm.start <- L.stocks.data$bm.bfg[1]


M.stocks.data <- bm.best.b |> collapse::fsubset(Year == 2017) |> collapse::fsubset(BFG == "M")
M.stocks <- unique(M.stocks.data$species)


#initialize starting objects
n.yrs.proj <- 10
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
start.res.df <- data.frame(sim = 0, year = 2017,
                           stock = c(L.stocks,M.stocks), 
                           bm.stock = c(L.stocks.data$bm.stock,M.stocks.data$bm.stock),
                           bm.fg = c(L.stocks.data$bm.bfg, M.stocks.data$bm.bfg),
                           bm.com = c(L.stocks.data$bm.eco, M.stocks.data$bm.eco),
                           removals = 0, lambda = 0,stock.K = 0, fg.K = 0, com.K = 0)

#initialize a data frame to hold results of simulation
res.df <- NULL

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
  bm.sim.K.L <- bm.sim.K.com$bm * c(prop.sim.K.L)
  
  #this is the (stand-by) value for the mediums (will be updated later after L dynamics are calculated)
  bm.sim.K.M <- bm.sim.K.com$bm * (1 - prop.sim.K.L)
  
  for (t in 1:n.yrs.proj){
    #Get stock-specific K values using a multinomial distribution:
    #get proportions of start L K that each stock uses with a multinomial distribution
    L.multinom.year <- rmultinom(1, bm.sim.K.L[t], prob = (L.stocks.data$prop.bm.stock.bfg))
    
    #if t = 1, store stock Ks in start value data frame
    if (t == 1) start.res.df$stock.K <- L.multinom.year
    
    #regardless of year, store values in this intermediate results data frame
      L.sim.K.stock <- data.frame(Year = t + 2016, 
                                  sim = i,
                                  stock = L.stocks, 
                                  bm.sim.K.stock = L.multinom.year,
                                  fg = "L")

    #Ks are now done!
    
    #now onto L lambdas:
    for(s in L.stocks)
    {
      # Reset samples
      #browser()
      stock.lambdas <- lambdas |> collapse::fsubset(common == s)
      #bm.ts.stock <- bm.final[bm.final$Stock == s & bm.final$Year %in% first.year:last.year,]  
      #a=s
      #browser()
      #res.ts[[s]] <- data.frame(net.bm = L.stocks.data$bm.stock[L.stocks.data$species == s],
                                #removals = NA, mgmt.update = NA,
                                #update.type = NA, current.u = NA, next.yrs.u = NA,
                                #Stock = s,sim= i,lambda = NA,Years=L.stocks.data$Year[1],
                                #fg = "Wrong",
                                #K.bm = NA)
      
      #this needs to be indexed according to time
      if (t == 1)bm.last.yr <- start.res.df$bm.stock[start.res.df$stock == s]
      
      if (t > 1)bm.last.yr <- res.df$bm.stock[res.df$stock == s 
                                                        & res.df$year == t+2016
                                                        & res.df$sim == i]
      
      # Sort out which of the years are low or high bm
      # I'm using 0.5 as the cut off, other options are valid (0.4 is my fav...)
      low.vs.high <- 0.4
      low.vs.high.bm <- low.vs.high * max(stock.lambdas$year.minus.one)
      # Have to drop the final year because we don't have a lambda estimate for the final year
      low.bm.years <- which(stock.lambdas$year.minus.one[-nrow(stock.lambdas)] < low.vs.high.bm)
      if(length(low.bm.years) == 0) low.years <- F else low.years <- T
      high.bm.years <- which(stock.lambdas$year.minus.one[-nrow(stock.lambdas)] >= low.vs.high.bm)
      # ANd get what our carrying capacity is at this moment
      cur.K <- L.sim.K.stock$bm.sim.K.stock[L.sim.K.stock$stock == s &
                                              L.sim.K.stock$Year == t + 2016 &
                                              L.sim.K.stock$sim == i]
      
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
          lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
          #if(is.na(lam.samp)) browser()
        } # end if(bm.last.yr < low.vs.high.bm) 
        
        if(bm.last.yr >= low.vs.high.bm & bm.last.yr < cur.K) 
        {
          lam.med <- median(stock.lambdas$lambda[high.bm.years],na.rm=T)
          lam.sd <- sd(log(stock.lambdas$lambda[high.bm.years]),na.rm=T)
          lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
          
        } # end if(bm.last.yr < low.vs.high.bm) 
        
      } # end if(method != "sample")
      #if(is.na(lam.samp)) browser()
      while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      # Final one, if we are above the K, we are just doing the high biomass scenario for now.
      # Solution, sample from the lambdas at high biomass, but only take lambdas that are <= 1
      if(bm.last.yr >= cur.K) 
      {
        lam.med <- median(stock.lambdas$lambda[high.bm.years],na.rm=T)
        lam.sd <- sd(log(stock.lambdas$lambda[high.bm.years]),na.rm=T)
        #going to switch to stretched beta dist
        #beta distribution bounded by minimum value in ts and maximum value in ts
        #lam.samp <- r_(1, median lambda, min = min observed, max = max observed)
        lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
        #if(is.na(lam.samp)) browser()
        while(lam.samp >1) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      }
      #while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      
      ex.curr <- 0.1
      #just make this a vector for x years from now on; same for mgmt
      removals <- bm.last.yr * ex.curr
      #apply removals to results
      tst.res <- lam.samp*(bm.last.yr - removals)
      
      if (t == 1) {
        start.res.df$removals[start.res.df$stock == s] <- removals
        start.res.df$lambda[start.res.df$stock == s] <- lam.samp
        
        L.sim.bm.stock.tmp <- data.frame(sim = i, year = t + 2016, stock = s,
                                       bm.stock = tst.res, bm.fg = NULL, bm.com = NULL,
                                       removals = NULL, lambda = NULL, stock.K = NULL,
                                       fg.K = NULL, com.K = NULL)
      } 
      
      if (t > 1){
        
        res.df.sim.tmp$removals[]
        res.df.stock.tmp <- data.frame(sim = i, year = t + 2016, stock = s,
                                       bm.stock = tst.res, bm.fg = NULL, bm.com = NULL,
                                       removals = NULL, lambda = NULL, stock.K = NULL,
                                       fg.K = NULL, com.K = NULL)
        
      }
      
        

      
      if (s == "Atlantic cod" && t == 1){
        L.sim.bm.stock.tmp[[t]] <- data.frame(Year = NULL, sim = NULL,
                                     species = NULL, bm.stock = NULL,
                                     lambda = NULL, fg = NULL, removals = NULL,
                                     K.stock = NULL,
                                     K.fg = NULL)
        L.sim.bm.stock.tmp[[t]] <- rbind(L.sim.bm.stock.tmp[[t]], data.frame(Year = t + 2016,
                                                           sim = i,
                                                           species = s, 
                                                           bm.stock = L.stocks.data$bm.stock[L.stocks.data$species == s],
                                                           lambda = NA,
                                                           fg = "L",
                                                           removals = NA,
                                                           K.stock = cur.K,
                                                           K.fg = bm.sim.K.L[t]))
        L.sim.bm.stock.tmp[[t]] <- rbind(L.sim.bm.stock.tmp[[t]], data.frame(Year = t + 2017,
                                                                             sim = i,
                                                                             species = s, 
                                                                             bm.stock = tst.res,
                                                                             lambda = lam.samp,
                                                                             fg = "L",
                                                                             removals = removals,
                                                                             K.stock = cur.K,
                                                                             K.fg = bm.sim.K.L[t]))
      }else{
        L.sim.bm.stock.tmp[[t]] <- rbind(L.sim.bm.stock.tmp[[t]], data.frame(Year = t + 2016,
                                                            sim = i,
                                                           species = s, 
                                                           bm.stock = tst.res,
                                                           lambda = lam.samp,
                                                           fg = "L",
                                                           removals = removals,
                                                           K.stock = cur.K,
                                                           K.fg = bm.sim.K.L[t]))
      }
      
    }#end stock K loop
    L.sim.bm.stock <- do.call("rbind", L.sim.bm.stock.tmp)
    
    #add up all larges
    #will split into LP and LB later on
    L.bm <- L.sim.bm.stock |> group_by(Year, fg) |> 
      summarize(L.bm = sum(bm.stock, na.rm = TRUE), .groups = "keep") |> 
      arrange(Year, fg) |>
      ungroup()
    
    #this is the difference between realized biomass of L and L K (theoretical size of L)
    if (t == 1) L.K.space <- as.vector(bm.sim.K.L)[t] - L.bm.start
    if (t > 1) L.K.space <- as.vector(bm.sim.K.L)[t] - as.vector(bm.sim.K.L)[t - 1]
    
    #adjust M K according to realized growth of Ls
    M.K.real <- c(bm.sim.K.M[t]) + L.K.space
    
    #get K values for mediums for YEAR 1
    #get proportions of start M K that each stock uses with a multinomial distribution
    M.multinom.year <- rmultinom(1, M.K.real[t], prob = (M.stocks.data$prop.bm.stock.bfg))
    M.sim.K.stock <- data.frame(Year = rep(t + 2016, length(M.multinom.year)), 
                                sim = rep(i, length(M.multinom.year)),
                                stock = M.stocks, 
                                bm.sim.K.stock = M.multinom.year,
                                fg = "M")
    
    #now onto M lambdas for year 1:
    
    for(s in M.stocks)
    {
      # Reset samples
      #browser()
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
      if (t == 1)bm.last.yr <- M.stocks.data$bm.stock[M.stocks.data$species == s]
      
      if (t > 1)bm.last.yr <- L.M.sim.bm.stock$bm.stock[L.M.sim.bm.stock$species == s 
                                                        & L.M.sim.bm.stock$Year == t+2016
                                                        & L.M.sim.bm.stock$sim == i]
      
      # Sort out which of the years are low or high bm
      # I'm using 0.5 as the cut off, other options are valid (0.4 is my fav...)
      low.vs.high <- 0.4
      low.vs.high.bm <- low.vs.high * max(stock.lambdas$year.minus.one)
      # Have to drop the final year because we don't have a lambda estimate for the final year
      low.bm.years <- which(stock.lambdas$year.minus.one[-nrow(stock.lambdas)] < low.vs.high.bm)
      if(length(low.bm.years) == 0) low.years <- F else low.years <- T
      high.bm.years <- which(stock.lambdas$year.minus.one[-nrow(stock.lambdas)] >= low.vs.high.bm)
      # ANd get what our carrying capacity is at this moment
      cur.K <- M.sim.K.stock$bm.sim.K.stock[M.sim.K.stock$stock == s &
                                              M.sim.K.stock$Year == t + 2016 &
                                              M.sim.K.stock$sim == i]
      
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
          lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
          #if(is.na(lam.samp)) browser()
        } # end if(bm.last.yr < low.vs.high.bm) 
        
        if(bm.last.yr >= low.vs.high.bm & bm.last.yr < cur.K) 
        {
          lam.med <- median(stock.lambdas$lambda[high.bm.years],na.rm=T)
          lam.sd <- sd(log(stock.lambdas$lambda[high.bm.years]),na.rm=T)
          lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
          
        } # end if(bm.last.yr < low.vs.high.bm) 
        
      } # end if(method != "sample")
      #if(is.na(lam.samp)) browser()
      while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      # Final one, if we are above the K, we are just doing the high biomass scenario for now.
      # Solution, sample from the lambdas at high biomass, but only take lambdas that are <= 1
      if(bm.last.yr >= cur.K) 
      {
        lam.med <- median(stock.lambdas$lambda[high.bm.years],na.rm=T)
        lam.sd <- sd(log(stock.lambdas$lambda[high.bm.years]),na.rm=T)
        #going to switch to stretched beta dist
        #beta distribution bounded by minimum value in ts and maximum value in ts
        #lam.samp <- r_(1, median lambda, min = min observed, max = max observed)
        lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
        #if(is.na(lam.samp)) browser()
        while(lam.samp >1) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      }
      #while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      
      ex.curr <- 0.1
      #just make this a vector for x years from now on; same for mgmt
      removals <- bm.last.yr * ex.curr
      #apply removals to results
      tst.res <- lam.samp*(bm.last.yr - removals)
      
      if (s == "Silver hake" && t == 1){
        M.sim.bm.stock.tmp[[t]] <- data.frame(Year = NULL, sim = NULL,
                                              species = NULL, bm.stock = NULL,
                                              lambda = NULL, fg = NULL, removals = NULL,
                                              K.stock = NULL,
                                              K.fg = NULL)
        M.sim.bm.stock.tmp[[t]] <- rbind(M.sim.bm.stock.tmp[[t]], data.frame(Year = t + 2016,
                                                                             sim = i,
                                                                             species = s, 
                                                                             bm.stock = M.stocks.data$bm.stock[M.stocks.data$species == s],
                                                                             lambda = lam.samp,
                                                                             fg = "M",
                                                                             removals = removals,
                                                                             K.stock = cur.K,
                                                                             K.fg = M.K.real[t]))
      }else{
        M.sim.bm.stock.tmp[[t]] <- rbind(M.sim.bm.stock.tmp[[t]], data.frame(Year = t + 2016,
                                                                             sim = i,
                                                                             species = s, bm.stock = tst.res,
                                                                             lambda = lam.samp,
                                                                             fg = "M",
                                                                             removals = removals,
                                                                             K.stock = cur.K,
                                                                             K.fg = M.K.real[t]))
      }
      
      
    }#end stock K loop
    #get everything up into one data frame
    M.sim.bm.stock <- do.call("rbind", M.sim.bm.stock.tmp)
    
    L.M.sim.bm.stock <- rbind(L.sim.bm.stock, M.sim.bm.stock)
  }
 
    
    
    
    
  
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

for (j in 1:nrow(it.L.Ks)){
  #in the first year of the projection...
  if (it.L.Ks$Years[j] == 1){
    sim.K.stock <- data.frame(Year = NULL, sim = NULL, stock = NULL, prop.L.stock = NULL)
    
    
  }
  
  
sim.stocks.K <- do.call("rbind",sim.K.stock)
sim.L.K <- do.call("rbind",bm.L.Ks)
sim.com.K <- do.call("rbind",bm.sim.K.com)
