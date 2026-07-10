## NEXT STEPS
# get time loop working (need to fix indexing so that you start with empty df
#and fill it as you go) - DONE!

#compare summary sim stats to input sim stats to see how well sim is capturing trends in input data
#V1 is no restrictions on lambdas
  #wolffish, pollock and longfin hake are tanking for no reason
#limit minimum and maximum lambdas (nothing less than 90% decline or tenfold increase for first pass)
#0.1 and 10
  #THIS IS V2 (runs about the same speed as V1)
#can let lambdas go 20% below/above historical min/max (some flexibility but not too wild)
  #THIS IS V3 (runs about the same speed as V1 & V2)


#Dynamics:
  #DK did K-space to biomass ratio
  #if TL doesn't use all its K-space, the proportion of K-space not used creates a direct penalty on the next TL
  #not enough support at lower TL, so higher TL doesn't do as well
  #use distance from ORIGINAL K
    #if you use adjusted K, then dynamics can flip

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
start.com.diff <- start.com.sim - mn.com.sim
#sd.com.sim <- sd(log(com.bm))
sd.com.sim <- sd(com.bm)

L.com.prop <- unique(bm.best.b$prop.eco.bfg[bm.best.b$BFG == "L"])
start.L.com.prop <- L.com.prop[length(L.com.prop)]
L.com.prop.cor <- pacf(L.com.prop)
L.com.prop.cor.lag.1 <- L.com.prop.cor$acf[1]
L.com.prop.cor.lag.2 <- L.com.prop.cor$acf[2]

L.com.prop.logit <- logit(L.com.prop)
start.L.com.prop.logit <- L.com.prop.logit[length(L.com.prop.logit)]


#This doesn't equal the sd of L.com.prop
mn.L.sim.logit <- mean(L.com.prop.logit)
med.L.sim.logit <- median(L.com.prop.logit)
sd.L.sim.logit <- sd(L.com.prop.logit)

start.L.com.prop.diff.logit <- start.L.com.prop.logit - mn.L.sim.logit

#starting proportions for each stock
L.stocks.data <- bm.best.b |> collapse::fsubset(Year == 2017) |> collapse::fsubset(BFG == "L")
L.stocks <- as.character(unique(L.stocks.data$species))
L.bm.start <- L.stocks.data$bm.bfg[1]


M.stocks.data <- bm.best.b |> collapse::fsubset(Year == 2017) |> collapse::fsubset(BFG == "M")
M.stocks <- unique(M.stocks.data$species)


#initialize starting objects
n.yrs.proj <- 30
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
res.L.stock.df <- NULL
res.M.stock.df <- NULL
prop.K.L <- NULL
com.sim <- NULL

n.sims <- 10

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

#initialize a list to hold results of simulation
res.yr.df <- NULL



#include first year biomasses into simulation results data frame
#res.df <- rbind(res.df, start.res.df)

#stocks with lower biomasses (10000 - 1000000) and few high biomass years (4-5 max) tend to tank during projections
  #increasing number of high biomass years by changing lambda cutoff seems to be helping
  #maybe need things to be stock-specific - median of historical lambdas?

for (i in 1:n.sims){
  #browser()
  res.yr.df.tmp <- NULL
  #STARTING WITH YEAR 1 OF PROJECTION:
  #go through entire first year (lambdas, K-space), then continue with other years
  #browser()
  #size of community is autocorrelated
  # [[i]] should go back on this later!!
  
  #this is fine
  #bm.sim.K.com <- data.frame(bm = c(arima.sim(model =list(ar = com.cor.lag.1),n = n.yrs.proj,n.start=1,start.innov = start.com.diff/com.cor.lag.1,
                                              #innov = c(1,rlnorm(n.yrs.proj-1,0,sd.com.sim))*mn.com.sim - mn.com.sim) + mn.com.sim),
                             #Years = 1:n.yrs.proj,sim = i)
  
  bm.sim.K.com <- data.frame(bm = c(arima.sim(model =list(ar = com.cor.lag.1),n = n.yrs.proj,n.start=1,start.innov = 0/com.cor.lag.1,
                                              innov = c(0,rnorm(n.yrs.proj-1,0,sd.com.sim/2))) + mn.com.sim),
                             Years = 1:n.yrs.proj,sim = i)
  
  #bm.sim.K.com <- data.frame(bm = c(arima.sim(model =list(ar = com.cor.lag.1),n = n.yrs.proj,n.start=1,start.innov = start.com.diff/com.cor.lag.1,
                                              #innov = c(0,rnorm(n.yrs.proj-1,0,sd.com.sim/2))) + mn.com.sim),
                             #Years = 1:n.yrs.proj,sim = i)
  #making sd/2 because observational error in survey data creates biologically unreasonable spikiness in ts
  #reducing sd eliminates some variability and helps account for this error
  #also prevents negative biomasses from happening
  #we know variability is too high because simulated biomass went to 0 (and below)
  
  #bm.sim.K.com <- data.frame(x = 1:n.yrs.proj, bm = median(com.bm))
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
  
  #2-lag version was starting from wrong place, switched back to 1 lag version
  #looks too variable
  
  #prop.sim.K.L <- rmultinom(1,bm.sim.K.com$bm, prob = median(L.com.prop))
  
  #reducing sd, changign to 1st order and switching to mean helps
  prop.sim.K.L <-inv.logit(arima.sim(model =list(ar = c(L.com.prop.cor.lag.1)),n = n.yrs.proj,
                                     n.start =1, start.innov = c(start.L.com.prop.diff.logit/L.com.prop.cor.lag.1), 
                                     innov = c(0,rnorm(n.yrs.proj-1,0,sd.L.sim.logit/2))) + mn.L.sim.logit)
  
  #same rationale for sd/2 as above
  
  #prop.sim.K.L <- inv.logit(arima.sim(model = list(ar = c(L.com.prop.cor.lag.1)), n = n.yrs.proj,
                            #n.start = 1, start.innov = start.L.com.prop.diff.logit/L.com.prop.cor.lag.1,
                            #innov = c(1,rlnorm(n.yrs.proj-1,0,sd.L.sim.logit))*med.L.sim.logit - med.L.sim.logit) + med.L.sim.logit)
  
  #ORIGINAL - NOT STARTING FROM 2017 POINT
  #prop.sim.K.L <-inv.logit(arima.sim(model =list(ar = c(L.com.prop.cor.lag.1,L.com.prop.cor.lag.2)),n = n.yrs.proj,
                                     #n.start =2, start.innov = c(start.L.com.prop.diff.logit/L.com.prop.cor.lag.1,
                                                                 #start.L.com.prop.diff.logit/L.com.prop.cor.lag.1), 
                                     #innov = c(0,rnorm(n.yrs.proj-1,0,sd.L.sim.logit)))+med.L.sim.logit)
  
  #prop.sim.K.L <- c(rep(min(in.prop.com.l$L.com.prop), n.yrs.proj))
  
  #this is the K value for the larges
  bm.sim.K.L <- bm.sim.K.com$bm * c(prop.sim.K.L)
  
  #this is the (stand-by) value for the mediums (will be updated later after L dynamics are calculated)
  bm.sim.K.M <- bm.sim.K.com$bm * (1 - c(prop.sim.K.L))
  
  for (t in 1:n.yrs.proj){
    
    #Get stock-specific K values using a multinomial distribution:
    #get proportions of start L K that each stock uses with a multinomial distribution
    
    if (t == 1){
      res.L.stock.df[[i]] <- data.frame(sim = i, year = 2017,
                                   stock = L.stocks, 
                                   fg = L.stocks.data$BFG,
                                   bm.stock = L.stocks.data$bm.stock,
                                   bm.fg = L.stocks.data$bm.bfg,
                                   prop.stock.fg.bm = L.stocks.data$prop.bm.stock.bfg,
                                   bm.com = L.stocks.data$bm.eco,
                                   removals = NA, lambda = NA,stock.K = NA, fg.K = NA, com.K = NA, multinom = NA)
      
      res.M.stock.df[[i]] <- data.frame(sim = i, year = 2017,
                                   stock = M.stocks, 
                                   fg = M.stocks.data$BFG,
                                   bm.stock = M.stocks.data$bm.stock,
                                   bm.fg = M.stocks.data$bm.bfg,
                                   prop.stock.fg.bm = M.stocks.data$prop.bm.stock.bfg,
                                   bm.com = M.stocks.data$bm.eco,
                                   removals = NA, lambda = NA,stock.K = NA, fg.K = NA, com.K = NA, multinom = NA)
     
#mn.prop.stock.fg.bm <- bm.best.b |> collapse::fsubset(Year %in% c(2000:2017)) |> collapse::fgroup_by(species) |> collapse::fsummarize(mn.prop.stock.fg.bm = mean(prop.bm.stock.bfg)) |> ungroup()

mn.prop.stock.fg.bm <- bm.best.b[bm.best.b$Year %in% c(2000:2017),]
mn.stock.fg <- data.frame(stock = c(L.stocks, M.stocks), mn.prop.fg = 0)

for (m in mn.stock.fg$stock){
  #browser()
  mn <- mean(mn.prop.stock.fg.bm$prop.bm.stock.bfg[mn.prop.stock.fg.bm$species == m])
  mn.stock.fg$mn.prop.fg[mn.stock.fg$stock == m] <- mn
}
      
#L.multinom.year <- rmultinom(1, bm.sim.K.L[t], 
                             #prob = (res.L.stock.df[[1]]$prop.stock.fg.bm[res.L.stock.df[[1]]$year == 2017]))

L.multinom.year <- rmultinom(1, bm.sim.K.L[t], 
                             prob = (mn.stock.fg$mn.prop.fg[mn.stock.fg$stock %in% L.stocks]))

    }
    
    if (t > 1) {
      L.multinom.year <- rmultinom(1, bm.sim.K.L[t], 
                                   prob = (res.L.stock.df[[i]]$prop.stock.fg.bm[res.L.stock.df[[i]]$year == t + 2015]))
      
    }

#L Ks are now done!
    
    #now onto L lambdas:
    count <- 0
   # browser()
    for(s in L.stocks)
    {
      # Reset samples
      #browser()
      #if (s == "Atlantic cod") browser()
      count <- count + 1
      stock.lambdas <- lambdas |> collapse::fsubset(common == s)
      #bm.ts.stock <- bm.final[bm.final$Stock == s & bm.final$Year %in% first.year:last.year,]  
      #a=s

      #this needs to be indexed according to time
      bm.last.yr <- res.L.stock.df[[i]]$bm.stock[res.L.stock.df[[i]]$stock == s 
                                              & res.L.stock.df[[i]]$year == t+2016
                                              & res.L.stock.df[[i]]$sim == i]
      
      # Sort out which of the years are low or high bm
      # I'm using 0.5 as the cut off, other options are valid (0.4 is my fav...)

      low.vs.high <- 1
      #low.vs.high <- quantile(stock.lambdas$lambda, probs = c(0.65), na.rm = T)
      #if (s == "Atlantic wolffish") low.vs.high <- 0.35 #GOOD
      #if (s == "Pollock") low.vs.high <- 0.15 #GOOD
      #if (s == "Thorny skate") low.vs.high <- 0.3 
      #if (s == "Atlantic cod") low.vs.high <- 0.3
      low.vs.high.bm <- low.vs.high * median(stock.lambdas$year.minus.one)
      #low.vs.high.bm <- low.vs.high * quantile(stock.lambdas$year.minus.one, probs = c(0.7), na.rm = T)
      #low.vs.high.bm <- low.vs.high * max(stock.lambdas$year.minus.one)
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
        
        if (bm.last.yr < quantile(stock.lambdas$year.minus.one, probs = c(0.5), na.rm = T)){
          lam.med <- median(stock.lambdas$lambda[stock.lambdas$lambda > 1], na.rm = T)
          lam.sd <- sd(log(stock.lambdas$lambda[stock.lambdas$lambda > 1]), na.rm = T)
          
          lam.samp <- rlnorm(1,log(lam.med),lam.sd)
        }
        
       #if(bm.last.yr > quantile(stock.lambdas$year.minus.one, probs = c(0.4), na.rm = T) & bm.last.yr < low.vs.high.bm) 
        #{
         # if(length(low.bm.years) >0)
          #{
           # lam.med <- median(stock.lambdas$lambda[low.bm.years],na.rm=T)
            #lam.sd <- sd(log(stock.lambdas$lambda[low.bm.years]),na.rm=T)
            #lam.med <- median(stock.lambdas$lambda,na.rm=T)
            #lam.sd <- sd(log(stock.lambdas$lambda),na.rm=T)

          # if(length(low.bm.years) == 1) lam.sd <- 0.2 # In case there is just one low biomass year
          #}
          #if(length(low.bm.years) ==0)
          #{
           # lam.med <- median(stock.lambdas$lambda,na.rm=T)
           # lam.sd <- sd(log(stock.lambdas$lambda),na.rm=T)
          #}
          #lam.samp <- rlnorm(1,log(lam.med),lam.sd)
          #if(is.na(lam.samp)) browser()
        #} # end if(bm.last.yr < low.vs.high.bm) 
        
        if(bm.last.yr >= low.vs.high.bm & bm.last.yr < cur.K) 
        {
          lam.med <- median(stock.lambdas$lambda,na.rm=T)
          lam.sd <- sd(log(stock.lambdas$lambda),na.rm=T)
          #lam.med <- median(stock.lambdas$lambda[high.bm.years],na.rm=T)
          #lam.sd <- sd(log(stock.lambdas$lambda[high.bm.years]),na.rm=T)
          lam.samp <- rlnorm(1,log(lam.med),lam.sd)
          
        } # end if(bm.last.yr < low.vs.high.bm) 
        
      } # end if(method != "sample")
      #if(is.na(lam.samp)) browser()
      while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.med),lam.sd)
      #while(lam.samp < 0.1 || lam.samp > 10) lam.samp <- rlnorm(1,log(lam.med),lam.sd)
      #while (lam.samp < (min(stock.lambdas$lambda)) || lam.samp > 1.5*(max(stock.lambdas$lambda))) lam.samp <- rlnorm(1,log(lam.med),lam.sd)

      # Final one, if we are above the K, we are just doing the high biomass scenario for now.
      # Solution, sample from the lambdas at high biomass, but only take lambdas that are <= 1
      if(bm.last.yr >= cur.K) 
      {
        lam.med <- median(stock.lambdas$lambda[high.bm.years],na.rm=T)
        lam.sd <- sd(log(stock.lambdas$lambda[high.bm.years]),na.rm=T)
        #going to switch to stretched beta dist
        #beta distribution bounded by minimum value in ts and maximum value in ts
        #lam.samp <- r_(1, median lambda, min = min observed, max = max observed)
        
        #lam.med <- median(stock.lambdas$lambda[stock.lambdas$lambda < 1], na.rm = T)
        #lam.sd <- sd(log(stock.lambdas$lambda[stock.lambdas$lambda < 1]), na.rm = T)
        lam.samp <- rlnorm(1,log(lam.med),lam.sd)
        #if(is.na(lam.samp)) browser()
        
        ## add restrictions on lambda by adding || statements here ##
        while(lam.samp > 1) lam.samp <- rlnorm(1,log(lam.med),lam.sd)
        #while(lam.samp < 0.1 || lam.samp > 10) lam.samp <- rlnorm(1,log(lam.med),lam.sd)
        #while (lam.samp < (min(stock.lambdas$lambda)) || lam.samp > 1.5*(max(stock.lambdas$lambda))) lam.samp <- rlnorm(1,log(lam.med),lam.sd)
      }
      #while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      
      ex.curr <- 0
      #just make this a vector for x years from now on; same for mgmt
      removals <- bm.last.yr * ex.curr
      #apply removals to results
      tst.res <- lam.samp*(bm.last.yr - removals)

      
      if (t == 1) {
        res.L.stock.df[[i]]$removals[res.L.stock.df[[i]]$stock == s] <- removals
        res.L.stock.df[[i]]$lambda[res.L.stock.df[[i]]$stock == s] <- lam.samp
        res.L.stock.df[[i]]$fg.K[res.L.stock.df[[i]]$stock == s] <- bm.sim.K.L[t]
        res.L.stock.df[[i]]$com.K[res.L.stock.df[[i]]$stock == s] <- bm.sim.K.com$bm[t]
        res.L.stock.df[[i]]$stock.K[res.L.stock.df[[i]]$stock == s] <- cur.K
        #res.L.stock.df[[i]]$prop.stock.fg.bm[res.L.stock.df[[i]]$stock == s] <- cur.K/bm.sim.K.L[t]
        res.L.stock.df[[i]]$multinom[res.L.stock.df[[i]]$stock == s] <- L.multinom.year[count]
        
        res.L.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "L",
                                         bm.stock = tst.res, bm.fg = NA, 
                                         prop.stock.fg.bm = res.L.stock.df,
                                         bm.com = NA,
                                         removals = NA, lambda = NA, stock.K = NA,
                                         fg.K = NA, com.K = NA, multinom = NA)
        
      }
      
      if (t > 1){
        #browser()
        res.L.stock.df[[i]]$removals[res.L.stock.df[[i]]$year == t + 2016 &
                                       res.L.stock.df[[i]]$stock == s] <- removals
        res.L.stock.df[[i]]$lambda[res.L.stock.df[[i]]$year == t + 2016 &
                                     res.L.stock.df[[i]]$stock == s] <- lam.samp
        res.L.stock.df[[i]]$stock.K[res.L.stock.df[[i]]$year == t + 2016 &
                                      res.L.stock.df[[i]]$stock == s] <- cur.K
        res.L.stock.df[[i]]$fg.K[res.L.stock.df[[i]]$year == t + 2016 &
                                   res.L.stock.df[[i]]$stock == s] <- bm.sim.K.L[t]
        res.L.stock.df[[i]]$com.K[res.L.stock.df[[i]]$year == t + 2016 &
                                    res.L.stock.df[[i]]$stock == s] <- bm.sim.K.com$bm[t]
        #res.L.stock.df[[i]]$prop.stock.fg.bm[res.L.stock.df[[i]]$year == t + 2016 &
                                   #res.L.stock.df[[i]]$stock == s] <- 7
        res.L.stock.df[[i]]$multinom[res.L.stock.df[[i]]$year == t + 2016 &
                                       res.L.stock.df[[i]]$stock == s] <- L.multinom.year[count]
        
        res.L.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "L",
                                         bm.stock = tst.res, bm.fg = NA,
                                         prop.stock.fg.bm = NA,
                                         bm.com = NA,
                                         removals = NA, lambda = NA, stock.K = NA,
                                         fg.K = NA, com.K = NA, multinom = NA)
        
        #in timestep t >= 2, add in bm.fg and bm.com 
      }
    
      #CAUTION
    res.L.stock.df[[i]] <- rbind(res.L.stock.df[[i]], res.L.stock.df.tmp)
    #reset temp filler vector
    res.L.stock.df.tmp <- NULL
      
    }#end stock K loop
    #update 2017 data with right simulation year

    #add up all larges

    

    #will split into LP and LB later on
    L.bm <- res.L.stock.df[[i]] |> group_by(year, fg) |> 
      summarize(L.bm = sum(bm.stock, na.rm = TRUE), .groups = "keep") |> 
      arrange(year, fg) |>
      ungroup()
    
    res.L.stock.df[[i]]$bm.fg[res.L.stock.df[[i]]$year == t + 2017] <- L.bm$L.bm[L.bm$year == t + 2017]
    
    
    
    #if (t == 1) {
      
      #res.L.stock.df[[i]]$prop.stock.fg.bm[res.L.stock.df[[i]]$year == 2017] <- 
      #res.L.stock.df[[i]]$stock.K[res.L.stock.df[[i]]$year == 2017]/res.L.stock.df[[i]]$fg.K[res.L.stock.df[[i]]$year == 2017]
      #}
    
    #if (t > 1){
      
    #res.L.stock.df[[i]]$prop.stock.fg.bm[res.L.stock.df[[i]]$year == t + 2016] <- 
      #res.L.stock.df[[i]]$stock.K[res.L.stock.df[[i]]$year == t + 2016]/res.L.stock.df[[i]]$fg.K[res.L.stock.df[[i]]$year == t + 2016]}
    
    print(paste("done larges", s, t, i))
    

    
    
    #this is the difference between realized biomass of L and L K (theoretical size of L)
    if (t == 1) L.K.space <- bm.sim.K.L[t] - start.res.L.df$bm.fg[1]
    if (t > 1) L.K.space <- bm.sim.K.L[t] - bm.sim.K.L[t - 1]
    
    #adjust M K according to realized growth of Ls
    #M.K.real <- bm.sim.K.M[t] + L.K.space
    M.K.real <- bm.sim.K.M[t]
    
    #get K values for mediums for YEAR 1
    #get proportions of start M K that each stock uses with a multinomial distribution
    
    if (t == 1){
      #M.multinom.year <- rmultinom(1, M.K.real, 
      #prob = (res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$year == 2017]))
      
      M.multinom.year <- rmultinom(1, M.K.real, 
      prob = (mn.stock.fg$mn.prop.fg[mn.stock.fg$stock %in% M.stocks]))
    }
    
    if (t > 1) {
      M.multinom.year <- rmultinom(1, M.K.real, 
                                   prob = (res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$year == t + 2015]))
      
    }
    
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
      bm.last.yr <- res.M.stock.df[[i]]$bm.stock[res.M.stock.df[[i]]$stock == s 
                                        & res.M.stock.df[[i]]$year == t+2016
                                        & res.M.stock.df[[i]]$sim == i]
      
      
      # Sort out which of the years are low or high bm
      # I'm using 0.5 as the cut off, other options are valid (0.4 is my fav...)
      #if (s %in% c("Redfish", "Silver hake")){
        #browser()
      #}
      low.vs.high <- 1
      #low.vs.high <- median(stock.lambdas$lambda)
      #if (s == "Longfin hake") low.vs.high <- 0.2 #GOOD
      #if (s == "Sea raven") low.vs.high <- 0.15
      #if (s == "Longhorn sculpin") low.vs.high <- 0.3 #GOOD
      #if (s == "Witch flounder") low.vs.high <- 0.3 #GOOD
      #low.vs.high.bm <- low.vs.high * quantile(stock.lambdas$year.minus.one, probs = c(0.7), na.rm = T)
      low.vs.high.bm <- low.vs.high * median(stock.lambdas$year.minus.one)
      #low.vs.high.bm <- low.vs.high * max(stock.lambdas$year.minus.one)
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
        if (bm.last.yr < quantile(stock.lambdas$year.minus.one, probs = c(0.5), na.rm = T)){
          lam.med <- median(stock.lambdas$lambda[stock.lambdas$lambda > 1], na.rm = T)
          lam.sd <- sd(log(stock.lambdas$lambda[stock.lambdas$lambda > 1]), na.rm = T)
          
          lam.samp <- rlnorm(1,log(lam.med),lam.sd)
        }
        
        # The fun way to do it is to do something multivariate! Note these are instantaneous now!!
        #if(bm.last.yr > quantile(stock.lambdas$year.minus.one, probs = c(0.4), na.rm = T) & bm.last.yr < low.vs.high.bm) 
        #{

          #if(length(low.bm.years) >0)
          #{
            #lam.med <- median(stock.lambdas$lambda[low.bm.years],na.rm=T)
            #lam.sd <- sd(log(stock.lambdas$lambda[low.bm.years]),na.rm=T)
            #lam.med <- median(stock.lambdas$lambda,na.rm=T)
            #lam.sd <- sd(log(stock.lambdas$lambda),na.rm=T)
            #if(length(low.bm.years) == 1) lam.sd <- 0.2 # In case there is just one low biomass year
          #}
          #if(length(low.bm.years) ==0)
          #{
            #lam.med <- median(stock.lambdas$lambda,na.rm=T)
            #lam.sd <- sd(log(stock.lambdas$lambda),na.rm=T)
          #}
          #lam.samp <- rlnorm(1,log(lam.med),lam.sd)
          #if(is.na(lam.samp)) browser()
        #} # end if(bm.last.yr < low.vs.high.bm) 
        
        if(bm.last.yr >= low.vs.high.bm & bm.last.yr < cur.K) 
        {
          lam.med <- median(stock.lambdas$lambda,na.rm=T)
          lam.sd <- sd(log(stock.lambdas$lambda),na.rm=T)
          #lam.med <- median(stock.lambdas$lambda[high.bm.years],na.rm=T)
          #lam.sd <- sd(log(stock.lambdas$lambda[high.bm.years]),na.rm=T)
          lam.samp <- rlnorm(1,log(lam.med),lam.sd)
          
        } # end if(bm.last.yr < low.vs.high.bm) 
        
      } # end if(method != "sample")
      #if(is.na(lam.samp)) browser()
      while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.med),lam.sd)
      #while(lam.samp < 0.1 || lam.samp > 10) lam.samp <- rlnorm(1,log(lam.med),lam.sd)
      #while (lam.samp < (min(stock.lambdas$lambda)) || lam.samp > 1.5*(max(stock.lambdas$lambda))) lam.samp <- rlnorm(1,log(lam.med),lam.sd)
      
      # Final one, if we are above the K, we are just doing the high biomass scenario for now.
      # Solution, sample from the lambdas at high biomass, but only take lambdas that are <= 1
      if(bm.last.yr >= cur.K) 
      {
        lam.med <- median(stock.lambdas$lambda[high.bm.years],na.rm=T)
        lam.sd <- sd(log(stock.lambdas$lambda[high.bm.years]),na.rm=T)
        #going to switch to stretched beta dist
        #beta distribution bounded by minimum value in ts and maximum value in ts
        #lam.samp <- r_(1, median lambda, min = min observed, max = max observed)
        
        #lam.med <- median(stock.lambdas$lambda[stock.lambdas$lambda < 1], na.rm = T)
        #lam.sd <- sd(log(stock.lambdas$lambda[stock.lambdas$lambda < 1]), na.rm = T)
        lam.samp <- rlnorm(1,log(lam.med),lam.sd)
        #if(is.na(lam.samp)) browser()
        while(lam.samp >1 ) lam.samp <- rlnorm(1,log(lam.med),lam.sd)
        #while(lam.samp < 0.1 || lam.samp > 10) lam.samp <- rlnorm(1,log(lam.med),lam.sd)
        #while (lam.samp < (min(stock.lambdas$lambda)) || lam.samp > 1.5*(max(stock.lambdas$lambda))) lam.samp <- rlnorm(1,log(lam.med),lam.sd)
      }
      #while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      
      ex.curr <- 0
      #just make this a vector for x years from now on; same for mgmt
      removals <- bm.last.yr * ex.curr
      #apply removals to results
      tst.res <- lam.samp*(bm.last.yr - removals)
      

      
      if (t == 1) {
        res.M.stock.df[[i]]$removals[res.M.stock.df[[i]]$stock == s] <- removals
        res.M.stock.df[[i]]$lambda[res.M.stock.df[[i]]$stock == s] <- lam.samp
        res.M.stock.df[[i]]$fg.K[res.M.stock.df[[i]]$stock == s] <- M.K.real
        res.M.stock.df[[i]]$com.K[res.M.stock.df[[i]]$stock == s] <- bm.sim.K.com$bm[t]
        res.M.stock.df[[i]]$stock.K[res.M.stock.df[[i]]$stock == s] <- cur.K
        #res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$stock == s] <- cur.K/M.K.real
        res.M.stock.df[[i]]$multinom[res.M.stock.df[[i]]$stock == s] <- M.multinom.year[count]
        
        res.M.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "M",
                                         bm.stock = tst.res, bm.fg = NA, 
                                         prop.stock.fg.bm = NA,
                                         bm.com = NA,
                                         removals = NA, lambda = NA, stock.K = NA,
                                         fg.K = NA, com.K = NA, multinom = NA)
      }
      
      if (t > 1){
        
        res.M.stock.df[[i]]$removals[res.M.stock.df[[i]]$year == t + 2016 &
                                       res.M.stock.df[[i]]$stock == s] <- removals
        res.M.stock.df[[i]]$lambda[res.M.stock.df[[i]]$year == t + 2016 &
                                       res.M.stock.df[[i]]$stock == s] <- lam.samp
        res.M.stock.df[[i]]$stock.K[res.M.stock.df[[i]]$year == t + 2016 &
                                       res.M.stock.df[[i]]$stock == s] <- cur.K
        res.M.stock.df[[i]]$fg.K[res.M.stock.df[[i]]$year == t + 2016 &
                                       res.M.stock.df[[i]]$stock == s] <- M.K.real
        res.M.stock.df[[i]]$com.K[res.M.stock.df[[i]]$year == t + 2016 &
                                       res.M.stock.df[[i]]$stock == s] <- bm.sim.K.com$bm[t]
        #res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$year == t + 2016 &
                                               #res.M.stock.df[[i]]$stock == s] <- cur.K/M.K.real
        res.M.stock.df[[i]]$multinom[res.M.stock.df[[i]]$stock == s &
                                       res.M.stock.df[[i]]$year == t + 2016] <- M.multinom.year[count]
      
        res.M.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "M",
                                         bm.stock = tst.res, bm.fg = NA,
                                         prop.stock.fg.bm = NA,
                                         bm.com = NA,
                                         removals = NA, lambda = NA, stock.K = NA,
                                         fg.K = NA, com.K = NA, multinom = NA)
      }

      res.M.stock.df[[i]] <- rbind(res.M.stock.df[[i]], res.M.stock.df.tmp)
      #reset temp filler vector
      res.M.stock.df.tmp <- NULL
      
    }#end stock K loop
    
    #sum up M biomass
    M.bm <- res.M.stock.df[[i]] |> collapse::fsubset(year == t + 2017) |> group_by(year, fg) |> 
      summarize(M.bm = sum(bm.stock, na.rm = TRUE), .groups = "keep") |> 
      arrange(year, fg) |>
      ungroup()
    
    res.M.stock.df[[i]]$bm.fg[res.M.stock.df[[i]]$year == t + 2017] <- M.bm$M.bm[M.bm$year == t + 2017]
    
    #if (t == 1) {res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$year == 2017] <- 
      #res.M.stock.df[[i]]$stock.K[res.M.stock.df[[i]]$year == 2017]/res.M.stock.df[[i]]$fg.K[res.M.stock.df[[i]]$year == 2017]}
    
    
    
    #if (t > 1) {res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$year == t + 2016] <- 
      #res.M.stock.df[[i]]$stock.K[res.M.stock.df[[i]]$year == t + 2016]/res.M.stock.df[[i]]$fg.K[res.M.stock.df[[i]]$year == t + 2016]}
    
    
    print("done mediums")
    
    #get community biomass for the year
    #browser()
    bm.com <- res.M.stock.df[[i]]$bm.fg[res.M.stock.df[[i]]$year == t + 2017][1] + 
      res.L.stock.df[[i]]$bm.fg[res.L.stock.df[[i]]$year == t + 2017][1]
    res.L.stock.df[[i]]$bm.com[res.L.stock.df[[i]]$year == t + 2017] <- bm.com
    res.M.stock.df[[i]]$bm.com[res.M.stock.df[[i]]$year == t + 2017] <- bm.com
    
    #if (t == 2) browser()
    
    #if (t == 2) browser()
    #results for year t
    #res.yr.df.tmp[[t]] <- rbind(res.L.stock.df[[i]][res.L.stock.df[[i]]$year == t + 2017,], 
                            #res.M.stock.df[[i]][res.M.stock.df[[i]]$year == t + 2017,])
    browser()
  } #end t loop
  #res.yr.df[[i]] <- do.call("rbind", res.yr.df.tmp)
  
  #res.yr.df[[i]] <- rbind(res.L.stock.df[[i]], res.M.stock.df[[i]])
  
prop.K.L[[i]] <- data.frame(sim = i, year = 1:30, prop.com.L.K = c(prop.sim.K.L))
com.sim[[i]] <- data.frame(sim = i, year = 1:30, com.K = c(bm.sim.K.com$bm))
}
res.L.stock <- do.call("rbind", res.L.stock.df)
res.M.stock <- do.call("rbind", res.M.stock.df)

prop.K.L <- do.call("rbind", prop.K.L)
prop.K.L$year <- prop.K.L$year + 2016

com.sim <- do.call("rbind", com.sim)
com.sim$year <- com.sim$year + 2016

res.sim.df <- rbind(res.L.stock, res.M.stock)
#get everything into one data frame

########################################################
#testing - shows individual sims
res.sim.df$sim <- as.factor(res.sim.df$sim)
p.lambdas$sim <- 0

stock.eco.K <- res.sim.df
stock.eco.K$stock.in.fg.K <- stock.eco.K$stock.K/stock.eco.K$fg.K
stock.eco.K <- na.omit(stock.eco.K)
bm.best.p <- bm.best.b
bm.best.p$sim <- as.factor(0)
colnames(bm.best.p)[colnames(bm.best.p) == "species"] <- "stock"


ggplot(stock.eco.K[stock.eco.K$sim == 2,]) +geom_line(aes(x=year, y=stock.in.fg.K, group=sim, colour = sim)) +
  #geom_line(data = bm.best.p, aes(x = Year, y = prop.bm.stock.bfg, group = sim), 
            #colour = "blue")+
  #geom_vline(xintercept = 2017, linetype = 2, colour = "blue") +
  #geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~stock) +  scale_y_log10(name="Proportion of fg K") + 
  guides(colour = guide_legend(nrow = 5))

sim.L.K <- res.sim.df[res.sim.df$fg == "M",]
sim.L.K <- na.omit(sim.L.K)

ggplot(sim.L.K) + geom_line(aes(x=year, y=fg.K, group=sim, colour = sim)) +
  geom_line(data = bm.best.p[bm.best.p$BFG == "M",], aes(x = Year, y = bm.bfg, group = sim), 
            colour = "blue")+
  geom_vline(xintercept = 2017, linetype = 2, colour = "blue") +
  #geom_hline(yintercept = 0, linetype = 2)
  scale_y_log10(name="Med K") + 
  guides(colour = guide_legend(nrow = 5))

#if we do input x, do we get x in the outputs? Good way to check sims
  #sim should return parameters
sim.pacf <- data.frame(sim = n.sims, com.lag.1 = rep(0,n.sims), l.lag.1 = rep(0, n.sims))
for (i in 1:n.sims){
  pacf.com <- pacf(unique(res.sim.df$com.K[res.sim.df$sim == i]), na.action = na.pass)
  pacf.l <- pacf(unique(res.sim.df$fg.K[res.sim.df$sim == i & res.sim.df$fg == "L"]), na.action = na.pass)
  sim.pacf$sim[i] <- i
  sim.pacf$com.lag.1[i] <- pacf.com$acf[1]
  sim.pacf$l.lag.1[i]<- pacf.l$acf[1]
}


in.prop.com.l <- data.frame(sim = 0, year = 2000:2017, L.com.prop = L.com.prop)

ggplot(prop.K.L) +geom_line(aes(x=year, y=prop.com.L.K, group=sim, colour = sim)) +
  geom_line(data = in.prop.com.l, aes(x = year, y = L.com.prop, group = sim), 
            colour = "blue")+
  geom_vline(xintercept = 2017, linetype = 2, colour = "blue") +
  #geom_hline(yintercept = 0, linetype = 2)
  #scale_y_log10(name="Biomass") + 
  guides(colour = guide_legend(nrow = 5))

in.prop.com <- data.frame(sim = 0, year = 2000:2017, com.K = com.bm)

ggplot(com.sim) +geom_line(aes(x=year, y=com.K, group=sim, colour = sim)) +
  #geom_line(data = in.prop.com, aes(x = year, y = com.K, group = sim), 
            #colour = "blue")+
  geom_hline(yintercept = mean(in.prop.com), colour = "blue") +
  geom_vline(xintercept = 2017, linetype = 2, colour = "blue") +
  #geom_hline(yintercept = 0, linetype = 2)
  #scale_y_log10(name="Biomass") + 
  guides(colour = guide_legend(nrow = 5))


########################################################
quants <- res.sim.df |>  collapse::fgroup_by(year,stock,fg) |> collapse::fsummarize(L.50 = quantile(bm.stock,probs=c(0.25),na.rm=T),
                                                                                    med = median(bm.stock,na.rm=T),
                                                                                    U.50 = quantile(bm.stock,probs=c(0.75),na.rm=T))

#comparing inputs stats to simulation stats to see if sim is capturing trends in data
bm.com.test <- res.sim.df[,c(1,2,8)]
bm.com.test <- unique(bm.com.test)
sim.com.bm <- bm.com.test
quants.sim.com.bm <- sim.com.bm |> collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(bm.com,probs=c(0.25),na.rm=T),
                                                                                     med = median(bm.com,na.rm=T),
                                                                                     U.50 = quantile(bm.com,probs=c(0.75),na.rm=T))
mn.in <- data.frame(stock = c(L.stocks, M.stocks), fg = c((rep("L", length(L.stocks))), (rep("M", length(M.stocks)))), mean = 0)
count <- 0

for (i in mn.in$stock){
  count <- count + 1
  stock.in <- p.lambdas$pop.bm[p.lambdas$stock == i]
  mn.stock.in <- mean(stock.in)
  mn.in$mean[count] <- mn.stock.in
}

p.lambdas <- lambdas
colnames(p.lambdas)[colnames(p.lambdas) == "common"] <- "stock"
p.lambdas$mn.in <- rep(0, nrow(p.lambdas))
for (i in 1:nrow(p.lambdas)){
  #browser()
  stock.p <- p.lambdas$stock[i]
  p.lambdas$mn.in[i] <- mn.in$mean[mn.in$stock == stock.p]
}

ggplot(quants[quants$fg == "L",]) +geom_line(aes(x=year, y=L.50, group=stock,), colour = "grey") +
  geom_line(data = quants[quants$fg == "L",], aes(x = year, y = U.50, group = stock), colour = "grey") +
  geom_line(data = quants[quants$fg == "L",], aes(x = year, y = med, group = stock)) +
  geom_line(data = p.lambdas[p.lambdas$FG %in% c("LP", "LB"),], aes(x = year, y = pop.bm, group = stock), 
            colour = "blue")+
  geom_line(data = p.lambdas[p.lambdas$FG %in% c("LP", "LB"),], aes(x = year, y = mn.in, group = stock)) +
  geom_vline(xintercept = 2017, linetype = 2, colour = "blue") +
  facet_wrap(~stock) +  scale_y_log10(name="Biomass") + 
  guides(colour = guide_legend(nrow = 5))

ggplot(quants[quants$fg == "M",]) +geom_line(aes(x=year, y=L.50, group=stock,), colour = "grey") +
  geom_line(data = quants[quants$fg == "M",], aes(x = year, y = U.50, group = stock), colour = "grey") +
  geom_line(data = quants[quants$fg == "M",], aes(x = year, y = med, group = stock)) +
  geom_line(data = p.lambdas[p.lambdas$FG %in% c("MP", "MBZ"),], aes(x = year, y = pop.bm, group = stock), 
            colour = "blue")+
  geom_line(data = p.lambdas[p.lambdas$FG %in% c("MP", "MBZ"),], aes(x = year, y = mn.in, group = stock)) +
  geom_vline(xintercept = 2017, linetype = 2, colour = "blue") +
  facet_wrap(~stock) +  scale_y_log10(name="Biomass") + 
  guides(colour = guide_legend(nrow = 5))


quants.sim.com.bm$L.50 <- quants.sim.com.bm$L.50/1000000
quants.sim.com.bm$U.50 <- quants.sim.com.bm$U.50/1000000
quants.sim.com.bm$med <- quants.sim.com.bm$med/1000000

input.com.bm <- data.frame(year = 2000:2017, bm = com.bm/1000000)

p.sim.in.com <- ggplot(quants.sim.com.bm) +geom_line(aes(x=year, y=L.50), colour = "grey") +
  geom_line(data = quants.sim.com.bm, aes(x = year, y = U.50), colour = "grey") +
  geom_line(data = quants.sim.com.bm, aes(x = year, y = med)) +
  geom_line(data = input.com.bm, aes(x = year, y = bm), 
            colour = "blue")+
  geom_hline(yintercept = median(input.com.bm$bm)) +
  geom_vline(xintercept = 2017, linetype = 2, colour = "blue") #+
#scale_y_log10(name="Biomass (millions)")
p.sim.in.com

########################################################