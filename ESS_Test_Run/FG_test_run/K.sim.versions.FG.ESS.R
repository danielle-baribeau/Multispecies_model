# This is a function that can run the three different candidate K sims I am considering for the ESS MPM #

# Notes
#1: version: # which version of the K sim you want to run
              # can be version 1 (TL logic), 
              # version 2 (t-d, LP controls size of ecosystem)
              # input "1" or "2"
              #defaults to version 2 (because I think it is more right than version 1)
#2: props     #data frame with acfs/ccfs of all biomasses/proportions needed in each version:
  #eco.bm    # ecosystem biomass (only actually used in 1; fix this)
  #lp.prop.eco   #proportion of ecosystem biomass made of LP (only used in 1 and 2)
  #lb.prop.lb.m  #proportion of LB + M biomass made up of LB (only used in 2)
  #mp.prop.m   #proportion of M biomass made up of MP (only used in 2)
  #mbz.prop.m  #proportion of M biomass made up of MBZ (only used in 3)
  #m.prop.eco  #proportion of ecosystem biomass made up of M (only used in 3)
  #lb.prop.l #proportion of L biomass made up of LB (only used in 3)


sim.K.fg<-function(version = 2,props = props,
                      repo.loc = "D:/GitHub/Multispecies_model",method = "not_sample")
  
#ISSUE - NOT PARTITIONING MEDIUMS IN ANY OF THESE VERSIONS! NON-SIGNIFICANT PACFs
{
  if (version == 1){
    #main issue with version 1 - no significant autocorrelation structure in K.cor!
    for(i in 1:n.sims) 
    {
      browser()
      # The ecosystem K, using the mean of the ecosystem with the correlation observed of the time series.
      # This starts the time series at the last value of the time series, then moves it to the mean value, bam!!  This will be done for each of these arima sims.
      sim.eco.bm[[i]] <- data.frame(bm = c(arima.sim(model =list(ar = K.cor$acf[1]),n = n.yrs.proj,n.start=1,start.innov = start.eco.diff/K.cor$acf[1],
                                                     innov = c(0,rnorm(n.yrs.proj-1,0,sd.eco.bm))) + mn.eco.bm),
                                    Years = 1:n.yrs.proj,sim = i) 
      #pacf(sim.eco.bm[[i]]$bm)
      #NO SIGNIFICANT CORRELATIONS HERE - DOES THIS MEAN MODEL IS WORKING PROPERLY?
      
      # So then from my simulated ecosystem I want each FG to get its cut of the biomass, 
      # FIX: I am using the AR2, but I know the start innovation is slightly incorrect, but it make almost no difference for the NS case so I'll stick with it
      # so probably should figure out how to specify that right as it just works by luck here I think, if the difference was larger
      # or correlations different it wouldn't do so well (e.g., it isn't nice for the stock level ones.)
      sim.tl.3.prop.bm <-inv.logit(mn.tl.3.logit + 
                                     arima.sim(model =list(ar = c(tl.3.prop.bm.lag.1,tl.3.prop.bm.lag.2)),n = n.yrs.proj,
                                               n.start =2, start.innov = c(start.tl.3.diff/tl.3.prop.bm.lag.1,start.tl.3.diff/tl.3.prop.bm.lag.1), 
                                               innov = c(0,rnorm(n.yrs.proj-1,0,sd.tl.3.logit))))
      
      bm.sim.3 <- sim.tl.3.prop.bm * sim.eco.bm[[i]]$bm
      # So this is what is left for 3 and 4
      bm.sim.4 <- sim.eco.bm[[i]]$bm - bm.sim.3
      #bm.left.4.5<- sim.eco.bm[[i]]$bm - bm.sim.3
      # So then we use the historical split between 4 and 5 can see 5 gets about 1/3-1-5 of 3
      # so then simulate this split
      #sim.tl.5.4.prop.bm <- inv.logit(mn.tl.4.5.logit + 
      # arima.sim(model =list(ar = tl.4.5.prop.bm.lag.1),n = n.yrs.proj,
      #  n.start =1, start.innov = c(start.tl.4.5.diff/tl.4.5.prop.bm.lag.1), 
      #  innov = c(0,rnorm(n.yrs.proj-1,0,sd.tl.3.logit))))
      # And now TL 5 gets this proportion of the 4 and 5 biomass
      #bm.sim.5 <- bm.left.4.5 * sim.tl.5.4.prop.bm
      # And TL4 gets the rest, and so the ecosystem biomass is a portion of the whole biomass
      #bm.sim.4 <- sim.eco.bm[[i]]$bm - bm.sim.3-bm.sim.5
      
      bm.trophic.Ks[[i]] <- data.frame(Years = rep(1:n.yrs.proj,2), sim =i,
                                       bm.tl = c(bm.sim.3,bm.sim.4),troph.cat = as.factor(sort(rep(c(3,4),n.yrs.proj))),
                                       bm.eco = rep(sim.eco.bm[[i]]$bm,2))
      bm.trophic.Ks[[i]]$prop.bm.tl <- bm.trophic.Ks[[i]]$bm.tl/bm.trophic.Ks[[i]]$bm.eco
      
      # OK, so now we have the trophic level K values simulated in a 'nice' way. Next how do we partition these to the stocks
      # Give each stock a proportion of the K in it's ecosystem based on their historical cuts of the K, and include the time series correlation in that.
      # I'm going to build in correlation to their K time series (this could 100% be fishery induced correlation), could also put in 
      # cross correlation for species with multiple stocks, but for now, let's just do the AR1/2 thing with this for the proportion of the trophic level 
      # biomass each stock gets.
      
      for(tl in troph.levels)
      {
        tl.stocks <- unique(bm.best$Stock[bm.best$troph.cat==tl])
        n.stock.tl <- length(tl.stocks)
        count =0
        for(s in tl.stocks)
        {
          count = count+1
          # Now get the time series for each stock...
          if(count == 1 ||  n.stock.tl != 2)
          {
            tmp.dat <- bm.best[bm.best$Stock ==s,]
            tmp.cor <- pacf(tmp.dat$prop.bm.tl,plot=F) # Get the correlation, use AR1 and AR2 but no more.
            tmp.cor.lag.1 <- tmp.cor$acf[1]
            #tmp.cor.lag.2 <- tmp.cor$acf[2]
            #tmp.beta <- estBetaParams(mean(tmp.dat$prop.bm.tl),sd(tmp.dat$prop.bm.tl)^2)
            # Logit tranform the proportions and do the ARIMA on the logits
            bm.logit <- logit(tmp.dat$prop.bm.stock.tl)
            start.bm.logit <- bm.logit[length(bm.logit)]
            mn.bm.logit <- mean(bm.logit)
            # Use the most recent bm on logit scale as the 'mean' value for the simulation
            mn.bm.logit <- start.bm.logit
            # And the standard deviation
            sd.bm.logit <- sd(bm.logit)
            diff.bm.logit <- start.bm.logit - mn.bm.logit
            
            # Then backtransform and everything will stay positive! Just using the AR1 term for these
            # FIX: SEE above comment for where I'm using the AR2, here using the AR2 would give some poor starting values
            # So I'm not comfy doing that (it works by luck in the above for the NS IMHO.)
            tmp.prop.bm <- c(inv.logit(arima.sim(model =list(ar = c(tmp.cor.lag.1)),
                                                 n.start = 1, start.innov = c(diff.bm.logit/tmp.cor.lag.1),
                                                 n = n.yrs.proj,innov = c(0,rnorm(n.yrs.proj-1,0,sd.bm.logit))) + mn.bm.logit))
            sim.Ks[[s]] <- data.frame(Years = 1:n.yrs.proj, sim = i,
                                      Stock = s, troph.cat = tl,
                                      bm.stock = tmp.prop.bm*bm.trophic.Ks[[i]]$bm.tl[bm.trophic.Ks[[i]]$troph.cat==tl])
          } # end the if(count == 1 ||  n.stock.tl != 2)
          # If there are only 2 stocks in a trophic level, then the second stock get the rest of the trophic levels biomass
          
          if(count == 2 & n.stock.tl == 2) 
          {
            sim.Ks[[s]] <-  data.frame(Years = 1:n.yrs.proj, sim=i,
                                       Stock = s, troph.cat = as.numeric(tl),
                                       bm.stock = (1-tmp.prop.bm)*bm.trophic.Ks[[i]]$bm.tl[bm.trophic.Ks[[i]]$troph.cat==tl])
          } # end the case of just 2 stocks
        } # end the stocks loop
      } # end the trophic level loop
      sim.K.stock[[i]] <- do.call("rbind",sim.Ks)
      
    } # end the simulation loop
  }#end of version 1 simulation
  
  if (version == 2){
    #start with LP
    for(i in 1:n.sims) 
    {
      browser()
      # The ecosystem K, using the mean of the ecosystem with the correlation observed of the time series.
      # This starts the time series at the last value of the time series, then moves it to the mean value, bam!!  This will be done for each of these arima sims.
      sim.eco.bm[[i]] <- data.frame(bm = c(arima.sim(model =list(ar = LP.cor$acf[1]),n = n.yrs.proj,n.start=1,start.innov = start.eco.diff/K.cor$acf[1],
                                                     innov = c(0,rnorm(n.yrs.proj-1,0,sd.eco.bm))) + mn.eco.bm),
                                    Years = 1:n.yrs.proj,sim = i) 
      #pacf(sim.eco.bm[[i]]$bm)
      #NO SIGNIFICANT CORRELATIONS HERE - DOES THIS MEAN MODEL IS WORKING PROPERLY?
      
      # So then from my simulated ecosystem I want each FG to get its cut of the biomass, 
      # FIX: I am using the AR2, but I know the start innovation is slightly incorrect, but it make almost no difference for the NS case so I'll stick with it
      # so probably should figure out how to specify that right as it just works by luck here I think, if the difference was larger
      # or correlations different it wouldn't do so well (e.g., it isn't nice for the stock level ones.)
      sim.tl.3.prop.bm <-inv.logit(mn.tl.3.logit + 
                                     arima.sim(model =list(ar = c(tl.3.prop.bm.lag.1,tl.3.prop.bm.lag.2)),n = n.yrs.proj,
                                               n.start =2, start.innov = c(start.tl.3.diff/tl.3.prop.bm.lag.1,start.tl.3.diff/tl.3.prop.bm.lag.1), 
                                               innov = c(0,rnorm(n.yrs.proj-1,0,sd.tl.3.logit))))
      
      bm.sim.3 <- sim.tl.3.prop.bm * sim.eco.bm[[i]]$bm
      # So this is what is left for 3 and 4
      bm.sim.4 <- sim.eco.bm[[i]]$bm - bm.sim.3
      #bm.left.4.5<- sim.eco.bm[[i]]$bm - bm.sim.3
      # So then we use the historical split between 4 and 5 can see 5 gets about 1/3-1-5 of 3
      # so then simulate this split
      #sim.tl.5.4.prop.bm <- inv.logit(mn.tl.4.5.logit + 
      # arima.sim(model =list(ar = tl.4.5.prop.bm.lag.1),n = n.yrs.proj,
      #  n.start =1, start.innov = c(start.tl.4.5.diff/tl.4.5.prop.bm.lag.1), 
      #  innov = c(0,rnorm(n.yrs.proj-1,0,sd.tl.3.logit))))
      # And now TL 5 gets this proportion of the 4 and 5 biomass
      #bm.sim.5 <- bm.left.4.5 * sim.tl.5.4.prop.bm
      # And TL4 gets the rest, and so the ecosystem biomass is a portion of the whole biomass
      #bm.sim.4 <- sim.eco.bm[[i]]$bm - bm.sim.3-bm.sim.5
      
      bm.trophic.Ks[[i]] <- data.frame(Years = rep(1:n.yrs.proj,2), sim =i,
                                       bm.tl = c(bm.sim.3,bm.sim.4),troph.cat = as.factor(sort(rep(c(3,4),n.yrs.proj))),
                                       bm.eco = rep(sim.eco.bm[[i]]$bm,2))
      bm.trophic.Ks[[i]]$prop.bm.tl <- bm.trophic.Ks[[i]]$bm.tl/bm.trophic.Ks[[i]]$bm.eco
      
      # OK, so now we have the trophic level K values simulated in a 'nice' way. Next how do we partition these to the stocks
      # Give each stock a proportion of the K in it's ecosystem based on their historical cuts of the K, and include the time series correlation in that.
      # I'm going to build in correlation to their K time series (this could 100% be fishery induced correlation), could also put in 
      # cross correlation for species with multiple stocks, but for now, let's just do the AR1/2 thing with this for the proportion of the trophic level 
      # biomass each stock gets.
      
      for(tl in troph.levels)
      {
        tl.stocks <- unique(bm.best$Stock[bm.best$troph.cat==tl])
        n.stock.tl <- length(tl.stocks)
        count =0
        for(s in tl.stocks)
        {
          count = count+1
          # Now get the time series for each stock...
          if(count == 1 ||  n.stock.tl != 2)
          {
            tmp.dat <- bm.best[bm.best$Stock ==s,]
            tmp.cor <- pacf(tmp.dat$prop.bm.tl,plot=F) # Get the correlation, use AR1 and AR2 but no more.
            tmp.cor.lag.1 <- tmp.cor$acf[1]
            #tmp.cor.lag.2 <- tmp.cor$acf[2]
            #tmp.beta <- estBetaParams(mean(tmp.dat$prop.bm.tl),sd(tmp.dat$prop.bm.tl)^2)
            # Logit tranform the proportions and do the ARIMA on the logits
            bm.logit <- logit(tmp.dat$prop.bm.stock.tl)
            start.bm.logit <- bm.logit[length(bm.logit)]
            mn.bm.logit <- mean(bm.logit)
            # Use the most recent bm on logit scale as the 'mean' value for the simulation
            mn.bm.logit <- start.bm.logit
            # And the standard deviation
            sd.bm.logit <- sd(bm.logit)
            diff.bm.logit <- start.bm.logit - mn.bm.logit
            
            # Then backtransform and everything will stay positive! Just using the AR1 term for these
            # FIX: SEE above comment for where I'm using the AR2, here using the AR2 would give some poor starting values
            # So I'm not comfy doing that (it works by luck in the above for the NS IMHO.)
            tmp.prop.bm <- c(inv.logit(arima.sim(model =list(ar = c(tmp.cor.lag.1)),
                                                 n.start = 1, start.innov = c(diff.bm.logit/tmp.cor.lag.1),
                                                 n = n.yrs.proj,innov = c(0,rnorm(n.yrs.proj-1,0,sd.bm.logit))) + mn.bm.logit))
            sim.Ks[[s]] <- data.frame(Years = 1:n.yrs.proj, sim = i,
                                      Stock = s, troph.cat = tl,
                                      bm.stock = tmp.prop.bm*bm.trophic.Ks[[i]]$bm.tl[bm.trophic.Ks[[i]]$troph.cat==tl])
          } # end the if(count == 1 ||  n.stock.tl != 2)
          # If there are only 2 stocks in a trophic level, then the second stock get the rest of the trophic levels biomass
          
          if(count == 2 & n.stock.tl == 2) 
          {
            sim.Ks[[s]] <-  data.frame(Years = 1:n.yrs.proj, sim=i,
                                       Stock = s, troph.cat = as.numeric(tl),
                                       bm.stock = (1-tmp.prop.bm)*bm.trophic.Ks[[i]]$bm.tl[bm.trophic.Ks[[i]]$troph.cat==tl])
          } # end the case of just 2 stocks
        } # end the stocks loop
      } # end the trophic level loop
      sim.K.stock[[i]] <- do.call("rbind",sim.Ks)
      
    } # end the simulation loop
  }
  
  sim.K.stocks <- do.call("rbind",sim.K.stock)
  sim.troph.K <- do.call("rbind",bm.trophic.Ks)
  sim.eco.K <- do.call("rbind",sim.eco.bm)
  # Wrap up the K time series for each simulation
  #sim.K.stocks$Species <- substr(sim.K.stocks$Stock,14,100)
}