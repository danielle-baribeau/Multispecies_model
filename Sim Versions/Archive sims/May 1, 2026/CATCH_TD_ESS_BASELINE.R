#in a thesis, chapter 1 is description of model/proof of concept

#paper 1: fishing
#can we support otter trawl fishery?
#paper 2: environment - universal degradation vs. CRIB?
#get trophic structure back to the way it was in the 1970s/1980s?
#what would have to happen to reverse the regime shift?
#no fishing
#TD-BU is best idea - "predator-prey" model!!!

#Next steps: figure out how to compare scenarios
#aggregated figs like DK
#compare different scenarios
#play around with env dev (3% per decade)
#get BU scenario working
#I don't need HCRs! (check GitHub for DK's example)
#integrate historical removals 
#put together an outline of each paper FIRST
#topic sentence of each paragraph I want in the paper
#do intro, methods, results this way


#"If fisheries were going to be expanded by a significant amount, what would happen to the system?"
#Keep removals equal to background historical fishing, and vary one stock at a time
#reinstatement of targeted fishing on the ESS?
#is there capacity to reinstate groundfisheries on the ESS?
#how long would it take to see impacts of fisheries changes on groundfish in the ESS?
#is there scope for multispecies groundfish fisheries to increase?
#fishing as normal in CHPs fishery? YES
#going out directed fishing in CHP - pollock and haddock managed under RP, incidentally catch the rest
#need catch composition of CHP (because this would catch the other species)
#THE PROBLEM
#direct for pollock or haddock - catching other stuff while you do that
#groundfish fishery has to land everything
#means we have the data for what was landed
#catch composition depends on what's in the water and the fishermen knowing what's in the water
#fishers change behaviour, but not gear, so can't tell what they're directing for
#gonna be extreme highs and lows in catch composition - could have different scenarios for incidental catch
#after collapse in the ESS, is the community resilient enough for target fishing on top species?
#CHP fishery
#don't have to care about halibut because this is a longline!
#if we introduced a mobile groundfish fishery, what would happen to the stocks?
#compare to no fishing? Or compare to ecosystem chapter
#make fishery paper (TD or BU) based on whichever is least resilient (unless literature says otherwise)
#taking mroe extreme scenario would be precautionary
#intro
#ESS community has completely restructured - pulled out a ton of biomass
#what is happening now on ESS? Snowcrab fishery would catch cod, etc.
#Could we expand Canadian fisheries using groundfish mobile gear?
#mostly concerned about changes in bycatch, not target (because we'll be managing according to the PP)
#it's a bycatch paper! Kind of :P

#discussion - possible extensions if catchability data was better etc.

#also need a reference point comparison for "bycatch" species
#develop LRPs for the bycatch species
#monitor when the fishery dips these species below; don't stop fishing
#"we can maintain all bycatch species above LRP if exploitation on haddock stays above x"

#second paper - does climate resilience depend on ecosystem structure?
#bottom-up versus top-down of existing model
#which one supports ecosystem
#no additional fishing (adding more fishing just increases unknowns)
#would we expect fundamentally different resilience depending on how it's structured?
#if you hammer larges, you've essentially created a bottom-up community - how does this community react?
#constant penalty based off of CRIB scores?
#make climate change penalty as simple as possible
#put overall penalty on ecosystem biomass, then maybe partition this to stocks based on CRIB later on
#compare TD, BU, TD + decline, BU + decline
#resilience of functional group size

#fishing as normal in silver hake fishery? (only in very specific places, too niche)


#NEXT STEPS:
#turn R off and on again, then re-run everything again
#if it runs, SAVE AS "BASELINE" THEN OPEN NEW FILE TO DO THE REST
#500 sims
#adding exploitation on haddock and pollock, with bycatch of cod (cod is under moratorium)
#mechanism to track removals
#have to fish a given exploitation rate of haddock, which will be associated with a certain proportion of catch of other species
#another multinomial? Exploitation rate across the group
#HCR on haddock and pollock

#not going to add removals to input biomass because survey will already account for this

## Dynamic predator-prey pools!
#this model is only designed for making small perturbations to the ecosystem, not 
#major restructuring of ecosystem
#"stock has pool of resources available to it - whether it uses this pool depends on lambda"

#TO - DO:

#turn lambda sampling into a linear regression (log lambda and biomass)
#need to add in variation (sd - stock-specific?)


#MEDIUMS TOO!
#also try running with a complete reset - Looks good! DONE
#also look at density dependence patterns - try different values
#plot lambda against biomass to see if there are any weird patterns
#try stop forcing lambda to be under 1 if we go above K?
#try with 50 sims
#linear regression on lambda versus biomass with some random sampling? (if they look linear)
#just plot these first for each stockcbn 
#multiply bm.sim.K.L by proportions of mean input time series for now, see if ex is breaking things
#biomass multiplication hack to prevent tiny proportions


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
#max.com.sim <- max(com.bm)
mn.com.sim <- mean(com.bm)
#med.com.sim <- median(com.bm)
#q.com.sim <- c(quantile(com.bm, probs = c(0.75), na.rm = T))
start.com.diff <- start.com.sim - mn.com.sim
#sd.com.sim <- sd(log(com.bm))
sd.com.sim <- sd(com.bm)/sqrt(length(com.bm))

L.com.prop <- bm.best.b[bm.best.b$BFG == "L",]
L.com.prop <- L.com.prop[,c(2,8)]
L.com.prop <- unique(L.com.prop)
L.com.prop <- L.com.prop$prop.eco.bfg
start.L.com.prop <- L.com.prop[length(L.com.prop)]
L.com.prop.cor <- pacf(L.com.prop)
L.com.prop.cor.lag.1 <- L.com.prop.cor$acf[1]
L.com.prop.cor.lag.2 <- L.com.prop.cor$acf[2]

L.com.prop.logit <- logit(L.com.prop)
start.L.com.prop.logit <- L.com.prop.logit[length(L.com.prop.logit)]


#This doesn't equal the sd of L.com.prop
mn.L.sim.logit <- mean(L.com.prop.logit)
#med.L.sim.logit <- median(L.com.prop.logit)
sd.L.sim.logit <- sd(L.com.prop.logit)/sqrt(length(L.com.prop.logit))
#max.L.sim.logit <- max(L.com.prop.logit)
#q.L.sim.logit <- c(quantile(L.com.prop.logit, probs = c(0.75), na.rm = T))

start.L.com.prop.diff.logit <- start.L.com.prop.logit - mn.L.sim.logit

#starting proportions for each stock
L.stocks.data <- bm.best.b |> collapse::fsubset(Year == 2017) |> collapse::fsubset(BFG == "L")
L.stocks <- as.character(unique(L.stocks.data$species))
L.bm.start <- L.stocks.data$bm.bfg[1]

LP.stocks <- as.character(unique(bm.best.m$species[bm.best.m$MFG == "LP"]))


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
M.K.temp <- NULL

n.sims <- 10

env.dev <- 0.99

lambdas <- do.call("rbind", stocks.lst)
#truncate to years we are looking at
lambdas <- lambdas |> collapse::fsubset(year %in% yrs)

#initialize a list to hold results of simulation
res.yr.df <- NULL

mn.prop.stock.fg.bm <- bm.best.b[bm.best.b$Year %in% c(2000:2017),]

#mn.prop.stock.fg.bm$new.prop.bm.stock.bfg <- mn.prop.stock.fg.bm$bm.stock/mn.prop.stock.fg.bm$bm.bfg
mn.stock.fg <- data.frame(stock = c(L.stocks, M.stocks), fg = "K", mn.bm = 0, mn.fg = 0, mn.prop.fg = 0)

#get mean fg biomasses
mn.fg.M <- mn.prop.stock.fg.bm[mn.prop.stock.fg.bm$BFG == "M",]
mn.fg.M <- unique(mn.fg.M$bm.bfg)
mn.fg.M <- mean(mn.fg.M)

mn.fg.L <- mn.prop.stock.fg.bm[mn.prop.stock.fg.bm$BFG == "L",]
mn.fg.L <- unique(mn.fg.L$bm.bfg)
mn.fg.L <- mean(mn.fg.L)

for (m in mn.stock.fg$stock){
  #browser()
  if (m %in% L.stocks) {
    mn.stock.fg$fg[mn.stock.fg$stock == m] <- "L"
    mn.stock.fg$mn.fg[mn.stock.fg$stock == m] <- mn.fg.L/10
  }
  
  if (m %in% M.stocks) {
    mn.stock.fg$fg[mn.stock.fg$stock == m] <- "M"
    mn.stock.fg$mn.fg[mn.stock.fg$stock == m] <- mn.fg.M/10
  }
  
  mn.stock.fg$mn.bm[mn.stock.fg$stock == m] <- mean(lambdas$pop.bm[lambdas$common == m]/10)
  
  mn.stock.fg$mn.prop.fg[mn.stock.fg$stock == m] <- mn.stock.fg$mn.bm[mn.stock.fg$stock == m]/mn.stock.fg$mn.fg[mn.stock.fg$stock == m]
}

p.lambdas <- lambdas
colnames(p.lambdas)[9] <- "stock"
pred.lam <- NULL
s.p.lam <- NULL
for (d in c(L.stocks, M.stocks)){
  s.lam <- p.lambdas[p.lambdas$stock == d,]
  max.bm <- median(s.lam$year.minus.one)
  s.lam$bm.prop <- s.lam$year.minus.one/max.bm
  mod.tmp <- lm(log(lambda) ~ bm.prop,data = s.lam)
  # Grabbing the model p value... this is decision is very arbitrary
  p.val <- summary(mod.tmp)$coefficients[2,4]
  # Making a prediction object for the stock
  pred.tmp <- data.frame(lambda = 0,bm.prop = seq(0,1,by=0.0001),Stock = d)
  
  # Now get the appropriate lambda at each proportion of biomass
  pred.tmp$lambda <- exp(predict(mod.tmp,newdata = pred.tmp))
  pred.tmp$sd <- summary(mod.tmp)$sigma #This is residual standard error!
  r.square <- summary(mod.tmp)$r.square
  pred.tmp$bm.prop <- round(pred.tmp$bm.prop,digits=4)
  
  pred.lam <- rbind(pred.lam, pred.tmp)
  s.p.lam <- rbind(s.p.lam, s.lam)
  
  #min.lam <- min(s.lam$lambda)
  print(paste(d))
  print(paste(summary(s.lam$lambda)))
  #print(paste(r.square))
  #print(paste(min.lam))
}


for (i in 1:n.sims){
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
  
  bm.sim.K.com <- data.frame(bm = c(arima.sim(model =list(ar = com.cor.lag.1),n = n.yrs.proj + 1,n.start=1,start.innov = start.com.diff/com.cor.lag.1,
                                              innov = c(0,rnorm(n.yrs.proj,0,sd.com.sim))) + mn.com.sim)[2:(n.yrs.proj + 1)],
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
  
  #1: this gives you exponential decline from starting x on each original value per year
  for (n in 1:nrow(bm.sim.K.com)){
    if (n == 1) bm.sim.K.com$env.dev[n] <- 1
    if (n == 2) bm.sim.K.com$env.dev[n] <- env.dev
    if (n > 2) bm.sim.K.com$env.dev[n] <- env.dev - ((1-env.dev)*(bm.sim.K.com$Years[n]-2))
    bm.sim.K.com$adj.bm[n] <- bm.sim.K.com$bm[n]*bm.sim.K.com$env.dev[n]
  }#end 1
  
  #2-lag version was starting from wrong place, switched back to 1 lag version
  #looks too variable
  
  #prop.sim.K.L <- rmultinom(1,bm.sim.K.com$bm, prob = median(L.com.prop))  #1: this gives you exponential decline from starting x on each original value per year
  
  
  #2: this gives you x% of each original value in the time series
  #for (n in 1:nrow(bm.sim.K.com)){
  #if (n == 30)browser()
  #bm.sim.K.com$env.dev[n] <- env.dev
  #bm.sim.K.com$adj.bm[n] <- bm.sim.K.com$bm[n]*bm.sim.K.com$env.dev[n]
  #}#end 2
  
  #3:this gives you x% of first original time series value in each subsequent year
  #for (n in 1:nrow(bm.sim.K.com)){
  #if (n == 1) {
  #bm.sim.K.com$env.dev[n] <- 1
  #bm.sim.K.com$adj.bm[n] <- bm.sim.K.com$bm[1]
  #}
  #if (n == 2){
  #bm.sim.K.com$env.dev[n] <- env.dev
  #bm.sim.K.com$adj.bm[n] <- bm.sim.K.com$bm[1]*bm.sim.K.com$env.dev[n]
  #}
  #if (n > 2){
  # bm.sim.K.com$env.dev[n] <- env.dev - ((1-env.dev)*(bm.sim.K.com$Years[n]-2))
  # bm.sim.K.com$adj.bm[n] <- bm.sim.K.com$bm[1]*bm.sim.K.com$env.dev[n]
  #}
  #}#end 3
  #this works, but gets us to 0 at year 21 - maybe with rates of decline lower than 5%?
  #works for 1-3% 
  
  #I don't think these are right
  
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
  
  #prop.sim.K.L <- inv.logit(arima.sim(model = list(ar = c(L.com.prop.cor.lag.1)), n = n.yrs.proj,
  #n.start = 1, start.innov = start.L.com.prop.diff.logit/L.com.prop.cor.lag.1,
  #innov = c(1,rlnorm(n.yrs.proj-1,0,sd.L.sim.logit))*med.L.sim.logit - med.L.sim.logit) + med.L.sim.logit)
  
  #ORIGINAL - NOT STARTING FROM 2017 POINT
  #prop.sim.K.L <-inv.logit(arima.sim(model =list(ar = c(L.com.prop.cor.lag.1,L.com.prop.cor.lag.2)),n = n.yrs.proj,
  #n.start =2, start.innov = c(start.L.com.prop.diff.logit/L.com.prop.cor.lag.1,
  #start.L.com.prop.diff.logit/L.com.prop.cor.lag.1), 
  #innov = c(0,rnorm(n.yrs.proj-1,0,sd.L.sim.logit)))+med.L.sim.logit)
  
  #prop.sim.K.L <- c(rep(min(in.prop.com.l$L.com.prop), n.yrs.proj))
  #browser()
  #this is the K value for the larges
  #bm.sim.K.L <- bm.sim.K.com$bm * c(prop.sim.K.L)
  
  
  #reducing sd, changign to 1st order and switching to mean helps
  prop.sim.K.L <-inv.logit(arima.sim(model =list(ar = c(L.com.prop.cor.lag.1)),n = n.yrs.proj + 1,
                                     n.start =1, start.innov = start.L.com.prop.diff.logit/L.com.prop.cor.lag.1, 
                                     innov = c(0,rnorm(n.yrs.proj,0,sd.L.sim.logit))) + mn.L.sim.logit)[2:(n.yrs.proj + 1)]
  
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
  
  bm.sim.K.L <- bm.sim.K.com$adj.bm * c(prop.sim.K.L)
  
  #bm.sim.K.L <- bm.sim.K.com$bm * c(prop.sim.K.L)
  
  #this is the (stand-by) value for the mediums (will be updated later after L dynamics are calculated)
  #bm.sim.K.M <- bm.sim.K.com$bm * (1 - c(prop.sim.K.L))
  bm.sim.K.M <- bm.sim.K.com$adj.bm * (1 - c(prop.sim.K.L))
  
  for (t in 1:n.yrs.proj){
    #browser()
    #if (t == 2)browser()
    
    #Get stock-specific K values using a multinomial distribution:
    #get proportions of start L K that each stock uses with a multinomial distribution
    
    if (t == 1){
      res.L.stock.df[[i]] <- data.frame(sim = i, year = 2017,
                                        stock = L.stocks, 
                                        fg = L.stocks.data$BFG,
                                        bm.stock = L.stocks.data$bm.stock,
                                        bm.fg = L.stocks.data$bm.bfg,
                                        prop.stock.fg.bm = NA,
                                        bm.com = L.stocks.data$bm.eco,
                                        removals = NA, lambda = NA,stock.K = NA, fg.K = NA, com.K = NA, multinom = NA,
                                        lam.med = NA)
      
      res.M.stock.df[[i]] <- data.frame(sim = i, year = 2017,
                                        stock = M.stocks, 
                                        fg = M.stocks.data$BFG,
                                        bm.stock = M.stocks.data$bm.stock,
                                        bm.fg = M.stocks.data$bm.bfg,
                                        prop.stock.fg.bm = NA,
                                        bm.com = M.stocks.data$bm.eco,
                                        removals = NA, lambda = NA,stock.K = NA, fg.K = NA, com.K = NA, multinom = NA,
                                        lam.med = NA)
      
      rm.bm.sim.K.L <- bm.sim.K.L[1]/10
      
      L.multinom.year <- rmultinom(1, rm.bm.sim.K.L, 
                                   prob = (mn.stock.fg$mn.prop.fg[mn.stock.fg$stock %in% L.stocks]))
      L.multinom.year <- L.multinom.year*10
      
    }
    
    if (t > 1) {
      res.L.bm.last.yr <- res.L.stock.df[[i]][res.L.stock.df[[i]]$year == t + 2015,]
      res.L.bm.last.yr <- res.L.bm.last.yr[,c(1:4,7,11:12)]
      #res.L.bm.last.yr$prop.stock.fg.bm <- res.L.stock.df[[i]]$prop.stock.fg.bm[res.L.stock.df[[i]]$year == t + 2016]
      res.L.bm.last.yr$new.prop.stock.fg.bm <- 0
      res.L.bm.last.yr$stock.K <- res.L.bm.last.yr$stock.K/10
      res.L.bm.last.yr$fg.K <- res.L.bm.last.yr$fg.K/10
      res.L.bm.last.yr$new.prop.stock.fg.bm <- res.L.bm.last.yr$stock.K/res.L.bm.last.yr$fg.K
      #2015 is needed here because you need to fill in the t - 2 K proportion slots (this is because year 1 has both 2017 and 2018)
      
      rm.bm.sim.K.L <- bm.sim.K.L[t]/10
      
      #if (rm.bm.sim.K.L > 1e9){
      #while (rm.bm.sim.K.L > 1e9){
      #rm.bm.sim.K.L/10
      #multicount.l <- multicount + 1
      #} 
      #L.multinom.year <- rmultinom(1, rm.bm.sim.K.L, 
      #                             prob = (mn.stock.fg$mn.prop.fg[mn.stock.fg$stock %in% L.stocks]))
      #L.multinom.year <- L.multinom.year*(10*multicount.l)
      #rm.bm.sim.K.L <- bm.sim.K.L[t]*50
      
      #}
      #if (rm.bm.sim.K.L < 1e9) {L.multinom.year <- rmultinom(1, rm.bm.sim.K.L, 
      # prob = (res.L.bm.last.yr$new.prop.stock.fg.bm))}
      #L.multinom.year <- L.multinom.year/50
      
      L.multinom.year <- rmultinom(1, rm.bm.sim.K.L, 
                                   prob = (res.L.bm.last.yr$new.prop.stock.fg.bm))
      L.multinom.year <- L.multinom.year*10
      
      
    }
    
    #L Ks are now done!
    
    #now onto L lambdas:
    #if (t == 2) browser()
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
                                                 & res.L.stock.df[[i]]$year == t+2016]
      
      # Sort out which of the years are low or high bm
      # I'm using 0.5 as the cut off, other options are valid (0.4 is my fav...)
      
      # Get what our carrying capacity is at this moment
      cur.K <- L.multinom.year[count]
      
      #browser()
      
      #########################################################################
      if (bm.last.yr <= cur.K){
        
        bm.p <- round(bm.last.yr/cur.K, digits = 4)
        #if (t > 1){
        # bm.p <- round(bm.last.yr/res.L.stock.df[[i]]$stock.K[res.L.stock.df[[i]]$stock == s &
        #res.L.stock.df[[i]]$year == t + 2015], digits = 4)
        # }
        lam.mn <- pred.lam$lambda[pred.lam$bm.prop == bm.p & pred.lam$Stock == s]
        #lam.sd <- 0.1
        lam.sd <- pred.lam$sd[pred.lam$bm.prop == bm.p & pred.lam$Stock == s]#/sqrt(length(pred.lam$sd[pred.lam$bm.prop == bm.p & pred.lam$Stock == s])))
        #already SE
        lam.samp <- rlnorm(1,log(lam.mn), lam.sd)
        #while(lam.samp > 10 || lam.samp < 0.1) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
        #if (bm.p < 0.01) while(lam.samp < 1 | lam.samp > 1.5*max(stock.lambdas$lambda)) lam.samp <- rlnorm(1, log(lam.mn), lam.sd)
        #if (bm.p >= 0.01)
        while(lam.samp > 1.5*max(stock.lambdas$lambda) | lam.samp < 0.9*min(stock.lambdas$lambda)) lam.samp <- rlnorm(1, log(lam.mn), lam.sd)
        #while(lam.samp > 1.5*max(stock.lambdas$lambda) | lam.samp < 0.9*min(stock.lambdas$lambda) | (lam.samp < 0.9*min(stock.lambdas$lambda) & lam.samp < 0.1)) lam.samp <- rlnorm(1, log(lam.mn), lam.sd)
        bm.p <- NULL
      }
      #getting lambda using a linear regression
      #lam.reg <- lm(lambda ~ year.minus.one, data = stock.lambdas)
      #intercept <- reg.res$intercept[reg.res$stock == s]
      #slope <- reg.res$slope[reg.res$stock == s]
      
      #if (bm.last.yr < cur.K){
      #lam.med <- (slope*bm.last.yr) + intercept
      #lam.med <- exp(lam.med)
      #lam.sd <- 0.2
      #lam.sd <- sd(log(stock.lambdas$lambda[high.bm.years]),na.rm=T)/2
      
      #lam.samp <- rlnorm(1,log(lam.med),lam.sd)
      #}
      
      # Final one, if we are above the K, we are just doing the high biomass scenario for now.
      # Solution, sample from the lambdas at high biomass, but only take lambdas that are <= 1
      if(bm.last.yr > cur.K) 
        
      {
        #lam.mn <- 1/median(pred.lam$lambda[pred.lam$Stock == s])
        #lam.sd <- (median((pred.lam$sd[pred.lam$Stock == s])/sqrt(length(median((pred.lam$sd[pred.lam$Stock == s]))))))#/2
        #lam.sd <- 0.1
        
        lam.mn <- cur.K/bm.last.yr
        lam.sd <- pred.lam$sd[pred.lam$Stock == s][1]
        #already SE
        
        lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
        print(paste("Above K", s, t, i))
        
        #lam.med <- median(stock.lambdas$lambda[stock.lambdas$lambda < 1])
        #lam.sd <- 0.2
        #lam.sd <- sd(log(stock.lambdas$lambda[high.bm.years]),na.rm=T)/2
        #lam.samp <- rlnorm(1,log(lam.med),lam.sd)
        
        ## add restrictions on lambda by adding || statements here ##
        if(lam.mn > 0.5) while(lam.samp > mean(stock.lambdas$lambda) | lam.samp < ((0.9*(min(stock.lambdas$lambda))))) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
        #if(lam.mn < 0.5) while(lam.samp > (0.5*mean(stock.lambdas$lambda)) | lam.samp < (0.5*(min(stock.lambdas$lambda)))) lam.samp <- rlnorm(1,log(lam.mn*0.8),lam.sd)
        if (lam.mn < 0.5) while(lam.samp > mean(stock.lambdas$lambda) | lam.samp <  ((0.5*(min(stock.lambdas$lambda))))) lam.samp <- rlnorm(1,log(0.8*lam.mn),lam.sd)
        #while(lam.samp > median(stock.lambdas$lambda) || lam.samp < 0.5*min(stock.lambdas$lambda)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
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
        res.L.stock.df[[i]]$fg.K[res.L.stock.df[[i]]$stock == s] <- sum(L.multinom.year)
        res.L.stock.df[[i]]$com.K[res.L.stock.df[[i]]$stock == s] <- bm.sim.K.com$adj.bm[t]
        res.L.stock.df[[i]]$stock.K[res.L.stock.df[[i]]$stock == s] <- cur.K
        res.L.stock.df[[i]]$multinom[res.L.stock.df[[i]]$stock == s] <- L.multinom.year[count]
        res.L.stock.df[[i]]$lam.med[res.L.stock.df[[i]]$stock == s] <- lam.mn
        res.L.stock.df[[i]]$prop.stock.fg.bm[res.L.stock.df[[i]]$stock == s] <- cur.K/sum(L.multinom.year)
        
        res.L.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "L",
                                         bm.stock = tst.res, bm.fg = NA, 
                                         prop.stock.fg.bm = NA,
                                         bm.com = NA,
                                         removals = NA, lambda = NA, stock.K = NA,
                                         fg.K = NA, com.K = NA, multinom = NA,
                                         lam.med = NA)
        
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
                                   res.L.stock.df[[i]]$stock == s] <- sum(L.multinom.year)
        res.L.stock.df[[i]]$com.K[res.L.stock.df[[i]]$year == t + 2016 &
                                    res.L.stock.df[[i]]$stock == s] <- bm.sim.K.com$adj.bm[t]
        res.L.stock.df[[i]]$multinom[res.L.stock.df[[i]]$year == t + 2016 &
                                       res.L.stock.df[[i]]$stock == s] <- L.multinom.year[count]
        res.L.stock.df[[i]]$lam.med[res.L.stock.df[[i]]$year == t + 2016 &
                                      res.L.stock.df[[i]]$stock == s] <- lam.mn
        res.L.stock.df[[i]]$prop.stock.fg.bm[res.L.stock.df[[i]]$year == t + 2016 &
                                               res.L.stock.df[[i]]$stock == s] <-cur.K/sum(L.multinom.year)
        
        res.L.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "L",
                                         bm.stock = tst.res, bm.fg = NA,
                                         prop.stock.fg.bm = NA,
                                         bm.com = NA,
                                         removals = NA, lambda = NA, stock.K = NA,
                                         fg.K = NA, com.K = NA, multinom = NA,
                                         lam.med = NA)
        
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
    #browser()
    L.bm <- sum(res.L.stock.df[[i]]$bm.stock[res.L.stock.df[[i]]$year == t + 2017], na.rm = T)
    res.L.stock.df[[i]]$bm.fg[res.L.stock.df[[i]]$year == t + 2017] <- L.bm
    
    LP.bm <- sum(res.L.stock.df[[i]]$bm.stock[res.L.stock.df[[i]]$stock %in% LP.stocks 
                                              & res.L.stock.df[[i]]$year == t + 2017])
    
    LP.K <- sum(res.L.stock.df[[i]]$stock.K[res.L.stock.df[[i]]$stock %in% LP.stocks 
                                            & res.L.stock.df[[i]]$year == t + 2016])
    LP.bm.prop.K <- LP.bm/LP.K
    
    
    
    #if (t == 1) {
    #res.L.stock.df[[i]]$prop.stock.fg.bm[res.L.stock.df[[i]]$year == 2017] <- 
    #res.L.stock.df[[i]]$stock.K[res.L.stock.df[[i]]$year == 2017]/res.L.stock.df[[i]]$fg.K[res.L.stock.df[[i]]$year == 2017]
    
    #res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$year == 2018] <- 
    #res.M.stock.df[[i]]$stock.K[res.M.stock.df[[i]]$year == 2017]/res.M.stock.df[[i]]$fg.K[res.M.stock.df[[i]]$year == 2017]
    #res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$year == t + 2017] <- 
    #res.M.stock.df[[i]]$stock.K[res.M.stock.df[[i]]$year == t + 2017]/res.M.stock.df[[i]]$fg.K[res.M.stock.df[[i]]$year == t + 2017]
    #}
    
    
    
    #if (t > 1) {
    #res.L.stock.df[[i]]$prop.stock.fg.bm[res.L.stock.df[[i]]$year == t + 2016] <- 
    #res.L.stock.df[[i]]$stock.K[res.L.stock.df[[i]]$year == t + 2016]/res.L.stock.df[[i]]$fg.K[res.L.stock.df[[i]]$year == t + 2016]
    #}
    
    
    #OLD VERSION 
    #L.bm <- res.L.stock.df[[i]] |> group_by(year, fg) |> 
    #summarize(L.bm = sum(bm.stock, na.rm = TRUE), .groups = "keep") |> 
    #arrange(year, fg) |>
    #ungroup()
    
    #LP.res <- res.L.stock.df[[i]] |> collapse::fsubset(stock %in% LP.stocks) |> group_by(year, fg) |> 
    #summarize(LP.bm = sum(bm.stock, na.rm = TRUE), LP.K = sum(stock.K, na.rm = TRUE), .groups = "keep") |> 
    #arrange(year, fg) |>
    # ungroup()
    #if (t == 1) LP.res$prop.bm.K <- 0
    # LP.res$prop.bm.K[LP.res$year == t + 2017] <- LP.res$LP.bm[LP.res$year == t+ 2017]/LP.res$LP.K[LP.res$year == t + 2016]
    
    
    #res.L.stock.df[[i]]$bm.fg[res.L.stock.df[[i]]$year == t + 2017] <- L.bm$L.bm[L.bm$year == t + 2017]
    
    
    
    #if (t == 1) {
    #res.L.stock.df[[i]]$prop.stock.fg.bm[res.L.stock.df[[i]]$year == 2017] <- 
    # res.L.stock.df[[i]]$stock.K[res.L.stock.df[[i]]$year == 2017]/res.L.stock.df[[i]]$fg.K[res.L.stock.df[[i]]$year == 2017]}
    
    #if (t > 1){
    
    # res.L.stock.df[[i]]$prop.stock.fg.bm[res.L.stock.df[[i]]$year == t + 2016] <- 
    # res.L.stock.df[[i]]$stock.K[res.L.stock.df[[i]]$year == t + 2016]/res.L.stock.df[[i]]$fg.K[res.L.stock.df[[i]]$year == t + 2016]}
    
    print(paste("done larges", s, t, i))
    
    
    #this is the ratio between realized biomass of LP and LP K (theoretical size of LP)
    #give ratio to mediums
    M.K.real <- bm.sim.K.M[t]/LP.bm.prop.K
    #M.K.real <- bm.sim.K.M[t]
    
    #if (M.K.real > bm.sim.K.com$bm[t]) {
    # M.K.real <- bm.sim.K.com$bm[t]
    #}
    
    ##########################
    #Old version
    #this is the difference between realized biomass of L and L K (theoretical size of L)
    #if (t == 1) L.K.space <- bm.sim.K.L[t] - start.res.L.df$bm.fg[1]
    #if (t > 1) L.K.space <- bm.sim.K.L[t] - bm.sim.K.L[t - 1]
    
    #adjust M K according to realized growth of Ls
    #M.K.real <- bm.sim.K.M[t] + L.K.space
    #M.K.real <- bm.sim.K.M[t]
    
    #get K values for mediums for YEAR 1
    #get proportions of start M K that each stock uses with a multinomial distribution
    
    if (t == 1){
      rm.bm.sim.K.L <- M.K.real/10
      
      M.multinom.year <- rmultinom(1, rm.bm.sim.K.L, 
                                   prob = (mn.stock.fg$mn.prop.fg[mn.stock.fg$stock %in% M.stocks]))
      M.multinom.year <- M.multinom.year*10
    }
    
    if (t > 1) {
      res.M.bm.last.yr <- res.M.stock.df[[i]][res.M.stock.df[[i]]$year == t + 2015,]
      res.M.bm.last.yr <- res.M.bm.last.yr[,c(1:4,7,11:12)]
      #res.M.bm.last.yr$prop.stock.fg.bm <- res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$year == t + 2016]
      res.M.bm.last.yr$new.prop.stock.fg.bm <- NA
      res.M.bm.last.yr$stock.K <- res.M.bm.last.yr$stock.K/10
      res.M.bm.last.yr$fg.K <- res.M.bm.last.yr$fg.K/10
      res.M.bm.last.yr$new.prop.stock.fg.bm <- res.M.bm.last.yr$stock.K/res.M.bm.last.yr$fg.K
      
      rm.bm.sim.K.M <- M.K.real/10
      #multicount.m <- 0
      
      #if (rm.bm.sim.K.M > 1e9){
      #while (rm.bm.sim.K.M > 1e9){
      #rm.bm.sim.K.M/10
      #multicount.m <- multicount + 1
      #} 
      #M.multinom.year <- rmultinom(1, rm.bm.sim.K.M, 
      #   prob = (mn.stock.fg$mn.prop.fg[mn.stock.fg$stock %in% M.stocks]))
      #M.multinom.year <- M.multinom.year*(10*multicount.m)
      #rm.bm.sim.K.M <- M.K.real*50
      
      #}
      #M.multinom.year <- rmultinom(1, M.K.real*50, 
      #   prob = (res.M.bm.last.yr$new.prop.stock.fg.bm))
      #M.multinom.year <- M.multinom.year/50
      
      M.multinom.year <- rmultinom(1, rm.bm.sim.K.M, 
                                   prob = (res.M.bm.last.yr$new.prop.stock.fg.bm))
      M.multinom.year <- M.multinom.year*10
    }
    
    M.K.used <- sum(M.multinom.year)
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
      cur.K <- M.multinom.year[count]
      #intercept <- reg.res$intercept[reg.res$stock == s]
      #slope <- reg.res$slope[reg.res$stock == s]
      
      if (bm.last.yr <= cur.K){
        bm.p <- round(bm.last.yr/cur.K, digits = 4)
        #if (t > 1){
        #bm.p <- round(bm.last.yr/res.L.stock.df[[i]]$stock.K[res.L.stock.df[[i]]$stock == s &
        #res.L.stock.df[[i]]$year == t + 2015], digits = 4)
        #}
        
        lam.mn <- pred.lam$lambda[pred.lam$bm.prop == bm.p & pred.lam$Stock == s]
        lam.sd <- pred.lam$sd[pred.lam$bm.prop == bm.p & pred.lam$Stock == s]#/sqrt(length(pred.lam$sd[pred.lam$bm.prop == bm.p & pred.lam$Stock == s])))#/2
        #already SE
        #lam.sd <- 0.1
        #lam.sd <- 0.1
        #lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
        
        #lam.med <- (slope*bm.last.yr) + intercept
        #lam.med <- exp(lam.med)
        #lam.sd <- 0.2
        #lam.sd <- sd(log(stock.lambdas$lambda[high.bm.years]),na.rm=T)/2
        
        lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
        #while(lam.samp > 10 || lam.samp < 0.1) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
        #if (bm.p < 0.1) while(lam.samp < 0.8) lam.samp <- rlnorm(1, log(lam.mn), lam.sd)
        #if (bm.p >= 0.1) 
        #if (bm.p < 0.01) while(lam.samp < 1 | lam.samp > 1.5*max(stock.lambdas$lambda)) lam.samp <- rlnorm(1, log(lam.mn), lam.sd)
        #if (bm.p >= 0.01)
        while(lam.samp > 1.5*max(stock.lambdas$lambda) | lam.samp < 0.9*min(stock.lambdas$lambda)) lam.samp <- rlnorm(1, log(lam.mn), lam.sd)
        bm.p <- NULL
      }
      
      # Final one, if we are above the K, we are just doing the high biomass scenario for now.
      # Solution, sample from the lambdas at high biomass, but only take lambdas that are <= 1
      if(bm.last.yr > cur.K) 
      {
        #lam.mn <- 1/median(pred.lam$lambda[pred.lam$Stock == s])
        #lam.sd <- (median((pred.lam$sd[pred.lam$Stock == s])/sqrt(length(median((pred.lam$sd[pred.lam$Stock == s]))))))#/2
        #lam.sd <- 0.1
        
        lam.mn <- cur.K/bm.last.yr
        lam.sd <- pred.lam$sd[pred.lam$Stock == s][1]
        ##already SE
        
        lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
        print(paste("Above K", s, t, i))
        
        #lam.med <- median(stock.lambdas$lambda[stock.lambdas$lambda < 1])
        #lam.sd <- 0.2
        #lam.sd <- sd(log(stock.lambdas$lambda[high.bm.years]),na.rm=T)/2
        #lam.samp <- rlnorm(1,log(lam.med),lam.sd)
        
        ## add restrictions on lambda by adding || statements here ##
        if(lam.mn > 0.5) while(lam.samp > mean(stock.lambdas$lambda) | lam.samp < ((0.9*(min(stock.lambdas$lambda))))) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
        #if(lam.mn < 0.5) while(lam.samp > (0.5*mean(stock.lambdas$lambda)) | lam.samp < (0.5*(min(stock.lambdas$lambda)))) lam.samp <- rlnorm(1,log(lam.mn*0.8),lam.sd)
        if (lam.mn < 0.5) while(lam.samp > mean(stock.lambdas$lambda) | lam.samp <  ((0.5*(min(stock.lambdas$lambda))))) lam.samp <- rlnorm(1,log(0.8*lam.mn),lam.sd)
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
        res.M.stock.df[[i]]$com.K[res.M.stock.df[[i]]$stock == s] <- bm.sim.K.com$adj.bm[t]
        res.M.stock.df[[i]]$stock.K[res.M.stock.df[[i]]$stock == s] <- cur.K
        res.M.stock.df[[i]]$multinom[res.M.stock.df[[i]]$stock == s] <- M.multinom.year[count]
        res.M.stock.df[[i]]$lam.med[res.M.stock.df[[i]]$stock == s] <- lam.mn
        res.M.stock.df[[i]]$fg.K[res.M.stock.df[[i]]$stock == s] <- sum(M.multinom.year)#M.K.real
        res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$stock == s] <- cur.K/sum(M.multinom.year)
        
        res.M.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "M",
                                         bm.stock = tst.res, bm.fg = NA, 
                                         prop.stock.fg.bm = NA,
                                         bm.com = NA,
                                         removals = NA, lambda = NA, stock.K = NA,
                                         fg.K = NA, com.K = NA, multinom = NA,
                                         lam.med = NA)
      }
      
      if (t > 1){
        
        res.M.stock.df[[i]]$removals[res.M.stock.df[[i]]$year == t + 2016 &
                                       res.M.stock.df[[i]]$stock == s] <- removals
        res.M.stock.df[[i]]$lambda[res.M.stock.df[[i]]$year == t + 2016 &
                                     res.M.stock.df[[i]]$stock == s] <- lam.samp
        res.M.stock.df[[i]]$stock.K[res.M.stock.df[[i]]$year == t + 2016 &
                                      res.M.stock.df[[i]]$stock == s] <- cur.K
        res.M.stock.df[[i]]$fg.K[res.M.stock.df[[i]]$year == t + 2016 &
                                   res.M.stock.df[[i]]$stock == s] <- sum(M.multinom.year)#M.K.real
        res.M.stock.df[[i]]$com.K[res.M.stock.df[[i]]$year == t + 2016 &
                                    res.M.stock.df[[i]]$stock == s] <- bm.sim.K.com$adj.bm[t]
        res.M.stock.df[[i]]$multinom[res.M.stock.df[[i]]$stock == s &
                                       res.M.stock.df[[i]]$year == t + 2016] <- M.multinom.year[count]
        res.M.stock.df[[i]]$lam.med[res.M.stock.df[[i]]$stock == s &
                                      res.M.stock.df[[i]]$year == t + 2016] <- lam.mn
        res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$year == t + 2016 &
                                               res.M.stock.df[[i]]$stock == s] <-cur.K/sum(M.multinom.year)
        
        res.M.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "M",
                                         bm.stock = tst.res, bm.fg = NA,
                                         prop.stock.fg.bm = NA,
                                         bm.com = NA,
                                         removals = NA, lambda = NA, stock.K = NA,
                                         fg.K = NA, com.K = NA, multinom = NA,
                                         lam.med = NA)
      }
      
      res.M.stock.df[[i]] <- rbind(res.M.stock.df[[i]], res.M.stock.df.tmp)
      #reset temp filler vector
      res.M.stock.df.tmp <- NULL
      
    }#end stock K loop
    
    #sum up M biomass
    M.bm <- sum(res.M.stock.df[[i]]$bm.stock[res.M.stock.df[[i]]$year == t + 2017], na.rm = T)
    res.M.stock.df[[i]]$bm.fg[res.M.stock.df[[i]]$year == t + 2017] <- M.bm
    
    #old version
    #M.bm <- res.M.stock.df[[i]] |> collapse::fsubset(year == t + 2017) |> group_by(year, fg) |> 
    # summarize(M.bm = sum(bm.stock, na.rm = TRUE), .groups = "keep") |> 
    #arrange(year, fg) |>
    #ungroup()
    
    #res.M.stock.df[[i]]$bm.fg[res.M.stock.df[[i]]$year == t + 2017] <- M.bm$M.bm[M.bm$year == t + 2017]
    
    #if (t == 1) {
    # res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$year == 2017] <- 
    #  res.M.stock.df[[i]]$stock.K[res.M.stock.df[[i]]$year == 2017]/res.M.stock.df[[i]]$fg.K[res.M.stock.df[[i]]$year == 2017]
    
    # res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$year == 2018] <- 
    #  res.M.stock.df[[i]]$stock.K[res.M.stock.df[[i]]$year == 2017]/res.M.stock.df[[i]]$fg.K[res.M.stock.df[[i]]$year == 2017]
    #res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$year == t + 2017] <- 
    #res.M.stock.df[[i]]$stock.K[res.M.stock.df[[i]]$year == t + 2017]/res.M.stock.df[[i]]$fg.K[res.M.stock.df[[i]]$year == t + 2017]
    #}
    
    
    
    #if (t > 1) {
    #res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$year == t + 2016] <- 
    #res.M.stock.df[[i]]$stock.K[res.M.stock.df[[i]]$year == t + 2016]/res.M.stock.df[[i]]$fg.K[res.M.stock.df[[i]]$year == t + 2016]
    #}
    
    
    print("done mediums")
    
    #get community biomass for the year
    new.bm.com <- L.bm + M.bm
    res.L.stock.df[[i]]$bm.com[res.L.stock.df[[i]]$year == t + 2017] <- new.bm.com
    res.M.stock.df[[i]]$bm.com[res.M.stock.df[[i]]$year == t + 2017] <- new.bm.com
    
    #if (t == 2) browser()
    
    #if (t == 2) browser()
    #results for year t
    #res.yr.df.tmp[[t]] <- rbind(res.L.stock.df[[i]][res.L.stock.df[[i]]$year == t + 2017,], 
    #res.M.stock.df[[i]][res.M.stock.df[[i]]$year == t + 2017,])
    
  } #end t loop
  #res.yr.df[[i]] <- do.call("rbind", res.yr.df.tmp)
  
  #res.yr.df[[i]] <- rbind(res.L.stock.df[[i]], res.M.stock.df[[i]])
  
  prop.K.L[[i]] <- data.frame(sim = i, year = 1:30, prop.com.L.K = c(prop.sim.K.L))
  com.sim[[i]] <- data.frame(sim = i, year = 1:30, com.K = c(bm.sim.K.com$adj.bm))
  #M.K.temp[[i]] <- data.frame(sim = i, year = 2018:2047, temp.M.K = c(bm.sim.K.M))
  L.bm <- NULL
  LP.bm. <- NULL
  M.bm <- NULL
}
res.L.stock <- do.call("rbind", res.L.stock.df)
res.M.stock <- do.call("rbind", res.M.stock.df)

prop.K.L <- do.call("rbind", prop.K.L)
prop.K.L$year <- prop.K.L$year + 2016

com.sim <- do.call("rbind", com.sim)
com.sim$year <- com.sim$year + 2016

#M.K.temp <- do.call("rbind", M.K.temp)

res.sim.df <- rbind(res.L.stock, res.M.stock)
#get everything into one data frame

per.dec <- NULL
for (d in c(L.stocks, M.stocks)){
  s.res <- na.omit(res.sim.df[res.sim.df$stock == d,])
  #print(paste(d))
  #print(paste(summary(s.res$lambda)))
  for (x in 1:n.sims){
    per.dec.tmp <- data.frame(stock = "N", sim = 0, bm.start = 0, bm.end = 0, per.dec = 0)
    s.sim <- s.res[s.res$sim == x,]
    s.sim.start <- s.sim$bm.stock[s.sim$year == 2017]
    s.sim.end <- s.sim$bm.stock[s.sim$year == 2046]
    s.sim.diff <- s.sim.start - s.sim.end
    s.sim.per <- (s.sim.diff/s.sim.start)*100
    s.sim.per <- s.sim.per*-1
    per.dec.tmp$stock <- d
    per.dec.tmp$sim <- x
    per.dec.tmp$bm.start <- s.sim.start
    per.dec.tmp$bm.end <- s.sim.end
    per.dec.tmp$per.dec <- s.sim.per
    per.dec <- rbind(per.dec, per.dec.tmp)
    per.dec.tmp <- data.frame(stock = "N", sim = 0, bm.start = 0, bm.end = 0, per.dec = 0)
  }
  print(paste(d))
  print(paste(summary(per.dec$per.dec[per.dec$stock == d])))
}

########################################################

########################################################
no.fish <- res.sim.df
fishL <-res.sim.df

q.no.fish <- no.fish[no.fish$fg == "M",]|>  collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(bm.fg,probs=c(0.25),na.rm=T),
                                                                                              med = median(bm.fg,na.rm=T),
                                                                                              U.50 = quantile(bm.fg,probs=c(0.75),na.rm=T))
q.fishL <- fishL[fishL$fg == "M",] |>  collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(bm.fg,probs=c(0.25),na.rm=T),
                                                                                         med = median(bm.fg,na.rm=T),
                                                                                         U.50 = quantile(bm.fg,probs=c(0.75),na.rm=T))

ggplot(q.no.fish) +geom_line(aes(x=year, y=L.50), colour = "grey") +
  geom_line(aes(x=year, y=med)) +
  geom_line(aes(x=year, y=U.50), colour = "grey") +
  geom_line(data = q.fishL, aes(x=year, y=L.50), colour = "lightblue") +
  geom_line(data = q.fishL, aes(x=year, y=med), colour = "darkblue") +
  geom_line(data = q.fishL, aes(x=year, y=U.50), colour = "lightblue") +
  labs(y = "Biomass") +
  guides(colour = guide_legend(nrow = 5))


ggplot(na.omit(res.sim.df[res.sim.df$sim == 2,])) + geom_line(aes(x = year, y = prop.stock.fg.bm, group = sim)) +
  facet_wrap(~stock, scale = "free", ncol = 5)

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

p.lambdas <- lambdas
colnames(p.lambdas)[colnames(p.lambdas) == "common"] <- "stock"
p.lambdas$mn.in <- rep(0, nrow(p.lambdas))

mn.plot <- NULL

for (i in mn.in$stock){
  count <- count + 1
  stock.in <- p.lambdas$pop.bm[p.lambdas$stock == i]
  mn.stock.in <- median(stock.in)
  mn.in$mean[count] <- mn.stock.in
}


for (i in 1:nrow(p.lambdas)){
  #browser()
  stock.p <- p.lambdas$stock[i]
  p.lambdas$mn.in[i] <- mn.in$mean[mn.in$stock == stock.p]
}

L.bm.hist <- data.frame(year = 2000:2017, L.bm = com.bm*L.com.prop)

#[res.sim.df$fg == "L" & res.sim.df$sim == 2,]
ggplot(na.omit(res.sim.df)) +
  geom_line(aes(x = year, y = fg.K/1000000, group = sim), colour = "lightblue") +
  geom_line(aes(x = year, y = bm.fg/1000000, group = sim)) +
  geom_line(data = L.bm.hist, aes(x = year, y = L.bm/1000000), 
            colour = "blue") +
  #geom_line(data = p.lambdas[p.lambdas$FG %in% c("LP", "LB"),], aes(x = year, y = mn.in/1000000, group = stock)) +
  geom_vline(xintercept = 2017, linetype = 2, colour = "blue") +
  geom_hline(yintercept = median(L.bm.hist$L.bm)/1000000, colour = "blue") +
  #geom_line(aes(x = year, y = fg.K/1000000, group = sim), colour = "purple") +
  labs(x = "Year", y = "Biomass (millions of kg)") +
  theme(axis.text.x = element_text(size=14, angle=45),
        axis.text.y = element_text(size=14))

#[ & res.sim.df$sim == 11,]
ggplot(na.omit(res.sim.df[res.sim.df$fg == "L" & res.sim.df$year != 2017,])) +
  geom_line(aes(x = year, y = stock.K/1000000, group = sim), colour = "lightblue") +
  #geom_line(aes(x = year, y = bm.stock/1000000, group = sim)) +
  geom_line(data = p.lambdas[p.lambdas$FG %in% c("LP", "LB"),], aes(x = year, y = pop.bm/1000000, group = stock), 
            colour = "blue")+
  geom_line(data = p.lambdas[p.lambdas$FG %in% c("LP", "LB"),], aes(x = year, y = mn.in/1000000, group = stock)) +
  geom_vline(xintercept = 2017, linetype = 2, colour = "blue") +  facet_wrap(~stock, ncol = 4, scale = "free") +  #scale_y_log10(name="Biomass") + 
  #geom_line(aes(x = year, y = fg.K/1000000, group = sim), colour = "purple") +
  labs(x = "Year", y = "Biomass (millions of kg)") +
  theme(axis.text.x = element_text(size=14, angle=45),
        axis.text.y = element_text(size=14)) +
  guides(colour = guide_legend(ncol = 5))

ggplot(na.omit(res.sim.df[res.sim.df$fg == "M"& res.sim.df$year != 2017 ,])) +
  geom_line(aes(x = year, y = stock.K/1000000, group = sim), colour = "lightblue") +
  #geom_line(aes(x = year, y = bm.stock/1000000, group = sim)) +
  geom_line(data = p.lambdas[p.lambdas$FG %in% c("MP", "MBZ"),], aes(x = year, y = pop.bm/1000000, group = stock), 
            colour = "blue")+
  geom_line(data = p.lambdas[p.lambdas$FG %in% c("MP", "MBZ"),], aes(x = year, y = mn.in/1000000, group = stock)) +
  geom_vline(xintercept = 2017, linetype = 2, colour = "blue") +  facet_wrap(~stock, ncol = 4, scale = "free") +  #scale_y_log10(name="Biomass") + 
  labs(x = "Year", y = "Biomass (millions of kg)") +
  theme(axis.text.x = element_text(size=14, angle=45),
        axis.text.y = element_text(size=14)) +
  guides(colour = guide_legend(ncol = 5))

ggplot(quants[quants$fg == "L",]) +
  #geom_line(data = na.omit(res.sim.df[res.sim.df$fg == "L",]), aes(x = year, y = stock.K/1000000, group = sim), colour = "lightblue") +
  geom_line(aes(x=year, y=L.50/1000000, group=stock,), colour = "grey") +
  geom_line(aes(x = year, y = U.50/1000000, group = stock), colour = "grey") +
  geom_line(aes(x = year, y = med/1000000, group = stock)) +
  geom_line(data = p.lambdas[p.lambdas$FG %in% c("LP", "LB"),], aes(x = year, y = pop.bm/1000000, group = stock), 
            colour = "blue")+
  #geom_line(data = na.omit(res.sim.df[res.sim.df$fg == "L",]), aes(x = year, y = stock.K/1000000, group = sim),
  #colour = "lightblue") +
  geom_line(data = p.lambdas[p.lambdas$FG %in% c("LP", "LB"),], aes(x = year, y = mn.in/1000000, group = stock)) +
  geom_vline(xintercept = 2017, linetype = 2, colour = "blue") +  facet_wrap(~stock, ncol = 4, scale = "free") +  #scale_y_log10(name="Biomass") + 
  labs(x = "Year", y = "Biomass (millions of kg)") +
  theme(axis.text.x = element_text(size=14, angle=45),
        axis.text.y = element_text(size=14)) +
  guides(colour = guide_legend(ncol = 5))

ggplot(quants[quants$fg == "M",]) +
  #geom_line(data = na.omit(res.sim.df[res.sim.df$fg == "M",]), aes(x = year, y = stock.K/1000000, group = sim), colour = "lightblue") +
  geom_line(aes(x=year, y=L.50/1000000, group=stock,), colour = "grey") +
  geom_line(aes(x = year, y = U.50/1000000, group = stock), colour = "grey") +
  geom_line(aes(x = year, y = med/1000000, group = stock)) +
  geom_line(data = p.lambdas[p.lambdas$FG %in% c("MP", "MBZ"),], aes(x = year, y = pop.bm/1000000, group = stock), 
            colour = "blue")+
  #geom_line(data = na.omit(res.sim.df[res.sim.df$fg == "M",]), aes(x = year, y = stock.K/1000000, group = sim),
  #colour = "lightblue") +
  geom_line(data = p.lambdas[p.lambdas$FG %in% c("MP", "MBZ"),], aes(x = year, y = mn.in/1000000, group = stock)) +
  geom_vline(xintercept = 2017, linetype = 2, colour = "blue") +
  facet_wrap(~stock, ncol = 4, scale = "free") +  #scale_y_log10(name="Biomass") + 
  labs(x = "Year", y = "Biomass (millions of kg)") +
  theme(axis.text.x = element_text(size=14, angle=45),
        axis.text.y = element_text(size=14)) +
  guides(colour = guide_legend(nrow = 5))


ggplot(na.omit(res.sim.df[res.sim.df$stock %in% LP.stocks,])) + geom_point(aes(x = year, y = (bm.stock/1000000)/(stock.K/1000000), group = sim), size = 0.5)+
  facet_wrap(~stock, ncol = 3, scale = "free") +
  geom_hline(yintercept = 1, colour = "blue") +
  labs(x = "Year", y = "Biomass/K (millions of kg)")

ggplot(na.omit(res.sim.df[res.sim.df$stock %in% M.stocks,])) + geom_point(aes(x = year, y = (bm.stock/1000000)/(stock.K/1000000), group = sim), size = 0.5)+
  facet_wrap(~stock, ncol = 4, scale = "free") +
  geom_hline(yintercept = 1, colour = "blue") +
  labs(x = "Year", y = "Biomass/K (millions of kg)")


LP.res <- res.sim.df[res.sim.df$stock %in% LP.stocks,]
LP.res$prop.stock.K.bm <- LP.res$bm.stock/LP.res$stock.K

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

#how close are mediums to their Ks?
quants <- res.sim.df |>  collapse::fgroup_by(year,stock,fg) |> collapse::fsummarize(L.50 = quantile(bm.stock,probs=c(0.25),na.rm=T),
                                                                                    med = median(bm.stock,na.rm=T),
                                                                                    U.50 = quantile(bm.stock,probs=c(0.75),na.rm=T))

quants.k <- res.sim.df |>  collapse::fgroup_by(year,stock,fg) |> collapse::fsummarize(L.50 = quantile(stock.K,probs=c(0.25),na.rm=T),
                                                                                      med = median(stock.K,na.rm=T),
                                                                                      U.50 = quantile(stock.K,probs=c(0.75),na.rm=T))


ggplot(na.omit(quants[quants$fg == "M",])) +geom_line(aes(x=year, y=L.50/1000000, group=stock,), colour = "grey") +
  geom_line(aes(x = year, y = U.50/1000000, group = stock), colour = "grey") +
  geom_line(aes(x = year, y = med/1000000, group = stock)) +
  geom_line(data = na.omit(quants.k[quants.k$fg == "M",]), aes(x=year, y=L.50/1000000, group=stock,), colour = "lightblue") +
  geom_line(data = na.omit(quants.k[quants.k$fg == "M",]), aes(x = year, y = U.50/1000000, group = stock), colour = "lightblue") +
  geom_line(data = na.omit(quants.k[quants.k$fg == "M",]), aes(x = year, y = med/1000000, group = stock), colour = "blue") +
  facet_wrap(~stock, ncol = 4, scale = "free")

quants.l <- res.sim.df |>  collapse::fgroup_by(year,stock,fg) |> collapse::fsummarize(L.50 = quantile(lambda,probs=c(0.25),na.rm=T),
                                                                                      med = median(lambda,na.rm=T),
                                                                                      U.50 = quantile(lambda,probs=c(0.75),na.rm=T))


ggplot(na.omit(quants[quants$fg == "M",])) +#geom_line(aes(x=year, y=L.50/1000000, group=stock,), colour = "grey") +
  #geom_line(aes(x = year, y = U.50/1000000, group = stock), colour = "grey") +
  geom_line(aes(x = year, y = med/1000000, group = stock)) +
  #geom_line(data = na.omit(quants.l[quants.l$fg == "M",]), aes(x=year, y=L.50, group=stock,), colour = "lightblue") +
  #geom_line(data = na.omit(quants.l[quants.l$fg == "M",]), aes(x = year, y = U.50, group = stock), colour = "lightblue") +
  geom_point(data = na.omit(quants.l[quants.l$fg == "M",]), aes(x = year, y = med, group = stock), colour = "blue") +
  geom_line(data = p.lambdas[p.lambdas$FG %in% c("MP", "MBZ"),], aes(x = year, y = mn.in/1000000, group = stock), colour = "grey") +
  facet_wrap(~stock, ncol = 4, scale = "free") + scale_y_log10(name = "lam and biomass")

ggplot(na.omit(res.sim.df[res.sim.df$fg == "M",])) + geom_point(aes(x = bm.stock/1000000, y = lambda, group = sim, colour = sim)) +
  facet_wrap(~stock, ncol = 4, scale = "free")

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
stock.eco.K$fg.in.com.K <- stock.eco.K$fg.K/stock.eco.K$com.K

fg.res.df <- res.sim.df[,c(1:2,3:4,6,8,12,13)]
fg.res.df <- unique(fg.res.df)
fg.res.df <- na.omit(fg.res.df)

#ggplot(M.K.temp) +  geom_line(aes(x = year, y = temp.M.K, group = sim), colour = "darkblue")

ggplot(fg.res.df) +   geom_line(data=fg.res.df, aes(x = year, y = com.K/1000000, group = sim)) +
  geom_line(aes(x=year, y=fg.K/1000000, group=sim), colour = "gold") +
  #geom_line(data = bm.best.p, aes(x = Year, y = prop.bm.stock.bfg, group = sim), 
  #colour = "blue")+
  geom_vline(xintercept = 2017, linetype = 2, colour = "blue") +
  #geom_line(data = M.K.temp, aes(x = year, y = temp.M.K/1000000, group = sim), colour = "darkblue") +
  geom_hline(yintercept = max(fg.res.df$com.K/1000000), linetype = 2, colour = "blue") +
  #geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~fg) +  scale_y_log10(name="K (millions)") #+ 
#guides(colour = guide_legend(nrow = 5))


ggplot(fg.res.df) +   geom_line(data=fg.res.df, aes(x = year, y = bm.com/1000000, group = sim)) +
  geom_line(aes(x=year, y=bm.fg/1000000, group=sim), colour = "gold") +
  #geom_line(data = bm.best.p, aes(x = Year, y = prop.bm.stock.bfg, group = sim), 
  #colour = "blue")+
  geom_vline(xintercept = 2017, linetype = 2, colour = "blue") +
  geom_hline(yintercept = max(fg.res.df$bm.com/1000000), linetype = 2, colour = "blue") +
  #geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~fg) +  scale_y_log10(name="Biomass (millions)") #+ 
#guides(colour = guide_legend(nrow = 5))

sim.L.K <- res.sim.df[res.sim.df$fg == "M",]
sim.L.K <- na.omit(sim.L.K)

ggplot(sim.L.K) + geom_line(aes(x=year, y=fg.K, group=sim, colour = sim)) +
  geom_line(data = bm.best.p[bm.best.p$BFG == "M",], aes(x = Year, y = bm.bfg, group = sim), 
            colour = "blue")+
  geom_vline(xintercept = 2017, linetype = 2, colour = "blue") +
  #geom_hline(yintercept = 0, linetype = 2)
  scale_y_log10(name="Med K") 

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

for (d in c(L.stocks, M.stocks)){
  s.res <- na.omit(res.sim.df[res.sim.df$stock == d,])
  min.bm <- min(s.res$bm.stock)
  min.K <- min(s.res$stock.K)
  print(paste(d))
  #print(paste("min bm:", min.bm))
  #print(paste("min.K:", min.K))
  print(paste(summary(s.res$lambda)))
}