#In this version, we are setting up environmental degradation based on CRIB scores of each species
#In CRIB, higher score = more vulnerable
#general process for CRIB penalties:
#Everything above 0.5 will get worse, everything below gets better
#1. Rescale CRIB to go from 0 to 1 (logit transform, > 0.5 = +, < 0.5 = -, 0.5 = 0)
# we had talked about multiplying rescaled CRIB scores by -1 to get right magnitude, but this doesn't
#seem to be necessary
#3. Make overall size of community RP (resource pool) decline (choose percentage, p)
#4. Partition community RP decline among stocks based on CRIB scores
#Adjusted RP = RP - ((p/100)*RP*CRIB)


#what inputs do we need?
#community biomass time series + pacf() on this time series

#proportion of community biomass taken up by larges time series + pacf() (logit-transformed)

#for all these, need starting value, mean, start-mean and sd
com.bm <- eco.tot.bm.best$bm.eco  
com.cor <- pacf(com.bm)
com.cor.lag.1 <- com.cor$acf[1]

start.com.sim <- com.bm[length(com.bm)]
#max.com.sim <- max(com.bm)
mn.com.sim <- mean(com.bm)
med.com.sim <- median(com.bm)
#q.com.sim <- c(quantile(com.bm, probs = c(0.75), na.rm = T))
start.com.diff <- start.com.sim - med.com.sim
#se.com.sim <- sd(log(com.bm))
se.com.sim <- sd(com.bm)/sqrt(length(com.bm))

L.com.prop <- bm.best.b[bm.best.b$BFG == "L",]
L.com.prop <- L.com.prop[,c("Year","prop.eco.bfg")]
L.com.prop <- unique(L.com.prop)
L.com.prop <- L.com.prop$prop.eco.bfg
start.L.com.prop <- L.com.prop[length(L.com.prop)]
L.com.prop.cor <- pacf(L.com.prop)
L.com.prop.cor.lag.1 <- L.com.prop.cor$acf[1]

L.com.prop.logit <- logit(L.com.prop)
start.L.com.prop.logit <- L.com.prop.logit[length(L.com.prop.logit)]

mn.L.sim.logit <- mean(L.com.prop.logit)
med.L.sim.logit <- median(L.com.prop.logit)
se.L.sim.logit <- sd(L.com.prop.logit)/sqrt(length(L.com.prop.logit))
#max.L.sim.logit <- max(L.com.prop.logit)
#q.L.sim.logit <- c(quantile(L.com.prop.logit, probs = c(0.75), na.rm = T))

start.L.com.prop.diff.logit <- start.L.com.prop.logit - med.L.sim.logit

#pulling 2017 data for all stocks
L.stocks.data <- bm.best.b |> collapse::fsubset(Year == 2017) |> collapse::fsubset(BFG == "L")
L.stocks <- as.character(unique(L.stocks.data$species))
#L.bm.start <- L.stocks.data$bm.bfg[1]

LP.stocks <- as.character(unique(bm.best.m$species[bm.best.m$MFG == "LP"]))


M.stocks.data <- bm.best.b |> collapse::fsubset(Year == 2017) |> collapse::fsubset(BFG == "M")
M.stocks <- unique(M.stocks.data$species)


#initialize starting objects
n.yrs.proj <- 100
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

n.sims <- 50

#this now means something different - max percentage you want stock pool to decline by
env.dev <- 10

lambdas <- do.call("rbind", stocks.lst)
#truncate to years we are looking at
lambdas <- lambdas |> collapse::fsubset(year %in% yrs)

#mn.prop.stock.fg.bm$new.prop.bm.stock.bfg <- mn.prop.stock.fg.bm$bm.stock/mn.prop.stock.fg.bm$bm.bfg
mn.stock.fg <- data.frame(stock = c(L.stocks, M.stocks), mn.bm = 0, mn.fg = 0, mn.prop.fg = 0)

#get mean fg biomasses
#change this to be sum of all median stock biomasses to get proportions to sum to 1
#mn.fg.M <- mn.prop.stock.fg.bm[mn.prop.stock.fg.bm$BFG == "M",]
#mn.fg.M <- unique(mn.fg.M$bm.bfg)
#mn.fg.M <- median(mn.fg.M)

#mn.fg.L <- mn.prop.stock.fg.bm[mn.prop.stock.fg.bm$BFG == "L",]
#mn.fg.L <- unique(mn.fg.L$bm.bfg)
#mn.fg.L <- median(mn.fg.L)
#NEW: get CRIB scores set up and associated with each stock
load("C:/Users/BARIBEAUD/Desktop/GitHub/Multispecies_model/ESS_Test_Run/Sim Versions/Current sims/CRIB_ESS_DBaribeau.RData")
#this gives you CRIB scores for all species except longfin hake in individual grid cells across ESS
#for now, going to average across grid cells (check DB notes to see if this was OK approach)
#NOTE: only have CRIB scores for Acadian redfish
#NOTE 2: NO CRIB for longfin hake


crib <- data.frame(common = c(L.stocks, M.stocks), latin = c("Gadus morhua", "Melanogrammus aeglefinus", "Urophycis tenuis",
                                                             "Pollachius virens", "Hippoglossoides platessoides", "Anarhichas lupus",
                                                             "Amblyraja radiata", "Lophius americanus", "Merluccius bilinearis",
                                                             "Sebastes fasciatus", "Glyptocephalus cynoglossus",
                                                             "Phycis chesteri", "Malacoraja senta", "Myoxocephalus octodecemspinosus",
                                                             "Hemitripterus americanus"),
                   med.crib = 0, med.sd = 0)

#get median CRIB score for each species
for (l in crib$latin) {
  #just putting in placeholder for longfin hake for now; need to discuss with H and D
  if (l == "Phycis chesteri") {
    crib$med.crib[crib$latin == l] <- NA
  }
  
  if (l != "Phycis chesteri"){
    vdata.tmp <- vdata[vdata$spname == l,]
    hist(vdata.tmp$vuln, main = paste((l)))
    #everyone looks right-skewed except for longhorn sculpin and sea raven (left skewed)
    #going to use lognormal sampling, same approach as lambda sampling
    #also going to use median instead of mean CRIB
    crib$med.crib[crib$latin == l] <- median(vdata.tmp$vuln)
    crib$med.sd[crib$latin == l] <- median(vdata.tmp$vulnsd)
  }
}#end of med CRIB loop
#checked cod and redfish; it works!

#make longfin hake crib the average of all other cribs
crib$med.crib[crib$common == "Longfin hake"] <- median(crib$med.crib, na.rm = T)



for (m in mn.stock.fg$stock){
  #browser()
  #mn.stock.fg$mn.bm[mn.stock.fg$stock == m] <- median(lambdas$pop.bm[lambdas$common == m]/10)
  mn.stock.fg$mn.bm[mn.stock.fg$stock == m] <- median(lambdas$full.bm[lambdas$common == m]/10)
}

mn.fg.L <- sum(mn.stock.fg$mn.bm[mn.stock.fg$stock %in% L.stocks])
mn.fg.M <- sum(mn.stock.fg$mn.bm[mn.stock.fg$stock %in% M.stocks])


for (n in mn.stock.fg$stock){
  mn.stock.fg$mn.fg[mn.stock.fg$stock == n] <- ifelse (n %in% L.stocks, mn.fg.L, mn.fg.M)
}

mn.stock.fg$mn.prop.fg <- mn.stock.fg$mn.bm/mn.stock.fg$mn.fg
p.lambdas <- lambdas
colnames(p.lambdas)[colnames(p.lambdas) == "common"] <- "stock"



for (i in 1:n.sims){
  #browser()
  #STARTING WITH YEAR 1 OF PROJECTION:
  #go through entire first year (lambdas, K-space), then continue with other years
  #NEW: sample CRIB values for the run
  crib.sim <- data.frame(common = crib$common, latin = crib$latin, samp.crib = crib$med.crib, logit.crib = 0)
  #for (l in crib.sim$latin){
    #going to set SD to 0.1 for now, until we discuss sampling issue
    #chose lognormal because of skewedness in data
    
    #crib.sim$samp.crib[crib.sim$latin == l] <- rlnorm(1, log(crib$med.crib[crib$latin == l]), 0.1)
    
    
    
    #tried this approach, but some of SDs were too big
    #TO DO: check this with 100 sims
    #crib.sim$samp.crib[crib.sim$latin == l] <- rlnorm(1, log(crib$med.crib[crib$latin == l]), crib$med.sd[crib$latin == l]/2)
    #while (crib.sim$samp.crib[crib.sim$latin == l] > 1 || crib.sim$samp.crib[crib.sim$latin == l] < 0) crib.sim$samp.crib[crib.sim$latin == l] <- rlnorm(1, log(crib$med.crib[crib$latin == l]), crib$med.sd[crib$latin == l]/2)
    
    #crib.sim$samp.crib[crib.sim$latin == l] <- rlnorm(1, log(crib$avg.crib[crib$latin == l]), crib$avg.sd[crib$latin == l])
    #crib.sim$samp.crib[crib.sim$latin == l] <- rnorm(1, crib$avg.crib[crib$latin == l], crib$avg.sd[crib$latin == l])
  #}#end of CRIB sampling loop
  
  #rescale CRIB - 2 steps
  #naturally distributed around 0.5, so can logit right away and not worry about rescaling
  #1. rescale by dividing by max value
  #crib.sim$resc.crib <- crib.sim$samp.crib/max(crib.sim$samp.crib)

  #get highest CRIB to be 0.99 so that logit works
  #for (p in 1:nrow(crib.sim)){
    #if (crib.sim$resc.crib[p] == 1) crib.sim$resc.crib[p] <- 0.99
  #}
  
  #2. rescale with logit transformation
  #centers distribution around 0; lets you look at winners (> 0.5) and losers (< 0.5)
  crib.sim$logit.crib <- logit(crib.sim$samp.crib)
  
  
  bm.sim.K.com <- data.frame(bm = c(arima.sim(model =list(ar = com.cor.lag.1),n = n.yrs.proj + 1,n.start=1,start.innov = start.com.diff/com.cor.lag.1,
                                              innov = c(0,rnorm(n.yrs.proj,0,se.com.sim))) + med.com.sim)[2:(n.yrs.proj + 1)],
                             Years = 1:n.yrs.proj,sim = i)
 
  #1: this gives you linear decline from starting x on each original value per year
  #bm.sim.K.com$env.dev <- 0
  #bm.sim.K.com$bm <- 0
  #for (x in 1:nrow(bm.sim.K.com)){
    #if (x == 1) bm.sim.K.com$env.dev[x] <- 1
    #if (x == 2) bm.sim.K.com$env.dev[x] <- env.dev
    #if (x > 2) bm.sim.K.com$env.dev[x] <- env.dev - ((1-env.dev)*(bm.sim.K.com$Years[x]-2))
    #bm.sim.K.com$bm[x] <- bm.sim.K.com$bm[x]*bm.sim.K.com$env.dev[x]
  #}#end 1
  
 
  #reducing sd, changing to 1st order and switching to mean helps
  prop.sim.K.L <-inv.logit(arima.sim(model =list(ar = c(L.com.prop.cor.lag.1)),n = n.yrs.proj + 1,
                                     n.start =1, start.innov = start.L.com.prop.diff.logit/L.com.prop.cor.lag.1, 
                                     innov = c(0,rnorm(n.yrs.proj,0,se.L.sim.logit))) + med.L.sim.logit)[2:(n.yrs.proj + 1)]
  
  
  bm.sim.K.L <- bm.sim.K.com$bm * c(prop.sim.K.L)
  
  #bm.sim.K.L <- bm.sim.K.com$bm * c(prop.sim.K.L)
  
  #this is the (stand-by) value for the mediums (will be updated later after L dynamics are calculated)
  #bm.sim.K.M <- bm.sim.K.com$bm * (1 - c(prop.sim.K.L))
  bm.sim.K.M <- bm.sim.K.com$bm * (1 - c(prop.sim.K.L))
  
  for (t in 1:n.yrs.proj){
    #browser()
    #if (t == 2)browser()
    
    #Get stock-specific K values using a multinomial distribution:
    #get proportions of start L K that each stock uses with a multinomial distribution
    
    if (t == 1){
      res.L.stock.df[[i]] <- data.frame(sim = i, year = L.stocks.data$Year,
                                        stock = L.stocks, 
                                        fg = L.stocks.data$BFG,
                                        bm.stock = L.stocks.data$bm.stock,
                                        bm.fg = L.stocks.data$bm.bfg,
                                        prop.stock.fg.bm = NA,
                                        bm.com = L.stocks.data$bm.eco,
                                        removals = NA, lambda = NA,stock.K = NA, fg.K = NA, com.K = NA, multinom = NA,
                                        lam.med = NA)
      
      res.M.stock.df[[i]] <- data.frame(sim = i, year = M.stocks.data$Year,
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
      res.L.bm.last.yr <- res.L.bm.last.yr[,c("sim", "year", "stock", "fg", "prop.stock.fg.bm", "stock.K", "fg.K")]
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
    #NEW: adjust stock Ks according to their CRIB scores
    #Adjusted RP = RP - ((p/100)*RP*CRIB)
    
    L.Ks <- data.frame(stock = L.stocks, org.crib = crib$med.crib[crib$common %in% L.stocks],
                       logit.crib = crib.sim$logit.crib[crib.sim$common %in% L.stocks], 
                       K = L.multinom.year, adj.K = 0)
    
    for (l in L.Ks$stock){
      crib.score <- L.Ks$org.crib[L.Ks$stock == l]
      K <- L.Ks$K[L.Ks$stock == l]
      L.Ks$adj.K[L.Ks$stock == l] <- K - ((env.dev/100)*K*crib.score)
    }
    
    L.Ks$per.rem <- L.Ks$adj.K/L.Ks$K
    #sum(L.Ks$adj.K)/sum(L.multinom.year)
    
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
      stock.lambdas <- na.omit(lambdas) |> collapse::fsubset(common == s)
      #bm.ts.stock <- bm.final[bm.final$Stock == s & bm.final$Year %in% first.year:last.year,]  
      #a=s
      
      #this needs to be indexed according to time
      bm.last.yr <- res.L.stock.df[[i]]$bm.stock[res.L.stock.df[[i]]$stock == s 
                                                 & res.L.stock.df[[i]]$year == t+2016]
      
      # Sort out which of the years are low or high bm
      # I'm using 0.5 as the cut off, other options are valid (0.4 is my fav...)
      
      # Get what our carrying capacity is at this moment
      cur.K <- L.Ks$adj.K[L.Ks$stock == s]
      
      #browser()
      
      #########################################################################
      if (bm.last.yr <= cur.K){
        
        bm.p <- round(bm.last.yr/cur.K, digits = 4)
        
        input.lam <- s.p.lam[s.p.lam$stock == s,]
        lam.samp <- sample(input.lam$lambda[input.lam$bm.prop <= 1], 1)

        bm.p <- NULL
      }

      if(bm.last.yr > cur.K) 

      {
        bm.p <- round(bm.last.yr/cur.K, digits = 4)
        
        input.lam <- s.p.lam[s.p.lam$stock == s,]
        lam.samp <- sample(input.lam$lambda[input.lam$bm.prop > 1], 1)
        
        bm.p <- NULL
        
      }
      #while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      
      #IN TOP DOWN, FISHING MEDIUMS DOES NOTHING TO LARGES!!!
      ex.curr <- 0
      #just make this a vector for x years from now on; same for mgmt
      removals <- bm.last.yr * ex.curr
      #apply removals to results
      tst.res <- lam.samp*(bm.last.yr - removals)
      
      
      if (t == 1) {
        res.L.stock.df[[i]]$removals[res.L.stock.df[[i]]$stock == s] <- removals
        res.L.stock.df[[i]]$lambda[res.L.stock.df[[i]]$stock == s] <- lam.samp
        res.L.stock.df[[i]]$fg.K[res.L.stock.df[[i]]$stock == s] <- sum(L.Ks$adj.K)
        res.L.stock.df[[i]]$com.K[res.L.stock.df[[i]]$stock == s] <- bm.sim.K.com$bm[t]
        res.L.stock.df[[i]]$stock.K[res.L.stock.df[[i]]$stock == s] <- cur.K
        res.L.stock.df[[i]]$multinom[res.L.stock.df[[i]]$stock == s] <- L.multinom.year[count]
        res.L.stock.df[[i]]$lam.med[res.L.stock.df[[i]]$stock == s] <- NA
        res.L.stock.df[[i]]$prop.stock.fg.bm[res.L.stock.df[[i]]$stock == s] <- cur.K/sum(L.Ks$adj.K)
        
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
                                   res.L.stock.df[[i]]$stock == s] <- sum(L.Ks$adj.K)
        res.L.stock.df[[i]]$com.K[res.L.stock.df[[i]]$year == t + 2016 &
                                    res.L.stock.df[[i]]$stock == s] <- bm.sim.K.com$bm[t]
        res.L.stock.df[[i]]$multinom[res.L.stock.df[[i]]$year == t + 2016 &
                                       res.L.stock.df[[i]]$stock == s] <- L.multinom.year[count]
        res.L.stock.df[[i]]$lam.med[res.L.stock.df[[i]]$year == t + 2016 &
                                      res.L.stock.df[[i]]$stock == s] <- NA
        res.L.stock.df[[i]]$prop.stock.fg.bm[res.L.stock.df[[i]]$year == t + 2016 &
                                               res.L.stock.df[[i]]$stock == s] <-cur.K/sum(L.Ks$adj.K)
        
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
   
    print(paste("done larges", s, t, i))
    
    
    #this is the ratio between realized biomass of LP and LP K (theoretical size of LP)
    #give ratio to mediums
    M.K.real <- bm.sim.K.M[t]/LP.bm.prop.K
    #M.K.real <- bm.sim.K.M[t]
    
    #if (M.K.real > bm.sim.K.com$bm[t]) {
    
    if (t == 1){
      rm.bm.sim.K.M <- M.K.real/10
      
      M.multinom.year <- rmultinom(1, rm.bm.sim.K.M, 
                                   prob = (mn.stock.fg$mn.prop.fg[mn.stock.fg$stock %in% M.stocks]))
      M.multinom.year <- M.multinom.year*10
    }
    
    if (t > 1) {
      res.M.bm.last.yr <- res.M.stock.df[[i]][res.M.stock.df[[i]]$year == t + 2015,]
      res.M.bm.last.yr <- res.M.bm.last.yr[,c("sim", "year", "stock", "fg", "prop.stock.fg.bm", "stock.K", "fg.K")]
      #res.M.bm.last.yr$prop.stock.fg.bm <- res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$year == t + 2016]
      res.M.bm.last.yr$new.prop.stock.fg.bm <- NA
      res.M.bm.last.yr$stock.K <- res.M.bm.last.yr$stock.K/10
      res.M.bm.last.yr$fg.K <- res.M.bm.last.yr$fg.K/10
      res.M.bm.last.yr$new.prop.stock.fg.bm <- res.M.bm.last.yr$stock.K/res.M.bm.last.yr$fg.K
      
      rm.bm.sim.K.M <- M.K.real/10
      
      M.multinom.year <- rmultinom(1, rm.bm.sim.K.M, 
                                   prob = (res.M.bm.last.yr$new.prop.stock.fg.bm))
      M.multinom.year <- M.multinom.year*10
    }
    #NEW: adjust stock Ks according to their CRIB scores
    #Adjusted RP = RP - ((p/100)*RP*CRIB)
    
    M.Ks <- data.frame(stock = M.stocks, org.crib = crib$med.crib[crib$common %in% M.stocks],
                       logit.crib = crib.sim$logit.crib[crib.sim$common %in% M.stocks], 
                       K = M.multinom.year, adj.K = 0)
    
    for (l in M.Ks$stock){
      crib.score <- M.Ks$org.crib[M.Ks$stock == l]
      K <- M.Ks$K[M.Ks$stock == l]
      M.Ks$adj.K[M.Ks$stock == l] <- K - ((env.dev/100)*K*crib.score)
    }
    
    M.Ks$per.rem <- M.Ks$adj.K/M.Ks$K
    #sum(M.Ks$adj.K)/sum(M.multinom.year)
    
    #now onto M lambdas for year 1:
    
    count <- 0
    for(s in M.stocks)
    {
      # Reset samples
      #browser()
      count <- count + 1
      stock.lambdas <- na.omit(lambdas)  |> collapse::fsubset(common == s)

      bm.last.yr <- res.M.stock.df[[i]]$bm.stock[res.M.stock.df[[i]]$stock == s 
                                                 & res.M.stock.df[[i]]$year == t+2016]
      
      
      cur.K <- M.Ks$adj.K[M.Ks$stock == s]
      #intercept <- reg.res$intercept[reg.res$stock == s]
      #slope <- reg.res$slope[reg.res$stock == s]
      
      if (bm.last.yr <= cur.K){
        bm.p <- round(bm.last.yr/cur.K, digits = 4)
        
        input.lam <- s.p.lam[s.p.lam$stock == s,]
        lam.samp <- sample(input.lam$lambda[input.lam$bm.prop <= 1], 1)
        bm.p <- NULL
      }
      
      if(bm.last.yr > cur.K) 
      {
        bm.p <- round(bm.last.yr/cur.K, digits = 4)
        
        input.lam <- s.p.lam[s.p.lam$stock == s,]
        lam.samp <- sample(input.lam$lambda[input.lam$bm.prop > 1], 1)
        
        bm.p <- NULL

      }
      
      ex.curr <- 0
      #just make this a vector for x years from now on; same for mgmt
      removals <- bm.last.yr * ex.curr
      #apply removals to results
      tst.res <- lam.samp*(bm.last.yr - removals)
      
      
      
      if (t == 1) {
        res.M.stock.df[[i]]$removals[res.M.stock.df[[i]]$stock == s] <- removals
        res.M.stock.df[[i]]$lambda[res.M.stock.df[[i]]$stock == s] <- lam.samp
        res.M.stock.df[[i]]$com.K[res.M.stock.df[[i]]$stock == s] <- bm.sim.K.com$bm[t]
        res.M.stock.df[[i]]$stock.K[res.M.stock.df[[i]]$stock == s] <- cur.K
        res.M.stock.df[[i]]$multinom[res.M.stock.df[[i]]$stock == s] <- M.multinom.year[count]
        res.M.stock.df[[i]]$lam.med[res.M.stock.df[[i]]$stock == s] <- NA
        res.M.stock.df[[i]]$fg.K[res.M.stock.df[[i]]$stock == s] <- sum(M.Ks$adj.K)
        res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$stock == s] <- cur.K/sum(M.Ks$adj.K)
        
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
                                   res.M.stock.df[[i]]$stock == s] <- sum(M.Ks$adj.K)#M.K.real
        res.M.stock.df[[i]]$com.K[res.M.stock.df[[i]]$year == t + 2016 &
                                    res.M.stock.df[[i]]$stock == s] <- bm.sim.K.com$bm[t]
        res.M.stock.df[[i]]$multinom[res.M.stock.df[[i]]$stock == s &
                                       res.M.stock.df[[i]]$year == t + 2016] <- M.multinom.year[count]
        res.M.stock.df[[i]]$lam.med[res.M.stock.df[[i]]$stock == s &
                                      res.M.stock.df[[i]]$year == t + 2016] <- NA
        res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$year == t + 2016 &
                                               res.M.stock.df[[i]]$stock == s] <-cur.K/sum(M.Ks$adj.K)
        
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

    
  } #end t loop
  
  
  prop.K.L[[i]] <- data.frame(sim = i, year = 1:n.yrs.proj, prop.com.L.K = c(prop.sim.K.L))
  com.sim[[i]] <- data.frame(sim = i, year = 1:n.yrs.proj, com.K = c(bm.sim.K.com$bm))
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
#ratio of biomass to K for FGs? Stocks?
ratio.fg <- na.omit(res.sim.df) |> group_by(fg, year, sim) |> reframe(ratio.bm.K = bm.fg/fg.K)
ratio.stock <- na.omit(res.sim.df) |> group_by(stock, fg, year, sim) |> reframe(ratio.bm.K = bm.stock/stock.K)

for (i in c(L.stocks, M.stocks)){
  stock.r <- ratio.stock[ratio.stock$stock == i,]
  print(paste(i))
  print(mean(stock.r$ratio.bm.K))
  print(median(stock.r$ratio.bm.K))
}

########################################################
#community biomass
res.L <- res.sim.df[res.sim.df$fg == "L",]
res.L <- res.L[,c("sim", "year", "fg.K")]
res.L <- unique(res.L)
res.L <- na.omit(res.L)

res.M <- unique(res.sim.df[res.sim.df$fg == "M",])
res.M <- res.M[,c("sim", "year", "fg.K")]
res.M <- unique(res.M)
res.M <- na.omit(res.M)

res.real <- data.frame(sim = res.L$sim, year = res.L$year, real.com.K = res.L$fg.K + res.M$fg.K)


com.q <- res.sim.df |>  group_by(year) |> summarize(L.50 = quantile(bm.com,probs=c(0.25),na.rm=T),
                                                    med = median(bm.com,na.rm=T),
                                                    U.50 = quantile(bm.com,probs=c(0.75),na.rm=T))

com.K.q <- res.real |>  group_by(year) |> summarize(L.50 = quantile(real.com.K,probs=c(0.25),na.rm=T),
                                                    med = median(real.com.K,na.rm=T),
                                                    U.50 = quantile(real.com.K,probs=c(0.75),na.rm=T))
hist.com.bm <- data.frame(year = 2000:2017, com.bm = com.bm)

ggplot(hist.com.bm) + geom_line(aes(x = year, y = com.bm/1000000)) +
  geom_vline(xintercept = 2017, linetype = 2) +
  geom_hline(yintercept = (median(hist.com.bm$com.bm))/1000000) +
  geom_line(data = com.K.q, aes(x = year, y = med/1000000), , colour = "lightblue") +
  geom_ribbon(data = com.K.q, aes(x = year, ymax = U.50/1000000, ymin = L.50/1000000), 
              alpha = 0.3, fill = "lightblue") +
  geom_hline(yintercept = (median(com.K.q$med))/1000000, colour = "lightblue") +
  geom_line(data = com.q, aes(x = year, y = med/1000000), , colour = "darkblue") +
  geom_ribbon(data = com.q, aes(x = year, ymax = U.50/1000000, ymin = L.50/1000000), 
              alpha = 0.1, fill = "darkblue") +
  geom_hline(yintercept = (median(com.q$med))/1000000, colour = "darkblue") +
  labs(x = "Year", y = "Community biomass (K in lightblue) \n(millions of kg)") +
  theme(axis.title = element_text(face = "bold", size = 20))

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
  stock.in <- p.lambdas$full.bm[p.lambdas$stock == i]
  mn.stock.in <- median(stock.in)
  mn.in$mean[count] <- mn.stock.in
}

for (i in 1:nrow(p.lambdas)){
  #browser()
  stock.p <- p.lambdas$stock[i]
  p.lambdas$mn.in[i] <- mn.in$mean[mn.in$stock == stock.p]
}

#colnames(p.lambdas)[colnames(p.lambdas) == "FG"] <- "fg"

for (x in 1:nrow(p.lambdas)){
  p.lambdas$FG[x] <- ifelse (p.lambdas$stock[x] %in% L.stocks, "L", "M")
}
colnames(p.lambdas)[colnames(p.lambdas) == "FG"] <- "fg"

ggplot(quants) +
  #geom_line(data = na.omit(res.sim.df[res.sim.df$fg == "L",]), aes(x = year, y = stock.K/1000000, group = sim), colour = "lightblue") +
  
  geom_ribbon(data = quants, aes(x = year, ymax = U.50/1000000, ymin = L.50/1000000), 
              alpha = 0.6, fill = "white") +
  geom_line(aes(x = year, y = med/1000000, group = stock), colour = "white") +
  geom_line(data = p.lambdas, aes(x = year, y = full.bm/1000000, group = stock), colour = "darkgrey") +
  #geom_line(data = na.omit(res.sim.df[res.sim.df$fg == "L",]), aes(x = year, y = stock.K/1000000, group = sim),
  #colour = "lightblue") +
  geom_line(data = p.lambdas, aes(x = year, y = mn.in/1000000, group = stock), colour = "white") +
  geom_vline(xintercept = 2017, linetype = 2, size = 0.9, colour = "darkgrey") +
  #scale_y_continuous(limits = c(0, max(p.lambdas$full.bm))) +
  #facet_wrap(~stock, ncol = 5, scale = "free") +
  facet_wrap(~interaction(stock, fg , sep = " FG: ",drop=T),scales='free_y') +
  labs(x = "Year", y = "Biomass (millions of kg)") +
  theme(axis.text.x = element_text(size=14, angle=45),
        axis.text.y = element_text(size=14),
        #plot.background = element_rect(fill = 'black'),
        panel.background =element_rect(fill = 'black') )
