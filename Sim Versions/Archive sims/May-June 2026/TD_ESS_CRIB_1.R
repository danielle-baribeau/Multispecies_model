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


#Set up inputs
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

n.sims <- 100

#this now means something different - max percentage you want stock pool to decline by
env.dev <- 1

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
    crib$med.crib[crib$latin == l] <- 0.5
    crib$med.sd[crib$latin == l] <- 0.1
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


#for each sim, going to sample from distribution to get sim CRIB value, then proceed with rescaling and logit
#ISSUE - some of the SDs are greater than 1, which gives unrealistic CRIB scores when sampling
  #maybe ask DB about better way to handle this?

lambdas <- do.call("rbind", stocks.lst)
#truncate to years we are looking at
lambdas <- lambdas |> collapse::fsubset(year %in% yrs)

#mn.prop.stock.fg.bm$new.prop.bm.stock.bfg <- mn.prop.stock.fg.bm$bm.stock/mn.prop.stock.fg.bm$bm.bfg
mn.stock.fg <- data.frame(stock = c(L.stocks, M.stocks), mn.bm = 0, mn.fg = 0, mn.prop.fg = 0)

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
pred.lam <- NULL
s.p.lam <- NULL
for (d in c(L.stocks, M.stocks)){
  s.lam <- na.omit(p.lambdas[p.lambdas$stock == d,])
  max.bm <- median(na.omit(s.lam$year.minus.one))
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
  print(paste(p.val))
  #print(paste(summary(s.lam$lambda)))
  #print(paste(r.square))
  #print(paste(min.lam))
}


for (i in 1:n.sims){
  #browser()
  #browser()
  #STARTING WITH YEAR 1 OF PROJECTION:
  #go through entire first year (lambdas, K-space), then continue with other years
  
  #NEW: sample CRIB values for the run
  crib.sim <- data.frame(common = crib$common, latin = crib$latin, samp.crib = 0, resc.crib = 0, logit.crib = 0)
  for (l in crib.sim$latin){
    #going to set SD to 0.1 for now, until we discuss sampling issue
    #chose lognormal because of skewedness in data
    
    crib.sim$samp.crib[crib.sim$latin == l] <- rlnorm(1, log(crib$med.crib[crib$latin == l]), 0.1)
    
    #tried this approach, but some of SDs were too big
      #TO DO: check this with 100 sims
    #crib.sim$samp.crib[crib.sim$latin == l] <- rlnorm(1, log(crib$med.crib[crib$latin == l]), crib$med.sd[crib$latin == l]/2)
    #while (crib.sim$samp.crib[crib.sim$latin == l] > 1 || crib.sim$samp.crib[crib.sim$latin == l] < 0) crib.sim$samp.crib[crib.sim$latin == l] <- rlnorm(1, log(crib$med.crib[crib$latin == l]), crib$med.sd[crib$latin == l]/2)
    
    #crib.sim$samp.crib[crib.sim$latin == l] <- rlnorm(1, log(crib$avg.crib[crib$latin == l]), crib$avg.sd[crib$latin == l])
    #crib.sim$samp.crib[crib.sim$latin == l] <- rnorm(1, crib$avg.crib[crib$latin == l], crib$avg.sd[crib$latin == l])
  }#end of CRIB sampling loop
  
  #rescale CRIB - 2 options
  #1. rescale by dividing by max value
    #in this scenario, CRIB penalty becomes punitive for all species except top CRIB score
  crib.sim$resc.crib <- crib.sim$samp.crib/max(crib.sim$samp.crib)
  
  #2. rescale with logit transformation
    #centers distribution around 0; lets you look at winners (> 0.5) and losers (< 0.5)
  crib.sim$logit.crib <- logit(crib.sim$samp.crib)
  
  #tested with both options:
    #1. 
  
  #NOTE: in my notes I wrote down doing both steps - but logit transforming a 1 is Infinity
    #misunderstanding?
  
  bm.sim.K.com <- data.frame(bm = c(arima.sim(model =list(ar = com.cor.lag.1),n = n.yrs.proj + 1,n.start=1,start.innov = start.com.diff/com.cor.lag.1,
                                              innov = c(0,rnorm(n.yrs.proj,0,se.com.sim))) + med.com.sim)[2:(n.yrs.proj + 1)],
                             Years = 1:n.yrs.proj,sim = i)
  #don't want first year because you're starting in 2017
  #SE makes simulated data bounds narrower than input data, but need to have constraints somewhere
    #mean is bang-on still
  
  #proportion of community taken up by larges
  prop.sim.K.L <-inv.logit(arima.sim(model =list(ar = c(L.com.prop.cor.lag.1)),n = n.yrs.proj + 1,
                                     n.start =1, start.innov = start.L.com.prop.diff.logit/L.com.prop.cor.lag.1, 
                                     innov = c(0,rnorm(n.yrs.proj,0,se.L.sim.logit))) + med.L.sim.logit)[2:(n.yrs.proj + 1)]
  
  bm.sim.K.L <- bm.sim.K.com$bm * c(prop.sim.K.L)
  
  #this is the (stand-by) value for the mediums (will be updated later after L dynamics are calculated)
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
                                        lam.med = NA, crib.sim = NA)
      
      res.M.stock.df[[i]] <- data.frame(sim = i, year = M.stocks.data$Year,
                                        stock = M.stocks, 
                                        fg = M.stocks.data$BFG,
                                        bm.stock = M.stocks.data$bm.stock,
                                        bm.fg = M.stocks.data$bm.bfg,
                                        prop.stock.fg.bm = NA,
                                        bm.com = M.stocks.data$bm.eco,
                                        removals = NA, lambda = NA,stock.K = NA, fg.K = NA, com.K = NA, multinom = NA,
                                        lam.med = NA, crib.sim = NA)
      
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
      
      L.multinom.year <- rmultinom(1, rm.bm.sim.K.L, 
                                   prob = (res.L.bm.last.yr$new.prop.stock.fg.bm))
      L.multinom.year <- L.multinom.year*10
      
      
    }
    
    #NEW: adjust stock Ks according to their CRIB scores
    #Adjusted RP = RP - ((p/100)*RP*CRIB)
    
    L.Ks <- data.frame(stock = L.stocks, org.crib = crib$med.crib[crib$common %in% L.stocks],
                       logit.crib = crib.sim$logit.crib[crib.sim$common %in% L.stocks], 
                       resc.crib = crib.sim$resc.crib[crib.sim$common %in% L.stocks],
                       K = L.multinom.year, adj.K = 0)
    
    for (l in L.Ks$stock){
      crib.score <- crib.sim$resc.crib[crib.sim$common == l]
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
      #cur.K <- L.multinom.year[count]
      cur.K <- L.Ks$adj.K[L.Ks$stock == s]
      
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
        lam.sd <- pred.lam$sd[pred.lam$bm.prop == bm.p & pred.lam$Stock == s]#this is SD
        lam.samp <- rlnorm(1,log(lam.mn), lam.sd)
        #while(lam.samp > 10 || lam.samp < 0.1) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
        #if (bm.p < 0.01) while(lam.samp < 1 | lam.samp > 1.5*max(stock.lambdas$lambda)) lam.samp <- rlnorm(1, log(lam.mn), lam.sd)
        #if (bm.p >= 0.01)
        while(lam.samp > 5 | lam.samp < 0.2) lam.samp <- rlnorm(1, log(lam.mn), lam.sd)
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
        if(lam.mn > 0.5) while(lam.samp > mean(stock.lambdas$lambda) | lam.samp < 0.2) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
        #if(lam.mn < 0.5) while(lam.samp > (0.5*mean(stock.lambdas$lambda)) | lam.samp < (0.5*(min(stock.lambdas$lambda)))) lam.samp <- rlnorm(1,log(lam.mn*0.8),lam.sd)
        if (lam.mn < 0.5) while(lam.samp > mean(stock.lambdas$lambda) | lam.samp < 0.2) lam.samp <- rlnorm(1,log(0.8*lam.mn),lam.sd)
          #this 
        
        #while(lam.samp > median(stock.lambdas$lambda) || lam.samp < 0.5*min(stock.lambdas$lambda)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
        #while(lam.samp < 0.1 || lam.samp > 10) lam.samp <- rlnorm(1,log(lam.med),lam.sd)
        #while (lam.samp < (min(stock.lambdas$lambda)) || lam.samp > 1.5*(max(stock.lambdas$lambda))) lam.samp <- rlnorm(1,log(lam.med),lam.sd)
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
        res.L.stock.df[[i]]$multinom[res.L.stock.df[[i]]$stock == s] <- L.Ks$adj.K[L.Ks$stock == s]
        res.L.stock.df[[i]]$lam.med[res.L.stock.df[[i]]$stock == s] <- lam.mn
        res.L.stock.df[[i]]$prop.stock.fg.bm[res.L.stock.df[[i]]$stock == s] <- cur.K/sum(L.Ks$adj.K)
        res.L.stock.df[[i]]$crib.sim[res.L.stock.df[[i]]$stock == s] <- crib.sim$samp.crib[crib.sim$common == s]
        
        res.L.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "L",
                                         bm.stock = tst.res, bm.fg = NA, 
                                         prop.stock.fg.bm = NA,
                                         bm.com = NA,
                                         removals = NA, lambda = NA, stock.K = NA,
                                         fg.K = NA, com.K = NA, multinom = NA,
                                         lam.med = NA, crib.sim = NA)
        
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
                                       res.L.stock.df[[i]]$stock == s] <- L.Ks$adj.K[L.Ks$stock == s]
        res.L.stock.df[[i]]$lam.med[res.L.stock.df[[i]]$year == t + 2016 &
                                      res.L.stock.df[[i]]$stock == s] <- lam.mn
        res.L.stock.df[[i]]$prop.stock.fg.bm[res.L.stock.df[[i]]$year == t + 2016 &
                                               res.L.stock.df[[i]]$stock == s] <-cur.K/sum(L.Ks$adj.K)
        res.L.stock.df[[i]]$crib.sim[res.L.stock.df[[i]]$stock == s] <- crib.sim$samp.crib[crib.sim$common == s]
        
        res.L.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "L",
                                         bm.stock = tst.res, bm.fg = NA,
                                         prop.stock.fg.bm = NA,
                                         bm.com = NA,
                                         removals = NA, lambda = NA, stock.K = NA,
                                         fg.K = NA, com.K = NA, multinom = NA,
                                         lam.med = NA, crib.sim = NA)
        
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
                       #logit.crib = crib.sim$logit.crib[crib.sim$common %in% M.stocks], 
                       resc.crib = crib.sim$resc.crib[crib.sim$common %in% M.stocks],
                       K = M.multinom.year, adj.K = 0)
    
    for (l in M.Ks$stock){
      crib.score <- crib.sim$resc.crib[crib.sim$common == l]
      K <- M.Ks$K[M.Ks$stock == l]
      M.Ks$adj.K[M.Ks$stock == l] <- K - ((env.dev/100)*K*crib.score)
    }
    
    M.Ks$per.rem <- M.Ks$adj.K/M.Ks$K
    sum(M.Ks$adj.K)/sum(M.multinom.year)
    
    
    #now onto M lambdas for year 1:
    
    count <- 0
    for(s in M.stocks)
    {
      # Reset samples
      #browser()
      count <- count + 1
      stock.lambdas <- na.omit(lambdas)  |> collapse::fsubset(common == s)
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
                                                 & res.M.stock.df[[i]]$year == t+2016]
      
      
      # Sort out which of the years are low or high bm
      # I'm using 0.5 as the cut off, other options are valid (0.4 is my fav...)
      #if (s %in% c("Redfish", "Silver hake")){
      #browser()
      #}
      #cur.K <- M.multinom.year[count]
      cur.K <- M.Ks$adj.K[M.Ks$stock == s]
      #intercept <- reg.res$intercept[reg.res$stock == s]
      #slope <- reg.res$slope[reg.res$stock == s]
      
      if (bm.last.yr <= cur.K){
        bm.p <- round(bm.last.yr/cur.K, digits = 4)
        #if (t > 1){
          #bm.p <- round(bm.last.yr/res.L.stock.df[[i]]$stock.K[res.L.stock.df[[i]]$stock == s &
                                                                 #res.L.stock.df[[i]]$year == t + 2015], digits = 4)
        #}

        lam.mn <- pred.lam$lambda[pred.lam$bm.prop == bm.p & pred.lam$Stock == s]
        lam.sd <- pred.lam$sd[pred.lam$bm.prop == bm.p & pred.lam$Stock == s]
        #THIS IS SD!!
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
        while(lam.samp > 5 | lam.samp < 0.2) lam.samp <- rlnorm(1, log(lam.mn), lam.sd)
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
        if(lam.mn > 0.5) while(lam.samp > mean(stock.lambdas$lambda) | lam.samp < 0.2) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
        #if(lam.mn < 0.5) while(lam.samp > (0.5*mean(stock.lambdas$lambda)) | lam.samp < (0.5*(min(stock.lambdas$lambda)))) lam.samp <- rlnorm(1,log(lam.mn*0.8),lam.sd)
        if (lam.mn < 0.5) while(lam.samp > mean(stock.lambdas$lambda) | lam.samp <  0.2) lam.samp <- rlnorm(1,log(0.8*lam.mn),lam.sd)
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
        res.M.stock.df[[i]]$com.K[res.M.stock.df[[i]]$stock == s] <- bm.sim.K.com$bm[t]
        res.M.stock.df[[i]]$stock.K[res.M.stock.df[[i]]$stock == s] <- cur.K
        res.M.stock.df[[i]]$multinom[res.M.stock.df[[i]]$stock == s] <- M.Ks$adj.K[M.Ks$stock == s]
        res.M.stock.df[[i]]$lam.med[res.M.stock.df[[i]]$stock == s] <- lam.mn
        res.M.stock.df[[i]]$fg.K[res.M.stock.df[[i]]$stock == s] <- sum(M.Ks$adj.K)
        res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$stock == s] <- cur.K/sum(M.Ks$adj.K)
        res.M.stock.df[[i]]$crib.sim[res.M.stock.df[[i]]$stock == s] <- crib.sim$samp.crib[crib.sim$common == s]
        
        res.M.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "M",
                                         bm.stock = tst.res, bm.fg = NA, 
                                         prop.stock.fg.bm = NA,
                                         bm.com = NA,
                                         removals = NA, lambda = NA, stock.K = NA,
                                         fg.K = NA, com.K = NA, multinom = NA,
                                         lam.med = NA, crib.sim = NA)
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
                                       res.M.stock.df[[i]]$year == t + 2016] <- M.Ks$adj.K[M.Ks$stock == s]
        res.M.stock.df[[i]]$lam.med[res.M.stock.df[[i]]$stock == s &
                                      res.M.stock.df[[i]]$year == t + 2016] <- lam.mn
        res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$year == t + 2016 &
                                               res.M.stock.df[[i]]$stock == s] <-cur.K/sum(M.Ks$adj.K)
        res.M.stock.df[[i]]$crib.sim[res.M.stock.df[[i]]$stock == s] <- crib.sim$samp.crib[crib.sim$common == s]
        
        res.M.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "M",
                                         bm.stock = tst.res, bm.fg = NA,
                                         prop.stock.fg.bm = NA,
                                         bm.com = NA,
                                         removals = NA, lambda = NA, stock.K = NA,
                                         fg.K = NA, com.K = NA, multinom = NA,
                                         lam.med = NA, crib.sim = NA)
      }
      
      res.M.stock.df[[i]] <- rbind(res.M.stock.df[[i]], res.M.stock.df.tmp)
      #reset temp filler vector
      res.M.stock.df.tmp <- NULL
      
    }#end stock K loop
    
    #sum up M biomass
    M.bm <- sum(res.M.stock.df[[i]]$bm.stock[res.M.stock.df[[i]]$year == t + 2017], na.rm = T)
    res.M.stock.df[[i]]$bm.fg[res.M.stock.df[[i]]$year == t + 2017] <- M.bm
    
    print("done mediums")
    
    #get community biomass for the year
    new.bm.com <- L.bm + M.bm
    res.L.stock.df[[i]]$bm.com[res.L.stock.df[[i]]$year == t + 2017] <- new.bm.com
    res.M.stock.df[[i]]$bm.com[res.M.stock.df[[i]]$year == t + 2017] <- new.bm.com

    
  } #end t loop
  
  
  prop.K.L[[i]] <- data.frame(sim = i, year = 1:30, prop.com.L.K = c(prop.sim.K.L))
  com.sim[[i]] <- data.frame(sim = i, year = 1:30, com.K = c(bm.sim.K.com$bm))
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

#plotting simulated crib scores
crib.res <- res.sim.df[,c("sim", "stock", "crib.sim")]
crib.res <- unique(crib.res)
crib.res <- na.omit(crib.res)
crib$sd.low <- crib$med.crib - crib$med.sd
crib$sd.high <- crib$med.crib + crib$med.sd

ggplot(crib.res) + geom_boxplot(aes(x = stock, y = crib.sim)) + 
  geom_point(data = crib, aes(x = common, y = med.crib), colour = "blue") +
  #geom_errorbar(data = crib, aes(x = common, y = med.crib, ymin = sd.low, ymax = sd.high), colour = "blue")+
  theme(axis.text.x = element_text(size = 14, angle=90),
        axis.text.y = element_text(size=14))

########################################################

########################################################


quants <- res.sim.df |>  collapse::fgroup_by(year,stock,fg) |> collapse::fsummarize(L.50 = quantile(bm.stock,probs=c(0.25),na.rm=T),
                                                                                    med = median(bm.stock,na.rm=T),
                                                                                    U.50 = quantile(bm.stock,probs=c(0.75),na.rm=T))
mn.in <- data.frame(stock = c(L.stocks, M.stocks), fg = c((rep("L", length(L.stocks))), (rep("M", length(M.stocks)))), mean = 0)
count <- 0

p.lambdas <- lambdas
colnames(p.lambdas)[colnames(p.lambdas) == "common"] <- "stock"
colnames(p.lambdas)[colnames(p.lambdas) == "FG"] <- "fg"
p.lambdas$mn.in <- rep(0, nrow(p.lambdas))
for (i in 1:nrow(p.lambdas)){
  if (grepl("M", p.lambdas$fg[i]) == TRUE) p.lambdas$fg[i] <- "M"
  if (grepl("L", p.lambdas$fg[i]) == TRUE) p.lambdas$fg[i] <- "L"
}

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

ggplot(quants) +
  geom_ribbon(data = quants, aes(x = year, ymax = U.50/1000000, ymin = L.50/1000000), alpha = 0.6, fill = "white") +
  geom_line(aes(x = year, y = med/1000000, group = stock), colour = "white") +
  geom_line(data = p.lambdas, aes(x = year, y = full.bm/1000000, group = stock), colour = "darkgrey") +
  geom_line(data = p.lambdas, aes(x = year, y = mn.in/1000000, group = stock), colour = "white") +
  geom_vline(xintercept = 2017, linetype = 2, size = 0.9, colour = "darkgrey") +
  #facet_wrap(~stock, scales = 'free_y') +
  facet_wrap(~interaction(stock, fg , sep = " FG: ",drop=T),scales='free_y') +
  labs(x = "Year", y = "Biomass (millions of kg)") +
  theme(axis.text.x = element_text(size=14, angle=45),
        axis.text.y = element_text(size=14),
        panel.background =element_rect(fill = 'black'))


#community biomass
com.quants <- res.sim.df |>  collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(bm.com,probs=c(0.25),na.rm=T),
                                                                                    med = median(bm.com,na.rm=T),
                                                                                    U.50 = quantile(bm.com,probs=c(0.75),na.rm=T))

med.com <- median(com.bm)
hist.com <- data.frame(year = 2000:2017, com.bm = com.bm)

ggplot(com.quants) +
  geom_ribbon(aes(x = year, ymax = U.50/1000000, ymin = L.50/1000000), alpha = 0.6, fill = "white") +
  geom_line(aes(x = year, y = med/1000000), colour = "white") +
  geom_line(data = hist.com, aes(x = year, y = com.bm/1000000), colour = "darkgrey") +
  geom_hline(yintercept = med.com/1000000, colour = "white") +
  geom_vline(xintercept = 2017, linetype = 2, size = 0.9, colour = "darkgrey") +
  labs(x = "Year", y = "Community Biomass (millions of kg)") +
  theme(axis.text.x = element_text(size=14, angle=45),
        axis.text.y = element_text(size=14),
        panel.background =element_rect(fill = 'black'))

