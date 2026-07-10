
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

n.sims <- 1

env.dev <- 1

lambdas <- do.call("rbind", stocks.lst)
#truncate to years we are looking at
lambdas <- lambdas |> collapse::fsubset(year %in% yrs)


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
  s.lam <- p.lambdas[p.lambdas$stock == d & p.lambdas$year != 2000,]
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
  print(paste(p.val))
  #print(paste(summary(s.lam$lambda)))
  #print(paste(r.square))
  #print(paste(min.lam))
}

#getting average historical removals for each stock
hist.rem <- data.frame(stock = c(L.stocks, M.stocks), avg.ex = 0)

lambdas$ex <- lambdas$removals/lambdas$pop.bm

for (c in c(L.stocks, M.stocks)){
  hist.rem$avg.ex[hist.rem$stock == c] <- mean(lambdas$ex[lambdas$common == c])
  if (hist.rem$avg.ex[hist.rem$stock == c] == 0) print(paste(c,"no removals"))
}



for (i in 1:n.sims){
 
  bm.sim.K.com <- data.frame(bm = c(arima.sim(model =list(ar = com.cor.lag.1),n = n.yrs.proj + 1,n.start=1,start.innov = start.com.diff/com.cor.lag.1,
                                              innov = c(0,rnorm(n.yrs.proj,0,se.com.sim))) + med.com.sim)[2:(n.yrs.proj + 1)],
                             Years = 1:n.yrs.proj,sim = i)
  
  #1: this gives you linear decline from starting x on each original value per year
  bm.sim.K.com$env.dev <- 0
  bm.sim.K.com$adj.bm <- 0
  for (x in 1:nrow(bm.sim.K.com)){
    if (x == 1) bm.sim.K.com$env.dev[x] <- 1
    if (x == 2) bm.sim.K.com$env.dev[x] <- env.dev
    if (x > 2) bm.sim.K.com$env.dev[x] <- env.dev - ((1-env.dev)*(bm.sim.K.com$Years[x]-2))
    bm.sim.K.com$adj.bm[x] <- bm.sim.K.com$bm[x]*bm.sim.K.com$env.dev[x]
    bm.sim.K.com$adj.bm <- c(bm.sim.K.com$adj.bm)
  }#end 1
  
  
  #reducing sd, changing to 1st order and switching to mean helps
  prop.sim.K.L <-inv.logit(arima.sim(model =list(ar = c(L.com.prop.cor.lag.1)),n = n.yrs.proj + 1,
                                     n.start =1, start.innov = start.L.com.prop.diff.logit/L.com.prop.cor.lag.1, 
                                     innov = c(0,rnorm(n.yrs.proj,0,se.L.sim.logit))) + med.L.sim.logit)[2:(n.yrs.proj + 1)]
  
  
  
  bm.sim.K.L <- bm.sim.K.com$adj.bm * c(prop.sim.K.L)

  bm.sim.K.M <- bm.sim.K.com$adj.bm * (1 - c(prop.sim.K.L))
  
  for (t in 1:n.yrs.proj){
   
    if (t == 1){
      res.L.stock.df[[i]] <- data.frame(sim = i, year = 2017,
                                        stock = L.stocks, 
                                        fg = L.stocks.data$BFG,
                                        bm.stock = L.stocks.data$bm.stock,
                                        bm.fg = L.stocks.data$bm.bfg,
                                        prop.stock.fg.bm = 0,
                                        bm.com = L.stocks.data$bm.eco,
                                        removals = 0, lambda = 0,stock.K = 0, fg.K = 0, com.K = 0, multinom = 0,
                                        ex.curr = 0)
      
      res.M.stock.df[[i]] <- data.frame(sim = i, year = 2017,
                                        stock = M.stocks, 
                                        fg = M.stocks.data$BFG,
                                        bm.stock = M.stocks.data$bm.stock,
                                        bm.fg = M.stocks.data$bm.bfg,
                                        prop.stock.fg.bm = 0,
                                        bm.com = M.stocks.data$bm.eco,
                                        removals = 0, lambda = 0,stock.K = 0, fg.K = 0, com.K = 0, multinom = 0,
                                        ex.curr = 0)
      
      rm.bm.sim.K.L <- bm.sim.K.L[1]/10
      
      L.multinom.year <- c(rmultinom(1, rm.bm.sim.K.L, 
                                     prob = (mn.stock.fg$mn.prop.fg[mn.stock.fg$stock %in% L.stocks])))
      L.multinom.year <- L.multinom.year*10
      
    }
    
    if (t > 1) {
      res.L.bm.last.yr <- res.L.stock.df[[i]][res.L.stock.df[[i]]$year == t + 2015,]
      res.L.bm.last.yr <- res.L.bm.last.yr[,c("sim", "year", "stock", "fg", "prop.stock.fg.bm", "stock.K", "fg.K")]
      res.L.bm.last.yr$new.prop.stock.fg.bm <- 0
      res.L.bm.last.yr$stock.K <- res.L.bm.last.yr$stock.K/10
      res.L.bm.last.yr$fg.K <- res.L.bm.last.yr$fg.K/10
      res.L.bm.last.yr$new.prop.stock.fg.bm <- res.L.bm.last.yr$stock.K/res.L.bm.last.yr$fg.K
      #2015 is needed here because you need to fill in the t - 2 K proportion slots (this is because year 1 has both 2017 and 2018)
      
      rm.bm.sim.K.L <- bm.sim.K.L[t]/10
      
      L.multinom.year <- c(rmultinom(1, rm.bm.sim.K.L, 
                                     prob = (res.L.bm.last.yr$new.prop.stock.fg.bm)))
      L.multinom.year <- L.multinom.year*10
      
      
    }
    
    count <- 0
    # browser()
    for(s in L.stocks)
    {
      
      count <- count + 1
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
      ex.curr <- 0.04
      #ex.curr <- hist.rem$avg.ex[hist.rem$stock == s]
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
        res.L.stock.df[[i]]$ex.curr[res.L.stock.df[[i]]$stock == s] <- ex.curr
        res.L.stock.df[[i]]$prop.stock.fg.bm[res.L.stock.df[[i]]$stock == s] <- cur.K/sum(L.multinom.year)
        
        res.L.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "L",
                                         bm.stock = tst.res, bm.fg = 0, 
                                         prop.stock.fg.bm = 0,
                                         bm.com = 0,
                                         removals = 0, lambda = 0, stock.K = 0,
                                         fg.K = 0, com.K = 0, multinom = 0,
                                         ex.curr = 0)
        
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
        res.L.stock.df[[i]]$ex.curr[res.L.stock.df[[i]]$year == t + 2016 &
                                      res.L.stock.df[[i]]$stock == s] <- ex.curr
        res.L.stock.df[[i]]$prop.stock.fg.bm[res.L.stock.df[[i]]$year == t + 2016 &
                                               res.L.stock.df[[i]]$stock == s] <-cur.K/sum(L.multinom.year)
        
        res.L.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "L",
                                         bm.stock = tst.res, bm.fg = 0,
                                         prop.stock.fg.bm = 0,
                                         bm.com = 0,
                                         removals = 0, lambda = 0, stock.K = 0,
                                         fg.K = 0, com.K = 0, multinom = 0,
                                         ex.curr = 0)
        
        #in timestep t >= 2, add in bm.fg and bm.com 
      }
      
      #CAUTION
      res.L.stock.df[[i]] <- rbind(res.L.stock.df[[i]], res.L.stock.df.tmp)
      #reset temp filler vector
      #res.L.stock.df.tmp <- NULL
      
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
    
    
    

    #give ratio to mediums
    M.K.real <- bm.sim.K.M[t]/LP.bm.prop.K

    if (t == 1){
      rm.bm.sim.K.M <- M.K.real/10
      
      M.multinom.year <- c(rmultinom(1, rm.bm.sim.K.M, 
                                     prob = (mn.stock.fg$mn.prop.fg[mn.stock.fg$stock %in% M.stocks])))
      M.multinom.year <- M.multinom.year*10
    }
    
    if (t > 1) {
      res.M.bm.last.yr <- res.M.stock.df[[i]][res.M.stock.df[[i]]$year == t + 2015,]
      res.M.bm.last.yr <- res.M.bm.last.yr[,c("sim", "year", "stock", "fg", "prop.stock.fg.bm", "stock.K", "fg.K")]
      res.M.bm.last.yr$new.prop.stock.fg.bm <- NA
      res.M.bm.last.yr$stock.K <- res.M.bm.last.yr$stock.K/10
      res.M.bm.last.yr$fg.K <- res.M.bm.last.yr$fg.K/10
      res.M.bm.last.yr$new.prop.stock.fg.bm <- res.M.bm.last.yr$stock.K/res.M.bm.last.yr$fg.K
      
      rm.bm.sim.K.M <- M.K.real/10

      M.multinom.year <- c(rmultinom(1, rm.bm.sim.K.M, 
                                     prob = (res.M.bm.last.yr$new.prop.stock.fg.bm)))
      M.multinom.year <- M.multinom.year*10
    }
    
    #now onto M lambdas for year 1:
    
    count <- 0
    for(s in M.stocks)
    {
      # Reset samples
      #browser()
      count <- count + 1

      bm.last.yr <- res.M.stock.df[[i]]$bm.stock[res.M.stock.df[[i]]$stock == s 
                                                 & res.M.stock.df[[i]]$year == t+2016]
      
      
    
      cur.K <- M.multinom.year[count]
    
      if (bm.last.yr <= cur.K){
        bm.p <- round(bm.last.yr/cur.K, digits = 4)
        
        input.lam <- s.p.lam[s.p.lam$stock == s,]
        lam.samp <- sample(input.lam$lambda[input.lam$bm.prop <= 1], 1)
        
     
        bm.p <- NULL
      }
      
      # Final one, if we are above the K, we are just doing the high biomass scenario for now.
      # Solution, sample from the lambdas at high biomass, but only take lambdas that are <= 1
      if(bm.last.yr > cur.K) 
      {
        bm.p <- round(bm.last.yr/cur.K, digits = 4)
        
        input.lam <- s.p.lam[s.p.lam$stock == s,]
        lam.samp <- sample(input.lam$lambda[input.lam$bm.prop > 1], 1)
        
        bm.p <- NULL
    
      }
      #while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      
      ex.curr <- 0.04
      #ex.curr <- hist.rem$avg.ex[hist.rem$stock == s]
      #just make this a vector for x years from now on; same for mgmt
      removals <- bm.last.yr * ex.curr
      if (is.numeric(removals) == FALSE) {
        print(paste("non-numeric removals", s, t))
        print(paste(str(removals)))
      }
      #apply removals to results
      tst.res <- lam.samp*(bm.last.yr - removals)
      
      
      
      if (t == 1) {
        res.M.stock.df[[i]]$removals[res.M.stock.df[[i]]$stock == s] <- removals
        res.M.stock.df[[i]]$lambda[res.M.stock.df[[i]]$stock == s] <- lam.samp
        res.M.stock.df[[i]]$com.K[res.M.stock.df[[i]]$stock == s] <- bm.sim.K.com$adj.bm[t]
        res.M.stock.df[[i]]$stock.K[res.M.stock.df[[i]]$stock == s] <- cur.K
        res.M.stock.df[[i]]$multinom[res.M.stock.df[[i]]$stock == s] <- M.multinom.year[count]
        res.M.stock.df[[i]]$ex.curr[res.M.stock.df[[i]]$stock == s] <- ex.curr
        res.M.stock.df[[i]]$fg.K[res.M.stock.df[[i]]$stock == s] <- sum(M.multinom.year)
        res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$stock == s] <- cur.K/sum(M.multinom.year)
        
        res.M.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "M",
                                         bm.stock = tst.res, bm.fg = 0, 
                                         prop.stock.fg.bm = 0,
                                         bm.com = 0,
                                         removals = 0, lambda = 0, stock.K = 0,
                                         fg.K = 0, com.K = 0, multinom = 0,
                                         ex.curr = 0)
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
        res.M.stock.df[[i]]$ex.curr[res.M.stock.df[[i]]$stock == s &
                                      res.M.stock.df[[i]]$year == t + 2016] <- ex.curr
        res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$year == t + 2016 &
                                               res.M.stock.df[[i]]$stock == s] <-cur.K/sum(M.multinom.year)
        
        res.M.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "M",
                                         bm.stock = tst.res, bm.fg = 0,
                                         prop.stock.fg.bm = 0,
                                         bm.com = 0,
                                         removals = 0, lambda = 0, stock.K = 0,
                                         fg.K = 0, com.K = 0, multinom = 0,
                                         ex.curr = 0)
      }
      
      res.M.stock.df[[i]] <- rbind(res.M.stock.df[[i]], res.M.stock.df.tmp)
      #reset temp filler vector
      #res.M.stock.df.tmp <- NULL
      
    }#end stock K loop
    
    #sum up M biomass
    M.bm <- sum(res.M.stock.df[[i]]$bm.stock[res.M.stock.df[[i]]$year == t + 2017], na.rm = T)
    res.M.stock.df[[i]]$bm.fg[res.M.stock.df[[i]]$year == t + 2017] <- M.bm
   
    new.bm.com <- L.bm + M.bm
    res.L.stock.df[[i]]$bm.com[res.L.stock.df[[i]]$year == t + 2017] <- new.bm.com
    res.M.stock.df[[i]]$bm.com[res.M.stock.df[[i]]$year == t + 2017] <- new.bm.com
    
    
  } #end t loop
  
  L.bm <- NULL
  LP.bm. <- NULL
  M.bm <- NULL
}
res.L.stock <- do.call("rbind", res.L.stock.df)
res.M.stock <- do.call("rbind", res.M.stock.df)

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


test <- res.sim.df
test$ex <- test$removals/test$bm.stock
unique(test$ex)
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
  geom_line(aes(x = year, y = med/1000000, group = stock), colour = "darkblue") +
  geom_ribbon(data = quants, aes(x = year, ymax = U.50/1000000, ymin = L.50/1000000), 
              alpha = 0.1, fill = "darkblue") +
  geom_line(data = p.lambdas, aes(x = year, y = full.bm/1000000, group = stock)) +
  #geom_line(data = na.omit(res.sim.df[res.sim.df$fg == "L",]), aes(x = year, y = stock.K/1000000, group = sim),
  #colour = "lightblue") +
  geom_line(data = p.lambdas, aes(x = year, y = mn.in/1000000, group = stock)) +
  geom_vline(xintercept = 2017, linetype = 2, size = 0.9, ) +
  #scale_y_continuous(limits = c(0, max(p.lambdas$full.bm))) +
  #facet_wrap(~stock, ncol = 5, scale = "free") +
  facet_wrap(~interaction(stock, fg , sep = " FG: ",drop=T),scales='free_y') +
  labs(x = "Year", y = "Biomass (millions of kg)") +
  theme(axis.text.x = element_text(size=14, angle=45),
        axis.text.y = element_text(size=14))

