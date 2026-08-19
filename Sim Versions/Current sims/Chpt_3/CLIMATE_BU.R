### Autocorrelation ###

#For ARIMAs, need:
#community biomass time series + pacf() on this time series
#proportion of community biomass taken up by larges time series + pacf() (logit-transformed)
#for both these, get starting value, mean, start-mean and sd
com.bm <- eco.tot.bm.best$bm.eco  
com.cor <- pacf(com.bm)
com.cor.lag.1 <- com.cor$acf[1]
com.cor.lag.2 <- com.cor$acf[2]

start.com.sim <- com.bm[length(com.bm)]
#max.com.sim <- max(com.bm)
mn.com.sim <- mean(com.bm)
med.com.sim <- median(com.bm)
#q.com.sim <- c(quantile(com.bm, probs = c(0.75), na.rm = T))
start.com.diff <- start.com.sim - med.com.sim
#sd.com.sim <- sd(log(com.bm))
sd.com.sim <- sd(com.bm)/sqrt(length(com.bm))

M.com.prop <- bm.best.b[bm.best.b$BFG == "M",]
M.com.prop <- M.com.prop[,c("Year","prop.eco.bfg")]
M.com.prop <- unique(M.com.prop)
M.com.prop <- M.com.prop$prop.eco.bfg
start.M.com.prop <- M.com.prop[length(M.com.prop)]
M.com.prop.cor <- pacf(M.com.prop)
M.com.prop.cor.lag.1 <- M.com.prop.cor$acf[1]
M.com.prop.cor.lag.2 <- M.com.prop.cor$acf[2]

M.com.prop.logit <- logit(M.com.prop)
start.M.com.prop.logit <- M.com.prop.logit[length(M.com.prop.logit)]

mn.M.sim.logit <- mean(M.com.prop.logit)
med.M.sim.logit <- median(M.com.prop.logit)
sd.M.sim.logit <- sd(M.com.prop.logit)/sqrt(length(M.com.prop.logit))

start.M.com.prop.diff.logit <- start.M.com.prop.logit - med.M.sim.logit

#starting proportions for each stock
L.stocks.data <- bm.best.b |> collapse::fsubset(Year == 2017) |> collapse::fsubset(BFG == "L")
L.stocks <- L.stocks.data$species
L.bm.start <- L.stocks.data$bm.bfg[1]

#LP.stocks <- as.character(unique(bm.best.m$species[bm.best.m$MFG == "LP"]))

M.stocks.data <- bm.best.b |> collapse::fsubset(Year == 2017) |> collapse::fsubset(BFG == "M")
M.stocks <- M.stocks.data$species


### Initialization ###
n.yrs.proj <- 100
n.sims <- 100
env.dev <- 0.99

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


lambdas <- do.call("rbind", stocks.lst)
#truncate to years we are looking at
lambdas <- lambdas |> collapse::fsubset(year %in% yrs)


### Historical proportions ###
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


### Lambda regression ###
p.lambdas <- lambdas
colnames(p.lambdas)[colnames(p.lambdas) == "common"] <- "stock"
pred.lam <- NULL
s.p.lam <- NULL
for (d in c(L.stocks, M.stocks)){
  s.lam <- na.omit(p.lambdas[p.lambdas$stock == d,])
  max.bm <- median(na.omit(s.lam$year.minus.one))
  s.lam$bm.prop <- s.lam$year.minus.one/max.bm
  mod.tmp <- lm(log(lambda) ~ bm.prop,data = s.lam)
  #Grabbing the model p value... this is decision is very arbitrary
  p.val <- summary(mod.tmp)$coefficients[2,4]
  #Making a prediction object for the stock
  pred.tmp <- data.frame(lambda = 0,bm.prop = seq(0,1,by=0.0001),Stock = d)
  
  #Now get the appropriate lambda at each proportion of biomass
  pred.tmp$lambda <- exp(predict(mod.tmp,newdata = pred.tmp))
  pred.tmp$sd <- summary(mod.tmp)$sigma #This is residual standard error!
  r.square <- summary(mod.tmp)$r.square
  pred.tmp$bm.prop <- round(pred.tmp$bm.prop,digits=4)
  
  pred.lam <- rbind(pred.lam, pred.tmp)
  s.p.lam <- rbind(s.p.lam, s.lam)
  
  min.lam <- min(s.lam$lambda)
  #print(paste(d))
  #print(paste(p.val))
  #print(paste(summary(s.lam$lambda)))
  #print(paste(r.square))
  #print(paste(min.lam))
}


### Historical exploitation rates ###
rem.sum$avg.u <- 0
rem.sum$med.u <- 0

lambdas$u <- lambdas$removals/lambdas$full.bm

for (c in rem.sum$common){
  rem.sum$avg.u[rem.sum$common == c] <- mean(lambdas$u[lambdas$common == c])
  rem.sum$med.u[rem.sum$common == c] <- median(lambdas$u[lambdas$common == c])
}


### Start simulation ###
set.seed(3)
for (i in 1:n.sims){
  res.yr.df.tmp <- NULL
  
  bm.sim.K.com <- data.frame(bm = c(arima.sim(model =list(ar = com.cor.lag.1),n = n.yrs.proj + 1,n.start=1,start.innov = start.com.diff/com.cor.lag.1,
                                              innov = c(0,rnorm(n.yrs.proj,0,sd.com.sim))) + med.com.sim)[2:(n.yrs.proj + 1)],
                             Years = 1:n.yrs.proj,sim = i)
  
  #Climate change: this gives you linear decline from starting env.dev on each original value per year
  for (n in 1:nrow(bm.sim.K.com)){
    if (n == 1) bm.sim.K.com$env.dev[n] <- 1
    if (n == 2) bm.sim.K.com$env.dev[n] <- env.dev
    if (n > 2) bm.sim.K.com$env.dev[n] <- env.dev - ((1-env.dev)*(bm.sim.K.com$Years[n]-2))
    bm.sim.K.com$adj.bm[n] <- bm.sim.K.com$bm[n]*bm.sim.K.com$env.dev[n]
  }#end 1
  
  prop.sim.K.M <-inv.logit(arima.sim(model =list(ar = c(M.com.prop.cor.lag.1)),n = n.yrs.proj + 1,
                                     n.start =1, start.innov = start.M.com.prop.diff.logit/M.com.prop.cor.lag.1, 
                                     innov = c(0,rnorm(n.yrs.proj,0,sd.M.sim.logit))) + med.M.sim.logit)[2:(n.yrs.proj + 1)]
  
  
  bm.sim.K.M <- bm.sim.K.com$adj.bm * c(prop.sim.K.M)
  
  #this is the (stand-by) value for the larges (will be updated later after M dynamics are calculated)
  bm.sim.K.L <- bm.sim.K.com$adj.bm * (1 - c(prop.sim.K.M))
  
  for (t in 1:n.yrs.proj){
    #Get stock-specific K values using a multinomial distribution:
    #proportions of start M K that each stock uses
    if (t == 1){
      res.L.stock.df[[i]] <- data.frame(sim = i, year = 2017,
                                        stock = L.stocks, 
                                        fg = L.stocks.data$BFG,
                                        bm.stock = L.stocks.data$bm.stock,
                                        bm.fg = L.stocks.data$bm.bfg,
                                        prop.stock.fg.bm = NA,
                                        bm.com = L.stocks.data$bm.eco,
                                        removals = NA, lambda = NA, stock.K = NA, 
                                        fg.K = NA, start.K = NA, com.K = NA, multinom = NA,
                                        ex.curr = NA)
      
      
      res.M.stock.df[[i]] <- data.frame(sim = i, year = 2017,
                                        stock = M.stocks, 
                                        fg = M.stocks.data$BFG,
                                        bm.stock = M.stocks.data$bm.stock,
                                        bm.fg = M.stocks.data$bm.bfg,
                                        prop.stock.fg.bm = NA,
                                        bm.com = M.stocks.data$bm.eco,
                                        removals = NA, lambda = NA, stock.K = NA, 
                                        fg.K = NA, start.K = NA, com.K = NA, multinom = NA,
                                        ex.curr = NA)
      
      rm.bm.sim.K.M <- bm.sim.K.M[1]/10
      
      M.multinom.year <- rmultinom(1, rm.bm.sim.K.M, 
                                   prob = (mn.stock.fg$mn.prop.fg[mn.stock.fg$stock %in% M.stocks]))
      M.multinom.year <- M.multinom.year*10
      
    }
    
    if (t > 1) {
      res.M.bm.last.yr <- res.M.stock.df[[i]][res.M.stock.df[[i]]$year == t + 2015,]
      res.M.bm.last.yr <- res.M.bm.last.yr[,c("sim", "year", "stock", "fg", "prop.stock.fg.bm", "stock.K", "fg.K")]
      res.M.bm.last.yr$new.prop.stock.fg.bm <- 0
      res.M.bm.last.yr$stock.K <- res.M.bm.last.yr$stock.K/10
      res.M.bm.last.yr$fg.K <- res.M.bm.last.yr$fg.K/10
      res.M.bm.last.yr$new.prop.stock.fg.bm <- res.M.bm.last.yr$stock.K/res.M.bm.last.yr$fg.K
      #2015 is needed here because you need to fill in the t - 2 K proportion slots (this is because year 1 has both 2017 and 2018)
      
      rm.bm.sim.K.M <- bm.sim.K.M[t]/10
      
      M.multinom.year <- rmultinom(1, rm.bm.sim.K.M, 
                                   prob = (res.M.bm.last.yr$new.prop.stock.fg.bm))
      M.multinom.year <- M.multinom.year*10
    }
    
    #M Ks are now done!
    
    #now onto M lambdas:
    count <- 0
    for(s in M.stocks)
    {
      # Reset samples
      count <- count + 1
      stock.lambdas <- na.omit(lambdas) |> collapse::fsubset(common == s)
      
      bm.last.yr <- res.M.stock.df[[i]]$bm.stock[res.M.stock.df[[i]]$stock == s 
                                                 & res.M.stock.df[[i]]$year == t+2016]
      
      # Get what our resource pool is at this moment
      cur.K <- M.multinom.year[count]
      
      if (bm.last.yr <= cur.K){
        #now just sampling from input lambdas based on stock resource use
        #see baseline file for previous sampling strategies
        bm.p <- round(bm.last.yr/cur.K, digits = 4)
        
        input.lam <- s.p.lam[s.p.lam$stock == s,]
        lam.samp <- sample(input.lam$lambda[input.lam$bm.prop <= 1], 1)
        
        bm.p <- NULL
      }
      
      if(bm.last.yr > cur.K) {
        bm.p <- round(bm.last.yr/cur.K, digits = 4)
        
        input.lam <- s.p.lam[s.p.lam$stock == s,]
        lam.samp <- sample(input.lam$lambda[input.lam$bm.prop > 1], 1)
        
        bm.p <- NULL
      }
      
      #Project realized biomass
      #Using historical exploitation rates
      ex.curr <- rem.sum$med.u[rem.sum$common == s]
      #ex.curr <- 0.1
      removals <- bm.last.yr * ex.curr
      
      #realized biomass:
      tst.res <- lam.samp*(bm.last.yr - removals)
      
      #store results
      if (t == 1) {
        res.M.stock.df[[i]]$removals[res.M.stock.df[[i]]$stock == s] <- removals
        res.M.stock.df[[i]]$lambda[res.M.stock.df[[i]]$stock == s] <- lam.samp
        res.M.stock.df[[i]]$start.K[res.M.stock.df[[i]]$stock == s] <- bm.sim.K.M[t]
        res.M.stock.df[[i]]$fg.K[res.M.stock.df[[i]]$stock == s] <- sum(M.multinom.year)
        res.M.stock.df[[i]]$com.K[res.M.stock.df[[i]]$stock == s] <- bm.sim.K.com$adj.bm[t]
        res.M.stock.df[[i]]$stock.K[res.M.stock.df[[i]]$stock == s] <- cur.K
        res.M.stock.df[[i]]$multinom[res.M.stock.df[[i]]$stock == s] <- M.multinom.year[count]
        res.M.stock.df[[i]]$ex.curr[res.M.stock.df[[i]]$stock == s] <- ex.curr
        res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$stock == s] <- cur.K/sum(M.multinom.year)
        
        res.M.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "M",
                                         bm.stock = tst.res, bm.fg = NA, 
                                         prop.stock.fg.bm = NA,
                                         bm.com = NA,
                                         removals = NA, lambda = NA, stock.K = NA,
                                         fg.K = NA, start.K = NA, com.K = NA, multinom = NA,
                                         ex.curr = NA)
        
      }
      
      if (t > 1){
        #browser()
        res.M.stock.df[[i]]$removals[res.M.stock.df[[i]]$year == t + 2016 &
                                       res.M.stock.df[[i]]$stock == s] <- removals
        res.M.stock.df[[i]]$lambda[res.M.stock.df[[i]]$year == t + 2016 &
                                     res.M.stock.df[[i]]$stock == s] <- lam.samp
        res.M.stock.df[[i]]$stock.K[res.M.stock.df[[i]]$year == t + 2016 &
                                      res.M.stock.df[[i]]$stock == s] <- cur.K
        res.M.stock.df[[i]]$fg.K[res.M.stock.df[[i]]$year == t + 2016 &
                                   res.M.stock.df[[i]]$stock == s] <- sum(M.multinom.year)
        res.M.stock.df[[i]]$start.K[res.M.stock.df[[i]]$year == t + 2016 &
                                      res.M.stock.df[[i]]$stock == s] <- bm.sim.K.M[t] 
        res.M.stock.df[[i]]$com.K[res.M.stock.df[[i]]$year == t + 2016 &
                                    res.M.stock.df[[i]]$stock == s] <- bm.sim.K.com$adj.bm[t]
        res.M.stock.df[[i]]$multinom[res.M.stock.df[[i]]$year == t + 2016 &
                                       res.M.stock.df[[i]]$stock == s] <- M.multinom.year[count]
        res.M.stock.df[[i]]$ex.curr[res.M.stock.df[[i]]$year == t + 2016 &
                                      res.M.stock.df[[i]]$stock == s] <- ex.curr
        res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$year == t + 2016 &
                                               res.M.stock.df[[i]]$stock == s] <-cur.K/sum(M.multinom.year)
        
        res.M.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "M",
                                         bm.stock = tst.res, bm.fg = NA,
                                         prop.stock.fg.bm = NA,
                                         bm.com = NA,
                                         removals = NA, lambda = NA, stock.K = NA,
                                         fg.K = NA, start.K = NA, com.K = NA, multinom = NA,
                                         ex.curr = NA)
      }
      
      res.M.stock.df[[i]] <- rbind(res.M.stock.df[[i]], res.M.stock.df.tmp)
      #reset temp filler vector
      res.M.stock.df.tmp <- NULL
      
    }#end stock dynamics loop
    
    #add up all mediums
    M.bm <- sum(res.M.stock.df[[i]]$bm.stock[res.M.stock.df[[i]]$year == t + 2017], na.rm = T)
    res.M.stock.df[[i]]$bm.fg[res.M.stock.df[[i]]$year == t + 2017] <- M.bm
    
    M.K <- sum(res.M.stock.df[[i]]$stock.K[res.M.stock.df[[i]]$year == t + 2016])
    
    #get proportion of resource pool that mediums used
    M.bm.prop.K <- M.bm/M.K
    
    print(paste("done mediums", s, t, i))
    
    #Adjust large resource pool
    L.K.real <- bm.sim.K.L[t]*M.bm.prop.K
    #L.K.real <- bm.sim.K.L[t]
    
    #Get resource pools for large stocks
    if (t == 1){
      rm.bm.sim.K.L <- L.K.real/10
      
      L.multinom.year <- rmultinom(1, rm.bm.sim.K.L, 
                                   prob = (mn.stock.fg$mn.prop.fg[mn.stock.fg$stock %in% L.stocks]))
      L.multinom.year <- L.multinom.year*10
    }
    
    if (t > 1) {
      res.L.bm.last.yr <- res.L.stock.df[[i]][res.L.stock.df[[i]]$year == t + 2015,]
      res.L.bm.last.yr <- res.L.bm.last.yr[,c("sim", "year", "stock", "fg", "prop.stock.fg.bm", "stock.K", "fg.K")]
      res.L.bm.last.yr$new.prop.stock.fg.bm <- NA
      res.L.bm.last.yr$stock.K <- res.L.bm.last.yr$stock.K/10
      res.L.bm.last.yr$fg.K <- res.L.bm.last.yr$fg.K/10
      res.L.bm.last.yr$new.prop.stock.fg.bm <- res.L.bm.last.yr$stock.K/res.L.bm.last.yr$fg.K
      
      rm.bm.sim.K.L <- L.K.real/10
      
      L.multinom.year <- rmultinom(1, rm.bm.sim.K.L, 
                                   prob = (res.L.bm.last.yr$new.prop.stock.fg.bm))
      L.multinom.year <- L.multinom.year*10
    }
    
    count <- 0
    for(s in L.stocks)
    {
      count <- count + 1
      stock.lambdas <- na.omit(lambdas) |> collapse::fsubset(common == s)
      
      bm.last.yr <- res.L.stock.df[[i]]$bm.stock[res.L.stock.df[[i]]$stock == s 
                                                 & res.L.stock.df[[i]]$year == t+2016
                                                 & res.L.stock.df[[i]]$sim == i]
      
      
      cur.K <- L.multinom.year[count]
      
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
      
      #Using historical exploitation rates
      ex.curr <- rem.sum$med.u[rem.sum$common == s]
      #ex.curr <- 0.1
      removals <- bm.last.yr * ex.curr
      
      #realized biomass:
      tst.res <- lam.samp*(bm.last.yr - removals)
      
      if (t == 1) {
        res.L.stock.df[[i]]$removals[res.L.stock.df[[i]]$stock == s] <- removals
        res.L.stock.df[[i]]$lambda[res.L.stock.df[[i]]$stock == s] <- lam.samp
        res.L.stock.df[[i]]$com.K[res.L.stock.df[[i]]$stock == s] <- bm.sim.K.com$adj.bm[t]
        res.L.stock.df[[i]]$stock.K[res.L.stock.df[[i]]$stock == s] <- cur.K
        res.L.stock.df[[i]]$multinom[res.L.stock.df[[i]]$stock == s] <- L.multinom.year[count]
        res.L.stock.df[[i]]$ex.curr[res.L.stock.df[[i]]$stock == s] <- ex.curr
        res.L.stock.df[[i]]$fg.K[res.L.stock.df[[i]]$stock == s] <- sum(L.multinom.year)
        res.L.stock.df[[i]]$start.K[res.L.stock.df[[i]]$stock == s] <- bm.sim.K.L[t]
        res.L.stock.df[[i]]$prop.stock.fg.bm[res.L.stock.df[[i]]$stock == s] <- cur.K/sum(L.multinom.year)
        
        res.L.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "L",
                                         bm.stock = tst.res, bm.fg = NA, 
                                         prop.stock.fg.bm = NA,
                                         bm.com = NA,
                                         removals = NA, lambda = NA, stock.K = NA,
                                         fg.K = NA, 
                                         start.K = NA,
                                         com.K = NA, multinom = NA,
                                         ex.curr = NA)
      }
      
      if (t > 1){
        
        res.L.stock.df[[i]]$removals[res.L.stock.df[[i]]$year == t + 2016 &
                                       res.L.stock.df[[i]]$stock == s] <- removals
        res.L.stock.df[[i]]$lambda[res.L.stock.df[[i]]$year == t + 2016 &
                                     res.L.stock.df[[i]]$stock == s] <- lam.samp
        res.L.stock.df[[i]]$stock.K[res.L.stock.df[[i]]$year == t + 2016 &
                                      res.L.stock.df[[i]]$stock == s] <- cur.K
        res.L.stock.df[[i]]$fg.K[res.L.stock.df[[i]]$year == t + 2016 &
                                   res.L.stock.df[[i]]$stock == s] <- sum(L.multinom.year)
        res.L.stock.df[[i]]$start.K[res.L.stock.df[[i]]$year == t + 2016 &
                                      res.L.stock.df[[i]]$stock == s] <- bm.sim.K.L[t]
        res.L.stock.df[[i]]$com.K[res.L.stock.df[[i]]$year == t + 2016 &
                                    res.L.stock.df[[i]]$stock == s] <- bm.sim.K.com$adj.bm[t]
        res.L.stock.df[[i]]$multinom[res.L.stock.df[[i]]$stock == s &
                                       res.L.stock.df[[i]]$year == t + 2016] <- L.multinom.year[count]
        res.L.stock.df[[i]]$ex.curr[res.L.stock.df[[i]]$stock == s &
                                      res.L.stock.df[[i]]$year == t + 2016] <- ex.curr
        res.L.stock.df[[i]]$prop.stock.fg.bm[res.L.stock.df[[i]]$year == t + 2016 &
                                               res.L.stock.df[[i]]$stock == s] <-cur.K/sum(L.multinom.year)
        
        res.L.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "L",
                                         bm.stock = tst.res, bm.fg = NA, 
                                         prop.stock.fg.bm = NA,
                                         bm.com = NA,
                                         removals = NA, lambda = NA, stock.K = NA,
                                         fg.K = NA, 
                                         start.K = NA,
                                         com.K = NA, multinom = NA,
                                         ex.curr = NA)
      }
      
      res.L.stock.df[[i]] <- rbind(res.L.stock.df[[i]], res.L.stock.df.tmp)
      #reset temp filler vector
      res.L.stock.df.tmp <- NULL
      
    }#end stock K loop
    
    #sum up L biomass
    L.bm <- sum(res.L.stock.df[[i]]$bm.stock[res.L.stock.df[[i]]$year == t + 2017], na.rm = T)
    res.L.stock.df[[i]]$bm.fg[res.L.stock.df[[i]]$year == t + 2017] <- L.bm
    
    print("done larges")
    
    #get community biomass for the year
    new.bm.com <- L.bm + M.bm
    res.L.stock.df[[i]]$bm.com[res.L.stock.df[[i]]$year == t + 2017] <- new.bm.com
    res.M.stock.df[[i]]$bm.com[res.M.stock.df[[i]]$year == t + 2017] <- new.bm.com
    
  } #end t loop
  L.bm <- NULL
  LP.bm. <- NULL
  M.bm <- NULL
}#end i loop

#unpack results
res.L.stock <- do.call("rbind", res.L.stock.df)
res.M.stock <- do.call("rbind", res.M.stock.df)

res.sim.df <- rbind(res.L.stock, res.M.stock)
#get everything into one data frame

#export results
setwd("C:/Users/BARIBEAUD/Desktop/GitHub/Multispecies_model/Sim Versions/Current sims/Chpt_3")
save(res.sim.df, file = "test BU CHPT 3.RData")

save(res.sim.df, file = "test_BU-Increased med u CHPT 2.RData")