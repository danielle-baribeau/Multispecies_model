### Autocorrelation ###

#For ARIMAs, need:
#community biomass time series + pacf() on this time series
#proportion of community biomass taken up by larges time series + pacf() (logit-transformed)
#for both these, get starting value, mean, start-mean and sd
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


### Initialization ###
n.yrs.proj <- 100
n.sims <- 500

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

#Plots
#{ggplot(na.omit(s.p.lam[s.p.lam$stock == "Atlantic cod",])) +
#geom_point(aes(x = bm.prop, y = log(lambda))) +
  #geom_smooth(aes(x = bm.prop, y = log(lambda)), method = "lm", colour = "darkorange") +
  #facet_wrap(~stock, scale = "free", ncol = 5) +
  #geom_hline(yintercept = 0) +
  #labs(x = "Biomass : Median historical biomass", y = "Log of growth rate") +
  #annotate("text", x = 3.5, y = 1.8, label = "p-value = 0.022") +
  #ggtitle("Atlantic cod") +
  #theme(axis.text.x = element_text(size = 25),
        #axis.text.y = element_text(size = 25),
        #axis.title.x = element_text(size = 30),
        #axis.title.y = element_text(size = 30),
        #title = element_text(size = 35))

#ggplot(na.omit(s.p.lam[s.p.lam$FG %in% c("MP", "MBZ"),])) +
  #geom_point(aes(x = bm.prop, y = log(lambda))) +
  #geom_smooth(aes(x = bm.prop, y = log(lambda)), method = "lm", colour = "darkorange") +
  #facet_wrap(~stock, scale = "free", ncol = 4) +
  #geom_hline(yintercept = 0) +
  #geom_vline(xintercept = 1) +
  #labs(x = "Biomass : Median historical biomass", y = "Log of growth rate")
#}


### Historical exploitation rates ###
rem.sum$avg.u <- 0
rem.sum$med.u <- 0

lambdas$u <- lambdas$removals/lambdas$full.bm

for (c in rem.sum$common){
  rem.sum$avg.u[rem.sum$common == c] <- mean(lambdas$u[lambdas$common == c])
  rem.sum$med.u[rem.sum$common == c] <- median(lambdas$u[lambdas$common == c])
}


### Start simulation ###
set.seed(2)
for (i in 1:n.sims){
  bm.sim.K.com <- data.frame(bm = c(arima.sim(model =list(ar = com.cor.lag.1),n = n.yrs.proj + 1,n.start=1,start.innov = start.com.diff/com.cor.lag.1,
                                              innov = c(0,rnorm(n.yrs.proj,0,se.com.sim))) + med.com.sim)[2:(n.yrs.proj + 1)],
                             Years = 1:n.yrs.proj,sim = i)

  prop.sim.K.L <-inv.logit(arima.sim(model =list(ar = c(L.com.prop.cor.lag.1)),n = n.yrs.proj + 1,
                                     n.start =1, start.innov = start.L.com.prop.diff.logit/L.com.prop.cor.lag.1, 
                                     innov = c(0,rnorm(n.yrs.proj,0,se.L.sim.logit))) + med.L.sim.logit)[2:(n.yrs.proj + 1)]
  

  bm.sim.K.L <- bm.sim.K.com$bm * c(prop.sim.K.L)
  
  #this is the (stand-by) value for the mediums (will be updated later after L dynamics are calculated)
  #bm.sim.K.M <- bm.sim.K.com$bm * (1 - c(prop.sim.K.L))
  bm.sim.K.M <- bm.sim.K.com$bm * (1 - c(prop.sim.K.L))
  
  for (t in 1:n.yrs.proj){
  #Get stock-specific resource pools using a multinomial distribution:
    #proportions of start M resource pool that each stock uses
    if (t == 1){
      res.L.stock.df[[i]] <- data.frame(sim = i, year = 2017,
                                        stock = L.stocks, 
                                        fg = L.stocks.data$BFG,
                                        bm.stock = L.stocks.data$bm.stock,
                                        bm.fg = L.stocks.data$bm.bfg,
                                        prop.stock.fg.bm = NA,
                                        bm.com = L.stocks.data$bm.eco,
                                        removals = NA, lambda = NA,stock.K = NA, fg.K = NA, 
                                        start.K = NA, com.K = NA, multinom = NA,
                                        ex.curr = NA)
      
      res.M.stock.df[[i]] <- data.frame(sim = i, year = 2017,
                                        stock = M.stocks, 
                                        fg = M.stocks.data$BFG,
                                        bm.stock = M.stocks.data$bm.stock,
                                        bm.fg = M.stocks.data$bm.bfg,
                                        prop.stock.fg.bm = NA,
                                        bm.com = M.stocks.data$bm.eco,
                                        removals = NA, lambda = NA, stock.K = NA, fg.K = NA,
                                        start.K = NA, com.K = NA, multinom = NA,
                                        ex.curr = NA)
      
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
    
    #L Ks are now done!
    
    #now onto L lambdas:
    count <- 0
    for(s in L.stocks)
    {
      #Reset samples
      count <- count + 1
      stock.lambdas <- na.omit(lambdas) |> collapse::fsubset(common == s)
      
      bm.last.yr <- res.L.stock.df[[i]]$bm.stock[res.L.stock.df[[i]]$stock == s 
                                                 & res.L.stock.df[[i]]$year == t+2016]

      
      # Get what our resource pool is at this moment
      cur.K <- L.multinom.year[count]
      
      #Sample lambda
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
    #Project realized biomass
      #Using historical exploitation rates
      ex.curr <- rem.sum$med.u[rem.sum$common == s]
      #ex.curr <- ifelse (s %in% LP.stocks, 0.1, rem.sum$avg.u[rem.sum$common == s])
      #ex.curr <- 0.1
      removals <- bm.last.yr * ex.curr
      
      #realized biomass:
      tst.res <- lam.samp*(bm.last.yr - removals)
      
      #store results
      if (t == 1) {
        res.L.stock.df[[i]]$removals[res.L.stock.df[[i]]$stock == s] <- removals
        res.L.stock.df[[i]]$lambda[res.L.stock.df[[i]]$stock == s] <- lam.samp
        res.L.stock.df[[i]]$fg.K[res.L.stock.df[[i]]$stock == s] <- sum(L.multinom.year)
        res.L.stock.df[[i]]$com.K[res.L.stock.df[[i]]$stock == s] <- bm.sim.K.com$bm[t]
        res.L.stock.df[[i]]$stock.K[res.L.stock.df[[i]]$stock == s] <- cur.K
        res.L.stock.df[[i]]$multinom[res.L.stock.df[[i]]$stock == s] <- L.multinom.year[count]
        res.L.stock.df[[i]]$ex.curr[res.L.stock.df[[i]]$stock == s] <- ex.curr
        res.L.stock.df[[i]]$prop.stock.fg.bm[res.L.stock.df[[i]]$stock == s] <- cur.K/sum(L.multinom.year)
        res.L.stock.df[[i]]$start.K[res.M.stock.df[[i]]$stock == s] <- bm.sim.K.L[t]
        
        res.L.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "L",
                                         bm.stock = tst.res, bm.fg = NA, 
                                         prop.stock.fg.bm = NA,
                                         bm.com = NA,
                                         removals = NA, lambda = NA, stock.K = NA,
                                         fg.K = NA, 
                                         start.K = NA, com.K = NA, multinom = NA,
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
        res.L.stock.df[[i]]$com.K[res.L.stock.df[[i]]$year == t + 2016 &
                                    res.L.stock.df[[i]]$stock == s] <- bm.sim.K.com$bm[t]
        res.L.stock.df[[i]]$multinom[res.L.stock.df[[i]]$year == t + 2016 &
                                       res.L.stock.df[[i]]$stock == s] <- L.multinom.year[count]
        res.L.stock.df[[i]]$ex.curr[res.L.stock.df[[i]]$year == t + 2016 &
                                      res.L.stock.df[[i]]$stock == s] <- ex.curr
        res.L.stock.df[[i]]$prop.stock.fg.bm[res.L.stock.df[[i]]$year == t + 2016 &
                                               res.L.stock.df[[i]]$stock == s] <-cur.K/sum(L.multinom.year)
        res.L.stock.df[[i]]$start.K[res.L.stock.df[[i]]$year == t + 2016 &
                                      res.L.stock.df[[i]]$stock == s] <- bm.sim.K.L[t]
        
        res.L.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "L",
                                         bm.stock = tst.res, bm.fg = NA,
                                         prop.stock.fg.bm = NA,
                                         bm.com = NA,
                                         removals = NA, lambda = NA, stock.K = NA,
                                         fg.K = NA,
                                         start.K = NA, com.K = NA, multinom = NA,
                                         ex.curr = NA)
      }
      
      res.L.stock.df[[i]] <- rbind(res.L.stock.df[[i]], res.L.stock.df.tmp)
      
    }#end stock dynamics loop
    
    #add up all larges
    L.bm <- sum(res.L.stock.df[[i]]$bm.stock[res.L.stock.df[[i]]$year == t + 2017], na.rm = T)
    res.L.stock.df[[i]]$bm.fg[res.L.stock.df[[i]]$year == t + 2017] <- L.bm
    
    #add up all large psicivores
    LP.bm <- sum(res.L.stock.df[[i]]$bm.stock[res.L.stock.df[[i]]$stock %in% LP.stocks 
                                              & res.L.stock.df[[i]]$year == t + 2017])
    
    #get large piscivore resource pool
    LP.K <- sum(res.L.stock.df[[i]]$stock.K[res.L.stock.df[[i]]$stock %in% LP.stocks 
                                            & res.L.stock.df[[i]]$year == t + 2016])
    #get proportion of resource pool that large piscivores used
    LP.bm.prop.K <- LP.bm/LP.K

    print(paste("done larges", s, t, i))
    
    #Adjust medium resource pool
    M.K.real <- bm.sim.K.M[t]/LP.bm.prop.K
    #M.K.real <- bm.sim.K.M[t]
    
    #Get resource pools for medium stocks
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
      count <- count + 1
      stock.lambdas <- na.omit(lambdas)  |> collapse::fsubset(common == s)

      bm.last.yr <- res.M.stock.df[[i]]$bm.stock[res.M.stock.df[[i]]$stock == s 
                                                 & res.M.stock.df[[i]]$year == t+2016]
      
      #get stock's resource pool
      cur.K <- M.multinom.year[count]
      
      #Sample lambda
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
      
      #Project realized biomass
      #Using historical exploitation rates
      ex.curr <- rem.sum$med.u[rem.sum$common == s]
      removals <- bm.last.yr * ex.curr
      
      #realized biomass:
      tst.res <- lam.samp*(bm.last.yr - removals)
      
    #store results
      if (t == 1) {
        res.M.stock.df[[i]]$removals[res.M.stock.df[[i]]$stock == s] <- removals
        res.M.stock.df[[i]]$lambda[res.M.stock.df[[i]]$stock == s] <- lam.samp
        res.M.stock.df[[i]]$com.K[res.M.stock.df[[i]]$stock == s] <- bm.sim.K.com$bm[t]
        res.M.stock.df[[i]]$stock.K[res.M.stock.df[[i]]$stock == s] <- cur.K
        res.M.stock.df[[i]]$multinom[res.M.stock.df[[i]]$stock == s] <- M.multinom.year[count]
        res.M.stock.df[[i]]$ex.curr[res.M.stock.df[[i]]$stock == s] <- ex.curr
        res.M.stock.df[[i]]$fg.K[res.M.stock.df[[i]]$stock == s] <- sum(M.multinom.year)
        res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$stock == s] <- cur.K/sum(M.multinom.year)
        res.M.stock.df[[i]]$start.K[res.M.stock.df[[i]]$stock == s] <- bm.sim.K.M[t]
        
        
        res.M.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "M",
                                         bm.stock = tst.res, bm.fg = NA, 
                                         prop.stock.fg.bm = NA,
                                         bm.com = NA,
                                         removals = NA, lambda = NA, stock.K = NA,
                                         fg.K = NA,
                                         start.K = NA, com.K = NA, multinom = NA,
                                         ex.curr = NA)
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
                                    res.M.stock.df[[i]]$stock == s] <- bm.sim.K.com$bm[t]
        res.M.stock.df[[i]]$multinom[res.M.stock.df[[i]]$stock == s &
                                       res.M.stock.df[[i]]$year == t + 2016] <- M.multinom.year[count]
        res.M.stock.df[[i]]$ex.curr[res.M.stock.df[[i]]$stock == s &
                                      res.M.stock.df[[i]]$year == t + 2016] <- ex.curr
        res.M.stock.df[[i]]$prop.stock.fg.bm[res.M.stock.df[[i]]$year == t + 2016 &
                                               res.M.stock.df[[i]]$stock == s] <-cur.K/sum(M.multinom.year)
        res.M.stock.df[[i]]$start.K[res.M.stock.df[[i]]$year == t + 2016 &
                                      res.M.stock.df[[i]]$stock == s] <- bm.sim.K.M[t] 
        
        res.M.stock.df.tmp <- data.frame(sim = i, year = t + 2017, stock = s,
                                         fg = "M",
                                         bm.stock = tst.res, bm.fg = NA,
                                         prop.stock.fg.bm = NA,
                                         bm.com = NA,
                                         removals = NA, lambda = NA, stock.K = NA,
                                         fg.K = NA,
                                         start.K = NA, com.K = NA, multinom = NA,
                                         ex.curr = NA)
      }
      
      res.M.stock.df[[i]] <- rbind(res.M.stock.df[[i]], res.M.stock.df.tmp)
      
    }#end stock dynamics loop
    
    #sum up M biomass
    M.bm <- sum(res.M.stock.df[[i]]$bm.stock[res.M.stock.df[[i]]$year == t + 2017], na.rm = T)
    res.M.stock.df[[i]]$bm.fg[res.M.stock.df[[i]]$year == t + 2017] <- M.bm
    

    #get community biomass for the year
    new.bm.com <- L.bm + M.bm
    res.L.stock.df[[i]]$bm.com[res.L.stock.df[[i]]$year == t + 2017] <- new.bm.com
    res.M.stock.df[[i]]$bm.com[res.M.stock.df[[i]]$year == t + 2017] <- new.bm.com

    
  } #end t loop

  L.bm <- NULL
  LP.bm. <- NULL
  M.bm <- NULL
} #end i loop

res.L.stock <- do.call("rbind", res.L.stock.df)
res.M.stock <- do.call("rbind", res.M.stock.df)

res.sim.df <- rbind(res.L.stock, res.M.stock)
#get everything into one data frame

#export results
setwd("C:/Users/BARIBEAUD/Desktop/GitHub/Multispecies_model/Sim Versions/Current sims/Chpt_2/median historical exploitation")
save(res.sim.df, file = "test_TD-Baseline med u CHPT 2.RData")

save(res.sim.df, file = "test_TD-Increased med u CHPT 2.RData")

fish.all.L <- res.sim.df

fish.LPs <- res.sim.df

no.fish <- res.sim.df



com.bm <- data.frame(year = 2000:2017, bm.com = eco.tot.bm.best$bm.eco)
#get all data into tonnes
com.bm$bm.com <- com.bm$bm.com*0.001
com.bm$sim = 0
med.com.bm <- mean(com.bm$bm.com)
med.com.bm.hist <- data.frame(year = 2000:2117, med = med.com.bm)
#get medians and quantiles for each scenario
com.L <- fish.all.L |>  collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(bm.com*0.001,probs=c(0.25),na.rm=T),
                                                                           med = median(bm.com*0.001,na.rm=T),
                                                                           U.50 = quantile(bm.com*0.001,probs=c(0.75),na.rm=T))

com.LP <- fish.LPs |>  collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(bm.com*0.001,probs=c(0.25),na.rm=T),
                                                                             med = median(bm.com*0.001,na.rm=T),
                                                                             U.50 = quantile(bm.com*0.001,probs=c(0.75),na.rm=T))

com.hist <- no.fish |>  collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(bm.com*0.001,probs=c(0.25),na.rm=T),
                                                                             med = median(bm.com*0.001,na.rm=T),
                                                                             U.50 = quantile(bm.com*0.001,probs=c(0.75),na.rm=T))

com.L$type <- "fish Ls"
com.LP$type <- "fish LPs"
com.hist$type <- "no fish"

com.sim <- rbind(com.L, com.LP, com.hist)

ggplot(com.bm) + geom_line(aes(x = year, y = bm.com/1000)) +
  #geom_vline(xintercept = 2017, linetype = 2) +
  geom_line(data = med.com.bm.hist, aes(x = year, y = med/1000), linetype = 2) +
  geom_line(data = com.sim, aes(x = year, y = med/1000, colour = type)) +
  #geom_ribbon(data = com.sim, aes(x = year, ymax = U.50/1000, ymin = L.50/1000, fill = type), 
  #alpha = 0.2) +
  scale_colour_manual(values = c("darkblue", "lightblue", "orange", "gold")) +
  scale_fill_manual(values = c("darkblue", "lightblue", "orange", "gold")) +
  labs(x = "Year", y = "Community biomass \n(thousands of tonnes)",
       colour = "Scenario") +
  #facet_wrap(~type, ncol = 2, scale = "free") +
  theme(axis.title = element_text(face = "bold", size = 20),
        legend.position = "top")

#export results
setwd("C:/Users/BARIBEAUD/Desktop/GitHub/Multispecies_model/Sim Versions/Current sims/Chpt_2")
save(res.sim.df, file = "TD-Baseline CHPT 2.RData")

save(res.sim.df, file = "TD-Increased CHPT 2.RData")

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
  geom_hline(yintercept = (mean(hist.com.bm$com.bm))/1000000) +
  geom_hline(yintercept = (median(hist.com.bm$com.bm))/1000000, linetype = "dashed") +
  geom_line(data = com.K.q, aes(x = year, y = med/1000000), , colour = "lightblue") +
  geom_ribbon(data = com.K.q, aes(x = year, ymax = U.50/1000000, ymin = L.50/1000000), 
              alpha = 0.3, fill = "lightblue") +
  geom_hline(yintercept = (median(com.K.q$med))/1000000, colour = "lightblue") +
  geom_line(data = com.q, aes(x = year, y = med/1000000), , colour = "darkblue") +
  geom_ribbon(data = com.q, aes(x = year, ymax = U.50/1000000, ymin = L.50/1000000), 
              alpha = 0.1, fill = "darkblue") +
  geom_hline(yintercept = (median(com.q$med))/1000000, colour = "darkblue") +
  labs(x = "Year", y = "Community biomass (K in lightblue) \n(millions of kg)") +
  theme(axis.title = element_text(face = "bold", size = 20)) +
  ggtitle("Top-down")



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

l <- ggplot(quants[quants$fg == "L",]) +
  #geom_line(data = na.omit(res.sim.df[res.sim.df$fg == "L",]), aes(x = year, y = stock.K/1000000, group = sim), colour = "lightblue") +
  geom_line(aes(x=year, y=L.50/1000000, group=stock,), colour = "grey") +
  geom_line(aes(x = year, y = U.50/1000000, group = stock), colour = "grey") +
  geom_line(aes(x = year, y = med/1000000, group = stock)) +
  geom_line(data = p.lambdas[p.lambdas$FG %in% c("LP", "LB"),], aes(x = year, y = full.bm/1000000, group = stock), 
            colour = "blue")+
  #geom_line(data = na.omit(res.sim.df[res.sim.df$fg == "L",]), aes(x = year, y = stock.K/1000000, group = sim),
  #colour = "lightblue") +
  geom_line(data = p.lambdas[p.lambdas$FG %in% c("LP", "LB"),], aes(x = year, y = mn.in/1000000, group = stock)) +
  geom_vline(xintercept = 2017, linetype = 2, colour = "blue") +  facet_wrap(~stock, ncol = 4, scale = "free") +  #scale_y_log10(name="Biomass") + 
  labs(x = "Year", y = "Biomass (millions of kg)") +
  theme(axis.text.x = element_text(size=14, angle=45),
        axis.text.y = element_text(size=14)) +
  guides(colour = guide_legend(ncol = 5))

m <- ggplot(quants[quants$fg == "M",]) +
  #geom_line(data = na.omit(res.sim.df[res.sim.df$fg == "M",]), aes(x = year, y = stock.K/1000000, group = sim), colour = "lightblue") +
  geom_line(aes(x=year, y=L.50/1000000, group=stock,), colour = "grey") +
  geom_line(aes(x = year, y = U.50/1000000, group = stock), colour = "grey") +
  geom_line(aes(x = year, y = med/1000000, group = stock)) +
  geom_line(data = p.lambdas[p.lambdas$FG %in% c("MP", "MBZ"),], aes(x = year, y = full.bm/1000000, group = stock), 
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

ggarrange(l,m, ncol = 1, label.x = "test", label.y = "test")

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