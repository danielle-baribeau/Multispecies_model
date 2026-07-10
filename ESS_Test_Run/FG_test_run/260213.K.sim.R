#NEW VERSION
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
  
  #biomass starting values
  L.bm <- bfg.bm.best |> collapse ::fsubset(BFG == "L")
  start.L.sim <- L.bm$bm.bfg[length(L.bm$bm.bfg)]
  
  LP.bm <- mfg.bm.best |> collapse ::fsubset(BFG == "LP")
  start.LP.sim <- LP.bm$bm.bfg[length(LP.bm$bm.bfg)]
  
  LB.bm <- mfg.bm.best |> collapse ::fsubset(BFG == "LB")
  start.LB.sim <- LB.bm$bm.bfg[length(LB.bm$bm.bfg)]
  mn.LB.sim <- start.LB.sim
  
  #get all stock starting values into 1 data frame
  L.stocks <- bm.best.b |> collapse ::fsubset(BFG == "L")
  L.stocks <- unique(L.stocks$Stock)
  n.L.stocks <- length(L.stocks)
  
  start.L.stock.sim <- data.frame(stock = c(rep(0, n.L.stocks)), 
                                  start.sim = c(rep(0, n.L.stocks)), 
                                  mn = c(rep(0, n.L.stocks)),
                                  sd = c(rep(0, n.L.stocks)))
  
  stock <- NULL
  for (s in 1:length(L.stocks)){
    #browser()
    stock <- bm.best.b |> collapse::fsubset(Stock == L.stocks[s])
    start.L.stock.sim$stock[s] <- L.stocks[s]
    start.L.stock.sim$start.sim[s] <- stock$bm.stock[length(stock$bm.stock)]
    start.L.stock.sim$med[s] <- median(stock$bm.stock)
    start.L.stock.sim$sd[s] <- sd(stock$bm.stock)
    stock <- NULL
  }#end of stock value loop
  #get difference between starting value and mean for each stock
  start.L.stock.sim$start.diff <- start.L.stock.sim$start.sim - start.L.stock.sim$med

#2000-2017
  #anything where number is a proportion = "prop"
  #anything where number is biomass = "bm"
  #if "K", then it is a K
  #no K = actual biomass
  #last element tells you where you are in the hierarchy
  
  n.yrs.proj <- 100 
  bm.sim.K.com <- NULL
  bm.L.Ks <- NULL
  arima.logit <- NULL
  sim.K.stock <- NULL
  bm.sim.K.stocks <- NULL
  sim.stocks.K <- NULL
  arima.logit.stocks <- NULL
  store.logit <- NULL
  
  n.sims <- 10
  
for (i in 1:n.sims){
  #browser()
  #size of community is autocorrelated
  bm.sim.K.com[[i]] <- data.frame(bm = c(arima.sim(model =list(ar = com.cor.lag.1),n = n.yrs.proj,n.start=1,start.innov = start.com.diff/com.cor.lag.1,
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
  #this is the K values for the larges
  bm.sim.K.L <- bm.sim.K.com[[i]]$bm * prop.sim.K.L
  
  #storing L K values:
  bm.L.Ks[[i]] <- data.frame(Years = rep(1:n.yrs.proj), sim =i,
                              bm.L = c(bm.sim.K.L),
                              fg = as.factor(sort(rep(c("L"),n.yrs.proj))),
                              bm.com = rep(bm.sim.K.com[[i]]$bm))
  bm.L.Ks[[i]]$prop.com.fg <- bm.L.Ks[[i]]$bm.L/bm.L.Ks[[i]]$bm.com
  
  #storing results of prop.sim.K.L arima
  if (i == 1){
    arima.logit <- data.frame(sim = rep(i, n.yrs.proj),
                              L.com = as.vector(prop.sim.K.L))
  }
  if (i > 1){
    arima.logit <- rbind(arima.logit, data.frame(sim = rep(i, n.yrs.proj),
                                                 L.com = as.vector(prop.sim.K.L)))
  }
  
  #now going to partition L K to all the L stocks:
  #set up object for storing ARIMA sim results
  stocks.arima.logit <- data.frame(stock = NA, sim = NA, stock.L = NA)
  #set up counter
  count = 0
  for (f in L.stocks){
    #browser()
    count = count + 1

    if(count == 1 || n.L.stocks != 2){
      #browser()
      tmp.dat <- bm.best.b |> collapse::fsubset(BFG == "L")
      tmp.dat <- tmp.dat |> collapse::fsubset(Stock == f)
      tmp.prop.L.cor <- pacf(tmp.dat$prop.bm.stock.bfg, plot = F) #using lags 1 and 2
      tmp.prop.L.cor.lag.1 <- tmp.prop.L.cor$acf[1]
      tmp.prop.L.cor.lag.2 <- tmp.prop.L.cor$acf[2]
      #logit-transform the proportions and do the ARIMA on the logits
      tmp.prop.L.logit <- logit(tmp.dat$prop.bm.stock.bfg)
      start.tmp.prop.L.logit <- tmp.prop.L.logit[length(tmp.prop.L.logit)]
      mn.tmp.prop.L.logit <- start.tmp.prop.L.logit
      sd.tmp.prop.L.logit <- sd(tmp.prop.L.logit)
      diff.tmp.prop.L.logit <- start.tmp.prop.L.logit - mn.tmp.prop.L.logit
      
      #now start the ARIMA
      prop.sim.K.L.stock <- inv.logit(arima.sim(model =list(ar = c(tmp.prop.L.cor.lag.1)),n = n.yrs.proj,
                                              n.start =1, start.innov = c(start.tmp.prop.L.logit/tmp.prop.L.cor.lag.1), 
                                              innov = c(0,rnorm(n.yrs.proj-1,0,sd.tmp.prop.L.logit))) +  mn.tmp.prop.L.logit)

      
      if (f == 10){
        store.logit <- data.frame(Year = c(seq(1, n.yrs.proj)), sim = c(rep(i, n.yrs.proj)), 
                                  stock = c(rep(f,n.yrs.proj)), prop.sim.K.L.stock =
                                    as.vector(prop.sim.K.L.stock))
      }
      if (f != 10){
        store.logit <- rbind(store.logit, data.frame(Year = c(seq(1, n.yrs.proj)), sim = c(rep(i, n.yrs.proj)), 
                                                     stock = c(rep(f,n.yrs.proj)), prop.sim.K.L.stock =
                                                       as.vector(prop.sim.K.L.stock)))
      }
      #get the K values for the stock in terms of biomass
      bm.sim.K.stock[[f]] <- data.frame(Years = 1:n.yrs.proj, sim = i,
                                Stock = f,
                                bm.stock = as.vector(prop.sim.K.L.stock)*bm.L.Ks[[i]]$bm.L)
      
    }

  }#end of stock K loop
  if (i == 1){
    arima.logit.stocks <- store.logit
  }
  if (i > 1) {
    arima.logit.stocks <- rbind(arima.logit.stocks, store.logit)
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
  sim.K.stock[[i]] <- do.call("rbind",bm.sim.K.stock)
}
  sim.stocks.K <- do.call("rbind",sim.K.stock)
  sim.L.K <- do.call("rbind",bm.L.Ks)
  sim.com.K <- do.call("rbind",bm.sim.K.com)

  eco.plt <- sim.com.K
  eco.plt$type <- rep("Simulated", nrow(eco.plt))
  eco.input <- data.frame (bm = eco.tot.bm.best$bm.eco, Years = eco.tot.bm.best$Year,
                           sim = c(rep(0, nrow(eco.tot.bm.best))), type = c(rep("Inputs", nrow(eco.tot.bm.best))))
  
  eco.bm.plt <- rbind(eco.plt, eco.input)
  eco.bm.plt$bm <- eco.bm.plt$bm/1000000
  eco.bm.plt$sim <- as.factor(eco.bm.plt$sim)
  
  sim.eco.bm.plt <- ggplot(eco.bm.plt, aes(x = sim, y = bm, fill = type)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white", color = "black") +
    labs(x="Simulation iteration", y = "Community biomass (millions)", fill = "Data type") +
    scale_fill_manual(values = c("darkgrey", "lightblue"))
  sim.eco.bm.plt
  
  sim.eco.bm.plt.2 <- ggplot(eco.bm.plt, aes(x = type, y = bm, fill = type)) +
    geom_violin() +
    geom_boxplot(width = 0.2, fill = "white", color = "black") +
    labs(x="Data type", y = "Community biomass (millions)") +
    scale_fill_manual(values = c("darkgrey", "lightblue")) +
    theme(legend.position = "none")
  sim.eco.bm.plt.2
  
  #repeat with logit arimas
  L.com.arima <- arima.logit
  L.com.arima$type <- c(rep("Simulated", nrow(L.com.arima)))
  L.com.input <- data.frame(sim = c(rep(0,length(L.com.prop))), L.com = L.com.prop, 
                            type = c(rep("Input", length(L.com.prop))))
  L.com.plt <- rbind(L.com.arima, L.com.input)
  L.com.plt$sim <- as.factor(L.com.plt$sim)
  
  sim.plt <- ggplot(L.com.plt, aes(x = sim, y = L.com, fill = type)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white", color = "black") +
    labs(x="Simulation iteration", y = "L/com", fill = "Data type") +
    scale_fill_manual(values = c("darkgrey", "lightblue"))
  sim.plt
  
  sim.plt.2 <- ggplot(L.com.plt, aes(x = type, y = L.com, fill = type)) +
    geom_violin() +
    geom_boxplot(width = 0.2, fill = "white", color = "black") +
    labs(x="Data type", y = "L/com") +
    scale_fill_manual(values = c("darkgrey", "lightblue")) +
    theme(legend.position = "none")
  sim.plt.2
  
  #repeat with stock arimas
  pdf(file = "C:/Users/BARIBEAUD/Desktop/GitHub/ESS_Test_Run/FG_test_run/260220/K sim/test.stock.arimas.K.pdf",
      height = 12, width = 24)
  for (i in L.stocks){
    #browser()
    stock.arima <- arima.logit.stocks |> collapse::fsubset(stock == i)
    stock.arima$type <- c(rep("Simulated", nrow(stock.arima)))
    stock.data <- bm.best.b |> collapse::fsubset(Stock == i)
    stock.inputs <- data.frame(Year = yrs, sim = c(rep(0, length(yrs))), stock = c(rep(i, length(yrs))),
                               prop.sim.K.L.stock = stock.data$prop.bm.stock.bfg,
                               type = c(rep("Input", length(yrs))))
    stock.plt <- rbind(stock.arima, stock.inputs)
    stock.plt$sim <- as.factor(stock.plt$sim)
    
    plot <- ggplot(stock.plt, aes(x = sim, y = prop.sim.K.L.stock, fill = type)) +
      geom_violin() +
      geom_boxplot(width = 0.1, fill = "white", color = "black") +
      labs(x="Simulation iteration", y = "stock/L", fill = "Data type") +
      scale_fill_manual(values = c("darkgrey", "lightblue")) +
      ggtitle(paste(as.character(i)))
    print(plot)

  }
  dev.off()
  
  pdf(file = "C:/Users/BARIBEAUD/Desktop/GitHub/ESS_Test_Run/FG_test_run/260220/K sim/stock.arimas.K.overall.pdf",
      height = 12, width = 12)
  for (i in L.stocks){
    #browser()
    stock.arima <- arima.logit.stocks |> collapse::fsubset(stock == i)
    stock.arima$type <- c(rep("Simulated", nrow(stock.arima)))
    stock.data <- bm.best.b |> collapse::fsubset(Stock == i)
    stock.inputs <- data.frame(Year = yrs, sim = c(rep(0, length(yrs))), stock = c(rep(i, length(yrs))),
                               prop.sim.K.L.stock = stock.data$prop.bm.stock.bfg,
                               type = c(rep("Input", length(yrs))))
    stock.plt <- rbind(stock.arima, stock.inputs)
    stock.plt$sim <- as.factor(stock.plt$sim)
    
    plot <- ggplot(stock.plt, aes(x = type, y = prop.sim.K.L.stock, fill = type)) +
      geom_violin() +
      geom_boxplot(width = 0.2, fill = "white", color = "black") +
      labs(x="Data type", y = "stock/L") +
      scale_fill_manual(values = c("darkgrey", "lightblue")) +
      theme(legend.position = "none") +
      ggtitle(paste(as.character(i)))
    print(plot)
    
  }
  dev.off()  
  
  monkfish <- bm.best.b |> collapse::fsubset(Stock == 400)
  pacf(monkfish$prop.bm.stock.bfg)
  #so significant autocorrelation in monkfish input time series
  sim1 <- arima.logit.stocks |> collapse::fsubset(sim == 1)
  year1 <- sim1|> collapse::fsubset(Year == 1)
  pie(year1$prop.sim.K.L.stock, year1$stock)
  year1.in <- bm.best.b |> collapse::fsubset(BFG == "L") |> collapse::fsubset(Year == 2000)
  pie(year1.in$prop.bm.stock.bfg, year1.in$Stock)
  
  
  #compare biomass from simulation to input data
  #ecosystem biomass
  #get quantiles
  q.sim.eco <- sim.eco.K |>  collapse::fgroup_by(Years) |> collapse::fsummarize(L.50 = quantile(bm,probs=c(0.25),na.rm=T),
                                                                                med = median(bm,na.rm=T),
                                                                                U.50 = quantile(bm,probs=c(0.75),na.rm=T)) |> ungroup()
  #make projection years into actual years
  for (i in 1:nrow(q.sim.eco)){
    q.sim.eco$Years[i] <- q.sim.eco$Years[i] + 2017
  }
  #set up input time series to match layout of projection quantiles
  q.eco.inputs <- data.frame(Years = eco.input$Years, type = eco.input$type, bm = com.bm)
  #put biomass into a readable format
  q.sim.eco$L.50 <- q.sim.eco$L.50/1000000
  q.sim.eco$med <- q.sim.eco$med/1000000
  q.sim.eco$U.50 <- q.sim.eco$U.50/1000000
  q.eco.inputs$bm <- q.eco.inputs$bm/1000000
  #plot quantiles and input time series
  q.eco.sim.plt <- ggplot(q.eco.inputs, aes(x = Years, y = bm)) +
    geom_line() +
    geom_line(data=q.sim.eco, aes(x=Years, y=L.50), linetype = 2, colour = "darkgrey") +
    geom_line(data=q.sim.eco, aes(x=Years, y=med), linetype = 2) +
    geom_line(data=q.sim.eco, aes(x=Years, y=U.50), linetype = 2, colour = "darkgrey") +
    scale_x_continuous(breaks = seq(2000, 2117, by = 20)) +
    labs(x = "Year", y = "Community K (biomass in millions)")
  q.eco.sim.plt
  
#functional groups
  q.sim.fg <- sim.fg.K |>  collapse::fsubset(fg == "L") |>
    collapse::fgroup_by(Years) |> collapse::fsummarize(L.50 = quantile(bm.fg,probs=c(0.25),na.rm=T),
                                                       med = median(bm.fg,na.rm=T),
                                                       U.50 = quantile(bm.fg,probs=c(0.75),na.rm=T)) |> ungroup()
  #make projection years into actual years
  for (i in 1:nrow(q.sim.fg)){
    q.sim.fg$Years[i] <- q.sim.fg$Years[i] + 2017
  }
  #set up input time series to match layout of projection quantiles
  fg.bm <- fg.bm.best |> collapse::fsubset(FG == "L")
  fg.inputs <- data.frame(Years = eco.input$Years, type = eco.input$type, bm = L.bm)
  #put biomass into a readable format
  q.sim.fg$L.50 <- q.sim.fg$L.50/1000000
  q.sim.fg$med <- q.sim.fg$med/1000000
  q.sim.fg$U.50 <- q.sim.fg$U.50/1000000
  fg.inputs$bm <- fg.inputs$bm/1000000
  #plot quantiles and input time series
  q.fg.sim.plt <- ggplot(fg.inputs, aes(x = Years, y = bm)) +
    geom_line() +
    geom_line(data=q.sim.fg, aes(x=Years, y=L.50), linetype = 2, colour = "darkgrey") +
    geom_line(data=q.sim.fg, aes(x=Years, y=med), linetype = 2) +
    geom_line(data=q.sim.fg, aes(x=Years, y=U.50), linetype = 2, colour = "darkgrey") +
    scale_x_continuous(breaks = seq(2000, 2027, by = 2)) +
    labs(x = "Year", y = "Large K (biomass in millions)")
  q.fg.sim.plt
  
  
  

