setwd("C:/Users/BARIBEAUD/Desktop/GitHub/Multispecies_model/Sim Versions/Current sims/Chpt_3")
load("test BU CHPT 3.RData")
test.dat <- res.sim.df

#summary of what I've tried
  #simpson diversity index - not much change showing up, because all stocks are winners
    #would take a lot to shift proportions, 1% annual env dev wasn't enough
    #good QC check? Illuminates characteristics of model
  #change point analysis - not suitable for autoregressive data, which this is :(
  #runs test - also assumes independence, but could try and reduce autocorrelation by binning
    #or acknowledge problem up front

#Simpson diversity index
library(vegan)
library(tidyr)
#summary across sims
quants <- res.sim.df |>  collapse::fgroup_by(year,stock,fg) |> collapse::fsummarize(L.50 = quantile(bm.stock,probs=c(0.25),na.rm=T),
                                                                                    med = median(bm.stock,na.rm=T),
                                                                                    U.50 = quantile(bm.stock,probs=c(0.75),na.rm=T))
#get data in right format
test.df <- quants[,c("year", "stock", "med")]
l.test.df <- test.df %>% 
  pivot_wider(names_from = "stock", values_from = "med")
l.test.m <- as.matrix(l.test.df)
simpson <- diversity(l.test.m, index = "simpson")
simpson.plot <- data.frame(year = unique(res.sim.df$year), simpson = simpson)
ggplot(simpson.plot) + geom_line(aes(x = year, y = simpson))
#input data simpson
input.df <- p.lambdas[,c("year", "stock", "full.bm")]
l.input.df <- input.df %>% pivot_wider(names_from = "stock", values_from = "full.bm")
l.input.m <- as.matrix(l.input.df)
input.simpson <- diversity(l.input.m, index = "simpson")
input.simpson.plot <- data.frame(year = c(unique(p.lambdas$year)),
                                 simpson = c(input.simpson))

ggplot(simpson.plot) + geom_line(aes(x = year, y = simpson), colour = "darkblue") +
  geom_hline(yintercept = median(simpson), colour = "darkblue") +
  geom_line(data = input.simpson.plot, aes(x = year, y = simpson)) +
  geom_hline(yintercept = median(input.simpson))

#runs test
install.packages("tseries")
library(tseries)
#setup - get interannual change into binary (+/-) format
com.q <- test.dat[test.dat$sim == 17,] |>  collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(bm.com*0.001,probs=c(0.25),na.rm=T),
                                                                           med = median(bm.com*0.001,na.rm=T),
                                                                           U.50 = quantile(bm.com*0.001,probs=c(0.75),na.rm=T))
com.q$med.minus.one <- c(NA, com.q$med[1:length(com.q$med) - 1])
com.q$med.diff <- com.q$med - com.q$med.minus.one
com.q$med.diff.cat <- NA
for(i in 2:nrow(com.q)){
  if (com.q$med.diff[i] < 1) com.q$med.diff.cat[i] <- "-"
  if (com.q$med.diff[i] > 1) com.q$med.diff.cat[i] <- "+"
  if (com.q$med.diff[i] == 0) com.q$med.diff.cat[i] <- "NULL"
}
com.run <- com.q$med.diff.cat[2:length(com.q$med.diff.cat)]
com.run <- as.factor(com.run)

#perform runs test
res <- runs.test(com.run)

#try subsetting data to every 10 years
com.10.q <- com.q[c(1,seq(11,101,10)),]
com.10.q$med.minus.one <- c(NA, com.10.q$med[1:length(com.10.q$med) - 1])
com.10.q$med.diff <- com.10.q$med - com.10.q$med.minus.one
com.10.q$med.diff.cat <- NA
for (i in 2:nrow(com.10.q)){
  if (com.10.q$med.diff[i] < 1) com.10.q$med.diff.cat[i] <- "-"
  if (com.10.q$med.diff[i] > 1) com.10.q$med.diff.cat[i] <- "+"
  if (com.10.q$med.diff[i] == 0) com.10.q$med.diff.cat[i] <- "NULL"
}

dw.model <- lm(com.10.q$med.diff[2:length(com.10.q$med.diff)] ~ com.10.q$year[2:length(com.10.q$year)])
durbinWatsonTest(com.10.q$med.diff[2:length(com.10.q$med.diff)])
test <- durbinWatsonTest(dw.model)
acf(com.10.q$med.diff[2:length(com.10.q$med.diff)])

library(lmtest)
dwtest(com.10.q$med.diff[2:length(com.10.q$med.diff)])

acf(com.10.q$med.diff[2:length(com.10.q$med.diff)])

res.10 <- runs.test(as.factor(com.10.q$med.diff.cat[2:length(com.10.q$med.diff.cat)]))

for (r in 2:length(com.10.q$med.diff.cat)){
  if (com.10.q$med.diff.cat[r] == "-"){
    potential.dec <- com.10.q$med.diff.cat[r:length(com.10.q$med.diff.cat)]
    plus <- 0
    minus <- 0
    for (i in 1:length(potential.dec)){
      if(potential.dec[i] == "+") plus <- plus + 1
      if (potential.dec[i] == "-") minus <- minus + 1
    }
    minus.percent <- minus/(minus + plus)
    
    if (minus.percent > 0.5){
      dec.decade = (r - 2)*10
      print(paste("decline point is", dec.decade))
      break
    }
  }#end of minus loop
}

#t-test for determining if changes are significant
#bin 10-yr data into sections
pre <- com.q[1:21,]$med
post <- com.q[22:100,]$med
boxplot(pre, post, names = c("pre", "post"))
#variances are not equal, so use welch's t-test (unequal sample sizes as well)
t.test(pre, post, var.equal = FALSE)
#statistically significant difference! Supports that change happened around 20 years in!

#what if you just do it on one sim at a time?
sim.10 <- NULL
sum.stat
for (i in unique(res.sim.df$sim)){
  #get data of interest
  sim.dat <- res.sim.df[res.sim.df$sim == i, c("sim", "year", "bm.com")]

  #set up % change between decades
  sim.dat <- sim.dat[c(1,seq(11,101,10)),]
  sim.dat$bm.minus.10 <- c(NA, sim.dat$bm.com[1:length(sim.dat$bm.com) - 1])
  sim.dat$bm.diff <- sim.dat$bm.com - sim.dat$bm.minus.10
  dw.model <- lm(sim.dat$bm.diff ~ sim.dat$year)
  #check for autocorrelation
  dw.test <- durbinWatsonTest(dw.model)

  #categorize % change as + or -
  sim.dat$bm.diff.cat <- NA
  for (i in 2:nrow(sim.dat)){
    if (sim.dat$bm.diff[i] < 1) sim.dat$bm.diff.cat[i] <- "-"
    if (sim.dat$bm.diff[i] > 1) sim.dat$bm.diff.cat[i] <- "+"
    if (sim.dat$bm.diff[i] == 0) sim.dat$bm.diff.cat[i] <- "NULL"
  }
  
  #do runs test
  r.test <- runs.test(as.factor(sim.dat$bm.diff.cat[2:length(sim.dat$bm.diff.cat)]))
  
  if (r.test$p-value >= 0.05){
    #ID decline point
    for (r in 2:length(sim.dat$bm.diff.cat)){
      if (sim.dat$bm.diff.cat[r] == "-"){
        potential.dec <- sim.dat$bm.diff.cat[r:length(sim.dat$bm.diff.cat)]
        plus <- 0
        minus <- 0
        for (i in 1:length(potential.dec)){
          if(potential.dec[i] == "+") plus <- plus + 1
          if (potential.dec[i] == "-") minus <- minus + 1
        }
        minus.percent <- minus/(minus + plus)
        
        if (minus.percent > 0.5){
          dec.decade <- (r - 2)*10
          print(paste("decline point is", dec.decade))
          break
        }
      }#end of minus loop
      
    }
  }else {dec.decade <- 0}
  
  
  
  #check if pre- and post- decline decade are significant
  decade <- 0
  for (x in sim.dat$bm.diff.cat){
    if (sim.dat$bm.diff.cat[x] == "+") decade <- decade + 10
    if (sim.dat$bm.diff.cat[x] == "-") break
  }
  
  
  
  
  #store all summary stats
  sum.stat.tmp <- data.frame(sim = i, auto.lag1 = dw.test$r,
                             dw.stat = dw.test$dw, dw.pval$p,
                             rt.stat = r.test$statistic, rt.pval = r.test$p-value)
  
  
  sim.10[[i]] <- sim.dat
  sim.dat <- NULL
}
 
