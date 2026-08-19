#July 23:
  #add up actual community resource pool for all stocks
  #make rest of results plots and write up outlines


#Making plots of multiple scenarios
library(ggplot2)
library(cowplot)

#Main text plots
  #Community biomass time series (4 lines on one plot)
  #FG biomass time series (4 lines per plot, 2 panels)
  #RP versus biomass time series for community (4 panels, 1 per scenario)
  #Stock-specific percent decline every 10 years (heatmap, 4 figures)
  #Violin of scenario lambdas - compare observed versus simulated

#Supplementary plots
  #Stock-specific:
    #biomass time series
    #raw sim biomasses
    #resource pool time series
    #raw resource pool biomasses

#### Load data and set things up ####

#load simulation outputs
setwd("C:/Users/BARIBEAUD/Desktop/GitHub/Multispecies_model/Sim Versions/Current sims/Chpt_2/median historical exploitation")
load("BU-Increased med u CHPT 2.RData")
bu.inc <- res.sim.df
load("BU-Baseline med u CHPT 2.RData")
bu.hist <- res.sim.df
load("TD-Baseline med u CHPT 2.RData")
td.hist <- res.sim.df
load("TD-Increased med u CHPT 2.RData")
td.inc <- res.sim.df

#add scenario labels
bu.inc$scenario <- "BU-Increased"
bu.hist$scenario <- "BU-Baseline"
td.hist$scenario <- "TD-Baseline"
td.inc$scenario <- "TD-Increased"

#store together
res.all <- rbind(bu.inc, bu.hist, td.hist, td.inc)

#get stock names
M.stocks <- unique(bu.inc$stock[bu.inc$fg == "M"])
L.stocks <- unique(bu.inc$stock[bu.inc$fg == "L"])
LB.stocks <- c("Haddock", "Atlantic wolffish", "Thorny skate")
LP.stocks <- c("Atlantic cod", "White hake", "Pollock", "American plaice", "Monkfish")

##### Community time series #####
  ##community removals
      #get historical removals
      r.com.hist <- p.lambdas |> group_by(year) |> summarize(tot.rem = sum(removals)) |> ungroup()
      #convert to tonnes
      r.com.hist$tot.rem <- r.com.hist$tot.rem*0.001
      med.rcom.hist <- median(r.com.hist$tot.rem)
      med.rcom.bm.hist <- data.frame(year = 2000:2117, med = med.rcom.hist)
      
    #get yearly removals per sim
  rcom.bu.inc <- bu.inc |>  collapse::fgroup_by(sim, year) |> collapse::fsummarize(tot.rem = sum(removals)) |> ungroup()
  rcom.bu.hist <- bu.hist |>  collapse::fgroup_by(sim, year) |> collapse::fsummarize(tot.rem = sum(removals)) |> ungroup()
  rcom.td.inc <- td.inc |>  collapse::fgroup_by(sim, year) |> collapse::fsummarize(tot.rem = sum(removals)) |> ungroup()
  rcom.td.hist <- td.hist |>  collapse::fgroup_by(sim, year) |> collapse::fsummarize(tot.rem = sum(removals)) |> ungroup()
    #summary stats across sims
  r.bu.inc <- rcom.bu.inc |>  collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(tot.rem*0.001,probs=c(0.25),na.rm=T),
                                                                                med = median(tot.rem*0.001,na.rm=T),
                                                                                U.50 = quantile(tot.rem*0.001,probs=c(0.75),na.rm=T))
  r.bu.hist <- rcom.bu.hist |>  collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(tot.rem*0.001,probs=c(0.25),na.rm=T),
                                                                                med = median(tot.rem*0.001,na.rm=T),
                                                                                U.50 = quantile(tot.rem*0.001,probs=c(0.75),na.rm=T)) 
  r.td.inc <- rcom.td.inc |>  collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(tot.rem*0.001,probs=c(0.25),na.rm=T),
                                                                                med = median(tot.rem*0.001,na.rm=T),
                                                                                U.50 = quantile(tot.rem*0.001,probs=c(0.75),na.rm=T))
  r.td.hist <- rcom.td.hist |>  collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(tot.rem*0.001,probs=c(0.25),na.rm=T),
                                                                                med = median(tot.rem*0.001,na.rm=T),
                                                                                U.50 = quantile(tot.rem*0.001,probs=c(0.75),na.rm=T))
  #add scenarios    
      r.bu.inc$type <- "BU-Increased"
      r.bu.hist$type <- "BU-Baseline"
      r.td.hist$type <- "TD-Baseline"
      r.td.inc$type <- "TD-Increased"
      
      r.sim <- rbind(r.bu.inc, r.bu.hist, r.td.hist, r.td.inc)
      #get right years (plotting so that year matches the removals used to estimate biomass in that year)
        #this way, we don't have 2 sets of 2017 removals on plot
      r.sim$year <- r.sim$year + 1
      
    #plot!
      ggplot(r.com.hist) + geom_line(aes(x = year, y = tot.rem/1000)) +
        #geom_vline(xintercept = 2017, linetype = 2) +
        geom_line(data = med.rcom.bm.hist, aes(x = year, y = med/1000), linetype = 2) +
        geom_line(data = r.sim, aes(x = year, y = med/1000, colour = type)) +
        #geom_ribbon(data = com.sim, aes(x = year, ymax = U.50/1000, ymin = L.50/1000, fill = type), 
        #alpha = 0.2) +
        scale_colour_manual(values = c("darkblue", "lightblue", "orange", "gold")) +
        scale_fill_manual(values = c("darkblue", "lightblue", "orange", "gold")) +
        labs(x = "Year", y = "Community removals \n(thousands of tonnes)",
             colour = "Scenario") +
        #facet_wrap(~type, ncol = 2, scale = "free") +
        theme(axis.title = element_text(face = "bold", size = 20),
              legend.position = "top")
      
 ##community biomass    
      #get inputs
    com.bm <- data.frame(year = 2000:2017, bm.com = eco.tot.bm.best$bm.eco)
    #get all data into tonnes
    com.bm$bm.com <- com.bm$bm.com*0.001
    com.bm$sim = 0
    med.com.bm <- median(com.bm$bm.com)
    med.com.bm.hist <- data.frame(year = 2000:2117, med = med.com.bm)
    #get medians and quantiles for each scenario
    com.bu.inc <- bu.inc |>  collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(bm.com*0.001,probs=c(0.25),na.rm=T),
                                                                              med = median(bm.com*0.001,na.rm=T),
                                                                              U.50 = quantile(bm.com*0.001,probs=c(0.75),na.rm=T))
    
    com.bu.hist <- bu.hist |>  collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(bm.com*0.001,probs=c(0.25),na.rm=T),
                                                                              med = median(bm.com*0.001,na.rm=T),
                                                                              U.50 = quantile(bm.com*0.001,probs=c(0.75),na.rm=T))
    
    com.td.hist <- td.hist |>  collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(bm.com*0.001,probs=c(0.25),na.rm=T),
                                                                                 med = median(bm.com*0.001,na.rm=T),
                                                                                 U.50 = quantile(bm.com*0.001,probs=c(0.75),na.rm=T))
    
    com.td.inc <- td.inc |>  collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(bm.com*0.001,probs=c(0.25),na.rm=T),
                                                                                 med = median(bm.com*0.001,na.rm=T),
                                                                                 U.50 = quantile(bm.com*0.001,probs=c(0.75),na.rm=T))
    #join into one data frame
    com.bu.inc$type <- "BU-Increased"
    com.bu.hist$type <- "BU-Baseline"
    com.td.hist$type <- "TD-Baseline"
    com.td.inc$type <- "TD-Increased"
    
    com.sim <- rbind(com.bu.inc, com.bu.hist, com.td.hist, com.td.inc)
    
    #plot!
    ggplot(com.bm) + geom_line(aes(x = year, y = bm.com/1000)) +
      #geom_vline(xintercept = 2017, linetype = 2) +
      geom_line(data = med.com.bm.hist, aes(x = year, y = med/1000), linetype = 2) +
      geom_line(data = com.sim, aes(x = year, y = med/1000, colour = type)) +
      #geom_ribbon(data = com.sim, aes(x = year, ymax = U.50/1000, ymin = L.50/1000, fill = type), 
                  #alpha = 0.2) +
      scale_colour_manual(values = c("darkblue", "lightblue", "orange", "gold")) +
      scale_fill_manual(values = c("darkblue", "lightblue", "orange", "gold")) +
      labs(x = "Year", y = "Community biomass \n(thousands of tonnes)",
           colour = "Scenario", fill = "Scenario") +
      #facet_wrap(~type, ncol = 2, scale = "free") +
      theme(axis.title = element_text(face = "bold", size = 20),
            legend.position = "top")
    
    
  ##community resource pools versus biomass
    #get yearly FG RPs per sim
    fgrp.bu.inc <- unique(bu.inc[,c(1,2,4,12,17)])
    fgrp.bu.hist <- unique(bu.hist[,c(1,2,4,12,17)])
    fgrp.td.inc <- unique(td.inc[,c(1,2,4,12,17)])
    fgrp.td.hist <- unique(td.hist[,c(1,2,4,12,17)])
    
    #get community RP by adding RPs for L and M together
    crp.bu.inc <- fgrp.bu.inc |>  collapse::fgroup_by(sim, year) |> collapse::fsummarize(com.rp = sum(fg.K)) |> ungroup()
    crp.bu.hist <- fgrp.bu.hist |>  collapse::fgroup_by(sim, year) |> collapse::fsummarize(com.rp = sum(fg.K)) |> ungroup()
    crp.td.inc <- fgrp.td.inc |>  collapse::fgroup_by(sim, year) |> collapse::fsummarize(com.rp = sum(fg.K)) |> ungroup()
    crp.td.hist <- fgrp.td.hist |>  collapse::fgroup_by(sim, year) |> collapse::fsummarize(com.rp = sum(fg.K)) |> ungroup()
    
    #summary stats across sims (in tonnes)
    com.rp.bu.inc <- crp.bu.inc |>  collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(com.rp*0.001,probs=c(0.25),na.rm=T),
                                                                               med = median(com.rp*0.001,na.rm=T),
                                                                               U.50 = quantile(com.rp*0.001,probs=c(0.75),na.rm=T))
    
    com.rp.bu.hist <- crp.bu.hist |>  collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(com.rp*0.001,probs=c(0.25),na.rm=T),
                                                                                      med = median(com.rp*0.001,na.rm=T),
                                                                                      U.50 = quantile(com.rp*0.001,probs=c(0.75),na.rm=T))
    
    com.rp.td.inc <- crp.td.inc |>  collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(com.rp*0.001,probs=c(0.25),na.rm=T),
                                                                                      med = median(com.rp*0.001,na.rm=T),
                                                                                      U.50 = quantile(com.rp*0.001,probs=c(0.75),na.rm=T))
    
    com.rp.td.hist <- crp.td.hist |>  collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(com.rp*0.001,probs=c(0.25),na.rm=T),
                                                                                      med = median(com.rp*0.001,na.rm=T),
                                                                                      U.50 = quantile(com.rp*0.001,probs=c(0.75),na.rm=T))
    #join into one data frame
    com.rp.bu.inc$type <- "BU-Increased"
    com.rp.bu.hist$type <- "BU-Baseline"
    com.rp.td.hist$type <- "TD-Baseline"
    com.rp.td.inc$type <- "TD-Increased"
    
    com.rp.sim <- rbind(com.rp.bu.inc, com.rp.bu.hist, com.rp.td.hist, com.rp.td.inc)
    
    com.rp.sim$bm <- "Resource pool"
    com.sim$bm <- "Realized biomass"
    com.rp.bm.sim <- rbind(com.sim, com.rp.sim)
    
    #look at proportion of RP used by community each year
    prop.com.sim <- com.sim[,c(1,3,5)]
    colnames(prop.com.sim)[colnames(prop.com.sim) == "med"] <- "com.bm"
    prop.com.rp.sim <- com.rp.sim[,c(1,3,5)]
    colnames(prop.com.rp.sim)[colnames(prop.com.rp.sim) == "med"] <- "rp.bm"
    prop.com.rp.bm.sim <- merge(prop.com.sim, prop.com.rp.sim, by = c("year", "type"))
    prop.com.rp.bm.sim$prop.rp <- prop.com.rp.bm.sim$com.bm/prop.com.rp.bm.sim$rp.bm
    
    #summarize proportion of RP used for each scenario
    sum.prop.rp.com <- prop.com.rp.bm.sim |> group_by(type) |> summarize(med.prop.rp = median(prop.rp, na.rm = T)) |> ungroup()
    
    #plot!
    ggplot(com.rp.bm.sim) +
      geom_line(data = com.bm, aes(x = year, y = bm.com/1000)) +
      geom_line(aes(x = year, y = med/1000, colour = bm)) +
      scale_colour_manual(values = c("#69757D", "#AFB3B3")) +
      geom_line(data = med.com.bm.hist, aes(x = year, y = med/1000), linetype = 2) +
      labs(x = "Year", y = "Community size \n(thousands of tonnes)",
           colour = "Type") +
      facet_wrap(~type, scale = "free") +
      theme(axis.title = element_text(face = "bold", size = 20),
            legend.position = "top")
      

    
    ##### Functional groups #####
  ##functional group biomass time series
    #get inputs
    fg.bm.hist <- unique(bm.best.m[,c(2,3,6)])
    #get all data into tonnes
    fg.bm.hist$bm.mfg <- fg.bm.hist$bm.mfg*0.001
    colnames(fg.bm.hist)[colnames(fg.bm.hist) == "MFG"] <- "fg"
    #get median historical bm
    med.fg.hist <- fg.bm.hist |> group_by(fg) |> summarize(med = median(bm.mfg, na.rm = T)) |> ungroup()
    fg <- c("LP", "LB", "M")
    med.fg.bm.hist <- NULL
    for (f in fg){
      med.fg.tmp <- data.frame(year = 2000:2117, med = med.fg.hist$med[med.fg.hist$fg == f], fg = f)
      med.fg.bm.hist[[f]] <- med.fg.tmp
    }#end of historical median setup loop
    med.fg.bm <- do.call("rbind", med.fg.bm.hist)
    
    #get LP functional group ID into results (doing removals and RP at same time for later)
    mfg.res.LB <- res.all[res.all$stock %in% c("Haddock", "Atlantic wolffish", "Thorny skate"),]
    bm.LB <- mfg.res.LB |> collapse::fgroup_by(scenario, sim, year) |> collapse::fsummarize(bm = sum(bm.stock, na.rm = T),
                                                                                            rem = sum(removals, na.rm = T),
                                                                                            rp = sum(stock.K, na.rm = T)) |> ungroup()
    bm.LB$fg <- "LB"
    
    mfg.res.LP <- res.all[res.all$stock %in% LP.stocks,]
    bm.LP <- mfg.res.LP |> collapse::fgroup_by(scenario, sim, year) |> collapse::fsummarize(bm = sum(bm.stock, na.rm = T),
                                                                                            rem = sum(removals, na.rm = T),
                                                                                            rp = sum(stock.K, na.rm = T)) |> ungroup()
    bm.LP$fg <- "LP"
    
    mfg.res.M <- res.all[res.all$fg == "M",]
    bm.M <- mfg.res.M |> collapse::fgroup_by(scenario, sim, year) |> collapse::fsummarize(bm = sum(bm.stock, na.rm = T),
                                                                                          rem = sum(removals, na.rm = T),
                                                                                          rp = sum(stock.K, na.rm = T)) |> ungroup()
    bm.M$fg <- "M"
    
    mfg.res <- rbind(bm.LP, bm.LB, bm.M)
    
    #get medians and quantiles for each fg in each scenario (in tonnes)
    q.mfg.res <- mfg.res |> collapse::fgroup_by(scenario, year, fg) |> collapse::fsummarize(L.50 = quantile(bm*0.001,probs=c(0.25),na.rm=T),
                                                                                            med = median(bm*0.001,na.rm=T),
                                                                                            U.50 = quantile(bm*0.001,probs=c(0.75),na.rm=T))
    
    #plot!
    ggplot(fg.bm.hist) + geom_line(aes(x = Year, y = bm.mfg/1000, group = fg)) +
      geom_line(data = med.fg.bm, aes(x = year, y = med/1000, group = fg), linetype = 2) +
      geom_line(data = q.mfg.res, aes(x = year, y = med/1000, colour = scenario)) +
      scale_colour_manual(values = c("darkblue", "lightblue", "orange", "gold")) +
      scale_fill_manual(values = c("darkblue", "lightblue", "orange", "gold")) +
      labs(x = "Year", y = "Biomass (thousands of tonnes)",
           colour = "Scenario") +
      scale_y_log10() +
      facet_wrap(~fg, ncol = 3) +
      theme(axis.title = element_text(face = "bold", size = 20),
            legend.position = "top")
    
    
    #test - area chart
    q.mfg.res$fg <- as.factor(q.mfg.res$fg)
    levels(q.mfg.res$fg) <- c("LP", "LB", "M")
    fg.bm.hist$fg <- as.factor(fg.bm.hist$fg)
    levels(fg.bm.hist$fg) <- c("LP", "LB", "M")
    
    ggplot(q.mfg.res, aes(x = year, y = med/1000, fill = fg)) +
      geom_area(alpha = 0.6) +
      geom_area(data = fg.bm.hist, aes(x = Year, y = bm.mfg/1000, fill = fg)) +
      geom_hline(yintercept = unique(med.com.bm.hist$med/1000), linetype = 2, size = 0.7) +
      #geom_vline(xintercept = 2017, linetype = 1, size = 0.7, colour = "darkgrey") +
      scale_fill_manual(values = c("orange", "gold", "lightblue")) +
      labs(x = "Year", y = "Biomass (thousands of tonnes)",
           fill = "Functional group") +
      facet_wrap(~scenario, scale = "free") +
      theme(panel.background = element_rect("white"),
            legend.position = "top")
    
    ##functional group removals
    #get input fg removals in tonnes
    rem.LB.stocks <- lambdas[lambdas$common %in% c("Haddock", "Atlantic wolffish", "Thorny skate"),]
    rem.LB.stocks$removals <- rem.LB.stocks$removals*0.001
    rem.LB.hist <- rem.LB.stocks |> group_by(year) |> summarize(fg.rem = sum(removals, na.rm = T)) |> ungroup()
    rem.LB.hist$fg <- "LB"
    
    rem.LP.stocks <- lambdas[lambdas$common %in% LP.stocks,]
    rem.LP.stocks$removals <- rem.LP.stocks$removals*0.001
    rem.LP.hist <- rem.LP.stocks |> group_by(year) |> summarize(fg.rem = sum(removals, na.rm = T)) |> ungroup()
    rem.LP.hist$fg <- "LP"
    
    rem.M.stocks <- lambdas[lambdas$common %in% M.stocks,]
    rem.M.stocks$removals <- rem.M.stocks$removals*0.001
    rem.M.hist <- rem.M.stocks |> group_by(year) |> summarize(fg.rem = sum(removals, na.rm = T)) |> ungroup()
    rem.M.hist$fg <- "M"
    
    fg.rem.hist <- rbind(rem.LB.hist, rem.LP.hist, rem.M.hist)
    
    #get median historical removals
    med.rem <- fg.rem.hist |> group_by(fg) |> summarize(med = median(fg.rem, na.rm = T)) |> ungroup()
    fg <- c("LP", "LB", "M")
    med.fg.rem.hist <- NULL
    for (f in fg){
      med.fg.tmp <- data.frame(year = 2000:2117, med = med.rem$med[med.fg.rem$fg == f], fg = f)
      med.fg.rem.hist[[f]] <- med.fg.tmp
    }#end of historical median setup loop
    med.fg.rem <- do.call("rbind", med.fg.rem.hist)
    
    #get medians and quantiles for each fg in each scenario (in tonnes) (got removals for each FG summarized in bm section)
    q.fg.rem <- mfg.res |> collapse::fgroup_by(scenario, year, fg) |> collapse::fsummarize(L.50 = quantile(rem*0.001,probs=c(0.25),na.rm=T),
                                                                                            med = median(rem*0.001,na.rm=T),
                                                                                            U.50 = quantile(rem*0.001,probs=c(0.75),na.rm=T))
    
    #get to right years
    q.fg.rem <- q.fg.rem[q.fg.rem$year != 2117,]
    q.fg.rem$year <- q.fg.rem$year + 1
    
    #plot!
    ggplot(fg.rem.hist) + geom_line(aes(x = year, y = fg.rem/1000, group = fg)) +
      geom_line(data = med.fg.rem, aes(x = year, y = med/1000, group = fg), linetype = 2) +
      geom_line(data = q.fg.rem, aes(x = year, y = med/1000, colour = scenario)) +
      scale_colour_manual(values = c("darkblue", "lightblue", "orange", "gold")) +
      scale_fill_manual(values = c("darkblue", "lightblue", "orange", "gold")) +
      labs(x = "Year", y = "Removals (thousands of tonnes)",
           colour = "Scenario") +
      scale_y_log10() +
      facet_wrap(~fg, ncol = 3) +
      theme(axis.title = element_text(face = "bold", size = 20),
            legend.position = "top")
    
    
  ##FG resource pools versus biomass [currently not using]
      #not super informative and can't get faceting right
    #get medians and quantiles for each fg in each scenario (in tonnes) (got RP and bm for each FG summarized in bm section)
    q.fg.bm <- q.mfg.res
    q.fg.bm$bm <- "Realized biomass"
    
    q.fg.rp <- mfg.res |> collapse::fgroup_by(scenario, year, fg) |> collapse::fsummarize(L.50 = quantile(rp*0.001,probs=c(0.25),na.rm=T),
                                                                                           med = median(rp*0.001,na.rm=T),
                                                                                           U.50 = quantile(rp*0.001,probs=c(0.75),na.rm=T))
    q.fg.rp$bm <- "Resource pool"

    fg.rp.bm.sim <- rbind(q.fg.rp, q.fg.bm)
    
    #look at proportion of RP used by community each year
    prop.fg.sim <- q.fg.bm[,c(1,2,3,5)]
    colnames(prop.fg.sim)[colnames(prop.fg.sim) == "med"] <- "fg.bm"
    prop.fg.rp.sim <- q.fg.rp[,c(1,2,3,5)]
    colnames(prop.fg.rp.sim)[colnames(prop.fg.rp.sim) == "med"] <- "rp.bm"
    prop.fg.rp.bm.sim <- merge(prop.fg.sim, prop.fg.rp.sim, by = c("year", "scenario", "fg"))
    prop.fg.rp.bm.sim$prop.rp <- prop.fg.rp.bm.sim$fg.bm/prop.fg.rp.bm.sim$rp.bm
    
    #summarize proportion of RP used for each scenario
    sum.prop.rp.fg <- prop.fg.rp.bm.sim |> collapse::fgroup_by(scenario, fg) |> collapse::fsummarize(med.prop.rp = median(prop.rp, na.rm = T)) |> ungroup()
    
    #plot!
    ggplot(fg.rp.bm.sim) +
      geom_line(data = fg.bm.hist, aes(x = Year, y = bm.mfg/1000, group = fg)) +
      geom_line(aes(x = year, y = med/1000, colour = bm)) +
      scale_colour_manual(values = c("#69757D", "#AFB3B3")) +
      geom_line(data = med.fg.bm, aes(x = year, y = med/1000, group = fg), linetype = 2) +
      labs(x = "Year", y = "Functional group size \n(thousands of tonnes)",
           colour = "Type") +
      facet_wrap(~interaction(scenario, fg , sep = " ",drop=T),scales='free_y') +
      theme(axis.title = element_text(face = "bold", size = 20),
            legend.position = "top")    
    
    ##### Stocks #####
    
    
### Stock biomasses/removals per scenario
    quants <- res.all |> collapse::fgroup_by(scenario, stock, year) |> collapse::fsummarize(L.50 = quantile(bm.stock*0.001,probs=c(0.25),na.rm=T),
                                                                                            med = median(bm.stock*0.001,na.rm=T),
                                                                                            U.50 = quantile(bm.stock*0.001,probs=c(0.75),na.rm=T))
    
    for (x in 1:nrow(quants)){
      if(quants$stock[x] %in% LB.stocks) quants$fg[x] <- "LB"
      if(quants$stock[x] %in% LP.stocks) quants$fg[x] <- "LP"
      if(quants$stock[x] %in% M.stocks) quants$fg[x] <- "M"
    }
    
    #set up historical data
    p.lambdas <- lambdas
    colnames(p.lambdas)[colnames(p.lambdas) == "common"] <- "stock"
    colnames(p.lambdas)[colnames(p.lambdas) == "FG"] <- "fg"
    
    for (x in 1:nrow(p.lambdas)){
      if(p.lambdas$stock[x] %in% LB.stocks) p.lambdas$fg[x] <- "LB"
      if(p.lambdas$stock[x] %in% LP.stocks) p.lambdas$fg[x] <- "LP"
      if(p.lambdas$stock[x] %in% M.stocks) p.lambdas$fg[x] <- "M"
    }
    


    #stock names and fgs
    stock.info <- data.frame(stock  = c(LP.stocks, LB.stocks, M.stocks), fg = NA)
    for (i in 1:nrow(stock.info)) {
      if(stock.info$stock[i] %in% LB.stocks) stock.info$fg[i] <- "LB"
      if(stock.info$stock[i] %in% LP.stocks) stock.info$fg[i] <- "LP"
      if(stock.info$stock[i] %in% M.stocks) stock.info$fg[i] <- "M"
    }
    
    #get historical median in plotting format
    med.in <- NULL
    count <- 0
    for (f in c(LP.stocks, LB.stocks, M.stocks)){
      count <- count + 1
      med.in.tmp <- data.frame(year= 2000:2116, stock = f, fg = stock.info$fg[stock.info$stock == f], med.bm = NA)
      med.in.tmp$med.bm <- median(p.lambdas$full.bm[p.lambdas$stock == f])
      med.in[[count]] <- med.in.tmp
      med.in.tmp <- NULL
    }
    med.in.res <- do.call("rbind", med.in)
    
    #plot trajectory per stocks
    ggplot(quants[quants$scenario == "TD-Increased",]) +
      geom_line(aes(x = year, y = med/1000, colour = stock)) +
      scale_y_log10() +
      labs(x = "Year", y = "Biomass (thousands of tonnes)", colour = "Stock") +
      ggtitle("Top-down; 10% fishing on large piscivores") +
      theme(axis.text.x = element_text(size=14, angle=45),
            axis.text.y = element_text(size=14),
            legend.position = "right")
    
    
    #plot stock biomass (CIs)
    ggplot(quants[quants$scenario == "TD-Increased" & quants$med != 0,]) +
      geom_ribbon(aes(x = year, ymax = U.50/100, ymin = L.50/100, group = stock),
                  fill = "gold", colour = "gold", alpha = 0.4, size = 0.1) +
      geom_line(data = p.lambdas[p.lambdas$removals != 0,], aes(x = year, y = (removals*0.001)/100, group = stock), size = 0.8, colour = "black") +
      geom_line(data = med.in.res[med.in.res$med.bm != 0,], aes(x = year, y = (med.bm*0.001)/100, group = stock), linetype = 2, size = 0.8, colour = "black") +
      geom_line(aes(x = year, y = med/100, group = stock), colour = "black", size = 0.8, alpha = 0.7) +
      #geom_vline(xintercept = 2017, linetype = 2, size = 0.9, ) +
      facet_wrap(~interaction(stock, factor(fg, levels = c("LP", "LB", "M")) , sep = " FG: ",drop=T),scales='free_y') +
      labs(x = "Year", y = "Removals (hundreds of tonnes)") +
      theme(axis.text.x = element_text(size=14, angle=45),
            axis.text.y = element_text(size=14))
    #save as 1372, 872
    
    #get proper fgs into res.all
    res.lp <- res.all[res.all$stock %in% LP.stocks,]
    res.lp$fg <- "LP"
    res.lb <- res.all[res.all$stock %in% LB.stocks,]
    res.lb$fg <- "LB"
    res.m <- res.all[res.all$stock %in% M.stocks,]
    res.m$fg <- "M"
    res.fg <- rbind(res.lp, res.lb, res.m)
    
    res.rem.fg <- na.omit(res.fg)

    #plot raw sim stock biomass
    ggplot(res.fg[res.fg$scenario == "TD-Baseline" & res.fg$year != 2117,]) +
      #geom_line(data = p.lambdas, aes(x = year, y = (full.bm*0.001)/1000, group = stock, colour = sim), size = 0.8, colour = "black") +
      geom_line(aes(x = year, y = (stock.K*0.001)/1000, group = sim, colour = sim), size = 0.8) +
      geom_line(data = med.in.res[med.in.res$year > 2016,], aes(x = year, y = (med.bm*0.001)/1000, group = stock), linetype = 2, size = 0.6, colour = "black") +
      scale_colour_gradient(low = "grey", high = "orange") +
      facet_wrap(~interaction(stock, factor(fg, levels = c("LP", "LB", "M")) , sep = " FG: ",drop=T),scales='free_y') +
      labs(x = "Year", y = "Resource pool (thousands of tonnes)") +
      theme(axis.text.x = element_text(size=14, angle=45),
            axis.text.y = element_text(size=14),
            legend.position = "none")
    
    #get biomass versus resource pool for all stocks and plot
    rp.quants <- res.all |> collapse::fgroup_by(scenario, stock, year) |> collapse::fsummarize(L.50 = quantile(stock.K*0.001,probs=c(0.25),na.rm=T),
                                                                                            med = median(stock.K*0.001,na.rm=T),
                                                                                            U.50 = quantile(stock.K*0.001,probs=c(0.75),na.rm=T))
    
    for (x in 1:nrow(rp.quants)){
      if(rp.quants$stock[x] %in% LB.stocks) rp.quants$fg[x] <- "LB"
      if(rp.quants$stock[x] %in% LP.stocks) rp.quants$fg[x] <- "LP"
      if(rp.quants$stock[x] %in% M.stocks) rp.quants$fg[x] <- "M"
    }
    
    #join rp and realized median biomasses together
    colnames(rp.quants)[colnames(rp.quants) == "med"] <- "rp.med"
    colnames(rp.quants)[colnames(rp.quants) == "U.50"] <- "rp.U.50"
    colnames(rp.quants)[colnames(rp.quants) == "L.50"] <- "rp.L.50"
    rp.bm.plot <- merge(quants, rp.quants, by = c("scenario", "stock", "year", "fg"))
    
    rp.bm.plot$prop.med <- rp.bm.plot$med/rp.bm.plot$rp.med
    rp.bm.plot$prop.U.50 <- rp.bm.plot$U.50/rp.bm.plot$rp.U.50
    rp.bm.plot$prop.L.50 <- rp.bm.plot$L.50/rp.bm.plot$rp.L.50
    
    #set up historical proportion
    hist.rp.bm <- s.p.lam
    hist.rp.bm$year <- hist.rp.bm$year - 1
    
    #make plot (summary)
    ggplot(rp.bm.plot) +
      geom_hline(yintercept = 1, colour = "black", size = 0.8, linetype = 1, alpha = 0.7) +
      geom_line(aes(x = year, y = prop.med, colour = scenario), size = 0.8) +
      #geom_ribbon(aes(x = year, ymax = U.50, ymin = L.50, colour = scenario),
                  #fill = "darkblue", colour = "darkblue", alpha = 0.3, size = 0.1) +
      #geom_line(data = hist.rp.bm, aes(x = year, y = bm.prop, group = stock), size = 0.8) +
      scale_colour_manual(name = "Scenario",
                          breaks = c("BU-Baseline", "TD-Baseline",
                                     "BU-Increased", "TD-Increased"),
                          values = c("darkblue", "orange",
                                     "lightblue", "gold")) +
      labs(x = "Year", y = "Proportion of resource pool biomass",
           colour = "Type") +
      facet_wrap(~stock, scale = "free") +
      theme(axis.title = element_text(face = "bold", size = 20),
            legend.position = "top")
    
    #plot raw sim stock rp
    res.sim <- res.fg[res.fg$year %in% 2017:2116,]
    ggplot(res.sim[res.sim$scenario == "Bu-Baseline",]) +
      #geom_line(data = p.lambdas, aes(x = year, y = (full.bm*0.001)/1000, group = stock, colour = sim), size = 0.8, colour = "black") +
      geom_line(aes(x = year, y = (stock.K*0.001)/1000, group = sim), size = 0.8) +
      #geom_line(data = med.in.res, aes(x = year, y = (med.bm*0.001)/1000, group = stock), linetype = 2, size = 0.6, colour = "black") +
      scale_colour_gradient(low = "darkgrey", high = "darkblue") +
      facet_wrap(fg, scale - "free") +
      #facet_wrap(~interaction(stock, factor(fg, levels = c("LP", "LB", "M")) , sep = " FG: ",drop=T),scales='free') +
      labs(x = "Year", y = "Biomass (thousands of tonnes)") +
      theme(axis.text.x = element_text(size=14, angle=45),
            axis.text.y = element_text(size=14),
            legend.position = "none")
    
#stock biomass vs resource pools
    #get yearly stock RPs per sim
    srp.res <- res.all[,c(1:5,11,17)]

    #summary stats across sims (in tonnes)
    q.srp <- srp.res |> collapse::fgroup_by(scenario, year, stock) |> collapse::fsummarize(L.50 = quantile(stock.K*0.001,probs=c(0.25),na.rm=T),
                                                                                           med = median(stock.K*0.001,na.rm=T),
                                                                                           U.50 = quantile(stock.K*0.001,probs=c(0.75),na.rm=T))
    
    q.sbm <- srp.res |> collapse::fgroup_by(scenario, year, stock) |> collapse::fsummarize(L.50 = quantile(bm.stock*0.001,probs=c(0.25),na.rm=T),
                                                                                           med = median(bm.stock*0.001,na.rm=T),
                                                                                           U.50 = quantile(bm.stock*0.001,probs=c(0.75),na.rm=T))
    
    #combine resource pool and biomass into one data frame
    q.srp$bm <- "Resource pool"
    q.sbm$bm <- "Realized biomass"
    
    s.rp.bm.sim <- rbind(q.srp, q.sbm)
    
    q.srp.prop <- q.srp[,c(1:3, 5)]
    colnames(q.srp.prop)[colnames(q.srp.prop) == "med"] <- "med.rp"
    q.sbm.prop <- q.sbm[,c(1:3,5)]
    colnames(q.sbm.prop)[colnames(q.sbm.prop) == "med"] <- "med.bm"
    q.s.prop <- merge(q.srp.prop, q.sbm.prop, by = c("scenario", "year", "stock"))
    q.s.prop$prop.rp <- q.s.prop$med.bm/q.s.prop$med.rp
    
    l.prop.plot <- s.p.lam
    l.prop.plot$scenario <- "Historical"
    l.prop.plot$med.bm.prop <- 0
    for (i in c(M.stocks, LP.stocks, LB.stocks)){
      l.prop.plot$med.bm.prop[l.prop.plot$stock == i] <- mean(l.prop.plot$bm.prop[l.prop.plot$stock == i])
    }
    
    #getting year indexing right
    l.prop.plot$year <- l.prop.plot$year - 1
    
    ggplot(q.s.prop[]) +
      geom_hline(yintercept = 1, colour = "darkgrey") +
      geom_line(data = l.prop.plot[], aes(x = year, y = bm.prop, colour = scenario)) +
      geom_line(aes(x = year, y = prop.rp, colour = scenario)) +
      scale_x_continuous(breaks = c(2000, 2050, 2100)) +
      scale_colour_manual(name = "Scenario",
                          breaks = c("Historical", "BU-Baseline", "TD-Baseline",
                                     "BU-Increased", "TD-Increased"),
                          values = c("black", "darkblue", "orange", "lightblue", "gold")) +
      labs(x = "Year", y = "Proportion of resource pool biomass",
           colour = "Type") +
      facet_wrap(~stock, scale = "free") +
      theme(axis.title = element_text(face = "bold", size = 20),
            legend.position = "top")
    
    ggplot() +
      geom_hline(yintercept = 1, colour = "black") +
      geom_line(data = l.prop.plot[l.prop.plot$stock %in% c(LP.stocks, LB.stocks),], aes(x = year, y = med.bm.prop, colour = scenario)) +
      geom_line(data = summary.means())
      scale_colour_manual(name = "Scenario",
                          breaks = c("Historical", "BU-Baseline", "TD-Baseline", "BU-Increased", "TD-Increased"),
                          values = c("darkgrey", "darkblue", "orange", "lightblue", "gold")) +
        labs(x = "Year", y = "Proportion of resource pool biomass",
             colour = "Type") +
        facet_wrap(~stock, scale = "free") +
        theme(axis.title = element_text(face = "bold", size = 20),
              legend.position = "top")
      
#what will go into function
    #call it ts.med(requires data frames with year and median of historical data and each scenario)
    ggplot(com.bm) + geom_line(aes(x = year, y = bm.com/1000, colour = "Historical")) +
      geom_vline(xintercept = 2017, linetype = 2) +
      geom_line(data = med.com.bm.hist, aes(x = year, y = med/1000, colour = "Historical")) +
      geom_line(data = com.bnull, aes(x = year, y = med/1000, colour = "BU-Neutral")) +
      geom_ribbon(data = com.bnull, aes(x = year, ymax = U.50/1000, ymin = L.50/1000), alpha = 0.1, fill = "#0070B2") +
      geom_line(data = med.com.tnull, aes(x = year, y = med/1000, colour = "TD-Neutral")) +
      geom_ribbon(data = com.tnull, aes(x = year, ymax = U.50/1000, ymin = L.50/1000), alpha = 0.1, fill = "#F4D166") +
      labs(x = "Year", y = "Community biomass (millions of kg)") +
      scale_color_manual(
        name = 'Scenarios',
        breaks = c('Historical', 'BU-Neutral', 'TD-Neutral'),
        values = c('Historical' = "black", 
                   'BU-Neutral' = "#0070B2", 
                   'TD-Neutral' = "#F4D166")) +
      theme(
            legend.title = element_text(size = 20), legend.text = element_text(size = 14),
            axis.title = element_blank())
    
    td <- ggplot(com.bm) + geom_line(aes(x = year, y = bm.com/1000, colour = "Historical")) +
      geom_vline(xintercept = 2017, linetype = 2) +
      geom_line(data = med.com.bm.hist, aes(x = year, y = med/1000, colour = "Historical")) +
      geom_line(data = com.tnull, aes(x = year, y = med/1000, colour = "TD-Neutral")) +
      geom_ribbon(data = com.tnull, aes(x = year, ymax = U.50/1000, ymin = L.50/1000), alpha = 0.1, fill = "#F4D166") +
      #labs(x = "Year", y = "Community biomass (thousands of tonnes)") +
      scale_color_manual(
        name = 'Scenarios',
        breaks = c('Historical', 'BU-Neutral', 'TD-Neutral'),
        values = c('Historical' = "black",
                   'TD-Neutral' = "#F4D166")) +
      theme(
            legend.title = element_text(size = 20), legend.text = element_text(size = 14),
            axis.title = element_blank(),
            legend.position = "none")
    td
    
    bu <- ggplot(com.bm) + geom_line(aes(x = year, y = bm.com/1000, colour = "Historical")) +
      geom_vline(xintercept = 2017, linetype = 2) +
      geom_line(data = med.com.bm.hist, aes(x = year, y = med/1000, colour = "Historical")) +
      geom_line(data = com.bnull, aes(x = year, y = med/1000, colour = "BU-Neutral")) +
      geom_ribbon(data = com.bnull, aes(x = year, ymax = U.50/1000, ymin = L.50/1000), alpha = 0.1, fill = "#0070B2") +
      #labs(x = "Year", y = "Community biomass (millions of kg)") +
      scale_color_manual(
        name = 'Scenarios',
        breaks = c('Historical', 'BU-Neutral'),
        values = c('Historical' = "black", 
                   'BU-Neutral' = "#0070B2")) +
      theme(#axis.text.x = element_text(size=14, angle=45),
            #axis.text.y = element_text(size=14),
            legend.title = element_text(size = 20), legend.text = element_text(size = 14),
            axis.title = element_blank(),
            legend.position = "none")
    bu

    temp <- plot_grid(td,bu, ncol = 1)
    
    final <- ggdraw(temp) + draw_label("Year", x = 0.5, y = 0, vjust = 0, angle = 0, size = 25) +
      draw_label("Community biomass (thousands of tonnes)", x = 0, y = 0.6, vjust = 1, angle = 90, size = 25)
    final
    
    legend <- get_legend(
      final + theme(legend.box.margin = margin(0,0,0,12))
    )
    

##### Lambda comparison for null scenarios #####
#goal? show that model can recreate historical conditions
#isolate input lambdas
  lam.in <- na.omit(p.lambdas[,c("stock","lambda")])
  lam.in$scenario <- "Historical"
  lam.in$fg <- "S"
  for (i in 1:nrow(lam.in)){
    if(lam.in$stock[i] %in% LB.stocks) lam.in$fg[i] <- "LB"
    if(lam.in$stock[i] %in% LP.stocks) lam.in$fg[i] <- "LP"
    if (lam.in$stock[i] %in% M.stocks) lam.in$fg[i] <- "M"
  }
  
  #get simulated lambdas
  lam.sim <- na.omit(res.all[,c("stock","lambda", "scenario")])
  
  #get fg info for simulated lambdas (need to do it this way because df is so big)
  lam.sim.lp <- lam.sim[lam.sim$stock %in% LP.stocks,]
  lam.sim.lp$fg <- "LP"
  lam.sim.m <- lam.sim[lam.sim$stock %in% M.stocks,]
  lam.sim.m$fg <- "M"
  lam.sim.lb <- lam.sim[lam.sim$stock %in% LB.stocks,]
  lam.sim.lb$fg <- "LB"
  lam.sim <- rbind(lam.sim.lp, lam.sim.m, lam.sim.lb)
  
  lam.plot <- rbind(lam.in, lam.sim)
  
  #lam.plot$stock <- as.factor(lam.plot$stock)
  #levels(lam.plot$stock) <- c("Haddock","Atlantic wolffish","Thorny skate",
                              #"Atlantic cod","White hake","Pollock","American plaice","Monkfish",
                              #"Silver hake","Redfish","Witch flounder","Longfin hake",
                              #"Smooth skate","Longhorn sculpin","Sea raven")
  
  #make boxplot
    ggplot(lam.plot) + 
    geom_boxplot(aes(x = lambda, y = stock, fill = factor(scenario, levels = c("BU-Increased", "BU-Baseline",
                                                                               "Historical", "TD-Baseline", "TD-Increased"))), alpha = 0.5,
                 outlier.colour = "white", outlier.alpha = 1,
                 colour = "white", linewidth = 0.5) +
    scale_fill_manual(values = c("lightblue","darkblue", "black", "orange", "gold")) +
    geom_vline(xintercept = 1, linetype = 2, colour = "white", linewidth = 0.5) +
      scale_y_discrete(labels = scales::label_wrap(10)) +
    labs(y = "Stock", fill = "Scenario") + scale_x_log10(name="λ (log scale)") +
    facet_wrap(~factor(fg, levels = c("M", "LB", "LP")), scale = "free") +
    #facet_wrap(~fg, scale = "free") +
    theme(panel.background = element_rect("darkgrey"),
          legend.position = "top",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 17))
    
    stat_summary(
      fun = mean,
      fun.min = function(x) mean(x) - sd(x),
      fun.max = function(x) mean(x) + sd(x),
      geom = "pointrange"
    ) +
      
  #violin option
      ggplot(lam.plot) + 
      geom_violin(aes(x = lambda, y = stock, fill = factor(scenario, levels = c("BU-Increased", "BU-Baseline",
                                                                                 "Historical", "TD-Baseline", "TD-Increased"))), alpha = 0.5,
                   colour = "white", linewidth = 0.5, quantiles_linetype = c(0.5)) +
      scale_fill_manual(values = c("lightblue","darkblue", "black", "orange", "gold")) +
      geom_vline(xintercept = 1, linetype = 2, colour = "white", linewidth = 0.5) +
      scale_y_discrete(labels = scales::label_wrap(10)) +
      labs(y = "Stock", fill = "Scenario") + scale_x_log10(name="λ (log scale)") +
      facet_wrap(~factor(fg, levels = c("M", "LB", "LP")), scale = "free") +
      #facet_wrap(~fg, scale = "free") +
      theme(panel.background = element_rect("darkgrey"),
            legend.position = "top",
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 17))


##### Percent change comparisons #####
#goal? Compare how results change across scenarios at specific timepoints in the simulation
#stock-specific
#1. Biomass
  #first, get percent change info for each scenario
  #ID intervals over which you want to look at per change
  #intervals <- c(seq(20,100,20))
  intervals <- c(100)
  sim.types <- c(unique(res.all$scenario))
  delta.res <- NULL
  delta.stock <- NULL
  delta.all <- NULL
  
  for (i in sim.types){
    #isolate data to specific type
    sim.res <- res.all[res.all$scenario == i,]
    for (d in c(unique(res.all$stock))){
      #isolate type data to specific stock
      sim.stock <- sim.res[sim.res$stock == d,]
    for(x in 1:n.sims){
      #isolate stock data to specific sim
      s.sim <- sim.stock[sim.stock$sim == x,]
      #get starting biomass (from 2017 + warmup period)
      s.sim.start <- s.sim$bm.stock[s.sim$year == 2017 + 5]
      #initialize temporary data frame to hold results
      delta.tmp <- data.frame(stock = d, sim = x, intervals = intervals, bm.start = s.sim.start,
                              bm.int = NA, delta = NA, type = i)
    for(t in intervals){
      #get biomass in interval year
      if (t == 100)s.sim.end <- median(s.sim$bm.stock[s.sim$year %in% 2107:2117])
      if (t == 15) s.sim.end <- median(s.sim$bm.stock[s.sim$year %in% 2022:2032])
      #difference between starting biomass and interval year biomass
      s.sim.diff <- s.sim.start - s.sim.end
      #get percent change
      delta <- (s.sim.diff/s.sim.start)*100
      #get correct magnitude
      delta <- delta*-1
      #store results in temporary data frame
      delta.tmp$bm.int[delta.tmp$intervals == t] <- s.sim.end
      delta.tmp$delta[delta.tmp$intervals == t] <- delta
      }#end of interval loop
      #add temporary data frame to main results data frame
      delta.res[[x]] <- delta.tmp
      #reset delta.tmp
      delta.tmp <- NULL
      print(paste("done sim", x))
    }#end of sim loop
      delta.res <- do.call("rbind", delta.res)
      delta.stock[[d]] <- delta.res
      delta.res <- NULL
      #add temporary data frame to main results data frame
      print(paste("done ", d))
  }#end of stock loop
    delta.stock <- do.call("rbind", delta.stock)
    delta.all[[i]] <- delta.stock
    delta.stock <- NULL
    print(paste("done scenario", i))
}#end of percent change loop
  delta.all <- do.call("rbind", delta.all)
  row.names(delta.all) <- c(seq(1,nrow(delta.all)))
  
  #med.delta <- delta.all |> collapse::fgroup_by(type, stock, intervals) |> collapse::fsummarize(med = median(delta, na.rm = T))
  q.delta <- delta.all |> collapse::fgroup_by(type, stock, intervals) |> collapse::fsummarize(L.50 = quantile(delta,probs=c(0.25),na.rm=T),
                                                                                              med = median(delta,na.rm=T),
                                                                                              U.50 = quantile(delta,probs=c(0.75),na.rm=T))
  
  #add functional group labels to med.delta
  q.delta$fg <- "S"
  for (i in 1:nrow(q.delta)){
    if(q.delta$stock[i] %in% LB.stocks) q.delta$fg[i] <- "LB"
    if(q.delta$stock[i] %in% LP.stocks) q.delta$fg[i] <- "LP"
    if(q.delta$stock[i] %in% M.stocks) q.delta$fg[i] <- "M"
  }#end fg loop for q.delta
  
  #general summary of percent change across scenarios
  median(q.delta$med[q.delta$intervals == 100 & q.delta$type == "BU-Baseline"])
  median(q.delta$med[q.delta$intervals == 100 & q.delta$type == "BU-Increased"])
  median(q.delta$med[q.delta$intervals == 100 & q.delta$type == "TD-Baseline"])
  median(q.delta$med[q.delta$intervals == 100 & q.delta$type == "TD-Increased"])
  
  ggplot(q.delta[q.delta$intervals == 100 & q.delta$type == "TD-Increased",]) +
    geom_bar(aes(x = med, y = factor(stock, levels = c(M.stocks, LB.stocks, LP.stocks)),
                 fill = factor(fg, levels = c("M", "LB", "LP"))),
             stat = "identity") + 
    #scale_y_discrete(labels = scales::label_wrap(10)) +
    scale_fill_manual(values = c("darkblue", "orange", "gold")) +
    labs(x = "Median percent change", y = "Stock", fill = "Functional group") +
    ggtitle("Top-down; 10% exploitation on LP") +
    theme(legend.position = "right",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 17),
          plot.title = element_text(hjust = 0.5))
    
    
  ggplot(q.delta[q.delta$intervals == 100,], aes(x = med, y = stock, 
            fill = factor(type, levels = c("BU-Increased", "BU-Baseline",
                                           "TD-Baseline","TD-Increased")))) +
    geom_bar(position = "dodge", stat = "identity", width = 0.8, 
             linewidth = 0.5) +
    #geom_text(data = q.delta[q.delta$intervals == 100,],
              #aes(x = med,
                 # y = stock,
                  #label = paste(round(med, 0), "%")),
              #position = "stack")+
    scale_y_discrete(labels = scales::label_wrap(10)) +
    geom_vline(xintercept = 0, colour = "black", linetype = 2, size = 0.4) +
    scale_fill_manual(values = c("lightblue", "darkblue", "orange", "gold")) +
    labs(x = "Percent change", y = "Stock", fill = "Scenario") +
    facet_wrap(~factor(fg, levels = c("M", "LB", "LP")), scale = "free") +
    theme(legend.position = "top",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 17))
  
  #add functional group labels to delta.all
  delta.all$fg <- "S"
  for (i in 1:nrow(delta.all)){
    if(delta.all$stock[i] %in% LB.stocks) delta.all$fg[i] <- "LB"
    if(delta.all$stock[i] %in% LP.stocks) delta.all$fg[i] <- "LP"
    if(delta.all$stock[i] %in% M.stocks) delta.all$fg[i] <- "M"
  }#end fg loop for delta.all
  
  ggplot(delta.all[delta.all$type %in% c("TD-Baseline","TD-Increased"),]) +
    geom_violin(aes(x = stock, y = delta, fill = intervals), trim = T) +
    #facet_wrap(~factor(fg, levels = c("M", "LB", "LP")), scale = "free") +
    theme(legend.position = "top",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 17))
    
  
    #heatmaps
    #heatmap with percent change at regular intervals
      ggplot(med.delta[med.delta$type == "TD-Baseline",]) + 
        geom_tile(aes(x = intervals, y = stock, fill = med)) +
        geom_text(aes(x = intervals, y = stock, label = round(med, digits = 0)), colour = "lightgrey",
                  size = 3.5) +
        scale_fill_viridis_c(option = "magma") +
        scale_x_continuous(breaks = seq(20,100,by = 20)) +
        labs(x = "Projection year", y = "Stock", fill = "Percent\nchange") +
        facet_wrap(~fg, scales = "free", ncol = 1) +
        theme(strip.text.x = element_blank())
      
      LB.stocks <- c("Haddock", "Atlantic wolffish", "Thorny skate")
      stock.colours <- data.frame(stock = levels(med.delta$stock),
                                  colour = NA)
      for (i in 1:nrow(stock.colours)){
        if(stock.colours$stock[i] %in% LP.stocks) stock.colours$colour[i] <- "darkblue"
        if(stock.colours$stock[i] %in% LB.stocks) stock.colours$colour[i] <- "darkgrey"
        if(stock.colours$stock[i] %in% M.stocks) stock.colours$colour[i] <- "black"
      }
    
      ggplot(med.delta[med.delta$type == "BU-Baseline",]) + 
        geom_tile(aes(x = intervals, y = stock, fill = med)) +
        geom_rect(aes(xmin = 10, xmax = 110, ymin = 8.5, ymax = 15.5),
                  fill = "transparent", colour = "white", size = 1) +
        scale_fill_gradient2(low = "darkblue", mid = "white", high = "gold", midpoint = 0,
                             name = "Percent change\nin biomass") +
        #scale_fill_viridis_c(option = "magma") +
        geom_rect(aes(xmin = 10, xmax = 110, ymin = 3.5, ymax = 8.5),
                  fill = "transparent", colour = "white", size = 1) +
        geom_rect(aes(xmin = 10, xmax = 110, ymin = 0.5, ymax = 3.5),
                  fill = "transparent", colour = "white", size = 1) +
        scale_x_continuous(breaks = seq(20,100,by = 20)) +
        labs(x = "Projection year", y = "Stock", fill = "Percent\nchange") +
        #facet_wrap(~fg, scales = "free", ncol = 1)
        theme(axis.text.y = element_text(colour = stock.colours$colour),
              panel.background = element_rect("darkgrey"))
      
      
      
      ggplot(med.delta[med.delta$stock %in% L.stocks,]) + geom_tile(aes(x = intervals, y = stock, fill = med)) +
        scale_fill_gradient(low = "darkgrey",
                            high = "lightgrey") +
        geom_text(aes(x = intervals, y = stock, label = round(med, digits = 1)), color = "black", size = 3) +
        facet_wrap(~type, scales = "free") +
        guides(fill = guide_colourbar(title = "Percent \nchange")) +
        theme(legend.position = "none")
      
#2. Removals
      intervals <- 99
      sim.types <- c(unique(res.all$scenario))
      delta.res <- NULL
      delta.stock <- NULL
      delta.all <- NULL
      
      for (i in sim.types){
        #isolate data to specific type
        sim.res <- res.all[res.all$scenario == i,]
        for (d in c(unique(res.all$stock))){
          #isolate type data to specific stock
          sim.stock <- sim.res[sim.res$stock == d,]
          for(x in 1:n.sims){
            #isolate stock data to specific sim
            s.sim <- sim.stock[sim.stock$sim == x,]
            #get starting biomass (from 2017 + warmup period)
            s.sim.start <- s.sim$removals[s.sim$year == 2017 + 5]
            #initialize temporary data frame to hold results
            delta.tmp <- data.frame(stock = d, sim = x, intervals = intervals, bm.start = s.sim.start,
                                    bm.int = NA, delta = NA, type = i)
            for(t in intervals){
              #get biomass in interval year
              #if (t == 99)
                s.sim.end <- median(s.sim$removals[s.sim$year %in% 2106:2116])
             
              #difference between starting biomass and interval year biomass
              s.sim.diff <- s.sim.start - s.sim.end
              #get percent change
              delta <- (s.sim.diff/s.sim.start)*100
              #get correct magnitude
              delta <- delta*-1
              #store results in temporary data frame
              delta.tmp$bm.int[delta.tmp$intervals == t] <- s.sim.end
              delta.tmp$delta[delta.tmp$intervals == t] <- delta
            }#end of interval loop
            #add temporary data frame to main results data frame
            delta.res[[x]] <- delta.tmp
            #reset delta.tmp
            delta.tmp <- NULL
            print(paste("done sim", x))
          }#end of sim loop
          delta.res <- do.call("rbind", delta.res)
          delta.stock[[d]] <- delta.res
          delta.res <- NULL
          #add temporary data frame to main results data frame
          print(paste("done ", d))
        }#end of stock loop
        delta.stock <- do.call("rbind", delta.stock)
        delta.all[[i]] <- delta.stock
        delta.stock <- NULL
        print(paste("done scenario", i))
      }#end of percent change loop
      delta.all <- do.call("rbind", delta.all)
      row.names(delta.all) <- c(seq(1,nrow(delta.all)))
      
      #med.delta <- delta.all |> collapse::fgroup_by(type, stock, intervals) |> collapse::fsummarize(med = median(delta, na.rm = T))
      q.delta <- delta.all |> collapse::fgroup_by(type, stock, intervals) |> collapse::fsummarize(L.50 = quantile(delta,probs=c(0.25),na.rm=T),
                                                                                                  med = median(delta,na.rm=T),
                                                                                                  U.50 = quantile(delta,probs=c(0.75),na.rm=T))
      
      #add functional group labels to med.delta
      q.delta$fg <- "S"
      for (i in 1:nrow(q.delta)){
        if(q.delta$stock[i] %in% LB.stocks) q.delta$fg[i] <- "LB"
        if(q.delta$stock[i] %in% LP.stocks) q.delta$fg[i] <- "LP"
        if(q.delta$stock[i] %in% M.stocks) q.delta$fg[i] <- "M"
      }#end fg loop for q.delta
      
      #add in 0s for unfished stocks
      for (i in 1:nrow(q.delta)){
        if (is.na(q.delta$med[i])) q.delta$med[i] <- 0
      }

      #general summary of percent change across scenarios
      median(q.delta$med[q.delta$intervals == 99 & q.delta$type == "BU-Baseline"])
      median(q.delta$med[q.delta$intervals == 99 & q.delta$type == "BU-Increased"])
      hist(q.delta$med[q.delta$intervals == 99 & q.delta$type == "TD-Baseline"])
      #skewed
      hist(q.delta$med[q.delta$intervals == 99 & q.delta$type == "TD-Increased"])
      #skewed
      
      ggplot(q.delta[q.delta$intervals == 99 & q.delta$stock != "Thorny skate",], aes(x = med, y = stock,
              fill = factor(type, levels = c("BU-Increased", "BU-Baseline",
                                             "TD-Baseline","TD-Increased")))) +
        geom_bar(position = "dodge", stat = "identity", width = 0.8, 
                 linewidth = 0.5) +
        scale_y_discrete(labels = scales::label_wrap(10)) +
        geom_vline(xintercept = 0, colour = "black", linetype = 2, size = 0.4) +
        scale_fill_manual(values = c("lightblue", "darkblue", "orange", "gold")) +
        labs(x = "Percent change", y = "Stock", fill = "Scenario") +
        facet_wrap(~factor(fg, levels = c("M", "LB", "LP")), scale = "free") +
        theme(legend.position = "top",
              axis.text = element_text(size = 14),
              axis.title = element_text(size = 17))