#Making plots of multiple scenarios
library(ggplot2)
library(cowplot)

#Paper 1: comparing (1) BU and (2) TD with (a) no environmental constraint and (b) environmental constraint
  #how do these scenarios compare at the community, functional group and stock-specific scales?
    #plots for resource pools, biomass and lambda

    #fishing at historical intensities
    #worth bringing in economic value of the catch for commercial species?
      #would need price per unit weight
      #could use this site: https://www.dfo-mpo.gc.ca/stats/commercial/sea-maritimes-eng.htm#user-guide
        #divide total catch by total value
        #broken down into some species

  #could have all 4 combos on same plot (like DK)
  #could also have 4 panels on one figure
  #could also do a heatmap - % of simulations that fall below mean historical biomass per year
  #could also do violin plots - distribution of biomass across simulations at short-term and long-term check-in points

  #going to try all 4 of these:

#read in data
b.env <- res.sim.df 
b.null <- res.sim.df
t.env <- res.sim.df
t.null <- res.sim.df

b.fish <- res.sim.df
t.fish <- res.sim.df

####################################################################
M.stocks <- unique(b.null$stock[b.null$fg == "M"])
L.stocks <- unique(b.null$stock[b.null$fg == "L"])

#group nulls together
bnull.type <- b.null
bnull.type$type <- "BU-Neutral"
tnull.type <- t.null
tnull.type$type <- "TD-Neutral"
bnull.type <- bnull.type[,-13]
res.neut <- rbind(bnull.type, tnull.type)

#group envs together
benv.type <- b.env
benv.type$type <- "BU-Env"
tenv.type <- t.env
tenv.type$type <- "TD-Env"
res.env <- rbind(benv.type, tenv.type)

#group fishs together
bfish.type <- b.fish
bfish.type$type <- "BU-Fish"
tfish.type <- t.fish
tfish.type$type <- "TD-Fish"
res.fish <- rbind(bfish.type, tfish.type)

#store all scenarios together with labels
res.all <- rbind(res.neut, res.env, res.fish)

##### Time series like DK for null scenarios #####
  #goal? show that model can recreate historical conditions
  #community biomass example (would also make these for FGs and stocks)
    {
      #get inputs
    com.bm <- data.frame(year = 2000:2017, bm.com = eco.tot.bm.best$bm.eco)
    #get all data into tonnes
    com.bm$bm.com <- com.bm$bm.com*0.001
    com.bm$sim = 0
    med.com.bm <- median(com.bm$bm.com)
    med.com.bm.hist <- data.frame(year = 2000:2047, med = med.com.bm)
    #get medians and quantiles for each scenario
    com.bnull <- b.null |>  collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(bm.com*0.001,probs=c(0.25),na.rm=T),
                                                                              med = median(bm.com*0.001,na.rm=T),
                                                                              U.50 = quantile(bm.com*0.001,probs=c(0.75),na.rm=T))
    
    com.tnull <- t.null |>  collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(bm.com*0.001,probs=c(0.25),na.rm=T),
                                                                              med = median(bm.com*0.001,na.rm=T),
                                                                              U.50 = quantile(bm.com*0.001,probs=c(0.75),na.rm=T))
    #join into one data frame
    com.tnull$type <- "TD-Neutral"
    com.bnull$type <- "BU-Neutral"
    
    com.sim <- rbind(com.tnull, com.bnull)
    com.sim <- com.bnull
    
    ggplot(com.bm) + geom_line(aes(x = year, y = bm.com/1000)) +
      geom_vline(xintercept = 2017, linetype = 2) +
      geom_line(data = med.com.bm.hist, aes(x = year, y = med/1000)) +
      geom_line(data = com.sim[com.sim$type == "BU-Neutral",], aes(x = year, y = med/1000), , colour = "darkblue") +
      geom_ribbon(data = com.sim[com.sim$type == "BU-Neutral",], aes(x = year, ymax = U.50/1000, ymin = L.50/1000), 
                  alpha = 0.1, fill = "darkblue") +
      #geom_line(data = com.sim[com.sim$type == "TD-Neutral",], aes(x = year, y = med/1000), , colour = "darkblue") +
      #geom_ribbon(data = com.sim[com.sim$type == "TD-Neutral",], aes(x = year, ymax = U.50/1000, ymin = L.50/1000), 
       #           alpha = 0.1, fill = "darkblue") +
      labs(x = "Year", y = "Community biomass \n(thousands of tonnes)") +
      #facet_wrap(~type, ncol = 2, scale = "free") +
      theme(axis.title = element_text(face = "bold", size = 20))
    
    b.null <- res.sim.df
    bnull.type <- b.null
    bnull.type$type <- "BU-Neutral"
    
    com.bm <- data.frame(year = 2000:2017, bm.com = eco.tot.bm.best$bm.eco)
    #get all data into tonnes
    com.bm$bm.com <- com.bm$bm.com*0.001
    com.bm$sim = 0
    med.com.bm <- median(com.bm$bm.com)
    med.com.bm.hist <- data.frame(year = 2000:2047, med = med.com.bm)
    #get medians and quantiles for each scenario
    
    ggplot(com.bm) + geom_line(aes(x = year, y = bm.com/1000, group = sim)) +
      geom_vline(xintercept = 2017, linetype = 2) +
      geom_line(data = bnull.type, aes(x = year, y = (bm.com*0.001)/1000, group = sim), , colour = "darkblue") +
      geom_line(data = med.com.bm.hist, aes(x = year, y = med/1000), colour = "gold") +
      labs(x = "Year", y = "Community biomass \n(thousands of tonnes)") +
      theme(axis.title = element_text(face = "bold", size = 20))
    
    
    
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
    }

##### Lambda comparison for null scenarios #####
#goal? show that model can recreate historical conditions
{#isolate input lambdas
  lam.in <- p.lambdas[,c("lambda","stock")]
  lam.in$type <- "Historical"
  #join historical lambdas with simulated lambdas
  lam.sim <- na.omit(res.fish[,c("stock","lambda","type")])
  lam.plot <- rbind(lam.in, lam.sim)
  
  #make boxplot
  ggplot(lam.plot) + 
    geom_boxplot(aes(x = lambda, y = stock, fill = as.factor(type))) +
    scale_fill_manual("Scenario", values = c("lightblue", "darkgrey", "gold")) +
    geom_vline(xintercept = 1, linetype = 1, colour = "white") +
    labs(y = "Stock") + scale_x_log10(name="Lambda")
  
  ggplot(lam.plot[lam.plot$stock %in% L.stocks,]) + 
    geom_hline(yintercept = 1, colour = "grey") +
    geom_boxplot(aes(x = type, y = lambda, fill = as.factor(type))) +
    scale_fill_manual("Scenario", values = c("lightblue", "darkgrey", "gold")) +
    #geom_hline(yintercept = 1, linetype = 2) +
    labs(x = "Scenario") +
    scale_y_log10(name = "Lambda") +
    facet_wrap(~stock, ncol = 4) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          strip.background = element_rect(fill = "lightgrey", colour = "black"))
  
  #lambda plots as violins?
  ggplot(lam.plot[lam.plot$stock %in% L.stocks,]) + 
    geom_hline(yintercept = 1, colour = "grey") +
    geom_violin(aes(x = type, y = lambda, fill = as.factor(type))) +
    geom_boxplot(aes(x = type, y = lambda), width = 0.2, fill = "white", color = "black") +
    scale_fill_manual("Scenario", values = c("lightblue", "darkgrey", "gold")) +
    #geom_hline(yintercept = 1, linetype = 2) +
    labs(x = "Scenario") +
    scale_y_log10(name = "Lambda") +
    facet_wrap(~stock, ncol = 4) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          strip.background = element_rect(fill = "lightgrey", colour = "black"))}

##### Percent change comparisons #####
#goal? Compare how results change across scenarios at specific timepoints in the simulation
{
  #first, get percent change info for each scenario
  #ID intervals over which you want to look at per change
  intervals <- c(seq(5,30,5))
  sim.types <- c(unique(res.all$type))
  count <- 0
  delta.res <- NULL
  delta.stock <- NULL
  delta.all <- NULL
  
  for (i in sim.types){
    #isolate data to specific type
    sim.res <- res.all[res.all$type == i,]
    for (d in c(unique(res.all$stock))){
      count <- count + 1
      #isolate type data to specific stock
      sim.stock <- sim.res[sim.res$stock == d,]
    for(x in 1:n.sims){
      #isolate stock data to specific sim
      s.sim <- sim.stock[sim.stock$sim == x,]
      #get starting biomass (from 2017)
      s.sim.start <- s.sim$bm.stock[s.sim$year == 2017]
      #initialize temporary data frame to hold results
      delta.tmp <- data.frame(stock = d, sim = x, intervals = intervals, bm.start = s.sim.start,
                              bm.int = 0, delta = 0, type = i)
    for(t in intervals){
      #get biomass in interval year
      s.sim.end <- s.sim$bm.stock[s.sim$year == t + 2017]
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
      delta.res[[x]] <- (delta.tmp)
      #reset delta.tmp
      delta.tmp <- NULL
    }#end of sim loop
      delta.res <- do.call("rbind", delta.res)
      delta.stock[[count]] <- (delta.res)
      delta.res <- NULL
      #add temporary data frame to main results data frame
  }#end of stock loop
    delta.stock <- do.call("rbind", delta.stock)
    delta.all[[i]] <- (delta.stock)
    delta.stock <- NULL
}#end of percent change loop
  delta.all <- do.call("rbind", delta.all)
  row.names(delta.all) <- c(seq(1,nrow(delta.all)))
  
  med.delta <- delta.all |> collapse::fgroup_by(type, stock, intervals) |> collapse::fsummarize(med = median(delta, na.rm = T))
  df1_sorted <- df1 |>

    #heatmaps
    {#heatmap with percent change at regular intervals
      ggplot(med.delta[med.delta$stock %in% L.stocks,]) + geom_raster(aes(x = intervals, y = stock, fill = med)) +
        scale_fill_viridis_c(option = "magma") +
        facet_wrap(~type, scales = "free")
      
      ggplot(med.delta[med.delta$stock %in% L.stocks,]) + geom_tile(aes(x = intervals, y = stock, fill = med)) +
        scale_fill_gradient(low = "darkgrey",
                            high = "lightgrey") +
        geom_text(aes(x = intervals, y = stock, label = round(med, digits = 1)), color = "black", size = 3) +
        facet_wrap(~type, scales = "free") +
        guides(fill = guide_colourbar(title = "Percent \nchange")) +
        theme(legend.position = "none")
      
      
      #test this package with similarity analysis
      #restructure data
      dist.stock = NULL
      dist.res = NULL
      dist.delta <- NULL
      s.num <- 0
      i.num <- 0
      
      for (i in sim.types){
        i.num <- i.num + 1
        for (d in c(L.stocks, M.stocks)){
          s.num <- s.num + 1
          dist.tmp <- data.frame(stock = d, type = i, "5" = 0, "10" = 0, "15" = 0,
                                 "20" = 0, "25" = 0, "30" = 0)
          t.count <- 0
          for (t in intervals){
            t.count <- t.count + 1
            dist.tmp[2+t.count] <- med.delta$med[med.delta$type == i & med.delta$stock == d & med.delta$intervals == t]
          }#end of intervals loop
          dist.stock[[s.num]] <- (dist.tmp)
          dist.tmp <- NULL
        }#end of stock loop
        dist.stock <- do.call("rbind", dist.stock)
        dist.delta[[i.num]] <- dist.stock
        dist.stock <- NULL
      }#end of restructuring loop
      dist.delta <- do.call("rbind", dist.delta)
      
      #install.packages("pheatmap")
      #library(pheatmap)
      delta.map <- dist.delta[dist.delta$type == "BU-Neutral",]
      names <- delta.map[,1]
      delta.map <- delta.map[,-c(1,2)]
      row.names(delta.map) <- names
      colnames(delta.map) <- as.character(intervals)
      delta.map <- as.matrix(delta.map)
      cols <- colorRampPalette(c("lightblue", "white", "orange"))(100)
      bnull.hm <- pheatmap(delta.map,
                           display_numbers=TRUE, 
                           scale = "row",
                           color = cols,
                           number_colour = "black",
                           border_colour = "black",
                           cluster_cols = FALSE,
                           cellwidth = 30,
                           cellheight = 30,
                           main = "BU-Neutral")
      
      
      delta.map <- dist.delta[dist.delta$type == "BU-Env",]
      names <- delta.map[,1]
      delta.map <- delta.map[,-c(1,2)]
      row.names(delta.map) <- names
      colnames(delta.map) <- as.character(intervals)
      delta.map <- as.matrix(delta.map)
      benv.hm <- pheatmap(delta.map,
                          display_numbers=TRUE, 
                          scale = "row",
                          color = cols,
                          number_colour = "black",
                          border_colour = "black",
                          cluster_cols = FALSE,
                          cellwidth = 30,
                          cellheight = 30,
                          main = "BU-Env")}
  
    #violin plots
  {
    ggplot(delta.all[delta.all$stock %in% L.stocks & delta.all$intervals %in% c(5, 30),]) + 
      geom_violin(aes(x = rev(as.character(intervals)), y = delta, fill = as.factor(type))) +
      #scale_x_continuous(limits = c(5, 30), breaks = c(5)) +
      facet_wrap(~stock, ncol = 4)
    
    ggplot(delta.all[delta.all$type %in% c("BU-Neutral", "BU-Env") 
                     & delta.all$intervals %in% c(5, 30) &
                       delta.all$stock %in% L.stocks,]) + 
      geom_violin(aes(x = stock, y = delta, fill = as.factor(intervals))) +
      labs(fill = "Years \nelapsed") +
      #scale_x_continuous(limits = c(5, 30), breaks = c(5)) +
      facet_wrap(~type, ncol = 1) +
      theme(axis.text.x = element_text(size=14, angle=90),
            axis.text.y = element_text(size=14))
    
    
    ggplot(lam.plot[lam.plot$stock %in% L.stocks,]) + 
      geom_hline(yintercept = 1, colour = "grey") +
      geom_violin(aes(x = type, y = lambda, fill = as.factor(type))) +
      geom_boxplot(aes(x = type, y = lambda), width = 0.2, fill = "white", color = "black") +
      scale_fill_manual("Scenario", values = c("lightblue", "darkgrey", "gold")) +
      #geom_hline(yintercept = 1, linetype = 2) +
      labs(x = "Scenario") +
      scale_y_log10(name = "Lambda") +
      facet_wrap(~stock, ncol = 4) +
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            strip.background = element_rect(fill = "lightgrey", colour = "black"))
  }
  

  
}

##### Fishing ########
#My main three functions in Github

#/Scripts/Population_dynamics_function.R

#/Scripts/NS_catch_function.R
#/Scripts/trophic_model_function_bottom_up_multi_test.R

#Also will want to look at the RUNME as I have the management scenarios loaded in this way.  
#In theory you should be able to run everything if you loaded the Repo (just don't overwrite yours!!)
#trophic_level_RUNME.R
 
