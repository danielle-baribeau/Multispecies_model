# Here we develop a multi-species model for the Eastern Scotian Shelf.


# Notes

#1: stocks:     # historical estimates, realized lambdas, and functional groups for each stock 
#in your ecosystem, should be a list with each stocks being a names list, 
#with the name of the list being a unique stock name.
#2 lambdas:     # Lambda estimates from survey data
#3: n.yrs.proj  # How many years into the future we are going to project the stocks
#4: n.sims      # The numbers of simulations to run, keeping low for testing...
#5: repo.loc    # Location of the Github repo, defaults to "D:/GitHub/Multispecies_model/"
#6: mgmt        #list of management plan information for each stock


trophic.mod<-function(stocks = NULL,lambdas= NULL,n.yrs.proj = 50, n.sims = 20,
                      mgmt = list(mgmt = mgmt.scen,er.mn = NULL,er.sd = NULL),
                      repo.loc = "D:/GitHub/Multispecies_model",method = "not_sample")
{
  stock.eco <- unique(stocks$code)
  
  library(tidyverse)
  library(GGally)
  library(cowplot)
  library(ggthemes)
  library(boot)
  library(ggplot2)
  # Set the base plot theme
  theme_set(theme_few(base_size = 22))
  options(scipen = 999)
  # Download the function to go from inla to sf
  funs <- c("https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/simple_Lotka_r.r",
            "https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/simple_forward_sim.r",
            "https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/forward_project.r"
            
  )
  # Now run through a quick loop to load each one, just be sure that your working directory is read/write!
  for(fun in funs) 
  {
    download.file(fun,destfile = basename(fun))
    source(paste0(getwd(),"/",basename(fun)))
    file.remove(paste0(getwd(),"/",basename(fun)))
  }
  
  
  ########################## Section 2 Parameters ########################## Section 2 Parameters ########################## Section 2 Parameters
  
  # So here we are working to get the 'ecosystem' carrying capacity by looking at the total biomass for the ESS stocks we have
  # data for over the period of time we have data for all the stocks.
  # So here we pull out the data we need to look at total abundance and total biomass in the system by year...
  bm.tst <- stocks
  colnames(bm.tst)[colnames(bm.tst) == "code"] <- "Stock"
  colnames(bm.tst)[colnames(bm.tst) == "year"] <- "Year"
  colnames(bm.tst)[colnames(bm.tst) == "common"] <- "Species"
  colnames(bm.tst)[colnames(bm.tst) == "pop.bm"] <- "bm"
  
  # Look at the biomass and abundance in the ecosystem
  # (FIX, about 1% of the catch biomasses are larger than the actual biomass observed, take a look
  # and make sure that there isn't something mis-aligned for one of the stocks.)
  bm.tot <- bm.tst |> collapse::fgroup_by(Stock,Year,FG,Species) |> 
    collapse::fsummarize(bm = sum(bm,na.rm=T))
  
  # The 'ecosystem' biomass and numbers
  eco.tot.bm <- bm.tot |> collapse::fgroup_by(Year) |> 
    collapse::fsummarize(bm.eco = sum(bm))
  #functional group biomass and numbers.
  fg.bm <- bm.tot |> collapse::fgroup_by(Year,FG) |> 
    collapse::fsummarize(bm.fg = sum(bm))
  
  # All the bm together
  fg.eco.bm <- left_join(fg.bm, eco.tot.bm, by = "Year")
  fg.eco.bm$prop.eco.fg <- fg.eco.bm$bm.fg/fg.eco.bm$bm.eco
  #proportions add up to 1!
  
  # So now take the bm.tot and merge that with the total biomass and the functional group biomass so we can
  # look at what the stock does within its FG.
  
  # Now we combine the ecosystem results with the stock biomass's
  bm.final <- left_join(bm.tot,fg.eco.bm,by=c("Year","FG"))
  names(bm.final) <- c("Stock","Year","FG","species","bm.stock","bm.fg",'bm.eco',
                       'prop.eco.fg')
  # Get the proportion of the total biomass each stock accounts for
  bm.final <- bm.final |> collapse::fmutate(prop.bm.stock.eco = bm.stock/bm.eco,
                                            prop.bm.stock.fg = bm.stock/bm.fg)
  #test <- bm.final |> collapse::fsubset(Year == 1985) |> collapse::fsubset(FG == "MBZ")
  #sum(test$prop.bm.stock.fg)
  #checked 1985 - proportions add up to 1!
  #bm.final is 540 rows long
  
  # Remove 0s from the data
  #bm.final <- bm.final[bm.final$bm.stock > 0,]
  #no 0s in the ESS data
  bm.final <- as.data.frame(bm.final)
  
  #the years we want to look at from these data
  yrs <- seq(2000,2017)
  bm.best <- bm.final |> collapse::fsubset(Year %in% yrs)
  eco.tot.bm.best <- eco.tot.bm |> collapse::fsubset(Year %in% yrs)
  
  #getting biomass in more readable format (millions)
  plt.bm.best <- bm.best
  plt.bm.best$bm.stock <- plt.bm.best$bm.stock/1000000
  plt.bm.best$bm.fg <- plt.bm.best$bm.fg/1000000
  plt.bm.best$bm.eco <- plt.bm.best$bm.eco/1000000
  plt.eco.tot.bm <- eco.tot.bm.best
  plt.eco.tot.bm$bm.eco <- plt.eco.tot.bm$bm.eco/1000000
  # Biomass by functional group over time
  bm.fg.plt <- ggplot(plt.bm.best) + geom_line(aes(x=Year,y=bm.fg,group=FG,color=FG)) + 
    scale_color_manual(values = c("blue","black","darkgrey","lightblue")) + scale_y_log10(name="Biomass (millions) (log scale)") +
    theme(legend.position = 'bottom') + labs(colour = "Functional groups") +
    geom_line(data=plt.eco.tot.bm, aes(x=Year, y=bm.eco), color='black', linetype = 2)
  bm.fg.plt
  #drop in 2015 LBs driven by ~50% decline in wolffish and haddock (thorny skate also declined, but
  #not by as much)
  save_plot(paste0(repo.loc,"/FG_test_run/Biomass_by_fun_group.png"),bm.fg.plt,base_height = 8,base_width = 11)
  
  # This is real good now...
  prop.bm.fg.plt <- ggplot(bm.best) + geom_line(aes(x=Year,y=prop.eco.fg,group=FG,color=FG)) + 
    scale_color_manual(values = c("blue","black","darkgrey","lightblue")) + 
    scale_y_continuous(name="Proportion of ecosystem biomass") +   theme(legend.position = 'bottom') +
    labs(colour = "Functional groups")
  prop.bm.fg.plt
  #checking to make sure the fg/eco bm always add up to 1:
  for (i in 1:length(yrs)){
    test.yr <- bm.best |> collapse::fsubset(Year == yrs[i])
    test.yr <- test.yr[order(test.yr$FG), ]
    for (j in 1:nrow(test.yr)){
      if (identical(test.yr$FG[j], test.yr$FG[j+1])){
        test.yr$FG[j] <- NA
      }#end of if loop
    }#end of nrow for loop
    test.yr <- na.omit(test.yr)
    prop.eco.tot <- sum(test.yr$prop.eco.fg)
    print(prop.eco.tot)
  }#end of prop check for loop
  #all good!
  save_plot(paste0(repo.loc,"/FG_test_run/Prop_eco_biomass_by_fun_group.png"),prop.bm.fg.plt,base_height = 8,base_width = 11)
  
  # The biomass for the ecosystem
  bm.eco.plt <- ggplot(plt.eco.tot.bm) + geom_line(aes(x=Year,y=bm.eco)) + 
    scale_y_continuous(name="Biomass (millions)",limits = c(0,NA))
  bm.eco.plt
  save_plot(paste0(repo.loc,"/FG_test_run/Biomass_ess_ecosystem.png"),bm.eco.plt,base_height = 8,base_width = 11)
  
  #WIDER SCALE:
  #look at proportions at broad functional group level(large and medium; will function like trophic categories)
  bm.broad <- stocks
  colnames(bm.broad)[colnames(bm.broad) == "code"] <- "Stock"
  colnames(bm.broad)[colnames(bm.broad) == "year"] <- "Year"
  colnames(bm.broad)[colnames(bm.broad) == "common"] <- "Species"
  colnames(bm.broad)[colnames(bm.broad) == "pop.bm"] <- "bm"
  #add broad functional groups
  bm.broad$BFG <- NULL
  for (i in 1:nrow(bm.broad)){
    if (grepl('M', bm.broad$FG[i]) == TRUE){
      bm.broad$BFG[i] <- "M"
    } 
    if (grepl('L', bm.broad$FG[i]) == TRUE){
      bm.broad$BFG[i] <- "L"
    } 
  }
  
  bm.tot.b <- bm.broad |> collapse::fgroup_by(Stock,Year,BFG,Species) |> 
    collapse::fsummarize(bm = sum(bm,na.rm=T))
  
  #functional group biomass and numbers.
  bfg.bm <- bm.tot.b |> collapse::fgroup_by(Year,BFG) |> 
    collapse::fsummarize(bm.bfg = sum(bm))
  
  # All the bm together
  bfg.eco.bm <- left_join(bfg.bm, eco.tot.bm, by = "Year")
  bfg.eco.bm$prop.eco.bfg <- bfg.eco.bm$bm.bfg/bfg.eco.bm$bm.eco
  #proportions add up to 1!
  
  # Now we combine the ecosystem results with the stock biomass's
  bm.final.b <- left_join(bm.tot.b,bfg.eco.bm,by=c("Year","BFG"))
  names(bm.final.b) <- c("Stock","Year","BFG","species","bm.stock","bm.bfg",'bm.eco',
                         'prop.eco.bfg')
  # Get the proportion of the total biomass each stock accounts for
  bm.final.b <- bm.final.b |> collapse::fmutate(prop.bm.stock.eco = bm.stock/bm.eco,
                                                prop.bm.stock.bfg = bm.stock/bm.bfg)
  
  bm.final.b <- as.data.frame(bm.final.b)
  bm.best.b <- bm.final.b |> collapse::fsubset(Year %in% yrs)
  
  #getting biomass in more readable format (millions)
  plt.bm.best.b <- bm.best.b
  plt.bm.best.b$bm.stock <- plt.bm.best.b$bm.stock/1000000
  plt.bm.best.b$bm.bfg <- plt.bm.best.b$bm.bfg/1000000
  plt.bm.best.b$bm.eco <- plt.bm.best.b$bm.eco/1000000
  
  #repeat above plots for broad functional groups
  # Biomass by broad functional group over time
  ggplot(bm.best.b[bm.best.b$BFG == "L",]) + geom_area(aes(x = Year, y = bm.bfg/1000000))
  
  
  geom_line(aes(x = Year, y = bm.eco/1000000), linetype = 2)+
    geom_line(data = bm.best.b[bm.best.b$BFG == "M",], aes(x = Year, y = bm.bfg/1000000))+
    geom_ribbon(aes(ymin = min(bm.bfg/1000000), ymax = max(bm.bfg/1000000)),fill = "lightblue") + 
    geom_area(data = bm.best.b[bm.best.b$BFG == "L",], aes(x=Year, y = bm.bfg/1000000), fill = "darkblue") +
    labs(x = Year, y = "Biomass (millions kg)")
  
  
  bm.bfg.plt <- ggplot(plt.bm.best.b) + geom_line(aes(x=Year,y=bm.bfg,group=BFG,color=BFG), linewidth = 1.5) + 
    scale_color_manual(values = c("blue","lightblue")) + #scale_y_log10(name="Biomass (millions) (log scale)") +
    labs(colour = "Functional groups") +
    geom_line(aes(x = Year, y = bm.eco), linetype = 2, linewidth = 1.5) +
    labs(x = "Year", y = "Biomass (millions kg)") +
    theme_classic() +
    theme(legend.position = 'bottom', 
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20)) 
  bm.bfg.plt
  save_plot(paste0(repo.loc,"/FG_test_run/Describing ecosystem/BFG/Biomass_by_bfun_group.png"),bm.bfg.plt,base_height = 8,base_width = 11)
  
  prop.bm.bfg.plt <- ggplot(plt.bm.best.b) + geom_line(aes(x=Year,y=prop.eco.bfg,group=BFG,color=BFG)) + 
    scale_color_manual(values = c("blue","lightblue")) + 
    scale_y_continuous(name="Proportion of ecosystem biomass") +   theme(legend.position = 'bottom') +
    labs(colour = "Broad functional groups")
  prop.bm.bfg.plt
  save_plot(paste0(repo.loc,"/FG_test_run/Describing ecosystem/BFG/Prop_eco_biomass_by_bfun_group.2000.png"),prop.bm.bfg.plt,base_height = 8,base_width = 11)
  
  #MEDIUM SCALE: split between LPs, LBs and Mediums (to match option 1 for FG K simulation ideas)
  bm.med <- stocks
  colnames(bm.med)[colnames(bm.med) == "code"] <- "Stock"
  colnames(bm.med)[colnames(bm.med) == "year"] <- "Year"
  colnames(bm.med)[colnames(bm.med) == "common"] <- "Species"
  colnames(bm.med)[colnames(bm.med) == "pop.bm"] <- "bm"
  #add medium/broad functional groups
  bm.med$BFG <- NULL
  
  for (i in 1:nrow(bm.med)){
    if (bm.med$FG[i] == "MP" || bm.med$FG[i] == "MBZ"){
      bm.med$BFG[i] <- "M"
    } 
    if (bm.med$FG[i] == "LP"){
      bm.med$BFG[i] <- "LP"
    } 
    if (bm.med$FG[i] == "LB"){
      bm.med$BFG[i] <- "LB"
    } 
  }
  
  
  
  bm.tot.m <- bm.med |> collapse::fgroup_by(Stock,Year,BFG,Species) |> 
    collapse::fsummarize(bm = sum(bm,na.rm=T))
  
  #functional group biomass and numbers.
  mfg.bm <- bm.tot.m |> collapse::fgroup_by(Year,BFG) |> 
    collapse::fsummarize(bm.mfg = sum(bm))
  
  # All the bm together
  mfg.eco.bm <- left_join(mfg.bm, eco.tot.bm, by = "Year")
  mfg.eco.bm$prop.eco.mfg <- mfg.eco.bm$bm.mfg/mfg.eco.bm$bm.eco
  #proportions add up to 1!
  
  # Now we combine the ecosystem results with the stock biomass's
  bm.final.m <- left_join(bm.tot.m,mfg.eco.bm,by=c("Year","BFG"))
  names(bm.final.m) <- c("Stock","Year","MFG","species","bm.stock","bm.mfg",'bm.eco',
                         'prop.eco.mfg')
  # Get the proportion of the total biomass each stock accounts for
  bm.final.m <- bm.final.m |> collapse::fmutate(prop.bm.stock.eco = bm.stock/bm.eco,
                                                prop.bm.stock.mfg = bm.stock/bm.mfg)
  
  bm.final.m <- as.data.frame(bm.final.m)
  bm.best.m <- bm.final.m |> collapse::fsubset(Year %in% yrs)
  
  #getting biomass in more readable format (millions)
  plt.bm.best.m <- bm.best.m
  plt.bm.best.m$bm.stock <- plt.bm.best.m$bm.stock/1000000
  plt.bm.best.m$bm.mfg <- plt.bm.best.m$bm.mfg/1000000
  plt.bm.best.m$bm.eco <- plt.bm.best.m$bm.eco/1000000
  
  +
    geom_line(data=plt.eco.tot.bm, aes(x=Year, y=bm.eco), color='black', linetype = 2)
  
  
  
  #repeat above plots for medium functional groups
  # Biomass by medium functional group over time
  bm.mfg.plt <- ggplot(plt.bm.best.m) + geom_line(aes(x=Year,y=bm.mfg,group=MFG,color=MFG)) + 
    scale_color_manual(values = c("black", "blue","lightblue")) + scale_y_log10(name="Biomass (millions) (log scale)") +
    theme(legend.position = 'bottom') + labs(colour = "Medium functional groups")
  bm.mfg.plt
  save_plot(paste0(repo.loc,"/FG_test_run/Biomass_by_mfun_group.png"),bm.mfg.plt,base_height = 8,base_width = 11)
  #drop in LB after 2009 - which species is responsible?
  #seems to be haddock doing super well in 2009 that is causing this (wolf and skate get lower as well)
  #LB biomass in 2010 is close to that of 2008
  #wolffish and skate stay low into 2012; haddock returns to 2008 size
  
  #plot LP and M
  LP <- mfg.bm.best |> collapse::fsubset(BFG == "LP")
  M <- mfg.bm.best |> collapse::fsubset(BFG == "M")
  
  ggplot(LP) + geom_line(aes(x = Year, y = bm.mfg)) +
    geom_line(data = M, aes(x = Year, y = bm.mfg), colour = "blue")
  
  ccf(LP$bm.mfg, M$bm.mfg)
  
  M.detrend <- lm()
  
  prop.bm.mfg.plt <- ggplot(plt.bm.best.m) + geom_line(aes(x=Year,y=prop.eco.mfg,group=MFG,color=MFG)) + 
    scale_color_manual(values = c("black", "blue","lightblue")) + 
    scale_y_continuous(name="Proportion of ecosystem biomass") +   theme(legend.position = 'bottom') +
    labs(colour = "Medium functional groups")
  prop.bm.mfg.plt
  save_plot(paste0(repo.loc,"/FG_test_run/Prop_eco_biomass_by_mfun_group.png"),prop.bm.mfg.plt,base_height = 8,base_width = 11)
  
  
  
  # So now we want to look at stock level within a functional group
  # add some colors...
  bm.best$color <- "black"
  bm.best$color[bm.best$species %in% c("Haddock", "Silver hake", "Atlantic cod", "Redfish")] <- "blue"
  bm.best$color[bm.best$species %in% c("Atlantic wolffish", "Longhorn sculpin", "White hake",
                                       "Witch flounder")] <- "lightblue"
  bm.best$color[bm.best$species %in% c("Thorny skate", "Sea raven", "Pollock", "Longfin hake")] <- "black"
  bm.best$color[bm.best$species %in% c("American plaice", "Smooth skate")] <- "gold"
  bm.best$color[bm.best$species %in% c("Monkfish")] <- "grey"
  
  # Put in Species + functional group
  bm.best$spec.fg <- paste(bm.best$species,"(",bm.best$FG,")")
  # Pull out meta data
  meta.dat <- bm.best |> dplyr::group_by(Stock,FG,species,color,spec.fg) |> filter(row_number() >= (n() ))
  meta.dat <- meta.dat[,c("Stock","FG","species","color","spec.fg")]
  
  colors <- distinct(bm.best, spec.fg, color)
  pal <- colors$color
  names(pal) <- colors$spec.fg
  # Another color thing
  colors2 <- distinct(bm.best, species, color)
  pal2 <- colors$color
  names(pal2) <- colors$species
  
  
  stock.prop.bm.plt <- ggplot(bm.best) + geom_line(aes(x=Year,y=prop.bm.stock.fg,group = Stock,color=spec.fg),linewidth=2) + 
    facet_wrap(~FG) + guides(colour = guide_legend(nrow = 15)) + theme(legend.position = 'right') +
    scale_y_log10(name= "Proportion of FG biomass",n.breaks=10) +
    scale_color_manual(values=pal)
  stock.prop.bm.plt
  save_plot(paste0(repo.loc,"/FG_test_run/Prop_Biomass_ess_by_stock.png"),stock.prop.bm.plt,base_height = 8,base_width = 15)
  
  stock.bm.plt <- ggplot(bm.best) + geom_line(aes(x=Year,y=bm.stock,group = Stock,color=spec.fg),linewidth=2) + 
    facet_wrap(~FG) +
    scale_y_log10(name = "Biomass",n.breaks=7) + theme(legend.position = 'right') +
    guides(colour = guide_legend(nrow = 15)) + scale_color_manual(values=pal)
  stock.bm.plt
  save_plot(paste0(repo.loc,"/FG_test_run/Biomass_ess_by_stock.png"),stock.bm.plt,base_height = 8,base_width = 15)
  
  ######################################next, stock level within broad functional groups
  test <- bm.broad[,c("Species", "FG", "BFG")]
  bm.best.b$color <- "black"
  bm.best.b$color[bm.best.b$species %in% c("Silver hake", "Atlantic cod")] <- "blue"
  bm.best.b$color[bm.best.b$species %in% c("Haddock", "Redfish")] <- "lightblue"
  bm.best.b$color[bm.best.b$species %in% c("White hake", "Witch flounder")] <- "black"
  bm.best.b$color[bm.best.b$species %in% c("Pollock", "Longfin hake")] <- "gold"
  bm.best.b$color[bm.best.b$species %in% c("American plaice", "Longhorn sculpin")] <- "lightgrey"
  bm.best.b$color[bm.best.b$species %in% c("Atlantic wolffish", "Sea raven")] <- "darkgrey"
  bm.best.b$color[bm.best.b$species %in% c("Thorny skate", "Smooth skate")] <- "darkblue"
  bm.best.b$color[bm.best.b$species %in% c("Monkfish")] <- "darkorchid"
  
  # Put in Species + functional group
  bm.best.b$spec.bfg <- paste(bm.best.b$species,"(",bm.best.b$BFG,")")
  colors <- distinct(bm.best.b, spec.bfg, color)
  pal <- colors$color
  names(pal) <- colors$spec.bfg
  
  bstock.prop.bm.plt <- ggplot(bm.best.b) + geom_line(aes(x=Year,y=prop.bm.stock.bfg,group = Stock,color=spec.bfg),linewidth=2) + 
    facet_wrap(~BFG) + guides(colour = guide_legend(nrow = 15)) + theme(legend.position = 'right') +
    scale_y_log10(name= "Proportion of broad FG biomass",n.breaks=10) +
    scale_color_manual(values=pal)
  bstock.prop.bm.plt
  save_plot(paste0(repo.loc,"/FG_test_run/Prop_Biomass_ess_by_stock_bfg.png"),bstock.prop.bm.plt,base_height = 8,base_width = 15)
  
  bstock.bm.plt <- ggplot(bm.best.b) + geom_line(aes(x=Year,y=bm.stock,group = Stock,color=spec.bfg),linewidth=2) + 
    facet_wrap(~BFG) +
    scale_y_log10(name = "Biomass",n.breaks=7) + theme(legend.position = 'right') +
    guides(colour = guide_legend(nrow = 15)) + scale_color_manual(values=pal)
  bstock.bm.plt
  save_plot(paste0(repo.loc,"/FG_test_run/Biomass_ess_by_stock_bfg.png"),bstock.bm.plt,base_height = 8,base_width = 15)
  
  ###################################### Lastly, stock level within medium functional groups
  bm.best.m$color <- "black"
  bm.best.m$color[bm.best.m$species %in% c("Haddock", "Atlantic cod", "Silver hake")] <- "blue"
  bm.best.m$color[bm.best.m$species %in% c("Atlantic wolffish", 
                                           "White hake", "Redfish")] <- "lightblue"
  bm.best.m$color[bm.best.m$species %in% c("Thorny skate", "Pollock", "Redfish")] <- "black"
  bm.best.m$color[bm.best.m$species %in% c("American plaice", "Witch flounder")] <- "gold"
  bm.best.m$color[bm.best.m$species %in% c("Monkfish", "Longfin hake")] <- "lightgrey"
  bm.best.m$color[bm.best.m$species %in% c("Smooth skate")] <- "darkblue"
  bm.best.m$color[bm.best.m$species %in% c("Longhorn sculpin")] <- "darkgrey"
  bm.best.m$color[bm.best.m$species %in% c("Sea raven")] <- "darkorchid"
  
  # Put in Species + functional group
  bm.best.m$spec.mfg <- paste(bm.best.m$species,"(",bm.best.m$MFG,")")
  colors <- distinct(bm.best.m, spec.mfg, color)
  pal <- colors$color
  names(pal) <- colors$spec.mfg
  
  mstock.prop.bm.plt <- ggplot(bm.best.m) + geom_line(aes(x=Year,y=prop.bm.stock.mfg,group = Stock,color=spec.mfg),linewidth=2) + 
    facet_wrap(~MFG) + guides(colour = guide_legend(nrow = 5)) + theme(legend.position = 'bottom') +
    scale_y_log10(name= "Proportion of medium FG biomass",n.breaks=10) +
    scale_color_manual(values=pal)
  mstock.prop.bm.plt
  save_plot(paste0(repo.loc,"/FG_test_run/Prop_Biomass_ess_by_stock_mfg.png"),mstock.prop.bm.plt,base_height = 8,base_width = 15)
  
  mstock.bm.plt <- ggplot(bm.best.m) + geom_line(aes(x=Year,y=bm.stock,group = Stock,color=spec.mfg),linewidth=2) + 
    facet_wrap(~MFG) +
    scale_y_log10(name = "Biomass",n.breaks=7) + theme(legend.position = 'bottom') +
    guides(colour = guide_legend(nrow = 5)) + scale_color_manual(values=pal)
  mstock.bm.plt
  save_plot(paste0(repo.loc,"/FG_test_run/Biomass_ess_by_stock_mfg.png"),mstock.bm.plt,base_height = 8,base_width = 15)
  
  #overall, it all looks reasonable! At least from 2000 onwards, mediums are more plentiful than larges
  #trophic cascade? Matches literature well
  
  # So Model 1: You're Basic
  # OK, so within a FG each stock has it's own carrying capacity, that is nested within the FG carrying capacity
  # so if the FG is below the carrying capacity each stock gets a bit of that K space for the logistic model. 
  # The percentage of the K-space they get is contingent on their historic % of the carrying capacity the stock has had.
  # Go with the logistic model too, but I need to build in some uncertainty to the logistic projection
  # Initially I'm thinking I'll do (this comment will be outdated by the time you, dear reader, are reading this)
  # Step 1: We have a total K for the ecosystem based on past K's, let it vary
  # Step 2: We partition that to large piscivores and benthivores, based on historic splits
  # Step 4: Get mediums by looking at ratio of mediums to piscivores over time
  # Step 5: get medium piscivores and benthivories/zoopiscivores based on historic splits within mediums
  # Step 3: We then partition FG K to each stock, again based on historic proportion of the K, I wonder how 
  #         this will work if a stock is over-fished, the others will be able to fill some of the K-space, but 
  #         probably not all of it?
  # Step 4: Run the logistic model with the K each stock gets apportioned and we have a FG multispecies model
  
  # I think this should work, if we over-fish a stock everyone gets a bit of the free K space (including the overfished stock)
  # probably means the FG K isn't entirely filled  which could cause problems
  # So I can build the ecosystem to have a K that is based on the observed ecosystem biomass history and portion that out
  # to each of the FGs appropriately, BUT, the population won't necessarily reach that K in any given year, but 
  # I guess it should come close. So for base model we have the ecosystem biomass as our K, and 
  # then we see if the model is able to get the populations to achieve that K. If we fish
  # a bunch of stocks too hard, we have the K, but it'll never reach it. So, assumption that is could be
  # a bit problematic, we assume the past ecosystem biomass is K for these stocks, but with this logic, if 
  # we overfish we won't reach K, so we are assuming that these stocks were not overfished in totality and thus
  # the historic B trend is K (but in reality K was probably > B).... that is if you believe in K in any way shape or form.
  
  
  
  eco.tot.bm.best <- eco.tot.bm |> collapse::fsubset(Year %in% yrs)
  fg.bm.best <- fg.bm |> collapse::fsubset(Year %in% yrs)
  bfg.bm.best <- bfg.bm |> collapse::fsubset(Year %in% yrs)
  mfg.bm.best <- mfg.bm |> collapse::fsubset(Year %in% yrs)
  
  
  ##### NEW STUFF FOR FISHING ##########################################
  ## Finding reference points for each stock (low and high thresholds) ##
  #using 40% and 80% of median historic biomass
  
  #put this into a new function
  #allow for 0.4 and 0.8 to be changed
  #also make it so you don't have to pick the median
  #function inputs that allow user to input anything, but make defaults median/0.4/0.8
  
  #find median historic biomass values for all stocks
  #initialize results objects
  median.hist.bm <- c(rep(0, length(stock.eco)))
  count.med <- 1
  #find median
  for (s in stock.eco){
    #subset bm.best data to stock of interest
    bm.best.stock <- bm.best |> collapse::fsubset(Stock == s)
    #find median
    median.bm.stock <- median(bm.best.stock$bm.stock)
    #store in results vector
    median.hist.bm[count.med] <- median.bm.stock
    #update counter
    count.med <- count.med + 1
  }#end of median loop
  
  #find threshold values for all stocks
  #initialize results objects
  min.threshold <- c(rep(0, length(eco.stocks)))
  max.threshold <- c(rep(0, length(eco.stocks)))
  count.thresh <- 1
  #find thresholds
  for (s in stock.eco) {
    min.threshold[count.thresh] <- 0.4*median.hist.bm[count.thresh]
    max.threshold[count.thresh] <- 0.8*median.hist.bm[count.thresh]
    count.thresh <- count.thresh + 1
  }#end of threshold loop
  
  #add thresholds into management plan data frame
  mgmt.scen$low.threshold <- min.threshold
  mgmt.scen$high.threshold <- max.threshold
  
  
  ## Finding equation of the line between low and high thresholds for each stock ##
  #initialize objects
  bm.intercept <- c(rep(0, length(stock.eco)))
  bm.slope <- c(rep(0, length(stock.eco)))
  count.eqn <- 1
  
  #find eqn components***NEW CHANGES HERE
  #do they need to have different slopes?
  for (s in 1:length(stock.eco)) {
    #browser()
    #make a data frame that holds both thresholds for the stock
    #get halfway mark between low and high thresholds (for u of 0.2) - need more than 2 points to get a
    #different slope
    half.bm <- min.threshold[s] + ((max.threshold[s] - min.threshold[s])/2)
    thresholds <- data.frame(u = c(0, 0.2, 0.4), 
                             bm = c(min.threshold[s], half.bm, max.threshold[s]))
    #run model
    u.bm.mod <- lm(u ~ bm, thresholds)
    #find and extract coefficients
    u.bm.coef <- coef(u.bm.mod)
    u.bm.intercept <- u.bm.coef[1]
    u.bm.slope <- u.bm.coef[2]
    #store all values in results vectors
    bm.intercept[s] <- u.bm.intercept
    bm.slope[s] <- u.bm.slope
  }#end of equation loop
  
  #add equation components into management plan data frame
  mgmt.scen$intercept <- bm.intercept
  mgmt.scen$slope <- bm.slope
  
  #put biomass into readable format
  mgmt.scen$low.threshold <- mgmt.scen$low.threshold/1000000
  
  ###################################################################
  
  # So everything will need to get wrapped up in a simulation loop
  # Initialize some things, or maybe no things
  res.ts <- NULL
  ts.unpack <- NULL
  count.stock <- 1
  for(j in 1:n.sims)
  {
    st.time <- Sys.time()
    
    for(t in 1:n.yrs.proj)
    {
      # Get some starting points. These are for the current year
      base.eco.K.tmp <- sim.eco.K |> collapse::fsubset(sim == j & Years ==t)
      base.tl.K.tmp <- sim.troph.K |> collapse::fsubset(sim == j & Years ==t)
      base.stock.K.tmp <- sim.K.stocks |> collapse::fsubset(sim == j & Years ==t)
      # Now get the stock biomass from last year.
      if(t ==1)
      {
        stock.bm.last <- init.stock.bm
        stock.bm.last <- stock.bm.last[order(stock.bm.last$troph.cat),]
        eco.bm.last <- init.eco.bm
        tl.bm.last <- init.tl.bm
      }
      # Then we'll need to get these from the model simulations.
      if(t > 1)
      {
        # Use the handy av.wgt data.frame I made above
        bm.stocks <- data.frame(bm = NA,meta.dat)
        for(s in stock.eco) bm.stocks$bm[bm.stocks$Stock == s] <- res.ts[[s]]$net.bm[res.ts[[s]]$Years == t-1]
        bm.stocks$bm <- bm.stocks$bm
        stock.bm.last <- bm.stocks
        eco.bm.last <- sum(bm.stocks$bm)
        tl.bm.last <- bm.stocks |> collapse::fgroup_by(troph.cat) |> collapse::fsummarise(bm.tl = sum(bm))
      }
      
      # Now we need to figure out what K space is available for each stock within the trophic level.
      # First is our Trophic level above the K we have for it.
      # So this is the K space available in a given trophic level in a year
      base.tl.K.tmp$prop.K.space <- base.tl.K.tmp$bm.tl/tl.bm.last$bm.tl
      # We can then adjust the stock K's by the available K space in each stock
      
      base.stock.K.tmp$K.space <- NA
      base.stock.K.tmp$K.space[base.stock.K.tmp$troph.cat ==3] <- base.stock.K.tmp$bm.stock[base.stock.K.tmp$troph.cat ==3] * 
        (base.tl.K.tmp$prop.K.space[base.tl.K.tmp$troph.cat ==3]-1)
      base.stock.K.tmp$K.space[base.stock.K.tmp$troph.cat ==4] <- base.stock.K.tmp$bm.stock[base.stock.K.tmp$troph.cat ==4] * 
        (base.tl.K.tmp$prop.K.space[base.tl.K.tmp$troph.cat ==4]-1)
      base.stock.K.tmp$K.space[base.stock.K.tmp$troph.cat ==5] <- base.stock.K.tmp$bm.stock[base.stock.K.tmp$troph.cat ==5] * 
        (base.tl.K.tmp$prop.K.space[base.tl.K.tmp$troph.cat ==5]-1)
      
      base.stock.K.tmp$adj.K <- base.stock.K.tmp$bm.stock + base.stock.K.tmp$K.space
      
      # So now I have Carrying Capacities that take up (or lose) any available K space.
      # Now we can convert these to numbers using the historic 'average weight' of the stocks, to avoid complication
      # I'm just using the average of the average weight for each stock...
      #base.stock.K.tmp <- left_join(base.stock.K.tmp,av.wgt,by=c("Stock","troph.cat"))
      # And now we can get a K in numbers....
      #base.stock.K.tmp$adj.K.num <- base.stock.K.tmp$adj.K/base.stock.K.tmp$mn.wgt
      # Since I have Years and sim recorded, I should just be able to recursivly rbind this...
      if(t ==1 & j == 1) 
      {
        base.stock.K <- base.stock.K.tmp
        base.tl.K <- base.tl.K.tmp
        base.eco.K <- base.eco.K.tmp
      } else {
        base.stock.K <- rbind(base.stock.K,base.stock.K.tmp)
        base.tl.K <- rbind(base.tl.K,base.tl.K.tmp)
        base.eco.K <- rbind(base.eco.K,base.eco.K.tmp)
      } # end the else...
      for(s in stock.eco)
      {
        # Reset samples
        #browser()
        count.stock <- 1
        stock.lambdas <- lambdas[[count.stock]] 
        tmp.bm.last <- stock.bm.last |> collapse::fsubset(Stock == s)
        tmp.stock.K <- base.stock.K.tmp |> collapse::fsubset(Stock == s)
        bm.ts.stock <- bm.final[bm.final$Stock == s & bm.final$Year %in% first.year:last.year,]  
        #a=s
        #browser()
        # Now get the final year bm
        if(t == 1) 
        { 
          
          bm.start <- bm.ts.stock$bm.stock[bm.ts.stock$Year == last.year]
          res.ts[[s]] <- data.frame(net.bm = bm.start,removals = NA, mgmt.update = NA,
                                    update.type = NA, current.u = NA, next.yrs.u = NA,
                                    Stock = s,sim= j,lambda = NA,Years=t-1,
                                    troph.cat = as.numeric(unique(bm.tot$troph.cat[bm.tot$Stock ==s])),
                                    K.bm = NA)
          
        } else{ bm.start <- res.ts[[s]]$net.bm[res.ts[[s]]$Years == t-1]}
        
        # Sort out which of the years are low or high bm
        # I'm using 0.5 as the cut off, other options are valid (0.4 is my fav...)
        low.vs.high <- 0.5
        low.vs.high.bm <- low.vs.high * max(bm.ts.stock$bm.stock)
        # Have to drop the final year because we don't have a lambda estimate for the final year
        low.bm.years <- which(bm.ts.stock$bm.stock[-nrow(bm.ts.stock)] < low.vs.high.bm)
        if(length(low.bm.years) == 0) low.years <- F else low.years <- T
        high.bm.years <- which(bm.ts.stock$bm.stock[-nrow(bm.ts.stock)] >= low.vs.high.bm)
        # ANd get what our carrying capacity is at this moment
        cur.K <- tmp.stock.K$adj.K
        
        # So first, get a sample 
        method <- "not_sample"
        if(method == 'sample')
        {
          # Pick one of these to sample if that's how we want to roll, if we have low biomass years (as
          # defined by our cut off low.vs.high)
          if(low.years == T)
          {
            if(bm.start < low.vs.high.bm) samp <- sample(low.bm.years,1)
            if(bm.start >= low.vs.high.bm) samp <- sample(high.bm.years,1)
          } # end If we have low years
          if(low.years == F) samp <- sample(high.bm.years,1)
          # The simple way to do it is just to sample from the natural mortality distribution
          # Now using the right lambda, go look at trends from the stocks that are declining to see what's up there.
          lam.samp <- stock.lambdas$lambda[samp] # Get the sample years.  
        } # end the sample method.
        
        # Or do it the fun way...
        if(method != "sample")
        {
          
          # The fun way to do it is to do something multivariate! Note these are instantaneous now!!
          if(bm.start < low.vs.high.bm) 
          {
            if(length(low.bm.years) >0)
            {
              lam.mn <- mean(stock.lambdas$lambda[low.bm.years],na.rm=T)
              lam.sd <- sd(log(stock.lambdas$lambda[low.bm.years]),na.rm=T)
              if(length(low.bm.years) == 1) lam.sd <- 0.2 # In case there is just one low biomass year
            }
            if(length(low.bm.years) ==0)
            {
              lam.mn <- mean(stock.lambdas$lambda,na.rm=T)
              lam.sd <- sd(log(stock.lambdas$lambda),na.rm=T)
            }
            lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
            #if(is.na(lam.samp)) browser()
          } # end if(bm.start < low.vs.high.bm) 
          
          if(bm.start >= low.vs.high.bm & bm.start < cur.K) 
          {
            lam.mn <- mean(stock.lambdas$lambda[high.bm.years],na.rm=T)
            lam.sd <- sd(log(stock.lambdas$lambda[high.bm.years]),na.rm=T)
            lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
            
          } # end if(bm.start < low.vs.high.bm) 
          
        } # end if(method != "sample")
        #if(is.na(lam.samp)) browser()
        while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
        # Final one, if we are above the K, we are just doing the high biomass scenario for now.
        # Solution, sample from the lambdas at high biomass, but only take lambdas that are <= 1
        if(bm.start >= cur.K) 
        {
          lam.mn <- mean(stock.lambdas$lambda[high.bm.years],na.rm=T)
          lam.sd <- sd(log(stock.lambdas$lambda[high.bm.years]),na.rm=T)
          #going to switch to stretched beta dist
          #beta distribution bounded by minimum value in ts and maximum value in ts
          #lam.samp <- r_(1, median lambda, min = min observed, max = max observed)
          lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
          #if(is.na(lam.samp)) browser()
          while(lam.samp >1) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
        }
        #while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
        
        ################################## New catch function stuff #########################################
        #make a stock-specific management plan
        mgmt.stock <- mgmt.scen |> collapse::fsubset(stock == s)
        #isolate stock position number/assessment interval from mgmt plan
        stock.num <- mgmt.stock$stock.num
        assessment <- mgmt.stock$assessment.interval
        #run catch function
        catch.tmp <- proj.catch.eqn(mgmt.stock)
        
        #store results of catch function
        tst.res <- catch.tmp$tst.res
        removals <- catch.tmp$removals
        change <- catch.tmp$update
        change.type <- catch.tmp$update.type
        tst.res <- catch.tmp$tst.res
        ex.rate.curr <- catch.tmp$ex.curr
        next.yrs.u <- catch.tmp$ex.next
        
        # Because I'm only doing this one year at a time, there's something in here I need to mess around with to get the output tidy...
        res.ts[[s]] <- rbind(res.ts[[s]] ,data.frame(net.bm = tst.res,removals = removals, mgmt.update = change,
                                                     update.type = change.type, current.u = ex.rate.curr, next.yrs.u = next.yrs.u,
                                                     Stock = s,sim= j,lambda = lam.samp, Years=t,
                                                     troph.cat = as.numeric(unique(bm.tot$troph.cat[bm.tot$Stock ==s])),
                                                     K.bm = tmp.stock.K$adj.K))
        
        #update management plan to hold u for next year as new current u
        mgmt.scen$ex.curr[mgmt.scen$stock == s] <- next.yrs.u
        
        count.stock <- count.stock + 1
      }#end stock loop
      #if this an assessment year, update default exploitation rate in original management plan for the stock
      
      # Unpack the results
      ts.unpack[[j]] <- do.call('rbind',res.ts)
      
      #ggplot(ts.unpack[[j]]) + geom_line(aes(x= Years,y=abund,group=Stock,color=Stock)) + facet_wrap(~troph.cat) + scale_y_log10()
      
    }# end t loop
    # Pop a note when done each simulation
    timer <- Sys.time() - st.time
    print(paste("Simulation ", j))
    print(signif(timer,digits=2))
    
    #reset management plan for next iteration
    mgmt.scen$ex.curr <- c(rep(0.1, length(eco.stocks)))
    
  }#end n.sims
  
  # Unpack all the results.
  ts.final <- do.call("rbind",ts.unpack)
  #store ts.final in data frame form
  ts.final.df <- ts.final
  #order according to stock, year and simulation iteration
  ts.final.df <- ts.final.df[order(ts.final.df$Stock), ]
  #Group results by stock into lists
  ts.final.tmp <- NULL
  for(s in 1:length(stock.eco)) ts.final.tmp[[stock.eco[s]]] <- ts.final[ts.final$Stock==stock.eco[s],]
  ts.final <- ts.final.tmp
  
  ################################## Edits stop here ###############################################
  #DANIELLE'S OUTPUTS
  #THESE SHOULD GO IN SEPARATE FUNCTIONS
  
  #truncate stock IDs so that they can be plotted more easily
  #keeping species name and area; removing assessor
  #just doing this manually; specific to North Sea - couldn't figure out a good way to make it general
  short.names <- eco.stocks
  
  #STOCK STATUS
  #heatmap showing the percentage of simulation iterations where stock is at a status of interest
  #each assessment year
  #use numeric "update type' results (a.k.a "status code") to get stock status
  #1 = critical
  #2 = cautious
  #3 = healthy
  
  #choose status code you want to investigate:
  status.code <- 3
  
  #1. Make data frame that holds assessment year, stock and update type for all stocks
  #intialize intermediary objects
  name <- NULL
  stock.res <- NULL
  #initialize object to hold results
  stock.status.data <- NULL
  for (f in 1:length(eco.stocks)){
    #browser()
    #pull out stock-specific ts.final results
    stock.res <- ts.final.df |> collapse::fsubset(Stock == eco.stocks[f])
    #remove all years that aren't assessment years from stock-specific results
    stock.res <- stock.res |> collapse::fsubset(mgmt.update != 0)
    #condense stock.status to columns of interest
    #stock.res <- stock.res[ , c('mgmt.update', 'Years', 'sim')]
    #create data frame to hold proportion of simulations with stock status of interest
    name <- short.names[f]
    stock.status <- data.frame(stock = c(rep(name, length(unique(stock.res$Years)))),
                               assessment.year = unique(stock.res$Years), 
                               prop.status.code = c(rep(2, length(unique(stock.res$Years))))) 
    #prop.status.code initialized as 2 because a
    #proportion can't be higher than 1; for QC checks
    #Now, go through stock specific results and count the number of status codes of interest for each 
    #assessment year
    #initialize loop and results objects
    a.yr <- NULL
    yr.res <- NULL
    count.status.code <- NULL
    #loop through new data frame and fill in prop.status.code
    for (y in 1:length(stock.status$assessment.year)){
      #set up counter to track results
      count.status.code <- 0
      #store specific assessment year
      a.yr <- stock.status$assessment.year[y]
      #subset results to the year of interest
      yr.res <- stock.res |> collapse::fsubset(Years %in% a.yr)
      #count how many times status code of interest shows up across the iterations of the year of interest
      for (n in 1:length(yr.res$mgmt.update)){
        if (yr.res$mgmt.update[n] == status.code){
          count.status.code <- count.status.code + 1
        }
      }#end of counting stock status code loop
      #fill in prop.status.code:
      stock.status$prop.status.code[y] <- count.status.code/n.sims
      
    }#end of "fill stock.status" loop
    
    #store stock-specific results in main results data frame (long format for heatmap)
    
    #attach stock-specific results to full results data frame
    stock.status.data <- rbind(stock.status.data, stock.status)
  }#end of stock status heatmap loop
  #checked yrs 3 and 18 for sand lance and haddock - it worked!
  
  #2. Make heatmap
  #get min and max years from stock.status.data
  hm.min <- min(stock.status.data$assessment.year)
  hm.max <- max(stock.status.data$assessment.year)
  hm.interval <- stock.status.data$assessment.year[2] - hm.min
  #NOTE - for this to work, assessment interval needs to be the same for all stocks
  #plot data
  hm.stock.status <- ggplot(stock.status.data, aes(assessment.year, as.character(stock))) +                 
    geom_tile(aes(fill = prop.status.code)) +
    scale_fill_gradient(low = "#ADD8E6", high = "black") +
    labs(x = "Assessment year\n(Number of years from start of simulation)", y = "Stock ID",
         fill = "Proportion of\niterations") +
    theme(legend.position = "top",
          legend.key.width = unit(1.5, 'cm')) +
    scale_x_continuous(breaks = seq(hm.min, hm.max, by = hm.interval))
  hm.stock.status    
  
  #Set up input data for biomass, exploitation rates, lambda and K.bm in the same data frame
  #similar to stock status, except all projection years are included
  #individual stocks:
  #reset intermediary object
  name <- NULL
  stock.res <- NULL
  #initialize object to hold results
  res.data <- NULL
  f <- 1
  for (f in 1:length(eco.stocks)){
    # browser()
    #pull out stock-specific ts.final results
    stock.res <- ts.final.df |> collapse::fsubset(Stock == eco.stocks[f])
    #add stock ID to stock.res
    name <- eco.stocks[f]
    #add stockID to stock.res
    stock.res$stock.ID <- c(rep(name, length(stock.res$net.bm)))
    #add stock.res onto results df
    res.data <- rbind(res.data, stock.res)
  }#end of stock-specific results df setup
  
  
  quants <- ts.final.df |>  collapse::fgroup_by(Years,Stock,troph.cat) |> collapse::fsummarize(L.50 = quantile(net.bm,probs=c(0.25),na.rm=T),
                                                                                               med = median(net.bm,na.rm=T),
                                                                                               U.50 = quantile(net.bm,probs=c(0.75),na.rm=T))#,
  
  
  #BIOMASS
  #option 1 - show historical thresholds on this graph
  pdf(file = 'C:/Users/BARIBEAUD/Desktop/GitHub/ESS_Test_Run/Results/cn.no.outliers.dog.scul/2000/No_Fishing_Results/stock_bm.2000.pdf',
      height = 12, width = 12)
  count.bm <- 1
  name <- NULL
  bm.plot <- NULL
  for(f in 1:length(eco.stocks)){
    #isolate stock-specific historical data
    s.hist.bm <- bm.best |> collapse::fsubset(Stock == eco.stocks[f])
    #get minimum, median and maximum historical biomass
    hist.min.bm <- min(s.hist.bm$bm.stock)
    hist.med.bm <- median(s.hist.bm$bm.stock)
    hist.max.bm <- max(s.hist.bm$bm.stock)
    #get stock-specific data from results df
    name <- eco.stocks[f]
    stock.data <- quants |> collapse::fsubset(Stock == name)
    #make plot
    bm.plot[[as.character(eco.stocks[f])]] <- ggplot(stock.data, aes(Years, L.50)) +
      geom_line(linetype = 2) +
      geom_line(aes(x = Years, y = med)) +
      geom_line(aes(x = Years, y = U.50), linetype = 2) +
      scale_x_continuous(breaks = seq(0, n.yrs.proj, by = 2)) +
      coord_cartesian(xlim = c(0,n.yrs.proj)) +
      scale_y_log10(name = "Projected biomass") +
      labs(x = "Assessment year\n(Number of years from start of simulation)", title = paste(name)) +
      geom_hline(yintercept = hist.min.bm, linetype = 2, colour = "#5E8FB9", linewidth = 2) +
      annotate("text", x=-0.3, y=(hist.min.bm), label="1", fontface="italic",
               colour = "black", size = 7) +
      geom_hline(yintercept = hist.med.bm, linetype = 2, colour = "#5E8FB9", linewidth = 2) +
      annotate("text", x=-0.3, y=(hist.med.bm), label="2", fontface="italic",
               colour = "black", size = 7) +
      geom_hline(yintercept = hist.max.bm, linetype = 2, colour = "#5E8FB9", linewidth = 2) +
      annotate("text", x=-0.3, y=(hist.max.bm), label="3", fontface="italic",
               colour = "black", size = 7) +
      theme(plot.title = element_text(hjust = 0.5))
    print(bm.plot[[as.character(eco.stocks[f])]])
  }#end of biomass pdf plot loop
  dev.off()
  
  #option 2 - show stock status zones (alternative to heatmap?)
  pdf(file = 'C:/Users/BARIBEAUD/Desktop/GitHub/ESS_Test_Run/Results/cn.no.outliers.dog.scul/status_stock_bm.2000.pdf',
      height = 12, width = 12)
  count.bm <- 1
  name <- NULL
  for(f in 1:length(eco.stocks)){
    #get LRP and USR for this stock
    s.lrp <- mgmt.scen$low.threshold[f]
    s.usr <- mgmt.scen$high.threshold[f]
    #get stock-specific data from results df
    name <- eco.stocks[f]
    stock.data <- quants |> collapse::fsubset(Stock == name)
    #make plot
    bm.plot[[as.character(f)]] <- ggplot(stock.data, aes(Years, L.50)) +
      geom_line(linetype = 2) +
      geom_line(aes(x = Years, y = med)) +
      geom_line(aes(x = Years, y = U.50), linetype = 2) +
      scale_x_continuous(breaks = seq(0, n.yrs.proj, by = 2)) +
      scale_y_log10(name = "Projected biomass") +
      labs(x = "Assessment year\n(Number of years from start of simulation)", title = paste(name)) +
      annotate("rect", xmin = 0, xmax =n.yrs.proj, ymin = min(stock.data$L.50), ymax = s.lrp,
               alpha = .2,fill = "#D55E00") +
      geom_hline(yintercept = s.lrp, linetype = 2, colour = "darkgrey") +
      annotate("text", x=-1.5, y=s.lrp, label="LRP", fontface="italic", colour = "black", size = 6) +
      annotate("rect", xmin = 0, xmax = n.yrs.proj, ymin = s.lrp, ymax = s.usr,
               alpha = .2,fill = "#F0E442") +
      geom_hline(yintercept = s.usr, linetype = 2, colour = "darkgrey") +
      annotate("text", x=-1.5, y=s.usr, label="USR", fontface="italic", colour = "black", size = 6) +
      annotate("rect", xmin = 0, xmax = n.yrs.proj, ymin = s.usr, ymax = max(stock.data$U.50),
               alpha = .2,fill = "#009E73") +
      theme(plot.title = element_text(hjust = 0.5))
    print(bm.plot[[as.character(f)]])
  }#end of biomass pdf plot loop
  dev.off()
  
  #option 3 - violins of biomass at 5 years and 50 years for each stock
  #could also apply this to lambda or u
  yrs <- c(2, 10)
  violin.res <- ts.final.df |> collapse::fsubset(Years %in% yrs)
  violin.res$Years <- as.factor(violin.res$Years)
  
  #add in stock index numbers; too difficult to fit it all on x-axis of violin plot
  count.index <- 1
  violin.res$stock.index <- c(rep(0, length(violin.res$stock.ID)))
  f <- 1
  for (f in 1:length(eco.stocks)){
    name <- eco.stocks[f]
    for (j in 1:length(violin.res$stock.index)){
      if (violin.res$Stock[j] == name){
        violin.res$stock.index[j] <- count.index
      }#end of index if loop
    }#end of violin.res index loop
    #update count index
    count.index <- count.index + 1
  }#end of stock index loop
  
  #Biomass violin:
  colours <- list('2' = '#F0E442', '10' = '#E69F00')
  v.bm <- ggplot(violin.res, aes(x=as.character(Stock), y=net.bm, fill = Years)) +
    geom_violin() +
    scale_y_log10(name = "Projected biomass") +
    labs(x = "Stock", fill = "Elapsed\nprojection years") +
    scale_fill_manual(values = colours) +
    theme(legend.position = "top",
          legend.key.width = unit(0.5, 'cm'))
  
  v.bm
  
  
  
  
  #EXPLOITATION RATES (same options as for biomass)
  #show historical thresholds on this graph (can't do this yet because I don't have historical
  #fishing data set up)
  
  quants.u <- ts.final.df |>  collapse::fgroup_by(Years,Stock,troph.cat) |> collapse::fsummarize(L.50 = quantile(current.u,probs=c(0.25),na.rm=T),
                                                                                                 med = median(current.u,na.rm=T),
                                                                                                 U.50 = quantile(current.u,probs=c(0.75),na.rm=T))
  
  pdf(file = 'C:/Users/BARIBEAUD/Desktop/GitHub/ESS_Test_Run/Results/cn.no.outliers.dog.scul/stock_u.2000.pdf',
      height = 12, width = 12)
  count.bm <- 1
  name <- NULL
  for(f in 1:length(eco.stocks)){
    #get stock-specific data from results df
    name <- eco.stocks[f]
    stock.data <- quants.u |> collapse::fsubset(Stock == name)
    #remove year 0s from dataset
    stock.data <- stock.data[complete.cases(stock.data), ]
    #make plot
    bm.plot[[as.character(eco.stocks[f])]] <- ggplot(stock.data, aes(Years, L.50)) +
      geom_line(linetype = 2) +
      geom_line(aes(x = Years, y = med)) +
      geom_line(aes(x = Years, y = U.50), linetype = 2) +
      scale_x_continuous(breaks = seq(0, n.yrs.proj, by = 2)) +
      labs(x = "Assessment year\n(Number of years from start of simulation)", 
           y = "Projected exploitation rate (u)", title = paste(name)) +
      theme(plot.title = element_text(hjust = 0.5))
    print(bm.plot[[as.character(eco.stocks[f])]])
  }#end of biomass pdf plot loop
  dev.off()
  
  
  #LAMBDA
  quants.l <- ts.final.df |>  collapse::fgroup_by(Years,Stock,troph.cat) |> collapse::fsummarize(L.50 = quantile(lambda,probs=c(0.25),na.rm=T),
                                                                                                 med = median(lambda,na.rm=T),
                                                                                                 U.50 = quantile(lambda,probs=c(0.75),na.rm=T))
  #show historical thresholds on this graph
  pdf(file = 'C:/Users/BARIBEAUD/Desktop/GitHub/ESS_Test_Run/Results/cn.no.outliers.dog.scul/stock_lam.83.pdf',
      height = 12, width = 12)
  count.bm <- 1
  name <- NULL
  for(f in 1:length(eco.stocks)){
    #isolate stock-specific historical data
    s.hist.lam <-  stocks |> collapse::fsubset(code == eco.stocks[f])
    #get minimum, median and maximum historical biomass
    hist.min.lam <- min(s.hist.lam$lambda)
    hist.med.lam <- median(s.hist.lam$lambda)
    hist.max.lam <- max(s.hist.lam$lambda)
    #get stock-specific data from results df
    name <- eco.stocks[f]
    stock.data <- quants.l |> collapse::fsubset(Stock == name)
    #remove year 0 rows from stock.data
    stock.data <- stock.data[complete.cases(stock.data), ]
    #make plot
    bm.plot[[as.character(eco.stocks[f])]] <- ggplot(stock.data, aes(Years, L.50)) +
      geom_line(linetype = 2) +
      geom_line(aes(x = Years, y = med)) +
      geom_line(aes(x = Years, y = U.50), linetype = 2) +
      scale_x_continuous(breaks = seq(1, n.yrs.proj, by = 2)) +
      labs(x = "Assessment year\n(Number of years from start of simulation)", 
           y = "Lambda", title = paste(name)) +
      geom_hline(yintercept = hist.min.lam, linetype = 2, colour = "#E35408", linewidth = 2) +
      annotate("text", x=-0.3, y=(hist.min.lam), label="1", fontface="italic",
               colour = "black", size = 9) +
      geom_hline(yintercept = hist.med.lam, linetype = 2, colour = "#E35408", linewidth = 2) +
      annotate("text", x=-0.3, y=(hist.med.lam), label="2", fontface="italic",
               colour = "black", size = 9) +
      geom_hline(yintercept = hist.max.lam, linetype = 2, colour = "#E35408", linewidth = 2) +
      annotate("text", x=-0.3, y=(hist.max.lam), label="3", fontface="italic",
               colour = "black", size = 9) +
      theme(plot.title = element_text(hjust = 0.5))
    print(bm.plot[[as.character(eco.stocks[f])]])
  }#end of lambda pdf plot loop
  dev.off()
  
  
  #ECONOMIC VALUE OF THE CATCH
  #would be same as option 1 for biomass
  
  #Going to repeat with trophic category and total ecosystem
  
  ##################################################################################################
  
  #ggplot(ts.final) + geom_line(aes(x= Years,y=abund,group=sim,color=sim)) + facet_wrap(~Stock) + scale_y_log10()
  
  
  #ts.final$fm <- ts.final$removals/ts.final$abund
  #av.wgt$troph.cat <- as.numeric(av.wgt$troph.cat)
  #ts.final <- left_join(ts.final,av.wgt,by=c("Stock","troph.cat"))
  #ts.final$biomass <- ts.final$abund*ts.final$mn.wgt
  #r.final <- do.call("rbind",r.unpack)
  
  quants <- ts.final |>  collapse::fgroup_by(Years,Stock,troph.cat) |> collapse::fsummarize(L.50 = quantile(bm,probs=c(0.25),na.rm=T),
                                                                                            med = median(bm,na.rm=T),
                                                                                            U.50 = quantile(bm,probs=c(0.75),na.rm=T))#,
  #fml.50 = quantile(fm,probs=c(0.25),na.rm=T),
  #fm = median(fm,na.rm=T),
  #fmu.50 = quantile(fm,probs=c(0.75),na.rm=T))
  
  ts.final <- left_join(ts.final,meta.dat,by = c("Stock","troph.cat"))
  quants <- left_join(quants,meta.dat,by = c("Stock","troph.cat"))
  
  
  # If happy save the 2 objects
  saveRDS(object = ts.final,file = paste0(repo.loc,"/Results/NS_projections_",n.sims,"_sims_",min(years),"_to_",max(years),
                                          "_time_series_projections.Rds"))
  
  saveRDS(object = quants,file = paste0(repo.loc,"/Results/NS_projections_",n.sims,"_sims_",min(years),"_to_",max(years),
                                        "_time_series_quantiles.Rds"))
  # 
  # saveRDS(object = r.final,file = paste0(repo.loc,"/Results/NS_projections_",n.sims,"_sims_",min(years),"_to_",max(years),
  #                                        "_r_projections.Rds"))
  
  
  # Two simple plots. 
  p.sims <- ggplot(ts.final ) + geom_line(aes(x=Years,y=bm,group = sim,color=sim),alpha=0.8) +
    facet_wrap(~Stock,scales = 'free_y') + 
    scale_x_continuous(breaks = seq(1,50,by=49),labels=c(2015,2065)) +
    scale_y_log10(name = "Biomass") + 
    theme(legend.position = 'none') 
  
  save_plot(paste0(repo.loc,"/Figures/biomass_trends.png"),p.sims,base_height = 12,base_width = 24)
  
  
  
  p.sims.quants <- ggplot(quants) + geom_line(aes(x=Years,y=med,group=Stock,color=spec.tl)) + 
    facet_wrap(~troph.cat,scales = 'free_y') +  scale_y_log10(name="Biomass") +   theme(legend.position = 'top') +
    guides(colour = guide_legend(nrow = 7)) +
    scale_x_continuous(breaks = seq(1,50,by=49),labels=c(2015,2065)) 
  #geom_ribbon(data=quants, aes(x=Years,ymax=U.50,ymin = L.50),alpha=0.5,fill='blue',color='blue') 
  save_plot(paste0(repo.loc,"/Figures/Quantile_abundance_trends.png"),p.sims.quants,base_height = 8,base_width = 16)
  
  
  colors <- distinct(bm.best, spec.tl, color)
  pal <- colors$color
  names(pal) <- colors$spec.tl
  
  p.sims.quants <- ggplot(quants) + geom_line(aes(x=Years,y=med,group=Stock,color=spec.tl),linewidth=2) + 
    facet_wrap(~troph.cat,scales = 'free_y') +  scale_y_log10(name="Biomass") +   theme(legend.position = 'top') +
    guides(colour = guide_legend(nrow = 5)) + scale_color_manual(values=pal) +
    scale_x_continuous(name="",breaks = seq(1,50,by=49),labels=c(2015,2065)) 
  #geom_ribbon(data=quants, aes(x=Years,ymax=U.50,ymin = L.50),alpha=0.5,fill='blue',color='blue') 
  save_plot(paste0(repo.loc,"/Figures/Quantile_biomass_trends.png"),p.sims.quants,base_height = 8,base_width = 16)
  
  
  
  # ggplot(bm.best) + geom_line(aes(x=Year,y=prop.bm.stock.tl,group = Stock,color=spec.tl),linewidth=2) + 
  #   facet_wrap(~troph.cat) + guides(colour = guide_legend(nrow = 5)) + theme(legend.position = 'top') +
  #   scale_y_log10(name= "Proportion of biomass",n.breaks=10) + scale_x_continuous(name="",labels = c(1990,2000,2010),breaks=c(1990,2000,2010))+
  #   scale_color_manual(values=pal)
  
  return(list(sim.quantiles = quants,
              sim.ts = ts.final,
              past.bm = bm.best,
              sim.K.stocks = sim.K.stocks,
              sim.troph.K = sim.troph.K,
              sim.eco.K = sim.eco.K))
} # end function