

dat.loc <- 'C:/Users/BARIBEAUD/Desktop/GitHub/ESS-Model-Setup/May 2026 Removals'
repo.loc <- 'C:/Users/BARIBEAUD/Desktop/GitHub/Multispecies_model/ESS_Test_Run'
load(file = paste0(dat.loc,"/removals_biomass_inputs.RData"))
# load in the FG model function
source(paste0(repo.loc,"/E_FG_model_function_testing_mgmt.R"))
#rename inputs
lam <- lam.ed.ts

# Get the right stocks, this is from our model runs
eco.stocks <- unique(lam$code)

# Get the functional groups
# FGs
eco.fg <- data.frame(common = c("Atlantic cod", "Haddock", "White hake", "Silver hake", "Pollock", 
                                 "Redfish","American plaice", "Witch flounder", "Atlantic wolffish", 
                                 "Longfin hake","Thorny skate", "Smooth skate", 
                                 "Longhorn sculpin", "Sea raven","Monkfish"),
                        code = eco.stocks,
                        FG = c("LP", "LB", "LP", "MP", "LP", "MBZ", "LP", "MBZ", "LB", "MBZ", "LB",
                               "MBZ", "MP", "MP", "LP"))
                        #LP = large piscivore
                        #LB = large benthivore
                        #MP = medium piscivore
                        #MBZ = medium benthivore/zoopiscivore


#add FGs to biomass inputs
stocks <- merge(lam, eco.fg, by = "code")

# Make it a list
stocks.lst <- NULL
for(s in 1:length(eco.stocks)) stocks.lst[[s]] <- stocks |> collapse::fsubset(code == eco.stocks[s])

#eco.lambdas <- NULL
#for(s in 1:length(eco.stocks)) eco.lambdas[[s]] <- lam |> collapse::fsubset(code == eco.stocks[s])



# Here we develop a multi-species model for the Eastern Scotian Shelf.

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
  #colnames(bm.tst)[colnames(bm.tst) == "pop.bm"] <- "bm"
  colnames(bm.tst)[colnames(bm.tst) == "full.bm"] <- "bm"
  
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

  #WIDER SCALE:
  #look at proportions at broad functional group level(large and medium; will function like trophic categories)
  bm.broad <- stocks
  colnames(bm.broad)[colnames(bm.broad) == "code"] <- "Stock"
  colnames(bm.broad)[colnames(bm.broad) == "year"] <- "Year"
  colnames(bm.broad)[colnames(bm.broad) == "common"] <- "Species"
  #colnames(bm.broad)[colnames(bm.broad) == "pop.bm"] <- "bm"
  colnames(bm.tst)[colnames(bm.tst) == "full.bm"] <- "bm"
  
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
    collapse::fsummarize(bm = sum(full.bm,na.rm=T))
  
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

  #MEDIUM SCALE: split between LPs, LBs and Mediums (to match option 1 for FG K simulation ideas)
  bm.med <- stocks
  colnames(bm.med)[colnames(bm.med) == "code"] <- "Stock"
  colnames(bm.med)[colnames(bm.med) == "year"] <- "Year"
  colnames(bm.med)[colnames(bm.med) == "common"] <- "Species"
  #colnames(bm.med)[colnames(bm.med) == "pop.bm"] <- "bm"
  colnames(bm.tst)[colnames(bm.tst) == "full.bm"] <- "bm"
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
    collapse::fsummarize(bm = sum(full.bm,na.rm=T))
  
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

  
  eco.tot.bm.best <- eco.tot.bm |> collapse::fsubset(Year %in% yrs)
  fg.bm.best <- fg.bm |> collapse::fsubset(Year %in% yrs)
  bfg.bm.best <- bfg.bm |> collapse::fsubset(Year %in% yrs)
  mfg.bm.best <- mfg.bm |> collapse::fsubset(Year %in% yrs)
  
###### Chapter 1 plots (872, 572) #####
#community biomass
  ggplot(eco.tot.bm.best) + geom_line(aes(x = Year, y = bm.eco/1000000)) +
    geom_hline(yintercept = median(eco.tot.bm.best$bm.eco/1000000), linetype = 2) +
    labs(x = "Year", y = "Community biomass (million kg)")
  
  
#functional group biomass
  ggplot(mfg.bm.best) + 
    geom_line(data = eco.tot.bm.best, aes(x = Year, y = bm.eco/1000000), linetype = 2) +
    geom_line(aes(x = Year, y = bm.mfg/1000000, group = BFG, colour = BFG)) +
    scale_colour_manual(values = c("#C26B51", "#E5A94E", "#8EB3D2" )) +
    labs(x = "Year", y = "Biomass (million kg)") +
    labs(colour = "Functional\ngroup")
  
  #functional group biomass 2
  ggplot(bfg.bm.best) + 
    geom_line(data = eco.tot.bm.best, aes(x = Year, y = bm.eco/1000000), linetype = 3) +
    geom_line(aes(x = Year, y = bm.bfg/1000000, group = BFG, colour = BFG)) +
    scale_colour_manual(values = c("#26456E", "#8EB3D2" )) +
    geom_hline(yintercept = median(bfg.bm.best$bm.bfg[bfg.bm.best$BFG == "L"]/1000000), colour = "#26456E", linetype = 2) +
    geom_hline(yintercept = median(bfg.bm.best$bm.bfg[bfg.bm.best$BFG == "M"]/1000000), colour = "#8EB3D2", linetype = 2) +
    labs(x = "Year", y = "Biomass (million kg)") +
    labs(colour = "Functional\ngroup")
  
  #proportion of stocks in functional groups
  l.prop.stat <- bm.best.b[bm.best.b$BFG == "L",] |> group_by(species) |> summarize(med.prop <- mean(prop.bm.stock.bfg))
  #mscale <- scales::seq_gradient_pal("#540046", "#D0E5EE", "Lab")(seq(0,1,length.out = 7))
  m.prop <- ggplot(bm.best.b[bm.best.b$BFG == "M",]) +
    geom_line(aes(x = Year, y = prop.bm.stock.bfg, group = species, colour = species)) +
    #scale_y_break(c(0.0225, 0.23)) +
    #scale_y_log10() +
    scale_colour_manual(values = c("#0054FF","#3299FF","#FF5500","#FF9932","#65CCFF",
                                   "#99DEFF","#FFEE99", "#CCFFFF")) +
    #scale_colour_viridis() +
    labs(x = "Year", y = "Proportion of medium biomass") +
    labs(colour = "Stock") +
    theme(panel.background = element_rect(fill = 'darkgrey')) +
    ggtitle("Regular y-axis")
  m.prop
  
  m.prop.log <- ggplot(bm.best.b[bm.best.b$BFG == "M",]) +
    geom_line(aes(x = Year, y = prop.bm.stock.bfg, group = species, colour = species)) +
    #scale_y_break(c(0.0225, 0.23)) +
    scale_y_log10() +
    scale_colour_manual(values = c("#0054FF","#3299FF","#FF5500","#FF9932","#65CCFF",
                                   "#99DEFF","#FFEE99", "#CCFFFF")) +
    #scale_colour_viridis() +
    labs(x = "Year", y = "Proportion of medium biomass (log scale)") +
    labs(colour = "Stock") +
    theme(panel.background = element_rect(fill = 'darkgrey')) +
    ggtitle("Adjusted y-axis")
  m.prop.log
  
  ggarrange(m.prop, m.prop.log, ncol = 1, common.legend = TRUE)
  
  l.prop <- ggplot(bm.best.b[bm.best.b$BFG == "L",]) +
    geom_line(aes(x = Year, y = prop.bm.stock.bfg, group = species, colour = species)) +
    #scale_y_break(c(0.0225, 0.23)) +
    #scale_y_log10() +
    scale_colour_manual(values = c("#0054FF","#3299FF","#FF5500","#FF9932","#65CCFF",
                                   "#99DEFF","#FFEE99", "#CCFFFF")) +
    #scale_colour_viridis_d(option = "plasma") +
    labs(x = "Year", y = "Proportion of large biomass") +
    labs(colour = "Stock") +
    theme(panel.background = element_rect(fill = 'darkgrey'),
          legend.position = "none") +
  ggtitle("Larges")
  l.prop
  
  lp.prop <- ggplot(bm.best.m[bm.best.m$MFG == "LP",]) +
    geom_line(aes(x = Year, y = prop.bm.stock.mfg, group = species, colour = species)) +
    #scale_colour_viridis_d(option = "plasma") +
    scale_colour_manual(values = c("#0054FF", "#3299FF","#65CCFF", "#99DEFF", "#CCFFFF"))+
    labs(x = "Year", y = "Proportion of large biomass", colour = "Stock") +
    theme(panel.background = element_rect(fill = 'darkgrey'),
          legend.position = "none") +
    ggtitle("Large piscivores")
  lp.prop
  
  lb.prop <- ggplot(bm.best.m[bm.best.m$MFG == "LB",]) +
    geom_line(aes(x = Year, y = prop.bm.stock.mfg, group = species, colour = species)) +
    #scale_colour_viridis_d(option = "plasma") +
    scale_colour_manual(values = c("#FF5500", "#FF9932","#FFEE99"))+
    labs(x = "Year", y = "Proportion of large biomass", colour = "Stock") +
    theme(panel.background = element_rect(fill = 'darkgrey'),
          legend.position = "none") +
    ggtitle("Large benthivores")
  lb.prop
  
  small.prop <- ggarrange(lp.prop, lb.prop, ncol = 1)
  small.prop
  
  ggarrange(l.prop, small.prop, common.legend = TRUE)

  
  #ggarrange(l.prop, lp.prop, lb.prop, ncol = 2, common.legend = TRUE)
  
  
  ggplot(bm.best.b[bm.best.b$BFG == "M",], aes(prop.bm.stock.bfg, Year)) +
    geom_col(orientation = "prop.bm.stock.bfg") +
    scale_x_break(c())
  
  #stocks
  l <- ggplot(bm.best.b[bm.best.b$BFG == "L",]) +
    geom_line(aes(x = Year, y = bm.stock/1000000)) +
    labs(x = "Year", y = "Stock biomass (million kg)") +
  facet_wrap(~species, scale = "free")
  
  m <- ggplot(bm.best.b[bm.best.b$BFG == "M",]) +
    geom_line(aes(x = Year, y = bm.stock/1000000)) +
    labs(x = "Year", y = "Stock biomass (million kg)") +
    facet_wrap(~species, scale = "free", ncol = 3)
  m
  
  median.bm <- NULL
  for (i in unique(bm.best.m$species)){
    dat <- bm.best.m[bm.best.m$species == i,]
    colnames(dat)[colnames(dat) == "BFG"] <- "MFG"
    dat$med.bm <- median(dat$bm.stock)
    median.bm[[i]] <- dat
  }
  med.bm <- do.call("rbind", median.bm)
  
  ggplot(bm.best.m) +
    geom_line(aes(x = Year, y = bm.stock/1000000)) +
    geom_smooth(method = "lm") +
    #geom_line(data = med.bm, aes(x = Year, y = med.bm/1000000, group = species), linetype = 2) +
    labs(x = "Year", y = "Stock biomass (million kg)") +
    facet_wrap(~interaction(species, MFG , sep = ", ",drop=T),scales='free_y')
  
  #proportion of stocks in community
  com.prop <- ggplot(bm.best.b) +
    geom_line(aes(x = Year, y = prop.bm.stock.eco, group = species, colour = species)) +
    scale_colour_viridis_d() +
    labs(x = "Year", y = "Proportion of community biomass") +
    labs(colour = "Stock") +
    theme(panel.background = element_rect(fill = 'darkgrey'))
  com.prop
  
  com.prop.log <- ggplot(bm.best.b) +
    geom_line(aes(x = Year, y = prop.bm.stock.eco, group = species, colour = species)) +
    scale_colour_viridis_d() +
    labs(x = "Year", y = "Proportion of community biomass (log scale)") +
    scale_y_log10() +
    labs(colour = "Stock") +
    theme(panel.background = element_rect(fill = 'darkgrey'))
  com.prop.log
  
  #lambdas
  lambdas <- do.call("rbind", stocks.lst)
  lambdas <- lambdas |> subset(year %in% yrs)
  lambdas <- na.omit(lambdas)
  for (i in 1:nrow(lambdas)){
    if (lambdas$FG[i] == "MBZ" || lambdas$FG[i] == "MP") lambdas$FG[i] <- "M"
  }
  
  for (i in unique(lambdas$common)){
    dat <- lambdas[lambdas$common == i,]
    print(paste(i))
    print(summary(dat$lambda))
  }
  
  ggplot(lambdas) + geom_point(aes(x = year, y = lambda, group = common)) +
    labs(x = "Year", y = "Realized rate of population growth (λ)") +
    geom_hline(yintercept = 1, linetype = 2) +
    facet_wrap(~interaction(common, FG , sep = ", ",drop=T),scales='free_y')
    
