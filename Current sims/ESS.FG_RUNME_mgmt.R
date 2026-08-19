# FIX For the MAESTRO work, lets simulate the ecosystem a fixed number of times, then run the population dynamics on each of these ecosystem
# simulations, instead of simulating the ecosystem for each of the population dynamics scenarios, for example, we can simulate 1,000 
# ecosystems.  Then, for the fishery dynamics, if we are testing 10 scenarios and simulating the population dynamics 100 times, that 
# all happens within a set of 1000 ecosystems scenarios (so each ecosystems would have 1000 population dynamics simulations run on it in this scenario)

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

#if using historical removals, pull these as well
load("C:/Users/BARIBEAUD/Desktop/GitHub/Multispecies_model/Sim Versions/Current sims/summary_hist_removals.RData")
#add in stock names
stock.names <- data.frame(common = c("Atlantic cod", "Haddock", "White hake", "Silver hake", "Pollock", 
                                     "Redfish","American plaice", "Witch flounder", "Atlantic wolffish", 
                                     "Longfin hake","Thorny skate", "Smooth skate", 
                                     "Longhorn sculpin", "Sea raven","Monkfish"),
                          code = eco.stocks)
rem.sum <- merge(rem.sum, stock.names, by = "code")


######################### Run this if you want to use the catch function ######################################
#TO-DO:
#New management plan options:
  #Using a function to fill in management scenario according to default settings
  #specify which scenario you want, and function will fill things in



#this management info will be used to inform catch projections in simulation
## using this in chapter 3 ##
  #user sets exploitation rate for first management cycle, then rate will be updated by the simulation in each assessment
mgmt.scen <- data.frame(stock = eco.stocks, stock.num = c(seq(1,length(eco.stocks))),
                        ex.curr = c(rep(0.1, length(eco.stocks))), u.sd = c(rep(0, length(eco.stocks))),
                        u.min = c(rep(0, length(eco.stocks))), u.max = c(rep(0.4, length(eco.stocks))),
                        assessment.interval = c(rep(3,length(eco.stocks))))

#spec-specific exploitation rates (DONE)
#go off fishery-affected lambdas; fishing has been minimal since 2000
  #projections - directed fishing on top of what's already happening
  #lambdas reflect a low level of bycatch that we can't characterize
#selectivity - do we care?
#get catchability into input time series biomass - matters for the Ks
#fishing scenarios

#to do before Dec 8th:
  #biomass - DONE
  #lambdas - DONE
  #do a run from 2000 onwards - DONE
  #fishing scenarios finalized - DONE

#climate? - implement a vector on the Ks - species-specific penalties (constant throughout projections) 
  #and the total ecosystem (annually varying)
  #does abundance factor heavily into CRIB components? Need to track this

#analysis of area fished post-collapse and effects of fishing in the future?
  
###############################################################################################################

#UPDATES COMPLETED
  #rnorm for u sampling has been added in - user can set sd for this ahead of time
  #u at LRP and USR now can get set by user

#trade-offs between number of functions you have and the understandability of your code


#UPDATES TO DO:
  #plan output plots - DONE
  #move eqn setup to a new function (see notes in trophic_model... file); make 0 and 0.4 default to u.min and u.max
  #stretched beta distribution - lambdas


#test <- trophic.mod(stocks = stock.lst,lambdas= eco.lambdas,n.sims=n.sims,
                   # mgmt = list(mgmt =mgmt.scen),
                    #n.yrs.proj= n.yrs.proj,repo.loc=repo.loc)

# Look at this...
#test$sim.ts

#FOR TESTING
#testing catch function
stocks = stocks
#lambdas= stocks.lst
n.sims=n.sims
#mgmt = list(mgmt =mgmt.scen)
n.yrs.proj= n.yrs.proj
repo.loc=repo.loc

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
