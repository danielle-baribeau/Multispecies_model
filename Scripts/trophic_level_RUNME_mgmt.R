# FIX For the MAESTRO work, lets simulate the ecosystem a fixed number of times, then run the population dynamics on each of these ecosystem
# simulations, instead of simulating the ecosystem for each of the population dynamics scenarios, for example, we can simulate 1,000 
# ecosystems.  Then, for the fishery dynamics, if we are testing 10 scenarios and simulating the population dynamics 100 times, that 
# all happens within a set of 1000 ecosystems scenarios (so each ecosystems would have 1000 population dynamics simulations run on it in this scenario)


n.yrs.proj <- 7 # How many years into the future we are going to project the stocks
n.sims <- 3 # The numbers of simulations to run, keeping low for testing...

dat.loc <- 'C:/Users/bowlbyh/Desktop/Graduate Studies/ICM'
repo.loc <- 'C:/Users/bowlbyh/Desktop/Graduate Studies/Multispecies_model'
load(file = paste0(dat.loc,"/Results/all_cleaned_forward_tune_summaries_no_age_corection_fec_nm.Rdata"))
load(file = paste0(dat.loc,"/Results/model_inputs_no_age_correction.Rdata"))
# load in the trophic model function
source(paste0(repo.loc,"/Scripts/trophic_model_function.R"))

# We've spelt the name of Turbot wrong and the name of arrow-tooth flounder. Tidy up of other names because of spaces and capitialization.
ASR.long$Genus[ASR.long$Genus == 'Scopthalmus'] <- "Scophthalmus" 
ASR.long$Genus[ASR.long$Genus == 'clupea'] <- "Clupea" 
ASR.long$Genus[ASR.long$Genus == ' Hippoglossoides'] <- "Hippoglossoides" 
ASR.long$Genus[ASR.long$Genus == 'Hippoglossoides '] <- "Hippoglossoides" 
ASR.long$Genus[ASR.long$Genus == 'Scomber '] <- "Scomber" 
ASR.long$Genus[ASR.long$Genus == 'Dicentrarchus '] <- "Dicentrarchus" 
ASR.long$Genus[ASR.long$Genus == 'Pollachius '] <- "Pollachius" 
ASR.long$Genus[ASR.long$Genus == 'Sardina '] <- "Sardina" 
ASR.long$Species[ASR.long$Species == 'Aeglefinus'] <- "aeglefinus" 
ASR.long$Species[ASR.long$Species == 'Chrysops'] <- "chrysops" 
ASR.long$Species[ASR.long$Species == ' harengus'] <- "harengus" 
ASR.long$Species[ASR.long$Species == 'stomais'] <- "stomias" 
ASR.long$Species[ASR.long$Species == 'Solea'] <- "solea" 
ASR.long$Gen.Spec <- paste(ASR.long$Genus,ASR.long$Species,sep=' ')


# Get the right stocks, this is from our model runs
Stocks <- names(for.tune.all)
eco.stocks <- Stocks[grep("NS",Stocks)]
# This one doesn't have the data we need
eco.stocks <- eco.stocks[eco.stocks != "ICES-WGHANSA_SP8abd_Sardina _pilchardus"]

# Get the trophic levels
# Trophic levels from Simon Jennings Paper in 2002 for NS...or fishbase if not in Jennings (e.g. lesser sand eel and Turbot)
eco.troph <- data.frame(Stock = eco.stocks,
                        Common = c("Herring","Lesser Sand eel","Sole","Atlantic cod",
                                   "Haddock","European plaice","Norway pout","Saithe",
                                   "Atlantic cod", "Whiting", "Sole", "European plaice","Turbot","Sole"),
                        TL = c(3.8,3.08,5.0,5.2,
                               4.7,4.5,4.2,4.6,
                               5.2,5.3,5.0,4.5,4.4,5.0))
# This gets some of the stock data we need.
stocks <- ASR.long |> collapse::fsubset(Stock %in% eco.stocks)
# Get the trophic levels for the stocks
stocks <- merge(stocks,eco.troph,by="Stock")
# Get the ages for the stocks



# Make it a list
stock.lst <- NULL
for(s in eco.stocks) stock.lst[[s]] <- stocks |> collapse::fsubset(Stock == s)

eco.lambdas <- NULL
for(s in names(for.tune.all)) if(s %in% eco.stocks) eco.lambdas[[s]] <- res.lambda.final[res.lambda.final$Stock ==s,]


######################### Run this if you want to use the catch function ######################################
#this management info will be used to inform catch projections in simulation
  #user sets exploitation rate for first management cycle, then rate will be updated by the simulation in each assessment
mgmt.plan <- data.frame(stock = eco.stocks, stock.num = c(seq(1,length(eco.stocks))),
                        ex.curr = c(rep(0.1, length(eco.stocks))),
                        assessment.interval = c(rep(3,length(eco.stocks))))
###############################################################################################################


test <- trophic.mod(stocks = stock.lst,lambdas= eco.lambdas,n.sims=n.sims,
                    mgmt = list(mgmt =mgmt.plan,er.mn = NULL,er.sd = NULL),
                    n.yrs.proj= n.yrs.proj,repo.loc=repo.loc)

# Look at this...
test$sim.ts

#testing catch function
stocks = stock.lst
lambdas= eco.lambdas
n.sims=n.sims
mgmt = list(mgmt =mgmt.plan,er.mn = NULL,er.sd = NULL)
n.yrs.proj= n.yrs.proj
repo.loc=repo.loc