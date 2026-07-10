## MPM FOR EASTERN SCOTIAN SHELF - FUNCTIONAL GROUPS ##

n.yrs.proj <- 10 # How many years into the future we are going to project the stocks
n.sims <- 10 # The numbers of simulations to run

dat.loc <- 'C:/Users/BARIBEAUD/Desktop/GitHub/ESS-Model-Setup/Feb2026_Data/No 5s'
repo.loc <- 'C:/Users/BARIBEAUD/Desktop/GitHub/ESS_Test_Run'
load(file = paste0(dat.loc,"/lambda_2024.RData"))
# load in the FG model function
source(paste0(repo.loc,"/ESS_FN.R"))
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

#this management info will be used to inform catch projections in simulation
#user sets exploitation rate for first management cycle, then rate will be updated by the simulation in each assessment
mgmt.scen <- data.frame(stock = eco.stocks, stock.num = c(seq(1,length(eco.stocks))),
                        ex.curr = c(rep(0.1, length(eco.stocks))), u.sd = c(rep(0, length(eco.stocks))),
                        u.min = c(rep(0, length(eco.stocks))), u.max = c(rep(0.4, length(eco.stocks))),
                        assessment.interval = c(rep(3,length(eco.stocks))))


#NOTE - SET SEED EVENTUALLY


#run sims
test <- trophic.mod(stocks = stock.lst,lambdas= eco.lambdas,n.sims=n.sims,
                    mgmt = list(mgmt =mgmt.scen),
                    n.yrs.proj= n.yrs.proj,repo.loc=repo.loc)

# Look at this...
test$sim.ts

#testing catch function
stocks = stocks
#lambdas= stocks.lst
n.sims=n.sims
mgmt = list(mgmt =mgmt.scen)
n.yrs.proj= n.yrs.proj
repo.loc=repo.loc
