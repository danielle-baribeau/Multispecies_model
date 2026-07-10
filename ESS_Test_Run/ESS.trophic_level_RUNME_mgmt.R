# FIX For the MAESTRO work, lets simulate the ecosystem a fixed number of times, then run the population dynamics on each of these ecosystem
# simulations, instead of simulating the ecosystem for each of the population dynamics scenarios, for example, we can simulate 1,000 
# ecosystems.  Then, for the fishery dynamics, if we are testing 10 scenarios and simulating the population dynamics 100 times, that 
# all happens within a set of 1000 ecosystems scenarios (so each ecosystems would have 1000 population dynamics simulations run on it in this scenario)


n.yrs.proj <- 10 # How many years into the future we are going to project the stocks
n.sims <- 10 # The numbers of simulations to run, keeping low for testing...

dat.loc <- 'C:/Users/BARIBEAUD/Desktop/GitHub/ESS_Test_Run/Inputs/cn.no.outliers.dog.scul'
repo.loc <- 'C:/Users/BARIBEAUD/Desktop/GitHub/ESS_Test_Run'
load(file = paste0(dat.loc,"/38_inputs_cn_ed.Rdata"))
# load in the trophic model function
source(paste0(repo.loc,"/E_trophic_model_function_testing_mgmt.R"))
#rename inputs
lam <- lam.ed.ts

# Get the right stocks, this is from our model runs
eco.stocks <- unique(lam$code)

# Get the trophic levels
# Trophic levels
eco.troph <- data.frame(common = c("Atlantic cod", "Haddock", "White hake", "Silver hake", "Pollock", 
                                 "Redfish","American plaice", "Witch flounder", "Atlantic wolffish", 
                                 "Longfin hake","Thorny skate", "Smooth skate", "Spiny dogfish", 
                                 "Longhorn sculpin", "Sea raven","Monkfish"),
                        code = eco.stocks,
                        TL = c(4.1, 4.0, 4.3, 4.5, 4.3, 4.2, 4.1,
                               3.2, 3.6, 3.2, 4.2, 3.5, 4.4, 3.6, 4.5, 4.5))
#add uncertainty into TL - randomly sample around absolute value
#eco.troph$TL.est <- NULL
#for (i in 1:nrow(eco.troph)) eco.troph$TL.est[i] <- rnorm(1, eco.troph$TL[i], 0.3)


#add trophic levels to biomass inputs
stocks <- merge(lam, eco.troph, by = "code")

# Make it a list
#stock.lst <- NULL
#for(s in 1:length(eco.stocks)) stock.lst[[s]] <- stocks |> collapse::fsubset(code == eco.stocks[s])

eco.lambdas <- NULL
for(s in 1:length(eco.stocks)) eco.lambdas[[s]] <- lam |> collapse::fsubset(code == eco.stocks[s])


######################### Run this if you want to use the catch function ######################################
#TO-DO:
#New management plan options:
  #Using a function to fill in management scenario according to default settings
  #specify which scenario you want, and function will fill things in



#this management info will be used to inform catch projections in simulation
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
  #do a run from 2000 onwards
  #fishing scenarios finalized

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
  #code default management scenarios - IP
  #plan output plots
  #move eqn setup to a new function (see notes in trophic_model... file); make 0 and 0.4 default to u.min and u.max
  #set up uncertainty in trophic level assignment
  #stretched beta distribution - lambdas


test <- trophic.mod(stocks = stock.lst,lambdas= eco.lambdas,n.sims=n.sims,
                    mgmt = list(mgmt =mgmt.scen),
                    n.yrs.proj= n.yrs.proj,repo.loc=repo.loc)

# Look at this...
test$sim.ts

#testing catch function
stocks = stocks
lambdas= eco.lambdas
n.sims=n.sims
mgmt = list(mgmt =mgmt.scen)
n.yrs.proj= n.yrs.proj
repo.loc=repo.loc
