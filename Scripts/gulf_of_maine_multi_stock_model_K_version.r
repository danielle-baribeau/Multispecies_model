# Here we develop a multi-species model for the North Sea.
#DB - January 2025: Running a test model for GOM/GB stocks (NAFO areas 5 and 6)

#################  Section 1 Loading #################  Section 1 Loading #################  Section 1 Loading  ###############################################
library(tidyverse)
library(GGally)
library(cowplot)
library(ggthemes)
# Set the base plot theme
theme_set(theme_few(base_size = 14))
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

#loop to load functions didn't work for me, so I sourced them instead for now
source("C:/Users/danxb/Desktop/Uni/Grad School Winter 2025/ICM/ICM/Scripts/functions/simple_Lotka_r.r")
source("C:/Users/danxb/Desktop/Uni/Grad School Winter 2025/ICM/ICM/Scripts/functions/simple_forward_sim.r")
source("C:/Users/danxb/Desktop/Uni/Grad School Winter 2025/ICM/ICM/Scripts/functions/forward_project.r")

#load(file = "D:/Github/ICM/Results/model_inputs.Rdata")
dat.loc <- 'C:/Users/danxb/Desktop/Uni/Grad School Winter 2025/ICM/ICM'
repo.loc <- "C:/Users/danxb/Desktop/Uni/Grad School Winter 2025/Multispecies_model"
#loc <- "C:/Users/Owner/Documents/Github/ICM"
load(file = paste0(dat.loc,"/Results/all_cleaned_forward_tune_summaries_fec_nm.Rdata"))
load(file = paste0(dat.loc,"/Results/model_inputs.Rdata"))

########################### End Section 1 Loading  ###################################### End Section 1 Loading ###############################################

########################## Section 2 Parameters ########################## Section 2 Parameters ########################## Section 2 Parameters

n.yrs.proj <- 10 # How many years into the future we are going to project the stocks
n.sims <- 1 # The numbers of simulations to run, keeping low for testing...

# Get the right stocks
#stocks for testing purposes
Stocks<-c("CERT-TRAC_GB_Melanogrammus_Aeglefinus")

#current process for getting target GOM stocks
#Stocks <- names(for.tune.all)
#Stocks <- Stocks[grep("NEFSC", Stocks)]
#Stocks <- c(Stocks, "CERT-TRAC_GB_Melanogrammus_Aeglefinus")
#Stocks <- Stocks[Stocks != "NEFSC-GARMIII_MA_Paralichthys_dentatus"]
#Stocks <- Stocks[Stocks != "NEFSC-GARMIII_SNE- MA_Limanda_ferruginea"]
#Stocks <- Stocks[Stocks != "NEFSC-SAW_USATL_Stenotomus_Chrysops"]
#Stocks <- Stocks[Stocks != "NEFSC_NWA_Scomber _scombrus"]


# So here we are working to get the 'ecosystem' carrying capacity by looking at the total biomass for the GOM stocks we have
# data for over the period of time we have data for all the stocks.
# So here we pull out the data we need to look at total abundance and total biomass in the system by year...
years.gom <- NULL
#holds all years we have for a given stock
vpa.gom <- NULL
#estimated abundance from ICM
bm.gom <- NULL
#biomass per age per year for stocks (unclear what num column is here)
num.gom <- NULL
#metadata on stock, and number of individuals per age per year
waa.gom <- NULL
#metadata on stock, and weight at age of individuals per year (?)
pnm.gom <- NULL
#estimated natural mortality per age from ICM (?)
rem.gom <- NULL
#removals per year for stocks (?)
mx.gom <- NULL
#estimated reproductive output per age from ICM (?)
am.gom <- NULL
#proportion mature at age for stocks
ages.gom <- NULL
#age classes for stocks (or actual ages?)

for(i in  Stocks)
{
  years.gom[[i]] <- years.tmp[[i]]
  #vpa.gom[[i]] <- vpa.tmp[[i]]
  ages.gom[[i]] <- ages.tmp[[i]]
  num.gom[[i]] <- ASR_long |> collapse::fsubset(Stock == i & type == "Num")
  num.gom[[i]] <- num.gom[[i]] |> collapse::fsubset(age != "tot")
  waa.gom[[i]] <- ASR_long|> collapse::fsubset(Stock == i & type == "WA")
  #if(i == "ICES-HAWG_NS_Ammodytes_dubius") waa.ns[[i]]$value <- waa.ns[[i]]$value/1000
  #above line is for NS stocks
  bm.gom[[i]] <- data.frame(Year = num.gom[[i]]$Year,Stock = num.gom[[i]]$Stock,age = num.gom[[i]]$age,
                           bm = num.gom[[i]]$value*waa.gom[[i]]$value,
                           num = num.gom[[i]]$value)
  pnm.gom[[i]] <- 1-exp(-for.tune.all[[i]]$nm.opt)
  mx.gom[[i]] <- for.tune.all[[i]]$fecund.opt
  vpa.gom[[i]] <- for.tune.all[[i]]$res$est.abund
  rem.gom[[i]] <- rem.tmp[[i]]
  rem.gom[[i]]$Stock <- i
  am.gom[[i]] <- am.tmp[[i]]
}
# Combine the biomass and abundance data into a dataframe
bm.tst <- do.call("rbind",bm.gom)
# Look at the biomass and abundance in the ecosystem
  #this removes extra metadata present in bm.tst
  #not familiar with what this does: |>
bm.tot <- bm.tst |> collapse::fgroup_by(Stock,Year) |> 
                    collapse::fsummarize(bm = sum(bm,na.rm=T),
                                         num = sum(num,na.rm=T))
# The 'ecosystem' biomass and numbers
  #combines bm and num of all stocks and reports this by year
eco.bm <- bm.tot |> collapse::fgroup_by(Year) |> 
                    collapse::fsummarize(num = sum(num),bm = sum(bm))


# Now we combine the ecosystem results with the stock biomass's
  #putting bm.tot and eco.bm into one data frame
bm.final <- left_join(bm.tot,eco.bm,by=c("Year"))
names(bm.final) <- c("Stock","Year","bm.stock","num.stock","num.total","bm.total")
# Get the proportion of the total biomass each stock accounts for
  #adding this information into bm.final
bm.final <- bm.final |> collapse::fmutate(bm.prop = bm.stock/bm.total,
                                       num.prop = num.stock/num.total)
bm.final <- bm.final[bm.final$bm.stock > 0,]
bm.final <- as.data.frame(bm.final)
# This gets the average weight of individuals in each stock, we'll need this later to get an approximate exploitation rate
  #adding this information into bm.final
bm.final$avg.weight <- bm.final$bm.stock/bm.final$num.stock

# Now we subset to the years we have data for all the stocks
what.year <- bm.final |> collapse::fgroup_by(Stock) |> collapse::fsummarize(min = min(Year),
                                                                      max = max(Year))
# The years we have data for all stocks
  #storing what.year information into usable format
first.year <- max(what.year$min)
last.year <- min(what.year$max)

# Now we subset the data to these years
bm.best <- bm.final |> collapse::fsubset(Year %in% first.year:last.year & Stock == Stocks[1], bm.total, num.total,Year) 

# Fix: This is not perfect way to get the past exploitation rates as the removals we have here are in numbers
# we'll need to think about whether using the avg.weight is appropriate it does get us 'close', but I think we can do better given the data we have.
  #call removals data
rem.tst <- do.call("rbind",rem.gom)
  #combining other necessary inputs with removals data
fm.dat <- left_join(bm.final,rem.tst,by=c("Stock",'Year'))
# This is where we go from numbers to a biomass and get an exploitation rate in biomass.
fm.dat$exploit <- (fm.dat$rem*fm.dat$avg.weight)/fm.dat$bm.stock

# # Autocorrelation in abundance and biomass time series for the 'ecosystem'
K.cor <- pacf(log(bm.best$num.total))
K.cor.bm <- pacf(log(bm.best$bm.total))
# for GOM stocks, ecosystem biomass is predictable from one year to the next
  #next year's biomass is highly related to current year's biomass
# visualizing the biomass time series
ggplot(bm.best,aes(x=Year,y=bm.total)) + geom_point() + geom_line()
# Numbers time series
ggplot(bm.best,aes(x=Year,y=num.total)) + geom_point() + geom_line()


# So now we can 'statistically represent the abundance and biomass carrying capacity, note we do this on the log scale
  #why do we do this on the log scale?
N.target.sum <- bm.best |> collapse::fsummarise(sd = sd(log(num.total)),
                                            mn = mean(log(num.total)),
                                            med = median(log(num.total)))
# Working on implementing using the biomass as K rather than nubmers...
bm.target.sum <- bm.best |>   collapse::fsummarise(sd = sd(log(bm.total)),
                                               mn = mean(log(bm.total)),
                                               med = median(log(bm.total)))


# So here is where the fun comes in, given we have the autocorrelation, mean and standard deviation of the carrying capacity time series
# We can make a time series with this same characteristics, which we can use to simulation the future carrying capacity of the system.
# We can represent the carrying capacity in numbers or in biomass.  
# FIX: The ICM runs in numbers, carrying capacity is probably better thought of in biomass, think of justifiable ways to get ICM in biomass
# or logical ways to convert result of ICM to biomass (we have the weight at age, but for now doing it this way)
K.devs <- NULL
#not sure what kind of data this is holding
cors <- NULL
K.sims <- NULL
#also not sure about this one; is this the values of K projected forwards?

#fills above empty vectors
for(i in 1:n.sims) 
{
   K.devs[[i]] <- as.data.frame(arima.sim(list(order = c(1,0,0), ar = K.cor$acf[1]), 
                                                                 n = n.yrs.proj,
                                                                 sd = N.target.sum$sd))
  #K.devs[[i]] <- data.frame(x = rep(0,n.yrs.proj)) # Lets see what happens when it's fixed at a mean level.
  K.sims[[i]] <- data.frame(year = 2015:2024,sim = i,num = exp(N.target.sum$med + as.numeric(K.devs[[i]]$x)))
}

#holds K.sims values to be plotted in code below
K.sims.4.plt <- data.frame(do.call('rbind',K.sims))

# Now does this work, lets visualize this...
# FIX: this is not bad, only issue I have is the starting point probably should be tied to the value in 2014, so will need to think about that
 #when running test for GOM, starting point overlaps with black line
# The simulation was also skewing a bit high, switching from using the mean to the median when getting K.sims took care of that.
ggplot(K.sims.4.plt,aes(x=year,y=num)) + 
                        geom_point(aes(group=sim,color=sim,fill=sim)) + geom_line(aes(group=sim,color=sim,fill=sim)) + 
                        geom_line(data=bm.best,aes(x=Year,y=num.total),color='black',linewidth=2)

# FIX: This is a simplification that we'll need to consider, we might want to make this
# be where we figure out how to partition the biomass by trophic level.  But the moment it's ok.

# For the moment, lets just cut the world up so each stock gets a fixed % of K based on their biomass
prop.stock <- bm.final |> collapse::fgroup_by(Stock) |> collapse::fsummarize(prop = median(num.prop))
prop.stock.bm <- bm.final |> collapse::fgroup_by(Stock) |> collapse::fsummarize(prop = median(bm.prop))

# Now we use the total carrying capacity and the proportion of the abundance of each stock gets to get a carrying capacity of each stock
# over the time series.
K.stock <- NULL
#K in abundance for each stock
fm.stock <- NULL
#fishing mortality at K for each stock (?)
tmp <- NULL
#temporary vector that holds precursor for K.stock value in for loop

for(i in 1:n.sims)
{
  tmp <- NULL
  for(s in Stocks)
  {
    # Pull out the proportion 
    prop <- prop.stock |> collapse::fsubset(Stock == s)
    tmp[[s]] <- prop$prop * as.numeric(K.sims[[i]]$num)
    # FIX: Really this should be in biomass too, but for now sticking with using abundance
    if(i == 1) fm.stock[[s]] <- fm.dat |> collapse::fsubset(Stock ==s) |> collapse::fsummarize(mn = median(exploit,na.rm=T),
                                                                                         sd = sd(log(exploit[exploit > 0]),na.rm=T))
  }
  K.stock[[i]] <- tmp
}

# Plots you are welcome to enjoy at your leisure
# ggplot(bm.best) + geom_line(aes(x=Year,y=bm.total),linewidth=2) + ylim(c(0,max(bm.final$bm.total))) + geom_hline(yintercept = mean(bm.best$bm.total))
# ggplot(bm.best) + geom_line(aes(x=Year,y=num.total)) + ylim(c(0,max(bm.final$num.total))) 
# ggplot(bm.best) + geom_line(aes(x=Year,y=bm.prop,group = Stock,color=Stock),linewidth=2) #+ ylim(c(0,0.5))
# ggplot(bm.best) + geom_bar(aes(x=Year,y=bm.prop,fill=Stock),position="stack",stat = 'identity') #+ ylim(c(0,0.5))
# prop.acfs <- NULL
# for(i in Stocks) prop.acfs[[i]] <- pacf(bm.best$bm.prop[bm.best$Stock == i])



############### Section 4 Multi-species model of Gulf of Maine ############### Section 4 Multi-species model of Gulf of Maine ############### Section 4 Multi-species model of Gulf of Maine

# Now we have our carrying capacity for each stock and we can get to business and running a model.
  #unsure how these relate to the simulation
tmp.mx <- NULL
tmp.nm <- NULL
tmp.mat <- NULL
tmp.age <- NULL
tmp.waa <- NULL
mx.dev <- NULL
nm.dev <- NULL
ts.unpack <- NULL
r.unpack <- NULL
Ks <- NULL

# Get the year range, going from the 'last' year to n.yrs.proj in the future, note this will go 1 year less than your intuition because
# we want n.yrs of data, i.e., 20 years is 2000 to 2019, not 2020... )
years <- (last.year+1):(last.year+n.yrs.proj)

# So everything will need to get wrapped up in a simulation loop
for(j in 1:n.sims)
{
  #start time of simulation is recorded here
  st.time <- Sys.time()
  # This is getting the K for the stock for a particular simulation
  Kss <- K.stock[[j]]
  
  res.ts <- NULL
  res.r <- NULL
  for(s in Stocks)
  {
    # Obtain a sample of the correct size from the historic natural mortality and fecundities, Note
    # That in this set up this picks the fecundity and natural mortality from the 'same' year, which is 'fine', but lots of options here.
    # FIX: We could do this sample as a proper statistical sample rather than picking exactly what has been observed in the past.
      #unclear as to what this line is doing, especially the [s]
    samp <- sample(nrow(mx.gom[[s]]),n.yrs.proj,replace=T)
      #I'm confused on what samp is holding; currently looks like this:
      #10 36 48 28 41 16 21 46  1 48 41 23  8 18 28 25 24 18 45 10 17 28 41 30 31
    # Historical natural mortality estimates (using our retrospective simulations)
    prop.nat.mort <- pnm.gom[[s]] 
    prop.nat.mort <- prop.nat.mort[samp,] # Get the sample years.
    
    # maturity ogive
    #age.mat <- am.ns[[s]]
    # Historical fecundity estimates (using our retrospective simulations)
    mx <- mx.gom[[s]] 
    mx <- mx[samp,] # Get the sample years.
    # Ages, everything is calculated back to 'age 0', needed so lotka optimization 'knows' the age it's dealing with
    ages <- 0:(ncol(prop.nat.mort)-1)
    # Carrying capacity for the stocks, currently in numbers, need to think about it in biomass terms IMHO, but not there yet.
    Ks <- Kss[[s]]
    # The number of individuals (using our retrospective simulations)
    vpa.gom  <- bm.final$num.stock[bm.final$Stock == s]
    # Now get the final year abundance
    N.start <- vpa.gom[length(vpa.gom)]
    
    # Run the simulations
    tst <- simp.for.sim(years= years,
                        nm = -(log(1-prop.nat.mort)),
                        ages = ages,
                        rems =  list(fm.stock[[s]]$mn,fm.stock[[s]]$sd), #fm,
                        fecund = mx,
                        N.start = N.start,
                        pop.model = 'bounded_exp', 
                        sim= "project",
                        n.sims = 1,
                        K = Ks,
                        repo='C:/Users/danxb/Desktop/Uni/Grad School Winter 2025/ICM/ICM')
    
# Tidy up the output, Abundance and population growth rates
res.ts[[s]] <- data.frame(tst$Pop[,-2],stock = s,sim= j)
res.r[[s]] <- data.frame(tst$r[,-3],stock=s,sim=j)

  } # end stock loop
  # Unpack the results
  ts.unpack[[j]] <- do.call('rbind',res.ts)
  r.unpack[[j]] <- do.call('rbind',res.r)
  
  # Pop a note when done each simulation
  timer <- Sys.time() - st.time
  print(paste("Simulation ", j))
  print(signif(timer,digits=2))
  
} # end n.sims

# Unpack all the results.
ts.final <- do.call("rbind",ts.unpack)
ts.final$fm <- ts.final$removals/ts.final$abund
r.final <- do.call("rbind",r.unpack)

quants <- ts.final |>  collapse::fgroup_by(years,stock) |> collapse::fsummarize(L.50 = quantile(abund,probs=c(0.25),na.rm=T),
                                                                          med = median(abund,na.rm=T),
                                                                          U.50 = quantile(abund,probs=c(0.75),na.rm=T),
                                                                          fml.50 = quantile(fm,probs=c(0.25),na.rm=T),
                                                                          fm = median(fm,na.rm=T),
                                                                          fmu.50 = quantile(fm,probs=c(0.75),na.rm=T))
# If happy save the 3 objects
saveRDS(object = ts.final,file = paste0(repo.loc,"/Results/NS_projections_",n.sims,"_sims_",min(years),"_to_",max(years),
                                         "_time_series_projections.Rds"))

saveRDS(object = quants,file = paste0(repo.loc,"/Results/NS_projections_",n.sims,"_sims_",min(years),"_to_",max(years),
                                      "_time_series_quantiles.Rds"))

saveRDS(object = r.final,file = paste0(repo.loc,"/Results/NS_projections_",n.sims,"_sims_",min(years),"_to_",max(years),
                                       "_r_projections.Rds"))


# Two simple plots. 
p.sims <- ggplot(ts.final ) + geom_line(aes(x=years,y=abund,group = sim,color=sim),alpha=0.8) +
  facet_wrap(~stock,scales = 'free_y') + ylim(c(0,NA)) + scale_x_continuous(breaks = seq(2010,2200,by=10)) 
#save_plot(paste0("D:/Github/ICM/Figures/NS_sims/NS_all_realizations_climate_starts_at_year_",c.effect,
#                 "_mx_decade_effect_",climate.mx.effect, "_nm_decade_effect_",climate.nm.effect,".png"),p.sims,base_height = 12,base_width = 20)

p.sims.quants <- ggplot(quants) + geom_line(aes(x=years,y=med)) + facet_wrap(~stock,scales = 'free_y') + ylim(c(0,NA)) +
  geom_ribbon(data=quants, aes(x=years,ymax=U.50,ymin = L.50),alpha=0.5,fill='blue',color='blue') 
#save_plot(paste0("D:/Github/ICM/Figures/NS_sims/NS_quantiles_climate_starts_at_year_",c.effect,
#                 "_mx_decade_effect_",climate.mx.effect, "_nm_decade_effect_",climate.nm.effect,".png"),p.sims.quants,base_height = 12,base_width = 20)



