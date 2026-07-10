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

# The correlation in the ecosystem biomass trend
K.cor <- pacf(eco.tot.bm.best$bm.eco) #NO SIG ACF STRUCTURE FOR 2000 - 2017

#what about 1995 - 2017? (still after cod moratorium in 1992, so regime shift has likely already been in progress)
yrs <- 2000:2017
#ecosystem biomass just barely AR(1) (not significant again in 1994; sig from 1995-1997(progressively less))
  #sig lost in 1989

#going to continue with 1995 - 2017

#what about autocorrelation in LP - strong correlation between LP and eco, so size of LP sets size of eco?
LP.cor <- pacf(mfg.bm.best$bm.mfg[mfg.bm.best$BFG=="LP"]) #yes! significant correlation at lag 2 (but not lag 1)
M.cor <- pacf(mfg.bm.best$bm.mfg[mfg.bm.best$BFG=="M"]) #no significant correlation for 2000, sig for 1995

# The cross correlation between the ecosystem biomass trend and the MFG biomasses
# All correlated, but m is strongest (matches lowest TL comment from NS code)
K.lp.cor <- ccf(eco.tot.bm.best$bm.eco,mfg.bm.best$bm.mfg[mfg.bm.best$BFG=="LP"]) 
  #cyclical cycle centered around lag 0; sig at lags 0 and -2 only
  #LP grows with eco in years 1-5, but decreases as eco grows in years farther away
K.lb.cor <- ccf(eco.tot.bm.best$bm.eco,mfg.bm.best$bm.mfg[mfg.bm.best$BFG=="LB"])
  #cyclical cycle centered around lag -3; sig for lags -4, -2 and -1 (not lag 0)
  #LB declines with eco around year -3; increases with eco in other years
K.M.cor <- ccf(eco.tot.bm.best$bm.eco,mfg.bm.best$bm.mfg[mfg.bm.best$BFG=="M"])
  #cyclical cycle centered around lag 0; sig at lags 0 and -1
  #M goes up with eco around lag 0; declines as eco goes up in further years
  #very similar structure to K.lp.cor
    #LP likely influences M

# Within FGs...
fg.lp.lb.cor <- ccf(mfg.bm.best$bm.mfg[mfg.bm.best$BFG=="LP"],mfg.bm.best$bm.mfg[mfg.bm.best$BFG=="LB"])
  #cyclical pattern around lag -3; sig at lag -1 only
  #LP declines as LB goes up around year -3; LP and LB are positively correlated in further years
  #very similar in structure to K.lb.cor
    #not sure how to interpret this ecologically
fg.lp.m.cor <- ccf(mfg.bm.best$bm.mfg[mfg.bm.best$BFG=="LP"],mfg.bm.best$bm.mfg[mfg.bm.best$BFG=="M"])
  #cyclical pattern around lag 0; sig at lags -1 and 0
  #LP and M are positively correlated around 0; negatively correlated in further years
  #very similar structure to both K.lp.cor and K.m.cor
    #LP is having strong influence on M
fg.lb.m.cor <- ccf(mfg.bm.best$bm.mfg[mfg.bm.best$BFG=="LB"],mfg.bm.best$bm.mfg[mfg.bm.best$BFG=="M"])
  #cyclical pattern around lag 3; sig at lags 1, 2 and 4
  #LB and M are negatively correlated around lag 3; become positive in later years
  #similar structure to K.lb.cor and fg.lp.lb.cor, just right-justified a few years
    #effects of LP on LB trickle down to M?

#also looking at correlation between MP and MBZ:
fg.mp.mbz.cor <- ccf(fg.bm.best$bm.fg[fg.bm.best$FG == "MP"], fg.bm.best$bm.fg[fg.bm.best$FG == "MBZ"])
  #cyclical pattern around lag 2; sig at lags 0-3
  #MBZ and MP are positively correlated around lag 2; become negatively correlated in later years
  #very similar to fg.lp.m, just right-justified a few years
    #makes ecological sense - LP affects all Ms, then MPs affect MBZs

#what about larges and specific smalls?
fg.mp.lp.cor <- ccf(fg.bm.best$bm.fg[fg.bm.best$FG == "MP"], fg.bm.best$bm.fg[fg.bm.best$FG == "LP"])
  #cyclical pattern around lag 1; sig at lag 1 only (not lag 0)
  #MP and LP are positively correlated near lag 1; become negatively correlated in further years
  #similar structure to K.lp.cor, K.m.cor, fg.lp.m.cor and fg.mp.mbz.cor
fg.mp.lb.cor <- ccf(fg.bm.best$bm.fg[fg.bm.best$FG == "MP"], fg.bm.best$bm.fg[fg.bm.best$FG == "LB"])
  #cyclical pattern around lag -2; sig at lags -1, -2 and -4
  #MP and LB are negatively correlated near lag -2; positively correlated as you get further from lag -2
  #similar structure to fg.lb.m.cor, K.lb.cor and fg.lp.lb.cor
fg.mbz.lp.cor <- ccf(fg.bm.best$bm.fg[fg.bm.best$FG == "MBZ"], fg.bm.best$bm.fg[fg.bm.best$FG == "LP"])
  #cyclical pattern around lag -2; sig at lags -2 and 6 (not 0)
  #pretty much the same as fg.mp.lb, just flipped (positive corr around lag -2, negative as you get further away)
  #similar structure to fg.mp.lp, jsut left-justified a bit
fg.mbz.lb.cor <- ccf(fg.bm.best$bm.fg[fg.bm.best$FG == "MBZ"], fg.bm.best$bm.fg[fg.bm.best$FG == "LB"])
  #cyclical pattern around lag -4; sig at lag -5 only
  #pretty much the same as fg.mp.lb
  #also almost mirror image of fg.mbz.lp

#larges and mediums overall?
eco.L.cor <- ccf(eco.tot.bm.best$bm.eco,bfg.bm.best$bm.bfg[bfg.bm.best$BFG=="L"])
  #cyclical pattern around lag 1; sig at lags 0 and -6
  #positively correlated near lag 1; negatively correlated in further years
  #very close to fg.mp.lp.cor, fg.mp.mbz.cor, K.lp.cor and K.m.cor
M.L.cor <- ccf(bfg.bm.best$bm.bfg[bfg.bm.best$BFG=="M"],bfg.bm.best$bm.bfg[bfg.bm.best$BFG=="L"])
  #cyclical pattern around lag 2 but NOT SIGNIFICANT! Closest is at lag -6
  #supports going with option 1 - using LP, LB and M as BFGs rather than L and M

#for all these cyclical patterns, shift seems to happen every 5-6 years

#looking at proportions:
fg.lp.lb.prop.cor <- ccf(bm.best.m$prop.eco.mfg[bm.best.m$MFG == "LP"], bm.best.m$prop.eco.mfg[bm.best.m$MFG == "LB"][1:n.years])
  #NO SIGNIFICANT CORRELATION! But cyclical pattern evey 3 years
LP.prop <- bm.best.m$prop.eco.mfg[bm.best.m$MFG == "LP"]
LP.prop <- unique(LP.prop)
M.prop <- bm.best.m$prop.eco.mfg[bm.best.m$MFG == "M"]
M.prop <- unique(M.prop)
prop <- data.frame(year = yrs, LP.prop = LP.prop, M.prop = M.prop)

ggplot(prop) + geom_line(aes(x = year, y = LP.prop)) +
  geom_line(data = prop, aes(x = year, y = M.prop), colour = "blue")

ccf(LP.prop, M.prop)

fg.lp.m.prop.cor <- ccf(bm.best.m$prop.eco.mfg[bm.best.m$MFG == "LP"], bm.best.m$prop.eco.mfg[bm.best.m$MFG == "M"])
  #cyclical pattern every 3 years; only significant at lags 0 (negative corr) and 3 (positive corr)
fg.lb.m.prop.cor <- ccf(bm.best.m$prop.eco.mfg[bm.best.m$MFG == "LB"], bm.best.m$prop.eco.mfg[bm.best.m$MFG == "M"][1:n.years])
  #cyclical pattern every 5-10 years; significant at lags -11, -10, (positive corr) -1, 0 and 1 (negative corr)

fg.lb.mbz.prop.cor <- ccf(bm.best$prop.eco.fg[bm.best$FG == "LB"], bm.best$prop.eco.fg[bm.best$FG == "MBZ"][1:n.years]) #~8-year cyclical pattern in positive-negative corr
  #cyclical pattern every 10 years; significant at lag -3 only (positive corr)
fg.lb.mp.prop.cor <- ccf(bm.best$prop.eco.fg[bm.best$FG == "LB"], bm.best$prop.eco.fg[bm.best$FG == "MP"][1:n.years]) #~6-year cyclical pattern in positive-neg
  #cyclical pattern every 5-10 years (similar to lb.mbz.prop); sig at lags -3, -1, 0 and 1
fg.lp.mbz.prop.cor <- ccf(bm.best$prop.eco.fg[bm.best$FG == "LP"], bm.best$prop.eco.fg[bm.best$FG == "MBZ"][1:n.years]) #only sig corr at lag of -3
  #cyclical pattern every 2ish years; sig at lag -3 only (positive corr)
fg.lp.mp.prop.cor <- ccf(bm.best$prop.eco.fg[bm.best$FG == "LP"], bm.best$prop.eco.fg[bm.best$FG == "MP"][1:n.years]) #sig only at lags 3 and 9 (positive - as Lp goes up, so does MP)
  #cyclical pattern every 3ish years; sig at lags 0 (negative corr) and 3 (positive corr)
fg.l.m.prop.cor <- ccf(bm.best.b$prop.eco.bfg[bm.best.b$BFG == "L"], bm.best.b$prop.eco.bfg[bm.best.b$BFG == "M"]) #10-year cyclical pattern in positive and negative corr
  #cyclical pattern every 10 years; sig at lags -15 to -5 (positive), -5 to 5 (negative) and 5 to 15 (positive)
    #lots of significance here!

#ccf.df <- data.frame(lag = seq(-12,12), prop.LP_LB = fg.lp.lb.prop.cor$acf, prop_LP_M = fg.lp.m.prop.cor$acf, 
#prop_LB_M = fg.lb.m.prop.cor$acf, prop_LB_MBZ = fg.lb.mbz.prop.cor$acf, prop..LB_MP = fg.lb.mp.prop.cor$acf, 
#prop_LP_MBZ = fg.lp.mbz.prop.cor$acf, prop_LP_MP = fg.lp.mp.prop.cor$acf,
#prop_M_L <- fg.l.m.prop.cor$acf)


#OK - simulation prep time. Pull out acfs/ccfs that we need for the model
#eco
eco.cor <- pacf(eco.tot.bm.best$bm.eco)
eco.cor.lag.1 <- eco.cor$acf[1]
eco.cor.lag.2 <- eco.cor$acf[2]

#LP/eco
LP.bm <- mfg.bm.best$bm.mfg[mfg.bm.best$BFG=="LP"]
LP.eco <- LP.bm/eco.tot.bm.best$bm.eco
LP.eco.cor <- pacf(LP.eco)
LP.eco.cor.lag.1 <- LP.eco.cor$acf[1]
LP.eco.cor.lag.2 <- LP.eco.cor$acf[2]

#LB/LB+M
  #first, get LB + M
  low.fgs <- c("M", "LB")
  LBM.bm <- mfg.bm.best |> collapse::fsubset(BFG %in% low.fgs) |> collapse::fgroup_by(Year) |> 
    collapse::fsummarize(bm = sum(bm.mfg))
  #next, divide LB by LBM
  LB.bm <- mfg.bm.best |> collapse::fsubset(BFG == "LB")
  LB.LBM <- LB.bm$bm.mfg/LBM.bm$bm
  #then pacf
  LB.LBM.cor <- pacf(LB.LBM)
  LB.LBM.cor.lag.1 <- LB.LBM.cor$acf[1]
  LB.LBM.cor.lag.2 <- LB.LBM.cor$acf[2]
  
#MP/M
bm.mp <- fg.bm.best|> collapse::fsubset(FG == "MP")
bm.tot.med <- mfg.bm.best |> collapse::fsubset(BFG == "M")
MP.M <- bm.mp$bm.fg/bm.tot.med$bm.mfg
MP.M.cor <- pacf(MP.M) 
MP.M.cor.lag.1 <- MP.M.cor$acf[1]
MP.M.cor.lag.2 <- MP.M.cor$acf[2]

mfg.levels <- sort(unique(bm.best.m$MFG))
fg.levels <- sort(unique(bm.best$FG))

sim.K.stock <- NULL
sim.Ks <- NULL
sim.eco.bm <- NULL
sim.LP.eco.prop.bm <- NULL
sim.LB.M.prop.bm <- NULL
sim.MP.M.prop.bm <- NULL
bm.mfg.Ks <- NULL
bm.fg.Ks <- NULL
arima.logit <- NULL

# Get necessary data on logit scale
LP.eco.logit <- logit(LP.eco)
LB.LBM.logit <- logit(LB.LBM)
MP.M.logit <- logit(MP.M)

# Starting values for the ecosystem and the proportions, logit needed for arima models with the proportions
#starting values for eco
start.eco.sim <- eco.tot.bm.best$bm.eco[length(eco.tot.bm.best$bm.eco)]

#starting values for LP/eco
start.LP.eco.sim <- LP.eco[length(LP.eco)]
start.LP.eco.sim.logit <- LP.eco.logit[length(LP.eco.logit)]
#starting values for LB/LB+M
start.LB.LBM.sim <- LB.LBM[length(LB.LBM)]
start.LB.LBM.sim.logit <- LB.LBM.logit[length(LB.LBM.logit)]
#starting values for MP/M
start.MP.M.sim <- MP.M[length(MP.M)]
start.MP.M.sim.logit <- MP.M.logit[length(MP.M.logit)]

# Mean values for LP/eco and LB/LB+M
#mn.LP.bm <- mean(mfg.bm.best$bm.mfg[mfg.bm.best$BFG=="LP"])
#mn.LP.eco <- mean(LP.eco)
#mn.LP.eco.logit <- mean(LP.eco.logit)
#mn.LB.LBM <- mean(LB.LBM)
#mn.LB.LBM.logit <- mean(LB.LBM.logit)


#or use most recent year as the mean...
mn.eco.bm <- start.eco.sim
mn.LP.eco <- start.LP.eco.sim
mn.LP.eco.logit <- start.LP.eco.sim.logit
mn.LB.LBM <- start.LB.LBM.sim
mn.LB.LBM.logit <- start.LB.LBM.sim.logit
mn.MP.M <- start.MP.M.sim
mn.MP.M.logit <- start.MP.M.sim.logit

# Difference between starting value an mean (if we use the most recent year as the 'mean', then this is 0)
start.eco.diff = start.eco.sim-mn.eco.bm
start.LP.eco.diff <- start.LP.eco.sim.logit - mn.LP.eco.logit
start.LB.LBM.diff <- start.LB.LBM.sim.logit - mn.LB.LBM.logit
start.MP.M.diff <- start.MP.M.sim.logit - mn.MP.M.logit

# the standard deviations (make them log scale because we are now sampling from a lognormal dist in the arima)
sd.eco.bm <- sd(log(eco.tot.bm.best$bm.eco))
#didn't log transform logits because of negative values - is this an issue?
sd.LP.eco.logit <- sd((LP.eco.logit))
sd.LB.LBM.logit <- sd((LB.LBM.logit))
sd.MP.M.logit <- sd((MP.M.logit))

#think this is an issue, because simulation isn't starting at value from last year of time series


#BAD VERSION:
for(i in 1:n.sims) 
{
  #browser()
 # The LP K, using the mean of the ecosystem with the correlation observed of the time series.
 # This starts the time series at the last value of the time series, then moves it to the mean value, bam!!  This will be done for each of these arima sims.
  sim.eco.bm[[i]] <- data.frame(bm = c(arima.sim(model =list(ar = eco.cor.lag.1, eco.cor.lag.2),n = n.yrs.proj,n.start=1,start.innov = start.eco.diff/eco.cor.lag.1,
                                                 start.eco.diff/eco.cor.lag.1,  #CHECK THAT THIS IS ALSO LAG 1
                                                 innov = c(1,rlnorm(n.yrs.proj-1,0,sd.eco.bm))) * mn.eco.bm),
                                Years = 1:n.yrs.proj,sim = i) #CHECK WHAT HAPPENS IF YOU PUT MEAN FIRST (LIKE BELOW SIM)
  #might get some super high #s with this new setup
  #also check to make sure we are centered around the median (mean might be skewed high but that is ok - problem with lognormal)
  #GETTING SOME ITERATIONS WITH NEGATIVE BIOMASS PREDICTIONS
                                            #start innov tells oyu how far awat from mean; divide by acf to get standardized distance
                                            #put both lags in as default
                                            #not able to get MA or I from a short time series - "accounted for first 2 AR terms"
  
  #summary of this ARIMA matches summary of input ecosystem biomass very well!
  #pacf(sim.eco.bm[[i]]$bm)

  # So then from my simulated ecosystem I want each FG to get its cut of the biomass
  #start with LP
  sim.LP.eco.prop.bm <-inv.logit(arima.sim(model =list(ar = c(LP.eco.cor.lag.1,LP.eco.cor.lag.2)),n = n.yrs.proj,
                                           n.start =2, start.innov = c(start.LP.eco.diff/LP.eco.cor.lag.1,start.LP.eco.diff/LP.eco.cor.lag.1), 
                                           innov = c(0,rnorm(n.yrs.proj-1,0,sd.LP.eco.logit))) + mn.LP.eco.logit)
  
  #again, very close to median of input data!!
  
  bm.sim.LP <- sim.LP.eco.prop.bm * sim.eco.bm[[i]]$bm
    #this is overshooting input data summary (median 79 mill in sim, 59 mill in data)
  
  # So this is what is left for LB and M
  bm.sim.LBM <- sim.eco.bm[[i]]$bm - bm.sim.LP
  # So then we use the historical split between LB and M (how much of LB + M taken up by LB)
   # so then simulate this split
  sim.LB.M.prop.bm <- inv.logit(arima.sim(model =list(ar = LB.LBM.cor.lag.1, LB.LBM.cor.lag.2),
                                    n = n.yrs.proj,n.start =2, 
                                    start.innov = c(start.LB.LBM.diff/LB.LBM.cor.lag.1,
                                                    start.LB.LBM.diff/LB.LBM.cor.lag.1), 
                                    innov = c(0,rnorm(n.yrs.proj-1,0,sd.LB.LBM.logit))) + mn.LB.LBM.logit)
  #summary is much higher than input data
  
  
  # And now LB gets this proportion of the LB and M biomass
  bm.sim.LB <- bm.sim.LBM * sim.LB.M.prop.bm
  # And M gets the rest, and so the ecosystem biomass is a portion of the whole biomass
  bm.sim.M <- sim.eco.bm[[i]]$bm - bm.sim.LP - bm.sim.LB
  
  #get split between MP and MBZ within M
  sim.MP.M.prop.bm <- inv.logit(arima.sim(model =list(ar = MP.M.cor.lag.1, MP.M.cor.lag.2),
                                            n = n.yrs.proj,n.start =2, 
                                            start.innov = c(start.MP.M.diff/MP.M.cor.lag.1,
                                                            start.MP.M.diff/MP.M.cor.lag.1), 
                                            innov = c(0,rnorm(n.yrs.proj-1,0,sd.MP.M.logit))) + mn.MP.M.logit )
  #get MP biomass
  bm.sim.MP <- bm.sim.M * sim.MP.M.prop.bm[[i]]
  #get MBZ biomass
  bm.sim.MBZ <- bm.sim.M - bm.sim.MP
  
  bm.fg.Ks[[i]] <- data.frame(Years = rep(1:n.yrs.proj,4), sim =i,
                                    bm.fg = c(bm.sim.LB,bm.sim.LP,bm.sim.MBZ,bm.sim.MP),
                                    fg = as.factor(sort(rep(c("LP","LB","MP","MBZ"),n.yrs.proj))),
                                    bm.eco = rep(sim.eco.bm[[i]]$bm,4))
  bm.fg.Ks[[i]]$prop.eco.fg <- bm.fg.Ks[[i]]$bm.fg/bm.fg.Ks[[i]]$bm.eco
  
  #this is just storing ARIMA results to compare with distribution of input values
  if (i == 1){
    arima.logit <- data.frame(sim = rep(i, n.yrs.proj),
                              LP.eco = as.vector(sim.LP.eco.prop.bm), 
                              LB.LBM = as.vector(sim.LB.M.prop.bm), 
                              MP.M = as.vector(sim.MP.M.prop.bm))
  }
  if (i > 1){
    arima.logit <- rbind (arima.logit, data.frame(sim = rep(i, n.yrs.proj),
                            LP.eco = as.vector(sim.LP.eco.prop.bm), 
                              LB.LBM = as.vector(sim.LB.M.prop.bm), 
                              MP.M = as.vector(sim.MP.M.prop.bm)))
  }
  
  # OK, so now we have the FG K values simulated in a 'nice' way. Next we partition these to the stocks
  # Give each stock a proportion of the K in it's ecosystem based on their historical cuts of the K, and include the time series correlation in that.
  # I'm going to build in correlation to their K time series (this could 100% be fishery induced correlation), could also put in 
  # cross correlation for species with multiple stocks, but for now, let's just do the AR1/2 thing with this for the proportion of the FG 
  # biomass each stock gets.
  
  for(f in fg.levels)
  {
    fg.stocks <- unique(bm.best$Stock[bm.best$FG==f])
    n.stock.fg <- length(fg.stocks)
    count =0
    for(s in fg.stocks)
    {
      count = count+1
      # Now get the time series for each stock...
      if(count == 1 ||  n.stock.fg != 2)
      {
        tmp.dat <- bm.best[bm.best$Stock ==s,]
        tmp.cor <- pacf(tmp.dat$prop.bm.stock.fg,plot=F) # Get the correlation, use AR1 and AR2 but no more.
        tmp.cor.lag.1 <- tmp.cor$acf[1]
        tmp.cor.lag.2 <- tmp.cor$acf[2]
        #tmp.beta <- estBetaParams(mean(tmp.dat$prop.bm.tl),sd(tmp.dat$prop.bm.tl)^2)
        # Logit tranform the proportions and do the ARIMA on the logits
        bm.logit <- logit(tmp.dat$prop.bm.stock.fg)
        start.bm.logit <- bm.logit[length(bm.logit)]
        mn.bm.logit <- mean(bm.logit)
        # Use the most recent bm on logit scale as the 'mean' value for the simulation
        mn.bm.logit <- start.bm.logit
        # And the standard deviation
        sd.bm.logit <- sd(bm.logit)
        diff.bm.logit <- start.bm.logit - mn.bm.logit
        
        # Then backtransform and everything will stay positive! Just using the AR1 term for these
        # FIX: SEE above comment for where I'm using the AR2, here using the AR2 would give some poor starting values
        # So I'm not comfy doing that (it works by luck in the above for the NS IMHO.)
        
        tmp.prop.bm <- c(inv.logit(arima.sim(model =list(ar = tmp.cor.lag.1, tmp.cor.lag.2),
                                             n.start = 1, start.innov = c(diff.bm.logit/tmp.cor.lag.1, diff.bm.logit/tmp.cor.lag.1),
                                             n = n.yrs.proj,innov = c(0,rnorm(n.yrs.proj-1,0,sd.bm.logit))) + mn.bm.logit))
        sim.Ks[[s]] <- data.frame(Years = 1:n.yrs.proj, sim = i,
                                       Stock = s, FG = f,
                                       bm.stock = tmp.prop.bm*bm.fg.Ks[[i]]$bm.fg[bm.fg.Ks[[i]]$fg==f])
      } # end the if(count == 1 ||  n.stock.tl != 2)
      # If there are only 2 stocks in a FG, then the second stock get the rest of the FG's biomass
      
      if(count == 2 & n.stock.fg == 2) 
      {
        sim.Ks[[s]] <-  data.frame(Years = 1:n.yrs.proj, sim=i,
                                   Stock = s, FG = f,
                                   bm.stock = (1-tmp.prop.bm)*bm.fg.Ks[[i]]$bm.fg[bm.fg.Ks[[i]]$fg==f])
      } # end the case of just 2 stocks
    } # end the stocks loop
  } # end the fg loop
  sim.K.stock[[i]] <- do.call("rbind",sim.Ks)
  
} # end the simulation loop

sim.K.stocks <- do.call("rbind",sim.K.stock)
sim.fg.K <- do.call("rbind",bm.fg.Ks)
sim.eco.K <- do.call("rbind",sim.eco.bm)

#compare simulations to input data
eco.plt <- sim.eco.K
eco.plt$type <- rep("Simulated", nrow(eco.plt))
eco.input <- data.frame (bm = eco.tot.bm.best$bm.eco, Years = eco.tot.bm.best$Year,
                            sim = c(rep(0, nrow(eco.tot.bm.best))), type = c(rep("Inputs", nrow(eco.tot.bm.best))))

eco.bm.plt <- rbind(eco.plt, eco.input)
eco.bm.plt$bm <- eco.bm.plt$bm/1000000
eco.bm.plt$sim <- as.factor(eco.bm.plt$sim)

sim.eco.bm.plt <- ggplot(eco.bm.plt, aes(x = sim, y = bm, fill = type)) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "white", color = "black") +
  labs(x="Simulation iteration", y = "Ecosystem biomass (millions)", fill = "Data type") +
  scale_fill_manual(values = c("darkgrey", "lightblue"))
sim.eco.bm.plt

sim.eco.bm.plt.2 <- ggplot(eco.bm.plt, aes(x = type, y = bm, fill = type)) +
  geom_violin() +
  geom_boxplot(width = 0.2, fill = "white", color = "black") +
  labs(x="Data type", y = "Ecosystem biomass (millions)") +
  scale_fill_manual(values = c("darkgrey", "lightblue")) +
  theme(legend.position = "none")
sim.eco.bm.plt.2

#repeat with logit arimas
LP.eco.arima <- arima.logit[,-c(3,4)]
LP.eco.arima$type <- c(rep("Simulated", nrow(LP.eco.arima)))
LP.eco.input <- data.frame(sim = c(rep(0,length(LP.eco))), LP.eco = LP.eco, type = c(rep("Input", length(LP.eco))))
LP.eco.plt <- rbind(LP.eco.arima, LP.eco.input)
LP.eco.plt$sim <- as.factor(LP.eco.plt$sim)

LB.LBM.arima <- arima.logit[,-c(2,4)]
LB.LBM.arima$type <- c(rep("Simulated", nrow(LB.LBM.arima)))
LB.LBM.input <- data.frame(sim = c(rep(0,length(LB.LBM))), LB.LBM = LB.LBM, type = c(rep("Input", length(LB.LBM))))
LB.LBM.plt <- rbind(LB.LBM.arima, LB.LBM.input)
LB.LBM.plt$sim <- as.factor(LB.LBM.plt$sim)

MP.M.arima <- arima.logit[,-c(2,3)]
MP.M.arima$type <- c(rep("Simulated", nrow(MP.M.arima)))
MP.M.input <- data.frame(sim = c(rep(0,length(MP.M))), MP.M = MP.M, type = c(rep("Input", length(MP.M))))
MP.M.plt <- rbind(MP.M.arima, MP.M.input)
MP.M.plt$sim <- as.factor(MP.M.plt$sim)

sim.plt <- ggplot(MP.M.plt, aes(x = sim, y = MP.M, fill = type)) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "white", color = "black") +
  labs(x="Simulation iteration", y = "MP/M", fill = "Data type") +
  scale_fill_manual(values = c("darkgrey", "lightblue"))
sim.plt

sim.plt.2 <- ggplot(MP.M.plt, aes(x = type, y = MP.M, fill = type)) +
  geom_violin() +
  geom_boxplot(width = 0.2, fill = "white", color = "black") +
  labs(x="Data type", y = "MP.M") +
  scale_fill_manual(values = c("darkgrey", "lightblue")) +
  theme(legend.position = "none")
sim.plt.2

#compare biomass from simulation to input data
#ecosystem biomass
#get quantiles
q.sim.eco <- sim.eco.K |>  collapse::fgroup_by(Years) |> collapse::fsummarize(L.50 = quantile(bm,probs=c(0.25),na.rm=T),
                                                                                               med = median(bm,na.rm=T),
                                                                                               U.50 = quantile(bm,probs=c(0.75),na.rm=T)) |> ungroup()
#make projection years into actual years
for (i in 1:nrow(q.sim.eco)){
  q.sim.eco$Years[i] <- q.sim.eco$Years[i] + 2017
}
#set up input time series to match layout of projection quantiles
q.eco.inputs <- data.frame(Years = eco.input$Years, type = eco.input$type, bm = eco.input$bm)
#put biomass into a readable format
q.sim.eco$L.50 <- q.sim.eco$L.50/1000000
q.sim.eco$med <- q.sim.eco$med/1000000
q.sim.eco$U.50 <- q.sim.eco$U.50/1000000
q.eco.inputs$bm <- q.eco.inputs$bm/1000000
#plot quantiles and input time series
q.eco.sim.plt <- ggplot(q.eco.inputs, aes(x = Years, y = bm)) +
  geom_line() +
  geom_line(data=q.sim.eco, aes(x=Years, y=L.50), linetype = 2, colour = "darkgrey") +
  geom_line(data=q.sim.eco, aes(x=Years, y=med), linetype = 2) +
  geom_line(data=q.sim.eco, aes(x=Years, y=U.50), linetype = 2, colour = "darkgrey") +
  scale_x_continuous(breaks = seq(2000, 2027, by = 2)) +
  labs(x = "Year", y = "Ecosystem K (biomass in millions)")
  q.eco.sim.plt

#FG biomass
  #get quantiles
  q.sim.fg <- sim.fg.K |>  collapse::fsubset(fg == "MBZ") |>
    collapse::fgroup_by(Years) |> collapse::fsummarize(L.50 = quantile(bm.fg,probs=c(0.25),na.rm=T),
                                                                                med = median(bm.fg,na.rm=T),
                                                                                U.50 = quantile(bm.fg,probs=c(0.75),na.rm=T)) |> ungroup()
  #make projection years into actual years
  for (i in 1:nrow(q.sim.fg)){
    q.sim.fg$Years[i] <- q.sim.fg$Years[i] + 2017
  }
  #set up input time series to match layout of projection quantiles
  fg.bm <- fg.bm.best |> collapse::fsubset(FG == "MBZ")
  fg.inputs <- data.frame(Years = yrs, type = eco.input$type, bm = fg.bm$bm.fg)
  #put biomass into a readable format
  q.sim.fg$L.50 <- q.sim.fg$L.50/1000000
  q.sim.fg$med <- q.sim.fg$med/1000000
  q.sim.fg$U.50 <- q.sim.fg$U.50/1000000
  fg.inputs$bm <- fg.inputs$bm/1000000
  #plot quantiles and input time series
  q.fg.sim.plt <- ggplot(fg.inputs, aes(x = Years, y = bm)) +
    geom_line() +
    geom_line(data=q.sim.fg, aes(x=Years, y=L.50), linetype = 2, colour = "darkgrey") +
    geom_line(data=q.sim.fg, aes(x=Years, y=med), linetype = 2) +
    geom_line(data=q.sim.fg, aes(x=Years, y=U.50), linetype = 2, colour = "darkgrey") +
    scale_x_continuous(breaks = seq(2000, 2027, by = 2)) +
    labs(x = "Year", y = "Medium benthi/zoopiscivore K (biomass in millions)")
  q.fg.sim.plt

# Wrap up the K time series for each simulation
#sim.K.stocks$Species <- substr(sim.K.stocks$Stock,14,100)
sim.K.stocks$Stock <- as.character(sim.K.stocks$Stock)
sim.K.stocks$bm.stock <- sim.K.stocks$bm.stock/1000000
sim.stock.K.plt <- ggplot(sim.K.stocks[sim.K.stocks$sim==1,]) + geom_line(aes(x=Years,y=bm.stock,group=Stock,color=Stock),linewidth=2) + 
                             facet_wrap(~FG) + scale_y_log10(name="Biomass (millions)") + theme(legend.position = 'right') +
                             guides(colour = guide_legend(nrow = 8))
sim.stock.K.plt
save_plot(filename = paste0(repo.loc,"/Figures/Simulation_stock_K.png"),sim.stock.K.plt,base_height = 8,base_width = 11)

sim.fg.K$bm.fg <- sim.fg.K$bm.fg/1000000
sim.tl.K.plt <- ggplot(sim.fg.K) + geom_line(aes(x=Years,y=bm.fg,group=as.factor(sim),color=as.factor(sim))) + 
                      facet_wrap(~fg) + theme(legend.position='none') + 
                      scale_y_log10(name="Biomass (millions)")
sim.tl.K.plt

save_plot(filename = paste0(repo.loc,"/Figures/Simulation_trophic_K.png"),sim.tl.K.plt,base_height = 8,base_width = 11)
sim.eco.K$bm <- sim.eco.K$bm/1000000
sim.eco.K.plt <- ggplot(sim.eco.K) + geom_line(aes(x=Years,y=bm,group=as.factor(sim),color=as.factor(sim))) +
                                 theme(legend.position = 'none') + labs(y = "Biomass (millions")
save_plot(filename = paste0(repo.loc,"/Figures/Simulation_eco_K.png"),sim.eco.K.plt,base_height = 8,base_width = 11)

# Comparing TL and ecosystem K going stock by stock with the trophic level and ecosystem K's that I originally made up
# And it's not perfect, but I think for a first pass this work, they keep the characteristics we want in terms of
# correlation and the K's are quite similar to the original ones. For TL3 it is perfect, for 4 and 5 it can be slightly off
# because I aimed to keep the trophic level having the correlation over focusing on getting the K exactly right
# There could be ways to do both I haven't thought of, but think this is ok for now.
 #tst <- sim.K.stocks |> collapse::fgroup_by(troph.cat,Years,sim) |> collapse::fsummarise(tot.bm = sum(bm.stock))
# tst2 <- left_join(tst,sim.troph.K,by=c("Years","troph.cat","sim"))
# tst2 <- tst2 |> collapse::fgroup_by(c('Years','sim')) |>collapse::fmutate(eco.bm.new = sum(tot.bm)) |> as.data.frame()
# # So they definitely differ stock by stock from the original trophic level splits, but most of the time
# # it is within 20% of the original
# tst2$per.diff <- 100*(tst2$tot.bm - tst2$bm.tl) / tst2$bm.tl
# hist(tst2$per.diff[tst2$troph.cat == 5])
# # They do retain the time series characteristics tho
# pacf(tst2$tot.bm[tst2$troph.cat == 3 & tst2$sim == 1])
# # Ecosystem is within 5% way more than 75% of the time
# summary(100*(tst2$eco.bm.new - tst2$bm.eco)/tst2$bm.eco)
# # And still has characteristics we want
# pacf(tst2$eco.bm.new[tst2$sim == 1][1:n.yrs.proj])
# # OK I think we're good to here, have simulated biomass time series :-)
#ggplot(sim.K.stocks |> collapse::fsubset(sim == 1)) + geom_line(aes(x=Years,y=bm.stock,group=Stock,color=troph.cat)) + scale_y_log10()
# ggplot(bm.best) + geom_line(aes(x=Year,y=bm.stock,group=Stock,color=troph.cat))+ scale_y_log10()

# Fix: This is not perfect way to get the past exploitation rates as the removals we have here are in numbers
# Give we have the database with the age specific removals and age specific weights, this
# should be tweaked to use that data. That said, for the moment this should be 'good enough' 
#rem.tst <- do.call("rbind",rem)
#fm.dat <- left_join(bm.best,rem.tst,by=c("Stock",'Year'))
# This is where we go from numbers to a biomass and get an exploitation rate in biomass.
#fm.dat$exploit <- (fm.dat$rem*fm.dat$avg.weight)/fm.dat$bm.stock
# Extract the LTR results.


############### Section 4 Multi-species model of North Sea ############### Section 4 Multi-species model of North Sea ############### Section 4 Multi-species model of North Sea

# Now we have our carrying capacity for each stock and we can get to business and running a model.

# Get the year range, going from the 'last' year to n.yrs.proj in the future, note this will go 1 year less than your intuition because
# we want n.yrs of data, i.e., 20 years is 2000 to 2019, not 2020... )
#browser()
years <- (last.year+1):(last.year+n.yrs.proj)
# Take the biomass data for the north sea and subset it to the years we have data
bm.mod.yrs <- bm.tst |> collapse::fsubset(Year %in% first.year:last.year)
for(s in stock.eco)
bm.start.year <- bm.mod.yrs |> collapse::fsubset(Year == last.year) |> 
                         collapse::fgroup_by(Stock,trophic,troph.cat,Species) |> 
                         collapse::fsummarise(bm.tot = sum(bm,na.rm=T))
# Get the initial ecosystem biomass..
init.eco.bm <- sum(bm.start.year$bm.tot)
init.tl.bm <- bm.start.year |> collapse::fgroup_by(troph.cat) |> collapse::fsummarise(bm.tl = sum(bm.tot))
init.stock.bm <- bm.start.year
# Get the average weight of the fish in the stocks so we can go from biomass to abundance for the model
# FIX: NOT SURE I NEED THIS ANYMORE This could definitely be done more sophisisticatedly!
#av.wgt <- bm.best |> collapse::fgroup_by(Stock,troph.cat) |> collapse::fsummarise(mn.wgt = mean(avg.weight,na.rm=T))
# FIX: NOT SURE I NEED THIS ANYMORE Let's try getting the most recent year weight to go from biomass to numbers as average may be somewhat misleading
# So here the idea is that the most recent years 
#av.wgt <- bm.best |> dplyr::group_by(Stock,troph.cat) |> filter(row_number() >= (n() ))
#av.wgt <- data.frame(Stock = av.wgt$Stock,troph.cat = as.numeric(av.wgt$troph.cat),mn.wgt = av.wgt$avg.weight)
# For some debugging, if still here you can delete I'm sure
#count = 0

##### DB - NEW STUFF FOR FISHING ##########################################
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