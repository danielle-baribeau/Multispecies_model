#Chapter 1 plots (872, 572)

#historical lambdas
for (i in 1:nrow(lambdas)){
  if (lambdas$common[i] %in% LP.stocks) lambdas$FG[i] <- "LP"
  if (lambdas$common[i] %in% LB.stocks) lambdas$FG[i] <- "LB"
  if (lambdas$common[i] %in% M.stocks) lambdas$FG[i] <- "M"
}
ggplot(lambdas) +
  geom_boxplot(aes(x = lambda, y = factor(common, levels = c(M.stocks, LB.stocks, LP.stocks)),
                   fill = factor(FG, levels = c("LP", "LB", "M"))),
               colour = "white", alpha = 0.7) +
  scale_x_log10() +
  #scale_fill_manual(values = c(""))
  scale_fill_viridis_d(option = "mako") +
  labs(x = "λ (log scale)", y = "Stock") +
  geom_vline(xintercept = 1, colour = "darkgrey") +
  labs(fill = "Functional group") +
  theme(panel.background = element_rect("darkgrey"),
        legend.position = "top")


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

