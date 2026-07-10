#are stocks close to their Ks?

res.bm.K <- res.sim.df
res.bm.K$ratio <- res.bm.K$bm.stock/res.bm.K$stock.K
res.bm.K <- na.omit(res.bm.K)
res.bm.K.quants <- res.bm.K |>  collapse::fgroup_by(year, fg, stock) |> collapse::fsummarize(L.50 = quantile(bm.K.stock,probs=c(0.25),na.rm=T),
                                                                                             med = median(bm.K.stock,na.rm=T),
                                                                                             U.50 = quantile(bm.K.stock,probs=c(0.75),na.rm=T))


bm.K.stock.plot <- ggplot(res.bm.K.quants[res.bm.K.quants$fg == "L",]) + geom_line(aes(x=year, y=L.50, group=stock), colour = "grey") +
  geom_line(aes(x=year, y=U.50, group=stock), colour = "grey") +
  geom_line(aes(x=year, y=med, group=stock)) +
  geom_hline(yintercept = 1, colour = "blue") +
  facet_wrap(~stock) +  labs(y="Biomass:K")
bm.K.stock.plot

bm.K.stock.plot <- ggplot(res.bm.K.quants[res.bm.K.quants$fg == "M",]) + geom_line(aes(x=year, y=L.50, group=stock), colour = "grey") +
  geom_line(aes(x=year, y=U.50, group=stock), colour = "grey") +
  geom_line(aes(x=year, y=med, group=stock)) +
  geom_hline(yintercept = 1, colour = "blue") +
  facet_wrap(~stock) +  labs(y="Biomass:K")
bm.K.stock.plot

bm.K.com.p <- ggplot(res.bm.K.quants) + geom_line(aes(x=year, y=L.50), colour = "grey") +
  geom_line(aes(x=year, y=U.50), colour = "grey") +
  geom_line(aes(x=year, y=med)) +
  geom_hline(yintercept = 1, colour = "blue") +
  labs(y = "Biomass:K")
bm.K.com.p


#comparing stock biomass to lambda values
res.bm.K <- res.sim.df
res.bm.K$bm.K.stock <- res.bm.K$bm.stock/res.bm.K$stock.K
res.bm.K <- na.omit(res.bm.K)

lambda.K.stock.plot <- ggplot(res.bm.K[res.bm.K$fg == "L",]) +   geom_boxplot(aes(x=lambda, y=bm.K.stock, group=stock)) +
  geom_vline(xintercept = 1, colour = "grey") +
  geom_hline(yintercept = 1, colour = "grey") +
  facet_wrap(~stock) +  labs(y = "Biomass:K") + coord_flip()
lambda.K.stock.plot

lambda.K.stock <- ggplot(res.bm.K[res.bm.K$fg == "L",]) + geom_boxplot(aes(x = lambda, y = bm.K.stock), group = stock) +
  facet_wrap(~stock)
lambda.K.stock


#boxplots of lambdas above/below K
res.bm.K <- res.sim.df
res.bm.K$ratio <- res.bm.K$bm.stock/res.bm.K$stock.K
res.bm.K <- na.omit(res.bm.K)
res.bm.K$test <- NA
for (i in 1:nrow(res.bm.K)){
  #browser()
  if(res.bm.K$ratio[i] > 1)res.bm.K$test[i] <- "Above K"
  if (res.bm.K$ratio[i] < 1) res.bm.K$test[i] <- "Below K"
  if (res.bm.K$ratio[i] == 1) res.bm.K$test[i] <- 0
}
res.bm.K$test <- as.factor(res.bm.K$test)

res.bm.K.test <- res.bm.K |> collapse::fgroup_by(test,stock) |> collapse::fsummarize(lam.mn = mean(lambda)) |> ungroup()

lambda.K.stock <- ggplot(res.bm.K, aes(x = stock, y = lambda, fill = test)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_boxplot() +
  #stat_summary(fun.y="mean", size = 0.1, colour = "darkblue") +
  #geom_point(data = res.bm.K.test, aes(x=stock, y=lam.mn), colour= "gold") +
  stat_summary(fun = "mean", geom = "point", shape = 8,
               size = 1, color = "gold") +
  scale_y_log10(name = "Lambda") +
  labs(fill = "Stock biomass") +
  scale_fill_manual(values = c("lightblue", "darkgrey")) +
  theme(axis.text.x = element_text(angle = 90))
lambda.K.stock

#boxplots of lambda checkpoints (0.4, 0.7)
res.bm.l <- res.sim.df
res.bm.l$ratio <- res.bm.K$bm.stock/res.bm.K$stock.K
res.bm.K <- na.omit(res.bm.K)
res.bm.K$test <- NA
for (i in 1:nrow(res.bm.K)){
  #browser()
  if(res.bm.K$ratio[i] > 1)res.bm.K$test[i] <- "Above K"
  if (res.bm.K$ratio[i] < 1) res.bm.K$test[i] <- "Below K"
  if (res.bm.K$ratio[i] == 1) res.bm.K$test[i] <- 0
}
res.bm.K$test <- as.factor(res.bm.K$test)

res.bm.K.test <- res.bm.K |> collapse::fgroup_by(test,stock) |> collapse::fsummarize(lam.mn = mean(lambda)) |> ungroup()

lambda.K.stock <- ggplot(res.bm.K, aes(x = stock, y = lambda, fill = test)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_boxplot() +
  #stat_summary(fun.y="mean", size = 0.1, colour = "darkblue") +
  #geom_point(data = res.bm.K.test, aes(x=stock, y=lam.mn), colour= "gold") +
  stat_summary(fun = "mean", geom = "point", shape = 8,
               size = 1, color = "gold") +
  scale_y_log10(name = "Lambda") +
  labs(fill = "Stock biomass") +
  scale_fill_manual(values = c("lightblue", "darkgrey")) +
  theme(axis.text.x = element_text(angle = 90))
lambda.K.stock

#The lower left should have points and the upper right shouldn't, so this is good . Can you turn those plots into boxplots showing above and below K. I really wanna see where the median lambda is below K (showing the mean would also be useful).  For Pollock I'm wondering if the stock is declining because everytime the stock goes above K it collapses (i.e., the resampling there is picking lambda's that are so small the stock is collapsing!



#large pop dynamics need to happen in separate place from mediums
#get down to lambdas with larges, then look at mediums

#next steps for larges
#Ks - from L to individual stocks
#population dynamics of each of the species
#compare population dynamics back to K (sum population dynamics to get how much larges are actually using)
#run all these steps in one go
#can add in fishing pressure, get everything set up - but just for larges
#ignore LP and LB for now
#once all the large stuff is working, think about deviates on ecosystem size (environmental change)
#THEN think about mediums

########################################################
quants <- res.sim.df |>  collapse::fgroup_by(year,stock,fg) |> collapse::fsummarize(L.50 = quantile(lambda,probs=c(0.25),na.rm=T),
                                                                                    med = median(lambda,na.rm=T),
                                                                                    U.50 = quantile(lambda,probs=c(0.75),na.rm=T))

#comparing inputs stats to simulation stats to see if sim is capturing trends in data
bm.com.test <- res.sim.df[,c(1,2,8)]
bm.com.test <- unique(bm.com.test)
sim.com.bm <- bm.com.test
quants.sim.com.bm <- sim.com.bm |> collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(bm.com,probs=c(0.25),na.rm=T),
                                                                                     med = median(bm.com,na.rm=T),
                                                                                     U.50 = quantile(bm.com,probs=c(0.75),na.rm=T))

p.sim.in.L <- ggplot(quants[quants$fg == "L",]) +geom_line(aes(x=year, y=L.50, group=stock,), colour = "grey") +
  geom_line(data = quants[quants$fg == "L",], aes(x = year, y = U.50, group = stock), colour = "grey") +
  geom_line(data = quants[quants$fg == "L",], aes(x = year, y = med, group = stock)) +
  geom_line(data = p.lambdas[p.lambdas$FG %in% c("LP", "LB"),], aes(x = year, y = lambda, group = stock), 
            colour = "blue")+
  geom_vline(xintercept = 2017, linetype = 2, colour = "blue") +
  facet_wrap(~stock) +  #scale_y_log10(name="Lambda") + 
  guides(colour = guide_legend(nrow = 5))
p.sim.in.L

ggplot(s.p.lam[s.p.lam$FG%in% c("LP", "LB"),]) + geom_point(aes(x = bm.prop, y = lambda, group = stock)) +
  geom_smooth(aes(x = bm.prop, y = lambda, group = stock), method = "lm") +
  geom_vline(xintercept = 1) +
  facet_wrap(~stock, scale = "free", ncol = 4)

ggplot(s.p.lam[s.p.lam$FG%in% c("MP", "MBZ"),]) + geom_point(aes(x = bm.prop, y = lambda, group = stock)) +
  geom_smooth(aes(x = bm.prop, y = lambda, group = stock), method = "lm") +
  geom_vline(xintercept = 1) +
  facet_wrap(~stock, scale = "free", ncol = 4)

ggplot(na.omit(res.sim.df[res.sim.df$stock == "Longfin hake"& res.sim.df$sim %in% 25:34,])) + geom_line(aes(x=year, y = lambda, group = sim),
                                                               colour = "lightblue") +
  geom_line(aes(x = year, y = (bm.stock/1000000)/(stock.K/1000000), group = sim)) +
  #geom_vline(xintercept = median(p.lambdas$pop.bm[p.lambdas$stock == "Smooth skate"], na.rm = T)) +
  geom_hline(yintercept = 1) +
  facet_wrap(~sim, scale = "free", ncol = 4)


ggplot(na.omit(res.sim.df[res.sim.df$stock == "Longfin hake"& res.sim.df$sim %in% 28,])) + geom_line(aes(x=year, y = lambda, group = sim),
                                                                                                     colour = "lightblue") +
  geom_line(aes(x = year, y = (bm.stock/1000000), group = sim)) +
  #geom_vline(xintercept = median(p.lambdas$pop.bm[p.lambdas$stock == "Smooth skate"], na.rm = T)) +
  geom_hline(yintercept = 1) +
  facet_wrap(~sim, scale = "free", ncol = 4)

ggplot(na.omit(res.sim.df[res.sim.df$stock == "Longfin hake"& res.sim.df$sim == 36,]))+
  geom_line(aes(x=year, y = bm.stock/100000, group = sim)) +
  geom_line(aes(x=year, y = stock.K/100000, group = sim), colour = "purple") +
  facet_wrap(~sim, scale = "free", ncol = 4)


View(res.sim.df[res.sim.df$stock == "Longfin hake" & res.sim.df$sim == 4,])

ggplot(na.omit(res.sim.df[res.sim.df$stock == "Longfin hake"& res.sim.df$sim == 3,])) + geom_line(aes(x=year, y = lambda, group = sim),
                                                                                                       colour = "lightblue") +
  geom_line(aes(x = year, y = (bm.stock/1000000)/(stock.K/1000000), group = sim)) +
  #geom_vline(xintercept = median(p.lambdas$pop.bm[p.lambdas$stock == "Smooth skate"], na.rm = T)) +
  geom_hline(yintercept = 1) +
  facet_wrap(~sim, scale = "free", ncol = 4)

ggplot(na.omit(res.sim.df[res.sim.df$fg == "L" & res.sim.df$sim == 7 & res.sim.df$year != 2017,])) + geom_text(aes(x=(bm.stock/1000000)/(stock.K/1000000), y = lambda, group = stock, label = substr(year, 3, 4)),
                                                                ) + #colour = "lightblue"
  geom_point(data = s.p.lam[s.p.lam$FG %in% c("LP", "LB"),], aes(x = bm.prop, y = lambda, group = stock)) +
  #geom_vline(xintercept = median(p.lambdas$pop.bm[p.lambdas$stock == "Smooth skate"], na.rm = T)) +
  geom_hline(yintercept = 1) +
 # geom_vline(xintercept = 0.5) +
  #scale_x_log10(name = "Biomass/K") + scale_y_log10(name = "Lambda") +
  facet_wrap(~stock, scale = "free", ncol = 4)

ggplot(na.omit(res.sim.df[res.sim.df$fg == "L" & res.sim.df$year != 2017,])) + geom_point(aes(x=(bm.stock/1000000)/(stock.K/1000000), y = lambda, group = stock), colour = "lightblue") + 
  #geom_point(data = s.p.lam[s.p.lam$FG %in% c("LP", "LB"),], aes(x = bm.prop, y = lambda, group = stock)) +
  #geom_vline(xintercept = median(p.lambdas$pop.bm[p.lambdas$stock == "Smooth skate"], na.rm = T)) +
  geom_hline(yintercept = 1) +
  # geom_vline(xintercept = 0.5) +
  scale_x_log10(name = "Biomass/K") + scale_y_log10(name = "Lambda") +
  facet_wrap(~stock, scale = "free", ncol = 4)


ggplot(na.omit(res.sim.df[res.sim.df$fg == "L",])) + geom_point(aes(x=(bm.stock/1000000)/(stock.K/1000000), y = lambda, group = stock, colour = year),
                                                                colour = "lightblue") + 
  geom_point(data = s.p.lam[s.p.lam$FG %in% c("LP", "LB"),], aes(x = bm.prop, y = lambda, group = stock)) +
  #geom_vline(xintercept = median(p.lambdas$pop.bm[p.lambdas$stock == "Smooth skate"], na.rm = T)) +
  geom_hline(yintercept = 1) +
  # geom_vline(xintercept = 0.5) +
  scale_x_log10(name = "Biomass/K") + scale_y_log10(name = "Lambda") +
  facet_wrap(~stock, scale = "free", ncol = 4)


ggplot(na.omit(res.sim.df[res.sim.df$fg == "L",])) + geom_line(aes(x=year, y = lambda, group = stock),
                                                                colour = "darkblue") +
  geom_line(aes(x=year, y = stock.K, group = stock), colour = "lightblue") +
  geom_line(data = s.p.lam[s.p.lam$FG %in% c("LP", "LB"),], aes(x = year, y = lambda, group = stock)) +
  geom_line(data = s.p.lam[s.p.lam$FG %in% c("LP", "LB"),], aes(x = year, y = year.minus.one, group = stock)) +
  #geom_vline(xintercept = median(p.lambdas$pop.bm[p.lambdas$stock == "Smooth skate"], na.rm = T)) +
  geom_hline(yintercept = 1) +
  scale_x_log10(name = "Year") + scale_y_log10(name = "Lambda and K") +
  facet_wrap(~stock, scale = "free", ncol = 4)

ggplot(na.omit(res.sim.df[res.sim.df$fg == "L",])) + geom_line(aes(x=year, y = bm.stock/stock.K, group = sim),
                                                               colour = "lightblue", alpha = 0.2) +
  geom_line(aes(x=year, y = lambda, group = sim), colour = "darkblue", alpha = 0.2) +
  geom_line(data = s.p.lam[s.p.lam$FG %in% c("LP", "LB"),], aes(x = year, y = lambda, group = stock), colour = "darkblue", alpha = 0.2) +
  geom_line(data = s.p.lam[s.p.lam$FG %in% c("LP", "LB"),], aes(x = year, y = bm.prop, group = stock), colour = "lightblue", alpha = 0.2) +
  #geom_vline(xintercept = median(p.lambdas$pop.bm[p.lambdas$stock == "Smooth skate"], na.rm = T)) +
  geom_hline(yintercept = 1) +
  scale_x_log10(name = "Year") + scale_y_log10(name = "Lambda and K") +
  facet_wrap(~stock, scale = "free", ncol = 4)

  
  test <- res.sim.df |> group_by(stock) |> subset(bm.stock/stock.K < 0.5) |> summarize(sd.lam = sd(lambda)) |> ungroup()
  test <- test |> collapse::fsubset(stock %in% L.stocks)
  test


ggplot(na.omit(res.sim.df[res.sim.df$fg == "M",])) + geom_point(aes(x=((bm.stock/1000000)/(stock.K/1000000)), y= lambda, group = stock),
                                         colour = "lightblue") +
  geom_point(data = s.p.lam[s.p.lam$FG %in% c("MP", "MBZ"),], aes(x = bm.prop, y = lambda, group = stock)) +
  #geom_vline(xintercept = median(p.lambdas$pop.bm[p.lambdas$stock == "Smooth skate"], na.rm = T)) +
  geom_hline(yintercept = 1) +
  scale_x_log10(name = "Biomass/K") + scale_y_log10(name = "Lambda") +
  facet_wrap(~stock, scale = "free", ncol = 4)

ggplot() + geom_point(data = lambdas[lambdas$FG %in% c("MP", "MBZ"),], aes(x = year.minus.one, y = lambda, group = common)) +
  facet_wrap(~common, scale = "free", ncol = 4)

ggplot(s.p.lam) + geom_point(data = s.p.lam, aes(x = bm.prop, y = log(lambda), group = stock)) +
  facet_wrap(~stock, scale = "free", ncol = 5)

colnames(pred.lam)[colnames(pred.lam) == "Stock"] <- "stock"
ggplot(pred.lam) + geom_point(aes(x = bm.prop, y = lambda, group = stock)) +
  geom_point(data = s.p.lam, aes(x = bm.prop, y = lambda, group = stock),
             colour = "lightblue") +
  facet_wrap(~stock, scale = "free", ncol = 5)


ggplot(na.omit(res.sim.df[res.sim.df$stock == "Witch flounder",])) + geom_point(aes(x=bm.stock/1000000, y = lambda),
                                         colour = "lightblue") +
  #geom_point(data = p.lambdas[p.lambdas$stock == "Witch flounder",], aes(x = year.minus.one/1000000, y = lambda)) +
  #geom_vline(xintercept = median(p.lambdas$pop.bm[p.lambdas$stock == "Smooth skate"], na.rm = T)) +
  #geom_hline(yintercept = 1) +
  scale_x_log10(name = "Biomass") + scale_y_log10(name = "Lambda") +
  facet_wrap(~stock, scale = "free", ncol = 5)


ggplot(na.omit(res.sim.df[res.sim.df$sim == 2,])) + geom_line(aes(x = year, y = lambda, group = sim)) +
  geom_hline(yintercept = 1) +
  facet_wrap(~stock, scale = "free", ncol = 5)

ggplot(na.omit(res.sim.df)) + geom_point(aes(x = bm.stock/1000000, y = lambda, group = stock, colour = sim)) +
  geom_hline(yintercept = 1) +
  labs(x = "Biomass (millions") +
  facet_wrap(~stock, scale = "free", ncol = 5)

ggplot(na.omit(res.sim.df[res.sim.df$sim == 1 & res.sim.df$fg == "M",])) + geom_line(aes(x = year, y = lambda, group = stock)) +
  geom_hline(yintercept = 1) +
  facet_wrap(~stock, scale = "free")

ggplot(p.lambdas[p.lambdas$FG %in% c("LP", "LB"),]) + geom_point(aes(x = year.minus.one, y = lambda, group = stock)) +
  geom_hline(yintercept = 1) +
  geom_smooth(aes(x = year.minus.one, y = lambda), method = "lm") +
  geom_line(aes(x = mn.in, y = lambda, group = stock)) +
  facet_wrap(~stock, scale = "free") + scale_y_log10(name = "l")

ggplot(p.lambdas[p.lambdas$FG %in% c("MP", "MBZ"),]) + geom_line(aes(x = year, y = lambda, group = stock)) +
  geom_hline(yintercept = 1) +
  geom_line(aes(x = mn.in, y = lambda, group = stock)) +
  facet_wrap(~stock, scale = "free") + scale_y_log10(name = "l")


p.sim.in.M <- ggplot(quants[quants$fg == "M",]) +geom_line(aes(x=year, y=L.50, group=stock,), colour = "grey") +
  geom_line(data = quants[quants$fg == "M",], aes(x = year, y = U.50, group = stock), colour = "grey") +
  geom_line(data = quants[quants$fg == "M",], aes(x = year, y = med, group = stock)) +
  geom_line(data = p.lambdas[p.lambdas$FG %in% c("MP", "MBZ"),], aes(x = year, y = pop.bm, group = stock), 
            colour = "blue")+
  geom_vline(xintercept = 2017, linetype = 2, colour = "blue") +
  facet_wrap(~stock) +  scale_y_log10(name="Biomass") + 
  guides(colour = guide_legend(nrow = 5))
p.sim.in.M

quants.sim.com.bm$L.50 <- quants.sim.com.bm$L.50/1000000
quants.sim.com.bm$U.50 <- quants.sim.com.bm$U.50/1000000
quants.sim.com.bm$med <- quants.sim.com.bm$med/1000000

input.com.bm <- data.frame(year = 2000:2017, bm = com.bm/1000000)

p.sim.in.com <- ggplot(quants.sim.com.bm) +geom_line(aes(x=year, y=L.50), colour = "grey") +
  geom_line(data = quants.sim.com.bm, aes(x = year, y = U.50), colour = "grey") +
  geom_line(data = quants.sim.com.bm, aes(x = year, y = med)) +
  geom_line(data = input.com.bm, aes(x = year, y = bm), 
            colour = "blue")+
  geom_vline(xintercept = 2017, linetype = 2, colour = "blue") #+
#scale_y_log10(name="Biomass (millions)")
p.sim.in.com

#comparing inputs stats to simulation stats to see if sim is capturing trends in data
bm.l.test <- res.sim.df[,c(1,2,4,6)]
#bm.l.test <- unique(bm.l.test)
sim.l.bm <- bm.l.test
quants.sim.l.bm <- sim.l.bm |> collapse::fgroup_by(year, fg) |> collapse::fsummarize(L.50 = quantile(bm.fg,probs=c(0.25),na.rm=T),
                                                                                     med = median(bm.fg,na.rm=T),
                                                                                     U.50 = quantile(bm.fg,probs=c(0.75),na.rm=T))
input.l.bm <- data.frame(year = 2000:2017, fg = bm.best.b$BFG, bm = bm.best.b$bm.bfg)

p.sim.in.l <- ggplot(quants.sim.l.bm) +geom_line(aes(x=year, y=L.50, group = fg), colour = "grey") +
  geom_line(data = quants.sim.l.bm, aes(x = year, y = U.50, group = fg), colour = "grey") +
  geom_line(data = quants.sim.l.bm, aes(x = year, y = med, group = fg)) +
  geom_line(data = input.l.bm, aes(x = year, y = bm, group = fg), 
            colour = "blue")+
  geom_vline(xintercept = 2017, linetype = 2, colour = "blue") +
  facet_wrap(~fg)
p.sim.in.l

p.dat <- res.sim.df[res.sim.df$stock == "Pollock",]
p.in <- bm.best.b[bm.best.b$species == "Pollock",]

quants.p.prop.fg <- p.dat |> collapse::fgroup_by(year) |> collapse::fsummarize(L.50 = quantile(prop.stock.fg.bm,probs=c(0.25),na.rm=T),
                                                                               med = median(prop.stock.fg.bm,na.rm=T),
                                                                               U.50 = quantile(prop.stock.fg.bm,probs=c(0.75),na.rm=T))

p.pollock <- ggplot(quants.p.prop.fg) + geom_line(aes (x = year, y = L.50), colour = "grey") +
  geom_line(data = quants.p.prop.fg, aes(x = year, y = U.50), colour = "grey") +
  geom_line(data = quants.p.prop.fg, aes(x = year, y = med)) +
  geom_line(data = p.in, aes(x = Year, y = prop.bm.stock.bfg), 
            colour = "blue")+
  geom_vline(xintercept = 2017, linetype = 2, colour = "blue") +
  labs(y = "Proportion of FG biomass")
p.pollock

#proportion of K assigned to each stock over the projection
quants.p.stocks <- na.omit(res.sim.df) |> collapse::fgroup_by(year, stock, fg) |> collapse::fsummarize(L.50 = quantile(stock.K,probs=c(0.25),na.rm=T),
                                                                                                       med = median(stock.K,na.rm=T),
                                                                                                       U.50 = quantile(stock.K,probs=c(0.75),na.rm=T))
bm.stocks.in <- bm.best.b
colnames(bm.stocks.in)[colnames(bm.stocks.in) == "species"] <- "stock"

p.stocks <- ggplot(quants.p.stocks) + geom_line(aes (x = year, y = L.50, group = stock), colour = "grey") +
  geom_line(data = quants.p.stocks, aes(x = year, y = U.50, group = stock), colour = "grey") +
  geom_line(data = quants.p.stocks, aes(x = year, y = med, group = stock)) +
  geom_line(data = quants, aes (x = year, y = L.50, group = stock), colour = "lightblue") +
  geom_line(data = quants, aes(x = year, y = U.50, group = stock), colour = "lightblue") +
  geom_line(data = quants, aes(x = year, y = med, group = stock), colour = "darkblue") +
  geom_line(data = bm.stocks.in, aes(x = Year, y = bm.stock, group = stock), 
            colour = "blue")+
  geom_vline(xintercept = 2017, linetype = 2, colour = "blue") +
  scale_y_log10(name = "Stock K & biomass") +
  facet_wrap(~stock)
p.stocks


######################################################## LAMBDAS
summary(lambdas$lambda[lambdas$common == "Pollock"])
summary(res.sim.df$lambda[res.sim.df$stock == "Pollock"])


quants.l <- res.sim.df |>  collapse::fgroup_by(year,stock,fg) |> collapse::fsummarize(L.50 = quantile(lambda,probs=c(0.25),na.rm=T),
                                                                                      med = median(lambda,na.rm=T),
                                                                                      U.50 = quantile(lambda,probs=c(0.75),na.rm=T))
quants.l <- na.omit(quants.l)

p.sim.in.L.l <- ggplot(quants.l[quants.l$fg == "L",]) +geom_line(aes(x=year, y=L.50, group=stock,), colour = "grey") +
  geom_line(data = quants.l[quants.l$fg == "L",], aes(x = year, y = U.50, group = stock), colour = "grey") +
  geom_line(data = quants.l[quants.l$fg == "L",], aes(x = year, y = med, group = stock)) +
  geom_line(data = p.lambdas[p.lambdas$FG %in% c("LP", "LB"),], aes(x = year, y = lambda, group = stock), 
            colour = "blue")+
  geom_vline(xintercept = 2017, linetype = 2, colour = "blue") +
  facet_wrap(~stock) +  labs(y = "Lambda")
p.sim.in.L.l

p.sim.in.M.l <- ggplot(quants.l[quants.l$fg == "M",]) +geom_line(aes(x=year, y=L.50, group=stock,), colour = "grey") +
  geom_line(data = quants.l[quants.l$fg == "M",], aes(x = year, y = U.50, group = stock), colour = "grey") +
  geom_line(data = quants.l[quants.l$fg == "M",], aes(x = year, y = med, group = stock)) +
  geom_line(data = p.lambdas[p.lambdas$FG %in% c("MP", "MBZ"),], aes(x = year, y = lambda, group = stock), 
            colour = "blue")+
  geom_vline(xintercept = 2017, linetype = 2, colour = "blue") +
  facet_wrap(~stock) +  labs(y = "Lambda")
p.sim.in.M.l


pollock.out <- res.sim.df[res.sim.df$stock == "Pollock" & res.sim.df$year %in% c(2017, 2018),]

ggplot(input.l.bm, aes(x = year, y = bm)) + geom_line() + geom_point(data = pollock.out, aes(x=year, y = bm.stock))

#low and high bm in input data
stock.names <- c(L.stocks, M.stocks)
low.high.lam <- NULL

for (s in c(L.stocks,M.stocks)){
  #browser()
  s.lam <- lambdas[lambdas$common == s,]
  low.high <- 1
  low.vs.high.bm <- low.high * median(s.lam$year.minus.one)
  s.lam$low.vs.high.bm <- low.vs.high.bm
  s.lam$low.or.high <- NA
  for (i in 1:nrow(s.lam)){
    #browser()
    s.lam$low.or.high[i] <- ifelse (s.lam$year.minus.one[i] >= low.vs.high.bm, "H", "L")
  }
  low.high.lam <- rbind(low.high.lam, s.lam)
}

low.high.lam$low.or.high <- as.factor(low.high.lam$low.or.high)

low.high.lam.med <- low.high.lam |> collapse::fgroup_by(common,low.or.high) |> collapse::fsummarize(med.lam = median(lambda)) |> ungroup()

low.high.L <- ggplot(low.high.lam[low.high.lam$FG %in% c("LP", "LB"),]) + 
  geom_point(aes(x = year.minus.one/1000000, y = lambda, group = common, colour = low.or.high)) +
  #geom_line(aes(x = year, y = low.vs.high.bm, group = common)) +
  geom_line(aes(x = low.vs.high.bm/1000000, y = lambda, group = common), colour = "blue") +
  facet_wrap(~common, scale = "free") +  #scale_x_log10(name="Biomass (millions of kg)") + 
  scale_colour_manual(values = c("black", "darkgrey"))
low.high.L

res.p.df.lam <- na.omit(res.sim.df[res.sim.df$stock == "Pollock",])
res.p.df.lam$low.or.high <- "N"

for (x in 1:nrow(res.p.df.lam)){
  med.test <- res.p.df.lam$bm.stock[x]/med.p
  if (med.test > 1){
    res.p.df.lam$low.or.high[x] <- "H"
  }
  if (med.test < 1){
    res.p.df.lam$low.or.high[x] <- "L"
  }
}

low.high.M <- ggplot(low.high.lam[low.high.lam$FG %in% c("MP", "MBZ"),]) + 
  geom_point(aes(x = year.minus.one/1000000, y = lambda, group = common, colour = low.or.high)) +
  #geom_line(aes(x = year, y = low.vs.high.bm, group = common)) +
  geom_line(aes(x = low.vs.high.bm/1000000, y = lambda, group = common), colour = "blue") +
  facet_wrap(~common, scale = "free") +  #scale_x_log10(name="Biomass (Millions of kg)") + 
  scale_colour_manual(values = c("black", "darkgrey"))
low.high.M


#################################################
#input biomasses - boxplot
bm.boxplot <- ggplot(lambdas, aes(x = common, y = year.minus.one)) +
  geom_boxplot() +
  scale_y_log10((name = "Biomass")) +
  stat_summary(fun = "mean", geom = "point", shape = 8,
               size = 1, color = "gold") +
  theme(axis.text.x = element_text(angle = 90))
bm.boxplot


#################################################
p.sims <- ggplot(res.sim.df[res.sim.df$sim == 6,]) + geom_line(aes(x=year,y=bm.stock,group = sim,color=sim),alpha=0.8) +
  facet_wrap(~stock,scales = 'free_y') + 
  geom_hline(yintercept = median(p.lambdas$pop.bm[p.lambdas$stock == "Pollock"], na.rm = T)) +
  scale_y_log10(name = "Biomass") + 
  theme(legend.position = 'none') 
p.sims

p.sims.quants <- ggplot(quants) + geom_line(aes(x=year,y=med,group=stock,color=stock),linewidth=2) + 
  geom_line(data = quants, aes(x = year, y = L.50, group = stock, colour = stock), linetype = 1) +
  geom_line(data = quants, aes(x = year, y = U.50, group = stock, colour = stock), linetype = 1) +
  facet_wrap(~stock) +  scale_y_log10(name="Biomass") +   theme(legend.position = 'none') +
  guides(colour = guide_legend(nrow = 5))
p.sims.quants

redfish.lam <- res.sim.df$lambda[res.sim.df$stock == "Redfish"]
redfish.dat <- res.sim.df[res.sim.df$stock == "Redfish",]
cod.dat <- res.sim.df[res.sim.df$stock == "Atlantic cod",]
ss.dat <- res.sim.df[res.sim.df$stock == "Smooth skate",]
ts.dat <- res.sim.df[res.sim.df$stock == "Thorny skate",]

wolf.dat <- res.sim.df[res.sim.df$stock == "Atlantic wolffish",]
wolf.lam <- lambdas[lambdas$code == 50,]
p.dat <- res.sim.df[res.sim.df$stock == "Pollock",]

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
q.eco.inputs <- data.frame(Years = eco.input$Years, type = eco.input$type, bm = com.bm)
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
  scale_x_continuous(breaks = seq(2000, 2117, by = 20)) +
  labs(x = "Year", y = "Community K (biomass in millions)")
q.eco.sim.plt

#functional groups
q.sim.fg <- sim.fg.K |>  collapse::fsubset(fg == "L") |>
  collapse::fgroup_by(Years) |> collapse::fsummarize(L.50 = quantile(bm.fg,probs=c(0.25),na.rm=T),
                                                     med = median(bm.fg,na.rm=T),
                                                     U.50 = quantile(bm.fg,probs=c(0.75),na.rm=T)) |> ungroup()
#make projection years into actual years
for (i in 1:nrow(q.sim.fg)){
  q.sim.fg$Years[i] <- q.sim.fg$Years[i] + 2017
}
#set up input time series to match layout of projection quantiles
fg.bm <- fg.bm.best |> collapse::fsubset(FG == "L")
fg.inputs <- data.frame(Years = eco.input$Years, type = eco.input$type, bm = L.bm)
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
  labs(x = "Year", y = "Large K (biomass in millions)")
q.fg.sim.plt

#good to look at lambdas versus biomass to figure out how many low/high bm years



##################################################
