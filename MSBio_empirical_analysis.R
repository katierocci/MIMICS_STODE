library(ggplot2)
library(tidyr)
library(dplyr)
library(car)
library(lmerTest)
library(emmeans)

rm(list = ls())

#load data
MicFG_soil <- read.csv("Example_simulations/Data/MicFG_soil.csv")

#format data
MFG_analysis <- MicFG_soil %>% mutate(log.vwc = log(vwc.avg)) %>% select(perc_decomp_T1, perc_decomp_T2, Copiotroph2, Oligotroph2, C_O, r_K, LIG_N, vwc.avg, log.vwc, MAT,site, plot) %>%
  pivot_longer(1:2, names_to = "time.point", values_to = "perc_decomp")

#standardize data
MFG_stdzd <- MFG_analysis %>% mutate_at(c('C_O', 'LIG_N', 'log.vwc'), ~(scale(.) %>% as.vector))

#calculate effect size
Obs_ES_mod <- lmer(perc_decomp ~ C_O + LIG_N + log.vwc + (1|site/plot), data = MFG_stdzd)
summary(Obs_ES_mod)
Obs_ES <- as.data.frame(fixef(Obs_ES_mod)) #fixed effects coefficients as effect size
Obs_ES$Vars <- rownames(Obs_ES)
colnames(Obs_ES)[1] <- "value"
Obs_ES <- Obs_ES[-1, ]
Obs_ES$mult <- ifelse(Obs_ES$value <0, -1, 1)
Obs_ES$rel_ES <- (abs(Obs_ES$value)/sum(abs(Obs_ES$value))) * 100 * Obs_ES$mult
Obs_ES$Vars <- factor(Obs_ES$Vars, levels = c("log.vwc", "LIG_N", "C_O"), labels = c("Soil moisture", "Litter quality", "Microbial community"))
ggplot(Obs_ES, aes(x=Vars, y=rel_ES)) + geom_bar(stat="identity", fill="black") + geom_text(aes(label=round(rel_ES, digits=1)), color="black", size=5, vjust=1) + geom_text(aes(label=round(rel_ES, digits=1)), color="black", size=5, vjust=-0.4) +
  theme_bw(base_size = 16) + ylab("Realtive effect size") + xlab("") + ylim(-40,63)#+ ggtitle("Effect size estimates for observations")
#checking collinearity
vif(Obs_ES_mod) #all under 1.5

#visually looking at effect size relationships
colorBlind7  <- c("#E69F00", "#56B4E9", "#009E73",
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
MFG_analysis <- MFG_analysis %>% mutate(site=factor(site, levels=c("TREE", "BART", "HARV", "GRSM", "SERC", "TALL", "LENO"))) #MAT order
#copiotroph:oligotroph
ggplot(MFG_analysis, aes(x=C_O, y=perc_decomp*100)) + geom_point(aes(color=site), size=3) + geom_smooth(method="lm", color="black") +
  xlab("Soil copiotroph:oligotroph ratio")+ylab("Percent litter remaining (%)")+ylim(0,100)+theme_bw(base_size = 16) +
  scale_color_manual(values=colorBlind7)
#soil moisture
ggplot(MFG_analysis, aes(x=log.vwc, y=perc_decomp*100)) + geom_point(aes(color=site), size=3) + geom_smooth(method="lm", color="black") +
  xlab("Log of soil volumetric water content (%)")+ylab("Percent litter remaining (%)")+ylim(0,100)+theme_bw(base_size = 16) +
  scale_color_manual(values=colorBlind7)
#litter lignin:N ratio
ggplot(MFG_analysis, aes(x=LIG_N, y=perc_decomp*100)) + geom_point(aes(color=site), size=3) + geom_smooth(method="lm", color="black") +
  xlab("Litter lignin:N ratio")+ylab("Percent litter remaining (%)")+ylim(0,100)+theme_bw(base_size = 16) +
  scale_color_manual(values=colorBlind7)
