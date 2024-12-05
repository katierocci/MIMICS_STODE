#using Derek's litterbag code to run MSBio litter decomp simulations

library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(rwa)
library(lmerTest)
library(randomForest)
library(car)
library(ranger)
library(caTools)
library(vip)


source("functions/MIMICS_sim_litterbag.R")
source("functions/MIMICS_calc_steady_state_pools.R")
source("functions/calc_Tpars.R")
source("Parameters/MIMICS_parameters_sandbox_20231129.R")
source("functions/RXEQ.R")
source("functions/MC_parameterization/set_parameter_defaults.R")



#-------------------------------
#Using MSBio data
#-------------------------------
#load site data
MSBio <- read.csv("Example_simulations/Data/Site_annual_clim_final.csv")
#match input data structure
#AGNPP (litterfall here) should be in grams Dry Weight (gDW) not gC! multiply by 2 here to remedy
MSBio2 <- MSBio %>% mutate(SITE = Site, ANPP = LITFALL_sum*2, TSOI = TSOI_mean, CLAY = PCT_CLAY_mean, GWC = H2OSOI_mean*100, W_SCALAR=W_SCALAR_mean) %>%
  select(SITE, ANPP, TSOI, CLAY, LIG_N, LIG_N_sp1, LIG_N_sp2, LIG_N_sp3, GWC, W_SCALAR, lci_SM_ratio, uci_SM_ratio) 
#fixing TALL ANPP
NEON_GPP <- read.csv("Example_simulations/Data/NEON_GPP.csv")
DailyInput <- read.csv("Example_simulations/Data/DailyInput.csv") %>% select(-MAT)
DailyInput$LITFALL[DailyInput$SITE == "TALL"] <- DailyInput$LITFALL[DailyInput$SITE == "TALL"]*0.60
DailyInput$ANPP[DailyInput$SITE == "TALL"] <- sum(DailyInput$LITFALL[DailyInput$SITE == "TALL"])
#replacing MSBio data with daily input sums and means to ensure comparable data between daily data and annual data
DI_sum <- DailyInput %>% select(SITE, ANPP, TSOI, W_SCALAR) %>% group_by(SITE, ANPP) %>% summarise(TSOI=mean(TSOI), W_SCALAR=mean(W_SCALAR))
MSBio3 <- MSBio2 %>% select(-ANPP, -TSOI, -W_SCALAR) %>% inner_join(DI_sum, by="SITE")
#filtering for only sites with microbial data to match observations
Mic_sites <- c("SERC","BART","TALL","TREE","LENO","HARV","GRSM")
MSBio_sites <- filter(MSBio3, SITE %in% Mic_sites)

#additional code for: (1) creating DailyInput file for all sites and (2) determining litterfall multiplier for TALL
#(1)
#daily data - change site name and MSBio2 row (e.g., 1=BART) to use different site daily input
# TREE_dailyinput <- read.csv("Example_simulations/Data/TREE_clim.csv")
#TREE_DI <- TREE_dailyinput %>% mutate(DAY=X, ANPP = rep(sum(LITFALL)*2,366), LITFALL=LITFALL*2, CLAY = rep(MSBio2[7,4], 366),
#                                       LIG_N = rep(MSBio2[7,5], 366), GWC = H2OSOI*100) %>% # , MAT=TBOT
#   select(DAY, ANPP, LITFALL, TSOI, CLAY, LIG_N, GWC, W_SCALAR) # MAT,
# DailyInput <- rbind(BART_DI, GRSM_DI, HARV_DI, LENO_DI, SERC_DI, TALL_DI, TREE_DI)
# DailyInput$SITE <- c(rep("BART", 366), rep("GRSM", 366), rep("HARV", 366), rep("LENO", 365), rep("SERC", 366), rep("TALL", 366), rep("TREE", 366))
# write.csv(DailyInput, "Example_simulations/Data/DailyInput.csv")

# TREE_comp <- data.frame(SITE="TREE", AGNPP_sum=sum(TREE_dailyinput$AGNPP), LF_sum=sum(TREE_dailyinput$LITFALL))
# LF_comp <- rbind(BART_comp, GRSM_comp, HARV_comp, LENO_comp, SERC_comp, TALL_comp, TREE_comp)
# write.csv(LF_comp, "LF_comp.csv")
#checking for NPP which should include above and belowground !
#BART_comp2 <- data.frame(SITE="BART", AGNPP_sum=sum(BART_dailyinput$NPP), LF_sum=sum(BART_dailyinput$LITFALL))
#LF_comp2 <- rbind(BART_comp2, GRSM_comp2, HARV_comp2, LENO_comp2, SERC_comp2, TALL_comp2, TREE_comp2)
#write.csv(LF_comp2, "LF_comp2.csv")
#LF_comp2 <- read.csv("LF_comp2.csv")
#(2)
# DI_40 <- DailyInput %>% filter(DAY==40)
# NPPComp <-DI_40%>% inner_join(NEON_GPP, by = "SITE")
# colorBlind7  <- c("#E69F00", "#56B4E9", "#009E73",
#                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 
# NPPComp <- NPPComp %>% mutate(SITE=factor(SITE, levels=c("TREE", "BART", "HARV", "GRSM", "SERC", "TALL", "LENO"))) #MAT order
# ggplot() + geom_point(data=NPPComp, aes(x=Annual.GPP, y= ANPP, color=SITE),size=4) +
#   geom_smooth(data=NPP_GoodSites, aes(x=Annual.GPP, y=ANPP), method='lm', color="black", se=FALSE)+ 
#   theme_bw(base_size = 16) + scale_color_manual(values=colorBlind7) + 
#   xlab(expression(paste("NEON GPP (gC m"^"-2"*" yr"^"-1"*")"))) +
#   ylab(expression(paste("CLM annual litterfall (gC m"^"-2"*" yr"^"-1"*")")))#TALL is too high and LENO is too low in comparison to NEON GPP
# NPP_GoodSites <- NPPComp %>% filter(SITE != "TALL")
# NPP_mod <- lm(ANPP~Annual.GPP, data=NPP_GoodSites)
# summary(NPP_mod) #510 + 0.18x #with derecho: 530 + 0.50x
# #determining multipliers for daily inputs -> new divided by old LITFALL sum estimates
# TALL_adj =  530 + 0.50*NEON_GPP[3,2]
# TALL_mult <- TALL_adj/DI_40[6,3] #0.60
# #checking if more aligned now
# NPPComp2 <-MSBio3 %>% inner_join(NEON_GPP, by = "SITE")
# ggplot(NPPComp2, aes(x=Annual.GPP, y= ANPP, color=SITE)) + geom_point(size=4) + theme_bw(base_size = 16)


#load in MSBio litter bag chemistry
### changed to fMET calculation in STODE script here!! Note that the two options are only somewhat related but less negatives in STODE equation
MSBio_BAGS <- MSBio_sites %>% select(SITE, LIG_N_sp1, LIG_N_sp2, LIG_N_sp3) %>% pivot_longer(2:4, names_to = "TYPE", values_to = "BAG_LIG_N")
MSBio_BAGS$CALC_MET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * (MSBio_BAGS$BAG_LIG_N))
MSBio_BAGS$CALC_MET[MSBio_BAGS$CALC_MET <0] = 0.01 #setting negatives to small number so 99% structural

BAG_init_size <- 100
BAGS <- MSBio_BAGS %>% select(SITE, TYPE, CALC_MET)
BAGS$BAG_LITm <- ((BAG_init_size * 1e3 / 1e4)/ depth) * BAGS$CALC_MET #g/m2 converted to mg/cm3
BAGS$BAG_LITs <- ((BAG_init_size * 1e3 / 1e4)/ depth) * (1-BAGS$CALC_MET) 


####
#run litterbag model 
####


#all sites and all litters with daily input looping through different soil moistures and litters as well
#prep the data
#create different soil moisture for steady state and daily input data
MSBio_sites_SM <- rbind(MSBio_sites, MSBio_sites, MSBio_sites)
MSBio_sites_SM <- MSBio_sites_SM %>% mutate(SM_type = c(rep("mean", 7), rep("max", 7), rep("min", 7))) %>% 
  mutate(W_SCALAR2 = case_when(SM_type == "mean" ~ W_SCALAR,
         SM_type == "max" ~ W_SCALAR*uci_SM_ratio,
         SM_type == "min" ~ W_SCALAR*lci_SM_ratio)) %>%
  mutate(W_SCALAR2 = case_when(W_SCALAR2>1~1, TRUE ~ W_SCALAR2)) %>%
  mutate(W_SCALAR = W_SCALAR2)
DailyInput_SM <- rbind(DailyInput, DailyInput, DailyInput)
SM_mult <- MSBio_sites %>% select(SITE, uci_SM_ratio, lci_SM_ratio)
DailyInput_SM <- DailyInput_SM %>% left_join(SM_mult, by="SITE") %>% mutate(SM_type = c(rep("mean", 2561), rep("max", 2561), rep("min", 2561))) %>% 
  mutate(W_SCALAR2 = case_when(SM_type == "mean" ~ W_SCALAR,
         SM_type == "max" ~ W_SCALAR *uci_SM_ratio,
         SM_type == "min" ~ W_SCALAR *lci_SM_ratio)) %>%
  mutate(W_SCALAR2 = case_when(W_SCALAR2>1~1, TRUE ~ W_SCALAR2)) %>%
  mutate(W_SCALAR = W_SCALAR2)


#just one parameter set (for default model)
BAGS_out_AllSites_SP = data.frame()
SM = c("mean", "max", "min")
for (SM_type2 in SM) {
  MSBio_sites_in <- filter(MSBio_sites_SM, SM_type==SM_type2)
  DailyInput_in <- filter(DailyInput_SM, SM_type==SM_type2)
  LQ = c("LIG_N_sp1", "LIG_N_sp2", "LIG_N_sp3")
  for (bag_type in LQ) {
    BAGS_mean <- filter(BAGS, TYPE==bag_type)
    MSBio_sites_in$LIG_N = MSBio_sites_in[[bag_type]]
    for (site in Mic_sites) {
      BAGS_input <- filter(BAGS_mean, SITE == site)
      forcing_input <- filter(MSBio_sites_in, SITE == site)
      daily_input <- filter(DailyInput_in, SITE == site)
      BO_DI <- MIMICS_LITBAG(forcing_df = forcing_input, litBAG = BAGS_input, dailyInput = daily_input, nspin_yrs=3, nspin_days=0, litadd_day=315, verbose=T)
      BAGS_out_AllSites_SP <- rbind(BAGS_out_AllSites_SP, BO_DI)
    }
  }
}


#multiple parameter sets (for calibrated model)
#load calibrated parameter sets
ES_Psets <- read.csv("ES_Psets_5000_NewInputs_ES.csv")
ES_Psets <- ES_Psets %>% filter(ID==175 | ID==176 | ID==190) 
BAGS_out_AllSites_Cal = data.frame()
Pset_ID <- ES_Psets$ID
for (i in Pset_ID) {
  ES_Pset_ID <- filter(ES_Psets, ID == i)
  print(i) #tracking pset
   tau_r <<- c(tau_r_default[1], tau_r_default[2] * ES_Pset_ID$Tau_r[1])
   beta <- beta_default * ES_Pset_ID$beta_x[1]
   vMOD <<- c(vMOD_default[1] * ES_Pset_ID$vMOD_m[1], vMOD_default[2] * ES_Pset_ID$vMOD_s[1], vMOD_default[3], vMOD_default[4] * ES_Pset_ID$vMOD_m[1], vMOD_default[5] * ES_Pset_ID$vMOD_s[1], vMOD_default[6])
  SM = c("mean", "max", "min")
  for (SM_type2 in SM) {
    MSBio_sites_in <- filter(MSBio_sites_SM, SM_type==SM_type2)
    DailyInput_in <- filter(DailyInput_SM, SM_type==SM_type2)
    LQ = c("LIG_N_sp1", "LIG_N_sp2", "LIG_N_sp3")
    for (bag_type in LQ) {
      BAGS_mean <- filter(BAGS, TYPE==bag_type)
      MSBio_sites_in$LIG_N = MSBio_sites_in[[bag_type]]
      for (site in Mic_sites) {
        BAGS_input <- filter(BAGS_mean, SITE == site)
        forcing_input <- filter(MSBio_sites_in, SITE == site)
        daily_input <- filter(DailyInput_in, SITE == site)
        BO_DI <- MIMICS_LITBAG(forcing_df = forcing_input, litBAG = BAGS_input, dailyInput = daily_input, nspin_yrs=3, nspin_days=0, litadd_day=315, verbose=T)
        BAGS_out_AllSites_Cal <- rbind(BAGS_out_AllSites_Cal,BO_DI)
      }
    }
  }
}


  
####
#plot output
####

#color palette for paper
colorBlind7  <- c("#E69F00", "#56B4E9", "#009E73",
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

#Formating observational data for comparing to field litter mass loss
Field_LML <- read.csv("Example_simulations/Data/Litter_decomp_all.csv")
LML_sum2 <- Field_LML  %>% group_by(site, time.point) %>% drop_na(percent.loss.litter) %>% filter(percent.loss.litter > 0) %>%
  summarize(mean.ML = mean(percent.loss.litter*100),
  n = n(),
  sd = sd(percent.loss.litter*100),
  SE = sd/sqrt(n),
  lci.ML = mean.ML - qt(1 - ((1 - 0.95) / 2), n - 1) * SE,
  uci.ML = mean.ML + qt(1 - ((1 - 0.95) / 2), n - 1) * SE,
  lci.ML.99 = mean.ML - qt(1 - ((1 - 0.99) / 2), n - 1) * SE,
  uci.ML.99 = mean.ML + qt(1 - ((1 - 0.99) / 2), n - 1) * SE,
  min.ML = min(percent.loss.litter*100),
  max.ML = max(percent.loss.litter*100),
  doy = mean(days_elapsed)) %>% mutate(doy=round(doy, digits=0)) %>%
filter(site %in% Mic_sites)

#add IDs for calibrated output
BAGS_out_AllSites_Cal <- BAGS_out_AllSites_Cal %>% mutate(ID=as.factor(rep(c(175, 176, 190), each=68985)))

#plotting format for Figure 2
BAGS_out_plot.SP <- BAGS_out_AllSites_SP %>% mutate(SITE.LT = paste(SITE, Litter_Type, sep=".")) %>% mutate(LIT_PerLoss = ((0.1 - (LITBAGm+LITBAGs))/0.1)*100)
BAGS_out_plot.Cal <- BAGS_out_AllSites_Cal %>% mutate(SITE.LT = paste(SITE, Litter_Type, sep=".")) %>% mutate(LIT_PerLoss = ((0.1 - (LITBAGm+LITBAGs))/0.1)*100)
#summary data for visualization
BO_plot_sum.SP <- BAGS_out_plot.SP %>% group_by(SITE,DAY) %>% summarise(mean=mean(LIT_PerLoss), min=min(LIT_PerLoss), max=max(LIT_PerLoss))
BO_plot_sum.Cal <- BAGS_out_plot.Cal %>% group_by(SITE,DAY) %>% summarise(mean=mean(LIT_PerLoss), min=min(LIT_PerLoss), max=max(LIT_PerLoss))
#reordering sites
BO_plot_sum.SP <- BO_plot_sum.SP %>% mutate(SITE=factor(SITE, levels=c("TREE", "BART", "HARV", "GRSM", "SERC", "TALL", "LENO"))) #MAT order
BO_plot_sum.Cal <- BO_plot_sum.Cal %>% mutate(SITE=factor(SITE, levels=c("TREE", "BART", "HARV", "GRSM", "SERC", "TALL", "LENO"))) #MAT order
LML_sum2 <- LML_sum2 %>% mutate(site=factor(site, levels=c("TREE", "BART", "HARV", "GRSM", "SERC", "TALL", "LENO"))) #MAT order
#default parameters
ggplot() +
  geom_ribbon(data=BO_plot_sum.SP, aes(y=100-mean, x=DAY-315, ymin = 100-min, ymax=100-max, group=SITE, fill=SITE, color=SITE), alpha = 0.3, size=0.5) +
  geom_point(data=LML_sum2, aes(y=100-mean.ML, x=doy, group=site, color=site), size = 3) +
  geom_errorbar(data=LML_sum2, aes(y=100-mean.ML, x=doy, ymin = 100-lci.ML, ymax = 100-uci.ML, group=site, color=site), width=0,linewidth=1) +
  ylab("Litter Bag C Remaining (%)") +
  xlab("Day") +
  xlim(0, 780) +
  scale_color_manual(values = colorBlind7) +
  scale_fill_manual(values = colorBlind7) +
  theme_bw(base_size = 20) 
#calibrated parameters
ggplot() +
  geom_ribbon(data=BO_plot_sum.Cal, aes(y=100-mean, x=DAY-315, ymin = 100-min, ymax=100-max, group=SITE, fill=SITE, color=SITE), alpha = 0.3, size=0.5) +
  geom_point(data=LML_sum2, aes(y=100-mean.ML, x=doy, group=site, color=site), size = 3) +
  geom_errorbar(data=LML_sum2, aes(y=100-mean.ML, x=doy, ymin = 100-lci.ML, ymax = 100-uci.ML, group=site, color=site), width=0,linewidth=1) +
  ylab("Litter Bag C Remaining (%)") +
  xlab("Day") +
  xlim(0, 780) +
  scale_color_manual(values = colorBlind7) +
  scale_fill_manual(values = colorBlind7) +
  theme_bw(base_size = 20) 

#comparing default and calibrated model (Figure 3)
#regressing SP vs calibrated
BAGS_out_AllSites_ES176 <- filter(BAGS_out_AllSites_Cal, ID==176) %>% select(-ID) 
BAGS_out_AllSites_ES175 <- filter(BAGS_out_AllSites_Cal, ID==175) %>% select(-ID)
BAGS_out_AllSites_ES190 <- filter(BAGS_out_AllSites_Cal, ID==190) %>% select(-ID)
#formating default
BO_initial_SP <- BAGS_out_AllSites_SP %>% filter(DAY==315) %>% mutate(MICrK = MICr/MICk, LITms = LITm/LITs) %>% pivot_longer(5:19, names_to = 'Pools', values_to = 'Carbon_SP') %>% 
  mutate(ID=1:length(Carbon_SP)) %>% select(Carbon_SP, ID)
#formatting calibrated
BO_initial_CE.175 <- BAGS_out_AllSites_ES175 %>% filter(DAY==315) %>% mutate(MICrK = MICr/MICk, LITms = LITm/LITs) %>% pivot_longer(5:19, names_to = 'Pools', values_to = 'Carbon_CE') %>% mutate(ID=1:length(Carbon_CE))
BO_initial_CE.176 <- BAGS_out_AllSites_ES176 %>% filter(DAY==315) %>% mutate(MICrK = MICr/MICk, LITms = LITm/LITs) %>% pivot_longer(5:19, names_to = 'Pools', values_to = 'Carbon_CE') %>% mutate(ID=1:length(Carbon_CE))
BO_initial_CE.190 <- BAGS_out_AllSites_ES190 %>% filter(DAY==315) %>% mutate(MICrK = MICr/MICk, LITms = LITm/LITs) %>% pivot_longer(5:19, names_to = 'Pools', values_to = 'Carbon_CE') %>% mutate(ID=1:length(Carbon_CE))
BAGS2 <- BAGS %>% mutate(SITE.LQ = paste(SITE, TYPE, sep=".")) %>% select(SITE.LQ, CALC_MET)
#combining default and calibrated output
BO_initial175 <- inner_join(BO_initial_SP, BO_initial_CE.175, by='ID') %>% mutate(SITE.LQ = paste(SITE, Litter_Type, sep=".")) %>% inner_join(BAGS2, by="SITE.LQ")
BO_initial176 <- inner_join(BO_initial_SP, BO_initial_CE.176, by='ID') %>% mutate(SITE.LQ = paste(SITE, Litter_Type, sep=".")) %>% inner_join(BAGS2, by="SITE.LQ")
BO_initial190 <- inner_join(BO_initial_SP, BO_initial_CE.190, by='ID') %>% mutate(SITE.LQ = paste(SITE, Litter_Type, sep=".")) %>% inner_join(BAGS2, by="SITE.LQ")
#combining all into one dataframe
BO_initial_sum <- rbind(BO_initial175, BO_initial176, BO_initial190) %>% group_by(SITE, Litter_Type, CALC_MET, SM_Type, Pools) %>% 
  summarise(sp.avg =mean(Carbon_SP), cal.avg =mean(Carbon_CE), n=n(),
  lci.sp = sp.avg - qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_SP)/sqrt(n)), uci.sp = sp.avg + qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_SP)/sqrt(n)),
  lci.cal = cal.avg - qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_CE)/sqrt(n)), uci.cal = cal.avg + qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_CE)/sqrt(n)))
#metabolic litter
LITm_pool <- c("LITm")
BO_initial_sum %>% filter(Pools %in% LITm_pool) %>% ggplot() + geom_point(aes(x=sp.avg*1000, y=cal.avg*1000, color=CALC_MET, shape=SM_Type), alpha=0.5, size=3) + 
  ylab(expression(paste("Steady state C (g C m"^"-2"*")"))) + xlab(expression(paste("Steady state C (g C m"^"-2"*")"))) +
  facet_grid(.~Pools) + geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) + theme(legend.position="none")
#structural litter
LITs_pool <- c("LITs")
BO_initial_sum %>% filter(Pools %in% LITs_pool) %>% ggplot() + geom_point(aes(x=sp.avg*1000, y=cal.avg*1000, color=CALC_MET, shape=SM_Type), alpha=0.5, size=3) + 
  ylab(expression(paste("Steady state C (g C m"^"-2"*")"))) + xlab(expression(paste("Steady state C (g C m"^"-2"*")"))) +
  facet_grid(.~Pools) + geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) + theme(legend.position="none")
#copiotroph:oligotroph ratio
MICrK_pool <- c("MICrK")
BO_initial_sum %>% filter(Pools %in% MICrK_pool) %>% ggplot() + geom_point(aes(x=sp.avg, y=cal.avg, color=CALC_MET, shape=SM_Type), alpha=0.5, size=3) + 
  ylab("Calibrated model") + xlab("Default model") + xlab("Copiotroph:oligotroph") + ylab("Copiotroph:oligotroph") +
  facet_grid(.~Pools) + geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) + theme(legend.position="none")
#decomposition rates
LB.i <- BAGS_out_AllSites_SP %>% filter(DAY==315) %>% mutate(SITE.LQ = paste(SITE, Litter_Type, sep=".")) %>% group_by(SITE.LQ) %>% 
                                                               summarise(LBm.i = mean(LITBAGm), LBs.i = mean(LITBAGs))
#formatting default
BO_Decomp_SP2 <- BAGS_out_AllSites_SP %>% mutate(SITE.LQ = paste(SITE, Litter_Type, sep=".")) %>% inner_join(LB.i, by="SITE.LQ") %>% filter(DAY>315) %>% 
  group_by(SITE,Litter_Type, SM_Type) %>% summarise(Drm = mean(Decomp_rate_rm), Drs = mean(Decomp_rate_rs), 
            Dkm = mean(Decomp_rate_km), Dks = mean(Decomp_rate_ks)) %>% 
  pivot_longer(4:7, names_to = 'Pools', values_to = 'Carbon_SP') %>% 
  mutate(combo = paste(SITE, Litter_Type, SM_Type, Pools, sep=".")) %>% select(Carbon_SP, combo)
#formatting calibrated
BO_Decomp_CE2.175 <- BAGS_out_AllSites_ES175 %>% filter(DAY>315) %>% mutate(SITE.LQ = paste(SITE, Litter_Type, sep=".")) %>% inner_join(LB.i, by="SITE.LQ") %>% group_by(SITE, Litter_Type, SM_Type) %>%
  summarise(Drm = mean(Decomp_rate_rm), Drs = mean(Decomp_rate_rs), Dkm = mean(Decomp_rate_km), Dks = mean(Decomp_rate_ks)) %>% 
  pivot_longer(4:7, names_to = 'Pools', values_to = 'Carbon_CE') %>% mutate(combo = paste(SITE, Litter_Type, SM_Type, Pools, sep="."))
BO_Decomp_CE2.176 <- BAGS_out_AllSites_ES176 %>% filter(DAY>315) %>% mutate(SITE.LQ = paste(SITE, Litter_Type, sep=".")) %>% inner_join(LB.i, by="SITE.LQ") %>% group_by(SITE, Litter_Type, SM_Type) %>%
  summarise(Drm = mean(Decomp_rate_rm), Drs = mean(Decomp_rate_rs), Dkm = mean(Decomp_rate_km), Dks = mean(Decomp_rate_ks)) %>% 
  pivot_longer(4:7, names_to = 'Pools', values_to = 'Carbon_CE') %>% mutate(combo = paste(SITE, Litter_Type, SM_Type, Pools, sep="."))
BO_Decomp_CE2.190 <- BAGS_out_AllSites_ES190 %>% filter(DAY>315) %>% mutate(SITE.LQ = paste(SITE, Litter_Type, sep=".")) %>% inner_join(LB.i, by="SITE.LQ") %>% group_by(SITE, Litter_Type, SM_Type) %>%
  summarise(Drm = mean(Decomp_rate_rm), Drs = mean(Decomp_rate_rs), Dkm = mean(Decomp_rate_km), Dks = mean(Decomp_rate_ks)) %>% 
  pivot_longer(4:7, names_to = 'Pools', values_to = 'Carbon_CE') %>% mutate(combo = paste(SITE, Litter_Type, SM_Type, Pools, sep="."))
#combining default and calibrated output
BO_decomp175.2 <- inner_join(BO_Decomp_SP2, BO_Decomp_CE2.175, by='combo') %>% mutate(SITE.LQ = paste(SITE.x, Litter_Type.x, sep=".")) %>% inner_join(BAGS2, by="SITE.LQ")
BO_decomp176.2 <- inner_join(BO_Decomp_SP2, BO_Decomp_CE2.176, by='combo') %>% mutate(SITE.LQ = paste(SITE.x, Litter_Type.x, sep=".")) %>% inner_join(BAGS2, by="SITE.LQ")
BO_decomp190.2 <- inner_join(BO_Decomp_SP2, BO_Decomp_CE2.190, by='combo') %>% mutate(SITE.LQ = paste(SITE.x, Litter_Type.x, sep=".")) %>% inner_join(BAGS2, by="SITE.LQ")
#combining all into one dataframe
BO_decomp_sum2 <- rbind(BO_decomp175.2, BO_decomp176.2, BO_decomp190.2) %>% group_by(SITE.x, Litter_Type.x, CALC_MET, SM_Type, Pools) %>% 
  summarise(sp.avg =mean(Carbon_SP), cal.avg =mean(Carbon_CE), n=n(),
            lci.sp = sp.avg - qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_SP)/sqrt(n)), uci.sp = sp.avg + qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_SP)/sqrt(n)),
            lci.cal = cal.avg - qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_CE)/sqrt(n)), uci.cal = cal.avg + qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_CE)/sqrt(n)))
Decomp_pools <- c("Drs", "Dkm") 
BO_decomp_sum2 %>% filter(Pools %in% Decomp_pools) %>% ggplot(aes(x=sp.avg*24000, y=cal.avg*24000)) + geom_point(aes(color=CALC_MET, shape=SM_Type), alpha=0.5, size=3) + facet_grid(.~Pools) + 
  geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) + theme(legend.position="none") +
  xlab(expression(paste("Decomposition flux (g C m"^"-2"*" d"^"-1"*")"))) + ylab(expression(paste("Decomposition flux (g C m"^"-2"*" d"^"-1"*")"))) #filter(Pools=="Dkm") %>%


#formatting data for model vs obs plots and effect size plots (Figures 1 and 2)
FieldData <- LML_sum2 %>% mutate(DAY=doy, SITE=site) %>% mutate(SITE.DAY=paste(SITE, DAY, sep=".")) %>% select(time.point,SITE.DAY, mean.ML, lci.ML, uci.ML)
LITi = 0.1
#default parameters
df_LML.SP <- BAGS_out_AllSites_SP %>% mutate(DAY.LitOut = DAY -314) %>% mutate(SITE.DAY=paste(SITE, DAY.LitOut, sep=".")) %>% 
  right_join(FieldData, by="SITE.DAY") %>% mutate(LIT_PerLoss = ((LITi - (LITBAGm+LITBAGs))/LITi)*100)
#calibrated parameters
df_LML.Cal <- BAGS_out_AllSites_Cal %>% mutate(DAY.LitOut = DAY -314) %>% mutate(SITE.DAY=paste(SITE, DAY.LitOut, sep=".")) %>% 
  right_join(FieldData, by="SITE.DAY") %>% mutate(LIT_PerLoss = ((LITi - (LITBAGm+LITBAGs))/LITi)*100)

#comparison of modeled and observed decomp at timepoint 1 and 2
#summarized values for plotting
#default
df_LML_sum.SP <- df_LML.SP %>% group_by(SITE, time.point) %>% summarise(mean.ML = mean(mean.ML), m.uci.ML = mean(uci.ML), m.lci.ML = mean(lci.ML), mean.LPL = mean(LIT_PerLoss),
                                                                       n=n(), SE = sd(LIT_PerLoss)/sqrt(n),
                                                                       min.LPL = min(LIT_PerLoss),
                                                                       max.LPL = max(LIT_PerLoss),
                                                                       lci.LPL = mean.LPL - qt(1 - ((1 - 0.95) / 2), n - 1) * SE,
                                                                       uci.LPL = mean.LPL + qt(1 - ((1 - 0.95) / 2), n - 1) * SE)
#calibrated
df_LML_sum.Cal <- df_LML.Cal %>% group_by(SITE, time.point) %>% summarise(mean.ML = mean(mean.ML), m.uci.ML = mean(uci.ML), m.lci.ML = mean(lci.ML), mean.LPL = mean(LIT_PerLoss),
                                                                        n=n(), SE = sd(LIT_PerLoss)/sqrt(n),
                                                                        min.LPL = min(LIT_PerLoss),
                                                                        max.LPL = max(LIT_PerLoss),
                                                                        lci.LPL = mean.LPL - qt(1 - ((1 - 0.95) / 2), n - 1) * SE,
                                                                        uci.LPL = mean.LPL + qt(1 - ((1 - 0.95) / 2), n - 1) * SE)
#summary stats and plots - change out df_LML_sum for default and calibrated versions of these
modelVobs <- lm(mean.LPL~mean.ML, data=df_LML_sum)
summary(modelVobs) #R2
sqrt(mean((df_LML_sum$mean.ML - df_LML_sum$mean.LPL)^2)) #RMSE
(1/length(df_LML_sum$n))*sum(df_LML_sum$mean.ML - df_LML_sum$mean.LPL) #bias
#reordering sites
df_LML_sum <- df_LML_sum %>% mutate(SITE=factor(SITE, levels=c("TREE", "BART", "HARV", "GRSM", "SERC", "TALL", "LENO"))) #MAT order
#plotting obs vs model
ggplot(df_LML_sum, aes(x=mean.ML, y=mean.LPL)) + geom_point(aes(color=SITE), size=4) + geom_smooth(method = "lm", color="black")  + xlim(0,80) + ylim(0,80) +
  geom_errorbar(aes(ymin=lci.LPL, ymax=uci.LPL, color=SITE), size=1) + geom_errorbarh(aes(xmin=m.lci.ML, xmax=m.uci.ML, color=SITE), size=1) +
  xlab("Observed litter percent C loss") + ylab("Modeled litter percent C loss") + geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) + 
  scale_color_manual(values = colorBlind7) + theme(legend.position="none")

#effect size estimates
#formatting data for statistical model
DI_means <- DailyInput_SM %>% mutate(SITE.SM = paste(SITE, SM_type, sep = ".")) %>% group_by(SITE.SM) %>% 
  summarise(W_SCALAR_mean=mean(W_SCALAR)) %>% select(SITE.SM, W_SCALAR_mean)
#default
MIC_init.SP <- BAGS_out_AllSites_SP %>% filter(DAY == 315) %>%mutate(MICrK.i =  MICr/MICk)%>% mutate(MICr.i = MICr)%>% mutate(MICK.i = MICk)%>%
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% select(SITE.SM.LQ, MICrK.i, MICr.i, MICK.i)
#calibrated
MIC_init.Cal.175 <- BAGS_out_AllSites_Cal %>% filter(ID == 175)%>% filter(DAY == 315) %>%mutate(MICrK.i =  MICr/MICk)%>% mutate(MICr.i = MICr)%>% mutate(MICK.i = MICk)%>%
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% select(SITE.SM.LQ, MICrK.i, MICr.i, MICK.i)
MIC_init.Cal.176 <- BAGS_out_AllSites_Cal %>% filter(ID == 176)%>% filter(DAY == 315) %>%mutate(MICrK.i =  MICr/MICk)%>% mutate(MICr.i = MICr)%>% mutate(MICK.i = MICk)%>%
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% select(SITE.SM.LQ, MICrK.i, MICr.i, MICK.i)
MIC_init.Cal.190 <- BAGS_out_AllSites_Cal %>% filter(ID == 190)%>% filter(DAY == 315) %>%mutate(MICrK.i =  MICr/MICk)%>% mutate(MICr.i = MICr)%>% mutate(MICK.i = MICk)%>%
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% select(SITE.SM.LQ, MICrK.i, MICr.i, MICK.i)
BAGS_LIGN <- MSBio_BAGS %>% mutate(SITE.LQ = paste(SITE, TYPE, sep = ".")) %>% select(SITE.LQ, BAG_LIG_N)
#bringing data together for statistical model
df_analysis.SP <- df_LML.SP %>% mutate(MICrK = MICr/MICk) %>% mutate(MIC=MICr+MICk) %>% mutate(SOC = SOMa+SOMc+SOMp) %>% 
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% mutate(SITE.SM = paste(SITE, SM_Type, sep = ".")) %>%
  mutate(SITE.LQ = paste(SITE, Litter_Type, sep = ".")) %>% inner_join(DI_means, by="SITE.SM") %>% inner_join(MIC_init.SP, by="SITE.SM.LQ") %>% 
  inner_join(BAGS_LIGN, by="SITE.LQ") %>% mutate(LQ.SM=paste(Litter_Type, SM_Type, sep="."))
df_analysis.Cal.175 <- df_LML.Cal %>% filter(ID == 175)%>% mutate(MICrK = MICr/MICk) %>% mutate(MIC=MICr+MICk) %>% mutate(SOC = SOMa+SOMc+SOMp) %>% 
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% mutate(SITE.SM = paste(SITE, SM_Type, sep = ".")) %>%
  mutate(SITE.LQ = paste(SITE, Litter_Type, sep = ".")) %>% inner_join(DI_means, by="SITE.SM") %>% inner_join(MIC_init.Cal.175, by="SITE.SM.LQ") %>% 
  inner_join(BAGS_LIGN, by="SITE.LQ") %>% mutate(LQ.SM=paste(Litter_Type, SM_Type, sep=".")) 
df_analysis.Cal.176 <- df_LML.Cal %>% filter(ID == 176)%>% mutate(MICrK = MICr/MICk) %>% mutate(MIC=MICr+MICk) %>% mutate(SOC = SOMa+SOMc+SOMp) %>% 
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% mutate(SITE.SM = paste(SITE, SM_Type, sep = ".")) %>%
  mutate(SITE.LQ = paste(SITE, Litter_Type, sep = ".")) %>% inner_join(DI_means, by="SITE.SM") %>% inner_join(MIC_init.Cal.176, by="SITE.SM.LQ") %>% 
  inner_join(BAGS_LIGN, by="SITE.LQ") %>% mutate(LQ.SM=paste(Litter_Type, SM_Type, sep=".")) 
df_analysis.Cal.190 <- df_LML.Cal %>% filter(ID == 190)%>% mutate(MICrK = MICr/MICk) %>% mutate(MIC=MICr+MICk) %>% mutate(SOC = SOMa+SOMc+SOMp) %>% 
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% mutate(SITE.SM = paste(SITE, SM_Type, sep = ".")) %>%
  mutate(SITE.LQ = paste(SITE, Litter_Type, sep = ".")) %>% inner_join(DI_means, by="SITE.SM") %>% inner_join(MIC_init.Cal.190, by="SITE.SM.LQ") %>% 
  inner_join(BAGS_LIGN, by="SITE.LQ") %>% mutate(LQ.SM=paste(Litter_Type, SM_Type, sep=".")) 
#logical checks - replace "df_analysis" with the respective df_analyses above to build a statistical model for the default parameters and each of the calibrated parameter sets (175, 176, 190)
df_check <- df_analysis %>% filter(MICrK > 0.01) %>%
  filter(MICrK < 100) %>%
  filter(MIC/SOC > 0.0001) %>%
  filter(MIC/SOC < 0.40) 
#statistical model
Obs_ES_mod <- lmer(LIT_PerLoss ~ scale(log(W_SCALAR_mean))+scale(BAG_LIG_N)+scale(MICrK.i)+ (1|SITE/LQ.SM), data = df_check)
#checking collinearity of models
vif(Obs_ES_mod) # for SP: all less than 5, up to 3.6
#for calibrated: (175) all less than 5, up to 3.6; (176) all less than 5, up to 4.7; (190) all less than 5, up to 4.0
Obs_ES <- as.data.frame(fixef(Obs_ES_mod)) #fixed effects coefficients as effect size
Obs_ES$Vars <- rownames(Obs_ES)
colnames(Obs_ES)[1] <- "value"
Obs_ES <- Obs_ES[-1, ]
Obs_ES$mult <- ifelse(Obs_ES$value <0, -1, 1)
Obs_ES$rel_ES <- (abs(Obs_ES$value)/sum(abs(Obs_ES$value))) * 100 * Obs_ES$mult
#rename Ob_ES to match parameter set and lines 399 and 400
Obs_ES_SP<- Obs_ES 
Obs_ES_Cal <- rbind(Obs_ES176, Obs_ES175, Obs_ES190) %>% group_by(Vars) %>% summarise(mean.rES = mean(rel_ES), sd.rES = sd(rel_ES)) %>% mutate(ID="Cal")
Obs_ES_SP2 <- Obs_ES_SP %>% mutate(mean.rES=rel_ES, sd.rES=0, ID="SP") %>% select(Vars, mean.rES, sd.rES, ID)
#creating data frame for observed effect sizes
Vars= c('scale(log(W_SCALAR_mean))','scale(BAG_LIG_N)', 'scale(MICrK.i)')
mean.rES = c(42.8, -34.2, -22.9)
sd.rES = c(0,0,0)
ID= c("Obs", "Obs", "Obs")
Obs_ES_Obs <- data.frame(Vars, mean.rES, sd.rES, ID)
#plotting observed, default, and calibrated on one plot
Obs_ES_Cal_points <- rbind(Obs_ES176, Obs_ES175, Obs_ES190) %>% mutate(ID="Cal") %>% select(Vars, rel_ES, ID) 
Obs_ES_SP_points <- Obs_ES_SP2 %>% mutate(rel_ES = c(NA,NA,NA)) %>% select(Vars, rel_ES, ID)
Obs_ES_Obs_points <- Obs_ES_Obs %>% mutate(rel_ES = c(NA,NA,NA)) %>% select(Vars, rel_ES, ID)
ES_points <- rbind(Obs_ES_Cal_points, Obs_ES_SP_points, Obs_ES_Obs_points)
Obs_ES_all <- rbind(Obs_ES_Obs, Obs_ES_SP2, Obs_ES_Cal) 
ggplot() + 
  geom_bar(data=Obs_ES_all, aes(x=factor(Vars, level=c('scale(log(W_SCALAR_mean))','scale(BAG_LIG_N)', 'scale(MICrK.i)'), labels = c("Soil moisture", "Lignin:N", "Copiotroph:oligotroph")),
               y=mean.rES, fill = factor(ID, level=c('Obs', 'SP', 'Cal'), labels = c("Observations", "Default", "Calibrated"))), stat="identity", position = position_dodge(), color="black", linewidth=1) +
  geom_point(data=ES_points, aes(x=factor(Vars, level=c('scale(log(W_SCALAR_mean))','scale(BAG_LIG_N)', 'scale(MICrK.i)'), labels = c("Soil moisture", "Lignin:N", "Copiotroph:oligotroph")),
                                    y=rel_ES, fill=factor(ID, level=c('Obs', 'SP', 'Cal'), labels = c("Observations", "Default", "Calibrated"))), 
             position = position_jitterdodge(dodge.width = 0.9), alpha=0.5, size=2) + #position_dodge(0.9)
  xlab("") + ylab("Relative effect size") + geom_abline(intercept=0, slope=0, color="black", linewidth=0.6) +
  theme_bw(base_size = 16)  + scale_fill_manual(name="Type", values = c(Observations = "black", Default = "white", Calibrated="darkgrey")) 
