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


##########
# MSBio runs
#############


#-------------------------------
#Using MSBio data
#-------------------------------
####
#load MSBio site and litter data and format to code structure
####

#load site data
MSBio <- read.csv("Example_simulations/Data/Site_annual_clim_validation.csv")
#match input data structure
#AGNPP should be in grams Dry Weight (gDW) not gC! multiply by 2 here to remedy
#switching AGNPP to LITFALL to match daily inputs!
MSBio2 <- MSBio %>% mutate(SITE = Site, ANPP = LITFALL_sum*2, TSOI = TSOI_mean, CLAY = PCT_CLAY_mean, GWC = H2OSOI_mean*100, W_SCALAR=W_SCALAR_mean) %>%
  select(SITE, ANPP, TSOI, CLAY, LIG_N, LIG_N_sp1, LIG_N_sp2, LIG_N_sp3, GWC, W_SCALAR, lci_SM_ratio, uci_SM_ratio) 
MSBio_sites <- MSBio2
#loading daily inputs 
DailyInput <- read.csv("Example_simulations/Data/DailyInput_validation.csv")
#replacing MSBio data with daily input sums and means to ensure comparable data between daily data and annual data
DI_sum <- DailyInput %>% select(SITE, ANPP, TSOI, W_SCALAR) %>% group_by(SITE, ANPP) %>% summarise(TSOI=mean(TSOI), W_SCALAR=mean(W_SCALAR))
MSBio3 <- MSBio2 %>% select(-ANPP, -TSOI, -W_SCALAR) %>% inner_join(DI_sum, by="SITE")


#additional code for creating daily input file
#daily data - change site name and MSBio2 row (1=MLBS, 2=SCBI, 3=UNDE) to use different site daily input
# UNDE_dailyinput <- read.csv("Example_simulations/Data/UNDE_clim.csv")
# UNDE_DI <- UNDE_dailyinput %>% mutate(DAY=X, ANPP = rep(sum(LITFALL)*2,366), LITFALL=LITFALL*2, CLAY = rep(MSBio2[3,4], 366),
#                                        LIG_N = rep(MSBio2[3,5], 366), GWC = H2OSOI*100, MAT=TBOT) %>%
#   select(DAY, ANPP, LITFALL, TSOI, MAT, CLAY, LIG_N, GWC, W_SCALAR)
# DailyInput <- rbind(MLBS_DI, SCBI_DI, UNDE_DI)
# DailyInput$SITE <- c(rep("MLBS", 366), rep("SCBI", 366), rep("UNDE", 366))
# write.csv(DailyInput, "Example_simulations/Data/DailyInput_validation.csv")



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
#below creates water scalar over 1 so maybe need to change all maxes where W_SCALAR over 1 is equal to 1? Mathematically, fine to go over 1....
MSBio_sites_SM <- MSBio_sites_SM %>% mutate(SM_type = c(rep("mean", 3), rep("max", 3), rep("min", 3))) %>% 
  mutate(W_SCALAR2 = case_when(SM_type == "mean" ~ W_SCALAR,
         SM_type == "max" ~ W_SCALAR*uci_SM_ratio,
         SM_type == "min" ~ W_SCALAR*lci_SM_ratio)) %>%
  mutate(W_SCALAR2 = case_when(W_SCALAR2>1~1, TRUE ~ W_SCALAR2)) %>%
  mutate(W_SCALAR = W_SCALAR2)
DailyInput_SM <- rbind(DailyInput, DailyInput, DailyInput)
SM_mult <- MSBio_sites %>% select(SITE, uci_SM_ratio, lci_SM_ratio)
DailyInput_SM <- DailyInput_SM %>% left_join(SM_mult, by="SITE") %>% mutate(SM_type = c(rep("mean", 1098), rep("max", 1098), rep("min", 1098))) %>% 
  mutate(W_SCALAR2 = case_when(SM_type == "mean" ~ W_SCALAR,
         SM_type == "max" ~ W_SCALAR *uci_SM_ratio,
         SM_type == "min" ~ W_SCALAR *lci_SM_ratio)) %>%
  mutate(W_SCALAR2 = case_when(W_SCALAR2>1~1, TRUE ~ W_SCALAR2)) %>%
  mutate(W_SCALAR = W_SCALAR2)


#just one parameter set
site_val <- c("MLBS", "SCBI", "UNDE")
BAGS_out_AllSites_SP = data.frame()
SM = c("mean", "max", "min")
for (SM_type2 in SM) {
  MSBio_sites_in <- filter(MSBio_sites_SM, SM_type==SM_type2)
  DailyInput_in <- filter(DailyInput_SM, SM_type==SM_type2)
  LQ = c("LIG_N_sp1", "LIG_N_sp2", "LIG_N_sp3")
  for (bag_type in LQ) {
    BAGS_mean <- filter(BAGS, TYPE==bag_type)
    MSBio_sites_in$LIG_N = MSBio_sites_in[[bag_type]]
    for (site in site_val) {
      BAGS_input <- filter(BAGS_mean, SITE == site)
      forcing_input <- filter(MSBio_sites_in, SITE == site)
      daily_input <- filter(DailyInput_in, SITE == site)
      BO_DI <- MIMICS_LITBAG(forcing_df = forcing_input, litBAG = BAGS_input, dailyInput = daily_input, nspin_yrs=3, nspin_days=0, litadd_day=315, verbose=T)
      BAGS_out_AllSites_SP <- rbind(BAGS_out_AllSites_SP, BO_DI)
    }
  }
}

#all litters and soil moistures
 site_val <- c("MLBS", "SCBI", "UNDE")
 ES_Psets <- read.csv("ES_Psets_5000_NewInputs_ES.csv")
 ES_Psets <- ES_Psets %>% filter(ID==175 | ID==176 | ID==190)
 BAGS_out_AllSites_ES_val = data.frame()
Pset_ID <- ES_Psets$ID
for (i in Pset_ID) {
  ES_Pset_ID <- filter(ES_Psets, ID == i)
  print(i) #tracking pset
   tau_r <<- c(tau_r_default[1], tau_r_default[2] * ES_Pset_ID$Tau_r[1])
   #tau_K <<- c(tau_K_default[1], tau_K_default[2] * ES_Pset_ID$Tau_K[1])
   #CUE <<- CUE_default * ES_Pset_ID$CUE_x[1]
   #vMOD <<- vMOD_default * ES_Pset_ID$vMOD_x[1]
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
      for (site in site_val) {
        BAGS_input <- filter(BAGS_mean, SITE == site)
        forcing_input <- filter(MSBio_sites_in, SITE == site)
        daily_input <- filter(DailyInput_in, SITE == site)
        BO_DI <- MIMICS_LITBAG(forcing_df = forcing_input, litBAG = BAGS_input, dailyInput = daily_input, nspin_yrs=3, nspin_days=0, litadd_day=315, verbose=T)
        BAGS_out_AllSites_ES_val <- rbind(BAGS_out_AllSites_ES_val,BO_DI)
      }
    }
  }
}





  
####
#plot output
####

colorBlind7  <- c("#E69F00", "#56B4E9", "#009E73",
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #yellow (LENO), blue (SERC), green (UNDE)

#Formating observational data for comparing to field litter mass loss
Field_LML <- read.csv("Example_simulations/Data/Litter_decomp_all.csv")
#Add Species to group_by to get species-specific summary
#note that because this is percent loss you do not need to convert to C!
#below doy is days elapsed!!
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
filter(site %in% site_val)


#comparison of baseline MIMICS to LML
#got rid of LITi calcs since all litters back to the same starting value
#LIT_init <- BAGS_out_AllSites_DI %>% filter(DAY == 10) %>% mutate(LITi = LITBAGm+LITBAGs) %>% 
#  mutate(SITE.LT = paste(SITE, Litter_Type, sep=".")) %>% select(SITE.LT, LITi)
#boxplot(LIT_init$LITi)
#BAGS_out_AllSites_ES_vMOD_noTk$ID <- as.factor(rep(1:12, each=68985)) #not lining up with numerically with what I would expect - not sure what's going on.... ES one worked so even weirder!
# #try counting each site and see if each site is repeated the same number of times!
BAGS_out_AllSites_ES_val$ID <- as.factor(rep(1:3, each=29565))
# BAGS_out_AllSites_test$ID <- as.factor(1:18)
# ID_test <- filter(BAGS_out_AllSites_test, DAY==50)
BAGS_out_plot <- BAGS_out_AllSites_SP %>% mutate(SITE.LT = paste(SITE, Litter_Type, sep=".")) %>% mutate(LIT_PerLoss = ((0.1 - (LITBAGm+LITBAGs))/0.1)*100)
#wide format for plotting
BAGS_out_wide = BAGS_out_plot %>% select(SITE, ID, Litter_Type, SM_Type, DAY, LIT_PerLoss) %>% #ID, 
  pivot_wider(names_from = Litter_Type, values_from = LIT_PerLoss)
#plotting - check!
#BAGS_out_wide_BART <- filter(BAGS_out_wide, SITE=='LENO')
ggplot() +
  geom_line(data=BAGS_out_wide, aes(y=100-LIG_N_sp1, x=DAY-315, group=SITE, color=SITE), linewidth=0.5, alpha=0.5) +
  geom_line(data=BAGS_out_wide, aes(y=100-LIG_N_sp2, x=DAY-315, group=SITE, color=SITE), linewidth=0.5, alpha=0.5) +
  geom_line(data=BAGS_out_wide, aes(y=100-LIG_N_sp3, x=DAY-315, group=SITE, color=SITE), linewidth=0.5, alpha=0.5) +
  #geom_ribbon(data=BAGS_out_wide, aes(y=100-LIG_N, x=DAY, ymin = 100-LIG_N_min, ymax=100-LIG_N_max, group=SITE, fill=SITE), alpha = 0.3) +
  geom_point(data=LML_sum2, aes(y=100-mean.ML, x=doy, group=site, color=site), size = 3) +
  geom_errorbar(data=LML_sum2, aes(y=100-mean.ML, x=doy, ymin = 100-lci.ML, ymax = 100-uci.ML, group=site, color=site), width=0,linewidth=1) +
  xlim(0, 780) +
  ylab("Litter Bag C Remaining (%)") +
  xlab("Days elapsed") +
  facet_wrap(.~SM_Type) +
  theme_bw(base_size = 20)
#summary data - works ok for visualization!
val_colors = c("#882255", "#999933", "#332288")
BO_plot_sum <- BAGS_out_plot  %>% group_by(SITE,DAY) %>% summarise(mean=mean(LIT_PerLoss), min=min(LIT_PerLoss), max=max(LIT_PerLoss)) #,SM_Type
tiff("MSBio_Fig3_SP.tiff", units="px", width=2000, height=1500, res=300)
ggplot() +
  geom_ribbon(data=BO_plot_sum, aes(y=100-mean, x=DAY-315, ymin = 100-min, ymax=100-max, group=SITE, fill=SITE, color=SITE), alpha = 0.3, size=0.7) +
  #geom_line(data=BO_plot_sum, aes(y=100-mean, x=DAY-315, group=SITE, color=SITE), linewidth=2, alpha = 0.3) +
  geom_point(data=LML_sum2, aes(y=100-mean.ML, x=doy, group=site, color=site), size = 3) +
  geom_errorbar(data=LML_sum2, aes(y=100-mean.ML, x=doy, ymin = 100-lci.ML, ymax = 100-uci.ML, group=site, color=site), width=0,linewidth=1) +
  ylab("Litter Bag C Remaining (%)") +
  xlab("Day") +
  xlim(0, 780) +
  theme_bw(base_size = 20) +
  scale_color_manual(values=val_colors) + scale_fill_manual(values=val_colors) #+
  #facet_wrap(.~SM_Type)
dev.off()
#seperate plotting of LITm and LITs
LIT_init <- BAGS_out_AllSites_SP %>% filter(DAY == 315) %>% mutate(LITm.i = LITBAGm) %>% mutate(LITs.i = LITBAGs) %>% 
  mutate(SITE.LT.SM = paste(SITE, Litter_Type, SM_Type,ID, sep=".")) %>% select(SITE.LT.SM, LITm.i, LITs.i)
BAGS_out_plot <- BAGS_out_AllSites_SP %>% mutate(SITE.LT.SM = paste(SITE, Litter_Type, SM_Type,ID, sep=".")) %>% inner_join(LIT_init, by="SITE.LT.SM") %>%
  mutate(LITm_PerLoss = ((LITm.i - (LITBAGm))/LITm.i)*100) %>% mutate(LITs_PerLoss = ((LITs.i - (LITBAGs))/LITs.i)*100) #%>%
  #filter(ID==20)
#plotting
BO_plot_sum.ms <- BAGS_out_plot %>% group_by(SITE,DAY) %>% summarise(mean.m=mean(LITm_PerLoss), min.m=min(LITm_PerLoss), max.m=max(LITm_PerLoss), 
                                                                             mean.s=mean(LITs_PerLoss), min.s=min(LITs_PerLoss), max.s=max(LITs_PerLoss)) %>%
  pivot_longer(3:8, names_to = 'LIT.ms', values_to = 'value') %>% ungroup() %>% mutate(LIT_type = ifelse(grepl("s", LIT.ms), "S", "M")) %>% 
  mutate(stat_type = rep(c("mean", "min", "max"), 15330)) %>% select(SITE, DAY, LIT_type, stat_type, value) %>% pivot_wider(names_from = stat_type, values_from = value)
#pivot wider based on new column that makes a structural column if LIT.ms contains "s"
ggplot() +
  geom_ribbon(data=BO_plot_sum.ms, aes(y=100-mean, x=DAY-315, ymin = 100-min, ymax=100-max, group=SITE, fill=SITE), alpha = 0.3) +
  #geom_line(data=BO_plot_sum, aes(y=100-mean, x=DAY-315, group=SITE, color=SITE), linewidth=2, alpha = 0.3) +
  geom_point(data=LML_sum2, aes(y=100-mean.ML, x=doy, group=site, color=site), size = 3) +
  geom_errorbar(data=LML_sum2, aes(y=100-mean.ML, x=doy, ymin = 100-lci.ML, ymax = 100-uci.ML, group=site, color=site), width=0,linewidth=1) +
  ylab("Litter Bag C Remaining (%)") +
  xlab("Day") +
  xlim(0, 780) +
  facet_wrap(.~LIT_type) +
  theme_bw(base_size = 20)


# #formating data for modelVobs, RWA, and effect size
#need to standardize and log this data before running it thru again
#MFG_stdzd <- MFG_analysis %>% mutate_at(c('C_O', 'LIG_N', 'log.vwc'), ~(scale(.) %>% as.vector))
FieldData <- LML_sum2 %>% mutate(DAY=doy, SITE=site) %>% mutate(SITE.DAY=paste(SITE, DAY, sep=".")) %>% select(time.point,SITE.DAY, mean.ML, lci.ML, uci.ML)
#LIT_init <- BAGS_out_AllSites_DI %>% filter(DAY == 10) %>% mutate(LITi = LITBAGm+LITBAGs) %>% select(SITE, LITi)
#boxplot(LIT_init$LITi)
#df <- BAGS_out_AllSites_DI %>% left_join(LIT_init, by = "SITE")
#BAGS_out_AllSites_ES_noTk$ID <- as.factor(rep(1:12, each=68985))
LITi = 0.1
df_LML <- BAGS_out_AllSites_ES_val %>% mutate(DAY.LitOut = DAY -314) %>% mutate(SITE.DAY=paste(SITE, DAY.LitOut, sep=".")) %>% 
  right_join(FieldData, by="SITE.DAY") %>% mutate(LIT_PerLoss = ((LITi - (LITBAGm+LITBAGs))/LITi)*100)
#write.csv#write.csv#write.csv(df_LML, 'df_LML_ES_vMOD_noTk_ParamUncert.csv')
#creating table with LML for each day
# df_LML_All <- BAGS_out_AllSites_var %>% filter(DAY>315) %>% mutate(LIT_PerLoss = ((LITi - (LITBAGm+LITBAGs))/LITi)*100)
# #filtering so each site only has days for which we have observations 
# T2_DAY <- df_LML %>% filter(time.point==2) %>% mutate(T2.DAY = DAY) %>% select(SITE, T2.DAY)
# T1_DAY <- df_LML %>% filter(time.point==1) %>% mutate(T1.DAY = DAY) %>% select(SITE, T1.DAY)
# df_LML_RI <- df_LML_All %>% inner_join(T2_DAY, by = "SITE") %>% inner_join(T1_DAY, by = "SITE") %>% filter(DAY <= T2.DAY+315) #accounts for 10 days of spinup before starting simulation
# df_LML_T1T2 <- df_LML_RI %>% filter(DAY==T2.DAY | DAY==T1.DAY)

#df_LML <- read.csv('df_LML_ES_noTk_ParamUncert.csv')

#comparison of modeled and observed decomp at timepoint 1 and 2
#site observed mean mass loss versus individual obs of modeled data
modelVobs <- lm(LIT_PerLoss~mean.ML, data=df_LML)
summary(modelVobs)
#RMSE
sqrt(mean((df_LML$mean.ML - df_LML$LIT_PerLoss)^2))
#bias
(1/length(df_LML$mean.ML))*sum(df_LML$mean.ML - df_LML$LIT_PerLoss)
#loop for muliple Psets
ID <- unique(df_LML$ID)
ES_fit <- data.frame()
for (i in ID) {
  df_LML_loop <- filter(df_LML, ID == i)
  df_LML_sum <- df_LML_loop %>% group_by(SITE, time.point) %>% summarise(mean.ML = mean(mean.ML), m.uci.ML = mean(uci.ML), m.lci.ML = mean(lci.ML), mean.LPL = mean(LIT_PerLoss),
                                                        n=n(), SE = sd(LIT_PerLoss)/sqrt(n),
                                                        min.LPL = min(LIT_PerLoss),
                                                        max.LPL = max(LIT_PerLoss),
                                                        lci.LPL = mean.LPL - qt(1 - ((1 - 0.95) / 2), n - 1) * SE,
                                                        uci.LPL = mean.LPL + qt(1 - ((1 - 0.95) / 2), n - 1) * SE)
  #site observed mean mass loss versus modeled mean mass loss
  modelVobs <- lm(mean.LPL~mean.ML, data=df_LML_sum)
  sum <- summary(modelVobs)
  aR2 <- sum$adj.r.squared
  #RMSE
  RMSE <- sqrt(mean((df_LML_sum$mean.ML - df_LML_sum$mean.LPL)^2))
  fit <- data.frame(i, aR2, RMSE)
  ES_fit <- rbind(ES_fit, fit)
}
df_LML_best <- filter(df_LML, ID == 4)
df_LML_sum <- df_LML %>% group_by(SITE, time.point) %>% summarise(mean.ML = mean(mean.ML), m.uci.ML = mean(uci.ML), m.lci.ML = mean(lci.ML), mean.LPL = mean(LIT_PerLoss),
                                                                       n=n(), SE = sd(LIT_PerLoss)/sqrt(n),
                                                                       min.LPL = min(LIT_PerLoss),
                                                                       max.LPL = max(LIT_PerLoss),
                                                                       lci.LPL = mean.LPL - qt(1 - ((1 - 0.95) / 2), n - 1) * SE,
                                                                       uci.LPL = mean.LPL + qt(1 - ((1 - 0.95) / 2), n - 1) * SE)
modelVobs <- lm(mean.LPL~mean.ML, data=df_LML_sum)
summary(modelVobs) #R2
sqrt(mean((df_LML_sum$mean.ML - df_LML_sum$mean.LPL)^2)) #RMSE
(1/length(df_LML_sum$n))*sum(df_LML_sum$mean.ML - df_LML_sum$mean.LPL) #bias
val_colors = c("#882255", "#999933", "#332288")
tiff("MSBio_Fig3_Cal_inset.tiff", units="px", width=1100, height=1030, res=300)
ggplot(df_LML_sum, aes(x=mean.ML, y=mean.LPL)) + geom_point(aes(color=SITE), size=4) + geom_smooth(method = "lm", color="black")  + xlim(0,80) + ylim(0,80) +
  geom_errorbar(aes(ymin=lci.LPL, ymax=uci.LPL, color=SITE), size=1) + geom_errorbarh(aes(xmin=m.lci.ML, xmax=m.uci.ML, color=SITE), size=1) +
  xlab("Observed litter percent C loss") + ylab("Modeled litter percent C loss") + geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) +
  scale_color_manual(values = val_colors) + theme(legend.position="none")
dev.off()

####
#RWA and effect size estimation - set up for multiple soil moistures and multiple litter qualities
####
#OPTION 1: preparing data for analysis over time
# DI_2y <- rbind(DailyInput, DailyInput)
# DI_2y$DAY2 <- c(0:365, 0:365, 0:365, 0:364, 0:365, 0:365, 0:365, 366:731, 366:731, 366:731, 364:728, 366:731, 366:731, 366:731) 
# DI_analysis <- DI_2y %>% inner_join(T2_DAY, by = "SITE") %>% filter(DAY2 >10 & DAY2 <= T2.DAY) %>% mutate(SITE.DAY = paste(SITE, DAY2, sep = "."))
# df_analysis <- df_LML_RI %>% mutate(MICrK = MICr/MICk) %>% mutate(MIC=MICr+MICk) %>% mutate(SOC = SOMa+SOMc+SOMp) %>% 
#   mutate(SITE.DAY = paste(SITE, DAY, sep = ".")) %>% inner_join(DI_analysis, by="SITE.DAY") #%>% filter(time.point==2)
#OPTION 2: mean daily input model output analysis 
DI_means <- DailyInput_SM %>% mutate(SITE.SM = paste(SITE, SM_type, sep = ".")) %>% group_by(SITE.SM) %>% 
  summarise(W_SCALAR_mean=mean(W_SCALAR), MAT_mean=mean(MAT)) %>% select(SITE.SM, W_SCALAR_mean, MAT_mean)
#initial MICrK
#BAGS_out_AllSites_ES$ID <- rep(1:18, each = 22995)
MIC_init <- BAGS_out_AllSites_ES_val %>% filter(ID == 1) %>% filter(DAY == 315) %>%mutate(MICrK.i =  MICr/MICk)%>%
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% select(SITE.SM.LQ, MICrK.i) #%>% filter(ID == 4)
#need bag means!
BAGS_LIGN <- MSBio_BAGS %>% mutate(SITE.LQ = paste(SITE, TYPE, sep = ".")) %>% select(SITE.LQ, BAG_LIG_N)
df_analysis <- df_LML %>% filter(ID == 1) %>% mutate(MICrK = MICr/MICk) %>% mutate(MIC=MICr+MICk) %>% mutate(SOC = SOMa+SOMc+SOMp) %>% 
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% mutate(SITE.SM = paste(SITE, SM_Type, sep = ".")) %>%
  mutate(SITE.LQ = paste(SITE, Litter_Type, sep = ".")) %>% inner_join(DI_means, by="SITE.SM") %>% 
  inner_join(MIC_init, by="SITE.SM.LQ") %>% inner_join(BAGS_LIGN, by="SITE.LQ") # %>% filter(ID == 3)
#logical checks
df_check <- df_analysis %>% filter(MICrK > 0.01) %>%
  filter(MICrK < 100) %>%
  filter(MIC/SOC > 0.0001) %>%
  filter(MIC/SOC < 0.40) 
#plot((df_check$MICr/df_check$MICk), df_check$BAG_LIG_N)
#relative weights analysis
df_check$log_WS <- log(df_check$W_SCALAR_mean)
# df_rwa <- df_check %>% select(LIT_PerLoss, log_WS, BAG_LIG_N, MICrK.i) %>% 
#   mutate_at(vars(c("log_WS", "BAG_LIG_N", "MICrK.i")), ~(scale(.) %>% as.vector))
# MIM_rwa <- rwa(df_rwa, "LIT_PerLoss", c("log_WS", "BAG_LIG_N", "MICrK.i"), applysigns = TRUE, plot = TRUE) #"MAT", 
# plot_rwa(MIM_rwa) 
#effect size
#df_ES <- df_check %>% mutate_at(vars(c("log_WS", "BAG_LIG_N", "MICrK.i")), ~(scale(.) %>% as.vector))
Obs_ES_mod <- lmer(LIT_PerLoss ~ scale(log(W_SCALAR_mean))+scale(BAG_LIG_N)+scale(MICrK.i)+ (1|SITE), data = df_check) #MAT #scale(BAG_LIG_N)+
#Obs_ES_mod2 <- lmer(LIT_PerLoss ~ log_WS+BAG_LIG_N+MICrK.i+ (1|SITE), data = df_ES)
summary(Obs_ES_mod)
Obs_ES <- as.data.frame(fixef(Obs_ES_mod)) #fixed effects coefficients as effect size
#Obs_ES <- as.data.frame(Obs_ES_mod$coefficients) #if not using fixed effects model do this
Obs_ES$Vars <- rownames(Obs_ES)
colnames(Obs_ES)[1] <- "value"
Obs_ES <- Obs_ES[-1, ]
Obs_ES$mult <- ifelse(Obs_ES$value <0, -1, 1)
Obs_ES$rel_ES <- (abs(Obs_ES$value)/sum(abs(Obs_ES$value))) * 100 * Obs_ES$mult
Obs_ES$Vars <- factor(Obs_ES$Vars, levels=c('scale(MICrK.i)','scale(BAG_LIG_N)', 'scale(log(W_SCALAR_mean))'),
                       labels=c('MICrK.i','BAG_LIG_N', 'log_WS'))
ggplot(Obs_ES, aes(x=Vars, y=rel_ES)) + geom_bar(stat="identity", fill="lightblue") + coord_flip() + 
  geom_text(aes(label=round(rel_ES, digits=1), vjust=1.5), size=5) +theme_bw(base_size = 16)
Obs_ES1<- Obs_ES
rbind(Obs_ES1, Obs_ES2, Obs_ES3) %>% group_by(Vars) %>% summarise(mean.rES = mean(rel_ES), sd.rES = sd(rel_ES)) %>%
  ggplot(aes(x=Vars, y=mean.rES)) + geom_bar(stat="identity", fill="lightblue", aes(group=Vars)) +  geom_errorbar(aes(ymax=mean.rES+sd.rES, ymin=mean.rES-sd.rES), width=0.1, size=1) + 
  coord_flip() + geom_text(aes(label=round(mean.rES, digits=1), vjust=1.5), size=5) +theme_bw(base_size = 16) #, hjust=1.1


#random forest
#test and training data - using 75 train-25 test split like in Georgiou et al., 2021
# MFG_rf <- as.data.frame(df_check %>% select(LIT_PerLoss, MAT, W_SCALAR, LIG_N, MICrK) %>% na.omit(.))
# split <- sample.split(MFG_rf, SplitRatio = 0.75)
# data_train <- subset(MFG_rf, split == "TRUE")
# data_test <- subset(MFG_rf, split == "FALSE")
# #run training data with ranger
# #MFG_analysis_NAo <- MFG_analysis %>% na.omit(.)
# ranger_train <- ranger(LIT_PerLoss ~ MAT + W_SCALAR + LIG_N + MICrK,
#                        data = data_train,
#                        importance = 'impurity',
#                        mtry = 1)
# print(ranger_train) #really similar to rF package
# pred_test2 <- predict(ranger_train, data = data_test)
# summary(r.rap <- data_test$LIT_PerLoss - pred_test2$predictions)
# (rmse.ra <- sqrt(sum(r.rap^2)/length(r.rap))) #11% - lower than rF
# plot(data_test$LIT_PerLoss ~ pred_test2$predictions, asp=1, pch=20, xlab="fitted", ylab="actual", main="Prediciton of Litter Decomposition, Ranger")
# grid(); abline(0,1)
# #comparing variable importance
# vip(ranger_train, title = "Ranger")
# #ensemble of random forests
# #function for getting random forest vips - removing train-test becasue using all data for this
# rf.vip <- function(df) {
#   MFG_rf <- as.data.frame(df %>% select(LIT_PerLoss, MAT, W_SCALAR, LIG_N, MICrK) %>% na.omit(.))
#   ranger_train <- ranger(LIT_PerLoss ~ MAT + W_SCALAR + LIG_N + MICrK,
#                        data = MFG_rf,
#                        importance = 'impurity',
#                        mtry = 1)
#   Obs_rf <- as.data.frame(ranger_train$variable.importance) #fixed effects coefficients as effect size
#   Obs_rf$Vars <- rownames(Obs_rf)
#   colnames(Obs_rf)[1] <- "value"
#   Obs_rf$rel_rf <- (abs(Obs_rf$value)/sum(abs(Obs_rf$value))) * 100
#   return(Obs_rf)
# }
# y = rf.vip(df_check)
# y %>% mutate(Vars=factor(Vars, levels = c("MICrK", "LIG_N", "MAT", "W_SCALAR"))) %>% ggplot(aes(x=Vars, y=rel_rf)) + geom_bar(stat="identity", fill="blue") + coord_flip() + geom_text(aes(label=round(rel_rf, digits=1)), color="red", size=7) +theme_bw(base_size = 16)
# #replicate function 10 times
# result <- t(replicate(10, rf.vip(df_check)))
# for (i in 1:10) {
#   rf_vip_out <- rbind(rf_vip_out, as.data.frame(result[i,]))
# }
# #assess variaiton in rf ensemble
# ggplot(rf_vip_out, aes(x=Vars, y=rel_rf, color=Vars)) + geom_boxplot() + theme_bw(base_size = 16)
# 
# 
# #pearson and spearman correlation
# df_cor <- df_check %>% select(LIT_PerLoss, MAT, W_SCALAR, LIG_N, MICrK)
# cor.p <- cor(df_cor, method="pearson")
# corrplot(cor.p, type = "upper") #aligned with RWA
# cor.s <- cor(df_cor, method="spearman")
# corrplot(cor.s, type = "lower")#aligned with RWA
# 
# #simple realtionships between variables
# ggplot(df_check, aes(x=W_SCALAR, y=LIT_PerLoss)) + geom_point(size=3, aes(color=SITE)) + geom_smooth(method="lm") + theme_bw(base_size = 16)
# ggplot(df_check, aes(x=MAT, y=LIT_PerLoss)) + geom_point(size=3, aes(color=SITE)) + geom_smooth(method="lm") + theme_bw(base_size = 16)
# ggplot(df_check, aes(x=LIG_N, y=LIT_PerLoss)) + geom_point(size=3, aes(color=SITE)) + geom_smooth(method="lm") + theme_bw(base_size = 16)
# ggplot(df_check, aes(x=MICrK.i, y=LIT_PerLoss)) + geom_point(size=3, aes(color=SITE)) + geom_smooth(method="lm") + theme_bw(base_size = 16)


#comparing spin up with mean vs individual litters
# BO_var <- filter(BAGS_out_AllSites_var, DAY == 315)
# BO_ss <- filter(BAGS_out_AllSites_ss, DAY == 315)
# ggplot(BO_var, aes(x=SITE, y=MICr/MICk, color=SITE)) +
#   geom_point(size=4) +
#   #geom_point(BAGS_out_AllSites_var, aes(x=DAY, y=MICr/MICk, color=SITE), shape=2) +
#   facet_grid(Litter_Type~SM_Type)
# ggplot(BO_var, aes(x=LITm/LITs, y=MICr/MICk, color=SITE)) +
#   geom_point(size=4) +
#   #geom_point(BAGS_out_AllSites_var, aes(x=DAY, y=MICr/MICk, color=SITE), shape=2) +
#   facet_grid(Litter_Type~SM_Type)
# ggplot(BO_ss, aes(x=LITm/LITs, y=MICr/MICk, color=SITE)) +
#   geom_point(size=4) +
#   #geom_point(BAGS_out_AllSites_var, aes(x=DAY, y=MICr/MICk, color=SITE), shape=2) +
#   facet_grid(Litter_Type~SM_Type)



###
#moisture function testing
###

#wide format MIMICS output for plotting
# BAGS_out_wide_fWmeth0 = BAGS_out %>% mutate(LITBAG = LITBAGm+LITBAGs) %>% select(SITE, Litter_Type, DAY, LITBAG) %>% pivot_wider(names_from = Litter_Type, values_from = LITBAG)
# BAGS_out_wide_fWmeth1 = BAGS_out_fWm1 %>% mutate(LITBAG = LITBAGm+LITBAGs) %>% select(SITE, Litter_Type, DAY, LITBAG) %>% pivot_wider(names_from = Litter_Type, values_from = LITBAG)
# BAGS_out_wide_fWmeth2 = BAGS_out_fWm2 %>% mutate(LITBAG = LITBAGm+LITBAGs) %>% select(SITE, Litter_Type, DAY, LITBAG) %>% pivot_wider(names_from = Litter_Type, values_from = LITBAG)
# BAGS_out_wide_fWmeth3 = BAGS_out_fWm3 %>% mutate(LITBAG = LITBAGm+LITBAGs) %>% select(SITE, Litter_Type, DAY, LITBAG) %>% pivot_wider(names_from = Litter_Type, values_from = LITBAG)
# 
# #wide format for fW effects on just Vmax or on both Vmax and tau
# BAGS_out_wide = BAGS_out_100y %>% mutate(LITBAG = LITBAGm+LITBAGs) %>% select(SITE, Litter_Type, DAY, LITBAG) %>% pivot_wider(names_from = Litter_Type, values_from = LITBAG)
# BAGS_out_wide_fW.tau = BAGS_out_fWt %>% mutate(LITBAG = LITBAGm+LITBAGs) %>% select(SITE, Litter_Type, DAY, LITBAG) %>% pivot_wider(names_from = Litter_Type, values_from = LITBAG)
# 
# #plot MIMICS output with different soil moisture and field litter mass loss together
# ggplot() +
#   #geom_line(data=BAGS_out_wide_fWmeth0, aes(y=(mean/0.1)*100, x=DAY, color ="fW=1"), linewidth=2, alpha=0.5) +
#   #geom_line(data=BAGS_out_wide_fWmeth1, aes(y=(mean/0.1)*100, x=DAY, color ="CORPSE"), linewidth=2, alpha=0.5) +
#   #geom_line(data=BAGS_out_wide_fWmeth2, aes(y=(mean/0.1)*100, x=DAY, color ="Calibrated"), linewidth=2, alpha=0.5) +
#   #geom_line(data=BAGS_out_wide_fWmeth3, aes(y=(mean/0.1)*100, x=DAY, color ="W_SCALAR"), linewidth=2, alpha=0.5) +
#   geom_line(data=BAGS_out_wide, aes(y=(mean/0.1)*100, x=DAY, color ="Default"), linewidth=2, alpha=0.5) +
#   #geom_line(data=BAGS_out_wide_fW.tau, aes(y=(mean/0.1)*100, x=DAY, color ="Tau Effects"), linewidth=2, alpha=0.5) +
#   #geom_ribbon(data=BAGS_out_wide, aes(y=(mean/0.1)*100, x=DAY, ymin = (lci/0.1)*100, ymax=(uci/0.1)*100), alpha = 0.3) +
#   #geom_point(data=LML_sum2, aes(y=(1-mean.ML)*100, x=doy+10), color = "#009E73", size = 3) +
#   #geom_errorbar(data=LML_sum2, aes(y=(1-mean.ML)*100, x=doy+10, ymin = (1-lci.ML)*100, ymax = (1-uci.ML)*100), width=0, color = "#009E73",linewidth=1) +
#   ylab("Litter Bag C Remaining (%)") +
#   xlab("Day") +
#   scale_color_manual(values = c("fW=1"="#E69F00", "CORPSE"="#56B4E9", "Calibrated"="#009E73", "W_SCALAR" = "#F0E442")) +
#   ggtitle("NPP option for turnover") +
#   theme_bw(base_size = 20)
# 
# 
# #plot MIMICS output with different soil moisture and microbial dynamics together
# ggplot() +
#   geom_line(data=BAGS_out_SERC_fW0, aes(y=MICr+MICk, x=DAY, color ="fW=1"), linewidth=2, alpha=0.3) +
#   geom_line(data=BAGS_out_SERC_fW1, aes(y=MICr+MICk, x=DAY, color ="CORPSE"), linewidth=2, alpha=0.3) +
#   geom_line(data=BAGS_out_SERC_fW2, aes(y=MICr+MICk, x=DAY, color ="Calibrated"), linewidth=2, alpha=0.3) +
#   geom_line(data=BAGS_out_SERC_fW3, aes(y=MICr+MICk, x=DAY, color ="W_SCALAR"), linewidth=2, alpha=0.3) +
#   #geom_line(data=BAGS_out, aes(y=MICr+MICk, x=DAY, color ="Default"), linewidth=2, alpha=0.5) +
#   #geom_line(data=BAGS_out_fWt, aes(y=MICr+MICk, x=DAY, color ="Tau Effects"), linewidth=2, alpha=0.5) +
#   ylab("microbial biomass") +
#   xlab("Day") +
#   #ylim(0,2) +
#   xlim(0,3650)+
#   scale_color_manual(values = c("fW=1"="#E69F00", "CORPSE"="#56B4E9", "Calibrated"="#009E73", "W_SCALAR" = "#F0E442")) +
#   ggtitle("beta option for turnover - SERC") +
#   theme_bw(base_size = 20)
# 
# #Mic biomass vs soil moisture
# BO_fWm0 <- BAGS_out_fWm0 %>% filter(DAY<366) %>% filter(Litter_Type=="mean") %>% left_join(SERC_DI, by="DAY")
# BO_fWm1 <- BAGS_out_fWm1 %>% filter(DAY<366) %>% filter(Litter_Type=="mean") %>% left_join(SERC_DI, by="DAY")
# BO_fWm2 <- BAGS_out_fWm2 %>% filter(DAY<366) %>% filter(Litter_Type=="mean") %>% left_join(SERC_DI, by="DAY")
# BO_fWm3 <- BAGS_out_fWm3 %>% filter(DAY<366) %>% filter(Litter_Type=="mean") %>% left_join(SERC_DI, by="DAY")
# ggplot() +
#   geom_line(data=BO_fWm0, aes(y=MICr+MICk, x=GWC, color = "fW=1"), linewidth=2, alpha=0.5) +
#   geom_line(data=BO_fWm1, aes(y=MICr+MICk, x=GWC, color = "CORPSE"), linewidth=2, alpha=0.5) +
#   geom_line(data=BO_fWm2, aes(y=MICr+MICk, x=GWC, color = "Calibrated"), linewidth=2, alpha=0.5)+
#   geom_line(data=BO_fWm3, aes(y=MICr+MICk, x=GWC, color = "W_SCALAR"), linewidth=2, alpha=0.5)+
#   ylim(0,0.7) +
#   ylab("Microbial biomass") +
#   xlab("GWC") +
#   scale_linetype_manual(values = c("fW=1"="E69F00", "CORPSE"="#56B4E9", "Calibrated"="#009E73", "W_SCALAR" = "#F0E442"))


###
#Daily vs. STODE
####

# SERC_daily <- rbind(SERC_DI, SERC_DI)
# SERC_daily$DAY <- 1:732
# BAGS_daily <- inner_join(BAGS_out_SERC_lowANPP, SERC_daily, by='DAY')
# BAGS_out_4y <- filter(BAGS_out_BART, DAY <1461)
# BAGS_BART_sum <- BAGS_out_BART %>% mutate(YEAR = c(rep(1, 365), rep(2, 365), rep(3, 365), rep(4, 365), rep(5, 365))) %>% group_by(YEAR) %>%
#   summarise(daily_mean_mic = mean(MICr+MICk), daily_mean_lit = mean(LITm+LITs), daily_mean_som = mean(SOMa+SOMc+SOMp))
# BAGS_out_BART <- BAGS_out_BART %>% mutate(YEAR = c(rep(1, 365), rep(2, 365), rep(3, 365), rep(4, 365), rep(5, 365))) %>%
#   inner_join(BAGS_BART_sum)
ggplot() +
  geom_line(data=BAGS_out_BART_daily, aes(y=MICr+MICk, x=DAY, color = "MIC", linetype ="daily"), linewidth=2, alpha=0.5) +
  geom_line(data=BAGS_out_BART_daily, aes(y=LITm+LITs, x=DAY, color = "Litter", linetype ="daily"), linewidth=2, alpha=0.5) +
  geom_line(data=BAGS_out_BART_daily, aes(y=SOMa+SOMc+SOMp, x=DAY, color = "SOM", linetype ="daily"), linewidth=2, alpha=0.5) +
  geom_line(data=BAGS_out_BART_SS, aes(y=MICr+MICk, x=DAY, color = "MIC", linetype ="steady state"), linewidth=2, alpha=0.5,) +
  geom_line(data=BAGS_out_BART_SS, aes(y=LITm+LITs, x=DAY, color = "Litter", linetype ="steady state"), linewidth=2, alpha=0.5) +
  geom_line(data=BAGS_out_BART_SS, aes(y=SOMa+SOMc+SOMp, x=DAY, color = "SOM", linetype ="steady state"), linewidth=2, alpha=0.5) +
  #geom_line(data=BAGS_out_BART, aes(y=daily_mean_mic, x=DAY, color = "MIC", linetype ="daily mean"), linewidth=2, alpha=0.5,) +
  #geom_line(data=BAGS_out_BART, aes(y=daily_mean_lit, x=DAY, color = "Litter", linetype ="daily mean"), linewidth=2, alpha=0.5) +
  #geom_line(data=BAGS_out_BART, aes(y=daily_mean_som, x=DAY, color = "SOM", linetype ="daily mean"), linewidth=2, alpha=0.5) +
  #geom_line(data=BAGS_daily, aes(y=ANPP/2, x=DAY, color = "ANPP"), linewidth=2, alpha=0.5) +
  #geom_line(data=BAGS_daily, aes(y=W_SCALAR, x=DAY, color = "W_SCALAR"), linewidth=2, alpha=0.5) +
  #scale_y_log10() +
  scale_color_manual(values = c("MIC"="#E69F00", "Litter"="#56B4E9", "SOM"="#009E73", "W_SCALAR" = "#F0E442", "ANPP" = "#CC79A7")) +
  scale_linetype_manual(values = c("daily"=1, "steady state"=2, "daily mean"=3)) + 
  ylab("C pools") +
  xlab("Day") +
  ggtitle("W_SCALAR moisture at BART") +
  theme_bw(base_size = 20)

###
#within site testing
####


#can varying LQ get the same variability as at the sites?
#within site
# BAGS_out_AllSites$Litter_Type = "mean"
# BAGS_out_AllSites_lci$Litter_Type = "lci"
# BAGS_out_AllSites_uci$Litter_Type = "uci"
# BAGS_out_AllSites <- rbind(BAGS_out_AllSites, BAGS_out_AllSites_lci, BAGS_out_AllSites_uci)
# WIsite <- c("LENO", "SERC", "UNDE")
#comparison of baseline MIMICS to LML and within site
# LIT_init <- BAGS_out_AllSites %>% filter(DAY == 10) %>% mutate(LITi = LITBAGm+LITBAGs) %>% mutate(SITE.LT = paste(SITE, Litter_Type, sep=".")) %>% select(SITE>LT, LITi)
# boxplot(LIT_init$LITi)
# BAGS_out_plot <- BAGS_out_AllSites %>% mutate(SITE.LT = paste(SITE, Litter_Type, sep=".")) %>% left_join(LIT_init, by = "SITE.LT") %>% mutate(LIT_PerLoss = ((LITi - (LITBAGm+LITBAGs))/LITi)*100)
# #within site
# BAGS_out_plot_withinsite <- BAGS_out_plot %>% filter(SITE %in% WIsite)
# LML_WIsite <- LML_sum2 %>% filter(site %in% WIsite)
# BAGS_out_wide = BAGS_out_plot_withinsite %>% filter(SITE=="LENO") %>% mutate(LITBAG = LITBAGm+LITBAGs) %>% select(SITE, Litter_Type, DAY, LITBAG) %>% 
#   pivot_wider(names_from = Litter_Type, values_from = LITBAG)
#plotting
# ggplot() +
#   geom_line(data=BAGS_out_plot, aes(y=100-LIT_PerLoss, x=DAY, group=SITE, color=SITE), linewidth=0.5) +
#   geom_ribbon(data=BAGS_out_wide, aes(y=100-mean, x=DAY, ymin = 100-lci, ymax=100-uci, alpha = 0.3)) +
#   geom_point(data=LML_WIsite, aes(y=100-mean.ML, x=doy+10, group=site, color=site), size = 3) +
#   geom_errorbar(data=LML_WIsite, aes(y=100-mean.ML, x=doy+10, ymin = 100-lci.ML, ymax = 100-uci.ML, group=site, color=site), width=0,linewidth=1) +
#   ylab("Litter Bag C Remaining (%)") +
#   xlab("Day") +
#   theme_bw(base_size = 20)
# #   theme_bw(base_size = 20)
# 
# #can varying soil moisture get the same variability as at the sites?
# #One Site at a time: wide format MIMICS output for plotting
# BAGS_out_wide = BAGS_out %>% mutate(LITBAG = LITBAGm+LITBAGs) %>% select(SITE, Litter_Type, DAY, LITBAG) %>% pivot_wider(names_from = Litter_Type, values_from = LITBAG)
# BAGS_out_wide_wet = BAGS_out %>% mutate(LITBAG = LITBAGm+LITBAGs) %>% select(SITE, Litter_Type, DAY, LITBAG) %>% pivot_wider(names_from = Litter_Type, values_from = LITBAG)
# BAGS_out_wide_dry = BAGS_out %>% mutate(LITBAG = LITBAGm+LITBAGs) %>% select(SITE, Litter_Type, DAY, LITBAG) %>% pivot_wider(names_from = Litter_Type, values_from = LITBAG)
#All sites together
# BAGS_out_AllSites$Mois = "avg"
# BAGS_out_AllSites_wet$Mois = "wet"
# BAGS_out_AllSites_dry$Mois = "dry"
# BAGS_out_AllSites <- rbind(BAGS_out_AllSites, BAGS_out_AllSites_wet, BAGS_out_AllSites_dry)
# WIsite <- c("LENO", "SERC", "UNDE")
# LIT_init <- BAGS_out_AllSites %>% filter(DAY == 10) %>% mutate(LITi = LITBAGm+LITBAGs) %>% mutate(SITE.MT = paste(SITE, Mois, sep=".")) %>% select(SITE.MT, LITi)
# boxplot(LIT_init$LITi)
# BAGS_out_plot <- BAGS_out_AllSites %>% mutate(SITE.MT = paste(SITE, Mois, sep=".")) %>% left_join(LIT_init, by = "SITE.MT") %>% mutate(LIT_PerLoss = ((LITi - (LITBAGm+LITBAGs))/LITi)*100)
# #within site
# BAGS_out_plot_withinsite <- BAGS_out_plot %>% filter(SITE %in% WIsite)
# LML_WIsite <- LML_sum2 %>% filter(site %in% WIsite)
# BAGS_out_wide = BAGS_out_plot_withinsite  %>% mutate(LITBAG = LITBAGm+LITBAGs) %>% select(SITE, Mois, DAY, LIT_PerLoss) %>%
#   pivot_wider(names_from = Mois, values_from = LIT_PerLoss)
# ggplot() +
#   #geom_line(data=BAGS_out_plot_withinsite, aes(y=100-LIT_PerLoss, x=DAY, group=SITE, color=SITE)) +
#   geom_ribbon(data=BAGS_out_wide, aes(y=100-avg, x=DAY, ymin = 100-dry, ymax=100-wet, fill=SITE), alpha=0.5) +
#   #geom_line(data=BAGS_out_wide_dry, aes(y=(mean/0.1)*100, x=DAY), linewidth=1.5, alpha=0.5, color = "#E69F00", linetype =1) +
#   #geom_line(data=BAGS_out_wide, aes(y=(mean/0.1)*100, x=DAY), linewidth=1.5, alpha=0.5, color = "#E69F00", linetype =2) +
#   #geom_line(data=BAGS_out_wide_mid, aes(y=(mean/0.1)*100, x=DAY), linewidth=1.5, alpha=0.5, color = "#E69F00", linetype =3) +
#   #geom_ribbon(data=BAGS_out_wide, aes(y=(mean/0.1)*100, x=DAY, ymin = (lci/0.1)*100, ymax=(uci/0.1)*100), alpha = 0.3) +
#   geom_point(data=LML_WIsite, aes(y=100-mean.ML, x=doy+10, color=site), size = 3) +
#   geom_errorbar(data=LML_WIsite, aes(y=100-mean.ML, x=doy+10, ymin = 100-lci.ML, ymax = 100-uci.ML, color=site), width=0, linewidth=1) +
#   ylab("Litter Bag C Remaining (%)") +
#   xlab("Day") +
#   labs(linetype="Parameter Set") +
#   ggtitle("Variation in litter decomposition with +/- 8% W_SCALAR") +
#   theme_bw(base_size = 20)

####
#plot output - comparing observed to MIMICS microbial community
####

#empirical microbe data
# MSBio_rK <- read.csv("MSBio_rK.csv")
# SERC_rK <- MSBio_rK %>% filter(site == "SERC" & time.point == 0) %>% mutate(rK = r/K) %>% mutate(CO = Copiotroph/Oligotroph)
# 
# #bringing data into one dataframe
# a <- data.frame(group = "model", value = (BAGS_out$MICr/BAGS_out$MICk))
# b <- data.frame(group = "obs_rK", value = SERC_rK$rK)
# c <- data.frame(group = "obs_CO", value = SERC_rK$CO)
# plot.data <- rbind(a,b,c)
# 
# #how does r:K change over time? very little...
# ggplot() + 
#   geom_line(data = BAGS_out, aes(x=DAY, y=(MICr/MICk)))
# #comparing r:K values in empirical and in MIMICS
# ggplot(plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot() +
#   xlab("Type of data") + ylab("r:K or C:O") + ggtitle("SERC Microbial Community Comparison") + 
#   theme_bw(base_size = 20) + theme(legend.position="none")
# #empirical data not directly related to LQ and related to moisture a little bit if you use r:K
# 




