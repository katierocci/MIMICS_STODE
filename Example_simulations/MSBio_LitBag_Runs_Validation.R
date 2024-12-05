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


####
#load MSBio site and litter data and format to code structure
####

MSBio <- read.csv("Example_simulations/Data/Site_annual_clim_validation.csv")
#match input data structure
#AGNPP should be in grams Dry Weight (gDW) not gC! multiply by 2 here to remedy
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


#just one parameter set (default paramters)
site_val <- c("MLBS", "SCBI", "UNDE")
BAGS_out_AllSites_SP_val = data.frame()
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
      BAGS_out_AllSites_SP_val <- rbind(BAGS_out_AllSites_SP_val, BO_DI)
    }
  }
}

#multiple parameter sets (calibrated)
 site_val <- c("MLBS", "SCBI", "UNDE")
 ES_Psets <- read.csv("ES_Psets_5000_NewInputs_ES.csv")
 ES_Psets <- ES_Psets %>% filter(ID==175 | ID==176 | ID==190)
 BAGS_out_AllSites_ES_val = data.frame()
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


#plotting observed litter mass loss vs modeled
BAGS_out_AllSites_ES_val$ID <- as.factor(rep(c(176, 175, 190), each=29565))
BAGS_out_plot.SP <- BAGS_out_AllSites_SP_val %>% mutate(SITE.LT = paste(SITE, Litter_Type, sep=".")) %>% mutate(LIT_PerLoss = ((0.1 - (LITBAGm+LITBAGs))/0.1)*100)
BAGS_out_plot.Cal <- BAGS_out_AllSites_ES_val %>% mutate(SITE.LT = paste(SITE, Litter_Type, sep=".")) %>% mutate(LIT_PerLoss = ((0.1 - (LITBAGm+LITBAGs))/0.1)*100)
#summary data for visualization
val_colors = c("#882255", "#999933", "#332288")
BO_plot_sum.SP <- BAGS_out_plot.SP  %>% group_by(SITE,DAY) %>% summarise(mean=mean(LIT_PerLoss), min=min(LIT_PerLoss), max=max(LIT_PerLoss)) 
BO_plot_sum.Cal <- BAGS_out_plot.Cal  %>% group_by(SITE,DAY) %>% summarise(mean=mean(LIT_PerLoss), min=min(LIT_PerLoss), max=max(LIT_PerLoss)) 
#default
ggplot() +
  geom_ribbon(data=BO_plot_sum.SP, aes(y=100-mean, x=DAY-315, ymin = 100-min, ymax=100-max, group=SITE, fill=SITE, color=SITE), alpha = 0.3, size=0.7) +
  geom_point(data=LML_sum2, aes(y=100-mean.ML, x=doy, group=site, color=site), size = 3) +
  geom_errorbar(data=LML_sum2, aes(y=100-mean.ML, x=doy, ymin = 100-lci.ML, ymax = 100-uci.ML, group=site, color=site), width=0,linewidth=1) +
  ylab("Litter Bag C Remaining (%)") +
  xlab("Day") +
  xlim(0, 780) +
  theme_bw(base_size = 20) +
  scale_color_manual(values=val_colors) + scale_fill_manual(values=val_colors) 
#calibrated
ggplot() +
  geom_ribbon(data=BO_plot_sum.Cal, aes(y=100-mean, x=DAY-315, ymin = 100-min, ymax=100-max, group=SITE, fill=SITE, color=SITE), alpha = 0.3, size=0.7) +
  geom_point(data=LML_sum2, aes(y=100-mean.ML, x=doy, group=site, color=site), size = 3) +
  geom_errorbar(data=LML_sum2, aes(y=100-mean.ML, x=doy, ymin = 100-lci.ML, ymax = 100-uci.ML, group=site, color=site), width=0,linewidth=1) +
  ylab("Litter Bag C Remaining (%)") +
  xlab("Day") +
  xlim(0, 780) +
  theme_bw(base_size = 20) +
  scale_color_manual(values=val_colors) + scale_fill_manual(values=val_colors) 


#formatting data for model vs obs and for effect sizes
FieldData <- LML_sum2 %>% mutate(DAY=doy, SITE=site) %>% mutate(SITE.DAY=paste(SITE, DAY, sep=".")) %>% select(time.point,SITE.DAY, mean.ML, lci.ML, uci.ML)
LITi = 0.1
#default
df_LML.SP <- BAGS_out_AllSites_SP_val %>% mutate(DAY.LitOut = DAY -314) %>% mutate(SITE.DAY=paste(SITE, DAY.LitOut, sep=".")) %>% 
  right_join(FieldData, by="SITE.DAY") %>% mutate(LIT_PerLoss = ((LITi - (LITBAGm+LITBAGs))/LITi)*100)
#calibrated
df_LML.Cal <- BAGS_out_AllSites_ES_val %>% mutate(DAY.LitOut = DAY -314) %>% mutate(SITE.DAY=paste(SITE, DAY.LitOut, sep=".")) %>% 
  right_join(FieldData, by="SITE.DAY") %>% mutate(LIT_PerLoss = ((LITi - (LITBAGm+LITBAGs))/LITi)*100)
#model vs obs
#default
df_LML_sum.SP <- df_LML.SP %>% group_by(SITE, time.point) %>% summarise(mean.ML = mean(mean.ML), m.uci.ML = mean(uci.ML), m.lci.ML = mean(lci.ML), mean.LPL = mean(LIT_PerLoss),
                                                                       n=n(), SE = sd(LIT_PerLoss)/sqrt(n),
                                                                       min.LPL = min(LIT_PerLoss),
                                                                       max.LPL = max(LIT_PerLoss),
                                                                       lci.LPL = mean.LPL - qt(1 - ((1 - 0.95) / 2), n - 1) * SE,
                                                                       uci.LPL = mean.LPL + qt(1 - ((1 - 0.95) / 2), n - 1) * SE)
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
val_colors = c("#882255", "#999933", "#332288")
#plotting
ggplot(df_LML_sum, aes(x=mean.ML, y=mean.LPL)) + geom_point(aes(color=SITE), size=4) + geom_smooth(method = "lm", color="black")  + xlim(0,80) + ylim(0,80) +
  geom_errorbar(aes(ymin=lci.LPL, ymax=uci.LPL, color=SITE), size=1) + geom_errorbarh(aes(xmin=m.lci.ML, xmax=m.uci.ML, color=SITE), size=1) +
  xlab("Observed litter percent C loss") + ylab("Modeled litter percent C loss") + geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) +
  scale_color_manual(values = val_colors) + theme(legend.position="none")


#effect size estimation
#formatting data for statistical model
DI_means <- DailyInput_SM %>% mutate(SITE.SM = paste(SITE, SM_type, sep = ".")) %>% group_by(SITE.SM) %>% 
  summarise(W_SCALAR_mean=mean(W_SCALAR), MAT_mean=mean(MAT)) %>% select(SITE.SM, W_SCALAR_mean, MAT_mean)
#default
MIC_init.SP <- BAGS_out_AllSites_SP_val %>% filter(DAY == 315) %>%mutate(MICrK.i =  MICr/MICk)%>% mutate(MICr.i = MICr)%>% mutate(MICK.i = MICk)%>%
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% select(SITE.SM.LQ, MICrK.i, MICr.i, MICK.i)
#calibrated
MIC_init.Cal.175 <- BAGS_out_AllSites_ES_val %>% filter(ID == 175)%>% filter(DAY == 315) %>%mutate(MICrK.i =  MICr/MICk)%>% mutate(MICr.i = MICr)%>% mutate(MICK.i = MICk)%>%
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% select(SITE.SM.LQ, MICrK.i, MICr.i, MICK.i)
MIC_init.Cal.176 <- BAGS_out_AllSites_ES_val %>% filter(ID == 176)%>% filter(DAY == 315) %>%mutate(MICrK.i =  MICr/MICk)%>% mutate(MICr.i = MICr)%>% mutate(MICK.i = MICk)%>%
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% select(SITE.SM.LQ, MICrK.i, MICr.i, MICK.i)
MIC_init.Cal.190 <- BAGS_out_AllSites_ES_val %>% filter(ID == 190)%>% filter(DAY == 315) %>%mutate(MICrK.i =  MICr/MICk)%>% mutate(MICr.i = MICr)%>% mutate(MICK.i = MICk)%>%
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% select(SITE.SM.LQ, MICrK.i, MICr.i, MICK.i)
BAGS_LIGN <- MSBio_BAGS %>% mutate(SITE.LQ = paste(SITE, TYPE, sep = ".")) %>% select(SITE.LQ, BAG_LIG_N)
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
#logical checks
df_check <- df_analysis %>% filter(MICrK > 0.01) %>%
  filter(MICrK < 100) %>%
  filter(MIC/SOC > 0.0001) %>%
  filter(MIC/SOC < 0.40) 
#logical checks - replace "df_analysis" with the respective df_analyses above to build a statistical model for the default parameters and each of the calibrated parameter sets (175, 176, 190)
#statistical model
Obs_ES_mod <- lmer(LIT_PerLoss ~ scale(log(W_SCALAR_mean))+scale(BAG_LIG_N)+scale(MICrK.i)+ (1|SITE), data = df_check)
Obs_ES <- as.data.frame(fixef(Obs_ES_mod)) 
Obs_ES$Vars <- rownames(Obs_ES)
colnames(Obs_ES)[1] <- "value"
Obs_ES <- Obs_ES[-1, ]
Obs_ES$mult <- ifelse(Obs_ES$value <0, -1, 1)
Obs_ES$rel_ES <- (abs(Obs_ES$value)/sum(abs(Obs_ES$value))) * 100 * Obs_ES$mult
Obs_ES$Vars <- factor(Obs_ES$Vars, levels=c('scale(MICrK.i)','scale(BAG_LIG_N)', 'scale(log(W_SCALAR_mean))'),
                       labels=c('MICrK.i','BAG_LIG_N', 'log_WS'))
Obs_ES$rel_ES
