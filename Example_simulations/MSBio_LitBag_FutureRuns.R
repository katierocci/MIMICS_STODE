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
library(ggnewscale)


source("functions/MIMICS_sim_litterbag.R")
source("functions/MIMICS_calc_steady_state_pools.R")
source("functions/calc_Tpars.R")
source("Parameters/MIMICS_parameters_sandbox_20231129.R")
source("functions/RXEQ.R")
source("functions/MC_parameterization/set_parameter_defaults.R")


##########
# load input data
#############

#loop for annual climate data for future
Mic_sites <- c("SERC","BART","TALL","TREE","LENO","HARV","GRSM")
years = c(2072, 2073, 2074) #choosing future years to model
SSP_AnnualClim = data.frame()
for (site in Mic_sites) {
  NEON_SSP <- read.csv(paste("Example_simulations/Data/",site,"_SSP370_anomalies.csv", sep=""))
  NEON_SSP <- NEON_SSP %>% mutate(SITE=site, LIG_N_sp1=MSBio_sites$LIG_N_sp1[MSBio_sites$SITE == site], LIG_N_sp2=MSBio_sites$LIG_N_sp2[MSBio_sites$SITE == site], 
                                  LIG_N_sp3=MSBio_sites$LIG_N_sp3[MSBio_sites$SITE == site], CLAY=MSBio_sites$CLAY[MSBio_sites$SITE == site], 
                                  lci_SM_ratio=MSBio_sites$lci_SM_ratio[MSBio_sites$SITE == site], uci_SM_ratio=MSBio_sites$uci_SM_ratio[MSBio_sites$SITE == site])
  NEON_SSP_AnnualClim <- NEON_SSP %>% filter(YEAR %in% years) %>% group_by(YEAR) %>%
    summarise(SITE= first(SITE), ANPP = sum(LITFALL)*2, TSOI = mean(TSOI),GWC = mean(H2OSOI*100), W_SCALAR=mean(W_SCALAR),
              LIG_N_sp1 = mean(LIG_N_sp1), LIG_N_sp2 = mean(LIG_N_sp2),LIG_N_sp3 = mean(LIG_N_sp3), CLAY=mean(CLAY),
              lci_SM_ratio = mean(lci_SM_ratio), uci_SM_ratio=mean(uci_SM_ratio)) %>%
    ungroup() %>% summarise(SITE= first(SITE), ANPP = mean(ANPP), TSOI = mean(TSOI),GWC = mean(GWC), W_SCALAR=mean(W_SCALAR),
                            LIG_N_sp1 = mean(LIG_N_sp1), LIG_N_sp2 = mean(LIG_N_sp2),LIG_N_sp3 = mean(LIG_N_sp3), CLAY=mean(CLAY),
                            lci_SM_ratio = mean(lci_SM_ratio), uci_SM_ratio=mean(uci_SM_ratio))
  SSP_AnnualClim <- rbind(SSP_AnnualClim, NEON_SSP_AnnualClim)
  }

#loading daily inputs
Mic_sites <- c("SERC","BART","TALL","TREE","LENO","HARV","GRSM")
years = c(2072, 2073, 2074) 
SSP_DailyInput = data.frame()
for (site in Mic_sites) {
  NEON_SSP <- read.csv(paste("Example_simulations/Data/",site,"_SSP370_anomalies.csv", sep=""))
  DailyInput_SSP <- NEON_SSP %>% filter(YEAR %in% years) %>% group_by(YEAR) %>% mutate(ANPP = sum(LITFALL)*2) %>% ungroup() %>%
    mutate(SITE=site, ANPP = ANPP, LITFALL=LITFALL*2,
    CLAY = rep(MSBio_sites$CLAY[MSBio_sites$SITE == site], 1095), 
    LIG_N = rep(MSBio_sites$LIG_N[MSBio_sites$SITE == site], 1095), 
    GWC = H2OSOI*100, MAT=TBOT) %>%
                  select(SITE, YEAR, DOY, ANPP, LITFALL, TSOI, MAT, CLAY, LIG_N, GWC, W_SCALAR) #%>% filter(YEAR %in% years)
  SSP_DailyInput <- rbind(SSP_DailyInput, DailyInput_SSP)
}

#modifying TALL ANPP as for transient runs
SSP_DailyInput$LITFALL[SSP_DailyInput$SITE == "TALL"] <- SSP_DailyInput$LITFALL[SSP_DailyInput$SITE == "TALL"]*0.60
TALL_ANPP_ByYear <- SSP_DailyInput %>% filter(YEAR %in% years) %>% filter(SITE=="TALL") %>% group_by(YEAR) %>% 
  mutate(ANPP = sum(LITFALL)) %>% filter(DOY==1) %>% select(YEAR, SITE, ANPP)
SSP_DailyInput$ANPP[SSP_DailyInput$SITE == "TALL" & SSP_DailyInput$YEAR == "2072"] <-as.numeric(TALL_ANPP_ByYear[1,3])
SSP_DailyInput$ANPP[SSP_DailyInput$SITE == "TALL" & SSP_DailyInput$YEAR == "2073"] <-as.numeric(TALL_ANPP_ByYear[2,3])
SSP_DailyInput$ANPP[SSP_DailyInput$SITE == "TALL" & SSP_DailyInput$YEAR == "2074"] <-as.numeric(TALL_ANPP_ByYear[3,3])
ANPP.TALL <- mean(TALL_ANPP_ByYear$ANPP) 
SSP_AnnualClim$ANPP[SSP_AnnualClim$SITE == "TALL"] <- ANPP.TALL


#load in MSBio litter bag chemistry
### changed to fMET calculation in STODE script here!! Note that the two options are only somewhat related but less negatives in STODE equation
MSBio4 <- read.csv("Example_simulations/Data/Site_annual_clim_final.csv") %>% mutate(SITE=Site)
MSBio_BAGS <- MSBio4 %>% select(SITE, LIG_N_sp1, LIG_N_sp2, LIG_N_sp3) %>% pivot_longer(2:4, names_to = "TYPE", values_to = "BAG_LIG_N")
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
SSP_AnnualClim_SM <- rbind(SSP_AnnualClim, SSP_AnnualClim, SSP_AnnualClim)
#below creates water scalar over 1 so maybe need to change all maxes where W_SCALAR over 1 is equal to 1? Mathematically, fine to go over 1....
SSP_AnnualClim_SM <- SSP_AnnualClim_SM %>% mutate(SM_type = c(rep("mean",7), rep("max",7), rep("min",7))) %>% 
  mutate(W_SCALAR2 = case_when(SM_type == "mean" ~ W_SCALAR,
                               SM_type == "max" ~ W_SCALAR*uci_SM_ratio,
                               SM_type == "min" ~ W_SCALAR*lci_SM_ratio)) %>%
  mutate(W_SCALAR2 = case_when(W_SCALAR2>1~1, TRUE ~ W_SCALAR2)) %>%
  mutate(W_SCALAR = W_SCALAR2)
SSP_DailyInput_SM <- rbind(SSP_DailyInput, SSP_DailyInput, SSP_DailyInput)
SM_mult <- SSP_AnnualClim %>% select(CLAY, uci_SM_ratio, lci_SM_ratio)
SSP_DailyInput_SM <- SSP_DailyInput_SM %>% left_join(SM_mult, by="CLAY") %>% mutate(SM_type = c(rep("mean", 7665), rep("max", 7665), rep("min", 7665))) %>% 
  mutate(W_SCALAR2 = case_when(SM_type == "mean" ~ W_SCALAR,
                               SM_type == "max" ~ W_SCALAR *uci_SM_ratio,
                               SM_type == "min" ~ W_SCALAR *lci_SM_ratio)) %>%
  mutate(W_SCALAR2 = case_when(W_SCALAR2>1~1, TRUE ~ W_SCALAR2)) %>%
  mutate(W_SCALAR = W_SCALAR2)


#default parameter set
BAGS_out_SP_2070s = data.frame()
SM = c("mean", "max", "min")
for (SM_type2 in SM) {
  MSBio_sites_in <- filter(SSP_AnnualClim_SM, SM_type==SM_type2)
  DailyInput_in <- filter(SSP_DailyInput_SM, SM_type==SM_type2)
  LQ = c("LIG_N_sp1", "LIG_N_sp2", "LIG_N_sp3")
  for (bag_type in LQ) {
    BAGS_mean <- filter(BAGS, TYPE==bag_type)
    MSBio_sites_in$LIG_N = MSBio_sites_in[[bag_type]]
    for (site in Mic_sites) {
      BAGS_input <- filter(BAGS_mean, SITE == site)
      forcing_input <- filter(MSBio_sites_in, SITE == site)
      daily_input <- filter(DailyInput_in, SITE == site)
      BO_DI <- MIMICS_LITBAG(forcing_df = forcing_input, litBAG = BAGS_input, dailyInput = daily_input, nspin_yrs=3, nspin_days=0, litadd_day=315, verbose=T)
      BAGS_out_SP_2070s <- rbind(BAGS_out_SP_2070s, BO_DI)
    }
  }
}

#multiple parameter sets (calibrated model)
ES_Psets <- read.csv("ES_Psets_5000_NewInputs_ES.csv")
BAGS_out_Cal_2070s = data.frame()
ES_Psets <- ES_Psets %>% filter(ID==176 | ID==175 | ID==190)
Pset_ID <- ES_Psets$ID
for (i in Pset_ID) {
  ES_Pset_ID <- filter(ES_Psets, ID == i)
  print(i) #tracking pset
  tau_r <<- c(tau_r_default[1], tau_r_default[2] * ES_Pset_ID$Tau_r[1])
  beta <- beta_default * ES_Pset_ID$beta_x[1]
  vMOD <<- c(vMOD_default[1] * ES_Pset_ID$vMOD_m[1], vMOD_default[2] * ES_Pset_ID$vMOD_s[1], vMOD_default[3], vMOD_default[4] * ES_Pset_ID$vMOD_m[1], vMOD_default[5] * ES_Pset_ID$vMOD_s[1], vMOD_default[6])
  SM = c("mean", "max", "min")
  for (SM_type2 in SM) {
    MSBio_sites_in <- filter(SSP_AnnualClim_SM, SM_type==SM_type2)
    DailyInput_in <- filter(SSP_DailyInput_SM, SM_type==SM_type2)
    LQ = c("LIG_N_sp1", "LIG_N_sp2", "LIG_N_sp3")
    for (bag_type in LQ) {
      BAGS_mean <- filter(BAGS, TYPE==bag_type)
      MSBio_sites_in$LIG_N = MSBio_sites_in[[bag_type]]
      for (site in Mic_sites) {
        BAGS_input <- filter(BAGS_mean, SITE == site)
        forcing_input <- filter(MSBio_sites_in, SITE == site)
        daily_input <- filter(DailyInput_in, SITE == site)
        BO_DI <- MIMICS_LITBAG(forcing_df = forcing_input, litBAG = BAGS_input, dailyInput = daily_input, nspin_yrs=3, nspin_days=0, litadd_day=315, verbose=T)
        BAGS_out_Cal_2070s <- rbind(BAGS_out_Cal_2070s,BO_DI)
      }
    }
  }
}



####
#plot output - ***REQUIRES RUNNING "MSBio_LitBag_Runs.R" TO THIS POINT IN CODE (line 177) FIRST***
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
  filter(site %in% Mic_sites) %>% mutate(SITE=site)


#adding  parameter set IDs
BAGS_out_Cal_2070s$ID <- as.factor(rep(c(176, 175, 190), each=68985))
BAGS_out_AllSites_Cal$ID <- as.factor(rep(c(176, 175, 190), each=68985))
#plotting format for each paramter set and time frame, where 1=sp_2072; 2=cal_2072; 3=sp_2022; 4=cal_2022
BAGS_out_plot1 <- BAGS_out_SP_2070s %>% mutate(LIT_PerLoss = ((0.1 - (LITBAGm+LITBAGs))/0.1)*100)
BAGS_out_plot2 <- BAGS_out_Cal_2070s %>% mutate(LIT_PerLoss = ((0.1 - (LITBAGm+LITBAGs))/0.1)*100)
BAGS_out_plot3 <- BAGS_out_AllSites_SP %>% mutate(LIT_PerLoss = ((0.1 - (LITBAGm+LITBAGs))/0.1)*100)
BAGS_out_plot4 <- BAGS_out_AllSites_Cal %>% mutate(LIT_PerLoss = ((0.1 - (LITBAGm+LITBAGs))/0.1)*100)
#summarizing data for visualization
BO_plot_sum1 <- BAGS_out_plot1 %>%group_by(SITE, DAY,SM_Type) %>% summarise(mean=mean(LIT_PerLoss), min=min(LIT_PerLoss), max=max(LIT_PerLoss))
BO_plot_sum2 <- BAGS_out_plot2 %>%group_by(SITE, DAY,SM_Type) %>% summarise(mean=mean(LIT_PerLoss), min=min(LIT_PerLoss), max=max(LIT_PerLoss))
BO_plot_sum3 <- BAGS_out_plot3 %>%group_by(SITE, DAY,SM_Type) %>% summarise(mean=mean(LIT_PerLoss), min=min(LIT_PerLoss), max=max(LIT_PerLoss))
BO_plot_sum4 <- BAGS_out_plot4 %>%group_by(SITE, DAY,SM_Type) %>% summarise(mean=mean(LIT_PerLoss), min=min(LIT_PerLoss), max=max(LIT_PerLoss))
#plotting decomposition trajectories over time (supplemental figure 8)
ggplot() +
  geom_ribbon(data=BO_plot_sum1, aes(y=100-mean, x=DAY-315, ymin = 100-min, ymax=100-max, group=SITE, fill="Default_future"), alpha = 0.3) +
  geom_ribbon(data=BO_plot_sum2, aes(y=100-mean, x=DAY-315, ymin = 100-min, ymax=100-max, group=SITE, fill="Calibrated_future"), alpha = 0.3) +
  geom_ribbon(data=BO_plot_sum3, aes(y=100-mean, x=DAY-315, ymin = 100-min, ymax=100-max, group=SITE, fill="Default_historical"), alpha = 0.3) +
  geom_ribbon(data=BO_plot_sum4, aes(y=100-mean, x=DAY-315, ymin = 100-min, ymax=100-max, group=SITE, fill="Calibrated_historical"), alpha = 0.3) +
  geom_point(data=LML_sum2, aes(y=100-mean.ML, x=doy, group=SITE), size = 3) +
  geom_errorbar(data=LML_sum2, aes(y=100-mean.ML, x=doy, ymin = 100-lci.ML, ymax = 100-uci.ML, group=SITE), width=0,linewidth=1) +
  ylab("Litter Bag C Remaining (%)") +
  xlab("Day") +
  xlim(0, 780) +
  facet_wrap(.~SITE) +
  theme_bw(base_size = 20) + 
  scale_fill_manual(name='Model',values=c(Default_historical="dodgerblue", Calibrated_historical="pink", Default_future="blue", Calibrated_future="red"))


#other pools 
BAGS_out_SP_2070s$ID = 1
BAGS_out_SP_2070s$type="sp_2070"
BAGS_out_Cal_2070s$type="cal_2070"
BAGS_out_AllSites_SP$ID = 1
BAGS_out_AllSites_SP$type="sp_2018"
BAGS_out_AllSites_Cal$type="cal_2018"
BO_all <- rbind(BAGS_out_SP_2070s, BAGS_out_Cal_2070s, BAGS_out_AllSites_SP, BAGS_out_AllSites_Cal)
BO_all <- BO_all %>% mutate(SITE=factor(SITE, levels=c("TREE", "BART", "HARV", "GRSM", "SERC", "TALL", "LENO"))) #MAT order
#differences between future and historical - means with point under for average time of litterbag deployment 
#metabolic litter
BO_dif  <- BO_all %>% filter(DAY>315 | DAY<964)%>% group_by(SITE, type) %>% summarise(LBm_avg=mean(LITBAGm)) %>% pivot_wider(names_from = "type", values_from = "LBm_avg") %>% mutate(sp=((sp_2070-sp_2018)/sp_2018)*100, cal=((cal_2070-cal_2018)/cal_2018)*100) #%>% mutate(LIT2MIC=(LITs+LITm)/(MICr+MICk)) 
BO_dif2 <- BO_all %>% select(SITE, Litter_Type, SM_Type, type, DAY, LITBAGm) %>% filter(DAY>315 | DAY<964) %>% group_by(SITE, Litter_Type, SM_Type, type) %>% 
  summarise(LBm = mean(LITBAGm)) %>% pivot_wider(names_from = "type", values_from = "LBm") %>% 
  mutate(sp=((sp_2070-sp_2018)/sp_2018)*100, cal=((cal_2070-cal_2018)/cal_2018)*100) 
ggplot() + geom_point(data=BO_dif2, aes(x=SITE, y=sp, group = SITE, colour = SITE, shape = "Default"), size=3) + 
  geom_point(data=BO_dif2, aes(x=SITE, y=cal, group = SITE, colour = SITE, shape = "Calibrated"), size=3) +
  geom_point(data=BO_dif, aes(x=SITE, y=sp, group = SITE, colour = SITE, shape = "Default"), size=8, alpha=0.6, stroke =2) +
  geom_point(data=BO_dif, aes(x=SITE, y=cal, group = SITE, colour = SITE, shape = "Calibrated"), size=8, alpha=0.6, stroke=2) +
  ylab("Percent difference between \n future and historical (%)") + scale_color_manual(values=colorBlind7) + theme_bw(base_size = 16) + 
  scale_shape_manual(name='Model', breaks=c('Default', 'Calibrated'), values=c('Default'=1, 'Calibrated'=16)) +ylim(-31, 8) + theme(legend.position = "none")
#structural litter
BO_dif  <- BO_all %>% filter(DAY>315 | DAY<964)%>% group_by(SITE, type) %>% summarise(LBs_avg=mean(LITBAGs)) %>% pivot_wider(names_from = "type", values_from = "LBs_avg") %>% mutate(sp=((sp_2070-sp_2018)/sp_2018)*100, cal=((cal_2070-cal_2018)/cal_2018)*100) #%>% mutate(LIT2MIC=(LITs+LITm)/(MICr+MICk)) 
BO_dif2 <- BO_all %>% select(SITE, Litter_Type, SM_Type, type, DAY, LITBAGs) %>% filter(DAY>315 | DAY<964) %>% group_by(SITE, Litter_Type, SM_Type, type) %>% 
  summarise(LBs = mean(LITBAGs)) %>% pivot_wider(names_from = "type", values_from = "LBs") %>% 
  mutate(sp=((sp_2070-sp_2018)/sp_2018)*100, cal=((cal_2070-cal_2018)/cal_2018)*100) 
ggplot() + geom_point(data=BO_dif2, aes(x=SITE, y=sp, group = SITE, colour = SITE, shape = "Default"), size=3) + 
  geom_point(data=BO_dif2, aes(x=SITE, y=cal, group = SITE, colour = SITE, shape = "Calibrated"), size=3) +
  geom_point(data=BO_dif, aes(x=SITE, y=sp, group = SITE, colour = SITE, shape = "Default"), size=8, alpha=0.6, stroke =2) +
  geom_point(data=BO_dif, aes(x=SITE, y=cal, group = SITE, colour = SITE, shape = "Calibrated"), size=8, alpha=0.6, stroke=2) +
  ylab("Percent difference between \n future and historical (%)") + scale_color_manual(values=colorBlind7) + theme_bw(base_size = 16) + 
  scale_shape_manual(name='Model', breaks=c('Default', 'Calibrated'), values=c('Default'=1, 'Calibrated'=16)) +ylim(-31, 8) + theme(legend.position = "none")

#loop for summary table of differences between future and historical (supplementtary table 2)
Pools <- c('LITBAGm', 'LITBAGs', 'LITm', 'LITs', 'MICr', 'MICk', 'SOMp', 'SOMc', 'SOMa', 'Decomp_rate_rm', 'Decomp_rate_rs', 'Decomp_rate_km', 'Decomp_rate_ks')
BO_time <- BO_all %>% filter(DAY>315 | DAY<964)
CC_difs_sum <- data.frame()
for (pool in Pools) {
  #pool <- sym(pool)
  pool <- noquote(pool)
  BO_dif_loop  <- BO_time %>% group_by(SITE, type) %>% summarise(avg=mean(get(pool))) %>% 
    pivot_wider(names_from = "type", values_from = "avg") %>% 
    mutate(sp=((sp_2070-sp_2018)/sp_2018)*100, cal=((cal_2070-cal_2018)/cal_2018)*100)
  BO_dif_loop$POOL <- pool
  CC_difs_sum = rbind(CC_difs_sum, BO_dif_loop)
}
#format table for paper
CC_difs_table <- CC_difs_sum %>% select(SITE, POOL, sp, cal) %>% pivot_longer(3:4, names_to = 'Model_type', values_to = 'PerDiff') %>%
  pivot_wider(names_from = POOL, values_from = PerDiff)


#comparing decomposition rates from 2022 to 2072 (for Figure 5)
FieldData <- LML_sum2 %>% mutate(DAY=doy, SITE=site) %>% mutate(SITE.DAY=paste(SITE, DAY, sep=".")) %>% select(time.point,SITE.DAY, mean.ML, lci.ML, uci.ML)
LITi = 0.1
#1=sp_2072; 2=cal_2072; 3=sp_2022; 4=cal_2022
BAGS_out_SP_2070s$ID = 1
BAGS_out_SP_2070s$type="sp_2072"
BAGS_out_Cal_2070s$type="cal_2072"
BAGS_out_AllSites_SP$ID = 1
BAGS_out_AllSites_SP$type="sp_2022"
BAGS_out_AllSites_Cal$type="cal_2022"
df_LML1 <- BAGS_out_SP_2070s %>% mutate(DAY.LitOut = DAY -314) %>% mutate(SITE.DAY=paste(SITE, DAY.LitOut, sep=".")) %>%
  right_join(FieldData, by="SITE.DAY") %>% mutate(LIT_PerLoss = ((LITi - (LITBAGm+LITBAGs))/LITi)*100)
df_LML2 <- BAGS_out_Cal_2070s %>% mutate(DAY.LitOut = DAY -314) %>% mutate(SITE.DAY=paste(SITE, DAY.LitOut, sep=".")) %>%
  right_join(FieldData, by="SITE.DAY") %>% mutate(LIT_PerLoss = ((LITi - (LITBAGm+LITBAGs))/LITi)*100)
df_LML3 <- BAGS_out_AllSites_SP %>% mutate(DAY.LitOut = DAY -314) %>% mutate(SITE.DAY=paste(SITE, DAY.LitOut, sep=".")) %>%
  right_join(FieldData, by="SITE.DAY") %>% mutate(LIT_PerLoss = ((LITi - (LITBAGm+LITBAGs))/LITi)*100)
df_LML4 <- BAGS_out_AllSites_Cal %>% mutate(DAY.LitOut = DAY -314) %>% mutate(SITE.DAY=paste(SITE, DAY.LitOut, sep=".")) %>%
  right_join(FieldData, by="SITE.DAY") %>% mutate(LIT_PerLoss = ((LITi - (LITBAGm+LITBAGs))/LITi)*100)
#1:1 plots
ID_175 <- c(1, 175) #IDs to keep in next lines
ID_176 <- c(1, 176) #IDs to keep in next lines
ID_190 <- c(1, 190) #IDs to keep in next lines
LML_121.175 <- rbind(df_LML1, df_LML2, df_LML3, df_LML4) %>% mutate(S.LT.SM =paste(SITE, Litter_Type, SM_Type, sep="."))  %>% select(SITE, Litter_Type, SM_Type, ID, S.LT.SM, time.point, type, LIT_PerLoss) %>% filter(ID %in% ID_175) %>% select(-ID) %>% pivot_wider(names_from = "type", values_from = LIT_PerLoss) %>% mutate(Pset_cal = "175")
LML_121.176 <- rbind(df_LML1, df_LML2, df_LML3, df_LML4) %>% mutate(S.LT.SM =paste(SITE, Litter_Type, SM_Type, sep="."))  %>% select(SITE, Litter_Type, SM_Type, ID, S.LT.SM, time.point, type, LIT_PerLoss) %>% filter(ID %in% ID_176) %>% select(-ID) %>% pivot_wider(names_from = "type", values_from = LIT_PerLoss) %>% mutate(Pset_cal = "176") 
LML_121.190 <- rbind(df_LML1, df_LML2, df_LML3, df_LML4) %>% mutate(S.LT.SM =paste(SITE, Litter_Type, SM_Type, sep="."))  %>% select(SITE, Litter_Type, SM_Type, ID, S.LT.SM, time.point, type, LIT_PerLoss) %>% filter(ID %in% ID_190) %>% select(-ID) %>% pivot_wider(names_from = "type", values_from = LIT_PerLoss) %>% mutate(Pset_cal = "190") 
LML_121 <- rbind(LML_121.175, LML_121.176, LML_121.190) %>% group_by(SITE, Litter_Type, SM_Type, time.point) %>% summarise(sp_2072.avg =mean(sp_2072), cal_2072.avg =mean(cal_2072), sp_2022.avg =mean(sp_2022), cal_2022.avg =mean(cal_2022), n=n(),
                                                                                                                           lci.sp72 = sp_2072.avg - qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(sp_2072)/sqrt(n)), uci.sp72 = sp_2072.avg + qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(sp_2072)/sqrt(n)),
                                                                                                                           lci.cal72 = cal_2072.avg - qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(cal_2072)/sqrt(n)), uci.cal72 = cal_2072.avg + qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(cal_2072)/sqrt(n)),
                                                                                                                           lci.sp22 = sp_2022.avg - qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(sp_2022)/sqrt(n)), uci.sp22 = sp_2022.avg + qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(sp_2022)/sqrt(n)),
                                                                                                                           lci.cal22 = cal_2022.avg - qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(cal_2022)/sqrt(n)), uci.cal22 = cal_2022.avg + qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(cal_2022)/sqrt(n)))
#ordering sites
LML_121 <- LML_121 %>% mutate(SITE=factor(SITE, levels=c("TREE", "BART", "HARV", "GRSM", "SERC", "TALL", "LENO"))) #MAT order
#three facets for coldest sites, mid sites, and warmest sites
cold <- c('TREE', 'BART', 'HARV')
mid <- c('GRSM', 'SERC')
LML_121_TG <- LML_121 %>% mutate(Temp_group = ifelse(SITE %in% cold, "6.8-8.9", ifelse(SITE %in% mid, "14.7-14.8", "18.1-19.2"))) %>% 
  mutate(Temp_group = factor(Temp_group, levels=c('6.8-8.9', '14.7-14.8', '18.1-19.2')))
ggplot(LML_121_TG) + geom_point(aes(x=sp_2022.avg, y=sp_2072.avg, shape = "Default", color=SITE), size=3, alpha=0.8) + 
  geom_point(aes(x=cal_2022.avg, y=cal_2072.avg, shape="Calibrated", color=SITE), size=3, alpha=0.8) + geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) + theme(legend.position = "none") +
  xlab("Historic litter mass loss (%)") + ylab("Future litter mass loss \n under climate change (%)") + facet_wrap(.~Temp_group) + scale_shape_manual(name='Model', breaks=c('Default', 'Calibrated'), values=c('Default'=1, 'Calibrated'=16)) +scale_color_manual(values = colorBlind7)
#relationships between differences in calibrated vs default responses to climate change and environmental drivers
MSBio_long <- MSBio_sites  %>% mutate(mean=W_SCALAR, min=W_SCALAR*lci_SM_ratio, max=W_SCALAR*uci_SM_ratio) %>% select(-LIG_N, -W_SCALAR) %>% pivot_longer(3:5, names_to = "Litter_Type", values_to = "LIG_N") %>%
  pivot_longer(8:10, names_to = "SM_Type", values_to = "W_SCALAR") %>% mutate(S.LT.ST = paste(SITE, Litter_Type, SM_Type, sep = "."))
lml_dif_drivers3 <- LML_121 %>% mutate(CC_dif = (cal_2072.avg-cal_2022.avg)-(sp_2072.avg-sp_2022.avg)) %>% mutate(S.LT.ST = paste(SITE, Litter_Type, SM_Type, sep = ".")) %>% select(SITE, Litter_Type,SM_Type, time.point, CC_dif, S.LT.ST) %>% inner_join(MSBio_long, by="S.LT.ST")
lml_dif_mod3 <- lmer(CC_dif ~ scale(CLAY) + scale(LIG_N) + scale(TSOI) + scale(lci_SM_ratio) + scale(log(W_SCALAR)) + (1|SITE.x/S.LT.ST), data=lml_dif_drivers3) # 
vif(lml_dif_mod3)
summary(lml_dif_mod3)
Anova(lml_dif_mod3, type=3) 
#sensitivity: all vars = lci and LIG:N; -CLAY: same with TSOI as marginal; -lci: clay and LIG:N; -LIG:N: just lci; -WS: same; -TSOI: same (TLDR not very sensitive to var inclusion)
#interactions: LIG_N*W_SCALAR, LIG_N*CLAY; LIG_N*TSOI - intersting that LIG:N interacts with almost every other variable - maybe also indicates site effect
lml_dif_drivers3 <- lml_dif_drivers3%>% mutate(SITE.x=factor(SITE.x, levels=c("TREE", "BART", "HARV", "GRSM", "SERC", "TALL", "LENO"))) #MAT order
#Litter lignin:N
ggplot(data=lml_dif_drivers3, aes(y=CC_dif, x=LIG_N)) + geom_point(aes(color=SITE.x), size=4, alpha=0.5) +
  geom_smooth(method = "lm", color="black") +scale_color_manual(values=colorBlind7) + theme_bw(base_size = 16) +
  ylab("Percent difference between \ncalibrated and default litter mass \nloss under climate change (%)") +xlab("Litter lignin:N") +theme(legend.position = "none")
#moisture varaiblitiy
ggplot(data=lml_dif_drivers3, aes(y=CC_dif, x=(uci_SM_ratio-1)*100)) + geom_point(aes(color=SITE.x), size=4, alpha=0.5) +
  geom_smooth(method = "lm", color="black") +scale_color_manual(values=colorBlind7) + theme_bw(base_size = 16) +
  ylab("Percent difference between \ncalibrated and default litter mass \nloss under climate change (%)") +xlab("Soil mositure varability (%)") +theme(legend.position = "none")
#adding change in climate vars as drivers
SITE.x = c("BART", "GRSM", "HARV", "LENO", "SERC", "TALL", "TREE")
NPP_CC = c(17, 12, 19, 11, 16, 14, 23)
TSOI_CC = c(26, 14, 35, 10, 14, 9, 37)
WS_CC = c(15, 3, 15, -1, 0, -9, -5)
CC_changes = data.frame(SITE.x, NPP_CC, TSOI_CC, WS_CC)
lml_dif_drivers4 <- lml_dif_drivers3 %>% inner_join(CC_changes, by="SITE.x")
lml_dif_mod4 <- lmer(CC_dif ~ scale(TSOI_CC) + scale(WS_CC) + scale(LIG_N) + scale(CLAY) +scale(lci_SM_ratio)+ (1|SITE.x/S.LT.ST), data=lml_dif_drivers4) #scale(NPP_CC) + 
vif(lml_dif_mod4) #high
Anova(lml_dif_mod4, type=3) 
#adding future values as drivers
SITE.x = c("BART", "GRSM", "HARV", "LENO", "SERC", "TALL", "TREE")
NPP_CC = c(716.5, 954.4, 758.6, 929.6, 865.2, 694.6, 636.6)
TSOI_CC = c(10.7, 16.7, 12.0, 21.0, 16.9, 19.7, 9.3)
WS_CC = c(0.63, 0.82, 0.76, 0.70, 0.62, 0.60, 0.42)
CC_changes = data.frame(SITE.x, NPP_CC, TSOI_CC, WS_CC)
lml_dif_drivers5 <- lml_dif_drivers3 %>% inner_join(CC_changes, by="SITE.x")
lml_dif_mod5 <- lmer(CC_dif ~ scale(TSOI_CC) + scale(log(WS_CC))+ scale(LIG_N) + scale(CLAY) +scale(lci_SM_ratio) + (1|SITE.x/S.LT.ST), data=lml_dif_drivers5) 
vif(lml_dif_mod5) 
Anova(lml_dif_mod5, type=3) 



####
#effect size estimation
####
#prepping data
DI_means <- SSP_DailyInput_SM %>% mutate(SITE.SM = paste(SITE, SM_type, sep = ".")) %>% group_by(SITE.SM) %>% 
  summarise(W_SCALAR_mean=mean(W_SCALAR), MAT_mean=mean(MAT)) %>% select(SITE.SM, W_SCALAR_mean, MAT_mean)
#default
MIC_init.SP <- BAGS_out_SP_2070s %>% filter(DAY == 315) %>%mutate(MICrK.i =  MICr/MICk)%>% mutate(MICr.i = MICr)%>% mutate(MICK.i = MICk)%>%
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% select(SITE.SM.LQ, MICrK.i, MICr.i, MICK.i)
#calibrated
MIC_init.Cal.175 <- BAGS_out_Cal_2070s %>% filter(ID == 175)%>% filter(DAY == 315) %>%mutate(MICrK.i =  MICr/MICk)%>% mutate(MICr.i = MICr)%>% mutate(MICK.i = MICk)%>%
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% select(SITE.SM.LQ, MICrK.i, MICr.i, MICK.i)
MIC_init.Cal.176 <- BAGS_out_Cal_2070s %>% filter(ID == 176)%>% filter(DAY == 315) %>%mutate(MICrK.i =  MICr/MICk)%>% mutate(MICr.i = MICr)%>% mutate(MICK.i = MICk)%>%
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% select(SITE.SM.LQ, MICrK.i, MICr.i, MICK.i)
MIC_init.Cal.190 <- BAGS_out_Cal_2070s %>% filter(ID == 190)%>% filter(DAY == 315) %>%mutate(MICrK.i =  MICr/MICk)%>% mutate(MICr.i = MICr)%>% mutate(MICK.i = MICk)%>%
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% select(SITE.SM.LQ, MICrK.i, MICr.i, MICK.i)
BAGS_LIGN <- MSBio_BAGS %>% mutate(SITE.LQ = paste(SITE, TYPE, sep = ".")) %>% select(SITE.LQ, BAG_LIG_N)
#bringing data together for statistical model
df_analysis.SP <- df_LML1 %>% mutate(MICrK = MICr/MICk) %>% mutate(MIC=MICr+MICk) %>% mutate(SOC = SOMa+SOMc+SOMp) %>% 
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% mutate(SITE.SM = paste(SITE, SM_Type, sep = ".")) %>%
  mutate(SITE.LQ = paste(SITE, Litter_Type, sep = ".")) %>% inner_join(DI_means, by="SITE.SM") %>% inner_join(MIC_init.SP, by="SITE.SM.LQ") %>% 
  inner_join(BAGS_LIGN, by="SITE.LQ") %>% mutate(LQ.SM=paste(Litter_Type, SM_Type, sep="."))
df_analysis.Cal.175 <- df_LML2 %>% filter(ID == 175)%>% mutate(MICrK = MICr/MICk) %>% mutate(MIC=MICr+MICk) %>% mutate(SOC = SOMa+SOMc+SOMp) %>% 
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% mutate(SITE.SM = paste(SITE, SM_Type, sep = ".")) %>%
  mutate(SITE.LQ = paste(SITE, Litter_Type, sep = ".")) %>% inner_join(DI_means, by="SITE.SM") %>% inner_join(MIC_init.Cal.175, by="SITE.SM.LQ") %>% 
  inner_join(BAGS_LIGN, by="SITE.LQ") %>% mutate(LQ.SM=paste(Litter_Type, SM_Type, sep=".")) 
df_analysis.Cal.176 <- df_LML2 %>% filter(ID == 176)%>% mutate(MICrK = MICr/MICk) %>% mutate(MIC=MICr+MICk) %>% mutate(SOC = SOMa+SOMc+SOMp) %>% 
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% mutate(SITE.SM = paste(SITE, SM_Type, sep = ".")) %>%
  mutate(SITE.LQ = paste(SITE, Litter_Type, sep = ".")) %>% inner_join(DI_means, by="SITE.SM") %>% inner_join(MIC_init.Cal.176, by="SITE.SM.LQ") %>% 
  inner_join(BAGS_LIGN, by="SITE.LQ") %>% mutate(LQ.SM=paste(Litter_Type, SM_Type, sep=".")) 
df_analysis.Cal.190 <- df_LML2 %>% filter(ID == 190)%>% mutate(MICrK = MICr/MICk) %>% mutate(MIC=MICr+MICk) %>% mutate(SOC = SOMa+SOMc+SOMp) %>% 
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% mutate(SITE.SM = paste(SITE, SM_Type, sep = ".")) %>%
  mutate(SITE.LQ = paste(SITE, Litter_Type, sep = ".")) %>% inner_join(DI_means, by="SITE.SM") %>% inner_join(MIC_init.Cal.190, by="SITE.SM.LQ") %>% 
  inner_join(BAGS_LIGN, by="SITE.LQ") %>% mutate(LQ.SM=paste(Litter_Type, SM_Type, sep=".")) 
#logical checks - replace "df_analysis" with the respective df_analyses above to build a statistical model for the default parameters and each of the calibrated parameter sets (175, 176, 190)
df_check <- df_analysis.Cal.190 %>% filter(MICrK > 0.01) %>%
  filter(MICrK < 100) %>%
  filter(MIC/SOC > 0.0001) %>%
  filter(MIC/SOC < 0.40) 
#statistical model
Obs_ES_mod <- lmer(LIT_PerLoss ~ scale(log(W_SCALAR_mean))+scale(BAG_LIG_N)+scale(MICrK.i)+ (1|SITE/LQ.SM), data = df_check)
Obs_ES <- as.data.frame(fixef(Obs_ES_mod)) #fixed effects coefficients as effect size
Obs_ES$Vars <- rownames(Obs_ES)
colnames(Obs_ES)[1] <- "value"
Obs_ES <- Obs_ES[-1, ]
Obs_ES$mult <- ifelse(Obs_ES$value <0, -1, 1)
Obs_ES$rel_ES <- (abs(Obs_ES$value)/sum(abs(Obs_ES$value))) * 100 * Obs_ES$mult
Obs_ES$Vars <- factor(Obs_ES$Vars, levels=c('scale(MICrK.i)','scale(BAG_LIG_N)', 'scale(log(W_SCALAR_mean))'),
                      labels=c('MICrK.i','BAG_LIG_N', 'log_WS'))
Obs_ES$rel_ES 