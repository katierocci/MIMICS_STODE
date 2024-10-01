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
# MSBio runs
#############


#-------------------------------
#Using MSBio data
#-------------------------------
####
#load MSBio site and litter data and format to code structure
####

#load baseline site data - needed?
# MSBio <- read.csv("Example_simulations/Data/Site_annual_clim_final.csv")
# #match input data structure
# #AGNPP should be in grams Dry Weight (gDW) not gC! multiply by 2 here to remedy
# #switching AGNPP to LITFALL to match daily inputs!
# #don't have gravimetric soil moisture, just volumetric, assuming a BD of 1g/cm3 makes them equivalent - could be bad assumption given this is BD of leaves
# MSBio2 <- MSBio %>% mutate(SITE = Site, ANPP = LITFALL_sum*2, TSOI = TSOI_mean, CLAY = PCT_CLAY_mean, GWC = H2OSOI_mean*100, W_SCALAR=W_SCALAR_mean) %>%
#   select(SITE, ANPP, TSOI, CLAY, LIG_N, LIG_N_sp1, LIG_N_sp2, LIG_N_sp3, GWC, W_SCALAR, lci_SM_ratio, uci_SM_ratio) 
# #fixing TALL and OSBS ANPP
# #NEON_GPP <- read.csv("Example_simulations/Data/NEON_GPP.csv")
# #MSBio3$ANPP[MSBio3$SITE == "TALL"] <- 510 + 0.41*NEON_GPP[9,2] #using relationship between NEON GPP and ANPP 
# #loading daily inputs and replacing TALL data to be more realistic
# DailyInput <- read.csv("Example_simulations/Data/DailyInput.csv") %>% select(-MAT)
# DailyInput$LITFALL[DailyInput$SITE == "TALL"] <- DailyInput$LITFALL[DailyInput$SITE == "TALL"]*0.60
# DailyInput$ANPP[DailyInput$SITE == "TALL"] <- sum(DailyInput$LITFALL[DailyInput$SITE == "TALL"])
# #replacing MSBio data with daily input sums and means to ensure comparable data between daily data and annual data
# DI_sum <- DailyInput %>% select(SITE, ANPP, TSOI, W_SCALAR) %>% group_by(SITE, ANPP) %>% summarise(TSOI=mean(TSOI), W_SCALAR=mean(W_SCALAR))
# MSBio3 <- MSBio2 %>% select(-ANPP, -TSOI, -W_SCALAR) %>% inner_join(DI_sum, by="SITE")
# #filtering for only sites with microbial data to match observations
# Mic_sites <- c("SERC","BART","TALL","TREE","LENO","HARV","GRSM")
# MSBio_sites <- filter(MSBio3, SITE %in% Mic_sites)


#loop for annual clim data for future
Mic_sites <- c("SERC","BART","TALL","TREE","LENO","HARV","GRSM")
years = c(2072, 2073, 2074) #choosing future years to model
SSP_AnnualClim = data.frame()
#need to do group and then ungroup in below code to make sure ANPP is annual and not over all years
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
#need to do group and then ungroup in below code to make sure ANPP is annual and not over all years
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

#checking all years
########################
#loading daily inputs
#TALL not corrected!!!
Mic_sites <- c("SERC","BART","TALL","TREE","LENO","HARV","GRSM")
years = c(2023:2099) 
SSP_DailyInput_AllYears = data.frame()
#need to do group and then ungroup in below code to make sure ANPP is annual and not over all years
for (site in Mic_sites) {
  NEON_SSP <- read.csv(paste("Example_simulations/Data/",site,"_SSP370_anomalies.csv", sep=""))
  DailyInput_SSP <- NEON_SSP %>% filter(YEAR %in% years) %>% group_by(YEAR) %>% mutate(ANPP = sum(LITFALL)*2) %>% ungroup() %>%
    mutate(SITE=site, ANPP = ANPP, LITFALL=LITFALL*2, #1095
           CLAY = rep(MSBio_sites$CLAY[MSBio_sites$SITE == site], 28105), #1095
           LIG_N = rep(MSBio_sites$LIG_N[MSBio_sites$SITE == site], 28105), #1095
           GWC = H2OSOI*100, MAT=TBOT) %>%
    select(SITE, YEAR, DOY, ANPP, LITFALL, TSOI, MAT, CLAY, LIG_N, GWC, W_SCALAR) #%>% filter(YEAR %in% years)
  SSP_DailyInput_AllYears <- rbind(SSP_DailyInput_AllYears, DailyInput_SSP)
}

ggplot(SSP_DailyInput_AllYears, aes(x=YEAR, y=ANPP, group = YEAR)) + geom_boxplot() +theme_bw() +facet_wrap(.~SITE)

#corrected TALL for selected years
colorBlind7  <- c("#E69F00", "#56B4E9", "#009E73",
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 
DI_sum <- DailyInput %>% group_by(SITE) %>% summarise(ANPP=mean(ANPP))
ggplot()+ geom_line(data=SSP_DailyInput, aes(x=YEAR, y=ANPP, color=SITE), size=4) + ylab(expression(paste("ANPP (g DW m"^"-2"*"yr"^"-1"*")"))) +
  geom_point(data=DI_sum, aes(x=2073, y=ANPP, color=SITE), size=5, alpha=0.7) +theme_bw(base_size = 16) + scale_color_manual(values=colorBlind7) #alpha doesn't work here bc the same points are overlid a bunch of times!
SSP_DI_Sum <- SSP_DailyInput %>% group_by(SITE) %>% summarise(ANPP=mean(ANPP)/2, W_SCALAR = mean(W_SCALAR), TSOI = mean(TSOI))
SSP_DI_Clim <- SSP_DailyInput %>% group_by(SITE, DOY) %>% summarise(LF=mean(LITFALL)/2, W_SCALAR = mean(W_SCALAR), TSOI = mean(TSOI)) 
ggplot() + geom_line(data = SSP_DI_Clim, aes(x=DOY, y=1), color="orange", linewidth=2, alpha=0.7) + ylim(0,4) +
  geom_line(data = DailyInput, aes(x=DAY, y=2), color="darkgreen", linewidth=2, alpha=0.5) + facet_wrap(.~factor(SITE, levels = c("TREE", "BART", "HARV", "GRSM", "SERC", "TALL", "LENO"))) + theme_bw(base_size = 16) 
############################################

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


#all sites and all litters with daily input looping through different soil moistures and litters as well
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

#species specific litter at steady state - multiple parameter sets
#315th day of the year is 11/11/21 which is the average day the litter was deployed
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
  filter(site %in% Mic_sites) %>% mutate(SITE=site)

#comparison of baseline MIMICS to LML
#adding Pset IDs
BAGS_out_Cal_2070s$ID <- as.factor(rep(c(176, 175, 190), each=68985))
BAGS_out_AllSites_Cal$ID <- as.factor(rep(c(176, 175, 190), each=68985))
#for 10 year runs
#BAGS_out_Cal_2070s$ID <- as.factor(rep(c(176, 175, 190), each=229950))
#BAGS_out_AllSites_Cal$ID <- as.factor(rep(c(176, 175, 190), each=229950))
# BAGS_out_AllSites_test$ID <- as.factor(1:18)
# ID_test <- filter(BAGS_out_AllSites_test, DAY==50)
BAGS_out_plot <- BAGS_out_SP_2070s %>% mutate(LIT_PerLoss = ((0.1 - (LITBAGm+LITBAGs))/0.1)*100)
#wide format for plotting
BAGS_out_wide = BAGS_out_plot %>% select(SITE, ID, Litter_Type, SM_Type, DAY, LIT_PerLoss) %>% #ID, 
  pivot_wider(names_from = Litter_Type, values_from = LIT_PerLoss)
#plotting - check!
#BAGS_out_wide_BART <- filter(BAGS_out_wide, SITE=='LENO')
ggplot() +
  geom_line(data=BAGS_out_wide%>% filter(SITE=="TREE"), aes(y=100-LIG_N_sp1, x=DAY-315, group=SITE, color=SITE), linewidth=0.5, alpha=0.5) +
  geom_line(data=BAGS_out_wide%>% filter(SITE=="TREE"), aes(y=100-LIG_N_sp2, x=DAY-315, group=SITE, color=SITE), linewidth=0.5, alpha=0.5) +
  geom_line(data=BAGS_out_wide%>% filter(SITE=="TREE"), aes(y=100-LIG_N_sp3, x=DAY-315, group=SITE, color=SITE), linewidth=0.5, alpha=0.5) +
  #geom_ribbon(data=BAGS_out_wide, aes(y=100-LIG_N, x=DAY, ymin = 100-LIG_N_min, ymax=100-LIG_N_max, group=SITE, fill=SITE), alpha = 0.3) +
  geom_point(data=LML_sum2 , aes(y=100-mean.ML, x=doy), size = 3) +
  geom_errorbar(data=LML_sum2 , aes(y=100-mean.ML, x=doy, ymin = 100-lci.ML, ymax = 100-uci.ML), width=0,linewidth=1) +
  xlim(0, 780) +
  ylab("Litter Bag C Remaining (%)") +
  xlab("Days elapsed") +
  facet_wrap(.~SM_Type) +
  theme_bw(base_size = 20)
#summary data - works ok for visualization!
#1=sp_2070; 2=cal_2070; 3=sp_2018; 4=cal_2018
LML_sum2 <- LML_sum2 %>% mutate(SITE=factor(SITE, levels=c("SERC","TREE", "BART", "GRSM", "LENO", "HARV", "TALL"))) #LIG:N order
BO_plot_sum1 <- BAGS_out_plot %>%group_by(SITE, DAY,SM_Type) %>% summarise(mean=mean(LIT_PerLoss), min=min(LIT_PerLoss), max=max(LIT_PerLoss))
BO_plot_sum1 <- BO_plot_sum1 %>% mutate(SITE=factor(SITE, levels=c("TREE", "TALL", "BART", "HARV", "SERC", "LENO", "GRSM"))) #ANPP order
ggplot() +
  geom_ribbon(data=BO_plot_sum1, aes(y=100-mean, x=DAY-315, ymin = 100-min, ymax=100-max, group=SITE, fill="Default_future"), alpha = 0.3) +
  geom_ribbon(data=BO_plot_sum2, aes(y=100-mean, x=DAY-315, ymin = 100-min, ymax=100-max, group=SITE, fill="Calibrated_future"), alpha = 0.3) +
  geom_ribbon(data=BO_plot_sum3, aes(y=100-mean, x=DAY-315, ymin = 100-min, ymax=100-max, group=SITE, fill="Default_historical"), alpha = 0.3) +
  geom_ribbon(data=BO_plot_sum4, aes(y=100-mean, x=DAY-315, ymin = 100-min, ymax=100-max, group=SITE, fill="Calibrated_historical"), alpha = 0.3) +
  #geom_line(data=BO_plot_sum, aes(y=100-mean, x=DAY-315, group=SITE, color=SITE), linewidth=2, alpha = 0.3) +
  geom_point(data=LML_sum2, aes(y=100-mean.ML, x=doy, group=SITE), size = 3) +
  geom_errorbar(data=LML_sum2, aes(y=100-mean.ML, x=doy, ymin = 100-lci.ML, ymax = 100-uci.ML, group=SITE), width=0,linewidth=1) +
  ylab("Litter Bag C Remaining (%)") +
  xlab("Day") +
  xlim(0, 780) +
  facet_wrap(.~SITE) +
  theme_bw(base_size = 20) + 
  scale_fill_manual(name='Model',values=c(Default_historical="dodgerblue", Calibrated_historical="pink", Default_future="blue", Calibrated_future="red"))

#10-year decomp values
BO_plot_sum <- rbind(BO_plot_sum1%>%mutate(type="sp_2072"), BO_plot_sum2%>%mutate(type="cal_2072"), 
                     BO_plot_sum3%>%mutate(type="sp_2022"), BO_plot_sum4%>%mutate(type="cal_2022")) %>%filter(DAY>314&SM_Type=="mean")
ggplot(BO_plot_sum, aes(x=SITE, y=100-mean,fill=type)) + geom_boxplot() +theme_bw(base_size = 16) +
  ylab("Litter mass \n remaining at 10 years") + scale_fill_manual(values = c("pink", "red", "dodgerblue", "blue"))
write.csv(BO_plot_sum, "BO_plot_sum_10yr.csv")
#asymtote decomp values
BO_plot_sum <- rbind(BO_plot_sum1%>%mutate(type="sp_2072"), BO_plot_sum2%>%mutate(type="cal_2072"), 
                     BO_plot_sum3%>%mutate(type="sp_2022"), BO_plot_sum4%>%mutate(type="cal_2022")) %>%filter(DAY==1000)
ggplot(BO_plot_sum, aes(x=SITE, y=100-mean,fill=type)) + geom_boxplot() +theme_bw(base_size = 16) +
  ylab("Litter mass \n remaining at 1000 days") + scale_fill_manual(values = c("pink", "red", "dodgerblue", "blue"))

#other pools besides litter decomp
BAGS_out_SP_2070s$ID = 1
BAGS_out_SP_2070s$type="sp_2070"
BAGS_out_Cal_2070s$type="cal_2070"
BAGS_out_AllSites_SP$ID = 1
BAGS_out_AllSites_SP$type="sp_2018"
BAGS_out_AllSites_Cal$type="cal_2018"
BO_all <- rbind(BAGS_out_SP_2070s, BAGS_out_Cal_2070s, BAGS_out_AllSites_SP, BAGS_out_AllSites_Cal)
BO_all <- BO_all %>% mutate(SITE=factor(SITE, levels=c("TREE", "BART", "HARV", "GRSM", "SERC", "TALL", "LENO"))) #MAT order
#BO_all <- BO_all %>% mutate(SITE=factor(SITE, levels=c("TREE", "TALL", "BART", "HARV", "SERC", "LENO", "GRSM"))) #ANPP order
ggplot() +
  geom_boxplot(data=BO_all, aes(y=(MICr/(MICr+MICk))/(MICk/(MICr+MICk)), x=SITE, group=interaction(type, SITE), fill = type))+
  scale_fill_manual(values=c("pink", "red", "dodgerblue", "blue")) +
  ylab("Litter") +
  xlab("SITE") +
  theme_bw(base_size = 20)
ggplot() +
  geom_boxplot(data=BO_all %>% filter(SITE=='TALL') %>% filter(DAY>601), aes(y=Decomp_rate_ks, x=Litter_Type, group=interaction(type, Litter_Type), fill = type))+
  scale_fill_manual(values=c("pink", "red", "dodgerblue", "blue")) +
  ylab("Litter") +
  xlab("SITE") +
  theme_bw(base_size = 20)
#differences between future and historical
BO_dif  <- BO_all%>% group_by(SITE, type) %>% summarise(LBs_avg=mean(LITBAGs)) %>% pivot_wider(names_from = "type", values_from = "LBs_avg") %>% mutate(sp=((sp_2070-sp_2018)/sp_2018)*100, cal=((cal_2070-cal_2018)/cal_2018)*100) #%>% mutate(LIT2MIC=(LITs+LITm)/(MICr+MICk)) 
#means
ggplot(data=BO_dif) +
  geom_point(data=BO_dif, aes(y=sp, x=SITE,color = SITE, shape="Default"), size=4, stroke=2)+
  geom_point(data=BO_dif, aes(y=cal, x=SITE,color = SITE, shape= "Calibrated"),  size=4, stroke=2)+
  ylab("Percent difference in mean \n value between 2072 and 2022") +
  xlab("SITE") +
  theme_bw(base_size = 20) +
  scale_shape_manual(name='Model', breaks=c('Default', 'Calibrated'), values=c('Default'=1, 'Calibrated'=16)) +
  scale_color_manual(values = colorBlind7)
#means and points under
BO_dif2 <- BO_all %>% select(SITE, Litter_Type, SM_Type, type, DAY, LITBAGs) %>% filter(DAY>315) %>% group_by(SITE, Litter_Type, SM_Type, type) %>% 
  summarise(LBs = mean(LITBAGs)) %>% pivot_wider(names_from = "type", values_from = "LBs") %>% 
  mutate(sp=((sp_2070-sp_2018)/sp_2018)*100, cal=((cal_2070-cal_2018)/cal_2018)*100) 
ggplot() + geom_point(data=BO_dif2, aes(x=SITE, y=sp, group = SITE, colour = SITE, shape = "Default"), size=3) + 
  geom_point(data=BO_dif2, aes(x=SITE, y=cal, group = SITE, colour = SITE, shape = "Calibrated"), size=3) +
  geom_point(data=BO_dif, aes(x=SITE, y=sp, group = SITE, colour = SITE, shape = "Default"), size=8, alpha=0.6, stroke =2) +
  geom_point(data=BO_dif, aes(x=SITE, y=cal, group = SITE, colour = SITE, shape = "Calibrated"), size=8, alpha=0.6, stroke=2) +
  ylab("Percent difference between \n future and historical (%)") + scale_color_manual(values=colorBlind7) + theme_bw(base_size = 16) + 
  scale_shape_manual(name='Model', breaks=c('Default', 'Calibrated'), values=c('Default'=1, 'Calibrated'=16)) 
#just TALL
BO_dif_TALL  <- BO_all %>% filter(SITE=='TALL')  %>% filter(DAY>315 | DAY<601) %>% group_by(Litter_Type, type) %>% summarise(MICk_avg=mean(MICk)) %>% pivot_wider(names_from = "type", values_from = "MICk_avg") %>% mutate(sp=((sp_2070-sp_2018)/sp_2018)*100, cal=((cal_2070-cal_2018)/cal_2018)*100) #%>% mutate(LIT2MIC=(LITs+LITm)/(MICr+MICk)) 
ggplot(data=BO_dif_TALL) +
  geom_point(aes(y=sp, x=Litter_Type, shape="Default"), size=4, stroke=2, color="#D55E00")+
  geom_point(aes(y=cal, x=Litter_Type, shape= "Calibrated"), size=4, stroke=2, color="#D55E00")+
  ylab("Percent difference in mean oligotorphic \n biomass between 2072 and 2022") +
  xlab("SITE") +
  theme_bw(base_size = 20) +
  scale_shape_manual(name='Model', breaks=c('Default', 'Calibrated'), values=c('Default'=1, 'Calibrated'=16))
#loop for summary table of differences between future and historical
Pools <- c('LITBAGm', 'LITBAGs', 'LITm', 'LITs', 'MICr', 'MICk', 'SOMp', 'SOMc', 'SOMa', 'Decomp_rate_rm', 'Decomp_rate_rs', 'Decomp_rate_km', 'Decomp_rate_ks')
CC_difs_sum <- data.frame()
for (pool in Pools) {
  #pool <- sym(pool)
  pool <- noquote(pool)
  BO_dif_loop  <- BO_all %>% group_by(SITE, type) %>% summarise(avg=mean(get(pool))) %>% 
    pivot_wider(names_from = "type", values_from = "avg") %>% 
    mutate(sp=((sp_2070-sp_2018)/sp_2018)*100, cal=((cal_2070-cal_2018)/cal_2018)*100)
  BO_dif_loop$POOL <- pool
  CC_difs_sum = rbind(CC_difs_sum, BO_dif_loop)
}
#format table for paper
CC_difs_table <- CC_difs_sum %>% select(SITE, POOL, sp, cal) %>% pivot_longer(3:4, names_to = 'Model_type', values_to = 'PerDiff') %>%
  pivot_wider(names_from = POOL, values_from = PerDiff)
write.csv(CC_difs_table, 'CC_difs_table.csv')

#comparing decomposition rates from 2022 to 2072
FieldData <- LML_sum2 %>% mutate(DAY=doy, SITE=site) %>% mutate(SITE.DAY=paste(SITE, DAY, sep=".")) %>% select(time.point,SITE.DAY, mean.ML, lci.ML, uci.ML)
LITi = 0.1
#1=sp_2070; 2=cal_2070; 3=sp_2018; 4=cal_2018
BAGS_out_SP_2070s$ID = 1
BAGS_out_SP_2070s$type="sp_2072"
BAGS_out_Cal_2070s$type="cal_2072"
BAGS_out_AllSites_SP$ID = 1
BAGS_out_AllSites_SP$type="sp_2022"
BAGS_out_AllSites_Cal$type="cal_2022"
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
#Litter_Type, SM_Type, 
#LML_121 <- LML_121 %>% mutate(SITE=factor(SITE, levels=c("TREE", "TALL", "BART", "HARV", "SERC", "LENO", "GRSM"))) #ANPP order
LML_121 <- LML_121 %>% mutate(SITE=factor(SITE, levels=c("TREE", "BART", "HARV", "GRSM", "SERC", "TALL", "LENO"))) #MAT order
#sp vs calibrated
#ggplot(LML_121) + geom_point(aes(x=sp_2070, y=cal_2070, shape = Litter_Type, color="Transient"), size=3, alpha=0.5) + geom_point(aes(x=sp_2018, y=cal_2018, shape=Litter_Type, color="SSP-370"), size=3, alpha=0.5) + geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) + 
#  xlab("Default model litter mass loss (%)") + ylab("Calibrated model litter mass loss (%)") + facet_wrap(.~SITE) + scale_color_manual(name='Type', breaks=c('Transient', 'SSP-370'), values=c('Transient'='red', 'SSP-370'='blue'))
#transient vs CC - one Pset as a time
ggplot(LML_121) + geom_point(aes(x=sp_2022.avg, y=sp_2072.avg, shape = as.factor(time.point), color="Default"), size=3, alpha=0.5) + geom_point(aes(x=cal_2022.avg, y=cal_2072.avg, shape=as.factor(time.point), color="Calibrated"), size=3, alpha=0.5) + geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) + 
  #geom_errorbar(aes(ymin=lci.LPL, ymax=uci.LPL, color=SITE), size=1) + geom_errorbarh(aes(xmin=m.lci.ML, xmax=m.uci.ML, color=SITE), size=1) + #hmm how to connec tthese??
  xlab("Transient litter mass loss (%)") + ylab("SSP-370 litter mass loss (%)") + facet_wrap(.~SITE) + scale_color_manual(name='Type', breaks=c('Default', 'Calibrated'), values=c('Default'='red', 'Calibrated'='blue'))
#transient vs CC - all three psets w/ errorbars
ggplot(LML_121) + geom_point(aes(x=sp_2022.avg, y=sp_2072.avg, shape = as.factor(time.point), color="Default"), size=3, alpha=0.5) + 
  geom_errorbar(aes(x=sp_2022.avg, y=sp_2072.avg, ymin=lci.sp72, ymax=uci.sp72, color='Default'), size=1) + geom_errorbarh(aes(x=sp_2022.avg, y=sp_2072.avg, xmin=lci.sp22, xmax=uci.sp22, color='Default'), size=1) +
  geom_point(aes(x=cal_2022.avg, y=cal_2072.avg, shape=as.factor(time.point), color="Calibrated"), size=3, alpha=0.5) + geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) + 
  geom_errorbar(aes(x=cal_2022.avg, y=cal_2072.avg, ymin=lci.cal72, ymax=uci.cal72, color='Calibrated'), size=1) + geom_errorbarh(aes(x=cal_2022.avg, y=cal_2072.avg, xmin=lci.cal22, xmax=uci.cal22, color='Calibrated'), size=1) +
  xlab("Transient litter mass loss (%)") + ylab("SSP-370 litter mass loss (%)") + facet_wrap(.~SITE) + scale_color_manual(name='Type', breaks=c('Default', 'Calibrated'), values=c('Default'='red', 'Calibrated'='blue'))
#transient vs CC - all three psets w/ just averages
ggplot(LML_121) + geom_point(aes(x=sp_2022.avg, y=sp_2072.avg, shape = "Default", color=SITE), size=3, alpha=0.8) + 
  geom_point(aes(x=cal_2022.avg, y=cal_2072.avg, shape="Calibrated", color=SITE), size=3, alpha=0.8) + geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) + 
  xlab("Historic litter mass loss (%)") + ylab("Future litter mass loss under climate change (%)") + facet_wrap(.~SITE) + scale_shape_manual(name='Model', breaks=c('Default', 'Calibrated'), values=c('Default'=1, 'Calibrated'=16)) +scale_color_manual(values = colorBlind7)
#all sites on one plot?
ggplot(LML_121) + geom_point(aes(x=sp_2022.avg, y=sp_2072.avg, shape = "Default", color=SITE), size=4, alpha=0.8) + 
  geom_point(aes(x=cal_2022.avg, y=cal_2072.avg, shape="Calibrated", color=SITE), size=4, alpha=0.8) + geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) + 
  xlab("Historic litter mass loss (%)") + ylab("Future litter mass loss under climate change (%)") + scale_shape_manual(name='Model', breaks=c('Default', 'Calibrated'), values=c('Default'=1, 'Calibrated'=16)) +scale_color_manual(values = colorBlind7)
#three facets for coldest sites, mid sites, and warmest sites
cold <- c('TREE', 'BART', 'HARV')
mid <- c('GRSM', 'SERC')
LML_121_TG <- LML_121 %>% mutate(Temp_group = ifelse(SITE %in% cold, "6.8-8.9", ifelse(SITE %in% mid, "14.7-14.8", "18.1-19.2"))) %>% 
  mutate(Temp_group = factor(Temp_group, levels=c('6.8-8.9', '14.7-14.8', '18.1-19.2')))
ggplot(LML_121_TG) + geom_point(aes(x=sp_2022.avg, y=sp_2072.avg, shape = "Default", color=SITE), size=3, alpha=0.8) + 
  geom_point(aes(x=cal_2022.avg, y=cal_2072.avg, shape="Calibrated", color=SITE), size=3, alpha=0.8) + geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) + 
  xlab("Historic litter mass loss (%)") + ylab("Future litter mass loss \n under climate change (%)") + facet_wrap(.~Temp_group) + scale_shape_manual(name='Model', breaks=c('Default', 'Calibrated'), values=c('Default'=1, 'Calibrated'=16)) +scale_color_manual(values = colorBlind7)
#showing LQ with color
LQ <- MSBio_BAGS %>% mutate(SITE.LT = paste(SITE, TYPE, sep=".")) %>% select(SITE.LT, BAG_LIG_N) %>% mutate(BAG_LIG_N = as.numeric(BAG_LIG_N)) %>% mutate(BAG_LIG_N.2 = as.numeric(BAG_LIG_N))
LML_121_LQ <- LML_121 %>% mutate(SITE.LT = paste(SITE, Litter_Type, sep=".")) %>% inner_join(LQ, by='SITE.LT') 
ggplot(LML_121_LQ) + geom_point(aes(x=sp_2022.avg, y=sp_2072.avg, shape = as.factor(time.point), color=BAG_LIG_N), size=3, alpha=0.5) + 
  scale_color_viridis_c(option = 'rocket', direction = -1, end = 0.7, begin=0.3) + new_scale_color() +
  geom_point(aes(x=cal_2022.avg, y=cal_2072.avg, shape=as.factor(time.point), color=BAG_LIG_N.2), size=3, alpha=0.5) + scale_color_viridis_c(option = 'mako', direction = -1, end = 0.7, begin=0.3) + #backwards direction so need to do end and begin as opposites
  geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) + 
  xlab("Transient litter mass loss (%)") + ylab("SSP3-7.0 litter mass loss (%)") + facet_wrap(.~SITE) 
#showing LQ with shapes
LML_121_shape <- LML_121 %>% mutate(LT.TP = paste(Litter_Type, time.point, sep = "."))
ggplot(LML_121_shape) + geom_point(aes(x=sp_2022.avg, y=sp_2072.avg, shape = LT.TP, color="Default"), size=3, alpha=0.5) + 
  geom_point(aes(x=cal_2022.avg, y=cal_2072.avg, shape=LT.TP, color="Calibrated"), size=3, alpha=0.5) + geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) + 
  xlab("Transient litter mass loss (%)") + ylab("SSP3-7.0 litter mass loss (%)") + facet_wrap(.~SITE) + 
  scale_color_manual(name='Type', breaks=c('Default', 'Calibrated'), values=c('Default'='red', 'Calibrated'='blue')) + scale_shape_manual(values = c(1, 19, 2, 17, 0, 15))
#plotting difference in LML under CC for each site and model types 
LML_wide <- LML_121 %>%  mutate(LPL_sp=sp_2072.avg-sp_2022.avg, LPL_cal=cal_2072.avg-cal_2022.avg) %>% filter(time.point==2) %>%
  select(SITE, Litter_Type, SM_Type, time.point, LPL_sp, LPL_cal)
LML_wide_sum <- LML_wide %>% group_by(SITE) %>% summarise(mean_LPL_sp=mean(LPL_sp), mean_LPL_cal=mean(LPL_cal), min_LPL_sp=min(LPL_sp), min_LPL_cal=min(LPL_cal), max_LPL_sp=max(LPL_sp), max_LPL_cal=max(LPL_cal))
#big and small points
ggplot() + geom_point(data=LML_wide, aes(x=SITE, y=LPL_sp, group = SITE, colour = SITE, shape = "Default"), size=3) + 
  geom_point(data=LML_wide, aes(x=SITE, y=LPL_cal, group = SITE, colour = SITE, shape = "Calibrated"), size=3) +
  geom_point(data=LML_wide_sum, aes(x=SITE, y=mean_LPL_sp, group = SITE, colour = SITE, shape = "Default"), size=8, alpha=0.6) +
  geom_point(data=LML_wide_sum, aes(x=SITE, y=mean_LPL_cal, group = SITE, colour = SITE, shape = "Calibrated"), size=8, alpha=0.6) +
  ylab("2072-2022 litter mass loss (%)") + scale_color_manual(values=colorBlind7) + theme_bw(base_size = 16) + 
  scale_shape_manual(name='Model', breaks=c('Default', 'Calibrated'), values=c('Default'=1, 'Calibrated'=16)) 
#just TALL big and small points
ggplot() + geom_point(data=LML_wide %>%filter(SITE=='TALL'), aes(x=Litter_Type, y=LPL_sp, group = Litter_Type, colour = Litter_Type, shape = "Default"), size=3) + 
  geom_point(data=LML_wide %>%filter(SITE=='TALL'), aes(x=Litter_Type, y=LPL_cal, group = Litter_Type, colour = Litter_Type, shape = "Calibrated"), size=3) +
  ylab("2072-2022 litter mass loss (%)") + scale_color_manual(values=colorBlind7) + theme_bw(base_size = 16) + 
  scale_shape_manual(name='Model', breaks=c('Default', 'Calibrated'), values=c('Default'=1, 'Calibrated'=16))
#points with min and max
ggplot() +
  geom_pointrange(data=LML_wide_sum, aes(x=SITE, y=mean_LPL_sp, ymin=min_LPL_sp, ymax=max_LPL_sp, group = SITE, colour = SITE, shape = "Default"), size=2, alpha=0.8) +
  geom_pointrange(data=LML_wide_sum, aes(x=SITE, y=mean_LPL_cal, ymin=min_LPL_cal, ymax=max_LPL_cal, group = SITE, colour = SITE, shape = "Calibrated"), size=2, alpha=0.8) +
  ylab("2072-2022 litter mass loss (%)") + scale_color_manual(values=colorBlind7) + theme_bw(base_size = 16) + 
  scale_shape_manual(name='Model', breaks=c('Default', 'Calibrated'), values=c('Default'=1, 'Calibrated'=16)) 
#points just means
ggplot() +
  geom_point(data=LML_wide_sum, aes(x=SITE, y=mean_LPL_sp, group = SITE, colour = SITE, shape = "Default"), size=10, alpha=0.8, stroke=2) +
  geom_point(data=LML_wide_sum, aes(x=SITE, y=mean_LPL_cal, group = SITE, colour = SITE, shape = "Calibrated"), size=10, alpha=0.8, stroke =2) +
  ylab(" Mean future-historic \n litter mass loss (%)") + scale_color_manual(values=colorBlind7) + theme_bw(base_size = 24) + 
  scale_shape_manual(name='Model', breaks=c('Default', 'Calibrated'), values=c('Default'=1, 'Calibrated'=16)) +ylim(0,20)
#violins
ggplot(LML_wide) + geom_boxplot(aes(x=SITE, y=LPL_sp, group = SITE, colour = SITE, fill = "Default"), size=1) + 
  geom_boxplot(aes(x=SITE, y=LPL_cal, group = SITE, colour = SITE, fill = "Calibrated"), size=1) +  
  ylab("2070-2018 litter mass loss (%)") + scale_color_manual(values=colorBlind7) + theme_bw(base_size = 16) + 
  scale_fill_manual(name='Model', breaks=c('Default', 'Calibrated'), values=c('Default'='white', 'Calibrated'='grey'))

#comparing means from 2070-2018
lml4<- df_LML4 %>% group_by(time.point, SITE) %>% summarise(LPL = mean(LIT_PerLoss), sd= sd(LIT_PerLoss)) %>% mutate(ID='cal_2022')
lml <- rbind(lml1, lml2, lml3, lml4)
#lml <- lml %>% mutate(SITE=factor(SITE, levels=c("TREE", "BART", "HARV", "GRSM", "SERC", "TALL", "LENO"))) #MAT order
#lml <- lml %>% mutate(SITE=factor(SITE, levels=c("TREE", "BART", "SERC", "TALL", "HARV", "LENO", "GRSM"))) #moisture order
#lml <- lml %>% mutate(SITE=factor(SITE, levels=c("SERC","TREE", "BART", "GRSM", "LENO", "HARV", "TALL"))) #LIG:N order
lml <- lml %>% mutate(SITE=factor(SITE, levels=c("TREE", "TALL", "BART", "HARV", "SERC", "LENO", "GRSM"))) #ANPP order
ggplot(lml, aes(x=as.factor(time.point), y=LPL)) + geom_point(aes(color=ID), size=4, alpha=0.8) + geom_errorbar(aes(ymin=LPL-sd, ymax = LPL+sd), width=0) +
  theme_bw(base_size = 16) + scale_color_manual(values=c("pink", "red", "dodgerblue","blue")) +facet_grid(.~SITE)
lml_wide <- lml %>%  pivot_wider(names_from = "ID", values_from = 3:4) %>% mutate(LPL_sp=LPL_sp_2072-LPL_sp_2022, LPL_cal=LPL_cal_2072-LPL_cal_2022,
                                                                                  sd_sp=sd_sp_2072+sd_sp_2022, sd_cal=sd_cal_2072+sd_cal_2022) %>%
  select(time.point, SITE, LPL_sp, LPL_cal, sd_sp, sd_cal) %>% pivot_longer(cols=3:6, names_to = c(".value","Pset"), names_sep = "_")
ggplot(lml_wide) + geom_point(aes(x=as.factor(time.point), y=LPL, group = Pset, colour = Pset), size=4, alpha=0.8) + 
  geom_errorbar(aes(x=as.factor(time.point), ymin = LPL-sd, ymax=LPL+sd, group = Pset, colour = Pset), width=0.2, alpha=0.8) + 
  ylab("2070-2018 litter mass loss (%)") + scale_color_manual(values=c("red", "blue")) + theme_bw(base_size = 16)+ facet_grid(.~SITE)
#drivers of climate change responses for calibrated an default models
#lml_CC_drivers <- lml_wide %>% inner_join(MSBio_sites, by="SITE")
#lml_CC_mod <- lmer(LPL ~ scale(CLAY)+scale(ANPP)+scale(TSOI)+scale(W_SCALAR)+scale(LIG_N)+(1|SITE), data=lml_CC_drivers %>% filter(Pset=="cal"))
#summary(lml_CC_mod)
#lml_CC_drivers_cal <- lml_CC_drivers %>% filter(Pset=="cal")
#lml_CC_drivers_sp <- lml_CC_drivers %>% filter(Pset=="sp")
#plot(lml_CC_drivers_cal$TSOI, lml_CC_drivers_cal$LPL)
#taking a look at realtionships
#plot(lml_dif_drivers$CLAY, lml_dif_drivers$LPL_dif)
lml_dif  <- lml %>%  pivot_wider(names_from = "ID", values_from = 3:4) %>% mutate(LPL_sp=LPL_sp_2072-LPL_sp_2022, LPL_cal=LPL_cal_2072-LPL_cal_2022,
                                                                                  sd_sp=sd_sp_2072+sd_sp_2022, sd_cal=sd_cal_2072+sd_cal_2022) %>%
  select(time.point, SITE, LPL_sp, LPL_cal, sd_sp, sd_cal) %>% mutate(LPL_dif=LPL_cal-LPL_sp, sd_dif=sd_cal+sd_sp)
ggplot(lml_dif) + geom_point(aes(x=as.factor(time.point), y=LPL_dif), size=4, alpha=0.8, color="purple") + 
  geom_errorbar(aes(x=as.factor(time.point), ymin = LPL_dif-sd_dif, ymax=LPL_dif+sd_dif), width=0.2, alpha=0.8, color="purple") + 
  ylab("Difference between calibrate and starting \n point litter mass loss under climate change (%)") + theme_bw(base_size = 10)+ facet_grid(.~SITE)
#drivers of differences between calibrated and default models
lml_dif_drivers <- lml_dif %>% inner_join(MSBio_sites, by="SITE")
lml_dif_mod <- lmer(LPL_dif ~ scale(CLAY)+scale(TSOI)+scale(W_SCALAR)+scale(LIG_N) + (1|SITE), data=lml_dif_drivers) #scale(ANPP)+
lml_dif_mod_lm <- lm(LPL_dif ~ scale(TSOI)+scale(W_SCALAR) + scale(LIG_N) +scale(CLAY), data=lml_dif_drivers) #dropping ANPP cause the least change in R2 # +scale(ANPP)
ANPP_mod <- lm(scale(ANPP) ~ scale(TSOI)+scale(W_SCALAR) + scale(LIG_N) +scale(CLAY), data=lml_dif_drivers) 
vif(lml_dif_mod)
summary(lml_dif_mod) 
#taking a look at realtionships
plot(lml_dif_drivers$ANPP, lml_dif_drivers$LPL_dif)
ggplot(data=lml_dif_drivers, aes(y=LPL_dif, x=scale(ANPP))) + geom_point() +geom_smooth(method = "lm")

#by LQ instead of site
lml1_LT <- df_LML1 %>% group_by(time.point, SITE, Litter_Type) %>% summarise(LPL = mean(LIT_PerLoss), sd= sd(LIT_PerLoss)) %>% mutate(ID='sp_2070')
lml_LT <- rbind(lml1_LT, lml2_LT, lml3_LT, lml4_LT)
SITE_LQ <- SSP_AnnualClim %>% select(SITE, LIG_N_sp1, LIG_N_sp2, LIG_N_sp3) %>% pivot_longer(2:4, names_to = 'Litter_Type', values_to = 'LIG_N') %>% mutate(SITE.LQ=paste(SITE, Litter_Type, sep = "."))
lml_LT_dif  <- lml_LT %>%  pivot_wider(names_from = "ID", values_from = 4:5) %>% mutate(LPL_sp=LPL_sp_2070-LPL_sp_2018, LPL_cal=LPL_cal_2070-LPL_cal_2018,
                                                                                  sd_sp=sd_sp_2070+sd_sp_2018, sd_cal=sd_cal_2070+sd_cal_2018) %>%
  select(time.point, SITE, Litter_Type, LPL_sp, LPL_cal, sd_sp, sd_cal) %>% mutate(LPL_dif=LPL_cal-LPL_sp, sd_dif=sd_cal+sd_sp) %>% mutate(SITE.LQ=paste(SITE, Litter_Type, sep = ".")) %>% 
  inner_join(SITE_LQ, by="SITE.LQ") %>% mutate(LQ_bin = ifelse(LIG_N<20, "<20", ifelse(LIG_N<30, "20-30", ifelse(LIG_N<40, "30-40", ">40")))) %>%
  mutate(LQ_bin=factor(LQ_bin, levels=c("<20","20-30", "30-40", ">40"))) 
ggplot(lml_LT_dif) + geom_point(aes(x=as.factor(time.point), y=LPL_dif, color=LIG_N), size=4) +
  geom_boxplot(aes(x=as.factor(time.point), y=LPL_dif), alpha=0.2, fill="black") +  scale_color_gradient(low="yellow", high="red") +
  ylab("Difference between calibrate and starting \n point litter mass loss under climate change (%)") + theme_bw(base_size = 16)+ facet_grid(.~LQ_bin)

#comparison of modeled and observed decomp at timepoint 1 and 2
#site observed mean mass loss versus individual obs of modeled data
modelVobs <- lm(LIT_PerLoss~mean.ML, data=df_LML)
summary(modelVobs)
#RMSE
sqrt(mean((df_LML$mean.ML - df_LML$LIT_PerLoss)^2))
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
df_LML_best <- filter(df_LML, ID == 20 | ID==3 | ID==5)
df_LML_sum <- df_LML_best %>% group_by(SITE, time.point) %>% summarise(mean.ML = mean(mean.ML), m.uci.ML = mean(uci.ML), m.lci.ML = mean(lci.ML), mean.LPL = mean(LIT_PerLoss),
                                                                       n=n(), SE = sd(LIT_PerLoss)/sqrt(n),
                                                                       min.LPL = min(LIT_PerLoss),
                                                                       max.LPL = max(LIT_PerLoss),
                                                                       lci.LPL = mean.LPL - qt(1 - ((1 - 0.95) / 2), n - 1) * SE,
                                                                       uci.LPL = mean.LPL + qt(1 - ((1 - 0.95) / 2), n - 1) * SE)
modelVobs <- lm(mean.LPL~mean.ML, data=df_LML_sum)
summary(modelVobs)
sqrt(mean((df_LML_sum$mean.ML - df_LML_sum$mean.LPL)^2))
ggplot(df_LML_sum, aes(x=mean.ML, y=mean.LPL)) + geom_point(aes(color=SITE), size=4) + geom_smooth(method = "lm", color="black")  + xlim(0,80) + ylim(0,80) +
  geom_errorbar(aes(ymin=lci.LPL, ymax=uci.LPL, color=SITE), size=1) + geom_errorbarh(aes(xmin=m.lci.ML, xmax=m.uci.ML, color=SITE), size=1) +
  xlab("Observed litter percent C loss") + ylab("Modeled litter percent C loss") + geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16)


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
DI_means <- SSP_DailyInput_SM %>% mutate(SITE.SM = paste(SITE, SM_type, sep = ".")) %>% group_by(SITE.SM) %>% 
  summarise(W_SCALAR_mean=mean(W_SCALAR), MAT_mean=mean(MAT)) %>% select(SITE.SM, W_SCALAR_mean, MAT_mean)
#initial MICrK
#BAGS_out_AllSites_ES$ID <- rep(1:18, each = 22995)
MIC_init <- BAGS_out_SP_2070s %>% filter(DAY == 315) %>%mutate(MICrK.i =  MICr/MICk)%>%
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% select(SITE.SM.LQ, MICrK.i) #%>% filter(ID == 4)
#need bag means!
BAGS_LIGN <- MSBio_BAGS %>% mutate(SITE.LQ = paste(SITE, TYPE, sep = ".")) %>% select(SITE.LQ, BAG_LIG_N)
df_analysis <- df_LML1 %>% mutate(MICrK = MICr/MICk) %>% mutate(MIC=MICr+MICk) %>% mutate(SOC = SOMa+SOMc+SOMp) %>% 
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% mutate(SITE.SM = paste(SITE, SM_Type, sep = ".")) %>%
  mutate(SITE.LQ = paste(SITE, Litter_Type, sep = ".")) %>% inner_join(DI_means, by="SITE.SM") %>% 
  inner_join(MIC_init, by="SITE.SM.LQ") %>% inner_join(BAGS_LIGN, by="SITE.LQ") %>% mutate(LQ.SM = paste(Litter_Type, SM_Type)) # %>% filter(ID == 3)
#logical checks
df_check <- df_analysis %>% filter(MICrK > 0.01) %>%
  filter(MICrK < 100) %>%
  filter(MIC/SOC > 0.0001) %>%
  filter(MIC/SOC < 0.40) 
#plot((df_check$MICr/df_check$MICk), df_check$BAG_LIG_N)
df_check$log_WS <- log(df_check$W_SCALAR_mean)
#effect size
#df_ES <- df_check %>% mutate_at(vars(c("log_WS", "BAG_LIG_N", "MICrK.i")), ~(scale(.) %>% as.vector))
Obs_ES_mod <- lmer(LIT_PerLoss ~ scale(log(W_SCALAR_mean))+scale(BAG_LIG_N)+scale(MICrK.i)+ (1|SITE/LQ.SM), data = df_check) #MAT #scale(BAG_LIG_N)+
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
ggplot(Obs_ES, aes(x=Vars, y=rel_ES)) + geom_bar(stat="identity", fill="darkolivegreen3") + coord_flip() + 
  geom_text(aes(label=round(rel_ES, digits=1), vjust=1.5), size=5) +theme_bw(base_size = 16)
Obs_ES175<- Obs_ES
rbind(Obs_ES175, Obs_ES176, Obs_ES190) %>% group_by(Vars) %>% summarise(mean.rES = mean(rel_ES), sd.rES = sd(rel_ES)) %>%
  ggplot(aes(x=Vars, y=mean.rES)) + geom_bar(stat="identity", fill="darkolivegreen3", aes(group=Vars)) +  geom_errorbar(aes(ymax=mean.rES+sd.rES, ymin=mean.rES-sd.rES), width=0.1, size=1) + 
  coord_flip() + geom_text(aes(label=round(mean.rES, digits=1), vjust=1.5), size=5) +theme_bw(base_size = 16) #, hjust=1.1

