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
MSBio <- read.csv("Example_simulations/Data/Site_annual_clim_final.csv")
#match input data structure
#AGNPP should be in grams Dry Weight (gDW) not gC! multiply by 2 here to remedy
#switching AGNPP to LITFALL to match daily inputs!
#don't have gravimetric soil moisture, just volumetric, assuming a BD of 1g/cm3 makes them equivalent - could be bad assumption given this is BD of leaves
MSBio2 <- MSBio %>% mutate(SITE = Site, ANPP = LITFALL_sum*2, TSOI = TSOI_mean, CLAY = PCT_CLAY_mean, GWC = H2OSOI_mean*100, W_SCALAR=W_SCALAR_mean) %>%
  select(SITE, ANPP, TSOI, CLAY, LIG_N, LIG_N_sp1, LIG_N_sp2, LIG_N_sp3, GWC, W_SCALAR, lci_SM_ratio, uci_SM_ratio) 
#fixing TALL and OSBS ANPP
NEON_GPP <- read.csv("Example_simulations/Data/NEON_GPP.csv")
#MSBio3$ANPP[MSBio3$SITE == "TALL"] <- 510 + 0.41*NEON_GPP[9,2] #using relationship between NEON GPP and ANPP 
#loading daily inputs and replacing TALL data to be more realistic
DailyInput <- read.csv("Example_simulations/Data/DailyInput.csv") %>% select(-MAT)
DailyInput$LITFALL[DailyInput$SITE == "TALL"] <- DailyInput$LITFALL[DailyInput$SITE == "TALL"]*0.60
DailyInput$ANPP[DailyInput$SITE == "TALL"] <- sum(DailyInput$LITFALL[DailyInput$SITE == "TALL"])
#replacing MSBio data with daily input sums and means to ensure comparable data between daily data and annual data
DI_sum <- DailyInput %>% select(SITE, ANPP, TSOI, W_SCALAR) %>% group_by(SITE, ANPP) %>% summarise(TSOI=mean(TSOI), W_SCALAR=mean(W_SCALAR))
MSBio3 <- MSBio2 %>% select(-ANPP, -TSOI, -W_SCALAR) %>% inner_join(DI_sum, by="SITE")
#filtering for only sites with microbial data to match observations
Mic_sites <- c("SERC","BART","TALL","TREE","LENO","HARV","GRSM")
MSBio_sites <- filter(MSBio3, SITE %in% Mic_sites)

#additional code for: (1) creating DailyInput file for all sites; (2) determining litterfall multiplier for TALL; (3) checking litterfall vs ANPP for inputs; (4) OLD: determining ANPP and LitFall multipliers when we were using both
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
#(3)
#need to start from raw data
# RawDaily = data.frame()
# sites_daily <- c("BART", "GRSM", "HARV", "LENO", "SERC", "TALL", "TREE")
# for (site in sites_daily) {
#   daily <-  read.csv(paste("Example_simulations/Data/",site,"_clim.csv", sep=""))
#   daily$SITE = site
#   RawDaily <- rbind(RawDaily, daily)  
# }
# ggplot(RawDaily %>%filter(SITE=='TALL'), aes(x=X, y=LITFALL*2)) + geom_point() +theme_bw() +xlab("DOY") #+ facet_grid(.~SITE)
# ggplot(RawDaily %>%filter(SITE=='TALL'), aes(x=X, y=AGNPP*2)) + geom_point() +theme_bw() +xlab("DOY") #+ facet_grid(.~SITE)
#(4)
# #replace TALL with more realistic data - LENO litterfall and ANPP are similar but TALL litterfall still very off
# test1 <- filter(DailyInput, DAY == 40)
# test2 <- inner_join(MSBio_sites, test1, by = "SITE")
# ggplot(test2, aes(x=ANPP.x, y=ANPP.y, color=SITE)) + geom_point(size=3) + theme_bw(base_size = 16) #comparing adjusted annual ANPP (ANPP.x) to sum of litterfall
# LF_GoodSites <- test2 %>% filter(SITE != "TALL")
# LF_mod <- lm(ANPP.y~ANPP.x, data=LF_GoodSites)
# summary(LF_mod) #-7.2 + 1.7x #derecho: 101 + 1.9x
# TALL_ANPP.adj <- 101 + 1.9*MSBio_sites[6,2]
# TALL_mult <- TALL_ANPP.adj/test1[6,2] #0.663 #new: 0.62



#load in MSBio litter bag chemistry
### changed to fMET calculation in STODE script here!! Note that the two options are only somewhat related but less negatives in STODE equation
MSBio_BAGS <- MSBio_sites %>% select(SITE, LIG_N_sp1, LIG_N_sp2, LIG_N_sp3) %>% pivot_longer(2:4, names_to = "TYPE", values_to = "BAG_LIG_N")
MSBio_BAGS$CALC_MET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * (MSBio_BAGS$BAG_LIG_N))
MSBio_BAGS$CALC_MET[MSBio_BAGS$CALC_MET <0] = 0.01 #setting negatives to small number so 99% structural

BAG_init_size <- 100
BAGS <- MSBio_BAGS %>% select(SITE, TYPE, CALC_MET)
BAGS$BAG_LITm <- ((BAG_init_size * 1e3 / 1e4)/ depth) * BAGS$CALC_MET #g/m2 converted to mg/cm3
BAGS$BAG_LITs <- ((BAG_init_size * 1e3 / 1e4)/ depth) * (1-BAGS$CALC_MET) 
#initial litter = 0.33 because of unit conversions here


####
#run litterbag model 
####


# #Individual site example (SERC - row 8)
# BAGS_TREE <- filter(BAGS_sites, Site == "TREE" & TYPE == "mean")
# BAGS_TREE <- BAGS_TREE[,2:5]
# BAGS_out_TREE_SS <- BAGS_TREE %>% split(1:nrow(BAGS_TREE)) %>% map(~ MIMICS_LITBAG(litBAG=.,
#                                                                            forcing_df=MSBio_sites[7,],
#                                                                            dailyInput = TREE_DI, 
#                                                                            nspin_yrs=2,
#                                                                            nspin_days=0,
#                                                                            litadd_day=10,
#                                                                            verbose=T)) %>% bind_rows()
# 
# #all sites and all litters
# BAGS_mean <- filter(BAGS_sites, TYPE=="mean")
# BAGS_input <- split(BAGS_mean, 1:nrow(BAGS_mean))
# forcing_input <- split(MSBio_sites, 1:nrow(MSBio_sites))
# BAGS_out_AllSites <- map2(forcing_input, BAGS_input, ~MIMICS_LITBAG(forcing_df = .x, litBAG = .y, nspin_yrs=2, nspin_days=0, litadd_day=10, verbose=T)) %>% bind_rows()
# 
# #all sites and all litters with daily input
# #switch to loop since there isn't a map function that can handle vectors and lists together (I don't think)
# BAGS_mean <- filter(BAGS, TYPE=="LIG_N")
# BAGS_out_AllSites_DI = data.frame()
# for (site in Mic_sites) {
#   BAGS_input <- filter(BAGS_mean, SITE == site)
#   forcing_input <- filter(MSBio_sites, SITE == site)
#   daily_input <- filter(DailyInput, SITE == site)
#   BO_DI <- MIMICS_LITBAG(forcing_df = forcing_input, litBAG = BAGS_input, dailyInput = daily_input, nspin_yrs=2, nspin_days=0, litadd_day=10, verbose=T) 
#   BAGS_out_AllSites_DI <- rbind(BAGS_out_AllSites_DI,BO_DI)
# }

#all sites and all litters with daily input looping through different soil moistures and litters as well
#prep the data
#create different soil moisture for steady state and daily input data
MSBio_sites_SM <- rbind(MSBio_sites, MSBio_sites, MSBio_sites)
#below creates water scalar over 1 so maybe need to change all maxes where W_SCALAR over 1 is equal to 1? Mathematically, fine to go over 1....
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

#mean litter at steady state
# BAGS_out_AllSites_var = data.frame()
# SM = c("mean", "max", "min")
# for (SM_type2 in SM) {
#   MSBio_sites_in <- filter(MSBio_sites_SM, SM_type==SM_type2)
#   DailyInput_in <- filter(DailyInput_SM, SM_type==SM_type2)
#   LQ = c("LIG_N_sp1", "LIG_N_sp2", "LIG_N_sp3")
#     for (bag_type in LQ) {
#       BAGS_mean <- filter(BAGS, TYPE==bag_type)
#         for (site in Mic_sites) {
#           BAGS_input <- filter(BAGS_mean, SITE == site)
#           forcing_input <- filter(MSBio_sites_in, SITE == site)
#           daily_input <- filter(DailyInput_in, SITE == site)
#           BO_DI <- MIMICS_LITBAG(forcing_df = forcing_input, litBAG = BAGS_input, dailyInput = daily_input, nspin_yrs=3, nspin_days=0, litadd_day=315, verbose=T)
#           BAGS_out_AllSites_var <- rbind(BAGS_out_AllSites_var,BO_DI)
#         }
#       }
# }

#just one parameter set - must do this one first to ensure a modifyied parameter sets isn't input to the starting point model
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


#species specific litter at steady state - multiple parameter sets
#315th day of the year is 11/11/21 which is the average day the litter was deployed
ES_Psets <- read.csv("ES_Psets_5000_NewInputs_ES.csv")
ES_Psets <- ES_Psets %>% filter(ID==175 | ID==176 | ID==190) 
BAGS_out_AllSites_Cal = data.frame()
Pset_ID <- ES_Psets$ID
for (i in Pset_ID) {
  ES_Pset_ID <- filter(ES_Psets, ID == i)
  print(i) #tracking pset
   tau_r <<- c(tau_r_default[1], tau_r_default[2] * ES_Pset_ID$Tau_r[1])
   #tau_K <<- c(tau_K_default[1], tau_K_default[2] * ES_Pset_ID$Tau_K[1])
   #CUE <<- CUE_default * ES_Pset_ID$CUE_x[1]
   #vMOD <<- vMOD_default * ES_Pset_ID$vMOD_x[1]
   beta <- beta_default * ES_Pset_ID$beta_x[1]
   #beta <- c(beta_default[1] * ES_Pset_ID$beta_r[1], beta_default[2] * ES_Pset_ID$beta_k[1])
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

#write.csv(BAGS_out_AllSites_ES, "BO_NI_vMODms.csv")
#write.csv(BAGS_out_AllSites_ES2, "BO_looseES.csv")

# #testing pset multiplier loop - working besides beta which cannot be changed for some reason....
# RWA_Psets <- read.csv('RWA_Psets_MinMax.csv')
# ES_Psets <- read.csv('ES_Psets_MinMax.csv')
# BAGS_out_AllSites_test = data.frame()
# forcing_BART=filter(MSBio_sites_SM, SITE == 'BART' & SM_type=='mean')
# BAGS_Sp1=filter(BAGS, SITE == 'BART' & TYPE=='LIG_N_sp1')
# daily_mean=filter(DailyInput_SM, SITE == 'BART' & SM_type=='mean')
# Pset_ID <- ES_Psets$ID
# for (i in Pset_ID) {
#   ES_Pset_ID <- filter(ES_Psets, ID == i)
#   print(i) #tracking pset
#   tau_r <<- c(tau_r_default[1], tau_r_default[2] * ES_Pset_ID$Tau_r[1])
#   tau_K <<- c(tau_K_default[1], tau_K_default[2] * ES_Pset_ID$Tau_K[1])
#   CUE <<- CUE_default * ES_Pset_ID$CUE_x[1]
#   vMOD <<- vMOD_default * ES_Pset_ID$vMOD_x[1]
#   beta <- beta_default * ES_Pset_ID$beta_x[1] #works without double-headed arrow, not sure why double-headed arrow isn't working
#   BO_DI <- MIMICS_LITBAG(forcing_df = forcing_BART, litBAG = BAGS_Sp1, dailyInput = daily_mean, nspin_yrs=1, nspin_days=0, litadd_day=10, verbose=T)
#   BAGS_out_AllSites_test <- rbind(BAGS_out_AllSites_test,BO_DI)
# }




# LQ = c("LIG_N_sp1", "LIG_N_sp2", "LIG_N_sp3")
# for (bag_type in LQ) {
#   BAGS_mean <- filter(BAGS, TYPE==bag_type)
#   MSBio_sites_in$LIG_N = MSBio_sites_in[[bag_type]]
#   print(MSBio_sites_in$LIG_N)
# }

# #just LQ
# BAGS_out_AllSites_var = data.frame()
# LQ = c("LIG_N", "LIG_N_max", "LIG_N_min")
# for (bag_type in LQ) {
#   BAGS_mean <- filter(BAGS, TYPE==bag_type)
#   for (site in Mic_sites) {
#     BAGS_input <- filter(BAGS_mean, SITE == site)
#     forcing_input <- filter(MSBio_sites, SITE == site)
#     daily_input <- filter(DailyInput, SITE == site)
#     BO_DI <- MIMICS_LITBAG(forcing_df = forcing_input, litBAG = BAGS_input, dailyInput = daily_input, nspin_yrs=2, nspin_days=0, litadd_day=10, verbose=T)
#     BAGS_out_AllSites_var <- rbind(BAGS_out_AllSites_var,BO_DI)
#   }
# }
# 
# 
# #just one site soil moisture loop
#   BAGS_out_AllSites_test = data.frame()
#   SM = c("mean", "max", "min")
#   for (SM_type2 in SM) {
#     MSBio_sites_in <- filter(MSBio_sites_SM, SM_type==SM_type2)
#     DailyInput_in <- filter(DailyInput_SM, SM_type==SM_type2)
#       BAGS_mean <- filter(BAGS, TYPE=="LIG_N")
#       BAGS_input <- filter(BAGS_mean, SITE == "BART")
#         forcing_input <- filter(MSBio_sites_in, SITE == "BART")
#         daily_input <- filter(DailyInput_in, SITE == "BART")
#         BO_DI <- MIMICS_LITBAG(forcing_df = forcing_input, litBAG = BAGS_input, dailyInput = daily_input, nspin_yrs=0.5, nspin_days=0, litadd_day=10, verbose=T)
#         BAGS_out_AllSites_test <- rbind(BAGS_out_AllSites_test,BO_DI)
#       }
# 
#   #just one site daily input
#   BAGS_mean <- filter(BAGS, TYPE=="LIG_N")
#   BAGS_input <- filter(BAGS_mean, SITE == "BART")
#   forcing_input <- filter(MSBio_sites, SITE == "BART")
#   daily_input <- filter(DailyInput, SITE ==  "BART")
#   BO_BART <- MIMICS_LITBAG(forcing_df = forcing_input, litBAG = BAGS_input, dailyInput = daily_input, nspin_yrs=2, nspin_days=0, litadd_day=10, verbose=T) 
#   
  
####
#plot output
####

colorBlind7  <- c("#E69F00", "#56B4E9", "#009E73",
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #yellow (LENO), blue (SERC), green (UNDE)
MSBio_palette <- c("#313695", "#4575B4", "#74ADD1", "#FDAE61", "#F46D43", "#D73027", "#A50026") #temperature color order

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
filter(site %in% Mic_sites)

#quick little MAOM-moisture test
#SM_comp <- MSBio_sites_SM %>% mutate(SITE.SM = paste(SITE, SM_type)) %>% select(SITE.SM, W_SCALAR)
#BAGS_out_AllSites_SP %>% mutate(SITE.SM = paste(SITE, SM_Type)) %>% inner_join(SM_comp, by="SITE.SM") %>% filter(DAY==315) %>% ggplot(aes(x=W_SCALAR, y=SOMp)) + geom_point(aes(color=SITE), size=4) +theme_bw() +geom_smooth()
#BAGS_out_AllSites_Cal %>% mutate(SITE.SM = paste(SITE, SM_Type)) %>% inner_join(SM_comp, by="SITE.SM") %>% filter(DAY==315) %>% ggplot(aes(x=W_SCALAR, y=SOMp)) + geom_point(aes(color=SITE), size=4) +theme_bw() +geom_smooth()
#SP_SM.SOMp <- BAGS_out_AllSites_SP %>% mutate(SITE.SM = paste(SITE, SM_Type)) %>% inner_join(SM_comp, by="SITE.SM") %>% filter(DAY==315) 
#Cal_SM.SOMp <- BAGS_out_AllSites_Cal %>% mutate(SITE.SM = paste(SITE, SM_Type)) %>% inner_join(SM_comp, by="SITE.SM") %>% filter(DAY==315) 
#cor.test(SP_SM.SOMp$W_SCALAR, SP_SM.SOMp$SOMp)
#cor.test(Cal_SM.SOMp$W_SCALAR, Cal_SM.SOMp$SOMp) #only slightly higher R

#comparison of baseline MIMICS to LML
#got rid of LITi calcs since all litters back to the same starting value
#LIT_init <- BAGS_out_AllSites_DI %>% filter(DAY == 10) %>% mutate(LITi = LITBAGm+LITBAGs) %>% 
#  mutate(SITE.LT = paste(SITE, Litter_Type, sep=".")) %>% select(SITE.LT, LITi)
#boxplot(LIT_init$LITi)
#BAGS_out_AllSites_ES_vMOD_noTk$ID <- as.factor(rep(1:12, each=68985)) #not lining up with numerically with what I would expect - not sure what's going on.... ES one worked so even weirder!
# #try counting each site and see if each site is repeated the same number of times!
#BO_ES_1 <- readRDS("Analysis/MC_Output/BAGS_out_AllSites_ES_5000_NewInputs_10_1.rds") %>% select(-4, -5, -6, -7) %>% mutate(ID=rep(1:100, each=68985))
#BO_ES_2 <- readRDS("Analysis/MC_Output/BAGS_out_AllSites_ES_5000_NewInputs_10_2.rds") %>% select(-4, -5, -6, -7) %>%  mutate(ID=rep(101:190, each=68985))
#BO_ES_3 <- readRDS("Analysis/MC_Output/BAGS_out_AllSites_ES_5000_NewInputs_3.rds") %>% select(-4, -5, -6, -7) %>%  mutate(ID=rep(201:300, each=68985))
#BO_ES_4 <- readRDS("Analysis/MC_Output/BAGS_out_AllSites_ES_5000_NewInputs_4.rds") %>% select(-4, -5, -6, -7) %>%  mutate(ID=rep(301:388, each=68985))
#BAGS_out_AllSites_ES <- rbind(BO_ES_1, BO_ES_2)
#BAGS_out_AllSites_ES_best <- BAGS_out_AllSites_ES %>% filter(ID==241 | ID== 282 | ID== 48)
BAGS_out_AllSites_Cal <- BAGS_out_AllSites_Cal %>% mutate(ID=as.factor(rep(c(175, 176, 190), each=68985)))
#BAGS_out_AllSites_ES_NoCUE2 <- BAGS_out_AllSites_ES_NoCUE[ , -c(4:7)] #removing duplicate column names
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
BO_plot_sum <- BAGS_out_plot %>% group_by(SITE,DAY) %>% summarise(mean=mean(LIT_PerLoss), min=min(LIT_PerLoss), max=max(LIT_PerLoss)) #,SM_Type
# filter(ID==215 | ID== 71 | ID==227 | ID==38 | ID==164 | ID==21) %>%
BO_plot_sum <- BO_plot_sum %>% mutate(SITE=factor(SITE, levels=c("TREE", "BART", "HARV", "GRSM", "SERC", "TALL", "LENO"))) #MAT order
LML_sum2 <- LML_sum2 %>% mutate(site=factor(site, levels=c("TREE", "BART", "HARV", "GRSM", "SERC", "TALL", "LENO"))) #MAT order
tiff("MSBio_Fig2_SP.tiff", units="px", width=2000, height=1500, res=300)
ggplot() +
  geom_ribbon(data=BO_plot_sum, aes(y=100-mean, x=DAY-315, ymin = 100-min, ymax=100-max, group=SITE, fill=SITE, color=SITE), alpha = 0.3, size=0.5) +
  #geom_line(data=BO_plot_sum, aes(y=100-mean, x=DAY-315, group=SITE, color=SITE), linewidth=2, alpha = 0.3) +
  geom_point(data=LML_sum2, aes(y=100-mean.ML, x=doy, group=site, color=site), size = 3) +
  geom_errorbar(data=LML_sum2, aes(y=100-mean.ML, x=doy, ymin = 100-lci.ML, ymax = 100-uci.ML, group=site, color=site), width=0,linewidth=1) +
  ylab("Litter Bag C Remaining (%)") +
  xlab("Day") +
  xlim(0, 780) +
  scale_color_manual(values = colorBlind7) +
  scale_fill_manual(values = colorBlind7) +
  theme_bw(base_size = 20) #+
  #facet_wrap(.~SM_Type)
dev.off()
#seperate plotting of LITm and LITs
LIT_init <- BAGS_out_AllSites_SP %>% filter(DAY == 315) %>% mutate(LITm.i = LITBAGm) %>% mutate(LITs.i = LITBAGs) %>% 
  mutate(SITE.LT.SM = paste(SITE, Litter_Type, SM_Type, sep=".")) %>% select(SITE.LT.SM, LITm.i, LITs.i)
BAGS_out_plot <- BAGS_out_AllSites_Cal %>% mutate(SITE.LT.SM = paste(SITE, Litter_Type, SM_Type, sep=".")) %>% inner_join(LIT_init, by="SITE.LT.SM") %>%
  mutate(LITm_PerLoss = ((LITm.i - (LITBAGm))/LITm.i)*100) %>% mutate(LITs_PerLoss = ((LITs.i - (LITBAGs))/LITs.i)*100) #%>% filter(ID==190)
#plotting
BO_plot_sum.ms <- BAGS_out_plot %>% group_by(SITE,DAY) %>% summarise(mean.m=mean(LITm_PerLoss), min.m=min(LITm_PerLoss), max.m=max(LITm_PerLoss), 
                                                                             mean.s=mean(LITs_PerLoss), min.s=min(LITs_PerLoss), max.s=max(LITs_PerLoss)) %>%
  pivot_longer(3:8, names_to = 'LIT.ms', values_to = 'value') %>% ungroup() %>% mutate(LIT_type = ifelse(grepl("s", LIT.ms), "S", "M")) %>% 
  mutate(stat_type = rep(c("mean", "min", "max"), 15330)) %>% select(SITE, DAY, LIT_type, stat_type, value) %>% pivot_wider(names_from = stat_type, values_from = value)
#pivot wider based on new column that makes a structural column if LIT.ms contains "s"
BO_plot_sum.ms <- BO_plot_sum.ms %>% mutate(SITE=factor(SITE, levels=c("TREE", "BART", "HARV", "GRSM", "SERC", "TALL", "LENO"))) #MAT order
ggplot() +
  geom_ribbon(data=BO_plot_sum.ms, aes(y=100-mean, x=DAY-315, ymin = 100-min, ymax=100-max, group=SITE, fill=SITE), alpha = 0.3) +
  #geom_line(data=BO_plot_sum, aes(y=100-mean, x=DAY-315, group=SITE, color=SITE), linewidth=2, alpha = 0.3) +
  geom_point(data=LML_sum2, aes(y=100-mean.ML, x=doy, group=site, color=site), size = 3) +
  geom_errorbar(data=LML_sum2, aes(y=100-mean.ML, x=doy, ymin = 100-lci.ML, ymax = 100-uci.ML, group=site, color=site), width=0,linewidth=1) +
  ylab("Litter Bag C Remaining (%)") +
  xlab("Day") +
  xlim(0, 780) +
  facet_wrap(.~LIT_type) +
  theme_bw(base_size = 20) + scale_fill_manual(values=colorBlind7)+ scale_color_manual(values=colorBlind7)


#looking at model output to see how calibrated differs from starting point model
#regressing SP vs calibrated
BAGS_out_AllSites_ES176 <- filter(BAGS_out_AllSites_Cal, ID==176) %>% select(-ID) #missing decomp rates! Need these for analyzing differences 
BAGS_out_AllSites_ES175 <- filter(BAGS_out_AllSites_Cal, ID==175) %>% select(-ID)
BAGS_out_AllSites_ES190 <- filter(BAGS_out_AllSites_Cal, ID==190) %>% select(-ID)
BO_initial_SP <- BAGS_out_AllSites_SP %>% filter(DAY==315) %>% mutate(MICrK = MICr/MICk, LITms = LITm/LITs) %>% pivot_longer(5:19, names_to = 'Pools', values_to = 'Carbon_SP') %>% 
  mutate(ID=1:length(Carbon_SP)) %>% select(Carbon_SP, ID)
BO_initial_CE <- BAGS_out_AllSites_ES175 %>% filter(DAY==315) %>% mutate(MICrK = MICr/MICk, LITms = LITm/LITs) %>% pivot_longer(5:19, names_to = 'Pools', values_to = 'Carbon_CE') %>% mutate(ID=1:length(Carbon_CE))
BAGS2 <- BAGS %>% mutate(SITE.LQ = paste(SITE, TYPE, sep=".")) %>% select(SITE.LQ, CALC_MET)
BO_initial175 <- inner_join(BO_initial_SP, BO_initial_CE, by='ID') %>% mutate(SITE.LQ = paste(SITE, Litter_Type, sep=".")) %>% inner_join(BAGS2, by="SITE.LQ")
low_pools <- c("LITm", "LITms", "MICk", "MICr")
BO_initial %>% filter(Pools %in% low_pools) %>% ggplot(aes(x=Carbon_SP, y=Carbon_CE)) + geom_point(aes(color=CALC_MET, shape=SM_Type), alpha=0.5, size=3) + facet_grid(.~Pools) + 
  geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16)
high_pools <- c("LITs", "SOMa", "SOMc", "SOMp") #"MICrK", 
BO_initial %>% filter(Pools %in% high_pools) %>% ggplot(aes(x=Carbon_SP, y=Carbon_CE)) + geom_point(aes(color=CALC_MET, shape=SM_Type), alpha=0.5, size=3) + facet_grid(.~Pools) + 
  geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16)
MICrK_pool <- c("MICrK") 
BO_initial190 %>% filter(Pools %in% MICrK_pool) %>% ggplot(aes(x=Carbon_SP, y=Carbon_CE)) + geom_point(aes(color=CALC_MET, shape=SM_Type), alpha=0.5, size=3) + facet_grid(.~Pools) + 
  geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) +ylim(0,3.6) +xlim(0,1.5)
#combining all into one plot
BO_initial_sum <- rbind(BO_initial175, BO_initial176, BO_initial190) %>% group_by(SITE, Litter_Type, CALC_MET, SM_Type, Pools) %>% 
  summarise(sp.avg =mean(Carbon_SP), cal.avg =mean(Carbon_CE), n=n(),
  lci.sp = sp.avg - qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_SP)/sqrt(n)), uci.sp = sp.avg + qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_SP)/sqrt(n)),
  lci.cal = cal.avg - qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_CE)/sqrt(n)), uci.cal = cal.avg + qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_CE)/sqrt(n)))
LITm_pool <- c("LITm")
tiff("MSBio_Fig4_LITm.tiff", units="px", width=950, height=1150, res=300)
BO_initial_sum %>% filter(Pools %in% LITm_pool) %>% ggplot() + geom_point(aes(x=sp.avg*1000, y=cal.avg*1000, color=CALC_MET, shape=SM_Type), alpha=0.5, size=3) + 
  ylab(expression(paste("Steady state C (g C m"^"-2"*")"))) + xlab(expression(paste("Steady state C (g C m"^"-2"*")"))) +
  facet_grid(.~Pools) + geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) + theme(legend.position="none")
  #geom_errorbar(aes(x=sp.avg, y=cal.avg, ymin=lci.cal, ymax=uci.cal)) +
dev.off()
high_pools <- c("LITs", "SOMa", "SOMc")
BO_initial_sum %>% filter(Pools %in% high_pools) %>% ggplot() + geom_point(aes(x=sp.avg, y=cal.avg, color=CALC_MET, shape=SM_Type), alpha=0.5, size=3) + 
  ylab(expression(paste("Calibrated model C (kg C m"^"-2"*")"))) + xlab(expression(paste("Default model C (kg C m"^"-2"*")"))) +
  facet_grid(.~Pools) + geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16)
MICrK_pool <- c("MICrK")
tiff("MSBio_Fig4_MICrK.tiff", units="px", width=850, height=1150, res=300)
BO_initial_sum %>% filter(Pools %in% MICrK_pool) %>% ggplot() + geom_point(aes(x=sp.avg, y=cal.avg, color=CALC_MET, shape=SM_Type), alpha=0.5, size=3) + 
  ylab("Calibrated model") + xlab("Default model") + xlab("Copiotroph:oligotroph") + ylab("Copiotroph:oligotroph") +
  facet_grid(.~Pools) + geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) + theme(legend.position="none")
dev.off()
#other pools
MIC_pools <- c("MICr", "MICk")
BO_initial_sum %>% filter(Pools %in% MIC_pools) %>% ggplot() + geom_point(aes(x=sp.avg, y=cal.avg, color=SITE, shape=SM_Type), alpha=0.5, size=3) + 
  ylab("Calibrated model") + xlab("Default model") + xlab("Copiotroph or oligotroph C") + ylab("Copiotroph or oligotroph C") +
  facet_grid(.~Pools) + geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16)
#decomposition rates
#scaled by microbial biomass or initial litter - checked and the same for any calibration
LB.i <- BAGS_out_AllSites_SP %>% filter(DAY==315) %>% mutate(SITE.LQ = paste(SITE, Litter_Type, sep=".")) %>% group_by(SITE.LQ) %>% 
                                                               summarise(LBm.i = mean(LITBAGm), LBs.i = mean(LITBAGs))
BO_Decomp_SP <- BAGS_out_AllSites_SP %>% mutate(SITE.LQ = paste(SITE, Litter_Type, sep=".")) %>% filter(DAY>315) %>%#inner_join(LB.i, by="SITE.LQ") %>% filter(DAY>315) %>% 
  group_by(SITE,Litter_Type, SM_Type) %>%  #inner_join(BAGS2, by="SITE.LQ") %>% 
  summarise(Drm = mean(Decomp_rate_rm/MICr), Drs = mean(Decomp_rate_rs/MICr), 
  Dkm = mean(Decomp_rate_km/MICk), Dks = mean(Decomp_rate_ks/MICk)) %>% 
  pivot_longer(4:7, names_to = 'Pools', values_to = 'Carbon_SP') %>% #ggplot(aes(x=DAY, y=Carbon_SP)) + geom_point(aes(color=CALC_MET)) + facet_grid(.~Pools)
  mutate(combo = paste(SITE, Litter_Type, SM_Type, Pools, sep=".")) %>% select(Carbon_SP, combo)
BO_Decomp_CE <- BAGS_out_AllSites_ES175 %>% filter(DAY>315) %>% #mutate(SITE.LQ = paste(SITE, Litter_Type, sep=".")) %>% inner_join(LB.i, by="SITE.LQ") %>% 
  group_by(SITE, Litter_Type, SM_Type) %>%
  summarise(Drm = mean(Decomp_rate_rm/MICr), Drs = mean(Decomp_rate_rs/MICr), Dkm = mean(Decomp_rate_km/MICk), Dks = mean(Decomp_rate_ks/MICk)) %>% 
  pivot_longer(4:7, names_to = 'Pools', values_to = 'Carbon_CE') %>% 
  mutate(combo = paste(SITE, Litter_Type, SM_Type, Pools, sep="."))
BAGS2 <- BAGS %>% mutate(SITE.LQ = paste(SITE, TYPE, sep=".")) %>% select(SITE.LQ, CALC_MET)
BO_decomp175 <- inner_join(BO_Decomp_SP, BO_Decomp_CE, by='combo') %>% mutate(SITE.LQ = paste(SITE.x, Litter_Type.x, sep=".")) %>% inner_join(BAGS2, by="SITE.LQ")
BO_decomp190 %>% filter(Pools=="Dks") %>% ggplot(aes(x=Carbon_SP, y=Carbon_CE)) + geom_point(aes(color=CALC_MET, shape=SM_Type), alpha=0.5, size=3) + facet_grid(.~Pools) + 
  geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) #filter(Pools=="Dkm") %>%
#means of BOs
BO_decomp_sum <- rbind(BO_decomp175, BO_decomp176, BO_decomp190) %>% group_by(SITE.x, Litter_Type.x, CALC_MET, SM_Type, Pools) %>% 
  summarise(sp.avg =mean(Carbon_SP), cal.avg =mean(Carbon_CE), n=n(),
            lci.sp = sp.avg - qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_SP)/sqrt(n)), uci.sp = sp.avg + qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_SP)/sqrt(n)),
            lci.cal = cal.avg - qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_CE)/sqrt(n)), uci.cal = cal.avg + qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_CE)/sqrt(n)))
BO_decomp_sum %>% filter(Pools=="Drs") %>% ggplot(aes(x=sp.avg, y=cal.avg)) + geom_point(aes(color=CALC_MET, shape=SM_Type), alpha=0.5, size=3) + facet_grid(.~Pools) + 
  geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) + theme(legend.position="none") +
  xlab(expression(paste("Specific decomposition (hr"^"-1"*")"))) + ylab(expression(paste("Specific decomposition (hr"^"-1"*")"))) #filter(Pools=="Dkm") %>%
#overall mean decomp rate
BO_decomp %>% group_by(SITE.x,Litter_Type.x, SM_Type) %>% summarise(decomp.SP = mean(Carbon_SP), decomp.CE = mean(Carbon_CE), fMET=mean(CALC_MET)) %>%
  ggplot(aes(x=decomp.SP, y=decomp.CE)) + geom_point(aes(color=fMET, shape=SM_Type), alpha=0.5, size=3) + 
  geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) 
#decomp rate not relativized
BO_Decomp_SP2 <- BAGS_out_AllSites_SP %>% mutate(SITE.LQ = paste(SITE, Litter_Type, sep=".")) %>% inner_join(LB.i, by="SITE.LQ") %>% filter(DAY>315) %>% 
  group_by(SITE,Litter_Type, SM_Type) %>% summarise(Drm = mean(Decomp_rate_rm), Drs = mean(Decomp_rate_rs), 
            Dkm = mean(Decomp_rate_km), Dks = mean(Decomp_rate_ks)) %>% 
  pivot_longer(4:7, names_to = 'Pools', values_to = 'Carbon_SP') %>% 
  mutate(combo = paste(SITE, Litter_Type, SM_Type, Pools, sep=".")) %>% select(Carbon_SP, combo)
BO_Decomp_CE2 <- BAGS_out_AllSites_ES190 %>% filter(DAY>315) %>% mutate(SITE.LQ = paste(SITE, Litter_Type, sep=".")) %>% inner_join(LB.i, by="SITE.LQ") %>% group_by(SITE, Litter_Type, SM_Type) %>%
  summarise(Drm = mean(Decomp_rate_rm), Drs = mean(Decomp_rate_rs), Dkm = mean(Decomp_rate_km), Dks = mean(Decomp_rate_ks)) %>% 
  pivot_longer(4:7, names_to = 'Pools', values_to = 'Carbon_CE') %>% 
  mutate(combo = paste(SITE, Litter_Type, SM_Type, Pools, sep="."))
BO_decomp190.2 <- inner_join(BO_Decomp_SP2, BO_Decomp_CE2, by='combo') %>% mutate(SITE.LQ = paste(SITE.x, Litter_Type.x, sep=".")) %>% inner_join(BAGS2, by="SITE.LQ")
BO_decomp2 %>% group_by(SITE.x,Litter_Type.x, SM_Type) %>% summarise(decomp.SP = mean(Carbon_SP), decomp.CE = mean(Carbon_CE), fMET=mean(CALC_MET)) %>%
  ggplot(aes(x=decomp.SP, y=decomp.CE)) + geom_point(aes(color=fMET), alpha=0.5, size=3) + #facet_grid(.~Pools) + 
  geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) #filter(Pools=="Dkm") %>%
#means of BO
BO_decomp_sum2 <- rbind(BO_decomp175.2, BO_decomp176.2, BO_decomp190.2) %>% group_by(SITE.x, Litter_Type.x, CALC_MET, SM_Type, Pools) %>% 
  summarise(sp.avg =mean(Carbon_SP), cal.avg =mean(Carbon_CE), n=n(),
            lci.sp = sp.avg - qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_SP)/sqrt(n)), uci.sp = sp.avg + qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_SP)/sqrt(n)),
            lci.cal = cal.avg - qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_CE)/sqrt(n)), uci.cal = cal.avg + qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_CE)/sqrt(n)))
Decomp_pools <- c("Drs", "Dkm") # "Drs", "Dkm", "Drm", "Dks"
tiff("MSBio_Fig4_decomp.tiff", units="px", width=1650, height=1150, res=300)
BO_decomp_sum2 %>% filter(Pools %in% Decomp_pools) %>% ggplot(aes(x=sp.avg*24000, y=cal.avg*24000)) + geom_point(aes(color=CALC_MET, shape=SM_Type), alpha=0.5, size=3) + facet_grid(.~Pools) + 
  geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) + theme(legend.position="none") +
  xlab(expression(paste("Decomposition flux (g C m"^"-2"*" d"^"-1"*")"))) + ylab(expression(paste("Decomposition flux (g C m"^"-2"*" d"^"-1"*")"))) #filter(Pools=="Dkm") %>%
dev.off()
#litterbag decomposition - work n this!!
BO_Decomp_SP3 <- BAGS_out_AllSites_SP %>% mutate(SITE.LQ = paste(SITE, Litter_Type, sep=".")) %>% inner_join(LB.i, by="SITE.LQ") %>% filter(DAY>315) %>% 
  group_by(SITE,Litter_Type, SM_Type) %>% summarise(LBm = mean(as.numeric(LITBAGm)), LBs = mean(as.numeric(LITBAGs))) %>% 
  pivot_longer(4:5, names_to = 'Pools', values_to = 'Carbon_SP') %>% 
  mutate(combo = paste(SITE, Litter_Type, SM_Type, Pools, sep=".")) %>% select(Carbon_SP, combo)
BO_Decomp_CE3 <- BAGS_out_AllSites_ES175 %>% filter(DAY>315) %>% mutate(SITE.LQ = paste(SITE, Litter_Type, sep=".")) %>% inner_join(LB.i, by="SITE.LQ") %>% group_by(SITE, Litter_Type, SM_Type) %>%
  summarise(LBm = mean(as.numeric(LITBAGm)), LBs = mean(as.numeric(LITBAGs))) %>% 
  pivot_longer(4:5, names_to = 'Pools', values_to = 'Carbon_CE') %>% 
  mutate(combo = paste(SITE, Litter_Type, SM_Type, Pools, sep="."))
BO_decomp175.3 <- inner_join(BO_Decomp_SP3, BO_Decomp_CE3, by='combo') %>% mutate(SITE.LQ = paste(SITE.x, Litter_Type.x, sep=".")) %>% inner_join(BAGS2, by="SITE.LQ")
#means of BO
BO_decomp_sum3 <- rbind(BO_decomp175.3, BO_decomp176.3, BO_decomp190.3) %>% group_by(SITE.x, Litter_Type.x, CALC_MET, SM_Type, Pools) %>% 
  summarise(sp.avg =mean(Carbon_SP), cal.avg =mean(Carbon_CE), n=n(),
            lci.sp = sp.avg - qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_SP)/sqrt(n)), uci.sp = sp.avg + qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_SP)/sqrt(n)),
            lci.cal = cal.avg - qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_CE)/sqrt(n)), uci.cal = cal.avg + qt(1 - ((1 - 0.95) / 2), n - 1) * (sd(Carbon_CE)/sqrt(n)))
LB_pools <- c("LBm", "LBs") # "Drs", "Dkm", "Drm", "Dks"
BO_decomp_sum3 %>% filter(Pools %in% LB_pools) %>% ggplot(aes(x=sp.avg, y=cal.avg)) + geom_point(aes(color=CALC_MET, shape=SM_Type), alpha=0.5, size=3) + facet_grid(.~Pools) + 
  geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) + theme(legend.position="none") +
  xlab(expression(paste("Litterbag C (kg C m"^"-2"*")"))) + ylab(expression(paste("Litterbag C (kg C m"^"-2"*")"))) #filter(Pools=="Dkm") %>%
#relationships between variables 
BAGS_out_AllSites_ES190 %>% filter(DAY==315) %>% mutate(SITE.LQ = paste(SITE, Litter_Type, sep=".")) %>% inner_join(BAGS2, by="SITE.LQ") %>%
  ggplot(aes(y=MICr/MICk, x=CALC_MET)) + geom_point(aes(color=LITm/LITs, shape=SM_Type), alpha=0.5, size=3) + ylim(0,3) #+ xlim(0.05,0.25)#color=CALC_MET, 
MSBio_SM <- MSBio_sites_SM %>% mutate(SITE.SM = paste(SITE, SM_type, sep="."))
BAGS_out_AllSites_ES190 %>% filter(DAY==315) %>% mutate(SITE.LQ = paste(SITE, Litter_Type, sep=".")) %>%  inner_join(BAGS2, by="SITE.LQ") %>% mutate(SITE.SM = paste(SITE, SM_Type, sep=".")) %>% 
  inner_join(MSBio_SM, by="SITE.SM") %>% mutate(fct=cut(CALC_MET, breaks=c(0,0.4,0.5, 0.75), labels=c("low fMET","mid fMET","high fMET"))) %>% 
  ggplot(aes(y=MICr/MICk, x=TSOI)) + geom_point(aes(color=W_SCALAR2), alpha=0.5, size=3) +ylim(0,3) + facet_grid(.~fct) #aes(color=CALC_MET, shape=SM_Type), 
BO_i_rK <- BO_initial %>% filter(Carbon_CE >0) %>% mutate(SITE.SM = paste(SITE, SM_Type, sep=".")) %>% inner_join(MSBio_SM, by="SITE.SM") %>% mutate(SITE = SITE.x) %>% filter(Pools == "MICrK")
BOi_CE <- lmer(Carbon_CE ~ CALC_MET*TSOI*W_SCALAR2 + (1|SITE), data=BO_i_rK)
summary(BOi_CE)
#Looking at the spread of an individual pool
#LITms
# BAGS_out_AllSites_SP %>% filter(DAY==315) %>% mutate(LITms = LITm/LITs) %>%
#   mutate(SITE.LQ = paste(SITE, Litter_Type, sep=".")) %>% inner_join(BAGS2, by="SITE.LQ") %>%
#   ggplot() + geom_point(aes(y=LITms, x=SITE, color=CALC_MET, shape=SM_Type), size=3, alpha=0.5) +
#   scale_x_discrete(guide = guide_axis(angle = 90)) + theme_bw(base_size = 16) #+ ylim(0.03,0.2)
# #MICrK
# BAGS_out_AllSites_SP%>% filter(DAY==315) %>% mutate(MICrk = MICr/MICk) %>%
#   mutate(SITE.LQ = paste(SITE, Litter_Type, sep=".")) %>% inner_join(BAGS2, by="SITE.LQ") %>%
#   ggplot() + geom_point(aes(y=MICrk, x=SITE, color=CALC_MET, shape=SM_Type), size=3, alpha=0.5) +
#   scale_x_discrete(guide = guide_axis(angle = 90)) + theme_bw(base_size = 16) #+ ylim(0.03,0.2)
#comparing seasonality of decomp
# temp_order = c('TREE', 'BART', 'HARV', 'GRSM', 'SERC', 'TALL', 'LENO')
# BAGS_out_AllSites_SP %>% filter(DAY > 314 & DAY < 679) %>%  group_by(SITE,DAY) %>% summarise(min.MICr=min(MICr), max.MICr=max(MICr),min.MICk=min(MICk), max.MICk=max(MICk),
#                                                                               min.LBm=min(LITBAGm), max.LBm=max(LITBAGm),min.LBs=min(LITBAGs), max.LBs=max(LITBAGs)) %>%
#   mutate(SITE = factor(SITE, levels=temp_order)) %>%
#   ggplot() + geom_ribbon(aes(ymin=min.MICr,ymax=max.MICr, x=DAY, fill ="MICr"), alpha=0.7) +
#   geom_ribbon(aes(ymin=min.MICk, ymax=max.MICk, x=DAY, fill ="MICk"), alpha=0.7) + geom_ribbon(aes(ymin=min.LBm, ymax=max.LBm, x=DAY, fill ="LITBAGm"), alpha=0.7) +
#   geom_ribbon(aes(ymin=min.LBs, ymax=max.LBs, x=DAY, fill ="LITBAGs"), alpha=0.7) + facet_grid(.~SITE) + ylim(0, 0.3)
# #variation in relationship between r:K and m:s over time
# BAGS_out_AllSites_CalES.13 %>% mutate(MICrK=MICr/MICk) %>% mutate(LITms = LITm/LITs) %>% group_by(DAY) %>% summarise(slope = lm(LITms ~ MICrK)$coefficients[2]) %>%
#   ggplot(aes(x=DAY, y=slope)) +geom_point() +theme_bw(base_size = 16)
# BO.MicRange.SP <- BAGS_out_AllSites_SP %>% filter(DAY > 314 & DAY < 679) %>%  group_by(SITE) %>% summarise(min.MICr=min(MICr), max.MICr=max(MICr),
#                                                                                             min.MICk=min(MICk), max.MICk=max(MICk))%>%
#   mutate(SITE = factor(SITE, levels=temp_order)) %>% mutate(MICr.range.SP = max.MICr-min.MICr, MICk.range.SP = max.MICk-min.MICk) %>%
#   select(SITE, MICr.range.SP, MICk.range.SP)
# BO.MicRange <-BAGS_out_AllSites_CalES %>% filter(DAY > 314 & DAY < 679) %>%  group_by(SITE) %>% summarise(min.MICr=min(MICr), max.MICr=max(MICr),
#                                                                                             min.MICk=min(MICk), max.MICk=max(MICk))%>%
#   mutate(SITE = factor(SITE, levels=temp_order)) %>% mutate(MICr.range.CE = max.MICr-min.MICr, MICk.range.CE = max.MICk-min.MICk) %>%
#   select(SITE, MICr.range.CE, MICk.range.CE) %>% inner_join(BO.MicRange.SP, by="SITE")
# ggplot(BO.MicRange) + geom_point(aes(y=MICr.range.SP, x=SITE, color="MICr", shape="SP"), size=3, alpha=0.7) + geom_point(aes(y=MICk.range.SP, x=SITE, color="MICk", shape="SP"), size=3, alpha=0.7) +
#   geom_point(aes(y=MICr.range.CE, x=SITE, color="MICr", shape="CE"), size=3, alpha=0.7) + geom_point(aes(y=MICk.range.CE, x=SITE, color="MICk", shape="CE"), size=3, alpha=0.7) +
#   ylab('Mic range') + theme_bw(base_size = 16)
#comparing decomposition rates between model runs for more interpreatble output
# BO_DecompRate_SP  <- BAGS_out_AllSites_SP %>% filter(DAY == 315 | DAY== 679) %>%  mutate(DAY.TYPE = ifelse(DAY==315, "start", "end")) %>% select(SITE, Litter_Type, SM_Type, DAY.TYPE, LITBAGs) %>%
#   pivot_wider(names_from = DAY.TYPE, values_from = LITBAGs) %>% mutate(LITs.rate_SP = (end-start)/315) %>% mutate(ID = paste(SITE, Litter_Type, SM_Type, sep=".")) %>% 
#   select(ID, LITs.rate_SP) #315, 679
# BO_DecompRate <- BAGS_out_AllSites_CalES.noTk.9 %>% filter(DAY == 315 | DAY== 679) %>%  mutate(DAY.TYPE = ifelse(DAY==315, "start", "end")) %>% select(SITE, Litter_Type, SM_Type, DAY.TYPE, LITBAGs) %>%
#   pivot_wider(names_from = DAY.TYPE, values_from = LITBAGs) %>% mutate(LITs.rate_Cal = (end-start)/315) %>% mutate(ID = paste(SITE, Litter_Type, SM_Type, sep=".")) %>% 
#   select(ID, SITE, Litter_Type, SM_Type, LITs.rate_Cal) %>% inner_join(BO_DecompRate_SP, by='ID') %>% 
#   mutate(SITE.LQ = paste(SITE, Litter_Type, sep=".")) %>% inner_join(BAGS2, by="SITE.LQ")
# site_shapes <- c(16, 15, 17, 18, 8, 3, 10)
# ggplot(BO_DecompRate, aes(x=LITs.rate_SP, y=LITs.rate_Cal)) + geom_point(aes(color=CALC_MET, shape=SITE), alpha=0.7, size=3) + 
#   geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) + scale_shape_manual(values=site_shapes)
# BO_initial.rk <- BO_initial %>% filter(Pools == "MICrK") %>% mutate(ID = as.character(paste(SITE, Litter_Type, SM_Type, sep=".")))
# BO_MB_SP <- BAGS_out_AllSites_SP %>% filter(DAY==315) %>% mutate(MIC.SP = MICr+MICk) %>% mutate(ID = paste(SITE, Litter_Type, SM_Type, sep = ".")) %>% select(ID, MIC.SP)
# BO_MB <- BAGS_out_AllSites_CalES.noTk.9 %>% filter(DAY==315) %>% mutate(MIC.CE = MICr+MICk) %>% mutate(ID = paste(SITE, Litter_Type, SM_Type, sep = ".")) %>% inner_join(BO_MB_SP, by="ID")  %>% select(ID, MIC.SP, MIC.CE)
# BO_DecompRate %>% inner_join(BO_MB, by="ID") %>% ggplot() + geom_point(aes(x=MIC.CE, y=LITs.rate_Cal, color=CALC_MET), size=3, alpha=0.7) + 
#   theme_bw(base_size = 16) + ylim(-0.0002, 0) +  xlim(0.05,0.25) #+ xlab("Initial MICr:MICK")

# #formating data for modelVobs, RWA, and effect size
#need to standardize and log this data before running it thru again
#MFG_stdzd <- MFG_analysis %>% mutate_at(c('C_O', 'LIG_N', 'log.vwc'), ~(scale(.) %>% as.vector))
FieldData <- LML_sum2 %>% mutate(DAY=doy, SITE=site) %>% mutate(SITE.DAY=paste(SITE, DAY, sep=".")) %>% select(time.point,SITE.DAY, mean.ML, lci.ML, uci.ML)
#LIT_init <- BAGS_out_AllSites_DI %>% filter(DAY == 10) %>% mutate(LITi = LITBAGm+LITBAGs) %>% select(SITE, LITi)
#boxplot(LIT_init$LITi)
#df <- BAGS_out_AllSites_DI %>% left_join(LIT_init, by = "SITE")
#BAGS_out_AllSites_ES_noTk$ID <- as.factor(rep(1:12, each=68985))
LITi = 0.1
df_LML <- BAGS_out_AllSites_Cal %>% mutate(DAY.LitOut = DAY -314) %>% mutate(SITE.DAY=paste(SITE, DAY.LitOut, sep=".")) %>% 
  right_join(FieldData, by="SITE.DAY") %>% mutate(LIT_PerLoss = ((LITi - (LITBAGm+LITBAGs))/LITi)*100)
#write.csv(df_LML, 'df_LML_ES_vMOD_noTk_ParamUncert.csv')
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
#loop for muliple Psets
#df_LML2 <- df_LML %>% inner_join(ES_Psets, by="ID")
#rn_filter = c(1:4500)
#df_LML2 <- df_LML2 %>% filter(run_num %in% rn_filter)
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
df_LML_best <- filter(df_LML, ID== 176 | ID==175 | ID==190| ID==108| ID==21) #ID== 176 | ID==175 | ID==190| 
#summarized values
df_LML_sum <- df_LML %>% group_by(SITE, time.point) %>% summarise(mean.ML = mean(mean.ML), m.uci.ML = mean(uci.ML), m.lci.ML = mean(lci.ML), mean.LPL = mean(LIT_PerLoss),
                                                                       n=n(), SE = sd(LIT_PerLoss)/sqrt(n),
                                                                       min.LPL = min(LIT_PerLoss),
                                                                       max.LPL = max(LIT_PerLoss),
                                                                       lci.LPL = mean.LPL - qt(1 - ((1 - 0.95) / 2), n - 1) * SE,
                                                                       uci.LPL = mean.LPL + qt(1 - ((1 - 0.95) / 2), n - 1) * SE)
modelVobs <- lm(mean.LPL~mean.ML, data=df_LML_sum)
summary(modelVobs)
sqrt(mean((df_LML_sum$mean.ML - df_LML_sum$mean.LPL)^2)) #RMSE
(1/length(df_LML_sum$n))*sum(df_LML_sum$mean.ML - df_LML_sum$mean.LPL) #bias
df_LML_sum <- df_LML_sum %>% mutate(SITE=factor(SITE, levels=c("TREE", "BART", "HARV", "GRSM", "SERC", "TALL", "LENO"))) #MAT order
tiff("MSBio_Fig2_Cal_inset.tiff", units="px", width=1100, height=1030, res=300)
ggplot(df_LML_sum, aes(x=mean.ML, y=mean.LPL)) + geom_point(aes(color=SITE), size=4) + geom_smooth(method = "lm", color="black")  + xlim(0,80) + ylim(0,80) +
  geom_errorbar(aes(ymin=lci.LPL, ymax=uci.LPL, color=SITE), size=1) + geom_errorbarh(aes(xmin=m.lci.ML, xmax=m.uci.ML, color=SITE), size=1) +
  xlab("Observed litter percent C loss") + ylab("Modeled litter percent C loss") + geom_abline(intercept=0, slope=1, linetype=2) + theme_bw(base_size = 16) + 
  scale_color_manual(values = colorBlind7) + theme(legend.position="none")
dev.off()

#checking paramter relationships of top Psets
ES_fit_best <- ES_fit %>% filter(RMSE>5.3 & RMSE < 5.8) %>% filter(i != 48) %>% filter(i != 107) %>% mutate(ID=i) %>% inner_join(ES_Psets_all, by='ID')
plot(ES_fit_best$beta_x, ES_fit_best$Tau_r)



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
  summarise(W_SCALAR_mean=mean(W_SCALAR)) %>% select(SITE.SM, W_SCALAR_mean)
#initial MICrK
#BAGS_out_AllSites_ES$ID <- rep(1:18, each = 22995)
#ID == 68 | ID==16 | ID== 64 | ID==70 | ID==52
MIC_init <- BAGS_out_AllSites_SP %>% filter(DAY == 315) %>%mutate(MICrK.i =  MICr/MICk)%>% mutate(MICr.i = MICr)%>% mutate(MICK.i = MICk)%>%
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% select(SITE.SM.LQ, MICrK.i, MICr.i, MICK.i) #%>% filter(ID == 175)
#need bag means!
BAGS_LIGN <- MSBio_BAGS %>% mutate(SITE.LQ = paste(SITE, TYPE, sep = ".")) %>% select(SITE.LQ, BAG_LIG_N)
df_analysis <- df_LML %>% mutate(MICrK = MICr/MICk) %>% mutate(MIC=MICr+MICk) %>% mutate(SOC = SOMa+SOMc+SOMp) %>% 
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% mutate(SITE.SM = paste(SITE, SM_Type, sep = ".")) %>%
  mutate(SITE.LQ = paste(SITE, Litter_Type, sep = ".")) %>% inner_join(DI_means, by="SITE.SM") %>% inner_join(MIC_init, by="SITE.SM.LQ") %>% 
  inner_join(BAGS_LIGN, by="SITE.LQ") %>% mutate(LQ.SM=paste(Litter_Type, SM_Type, sep=".")) # %>% filter(ID == 175)
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
Obs_ES_mod <- lmer(LIT_PerLoss ~ scale(log(W_SCALAR_mean))+scale(BAG_LIG_N)+scale(MICrK.i)+ (1|SITE/LQ.SM), data = df_check) #MAT #scale(BAG_LIG_N)+
#Obs_ES_mod <- lmer(LIT_PerLoss ~ scale(log(W_SCALAR_mean))+scale(BAG_LIG_N)+scale(MICr.i)+scale(MICK.i)+ (1|SITE/LQ.SM), data = df_check) #MAT #scale(BAG_LIG_N)+
#Obs_ES_mod2 <- lmer(LIT_PerLoss ~ log_WS+BAG_LIG_N+MICrK.i+ (1|SITE), data = df_ES)
summary(Obs_ES_mod)
Obs_ES <- as.data.frame(fixef(Obs_ES_mod)) #fixed effects coefficients as effect size
#Obs_ES <- as.data.frame(Obs_ES_mod$coefficients) #if not using fixed effects model do this
Obs_ES$Vars <- rownames(Obs_ES)
colnames(Obs_ES)[1] <- "value"
Obs_ES <- Obs_ES[-1, ]
Obs_ES$mult <- ifelse(Obs_ES$value <0, -1, 1)
Obs_ES$rel_ES <- (abs(Obs_ES$value)/sum(abs(Obs_ES$value))) * 100 * Obs_ES$mult
#Obs_ES$Vars <- factor(Obs_ES$Vars, levels=c('scale(MICrK.i)','scale(BAG_LIG_N)', 'scale(log(W_SCALAR_mean))'),
#                       labels=c('MICrK.i','BAG_LIG_N', 'log_WS'))
ggplot(Obs_ES, aes(x=Vars, y=rel_ES)) + geom_bar(stat="identity", fill="darkolivegreen3") + coord_flip() + 
  geom_text(aes(label=round(rel_ES, digits=1), vjust=1.5), size=5) +theme_bw(base_size = 16)
Obs_ES_SP<- Obs_ES 
rbind(Obs_ES176, Obs_ES175, Obs_ES190) %>% group_by(Vars) %>% summarise(mean.rES = mean(rel_ES), sd.rES = sd(rel_ES)) %>%
  ggplot(aes(x=Vars, y=mean.rES)) + geom_bar(stat="identity", fill="darkolivegreen3", aes(group=Vars)) +  geom_errorbar(aes(ymax=mean.rES+sd.rES, ymin=mean.rES-sd.rES), width=0.1, size=1) + 
  coord_flip() + geom_text(aes(label=round(mean.rES, digits=1), vjust=1.5), size=5) +theme_bw(base_size = 16) #, hjust=1.1
Obs_ES_Cal <- rbind(Obs_ES176, Obs_ES175, Obs_ES190) %>% group_by(Vars) %>% summarise(mean.rES = mean(rel_ES), sd.rES = sd(rel_ES)) %>% mutate(ID="Cal")
Obs_ES_SP2 <- Obs_ES_SP %>% mutate(mean.rES=rel_ES, sd.rES=0, ID="SP") %>% select(Vars, mean.rES, sd.rES, ID)
#creating data frame for observed effect sizes
Vars= c('scale(log(W_SCALAR_mean))','scale(BAG_LIG_N)', 'scale(MICrK.i)')
mean.rES = c(42.8, -34.2, -22.9)
sd.rES = c(0,0,0)
ID= c("Obs", "Obs", "Obs")
Obs_ES_Obs <- data.frame(Vars, mean.rES, sd.rES, ID)
tiff("MSBio_Fig1.tiff", units="px", width=2000, height=1500, res=300)
rbind(Obs_ES_Obs, Obs_ES_SP2, Obs_ES_Cal) %>% ggplot() + 
  geom_bar(aes(x=factor(Vars, level=c('scale(log(W_SCALAR_mean))','scale(BAG_LIG_N)', 'scale(MICrK.i)'), labels = c("Soil moisture", "Lignin:N", "Copiotroph:oligotroph")),
        y=mean.rES, fill = factor(ID, level=c('Obs', 'SP', 'Cal'), labels = c("Observations", "Default", "Calibrated"))), stat="identity", position = position_dodge(), color="black", linewidth=1) +
  geom_errorbar(aes(x=factor(Vars, level=c('scale(log(W_SCALAR_mean))','scale(BAG_LIG_N)', 'scale(MICrK.i)'), labels = c("Soil moisture", "Lignin:N", "Copiotroph:oligotroph")),
                    y=mean.rES, ymax=mean.rES+sd.rES, ymin=mean.rES-sd.rES, 
                    group=factor(ID, level=c('Obs', 'SP', 'Cal'), labels = c("Observations", "Default", "Calibrated"))), width=0.1, size=1, position = position_dodge(0.9)) +
  xlab("") + ylab("Relative effect size") + geom_abline(intercept=0, slope=0, color="black", linewidth=0.6) +
  theme_bw(base_size = 16)  + scale_fill_manual(name="Type", values = c(Observations = "black", Default = "white", Calibrated="darkgrey")) 
dev.off()

#with points for calibrated
Obs_ES_Cal_points <- rbind(Obs_ES176, Obs_ES175, Obs_ES190) %>% mutate(ID="Cal") %>% select(Vars, rel_ES, ID) 
Obs_ES_SP_points <- Obs_ES_SP2 %>% mutate(rel_ES = c(NA,NA,NA)) %>% select(Vars, rel_ES, ID)
Obs_ES_Obs_points <- Obs_ES_Obs %>% mutate(rel_ES = c(NA,NA,NA)) %>% select(Vars, rel_ES, ID)
ES_points <- rbind(Obs_ES_Cal_points, Obs_ES_SP_points, Obs_ES_Obs_points)
Obs_ES_all <- rbind(Obs_ES_Obs, Obs_ES_SP2, Obs_ES_Cal) 
tiff("MSBio_Fig1.tiff", units="px", width=2000, height=1500, res=300)
ggplot() + 
  geom_bar(data=Obs_ES_all, aes(x=factor(Vars, level=c('scale(log(W_SCALAR_mean))','scale(BAG_LIG_N)', 'scale(MICrK.i)'), labels = c("Soil moisture", "Lignin:N", "Copiotroph:oligotroph")),
               y=mean.rES, fill = factor(ID, level=c('Obs', 'SP', 'Cal'), labels = c("Observations", "Default", "Calibrated"))), stat="identity", position = position_dodge(), color="black", linewidth=1) +
  #geom_errorbar(data=Obs_ES_all, aes(x=factor(Vars, level=c('scale(log(W_SCALAR_mean))','scale(BAG_LIG_N)', 'scale(MICrK.i)'), labels = c("Soil moisture", "Lignin:N", "Copiotroph:oligotroph")),
  #                  y=mean.rES, ymax=mean.rES+sd.rES, ymin=mean.rES-sd.rES, 
  #                  group=factor(ID, level=c('Obs', 'SP', 'Cal'), labels = c("Observations", "Default", "Calibrated"))), width=0.1, size=1, position = position_dodge(0.9)) +
  geom_point(data=ES_points, aes(x=factor(Vars, level=c('scale(log(W_SCALAR_mean))','scale(BAG_LIG_N)', 'scale(MICrK.i)'), labels = c("Soil moisture", "Lignin:N", "Copiotroph:oligotroph")),
                                    y=rel_ES, fill=factor(ID, level=c('Obs', 'SP', 'Cal'), labels = c("Observations", "Default", "Calibrated"))), 
             position = position_jitterdodge(dodge.width = 0.9), alpha=0.5, size=2) + #position_dodge(0.9)
  xlab("") + ylab("Relative effect size") + geom_abline(intercept=0, slope=0, color="black", linewidth=0.6) +
  theme_bw(base_size = 16)  + scale_fill_manual(name="Type", values = c(Observations = "black", Default = "white", Calibrated="darkgrey")) 
dev.off()

ggplot()+
geom_point(data=ES_points, aes(x=factor(Vars, level=c('scale(log(W_SCALAR_mean))','scale(BAG_LIG_N)', 'scale(MICrK.i)'), labels = c("Soil moisture", "Lignin:N", "Copiotroph:oligotroph")),
                          y=rel_ES, 
                          group=factor(ID, level=c('Obs', 'SP', 'Cal'), labels = c("Observations", "Default", "Calibrated"))), 
           position = position_dodge())
#collineraity?
vif(Obs_ES_mod) # for SP: all less than 5, up to 3.6
#for calibrated: (175) all less than 5, up to 3.6; (176) all less than 5, up to 4.7; (190) all less than 5, up to 4.0


#random forest#random forestposition_dodge2()
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




