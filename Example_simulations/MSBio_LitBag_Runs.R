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
#remove bottom rows if not using species specific litter
MSBio <- MSBio[1:11,]
#match input data structure
#AGNPP should be in grams Dry Weight (gDW) not gC! multiply by 2 here to remedy
#don't have gravimetric soil moisture, just volumetric, assuming a BD of 1g/cm3 makes them equivalent - could be bad assumption given this is BD of leaves
MSBio2 <- MSBio %>% mutate(SITE = Site, ANPP = AGNPP_sum*2, TSOI = TSOI_mean, CLAY = PCT_CLAY_mean, GWC = H2OSOI_mean*100, W_SCALAR=W_SCALAR_mean) %>%
  select(SITE, ANPP, TSOI, CLAY, LIG_N, LIG_N_sp1, LIG_N_sp2, LIG_N_sp3, GWC, W_SCALAR, lci_SM_ratio, uci_SM_ratio) 
#fixing TALL and OSBS ANPP
NEON_GPP <- read.csv("Example_simulations/Data/NEON_GPP.csv")
MSBio3 <- MSBio2
MSBio3$ANPP[MSBio3$SITE == "OSBS"] <- 547 + 0.18*NEON_GPP[6,2] #using relationship between NEON GPP and ANPP
MSBio3$ANPP[MSBio3$SITE == "TALL"] <- 547 + 0.18*NEON_GPP[9,2] #using relationship between NEON GPP and ANPP 
#filtering for only sites with microbial data to match observations
Mic_sites <- c("SERC","BART","TALL","TREE","LENO","HARV","GRSM")
MSBio_sites <- filter(MSBio3, SITE %in% Mic_sites)
#loading daily inputs and replacing TALL data to be more realistic
DailyInput <- read.csv("Example_simulations/Data/DailyInput.csv")
DailyInput$LITFALL[DailyInput$SITE == "TALL"] <- DailyInput$LITFALL[DailyInput$SITE == "TALL"]*0.663
DailyInput$ANPP[DailyInput$SITE == "TALL"] <- sum(DailyInput$LITFALL[DailyInput$SITE == "TALL"])

#checking relationship between soil moisture and W_Scalar
# ggplot(DailyInput, aes(x=GWC, y=W_SCALAR, color=SITE, group=SITE)) + geom_point() +geom_smooth(method = "loess")
# #mostly hump shaped curve but some sites don't reach downward part of curve and TREE is just mega weird
# xx<-seq(0,80,0.0312378)
# ggplot(DailyInput, aes(x=GWC, y=W_SCALAR)) + geom_point() +geom_smooth(method = "lm", formula = y~x + I(x^2), color="red")+
#   #lines(x=xx,y=test$coef[3]*xx^2+test$coef[2]*xx+test$coef[1]) #not working! Try again Monday!
# geom_line(aes(x=xx,y=test$coef[3]*xx^2+test$coef[2]*xx+test$coef[1]))
#   #, col='blue
# #if you include all sites there is broadly a hump shaped curve
# y=DailyInput$W_SCALAR
# x=DailyInput$GWC
# test <- lm(y~x + I(x^2))
# summary(test)
# test$coef

#additional code for: (1) determining ANPP multipliers for RALL and OSBS; (2) creating DailyInput file for all sites; (3) determining litterfall multiplier for TALL
#(1)
# NPPComp <-MSBio %>% mutate(SITE=Site) %>% inner_join(NEON_GPP, by = "SITE") 
# ggplot(NPPComp, aes(x=Annual.GPP, y= AGNPP_sum*2, color=SITE)) + geom_point(size=4) + theme_bw(base_size = 16)
# NPP_GoodSites <- NPPComp %>% filter(SITE != "TALL") %>% filter(SITE !="OSBS")
# NPP_mod <- lm((AGNPP_sum*2)~Annual.GPP, data=NPP_GoodSites)
# summary(NPP_mod) #547 + 0.18x
#(2)
#daily data - change site name and MSBio2 row (1=BART, 8=SERC) to use different site daily input
# BART_dailyinput <- read.csv("Example_simulations/Data/BART_clim.csv")
# BART_DI <- BART_dailyinput %>% mutate(DAY=X, ANPP = rep(sum(LITFALL)*2,366), LITFALL=LITFALL*2, CLAY = rep(MSBio2[1,4], 366), 
#                                        LIG_N = rep(MSBio2[1,9], 366), GWC = H2OSOI*100, MAT=TBOT) %>%
#   select(DAY, ANPP, LITFALL, TSOI, MAT, CLAY, LIG_N, GWC, W_SCALAR) 
# DailyInput <- rbind(BART_DI, GRSM_DI, HARV_DI, LENO_DI, SERC_DI, TALL_DI, TREE_DI)
# DailyInput$SITE <- c(rep("BART", 366), rep("GRSM", 366), rep("HARV", 366), rep("LENO", 365), rep("SERC", 366), rep("TALL", 366), rep("TREE", 366))
# write.csv(DailyInput, "Example_simulations/Data/DailyInput.csv")
#(3)
#replace TALL with more realistic data
# test1 <- filter(DailyInput, DAY == 40)
# test2 <- inner_join(MSBio_sites, test1, by = "SITE")
# ggplot(test2, aes(x=ANPP.x, y=ANPP.y, color=SITE)) + geom_point(size=3) + theme_bw(base_size = 16)
# LF_GoodSites <- test2 %>% filter(SITE != "TALL")
# LF_mod <- lm(ANPP.y~ANPP.x, data=LF_GoodSites)
# summary(LF_mod) #-7.2 + 1.7x
# TALL_ANPP.adj <- -7.2 + 1.7*MSBio_sites[6,2]
# TALL_mult <- TALL_ANPP.adj/test1[6,3] #0.663



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


#species specific litter at steady state
#315th day of the year is 11/11/21 which is the average day the litter was deployed
BAGS_out_AllSites_ss = data.frame()
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
      BAGS_out_AllSites_ss <- rbind(BAGS_out_AllSites_ss,BO_DI)
    }
  }
}

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


#comparison of baseline MIMICS to LML
#got rid of LITi calcs since all litters back to the same starting value
#LIT_init <- BAGS_out_AllSites_DI %>% filter(DAY == 10) %>% mutate(LITi = LITBAGm+LITBAGs) %>% 
#  mutate(SITE.LT = paste(SITE, Litter_Type, sep=".")) %>% select(SITE.LT, LITi)
#boxplot(LIT_init$LITi)
BAGS_out_plot <- BAGS_out_AllSites_ss  %>% mutate(SITE.LT = paste(SITE, Litter_Type, sep=".")) %>% mutate(LIT_PerLoss = ((0.1 - (LITBAGm+LITBAGs))/0.1)*100)
#wide format for plotting
BAGS_out_wide = BAGS_out_plot %>% select(SITE, Litter_Type, SM_Type, DAY, LIT_PerLoss) %>% 
  pivot_wider(names_from = Litter_Type, values_from = LIT_PerLoss)
#plotting - check!
ggplot() +
  geom_line(data=BAGS_out_wide, aes(y=100-LIG_N_sp1, x=DAY-315, group=SITE, color=SITE), linewidth=1.5, alpha=0.5) +
  geom_line(data=BAGS_out_wide, aes(y=100-LIG_N_sp2, x=DAY-315, group=SITE, color=SITE), linewidth=1.5, alpha=0.5) +
  geom_line(data=BAGS_out_wide, aes(y=100-LIG_N_sp3, x=DAY-315, group=SITE, color=SITE), linewidth=1.5, alpha=0.5) +
  #geom_ribbon(data=BAGS_out_wide, aes(y=100-LIG_N, x=DAY, ymin = 100-LIG_N_min, ymax=100-LIG_N_max, group=SITE, fill=SITE), alpha = 0.3) +
  geom_point(data=LML_sum2, aes(y=100-mean.ML, x=doy, group=site, color=site), size = 3) +
  geom_errorbar(data=LML_sum2, aes(y=100-mean.ML, x=doy, ymin = 100-lci.ML, ymax = 100-uci.ML, group=site, color=site), width=0,linewidth=1) +
  xlim(0, 780) +
  ylab("Litter Bag C Remaining (%)") +
  xlab("Days elapsed") +
  facet_wrap(.~SM_Type) +
  theme_bw(base_size = 20)
#summary data - works ok for visualization!
BO_plot_sum <- BAGS_out_plot %>% group_by(SITE,DAY,SM_Type) %>% summarise(mean=mean(LIT_PerLoss), min=min(LIT_PerLoss), max=max(LIT_PerLoss))
ggplot() +
  geom_ribbon(data=BO_plot_sum, aes(y=100-mean, x=DAY-315, ymin = 100-min, ymax=100-max, group=SITE, fill=SITE), alpha = 0.3) +
  geom_point(data=LML_sum2, aes(y=100-mean.ML, x=doy, group=site, color=site), size = 3) +
  geom_errorbar(data=LML_sum2, aes(y=100-mean.ML, x=doy, ymin = 100-lci.ML, ymax = 100-uci.ML, group=site, color=site), width=0,linewidth=1) +
  ylab("Litter Bag C Remaining (%)") +
  xlab("Day") +
  xlim(0, 780) +
  facet_wrap(.~SM_Type) +
  theme_bw(base_size = 20)


#comparing spin up with mean vs individual litters
BO_var <- filter(BAGS_out_AllSites_var, DAY == 315)
BO_ss <- filter(BAGS_out_AllSites_ss, DAY == 315)
ggplot(BO_var, aes(x=SITE, y=MICr/MICk, color=SITE)) +
  geom_point(size=4) +
  #geom_point(BAGS_out_AllSites_var, aes(x=DAY, y=MICr/MICk, color=SITE), shape=2) +
  facet_grid(Litter_Type~SM_Type)
ggplot(BO_var, aes(x=LITm/LITs, y=MICr/MICk, color=SITE)) +
  geom_point(size=4) +
  #geom_point(BAGS_out_AllSites_var, aes(x=DAY, y=MICr/MICk, color=SITE), shape=2) +
  facet_grid(Litter_Type~SM_Type)
ggplot(BO_ss, aes(x=LITm/LITs, y=MICr/MICk, color=SITE)) +
  geom_point(size=4) +
  #geom_point(BAGS_out_AllSites_var, aes(x=DAY, y=MICr/MICk, color=SITE), shape=2) +
  facet_grid(Litter_Type~SM_Type)


# #formating data for modelVobs, RWA, and effect size
#need to standardize and log this data before running it thru again
#MFG_stdzd <- MFG_analysis %>% mutate_at(c('C_O', 'LIG_N', 'log.vwc'), ~(scale(.) %>% as.vector))
FieldData <- LML_sum2 %>% mutate(DAY=doy, SITE=site) %>% mutate(SITE.DAY=paste(SITE, DAY, sep=".")) %>% select(time.point,SITE.DAY, mean.ML)
#LIT_init <- BAGS_out_AllSites_DI %>% filter(DAY == 10) %>% mutate(LITi = LITBAGm+LITBAGs) %>% select(SITE, LITi)
#boxplot(LIT_init$LITi)
#df <- BAGS_out_AllSites_DI %>% left_join(LIT_init, by = "SITE")
LITi = 0.1
df_LML <- BAGS_out_AllSites_ss %>% mutate(DAY.LitOut = DAY -314) %>% mutate(SITE.DAY=paste(SITE, DAY.LitOut, sep=".")) %>% 
  right_join(FieldData, by="SITE.DAY") %>% mutate(LIT_PerLoss = ((LITi - (LITBAGm+LITBAGs))/LITi)*100)
#creating table with LML for each day
# df_LML_All <- BAGS_out_AllSites_var %>% filter(DAY>315) %>% mutate(LIT_PerLoss = ((LITi - (LITBAGm+LITBAGs))/LITi)*100)
# #filtering so each site only has days for which we have observations 
# T2_DAY <- df_LML %>% filter(time.point==2) %>% mutate(T2.DAY = DAY) %>% select(SITE, T2.DAY)
# T1_DAY <- df_LML %>% filter(time.point==1) %>% mutate(T1.DAY = DAY) %>% select(SITE, T1.DAY)
# df_LML_RI <- df_LML_All %>% inner_join(T2_DAY, by = "SITE") %>% inner_join(T1_DAY, by = "SITE") %>% filter(DAY <= T2.DAY+315) #accounts for 10 days of spinup before starting simulation
# df_LML_T1T2 <- df_LML_RI %>% filter(DAY==T2.DAY | DAY==T1.DAY)

#comparison of modeled and observed decomp at timepoint 2
modelVobs <- lm(LIT_PerLoss~mean.ML, data=df_LML)
summary(modelVobs)
sqrt(mean((df_LML$mean.ML - df_LML$LIT_PerLoss)^2))
ggplot(df_LML, aes(x=mean.ML, y=LIT_PerLoss)) + geom_point(aes(color=SITE), size=4) + geom_smooth(method = "lm", color="black")  + xlim(0,80) + ylim(0,80) +
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
DI_means <- DailyInput_SM %>% mutate(SITE.SM = paste(SITE, SM_type, sep = ".")) %>% group_by(SITE.SM) %>% 
  summarise(W_SCALAR=mean(W_SCALAR), MAT=mean(MAT)) %>% select(SITE.SM, W_SCALAR, MAT)
#initial MICrK
MIC_init <- BAGS_out_AllSites_ss %>% filter(DAY == 315) %>% mutate(MICrK.i =  MICr/MICk)%>%
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% select(SITE.SM.LQ, MICrK.i) 
#need bag means!
BAGS_LIGN <- MSBio_BAGS %>% mutate(SITE.LQ = paste(SITE, TYPE, sep = ".")) %>% select(SITE.LQ, BAG_LIG_N)
df_analysis <- df_LML %>% mutate(MICrK = MICr/MICk) %>% mutate(MIC=MICr+MICk) %>% mutate(SOC = SOMa+SOMc+SOMp) %>% 
  mutate(SITE.SM.LQ = paste(SITE, SM_Type, Litter_Type, sep = ".")) %>% mutate(SITE.SM = paste(SITE, SM_Type, sep = ".")) %>%
  mutate(SITE.LQ = paste(SITE, Litter_Type, sep = ".")) %>% inner_join(DI_means, by="SITE.SM") %>% 
  inner_join(MIC_init, by="SITE.SM.LQ") %>% inner_join(BAGS_LIGN, by="SITE.LQ")
#logical checks
df_check <- df_analysis %>% filter(MICrK > 0.01) %>%
  filter(MICrK < 100) %>%
  filter(MIC/SOC > 0.0001) %>%
  filter(MIC/SOC < 0.40) 
plot((df_check$MICr/df_check$MICk), df_check$BAG_LIG_N)
#relative weights analysis
df_check$log_WS <- log(df_check$W_SCALAR)
df_rwa <- df_check %>% select(LIT_PerLoss, log_WS, BAG_LIG_N, MICrK.i) %>% 
  mutate_at(vars(c("log_WS", "BAG_LIG_N", "MICrK.i")), ~(scale(.) %>% as.vector))
MIM_rwa <- rwa(df_rwa, "LIT_PerLoss", c("log_WS", "BAG_LIG_N", "MICrK.i"), applysigns = TRUE, plot = TRUE) #"MAT", 
plot_rwa(MIM_rwa) 
#effect size
Obs_ES_mod <- lmer(LIT_PerLoss ~ scale(log(W_SCALAR))+scale(BAG_LIG_N)+scale(MICrK.i)+ (1|SITE), data = df_check) #MAT #scale(BAG_LIG_N)+
summary(Obs_ES_mod)
Obs_ES <- as.data.frame(fixef(Obs_ES_mod)) #fixed effects coefficients as effect size
#Obs_ES <- as.data.frame(Obs_ES_mod$coefficients) #if not using fixed effects model do this
Obs_ES$Vars <- rownames(Obs_ES)
colnames(Obs_ES)[1] <- "value"
Obs_ES <- Obs_ES[-1, ]
Obs_ES$mult <- ifelse(Obs_ES$value <0, -1, 1)
Obs_ES$rel_ES <- (abs(Obs_ES$value)/sum(abs(Obs_ES$value))) * 100 * Obs_ES$mult
ggplot(Obs_ES, aes(x=Vars, y=rel_ES)) + geom_bar(stat="identity", fill="blue") + coord_flip() + geom_text(aes(label=round(rel_ES, digits=1)), color="red", size=5) +theme_bw(base_size = 16)


#random forest
#test and training data - using 75 train-25 test split like in Georgiou et al., 2021
MFG_rf <- as.data.frame(df_check %>% select(LIT_PerLoss, MAT, W_SCALAR, LIG_N, MICrK) %>% na.omit(.))
split <- sample.split(MFG_rf, SplitRatio = 0.75)
data_train <- subset(MFG_rf, split == "TRUE")
data_test <- subset(MFG_rf, split == "FALSE")
#run training data with ranger
#MFG_analysis_NAo <- MFG_analysis %>% na.omit(.)
ranger_train <- ranger(LIT_PerLoss ~ MAT + W_SCALAR + LIG_N + MICrK,
                       data = data_train,
                       importance = 'impurity',
                       mtry = 1)
print(ranger_train) #really similar to rF package
pred_test2 <- predict(ranger_train, data = data_test)
summary(r.rap <- data_test$LIT_PerLoss - pred_test2$predictions)
(rmse.ra <- sqrt(sum(r.rap^2)/length(r.rap))) #11% - lower than rF
plot(data_test$LIT_PerLoss ~ pred_test2$predictions, asp=1, pch=20, xlab="fitted", ylab="actual", main="Prediciton of Litter Decomposition, Ranger")
grid(); abline(0,1)
#comparing variable importance
vip(ranger_train, title = "Ranger")
#ensemble of random forests
#function for getting random forest vips - removing train-test becasue using all data for this
rf.vip <- function(df) {
  MFG_rf <- as.data.frame(df %>% select(LIT_PerLoss, MAT, W_SCALAR, LIG_N, MICrK) %>% na.omit(.))
  ranger_train <- ranger(LIT_PerLoss ~ MAT + W_SCALAR + LIG_N + MICrK,
                       data = MFG_rf,
                       importance = 'impurity',
                       mtry = 1)
  Obs_rf <- as.data.frame(ranger_train$variable.importance) #fixed effects coefficients as effect size
  Obs_rf$Vars <- rownames(Obs_rf)
  colnames(Obs_rf)[1] <- "value"
  Obs_rf$rel_rf <- (abs(Obs_rf$value)/sum(abs(Obs_rf$value))) * 100
  return(Obs_rf)
}
y = rf.vip(df_check)
y %>% mutate(Vars=factor(Vars, levels = c("MICrK", "LIG_N", "MAT", "W_SCALAR"))) %>% ggplot(aes(x=Vars, y=rel_rf)) + geom_bar(stat="identity", fill="blue") + coord_flip() + geom_text(aes(label=round(rel_rf, digits=1)), color="red", size=7) +theme_bw(base_size = 16)
#replicate function 10 times
result <- t(replicate(10, rf.vip(df_check)))
for (i in 1:10) {
  rf_vip_out <- rbind(rf_vip_out, as.data.frame(result[i,]))
}
#assess variaiton in rf ensemble
ggplot(rf_vip_out, aes(x=Vars, y=rel_rf, color=Vars)) + geom_boxplot() + theme_bw(base_size = 16)


#pearson and spearman correlation
df_cor <- df_check %>% select(LIT_PerLoss, MAT, W_SCALAR, LIG_N, MICrK)
cor.p <- cor(df_cor, method="pearson")
corrplot(cor.p, type = "upper") #aligned with RWA
cor.s <- cor(df_cor, method="spearman")
corrplot(cor.s, type = "lower")#aligned with RWA

#simple realtionships between variables
ggplot(df_check, aes(x=W_SCALAR, y=LIT_PerLoss)) + geom_point(size=3, aes(color=SITE)) + geom_smooth(method="lm") + theme_bw(base_size = 16)
ggplot(df_check, aes(x=MAT, y=LIT_PerLoss)) + geom_point(size=3, aes(color=SITE)) + geom_smooth(method="lm") + theme_bw(base_size = 16)
ggplot(df_check, aes(x=LIG_N, y=LIT_PerLoss)) + geom_point(size=3, aes(color=SITE)) + geom_smooth(method="lm") + theme_bw(base_size = 16)
ggplot(df_check, aes(x=MICrK.i, y=LIT_PerLoss)) + geom_point(size=3, aes(color=SITE)) + geom_smooth(method="lm") + theme_bw(base_size = 16)

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




