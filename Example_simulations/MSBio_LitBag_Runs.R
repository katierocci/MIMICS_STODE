#using Derek's litterbag code to run MSBio litter decomp simulations
rm(list = ls())

library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)

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
MSBio <- read.csv("Example_simulations/Data/Site_annual_clim.csv")
#match input data strucutre
#AGNPP should be in gDW!! multiply by 2 here to remedy
#don't have gravimetric soil moisture, just volumetric, assuming a BD of 1g/cm3 makes them equivalent - could be bad assumption given this is BD of leaves
MSBio2 <- MSBio %>% mutate(SITE = Site, ANPP = AGNPP_sum*2, TSOI = TSOI_mean, CLAY = PCT_CLAY_mean, lig_N = LIG_N, GWC = H2OSOI_mean*100, W_SCALAR=W_SCALAR_mean) %>%
  select(SITE, ANPP, TSOI, CLAY, LIG, C, N, CN, LIG_N, GWC, W_SCALAR) 
#daily data - change site name and MSBio2 row (1=BART, 8=SERC) to use different site daily input
BART_dailyinput <- read.csv("Example_simulations/Data//BART_clim.csv")

#check sum of daily fluxes
sum(BART_dailyinput$LITFALL)
plot(BART_dailyinput$LITFALL)
lines(BART_dailyinput$AGNPP)
2*sum(BART_dailyinput$AGNPP)
sum(BART_dailyinput$NPP)
sum(BART_dailyinput$LITFALL)
plot(BART_dailyinput$TSOI/10)
lines(BART_dailyinput$W_SCALAR)
lines(BART_dailyinput$H2OSOI)
lines(BART_dailyinput$TSOI/10)

BART_DI <- BART_dailyinput %>% mutate(DAY=X, ANPP = AGNPP*2, CLAY = rep(MSBio2[1,4], 366), LIG_N = rep(MSBio2[1,9], 366), GWC = H2OSOI*100) %>%
  select(DAY, ANPP, TSOI, CLAY, LIG_N, GWC, W_SCALAR) 

# two test dataframes with daily and annual results
BART_DI2 <- BART_dailyinput %>% mutate(SITE= rep("BART", 366),DAY=X, ANPP = rep(sum(LITFALL)*2,366), LITFALL=LITFALL*2, CLAY = rep(MSBio2[1,4], 366), 
                                       LIG_N = rep(MSBio2[1,9], 366), GWC = H2OSOI*100) %>%
  select(SITE, DAY, ANPP, LITFALL, TSOI, CLAY, LIG_N, GWC, W_SCALAR) 

BART_DI2_mean <- BART_DI2  %>% mutate( SITE= "BART", ANPP = mean(ANPP), TSOI = mean(TSOI), CLAY = mean(CLAY), 
                                       LIG_N= mean(LIG_N), GWC=mean(GWC), W_SCALAR=mean(W_SCALAR)) %>%
  select(SITE, ANPP, TSOI, CLAY, LIG_N, GWC, W_SCALAR) 
BART_DI2_mean = BART_DI2_mean[1,]

vmax_DI1 = exp(BART_DI2$TSOI * 0.063 + 5.47) * 0.000008 * BART_DI2$W_SCALAR *10
plot(vmax_DI1)
abline(h=mean(vmax_DI1),lt=2)
vmax_mean = exp(BART_DI2_mean$TSOI * 0.063 + 5.47) * 0.000008 * BART_DI2_mean$W_SCALAR *10
abline(h=vmax_mean,lw=2)
#Option 1: MSBio litter bags with just variation in NEON litter (not separated by species)
MSBio_BAGS <- read.csv("Example_simulations/Data/NEON_MSB_LitVars.csv")

#Option2: MSBio litter bags with leaf and litter chemistry combined
# MSBio_BAGS <- read.csv("NEON_MSB_LeafChem.csv")
# MSBio_BAGS <- MSBio_BAGS[,2:8]
# #rename to match input
# MSBio_BAGS2 <- MSBio_BAGS %>% mutate(Site = siteID, TYPE = taxonID, BAG_LIG = leaflig, BAG_N = leafN, BAG_CN = leafCN) %>%
#   select(Site, TYPE, BAG_LIG, BAG_N, BAG_CN)
# #add combined chem
# COMBO_BAGS <- data.frame(Site = c("BART", 'GRSM', 'HARV', 'LENO', 'MLBS', 'OSBS', 'SCBI', 'SERC', 'TALL', 'TREE', 'UNDE'),
#                          TYPE = rep("COMBO", 11),
#                          BAG_LIG = MSBio2$LIG,
#                          BAG_N = MSBio2$N,
#                          BAG_CN = MSBio2$CN)
# MSBio_BAGS3 <- rbind(MSBio_BAGS2, COMBO_BAGS)


### changed to fMET calculation in STODE script here!! Note that the two options are only somewhat related but less negatives in STODE equation

MSBio_BAGS$CALC_MET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * (MSBio_BAGS$BAG_LIG/MSBio_BAGS$BAG_N))
MSBio_BAGS$CALC_MET[MSBio_BAGS$CALC_MET <0] = 0 #setting negatives to zero - might want to reconsider this for future runs, all strucutral seems pretty unlikely
#MSBio_BAGS$CALC_N <- (1 / MSBio_BAGS$BAG_CN) / 2.5 * 100 #why calculating from CN and not N directly?
#MSBio_BAGS$CALC_MET2 <- 0.85 - 0.013 * MSBio_BAGS$BAG_LIG/MSBio_BAGS$CALC_N #calculate fMET

BAG_init_size <- 100
BAGS <- MSBio_BAGS %>% select(Site, TYPE, CALC_MET)
BAGS$BAG_LITm <- ((BAG_init_size * 1e3 / 1e4)/ depth) * BAGS$CALC_MET
BAGS$BAG_LITs <- ((BAG_init_size * 1e3 / 1e4)/ depth) * (1-BAGS$CALC_MET) #initial litter = 0.1 because of unit conversions here 


####
#run litterbag model 
####

#Individual litters decomposing too fast
#try N re-absorption
#try min and max of litter COMBO values
#compare to litter mass loss in observations

#SERC and all SERC litters (1) for 2 years
BAGS_SERC <- filter(BAGS, Site == "SERC" & TYPE == "mean")
BAGS_SERC <- BAGS_SERC[,2:5]
BAGS_out_SERC_SS <- BAGS_SERC %>% split(1:nrow(BAGS_SERC)) %>% map(~ MIMICS_LITBAG(litBAG=.,
                                                                           forcing_df=MSBio2[8,],
                                                                           #dailyInput = SERC_DI, 
                                                                           nspin_yrs=5,
                                                                           nspin_days=0,
                                                                           litadd_day=10,
                                                                           verbose=T)) %>% bind_rows()


#LENO and all LENO litters for 2 years
# BAGS_LENO <- filter(BAGS, Site == "LENO")
# BAGS_LENO <- BAGS_LENO[,2:5]
# MSBio_LENO_lsm <- MSBio2
# MSBio_LENO_lsm[4, 10] = 32 #from empirical data - LCI
# BAGS_out <- BAGS_LENO %>% split(1:nrow(BAGS_LENO)) %>% map(~ MIMICS_LITBAG(litBAG=.,
#                                                                            forcing_df=MSBio2[4,],
#                                                                            dailyInput = LENO_DI,
#                                                                            nspin_yrs=2,
#                                                                            nspin_days=0,
#                                                                            litadd_day=10,
#                                                                            verbose=T)) %>% bind_rows()

#UNDE and all UNDE litters (11) for 2 years
# BAGS_UNDE <- filter(BAGS, Site == "UNDE")
# BAGS_UNDE <- BAGS_UNDE[,2:5]
# BAGS_out <- BAGS_UNDE %>% split(1:nrow(BAGS_UNDE)) %>% map(~ MIMICS_LITBAG(litBAG=.,
#                                                                            forcing_df=MSBio2[11,],
#                                                                            nspin_yrs=2,
#                                                                            nspin_days=0,
#                                                                            litadd_day=10,
#                                                                            verbose=T)) %>% bind_rows()

#BART and all BART litters (1) for 2 years
BAGS_BART <- filter(BAGS, Site == "BART" & TYPE == "mean")
BAGS_BART <- BAGS_BART[,2:5]
BAGS_out_BART_SS <- BAGS_BART %>% split(1:nrow(BAGS_BART)) %>% map(~ MIMICS_LITBAG(litBAG=.,
                                                                           forcing_df=MSBio2[1,],
                                                                           #dailyInput = BART_DI,
                                                                           nspin_yrs=5,
                                                                           nspin_days=0,
                                                                           litadd_day=10,
                                                                           verbose=T)) %>% bind_rows()
ss_testM = MIMICS_SS(BART_DI2_mean[1,])
ss_testM[[2]]$VMAX
ss_testM[[2]]$KM
ss_testM[[2]]$I

ss_test = MIMICS_SS(BART_DI2)
ss_test[[2]]$VMAX
ss_test[[2]]$KM
ss_test[[2]]$I

ss_testM[[1]] - ss_test[[1]]

BAGS_BART <- filter(BAGS, Site == "BART" & TYPE == "mean")
BAGS_BART <- BAGS_BART[,2:5]
BAGS_out_BART_day1 <- BAGS_BART %>% split(1:nrow(BAGS_BART)) %>% map(~ MIMICS_LITBAG(litBAG=.,
                                                                                   forcing_df=BART_DI2,
                                                                                   dailyInput = BART_DI2,
                                                                                   nspin_yrs=20,
                                                                                   nspin_days=0,
                                                                                   litadd_day=10,
                                                                                   verbose=T)) %>% bind_rows()

# WRW tests for new code:
BAGS_BART2 <- filter(BAGS, Site == "BART" & TYPE == "mean")
BAGS_BART2 <- BAGS_BART2[,2:5]
BAGS_out_BART_daily <- BAGS_BART2 %>% split(1:nrow(BAGS_BART2)) %>% map(~ MIMICS_LITBAG(litBAG=.,
                                                                                   forcing_df=BART_DI2_mean[1,],
                                                                                   dailyInput = BART_DI2,
                                                                                   nspin_yrs=5,
                                                                                   nspin_days=0,
                                                                                   litadd_day=10,
                                                                                   verbose=T)) %>% bind_rows()

BAGS_out_BART_SS2 <- BAGS_BART2 %>% split(1:nrow(BAGS_BART2)) %>% map(~ MIMICS_LITBAG(litBAG=.,
                                                                                        forcing_df=BART_DI2,
                                                                                        #dailyInput = BART_DI2,
                                                                                        nspin_yrs=20,
                                                                                        nspin_days=0,
                                                                                        litadd_day=10,
                                                                                        verbose=T)) %>% bind_rows()
#all sites and all litters
# BAGS_mean <- filter(BAGS, TYPE == "mean")
# BAGS_input <- split(BAGS_mean, 1:nrow(BAGS_mean))
# forcing_input <- split(MSBio2, 1:nrow(MSBio2))
# BAGS_out_AllSites <- map2(forcing_input, BAGS_input, ~MIMICS_LITBAG(.x, .y, nspin_yrs=2, nspin_days=0, litadd_day=10, verbose=T)) %>% bind_rows()

####
#plot output
####

colorBlind7  <- c("#E69F00", "#56B4E9", "#009E73",
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #yellow (LENO), blue (SERC), green (UNDE)

#Formating observational data for comparing to field litter mass loss
Field_LML <- read.csv("Example_simulations/Data/Litter_decomp_all.csv")
#Add Species to group_by to get species-specific summary
LML_sum2 <- Field_LML %>% filter(site == 'BART') %>% group_by(time.point) %>% drop_na(percent.loss.litter) %>% summarize(mean.ML = mean(percent.loss.litter/2),
                                                                                                n = n(),
                                                                                                sd = sd(percent.loss.litter/2),
                                                                                                SE = sd/sqrt(n),
                                                                                                lci.ML = mean.ML - qt(1 - ((1 - 0.95) / 2), n - 1) * SE,
                                                                                                uci.ML = mean.ML + qt(1 - ((1 - 0.95) / 2), n - 1) * SE,
                                                                                                doy = mean(days_elapsed))


###
#Daily vs. STODE
####
daily = BAGS_out_BART_day1
SS = BAGS_out_BART_SS2
ggplot() +
  geom_line(data=daily, aes(y=LITm+LITs, x=DAY/365, color = "Litter", linetype ="daily"), linewidth=1, alpha=0.5) +
  geom_line(data=daily, aes(y=MICr+MICk, x=DAY/365, color = "MIC", linetype ="daily"), linewidth=1, alpha=0.5) +
  geom_line(data=daily, aes(y=SOMa+SOMc+SOMp, x=DAY/365, color = "SOM", linetype ="daily"), linewidth=1, alpha=0.5) +
  geom_line(data=SS, aes(y=MICr+MICk, x=DAY/365, color = "MIC", linetype ="steady state"), linewidth=2, alpha=0.5,) +
  geom_line(data=SS, aes(y=LITm+LITs, x=DAY/365, color = "Litter", linetype ="steady state"), linewidth=2, alpha=0.5) +
  geom_line(data=SS, aes(y=SOMa+SOMc+SOMp, x=DAY/365, color = "SOM", linetype ="steady state"), linewidth=2, alpha=0.5) +
  scale_color_manual(values = c("MIC"="#E69F00", "Litter"="#56B4E9", "SOM"="#009E73", "W_SCALAR" = "#F0E442", "ANPP" = "#CC79A7")) +
  scale_linetype_manual(values = c("daily"=1, "steady state"=2, "daily mean"=3)) + 
  ylab("C pools") +
  xlab("Year") +
  ggtitle("W_SCALAR moisture, beta at BART") +
  theme_bw(base_size = 20)


ggplot() +
  geom_line(data=daily, aes(y=LITm, x=DAY/365, color = "LITm", linetype ="daily"), linewidth=1, alpha=0.5) +
  geom_line(data=daily, aes(y=LITs, x=DAY/365, color = "LITs", linetype ="daily"), linewidth=1, alpha=0.5) +
  geom_line(data=daily, aes(y=MICr, x=DAY/365, color = "MICr", linetype ="daily"), linewidth=1, alpha=0.5) +
  geom_line(data=daily, aes(y=MICk, x=DAY/365, color = "MICk", linetype ="daily"), linewidth=1, alpha=0.5) +
  geom_line(data=SS, aes(y=LITm, x=DAY/365, color = "LITm", linetype ="daily"), linewidth=2, alpha=0.9) +
  geom_line(data=SS, aes(y=LITs, x=DAY/365, color = "LITs", linetype ="daily"), linewidth=2, alpha=0.9) +
  geom_line(data=SS, aes(y=MICr, x=DAY/365, color = "MICr", linetype ="daily"), linewidth=2, alpha=0.9) +
  geom_line(data=SS, aes(y=MICk, x=DAY/365, color = "MICk", linetype ="daily"), linewidth=2, alpha=0.9) +
  scale_y_log10() +
  scale_color_manual(values = c("MICr"="#E69F00", "MICk"="#56B4E9", "LITm"="#009E73", "LITs" = "#F0E442", "ANPP" = "#CC79A7")) +
  scale_linetype_manual(values = c("daily"=1, "steady state"=2, "daily mean"=3)) + 
  ylab("C pools") +
  xlab("year") +
  ggtitle("W_SCALAR moisture, NPP at BART") +
  theme_bw(base_size = 20)

ggplot() +
  geom_line(data=daily, aes(y=SOMa, x=DAY/365, color = "SOMa", linetype ="daily"), linewidth=1, alpha=0.5) +
  geom_line(data=daily, aes(y=SOMc, x=DAY/365, color = "SOMc", linetype ="daily"), linewidth=1, alpha=0.5) +
  geom_line(data=daily, aes(y=SOMp, x=DAY/365, color = "SOMp", linetype ="daily"), linewidth=1, alpha=0.5) +
  geom_line(data=SS, aes(y=SOMa, x=DAY/365, color = "SOMa", linetype ="daily"), linewidth=2, alpha=0.9) +
  geom_line(data=SS, aes(y=SOMc, x=DAY/365, color = "SOMc", linetype ="daily"), linewidth=2, alpha=0.9) +
  geom_line(data=SS, aes(y=SOMp, x=DAY/365, color = "SOMp", linetype ="daily"), linewidth=2, alpha=0.9) +
  scale_color_manual(values = c("SOMa"="#E69F00", "SOMc"="#56B4E9", "SOMp"="#009E73", "LITs" = "#F0E442", "ANPP" = "#CC79A7")) +
  scale_linetype_manual(values = c("daily"=1, "steady state"=2, "daily mean"=3)) + 
  ylab("C pools") +
  xlab("Year") +
  ggtitle("W_SCALAR moisture, NPP at BART") +
  theme_bw(base_size = 20)

