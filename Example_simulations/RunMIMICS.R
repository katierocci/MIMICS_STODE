#load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)
library(purrr)
library(Metrics)

#load model components
source("functions/RXEQ.R")
source("Parameters/MIMICS_parameters_sandbox_20231129.R")
source("functions/calc_Tpars.R")
source("functions/MIMICS_calc_steady_state_pools.R")


##load forcing data
#long term ecological research site data
forcing_data2 <- read.csv("example_simulations/data/LTER_SITE_vwc.csv", as.is=T)
#colnames(forcing_data)[1] <- "Set"  # Fix for first column name import error
#SSA data
forcing_data <- read.csv("example_simulations/data/afsis_ref_updated4_MIMICS.csv", as.is=T)
forcing_data<- forcing_data[,-1]
colnames(forcing_data)[1] <- "Set"  # Fix for first column name import error
#match MIMICS names
forcing_data$ANPP <- forcing_data$NPP.gC.m2.yr
forcing_data$CLAY <- forcing_data$Clay_2um #Rose is using clay < 63 um
forcing_data$lig_N <- forcing_data$LIG_N
forcing_data$TSOI <- forcing_data$SoilTMP_C
forcing_data$theta_liq <- forcing_data$SoilMoi_m3m3 #this has been updated to be more correct
forcing_data <- forcing_data %>% drop_na(ANPP, CLAY, lig_N, TSOI, theta_liq)



#run the model for all sites
MIMruns <- forcing_data %>% split(1:nrow(forcing_data)) %>% map(~ MIMICS_SS(df=.))
MIMICS_ss_AllSites <- lapply(MIMruns, MIMICS_SS_format) %>% bind_rows()
MIMICS_ss_AllSites2 <- MIMICS_ss_AllSites[, c(1,124:142)] #just MIMICS inputs and outputs
MIMICS_ss_AllSites2$LIG_N <- MIMICS_ss_AllSites2$LIG_N...126

#plot data
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888",
                             "#E69F00", "#D55E00")

#look at pool variation pool by pool
hist(MIMICS_ss_AllSites2$SOMp)
hist(MIMICS_ss_AllSites2$MICk)

#Compare pools sizes and variation
MIM_ss_long <- MIMICS_ss_AllSites2 %>% select(SET,LITm, LITs, MICr, MICk, SOMa, SOMc, SOMp) %>% pivot_longer(cols = 2:8, names_to = 'Pool', values_to = 'C_content')
ggplot(MIM_ss_long) + geom_boxplot(aes(x=Pool, y=C_content, fill=Pool, group = Pool)) +ylab("Pool C content 0-30cm (mgC/cm3")+ xlab("Pool") + 
  scale_fill_manual(values=safe_colorblind_palette, guide="none") + theme_bw(base_size = 14)

#compare inputs and outputs
ggplot(MIMICS_ss_AllSites2) + geom_point(aes(x=CLAY, y=SOMp, colour = LIG_N), size=4) +ylab("SOMp 0-30cm (mgC/cm3)")+ xlab("Clay") +  theme_bw(base_size = 14)
