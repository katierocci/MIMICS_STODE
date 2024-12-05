### MIMICS MC for litterbag simulations


########################################
# Load R packages
########################################
library(dplyr)
library(rootSolve)
library(purrr)
library(furrr)
library(tidyr)


########################################
# Load MIMICS data and ftns
########################################
source("Parameters/MIMICS_parameters_sandbox_20231129.R") #Sets the initial parameters
source("functions/MIMICS_sim_litterbag.R")
source("functions/MC_parameterization/MIMICS_LitBag_repeat.R")
source("functions/MC_parameterization/set_parameter_defaults.R")

########################################
# Load forcing data
########################################

#load site data
MSBio <- read.csv("Example_simulations/Data/Site_annual_clim_final.csv")
#match input data structure
#AGNPP should be in grams Dry Weight (gDW) not gC! multiply by 2 here to remedy
MSBio2 <- MSBio %>% mutate(SITE = Site, ANPP = LITFALL_sum*2, TSOI = TSOI_mean, CLAY = PCT_CLAY_mean, GWC = H2OSOI_mean*100, W_SCALAR=W_SCALAR_mean) %>%
  select(SITE, ANPP, TSOI, CLAY, LIG_N, LIG_N_sp1, LIG_N_sp2, LIG_N_sp3, GWC, W_SCALAR, lci_SM_ratio, uci_SM_ratio) 
#fixing TALL ANPP
DailyInput <- read.csv("Example_simulations/Data/DailyInput.csv") %>% select(-MAT)
DailyInput$LITFALL[DailyInput$SITE == "TALL"] <- DailyInput$LITFALL[DailyInput$SITE == "TALL"]*0.60
DailyInput$ANPP[DailyInput$SITE == "TALL"] <- sum(DailyInput$LITFALL[DailyInput$SITE == "TALL"])
#replacing MSBio data with daily input sums and means to ensure comparable data between daily data and annual data
DI_sum <- DailyInput %>% select(SITE, ANPP, TSOI, W_SCALAR) %>% group_by(SITE, ANPP) %>% summarise(TSOI=mean(TSOI), W_SCALAR=mean(W_SCALAR))
MSBio3 <- MSBio2 %>% select(-ANPP, -TSOI, -W_SCALAR) %>% inner_join(DI_sum, by="SITE")
#filtering for only sites with microbial data to match observations
Mic_sites <- c("SERC","BART","TALL","TREE","LENO","HARV","GRSM")
data <- filter(MSBio3, SITE %in% Mic_sites)


#formatting data for multiple soil moisture types
data_SM <- rbind(data, data, data)
data_SM <- data_SM %>% mutate(SM_type = c(rep("mean", 7), rep("max", 7), rep("min", 7))) %>% 
  mutate(W_SCALAR2 = case_when(SM_type == "mean" ~ W_SCALAR,
                               SM_type == "max" ~ W_SCALAR*uci_SM_ratio,
                               SM_type == "min" ~ W_SCALAR*lci_SM_ratio)) %>%
  mutate(W_SCALAR2 = case_when(W_SCALAR2>1~1, TRUE ~ W_SCALAR2)) %>%
  mutate(W_SCALAR = W_SCALAR2)
DailyInput_SM <- rbind(DailyInput, DailyInput, DailyInput)
SM_mult <- data %>% select(SITE, uci_SM_ratio, lci_SM_ratio)
DailyInput_SM <- DailyInput_SM %>% left_join(SM_mult, by="SITE") %>% mutate(SM_type = c(rep("mean", 2561), rep("max", 2561), rep("min", 2561))) %>% 
  mutate(W_SCALAR2 = case_when(SM_type == "mean" ~ W_SCALAR,
                               SM_type == "max" ~ W_SCALAR *uci_SM_ratio,
                               SM_type == "min" ~ W_SCALAR *lci_SM_ratio)) %>%
  mutate(W_SCALAR2 = case_when(W_SCALAR2>1~1, TRUE ~ W_SCALAR2)) %>%
  mutate(W_SCALAR = W_SCALAR2)


#load in MSBio litter bag chemistry
MSBio_BAGS <- data %>% select(SITE, LIG_N_sp1, LIG_N_sp2, LIG_N_sp3) %>% pivot_longer(2:4, names_to = "TYPE", values_to = "BAG_LIG_N")
MSBio_BAGS$CALC_MET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * (MSBio_BAGS$BAG_LIG_N))
MSBio_BAGS$CALC_MET[MSBio_BAGS$CALC_MET <0] = 0.01 #setting negatives to small number so 99% structural

BAG_init_size <- 100
BAGS <- MSBio_BAGS %>% select(SITE, TYPE, CALC_MET)
BAGS$BAG_LITm <- ((BAG_init_size * 1e3 / 1e4)/ depth) * BAGS$CALC_MET #g/m2 converted to mg/cm3
BAGS$BAG_LITs <- ((BAG_init_size * 1e3 / 1e4)/ depth) * (1-BAGS$CALC_MET) 
BAGS_mean <- filter(BAGS, TYPE == "LIG_N_sp1")


####################################
# Use the brute force MIMICS ftn
####################################

# Set desired number of random parameter runs
MIM_runs <- 5

### Create random parameter dataframe
## Parameter range informed by range observed over 10+ MCMC analysis results
rand_params <- data.frame(
  Tau_r = runif(MIM_runs, 0.3, 2), #original range: 0.3,3; secondary range: 0.3, 2
  vMOD_m = runif(MIM_runs, 0.5, 2), #original range: 0.5,2
  vMOD_s = runif(MIM_runs, 0.5, 2), #original range: 0.5,2
  beta_x = runif(MIM_runs, 0.67 ,1.33) #varies from 0.67 to 1.33 since default is 1.5 and can only vary between 1 and 2
)

rand_params$run_num <- seq(1,MIM_runs,1)


# Set number of cores to use
no_cores <- availableCores() - 1
plan(multisession, gc = TRUE, workers = no_cores)

# Run MIMICS!
print(paste0("Starting ", MIM_runs, " runs"))
print(paste0("Start time: ", Sys.time()))

start_time <- Sys.time()


#DAILY INPUT WITH MULTIPLE LQ AND SOIL MOISTURE
MC_MIMICS <- rand_params %>% split(1:nrow(rand_params)) %>% future_map(~MIMrepeat(forcing_df = data_SM, litBAG = BAGS, dailyInput = DailyInput_SM, rparams = .), .progress=TRUE) %>%
  bind_rows()


wall_time <- Sys.time() - start_time
print(paste0("Wall time: ", as.character(wall_time)))


# Release CPU cores
plan(sequential)
nbrOfWorkers()

# Clean up memory
gc()

#####################################


##########################################
# Save MC output data - for use with computing clusters
##########################################
#saveRDS(MC_MIMICS, paste0("temp/MSBio_MC_", as.character(MIM_runs), "_", format(Sys.time(), "%Y%m%d_%H%M%S_"),  ".rds"))
