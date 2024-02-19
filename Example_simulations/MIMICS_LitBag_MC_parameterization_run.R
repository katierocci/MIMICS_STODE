### MIMICS MC for litterbag simulations


########################################
# Load R packages
########################################
library(dplyr)
library(rootSolve)
library(purrr)
library(furrr)


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
MSBio <- read.csv("Example_simulations/Data/Site_annual_clim.csv")
#match input data strucutre
#AGNPP should be in gDW!! multiply by 2 here to remedy
#don't have gravimetric soil moisture, just volumetric, assuming a BD of 1g/cm3 makes them equivalent - could be bad assumption given this is BD of leaves
data <- MSBio %>% mutate(SITE = Site, ANPP = AGNPP_sum*2, TSOI = TSOI_mean, CLAY = PCT_CLAY_mean, lig_N = LIG_N, GWC = H2OSOI_mean*100, W_SCALAR=W_SCALAR_mean) %>%
  select(SITE, ANPP, TSOI, CLAY, LIG, C, N, CN, LIG_N, GWC, W_SCALAR) 

#MSBio litter bags based variation in NEON litter (not separated by species)
MSBio_BAGS <- read.csv("Example_simulations/Data/NEON_MSB_LitVars.csv")
MSBio_BAGS$CALC_MET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * (MSBio_BAGS$BAG_LIG/MSBio_BAGS$BAG_N)) #gives some negative values...
MSBio_BAGS$CALC_MET[MSBio_BAGS$CALC_MET <0] = 0 #setting negatives to zero - might want to reconsider this for future runs, all strucutral seems pretty unlikely
BAG_init_size <- 100
BAGS <- MSBio_BAGS %>% select(Site, TYPE, CALC_MET)
BAGS$BAG_LITm <- ((BAG_init_size * 1e3 / 1e4)/ depth) * BAGS$CALC_MET
BAGS$BAG_LITs <- ((BAG_init_size * 1e3 / 1e4)/ depth) * (1-BAGS$CALC_MET) #initial litter = 0.1 because of unit conversions here
#each litter with each site
BAGS_mean <- BAGS %>% filter(TYPE == "mean")
#just one litter type
#BAGS_BART <- BAGS %>% filter(Site == "BART" & TYPE == "mean")


####################################
# Use the brute force MIMICS ftn
####################################

# Set desired number of random parameter runs
MIM_runs <- 100

### Create random parameter dataframe
## Parameter range informed by range observed over 10+ MCMC analysis results
rand_params <- data.frame( Tau_x = runif(MIM_runs, 0.3, 3),
  # Tau_r = runif(MIM_runs, 0.3, 3),
  # Tau_k = runif(MIM_runs, 0.3, 3)#,
  CUE_x = runif(MIM_runs, 0.5, 1.4), 
  # CUE_r = runif(MIM_runs, 0.5, 1.4),
  # CUE_k = runif(MIM_runs, 0.5, 1.4)
  #Vslope_x = runif(MIM_runs, 0.5, 2)#,
  #Vint_x = runif(MIM_runs, 0.8, 1.3)
  #Kslope_x = runif(MIM_runs, 0.5, 2),
  #Kint_x = runif(MIM_runs, 0.5, 2)#,
  vMOD_x = runif(MIM_runs, 0.5, 2),
  kMOD_x = runif(MIM_runs, 0.5, 2)  
)

rand_params$run_num <- seq(1,MIM_runs,1)

# must change variable in LitBag_MIMICS_repeat as well

# Set number of cores to use
no_cores <- availableCores() - 1
plan(multisession, gc = TRUE, workers = no_cores)

# Run MIMICS!

print(paste0("Starting ", MIM_runs, " runs"))
print(paste0("Start time: ", Sys.time()))

start_time <- Sys.time()


#below should take each row of the random parameters (a single run number) and run the MIMrepeat function for each row
#if only doing a single litter type, change mapping function in LitBag_MIMICS_repeat
MC_MIMICS <- rand_params %>% split(1:nrow(rand_params)) %>% future_map(~MIMrepeat(forcing_df = data, litBAG = BAGS_mean, rparams = .), .progress=TRUE) %>% 
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
