library(dplyr)
library(ggplot2)
library(DT)
library(purrr)
library(Metrics)

## Set working drive
setwd("C:/GitHub/MIMICS_STODE")

#----------------------------------------
# Load necessary model components
#----------------------------------------

# Load reverse Michaelis-Menton MIMICS function
source("functions/RXEQ.R")

# Load function to calculate model input variables
source("functions/calc_Tpars.R")

# Load function to calc MIMICS steady state pools
source("functions/MIMICS_calc_steady_state_pools.R")

# Load function to simulate litter bag decomposition
source("functions/MIMICS_sim_litterbag.R")

# Set MIMICS parameters
source("parameters/MIMICS_parameters_sandbox_20231129.R")

#--------------------------------------
# Datasets for example simulations
#--------------------------------------

# Forcing data for LTER sites
LTER <- read.csv("example_simulations/data/LTER_SITE_1.csv", as.is=T)

# Load LiDET litter bag data 
source("Utilities/load_LiDET_litterbags.R")
  # See the LiDET_BAGS' & 'BAGS' dataframes

#---------------------------------------------
# Example: Run single litter bag decomp simulation
#---------------------------------------------

# Use an LTER site
#forcing_input <- LTER[6,]

# Build your own example
forcing_input <- data.frame(Site = "TEST SITE",
                            ANPP = 744,
                            TSOI = 15,
                            MAT = 25,
                            CLAY = 15,
                            LIG = 21,
                            N = 1.02,
                            CN = 49,
                            LIG_N = 20.588,
                            GWC = 35)

# Lets build a daily temperature curve
#------------------------------------------
source("Utilities/temp_curve_builder.R")
# Create the temperature vector
daily_TSOI <- generate_temp_curve(winter_low = 5, summer_high = 25, mean_temp = 15)

# Now build out the daily inputs, keeping all same except TSOI
dailyInputs <- data.frame(ANPP = rep(forcing_input$ANPP, 365), 
                          CLAY = rep(forcing_input$CLAY, 365), 
                          TSOI = daily_TSOI,
                          MAT = rep(forcing_input$MAT, 365), 
                          LIG_N = rep(forcing_input$LIG_N, 365),
                          CN = rep(forcing_input$CN, 365),  
                          LIG = rep(forcing_input$LIG, 365),
                          GWC = rep(forcing_input$GWC, 365), 
                          theta_frzn = rep(0, 365))

# Check steady state pool values
ssMIMOUT <- MIMICS_SS(forcing_input) %>% MIMICS_SS_format()

# Select a litter type from the LiDET example litter bags
BAG_input <- BAGS[6,]

# Run litterbag decomp simulation
BAG_out <- MIMICS_LITBAG(forcing_input, BAG_input, dailyInput=dailyInputs, loop_dailyInput=TRUE, nspin_yrs=10, nspin_days=150, litadd_day=125, verbose=T)

# Calc MIMICS C pools
BAG_out$MIMSOC <- rowSums(BAG_out[,6:12])
BAG_out$MIMMIC <- rowSums(BAG_out[,8:9])  
BAG_out$MIMLIT <- rowSums(BAG_out[,6:7])  
BAG_out$MICtoSOC <- BAG_out$MIMMIC/BAG_out$MIMSOC

# Plot litterbag decomp over time
ggplot(BAG_out) + 
  geom_line(aes(y=LITBAGm, x=DAY, colour="LITm"), size=1, linetype="dashed") +
  geom_line(aes(y=LITBAGs, x=DAY, colour="LITs"), size=1, linetype="dashed") +
  geom_line(aes(y=LITBAGs+LITBAGm, x=DAY, colour="Total"), size=1) +
  scale_colour_manual(values = c("Dark Green","Orange","Black")) +
  ylab("Litter Bag C") +
  xlab("Day") +
  labs(title = "MIMICS Litter Bag Simulation",
       subtitle = paste0("Site: ", forcing_input$Site, "\nLitter type: ", BAG_input[1]),
       colour = "LIT Pool") +
  theme_bw()


