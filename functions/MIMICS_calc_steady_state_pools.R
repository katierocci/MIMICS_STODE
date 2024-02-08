# Required libraries
library(rootSolve)

###############################################
# MIMICS single point function
#> Input dataframe ("df") must contain columns: SITE, ANPP, fCLAY, TSOI, LIG_N or LIG + CN // Optional: MAT, GWC, W_SCALAR 
#>> With output required for running forward
###############################################
MIMICS_SS <- function(df){
  
  #DEBUG
  #df = LTER[1,]
  
  # Convert all column names to upper case
  colnames(df) <- toupper(colnames(df))
  
  ### Setup a var to collect run notes
  note <- ""
  
  ### Bring in forcing ANPP value, convert gDW to gC
  ANPP <- df$ANPP/2
  LITFALL <- df$LITFALL/2
  
  ### Bring in CLAY value, convert from percent to decimal
  fCLAY <- df$CLAY/100
  
  ### Bring in TSOI value
  TSOI <- df$TSOI
  
  ### Bring in lig:N forcing data
  LIG_N <- df$LIG_N
  LIG <- df$LIG
  CN <- df$CN

  # Bring in soil moisture information
  theta_liq  <- df$GWC/100  # GWC = Gravimetric water content
  theta_frzn <- 0           # Not used here. TODO Needs validation. 

  ### Bring in mean annual temperature data
  MAT <- df$MAT  
  
  # Use TSOI if 'MAT' column is not in forcing data
  if(is.null(MAT)){
    MAT <- TSOI
  }
  
  #Bring in W_SCALAR if present
  W_SCALAR = df$W_SCALAR

  ############################################################
  # MIMICS MODEL CODE STARTS HERE
  ############################################################
  # function calculates fMET with LIG_N if provided in input data.
  Tpars <- calc_Tpars_Conly(ANPP=ANPP, fCLAY=fCLAY, TSOI=TSOI, MAT=MAT,     
                            CN=CN, LIG=LIG, LIG_N=LIG_N,
                            theta_liq=theta_liq, theta_frzn=theta_frzn, 
                            W_SCALAR = W_SCALAR, LITFALL=LITFALL) 
  
  # Create arrays to hold output
  lit     <- Tpars$I
  mic     <- Tpars$I
  som     <- rep(NA, 3) 
  som[1]  <- Tpars$I[1]
  som[2]  <- Tpars$I[2]
  som[3]  <- Tpars$I[1] 
  CO2     <- rep(0, 2) 

  Ty    <- c(LIT_1 = lit[1], LIT_2 = lit[2], 
              MIC_1 = mic[1], MIC_2 = mic[2], 
              SOM_1 = som[1], SOM_2 = som[2], SOM_3 = som[3])
  
  # ## Set global parameters to ensure passing of variables to stode function
  # .GlobalEnv$VMAX <- Tpars$VMAX
  # .GlobalEnv$KM <- Tpars$KM
  # .GlobalEnv$fPHYS <- Tpars$fPHYS
  # .GlobalEnv$fCHEM <- Tpars$fCHEM
  # .GlobalEnv$fAVAI <- Tpars$fAVAI
  # .GlobalEnv$I <- Tpars$I
  # .GlobalEnv$tau <- Tpars$tau
  # .GlobalEnv$LITmin <- Tpars$LITmin
  # .GlobalEnv$SOMmin <- Tpars$SOMmin
  # .GlobalEnv$MICtrn <- Tpars$MICtrn
  # .GlobalEnv$desorb <- Tpars$desorb
  # .GlobalEnv$DEsorb <- Tpars$DEsorb
  # .GlobalEnv$OXIDAT <- Tpars$OXIDAT
  # .GlobalEnv$beta <- Tpars$beta
  # 
  # # ------------RUN THE MODEL-------------
  
  test  <- stode(y = Ty, time = 1e7, fun = RXEQ, parms = Tpars, positive = TRUE)
  
  ### Calc and get MIMICS output 
  MIMLIT    <- (test[[1]][[1]]+test[[1]][[2]]) * depth *1e4 / 1e6 #convert kgC/m2 from mgC/cm3 (0-30 cm) 
  MIMMIC    <- (test[[1]][[3]]+test[[1]][[4]]) * depth *1e4 / 1e6
  MIM_CO    <-  test[[1]][[3]]/test[[1]][[4]]
  MIMSOC    <- sum(test[[1]]) * depth * 1e4 / 1e6   
  
  table <- as.numeric(test[[1]])
  
  MIMout <- list()
  MIMout[[1]] <- as.numeric(test[[1]]) * depth *1e4 / 1e6   # MIMICS steady state C pools
  MIMout[[2]] <- Tpars                                      # MIMICS parameters used for/from simulation
  MIMout[[3]] <- df[1,]                                         # Forcing dataframe
  MIMout[[4]] <- as.numeric(test[[2]]) * depth *1e4 / 1e6   # CO2-r & CO2-k pools 
  
  return(MIMout)
}


#################################################
# Helper ftn to summarize MIMICS_SS() output
#################################################

MIMICS_SS_format <- function(MIMICS_SS_output) {
  
  MIMout_C_pools <- as.data.frame(t(as.data.frame(MIMICS_SS_output[[1]]))) 
  colnames(MIMout_C_pools) <- c("LITm", "LITs", "MICr", "MICk", "SOMp", "SOMc", "SOMa")
 
  MIMout_C_pools$MIMSOC <- rowSums(MIMout_C_pools[,1:7])
  MIMout_C_pools$MIMLIT <- MIMout_C_pools$LITm + MIMout_C_pools$LITs    
  MIMout_C_pools$MIMMIC <- MIMout_C_pools$MICr + MIMout_C_pools$MICk
  MIMout_C_pools$MIC_ratio <- MIMout_C_pools$MICr/MIMout_C_pools$MICk
  
  MIMout_CO2_pools <- as.data.frame(t(as.data.frame(MIMICS_SS_output[[4]])))
  colnames(MIMout_CO2_pools) <- c("CO2r", "CO2k")
  
  SOMp_desorb_yr <- as.numeric(MIMICS_SS_output[[2]]['desorb']) * 24 * 365
  MIMout_C_pools$SOMp_Turnover_yrs <- MIMout_C_pools$SOMp/SOMp_desorb_yr
  
  MIMout_single_tbl <- cbind(as.data.frame(MIMICS_SS_output[[3]]), MIMout_C_pools %>%
                               select(MIMSOC, MIMLIT, MIMMIC, MIC_ratio, LITm, LITs, MICr, 
                                      MICk, SOMa, SOMc, SOMp, SOMp_Turnover_yrs),
                             MIMout_CO2_pools)
  return(MIMout_single_tbl)
}


########################################################################
# Ftn that combines MIMICS_SS with format and returns pools dataframe
########################################################################

MIMICS_SS_pools <- function(df) {return(MIMICS_SS(df) %>% MIMICS_SS_format())}
