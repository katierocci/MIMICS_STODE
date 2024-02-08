
MIMICS_LITBAG <- function(forcing_df, litBAG, dailyInput=NA, loop_dailyInput=TRUE, nspin_yrs=10, nspin_days=200, litadd_day=143, verbose=T){
  
  #DEBUG
  # forcing_df <- LTER[6,]
  # litBAG <- BAGS[1,]
  # nspin_yrs <- 10
  # nspin_days <- 200
  # litadd_day <- 143
  # dailyInput <- dailyInputs
  # loop_dailyInput = TRUE
  # verbose = TRUE
  
  # Convert forcing data column names to uppercase
  colnames(forcing_df) <- toupper(colnames(forcing_df))
  
  if(verbose){
    print("-------------------------------------------------------")
    print(paste0("Starting ", forcing_df$SITE[1], " - ", litBAG[1]))
    print("-------------------------------------------------------")
  }
  # drop litterfall if proivded
  forcing_df_ss = forcing_df %>%
    select(ANPP, TSOI, CLAY, LIG_N, GWC, W_SCALAR) 
  # Get MIMICS steady_state output
  MIMss <- MIMICS_SS(forcing_df)

  # Create dataframe to store MIMICS pools over timesteps
  MIMfwd = MIMss[[1]]
  MIMfwd = (MIMfwd / depth) / (1e4 / 1e6)  #convert from gC m-2 to mgC cm-3
  
  # Get Tpars from ss simulation  
  Tpars = MIMss[[2]]

  #Init arrays to store daily output data
  nday   <- 365 * nspin_yrs + nspin_days
  day    <- seq(1,nday,1)
  year   <- (day-litadd_day)/365
  doy    <- 1
    
  LIT    <- array(NA, dim = c(2,nday))
  LITBAG <- array(NA, dim = c(2,nday))
  MIC    <- array(NA, dim = c(2,nday))
  SOM    <- array(NA, dim = c(3,nday))
  LITbag_1 = 0. # initial litter bag 
  LITbag_2 = 0

  sim_year = 0
  i = 1

  for (d in 1:nday)  { 
    # For recalc of Tpars from daily forcing data
    if (all(!is.na(dailyInput))) {
      
      # Loop daily input (i.e., repeat 365 rows in daily input)
      if(loop_dailyInput){
        input_doy = (d - 1) %% nrow(dailyInput) + 1        
      } else {input_doy <- d}

      # Set daily Tpars from "dailyInput" dataframe (add in the function arguments)
        Tpars_mod = calc_Tpars_Conly(ANPP = dailyInput$ANPP[input_doy]/2,
                             fCLAY = dailyInput$CLAY[input_doy]/100, 
                             TSOI = dailyInput$TSOI[input_doy], 
                             MAT = dailyInput$MAT[input_doy], 
                             LIG_N = dailyInput$LIG_N[input_doy],
                             CN = dailyInput$CN[input_doy],  # Only needed if LIG_N not supplied 
                             LIG = dailyInput$LIG[input_doy],  # Only needed if LIG_N not supplied 
                             theta_liq = dailyInput$GWC[input_doy]/100, 
                             theta_frzn = 0,  # Change to column name if frozen water content is available
                             W_SCALAR = dailyInput$W_SCALAR[input_doy],
                             LITFALL=dailyInput$LITFALL[input_doy]/2) 
    } else {
      # Use ss Tpars (i.e., same forcing variables for each sim day)
      Tpars_mod = Tpars
    }  
   
    # Run simulation at hourly timestep
    for (h in 1:24)   {
      Ty <- c( LIT_1 = MIMfwd[1], LIT_2 = MIMfwd[2], 
                 MIC_1 = MIMfwd[3], MIC_2 = MIMfwd[4], 
                 SOM_1 = MIMfwd[5], SOM_2 = MIMfwd[6], 
                 SOM_3 = MIMfwd[7])
      
      # Run MIMICS simulation step
      step = RXEQ(t=NA, y=Ty, pars=Tpars_mod)
      
      # Litter bag decomp simulation step
      LITbag  <- rep(NA,4)
      LITbag[1] <- Ty[[3]] * Tpars_mod$VMAX[1] * LITbag_1 / (Tpars_mod$KM[1] + Ty[[3]])   #MIC_1 mineralization of METABOLIC litter
      LITbag[2] <- Ty[[3]] * Tpars_mod$VMAX[2] * LITbag_2 / (Tpars_mod$KM[2] + Ty[[3]])   #MIC_1 mineralization of STRUC litter
      LITbag[3] <- Ty[[4]] * Tpars_mod$VMAX[4] * LITbag_1 / (Tpars_mod$KM[4] + Ty[[4]])   #mineralization of MET litter
      LITbag[4] <- Ty[[4]] * Tpars_mod$VMAX[5] * LITbag_2 / (Tpars_mod$KM[5] + Ty[[4]])   #mineralization of SRUCTURAL litter

      # Update MIMICS pools
      #---------------------------------------------
      MIMfwd = MIMfwd + unlist(step[[1]]) # MIMICS pools
      
      LITbag_1 <- LITbag_1 - LITbag[1] - LITbag[3] # Litter bag pools
      LITbag_2 <- LITbag_2 - LITbag[2] - LITbag[4] 
      

      # add litter on correct day
      if (d == litadd_day)   {
        if (h == 24)  {
          LITbag_1 <- LITbag_1 + as.numeric(litBAG[3]) # Metabolic pool add
          LITbag_2 <- LITbag_2 + as.numeric(litBAG[4]) # Structural pool add
        }
      }
          
      #write out daily results
      if (h == 24) {
        LIT[1,d] <- MIMfwd[1]
        LIT[2,d] <- MIMfwd[2]
        LITBAG[1,d]  <- LITbag_1
        LITBAG[2,d]  <- LITbag_2
        MIC[1,d] <- MIMfwd[3]
        MIC[2,d] <- MIMfwd[4]
        SOM[1,d] <- MIMfwd[5]
        SOM[2,d] <- MIMfwd[6]
        SOM[3,d] <- MIMfwd[7]
          
        #advance day of year counter
        if (doy == 365) {
          doy <- 1
          sim_year = sim_year + 1
          if(verbose){
              print(paste0("Finished MIMICS simulation year ", sim_year))
          }
        } else {
          doy <- doy + 1
        }                         
      } # close h=24 loop	   						
    }		#close hour loop
  }		#close daily loop

  # Compile and format output
  LITBAG_out <- rbind(as.data.frame(LITBAG), 
                      as.data.frame(LIT), 
                      as.data.frame(MIC),
                      as.data.frame(SOM))
  
  LITBAG_out <- LITBAG_out * depth * 1e4 / 1e6 # Soil depth and unit conversion 
  LITBAG_out <- as.data.frame(t(LITBAG_out))
  colnames(LITBAG_out) <- c("LITBAGm", "LITBAGs", "LITm", "LITs", "MICr", "MICk", "SOMp", "SOMc", "SOMa")
  LITBAG_out <- cbind(data.frame(SITE = forcing_df$SITE[1],
                                 Litter_Type = as.character(litBAG[1]),
                                 DAY=seq(1:nrow(LITBAG_out))), LITBAG_out)
  
  return(LITBAG_out)
}
