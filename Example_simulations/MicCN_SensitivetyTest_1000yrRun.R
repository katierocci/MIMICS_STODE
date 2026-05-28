## Set working drive
#setwd("C:/github/MIMICS_STODE")

#Libraries
library(rootSolve)
library(boot)
library(dplyr)
library(purrr)
library(ggplot2)
library(Metrics)
library(deSolve)

Start.time = Sys.time() # took 11 hours!!!!

# Bring in RXEQ function
source("CN_RXEQ.R")
#source("C:/github/MIMICS_HiRes/MIMICS_ftns/RXEQ_ftn.R")

########################################
# Set MIMICS default parameters
########################################
Vslope  <- rep(0.063, 6)
Vint    <- rep(5.47, 6)
aV      <- rep(0.000008, 6)  
Kslope  <- rep(c(0.025, 0.035, 0.025),2)
Kint    <- rep(3.19, 6)
aK      <- rep(10, 6)
vMOD    <- c(10, 2, 10, 3, 3, 2)
kMOD    <- c(8, 2, 4, 2, 4, 6)
KO      <- c(6, 6)
CUE     <- c(0.55, 0.25, 0.75, 0.35)
tau_r   <- c(0.00052, 0.3)
tau_K   <- c(0.00024, 0.1)
Tau_MOD <- c(100, 0.8, 1.2, 2)
fPHYS_r <- c(0.3, 1.3)
fPHYS_K <- c(0.2, 0.8)
fCHEM_r <- c(0.1, -3, 1)
fCHEM_K <- c(0.3, -3, 1)
fSOM_p  <- c(0.000015, -1.5)
PHYS_scalar <- c(2, -2, NA, NA, NA, NA)
FI      <- c(0.05, 0.05)
fmet_p <- c(1, 0.85, 0.013)
depth <- 30 # set soil depth

#Set default multipliers
Tau_MULT = 1
desorb_MULT = 1
fPHYS_MULT = 1

# N parameters
NUE <<- rep(0.85, 4) #was set to 0.9 - 0.85 from Kyker-Snowman
CN_m        <<- 15
CN_r        <<- 6
print(CN_r)
CN_K        <<- 10
print(CN_K)
Nleak <<- 0.2 #0.2
densDep <<- 1 #beta off 

###############################################
#run stode to get starting point for model
###############################################

#LTER input data
data = read.csv("LTER_SITE_1.csv")
df <- data[c(2,6,14),]

#loop for sites
all_site_output <- list()
for (lter in 1:nrow(df)) {
  df_site <- df[lter, ]
  
  #print(df_site)
  
  cat("Running site:", df_site$Site, "\n")

  #
  Site = df_site$Site
  ANPP = df_site$ANPP/2
  TSOI = df_site$MAT
  fCLAY = df_site$CLAY2/100
  lig_N = df_site$LIG/df_site$N
  fMET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * lig_N)
  CN_s        <<- (df_site$CN-CN_m*fMET)/(1-fMET)
  # df2 <- df %>% mutate(lig_N = LIG/N, fMET = fmet_p[1] * (fmet_p[2] - fmet_p[3] * lig_N)) %>%
  #   dplyr::select(Site, ANPP, MAT, SAND2, CLAY2, LIG, CN, lig_N, fMET)
  # write.csv(df2, "LTER_SensitivitySites.csv")
  
  #make your own input data
  # Site = 'ParamTraits_test'
  # ANPP = 500
  # TSOI = 20
  # fCLAY = 0.33
  # lig_N = 25
  # fMET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * lig_N)
  # CN_s        <<- (50-CN_m*fMET)/(1-fMET)
  
  ############################################################
  # MIMICS MODEL CODE STARTS HERE
  ############################################################
    
  # Calc litter input rate
  EST_LIT <- (ANPP / (365*24)) * 1e3 / 1e4
  #print(EST_LIT)# gC/m2/h (from gC/m2/y) then mgC/cm2/h(from gC/m2/h) 
  
  # ------------ caclulate parameters ---------------
  Vmax     <- exp(TSOI * Vslope + Vint) * aV 
  Km       <- exp(TSOI * Kslope + Kint) * aK
  
  #ANPP strongly correlated with MAP
  Tau_MOD1 <- sqrt(ANPP/Tau_MOD[1])         
  Tau_MOD2 <- Tau_MOD[4]                        
  Tau_MOD1[Tau_MOD1 < Tau_MOD[2]] <- Tau_MOD[2]
  Tau_MOD1[Tau_MOD1 > Tau_MOD[3]] <- Tau_MOD[3] 
  
  tau <- c(tau_r[1]*exp(tau_r[2]*fMET), 
           tau_K[1]*exp(tau_K[2]*fMET))   
  tau <- tau * Tau_MOD1 * Tau_MOD2 * Tau_MULT 
  
  fPHYS    <- c(fPHYS_r[1] * exp(fPHYS_r[2]*fCLAY), 
                fPHYS_K[1] * exp(fPHYS_K[2]*fCLAY)) 	            
  fCHEM    <- c(fCHEM_r[1] * exp(fCHEM_r[2]*fMET) * fCHEM_r[3], 
                fCHEM_K[1] * exp(fCHEM_K[2]*fMET) * fCHEM_K[3]) 	
  fAVAI    <- 1 - (fPHYS + fCHEM)
  
  desorb   <- fSOM_p[1] * exp(fSOM_p[2]*(fCLAY))                  
  
  desorb <- desorb * desorb_MULT
  fPHYS <- fPHYS * fPHYS_MULT
  
  pSCALAR  <- PHYS_scalar[1] * exp(PHYS_scalar[2]*(sqrt(fCLAY)))  #Scalar for texture effects on SOMp
  
  v_MOD    <- vMOD  
  k_MOD    <- kMOD 
  k_MOD[3] <- k_MOD[3] * pSCALAR    
  k_MOD[6] <- k_MOD[6] * pSCALAR    
  
  VMAX     <- Vmax * v_MOD 
  KM       <- Km / k_MOD
  
  #----------initialize pools---------------
  I       <- array(NA, dim=2)             
  I[1]    <- (EST_LIT / depth) * fMET     
  I[2]    <- (EST_LIT / depth) * (1-fMET)
  Inputs <<- I
  lit     <- I   
  mic     <- I  
  som     <- rep(NA, 3) 
  som[1]  <- I[1]
  som[2]  <- I[2]
  som[3]  <- I[1] 
  LITmin  <- rep(NA, dim=4)
  MICtrn  <- c(NA,NA,NA,NA,NA,NA)
  SOMmin  <- rep(NA, dim=2)
  DEsorb  <- rep(NA, dim=1)
  OXIDAT  <- rep(NA, dim=1)
  
  LIT_1_N  <<- 1e-4
  LIT_2_N  <<- 1e-4
  MIC_1_N  <<- 1e-4
  MIC_2_N  <<- 1e-4
  SOM_1_N  <<- 1e-4
  SOM_2_N  <<- 1e-4
  SOM_3_N  <<- 1e-4
  DIN      <<- 1e-4
  
  LITminN   <<- array(NA, dim=4)
  MICtrnN   <<- array(NA, dim=6)
  SOMminN   <<- array(NA, dim=2)
  DEsorbN   <<- array(NA, dim=1)
  OXIDATN   <<- array(NA, dim=1)
  
  DINup     <<- array(NA, dim=2)
  Overflow  <<- array(NA, dim=2)
  Nspill    <<- array(NA, dim=2)
  CNup      <<- array(NA, dim=2)
  upMIC_1   <<-  array(NA, dim=1)
  upMIC_1_N <<-  array(NA, dim=1)
  upMIC_2   <<-  array(NA, dim=1)
  upMIC_2_N <<-  array(NA, dim=1)
  
  
  Tpars <<- c( Inputs = I, VMAX = VMAX, KM = KM, CUE = CUE, 
               fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
               tau = tau, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
               desorb = desorb, DEsorb = DEsorb, OXIDAT = OXIDAT, KO = KO,
               LITminN = LITminN, SOMminN = SOMminN, MICtrnN = MICtrnN,
               DEsorbN = DEsorbN, OXIDATN = OXIDATN, densDep=densDep,
               CNup=CNup, DINup=DINup, Nspill=Nspill, Overflow=Overflow, 
               upMIC_1=upMIC_1, upMIC_1_N=upMIC_1_N,
               upMIC_2=upMIC_2, upMIC_2_N=upMIC_2_N,
               NUE=NUE, CN_m=CN_m, CN_s=CN_s, CN_r=CN_r, CN_K=CN_K, Nleak=Nleak)
  
  Ty    <<- c(LIT_1 = lit[1], LIT_2 = lit[2], 
              MIC_1 = mic[1], MIC_2 = mic[2], 
              SOM_1 = som[1], SOM_2 = som[2], SOM_3 = som[3],
              LIT_1_N = LIT_1_N, LIT_2_N = LIT_2_N, 
              MIC_1_N = MIC_1_N, MIC_2_N = MIC_2_N, 
              SOM_1_N = SOM_1_N, SOM_2_N = SOM_2_N, SOM_3_N = SOM_3_N,
              DIN = DIN)
  
  MIMss  <<- stode(y = Ty, time = 1e7, fun = CN_RXEQ, parms = Tpars, positive = TRUE)
  
  ############################
  #transient runs
  ###########################
  
  
  #######
  #C:N applied to be averaged out for copios and oligos
  #######
  
  #setup for dataframe and index collection
  MIM_output_test <- data_frame()
  
  #average microbial C:N - range plus MIMICS default (8)
  MicCN2 = c(6, 10, 13, 17, 45, 8)
  
  for (b in MicCN2) {
    print(b)
    ############
    # Update MIMICS C:N for r and k
    ##############
    CN_r        <<- 0.75*b #6
    #print(CN_r)
    CN_K        <<- 1.25*b #10
    #print(CN_K)
    
    ###################
    # MIMICS single point function
    ###################
    
    #make your own input data
    Site = df_site$Site
    ANPP = df_site$ANPP/2
    TSOI = df_site$MAT
    fCLAY = df_site$CLAY2/100
    lig_N = df_site$LIG/df_site$N
    fMET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * lig_N)
    CN_s        <- (df_site$CN-CN_m*fMET)/(1-fMET)
    nday = 1000*365 #number of days to run simulation for
    
    ############################################################
    # MIMICS MODEL CODE STARTS HERE
    ############################################################
    
    # Calc litter input rate
    EST_LIT <- (ANPP / (365*24)) * 1e3 / 1e4
    #print(EST_LIT)# gC/m2/h (from gC/m2/y) then mgC/cm2/h(from gC/m2/h) 
    
    # ------------ calculate parameters ---------------
    Vmax     <- exp(TSOI * Vslope + Vint) * aV 
    Km       <- exp(TSOI * Kslope + Kint) * aK
    
    #ANPP strongly correlated with MAP
    Tau_MOD1 <- sqrt(ANPP/Tau_MOD[1])         
    Tau_MOD2 <- Tau_MOD[4]                        
    Tau_MOD1[Tau_MOD1 < Tau_MOD[2]] <- Tau_MOD[2]
    Tau_MOD1[Tau_MOD1 > Tau_MOD[3]] <- Tau_MOD[3] 
    
    tau <- c(tau_r[1]*exp(tau_r[2]*fMET), 
             tau_K[1]*exp(tau_K[2]*fMET))   
    tau <- tau * Tau_MOD1 * Tau_MOD2 * Tau_MULT 
    #print(tau) #0.0015 per hour and 0.00061 per hour -> 13 and 5.3 yr-1 0.08 and 0.19 residence time 
    
    fPHYS    <- c(fPHYS_r[1] * exp(fPHYS_r[2]*fCLAY), 
                  fPHYS_K[1] * exp(fPHYS_K[2]*fCLAY)) 	            
    fCHEM    <- c(fCHEM_r[1] * exp(fCHEM_r[2]*fMET) * fCHEM_r[3], 
                  fCHEM_K[1] * exp(fCHEM_K[2]*fMET) * fCHEM_K[3]) 	
    fAVAI    <- 1 - (fPHYS + fCHEM)
    
    desorb   <- fSOM_p[1] * exp(fSOM_p[2]*(fCLAY))                  
    
    desorb <- desorb * desorb_MULT
    fPHYS <- fPHYS * fPHYS_MULT
    
    pSCALAR  <- PHYS_scalar[1] * exp(PHYS_scalar[2]*(sqrt(fCLAY)))  #Scalar for texture effects on SOMp
    
    v_MOD    <- vMOD  
    k_MOD    <- kMOD 
    k_MOD[3] <- k_MOD[3] * pSCALAR    
    k_MOD[6] <- k_MOD[6] * pSCALAR    
    
    VMAX     <- Vmax * v_MOD 
    KM       <- Km / k_MOD
    
    #----------initialize pools---------------
    MIMfwd = MIMss[[1]] #using stode run above to initiate
    MIMfwd = (MIMfwd / depth) / (1e4 / 1e6)  #convert from gC m-2 to mgC cm-3
    
    # Get Tpars from ss simulation  
    #Tpars = MIMss[[2]]
    
    #Init arrays to store daily output data
    day    <- seq(1,nday,1)
    year   <- day/365
    doy    <- 1
    
    out_i <- 0L
    n_out <- nday %/% (10 * 365)
    LIT   <- array(NA, dim = c(2, n_out))
    MIC   <- array(NA, dim = c(2, n_out))
    SOM   <- array(NA, dim = c(3, n_out))
    LIT_N <- array(NA, dim = c(2, n_out))
    MIC_N <- array(NA, dim = c(2, n_out))
    SOM_N <- array(NA, dim = c(3, n_out))

    
    sim_year = 0
    i = 1
    
    I       <- array(NA, dim=2)             
    I[1]    <- (EST_LIT / depth) * fMET     
    I[2]    <- (EST_LIT / depth) * (1-fMET)
    
    for (d in 1:nday)  { 
      Tpars <<- c( Inputs = I, VMAX = VMAX, KM = KM, CUE = CUE, 
                 fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
                 tau = tau, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
                 desorb = desorb, DEsorb = DEsorb, OXIDAT = OXIDAT, KO = KO,
                 LITminN = LITminN, SOMminN = SOMminN, MICtrnN = MICtrnN,
                 DEsorbN = DEsorbN, OXIDATN = OXIDATN, densDep=densDep,
                 CNup=CNup, DINup=DINup, Nspill=Nspill, Overflow=Overflow, 
                 upMIC_1=upMIC_1, upMIC_1_N=upMIC_1_N,
                 upMIC_2=upMIC_2, upMIC_2_N=upMIC_2_N,
                 NUE=NUE, CN_m=CN_m, CN_s=CN_s, CN_r=CN_r, CN_K=CN_K, Nleak=Nleak)
    
    # Run simulation at hourly timestep
      for (h in 1:24)   {
        #start from steady state
        Ty <- c(LIT_1 = MIMfwd[[1]], LIT_2 = MIMfwd[[2]],
               MIC_1 = MIMfwd[[3]], MIC_2 = MIMfwd[[4]],
               SOM_1 = MIMfwd[[5]], SOM_2 = MIMfwd[[6]],
               SOM_3 = MIMfwd[[7]],
               LIT_1_N = MIMfwd[[8]], LIT_2_N = MIMfwd[[9]],
               MIC_1_N = MIMfwd[[10]], MIC_2_N = MIMfwd[[11]],
               SOM_1_N = MIMfwd[[12]], SOM_2_N = MIMfwd[[13]],
               SOM_3_N = MIMfwd[[14]], DIN = MIMfwd[[15]])
        #start from 1
      #   Ty <- c(LIT_1 = 1, LIT_2 = 1, 
      #           MIC_1 = 1, MIC_2 = 1, 
      #           SOM_1 =1, SOM_2 = 1, 
      #           SOM_3 = 1,
      #           LIT_1_N = 1, LIT_2_N = 1, 
      #           MIC_1_N = 1, MIC_2_N = 1, 
      #           SOM_1_N = 1, SOM_2_N = 1, 
      #           SOM_3_N = 1, DIN = 1)
      # #print(Ty)
        
      # Run MIMICS simulation step
        step = CN_RXEQ(t=NA, y=Ty, pars=Tpars)
      
      
      
      # Update MIMICS pools
      #---------------------------------------------
        MIMfwd = MIMfwd + unlist(step[[1]]) # MIMICS pools
      
      #write out 10-year results
        if (h == 24) {
          
          # only save every 10 years
          if (d %% (10 * 365) == 0) {
            out_i          <- out_i + 1L
            LIT[1, out_i]   <- MIMfwd[1]
            LIT[2, out_i]   <- MIMfwd[2]
            MIC[1, out_i]   <- MIMfwd[3]
            MIC[2, out_i]   <- MIMfwd[4]
            SOM[1, out_i]   <- MIMfwd[5]
            SOM[2, out_i]   <- MIMfwd[6]
            SOM[3, out_i]   <- MIMfwd[7]
            LIT_N[1, out_i] <- MIMfwd[8]
            LIT_N[2, out_i] <- MIMfwd[9]
            MIC_N[1, out_i] <- MIMfwd[10]
            MIC_N[2, out_i] <- MIMfwd[11]
            SOM_N[1, out_i] <- MIMfwd[12]
            SOM_N[2, out_i] <- MIMfwd[13]
            SOM_N[3, out_i] <- MIMfwd[14]
          }
          
          if (doy == 365) {
            doy <- 1
            sim_year = sim_year + 1
            print(paste0("Finished MIMICS simulation year ", sim_year))
          } else {
            doy <- doy + 1
          }
        } # close h=24 loop   						
      }		#close hour loop
    }		#close daily loop
    
    MIM_out <- rbind(as.data.frame(LIT), 
                        as.data.frame(MIC),
                        as.data.frame(SOM),
                     as.data.frame(LIT_N), 
                     as.data.frame(MIC_N),
                     as.data.frame(SOM_N))
    
    MIM_out <- MIM_out * depth * 1e4 / 1e6  #mg/cm3 converted kg/m2
    MIM_out <- as.data.frame(t(MIM_out))
    colnames(MIM_out) <- c("LITm", "LITs", "MICr", "MICk", "SOMp", "SOMc", "SOMa","LITm_N", "LITs_N", "MICr_N", "MICk_N", "SOMp_N", "SOMc_N", "SOMa_N")
    #MIM_out2 <- cbind(data.frame(CN_val = b, Site = df_site$Site, DAY=seq(1:nrow(MIM_out))), MIM_out)
    MIM_out2 <- cbind(data.frame(CN_val = b, 
                                 Site   = df_site$Site, 
                                 YEAR   = seq(10, nday/365, by = 10),
                                 DAY    = seq(10*365, nday, by = 10*365)), MIM_out)
    MIM_output_test <- rbind(MIM_output_test, MIM_out2)
  } #close C:N loop
  all_site_output[[lter]] <- MIM_output_test # store each site's result
}#close site loop
final_output <- do.call(rbind, all_site_output) #working!

#write out CSV
write.csv(final_output,"MIM_Sens_CN_051526.csv")
 
# MIM_out_CN <- MIM_output_test %>% mutate(soil_CN = (SOMc+SOMp+SOMa+MICr+MICk)/(SOMc_N+SOMp_N+SOMa_N+MICr_N+MICk_N))
# write.csv(MIM_out_CN, "MIM_out_CN_04212026.csv")
# 
# MIM_out_CN <- read.csv("MIM_out_CN.csv")
# ggplot(MIM_out_CN) + geom_line(aes(x=DAY, y=soil_CN, group = as.factor(CN_val), color=as.factor(CN_val))) +theme_bw(base_size = 16) 
# ggplot(MIM_out_CN) + geom_line(aes(x=DAY, y=(SOMc+SOMp+SOMa+MICr+MICk), group = as.factor(CN_val), color=as.factor(CN_val))) +theme_bw(base_size = 16) 
# ggplot(MIM_out_CN) + geom_line(aes(x=DAY, y=(MICr/MICk), group = as.factor(CN_val), color=as.factor(CN_val))) +theme_bw(base_size = 16) 
# 
# MIM_out_CN_end <- MIM_out_CN %>% filter(DAY==365000)
# write.csv(MIM_out_CN_end, "MIM_out_CN_1000.csv")
#  
# 
# #soil C
# ggplot(MIM_out_CN_end) + geom_point(aes(x=CN_val, y=(MICr+MICk+SOMa+SOMc+SOMp)), size=4) + theme_bw(base_size = 16) + ylab("Soil C")
# #soil N
# ggplot(MIM_out_CN_end) + geom_point(aes(x=CN_val, y=(SOMa_N+SOMc_N+SOMp_N)), size=4) + theme_bw(base_size = 16) + ylab("Soil N in all three pools")
# #soil C:N
# ggplot(MIM_out_CN_end) + geom_point(aes(x=CN_val, y=soil_CN), size=4) + theme_bw(base_size = 16)+ ylab("Soil C:N")
# #MICr:MICk
# ggplot(MIM_out_CN_end) + geom_point(aes(x=CN_val, y=(MICr/MICk)), size=4) + theme_bw(base_size = 16)+ ylab("MICr:MICk")

####################################################################################################################


####
#Carbon USe Efficiency (CUE)
####

########################################
# Set MIMICS default parameters
########################################
Vslope  <- rep(0.063, 6)
Vint    <- rep(5.47, 6)
aV      <- rep(0.000008, 6)  
Kslope  <- rep(c(0.025, 0.035, 0.025),2)
Kint    <- rep(3.19, 6)
aK      <- rep(10, 6)
vMOD    <- c(10, 2, 10, 3, 3, 2)
kMOD    <- c(8, 2, 4, 2, 4, 6)
KO      <- c(6, 6)
CUE     <- c(0.55, 0.25, 0.75, 0.35)
tau_r   <- c(0.00052, 0.3)
tau_K   <- c(0.00024, 0.1)
Tau_MOD <- c(100, 0.8, 1.2, 2)
fPHYS_r <- c(0.3, 1.3)
fPHYS_K <- c(0.2, 0.8)
fCHEM_r <- c(0.1, -3, 1)
fCHEM_K <- c(0.3, -3, 1)
fSOM_p  <- c(0.000015, -1.5)
PHYS_scalar <- c(2, -2, NA, NA, NA, NA)
FI      <- c(0.05, 0.05)
fmet_p <- c(1, 0.85, 0.013)
depth <- 30 # set soil depth

#Set default multipliers
Tau_MULT = 1
desorb_MULT = 1
fPHYS_MULT = 1

# N parameters
NUE <<- rep(0.85, 4) #was set to 0.9 - 0.85 from Kyker-Snowman
CN_m        <<- 15
CN_r        <<- 6
print(CN_r)
CN_K        <<- 10
print(CN_K)
Nleak <<- 0.2 #0.2
densDep <<- 1 #beta off 

###############################################
#run stode to get starting point for model
###############################################

#LTER input data
data = read.csv("LTER_SITE_1.csv")
df <- data[c(2,6,14),]

#loop for sites
all_site_output_CUE <- list()
for (lter in 1:nrow(df)) {
  df_site <- df[lter, ]
  
  #print(df_site)
  
  cat("Running site:", df_site$Site, "\n")
  
  #
  Site = df_site$Site
  ANPP = df_site$ANPP/2
  TSOI = df_site$MAT
  fCLAY = df_site$CLAY2/100
  lig_N = df_site$LIG/df_site$N
  fMET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * lig_N)
  CN_s        <<- (df_site$CN-CN_m*fMET)/(1-fMET)
  # df2 <- df %>% mutate(lig_N = LIG/N, fMET = fmet_p[1] * (fmet_p[2] - fmet_p[3] * lig_N)) %>%
  #   dplyr::select(Site, ANPP, MAT, SAND2, CLAY2, LIG, CN, lig_N, fMET)
  # write.csv(df2, "LTER_SensitivitySites.csv")
  
  #make your own input data
  # Site = 'ParamTraits_test'
  # ANPP = 500
  # TSOI = 20
  # fCLAY = 0.33
  # lig_N = 25
  # fMET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * lig_N)
  # CN_s        <<- (50-CN_m*fMET)/(1-fMET)
  
  ############################################################
  # MIMICS MODEL CODE STARTS HERE
  ############################################################
  
  # Calc litter input rate
  EST_LIT <- (ANPP / (365*24)) * 1e3 / 1e4
  #print(EST_LIT)# gC/m2/h (from gC/m2/y) then mgC/cm2/h(from gC/m2/h) 
  
  # ------------ caclulate parameters ---------------
  Vmax     <- exp(TSOI * Vslope + Vint) * aV 
  Km       <- exp(TSOI * Kslope + Kint) * aK
  
  #ANPP strongly correlated with MAP
  Tau_MOD1 <- sqrt(ANPP/Tau_MOD[1])         
  Tau_MOD2 <- Tau_MOD[4]                        
  Tau_MOD1[Tau_MOD1 < Tau_MOD[2]] <- Tau_MOD[2]
  Tau_MOD1[Tau_MOD1 > Tau_MOD[3]] <- Tau_MOD[3] 
  
  tau <- c(tau_r[1]*exp(tau_r[2]*fMET), 
           tau_K[1]*exp(tau_K[2]*fMET))   
  tau <- tau * Tau_MOD1 * Tau_MOD2 * Tau_MULT 
  
  fPHYS    <- c(fPHYS_r[1] * exp(fPHYS_r[2]*fCLAY), 
                fPHYS_K[1] * exp(fPHYS_K[2]*fCLAY)) 	            
  fCHEM    <- c(fCHEM_r[1] * exp(fCHEM_r[2]*fMET) * fCHEM_r[3], 
                fCHEM_K[1] * exp(fCHEM_K[2]*fMET) * fCHEM_K[3]) 	
  fAVAI    <- 1 - (fPHYS + fCHEM)
  
  desorb   <- fSOM_p[1] * exp(fSOM_p[2]*(fCLAY))                  
  
  desorb <- desorb * desorb_MULT
  fPHYS <- fPHYS * fPHYS_MULT
  
  pSCALAR  <- PHYS_scalar[1] * exp(PHYS_scalar[2]*(sqrt(fCLAY)))  #Scalar for texture effects on SOMp
  
  v_MOD    <- vMOD  
  k_MOD    <- kMOD 
  k_MOD[3] <- k_MOD[3] * pSCALAR    
  k_MOD[6] <- k_MOD[6] * pSCALAR    
  
  VMAX     <- Vmax * v_MOD 
  KM       <- Km / k_MOD
  
  #----------initialize pools---------------
  I       <- array(NA, dim=2)             
  I[1]    <- (EST_LIT / depth) * fMET     
  I[2]    <- (EST_LIT / depth) * (1-fMET)
  Inputs <<- I
  lit     <- I   
  mic     <- I  
  som     <- rep(NA, 3) 
  som[1]  <- I[1]
  som[2]  <- I[2]
  som[3]  <- I[1] 
  LITmin  <- rep(NA, dim=4)
  MICtrn  <- c(NA,NA,NA,NA,NA,NA)
  SOMmin  <- rep(NA, dim=2)
  DEsorb  <- rep(NA, dim=1)
  OXIDAT  <- rep(NA, dim=1)
  
  LIT_1_N  <<- 1e-4
  LIT_2_N  <<- 1e-4
  MIC_1_N  <<- 1e-4
  MIC_2_N  <<- 1e-4
  SOM_1_N  <<- 1e-4
  SOM_2_N  <<- 1e-4
  SOM_3_N  <<- 1e-4
  DIN      <<- 1e-4
  
  LITminN   <<- array(NA, dim=4)
  MICtrnN   <<- array(NA, dim=6)
  SOMminN   <<- array(NA, dim=2)
  DEsorbN   <<- array(NA, dim=1)
  OXIDATN   <<- array(NA, dim=1)
  
  DINup     <<- array(NA, dim=2)
  Overflow  <<- array(NA, dim=2)
  Nspill    <<- array(NA, dim=2)
  CNup      <<- array(NA, dim=2)
  upMIC_1   <<-  array(NA, dim=1)
  upMIC_1_N <<-  array(NA, dim=1)
  upMIC_2   <<-  array(NA, dim=1)
  upMIC_2_N <<-  array(NA, dim=1)
  
  
  Tpars <<- c( Inputs = I, VMAX = VMAX, KM = KM, CUE = CUE, 
               fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
               tau = tau, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
               desorb = desorb, DEsorb = DEsorb, OXIDAT = OXIDAT, KO = KO,
               LITminN = LITminN, SOMminN = SOMminN, MICtrnN = MICtrnN,
               DEsorbN = DEsorbN, OXIDATN = OXIDATN, densDep=densDep,
               CNup=CNup, DINup=DINup, Nspill=Nspill, Overflow=Overflow, 
               upMIC_1=upMIC_1, upMIC_1_N=upMIC_1_N,
               upMIC_2=upMIC_2, upMIC_2_N=upMIC_2_N,
               NUE=NUE, CN_m=CN_m, CN_s=CN_s, CN_r=CN_r, CN_K=CN_K, Nleak=Nleak)
  
  Ty    <<- c(LIT_1 = lit[1], LIT_2 = lit[2], 
              MIC_1 = mic[1], MIC_2 = mic[2], 
              SOM_1 = som[1], SOM_2 = som[2], SOM_3 = som[3],
              LIT_1_N = LIT_1_N, LIT_2_N = LIT_2_N, 
              MIC_1_N = MIC_1_N, MIC_2_N = MIC_2_N, 
              SOM_1_N = SOM_1_N, SOM_2_N = SOM_2_N, SOM_3_N = SOM_3_N,
              DIN = DIN)
  
  MIMss  <<- stode(y = Ty, time = 1e7, fun = CN_RXEQ, parms = Tpars, positive = TRUE)
  
  ############################
  #transient runs
  ###########################
  
  
  #######
  #C:N applied to be averaged out for copios and oligos
  #######
  
  #setup for dataframe and index collection
  MIM_output_CUE <- data_frame()
  
  #CUE defaults
  #0.55, 0.25, 0.75, 0.35
  
  #setup for dataframe and index collection
  MIM_output_CUE <- data_frame()
  
  #average microbial turnover time 
  #varying from residence times of 0.05 to 0.5 years or TTs of 20 to 2 yr-1 or 0.0023 to 0.00023, in amounts of 0.001035
  #plus defaults (0.55, 0.75), plus average of defaults (0.65,0.65)
  CUEs = list(c(0.37, 0.37), c(0.55,0.55), c(0.66,0.66), c(0.72,0.72), c(0.78, 0.78), c(0.55, 0.75), c(0.65, 0.65))
  
  for (b in CUEs) {
    print(b)
    ############
    # Update MIMICS CUE
    ##############
    
    CUE <<- c(b[1], 0.25, b[2], 0.35) #updated for metabolic litter and SOMa specifically
    
    
    ###################
    # MIMICS single point function
    ###################
    
    #make your own input data
    Site = df_site$Site
    ANPP = df_site$ANPP/2
    TSOI = df_site$MAT
    fCLAY = df_site$CLAY2/100
    lig_N = df_site$LIG/df_site$N
    fMET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * lig_N)
    CN_s        <- (df_site$CN-CN_m*fMET)/(1-fMET)
    nday = 1000*365 #number of days to run simulation for
    
    ############################################################
    # MIMICS MODEL CODE STARTS HERE
    ############################################################
    
    # Calc litter input rate
    EST_LIT <- (ANPP / (365*24)) * 1e3 / 1e4
    #print(EST_LIT)# gC/m2/h (from gC/m2/y) then mgC/cm2/h(from gC/m2/h) 
    
    # ------------ calculate parameters ---------------
    Vmax     <- exp(TSOI * Vslope + Vint) * aV 
    Km       <- exp(TSOI * Kslope + Kint) * aK
    
    #ANPP strongly correlated with MAP
    Tau_MOD1 <- sqrt(ANPP/Tau_MOD[1])         
    Tau_MOD2 <- Tau_MOD[4]                        
    Tau_MOD1[Tau_MOD1 < Tau_MOD[2]] <- Tau_MOD[2]
    Tau_MOD1[Tau_MOD1 > Tau_MOD[3]] <- Tau_MOD[3] 
    
    #ANPP strongly correlated with MAP
    Tau_MOD1 <- sqrt(ANPP/Tau_MOD[1])         
    Tau_MOD2 <- Tau_MOD[4]                        
    Tau_MOD1[Tau_MOD1 < Tau_MOD[2]] <- Tau_MOD[2]
    Tau_MOD1[Tau_MOD1 > Tau_MOD[3]] <- Tau_MOD[3] 
    
    tau <- c(tau_r[1]*exp(tau_r[2]*fMET), 
             tau_K[1]*exp(tau_K[2]*fMET))   
    tau <- tau * Tau_MOD1 * Tau_MOD2 * Tau_MULT 
    #print(tau) #0.0015 per hour and 0.00061 per hour -> 13 and 5.3 yr-1 0.08 and 0.19 residence time 
    
    fPHYS    <- c(fPHYS_r[1] * exp(fPHYS_r[2]*fCLAY), 
                  fPHYS_K[1] * exp(fPHYS_K[2]*fCLAY)) 	            
    fCHEM    <- c(fCHEM_r[1] * exp(fCHEM_r[2]*fMET) * fCHEM_r[3], 
                  fCHEM_K[1] * exp(fCHEM_K[2]*fMET) * fCHEM_K[3]) 	
    fAVAI    <- 1 - (fPHYS + fCHEM)
    
    desorb   <- fSOM_p[1] * exp(fSOM_p[2]*(fCLAY))                  
    
    desorb <- desorb * desorb_MULT
    fPHYS <- fPHYS * fPHYS_MULT
    
    pSCALAR  <- PHYS_scalar[1] * exp(PHYS_scalar[2]*(sqrt(fCLAY)))  #Scalar for texture effects on SOMp
    
    v_MOD    <- vMOD  
    k_MOD    <- kMOD 
    k_MOD[3] <- k_MOD[3] * pSCALAR    
    k_MOD[6] <- k_MOD[6] * pSCALAR    
    
    VMAX     <- Vmax * v_MOD 
    KM       <- Km / k_MOD
    
    #----------initialize pools---------------
    MIMfwd = MIMss[[1]] #using stode run above to initiate
    MIMfwd = (MIMfwd / depth) / (1e4 / 1e6)  #convert from gC m-2 to mgC cm-3
    
    # Get Tpars from ss simulation  
    #Tpars = MIMss[[2]]
    
    #Init arrays to store daily output data
    day    <- seq(1,nday,1)
    year   <- day/365
    doy    <- 1
    
    out_i <- 0L
    n_out <- nday %/% (10 * 365)
    LIT   <- array(NA, dim = c(2, n_out))
    MIC   <- array(NA, dim = c(2, n_out))
    SOM   <- array(NA, dim = c(3, n_out))
    LIT_N <- array(NA, dim = c(2, n_out))
    MIC_N <- array(NA, dim = c(2, n_out))
    SOM_N <- array(NA, dim = c(3, n_out))
    
    sim_year = 0
    i = 1
    
    I       <- array(NA, dim=2)             
    I[1]    <- (EST_LIT / depth) * fMET     
    I[2]    <- (EST_LIT / depth) * (1-fMET)
    
    print(CUE)
    
    for (d in 1:nday)  { 
      Tpars <<- c( Inputs = I, VMAX = VMAX, KM = KM, CUE = CUE, 
                   fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
                   tau = tau, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
                   desorb = desorb, DEsorb = DEsorb, OXIDAT = OXIDAT, KO = KO,
                   LITminN = LITminN, SOMminN = SOMminN, MICtrnN = MICtrnN,
                   DEsorbN = DEsorbN, OXIDATN = OXIDATN, densDep=densDep,
                   CNup=CNup, DINup=DINup, Nspill=Nspill, Overflow=Overflow, 
                   upMIC_1=upMIC_1, upMIC_1_N=upMIC_1_N,
                   upMIC_2=upMIC_2, upMIC_2_N=upMIC_2_N,
                   NUE=NUE, CN_m=CN_m, CN_s=CN_s, CN_r=CN_r, CN_K=CN_K, Nleak=Nleak)
      
      # Run simulation at hourly timestep
      for (h in 1:24)   {
        #start from steady state
        Ty <- c(LIT_1 = MIMfwd[[1]], LIT_2 = MIMfwd[[2]],
                MIC_1 = MIMfwd[[3]], MIC_2 = MIMfwd[[4]],
                SOM_1 = MIMfwd[[5]], SOM_2 = MIMfwd[[6]],
                SOM_3 = MIMfwd[[7]],
                LIT_1_N = MIMfwd[[8]], LIT_2_N = MIMfwd[[9]],
                MIC_1_N = MIMfwd[[10]], MIC_2_N = MIMfwd[[11]],
                SOM_1_N = MIMfwd[[12]], SOM_2_N = MIMfwd[[13]],
                SOM_3_N = MIMfwd[[14]], DIN = MIMfwd[[15]])
        #start from 1
        #   Ty <- c(LIT_1 = 1, LIT_2 = 1, 
        #           MIC_1 = 1, MIC_2 = 1, 
        #           SOM_1 =1, SOM_2 = 1, 
        #           SOM_3 = 1,
        #           LIT_1_N = 1, LIT_2_N = 1, 
        #           MIC_1_N = 1, MIC_2_N = 1, 
        #           SOM_1_N = 1, SOM_2_N = 1, 
        #           SOM_3_N = 1, DIN = 1)
        # #print(Ty)
        
        # Run MIMICS simulation step
        step = CN_RXEQ(t=NA, y=Ty, pars=Tpars)
        
        
        
        # Update MIMICS pools
        #---------------------------------------------
        MIMfwd = MIMfwd + unlist(step[[1]]) # MIMICS pools
        
        #write out 10-year results
        if (h == 24) {
          
          # only save every 10 years
          if (d %% (10 * 365) == 0) {
            out_i          <- out_i + 1L
            LIT[1, out_i]   <- MIMfwd[1]
            LIT[2, out_i]   <- MIMfwd[2]
            MIC[1, out_i]   <- MIMfwd[3]
            MIC[2, out_i]   <- MIMfwd[4]
            SOM[1, out_i]   <- MIMfwd[5]
            SOM[2, out_i]   <- MIMfwd[6]
            SOM[3, out_i]   <- MIMfwd[7]
            LIT_N[1, out_i] <- MIMfwd[8]
            LIT_N[2, out_i] <- MIMfwd[9]
            MIC_N[1, out_i] <- MIMfwd[10]
            MIC_N[2, out_i] <- MIMfwd[11]
            SOM_N[1, out_i] <- MIMfwd[12]
            SOM_N[2, out_i] <- MIMfwd[13]
            SOM_N[3, out_i] <- MIMfwd[14]
          }
          
          if (doy == 365) {
            doy <- 1
            sim_year = sim_year + 1
            print(paste0("Finished MIMICS simulation year ", sim_year))
          } else {
            doy <- doy + 1
          }
        } # close h=24 loop   						
      }		#close hour loop
    }		#close daily loop
    
    MIM_out <- rbind(as.data.frame(LIT), 
                     as.data.frame(MIC),
                     as.data.frame(SOM),
                     as.data.frame(LIT_N), 
                     as.data.frame(MIC_N),
                     as.data.frame(SOM_N))
    
    MIM_out <- MIM_out * depth * 1e4 / 1e6  #mg/cm3 converted kg/m2
    MIM_out <- as.data.frame(t(MIM_out))
    colnames(MIM_out) <- c("LITm", "LITs", "MICr", "MICk", "SOMp", "SOMc", "SOMa","LITm_N", "LITs_N", "MICr_N", "MICk_N", "SOMp_N", "SOMc_N", "SOMa_N")
    #MIM_out2 <- cbind(data.frame(CN_val = b, Site = df_site$Site, DAY=seq(1:nrow(MIM_out))), MIM_out)
    print(b)
    MIM_out2 <- cbind(data.frame(CUE_val = b[2], 
                                 Site   = df_site$Site, 
                                 YEAR   = seq(10, nday/365, by = 10),
                                 DAY    = seq(10*365, nday, by = 10*365)), MIM_out)
    MIM_output_CUE <- rbind(MIM_output_CUE, MIM_out2)
  } #close C:N loop
  all_site_output_CUE[[lter]] <- MIM_output_CUE # store each site's result
}#close site loop
final_output_CUE <- do.call(rbind, all_site_output_CUE) #working!

write.csv(final_output_CUE, "MIM_Sens_CUE_052726.csv") #labeling on output is causing issues

End.time = Sys.time()
total.time = Start.time - End.time


#####
#check on test data
#####
final_output_check <- final_output %>% mutate(MICrK = MICr/MICk, MIC = MICr+MICk, SOM = MICr+MICk+SOMa+SOMp+SOMc,
                                              SOM_N = MICr_N+MICk_N+SOMa_N+SOMc_N+SOMp_N, soilCN = SOM/SOM_N,
                                              MIC_SOM = MIC/SOM*100)


