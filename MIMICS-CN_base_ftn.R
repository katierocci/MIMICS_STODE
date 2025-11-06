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
CN_K        <<- 10
Nleak <<- 0.2 #0.2
densDep <<- 1 #beta off 

###########################################
# MIMICS single point function
###########################################

#LTER input data
# data = read.csv("LTER_SITE_1.csv")
# df <- data[1,] 
# 
# Site = df$Site
# ANPP = df$ANPP/2 
# TSOI = df$MAT
# fCLAY = df$CLAY2/100
# lig_N = (df$LIG/100)/(1/(df$CN/2.5))
# fMET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * lig_N) 
# CN_s        <<- (df$CN-CN_m*fMET)/(1-fMET)

#make your own input data
Site = 'ParamTraits_test'
ANPP = 500 
TSOI = 20
fCLAY = 0.33
lig_N = 25
fMET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * lig_N) 
CN_s        <<- (50-CN_m*fMET)/(1-fMET)

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

test  <<- stode(y = Ty, time = 1e7, fun = CN_RXEQ, parms = Tpars, positive = TRUE)

### Calc and get MIMICS output 
MIMLIT    <- (test[[1]][[1]]+test[[1]][[2]])  * depth *1e4 / 1e6 #convert kgC/m2 from mgC/cm3 (0-30 cm) 
MIMMIC    <- (test[[1]][[3]]+test[[1]][[4]])  * depth *1e4 / 1e6
MIM_CO    <-  test[[1]][[3]]/test[[1]][[4]]
MIMSOC    <- sum(test[[1]])  * depth *1e4 / 1e6   

# C only version results:
# LIT_1       LIT_2       MIC_1       MIC_2       SOM_1       SOM_2       SOM_3 
# 1.656909890 4.669793042 0.025965813 0.004162176 0.772933300 3.266147911 3.854644406 

test[[1]]
output <- t(as.data.frame(test[[1]]))
attributes(test)


##############################
# ODE SOLVER
##############################
 
# end = 3e6
# step = 1e2
# time_steps = seq(1,end,step)
# n_steps = length(time_steps)
# 
# start_time <- Sys.time()
# 
# t_sol_ss = ode(
#   y = Ty,
#   times = time_steps,
#   func = CN_RXEQ,
#   parms = c(),
#   method = "ode45",
#   maxsteps = 1000000)
# 
# print(paste0("Task time: ", Sys.time() - start_time))
# gc()
# 
# sol <- as.data.frame(t_sol_ss)
# 
# ggplot(sol, aes(x=time, y=SOM_1)) + geom_point() + geom_line()


########################################
# SS forward
#######################################
source("calc_Tpars.R")

# Grab steady state pools from stode output above
ss = test[[1]]

# Create dataframe to store MIMICS pools over timesteps
 # Init first row = steady state
MIMfwd = t(as.data.frame(ss)) 

# Set number of days to sim forward
sim_days = 365*2

# Hourly model loop
for(i in 2:(sim_days*24)){ # interval = hour
  #print(i) #progress tracker
  
  # Recalc Tpars here, calls ftn in calc_Tpars.R
  #---------------------------------------------
  # e.g. 1 C warming over simulation period
  total_warming <- 0
  hr_warming <- total_warming/(sim_days*24) 
  TSOI_adj = -7+(hr_warming*(i-1))
  Tpars_mod = calc_Tpars(TSOI = TSOI_adj, ANPP = 141, CLAY = 5, CN =36.49635, LIG = 16.6) #>> Same example site input as used for stode
  
  # Update MIMICS pools
  #---------------------------------------------
  step = CN_RXEQ(t=NA, y=MIMfwd[i-1,], pars=Tpars_mod)
  pools_update = MIMfwd[i-1,] + unlist(step)
  MIMfwd = rbind(MIMfwd, t(as.data.frame(pools_update)))
}

# Steady state forward hourly data
ss_forward_output = as.data.frame(MIMfwd)

# example pool change over simulation period
ggplot(ss_forward_output, aes(y=SOM_3_N, x=1:(sim_days*24)/24)) + geom_point() +
  xlab("Days") +
  theme_bw()

