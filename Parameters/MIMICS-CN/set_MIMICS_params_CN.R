########################################
# Set MIMICS default parameters 
# used for MIMICS-CN code in Rocci et al. 2023, SBB
########################################
Vslope  <- rep(0.063, 6)
Vint    <- rep(5.47, 6)
aV      <- rep(0.000008, 6)  
Kslope  <- rep(c(0.025, 0.035, 0.025),2)
Kint    <- rep(3.19, 6)
aK      <- rep(10, 6)
vMOD    <- c(10, 2, 10, 3, 3, 2)
kMOD    <- c(8, 2, 4, 2, 4, 6)
KO      <- c(5, 5) #c(6, 6)
CUE     <- c(0.55, 0.25, 0.75, 0.35)
tau_r   <- c(0.00052, 0.3)
tau_K   <- c(0.00024, 0.1)
Tau_MOD <- c(100, 0.6, 1.3, 1) # Testbed values here c(100, 0.8, 1.2, 2)
fPHYS_r <- c(0.03, 1.3)  # Testbed values here c(0.3, 1.3)
fPHYS_K <- c(0.02, 0.8)  # Testbed values here c(0.2, 0.8)
fCHEM_r <- c(0.1, -3, 3) # Testbed values here c(0.1, -3, 1)
fCHEM_K <- c(0.3, -3, 3) # Testbed values here c(0.3, -3, 1)
fSOM_p  <- c(1.05e-6, -2) # testbed values here c(0.000015, -1.5)
PHYS_scalar <- c(3, -2, NA, NA, NA, NA) # Tesetebed c(2, -2, NA, NA, NA, NA)
FI     <- c(0.01, 0.3) # Testbed values here  c(0.05, 0.05)
fmet_p <- c(0.6, 0.85, 0.013)  # reduced first number to lower values
depth <- 30 # set soil depth
CN_m     <- 15
CN_r     <- 6
CN_K     <- 10
cnModNum <- 0.4
densDep <- 1
h2y        <- 24*365
MICROtoECO <- depth *1e4 / 1e6 #convert kgC/m2 from mgC/cm3 (0-30 cm)
              #depth * 1e4 * 1e-3  # mgC/cm3 to g/m2

#Set default multipliers
Tau_MULT = 1
desorb_MULT = 1
fPHYS_MULT = 1

# N parameters
NUE <<- rep(0.9, 4) #0.85

Nleak <<- 0.2 #0.2