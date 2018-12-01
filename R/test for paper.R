rm(list=ls()) # delete everything

# ptm_start <- proc.time()

path <- c("C:/Users/tonyd/OneDrive - University of Wollongong/R/iDAD model/output/")
sample_name <- c("Hobbit_MH2T")
# enter file name with data - create dataframe df with data from file
df <- read.csv("C:/Users/tonyd/OneDrive - University of Wollongong/R/iDAD model/input/Hobbit_MH2T_for_iDAD.csv") # create dataframe df with data from file

# nb of iterations for the simulation (start with 1; see below)
nbit <- 100

# enter value for fsum_target (squared difference); 
fsum_target <- 0.01
# start with low value (e.g. 0.05) and 1 iteration (nbit <- 1)
# if script runs rapidly, try a lower value for fsum_target
# if script takes too long, try a higher value for fsum_target
# once happy with the value and pace of script, adjust the min and max values
# for U480 (by looking at the calculated U48_0_final), K (by looking at K_final) and T (by looking at T_final)
# increase the number of iterations above (nbit) and run the script. Increase nbit to get a better error

# thickness of the sample (cm)
l <- 5.35 # Hobbit_1-1T: 3.5 cm; Hobbit_MH2T: 5.35 cm
# uranium concentration at the surface (ppm)
U_0 <- 25 # Hobbit_1-1T: 15 ppm; Hobbit_MH2T: 25 ppm
# min and max value for K (cm2/s)
K_min <- 1e-13
K_max <- 1e-11
# min and max value for the (234U/238U) activity ratio at the surface
U48_0_min <- 1.2 # Hobbit_1-1T: 1.3; Hobbit_MH2T: 1.25
U48_0_max <- 1.3 # Hobbit_1-1T: 1.4; Hobbit_MH2T: 1.3
# min and max value for the age (yr)
T_min <- 5e3 # Hobbit_1-1T: 50e3; Hobbit_MH2T: 1e3
T_max <- 10e3 # Hobbit_1-1T: 100e3; Hobbit_MH2T: 20e3

source('iDAD Monte Carlo.R')
source("iDAD_direct.R")
source("graphs&save.R")
setwd("C:/Users/tonyd/OneDrive - University of Wollongong/R/iDAD model/")
