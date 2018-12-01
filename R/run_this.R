rm(list=ls()) # delete everything

# ptm_start <- proc.time()

# path where data are
path <- "C:/Users/tonyd/OneDrive - University of Wollongong/LAB/MC-ICPMS laser ablation/20180321 Amy Tasmania"
# enter sample name as in the data file
sample_name <- "A"

# enter position of first (1) and last (2) surfaces (in um) 
x_s1 <- 38.06*1000
y_s1 <- 14.1*1000
x_s2 <- 38*1000
y_s2 <- 15.2*1000

# nb of iterations for the simulation (start with 1; see below)
nbit <- 1000

# enter value for squared difference; start with low value and 1 iteration (nbit <- 1)
# if script runs rapidly, try a lower value for fsum_target
# if script takes too long, try a higher value for fsum_target
# once happy with the value and pace of script, change number of iterations above (nbit)
fsum_target <- 5
# thickness of the sample (cm)
l <- 0.1
# uranium concentration at the surface (ppm)
U_0 <- 2
# min and max value for K (cm2/s)
K_min <- 1e-14
K_max <- 1e-11
# min and max value for the (234U/238U) activity ratio at the surface
U48_0_min <- 1.1
U48_0_max <- 2.0
# min and max value for the age (yr)
T_min <- 3e3
T_max <- 800e3

setwd(path)
# compute table with relative distances and activity ratios
source("C:/Users/tonyd/OneDrive - University of Wollongong/R/MC-ICPMS laser ablation/iolite_to_iDAD.R")

# enter file name with data - create dataframe df with data from file
df <- read.csv(paste(sample_name,"_for_iDAD.csv",sep = "")) # create dataframe df with data from file

setwd("C:/Users/tonyd/OneDrive - University of Wollongong/R/iDAD model") # set working directory
source('iDAD Monte Carlo.R')
source("iDAD_direct.R")
source("graphs&save.R")