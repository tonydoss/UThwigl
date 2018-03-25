# enter sample name as in the data file
sample_name <- "PB_T1"

# enter coordinates of the sample outer and inner surfaces (in um).
# (1) refers to the surface closest to the first analysis,
# (2) to the surface closest to the last analysis.
x_s1 <- 17.762*1000
y_s1 <- 78.807*1000
x_s2 <- 18.897*1000
y_s2 <- 78.360*1000

# nb of iterations for the simulation (start with 1; see below)
nbit <- 100

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
T_max <- 8e3

# compute table with relative distances and activity ratios
source("R/iolite_to_iDAD.R")

# enter file name with data - create dataframe df with data from file
df <- read.csv(paste(sample_name,"_for_iDAD.csv",sep = "")) # create dataframe df with data from file

setwd("../")
source('R/iDAD Monte Carlo.R')
source("R/iDAD_direct.R")
source("R/graphs&save.R")
