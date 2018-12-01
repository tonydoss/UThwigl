l <- l/10/2

l238 <- 0.1551e-9/(365.25*24*3600)
l234 <- 2.826e-6/(365.25*24*3600)
l232 <- 4.948e-11/(365.25*24*3600)
l230 <- 9.158e-6/(365.25*24*3600)

length_series <- 10

x_vec <- df$iDAD.position

T_sol <- vector(mode="numeric", length=nbit)
U48_0_sol <- vector(mode="numeric", length=nbit)
K_sol <- vector(mode="numeric", length=nbit)

U48calc_sol <- matrix(, nrow = nbit, ncol = length(x_vec))
Th0U8calc_sol <- matrix(, nrow = nbit, ncol = length(x_vec))

series238 <- vector(mode="numeric", length=length_series+1)
beta234 <- vector(mode="numeric", length=length_series+1)
beta230 <- vector(mode="numeric", length=length_series+1)
series234 <- vector(mode="numeric", length=length_series+1)
gamma_model <- vector(mode="numeric", length=length_series+1)
series230 <- vector(mode="numeric", length=length_series+1)
series230_0 <- vector(mode="numeric", length=length_series+1)
beta234_0 <- vector(mode="numeric", length=length_series+1)
beta230_0 <- vector(mode="numeric", length=length_series+1)
gamma_model_0 <- vector(mode="numeric", length=length_series+1)

A1 <- vector(mode="numeric", length=1)
DA2 <- vector(mode="numeric", length=1)
A2 <- vector(mode="numeric", length=1)
A3 <- vector(mode="numeric", length=1)
f <- vector(mode="numeric", length=1)

sum_sq <- vector(mode="numeric", length=length(x_vec))
U48calc <- vector(mode="numeric", length=length(x_vec))
Th0U8calc <- vector(mode="numeric", length=length(x_vec))
Th0U4calc <- vector(mode="numeric", length=length(x_vec))
U48obs <- vector(mode="numeric", length=length(x_vec))
Th0U8obs <- vector(mode="numeric", length=length(x_vec))

R48_min <- df$U234_U238_CORR - df$U234_U238_CORR_Int2SE
R48_max <- df$U234_U238_CORR + df$U234_U238_CORR_Int2SE
R08_min <- df$Th230_U238_CORR - df$Th230_U238_CORR_Int2SE
R08_max <- df$Th230_U238_CORR + df$Th230_U238_CORR_Int2SE

counter <- 0

# repeat simulation 'nbit' number of times for a given sample
while (counter < nbit){
  K <- runif(1, K_min, K_max)
  U48_0 <- runif(1, U48_0_min, U48_0_max)
  T <- runif(1, T_min, T_max)
  for (ii in 1:length(x_vec)){
    U48obs[ii] <- runif(1, R48_min[ii], R48_max[ii])
    Th0U8obs[ii] <- runif(1, R08_min[ii], R08_max[ii])
  }
  
  t <- T*(365.25*24*3600)
  A1_0 <- U_0*1e-6/238*6.02e23*l238 # (238U) at the surface of the bone (disintegrations per second)
  A2_0 <- U48_0*A1_0
  DA2_0 <- A2_0 - A1_0 # (234U) excess at the surface of the bone (disintegrations per second)
  
  i <- 0
  
  for (x in x_vec){
    # for (x in -0.5){
    i <- i+1
    
    for (n in 0:length_series){
      series238[n+1] <- (-1)^n/(2*n + 1)*exp(-K*((2*n + 1)^2)*pi^2*t/(4*l^2))*cos((2*n + 1)/2*pi*x/l)
      beta234[n+1] <- 1 + 4*l234*(l^2)/((2*n + 1)^2)*pi^2*K
      beta230[n+1] <- l230 - l234 - K*(2*n + 1)^2*pi^2/(4*l^2)
      gamma_model[n+1] <- l230*(1/(beta230[n+1] + l234) + (U48_0 - 1)*exp(-l234*t)/(beta234[n+1]*beta230[n+1]))
      series234[n+1] <- (-1)^n/((2*n + 1)*beta234[n+1])*exp(-l234*t-K*(2*n + 1)^2*pi^2*t/(4*l^2))*cos((2*n + 1)/2*pi*x/l)
      series230[n+1] <- (-1)^n/(2*n + 1)*gamma_model[n+1]*exp(-K*(2*n + 1)^2*pi^2*t/(4*l^2))*cos((2*n + 1)/2*pi*x/l)
    }
    
    sum_series238 <- sum(series238)
    sum_series234 <- sum(series234)
    sum_series230 <- sum(series230)
    
    A1[i] <- A1_0*(1 - 4/pi*sum_series238)
    DA2[i] <- DA2_0*(cosh(x*(l234/K)^0.5)/(cosh(l*l234/K)^0.5) - 4/pi*sum_series234)
    
    t0 <- 0
    
    for (n in 0:length_series){
      beta234_0[n+1] <- 1 + 4*l234*(l^2)/((2*n + 1)^2)*pi^2*K
      beta230_0[n+1] <- l230 - l234 - K*(2*n + 1)^2*pi^2/(4*l^2)
      gamma_model_0[n+1] <- l230*(1/(beta230[n+1] + l234) + (U48_0 - 1)*exp(-l234*t0)/(beta234[n+1]*beta230[n+1]))
      series230_0[n+1] <- (-1)^n/(2*n + 1)*gamma_model[n+1]*exp(-K*(2*n + 1)^2*pi^2*t0/(4*l^2))*cos((2*n + 1)/2*pi*x/l)
    }
    sum_series230_0 <- sum(series230_0)
    f[i] <- -A1_0*(1 + (U48_0 - 1)*cosh(x*(l234/K)^0.5)/(cosh(l*(l234/K)^0.5) - 4/pi*sum_series230_0))*exp(-l230*t)
    A3[i] <- f[i] + A1_0*(1 + (U48_0 - 1)*cosh(x*(l234/K)^0.5)/(cosh(l*(l234/K)^0.5) - 4/pi*sum_series230))
    A2[i] <- DA2[i] + A1[i]
    U48calc[i] <- A2[i]/A1[i]
    Th0U8calc[i] <- A3[i]/A1[i]
    Th0U4calc[i] <- A3[i]/A2[i]
  }
  
  for (z in 1:length(x_vec)){
    # for (z in 1){
    sum_sq[z] <- (U48calc[z] - U48obs[z])^2 + (Th0U8calc[z] - Th0U8obs[z])^2
  }
  
  fsum <- sum(sum_sq) 
  
  if(fsum < fsum_target){
    counter <- counter+1
    
    U48calc_sol[counter,] <- U48calc
    Th0U8calc_sol[counter,] <- Th0U8calc
    T_sol[counter] <- T
    U48_0_sol[counter] <- U48_0
    K_sol[counter] <- K
  }
}

diff <- T_sol - median(T_sol)
T_final <- T_sol[diff == min(abs(diff))]
K_final <- K_sol[diff == min(abs(diff))]
U48_0_final <- U48_0_sol[diff == min(abs(diff))]





