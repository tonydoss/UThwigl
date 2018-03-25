
i <- 0

K <- K_final
t <- T_final*(365.25*24*3600)
A2_0 <- U48_0_final*A1_0
DA2_0 <- A2_0 - A1_0 # (234U) excess at the surface of the bone (disintegrations per second)

U48calc_final <- vector(mode="numeric", length=length(x_vec))
Th0U8calc_final <- vector(mode="numeric", length=length(x_vec))
Th0U4calc_final <- vector(mode="numeric", length=length(x_vec))

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
  U48calc_final[i] <- A2[i]/A1[i]
  Th0U8calc_final[i] <- A3[i]/A1[i]
  Th0U4calc_final[i] <- A3[i]/A2[i]
}


df <- cbind(df, U48calc_final, Th0U8calc_final)