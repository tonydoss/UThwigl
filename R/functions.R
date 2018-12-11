



#' Title
#'
#' @param input_data 
#' @param nbit 
#' @param fsum_target 
#' @param U48_0_min 
#' @param U48_0_max 
#' @param l 
#' @param U_0 
#' @param K_min 
#' @param K_max 
#' @param T_min 
#' @param T_max 
#' @param print_summary 
#'
#' @return
#' @export
#'
#' @examples
sambridge_et_al <- function(input_data,
                            nbit = 1,
                            fsum_target = 0.01,
                            U48_0_min = 1.265, # Hobbit_1-1T: 1.3; Hobbit_MH2T: 1.265
                            U48_0_max = 1.275, # Hobbit_1-1T: 1.4; Hobbit_MH2T: 1.275
                            l = 5.35, # Hobbit_1-1T: 3.5 cm; Hobbit_MH2T: 5.35 cm
                            U_0 = 25, # Hobbit_1-1T: 15 ppm; Hobbit_MH2T: 25 ppm
                            K_min = 1e-13,
                            K_max = 1e-11,
                            T_min = 1e3, # Hobbit_1-1T: 50e3; Hobbit_MH2T: 1e3
                            T_max = 20e3, # Hobbit_1-1T: 100e3; Hobbit_MH2T: 20e3
                            print_summary = TRUE
){
  

  
  # check that the input data frame has the columns with the right names
  # ...
  # names(input_data)
  
  
  

# from iDAD Monte Carlo.R -------------------------------------------------


  l <- l/10/2
  
  l238 <- 0.1551e-9/(365.25*24*3600)
  l234 <- 2.826e-6/(365.25*24*3600)
  l232 <- 4.948e-11/(365.25*24*3600)
  l230 <- 9.158e-6/(365.25*24*3600)
  
  length_series <- 10
  
  x_vec <- input_data$iDAD.position
  
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
  
  R48_min <- input_data$U234_U238_CORR -  input_data$U234_U238_CORR_Int2SE
  R48_max <- input_data$U234_U238_CORR +  input_data$U234_U238_CORR_Int2SE
  R08_min <- input_data$Th230_U238_CORR - input_data$Th230_U238_CORR_Int2SE
  R08_max <- input_data$Th230_U238_CORR + input_data$Th230_U238_CORR_Int2SE
  
  counter <- 0
  
  # repeat simulation 'nbit' number of times for a given sample
  # this takes a long time...
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
  } # end of the while loop
  
  diff <- T_sol - median(T_sol)
  T_final <- T_sol[diff == min(abs(diff))]
  K_final <- K_sol[diff == min(abs(diff))]
  U48_0_final <- U48_0_sol[diff == min(abs(diff))]
  
  # end of that script
  
# from iDAD_direct.R ------------------------------------------------------
  
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
  
  
  output_data <- cbind(input_data, U48calc_final, Th0U8calc_final)
  
# from graphs&save.R -------------------------------------------  
  results <-
    as.data.frame(
      cbind(
        T_final / 1000,
        (quantile(T_sol, .67) - T_final) / 1000,
        (T_final - quantile(T_sol, .33)) / 1000,
        U48_0_final,
        U48_0_max - U48_0_final,
        U48_0_final - U48_0_min
      )
    )
  
  colnames(results) <-
    c(
      "Age (ka)",
      "Age 67% quantile (ka)" ,
      "Age 33% quantile (ka)",
      "U234_U238_0",
      "U234_U238_0 67% quantile" ,
      "U234_U238_0 33% quantile"
    )
  rownames(results) <- c("Results")
  
  calc_ratios <- as.data.frame(cbind(U48calc_final, Th0U8calc_final))
  colnames(calc_ratios) <- c("U234_U238_CALC", "Th230_U238_CALC")
  
  
if(print_summary) { 
  print(paste(
    "Age: ",
    round(T_final / 1000, digits = 1),
    " +",
    round((quantile(T_sol, .67) - T_final) / 1000, digits = 1),
    "/-",
    round((T_final - quantile(T_sol, .33)) / 1000, digits = 1),
    " ka",
    sep = ""
  ))
} else {
  # don't print anything
}
  
  
  
# collect the output into a list ------------------------------------
  return(list(results = results,
              diff = diff,
              T_final = T_final,
              K_final = K_final,
              U48_0_final = U48_0_final,
              output_data = output_data,
              calc_ratios = calc_ratios))

}


# testing...
# sambridge_et_al(Hobbit_MH2T_for_iDAD, nbit = 10)
