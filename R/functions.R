
globalVariables(c("iDAD.position",
                  "Th230_U238_CORR",
                  "Th230_U238_CORR_Int2SE",
                  "Th230_U238_CALC",
                  "U_ppm_Int2SE",
                  "U234_U238_CORR",
                  "U234_U238_CORR_Int2SE",
                  "U234_U238_CALC",
                  "T_sol",
                  "U_ppm",
                  "U_ppm_Int2SE"
                  

))

#' Set a custom message the the user appears when they attach the pkg
#' 
#' @param libname Nothing to do
#' @param pkgname Nothing to do
.onAttach <- function(libname, pkgname) {
  
  # show the citation, but how to get a citation as a chr?: 
  
  packageStartupMessage("To cite the iDADwigl package use:

  Dosseto, A. and B. Marwick, (2019) iDADwigl: An R package
  for open-system uranium-thorium dating Quaternary
  Geochronology 0:000--000, http://doi.org/10.17605/OSF.IO/D5P7S

A BibTeX entry for LaTeX users can be obtained with
'bibentry('iDADwigl')'.

As iDADwigl is continually evolving, you may want to cite
its version number. Find it with 'help(package=iDADwigl)'.
                        ")
 
}


#' iDADwigl
#'
#' @param input_data A data frame with the required columns
#' @param nbit  The number of iterations
#' @param fsum_target The sum of the squared differences between the calculated and observed activity ratios.
#' @param U48_0_min The minimum value allowed for the (^234^U/^238^U) activity ratio at the surface of the sample
#' @param U48_0_max The maximum value allowed for the (^234^U/^238^U) activity ratio at the surface of the sample
#' @param l The thickness of the sample in centimeters
#' @param U_0 The uranium concentration at the surface in ppm
#' @param K_min The minimum value allowed for the uranium diffusion coefficient (in cm^2^/s)
#' @param K_max The maximum value allowed for the uranium diffusion coefficient (in cm^2^/s)
#' @param T_min The minimum value for the age of the specimen (yr)
#' @param T_max The maximum value for the age of the specimen (yr)
#' @param print_summary Print a summary of the output to the console? Default is TRUE
#' @param with_plots Display a panel of plots of the output? Default is TRUE
#' 
#' @importFrom cowplot plot_grid
#' @importFrom stats median quantile runif 
#' @importFrom utils ? 
#'
#' @return A list of results
#' @export
#'
#' @examples
#' 
#' data("Hobbit_MH2T_for_iDAD")
#' output <- iDADwigl(Hobbit_MH2T_for_iDAD,
#' nbit = 1,
#' fsum_target = 0.01,
#' U48_0_min = 1.265, # Hobbit_1-1T: 1.3; Hobbit_MH2T: 1.265
#' U48_0_max = 1.275, # Hobbit_1-1T: 1.4; Hobbit_MH2T: 1.275
#' l = 5.35, # Hobbit_1-1T: 3.5 cm; Hobbit_MH2T: 5.35 cm
#' U_0 = 25, # Hobbit_1-1T: 15 ppm; Hobbit_MH2T: 25 ppm
#' K_min = 1e-13,
#' K_max = 1e-11,
#' T_min = 1e3, # Hobbit_1-1T: 50e3; Hobbit_MH2T: 1e3
#' T_max = 20e3, # Hobbit_1-1T: 100e3; Hobbit_MH2T: 20e3
#' print_summary = TRUE,
#' with_plots = TRUE)
#' 

iDADwigl <- function(input_data,
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
                     print_summary = TRUE,
                     with_plots = TRUE
                     
                            
                     
                     
){
  
  
  # check that the input data frame has the columns with the right names
  
  col_names_we_need <-  c("iDAD.position", "U234_U238_CORR", "U234_U238_CORR_Int2SE", "Th230_U238_CORR", "Th230_U238_CORR_Int2SE")
  
  if(all(col_names_we_need %in% colnames(input_data)))
  {
    message("All required columns are present in the input data. \n");
  } else {
    ?iDADwigl::iDADwigl
    stop("\nThe input data frame does not contain the necessary columns, or the columns are not named correctly. Please check the documentation for details of the required column names, update the column names using the `names()` function, and try again.\n")
  }
  

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
  
  T_sol_df = as.data.frame(T_sol)  
  
  
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
  colnames(output_data)[  c( (ncol(output_data)-1), ncol(output_data) ) ] <- 
    c("U234_U238_CALC", "Th230_U238_CALC")
  
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
  output <- (list(results = results,
                  diff = diff,
                  T_final = T_final,
                  K_final = K_final,
                  T_sol = T_sol_df,
                  U48_0_final = U48_0_final,
                  output_data = output_data))
  
# plot or not? ------------------------------------
if(with_plots){
  # draw plots in a panel
  T_sol_plot_output <- T_sol_plot(output)
  u_conc_profile_plot_output <- u_conc_profile_plot(output)
  u234_u238_ratio_plot_output <- u234_u238_ratio_plot(output)
  th230_u238_ratio_plot_output <- th230_u238_ratio_plot(output)
  
  p1 <-
    cowplot::plot_grid(T_sol_plot_output,
              u_conc_profile_plot_output,
              u234_u238_ratio_plot_output,
              th230_u238_ratio_plot_output,
              labels = "AUTO",
              ncol = 2)
  
  print(p1)
  
}else {
  # don't plot anything
}
  
  
return(output)

}


#--------------------------------------------------------------------
# functions to draw plots with the output

#' Histogram of the solution ages
#' 
#' @param output Output from the `iDADwigl()` function
#' @param big_size Size of the main text on the plot, default is 10
#' @param less_big_size Size of the minor text on the plot, default is 8
#' @param point_size Size of the data points on the plot, default is 2
#' @param digits Number of digits to round the numeric output displayed as the plot title, default is 1
#' @import ggplot2
#' @export

T_sol_plot <- function(output,
                       big_size = 10,
                       less_big_size = 8,
                       point_size = 2,
                       digits = 1){
  
  theme_plots <-
    theme(
      plot.title = element_text(size = big_size, hjust = 0.5),
      legend.title = element_blank(),
      axis.title.y = element_text(size = big_size),
      axis.title.x = element_text(size = big_size),
      axis.text.x = element_text(size = less_big_size),
      axis.text.y = element_text(size = less_big_size)
    ) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  ggplot(output$T_sol,
         aes(T_sol)) +
  geom_histogram(binwidth = 500,
                 fill = "white",
                 color = "black") +
    ggtitle(paste(
    "Age: ",
    round(output$T_final / 1000, digits = digits),
    " +",
    round((quantile(output$T_sol$T_sol, .67) - output$T_final) / 1000, digits = digits),
    "/-",
    round((output$T_final - quantile(output$T_sol$T_sol, .33)) / 1000, digits = digits),
    " ka",
    sep = ""
  )) +  theme_plots 

}


#' Uranium concentration profile for transect
#' 
#' 
#' @param output Output from the `iDADwigl()` function
#' @param big_size Size of the main text on the plot, default is 10
#' @param less_big_size Size of the minor text on the plot, default is 8
#' @param point_size Size of the data points on the plot, default is 2
#' @param digits Number of digits to round the numeric output displayed as the plot title, default is 1
#' @import ggplot2
#' @export

u_conc_profile_plot <- function(output,
           big_size = 10,
           less_big_size = 8,
           point_size = 2,
           digits = 1){
    
    theme_plots <-
      theme(
        plot.title = element_text(size = big_size, hjust = 0.5),
        legend.title = element_blank(),
        axis.title.y = element_text(size = big_size),
        axis.title.x = element_text(size = big_size),
        axis.text.x = element_text(size = less_big_size),
        axis.text.y = element_text(size = less_big_size)
      ) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
  ggplot(output$output_data) +
  geom_errorbar(aes(x = iDAD.position,
                    ymax = U_ppm + U_ppm_Int2SE,
                    ymin = U_ppm - U_ppm_Int2SE, width = 0.02)) + # plot error bars
  geom_point(aes(iDAD.position, U_ppm),
             color = "blue",
             size = point_size) +
  ylab("U (ppm)") +
  xlab("Relative distance from center") +
  theme_plots
}

#' Calculated (red) and observed (blue) (^234^U/^238^U) activity ratios for transect
#' 
#' 
#' @param output Output from the `iDADwigl()` function
#' @param big_size Size of the main text on the plot, default is 10
#' @param less_big_size Size of the minor text on the plot, default is 8
#' @param point_size Size of the data points on the plot, default is 2
#' @param digits Number of digits to round the numeric output displayed as the plot title, default is 1
#' @import ggplot2
#' @export

u234_u238_ratio_plot <- function(output,
                                 big_size = 10,
                                 less_big_size = 8,
                                 point_size = 2,
                                 digits = 1){
  
  theme_plots <-
    theme(
      plot.title = element_text(size = big_size, hjust = 0.5),
      legend.title = element_blank(),
      axis.title.y = element_text(size = big_size),
      axis.title.x = element_text(size = big_size),
      axis.text.x = element_text(size = less_big_size),
      axis.text.y = element_text(size = less_big_size)
    ) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  ggplot(output$output_data) +
  geom_errorbar(aes(x = iDAD.position,
                    ymax = U234_U238_CORR + U234_U238_CORR_Int2SE,
                    ymin = U234_U238_CORR - U234_U238_CORR_Int2SE,
                    width = 0.02)) + # plot error bars
  geom_point(aes(iDAD.position,
                 U234_U238_CORR),
             color = "blue",
             size = point_size) +
  geom_point(aes(iDAD.position,
                 U234_U238_CALC),
             color = "red",
             size = point_size) +
  ylab(expression("("^234 * "U/"^238 * "U)")) +
  xlab("Relative distance from center") +
  theme_plots
}

#' Calculated (red) and observed (blue) (^230^Th/^238^U) activity ratios for transect 
#' 
#' @param output Output from the `iDADwigl()` function
#' @param big_size Size of the main text on the plot, default is 10
#' @param less_big_size Size of the minor text on the plot, default is 8
#' @param point_size Size of the data points on the plot, default is 2
#' @param digits Number of digits to round the numeric output displayed as the plot title, default is 1
#' 
#' @import ggplot2
#' @export

th230_u238_ratio_plot <-  function(output,
         big_size = 10,
         less_big_size = 8,
         point_size = 2,
         digits = 1){
  
  theme_plots <-
    theme(
      plot.title = element_text(size = big_size, hjust = 0.5),
      legend.title = element_blank(),
      axis.title.y = element_text(size = big_size),
      axis.title.x = element_text(size = big_size),
      axis.text.x = element_text(size = less_big_size),
      axis.text.y = element_text(size = less_big_size)
    ) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  ggplot(output$output_data) +
  geom_errorbar(
    aes(
      x = iDAD.position,
      ymax = Th230_U238_CORR + Th230_U238_CORR_Int2SE,
      ymin = Th230_U238_CORR - Th230_U238_CORR_Int2SE,
      width = 0.02
    )
  ) + # plot error bars
  geom_point(aes(iDAD.position, Th230_U238_CORR),
             color = 'blue',
             size = point_size) +
  geom_point(aes(iDAD.position, Th230_U238_CALC),
             color = 'red',
             size = point_size) +
  ylab(expression("(" ^ 230 * "Th/" ^ 238 * "U)")) + xlab("Relative distance from center") +
  theme_plots
}

# document code objects ----------------------------------------

#' Data from transect 1 for Homo floresiensis ulna LB1/52
#'
#' A dataset containing the U-Th measurement values for H. floresiensis.
#'
#' @format A data frame with 6 rows and 8 variables:
#' \describe{
#'   \item{iDAD.position}{...}
#'   \item{U234_U238_CORR}{...}
#'   \item{U234_U238_CORR_Int2SE}{...}
#'   \item{iDAD.position.1}{...}
#'   \item{Th230_U238_CORR}{...}
#'   \item{Th230_U238_CORR_Int2SE}{...}
#'   \item{U_ppm}{...}
#'   \item{U_ppm_Int2SE}{...}
#'   ...
#' }
#' @source \url{http://dx.doi.org/10.1038/nature17179}
"Hobbit_1_1T_for_iDAD"

#' Data transect 2 for modern human femur 132A/LB/27D/03
#'
#' A dataset containing the U-Th measurement values for a modern human
#'
#' @format A data frame with 6 rows and 8 variables:
#' \describe{
#'   \item{iDAD.position}{...}
#'   \item{U234_U238_CORR}{...}
#'   \item{U234_U238_CORR_Int2SE}{...}
#'   \item{iDAD.position.1}{...}
#'   \item{Th230_U238_CORR}{...}
#'   \item{Th230_U238_CORR_Int2SE}{...}
#'   \item{U_ppm}{...}
#'   \item{U_ppm_Int2SE}{...}
#'   ...
#' }
#' @source \url{http://dx.doi.org/10.1038/nature17179}
"Hobbit_MH2T_for_iDAD"


