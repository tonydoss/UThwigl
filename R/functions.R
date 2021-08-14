
globalVariables(c("iDAD.position",
                  "Th230_U238",
                  "Th230_U238_2SE",
                  "Th230_U238_CALC",
                  "U_ppm_2SE",
                  "U234_U238",
                  "U234_U238_2SE",
                  "U234_U238_CALC",
                  "T_sol",
                  "U_ppm",
                  "U_ppm_2SE",
                  "ID",
                  "Age (ka)", 
                  "Age 2sd", 
                  "[234U/238U]i", 
                  "Ratio 2sd",
                  "Sample ID"
))

#' Set a custom message the the user appears when they attach the pkg
#' 
#' @param libname Nothing to do
#' @param pkgname Nothing to do
.onAttach <- function(libname, pkgname) {
  
  # show the citation, but how to get a citation as a chr?: 
  
  packageStartupMessage("To cite the UThwigl package use:

  Dosseto, A. and B. Marwick, (2021) UThwigl: An R package
  for closed- and open-system uranium-thorium dating Quaternary
  Geochronology 0:000--000, http://doi.org/10.17605/OSF.IO/D5P7S

A BibTeX entry for LaTeX users can be obtained with
'bibentry('UThwigl')'.

As UThwigl is continually evolving, you may want to cite
its version number. Find it with 'help(package=UThwigl)'.
                        ")
 
}


#' osUTh
#'
#' @param input_data A data frame with the required columns
#' @param nbit  The number of iterations
#' @param fsum_target The sum of the squared differences between the calculated and observed activity ratios.
#' @param U48_0_min The minimum value allowed for the (^234^U/^238^U) activity ratio at the surface of the sample
#' @param U48_0_max The maximum value allowed for the (^234^U/^238^U) activity ratio at the surface of the sample
#' @param U_0 The uranium concentration at the surface in ppm
#' @param K_min The minimum value allowed for the uranium diffusion coefficient (in cm^2^/s)
#' @param K_max The maximum value allowed for the uranium diffusion coefficient (in cm^2^/s)
#' @param T_min The minimum value for the age of the specimen (yr)
#' @param T_max The maximum value for the age of the specimen (yr)
#' @param print_age Print a summary of the output to the console? Default is TRUE
#' @param with_plots Display a panel of plots of the output? Default is TRUE
#' @param save_plots Save plots as a png file to the current working directory? Default: TRUE
#' @param save_output Save output data as a CSV file to the current working directory? Default: FALSE
#' 
#' 
#' @importFrom cowplot plot_grid
#' @importFrom stats median quantile runif 
#' @importFrom utils ? write.csv
#' @importFrom grDevices dev.off png
#'
#' @return A list of results
#' @export
#'
#' @examples
#' data("Hobbit_MH2T_for_iDAD")
#' output <- osUTh(Hobbit_MH2T_for_iDAD,
#' nbit = 1,
#' fsum_target = 0.01,
#' U48_0_min = 1.265,
#' U48_0_max = 1.275,
#' U_0 = 25,
#' K_min = 1e-13,
#' K_max = 1e-11,
#' T_min = 1e3,
#' T_max = 20e3,
#' print_age = TRUE,
#' with_plots = FALSE,
#' save_plots = FALSE)
#' 

osUTh <- function(input_data,
                     nbit = 1000,
                     fsum_target = 0.01,
                     U48_0_min = 1.265, # Hobbit_1-1T: 1.3; Hobbit_MH2T: 1.265
                     U48_0_max = 1.275, # Hobbit_1-1T: 1.4; Hobbit_MH2T: 1.275
                     U_0 = 25, # Hobbit_1-1T: 15 ppm; Hobbit_MH2T: 25 ppm
                     K_min = 1e-13,
                     K_max = 1e-11,
                     T_min = 1e3, # Hobbit_1-1T: 50e3; Hobbit_MH2T: 1e3
                     T_max = 20e3, # Hobbit_1-1T: 100e3; Hobbit_MH2T: 20e3
                     print_age = TRUE,
                     with_plots = TRUE,
                     save_plots = FALSE,
                     save_output = FALSE
                  ){
  
  
  # check that the input data frame has the columns with the right names
  
  col_names_we_need <-  c("x", "y", "Comments", "U234_U238", "U234_U238_2SE", "Th230_U238", "Th230_U238_2SE")
  
  if (all(col_names_we_need %in% colnames(input_data)))
  {
    message("All required columns are present in the input data. \n");
  } else {
    ?UThwigl::osUTh
    stop(paste0("\nThe input data frame does not contain the necessary columns, or the columns are not named correctly. 
         Please check the documentation for details of the required column names, 
         update the column names using the `names()` function, and try again.\n
         Your data has these column names: ", colnames(input_data), "\n
         This function requires data with these column names: ", col_names_we_need,
         "\n"))
  }
  


# calculate iDAD positions and sample thickness ---------------------------

  # Coordinates of outer and inner surfaces
  x_s1 <- input_data$x[grepl("outer surface", input_data$Comments)]
  y_s1 <- input_data$y[grepl("outer surface", input_data$Comments)]
  x_s2 <- input_data$x[grepl("inner surface", input_data$Comments)]
  y_s2 <- input_data$y[grepl("inner surface", input_data$Comments)]
  
  if (abs(x_s2 - x_s1) > abs(y_s2 - y_s1)) {
    input_data$iDAD.position <- (2/(x_s2 - x_s1)*(input_data$x - x_s2) + 1)
  } else {
    input_data$iDAD.position <- (2/(y_s2 - y_s1)*(input_data$y - y_s2) + 1)
  }
  
  # Calculate sample thickness (in cm)
  l <- sqrt((x_s2 - x_s1)^2 + (y_s2 - y_s1)^2)/10
  
  # Remove rows that don't have data (rows with coordinates of outer and inner surfaces)
  input_data <- input_data[!is.na(input_data$U234_U238),]
  # Remove column "Comments" (we don't need it anymore)
  input_data <-  input_data[ , ! names(input_data) %in% c("Comments")]
  
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
  
  R48_min <- input_data$U234_U238 -  input_data$U234_U238_2SE
  R48_max <- input_data$U234_U238 +  input_data$U234_U238_2SE
  R08_min <- input_data$Th230_U238 - input_data$Th230_U238_2SE
  R08_max <- input_data$Th230_U238 + input_data$Th230_U238_2SE
  
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
        (quantile(T_sol, .841) - T_final) / 1000,
        (T_final - quantile(T_sol, .159)) / 1000,
        U48_0_final,
        quantile(U48_0_sol, .841) - U48_0_final,
        U48_0_final - quantile(U48_0_sol, .159)
      )
    )
  
  colnames(results) <-
    c(
      "Age (ka)",
      "Age +1SD (ka)" ,
      "Age -1SD (ka)",
      "U234_U238_0",
      "U234_U238_0 +1SD" ,
      "U234_U238_0 -1SD"
    )
  rownames(results) <- c("Results")
  
  
if(print_age) { 
  print(paste(
    "Age: ",
    round(T_final / 1000, digits = 1),
    " +",
    round((quantile(T_sol, .841) - T_final) / 1000, digits = 1),
    "/-",
    round((T_final - quantile(T_sol, .159)) / 1000, digits = 1),
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
                  output_data = output_data,
                  plots = NULL))
  
  
  # for unique file names
  st <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  
# plot or not? ------------------------------------
if(with_plots){
  message("Drawing plots...")
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
  
  output$plots <-  p1
  message("Done.")
  
  if(save_plots){
    
    plot_file_name <- paste0("osUTh-plots-", st, ".png")
    message(paste0("Saving plots to ", getwd(), ", look for '", plot_file_name, "'" ))
    png(plot_file_name, width = 15, height = 15, units = "cm", res = 300 )
    print(output$plots)
    dev.off()
    
  } else {
    # don't save plots
  }
  
}else {
  # don't plot anything
}
  
  if(save_output){
    
    filename <-  paste0("osUTh-output-", st, ".csv")
    write.csv(output$results, filename)
    message(paste0("Saving output to ", getwd(), ", look for '", filename, "'"))
    
    
  } else {
    # don't save anything
  }
  
  
return(output)

}


#--------------------------------------------------------------------
# functions to draw plots with the output

#' Histogram of the solution ages
#' 
#' @param output Output from the `osUTh()` function
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
    round((quantile(output$T_sol$T_sol, .841) - output$T_final) / 1000, digits = digits),
    "/-",
    round((output$T_final - quantile(output$T_sol$T_sol, .159)) / 1000, digits = digits),
    " ka",
    sep = ""
  )) +  theme_plots 

}


#' Uranium concentration profile for transect
#' 
#' 
#' @param output Output from the `osUTh()` function
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
                    ymax = U_ppm + U_ppm_2SE,
                    ymin = U_ppm - U_ppm_2SE, width = 0.02)) + # plot error bars
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
#' @param output Output from the `osUTh()` function
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
                    ymax = U234_U238 + U234_U238_2SE,
                    ymin = U234_U238 - U234_U238_2SE,
                    width = 0.02)) + # plot error bars
  geom_point(aes(iDAD.position,
                 U234_U238),
             color = "blue",
             size = point_size) +
  geom_point(aes(iDAD.position,
                 U234_U238_CALC),
             color = "red",
             size = point_size) +
  ylab(expression("["^234 * "U/"^238 * "U]")) +
  xlab("Relative distance from center") +
  theme_plots
}

#' Calculated (red) and observed (blue) (^230^Th/^238^U) activity ratios for transect 
#' 
#' @param output Output from the `osUTh()` function
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
      ymax = Th230_U238 + Th230_U238_2SE,
      ymin = Th230_U238 - Th230_U238_2SE,
      width = 0.02
    )
  ) + # plot error bars
  geom_point(aes(iDAD.position, Th230_U238),
             color = 'blue',
             size = point_size) +
  geom_point(aes(iDAD.position, Th230_U238_CALC),
             color = 'red',
             size = point_size) +
  ylab(expression("[" ^ 230 * "Th/" ^ 238 * "U]")) + xlab("Relative distance from center") +
  theme_plots
}


#--------------------------------------------------------------------
# function for closed-system dating

#' csUTh
#'
#' csUTh calculates closed-system Th-230/U ages, including detrital correction.
#'
#'
#' @param input_data Input data frame. The following columns need to be present in this data frame, with these exact names: Sample_ID, U234_U238, U234_U238_2SE, Th230_U238, Th230_U238_2SE, Th232_U238, Th232_U238_2SE.
#' @param sample_name Name of the sample to calculate closed-system ages for. The string entered must match characters for the chosen sample in the column 'Sample_ID' of the data file. Default: 'MK16'.
#' @param nbitchoice Number of iterations in the model. Recommended to have at least 100. Default: 100.
#' @param detcorrectionchoice Do a detrital correction? Enter TRUE for yes, or FALSE for no. Default: TRUE
#' @param R28det (232Th/238U) activity ratio of the detritus. Default: 0.8
#' @param R28det_err Error on the (232Th/238U) activity ratio of the detritus. Default: 0.4
#' @param R08det (230Th/238U) activity ratio of the detritus. Default: 1
#' @param R08det_err Error on the (230Th/238U) activity ratio of the detritus. Default: 0.05
#' @param R48det (234U/238U) activity ratio of the detritus. Default: 1
#' @param R48det_err Error on the (234U/238U) activity ratio of the detritus. Default: 0.02
#' @param keepfiltereddata Save filtered data on which an outlier test was performed? Only recommended if all analyses of a same sample are supposed to give the same age. Enter TRUE for yes, or FALSE for no. Default: FALSE
#' @param print_age Print a summary of the output to the console? Default: TRUE
#' @param with_plots Draw plots? Default: TRUE
#' @param save_plots Save plots as a png file to the current working directory? Default: TRUE
#' @param save_output Save output data as a CSV file to the current working directory? Default: FALSE
#' 
#' @import deSolve ggplot2
#' @importFrom stats IQR optim sd
#' @importFrom grDevices dev.off png
#' @importFrom utils write.csv
#' 
#' @examples 
#' data("Pan2018")
#' # Solve for sample YP003
#' output <- csUTh(Pan2018,
#' sample_name = 'YP003',
#' nbitchoice = 100,
#' detcorrectionchoice = TRUE,
#' keepfiltereddata = FALSE,
#' print_age = TRUE,
#' with_plots = FALSE,
#' save_plots = FALSE)
#' 
#' @export

csUTh <- function(input_data, 
                  sample_name = 'YP003',
                  nbitchoice = 1000,
                  detcorrectionchoice = TRUE,
                  R28det = 0.8,
                  R28det_err = 0.4,
                  R08det = 1,
                  R08det_err = 0.05,
                  R48det = 1,
                  R48det_err = 0.02,
                  keepfiltereddata = FALSE,
                  print_age = TRUE,
                  with_plots = TRUE,
                  save_plots = TRUE,
                  save_output = FALSE
){
  
  # check that the input data frame has the columns with the right names
  
  
  col_names_we_absolutely_need <-  c("Sample_ID", "U234_U238", "U234_U238_2SE", "Th230_U238", "Th230_U238_2SE")
  col_names_we_need_either <- c("Th232_U238", "Th232_U238_2SE")
  col_names_we_need_or <- c("Th230_Th232", "Th230_Th232_2SE")
  
  if(all(col_names_we_absolutely_need %in% colnames(input_data)) & 
     (all(col_names_we_need_either %in% colnames(input_data)) | 
      all(col_names_we_need_or %in% colnames(input_data))))
  {
    message("All required columns are present in the input data. \n");
  } else {
    ?UThwigl::csUTh
    stop("\nThe input data frame does not contain the necessary columns, or the columns are not named correctly. Please check the documentation for details of the required column names, update the column names using the `names()` function, and try again.\n")
  }
  
  # calculate (232Th/238U) activity ratios and their error, if needed
  if("Th230_Th232" %in% colnames(input_data)){
    input_data$Th232_U238 <- input_data$Th230_U238/input_data$Th230_Th232
    input_data$Th232_U238_2SE <- input_data$Th232_U238*
      sqrt((input_data$Th230_U238_2SE/input_data$Th230_U238)^2)
  } 
  
  l234 <- 2.8262e-6 # 234U decay constant (a-1)
  l230 <- 9.1577e-6 # 230Th decay constant (a-1)
  
  # nb of times optimisation is repeated (for each sample)
  nbit <- nbitchoice
  
  lowerbound <- c(2, 01.0) # lower bound values for age (log10(yr)) and initial (234U/238U)
  upperbound <- c(6, 10.0) # upper bound values for age (log10(yr)) and initial (234U/238U)
  
  # use detrital correction (TRUE) or not (FALSE)
  detcorrection <- detcorrectionchoice
  
  # create dataframe with data only for samples to solve
  data <- subset(input_data, (grepl(sample_name, input_data$Sample_ID)))
  # number of samples to solve
  number_sampletosolve <- nrow(data)
  
  # parameters for detrial correction
  data$B <- (R08det - data$Th230_U238)/(R28det - data$Th232_U238)
  data$b <- (R48det - data$U234_U238)/(R28det - data$Th232_U238)
  data$r1 <- data$Th232_U238/(R28det - data$Th232_U238)
  data$r2 <- R28det/(R28det - data$Th232_U238)
  
  # detrital-corrected ratios
  data$U234_U238_DET <- data$U234_U238 - data$b*data$Th232_U238
  data$U234_U238_DET_ERR <- sqrt(data$b^2*(data$r2^2*data$Th232_U238_2SE^2 +
                                                  data$r1^2*R28det_err^2) +
                                        data$r2^2*data$U234_U238_2SE +
                                        data$r1^2*R48det_err^2)
  data$Th230_U238_DET <- data$Th230_U238 - data$B*data$Th232_U238
  data$Th230_U238_DET_ERR <- sqrt(data$B^2*(data$r2^2*data$Th232_U238_2SE^2 +
                                                   data$r1^2*R28det_err^2) +
                                         data$r2^2*data$Th230_U238_2SE +
                                         data$r1^2*R08det_err^2)
  
  # create vectors
  time_results <- vector(mode="numeric", length=number_sampletosolve)
  err_time_results <- vector(mode="numeric", length=number_sampletosolve)
  R48i_results <- vector(mode="numeric", length=number_sampletosolve)
  err_R48i_results <- vector(mode="numeric", length=number_sampletosolve)
  time2sd_results <- vector(mode="numeric", length=number_sampletosolve)
  R48i2sd_results <- vector(mode="numeric", length=number_sampletosolve)
  
  # repeat loop for each sample
  for (count in 1:number_sampletosolve){
    if (detcorrection == 'Y'){
      U48meas <- data$U234_U238_DET[count]
      Th0U8meas <- data$Th230_U238_DET[count]
      err_R08 <- data$Th230_U238_DET_ERR[count]
      err_R48 <- data$U234_U238_DET_ERR[count]
    } else {
      U48meas <- data$U234_U238[count]
      Th0U8meas <- data$Th230_U238[count]
      err_R08 <- data$Th230_U238_2SE[count]
      err_R48 <- data$U234_U238_2SE[count]
    }
    # create list for optimisation results
    sol=list()
    # create vectors
    U48calc <- vector(mode="numeric", length=nbit)
    Th0U8calc <- vector(mode="numeric", length=nbit)
    time <- vector(mode="numeric", length=nbit)
    R48i <- vector(mode="numeric", length=nbit)
    time_2sd <- vector(mode="numeric", length=nbit)
    R48i_2sd <- vector(mode="numeric", length=nbit)
    
    # repeat optimisation 'nbit' number of times for a given sample
    for (i in 1:nbit){
      # pick (234U/238U) and (230Th/238U) within range of measured values
      U48target <- runif(1, U48meas - err_R48, U48meas + err_R48)
      Th0U8target <- runif(1, Th0U8meas - err_R08, Th0U8meas + err_R08)
      
      # start optimisation with random age and initial (234U/23U) taken from the range of values allowed
      init_time <- runif(1, lowerbound[1], upperbound[1])
      init_R48i <- runif(1, lowerbound[2], upperbound[2])
      paraminit <- c(init_time, init_R48i)
      
      # function to minimise
      funmin <- function(x) {
        t <- 10^x[1] # age in yr
        U48i <- x[2] # intial (234U/238U)
        
        U48calc <- 1 + (U48i - 1)*exp(-l234*t) # (234U/238U)
        Th0U8calc <- 1 - exp(-l230*t) + (U48calc - 1)*(l230/(l230 - l234))*(1 - exp((l234 - l230)*t))
        
        fmin <- sum((U48calc - U48target)^2 + (Th0U8calc - Th0U8target)^2 ) # function to minimise
      }
      
      # optimisation
      sol <- optim(paraminit, funmin, method = "L-BFGS-B",
                   lower = lowerbound, upper = upperbound, control = list(factr = 1e-8))
      # store calculated age, initial (234U/23U) and calculated activity ratios for each optimisation
      time[i] <- 10^sol$par[1]
      R48i[i] <- sol$par[2]
      U48calc[i] <- 1 + (R48i[i] - 1)*exp(-l234*time[i]) # (234U/238U)
      Th0U8calc[i] <- 1 - exp(-l230*time[i]) - (U48calc[i] - 1)*(l230/(l234 - l230))*(1 - exp((l234 - l230)*time[i]))
    }
    
    # store results from all optimisations for a given sample
    results <- as.data.frame(cbind(time, R48i, U48calc, Th0U8calc))
    # take the median of all ages and initial (234U/23U)
    time_median <- median(results$time)
    time_2sd <- 2*sd(results$time)
    R48i_median <- median(results$R48i)
    R48i_2sd <- 2*sd(results$R48i)

    # store age, error on age and initial (234U/23U) for each sample
    time_results[count] <- time_median
    time2sd_results[count] <- time_2sd
    R48i_results[count] <- R48i_median
    R48i2sd_results[count] <- R48i_2sd
  }
  
  final_results <- as.data.frame(cbind(
                                       data,
                                       round(time_results/1000,3), round(time2sd_results/1000,3),
                                       round(R48i_results,3), round(R48i2sd_results,3)))
  final_results$Sample_ID <- data$Sample_ID
  colnames(final_results)[1] <- "Sample ID"
  colnames(final_results)[(ncol(final_results)-3):ncol(final_results)] <- c("Age (ka)", "Age 2sd", "[234U/238U]i", "Ratio 2sd")
  
  remove_outliers <- function(x, na.rm = TRUE, ...) {
    qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
    H <- 1.5 * IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- NA
    y[x > (qnt[2] + H)] <- NA
    y
  }
  
  final_results_filtered <- final_results
  final_results_filtered$`Age (ka)` <- remove_outliers(final_results$`Age (ka)`)
  final_results_filtered <- final_results_filtered[!is.na(final_results_filtered$`Age (ka)`),]
  
  if (!keepfiltereddata){
    plotdata <- final_results
  } else if (keepfiltereddata){
    plotdata <- final_results_filtered
  }
  
  # simplify output
  output <- list(results = NULL, plots = NULL)
  if("Th230_Th232" %in% colnames(input_data)){
    output$results <- plotdata[,c(1, 18:21)]}else{
      output$results <- plotdata[,c(1, 16:19)]
    }
  
  # plot initial (234U/238U)
  
  # draw plots
  st <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  
  # plot or not? ------------------------------------
  if(with_plots){
    message("Drawing plots...")
    
    p2 <- initial_234U_238U_plot(output) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    
    # plot ages
    p1 <- ages_plot(output) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    
    # draw plots in a panel
  p3 <- cowplot::plot_grid(p1, p2, ncol = 2)
  
  output$plots <-  p3
  
  # for unique file names
  
  
  if(save_plots){
    
    plot_file_name <- paste0("csUTh-plots-", st, ".png")
    message(paste0("Saving plots to ", getwd(), ", look for '", plot_file_name, "'"))
    png(plot_file_name, width = 15, height = 10, units = "cm", res = 300 )
    print(output$plots)
    dev.off()
    
  } else {
    # don't save plots
  }
  
  message("Done.")
  
  } else {
    # don't draw plots
  }
  
  output$print_age <- paste('Mean age: ',round(mean(output$results$`Age (ka)`, na.rm = TRUE),1),
                         '+/-', round(2*sd(output$results$`Age (ka)`, na.rm = TRUE)/
                                        sqrt(length(output$results$`Age (ka)`)), 1), ' ka')

  
  if(print_age){
    
  print(output$print_age)
    
  } else {
    # don't print anything
  }
  
  if(save_output){
    
    filename <-  paste0("csUTh-output-", st, ".csv")
    write.csv(output$results, filename)
    message(paste0("Saving output to ", getwd(), ", look for '", filename, "'"))
              
    
  } else {
    # don't save anything
  }
  
  # return results
  return(output)
}



#' Initial (^234^U/^238^U) plot
#' 
#' @param output Output from the `csUTh()` function
#' @param big_size Size of the main text on the plot, default is 10
#' @param less_big_size Size of the minor text on the plot, default is 8
#' @param point_size Size of the data points on the plot, default is 5
#' @import ggplot2
#' @export

initial_234U_238U_plot <- function(output,
                                 big_size = 10,
                                 less_big_size = 8,
                                 point_size = 5){
  
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
  
  ggplot(output$results, aes(`Sample ID`, `[234U/238U]i`)) + # plot ages
  geom_errorbar(aes(ymin = (`[234U/238U]i` - `Ratio 2sd`),
                    ymax = (`[234U/238U]i` + `Ratio 2sd`)), 
                width=0.1) + # plot error bars
  geom_point(size=point_size) + # plot points
  xlab("Sample ID") + # x axis label
  ylab(expression("Initial ["^234*"U/"^238*"U]")) + # y axis label
    theme_plots
}


#' Ages plot for closed-system analysis
#' 
#' 
#' @param output Output from the `csUTh()` function
#' @param big_size Size of the main text on the plot, default is 10
#' @param less_big_size Size of the minor text on the plot, default is 8
#' @param point_size Size of the data points on the plot, default is 5
#' @import ggplot2
#' @export
ages_plot <- function(output,
                   big_size = 10,
                   less_big_size = 8,
                   point_size = 5){
  
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
  
  ggplot(output$results, aes(`Sample ID`, `Age (ka)`)) + # plot ages
  geom_errorbar(aes(ymin = (`Age (ka)` - `Age 2sd`),
                    ymax = (`Age (ka)` + `Age 2sd`)), 
                width=0.1) + # plot error bars
  geom_point(size=point_size) + # plot points
  xlab("Sample ID") + # x axis label
  ylab("Age (ka)")  + # y axis label
  theme_plots
}


