# iDADwigl

iDADwigl is package for open-system uranium-thorium dating of geological and archaeological samples.
It is based on the diffusion-adsorption-decay (DAD) model of Sambridge et al. (2012), which allows for advective and diffusive transport of uranium and thorium isotopes, while including synchronous radioactive decay.

The key function of the package requires a data frame with the following column names: 

- `iDAD.position`
- `U234_U238_CORR`
- `U234_U238_CORR_Int2SE`
- `iDAD.position.1`
- `Th230_U238_CORR`
- `Th230_U238_CORR_Int2SE`
- `U_ppm`
- `U_ppm_Int2SE`        

This data frame can be created by importing an Excel or CSV file into the R environment using a generic function such as `read.csv` or `read_excel` from the `readxl` package. To help with preparing data for input into our function, two examples of input files are included. The code chunk below shows how to access one of the two example CSV files included in the package, and how to read it into the R environment.


```{r}
# get the path the one of the CSV files included in the package
path_to_included_example_csv_file <-
  system.file("extdata",
            "input/Hobbit_MH2T_for_iDAD.csv",
            package = "iDADwigl",
            mustWork = TRUE)

# read in the example CSV file included in the package
Hobbit_MH2T_for_iDAD <-
  read.csv(path_to_included_example_csv_file)
```

The columns `iDAD.position`, `U234_U238_CORR`, `U234_U238_CORR_Int2SE`, `Th230_U238_CORR` and `Th230_U238_CORR_Int2SE` must be present in the input data frame with these exact names for the model to function. The `iDADwigl()` function will check if the input data frame has these columns, and will stop with an error message if it does not find these columns. The `names` function can be used to update column names of a data frame to ensure they match the names that the model function requires. The order of the columns in the data frame is not important. 

The `iDAD.position` column corresponds to the coordinates of the (<sup>234</sup>U/<sup>238</sup>U) analyses, which as indicated above take values between -1 and 1. The second `iDAD.position.1` column is used if the coordinates of the (<sup>230</sup>Th/<sup>238</sup>U) analyses are different from those of the (<sup>234</sup>U/<sup>238</sup>U) analyses. 

Columns `U234_U238_CORR` and `U234_U238_CORR_Int2SE` are the (<sup>234</sup>U/<sup>238</sup>U) activity ratios and their 2\sigma errors. Columns `Th230_U238_CORR` and `Th230_U238_CORR_Int2SE` are the (<sup>230</sup>Th/<sup>238</sup>U) activity ratios and their 2\sigma errors. 

Columns `U_ppm` and `U_ppm_Int2SE` and calculated uranium concentrations (in ppm) and their 2\sigma errors. Uranium concentrations are not necessary for the model and only used for display of the U concentration profile in a figure.

# Details of the input parameters

Our key function, `iDADwigl()` has several arguments that need to be set before we can get meaningful results.

`nbit` is the number of iterations. For the first run, set to 1.

`fsum_target` is the sum of the squared differences between the calculated and observed activity ratios. Give it a low value to start with (e.g. 0.05). If script takes too long, try a higher value for fsum_target.

`U48_0_min` and `U48_0_max` are the minimum and maximum values allowed for the (<sup>234</sup>U/<sup>238</sup>U) activity ratio at the surface of the sample. Since (<sup>234</sup>U/<sup>238</sup>U) does not vary greatly over the time period generally studied, the values measured near the surface of the sample can be used as a guide. These values can be adjusted if the model fit to the data is not optimal. For Hobbit_1-1T they are taken to be 1.3 and 1.4, and for Hobbit_MH2T, 1.265 and 1.275, respectively.

`l` is the thickness of the sample in centimeters. For Hobbit_1-1T it is 3.5 cm, for Hobbit_MH2T it is 5.35 cm

`U_0` is the uranium concentration at the surface in ppm. This value does not significantly affect the model results and values from analyses near either surface of the sample can be used as a guide. For Hobbit_1-1T it is taken to be 15 ppm; for Hobbit_MH2T, 25 ppm.

`K_min` and `K_max` are the minimum and maximum values allowed for the uranium diffusion coefficient (in cm<sup>2</sup>/s). Values between 10<sup>-13</sup> and 1<sup>-11</sup> cm<sup>2</sup>/s are generally appropriate.

`T_min` and `T_max` are the minimum and maximum values for the age of the specimen (yr). If there is no estimated knowledge of the sample age, the range of values can be 1,000 to 500,000 yr and adjusted later. For Hobbit_1-1T, in the final model run, they are taken to be 50,000 and 100,000 yr, and for Hobbit_MH2T, 1,000 and 20,000 yr, respectively.

After setting the `U480` maximum and minimum values, run the function and adjust these min and max values by looking at the calculated `U48_0_final`, `K_final`, and `T_final`. Adjust `T_min` and `T_max` using first estimates of the age. As you iterate, increase the `nbit` value to reduce the error.

# How to run the model

To run the model, start `R` and install the package from GitHub. There are many ways to do this, one simple method is shown in the line below. This only needs to be done once per computer.

`source("https://install-github.me/tonydoss/iDADwigl")`

To run the function, attach the package and then run `iDADwigl()`, specifying the input data frame and the input parameters as described above. The code block below shows a quick example that will execute in less than five seconds on a typical 2.3 GHz Intel Core i5 laptop:  

```
library(iDADwigl)
output <- iDADwigl(Hobbit_MH2T_for_iDAD,
                   nbit = 500,
                   fsum_target = 0.01,
                   U48_0_min = 1.265, 
                   U48_0_max = 1.275, 
                   l = 5.35, 
                   U_0 = 25, 
                   K_min = 1e-13,
                   K_max = 1e-11,
                   T_min = 1e3, 
                   T_max = 20e3)
```

When run on the R console, this function will print a confirmation that the input data frame has the required columns, and print the resulting age value with an error reported as the 67% and 33% quantiles, as follows:

```
All required columns are present in the input data 👍
[1] "Age: 7 +0.6/-0.7 ka"
```

In a typical analysis we will explore the model fit by changing the range of allowed values for the (<sup>234</sup>U/<sup>238</sup>U) ratio at the surface and the age of the sample, and by decreasing the `fsum_target`. Once we obtain a satisfying fit (by visual inspection of the produced figures), we would increase `nbit` to a higher value (e.g. 1000) and run the model again.

## Inspecting the model's output

The model computes a Monte Carlo simulation where age of the sample, U diffusion coefficient and (<sup>234</sup>U/<sup>238</sup>U) ratio at the surface of the sample are taken randomly within the range of values allowed. Results are only kept if the calculated sum of the squared differences between the calculated and observed activity ratios is less than the value set in `fsum_target`. If this is the case, the calculated ratios and the set of solutions for age of the sample, U diffusion coefficient and (<sup>234</sup>U/<sup>238</sup>U) ratio at the surface of the sample are saved. 
The model stops once the number of sets of solutions reaches `nbit`.

The final calculated age `T_final` (in yr), U diffusion coefficient `K_final` (in cm<sup>2</sup>/s) and (<sup>234</sup>U/<sup>238</sup>U) ratio at the surface of the sample `U48_0_final` are the set of solutions where the solution age is the closest to the median age of the population of solutions. 
The uncertainty on each output paramter is calculated as the 67% and 33% quantiles of the population of solution sets.

`T_final`, `K_final` and `U48_0_final` are included in the model's output, along with their uncertainties. 
The function also includes a one-row data frame summarising the age.

The last item in the output is a copy of the input data with two additional columns, the calculated activity ratios, (<sup>234</sup>U/<sup>238</sup>U) and (<sup>230</sup>Th/<sup>238</sup>U), for each measurement location on the sample. 

## Visualising the model's output

`iDADwigl()` returns several figures useful for visualisation of the model results along with the data:
(i) a histogram of the solution ages, 
(ii) the U concentrations in the sample as a function of the relative distance from the center,
(iii) the measured (in blue) and modelled (in red) (<sup>234</sup>U/<sup>238</sup>U) activity ratios as a function of the relative distance from the center, and
(iv) the measured (in blue) and modelled (in red) (<sup>230</sup>Th/<sup>238</sup>U) activity ratios as a function of the relative distance from the center.

