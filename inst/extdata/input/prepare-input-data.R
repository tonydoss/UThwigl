

# assuming the working directory is the directory containing this file
Hobbit_1_1T_for_iDAD <- read.csv("Hobbit_1-1T_for_iDAD.csv")
Hobbit_MH2T_for_iDAD <- read.csv("Hobbit_MH2T_for_iDAD.csv")
iolite_export <- read.csv("example_data_csUTh.csv")


# create pkg data objects for easy re-use
usethis::use_data(Hobbit_1_1T_for_iDAD, Hobbit_1_1T_for_iDAD, overwrite = TRUE)
usethis::use_data(Hobbit_MH2T_for_iDAD, Hobbit_MH2T_for_iDAD, overwrite = TRUE)
usethis::use_data(iolite_export, iolite_export, overwrite = TRUE)
