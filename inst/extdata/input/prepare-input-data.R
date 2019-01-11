

# assuming the working directory is the directory containing this file
Hobbit_1_1T_for_iDAD <- read.csv("Hobbit_1-1T_for_iDAD.csv")
Hobbit_MH2T_for_iDAD <- read.csv("Hobbit_MH2T_for_iDAD.csv")


# create pkg data objects for easy re-use
usethis::use_data(Hobbit_1_1T_for_iDAD, Hobbit_1_1T_for_iDAD)
usethis::use_data(Hobbit_MH2T_for_iDAD, Hobbit_MH2T_for_iDAD)
