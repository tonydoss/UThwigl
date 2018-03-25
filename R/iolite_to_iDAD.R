setwd("input")
laser_experiment <- read.csv("laser coordinates.csv") # create dataframe df with data from file

iolite_results <- read.table("data.txt",
                             header = TRUE, sep = "\t", comment.char = "")

file_name_to_save <- paste(sample_name,"for_iDAD.csv", sep = "_")

laser_data_sample <- subset(laser_experiment, (grepl(sample_name, Name) & !grepl('MK', Name)))

iolite_data_sample <- subset(iolite_results, (grepl(sample_name, X) & !grepl('MK', X)))

if ((x_s2 - x_s1) > (y_s2 - y_s1)) {
  iDAD_pos <- as.data.frame(2/(x_s2 - x_s1)*(laser_data_sample$X - x_s2) + 1)
} else {
  iDAD_pos <- as.data.frame(2/(y_s2 - y_s1)*(laser_data_sample$Y - y_s2) + 1)
}

colnames(iDAD_pos) <- c("iDAD position")

df_for_iDAD <- cbind(iDAD_pos, iolite_data_sample[, c(31, 32)],
                     iDAD_pos, iolite_data_sample[, c(33, 34)], iolite_data_sample[, c(29, 30)])

setwd("../output/")
write.table(df_for_iDAD, file = file_name_to_save, sep = ",", row.names = F)
# save.image(file=paste(file_name_to_save,".RData"))
#
# library(ggplot2)
#
# ggplot(df_for_iDAD, aes(iDAD_pos, U_ppm)) + geom_point(size=5) + ggtitle(sample_name) +
#   geom_errorbar(aes(ymin=U_ppm-U_ppm_Int2SE, ymax=U_ppm+U_ppm_Int2SE), width=.05) +
#   ylab("U (ppm)") + xlab("Relative distance from center")
# ggsave(paste("U conc_",sample_name,".pdf", sep = ""))
#
# ggplot(df_for_iDAD, aes(iDAD_pos, U234_U238_CORR)) + geom_point(size=5) + ggtitle(sample_name) +
#   geom_errorbar(aes(ymin=U234_U238_CORR-U234_U238_CORR_Int2SE, ymax=U234_U238_CORR+U234_U238_CORR_Int2SE), width=.05) +
#   ylab(expression("("^234*"U/"^238*"U)")) + xlab("Relative distance from center")
# ggsave(paste("R48_",sample_name,".pdf", sep = ""))
#
# ggplot(df_for_iDAD, aes(iDAD_pos, Th230_U238_CORR)) + geom_point(size=5) + ggtitle(sample_name) +
#   geom_errorbar(aes(ymin=Th230_U238_CORR-Th230_U238_CORR_Int2SE, ymax=Th230_U238_CORR+Th230_U238_CORR_Int2SE), width=.05) +
#   ylab(expression("("^230*"Th/"^238*"U)")) + xlab("Relative distance from center")
# ggsave(paste("R08_",sample_name,".pdf", sep = ""))

