library(ggplot2)


theme_plots <- theme(plot.title = element_text(size = (18), hjust = 0.5),
                     legend.title=element_blank(),
                     axis.title.y = element_text(size=18),
                     axis.title.x = element_text(size=18),
                     axis.text.x=element_text(size=14),
                     axis.text.y=element_text(size=14)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

df_T_sol <- as.data.frame(T_sol)

ggplot(df_T_sol, aes(T_sol)) + geom_histogram(binwidth = 500, fill="white", color="black") + theme_plots

ggplot(df) +
  geom_errorbar(aes(x= iDAD.position, ymax = Th230_U238_CORR+Th230_U238_CORR_Int2SE, ymin = Th230_U238_CORR-Th230_U238_CORR_Int2SE, width=0.02)) + # plot error bars
  geom_point(aes(iDAD.position,Th230_U238_CORR), color='blue', size=5) +
  geom_point(aes(iDAD.position,Th0U8calc_final), color='red', size=5) +
  ylab(expression("("^230*"Th/"^238*"U)")) + xlab("Relative distance from center")+
  theme_plots + ggtitle(paste("Age: ",round(T_final/1000, digits=1)," +",
                              round((quantile(T_sol,.95)-T_final)/1000, digits = 1),"/-",
                              round((T_final - quantile(T_sol,.05))/1000, digits = 1)," ka", sep = ""))

setwd("output/")
ggsave(paste("R08_",sample_name,"_DAD.png", sep = ""), width = 12, height = 10, units = "cm", dpi = 300,)

ggplot(df) +
  geom_errorbar(aes(x= iDAD.position, ymax = U_ppm+U_ppm_Int2SE, ymin = U_ppm-U_ppm_Int2SE, width=0.02)) + # plot error bars
  geom_point(aes(iDAD.position,U_ppm), color='blue', size=5) +
  ylab("U (ppm)") + xlab("Relative distance from center") +
  theme_plots

ggsave(paste("U_",sample_name,"_DAD.png", sep = ""), width = 12, height = 10, units = "cm", dpi = 300,)

ggplot(df) +
  geom_errorbar(aes(x= iDAD.position, ymax = U234_U238_CORR+U234_U238_CORR_Int2SE, ymin = U234_U238_CORR-U234_U238_CORR_Int2SE, width=0.02)) + # plot error bars
  geom_point(aes(iDAD.position,U234_U238_CORR), color='blue', size=5) +
  geom_point(aes(iDAD.position,U48calc_final), color='red', size=5) +
  ylab(expression("("^234*"U/"^238*"U)")) + xlab("Relative distance from center") +
  theme_plots + ggtitle(paste("Age: ",round(T_final/1000, digits=1)," +",
                              round((quantile(T_sol,.95)-T_final)/1000, digits = 1),"/-",
                              round((T_final - quantile(T_sol,.05))/1000, digits = 1)," ka", sep = ""))

ggsave(paste("R48_",sample_name,"_DAD.png", sep = ""), width = 12, height = 10, units = "cm", dpi = 300,)

save.image(paste(sample_name,".RData", sep = ""))

print(paste("Age: ",round(T_final/1000, digits=1)," +",
            round((quantile(T_sol,.95)-T_final)/1000, digits = 1),"/-",
            round((T_final - quantile(T_sol,.05))/1000, digits = 1)," ka", sep = ""))

setwd("../")
