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
                              round((quantile(T_sol,.67)-T_final)/1000, digits = 1),"/-",
                              round((T_final - quantile(T_sol,.33))/1000, digits = 1)," ka", sep = ""))

setwd(paste(path_wd,"output/", sep = ""))

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
                              round((quantile(T_sol,.67)-T_final)/1000, digits = 1),"/-",
                              round((T_final - quantile(T_sol,.33))/1000, digits = 1)," ka", sep = ""))

ggsave(paste("R48_",sample_name,"_DAD.png", sep = ""), width = 12, height = 10, units = "cm", dpi = 300,)

save.image(paste(sample_name,".RData", sep = ""))

results <- as.data.frame(cbind(T_final/1000, (quantile(T_sol,.67)-T_final)/1000, (T_final - quantile(T_sol,.33))/1000, 
                         U48_0_final, U48_0_max - U48_0_final, U48_0_final - U48_0_min))

colnames(results) <- c("Age (ka)", "Age 67% quantile (ka)" ,"Age 33% quantile (ka)", 
                       "U234_U238_0", "U234_U238_0 67% quantile" ,"U234_U238_0 33% quantile")
rownames(results) <- c("Results")

write.table(results, file = paste(sample_name,"_model_results.csv", sep = ""), sep = ",", row.names = F)

calc_ratios <- as.data.frame(cbind(U48calc_final, Th0U8calc_final))
colnames(calc_ratios) <- c("U234_U238_CALC", "Th230_U238_CALC")

write.table(calc_ratios, file = paste(sample_name,"_calc_ratios.csv", sep = ""), sep = ",", row.names = F)

print(paste("Age: ",round(T_final/1000, digits=1)," +",
            round((quantile(T_sol,.67) - T_final)/1000, digits = 1),"/-",
            round((T_final - quantile(T_sol,.33))/1000, digits = 1)," ka", sep = ""))

setwd(path_wd)
