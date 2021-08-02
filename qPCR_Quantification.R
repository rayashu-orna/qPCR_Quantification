#!/usr/bin/env Rscript

# fill in the directory address where your data files are
setwd("full/directory/to/your/data/file")


# args <- commandArgs(trailingOnly=TRUE)

# read in qPCR Ct values
# fill in the name of your cleaned up Ct result table
Ct <- data.table::fread("Cleaned_up_Ct_table.txt")

# library size
#lib_size <- as.numeric(args[2])


# make standard curve
# fill in the name of your standard concentration table
Stn <- as.data.frame(Ct)[Ct$SampleType == "STN",]
Stn_concentration <- data.table::fread("STNs.txt")
Stn <- merge(Stn, Stn_concentration, "Sample")
Stn$Log10pM <- log10(Stn$pM)

Stn_fit <- lm(CT ~ Log10pM, data = Stn)
E <- scales::percent_format(accuracy = 0.01)(10^(-1/Stn_fit$coefficients[2]) -1)



# calculate concentrations of test libraries
Ct_test <- data.frame(Ct)[Ct$SampleType == "Test",]

Ct_test$pM  <- 10^((Ct_lib$CT-Stn_fit$coefficients[1])/Stn_fit$coefficients[2])
Ct_test$Log10pM <- log10(Ct_test$pM)
#Ct_lib$nM_undiluted <- Ct_lib$pM * Ct_lib$Dilution / 1000
#Ct_lib$nM_undil_sizeadj <- Ct_lib$nM_undiluted * 399 / lib_size




# plot test sample Ct-concentration against the standard curve
library(ggplot2)
library(ggpmisc)

p <- ggplot(Stn, aes(Log10pM, CT)) +
  geom_smooth(method = "lm", formula = y~x, se = FALSE, colour = "#de3259", size = 0.75) +
  geom_point(shape = 21, colour = "#402d84", size = 1.5, stroke = 0.75) +
  geom_point(data = Ct_test, shape = 3, aes(color = Sample)) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(stat(eq.label), stat(rr.label), sep = "*\", \"*")), 
               parse = TRUE,
               coef.digits = 4,
               rr.digits = 5,
               label.x = "right",
               size = 2) +
  annotate("text", lab = paste("Efficiency = ", E), 
           x = max(Stn$Log10pM), 
           y = max(Stn$CT) - 2.5,
           hjust = "inward",
           size = 2)

pdf("standard_sample.pdf", height = 1.8, width = 3)
  theme_set(theme_bw(base_size = 6))
  print(p)
dev.off()





library(dplyr)
Ct_test %>%
  group_by(Sample) %>%
  summarize(pM.mean = mean(pM))  -> Ct_test

#Ct_lib[,2:3] <- round(Ct_lib[,2:3], digits = 3)

write.table(Ct_test, "Sample.Conc.txt", quote = F, row.names = F)



