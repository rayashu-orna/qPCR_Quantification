#!/usr/bin/env Rscript

setwd("C:/Users/IanLanger/OneDrive - Orna Therapeutics/Desktop/qPCR/qPCR data/072721_2")



# args <- commandArgs(trailingOnly=TRUE)

# read in qPCR Ct values
Ct <- data.table::fread("WellResult_073021.txt")

# library size
#lib_size <- as.numeric(args[2])


# make standard curve
Stn <- as.data.frame(Ct)[Ct$SampleType == "STN",]
Stn_concentration <- data.table::fread("STNs.txt")
Stn <- merge(Stn, Stn_concentration, "Sample")
Stn$Log10pM <- log10(Stn$pM)

Stn_fit <- lm(Ct ~ Log10pM, data = Stn)
E <- scales::percent_format(accuracy = 0.01)(10^(-1/Stn_fit$coefficients[2]) -1)

library(ggplot2)
library(ggpmisc)

p <- ggplot(Stn, aes(Log10pM, Ct)) +
      geom_smooth(method = "lm", formula = y~x, se = FALSE, colour = "#de3259", size = 0.75) +
      geom_point(shape = 21, colour = "#402d84", size = 1.5, stroke = 0.75) +
      stat_poly_eq(formula = y ~ x, 
                   aes(label = paste(stat(eq.label), stat(rr.label), sep = "*\", \"*")), 
                   parse = TRUE,
                   coef.digits = 4,
                   rr.digits = 5,
                   label.x = "right",
                   size = 2) +
      annotate("text", lab = paste("Efficiency = ", E), 
               x = max(Stn$Log10pM), 
               y = max(Stn$Ct) - 1.5,
               hjust = "inward",
               size = 2)

pdf("standard_curve.pdf", height = 1.8, width = 2.5)
  theme_set(theme_bw(base_size = 6))
  print(p)
dev.off()


# calculate concentrations of test libraries
Ct_lib <- data.frame(Ct)[Ct$SampleType == "library",]

Ct_lib$pM  <- 10^((Ct_lib$Ct-Stn_fit$coefficients[1])/Stn_fit$coefficients[2])
#Ct_lib$nM_undiluted <- Ct_lib$pM * Ct_lib$Dilution / 1000
#Ct_lib$nM_undil_sizeadj <- Ct_lib$nM_undiluted * 399 / lib_size

library(dplyr)
Ct_lib %>%
  group_by(Sample) %>%
  summarize(pM.mean = mean(pM))  -> Ct_lib

#Ct_lib[,2:3] <- round(Ct_lib[,2:3], digits = 3)

write.table(Ct_lib, "Sample.Conc.tsv", quote = F, row.names = F)
