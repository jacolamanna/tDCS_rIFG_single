#Carico le librerie necessarie 
library(readxl)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggsignif)
library(ggplot2)
library(ggpubr)
library(broom)
library(AICcmodavg)
library(rstatix)
library(onewaytests)
library(AID)
library(MASS)
library(energy)
library(agop)
library(car)
library(Matrix)
library(lme4)
library(TestCor)
library(Hmisc)
library(RcmdrMisc)
library(cocor)
library(corrplot)
library(stats) 
library(FactoMineR) 
library(superb)
library(psych)

#Loading Dataset
data_Controlli_TS_LME <- read_excel("Desktop/TS OCD & NIBS/NIBS & Inhibitory Control - UniSR/Preliminar Analysis Plots & Datasets/Dataset_EXP_tDCS_PhD_LME.xlsx", col_names = TRUE)
summary(data_Controlli_TS_LME)
data_Controlli_TS_LME <- as.data.frame(data_Controlli_TS_LME)

#Data Preparation
if("k_Low_Post" %in% colnames(data_Controlli_TS_LME)) {
  data_Controlli_TS_LME$log_k_Low_Post <- log(data_Controlli_TS_LME$k_Low_Post)
} else {
  stop("La colonna 'k_Low_Post' non esiste nel dataframe.")
}
summary(data_Controlli_TS_LME)

if("k_High_Post" %in% colnames(data_Controlli_TS_LME)) {
  data_Controlli_TS_LME$log_k_High_Post <- log(data_Controlli_TS_LME$k_High_Post)
} else {
  stop("La colonna 'k_High_Post' non esiste nel dataframe.")
}
summary(data_Controlli_TS_LME)

if("k_High_Pre" %in% colnames(data_Controlli_TS_LME)) {
  data_Controlli_TS_LME$log_k_High_Pre <- log(data_Controlli_TS_LME$k_High_Pre)
} else {
  stop("La colonna 'k_High_Pre' non esiste nel dataframe.")
}
summary(data_Controlli_TS_LME)

if("k_Low_Pre" %in% colnames(data_Controlli_TS_LME)) {
  data_Controlli_TS_LME$log_k_Low_Pre <- log(data_Controlli_TS_LME$k_Low_Pre)
} else {
  stop("La colonna 'k_Low_Pre' non esiste nel dataframe.")
}
summary(data_Controlli_TS_LME)

data_Controlli_TS_LME <- data_Controlli_TS_LME %>% 
  mutate(Go_noGo_Mean_RT_GoTrials = (Go_noGo_Mean_RT_GoTrials_P + Go_noGo_Mean_RT_GoTrials_R)/2,
         Go_noGo_Mean_RT_noGoTrials = (Go_noGo_Mean_RT_noGOTrials_P + Go_noGo_Mean_RT_noGOTrials_R)/2,
         Go_noGo_Mean_RT = (Go_noGo_Mean_RT_GoTrials_P + Go_noGo_Mean_RT_GoTrials_R + Go_noGo_Mean_RT_noGOTrials_P + Go_noGo_Mean_RT_noGOTrials_R)/4,
         Go_no_Go_Mean_FA = (FA_Block_1 + FA_Block_2)/2,
         Go_noGo_Mean_CR = (CR_Block_1 + CR_Block_2)/2,
         Go_noGO_Mean_d_prime = (d_prime_B_1 + d_prime_B_2)/2)
summary(data_Controlli_TS_LME)

#Corrplots
data_Controlli_TS_LME_Real <- data_Controlli_TS_LME %>%
  dplyr::filter(Protocol == "Real")
data_Controlli_TS_LME_Real_f <- data_Controlli_TS_LME_Real  %>%
  dplyr::select(AUC_Low_Post, AUC_High_Post, log_k_Low_Post, log_k_High_Post, Go_noGo_Mean_RT_GoTrials, Go_no_Go_Mean_FA, SSRT) %>%
  rename(
    `AUC Low` = `AUC_Low_Post`,
    `AUC High` = `AUC_High_Post`,
    `Log(k) Low` = `log_k_Low_Post`,  
    `Log(k) High` = `log_k_High_Post`,                    
    `RT Go Trials` = `Go_noGo_Mean_RT_GoTrials`,
    `FA Rate` = `Go_no_Go_Mean_FA`,
    `SSRT` = `SSRT`)

pairs.panels(data_Controlli_TS_LME_Real_f,
             smooth = TRUE,      # If TRUE, draws loess smooths
             scale = FALSE,      # If TRUE, scales the correlation text font
             density = TRUE,     # If TRUE, adds density plots and histograms
             ellipses = FALSE,   # If TRUE, draws ellipses
             method = "spearman", # Correlation method (also "pearson" or "kendall")
             pch = 21,           # pch symbol
             lm = TRUE,          # If TRUE, plots linear fit rather than the LOESS (smoothed) fit
             cor = TRUE,         # If TRUE, reports correlations
             jiggle = FALSE,     # If TRUE, data points are jittered
             factor = 2,         # Jittering factor
             hist.col = "#F8766D",       # Histograms color
             stars = TRUE,       # If TRUE, adds significance level with stars
             ci = FALSE,         # If TRUE, adds confidence intervals
             cex.cor = 2.5,
             cex.labels = 1.5)      

res <- cor.test(data_Controlli_TS_LME_Real$`log_k_Low_Post`, data_Controlli_TS_LME_Real$`log_k_High_Post`, 
                method = "spearman")
res

data_Controlli_TS_LME_Sham <- data_Controlli_TS_LME %>%
  dplyr::filter(Protocol == "Sham")
data_Controlli_TS_LME_Sham_f <- data_Controlli_TS_LME_Sham  %>%
  dplyr::select(AUC_Low_Post, AUC_High_Post, log_k_Low_Post, log_k_High_Post, Go_noGo_Mean_RT_GoTrials, Go_no_Go_Mean_FA, SSRT) %>%
  rename(
    `AUC Low` = `AUC_Low_Post`,
    `AUC High` = `AUC_High_Post`,
    `Log(k) Low` = `log_k_Low_Post`,  
    `Log(k) High` = `log_k_High_Post`,                    
    `RT Go Trials` = `Go_noGo_Mean_RT_GoTrials`,
    `FA Rate` = `Go_no_Go_Mean_FA`,
    `SSRT` = `SSRT`)

pairs.panels(data_Controlli_TS_LME_Sham_f,
             smooth = TRUE,      # If TRUE, draws loess smooths
             scale = FALSE,      # If TRUE, scales the correlation text font
             density = TRUE,     # If TRUE, adds density plots and histograms
             ellipses = FALSE,   # If TRUE, draws ellipses
             method = "spearman", # Correlation method (also "pearson" or "kendall")
             pch = 21,           # pch symbol
             lm = TRUE,          # If TRUE, plots linear fit rather than the LOESS (smoothed) fit
             cor = TRUE,         # If TRUE, reports correlations
             jiggle = FALSE,     # If TRUE, data points are jittered
             factor = 2,         # Jittering factor
             hist.col = "#00BFC4",       # Histograms color
             stars = TRUE,       # If TRUE, adds significance level with stars
             ci = FALSE,         # If TRUE, adds confidence intervals
             cex.cor = 2.5,
             cex.labels = 1.5)      

res <- cor.test(data_Controlli_TS_LME_Sham$`log_k_Low_Post`, data_Controlli_TS_LME_Sham$`log_k_High_Post`, 
                method = "spearman")
res