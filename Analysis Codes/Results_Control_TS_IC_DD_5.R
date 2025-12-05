#Loading Libraries 
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
library(pbkrtest)
library(emmeans)
library(lmerTest)

#Loading Dataset
data_Controlli_TS_LME_Ext <- read_excel("Desktop/TS OCD & NIBS/NIBS & Inhibitory Control - UniSR/Preliminar Analysis Plots & Datasets/Dataset_EXP_tDCS_PhD_LME_Ext.xlsx", col_names = TRUE)
summary(data_Controlli_TS_LME_Ext)
data_Controlli_TS_LME_Ext <- as.data.frame(data_Controlli_TS_LME_Ext)

#Data Preparation
if("k_Post" %in% colnames(data_Controlli_TS_LME_Ext)) {
  data_Controlli_TS_LME_Ext$log_k_Post <- log(data_Controlli_TS_LME_Ext$k_Post)
} else {
  stop("La colonna 'k_Post' non esiste nel dataframe.")
}
summary(data_Controlli_TS_LME_Ext)

if("k_Pre" %in% colnames(data_Controlli_TS_LME_Ext)) {
  data_Controlli_TS_LME_Ext$log_k_Pre <- log(data_Controlli_TS_LME_Ext$k_Pre)
} else {
  stop("La colonna 'k_Pre' non esiste nel dataframe.")
}
summary(data_Controlli_TS_LME_Ext)

data_Controlli_TS_LME_Ext <- data_Controlli_TS_LME_Ext %>% 
  mutate(Go_noGo_Mean_RT_GoTrials = (Go_noGo_Mean_RT_GoTrials_P + Go_noGo_Mean_RT_GoTrials_R)/2,
         Go_noGo_Mean_RT_noGoTrials = (Go_noGo_Mean_RT_noGOTrials_P + Go_noGo_Mean_RT_noGOTrials_R)/2,
         Go_noGo_Mean_RT = (Go_noGo_Mean_RT_GoTrials_P + Go_noGo_Mean_RT_GoTrials_R + Go_noGo_Mean_RT_noGOTrials_P + Go_noGo_Mean_RT_noGOTrials_R)/4,
         Go_no_Go_Mean_FA = (FA_Block_1 + FA_Block_2)/2,
         Go_noGo_Mean_CR = (CR_Block_1 + CR_Block_2)/2,
         Go_noGO_Mean_d_prime = (d_prime_B_1 + d_prime_B_2)/2)
summary(data_Controlli_TS_LME_Ext)

##### LME SSRT #####
mod2 <- lmer(SSRT ~ 1 + Protocol + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + log_k_Pre + (1 | Id), data = data_Controlli_TS_LME_Ext, REML = FALSE)
mod1 <- lmer(SSRT ~ 1 + Protocol + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + log_k_Pre + Protocol_ord + Protocol:(MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + log_k_Pre) + (1 | Id), data = data_Controlli_TS_LME_Ext, REML = FALSE)

Anova(mod2, type = 'II')
Anova(mod2, type = 'II', test="F")     #,REML = TRUE in mod2
qqp(data_Controlli_TS_LME_Ext$SSRT)
qqp(resid(mod2), "norm")

Anova(mod1, type = 'III')
qqp(data_Controlli_TS_LME_Ext$SSRT)
qqp(resid(mod1), "norm")

#Kenward-Rogers
anova(mod2,mod1)
KRmodcomp(mod2,mod1)

emmeans_result <- emmeans(mod2, ~ Protocol)
print(emmeans_result)

aov_1 <- Anova(mod2, type = 'II')
pvals <- aov_1$`Pr(>Chisq)`
pvals_corr <- p.adjust(pvals, method = "fdr")
aov_1$`Pr(>Chisq)_corr_FDR` <- pvals_corr
aov_1

#Add Fitted data to dataset
data_Controlli_TS_LME_Ext$Fitted_SSRT <- fitted(mod2)
summary(data_Controlli_TS_LME_Ext)
data_Controlli_TS_LME_Ext_Filtered <- data_Controlli_TS_LME_Ext %>%
  filter(Reward == "High")

#Boxplot with Means tDCS ~ SSRT
P1_F <- ggplot(data_Controlli_TS_LME_Ext_Filtered, aes(x = Protocol, y = Fitted_SSRT, color = Protocol)) +
  geom_point(position = position_jitter(width = 0.2)) +
  geom_boxplot(color = "black", fill = NA, width = 0.2, outlier.shape = NA, position = position_dodge(width = 0.5), fatten = NULL) +
  stat_summary(aes(x = Protocol, y = Fitted_SSRT), fun = mean, geom = "crossbar", width = 0.2, color = "black", linetype = "solid", linewidth = 0.3, position = position_dodge(width = 0.5)) +
  labs(x = NULL, y = "Stop-Signal Reaction Time")  +
  theme_classic() +
  scale_y_continuous(limits = c(170, 295), breaks = seq(170, 295, by = 20)) +
  showSignificance( c(1,2), 293, -3.0) +
  annotate("text", x = mean(c(1,2)), y = 295, label = "***", size = 15) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),  
        legend.title = element_text(size = 18),
        legend.position = "top")

#MSAS - SSRT (NON UTILIZZATO NELL'ARTICOLO)
mod2 <- lmer(SSRT ~ 1 + Protocol + MSAS_18_sr + MSAS_18_cd + MSAS_18_so + MSAS_18_m + log_k_Pre + (1 | Id), data = data_Controlli_TS_LME_Ext, REML = FALSE)
mod1 <- lmer(SSRT ~ 1 + Protocol + MSAS_18_sr + MSAS_18_cd + MSAS_18_so + MSAS_18_m + log_k_Pre + Protocol:(MSAS_18_sr + MSAS_18_cd + MSAS_18_so + MSAS_18_m + log_k_Pre) + (1 | Id), data = data_Controlli_TS_LME_Ext, REML = FALSE)

Anova(mod2, type = 'III')
qqp(data_Controlli_TS_LME_Ext$SSRT)
qqp(resid(mod2), "norm")

Anova(mod1, type = 'III')
qqp(data_Controlli_TS_LME_Ext$SSRT)
qqp(resid(mod1), "norm")

#Kenward-Rogers
anova(mod2,mod1)
KRmodcomp(mod2,mod1)

##### LME GONOGO - RT GO TRIALS #####
mod2 <- lmer(Go_noGo_Mean_RT_GoTrials ~ 1 + Protocol + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + log_k_Pre + (1 | Id), data = data_Controlli_TS_LME_Ext, REML = FALSE)
mod1 <- lmer(Go_noGo_Mean_RT_GoTrials ~ 1 + Protocol + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + log_k_Pre + Protocol_ord + Protocol:(MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + log_k_Pre) + (1 | Id), data = data_Controlli_TS_LME_Ext, REML = FALSE)

Anova(mod2, type = 'III')
Anova(mod2, type = 'III', test="F")
qqp(data_Controlli_TS_LME_Ext$Go_noGo_Mean_RT_GoTrials)
qqp(resid(mod2), "norm")

Anova(mod1, type = 'III')
Anova(mod1, type = 'III', test="F")
qqp(data_Controlli_TS_LME_Ext$Go_noGo_Mean_RT_GoTrials)
qqp(resid(mod1), "norm")

#Kenward-Rogers
anova(mod2,mod1)
KRmodcomp(mod2,mod1)

emmeans_result <- emmeans(mod1, ~ Protocol)
print(emmeans_result)

aov_1 <- Anova(mod1, type = 'III')
pvals <- aov_1$`Pr(>Chisq)`
pvals_corr <- p.adjust(pvals, method = "fdr")
aov_1$`Pr(>Chisq)_corr_FDR` <- pvals_corr
aov_1

#Add Fitted data to dataset
data_Controlli_TS_LME_Ext$Fitted_Go_noGo_Mean_RT_GoTrials <- fitted(mod1)
summary(data_Controlli_TS_LME_Ext)
data_Controlli_TS_LME_Ext_Filtered <- data_Controlli_TS_LME_Ext %>%
  filter(Reward == "High")

#Boxplot with Means tDCS ~ RT GO TRIALS
P2_F <- ggplot(data_Controlli_TS_LME_Ext_Filtered, aes(x = Protocol, y = Fitted_Go_noGo_Mean_RT_GoTrials, color = Protocol)) +
  geom_point(position = position_jitter(width = 0.2)) +
  geom_boxplot(color = "black", fill = NA, width = 0.2, outlier.shape = NA, position = position_dodge(width = 0.5), fatten = NULL) +
  stat_summary(aes(x = Protocol, y = Fitted_Go_noGo_Mean_RT_GoTrials), fun = mean, geom = "crossbar", width = 0.2, color = "black", linetype = "solid", linewidth = 0.3, position = position_dodge(width = 0.5)) +
  labs(x = NULL, y = "Reaction Time Go Trials") +
  theme_classic() +
  scale_y_continuous(limits = c(250, 510), breaks = seq(250, 510, by = 50)) +
  showSignificance( c(1,2), 500, -5.5) +
  annotate("text", x = mean(c(1,2)), y = 510, label = "ns", size = 6) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),  
        legend.title = element_text(size = 18),
        legend.position = "top")

#MSAS - RT GO TRIALS (NON UTILIZZATO NELL'ARTICOLO)
mod2 <- lmer(Go_noGo_Mean_RT_GoTrials ~ 1 + Protocol + MSAS_18_sr + MSAS_18_cd + MSAS_18_so + MSAS_18_m + log_k_Pre + (1 | Id), data = data_Controlli_TS_LME_Ext, REML = FALSE)
mod1 <- lmer(Go_noGo_Mean_RT_GoTrials ~ 1 + Protocol + MSAS_18_sr + MSAS_18_cd + MSAS_18_so + MSAS_18_m + log_k_Pre + Protocol:(MSAS_18_sr + MSAS_18_cd + MSAS_18_so + MSAS_18_m + log_k_Pre) + (1 | Id), data = data_Controlli_TS_LME_Ext, REML = FALSE)

Anova(mod2, type = 'III')
qqp(data_Controlli_TS_LME_Ext$Go_noGo_Mean_RT_GoTrials)
qqp(resid(mod2), "norm")

Anova(mod1, type = 'III')
qqp(data_Controlli_TS_LME_Ext$Go_noGo_Mean_RT_GoTrials)
qqp(resid(mod1), "norm")

#Kenward-Rogers
anova(mod2,mod1)
KRmodcomp(mod2,mod1)

##### LME GONOGO - FA RATE #####
mod2 <- lmer(Go_no_Go_Mean_FA ~ 1 + Protocol + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + log_k_Pre + (1 | Id), data = data_Controlli_TS_LME_Ext, REML = FALSE)
mod1 <- lmer(Go_no_Go_Mean_FA ~ 1 + Protocol + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + log_k_Pre + Protocol_ord + Protocol:(MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + log_k_Pre) + (1 | Id), data = data_Controlli_TS_LME_Ext, REML = FALSE)

Anova(mod2, type = 'III')
qqp(data_Controlli_TS_LME_Ext$Go_no_Go_Mean_FA)
qqp(resid(mod2), "norm")

Anova(mod1, type = 'III')
Anova(mod1, type = 'III', test="F")
qqp(data_Controlli_TS_LME_Ext$Go_no_Go_Mean_FA)
qqp(resid(mod1), "norm")

#Kenward-Rogers
anova(mod2,mod1)
KRmodcomp(mod2,mod1)

emmeans_result <- emmeans(mod1, ~ Protocol)
print(emmeans_result)

aov_1 <- Anova(mod1, type = 'III')
pvals <- aov_1$`Pr(>Chisq)`
pvals_corr <- p.adjust(pvals, method = "fdr")
aov_1$`Pr(>Chisq)_corr_FDR` <- pvals_corr
aov_1

#Add Fitted data to dataset
data_Controlli_TS_LME_Ext$Fitted_Go_no_Go_Mean_FA <- fitted(mod1)
summary(data_Controlli_TS_LME_Ext)
data_Controlli_TS_LME_Ext_Filtered <- data_Controlli_TS_LME_Ext %>%
  filter(Reward == "High")

P3_F <- ggplot(data_Controlli_TS_LME_Ext_Filtered, aes(x = Protocol, y = Fitted_Go_no_Go_Mean_FA, color = Protocol)) +
  geom_point(position = position_jitter(width = 0.2)) +
  geom_boxplot(color = "black", fill = NA, width = 0.2, outlier.shape = NA, position = position_dodge(width = 0.5), fatten = NULL) +
  stat_summary(aes(x = Protocol, y = Fitted_Go_no_Go_Mean_FA), fun = mean, geom = "crossbar", width = 0.2, color = "black", linetype = "solid", linewidth = 0.3, position = position_dodge(width = 0.5)) +
  labs(x = NULL, y = "False Alarm Rate") +
  theme_classic() +
  scale_y_continuous(limits = c(0.0, 0.19), breaks = seq(0.0, 0.19, by = 0.03)) +
  showSignificance( c(1,2), 0.18, -0.004) +
  annotate("text", x = mean(c(1,2)), y = 0.188, label = "ns", size = 6) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),  
        legend.title = element_text(size = 18),
        legend.position = "top")

#MSAS - FA Rate (NON UTILIZZATO NELL'ARTICOLO)
mod2 <- lmer(Go_no_Go_Mean_FA ~ 1 + Protocol + MSAS_18_sr + MSAS_18_cd + MSAS_18_so + MSAS_18_m + log_k_Pre + (1 | Id), data = data_Controlli_TS_LME_Ext, REML = FALSE)
mod1 <- lmer(Go_no_Go_Mean_FA ~ 1 + Protocol + MSAS_18_sr + MSAS_18_cd + MSAS_18_so + MSAS_18_m + log_k_Pre + Protocol:(MSAS_18_sr + MSAS_18_cd + MSAS_18_so + MSAS_18_m + log_k_Pre) + (1 | Id), data = data_Controlli_TS_LME_Ext, REML = FALSE)

Anova(mod2, type = 'III')
qqp(data_Controlli_TS_LME_Ext$Go_no_Go_Mean_FA)
qqp(resid(mod2), "norm")

Anova(mod1, type = 'III')
qqp(data_Controlli_TS_LME_Ext$Go_no_Go_Mean_FA)
qqp(resid(mod1), "norm")

#Kenward-Rogers
anova(mod2,mod1)
KRmodcomp(mod2,mod1)

##### LME GONOGO - ACCURACY #####
mod2 <- lmer(Go_noGo_Mean_Accuracy ~ 1 + Protocol + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + log_k_Pre + (1 | Id), data = data_Controlli_TS_LME_Ext, REML = FALSE)
mod1 <- lmer(Go_noGo_Mean_Accuracy ~ 1 + Protocol + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + log_k_Pre + Protocol_ord + Protocol:(MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + log_k_Pre) + (1 | Id), data = data_Controlli_TS_LME_Ext, REML = TRUE)

Anova(mod2, type = 'III')
qqp(data_Controlli_TS_LME_Ext$Go_no_Go_Mean_FA)
qqp(resid(mod2), "norm")

Anova(mod1, type = 'III')
Anova(mod1, type = 'III', test="F")
qqp(data_Controlli_TS_LME_Ext$Go_no_Go_Mean_FA)
qqp(resid(mod1), "norm")

#Kenward-Rogers
anova(mod2,mod1)
KRmodcomp(mod2,mod1)

emmeans_result <- emmeans(mod1, ~ Protocol)
print(emmeans_result)

aov_1 <- Anova(mod1, type = 'III')
pvals <- aov_1$`Pr(>Chisq)`
pvals_corr <- p.adjust(pvals, method = "fdr")
aov_1$`Pr(>Chisq)_corr_FDR` <- pvals_corr
aov_1

#Add Fitted data to dataset
data_Controlli_TS_LME_Ext$Fitted_Go_noGo_Mean_Accuracy <- fitted(mod1)
summary(data_Controlli_TS_LME_Ext)
data_Controlli_TS_LME_Ext_Filtered <- data_Controlli_TS_LME_Ext %>%
  filter(Reward == "High")

P4_F <- ggplot(data_Controlli_TS_LME_Ext_Filtered, aes(x = Protocol, y = Fitted_Go_noGo_Mean_Accuracy, color = Protocol)) +
  geom_point(position = position_jitter(width = 0.2)) +
  geom_boxplot(color = "black", fill = NA, width = 0.2, outlier.shape = NA, position = position_dodge(width = 0.5), fatten = NULL) +
  stat_summary(aes(x = Protocol, y = Fitted_Go_noGo_Mean_Accuracy), fun = mean, geom = "crossbar", width = 0.2, color = "black", linetype = "solid", linewidth = 0.3, position = position_dodge(width = 0.5)) +
  labs(x = NULL, y = "Accuracy") +
  theme_classic() +
  scale_y_continuous(limits = c(0.75, 1.005), breaks = seq(0.7, 1.005, by = 0.05)) +
  showSignificance( c(1,2), 1.0, -0.005) +
  annotate("text", x = mean(c(1,2)), y = 1.005, label = "*", size = 15) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),  
        legend.title = element_text(size = 18),
        legend.position = "top")

#MSAS - Accuracy (NON UTILIZZATO NELL'ARTICOLO)
mod2 <- lmer(Go_noGo_Mean_Accuracy ~ 1 + Protocol + MSAS_18_sr + MSAS_18_cd + MSAS_18_so + MSAS_18_m + log_k_Pre + (1 | Id), data = data_Controlli_TS_LME_Ext, REML = FALSE)
mod1 <- lmer(Go_noGo_Mean_Accuracy ~ 1 + Protocol + MSAS_18_sr + MSAS_18_cd + MSAS_18_so + MSAS_18_m + log_k_Pre + Protocol:(MSAS_18_sr + MSAS_18_cd + MSAS_18_so + MSAS_18_m + log_k_Pre) + (1 | Id), data = data_Controlli_TS_LME_Ext, REML = FALSE)

Anova(mod2, type = 'III')
qqp(data_Controlli_TS_LME_Ext$Go_no_Go_Mean_FA)
qqp(resid(mod2), "norm")

Anova(mod1, type = 'III')
qqp(data_Controlli_TS_LME_Ext$Go_no_Go_Mean_FA)
qqp(resid(mod1), "norm")

#Kenward-Rogers
anova(mod2,mod1)
KRmodcomp(mod2,mod1)

emmeans_result <- emmeans(mod1, ~ Protocol)
print(emmeans_result)

#Compongo Figura 1 (a,b,c,d)
#P1_F <- P1_F + theme(legend.position = "none", aspect.ratio = 1)
#P3_F <- P3_F + theme(legend.position = "none", aspect.ratio = 1)
library("cowplot")
ggdraw() +
  draw_plot(P1_F, x = 0, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(P2_F, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(P3_F, x = 0, y = 0, width = 0.5, height = 0.5) +
  draw_plot(P4_F, x = 0.5, y = 0, width = 0.5, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 30,
                  x = c(0, 0.5, 0, 0.5), 
                  y = c(1, 1, 0.5, 0.5))

##### LME DD PRE-POST #####
data_Controlli_TS_LME_DD_Ext <- read_excel("Desktop/TS OCD & NIBS/NIBS & Inhibitory Control - UniSR/Preliminar Analysis Plots & Datasets/Dataset_EXP_tDCS_PhD_LME_Ext_PRE_POST.xlsx", col_names = TRUE)
summary(data_Controlli_TS_LME_DD_Ext)
data_Controlli_TS_LME_DD_Ext <- as.data.frame(data_Controlli_TS_LME_DD_Ext)

#Data Preparation
if("k" %in% colnames(data_Controlli_TS_LME_DD_Ext)) {
  data_Controlli_TS_LME_DD_Ext$log_k <- log(data_Controlli_TS_LME_DD_Ext$k)
} else {
  stop("La colonna 'k' non esiste nel dataframe.")
}
summary(data_Controlli_TS_LME_DD_Ext)

data_Controlli_TS_LME_DD_Ext <- data_Controlli_TS_LME_DD_Ext %>% 
  mutate(Go_noGo_Mean_RT_GoTrials = (Go_noGo_Mean_RT_GoTrials_P + Go_noGo_Mean_RT_GoTrials_R)/2,
         Go_noGo_Mean_RT_noGoTrials = (Go_noGo_Mean_RT_noGOTrials_P + Go_noGo_Mean_RT_noGOTrials_R)/2,
         Go_noGo_Mean_RT = (Go_noGo_Mean_RT_GoTrials_P + Go_noGo_Mean_RT_GoTrials_R + Go_noGo_Mean_RT_noGOTrials_P + Go_noGo_Mean_RT_noGOTrials_R)/4,
         Go_no_Go_Mean_FA = (FA_Block_1 + FA_Block_2)/2,
         Go_noGo_Mean_CR = (CR_Block_1 + CR_Block_2)/2,
         Go_noGO_Mean_d_prime = (d_prime_B_1 + d_prime_B_2)/2)
summary(data_Controlli_TS_LME_DD_Ext)

subset_low <- data_Controlli_TS_LME_DD_Ext %>% 
  filter(Reward == "Low")

subset_high <- data_Controlli_TS_LME_DD_Ext %>% 
  filter(Reward == "High")

mod2 <- lmer(log_k ~ 1 + Protocol*Time_DD + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + (1|Id), data = subset_low,
                   REML = FALSE)
mod1 <- lmer(log_k ~ 1 + Protocol*Time_DD + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + Protocol_ord + Protocol:(MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc) + (1 | Id), data = subset_low,
             REML = FALSE)

Anova(mod2, type = 'III')
qqp(subset_low$log_k)
qqp(resid(mod2), "norm")

Anova(mod1, type = 'III')
Anova(mod1, type = 'III', test="F")
qqp(subset_low$log_k)
qqp(resid(mod1), "norm")

#Kenward-Rogers
anova(mod2,mod1)
KRmodcomp(mod2,mod1)

aov_1 <- Anova(mod1, type = 'III')
pvals <- aov_1$`Pr(>Chisq)`
pvals_corr <- p.adjust(pvals, method = "fdr")
aov_1$`Pr(>Chisq)_corr_FDR` <- pvals_corr
aov_1

#Add Fitted data to dataset
subset_low$Fitted_log_k <- fitted(mod2)
summary(subset_low)
#Qui avevo provato a fare un boxplot con i dati fitatti ma non sono riuscito
#Create new datset with Fitted Data
write_xlsx(subset_low, path = "Desktop/subset_low.xlsx")
data_Controlli_TS_LME_DD_Low <- read_excel("Desktop/subset_low_G.xlsx", col_names = TRUE)
data_Controlli_TS_LME_DD_Low <- as.data.frame(data_Controlli_TS_LME_DD_Low)
data_Controlli_TS_LME_DD_Low$log_k_diff <- data_Controlli_TS_LME_DD_Low$F_log_Post - data_Controlli_TS_LME_DD_Low$F_log_Pre
summary(data_Controlli_TS_LME_DD_Low)
data_Controlli_TS_LME_DD_Low$Protocol <- as.factor(data_Controlli_TS_LME_DD_Low$Protocol)

P5_FF <- ggplot(data_Controlli_TS_LME_DD_Low, aes(x = Protocol, y = log_k_diff, color = Protocol)) +
  geom_point(aes(x = Protocol, y = log_k_diff, color = Protocol), position = position_jitter(width = 0.05)) +
  geom_boxplot(aes(x = Protocol, y = log_k_diff), color = "black", fill = NA, width = 0.2, outlier.shape = NA, position = position_dodge(width = 0.5), fatten = NULL) +
  stat_summary(aes(x = Protocol, y = log_k_diff), fun = mean, geom = "crossbar", width = 0.2, color = "black", linetype = "solid", linewidth = 0.3, position = position_dodge(width = 0.5)) +
  labs(x = NULL, y = "Log(k) Post - Log(k) Pre") +
  theme_classic() +
  scale_y_continuous(limits = c(-0.9, 0.65), breaks = seq(-0.9, 0.65, by = 0.3), labels = scales::number_format(accuracy = 0.1)) +
  showSignificance( c(1,2), 0.6, -0.05) +
  annotate("text", x = mean(c(1,2)), y = 0.63, label = "*", size = 15) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),  
        legend.title = element_text(size = 18),
        legend.position = "top")

mod2 <- lmer(log_k ~ 1 + Protocol*Time_DD + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + (1|Id), data = subset_high,
             REML = FALSE)
mod1 <- lmer(log_k ~ 1 + Protocol*Time_DD + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + Protocol_ord + Protocol:(MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc) + (1 | Id), data = subset_high,
             REML = FALSE)

Anova(mod2, type = 'III')
qqp(subset_low$log_k)
qqp(resid(mod2), "norm")

Anova(mod1, type = 'III')
Anova(mod1, type = 'III', test="F")
qqp(subset_low$log_k)
qqp(resid(mod1), "norm")

#Kenward-Rogers
anova(mod2,mod1)
KRmodcomp(mod2,mod1)

aov_1 <- Anova(mod1, type = 'III')
pvals <- aov_1$`Pr(>Chisq)`
pvals_corr <- p.adjust(pvals, method = "fdr")
aov_1$`Pr(>Chisq)_corr_FDR` <- pvals_corr
aov_1

#Add Fitted data to dataset
subset_high$Fitted_log_k <- fitted(mod2)
summary(subset_high)
#Qui avevo provato a fare un boxplot con i dati fitatti ma non sono riuscito
#Create new datset with Fitted Data
write_xlsx(subset_high, path = "Desktop/subset_high.xlsx")
data_Controlli_TS_LME_DD_High <- read_excel("Desktop/subset_high_G.xlsx", col_names = TRUE)
data_Controlli_TS_LME_DD_High <- as.data.frame(data_Controlli_TS_LME_DD_High)
data_Controlli_TS_LME_DD_High$log_k_diff <- data_Controlli_TS_LME_DD_High$F_log_Post - data_Controlli_TS_LME_DD_High$F_log_Pre
summary(data_Controlli_TS_LME_DD_High)
data_Controlli_TS_LME_DD_High$Protocol <- as.factor(data_Controlli_TS_LME_DD_High$Protocol)

P6_FF <- ggplot(data_Controlli_TS_LME_DD_High, aes(x = Protocol, y = log_k_diff, color = Protocol)) +
  geom_point(aes(x = Protocol, y = log_k_diff, color = Protocol), position = position_jitter(width = 0.05)) +
  geom_boxplot(aes(x = Protocol, y = log_k_diff), color = "black", fill = NA, width = 0.2, outlier.shape = NA, position = position_dodge(width = 0.5), fatten = NULL) +
  stat_summary(aes(x = Protocol, y = log_k_diff), fun = mean, geom = "crossbar", width = 0.2, color = "black", linetype = "solid", linewidth = 0.3, position = position_dodge(width = 0.5)) +
  labs(x = NULL, y = "Log(k) Post - Log(k) Pre") +
  theme_classic() +
  scale_y_continuous(limits = c(-0.8, 0.25), breaks = seq(-0.8, 0.25, by = 0.2), labels = scales::number_format(accuracy = 0.1)) +
  showSignificance( c(1,2), 0.2, -0.035) +
  annotate("text", x = mean(c(1,2)), y = 0.24, label = "ns", size = 5) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),  
        legend.title = element_text(size = 18),
        legend.position = "top")

#Creo Figura 3 (a, b, c, d)
library("cowplot")
P5_FF <- P5_FF + theme(legend.position = "none")
P6_FF <- P6_FF + theme(legend.position = "none")
ggdraw() +
  draw_plot(PDD1_F, x = 0, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(PDD2_F, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(P6_FF, x = 0, y = 0, width = 0.5, height = 0.5) +
  draw_plot(P5_FF, x = 0.5, y = 0, width = 0.5, height = 0.5) +
  draw_plot_label(label = c("a", "b", "c", "d"), size = 30,
                  x = c(0, 0.5, 0, 0.5), 
                  y = c(1, 1, 0.5, 0.5))


#MSAS Models (NON UTILIZZATI NELL'ARTICOLO)
mod2 <- lmer(log_k ~ 1 + Protocol*Time_DD + MSAS_18_sr + MSAS_18_cd + MSAS_18_so + MSAS_18_m + (1|Id), data = subset_low,
             REML = FALSE)
mod1 <- lmer(log_k ~ 1 + Protocol*Time_DD + MSAS_18_sr + MSAS_18_cd + MSAS_18_so + MSAS_18_m + Protocol:(MSAS_18_sr + MSAS_18_cd + MSAS_18_so + MSAS_18_m) + (1 | Id), data = subset_low,
             REML = FALSE)

Anova(mod2, type = 'III')
qqp(subset_low$log_k)
qqp(resid(mod2), "norm")

Anova(mod1, type = 'III')
qqp(subset_low$log_k)
qqp(resid(mod1), "norm")

#Kenward-Rogers
anova(mod2,mod1)
KRmodcomp(mod2,mod1)

mod2 <- lmer(log_k ~ 1 + Protocol*Time_DD + MSAS_18_sr + MSAS_18_cd + MSAS_18_so + MSAS_18_m + (1|Id), data = subset_high,
             REML = FALSE)
mod1 <- lmer(log_k ~ 1 + Protocol*Time_DD + MSAS_18_sr + MSAS_18_cd + MSAS_18_so + MSAS_18_m + Protocol:(MSAS_18_sr + MSAS_18_cd + MSAS_18_so + MSAS_18_m) + (1 | Id), data = subset_high,
             REML = FALSE)

Anova(mod2, type = 'III')
qqp(subset_low$log_k)
qqp(resid(mod2), "norm")

Anova(mod1, type = 'III')
qqp(subset_low$log_k)
qqp(resid(mod1), "norm")

#Kenward-Rogers
anova(mod2,mod1)
KRmodcomp(mod2,mod1)


### DD Boxplot Raw Data ###
data_Controlli_TS_LME_Ext <- read_excel("Desktop/TS OCD & NIBS/NIBS & Inhibitory Control - UniSR/Preliminar Analysis Plots & Datasets/Dataset_EXP_tDCS_PhD_LME_Ext.xlsx", col_names = TRUE)
summary(data_Controlli_TS_LME_Ext)
data_Controlli_TS_LME_Ext <- as.data.frame(data_Controlli_TS_LME_Ext)

#Data Preparation
if("k_Post" %in% colnames(data_Controlli_TS_LME_Ext)) {
  data_Controlli_TS_LME_Ext$log_k_Post <- log(data_Controlli_TS_LME_Ext$k_Post)
} else {
  stop("La colonna 'k_Post' non esiste nel dataframe.")
}
summary(data_Controlli_TS_LME_Ext)

if("k_Pre" %in% colnames(data_Controlli_TS_LME_Ext)) {
  data_Controlli_TS_LME_Ext$log_k_Pre <- log(data_Controlli_TS_LME_Ext$k_Pre)
} else {
  stop("La colonna 'k_Pre' non esiste nel dataframe.")
}
summary(data_Controlli_TS_LME_Ext)

data_Controlli_TS_LME_Ext <- data_Controlli_TS_LME_Ext %>% 
  mutate(Go_noGo_Mean_RT_GoTrials = (Go_noGo_Mean_RT_GoTrials_P + Go_noGo_Mean_RT_GoTrials_R)/2,
         Go_noGo_Mean_RT_noGoTrials = (Go_noGo_Mean_RT_noGOTrials_P + Go_noGo_Mean_RT_noGOTrials_R)/2,
         Go_noGo_Mean_RT = (Go_noGo_Mean_RT_GoTrials_P + Go_noGo_Mean_RT_GoTrials_R + Go_noGo_Mean_RT_noGOTrials_P + Go_noGo_Mean_RT_noGOTrials_R)/4,
         Go_no_Go_Mean_FA = (FA_Block_1 + FA_Block_2)/2,
         Go_noGo_Mean_CR = (CR_Block_1 + CR_Block_2)/2,
         Go_noGO_Mean_d_prime = (d_prime_B_1 + d_prime_B_2)/2)
summary(data_Controlli_TS_LME_Ext)

subset_low <- data_Controlli_TS_LME_Ext %>% 
  filter(Reward == "Low")

subset_high <- data_Controlli_TS_LME_Ext %>% 
  filter(Reward == "High")

subset_low$log_k_diff <- subset_low$log_k_Post - subset_low$log_k_Pre
summary(subset_low)
subset_low$Protocol <- as.factor(subset_low$Protocol)

subset_high$log_k_diff <- subset_high$log_k_Post - subset_high$log_k_Pre
summary(subset_high)
subset_high$Protocol <- as.factor(subset_high$Protocol)

P5_FF_Gr <- ggplot(subset_low, aes(x = Protocol, y = log_k_diff, color = Protocol)) +
  geom_point(aes(x = Protocol, y = log_k_diff, color = Protocol), position = position_jitter(width = 0.2)) +
  geom_boxplot(aes(x = Protocol, y = log_k_diff), color = "black", fill = NA, width = 0.2, outlier.shape = NA, position = position_dodge(width = 0.5), fatten = NULL) +
  stat_summary(aes(x = Protocol, y = log_k_diff), fun = mean, geom = "crossbar", width = 0.2, color = "black", linetype = "solid", linewidth = 0.3, position = position_dodge(width = 0.5)) +
  labs(x = NULL, y = "Log(k) Post - Log(k) Pre") +
  theme_classic() +
  scale_y_continuous(limits = c(-3.5, 3.5), breaks = seq(-4.0, 3.5, by = 1.0), labels = scales::number_format(accuracy = 0.1)) +
  showSignificance( c(1,2), 3.0, -0.20) +
  annotate("text", x = mean(c(1,2)), y = 3.2, label = "*", size = 15) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),  
        legend.title = element_text(size = 18),
        legend.position = "top")

P6_FF_Gr <- ggplot(subset_high, aes(x = Protocol, y = log_k_diff, color = Protocol)) +
  geom_point(aes(x = Protocol, y = log_k_diff, color = Protocol), position = position_jitter(width = 0.2)) +
  geom_boxplot(aes(x = Protocol, y = log_k_diff), color = "black", fill = NA, width = 0.2, outlier.shape = NA, position = position_dodge(width = 0.5), fatten = NULL) +
  stat_summary(aes(x = Protocol, y = log_k_diff), fun = mean, geom = "crossbar", width = 0.2, color = "black", linetype = "solid", linewidth = 0.3, position = position_dodge(width = 0.5)) +
  labs(x = NULL, y = "Log(k) Post - Log(k) Pre") +
  theme_classic() +
  scale_y_continuous(limits = c(-3.5, 3.5), breaks = seq(-4.0, 3.5, by = 1.0), labels = scales::number_format(accuracy = 0.1)) +
  showSignificance( c(1,2), 3.0, -0.20) +
  annotate("text", x = mean(c(1,2)), y = 3.3, label = "ns", size = 6) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),  
        legend.title = element_text(size = 18),
        legend.position = "top")

#Creo Figura 3 (a, b, c, d)
library("cowplot")
P5_FF_Gr <- P5_FF_Gr + theme(legend.position = "none")
P6_FF_Gr <- P6_FF_Gr + theme(legend.position = "none")
ggdraw() +
  draw_plot(PDD1_F, x = 0, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(PDD2_F, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(P6_FF_Gr, x = 0, y = 0, width = 0.5, height = 0.5) +
  draw_plot(P5_FF_Gr, x = 0.5, y = 0, width = 0.5, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 30,
                  x = c(0, 0.5, 0, 0.5), 
                  y = c(1, 1, 0.5, 0.5))

### DA QUA IN POI NON CONSIDERARE ###
#Barplot Interaction
mean_data <- subset_low %>%
  group_by(Time_DD, Protocol) %>%
  summarise(
    mean_log_k = mean(Fitted_log_k),
    sem_log_k = sd(Fitted_log_k) / sqrt(n()),
    sd_log_k = sd(Fitted_log_k)) 

mean_data$Time_DD <- as.factor(mean_data$Time_DD)
F5_F <- ggplot(mean_data, aes(x = Time_DD, y = mean_log_k, fill = Protocol)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", size = 0.7) +
  geom_errorbar(aes(ymin = mean_log_k - sem_log_k, ymax = mean_log_k + sem_log_k),
                width = 0.2, size = 0.8, position = position_dodge(width = 0.8)) +
  geom_point(data = subset_low, aes(x = Time_DD, y = Fitted_log_k, group = Protocol), 
             position = position_dodge(width = 1.3), size = 2, shape = 21, color = "black", fill = "white") +
  labs(x = NULL,
       y = "Log(k)") +
  theme_classic() +
  scale_x_discrete(labels = c("Pre", "Post")) +
  scale_y_continuous() +
  #annotate("text", x = "1", y = 1, label = "'", size = 10) + 
  #annotate("text", x = "2", y = 1, label = "*", size = 10) +
  #annotate("text", x = "3", y = 1, label = "*", size = 10) +
  #annotate("text", x = "4", y = 1, label = "*", size = 10) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),  
        legend.title = element_text(size = 18),
        legend.position = "top")

LME_DD_High_long <- lmer(log_k ~ 1 + Protocol*Time_DD + Protocol*MCQ_30_tot + Protocol*MSAS_18_tot + (1|Id), data = subset_high,
                         REML = FALSE)
Anova(LME_DD_High_long, type = 'III')
qqp(subset_high$log_k)
qqp(resid(LME_DD_High_long), "norm")

LME_DD_High_short <- lmer(log_k ~ 1 + Protocol*Time_DD + (1|Id), data = subset_high,
                          REML = FALSE)
Anova(LME_DD_High_short, type = 'III')
qqp(subset_high$log_k)
qqp(resid(LME_DD_High_short), "norm")
#Kenward-Rogers
anova(LME_DD_High_short,LME_DD_High_long)
KRmodcomp(LME_DD_High_short,LME_DD_High_long)

#CONCL. Anche in questo caso meglio il modello ristretto con salmente l'interazione Protocol*Time (Sign.: Time)


#####
#### Loading Data For DD Analisi ####
data_Controlli_TS_LME_Ext <- read_excel("Desktop/TS OCD & NIBS/NIBS & Inhibitory Control - UniSR/Preliminar Analysis Plots & Datasets/Dataset_EXP_tDCS_PhD_LME_Ext.xlsx", col_names = TRUE)
summary(data_Controlli_TS_LME_Ext)
data_Controlli_TS_LME_Ext <- as.data.frame(data_Controlli_TS_LME_Ext)

#Data Preparation
if("k_Post" %in% colnames(data_Controlli_TS_LME_Ext)) {
  data_Controlli_TS_LME_Ext$log_k_Post <- log(data_Controlli_TS_LME_Ext$k_Post)
} else {
  stop("La colonna 'k_Post' non esiste nel dataframe.")
}
summary(data_Controlli_TS_LME_Ext)

if("k_Pre" %in% colnames(data_Controlli_TS_LME_Ext)) {
  data_Controlli_TS_LME_Ext$log_k_Pre <- log(data_Controlli_TS_LME_Ext$k_Pre)
} else {
  stop("La colonna 'k_Pre' non esiste nel dataframe.")
}
summary(data_Controlli_TS_LME_Ext)

data_Controlli_TS_LME_Ext <- data_Controlli_TS_LME_Ext %>% 
  mutate(Go_noGo_Mean_RT_GoTrials = (Go_noGo_Mean_RT_GoTrials_P + Go_noGo_Mean_RT_GoTrials_R)/2,
         Go_noGo_Mean_RT_noGoTrials = (Go_noGo_Mean_RT_noGOTrials_P + Go_noGo_Mean_RT_noGOTrials_R)/2,
         Go_noGo_Mean_RT = (Go_noGo_Mean_RT_GoTrials_P + Go_noGo_Mean_RT_GoTrials_R + Go_noGo_Mean_RT_noGOTrials_P + Go_noGo_Mean_RT_noGOTrials_R)/4,
         Go_no_Go_Mean_FA = (FA_Block_1 + FA_Block_2)/2,
         Go_noGo_Mean_CR = (CR_Block_1 + CR_Block_2)/2,
         Go_noGO_Mean_d_prime = (d_prime_B_1 + d_prime_B_2)/2)
summary(data_Controlli_TS_LME_Ext)
data_Controlli_TS_LME_Ext$log_k_diff <- data_Controlli_TS_LME_Ext$log_k_Post - data_Controlli_TS_LME_Ext$log_k_Pre
data_low_reward <- data_Controlli_TS_LME_Ext %>%
  dplyr::filter(Reward == "Low")

P5_F <- ggplot(data_low_reward, aes(x = Protocol, y = log_k_diff, color = Protocol)) +
  geom_point(aes(x = Protocol, y = log_k_diff, color = Protocol), position = position_jitter(width = 0.2)) +
  geom_boxplot(color = "black", fill = NA, width = 0.2, outlier.shape = NA, position = position_dodge(width = 0.5), fatten = NULL) +
  stat_summary(aes(x = Protocol, y = log_k_diff), fun = mean, geom = "crossbar", width = 0.2, color = "black", linetype = "solid", linewidth = 0.3, position = position_dodge(width = 0.5)) +
  labs(x = NULL, y = "Log(k) Post - Log(k) Pre") +
  theme_classic() +
  scale_y_continuous(limits = c(-4.0, 3.2), breaks = seq(-4.0, 3.2, by = 1.0), labels = scales::number_format(accuracy = 0.1)) +
  showSignificance( c(1,2), 3.0, -0.2) +
  annotate("text", x = mean(c(1,2)), y = 3.1, label = "*", size = 15) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),  
        legend.title = element_text(size = 18),
        legend.position = "top")

data_high_reward <- data_Controlli_TS_LME_Ext %>%
  dplyr::filter(Reward == "High")

P6_F <- ggplot(data_high_reward, aes(x = Protocol, y = log_k_diff, color = Protocol)) +
  geom_point(aes(x = Protocol, y = log_k_diff, color = Protocol), position = position_jitter(width = 0.2)) +
  geom_boxplot(color = "black", fill = NA, width = 0.2, outlier.shape = NA, position = position_dodge(width = 0.5), fatten = NULL) +
  stat_summary(aes(x = Protocol, y = log_k_diff), fun = mean, geom = "crossbar", width = 0.2, color = "black", linetype = "solid", linewidth = 0.3, position = position_dodge(width = 0.5)) +
  labs(x = NULL, y = "Log(k) Post - Log(k) Pre") +
  theme_classic() +
  scale_y_continuous(limits = c(-4.0, 3.3), breaks = seq(-4.0, 3.3, by = 1.0), labels = scales::number_format(accuracy = 0.1)) +
  showSignificance( c(1,2), 3.0, -0.2) +
  annotate("text", x = mean(c(1,2)), y = 3.3, label = "ns", size = 5) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),  
        legend.title = element_text(size = 18),
        legend.position = "top")
