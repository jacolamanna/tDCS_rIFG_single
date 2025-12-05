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

LME_log_k_Pooled <- lmer(log_k_Post ~ 1 + Protocol*Reward + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + (1|Id) + (1|Protocol_ord), data = data_Controlli_TS_LME_Ext,
                         REML = FALSE)
Anova(LME_log_k_Pooled)
LME_AUC_Pooled <- lmer(AUC_Post ~ 1 + Protocol*Reward + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + (1|Id) + (1|Protocol_ord), data = data_Controlli_TS_LME_Ext,
                         REML = FALSE)
Anova(LME_AUC_Pooled)
#No tDCS effect (solo Reward), provo a fittare LME su High e Low separatamente

#Loading Datset
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

#LME tDCS ~ AUC High/Low
LME_AUC_High <- lmer(AUC_High_Post ~ 1 + Protocol + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + + MCQ_30_nc + MCQ_30_csc + (1|Id) + (1|Protocol_ord), data = data_Controlli_TS_LME,
                    REML = FALSE)
Anova(LME_AUC_High)
#Significant predictors: NONE (X2(1): 0.0056), p=0.9405)
qqp(data_Controlli_TS_LME$AUC_High_Post)
qqp(resid(LME_AUC_High), "norm")

LME_AUC_Low <- lmer(AUC_Low_Post ~ 1 + Protocol + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + + MCQ_30_nc + MCQ_30_csc + (1|Id) + (1|Protocol_ord), data = data_Controlli_TS_LME,
                     REML = FALSE)
Anova(LME_AUC_Low)
#Significant predictors: NONE (X2(1): 0.6555), p=0.4181)
qqp(data_Controlli_TS_LME$AUC_Low_Post)
qqp(resid(LME_AUC_Low), "norm")
#Non trovo effetto della tDCS su AUC, provo con i lok(k)

#LME tDCS ~ log(k) High/Low
LME_log_k_High <- lmer(log_k_High_Post ~ 1 + Protocol + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + (1|Id) + (1|Protocol_ord), data = data_Controlli_TS_LME,
                   REML = FALSE)
Anova(LME_log_k_High)
#Significant predictors: NONE (tDCS: X(1): 0.0826, p=0.7739)
qqp(data_Controlli_TS_LME$log_k_High_Post)
qqp(resid(LME_log_k_High), "norm")

LME_log_k_High_short <- lmer(log_k_High_Post ~ 1 + Protocol + (1|Id) + (1|Protocol_ord), data = data_Controlli_TS_LME,
                       REML = FALSE)
Anova(LME_log_k_High_short)
anova(LME_log_k_High_short)
#Significant predictors: NONE (X(1): 0.6757, p=0.4111)
qqp(data_Controlli_TS_LME$log_k_High_Post)
qqp(resid(LME_log_k_High_short), "norm")
#Kenward-Roger for LME comparison
anova(LME_log_k_High,LME_log_k_High_short)
KRmodcomp(LME_log_k_High,LME_log_k_High_short)
#No difference: X2(5): 2.9245, p = 0.7116
#Add Fitted data to dataset
data_Controlli_TS_LME$Fitted_log_k_High <- fitted(LME_log_k_High_short)
summary(data_Controlli_TS_LME)

LME_log_k_Low <- lmer(log_k_Low_Post ~ 1 + Protocol + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + (1|Id) + (1|Protocol_ord), data = data_Controlli_TS_LME,
                       REML = FALSE)
Anova(LME_log_k_Low)
#Significant predictors: NONE (tDCS: X2(1): 3.1585, p=0.07553)
qqp(data_Controlli_TS_LME$log_k_Low_Post)
qqp(resid(LME_log_k_Low), "norm")

LME_log_k_Low_short <- lmer(log_k_Low_Post ~ 1 + Protocol + (1|Id) + (1|Protocol_ord), data = data_Controlli_TS_LME,
                      REML = FALSE)
Anova(LME_log_k_Low_short)
anova(LME_log_k_Low_short)
#Significant predictors: tDCS (X2(1): 4.3453, p=0.03711)
qqp(data_Controlli_TS_LME$log_k_Low_Post)
qqp(resid(LME_log_k_Low_short), "norm")
#Kenward-Roger for LME comparison
anova(LME_log_k_Low,LME_log_k_Low_short)
KRmodcomp(LME_log_k_Low,LME_log_k_Low_short)
#No difference: X2(5): 1.0907, p = 0.9549
#Add Fitted data to dataset
data_Controlli_TS_LME$Fitted_log_k_Low <- fitted(LME_log_k_Low_short)
summary(data_Controlli_TS_LME)

#Plot log(k) Low e log(k) High
PDD3_F <- ggplot(data_Controlli_TS_LME, aes(x = Protocol, y = Fitted_log_k_High, color = Protocol)) +
  geom_point(position = position_jitter(width = 0.2)) +
  geom_boxplot(color = "black", fill = NA, width = 0.2, outlier.shape = NA, position = position_dodge(width = 0.5), fatten = NULL) +
  stat_summary(aes(x = Protocol, y = Fitted_log_k_High), fun = mean, geom = "crossbar", width = 0.2, color = "black", linetype = "solid", linewidth = 0.3, position = position_dodge(width = 0.5)) +
  labs(x = NULL, y = "Log(k)") +
  theme_classic() +
  scale_y_continuous(limits = c(-9.5, 3.1), breaks = seq(-9.5, 3.1, by = 2)) +
  showSignificance( c(1,2), 2.5, -0.4) +
  annotate("text", x = mean(c(1,2)), y = 3.1, label = "ns", size = 6) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),  
        legend.title = element_text(size = 18),
        legend.position = "top")

PDD4_F <- ggplot(data_Controlli_TS_LME, aes(x = Protocol, y = Fitted_log_k_Low, color = Protocol)) +
  geom_point(position = position_jitter(width = 0.2)) +
  geom_boxplot(color = "black", fill = NA, width = 0.2, outlier.shape = NA, position = position_dodge(width = 0.5), fatten = NULL) +
  stat_summary(aes(x = Protocol, y = Fitted_log_k_Low), fun = mean, geom = "crossbar", width = 0.2, color = "black", linetype = "solid", linewidth = 0.3, position = position_dodge(width = 0.5)) +
  labs(x = NULL, y = "Log(k)") +
  theme_classic() +
  scale_y_continuous(limits = c(-9.5, 2.0), breaks = seq(-9.5, 2.0, by = 2)) +
  showSignificance( c(1,2), 1.3, -0.4) +
  annotate("text", x = mean(c(1,2)), y = 1.5, label = "*", size = 15) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),  
        legend.title = element_text(size = 18),
        legend.position = "top")

#Plot IP Low
data_IP_Controlli_TS <- read_excel("Desktop/TS OCD & NIBS/NIBS & Inhibitory Control - UniSR/Preliminar Analysis Plots & Datasets/IPs_dataset_LOW.xlsx", col_names = TRUE)
summary(data_IP_Controlli_TS)
data_IP_Controlli_TS <- as.data.frame(data_IP_Controlli_TS)

data_IP_Controlli_TS_lomg <- data_IP_Controlli_TS %>%
  pivot_longer(cols = c(`0`, `1`, `6`, `12`, `24`, `60`, `120`), 
               names_to = "Delay", 
               values_to = "Score")
data_IP_Controlli_TS_lomg$Delay <- as.numeric(data_IP_Controlli_TS_lomg$Delay)
data_IP_Controlli_TS_lomg_summary <- data_IP_Controlli_TS_lomg %>%
  group_by(Protocol, Delay) %>%
  summarise(Mean = mean(Score),
            SD = sd(Score),
            N = n(),
            SE = SD / sqrt(N))

PDD2_F <- ggplot(data_IP_Controlli_TS_lomg_summary, aes(x = Delay, y = Mean, color = Protocol)) +
  geom_point() +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 1) +
  labs(x = "Delay (Months)",
       y = "Indifference Points (Euros)") +
  scale_x_continuous(limits = c(0, 123), breaks = seq(0, 120, by = 20)) +
  scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, by = 100)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),  
        legend.title = element_text(size = 18),
        legend.position = "top")

Real <- subset(data_IP_Controlli_TS,  Protocol == "Real", `120`,
                 drop = TRUE)
Sham <- subset(data_IP_Controlli_TS,  Protocol == "Sham", `120`,
                drop = TRUE)
res <- wilcox.test(Real, Sham, paired = TRUE)
res

#Plot IP High
data_IP_High_Controlli_TS <- read_excel("Desktop/TS OCD & NIBS/NIBS & Inhibitory Control - UniSR/Preliminar Analysis Plots & Datasets/IPs_dataset_HIGH.xlsx", col_names = TRUE)
summary(data_IP_High_Controlli_TS)
data_IP_High_Controlli_TS <- as.data.frame(data_IP_High_Controlli_TS)

data_IP_High_Controlli_TS_long <- data_IP_High_Controlli_TS %>%
  pivot_longer(cols = c(`0`, `1`, `6`, `12`, `24`, `60`, `120`), 
               names_to = "Delay", 
               values_to = "Score")
data_IP_High_Controlli_TS_long$Delay <- as.numeric(data_IP_High_Controlli_TS_long$Delay)
data_IP_High_Controlli_TS_long_summary <- data_IP_High_Controlli_TS_long %>%
  group_by(Protocol, Delay) %>%
  summarise(Mean = mean(Score),
            SD = sd(Score),
            N = n(),
            SE = SD / sqrt(N))

PDD1_F <- ggplot(data_IP_High_Controlli_TS_long_summary, aes(x = Delay, y = Mean, color = Protocol)) +
  geom_point() +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 1) +
  labs(x = "Delay (Months)",
       y = "Indifference Points (Euros)") +
  scale_x_continuous(limits = c(0, 123), breaks = seq(0, 120, by = 20)) +
  scale_y_continuous(limits = c(0, 10000), breaks = seq(0, 10000, by = 2000)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),  
        legend.title = element_text(size = 18),
        legend.position = "top")

Real <- subset(data_IP_High_Controlli_TS,  Protocol == "Real", `120`,
               drop = TRUE)
Sham <- subset(data_IP_High_Controlli_TS,  Protocol == "Sham", `120`,
               drop = TRUE)
res <- wilcox.test(Real, Sham, paired = TRUE)
res

#Creo Figura 2 (a, b, c, d)
library("cowplot")
PDD3_F <- PDD3_F + theme(legend.position = "none")
PDD4_F <- PDD4_F + theme(legend.position = "none")
ggdraw() +
  draw_plot(PDD1_F, x = 0, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(PDD2_F, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(PDD3_F, x = 0, y = 0, width = 0.5, height = 0.5) +
  draw_plot(PDD4_F, x = 0.5, y = 0, width = 0.5, height = 0.5) +
  draw_plot_label(label = c("a", "b", "c", "d"), size = 30,
                  x = c(0, 0.5, 0, 0.5), 
                  y = c(1, 1, 0.5, 0.5))
