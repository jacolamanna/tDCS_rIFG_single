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

### LME - Inhibitory Control - SSRT ####
LME_InC_SSRT <- lmer(SSRT ~ 1 + Protocol*MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + (1|Id) + (1|Protocol_ord), data = data_Controlli_TS_LME,
                     REML = FALSE)
Anova(LME_InC_SSRT)
anova(LME_InC_SSRT)
#Significant predictors: Protocol (X2(1): 6.0310, p=0.01406), MCQ_30_POS (X2(1): 6.0255, p=0.01410), MCQ_30_NC (X2(1): 4.0488, p=0.04420)
qqp(data_Controlli_TS_LME$SSRT)
qqp(resid(LME_InC_SSRT), "norm")

LME_InC_SSRT_l <- lmer(SSRT ~ 1 + Protocol*MCQ_30_pos + Protocol*MCQ_30_neg + Protocol*MCQ_30_cc + Protocol*MCQ_30_nc + Protocol*MCQ_30_csc + Protocol*AUC_Low_Pre + Protocol*AUC_High_Pre + (1|Id) + (1|Protocol_ord), data = data_Controlli_TS_LME,
                     REML = FALSE)
Anova(LME_InC_SSRT_l)
anova(LME_InC_SSRT_l)
anova(LME_InC_SSRT_l,LME_InC_SSRT_sign)
KRmodcomp(LME_InC_SSRT_l,LME_InC_SSRT_sign)

LME_InC_SSRT_short <- lmer(SSRT ~ 1 + Protocol + (1|Id) + (1|Protocol_ord), data = data_Controlli_TS_LME,
                           REML = FALSE)
Anova(LME_InC_SSRT_short)
anova(LME_InC_SSRT_short)
#Significant predictors: Protocol (X2(1): 9.794, p=0.001751), MCQ30_POS, MCQ_30_NC
qqp(data_Controlli_TS_LME$SSRT)
qqp(resid(LME_InC_SSRT_short), "norm")
#Kenward-Roger for LME comparison
anova(LME_InC_SSRT,LME_InC_SSRT_short)
KRmodcomp(LME_InC_SSRT,LME_InC_SSRT_short)
#No difference: X2(5): 8.5585, p = 0.1280

LME_InC_SSRT_sign <- lmer(SSRT ~ 1 + Protocol*MCQ_30_pos + (1|Id) + (1|Protocol_ord), data = data_Controlli_TS_LME,
                           REML = FALSE)
Anova(LME_InC_SSRT_sign)
anova(LME_InC_SSRT_sign)
#Significant predictors: Protocol (X2(1): 8.4654, p=0.00362), MCQ_30_POS (X2(1): 3.9342, p=0.04731), no interaction
qqp(data_Controlli_TS_LME$SSRT)
qqp(resid(LME_InC_SSRT_sign), "norm")
#Kenward-Roger for LME comparison
anova(LME_InC_SSRT,LME_InC_SSRT_sign)
KRmodcomp(LME_InC_SSRT,LME_InC_SSRT_sign)
#No difference: X2(3): 4.8137, p = 0.1860

#Add Fitted data to dataset
data_Controlli_TS_LME$Fitted_SSRT <- fitted(LME_InC_SSRT)
summary(data_Controlli_TS_LME)

#Boxplot with Means tDCS ~ SSRT
P1_F <- ggplot(data_Controlli_TS_LME, aes(x = Protocol, y = Fitted_SSRT, color = Protocol)) +
  geom_point(position = position_jitter(width = 0.2)) +
  geom_boxplot(color = "black", fill = NA, width = 0.2, outlier.shape = NA, position = position_dodge(width = 0.5), fatten = NULL) +
  stat_summary(aes(x = Protocol, y = Fitted_SSRT), fun = mean, geom = "crossbar", width = 0.2, color = "black", linetype = "solid", linewidth = 0.3, position = position_dodge(width = 0.5)) +
  labs(x = NULL, y = "Stop-Signal Reaction Time")  +
  theme_classic() +
  scale_y_continuous(limits = c(190, 272), breaks = seq(190, 272, by = 15)) +
  showSignificance( c(1,2), 270, -2.0) +
  annotate("text", x = mean(c(1,2)), y = 271, label = "**", size = 15) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),  
        legend.title = element_text(size = 18),
        legend.position = "top")

P1_Fb <- ggplot(data_Controlli_TS_LME) + 
  geom_point(aes(x=MCQ_30_pos, y=Fitted_SSRT, color = Protocol))+
  geom_smooth(aes(x=MCQ_30_pos, y=Fitted_SSRT), method=lm, se=FALSE, linetype="dashed", size=0.2, color = "black") +
  theme_classic() +
  labs(x = "Positive Beliefs about Thinking", y = "Stop-Signal Reaction Time") + 
  scale_y_continuous(limits = c(190, 272), breaks = seq(190, 272, by = 15)) +
  scale_x_continuous(limits = c(5, 25), breaks = seq(5, 25, by = 5)) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),  
        legend.title = element_text(size = 18),
        legend.position = "top") 

### LME - Inhibitory Control - GoTrials_RT e FA ####
LME_InC_GoNoGo <- lmer(Go_noGo_Mean_RT_GoTrials ~ 1 + Protocol*MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + (1|Id) + (1|Protocol_ord), data = data_Controlli_TS_LME,
                       REML = FALSE)
Anova(LME_InC_GoNoGo)
anova(LME_InC_GoNoGo)
#Significant predictors: MCQ_30_NC (X2(1): 4.6204, p=0.03159)
qqp(data_Controlli_TS_LME$Go_noGo_Mean_RT_GoTrials)
qqp(resid(LME_InC_GoNoGo), "norm")

LME_InC_GoNoGo_l <- lmer(Go_noGo_Mean_RT_GoTrials ~ 1 + Protocol*MCQ_30_pos + Protocol*MCQ_30_neg + Protocol*MCQ_30_cc + Protocol*MCQ_30_nc + Protocol*MCQ_30_csc + Protocol*AUC_Low_Pre + Protocol*AUC_High_Pre + (1|Id) + (1|Protocol_ord), data = data_Controlli_TS_LME,
                       REML = FALSE)
Anova(LME_InC_GoNoGo_l)
anova(LME_InC_GoNoGo_l)
anova(LME_InC_GoNoGo_l,LME_InC_GoNoGo_short)
KRmodcomp(LME_InC_SSRT_l,LME_InC_GoNoGo_short)

LME_InC_GoNoGo_short <- lmer(Go_noGo_Mean_RT_GoTrials ~ 1 + Protocol*MCQ_30_pos + (1|Id) + (1|Protocol_ord), data = data_Controlli_TS_LME,
                           REML = FALSE)
Anova(LME_InC_GoNoGo_short)
anova(LME_InC_GoNoGo_short)
#Significant predictors: NONE (X2(1): 0.2319; p=0.6301)
qqp(data_Controlli_TS_LME$Go_noGo_Mean_RT_GoTrials)
qqp(resid(LME_InC_GoNoGo_short), "norm")
#Kenward-Roger for LME comparison
anova(LME_InC_GoNoGo,LME_InC_GoNoGo_short)
KRmodcomp(LME_InC_GoNoGo,LME_InC_GoNoGo_short)
#No difference: X2(5): 4.9049, p = 0.4276

#Add Fitted data to dataset
data_Controlli_TS_LME$Fitted_Go_noGo_Mean_RT_GoTrials <- fitted(LME_InC_GoNoGo)
summary(data_Controlli_TS_LME)

#Boxplot with Means tDCS ~ RT_GoTrials
P2_F <- ggplot(data_Controlli_TS_LME, aes(x = Protocol, y = Fitted_Go_noGo_Mean_RT_GoTrials, color = Protocol)) +
  geom_point(position = position_jitter(width = 0.2)) +
  geom_boxplot(color = "black", fill = NA, width = 0.2, outlier.shape = NA, position = position_dodge(width = 0.5), fatten = NULL) +
  stat_summary(aes(x = Protocol, y = Fitted_Go_noGo_Mean_RT_GoTrials), fun = mean, geom = "crossbar", width = 0.2, color = "black", linetype = "solid", linewidth = 0.3, position = position_dodge(width = 0.5)) +
  labs(x = NULL, y = "Reaction Time Go Trials") +
  theme_classic() +
  scale_y_continuous(limits = c(365, 450), breaks = seq(365, 450, by = 15)) +
  showSignificance( c(1,2), 446, -2.0) +
  annotate("text", x = mean(c(1,2)), y = 450, label = "ns", size = 6) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),  
        legend.title = element_text(size = 18),
        legend.position = "top")

#FA
LME_InC_GoNoGo_FA <- lmer(Go_no_Go_Mean_FA ~ 1 + Protocol*MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + (1|Id) + (1|Protocol_ord), data = data_Controlli_TS_LME,
                       REML = FALSE)
Anova(LME_InC_GoNoGo_FA)
anova(LME_InC_GoNoGo_FA)
#Significant predictors: NONE 
qqp(data_Controlli_TS_LME$Go_no_Go_Mean_FA)
qqp(resid(LME_InC_GoNoGo_FA), "norm")

LME_InC_GoNoGo_FA_l <- lmer(Go_no_Go_Mean_FA ~ 1 + Protocol*MCQ_30_pos + Protocol*MCQ_30_neg + Protocol*MCQ_30_cc + Protocol*MCQ_30_nc + Protocol*MCQ_30_csc + Protocol*AUC_Low_Pre + Protocol*AUC_High_Pre + (1|Id) + (1|Protocol_ord), data = data_Controlli_TS_LME,
                         REML = FALSE)
Anova(LME_InC_GoNoGo_FA_l)
anova(LME_InC_GoNoGo_FA_l)
anova(LME_InC_GoNoGo_FA_l,LME_InC_GoNoGo_FA_short)
KRmodcomp(LME_InC_GoNoGo_FA_l,LME_InC_GoNoGo_FA_short)

LME_InC_GoNoGo_FA_short <- lmer(Go_no_Go_Mean_FA ~ 1 + Protocol*MCQ_30_pos + (1|Id) + (1|Protocol_ord), data = data_Controlli_TS_LME,
                          REML = FALSE)
Anova(LME_InC_GoNoGo_FA_short)
anova(LME_InC_GoNoGo_FA_short)
#Significant predictors: NONE (X(2): 1.6189, p=0.2033)
qqp(data_Controlli_TS_LME$Go_no_Go_Mean_FA)
qqp(resid(LME_InC_GoNoGo_FA_short), "norm")
#Kenward-Roger for LME comparison
anova(LME_InC_GoNoGo_FA,LME_InC_GoNoGo_FA_short)
KRmodcomp(LME_InC_GoNoGo_FA,LME_InC_GoNoGo_FA_short)
#No difference: X2(5): 2.6688, p = 0.7509

LME_InC_GoNoGo_FA_sign <- lmer(Go_no_Go_Mean_FA ~ 1 + Protocol*MCQ_30_pos + (1|Id) + (1|Protocol_ord), data = data_Controlli_TS_LME,
                          REML = FALSE)
Anova(LME_InC_GoNoGo_FA_sign)
anova(LME_InC_GoNoGo_FA_sign)
#Significant predictors: NONE (X(2): 1.7759, p=0.1827)
qqp(data_Controlli_TS_LME$Go_no_Go_Mean_FA)
qqp(resid(LME_InC_GoNoGo_FA_sign), "norm")
#Kenward-Roger for LME comparison
anova(LME_InC_GoNoGo_FA,LME_InC_GoNoGo_FA_sign)
KRmodcomp(LME_InC_GoNoGo_FA,LME_InC_GoNoGo_FA_sign)
#No difference: X2(3): 0.9168, p = 0.8214

#Add Fitted data to dataset
data_Controlli_TS_LME$Fitted_Go_no_Go_Mean_FA <- fitted(LME_InC_GoNoGo_FA_short)
summary(data_Controlli_TS_LME)

P3_F <- ggplot(data_Controlli_TS_LME, aes(x = Protocol, y = Fitted_Go_no_Go_Mean_FA, color = Protocol)) +
  geom_point(position = position_jitter(width = 0.2)) +
  geom_boxplot(color = "black", fill = NA, width = 0.2, outlier.shape = NA, position = position_dodge(width = 0.5), fatten = NULL) +
  stat_summary(aes(x = Protocol, y = Fitted_Go_no_Go_Mean_FA), fun = mean, geom = "crossbar", width = 0.2, color = "black", linetype = "solid", linewidth = 0.3, position = position_dodge(width = 0.5)) +
  labs(x = NULL, y = "False Alarm Rate") +
  theme_classic() +
  scale_y_continuous(limits = c(0.0, 0.17), breaks = seq(0.0, 0.17, by = 0.03)) +
  showSignificance( c(1,2), 0.162, -0.005) +
  annotate("text", x = mean(c(1,2)), y = 0.17, label = "ns", size = 6) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),  
        legend.title = element_text(size = 18),
        legend.position = "top")

#Compongo Figura 1 (a,b,c)
P1_F <- P1_F + theme(legend.position = "none", aspect.ratio = 1)
P1_Fb <- P1_Fb + theme(aspect.ratio = 1)
P3_F <- P3_F + theme(legend.position = "none", aspect.ratio = 1)
library("cowplot")
ggdraw() +
  draw_plot(P1_F, x = 0.0, y = 0.15, width = 1/3, height = 0.85) +
  draw_plot(P1_Fb, x = 1/3, y = 0.15, width = 1/3, height = 0.85) +
  draw_plot(P3_F, x = 2/3, y = 0.15, width = 1/3, height = 0.85) +
  draw_plot_label(label = c("a", "b", "c"), size = 30,
                  x = c(0, 1/3, 2/3), 
                  y = c(1.0, 1.0, 1.0))

#Compongo Figura 1 (a,b)
ggdraw() +
  draw_plot(P1_F, x = 0.0, y = 0.0, width = 0.5, height = 1.0) +
  draw_plot(P3_F, x = 0.5, y = 0.0, width = 0.5, height = 1.0) +
  draw_plot_label(label = c("a", "b"), size = 30,
                  x = c(0, 0.5), 
                  y = c(1.0, 1.0))

#Compongo Figura 1 (a,b,c,d)
P1_F <- P1_F + theme(legend.position = "none")
P2_F <- P2_F + theme(legend.position = "none")
P3_F <- P3_F + theme(legend.position = "none")
ggdraw() +
  draw_plot(P1_F, x = 0, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(P1_Fb, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(P2_F, x = 0, y = 0, width = 0.5, height = 0.5) +
  draw_plot(P3_F, x = 0.5, y = 0, width = 0.5, height = 0.5) +
  draw_plot_label(label = c("a", "b", "c", "d"), size = 30,
                  x = c(0, 0.5, 0, 0.5), 
                  y = c(1, 1, 0.5, 0.5))