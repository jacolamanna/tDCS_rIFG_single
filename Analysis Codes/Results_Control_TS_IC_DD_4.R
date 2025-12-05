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

##### LME SSRT #####
LME_SSRT_long <- lmer(SSRT ~ 1 + Protocol*MCQ_30_pos + Protocol*MCQ_30_neg + Protocol*MCQ_30_cc + Protocol*MCQ_30_nc + Protocol*MCQ_30_csc + Protocol*log_k_Pre + (1|Id), data = data_Controlli_TS_LME_Ext,
                 REML = FALSE)
Anova(LME_SSRT_long, type = 'III')
qqp(data_Controlli_TS_LME_Ext$SSRT)
qqp(resid(LME_SSRT_long), "norm")

LME_SSRT <- lmer(SSRT ~ 1 + Protocol + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + log_k_Pre + (1|Id), data = data_Controlli_TS_LME_Ext,
                 REML = FALSE)
Anova(LME_SSRT, type = 'III')
Anova(LME_SSRT, type = 'II')
qqp(data_Controlli_TS_LME_Ext$SSRT)
qqp(resid(LME_SSRT), "norm")
#Kenward-Rogers
anova(LME_SSRT,LME_SSRT_long)
KRmodcomp(LME_SSRT,LME_SSRT_long)

# CONCL. Meglio il modello con solo gli effetti principali (Sign.: Protocol, MCQ30_POS, MCQ_30_NEG e MCQ30_NC)
#        Nel modello con le interazioni non ho l'effetto prinicipale del protocollo, ma un'interazione significativa tra Protocol e MCQ30_NC (significativi anche main effects di MCQ30_NEG e MCQ_30_CC)

LME_SSRT_long_2 <- lmer(SSRT ~ 1 + Protocol*MSAS_18_sr + Protocol*MSAS_18_cd + Protocol*MSAS_18_so + Protocol*MSAS_18_m + (1|Id), data = data_Controlli_TS_LME_Ext,
                      REML = FALSE)
Anova(LME_SSRT_long_2, type = 'III')
qqp(data_Controlli_TS_LME_Ext$SSRT)
qqp(resid(LME_SSRT_long_2), "norm")

LME_SSRT_2 <- lmer(SSRT ~ 1 + Protocol + MSAS_18_sr + MSAS_18_cd + MSAS_18_so + MSAS_18_m + (1|Id), data = data_Controlli_TS_LME_Ext,
                        REML = FALSE)
Anova(LME_SSRT_2, type = 'III')
Anova(LME_SSRT_2, type = 'II')
qqp(data_Controlli_TS_LME_Ext$SSRT)
qqp(resid(LME_SSRT_2), "norm")
#Kenward-Rogers
anova(LME_SSRT_2,LME_SSRT_long_2)
KRmodcomp(LME_SSRT_2,LME_SSRT_long_2)

# CONCL. Leggermente meglio il modello con le interazioni. Non ho main effect del protocollo ma interazioni significative tra Protocol e MSAS18_SR/MSAS18_CD/MSAS18_SO (significativi anche main effect di MSAS18_SR e MSAS18_SO)

LME_SSRT_long_3a <- lmer(SSRT ~ 1 + Protocol*MCQ_30_pos + Protocol*MCQ_30_neg + Protocol*MCQ_30_cc + Protocol*MCQ_30_nc + Protocol*MCQ_30_csc + Protocol*MSAS_18_sr + Protocol*MSAS_18_cd + Protocol*MSAS_18_so + Protocol*MSAS_18_m + Protocol*log_k_Pre + (1|Id), data = data_Controlli_TS_LME_Ext,
                         REML = FALSE)
Anova(LME_SSRT_long_3a, type = 'III')
qqp(data_Controlli_TS_LME_Ext$SSRT)
qqp(resid(LME_SSRT_long_3a), "norm")

LME_SSRT_3a <- lmer(SSRT ~ 1 + Protocol + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + MSAS_18_sr + MSAS_18_cd + MSAS_18_so + MSAS_18_m + log_k_Pre + (1|Id), data = data_Controlli_TS_LME_Ext,
                         REML = FALSE)
Anova(LME_SSRT_3a, type = 'III')
Anova(LME_SSRT_3a, type = 'II')
qqp(data_Controlli_TS_LME_Ext$SSRT)
qqp(resid(LME_SSRT_long_3a), "norm")
#Kenward-Rogers
anova(LME_SSRT_3a,LME_SSRT_long_3a)
KRmodcomp(LME_SSRT_3a,LME_SSRT_long_3a)

#CONCL. Utilizzando le sottoscale sia di MCQ30 che MSAS risulta migliore il modello con solamente gli effetti principali (Sign.:Protocol, MCQ30_POS, MCQ30_NEG, MCQ30_NC, MSAS18_SR)
#       Nel modello con le interazioni non ho main effect del protocollo, ma un'interazione significativa tra Protocol e MSAS18_CD (Significativi anche i main effect di MCQ30_NEG, MCQ30_CC, MSAS18_SR, MSAS18_CD)

LME_SSRT_long_3b <- lmer(SSRT ~ 1 + Protocol*MSAS_18_tot + Protocol*MCQ_30_tot + Protocol*log_k_Pre + (1|Id), data = data_Controlli_TS_LME_Ext,
                        REML = FALSE)
Anova(LME_SSRT_long_3b, type = 'III')
qqp(data_Controlli_TS_LME_Ext$SSRT)
qqp(resid(LME_SSRT_long_3b), "norm")

LME_SSRT_3b <- lmer(SSRT ~ 1 + Protocol + MSAS_18_tot + MCQ_30_tot + log_k_Pre + (1|Id), data = data_Controlli_TS_LME_Ext,
                        REML = FALSE)
Anova(LME_SSRT_3b, type = 'III')
Anova(LME_SSRT_3b, type = 'II')
qqp(data_Controlli_TS_LME_Ext$SSRT)
qqp(resid(LME_SSRT_3b), "norm")
#Kenward-Rogers
anova(LME_SSRT_3b,LME_SSRT_long_3b)
KRmodcomp(LME_SSRT_3b,LME_SSRT_long_3b)

#CONCL. Utilizzando i totali di MCQ30 e MSAS18 risulta migliore il modello con solamente gli effetti principali (Sign.:Protocol)
#       Nel modello con le interazioni non ho main effect del protocollo e nessuna interazione significativa (Significativo solo il main effect di MCQ30_tot)

###     Forse la cosa migliore sarebbe testare i due modelli separati (i.e., solo main effect e solo interactions) senza confrontrali con Kenward-Rogers

##### LME Go-noGO #####
LME_RT_Go_long <- lmer(Go_noGo_Mean_RT_GoTrials ~ 1 + Protocol*MCQ_30_pos + Protocol*MCQ_30_neg + Protocol*MCQ_30_cc + Protocol*MCQ_30_nc + Protocol*MCQ_30_csc + Protocol*log_k_Pre + (1|Id), data = data_Controlli_TS_LME_Ext,
                       REML = FALSE)
Anova(LME_RT_Go_long, type = 'III')
qqp(data_Controlli_TS_LME_Ext$Go_noGo_Mean_RT_GoTrials)
qqp(resid(LME_RT_Go_long), "norm")

LME_RT_Go <- lmer(Go_noGo_Mean_RT_GoTrials ~ 1 + Protocol + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + log_k_Pre + (1|Id), data = data_Controlli_TS_LME_Ext,
                 REML = FALSE)
Anova(LME_RT_Go, type = 'III')
Anova(LME_RT_Go, type = 'II')
qqp(data_Controlli_TS_LME_Ext$Go_noGo_Mean_RT_GoTrials)
qqp(resid(LME_RT_Go), "norm")
#Kenward-Rogers
anova(LME_RT_Go,LME_RT_Go_long)
KRmodcomp(LME_RT_Go,LME_RT_Go_long)

LME_RT_Go_long_3b <- lmer(Go_noGo_Mean_RT_GoTrials ~ 1 + Protocol*MSAS_18_tot + Protocol*MCQ_30_tot + Protocol*log_k_Pre + (1|Id), data = data_Controlli_TS_LME_Ext,
                         REML = FALSE)
Anova(LME_RT_Go_long_3b, type = 'III')
qqp(data_Controlli_TS_LME_Ext$Go_noGo_Mean_RT_GoTrials)
qqp(resid(LME_RT_Go_long_3b), "norm")

LME_RT_Go_3b <- lmer(Go_noGo_Mean_RT_GoTrials ~ 1 + Protocol + MSAS_18_tot + MCQ_30_tot + log_k_Pre + (1|Id), data = data_Controlli_TS_LME_Ext,
                    REML = FALSE)
Anova(LME_RT_Go_3b, type = 'III')
Anova(LME_RT_Go_3b, type = 'II')
qqp(data_Controlli_TS_LME_Ext$Go_noGo_Mean_RT_GoTrials)
qqp(resid(LME_RT_Go_3b), "norm")
#Kenward-Rogers
anova(LME_RT_Go_3b,LME_RT_Go_long_3b)
KRmodcomp(LME_RT_Go_3b,LME_RT_Go_long_3b)

LME_FA_long <- lmer(Go_no_Go_Mean_FA ~ 1 + Protocol*MCQ_30_pos + Protocol*MCQ_30_neg + Protocol*MCQ_30_cc + Protocol*MCQ_30_nc + Protocol*MCQ_30_csc + Protocol*log_k_Pre + (1|Id), data = data_Controlli_TS_LME_Ext,
                    REML = FALSE)
Anova(LME_FA_long, type = 'III')
qqp(data_Controlli_TS_LME_Ext$Go_no_Go_Mean_FA)
qqp(resid(LME_FA_long), "norm")

LME_FA <- lmer(Go_no_Go_Mean_FA ~ 1 + Protocol + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + log_k_Pre + (1|Id), data = data_Controlli_TS_LME_Ext,
                  REML = FALSE)
Anova(LME_FA, type = 'III')
Anova(LME_FA, type = 'II')
qqp(data_Controlli_TS_LME_Ext$Go_no_Go_Mean_FA)
qqp(resid(LME_FA), "norm")
#Kenward-Rogers
anova(LME_FA,LME_FA_long)
KRmodcomp(LME_FA,LME_FA_long)

LME_FA_long_3b <- lmer(Go_no_Go_Mean_FA ~ 1 + Protocol*MCQ_30_tot + Protocol*MSAS_18_tot + Protocol*log_k_Pre + (1|Id), data = data_Controlli_TS_LME_Ext,
                    REML = FALSE)
Anova(LME_FA_long_3b, type = 'III')
qqp(data_Controlli_TS_LME_Ext$Go_no_Go_Mean_FA)
qqp(resid(LME_FA_long_3b), "norm")

LME_FA_3b <- lmer(Go_no_Go_Mean_FA ~ 1 + Protocol + MCQ_30_tot + MSAS_18_tot + log_k_Pre + (1|Id), data = data_Controlli_TS_LME_Ext,
                       REML = FALSE)
Anova(LME_FA_3b, type = 'III')
Anova(LME_FA_3b, type = 'II')
qqp(data_Controlli_TS_LME_Ext$Go_no_Go_Mean_FA)
qqp(resid(LME_FA_3b), "norm")
#Kenward-Rogers
anova(LME_FA_3b,LME_FA_long_3b)
KRmodcomp(LME_FA_3b,LME_FA_long_3b)

##### LME DD #####
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

LME_DD_Low <- lmer(log_k ~ 1 + Protocol*Time_DD + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + (1|Id), data = subset_low,
               REML = FALSE)
Anova(LME_DD_Low, type = 'III')
Anova(LME_DD_Low, type = 'II')
qqp(subset_low$log_k)
qqp(resid(LME_DD_Low), "norm")

LME_DD_Low_long <- lmer(log_k ~ 1 + Protocol*Time_DD + Protocol*MCQ_30_tot + Protocol*MSAS_18_tot + (1|Id), data = subset_low,
                   REML = FALSE)
Anova(LME_DD_Low_long, type = 'III')
qqp(subset_low$log_k)
qqp(resid(LME_DD_Low_long), "norm")

LME_DD_Low_short <- lmer(log_k ~ 1 + Protocol*Time_DD + (1|Id), data = subset_low,
                   REML = FALSE)
Anova(LME_DD_Low_short, type = 'III')
qqp(subset_low$log_k)
qqp(resid(LME_DD_Low_short), "norm")
#Kenward-Rogers
anova(LME_DD_Low_short,LME_DD_Low_long)
KRmodcomp(LME_DD_Low_short,LME_DD_Low_long)

#CONCL. Meglio il modello ristretto con salmente l'interazione Protocol*Time (Sign.: Protocol, Time, Protocol*Time)

LME_DD_High <- lmer(log_k ~ 1 + Protocol*Time_DD + MCQ_30_pos + MCQ_30_neg + MCQ_30_cc + MCQ_30_nc + MCQ_30_csc + (1|Id), data = subset_high,
                   REML = FALSE)
Anova(LME_DD_High, type = 'III')
qqp(subset_high$log_k)
qqp(resid(LME_DD_High), "norm")

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
