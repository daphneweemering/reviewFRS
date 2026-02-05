################################################################################
##                                                                            ##
##                                SCENARIOS                                   ##
##                                                                            ##
################################################################################
source("4. sim.JM.R")

library(MASS)
library(lme4)
library(JM)
library(data.table)
library(zoo)
library(mmrm)
library(survival)
library(parallel)
library(rstatix)
library(emmeans)
library(coin)
library(mice)
library(pbmcapply)
library(tidyverse)

set.seed(3239480)

X <- pbmclapply(1:10000, function(i){
  # Get the simulated data _____________________________________________________
  D <- JMD(N = 500, Nm = 7, Nm.pre = 4, B2 = 0, B4 = 0, g1 = 0, FUT = 12, DO.TRT = 0.1, DO.PLB = 0.1)
  
  # Fit the models _____________________________________________________________
  m1 <- lm(CFB_6 ~ TRT, data = D$LOCF)
  p1 <- drop1(m1, test = "Chisq")[2, 5]
  
  m2 <- survdiff(Surv(STIME2, EVENT2) ~ TRT, data = D$WD.AGG)
  p2 <- m2$pvalue
  
  m3 <- mmrm(CFB ~ TRT + VISIT + TRT:VISIT + us(VISIT|ID), data = D$LD[!D$LD$VISIT == 0, ])
  p3 <- as.data.frame(contrast(emmeans(m3, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$p.value
  
  m4 <- anova_test(data = D$LD[!is.na(D$LD$CFB), ], dv = CFB, wid = ID, within = VISIT, between = TRT, type = 3)
  p4 <- m4$`Sphericity Corrections`$`p[GG]`[2] 
  
  m5 <- lmer(TOT ~ TIME + TRT + TRT:TIME + (TIME|ID), data = D$LD, control = lmerControl(optimizer = "nlminbwrap"))
  p5 <- drop1(m5, test = "Chisq")[2, 4]
  
  m6 <- lm(RANK ~ TRT + TOT_0, data = D$WD.JR)
  p6 <- drop1(m6, test = "Chisq")[2, 5]
  
  m7 <- anova_test(data = D$LOCF2, dv = CFB, wid = ID, within = VISIT, between = TRT, type = 3)
  p7 <- m7$`Sphericity Corrections`$`p[GG]`[2] 
  
  m8 <- lm(TOT_6 ~ TRT, data = D$CC)
  p8 <- drop1(m8, test = "Chisq")[2, 5]
  
  m9 <- lm(RANK.genge ~ TRT + TOT_0, data = D$WD.JR)
  p9 <- drop1(m9, test = "Chisq")[2, 5]
  
  m10 <- coin::wilcox_test(CFB_6 ~ as.factor(TRT), data = D$ZERO, distribution = "exact")
  p10 <- coin::pvalue(m10)
  
  m11 <- coin::wilcox_test(CFB_6 ~ as.factor(TRT), data = D$WD.KAJI, distribution = "exact")
  p11 <- coin::pvalue(m11)
  
  m12 <- lm(CFB_6 ~ TRT, data = D$WD.KAUF)
  p12 <- drop1(m12, test = "Chisq")[2, 5]
  
  m13 <- lmer(CFB ~ TIME:TRT + BSLN:TIME + (TIME|ID), data = D$LD[!D$LD$VISIT == 0, ])
  p13 <- as.data.frame(contrast(emmeans(m13, ~ TRT | TIME, at = list(TIME = 12)), method = "pairwise"))$p.value
  
  m14 <- coin::wilcox_test(I((TOT_6 - TOT_0)/TIME_6) ~ as.factor(TRT), data = D$CC, distribution = "exact")
  p14 <- coin::pvalue(m14)
  
  m15 <- lmer(TOT ~ TIME + TRT + TIME:TRT + BSLN + (TIME|ID), data = D$LD, control = lmerControl(optimizer = "nlminbwrap"))
  p15 <- drop1(m15, test = "Chisq")[3, 4]
  
  m16 <- lm(SUM ~ TRT + TOT_0, data = D$WD.JR)
  p16 <- drop1(m16, test = "Chisq")[2, 5]
  
  m17 <- coin::wilcox_test(CFB_6 ~ as.factor(TRT), data = D$CC, distribution = "exact")
  p17 <- coin::pvalue(m17)
  
  m18 <- lm(CFB_6 ~ TRT + TOT_0, data = D$LOCF.MORA)
  p18 <- drop1(m18, test = "Chisq")[2, 5]
  
  m19 <- lm(CFB_6 ~ TRT, data = D$CC)
  p19 <- drop1(m19, test = "Chisq")[2, 5]
  
  m20 <- lmer(TOT ~ TRT + TIME + (1|ID), data = D$LD)
  p20 <- drop1(m20, test = "Chisq")[2, 4]
  
  m21 <- lmer(TOT ~ TIME + TRT + TIME:TRT + BSLN + (1|ID), data = D$LD)
  p21 <- drop1(m21, test = "Chisq")[3, 4]
  
  m22 <- mmrm(CFB ~ TRT + VISIT + TRT:VISIT + BSLN + BSLN:VISIT + us(VISIT|ID), data = D$LD[!D$LD$VISIT == 0, ])
  p22 <- as.data.frame(contrast(emmeans(m22, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$p.value
  
  m23 <- mmrm(CFB ~ TRT:VISIT + VISIT + BSLN + us(VISIT|ID), data = D$LD[!D$LD$VISIT == 0, ])
  p23 <- as.data.frame(contrast(emmeans(m23, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$p.value
  
  m24 <- mmrm(TOT ~ VISIT + TRT:VISIT + us(VISIT|ID), data = D$LD, reml = F)
  m24.1 <- mmrm(TOT ~ VISIT + us(VISIT|ID), data = D$LD, reml = F)
  p24 <- as.data.frame(contrast(emmeans(m24, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$p.value
  
  m25 <- mmrm(CFB ~ TRT:VISIT + VISIT + BSLN + us(VISIT|ID), data = D$LOCF2[!D$LOCF2$VISIT == 0, ])
  p25 <- as.data.frame(contrast(emmeans(m25, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$p.value
  
  m26 <- lm(CFB_6 ~ TRT + I(TOT_0 >= 37), data = D$CC)
  p26 <- drop1(m26, test = "Chisq")[2, 5]
  
  m27 <- lmer(SLP ~ TRT + (1|ID), data = D$LOCF.BERRY)
  p27 <- drop1(m27, test = "Chisq")[2, 4]
  
  m28 <- glm(RESP ~ TRT + BSLN, family = binomial(link = "logit"), data = D$RESP.CUD)
  p28 <- drop1(m28, ~ TRT, test = "Chisq")[2, 5]
  
  m29 <- glm(RESP ~ TRT, family = binomial(link = "logit"), data = D$RESP.ELIA)
  p29 <- drop1(m29, test = "Chisq")[2, 5]
  
  m30 <- lapply(1:5, function(i){
    lm(RANK ~ TRT + TOT_0, data = D$WD.MILL[D$WD.MILL$.imp == i, ])
  })
  estimates <- pool(m30)
  p30 <- summary(estimates)[2, 6]
  
  return(list(c(mget(paste0("p", 1:30)))))
}, mc.cores = 3)
