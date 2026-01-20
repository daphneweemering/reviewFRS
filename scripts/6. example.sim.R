# Because the Ceftriaxone ALS trial data is not publicly available, we created
# simulated data based on the joint modeling framework. We then fit the different
# models (without covariates) to the simulated data as an example. 

# get code to simulate data based on the joint modeling framework
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

# Simulate data with:
# - 500 patients (N = 500)
# - 7 trial measurements (Nm = 7)
# - 4 pre-baseline measurements (Nm.pre = 4)
# - No treatment effect on function or survival (B2 = 0, B4 = 0)
# - 12 months of follow-up (FUT = 12)
# - Dropout of 10% 
D <- JMD(N = 500, Nm = 7, Nm.pre = 4, B2 = 0, B4 = 0, g1 = 0, FUT = 12, DO.TRT = 0.1, DO.PLB = 0.1) 

# function to extract treatment effect estimate
est <- function(x){
  if (class(x)[1] == "lm"){
    e <- unname(coefficients(x)["TRT"])
  } else if (class(x)[1] == "mmrm"){
    e <- (as.data.frame(contrast(emmeans(x, ~ TRT), method = "pairwise"))$estimate)*-1
  } else if (class(x)[1] == "lmerMod"){
    e <- unname(fixef(x)["TIME:TRT"])
  }
  return(e)
}

# fit the different models (same models as in "5.models.JM.R")
m1 <- lm(CFB_6 ~ TRT, data = D$LOCF)
e1 <- est(m1)
p1 <- drop1(m1, test = "Chisq")[2, 5]

m2 <- survdiff(Surv(STIME2, EVENT2) ~ TRT, data = D$WD.AGG)
e2 <- (m2$obs[2] - m2$exp[2]) * -1 # *-1 so that positive e2 means TRT == 1 was better off. 
p2 <- m2$pvalue

m3 <- mmrm(CFB ~ TRT + VISIT + TRT:VISIT + us(VISIT|ID), data = D$LD[!D$LD$VISIT == 0, ])
e3 <- (as.data.frame(contrast(emmeans(m3, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$estimate)*-1
p3 <- as.data.frame(contrast(emmeans(m3, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$p.value

m4 <- anova_test(data = D$LD[!is.na(D$LD$CFB), ], dv = CFB, wid = ID, within = VISIT, between = TRT, type = 3)
e4 <- (as.data.frame(contrast(emmeans(lmer(CFB ~ TRT*VISIT + (1 | ID), data = D$LD[!is.na(D$LD$CFB), ]), ~ TRT | VISIT, at = list(VISIT = "6")), 
                              method = "pairwise"))$estimate)*-1
p4 <- m4$`Sphericity Corrections`$`p[GG]`[2] 

m5 <- lmer(TOT ~ TIME + TRT + TRT:TIME + (TIME|ID), data = D$LD, control = lmerControl(optimizer = "nlminbwrap"))
e5 <- est(m5)
p5 <- drop1(m5, test = "Chisq")[2, 4]

m6 <- lm(RANK ~ TRT + TOT_0, data = D$WD.JR)
e6 <- est(m6)
p6 <- drop1(m6, test = "Chisq")[2, 5]

m7 <- anova_test(data = D$LOCF2, dv = CFB, wid = ID, within = VISIT, between = TRT, type = 3)
e7 <- (as.data.frame(contrast(emmeans(lmer(CFB ~ TRT*VISIT + (1 | ID), data = D$LOCF2[!is.na(D$LOCF2$CFB), ]), ~ TRT | VISIT, at = list(VISIT = "6")), 
                              method = "pairwise"))$estimate)*-1
p7 <- m7$`Sphericity Corrections`$`p[GG]`[2] 

m8 <- lm(TOT_6 ~ TRT, data = D$CC)
e8 <- est(m8)
p8 <- drop1(m8, test = "Chisq")[2, 5]

m9 <- lm(RANK.genge ~ TRT + TOT_0, data = D$WD.JR)
e9 <- est(m9)
p9 <- drop1(m9, test = "Chisq")[2, 5]

m10 <- coin::wilcox_test(CFB_6 ~ as.factor(TRT), data = D$ZERO, distribution = "exact")
e10 <- m10@statistic@teststatistic
p10 <- coin::pvalue(m10)

m11 <- coin::wilcox_test(CFB_6 ~ as.factor(TRT), data = D$WD.KAJI, distribution = "exact")
e11 <- m11@statistic@teststatistic
p11 <- coin::pvalue(m11)

m12 <- lm(CFB_6 ~ TRT, data = D$WD.KAUF)
e12 <- est(m12)
p12 <- drop1(m12, test = "Chisq")[2, 5]

m13 <- lmer(CFB ~ TIME:TRT + BSLN:TIME + (TIME|ID), data = D$LD[!D$LD$VISIT == 0, ])
e13 <- (as.data.frame(contrast(emmeans(m13, ~ TRT | TIME, at = list(TIME = 12)), method = "pairwise"))$estimate)*-1
p13 <- as.data.frame(contrast(emmeans(m13, ~ TRT | TIME, at = list(TIME = 12)), method = "pairwise"))$p.value

m14 <- coin::wilcox_test(I((TOT_6 - TOT_0)/TIME_6) ~ as.factor(TRT), data = D$CC, distribution = "exact")
e14 <- m14@statistic@teststatistic
p14 <- coin::pvalue(m14)

m15 <- lmer(TOT ~ TIME + TRT + TIME:TRT + BSLN + (TIME|ID), data = D$LD, control = lmerControl(optimizer = "nlminbwrap"))
e15 <- est(m15)
p15 <- drop1(m15, test = "Chisq")[3, 4]

m16 <- lm(SUM ~ TRT + TOT_0, data = D$WD.JR)
e16 <- est(m16)
p16 <- drop1(m16, test = "Chisq")[2, 5]

m17 <- coin::wilcox_test(CFB_6 ~ as.factor(TRT), data = D$CC, distribution = "exact")
e17 <- m17@statistic@teststatistic
p17 <- coin::pvalue(m17)

m18 <- lm(CFB_6 ~ TRT + TOT_0, data = D$LOCF.MORA)
e18 <- est(m18)
p18 <- drop1(m18, test = "Chisq")[2, 5]

m19 <- lm(CFB_6 ~ TRT, data = D$CC)
e19 <- est(m19)
p19 <- drop1(m19, test = "Chisq")[2, 5]

m20 <- lmer(TOT ~ TRT + TIME + (1|ID), data = D$LD)
e20 <- unname(fixef(m20)["TRT"])
p20 <- drop1(m20, test = "Chisq")[2, 4]

m21 <- lmer(TOT ~ TIME + TRT + TIME:TRT + BSLN + (1|ID), data = D$LD)
e21 <- est(m21)
p21 <- drop1(m21, test = "Chisq")[3, 4]

m22 <- mmrm(CFB ~ TRT + VISIT + TRT:VISIT + BSLN + BSLN:VISIT + us(VISIT|ID), data = D$LD[!D$LD$VISIT == 0, ])
e22 <- (as.data.frame(contrast(emmeans(m22, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$estimate)*-1
p22 <- as.data.frame(contrast(emmeans(m22, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$p.value

m23 <- mmrm(CFB ~ TRT:VISIT + VISIT + BSLN + us(VISIT|ID), data = D$LD[!D$LD$VISIT == 0, ])
e23 <- (as.data.frame(contrast(emmeans(m23, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$estimate)*-1
p23 <- as.data.frame(contrast(emmeans(m23, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$p.value

m24 <- mmrm(TOT ~ VISIT + TRT:VISIT + us(VISIT|ID), data = D$LD, reml = F)
m24.1 <- mmrm(TOT ~ VISIT + us(VISIT|ID), data = D$LD, reml = F)
e24 <- (as.data.frame(contrast(emmeans(m24, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$estimate)*-1
p24 <- as.data.frame(contrast(emmeans(m24, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$p.value

m25 <- mmrm(CFB ~ TRT:VISIT + VISIT + BSLN + us(VISIT|ID), data = D$LOCF2[!D$LOCF2$VISIT == 0, ])
e25 <- (as.data.frame(contrast(emmeans(m25, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$estimate)*-1
p25 <- as.data.frame(contrast(emmeans(m25, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$p.value

m26 <- lm(CFB_6 ~ TRT + I(TOT_0 >= 37), data = D$CC)
e26 <- est(m26)
p26 <- drop1(m26, test = "Chisq")[2, 5]

m27 <- lmer(SLP ~ TRT + (1|ID), data = D$LOCF.BERRY)
e27 <- unname(fixef(m27)["TRT"])
p27 <- drop1(m27, test = "Chisq")[2, 4]

m28 <- glm(RESP ~ TRT + BSLN, family = binomial(link = "logit"), data = D$RESP.CUD)
e28 <- unname(coefficients(m28)["TRT"])
p28 <- drop1(m28, ~ TRT, test = "Chisq")[2, 5]

m29 <- glm(RESP ~ TRT, family = binomial(link = "logit"), data = D$RESP.ELIA)
e29 <- unname(coefficients(m29)["TRT"])
p29 <- drop1(m29, test = "Chisq")[2, 5]

m30 <- lapply(1:5, function(i){
  lm(RANK ~ TRT + TOT_0, data = D$WD.MILL[D$WD.MILL$.imp == i, ])
})
estimates <- pool(m30)
e30 <- summary(estimates)[2, 2]
p30 <- summary(estimates)[2, 6]

