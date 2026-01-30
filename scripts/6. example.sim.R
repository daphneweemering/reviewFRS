##### Because the Ceftriaxone ALS trial data is not publicly available, we created    #####
##### simulated data based on the joint modeling framework. We then fit the different #####
##### models (without covariates) to the simulated data as an example. Effect size    #####
##### and p-values are extracted.                                                     #####
###########################################################################################
 
# 1. Source the simulation code in order to get the data files.
source("4. sim.JM.R")

# 2. Simulate data and assign all the data files to D (list).  
D <- JMD(N = 500, Nm = 7, Nm.pre = 4, B2 = 0, B4 = 0, g1 = 0, FUT = 12, DO.TRT = 0.1, DO.PLB = 0.1) 
## Simulate data with:
# - 500 patients (N = 500)
# - 7 trial measurements (Nm = 7)
# - 4 pre-baseline measurements (Nm.pre = 4)
# - No treatment effect on function or survival (B2 = 0, B4 = 0)
# - 12 months of follow-up (FUT = 12)
# - Dropout of 10% 

## D CONTAINS: 
# - LD: long data file.
# - WD: wide data file.
# - LOCF: wide data file with last observation carried forward [LOCF] for all missing values.
# - ZERO: wide data file with zero score imputation for all missing values.
# - CC: wide data file with only complete observations. All rows with missing values are removed.
# - LOCF2: long data file with LOCF for all missing values.
# - WD.AGG: adjustments for the analysis method of Aggarwal [2010]. Includes STIME2 indicating time to death or 6-point decrease in ALSFRS-R.
# - WD.JR: adjustments for joint rank analyses so that they include the rank score and the sum score.
# - WD.MILL: adjustments for the joint rank analysis of Miller [2022], including multiple imputed datasets
# - LOCF.MORA: adjustments for the analysis method of Mora [2020]. Includes zero imputation for death-related missingness, and LOCF for non-death-related missingness.
# - WD.KAJI: adjustments for the analysis method of Kaji [2019]. Includes zero/worst score imputation for death-related missingness, and LOCF for non-death-related missingness.
# - WD.KAUF: adjustments for the analysis method of Kaufmann [2009]. Includes zero imputation for death-related missingness, and neirest-neighbor worst score imputation for non-death-related missingness.
# - LOCF.BERRY: adjustments for the analysis method of Berry [2019]. Includes pre-treatment slope and post-treatment slope.
# - RESP.ELIA: adjustments for the analysis method of Elia [2015]. Includes pre-treatment slope and post-treatment slope. Missing values are imputed with LOCF.
# - RESP.CUD: adjustments for the analysis method of Cudkowicz [2022]. Includes pre-treatment slope and post-treatment slope. Patients with death-related missingness were marked as non-responder, for other missingness LOCF was used.


# Function to extract treatment effect estimate
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

# Fit the models (same models as in "5.models.JM.R")

# For each model, we provide the following information:
# Outcome; analysis approach; missing data strategy death-events; missing data strategy
# non-death events; baseline adjustment 

# Change from baseline; linear regression; LOCF; LOCF; no baseline adjustment
m1 <- lm(CFB_6 ~ TRT, data = D$LOCF)
e1 <- est(m1)
p1 <- drop1(m1, test = "Chisq")[2, 5]

# Time-to-event; logrank test; event; censored; no baseline adjustment
m2 <- survdiff(Surv(STIME2, EVENT2) ~ TRT, data = D$WD.AGG)
e2 <- (m2$obs[2] - m2$exp[2]) * -1 # *-1 so that positive e2 means TRT == 1 was better off. 
p2 <- m2$pvalue

# Change from baseline; MMRM; direct modeling; direct modeling; no baseline adjustment
m3 <- mmrm(CFB ~ TRT + VISIT + TRT:VISIT + us(VISIT|ID), data = D$LD[!D$LD$VISIT == 0, ])
e3 <- (as.data.frame(contrast(emmeans(m3, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$estimate)*-1
p3 <- as.data.frame(contrast(emmeans(m3, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$p.value

# Change from baseline; repeated measures ANOVA; complete case; complete case; no baseline adjustment
m4 <- anova_test(data = D$LD[!is.na(D$LD$CFB), ], dv = CFB, wid = ID, within = VISIT, between = TRT, type = 3)
e4 <- (as.data.frame(contrast(emmeans(lmer(CFB ~ TRT*VISIT + (1 | ID), data = D$LD[!is.na(D$LD$CFB), ]), ~ TRT | VISIT, at = list(VISIT = "6")), 
                              method = "pairwise"))$estimate)*-1
p4 <- m4$`Sphericity Corrections`$`p[GG]`[2] 

# Progression rate; mixed effects model; direct modeling; direct modeling; no baseline adjustment
m5 <- lmer(TOT ~ TIME + TRT + TRT:TIME + (TIME|ID), data = D$LD, control = lmerControl(optimizer = "nlminbwrap"))
e5 <- est(m5)
p5 <- drop1(m5, test = "Chisq")[2, 4]

# Rank; joint rank analysis; joint rank analysis; comparison at last common visit; baseline adjustment
m6 <- lm(RANK ~ TRT + TOT_0, data = D$WD.JR)
e6 <- est(m6)
p6 <- drop1(m6, test = "Chisq")[2, 5]

# Change from baseline; repeated measures ANOVA; LOCF; LOCF; no baseline adjustment
m7 <- anova_test(data = D$LOCF2, dv = CFB, wid = ID, within = VISIT, between = TRT, type = 3)
e7 <- (as.data.frame(contrast(emmeans(lmer(CFB ~ TRT*VISIT + (1 | ID), data = D$LOCF2[!is.na(D$LOCF2$CFB), ]), ~ TRT | VISIT, at = list(VISIT = "6")), 
                              method = "pairwise"))$estimate)*-1
p7 <- m7$`Sphericity Corrections`$`p[GG]`[2] 

# Total score; linear regression; complete case; complete case; no baseline adjustment
m8 <- lm(TOT_6 ~ TRT, data = D$CC)
e8 <- est(m8)
p8 <- drop1(m8, test = "Chisq")[2, 5]

# Rank; joint rank analysis; joint rank analysis; slope; baseline adjustment
m9 <- lm(RANK.genge ~ TRT + TOT_0, data = D$WD.JR)
e9 <- est(m9)
p9 <- drop1(m9, test = "Chisq")[2, 5]

# Change from baseline; non-parametric test; zero imputation; zero imputation; no baseline adjustment
m10 <- coin::wilcox_test(CFB_6 ~ as.factor(TRT), data = D$ZERO, distribution = "exact")
e10 <- m10@statistic@teststatistic
p10 <- coin::pvalue(m10)

# Change from baseline; non-parametric test; zero imputation; LOCF; no baseline adjustment
m11 <- coin::wilcox_test(CFB_6 ~ as.factor(TRT), data = D$WD.KAJI, distribution = "exact")
e11 <- m11@statistic@teststatistic
p11 <- coin::pvalue(m11)

# Change from baseline; linear regression; zero imputation; worst score imputation; no baseline adjustment
m12 <- lm(CFB_6 ~ TRT, data = D$WD.KAUF)
e12 <- est(m12)
p12 <- drop1(m12, test = "Chisq")[2, 5]

# Change from baseline; mixed effects model; direct modeling; direct modeling; no baseline adjustment (interaction) 
# (looks at adjusted mean difference at the end of the study instead of progression rates)
m13 <- lmer(CFB ~ TIME:TRT + BSLN:TIME + (TIME|ID), data = D$LD[!D$LD$VISIT == 0, ])
e13 <- (as.data.frame(contrast(emmeans(m13, ~ TRT | TIME, at = list(TIME = 12)), method = "pairwise"))$estimate)*-1
p13 <- as.data.frame(contrast(emmeans(m13, ~ TRT | TIME, at = list(TIME = 12)), method = "pairwise"))$p.value

# Progression rate; non-parametric test; complete case; complete case; no baseline adjustment
m14 <- coin::wilcox_test(I((TOT_6 - TOT_0)/TIME_6) ~ as.factor(TRT), data = D$CC, distribution = "exact")
e14 <- m14@statistic@teststatistic
p14 <- coin::pvalue(m14)

# Progression rate; mixed effects model; direct modeling; direct modeling; baseline adjustment
m15 <- lmer(TOT ~ TIME + TRT + TIME:TRT + BSLN + (TIME|ID), data = D$LD, control = lmerControl(optimizer = "nlminbwrap"))
e15 <- est(m15)
p15 <- drop1(m15, test = "Chisq")[3, 4]

# Sum score; joint rank analysis; joint rank analysis; comparison at last common visit; baseline adjustment
m16 <- lm(SUM ~ TRT + TOT_0, data = D$WD.JR)
e16 <- est(m16)
p16 <- drop1(m16, test = "Chisq")[2, 5]

# Change from baseline; non-parametric test; complete case; complete case; no baseline adjustment
m17 <- coin::wilcox_test(CFB_6 ~ as.factor(TRT), data = D$CC, distribution = "exact")
e17 <- m17@statistic@teststatistic
p17 <- coin::pvalue(m17)

# Change from baseline; linear regression; zero score imputation; LOCF; baseline adjustment
m18 <- lm(CFB_6 ~ TRT + TOT_0, data = D$LOCF.MORA)
e18 <- est(m18)
p18 <- drop1(m18, test = "Chisq")[2, 5]

# Change from baseline; linear regression; complete case; complete case; no baseline adjustment
m19 <- lm(CFB_6 ~ TRT, data = D$CC)
e19 <- est(m19)
p19 <- drop1(m19, test = "Chisq")[2, 5]

# Progression rate; random-intercepts model; direct modeling; direct modeling; no baseline adjustment
m20 <- lmer(TOT ~ TRT + TIME + (1|ID), data = D$LD)
e20 <- unname(fixef(m20)["TRT"])
p20 <- drop1(m20, test = "Chisq")[2, 4]

# Progression rate; random-intercepts model; direct modeling; direct modeling; baseline adjustment
m21 <- lmer(TOT ~ TIME + TRT + TIME:TRT + BSLN + (1|ID), data = D$LD)
e21 <- est(m21)
p21 <- drop1(m21, test = "Chisq")[3, 4]

# Change from baseline; MMRM; direct modeling; direct modeling; baseline adjustment
m22 <- mmrm(CFB ~ TRT + VISIT + TRT:VISIT + BSLN + BSLN:VISIT + us(VISIT|ID), data = D$LD[!D$LD$VISIT == 0, ])
e22 <- (as.data.frame(contrast(emmeans(m22, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$estimate)*-1
p22 <- as.data.frame(contrast(emmeans(m22, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$p.value

# Change from baseline; MMRM; direct modeling; direct modeling; baseline adjustment
m23 <- mmrm(CFB ~ TRT:VISIT + VISIT + BSLN + us(VISIT|ID), data = D$LD[!D$LD$VISIT == 0, ])
e23 <- (as.data.frame(contrast(emmeans(m23, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$estimate)*-1
p23 <- as.data.frame(contrast(emmeans(m23, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$p.value

# Progression rate; MMRM; direct modeling; direct modeling; no baseline adjustment
m24 <- mmrm(TOT ~ VISIT + TRT:VISIT + us(VISIT|ID), data = D$LD, reml = F)
m24.1 <- mmrm(TOT ~ VISIT + us(VISIT|ID), data = D$LD, reml = F)
e24 <- (as.data.frame(contrast(emmeans(m24, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$estimate)*-1
p24 <- as.data.frame(contrast(emmeans(m24, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$p.value

# Change from baseline; MMRM; direct modeling; direct modeling; baseline adjustment
m25 <- mmrm(CFB ~ TRT:VISIT + VISIT + BSLN + us(VISIT|ID), data = D$LOCF2[!D$LOCF2$VISIT == 0, ])
e25 <- (as.data.frame(contrast(emmeans(m25, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$estimate)*-1
p25 <- as.data.frame(contrast(emmeans(m25, ~ TRT | VISIT, at = list(VISIT = "6")), method = "pairwise"))$p.value

# Change from baseline; linear regression; complete case; complete case; baseline adjustment
m26 <- lm(CFB_6 ~ TRT + I(TOT_0 >= 37), data = D$CC)
e26 <- est(m26)
p26 <- drop1(m26, test = "Chisq")[2, 5]

# Progression rate; two-step linear regression; available case analysis; available case analysis; no baseline adjustment
m27 <- lmer(SLP ~ TRT + (1|ID), data = D$LOCF.BERRY)
e27 <- unname(fixef(m27)["TRT"])
p27 <- drop1(m27, test = "Chisq")[2, 4]

# Proportion responders; logistic regression; non-responder; LOCF; baseline adjustment
m28 <- glm(RESP ~ TRT + BSLN, family = binomial(link = "logit"), data = D$RESP.CUD)
e28 <- unname(coefficients(m28)["TRT"])
p28 <- drop1(m28, ~ TRT, test = "Chisq")[2, 5]

# Proportion responders; logistic regression; LOCF; LOCF; no baseline adjustment
m29 <- glm(RESP ~ TRT, family = binomial(link = "logit"), data = D$RESP.ELIA)
e29 <- unname(coefficients(m29)["TRT"])
p29 <- drop1(m29, test = "Chisq")[2, 5]

# Rank; joint rank analysis; joint rank analysis; multiple imputation; baseline adjustment
m30 <- lapply(1:5, function(i){
  lm(RANK ~ TRT + TOT_0, data = D$WD.MILL[D$WD.MILL$.imp == i, ])
})
estimates <- pool(m30)
e30 <- summary(estimates)[2, 2]
p30 <- summary(estimates)[2, 6]

