##### Because the Ceftriaxone ALS trial data is not publicly available, we created    #####
##### simulated data based on the joint modeling framework. We then fit the different #####
##### models (without covariates) to the simulated data as an example. Effect size    #####
##### and p-values are extracted.                                                     #####
###########################################################################################

# :::::::::::::::::::::::::::::::::::: DATA :::::::::::::::::::::::::::::::::::#
# 1. Source the simulation code in order to get the data files.
source("4. sim.JM.R")


# 2. Simulate data and assign all the data files to D (list).  
set.seed(101) # seed for reproducible results (can be changed).

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


# ::::::::::::::::::::::::::::::::: ESTIMATES :::::::::::::::::::::::::::::::::#
# est(): Extract effect estimate, p-value, and z-score from a model

# Input:
#   x - a fitted model object (various types supported, see below)

# Output:
#   A list with:
#     "effect estimate" = estimated effect of TRT (or interaction)
#     "p-value" = p-value for TRT effect
#     "z-score" = z-score with sign matching effect direction

est <- function(x){
  # Linear regression: TRT main effect
  if (class(x)[1] == "lm"){
    e <- unname(coefficients(x)["TRT"])
    p <- drop1(x, ~TRT, test = "Chisq")[2, 5]
    z <- sign(e) * abs(qnorm(p/2, lower.tail = F))
    
    # Linear mixed-effects model with interaction TIME:TRT (specific model_id m13). 
    # Looks at treatment difference at month 12
  } else if (identical(attr(x, "model_id"), "m13")){
    e <- (as.data.frame(contrast(emmeans(x, ~ TRT | TIME, at = list(TIME = 12)), method = "revpairwise"))$estimate)
    p <- (as.data.frame(contrast(emmeans(x, ~ TRT | TIME, at = list(TIME = 12)), method = "revpairwise"))$p.value)
    z <- sign(e) * abs(qnorm(p/2, lower.tail = F))
    
    # MMRM estimating progression rate (model_id m24)
  } else if (identical(attr(x, "model_id"), "m24")){
    e <- summary(x)$coefficients["VISIT2:TRT", "Estimate"]
    p <- summary(x)$coefficients["VISIT2:TRT", "Pr(>|t|)"]
    z <- sign(e) * abs(qnorm(p/2, lower.tail = F))
    
    # Pooled estimates from multiple imputation (model_id m30)
  } else if (identical(attr(x, "model_id"), "m30")){
    estimates <- pool(x)
    e <- summary(estimates)[2, 2]
    p <- summary(estimates)[2, 6]
    z <- sign(e) * abs(qnorm(p/2, lower.tail = F))
    
    # Mixed model repeated measures (mmrm)
  } else if (class(x)[1] == "mmrm"){
    e <- as.data.frame(contrast(emmeans(x, ~ TRT | VISIT, at = list(VISIT = "6")), method = "revpairwise"))$estimate
    p <- as.data.frame(contrast(emmeans(x, ~ TRT | VISIT, at = list(VISIT = "6")), method = "revpairwise"))$p.value
    z <- sign(e) * abs(qnorm(p/2, lower.tail = F))
    
    # Linear mixed model without TIME:TRT interaction, just TRT main effect
  } else if (class(x)[1] == "lmerMod" && !"TIME:TRT" %in% names(fixef(x))){
    e <- unname(fixef(x)["TRT"])
    p <- drop1(x, ~TRT, test = "Chisq")[2, 4]
    z <- sign(e) * abs(qnorm(p/2, lower.tail = F))
    
    # Linear mixed model with TIME:TRT interaction
  } else if (class(x)[1] == "lmerMod"){
    e <- unname(fixef(x)["TIME:TRT"])
    p <- drop1(x, ~TIME:TRT, test = "Chisq")[2, 4]
    z <- sign(e) * abs(qnorm(p/2, lower.tail = F))
    
    # Survival analysis: logrank test
  } else if (class(x)[1] == "survdiff"){
    e <- summary(coxph(Surv(STIME2, EVENT2) ~ TRT, data = D$WD.AGG))$coefficients[1, 2]
    p <- x$pvalue
    z <- ifelse(e <= 1, qnorm(p/2, lower.tail = F), qnorm(p/2, lower.tail = T))
    
    # Repeated measures ANOVA with complete case (model_id m4)
  } else if (class(x)[1] == "anova_test" && identical(attr(x, "model_id"), "m4")){
    e <- suppressMessages(as.data.frame(contrast(emmeans(lmer(CFB ~ TRT*VISIT + (1 | ID), data = D$LD[!is.na(D$LD$CFB), ]), ~ TRT | VISIT, at = list(VISIT = "6")), 
                                method = "revpairwise"))$estimate)
    p <- x$`Sphericity Corrections`$`p[GG]`[2] 
    z <- sign(e) * abs(qnorm(p/2, , lower.tail = F))
    
    # Repeated measures ANOVA with LOCF (model_id m7)
  } else if (class(x)[1] == "anova_test" && identical(attr(x, "model_id"), "m7")){
    e <- suppressMessages(as.data.frame(contrast(emmeans(lmer(CFB ~ TRT*VISIT + (1 | ID), data = D$LOCF2[!is.na(D$LOCF2$CFB), ]), ~ TRT | VISIT, at = list(VISIT = "6")), 
                                                 method = "revpairwise"))$estimate)
    p <- x$`Sphericity Corrections`$`p[GG]`[2] 
    z <- sign(e) * abs(qnorm(p/2, , lower.tail = F))
    
    # Non-parametric Wilcoxon/Mann-Whitney test
  } else if (class(x)[1] == "ScalarIndependenceTest"){
    e <- statistic(x)*-1
    p <- coin::pvalue(x)
    z <- sign(e) * abs(qnorm(p/2, lower.tail = F))
    
    # Logistic regression (binary outcome)
  } else if (class(x)[1] == "glm"){
    e <- exp(unname(coefficients(x)["TRT"]))
    p <- drop1(x, ~TRT, test = "Chisq")[2, 5]
    z <- ifelse(e > 1, qnorm(p/2, lower.tail = F), qnorm(p/2, lower.tail = T)) 
    
  }
  
  return(list("effect estimate" = e, "p-value" = p, "z-score" = z))
}


# ::::::::::::::::::::::::::::::::::: MODELS ::::::::::::::::::::::::::::::::::#
# Fit the models (same models as in "5.models.JM.R")

# For each model, we provide the following information:
# Outcome; analysis approach; missing data strategy death-events; missing data strategy
# non-death events; baseline adjustment 

# Furthermore, it is possible to use the est() function to get the effect estimate, 
# p-value and z-score. This is done for m1, but can be done for all models.

# Change from baseline; linear regression; LOCF; LOCF; no baseline adjustment
m1 <- lm(CFB_6 ~ TRT, data = D$LOCF)
est(m1)

# Time-to-event; logrank test; event; censored; no baseline adjustment
m2 <- survdiff(Surv(STIME2, EVENT2) ~ TRT, data = D$WD.AGG)


# Change from baseline; MMRM; direct modeling; direct modeling; no baseline adjustment
m3 <- mmrm(CFB ~ TRT + VISIT + TRT:VISIT + us(VISIT|ID), data = D$LD[!D$LD$VISIT == 0, ])


# Change from baseline; repeated measures ANOVA; complete case; complete case; no baseline adjustment
m4 <- anova_test(data = D$LD[!is.na(D$LD$CFB), ], dv = CFB, wid = ID, within = VISIT, between = TRT, type = 3)
attr(m4, "model_id") <- "m4"

# Progression rate; mixed effects model; direct modeling; direct modeling; no baseline adjustment
m5 <- lmer(TOT ~ TIME + TRT + TRT:TIME + (TIME|ID), data = D$LD, control = lmerControl(optimizer = "nlminbwrap"))


# Rank; joint rank analysis; joint rank analysis; comparison at last common visit; baseline adjustment
m6 <- lm(RANK ~ TRT + TOT_0, data = D$WD.JR)


# Change from baseline; repeated measures ANOVA; LOCF; LOCF; no baseline adjustment
m7 <- anova_test(data = D$LOCF2, dv = CFB, wid = ID, within = VISIT, between = TRT, type = 3)
attr(m7, "model_id") <- "m7"

# Total score; linear regression; complete case; complete case; no baseline adjustment
m8 <- lm(TOT_6 ~ TRT, data = D$CC)


# Rank; joint rank analysis; joint rank analysis; slope; baseline adjustment
m9 <- lm(RANK.genge ~ TRT + TOT_0, data = D$WD.JR)


# Change from baseline; non-parametric test; zero imputation; zero imputation; no baseline adjustment
m10 <- coin::wilcox_test(CFB_6 ~ as.factor(TRT), data = D$ZERO, distribution = "exact")


# Change from baseline; non-parametric test; zero imputation; LOCF; no baseline adjustment
m11 <- coin::wilcox_test(CFB_6 ~ as.factor(TRT), data = D$WD.KAJI, distribution = "exact")


# Change from baseline; linear regression; zero imputation; worst score imputation; no baseline adjustment
m12 <- lm(CFB_6 ~ TRT, data = D$WD.KAUF)


# Change from baseline; mixed effects model; direct modeling; direct modeling; no baseline adjustment (interaction) 
# (looks at adjusted mean difference at the end of the study instead of progression rates)
m13 <- lmer(CFB ~ TIME:TRT + BSLN:TIME + (TIME|ID), data = D$LD[!D$LD$VISIT == 0, ])
attr(m13, "model_id") <- "m13"

# Progression rate; non-parametric test; complete case; complete case; no baseline adjustment
m14 <- coin::wilcox_test(I((TOT_6 - TOT_0)/TIME_6) ~ as.factor(TRT), data = D$CC, distribution = "exact")


# Progression rate; mixed effects model; direct modeling; direct modeling; baseline adjustment
m15 <- lmer(TOT ~ TIME + TRT + TIME:TRT + BSLN + (TIME|ID), data = D$LD, control = lmerControl(optimizer = "nlminbwrap"))


# Sum score; joint rank analysis; joint rank analysis; comparison at last common visit; baseline adjustment
m16 <- lm(SUM ~ TRT + TOT_0, data = D$WD.JR)


# Change from baseline; non-parametric test; complete case; complete case; no baseline adjustment
m17 <- coin::wilcox_test(CFB_6 ~ as.factor(TRT), data = D$CC, distribution = "exact")


# Change from baseline; linear regression; zero score imputation; LOCF; baseline adjustment
m18 <- lm(CFB_6 ~ TRT + TOT_0, data = D$LOCF.MORA)


# Change from baseline; linear regression; complete case; complete case; no baseline adjustment
m19 <- lm(CFB_6 ~ TRT, data = D$CC)


# Progression rate; random-intercepts model; direct modeling; direct modeling; no baseline adjustment
m20 <- lmer(TOT ~ TRT + TIME + (1|ID), data = D$LD)


# Progression rate; random-intercepts model; direct modeling; direct modeling; baseline adjustment
m21 <- lmer(TOT ~ TIME + TRT + TIME:TRT + BSLN + (1|ID), data = D$LD)


# Change from baseline; MMRM; direct modeling; direct modeling; baseline adjustment
m22 <- mmrm(CFB ~ TRT + VISIT + TRT:VISIT + BSLN + BSLN:VISIT + us(VISIT|ID), data = D$LD[!D$LD$VISIT == 0, ])


# Change from baseline; MMRM; direct modeling; direct modeling; baseline adjustment
m23 <- mmrm(CFB ~ TRT:VISIT + VISIT + BSLN + us(VISIT|ID), data = D$LD[!D$LD$VISIT == 0, ])


# Progression rate; MMRM; direct modeling; direct modeling; no baseline adjustment
D$LD$VISIT2 <- c(0, 1, 2, 3, 4, 5, 6)[match(D$LD$VISIT, c("0","1","2","3", "4", "5", "6"))]
m24 <- mmrm(TOT ~ VISIT2 + TRT:VISIT2 + us(VISIT|ID), data = D$LD, reml = F)
attr(m24, "model_id") <- "m24"

# Change from baseline; MMRM; direct modeling; direct modeling; baseline adjustment
m25 <- mmrm(CFB ~ TRT:VISIT + VISIT + BSLN + us(VISIT|ID), data = D$LOCF2[!D$LOCF2$VISIT == 0, ])


# Change from baseline; linear regression; complete case; complete case; baseline adjustment
m26 <- lm(CFB_6 ~ TRT + I(TOT_0 >= 37), data = D$CC)


# Progression rate; two-step linear regression; available case analysis; available case analysis; no baseline adjustment
m27 <- lmer(SLP ~ TRT + (1|ID), data = D$LOCF.BERRY)


# Proportion responders; logistic regression; non-responder; LOCF; baseline adjustment
m28 <- glm(RESP ~ TRT + BSLN, family = binomial(link = "logit"), data = D$RESP.CUD)


# Proportion responders; logistic regression; LOCF; LOCF; no baseline adjustment
m29 <- glm(RESP ~ TRT, family = binomial(link = "logit"), data = D$RESP.ELIA)


# Rank; joint rank analysis; joint rank analysis; multiple imputation; baseline adjustment
m30 <- lapply(1:5, function(i){
  lm(RANK ~ TRT + TOT_0, data = D$WD.MILL[D$WD.MILL$.imp == i, ])
})
attr(m30, "model_id") <- "m30"


# ::::::::::::::::::::::::::::::::::: FIGURE ::::::::::::::::::::::::::::::::::#
# collect z-scores
Z <- sapply(1:30, function(i) est(get(paste0("m", i)))$`z-score`)

# put everything together
results <- data.frame(model = 1:30, zscore = Z)

# create figure
par(mar = c(5, 7, 4, 2))
results <- results[order(results$zscore), ]
plot(NULL, xlab = "Z-score", ylab = "", axes = F, xlim = c(-3, 3), ylim = c(1, 30), frame.plot = T)
rect(xleft = -3.5, ybottom = -0.1, xright = -1.96, ytop = 31.1, col = "grey90", border = "grey90")
rect(xleft = 1.96, ybottom = -0.1, xright = 3.5, ytop = 31.1, col = "grey90", border = "grey90")
rect(xleft = results$zscore, ybottom = 1:nrow(results) - 0.08, xright = 0, ytop = 1:nrow(results) + 0.08, col = "grey50", border = "grey50")
points(results$zscore, 1:nrow(results), pch = 19)
axis(1, at = seq(-3, 3, by = 1))
segments(x0 = 0, y0 = -0.1, x1 = 0, y1 = 31.1, lty = 2)

text(x = -3.4, y = 1:30, paste0("Model ", results$model), adj = 1, xpd = T)

