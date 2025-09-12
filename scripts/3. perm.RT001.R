set.seed(3239480)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
suppressWarnings(source("2. studyadj.R")) 

# libraries
library(mmrm)
library(lme4)
library(survival)
library(parallel)
library(rstatix)
library(emmeans)
library(coin)
library(mice)
library(pbmcapply)


TE <- 0
TRT.cenrate <- 1
s <- 10000

X <- pbmclapply(1:s, function(i){
  
  # - ################################ DATA ################################ - #
  D <- D.ADJ.RT001(TE = TE, D = "RT001", TRT.cenrate = TRT.cenrate)
  
  D1 <- list(WD = D$WD, LOCF = D$LOCF, ZERO = D$ZERO, COMP = D$COMP,
             WD.AGG = D$WD.AGG, WD.JR = D$WD.JR, WD.KAUF = D$WD.KAUF, WD.MILL = D$WD.MILL)
  
  D2 <- list(LD = D$LD, LD2 = D$LD2, LOCF2 = D$LOCF2)

  # - ############################### MODELS ############################### - #
  # 1. -- Abe - 2014 - edaravone - PMID: 25286015 ------------------------------
  m1 <- lm(V3.CFB ~ TRT + DELFRS + ONSET + RILUSE, data = D1$LOCF)
  p1 <- drop1(m1, ~ TRT, test = "Chisq")[2, 5]
  
  # 2. -- Abe - 2017 - edaravone - PMID: 28522181 (larger 2017 trial) ----------
  m2 <- lm(V3.CFB ~ TRT + DELFRS + I(AGE >= 65), data = D1$LOCF)
  p2 <- drop1(m2, ~ TRT, test = "Chisq")[2, 5]
  
  # 3. -- Abe - 2017 - edaravone - PMID: 28872915 (smaller 2017 trial) ---------
  m3 <- lm(V3.CFB ~ TRT + DELFRS, data = D1$LOCF)
  p3 <- drop1(m3, ~ TRT, test = "Chisq")[2, 5]
  
  # 4. -- Aggarwal - 2010 - lithium - PMID: 20363190 ---------------------------
  m4 <- survdiff(Surv(STIME, EVENT) ~ TRT, data = D1$WD.AGG)
  p4 <- m4$pvalue
  
  # 5. -- Aizawa - 2022 - perampanel - PMID: 34191081 --------------------------
  m5 <- mmrm(CFB ~ TRT + VISIT + TRT:VISIT + us(VISIT|ID), data = D2$LD[!D2$LD$VISIT == 0, ])
  p5 <- unname(unlist(summary(emmeans(m5, specs = pairwise ~ TRT, at = list(VISIT = c("3")))$contrast)[6]))
  
  # 6. -- Amirzagar - 2015 - G-CSF - PMID: 25851895 ----------------------------
  m6 <- mmrm(CFB ~ TRT + VISIT + TRT:VISIT + us(VISIT|ID), data = D2$LD[!D2$LD$VISIT == 0, ])
  p6 <- unname(unlist(summary(emmeans(m6, specs = pairwise ~ TRT, at = list(VISIT = c("3")))$contrast)[6]))
  
  # 7. -- Berry - 2019 - MSC-NTF - PMID: 31740545 ------------------------------
  # not included (no pre-treatment slope)
  
  # 8. -- Boll - 2022 - valproate/lithium - PMID: 36049647 ---------------------
  m8 <- anova_test(data = D2$LD, dv = CFB, wid = ID, within = VISIT, between = TRT, type = 3)
  p8 <- m8$ANOVA$p[3]
  p8.1 <- m8$`Sphericity Corrections`$`p[GG]`[2] # Greenhouse-Geiser sphericity correction
  p8.2 <- m8$`Sphericity Corrections`$`p[HF]`[2] # Huynh-Feldt sphericity correction
  
  # 9. -- Cudkowicz - 2011 - dexpramipexole - PMID: 22101764 -------------------
  m9 <- lmer(TOT ~ TIME + TRT + TRT:TIME + (TIME|ID), data = D2$LD)
  p9 <- drop1(m9, ~ TIME:TRT, test = "Chisq")[2, 4]
  
  # 10. -- Cudkowicz - 2013 - dexpramipexole - PMID: 24067398 ------------------
  m10 <- lm(RANK ~ TRT + V0.TOT + DISDUR + ONSET + RILUSE, data = D1$WD.JR)
  p10 <- drop1(m10, ~ TRT, test = "Chisq")[2, 5]
  
  # 11. -- Cudkowicz - 2014 - ceftriaxone - PMID: 25297012 ---------------------
  m11 <- lmer(TOT ~ TIME + TRT + TRT:TIME + (TIME|ID), data = D2$LD)
  p11 <- drop1(m11, ~ TIME:TRT, test = "Chisq")[2, 4]
  
  # 12. -- Cudkowicz - 2022 - MSC-NTF - PMID: 34890069 -------------------------
  # not included (no pre-treatment slope)
  
  # 13. -- De Carvalho - 2010 - memantine - PMID: 20565333 ---------------------
  m13 <- anova_test(data = D2$LOCF2, dv = CFB, wid = ID, within = VISIT, between = TRT, type = 3)
  p13 <- m13$ANOVA$p[3]
  p13.1 <- m13$`Sphericity Corrections`$`p[GG]`[2] # Greenhouse-Geiser sphericity correction
  p13.2 <- m13$`Sphericity Corrections`$`p[HF]`[2] # Huynh-Feldt sphericity correction
  
  # 14. -- De la Rubia - 2019 - EH301 - PMID: 30668199 -------------------------
  m14 <- lm(V3.TOT ~ TRT, data = D1$COMP)
  p14 <- drop1(m14, test = "Chisq")[2, 5]
  
  # 15. -- Elia - 2015 - tauroursodeoxycholic acid - PMID: 25664595 ------------
  # not included (no pre-treatment slope)
  
  # 16. -- Genge - 2023 - ravulizumab - PMID: 37695623 -------------------------
  m16 <- lm(RANK.genge ~ TRT + AGE + SEX + V0.TOT + SVCP + DISDUR + ONSET + RILUSE, 
            data = D1$WD.JR)
  p16 <- drop1(m16, ~ TRT, test = "Chisq")[2, 5]
  
  # 17. -- Gordon - 2007 - minocycline - PMID: 17980667 ------------------------
  m17 <- lmer(TOT ~ TIME + TRT + TIME:TRT + (TIME|ID), data = D2$LD)
  p17 <- drop1(m17, ~ TIME:TRT, test = "Chisq")[2, 4]
  
  # 18. -- Juntas-Morales - 2020 - MD1003 - PMID: 32140672 ---------------------
  m18 <- coin::wilcox_test(V3.CFB ~ as.factor(TRT), data = D1$ZERO, distribution = "exact")
  p18 <- coin::pvalue(m18)
  
  # 19. -- Kaji - 2019 - methylcobalamin - PMID: 30636701 ----------------------
  m19 <- coin::wilcox_test(V3.CFB ~ as.factor(TRT), data = D1$LOCF, distribution = "exact")
  p19 <- coin::pvalue(m19)
  
  # 20. -- Kaufmann - 2009 - COQ10 - PMID: 19743457 ----------------------------
  m20 <- lm(V3.CFB ~ TRT, data = D1$WD.KAUF)
  p20 <- drop1(m20, test = "Chisq")[2, 5]
  
  # 21. -- Kim - 2023 - mecasin - PMID: 37257710 -------------------------------
  m21 <- lmer(CFB ~ TIME:TRT + BSLN:TIME + DISDUR:TIME + (TIME|ID), data = D2$LD[!D2$LD$VISIT == 0, ])
  p21 <- summary(contrast(emmeans(m21, ~ TRT, at = list(TIME = 6)), method = "pairwise"))$p
  
  # 22. -- Levine - 2012 - pioglitazone HCI-tretinoin - PMID: 22830016 ---------
  m22 <- coin::wilcox_test(I((V3.TOT - V0.TOT)/V3.TIME) ~ as.factor(TRT), data = D1$COMP, distribution = "exact")
  p22 <- coin::pvalue(m22)
  
  # 23. -- Meininger - 2009 - glatiramer - PMID: 19922128 ----------------------
  m23 <- lmer(TOT ~ TIME + TRT + TIME:TRT + AGE + ONSET + DISDUR + BSLN + RILUSE + 
                CNTRY + (TIME|ID), data = D2$LD, control = lmerControl(optimizer = "nlminbwrap"))
  p23 <- drop1(m23, ~ TIME:TRT, test = "Chisq")[2, 4]
  
  # 24. -- Meininger - 2017 - ozanezumab - PMID: 28139349 ----------------------
  m24 <- lm(SUM ~ TRT + V0.TOT + RILUSE + CNTRY, data = D1$WD.JR)
  p24 <- drop1(m24, ~ TRT, test = "Chisq")[2, 5]
  
  # 25. -- Miller - 2007 - TCH346 - PMID: 17709710 -----------------------------
  m25 <- lmer(TOT ~ TIME + TRT + TIME:TRT + (TIME|ID), data = D2$LD)
  p25 <- drop1(m25, ~ TIME:TRT, test = "Chisq")[2, 4]
  
  # 26. -- Miller - 2015 - NP001 - PMID: 25884010 ------------------------------
  m26 <- lmer(TOT ~ TIME + TIME:TRT + AGE + SEX +  DISDUR + RILUSE + ONSET + 
                BSLN + SVCL + (TIME|ID), data = D2$LD,
              control = lmerControl(optimizer = "nlminbwrap"))
  p26 <- drop1(m26, ~ TIME:TRT, test = "Chisq")[2, 4]
  
  # 27. -- Miller - 2022 - tofersen - PMID: 36129998 ---------------------------
  m27 <- lapply(1:5, function(i){
    lm(RANK ~ TRT + DISDUR + V0.TOT + RILUSE, data = D1$WD.MILL[D1$WD.MILL$.imp == i, ])
  })
  
  estimates <- pool(m27)
  p27 <- summary(estimates)[2, 6]
  
  # 28. -- Miller - 2022 - NP001 - PMID: 35098554 ------------------------------
  m28 <- coin::wilcox_test(V3.CFB ~ as.factor(TRT), data = D1$COMP, distribution = "exact")
  p28 <- coin::pvalue(m28)
  
  # 29. -- Mora - 2020 - masitinib - PMID: 31280619 ----------------------------
  m29 <- lm(V3.CFB ~ TRT + ONSET + V0.TOT + AGE + CNTRY + DELFRS, data = D1$LOCF)
  p29 <- drop1(m29, ~ TRT, test = "Chisq")[2, 5]
  
  # 30. -- Morimoto - 2023 - ropinirole - PMID: 37267913 -----------------------
  m30 <- mmrm(CFB ~ TRT + VISIT + TRT:VISIT + us(VISIT|ID), data = D2$LD[!D2$LD$VISIT == 0, ])
  p30 <- unname(unlist(summary(emmeans(m30, specs = pairwise ~ TRT, at = list(VISIT = c("3")))$contrast)[6]))
  
  # 31. -- Nagata - 2016 - bromocriptine mesylate - PMID: 26910108 -------------
  m31 <- lm(V3.CFB ~ TRT + AGE, data = D1$LOCF)
  p31 <- drop1(m31, ~ TRT, test = "Chisq")[2, 5]
  
  # 32. -- Nefussy - 2010 - G-CSF - PMID: 19449238 -----------------------------
  m32 <- lm(V3.CFB ~ TRT, data = D1$COMP)
  p32 <- drop1(m32, test = "Chisq")[2, 5]
  
  # 33. -- Oh - 2018 - BM-MSC - PMID: 30048006 ---------------------------------
  m33 <- lm(V3.CFB ~ TRT, data = D1$COMP)
  p33 <- drop1(m33, test = "Chisq")[2, 5]
  
  # 34. -- Oki - 2022 - methylcobalamin - PMID: 35532908 -----------------------
  m34 <- mmrm(CFB ~ TRT + VISIT + TRT:VISIT + CNTRY + ONSET + DISDUR + SVCP + 
                us(VISIT|ID), data = D2$LD[!D2$LD$VISIT == 0, ])
  p34 <- unname(unlist(summary(emmeans(m34, specs = pairwise ~ TRT, at = list(VISIT = c("3")))$contrast)[6]))
  
  # 35. -- Paganoni - 2020 - sodium phenylbutyrate-taurursodiol - PMID: 32877582
  m35 <- lmer(TOT ~ TIME + TIME:TRT + AGE:TIME + DELFRS:TIME + (TIME|ID), 
              data = D2$LD)
  p35 <- drop1(m35, ~ TIME:TRT, test = "Chisq")[2, 4]
  
  # 36. -- Park - 2015 - flecainide - PMID: 26844270 ---------------------------
  m36 <- lmer(TOT ~ TRT + TIME + DELFRS + (1|ID), data = D2$LD)
  p36 <- drop1(m36, ~ TRT, test = "Chisq")[2, 4]
  
  # 37. -- Scelsa - 2005 - indinavir - PMID: 15824372 --------------------------
  m37 <- lmer(TOT ~ TIME + TRT + TIME:TRT + BSLN + (1|ID), data = D2$LD)
  p37 <- drop1(m37, ~ TIME:TRT, test = "Chisq")[2, 4]
  
  # 38. -- Schröder - 2022 - herbal - PMID: 36324375 ---------------------------
  m38 <- lm(V3.CFB ~ TRT, data = D1$LOCF)
  p38 <- drop1(m38, test = "Chisq")[2, 5]
  
  # 39. -- Shefner - 2016 - tirasemtiv - PMID: 26982815 ------------------------
  m39 <- mmrm(CFB ~ TRT + VISIT + TRT:VISIT + BSLN + BSLN:VISIT + RILUSE + CNTRY +
                us(VISIT|ID), data = D2$LD[!D2$LD$VISIT == 0, ])
  p39 <- unname(unlist(summary(emmeans(m39, specs = pairwise ~ TRT, at = list(VISIT = c("3")))$contrast)[6]))
  
  # 40. -- Shibuya - 2015 - mexiletine - PMID: 25960085 ------------------------
  m40 <- lm(V3.CFB ~ TRT + I(V0.TOT >= 37) + I(AGE <= 65), data = D1$COMP)
  p40 <- drop1(m40, ~ TRT, test = "Chisq")[2, 5]
  
  # 41. -- Statland - 2019 - rasagiline - PMID: 30192007 -----------------------
  m41 <- lmer(TOT ~ TIME + TIME:TRT + DISDUR + RILUSE + (TIME|ID), data = D2$LD)
  p41 <- drop1(m41, ~ TIME:TRT, test = "Chisq")[2, 4]
  
  # 42. -- Vucic - 2021 - dimethyl fumerate - PMID: 34477330 -------------------
  m42 <- mmrm(CFB ~ VISIT + TRT:VISIT + BSLN + us(VISIT|ID), data = D2$LOCF2[!D2$LOCF2$VISIT == 0, ])
  p42 <- unname(unlist(summary(emmeans(m42, specs = pairwise ~ TRT, at = list(VISIT = c("3")))$contrast)[6]))
  
  # 43. -- Walk - 2023 - inosine - PMID: 36840949 ------------------------------
  m43 <- mmrm(TOT ~ VISIT + TRT:VISIT + us(VISIT|ID), data = D2$LD, reml = F)
  m43.1 <- mmrm(TOT ~ VISIT + us(VISIT|ID), data = D2$LD, reml = F)
  
  p43 <- -2*(logLik(m43.1) - logLik(m43))
  p43 <- pchisq(p43, df = 4, lower.tail = F)
  
  # 44. -- Weemering - 2023 - RT001 - PMID: 37550954 ---------------------------
  m44 <- mmrm(CFB ~ VISIT + TRT:VISIT + BSLN + CNTRY + TRICALS:VISIT + us(VISIT|ID), 
              data = D2$LD[!D2$LD$VISIT == 0, ])
  p44 <- unname(unlist(summary(emmeans(m44, specs = pairwise ~ TRT, at = list(VISIT = c("3")))$contrast)[6]))
  
  # 45. -- Weiss - 2016 - mexiletine - PMID: 26911633 --------------------------
  m45 <- lmer(CFB ~ TIME + TIME:TRT + ONSET + DISDUR + RILUSE + ONSET:TIME + 
                DISDUR:TIME + RILUSE:TIME + (TIME|ID), data = D2$LD[!D2$LD$VISIT == 0, ])
  p45 <- drop1(m45, ~ TIME:TRT, test = "Chisq")[2, 4]
  
  return(mget(paste0("p", c(1:6, 8:11, 13:14, 16:45))))
}, mc.cores = 3) 

