library(haven)
library(tidyr)
library(zoo)

setwd("/Users/dweemeri/surfdrive - Weemering, D.N. (Daphne)@surfdrive.surf.nl/data/ceftriaxone/Datasets")
A <- haven::read_sas("alsf.sas7bdat")

# TIME
A <- merge(A, A[A$Visit_Name == "Screening", c("study_id", "age_at_DATE_PERFORMED")],
           by = "study_id")
names(A)[names(A) == "age_at_DATE_PERFORMED.y"] <- "age_at_SCREENING"
A$TIME <- ifelse(is.na(A$age_at_DATE_PERFORMED.x) == F, (A$age_at_DATE_PERFORMED - A$age_at_SCREENING)*12,
                 (A$age_at_VISIT_DATE - A$age_at_SCREENING)*12)
A <- A[order(A$study_id, A$TIME), ]

# Keep only measurements up until one year of treatment to match phase 3 trial
A <- A[A$TIME <= 12.5, ]

# VISIT
A$VISIT <- ifelse(grepl("week|Week", A$Visit_Name), gsub("[^0-9]", "", A$Visit_Name), A$Visit_Name)
A$VISIT <- ifelse(A$VISIT %in% as.character(c(4, seq(16, 152, by = 8))),
                  paste0("V", match(A$VISIT, as.character(c(4, seq(16, 152, by = 8))))),
                  ifelse(A$Visit_Name == "Screening", "V0", A$Visit_Name))

# remove other visits for mmrm analyses
A[which(A$VISIT == "ALSFRS during week 36"), ]$VISIT <- NA
A[which(A$VISIT == "ALSFRS"), ]$VISIT <- NA
A[which(A$VISIT == "Unscheduled ALS-FRS"), ]$VISIT <- NA
A[which(A$VISIT == "Unscheduled ALSFRS"), ]$VISIT <- "V4"
A[which(A$VISIT == "ALSFRS 8 WK POST STUDY DISCONTINUATION"), ]$VISIT <- NA
A[which(A$VISIT == "End of treatment visit"), ]$VISIT <- "V2"
A[which(A$VISIT == "ITT-Wk 24"), ]$VISIT <- "V3"
A[which(A$VISIT == "Evaluations at Week 52 9-7-2011 (missed Week 48 visit)"), ]$VISIT <- "V6"
A[which(A$VISIT == "Intent to Treat week 28"), ]$VISIT <- "V3"
A[which(A$VISIT == "Week 44 Intent to treat"), ]$VISIT <- "V6"
A[which(A$VISIT == "ITT (Wk 16)"), ]$VISIT <- "V2"
A[which(A$VISIT == "ITT Outcomes Testing 06-03-10"), ]$VISIT <- "V3"
A[which(A$VISIT == "ITT Outcomes Testing 09-02-2010"), ]$VISIT <- "V4"
A[which(A$VISIT == "ITT Outcomes Testing 12-02-2010"), ]$VISIT <- "V6"
A[which(A$VISIT == "ITT Wk 36"), ]$VISIT <- "V5"
A[which(A$VISIT == "ITT-9-3-2010"), ]$VISIT <- "V6"
A[which(A$VISIT == "Outcomes Visit 1"), ]$VISIT <- "V3"
A[which(A$VISIT == "Outcomes visit 2"), ]$VISIT <- NA
A[which(A$VISIT == "Week 36 Outcomes Visit"), ]$VISIT <- "V4"
A[which(A$VISIT == "Repeat screening assessments"), ]$VISIT <- NA
A[which(A$VISIT == "UltraSound at Final Visit"), ]$VISIT <- NA
A[which(A$VISIT == "Unscheduled 4 - ALSFRS-R 4/12/2011"), ]$VISIT <- "V4"
A[which(A$VISIT == "Unscheduled - A:SFRS 09/21/2011"), ]$VISIT <- "V6"
A[which(A$VISIT == "Unscheduled 1-ALSFRS 03/12/2012"), ]$VISIT <- NA
A[which(A$VISIT == "Unscheduled - ALSFRS 05/21/2012"), ]$VISIT <- "V3"
A[which(A$VISIT == "Unscheduled ALSFRS 7/13/2012"), ]$VISIT <- "V4"
A[which(A$VISIT == "Unscheduled 8/21/2012"), ]$VISIT <- "V5"
A[which(A$VISIT == "Unscheduled - Final Study Visit"), ]$VISIT <- "V2"
A[which(A$VISIT == "Unscheduled 1 - ALSFRS 2/23/2012"), ]$VISIT <- "V5"
A[which(A$VISIT == "Unscheduled 1 - Vital Capacity"), ]$VISIT <- NA
A[which(A$VISIT == "Unscheduled 1 12/15/2010"), ]$VISIT <- NA
A[which(A$VISIT == "Unscheduled 14"), ]$VISIT <- "V4"
A[which(A$VISIT == "Unscheduled 13"), ]$VISIT <- "V5"
A[which(A$VISIT == "Unscheduled 2 - Abdominal Ultrasound"), ]$VISIT <- NA
A[which(A$VISIT == "Unscheduled 2, week 40 Vital Capacity"), ]$VISIT <- NA
A[which(A$VISIT == "Unscheduled 3 lab work only"), ]$VISIT <- NA
A[which(A$VISIT == "Unscheduled 4"), ]$VISIT <- NA
A[which(A$VISIT == "Unscheduled abd ultrasound"), ]$VISIT <- NA
A[which(A$VISIT == "Unscheduled Final Study Visit ALSFRS-R 8/21/2012"), ]$VISIT <- "V6"
A[which(A$VISIT == "Unscheduled re-test chemistry panel"), ]$VISIT <- NA
A[which(A$VISIT == "Unscheduled screening labs"), ]$VISIT <- NA
A[which(A$VISIT == "Unscheduled screening outside labs"), ]$VISIT <- NA
A[which(A$VISIT == "Unscheduled Ultrasound"), ]$VISIT <- NA
A[which(A$VISIT == "Unscheduled visit Sept. ALS clinic 9/7/07"), ]$VISIT <- NA
A[which(A$VISIT == "Unscheduled Week 28 ALSRS-R"), ]$VISIT <- NA
A[which(A$VISIT == "Unscheduled week 36"), ]$VISIT <- NA
A[which(A$VISIT == "V7" & A$study_id == "79"), ]$VISIT <- "V6"
A[which(A$VISIT == "V7" & A$study_id == "156"), ]$VISIT <- NA
A[which(A$VISIT == "V7" & A$study_id == "243"), ]$VISIT <- "V6"
A[which(A$VISIT == "Week 20 additional tests"), ]$VISIT <- NA
A[which(A$VISIT == "week 20 unscheduled ultrasound"), ]$VISIT <- NA
A[which(A$VISIT == "Week 36 abdominal ultrasound"), ]$VISIT <- NA
A[which(A$VISIT == "Week 40 - ALSFRS 7/25/2011"), ]$VISIT <- "V5"
A[which(A$VISIT == "Week 44 (ITT Outcomes)"), ]$VISIT <- "V6"
A[which(A$VISIT == "week 44 unscheduled ultrasound"), ]$VISIT <- NA
A[which(A$VISIT == "wk 32 assessments done at wk 36 visit"), ]$VISIT <- NA
A[which(A$VISIT == "V5" & A$study_id == "94"), ]$VISIT <- "V6"
A[which(A$VISIT == "Wk 36 (ALSFRS/ALSQoLQ/FVC/HHD)"), ]$VISIT <- "V5"
A[which(A$VISIT == "Wk 40 ALSFRS-R"), ]$VISIT <- NA

A$VISIT <- ifelse(A$VISIT == "Final Study Visit" & A$TIME > 3 & A$TIME <= 5, "V2", A$VISIT)
A$VISIT <- ifelse(A$VISIT == "Final Study Visit" & A$TIME > 5 & A$TIME <= 7, "V3", A$VISIT)
A$VISIT <- ifelse(A$VISIT == "Final Study Visit" & A$TIME > 7 & A$TIME <= 9, "V4", A$VISIT)
A$VISIT <- ifelse(A$VISIT == "Final Study Visit" & A$TIME > 9 & A$TIME <= 11, "V5", A$VISIT)
A$VISIT <- ifelse(A$VISIT == "Final Study Visit" & A$TIME > 11, "V6", A$VISIT)
A$VISIT <- ifelse(A$VISIT == "Final Study Visit" & A$TIME > 2.2 & A$TIME <= 2.5, "V2", A$VISIT)

A[which((A$VISIT == "Unscheduled 1" | A$VISIT == "Unscheduled 2" | A$VISIT == "Unscheduled 3") & is.na(A$als_total) == T), ]$VISIT <- NA
A$VISIT <- ifelse((A$VISIT == "Unscheduled 1" | A$VISIT == "Unscheduled 2" | A$VISIT == "Unscheduled 3") & A$TIME <= 3  & !ave(A$VISIT == "V1", A$study_id, FUN = any),
                    "V1", A$VISIT)
A$VISIT <- ifelse((A$VISIT == "Unscheduled 1" | A$VISIT == "Unscheduled 2" | A$VISIT == "Unscheduled 3") & A$TIME > 3 & A$TIME <= 5 & !ave(A$VISIT == "V2", A$study_id, FUN = any),
                    "V2", A$VISIT)
A$VISIT <- ifelse((A$VISIT == "Unscheduled 1" | A$VISIT == "Unscheduled 2" | A$VISIT == "Unscheduled 3") & A$TIME > 5 & A$TIME <= 7 & !ave(A$VISIT == "V3", A$study_id, FUN = any),
                    "V3", A$VISIT)
A$VISIT <- ifelse((A$VISIT == "Unscheduled 1" | A$VISIT == "Unscheduled 2" | A$VISIT == "Unscheduled 3") & A$TIME > 7 & A$TIME <= 9 & !ave(A$VISIT == "V4", A$study_id, FUN = any),
                    "V4", A$VISIT)
A$VISIT <- ifelse((A$VISIT == "Unscheduled 1" | A$VISIT == "Unscheduled 2" | A$VISIT == "Unscheduled 3") & A$TIME > 9 & A$TIME <= 11 & !ave(A$VISIT == "V5", A$study_id, FUN = any),
                    "V5", A$VISIT)
A$VISIT <- ifelse((A$VISIT == "Unscheduled 1" | A$VISIT == "Unscheduled 2" | A$VISIT == "Unscheduled 3") & A$TIME > 11 & !ave(A$VISIT == "V6", A$study_id, FUN = any),
                    "V6", A$VISIT)

A[which(A$VISIT == "Unscheduled 1"), ]$VISIT <- NA
A[which(A$VISIT == "Unscheduled 2"), ]$VISIT <- NA
A[which(A$VISIT == "Unscheduled 3"), ]$VISIT <- NA

# double visits
A <- A[!(A$study_id == "59" & A$VISIT == "V2" & A$TIME > 4), ]
A[(A$study_id == "62" & A$VISIT == "V5" & A$TIME > 10), ]$VISIT <- "V6"
A <- A[!(A$study_id == "66" & A$VISIT == "V2" & A$TIME < 3), ]
A <- A[-which(A$study_id == unique(A[(A$study_id == "120" & A$VISIT == "V5"), ])[1, ]$study_id &
                A$VISIT == unique(A[(A$study_id == "120" & A$VISIT == "V5"), ])[1, ]$VISIT &
                A$TIME == unique(A[(A$study_id == "120" & A$VISIT == "V5"), ])[1, ]$TIME)[2], ]
A[(A$study_id == "172" & A$VISIT == "V5" & A$TIME > 10), ]$VISIT <- "V6"
A <- A[!(A$study_id == "204" & A$VISIT == "V6"), ]
A <- A[!(A$study_id == "204" & A$VISIT == "V1"), ]
A[(A$study_id == "204" & A$VISIT == "V5" & A$TIME > 10), ]$VISIT <- "V6"
A <- A[!(A$study_id == "205" & A$VISIT == "V6" & A$TIME < 12), ]
A <- A[!(A$study_id == "218" & A$VISIT == "V4" & A$TIME < 8), ]
A <- A[!(A$study_id == "296" & A$VISIT == "V2" & A$TIME < 3), ]
A <- A[!(A$study_id == "303" & A$VISIT == "V2" & A$TIME > 4), ]
A[(A$study_id == "340" & A$VISIT == "V5" & A$TIME > 10), ]$VISIT <- "V6"
A <- A[!(A$study_id == "390" & A$VISIT == "V4" & A$TIME < 8), ]
A[(A$study_id == "394" & A$VISIT == "V5" & A$TIME > 10), ]$VISIT <- "V6"
A[(A$study_id == "412" & A$VISIT == "V5" & A$TIME > 10), ]$VISIT <- "V6"
A[(A$study_id == "413" & A$VISIT == "V5" & A$TIME > 10), ]$VISIT <- "V6"
A[(A$study_id == "462" & A$VISIT == "V5" & A$TIME > 10), ]$VISIT <- "V6"

A[(A$study_id == "188" & is.na(A$VISIT) & A$TIME > 5 & A$TIME < 7), ]$VISIT <- "V3"
A[(A$study_id == "188" & is.na(A$VISIT) & A$TIME > 10), ]$VISIT <- "V6"
A <- A[!(A$study_id == "188" & is.na(A$VISIT)), ]

A[(A$study_id == "277" & is.na(A$VISIT) & A$TIME > 7 & A$TIME < 8), ]$VISIT <- "V4"
A <- A[!(A$study_id == "277" & is.na(A$VISIT)), ]

A[(A$study_id == "329" & is.na(A$VISIT) & A$TIME > 5.9 & A$TIME < 6.1), ]$VISIT <- "V3"
A[(A$study_id == "329" & A$TIME > 7.1 & A$TIME < 7.3), ]$VISIT <- "V4"
A[(A$study_id == "329" & A$TIME > 9.5 & A$TIME < 9.7), ]$VISIT <- "V5"
A <- A[!(A$study_id == "329" & is.na(A$VISIT)), ]

# TRT
B <- haven::read_sas("treatment.sas7bdat")
B$TRT <- ifelse(B$treatment == "Active", 1, 0)
A <- merge(A, B[, c("study_id", "TRT")], by = "study_id")

# DEMOGRAPHICS
C <- haven::read_sas("demo.sas7bdat")
A <- merge(A, C[, c("study_id", "sex", "racial_categories")], by = "study_id" )
names(A)[match(c("sex", "racial_categories", "age_at_BASELINE"), names(A))] <- c("SEX", "RACE", "AGE")
A$RACE <- ifelse(A$RACE == 0, "WHITE",
                 ifelse(A$RACE == 1, "BLACK",
                        ifelse(A$RACE == 2, "ASIAN",
                               ifelse(A$RACE == 3, "HAWAI",
                                      ifelse(A$RACE == 4, "INDIAN/ALASKA",
                                             ifelse(A$RACE == 5, "MORE", "UNKN"))))))

# RILUZOLE, ONSET, FAMILIAL, DISDUR, DIAGDELAY
E <- haven::read_sas("mhas.sas7bdat")
A <- merge(A, E[, c("study_id", "riluzole", "siteofonset", "als_fam_hist",
                    "age_at_ALS_DATE_SYM", "age_at_DATE_OF_DIAG")], by = "study_id")
A$siteofonset <- ifelse(A$siteofonset == 0, "LIMB",
                        ifelse(A$siteofonset == 1, "BULB", "BOTH"))
A$als_fam_hist <- ifelse(A$als_fam_hist == 2, "FAM",
                         ifelse(A$als_fam_hist == 1, "SPOR", "UNKN"))
A$DISDUR <- (A$AGE - A$age_at_ALS_DATE_SYM)*12
A$DIAGDELAY <- (A$age_at_DATE_OF_DIAG - A$age_at_ALS_DATE_SYM)*12

A$riluzole <- ifelse(A$riluzole == 2, 1, 0)
L <- haven::read_sas("cm.sas7bdat")
L <- L[L$medication == "rilutek" | L$medication == "Rilutek" | L$medication == "riluzole" |
         L$medication == "Riluzole" | L$medication == "RILUZOLE", ]
L <- L[!L$concomitant_medications == 1, ]
L <- L[!duplicated(L$study_id), ]
A <- merge(A, L[, c("study_id", "medication")], by = "study_id", all.x = T)
A$riluzole <- ifelse(is.na(A$riluzole) == T & is.na(A$medication) == F, 1, A$riluzole)
A$riluzole <- ifelse(is.na(A$riluzole) == T & is.na(A$medication) == T, 0, A$riluzole)


# ESCORIAL
G <- haven::read_sas("als.sas7bdat")
A <- merge(A, G[, c("study_id", "criteria")], by = "study_id")
A$criteria <- ifelse(A$criteria == 1, "POS",
                     ifelse(A$criteria == 2, "PROB_LAB",
                            ifelse(A$criteria == 3, "PROB", "DEF")))

# CENSORING
H <- haven::read_sas("mort.sas7bdat")
A <- merge(A, H[, c("study_id", "did_die", "age_at_DATE_DEATH")], by = "study_id", all.x = T)
names(A)[names(A) == "did_die"] <- "CENSOR"
A$CENSOR <- ifelse(A$CENSOR == 2, 1, 0)
A$STIME <- ifelse(A$CENSOR == 1, (A$age_at_DATE_DEATH - A$AGE)*12,
                  ave(A$TIME, A$study_id, FUN = function(x) max(x, na.rm = TRUE)))
A$STIME2 <- ifelse(A$STIME > 12 & A$CENSOR == 1, 12, A$STIME)
A$CENSOR2 <- ifelse(A$STIME > 12 & A$CENSOR == 1, 0, A$CENSOR)
A$STIME <- A$STIME2
A$CENSOR <- A$CENSOR2

I <- haven::read_sas("disp.sas7bdat")
A <- merge(A, I[, c("study_id", "subj_final_disposition", "age_at_DATE_CONSENT_WITHDRAWAL",
                    "age_at_DATE_OF_LAST_VISIT")], by = "study_id", all.x = T)

A$STIME <- ifelse((A$subj_final_disposition == 3 | A$subj_final_disposition == 4) & is.na(A$age_at_DATE_CONSENT_WITHDRAWAL) == F,
                   (A$age_at_DATE_CONSENT_WITHDRAWAL - A$AGE)*12,
                   ifelse((A$subj_final_disposition == 3 | A$subj_final_disposition == 4) & is.na(A$age_at_DATE_OF_LAST_VISIT) == F,
                          (A$age_at_DATE_OF_LAST_VISIT - A$AGE)*12, A$STIME))
A$STIME <- ifelse(A$STIME > 12, 12, A$STIME)
A[A$study_id == "258", ]$STIME <- 12
A[A$study_id == "258", ]$CENSOR <- 0
A[A$study_id == "294", ]$CENSOR <- 0
A[A$study_id == "346", ]$CENSOR <- 0
A[A$study_id == "349", ]$CENSOR <- 0
A[A$study_id == "383", ]$CENSOR <- 0
A[A$study_id == "397", ]$CENSOR <- 0
# ID 298 does not have information on survival, should be removed from survival dataset

# FVC
K <- haven::read_sas("pfvc.sas7bdat")
FVC <- data.frame(study_id = unique(K$study_id),
                  FVCL = apply(K[K$Visit_Name == "Screening", c("vc1_subj", "vc2_subj", "vc3_subj")], 1, mean, na.rm = T),
                  FVCP = apply(K[K$Visit_Name == "Screening", c("vc1_predicted", "vc2_predicted", "vc3_predicted")], 1, mean, na.rm = T))
A <- merge(A, FVC, by = "study_id")

# these IDs have two or more TIME == 0 measurements, remove the TIME == 0 that are not right
A <- A[!(A$study_id == "221" & is.na(A$VISIT)), ]
A <- A[!(A$study_id == "257" & is.na(A$VISIT)), ]
A <- A[!(A$study_id == "29" & A$VISIT == "V1"), ]
A <- A[!(A$study_id == "151" & A$VISIT == "V1"), ]
A <- A[!(A$study_id == "193" & A$VISIT == "V1"), ]

# cleaned data
D <- data.frame(ID = A$study_id, TRT = A$TRT, TIME = A$TIME, VISIT = A$VISIT,
                TOT = A$als_total,
                CFB = A$als_total - tapply(A$als_total[A$TIME == 0], unique(A$study_id[A$TIME == 0]), mean)[as.character(A$study_id)],
                SEX = A$SEX, AGE = A$AGE, RACE = A$RACE, RILUSE = A$riluzole,
                ONSET = A$siteofonset, FAMILIAL = A$als_fam_hist, DISDUR = A$DISDUR,
                DIAGDELAY = A$DIAGDELAY, ESCORIAL = A$criteria, CENSOR = A$CENSOR,
                STIME = A$STIME, FVCL = A$FVCL, FVCP = A$FVCP, AGE_SYMP = A$age_at_ALS_DATE_SYM)
D <- D[order(D$ID, D$TIME), ]

# data formats
LD <- D[!is.na(D$TOT), ]
LD$VISIT <- as.numeric(gsub("V", "", LD$VISIT))
LD[c("ID", "VISIT")] <- lapply(LD[c("ID", "VISIT")], as.factor)

# DELFRS
tmpLD <- LD[LD$VISIT == "0", c("ID", "TOT", "DISDUR")]
tmpLD <- tmpLD[!is.na(tmpLD$ID), ]
tmpLD$DELFRS <- (tmpLD$TOT - 48) / tmpLD$DISDUR
LD <- merge(LD, tmpLD[, c("ID", "DELFRS")], by = "ID", all.x = TRUE)

# BASELINE
tmpLD <- LD[LD$VISIT == "0", c("ID", "TOT")]
tmpLD <- tmpLD[!is.na(tmpLD$ID), ]
names(tmpLD)[names(tmpLD) == "TOT"] <- "BSLN"
LD <- merge(LD, tmpLD[, c("ID", "BSLN")], by = "ID", all.x = T)
LD <- LD[order(LD$ID, LD$TIME), ]

# TRICALS RISK PROFILE
tmpLD <- LD[LD$VISIT == "0", c("ID", "FVCP", "DIAGDELAY", "DELFRS", "AGE_SYMP",
                               "ONSET", "ESCORIAL")]
tmpLD <- tmpLD[!is.na(tmpLD$ID), ]
tmpLD$ONSET <- ifelse(tmpLD$ONSET == "BULB", 1, 0)
tmpLD$ESCORIAL <- ifelse(tmpLD$ESCORIAL == "DEF", 1, 0)

tmpLD$TRICALS <- (0.474 * ((tmpLD$FVCP/100)^-1) + ((tmpLD$FVCP/100)^-0.5)) -
  (2.376 * ((tmpLD$DIAGDELAY/10)^-0.5) + log(tmpLD$DIAGDELAY/10)) -
  (1.839 * ((-tmpLD$DELFRS + 0.1)^-0.5)) - (0.264 * (tmpLD$AGE_SYMP/100)^-2) +
  (0.271 * tmpLD$ONSET) + (0.238 * tmpLD$ESCORIAL)

LD <- merge(LD, tmpLD[, c("ID", "TRICALS")], by = "ID", all.x = T)

# LD where rows with is.na(VISIT)
LD2 <- LD[!is.na(LD$VISIT), ]

WD <- LD %>%
  pivot_wider(
    id_cols = c(ID, TRT, SEX, AGE, RACE, RILUSE, ONSET, FAMILIAL, DISDUR, DIAGDELAY,
                DELFRS, ESCORIAL, FVCL, FVCP, CENSOR, STIME),
    names_from = VISIT,
    values_from = c(TOT, CFB, TIME),
    values_fill = NA,
    names_glue = "V{VISIT}.{.value}"
  )

WD <- WD[, !names(WD) %in% c("NA.TOT", "V0.CFB", "NA.CFB", "NA.TIME")]
WD <- WD[, c("ID", "TRT", "SEX", "AGE", "RACE", "RILUSE", "ONSET", "FAMILIAL",
             "DISDUR", "DIAGDELAY", "DELFRS", "ESCORIAL", "FVCL", "FVCP", "CENSOR",
             "STIME", paste0("V", 0:6, ".TOT"), paste0("V", 1:6, ".CFB"),
            paste0("V", 0:6, ".TIME"))]

LOCF <- WD[, which(names(WD) == "V0.TOT"):which(names(WD) == "V6.CFB")]
LOCF <- t(apply(LOCF, 1, zoo::na.locf))
LOCF <- cbind(WD[, which(names(WD) == "ID"):which(names(WD) == "STIME")], LOCF,
              WD[, which(names(WD) == "V0.TIME"):which(names(WD) == "V6.TIME")])

# LOCF long
LOCF2 <- LD
LOCF2 <- merge(expand.grid(ID = unique(LD$ID), VISIT = 0:6), LD,
               by = c("ID", "VISIT"), all.x = TRUE)
LOCF2 <- LOCF2[order(LOCF2$ID, LOCF2$VISIT), ]

for (var in setdiff(names(LD), c("TOT", "CFB", "TIME", "VISIT"))) {
  LOCF2[[var]] <- ave(LOCF2[[var]], LOCF2$ID, FUN = function(x) {
    fv <- x[!is.na(x)][1]
    ifelse(is.na(x), fv, x)
  })
}

# LOCF TOT
LOCF2$TOT <- ave(LOCF2$TOT, LOCF2$ID, FUN = function(x) {
  na.locf(x, na.rm = FALSE)
})

# LOCF CFB
LOCF2$CFB <- ave(LOCF2$CFB, LOCF2$ID, FUN = function(x) {
  na.locf(x, na.rm = FALSE)
})

# TIME
LOCF2$TIME <- ifelse(is.na(LOCF2$TIME),
                     tapply(LOCF2$TIME, LOCF2$VISIT, mean, na.rm = TRUE)[as.character(LOCF2$VISIT)], LOCF2$TIME)

LOCF2[c("ID", "VISIT")] <- lapply(LOCF2[c("ID", "VISIT")], as.factor)

ZERO <- WD
ZERO[, which(names(ZERO) == "V0.TOT"):which(names(ZERO) == "V6.TOT")] <-
  lapply(ZERO[, which(names(ZERO) == "V0.TOT"):which(names(ZERO) == "V6.TOT")],
         function(x) {x[is.na(x)] <- 0; x})
ZERO[, which(names(ZERO) == "V1.CFB"):which(names(ZERO) == "V6.CFB")] <-
  ZERO[, c(paste0("V", 1:6, ".TOT"))] - ZERO$V0.TOT

COMP <- WD[complete.cases(WD[, paste0("V", 1:6, ".TOT")]), ]

A <- A[order(A$study_id, A$TIME), ]
COMP <- COMP[order(COMP$ID), ]
WD <- WD[order(WD$ID), ]
ZERO <- ZERO[order(ZERO$ID), ]
LOCF <- LOCF[order(LOCF$ID), ]

LD <- LD[order(LD$ID, LD$TIME), ]
LD2 <- LD2[order(LD2$ID, LD2$TIME), ]
LOCF2 <- LOCF2[order(LOCF2$ID, LOCF2$TIME), ]

rm(list = setdiff(ls(), c(lsf.str(), c("WD", "LD", "LD2", "LOCF", "LOCF2", "ZERO", "COMP"), "TE", "TRT.cenrate")))


# sum(rowSums(is.na(WD[, c("V1.CFB", "V2.CFB", "V3.CFB", "V4.CFB", "V5.CFB", "V6.CFB")])) > 0)
# prop.table(table(WD$CENSOR))*100
# prop.table(table((WD$CENSOR == 0 & is.na(WD$V6.CFB))))*100
# tapply((is.na(WD$V6.CFB) & WD$CENSOR == 0), WD$TRT, function(x) prop.table(table(x))*100)
# tapply(WD$CENSOR, WD$TRT, function(x) prop.table(table(x))*100)

# #### 04/07
# WD$DO <- ifelse(is.na(WD$V6.TOT), 1, 0)
# 
# m <- glm(DO ~ V0.TOT, family = binomial(link = "logit"), data = WD)
# summary(m)
# 
# #### 15/07
# library(lme4)
# m <- lmer(TOT ~ TIME + (TIME|ID), data = LD)
# VarCorr(m, comp = c("Variance"))
# 
# print(vc <- VarCorr(m))
# print(vc, comp = c("Variance"))
# 
# ## 21/07
# LD <- LD[!LD$ID == "298", ]
# LL.RE <- nlme::lme(fixed = TOT ~ TIME + TIME:TRT, random = ~ TIME | ID, data = LD)
# LQ.RE <- nlme::lme(fixed = TOT ~ TIME + TIME:TRT, random = ~ TIME + I(TIME^2) | ID, data = LD)
# QL.RE <- nlme::lme(fixed = TOT ~ TIME + TIME:TRT + I(TIME^2) + I(TIME^2):TRT, random = ~ TIME | ID, data = LD)
# QQ.RE <- nlme::lme(fixed = TOT ~ TIME + TIME:TRT + I(TIME^2) + I(TIME^2):TRT, random = ~ TIME + I(TIME^2) | ID, data = LD)
# 
# 
# DM <- function(m, eff = c("TRT", "PBO", "BOTH")){
#   beta <- fixef(m)
# 
#   if (eff == "TRT"){
#     contrast <- if(length(beta) == 3) c(0, 12, 12) else c(0, 12, 144, 12, 144)
#   } else if (eff == "PBO"){
#     contrast <- if(length(beta) == 3) c(0, 12, 0) else c(0, 12, 144, 0, 0)
#   } else if (eff == "BOTH"){
#     contrast <- if(length(beta) == 3) c(0, 0, 12) else c(0, 0, 0, 12, 144)
#   }
# 
#   vcv <- if (class(m) == "lme") vcov(m) else vcov(m)[1:length(contrast), 1:length(contrast)]
# 
#   est <- sum(contrast * beta)
#   se <- as.vector(sqrt(t(contrast) %*% vcv %*% contrast))
#   tval <- est / se
#   df <- 2500
# 
#   p <- 2 * pt(-abs(tval), df)
# 
#   return(list(est, p))
# }
# 
# 
# library(emmeans)
# MMRM <- mmrm::mmrm(CFB ~ VISIT + VISIT:TRT + us(VISIT|ID), data = LD2[!LD2$VISIT == 0, ])
# em <- emmeans(MMRM, specs = pairwise ~ TRT, at = list(VISIT = c("6")))$contrast
# 
# library(survival)
# cox <- coxph(Surv(STIME, CENSOR) ~ TRT, x = TRUE, data = WD)
# 
# 
# WD$TIME <- apply(WD[, c("V1.TIME", "V2.TIME", "V3.TIME", "V4.TIME", "V5.TIME", "V6.TIME")], 1, function(row) {
#   nna <- na.omit(row)
#   if (length(nna) == 0) {
#     return(0)
#   } else {
#     return(tail(nna, 1))
#   }
# })
# 
# LD <- LD[LD$ID %in% WD$ID, ]
# WD <- WD[WD$ID %in% LD$ID, ]
# 
# # alpha == value
# library(JM)
# jm1 <- jointModel(lmeObject = LL.RE, survObject = cox, timeVar = "TIME", parameterization = "value", method = "weibull-PH-aGH")
# jm2 <- jointModel(lmeObject = LQ.RE, survObject = cox, timeVar = "TIME", parameterization = "value", method = "weibull-PH-aGH")
# jm3 <- jointModel(lmeObject = QL.RE, survObject = cox, timeVar = "TIME", parameterization = "value", method = "weibull-PH-aGH")
# jm4 <- jointModel(lmeObject = QQ.RE, survObject = cox, timeVar = "TIME", parameterization = "value", method = "weibull-PH-aGH")
# 
# # alpha == slope
# dform5 <- list(fixed = ~ 1 + TRT, indFixed = 2:3,
#               random = ~ 1, indRandom = 2)
# jm5 <- jointModel(lmeObject = LL.RE, survObject = cox, timeVar = "TIME", parameterization = "slope", method = "weibull-PH-aGH", derivForm = dform5)
# 
# dform6 <- list(fixed = ~ 1 + TRT, indFixed = 2:3,
#               random = ~ 1 + I(2*TIME), indRandom = 2:3)
# jm6 <- jointModel(lmeObject = LQ.RE, survObject = cox, timeVar = "TIME", parameterization = "slope", method = "weibull-PH-aGH", derivForm = dform6)
# 
# dform7 <- list(fixed = ~ 1 + TRT + I(2*TIME) + I(2*TIME*TRT), indFixed = 2:5,
#               random = ~ 1, indRandom = 2)
# jm7 <- jointModel(lmeObject = QL.RE, survObject = cox, timeVar = "TIME", parameterization = "slope", method = "weibull-PH-aGH", derivForm = dform7)
# 
# dform8 <- list(fixed = ~ 1 + TRT + I(2*TIME) + I(2*TIME*TRT), indFixed = 2:5,
#               random = ~ 1 + I(2*TIME), indRandom = 2:3)
# jm8 <- jointModel(lmeObject = QQ.RE, survObject = cox, timeVar = "TIME", parameterization = "slope", method = "weibull-PH-aGH", derivForm = dform8,
#                   control = list(knots = c(4, 8, 11)))
# 
# jm9 <- jointModel(lmeObject = LL.RE, survObject = cox, timeVar = "TIME", parameterization = "both", method = "weibull-PH-aGH", derivForm = dform5)
# jm10 <- jointModel(lmeObject = LQ.RE, survObject = cox, timeVar = "TIME", parameterization = "both", method = "weibull-PH-aGH", derivForm = dform6)
# jm11 <- jointModel(lmeObject = QL.RE, survObject = cox, timeVar = "TIME", parameterization = "both", method = "weibull-PH-aGH", derivForm = dform7)
# jm12 <- jointModel(lmeObject = QQ.RE, survObject = cox, timeVar = "TIME", parameterization = "both", method = "weibull-PH-aGH", derivForm = dform8)
# 
# 
# # # plot 
# PBO <- sapply(c(unname(mget(paste0("jm", 1:12)))), function(x) DM(x, "PBO")[[1]])
# TRT <- sapply(c(unname(mget(paste0("jm", 1:12)))), function(x) DM(x, "TRT")[[1]])
# y <- TRT - PBO
# 
# # add space between different association structures
# spc <- function(vec) {
#   s <- rep(NA, length(vec) + floor(length(vec) / 4))
#   id <- seq_along(vec) + floor((seq_along(vec) - 1) / 4)
#   s[id] <- vec
#   return(s)
# }
# 
# PBO <- spc(PBO)
# TRT <- spc(TRT)
# y <- spc(round(y, 2))
# 
# # D1 <- data.frame(MT = rep(c(rep("Linear", 2), rep("Quadratic", 2)), times = 3),
# #                  IT = rep(c("Linear", "Quadratic"), 6),
# #                  A  = rep(c("Value", "Slope", "Value + slope"), each = 4))
# 
# D1 <- data.frame(TTRJ  = paste0(rep(c(rep("Linear", 2), rep("Quadratic", 2)), times = 3), "   +   ", rep(c("Linear", "Quadratic"), 6)), 
#                  A  = rep(c("Value", "Slope", "Value + slope"), each = 4))
# 
# spc.df <- function(df, every = 4) {
#   n <- nrow(df)
#   nr <- matrix(NA, nrow = length(seq(every, n, by = every)), ncol = ncol(df))
#   colnames(nr) <- names(df)
#   
#   splt <- split(df, ceiling(seq_len(n) / every))
#   ndf <- do.call(rbind, lapply(splt, function(block) {
#     rbind(block, nr[1, , drop = FALSE])
#   }))
#   
#   # Remove the last NA row if it was added at the end unnecessarily
#   if (all(is.na(ndf[nrow(ndf), ]))) {
#     ndf <- ndf[-nrow(ndf), ]
#   }
#   
#   return(ndf)
# }
# 
# D1 <- spc.df(D1)
# 
# 
# 
# pdf("/Users/dweemeri/surfdrive - Weemering, D.N. (Daphne)@surfdrive.surf.nl/robustJM/robustJM/figures/FIG1.pdf", width = 8, height = 5.5)
# layout(matrix(1:3, nrow = 1), widths = c(0.87, 2, 1))
# 
# par(mar = c(3.5, 1, 0.2, 0))
# plot(NULL, xlim = c(0, 25), ylim = c(1, length(PBO)), axes = F, frame.plot = F, xlab = "", ylab = "")
# 
# 
# text(x = 1.2, y = 1:nrow(D1), D1$TTRJ, adj = 0, cex = 1.2, xpd = T)
# text(x = 1.2, y = nrow(D1) + 1, "Mean trajectory +\nindividual trajectory", adj = 0, cex = 1.2, font = 2)
# # text(x = 1.2, y = nrow(D1) + 1, "Mean\ntrajectory", adj = 0, cex = 1.2, font = 2, )
# # text(x = 12, y = nrow(D1) + 1, " + ", adj = 0, cex = 1.2, font = 2)
# # text(x = 15, y = nrow(D1) + 1, "Individual\ntrajectory", adj = 0, cex = 1.2, font = 2)
# 
# 
# text(x = -1, y = 2.5, expression(alpha * italic(" = value")), srt = 90, cex = 1.2, xpd = T)
# text(x = -1, y = 7.5, expression(alpha * italic(" = slope")), srt = 90, cex = 1.2, xpd = T)
# text(x = -1, y = 12.5, expression(alpha * italic(" = value & slope")), srt = 90, cex = 1.2, xpd = T)
# 
# # first plot
# par(mar = c(3.5, 0.5, 0.2, 0), mgp = c(2, 0.7, 0))
# plot(NULL, xlim = c(-15.6, -13.4), ylim = c(1, length(PBO)),
#      xlab = "",
#      ylab = "", axes = FALSE, cex.lab = 1.2)
# 
# segments(x0 = PBO, y0 = 1:length(PBO), x1 = TRT, y1 = 1:length(PBO), lty = 3, col = "grey60")
# rect(xleft = TRT + 0.1, ybottom = 1:length(PBO) - 0.05, xright = -10, ytop = 1:length(PBO) + 0.05, col = "grey90", border = "grey90")
# 
# points(PBO, 1:length(PBO), pch = 22, cex = 2.2, col = "royalblue3", bg = "#95B4CC")
# points(TRT, 1:length(TRT), pch = 22, cex = 2.2, col = "maroon", bg = "#FF6961")
# text(x = TRT + 0.25, y = 1:length(PBO), labels = y, cex = 1)
# 
# axis(1, c(-15.6, -13.4), labels = rep("", 2), lwd = 1, lwd.ticks = 0)
# axis(1, seq(-15.5, -13.5, by = 0.5), cex.axis = 1.1, lwd = 0, lwd.ticks = 1)
# mtext("ALSFRS\u00adR at end of follow\u00adup", side = 1,
#       line = 2, at = -14.5, cex = 0.85)
# 
# 
# # second plot
# par(mar = c(3.5, 0, 0.2, 0.5))
# plot(y, 1:length(PBO), axes = F, xlim = c(0.93, 1.37), xlab = "",
#      ylab = "", cex.lab = 1.1)
# 
# rect(xleft = 0.7, ybottom = 1:length(y) - 0.05, xright = y, ytop = 1:length(y) + 0.05, col = "grey90", border = "grey90")
# points(y, 1:length(PBO), pch = 21, col = "grey40", bg = "grey60", cex = 2.2)
# axis(1, c(0.93, 1.37), labels = rep("", 2), lwd = 1, lwd.ticks = 0)
# axis(1, seq(0.95, 1.35, by = 0.2), cex.axis = 1.1, lwd = 0, lwd.ticks = 1)
# mtext("Treatment effect", side = 1, line = 2, at = 1.15, cex = 0.85)
# 
# 
# dev.off()
# 
# 
# 
# 
