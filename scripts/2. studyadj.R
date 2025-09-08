# adjustments to the datasets based on study-specifics
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("1. TE.R")
# source("TEMP.TE.R")

library(mice)
library(data.table)


  ############################################################################
  ############ study-based data adjustments for the ceftriaxone dataset ######
  ############################################################################

D.ADJ.CEF <- function(TE, D = "CEF", TRT.cenrate){ # TRT.cenrate <- 1 == doubled dropout; 0 is original data
  D <- D.TE(TE = TE, D = D, TRT.cenrate)
  
  COMP <- D$COMP
  LD <- D$LD
  LOCF <- D$LOCF
  LOCF2 <- D$LOCF2
  WD <- D$WD
  ZERO <- D$ZERO
  
  # 4. -- Aggarwal - 2010 - lithium - PMID: 20363190 -----------------------------
  WD.AGG <- WD[!WD$ID == "298", ]

  # CENSOR2 for survival or ALSFRS-R <= -6
  WD.AGG$CENSOR2 <- apply(WD.AGG, 1, function(row) {
    cns <- as.numeric(row["CENSOR"]) == 1
    cns2 <- any(as.numeric(row[c("V1.CFB", "V2.CFB", "V3.CFB", "V4.CFB", "V5.CFB", "V6.CFB")]) <= -6, na.rm = TRUE)
    if (cns | cns2) 1 else 0
  })

  # STIME2 with the adjusted conditions
  WD.AGG$STIME2 <- apply(WD.AGG, 1, function(row) {
    idx <- which(as.numeric(row[c("V1.CFB", "V2.CFB", "V3.CFB", "V4.CFB", "V5.CFB", "V6.CFB")]) <= -6)

    if (as.numeric(row["CENSOR"]) == 1) {
      # condition 1: Someone died (CENSOR == 1)
      return(as.numeric(row["STIME"]))
    } else if (!as.numeric(row["CENSOR"]) == 1 && !as.numeric(row["CENSOR2"]) == 1) {
      # condition 2: Someone did not die (CENSOR == 0) and did not experience a decrease in ALSFRS-R <= -6 (CENSOR2 == 0)
      return(as.numeric(row["STIME"]))
    } else if (!as.numeric(row["CENSOR"]) == 1 && as.numeric(row["CENSOR2"]) == 1 && length(idx) > 0) {
      # condition 3: Someone did not die (CENSOR == 0), but experienced a decrease in ALSFRS-R <= -6 (CENSOR2 == 1)
      # assign the time corresponding to the last decrease in ALSFRS-R
      return(as.numeric(row[c("V1.TIME", "V2.TIME", "V3.TIME", "V4.TIME", "V5.TIME", "V6.TIME")])[max(idx)])
    } else {
      return(12)
    }
  })


  # # 7. -- Berry - 2019 - MSC-NTF - PMID: 31740545 --------------------------------
  # LD.BER <- LD[LD$ID %in% names(table(LD$ID)[table(LD$ID) > 1]), ]
  # POSTSLP <- lapply(unique(LD.BER$ID), function(i){
  #   POSTSLP <- summary(lm(TOT ~ TIME, data = subset(LD, LD$ID == i)))$coefficients[2, 1]
  # })
  # 
  # LD.BER$ID <- droplevels(LD.BER$ID)
  # LD.BER <- data.frame(ID = as.factor(rep(names(table(LD$ID)[table(LD$ID) > 1]), 2)),
  #                       SLP = c(tapply(LD.BER$DELFRS, LD.BER$ID, function(x) x[1]), do.call("rbind", POSTSLP)),
  #                       PER = rep(c(0, 1), each = length(c(names(table(LD$ID)[table(LD$ID) > 1])))))
  # LD.BER$TRT <- LD$TRT[match(LD.BER$ID, LD$ID)]


  # 10. -- Cudkowicz - 2013 - dexpramipexole - PMID: 24067398 --------------------
  # 16. -- Genge - 2023 - ravulizumab - PMID: 37695623 -------------------------
  # 24. -- Meininger - 2017 - ozanezumab - PMID: 28139349 ------------------------
  WD.JR <- as.data.table(WD[!is.na(WD$CENSOR), c("ID", "CENSOR", "STIME", "TRT")])
  LD.JR <- as.data.table(LD[!LD$ID == 298, c("ID", "VISIT", "CFB")])

  CFB <- as.data.table(WD[, c("ID", paste0("V", 1:6, ".CFB"))])
  CFB[, lastCFB := apply(.SD, 1, function(x) {
    last_val <- tail(na.omit(x), 1)
    if (length(last_val) == 0) NA_real_ else last_val
  }), .SDcols = patterns("^V\\d+\\.CFB")]

  WD.JR <- CFB[, .(ID, lastCFB)][WD.JR, on = "ID"]
  WD.JR[, c("SLOPE", "lastCFB") := .(lastCFB/STIME, NULL)]

  # change censor [0 == died; 1 == unknown] & change VISIT to numeric
  WD.JR[, CENSOR := 1 - CENSOR]
  LD.JR[, VISIT := as.numeric(VISIT)]

  # merge TTE (wide) and long datasets
  M1 <- LD.JR[WD.JR, on = "ID"]

  # cross join the combined data with itself on VISIT
  M2 <- M1[M1, on = "VISIT", allow.cartesian = T]

  # find out last shared visit for each pair of subjects
  maxvis <- data.table::copy(M2)[, .(x = max(VISIT)), by = c("ID", "i.ID")]

  # merge on the maximum visit and filter to combinations of different subjects
  # at their maximum visit
  M3 <- M2[maxvis, on = c("ID", "i.ID")][ID != i.ID & VISIT == x]
  data.table::setnames(M3, new = c("ID.x", "VISIT", "CFB.x", "CENSOR.x", "STIME.x", "TRT.x", "SLOPE.x",
                                   "ID.y", "CFB.y", "CENSOR.y", "STIME.y", "TRT.y", "SLOPE.y", "x"))

  # The method Cudkowicz (2013) and Meiniger (2017) used
  M3$score <- # If both die compare death times
    if_else(M3$CENSOR.x == 0 & M3$CENSOR.y == 0, sign(M3$STIME.x - M3$STIME.y),
            # If both survive compare function at last shared visit or a tie if
            # a comparison is not available
            if_else(M3$CENSOR.x == 1 & M3$CENSOR.y == 1,
                    if_else(!is.na(M3$CFB.x) & !is.na(M3$CFB.y),
                            sign(M3$CFB.x - M3$CFB.y),
                            0),
                    # If x died and y did not
                    if_else(M3$CENSOR.x == 0 & M3$CENSOR.y == 1,
                            # see if y has follow up after x death time
                            if_else(M3$STIME.x <= M3$STIME.y,
                                    # x loses if so
                                    -1,
                                    # otherwise return to comparison on function
                                    if_else(!is.na(M3$CFB.x) & !is.na(M3$CFB.y),
                                            sign(M3$CFB.x - M3$CFB.y),
                                            0)),
                            # If y died and x did not
                            if_else(M3$CENSOR.x == 1 & M3$CENSOR.y == 0,
                                    # see if x had followup after y death time
                                    if_else(M3$STIME.x >= M3$STIME.y,
                                            # x wins if so
                                            1,
                                            # otherwise return to comparison on function
                                            if_else(!is.na(M3$CFB.x) & !is.na(M3$CFB.y),
                                                    sign(M3$CFB.x - M3$CFB.y),
                                                    0)
                                    ),
                                    # This NA should be impossible because it
                                    # would require a value for CENSOR other than
                                    # 0 or 1
                                    NA)
                    )
            )
    )

  # The method Genge (2023) used
  M3$score.genge <-
    # If both die compare survival times
    if_else(M3$CENSOR.x == 0 & M3$CENSOR.y == 0, sign(M3$STIME.x - M3$STIME.y),
            # If x died and y did not
            if_else(M3$CENSOR.x == 0 & M3$CENSOR.y == 1, -1,
                    # If y died and x did not
                    if_else(M3$CENSOR.x == 1 & M3$CENSOR.y == 0, 1,
                            # If both are alive
                            if_else(M3$CENSOR.x == 1 & M3$CENSOR.y == 1,
                                    if_else(!is.na(M3$SLOPE.x) & !is.na(M3$SLOPE.y),
                                            sign(M3$SLOPE.x - M3$SLOPE.y), 0), 0))))


  # sum final scores and obtain ranks for Cudkowicz & Meiniger
  finscore <- M3[, .(x = sum(score), x2 = sum(score.genge)), by = ID.x]
  data.table::setnames(finscore, "ID.x", "ID")
  data.table::setorder(finscore, ID)
  finscore[, c("RANK", "RANK.genge") := .(rank(x), rank(x2))]
  finscore <- finscore[WD[,c("ID","TRT", "AGE", "SEX", "V0.TOT", "FVCP", "DISDUR",
                             "ONSET", "RILUSE")], on = "ID", nomatch = 0]
  setnames(finscore, "x", "SUM")
  setnames(finscore, "x2", "SUM.genge")
  WD.JR <- as.data.frame(finscore)


  # # 12. -- Cudkowicz - 2022 - MSC-NTF - PMID: 34890069 ---------------------------
  # # 15. -- Elia - 2015 - tauroursodeoxycholic acid - PMID: 25664595 ------------
  # LOCF2 <- as.data.table(LOCF2)
  # 
  # # obtain pre-slopes by extracting individual progression rates
  # m <- lmer(TOT ~ TIME + (TIME|ID), data = LD)
  # PRE <- data.table(ID = unique(LOCF2$ID),
  #                   PRE = fixef(m)["TIME"] + ranef(m)$ID$TIME)
  # 
  # RESP <- LOCF2[, .(POST = coef(lm(TOT ~ TIME))[2],
  #                   ID = ID, TRT = TRT, BSLN = BSLN, DISDUR = DISDUR,
  #                   ONSET = ONSET, RILUSE = RILUSE, CENSOR = CENSOR), by = ID]
  # 
  # RESP <- RESP[PRE, on = "ID"]
  # RESP <- RESP[, .SD[1], by = ID]
  # LOCF.RESP <- RESP[, c("RESP.ELIA", "RESP.CUD") := .(if_else(((POST - PRE)/PRE) <= TE, 1, 0),
  #                  if_else((CENSOR == 1 | ((POST - PRE)/PRE) <= TE), 1, 0)),
  #              by = ID]


  # 19. -- Kaji - 2019 - methylcobalamin - PMID: 30636701 ----------------------
  WD.KAJI <- WD
  WD.KAJI$V6.CFB <- ifelse((is.na(WD.KAJI$V6.CFB) & WD.KAJI$CENSOR == 1), -999, WD.KAJI$V6.CFB)
  WD.KAJI$V6.CFB <- apply(WD.KAJI[, grep("^V[0-9]+\\.CFB$", names(WD.KAJI))], 1, function(row) {
    lv <- tail(na.omit(row), 1)
    if (length(lv) == 0) NA else lv
  })

  WD.KAJI[is.na(WD.KAJI$V6.CFB), ]$V6.CFB <- 0 # two cases with no FU data


  # 20. -- Kaufmann - 2009 - COQ10 - PMID: 19743457 ----------------------------
  WD.KAUF <- WD
  
  # "neirest neighbour" worst score imputation
  nn <- function(i, d) {
    lrow <- d[i, ] # select row for each ID with missing data on V6.CFB
    prox <- d[d$TRT == lrow$TRT & !is.na(d$V6.CFB), ] # all IDs in same TRT with no missing data on V6.CFB (prox group)
    prox$diff <- abs(prox$V0.TOT - lrow$V0.TOT) # difference between ID and all IDs in the prox group

    return(min(prox[order(prox$diff), ][1:5, ]$V6.CFB)) # return V6.CFB of the five patients with BSLN closest to ID and pick the biggest decline
  }

  for (i in which(is.na(WD.KAUF$V6.CFB) & WD.KAUF$CENSOR != 1)) {
    WD.KAUF$V6.CFB[i] <- nn(i, WD.KAUF)
  }

  # impute zero for deceased patients with missing data
  WD.KAUF$V6.TOT <- ifelse((is.na(WD.KAUF$V6.TOT) & WD.KAUF$CENSOR == 1), 0, WD.KAUF$V6.TOT)
  WD.KAUF$V6.CFB <- ifelse((is.na(WD.KAUF$V6.CFB) & WD.KAUF$CENSOR == 1), WD.KAUF$V6.TOT - WD.KAUF$V0.TOT, WD.KAUF$V6.CFB)
  # WD.KAUF$V6.CFB <- ifelse((is.na(WD.KAUF$V6.CFB) & WD.KAUF$CENSOR == 1), 0, WD.KAUF$V6.CFB)
  WD.KAUF <- WD.KAUF[!is.na(WD.KAUF$V6.CFB), ]


  # 27. -- Miller - 2022 - tofersen - PMID: 36129998 ---------------------------
  predM <- make.predictorMatrix(WD)

  # imputations based on TRT, RILUSE and BSLN score
  predM[, !colnames(predM) %in% c("TRT", "RILUSE", "V0.TOT")] <- 0

  # only impute V6.TOT using pmm
  meth <- make.method(WD)
  meth[!names(meth) %in% "V6.TOT"] <- ""

  # impute datasets using TRT, RILUSE
  imp <- mice(WD, m = 5, maxit = 5, method = meth, predictorMatrix = predM)
  WD.MILL <- complete(imp, "long")

  # CFB
  WD.MILL <- as.data.table(WD.MILL)
  WD.MILL[is.na(V6.CFB), V6.CFB := V6.TOT - V0.TOT, by = .imp]

  # split WD.MILL by .imp
  WD.list <- split(WD.MILL, by = ".imp")

  # process each subset separately
  WD.list <- lapply(WD.list, function(dt) {

    # change censor [0 == died; 1 == unknown] & remove imputed data for patients who died
    dt[, CENSOR := 1 - CENSOR]
    dt[CENSOR == 0 & is.na(V6.TIME), V6.CFB := NA]

    dt <- dt[!is.na(CENSOR), .(.imp, ID, CENSOR, STIME, TRT, V6.CFB)]
    LD.MILL <- as.data.table(LOCF2[!LOCF2$ID == 298, c("ID", "VISIT")])

    M1 <- LD.MILL[dt, on = "ID"]
    M2 <- M1[M1, on = "VISIT", allow.cartesian = T][ID != i.ID & VISIT == 6]
    M2[, i..imp := NULL]

    data.table::setnames(M2, new = c("ID.x", "VISIT", "imp", "CENSOR.x", "STIME.x", "TRT.x", "CFB.x",
                                     "ID.y", "CENSOR.y", "STIME.y", "TRT.y", "CFB.y"))

    M2$score <-
      # If both die compare survival times
      if_else(M2$CENSOR.x == 0 & M2$CENSOR.y == 0, sign(M2$STIME.x - M2$STIME.y),
              # If x died and y did not
              if_else(M2$CENSOR.x == 0 & M2$CENSOR.y == 1,
                      if_else(is.na(M2$CFB.x) == T, -1, sign(M2$CFB.x - M2$CFB.y)),
                      # If y died and x did not
                      if_else(M2$CENSOR.x == 1 & M2$CENSOR.y == 0,
                              if_else(is.na(M2$CFB.y) == T, 1, sign(M2$CFB.x - M2$CFB.y)),
                              # If both are alive
                              if_else(M2$CENSOR.x == 1 & M2$CENSOR.y == 1,
                                      if_else(!is.na(M2$CFB.x) & !is.na(M2$CFB.y),
                                              sign(M2$CFB.x - M2$CFB.y), 0), 0))))

    # sum final scores
    finscore <- M2[, .(x = sum(score)), by = ID.x]
    data.table::setnames(finscore, "ID.x", "ID")
    data.table::setorder(finscore, ID)
    finscore[, RANK := rank(x)]
    finscore <- finscore[WD[,c("ID","TRT", "DISDUR", "V0.TOT", "RILUSE")], on = "ID", nomatch = 0]
    setnames(finscore, "x", "SUM")
    WD.MILL <- as.data.frame(finscore)

    return(WD.MILL)

  })

  WD.MILL <- rbindlist(WD.list)
  WD.MILL[, .imp := 1:5, by = ID]

  # 29. -- Mora - 2020 - masitinib - PMID: 31280619 ----------------------------
  LOCF.MORA <- WD
  LOCF.MORA$V6.TOT <- ifelse((is.na(LOCF.MORA$V6.TOT) & LOCF.MORA$CENSOR == 1), 0, LOCF.MORA$V6.TOT)
  LOCF.MORA$V6.CFB <- LOCF.MORA$V6.TOT - LOCF.MORA$V0.TOT
  
  LOCF.MORA.cols <- LOCF.MORA[, which(names(LOCF.MORA) == "V0.TOT"):which(names(LOCF.MORA) == "V6.CFB")]
  LOCF.MORA <- as.data.frame(t(apply(LOCF.MORA.cols, 1, function(x) na.locf(x, na.rm = FALSE))))
  LOCF.MORA <- cbind(WD[, which(names(WD) == "ID"):which(names(WD) == "V6.TIME")], LOCF.MORA)

  # rm(list = setdiff(ls(), c(lsf.str(), c("WD", "LD", "LOCF", "LOCF2", "ZERO",
  #                                        "COMP", "WD.AGG", "LD.BER", "WD.JR", "LOCF.RESP",
  #                                        "WD.KAJI", "WD.KAUF", "WD.MILL", "LOCF.MORA"), "TE", "TRT.cenrate")))
  # 
  # list(WD = WD, LD = LD, LOCF = LOCF, LOCF2 = LOCF2, ZERO = ZERO,
  #      COMP = COMP, WD.AGG = WD.AGG, LD.BER = LD.BER, WD.JR = WD.JR,
  #      LOCF.RESP = LOCF.RESP, WD.KAJI = WD.KAJI, WD.KAUF = WD.KAUF,
  #      WD.MILL = WD.MILL, LOCF.MORA = LOCF.MORA)
  
  rm(list = setdiff(ls(), c(lsf.str(), c("WD", "LD", "LOCF", "LOCF2", "ZERO",
                                         "COMP", "WD.AGG", "WD.JR", "WD.KAJI", 
                                         "WD.KAUF", "WD.MILL", "LOCF.MORA"), "TE", "TRT.cenrate")))
  
  list(WD = WD, LD = LD, LOCF = LOCF, LOCF2 = LOCF2, ZERO = ZERO, COMP = COMP, 
       WD.AGG = WD.AGG, WD.JR = WD.JR, WD.KAJI = WD.KAJI, WD.KAUF = WD.KAUF,
       WD.MILL = WD.MILL, LOCF.MORA = LOCF.MORA)

}







 ############################################################################
 ############ study-based data adjustments for the RT001 dataset ############
 ############################################################################

D.ADJ.RT001 <- function(TE, D = "RT001", TRT.cenrate){
  D <- D.TE(TE = TE, D = D, TRT.cenrate)

  COMP <- D$COMP
  LD <- D$LD
  LOCF <- D$LOCF
  LOCF2 <- D$LOCF2
  WD <- D$WD
  ZERO <- D$ZERO

  # 4. -- Aggarwal - 2010 - lithium - PMID: 20363190 ---------------------------
  WD.AGG <- WD
  WD.AGG$EVENT <- ifelse(WD.AGG$V1.CFB <= -6 | WD.AGG$V2.CFB <= -6 | WD.AGG$V3.CFB <= -6, 1, 0)
  WD.AGG$EVENT <- ifelse(is.na(WD.AGG$EVENT) == T, 0, WD.AGG$EVENT)
  WD.AGG$STIME <- ifelse(WD.AGG$V1.CFB <= -6, WD.AGG$V1.TIME,
                         ifelse(WD.AGG$V2.CFB <= -6, WD.AGG$V2.TIME,
                                ifelse(WD.AGG$V3.CFB <= -6, WD.AGG$V3.TIME, 7)))

  # censored IDs
  totc <- c("V0.TOT", "V1.TOT", "V2.TOT", "V3.TOT")
  timec <- c("V0.TIME", "V1.TIME", "V2.TIME", "V3.TIME")
  
  for (i in which(is.na(WD.AGG$STIME) & apply(WD.AGG[totc], 1, function(x) any(is.na(x))))) {
    l <- max(which(!is.na(unlist(WD.AGG[i, totc]))))
    WD.AGG$STIME[i] <- WD.AGG[i, timec[l]]
  }

  WD.AGG$STIME <- unlist(WD.AGG$STIME)

  # 10. -- Cudkowicz - 2013 - dexpramipexole - PMID: 24067398 --------------------
  # 16. -- Genge - 2023 - ravulizumab - PMID: 37695623 ---------------------------
  # 24. -- Meininger - 2017 - ozanezumab - PMID: 28139349 ------------------------
  # add a time variable
  LD <- setDT(LD)
  STIME <- LD[, .(STIME = max(TIME, na.rm = TRUE)), by = ID]
  WD <- merge(WD, STIME, by = "ID", all.x = TRUE)
  WD$STIME <- ifelse(is.na(WD$V3.CFB) == F, 7, WD$STIME)
  
  WD.JR <- as.data.table(WD[, c("ID", "TRT", "STIME")])
  LD.JR <- as.data.table(LD[, c("ID", "VISIT", "CFB")])
  
  CFB <- as.data.table(WD[, c("ID", paste0("V", 1:3, ".CFB"))])
  CFB[, lastCFB := apply(.SD, 1, function(x) {
    last_val <- tail(na.omit(x), 1)
    if (length(last_val) == 0) NA_real_ else last_val
  }), .SDcols = patterns("^V\\d+\\.CFB")]
  CFB$lastCFB <- ifelse(is.na(CFB$lastCFB), 0, CFB$lastCFB)
  
  WD.JR <- CFB[, .(ID, lastCFB)][WD.JR, on = "ID"]
  WD.JR[, c("SLOPE", "lastCFB") := .(lastCFB/STIME, NULL)]
  
  # change VISIT to numeric
  LD.JR[, VISIT := as.numeric(VISIT)]
  
  # merge TTE (wide) and long datasets
  M1 <- LD.JR[WD.JR, on = "ID"]
  
  # cross join the combined data with itself on VISIT
  M2 <- M1[M1, on = "VISIT", allow.cartesian = T]
  
  # find out last shared visit for each pair of subjects
  maxvis <- data.table::copy(M2)[, .(x = max(VISIT)), by = c("ID", "i.ID")]
  
  # merge on the maximum visit and filter to combinations of different subjects
  # at their maximum visit
  M3 <- M2[maxvis, on = c("ID", "i.ID")][ID != i.ID & VISIT == x]
  data.table::setnames(M3, new = c("ID.x", "VISIT", "CFB.x", "TRT.x", "STIME.x", "SLOPE.x",
                                   "ID.y", "CFB.y", "TRT.y", "STIME.y", "SLOPE.y", "x"))
  
  # The method Cudkowicz (2013) and Meiniger (2017) used
  M3$score <- if_else(!is.na(M3$CFB.x) & !is.na(M3$CFB.y), sign(M3$CFB.x - M3$CFB.y), 0)
  M3$score.genge <- if_else(!is.na(M3$SLOPE.x) & !is.na(M3$SLOPE.y), sign(M3$SLOPE.x - M3$SLOPE.y), 0)
  
  # sum final scores and obtain ranks for Cudkowicz & Meiniger
  finscore <- M3[, .(x = sum(score), x2 = sum(score.genge)), by = ID.x]
  data.table::setnames(finscore, "ID.x", "ID")
  data.table::setorder(finscore, ID)
  finscore[, c("RANK", "RANK.genge") := .(rank(x), rank(x2))]
  finscore <- finscore[WD[,c("ID","TRT", "AGE", "SEX", "V0.TOT", "SVCP", "DISDUR",
                             "ONSET", "RILUSE", "CNTRY")], on = "ID", nomatch = 0]
  setnames(finscore, "x", "SUM")
  setnames(finscore, "x2", "SUM.genge")
  WD.JR <- as.data.frame(finscore)


  # 20. -- Kaufmann - 2009 - COQ10 - PMID: 19743457 ----------------------------
  WD.KAUF <- WD
  
  # "neirest neighbour" worst score imputation
  nn <- function(i, d) {
    lrow <- d[i, ] # select row for each ID with missing data on V3.CFB
    prox <- d[d$TRT == lrow$TRT & !is.na(d$V3.CFB), ] # all IDs in same TRT with no missing data on V3.CFB (prox group)
    prox$diff <- abs(prox$V0.TOT - lrow$V0.TOT) # difference between ID and all IDs in the prox group
    
    return(min(prox[order(prox$diff), ][1:5, ]$V3.CFB)) # return V6.CFB of the five patients with BSLN closest to ID and pick the biggest decline
  }
  
  for (i in which(is.na(WD.KAUF$V3.CFB))) {
    WD.KAUF$V3.CFB[i] <- nn(i, WD.KAUF)
  }


  # 27. -- Miller - 2022 - tofersen - PMID: 36129998 ---------------------------
  predM <- make.predictorMatrix(WD)
  
  # imputations based on TRT, RILUSE and BSLN score
  predM[, !colnames(predM) %in% c("TRT", "RILUSE", "V0.TOT")] <- 0
  
  # only impute V6.TOT using pmm
  meth <- make.method(WD)
  meth[!names(meth) %in% "V3.TOT"] <- ""
  
  # impute datasets using TRT, RILUSE
  imp <- mice(WD, m = 5, maxit = 5, method = meth, predictorMatrix = predM)
  WD.MILL <- complete(imp, "long")
  
  # CFB
  WD.MILL <- as.data.table(WD.MILL)
  WD.MILL[is.na(V3.CFB), V3.CFB := V3.TOT - V0.TOT, by = .imp]
  
  # split WD.MILL by .imp
  WD.list <- split(WD.MILL, by = ".imp")
  
  # process each subset separately
  WD.list <- lapply(WD.list, function(dt) {
    
    dt <- dt[, .(.imp, ID, STIME, TRT, V3.CFB)]
    LD.MILL <- as.data.table(LOCF2[, c("ID", "VISIT")])
    
    M1 <- LD.MILL[dt, on = "ID"]
    M2 <- M1[M1, on = "VISIT", allow.cartesian = T][ID != i.ID & VISIT == 3]
    M2[, i..imp := NULL]
    
    data.table::setnames(M2, new = c("ID.x", "VISIT", "imp", "STIME.x", "TRT.x", "CFB.x",
                                     "ID.y", "STIME.y", "TRT.y", "CFB.y"))
    
    M2$score <- if_else(!is.na(M2$CFB.x) & !is.na(M2$CFB.y), sign(M2$CFB.x - M2$CFB.y), 0)
    
    # sum final scores
    finscore <- M2[, .(x = sum(score)), by = ID.x]
    data.table::setnames(finscore, "ID.x", "ID")
    data.table::setorder(finscore, ID)
    finscore[, RANK := rank(x)]
    finscore <- finscore[WD[,c("ID","TRT", "DISDUR", "V0.TOT", "RILUSE")], on = "ID", nomatch = 0]
    setnames(finscore, "x", "SUM")
    WD.MILL <- as.data.frame(finscore)
    
    return(WD.MILL)
    
  })
  
  WD.MILL <- rbindlist(WD.list)
  WD.MILL[, .imp := 1:5, by = ID]


  rm(list = setdiff(ls(), c(lsf.str(), c("WD", "LD", "LOCF", "LOCF2", "ZERO",
                                         "COMP", "WD.AGG", "WD.JR", "WD.KAUF", 
                                         "WD.MILL"), "TE", "TRT.cenrate")))

  list(WD = WD, LD = LD, LOCF = LOCF, LOCF2 = LOCF2, ZERO = ZERO, COMP = COMP, 
       WD.AGG = WD.AGG, WD.JR = WD.JR, WD.KAUF = WD.KAUF, WD.MILL = WD.MILL)

}
