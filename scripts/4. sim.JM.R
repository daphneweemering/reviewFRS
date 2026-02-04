################################################################################
##                                                                            ##
##                              DATA SIMULATION                               ##
##                                                                            ##
################################################################################

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

local({
# 1. functions required for simulating survival times --------------------------
# longitudinal process, current value and slope
mi1 <- function (t, B0 = 39.24208, B1 = -14.329175, B3 = 1.30610, u0i = u0i, u1i = u1i, u2i = u2i){                    
  B0 + u0i + B1*t + B3*t^2 + u1i*t + u2i*t^2
}

mi2 <- function(t, B1 = -14.329175, B3 = 1.30610, u1i = u1i, u2i = u2i){
  B1 + 2*B3*t + u1i + 2*u2i*t
}

# hazard function
h <- function (t,                           
               p = 1.0502200, g = 0.2449671,              
               TRT = 0, g1 = 0,          
               a1 = -0.0879232, a2 = -0.0486237,             
               B0 = 39.24208, B1 = -14.329175, B3 = 1.30610,
               u0i = u0i, u1i = u1i, u2i = u2i) {  
  
  p * exp (g + (g1 * TRT) + (a1 * mi1(t, B0 = B0, B1 = B1, B3 = B3, u0i = u0i, u1i = u1i,
                                      u2i = u2i)) + 
             (a2 * mi2 (t, B1 = B1, B3 = B3, u1i = u1i, u2i = u2i))) * (t^(p-1))
}

# survival probability for time t
SURV <- function (t, ...){
  exp (-integrate (h, 1e-05, t, ...)$value)
}

# inverse survival function
invS <- function (t, u,
                  p = p, g = g,
                  TRT = TRT, g1 = g1, 
                  a1 = a1, a2 = a2, B0 = B0, B1 = B1, B3 = B3,
                  u0i = u0i, u1i = u1i, u2i = u2i){
  SURV (t, p = p, TRT = TRT, g = g, g1 = g1, B0 = B0, B1 = B1, B3 = B3, a1 = a1, 
        a2 = a2, u0i = u0i, u1i = u1i, u2i = u2i)  -  u
}


# 2. simulate longitudinal and time-to-event data ------------------------------
vc <- matrix(c(54.1809471, 2.67757780, -0.201746783,
               2.67757780, 1.72379573, -0.078242505, 
               -0.2017468, -0.07824251, 0.005542269), nrow = 3, ncol = 3)

JMD <- function(N = 500, Nm = 7, Nm.pre, B2 = 0, B4 = 0, g1 = 0, FUT = 12, DO.TRT = 0.2, 
                DO.PLB = 0.2){ 
  
  # coefficients
  Us <- as.data.frame(mvrnorm(n = N, mu = c(0, 0, 0), Sigma = vc))
  
  colnames(Us) <- c("u0i", "u1i", "u2i")
  
  # 2:1 randomization
  D <- cbind(data.frame(ID  = seq(1:N), 
                        TRT = rep(0:1, times = c(round(1/3*N), N - round(1/3*N)))), Us)
  
  D$B0 <- 37.30083837
  D$B1 <- -1.26686491
  D$B3 <- 0.01598739
  
  D$u0i <- Us$u0i
  D$u1i <- Us$u1i
  D$u2i <- Us$u2i
  D$e   <- 4.736042
  
  D$p <- 1.0502200
  D$g <- -1.571722
  D$a1 <- -0.1088009 
  D$a2 <- -0.5024612
  
  D$g1 <- g1
  
  # numerical root finding to find survival times
  D$u <- runif(N)
  D$eT <- sapply(1:nrow(D), function(i) {
    L <- try(uniroot (invS, c(0, FUT), u = D[i, ]$u,  p = D[i, ]$p, g = D[i, ]$g,
                      TRT = D[i, ]$TRT, g1 = D[i, ]$g1, a1 = D[i, ]$a1, a2 = D[i, ]$a2,
                      B0 = D[i, ]$B0, B1 = D[i, ]$B1, B3 = D[i, ]$B3, u0i = D[i, ]$u0i,
                      u1i = D[i, ]$u1i, u2i = D[i, ]$u2i)$root, silent = T)
    if(class(L) == "try-error"){FUT}else{L}
  })

  ## random dropout
  D$cT <- rep(FUT, N)  # default censoring times (set to FUT)
  cenrates <- c("0" = DO.PLB, "1" = DO.TRT)  # cenrates per group
  n.all <- split(D$ID, D$TRT)  # split by TRT
  
  # select dropouts for each group 
  ci <- unlist(lapply(names(n.all), function(i) {
    sample(n.all[[i]], size = round(length(n.all[[i]]) * cenrates[i]), replace = FALSE)
  }))
  
  # left-skewed dropout times
  D$cT[ci] <- FUT * rbeta(length(ci), 3.5, 1.5)
  

  # longitudinal dataset with observed outcome Yi
  DL <- data.frame(ID    = rep(1:N, each = Nm + Nm.pre),
                   VISIT = rep(c(seq(-Nm.pre, -1, by = 1), seq(0, Nm-1)), times = N),
                   TIME  = rep(c(seq(-Nm.pre, -1, by = 1), seq(0, FUT, length.out = Nm)), times = N) + rnorm(N*(Nm + Nm.pre), 0, 5/365))
  DL$TIME[DL$VISIT == 0] <- 0
  
  DL <- merge(DL, D, by = "ID")
  
  DL$TOT <- with(DL, B0 + u0i +
                   # pre-treatment slope (same for all)
                   B1*TIME*(TIME < 0) + B3*TIME^2*(TIME < 0) +
                   # post-treatment slope (TE modifies slope)
                   B1*(1 - B2*TRT)*TIME*(TIME >= 0) +
                   B3*(1 - B4*TRT)*TIME^2*(TIME >= 0) +
                   # random slopes & quadratic terms & res. error
                   u1i*TIME + u2i*TIME^2 + rnorm(n = nrow(DL), mean = 0, sd = sqrt(e)))
  
  # get survival time given censoring
  D <- merge(D, DL[DL$VISIT == Nm-1, c("ID", "TIME")], by = "ID")
  
  D$EVENT <- ifelse((D$eT < D$cT) & (D$eT < D$TIME), 1, 0)
  D$STIME <- ifelse(D$EVENT == 1, D$eT, ifelse(D$cT < D$TIME, D$cT, D$TIME))
  
  DL <- merge(DL, D[, c("ID" , "EVENT","STIME")], by = "ID")
  
  DL$STATUS <- ifelse(((DL$TIME > DL$STIME) & DL$STIME < 12), 1, 0)
  
  # clean the data
  D <- D[order(D$ID), ]
  D$ID <- as.factor(as.integer(D$ID))
  
  DL <- DL[!DL$STATUS == 1, ] 
  DL <- DL[order(DL$ID, DL$TIME), ] 
  DL$ID <- as.factor(as.integer(DL$ID))
  D <- D[!(D$ID %in% setdiff(1:N, DL$ID)), ] 
  
  # CFB
  BSLN <- DL[DL$TIME == 0, c("ID", "TOT")]
  names(BSLN)[2] <- "BSLN"
  DL <- merge(DL, BSLN, by = "ID", all.x = TRUE)
  DL$CFB <- DL$TOT - DL$BSLN
  
  # reshape from long to wide format
  setDT(DL); setDT(D)
  D1 <- dcast(DL, ID ~ VISIT, value.var = c("TIME", "TOT", "CFB"))
  D <- D[D1, on = "ID"]

  cols_to_remove <- grep("^(TIME_-|TOT_-|CFB_-)", names(D), value = TRUE)
  D <- D[, !..cols_to_remove]
  
  # datasets with LOCF, CC, ZERO
  LOCF <- D[, which(names(D) == "TOT_0"):which(names(D) == "CFB_6")]
  LOCF <- t(apply(LOCF, 1, zoo::na.locf))
  LOCF <- data.frame(D[, which(names(D) == "ID"):which(names(D) == "TIME_6")], 
                     LOCF)
  
  CC <- D[complete.cases(D[, paste0("TOT_", 1:6)]), ]
  
  ZERO <- copy(D)
  ZERO[, (which(names(ZERO) == "TOT_0"):which(names(ZERO) == "TOT_6")) := 
         lapply(.SD, function(x) fifelse(is.na(x), 0, x)),
       .SDcols = which(names(ZERO) == "TOT_0"):which(names(ZERO) == "TOT_6")]
  
  ZERO[, grep("^CFB_\\d+$", names(ZERO), value = TRUE) := 
         .SD - TOT_0, 
       .SDcols = grep("^TOT_\\d+$", names(ZERO), value = TRUE)]
  
  # LOCF for the long data
  setDT(LOCF)
  LOCF2 <- melt(LOCF, id.vars = c("ID", "TRT", "EVENT", "STIME", "TOT_0"),
                measure.vars = list(
                  TIME = grep("^TIME_\\d+$", names(D), value = TRUE),
                  TOT = grep("^TOT_\\d+$", names(D), value = TRUE),
                  CFB = grep("^CFB_\\d+$", names(D), value = TRUE)),
                variable.name = "VISIT")
  setnames(LOCF2, "TOT_0", "BSLN")
  levels(LOCF2$VISIT) <- 0:6
  LOCF2 <- LOCF2[order(LOCF2$ID), ]
  
  LOCF3 <- rbind(DL[DL$VISIT < 0, c("ID", "TRT", "EVENT", "STIME", "BSLN", "VISIT", "TIME", "TOT", "CFB")], LOCF2)
  LOCF3 <- LOCF3[order(LOCF3$ID, LOCF3$VISIT), ]
  
  # add time values for the imputed TOT score rows
  LOCF3[is.na(TIME), TIME := (as.numeric(as.character(VISIT))*2) + rnorm(.N, 0, 5/365)]
  
  
  # remove CFB at V0 for the wide data.frames
  D <- D[, !c("CFB_0", "u0i", "u1i", "u2i", "B0", "B1", "B3", "e", "p", "g", "a1", "a2", "g1", "u")]
  LOCF <- LOCF[, !c("CFB_0", "u0i", "u1i", "u2i", "B0", "B1", "B3", "e", "p", "g", "a1", "a2", "g1", "u")]
  CC <- CC[, !c("CFB_0", "u0i", "u1i", "u2i", "B0", "B1", "B3", "e", "p", "g", "a1", "a2", "g1", "u")]
  ZERO <- ZERO[, !"CFB_0"]
  DL2 <- DL
  DL <- DL[DL$TIME >= 0, ]
  
  # . Study-based adjustments --------------------------------------------------
  # Aggarwal (2010) ____________________________________________________________
  WD.AGG <- copy(D)
  setDT(WD.AGG)
  WD.AGG[, EVENT2 := as.integer(
    EVENT == 1 | rowSums(do.call(cbind, lapply(.SD, function(x) x <= -6)), na.rm = TRUE) > 0
  ), .SDcols = paste0("CFB_", 1:6)]
  
  # VISIT where CFB <= -6
  WD.AGG[, idx := min(which(unlist(.SD) <= -6)), by = seq_len(nrow(WD.AGG)), .SDcols = paste0("CFB_", 1:6)]
  WD.AGG$idx <- ifelse(WD.AGG$idx == "Inf", NA, WD.AGG$idx)
  
  # TIME corresponding to first CFB_FRS <= -6
  tc <- paste0("TIME_", 1:6)
  WD.AGG[!is.na(idx) & idx > 0, TIME.FRS := unlist(Map(function(row, col) row[[col]],
                                                       split(.SD, seq_len(nrow(.SD))),
                                                       tc[idx])), .SDcols = tc]
  
  # fill in STIME2
  WD.AGG[, STIME2 := fcase(
    EVENT == 1, STIME,  # condition 1
    EVENT == 0 & EVENT2 == 0, STIME,  # condition 2
    EVENT == 0 & EVENT2 == 1 & !is.na(TIME.FRS), TIME.FRS,  # condition 3
    default = 12  # fallback
  )]
  
  
  # joint rank analyses ________________________________________________________
  WD.JR <- as.data.table(D[!is.na(D$EVENT), c("ID", "EVENT", "STIME", "TRT")])
  LD.JR <- as.data.table(DL[, c("ID", "VISIT", "CFB")])
  
  CFB <- as.data.table(D[, c("ID", paste0("CFB_", 1:6))])
  CFB[, lastCFB := apply(.SD, 1, function(x) {
    last_val <- tail(na.omit(x), 1)
    if (length(last_val) == 0) NA_real_ else last_val
  }), .SDcols = patterns("^CFB_\\d+$")]
  
  WD.JR <- CFB[, .(ID, lastCFB)][WD.JR, on = "ID"]
  WD.JR[, c("SLOPE", "lastCFB") := .(lastCFB/STIME, NULL)]
  # WD.JR[, SLOPE := lastCFB/STIME]
  
  # change event [0 == died; 1 == unknown] & change VISIT to numeric
  WD.JR[, EVENT := 1 - EVENT]
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
  finscore <- finscore[D[,c("ID","TRT", "TOT_0")], on = "ID", nomatch = 0]
  setnames(finscore, "x", "SUM")
  setnames(finscore, "x2", "SUM.genge")
  WD.JR <- as.data.frame(finscore)
  
  
  # 27. -- Miller - 2022 - tofersen - PMID: 36129998 ---------------------------
  WD <- copy(D)
  predM <- make.predictorMatrix(WD)
  
  # imputations based on TRT, RILUSE and BSLN score
  predM[, !colnames(predM) %in% c("TRT", "TOT_0")] <- 0
  
  # only impute V6.TOT using pmm
  meth <- make.method(WD)
  meth[!names(meth) %in% "TOT_6"] <- ""
  
  # impute datasets using TRT, RILUSE
  imp <- mice(WD, m = 5, maxit = 5, method = meth, predictorMatrix = predM)
  WD.MILL <- complete(imp, "long")
  
  # CFB
  WD.MILL <- as.data.table(WD.MILL)
  WD.MILL[is.na(CFB_6), CFB_6 := TOT_6 - TOT_0, by = .imp]
  
  # split WD.MILL by .imp
  WD.list <- split(WD.MILL, by = ".imp")
  
  # process each subset separately
  WD.list <- lapply(WD.list, function(dt) {
    
    # change censor [0 == died; 1 == unknown] & remove imputed data for patients who died
    dt[, EVENT := 1 - EVENT]
    dt[EVENT == 0 & is.na(TIME_6), CFB_6 := NA]
    
    dt <- dt[!is.na(EVENT), .(.imp, ID, EVENT, STIME, TRT, CFB_6)]
    LD.MILL <- as.data.table(LOCF2[, c("ID", "VISIT")])
    
    M1 <- LD.MILL[dt, on = "ID"]
    M2 <- M1[M1, on = "VISIT", allow.cartesian = T][ID != i.ID & VISIT == 6]
    M2[, i..imp := NULL]
    
    data.table::setnames(M2, new = c("ID.x", "VISIT", "imp", "EVENT.x", "STIME.x", "TRT.x", "CFB.x",
                                     "ID.y", "EVENT.y", "STIME.y", "TRT.y", "CFB.y"))
    
    M2$score <-
      # If both die compare survival times
      if_else(M2$EVENT.x == 0 & M2$EVENT.y == 0, sign(M2$STIME.x - M2$STIME.y),
              # If x died and y did not
              if_else(M2$EVENT.x == 0 & M2$EVENT.y == 1,
                      if_else(is.na(M2$CFB.x) == T, -1, sign(M2$CFB.x - M2$CFB.y)),
                      # If y died and x did not
                      if_else(M2$EVENT.x == 1 & M2$EVENT.y == 0,
                              if_else(is.na(M2$CFB.y) == T, 1, sign(M2$CFB.x - M2$CFB.y)),
                              # If both are alive
                              if_else(M2$EVENT.x == 1 & M2$EVENT.y == 1,
                                      if_else(!is.na(M2$CFB.x) & !is.na(M2$CFB.y),
                                              sign(M2$CFB.x - M2$CFB.y), 0), 0))))
    
    # sum final scores
    finscore <- M2[, .(x = sum(score)), by = ID.x]
    data.table::setnames(finscore, "ID.x", "ID")
    data.table::setorder(finscore, ID)
    finscore[, RANK := rank(x)]
    finscore <- finscore[WD[,c("ID","TRT", "TOT_0")], on = "ID", nomatch = 0]
    setnames(finscore, "x", "SUM")
    WD.MILL <- as.data.frame(finscore)
    
    return(WD.MILL)
    
  })
  
  WD.MILL <- rbindlist(WD.list)
  WD.MILL[, .imp := 1:5, by = ID]
  
  
  # Mora (2020) ________________________________________________________________
  LOCF.MORA <- copy(D)
  LOCF.MORA$TOT_6 <- ifelse((is.na(LOCF.MORA$TOT_6) & LOCF.MORA$EVENT == 1), 0, LOCF.MORA$TOT_6)
  LOCF.MORA[, paste0("TOT_", 0:6) := as.data.table(t(apply(.SD, 1, function(row) {
    out <- row
    for (i in seq_along(row)) {
      if (is.na(out[i]) && i > 1) out[i] <- out[i - 1]
    }
    out
  }))), .SDcols = paste0("TOT_", 0:6)]
  
  LOCF.MORA$CFB_6 <- LOCF.MORA$TOT_6 - LOCF.MORA$TOT_0
  
  # Kaji (2019) ________________________________________________________________
  WD.KAJI <- copy(D)
  WD.KAJI$CFB_6 <- ifelse((is.na(WD.KAJI$CFB_6) & WD.KAJI$EVENT == 1), -999, WD.KAJI$CFB_6)
  
  
  WD.KAJI$CFB_6 <- apply(WD.KAJI[, grep("^CFB_[0-9]+$", names(WD.KAJI), value = TRUE), with = FALSE], 1, function(row) {
    lv <- tail(na.omit(row), 1)
    if (length(lv) == 0) NA else lv
  })
  
  WD.KAJI[is.na(WD.KAJI$CFB_6), ]$CFB_6 <- 0 # cases with no FU data
  
  # Kaufmann (2009) ____________________________________________________________
  WD.KAUF <- copy(D)
  
  # "neirest neighbour" worst score imputation
  nn <- function(i, d) {
    lrow <- d[i, ] 
    prox <- d[d$TRT == lrow$TRT & !is.na(d$CFB_6), ]
    prox$diff <- abs(prox$TOT_0 - lrow$TOT_0) 
    
    return(min(prox[order(prox$diff), ][1:5, ]$CFB_6)) 
  }
  
  for (i in which(is.na(WD.KAUF$CFB_6) & WD.KAUF$EVENT != 1)) {
    WD.KAUF$CFB_6[i] <- nn(i, WD.KAUF)
  }
  
  # impute zero for deceased patients with missing data
  WD.KAUF$TOT_6 <- ifelse((is.na(WD.KAUF$TOT_6) & WD.KAUF$EVENT == 1), 0, WD.KAUF$TOT_6)
  WD.KAUF$CFB_6 <- ifelse((is.na(WD.KAUF$CFB_6) & WD.KAUF$EVENT == 1), WD.KAUF$TOT_6 - WD.KAUF$TOT_0, WD.KAUF$CFB_6)
  WD.KAUF <- WD.KAUF[!is.na(WD.KAUF$CFB_6), ]
  
  # Berry (2023) _______________________________________________________________
  # select valid IDs with >1 post-baseline observation
  vID <- DL2[TIME >= 0, .N, by = ID][N > 1, ID]
  
  # calculate slopes per ID
  LOCF.BERRY <- DL2[ID %in% vID, {
    POSTSLP <- summary(lm(TOT ~ TIME, data = .SD[TIME >= 0]))$coefficients[2, 1]
    PRESLP  <- summary(lm(TOT ~ TIME, data = .SD[TIME < 0]))$coefficients[2, 1]
    
    .(SLP = c(PRESLP, POSTSLP),
      PERIOD = c(0, 1))  # 0 = pre, 1 = post
  }, by = ID]
  
  LOCF.BERRY <- merge(LOCF.BERRY, unique(DL2[, .(ID, TRT)]), by = "ID", all.x = TRUE)
  LOCF.BERRY <- LOCF.BERRY[order(LOCF.BERRY$ID), ]
  
  # Cudkowicz (2022) & Elia (2015) _____________________________________________
  # compute PRE, POST, DELTA from model
  R <- function(dt) {
    dt[, {
      mod <- lm(TOT ~ TIME * I(TIME >= 0), data = .SD)
      c <- coef(mod)
      DELTA <- c["TIME:I(TIME >= 0)TRUE"]
      .(PRE = c["TIME"], POST = c["TIME"] + DELTA, DELTA = DELTA)
    }, by = ID]
  }
  
  RESP.ELIA <- R(LOCF3)
  RESP.CUD <- R(LOCF3)
  
  RESP.ELIA <- merge(RESP.ELIA, unique(LOCF3[, c("ID", "TRT", "BSLN", "EVENT")]), by = "ID", all.x = TRUE)
  RESP.CUD <- merge(RESP.CUD, unique(LOCF3[, c("ID", "TRT", "BSLN", "EVENT")]), by = "ID", all.x = TRUE)
  
  # responder status
  RESP.ELIA[, RESP := as.integer(DELTA >= B2)]
  RESP.CUD[, RESP := as.integer(!(DELTA < B2 | EVENT == 1))]

  LD <- DL[, !c("u0i", "u1i", "u2i", "B0", "B1", "B3", "e", "p", "g", "a1", "a2", "g1", "u", "eT", "cT")]
  WD <- D[, !c("eT", "cT")]
  
  LD$VISIT <- as.factor(LD$VISIT)
  
  return(list(LD = LD, WD = WD, LOCF = LOCF, ZERO = ZERO, CC = CC, LOCF2 = LOCF2,
              WD.AGG = WD.AGG, WD.JR = WD.JR, WD.MILL = WD.MILL, LOCF.MORA = LOCF.MORA, 
              WD.KAJI = WD.KAJI, WD.KAUF = WD.KAUF, LOCF.BERRY = LOCF.BERRY, 
              RESP.ELIA = RESP.ELIA, RESP.CUD = RESP.CUD))
  
}

assign("JMD", JMD, envir = .GlobalEnv)
})


