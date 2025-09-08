library(zoo)
library(dplyr)
library(data.table)

D.TE <- function(TE, D = c("RT001", "CEF"), TRT.cenrate) {
  # D <- "RT001" # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! REMOVE
  
  MAX.VISIT <- ifelse(D == "RT001", 3, 6)
  
  if(D == "RT001"){
    source("/Users/dweemeri/surfdrive - Weemering, D.N. (Daphne)@surfdrive.surf.nl/reviewFRS/reviewFRS/scripts/0. clean.RT001.R")
  } else{
    source("/Users/dweemeri/surfdrive - Weemering, D.N. (Daphne)@surfdrive.surf.nl/reviewFRS/reviewFRS/scripts/0. clean.CEF.R")
  }
  
  # D <- "RT001" # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! REMOVE
  # TE <- 0 # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! REMOVE
  
  
  # 1. PERMUTE TREATMENT ASSIGNMENT + TRT EFFECT (optional) ____________________
  # permute TRT 
  ptrt <- sample(LD$TRT[!duplicated(LD$ID)])
  for (j in seq_along(unique(LD$ID))) {
    LD$TRT[LD$ID == unique(LD$ID)[j]] <- ptrt[j]
  }
  
  # long dataset
  TRT <- lapply(unique(LD$ID), function(i) {
    d <- LD[LD$ID == i, ]
    if (all(d$TRT == 1)) {
      i_TE <- rnorm(1, mean = TE, sd = TE/2)
      sq <- cumsum(c(0, diff(d$TIME) * i_TE))
      # sq <- cumsum(c(0, diff(d$TIME) * TE)) --> original
      add <- floor(sq) + rbinom(length(d$TIME), size = 1, prob = sq - floor(sq))
      d$TOT + add
    } else d$TOT
  })
  
  LD$TOT <- unlist(TRT)
  LD$CFB <- LD$TOT - tapply(LD$TOT[LD$TIME == 0], unique(LD$ID[LD$TIME == 0]), mean)[as.character(LD$ID)]
  
  LD <- LD[!is.na(LD$VISIT), ]
  
  
  # 2. UNBALANCED DROPOUT ______________________________________________________
  # TRT.cenrate <- 0.2 # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! REMOVE
  
  if(TRT.cenrate > 0){
    LD <- setDT(LD)
    LD[, CFB := as.numeric(CFB)]
    tWD <- dcast(LD, ID ~ VISIT, value.var = c("TIME", "TOT", "CFB"))
    
    if(D == "CEF"){
      #. CEFTRIAXONE DATA ------------------------------------------------------
      tWD <- merge(tWD, unique(LD[, c("ID", "TRT", "CENSOR")]), by = "ID", all.x = TRUE)
      
      # 1. define the at-risk group (alive at end, TRT == 1)
      AR <- which(tWD$TRT == 1 & tWD$CENSOR == 0)
      
      # 2. get observed dropout rate
      obsP <- unname(prop.table(table(tWD[tWD$TRT == 1, ]$CENSOR == 0 & is.na(tWD[tWD$TRT == 1, ]$CFB_6)))[2])
      
      # 3. new dropout patterns and times
      tWD$DO <- ifelse(tWD$TRT == 1 & tWD$CENSOR == 0, rbinom(length(AR), size = 1, prob = (obsP*2)), 0)
      tWD$DO.t <- round(6.5 * rbeta(513, 4, 1.5)) # 6.5 bc 6 visits with a little margin for rounding
      
      # 4. create dropout by setting CFB and TOT NA based on DO.t
      cfb_cols <- grep("^CFB_", names(tWD), value = TRUE)
      tot_cols <- grep("^TOT_", names(tWD), value = TRUE)
      
      for (i in which(tWD$TRT == 1 & tWD$CENSOR == 0 & !is.na(tWD$CFB_6) & tWD$DO == 1)) {
        do.t <- tWD$DO.t[i]
        
        # determine which CFB and TOT columns to remove based on DO.t
        if (do.t == 1) {
          RMC <- c(cfb_cols[2:7], tot_cols[2:7])  # CFB_1 to CFB_6 and TOT_1 to TOT_6
        } else if (do.t >= 2 && do.t <= 6) {
          RMC <- c(cfb_cols[which(grepl(paste0("CFB_", do.t), cfb_cols)):length(cfb_cols)],
                   tot_cols[which(grepl(paste0("TOT_", do.t), tot_cols)):length(tot_cols)])
        } else if (do.t == 6) {
          RMC <- c(cfb_cols[which(grepl("CFB_6", cfb_cols))], tot_cols[which(grepl("TOT_6", tot_cols))])
        } else {
          RMC <- character(0)
        }
        
        # set those columns to NA for the current row
        tWD[i, (RMC) := NA]
      }
      
      # reshape tWD back to long format using melt
      tLD <- melt(tWD, id.vars = c("ID", "TRT", "TOT_0"),
                  measure.vars = list(
                    TIME = grep("^TIME_\\d+$", names(tWD), value = TRUE),
                    TOT = grep("^TOT_\\d+$", names(tWD), value = TRUE),
                    CFB = grep("^CFB_\\d+$", names(tWD), value = TRUE)),
                  variable.name = "VISIT")
      levels(tLD$VISIT) <- 0:6
      setnames(tLD, "TOT_0", "BSLN")
      tLD <- tLD[order(tLD$ID), ]
      tLD <- tLD[!is.na(tLD$TOT), ]
      LD[, c("TRT", "BSLN", "VISIT", "TIME", "TOT", "CFB") := NULL]
      LD <- merge(tLD, LD[!duplicated(LD$ID), ], by = "ID")
      
    } else {
      #. RT001 DATA ------------------------------------------------------------
      tWD <- merge(tWD, unique(LD[, c("ID", "TRT")]), by = "ID", all.x = TRUE)
      
      # 1. define the at-risk group (alive at end, TRT == 1)
      AR <- which(tWD$TRT == 1)
      
      # 2. get observed dropout rate
      obsP <- unname(prop.table(table(is.na(tWD[tWD$TRT == 1, ]$CFB_3)))[2])
      
      # 3. new dropout patterns and times
      tWD$DO <- ifelse(tWD$TRT == 1, rbinom(length(AR), size = 1, prob = obsP*2), 0)
      tWD$DO.t <- round(3.5 * rbeta(43, 4, 1.5)) # 3.5 bc 3 visits with a little margin for rounding
      
      # 4. create dropout by setting CFB and TOT NA based on DO.t
      cfb_cols <- grep("^CFB_", names(tWD), value = TRUE)
      tot_cols <- grep("^TOT_", names(tWD), value = TRUE)
      
      for (i in which(tWD$TRT == 1 & !is.na(tWD$CFB_3) & tWD$DO == 1)) {
        do.t <- tWD$DO.t[i]
        
        # determine which CFB and TOT columns to remove based on DO.t
        if (do.t == 1) {
          RMC <- c(cfb_cols[2:4], tot_cols[2:4])  # CFB_1 to CFB_3 and TOT_1 to TOT_3
        } else if (do.t == 2) {
          RMC <- RMC <- c(cfb_cols[3:4], tot_cols[3:4])  # CFB_2 to CFB_3 and TOT_2 to TOT_3
        } else if (do.t == 3) {
          RMC <- RMC <- c(cfb_cols[4], tot_cols[4])  # CFB_3 and TOT_3
        } else {
          RMC <- character(0)
        }
        
        # set those columns to NA for the current row
        tWD[i, (RMC) := NA]
      }
      
      # reshape tWD back to long format using melt
      tLD <- melt(tWD, id.vars = c("ID", "TRT", "TOT_0"),
                  measure.vars = list(
                    TIME = grep("^TIME_\\d+$", names(tWD), value = TRUE),
                    TOT = grep("^TOT_\\d+$", names(tWD), value = TRUE),
                    CFB = grep("^CFB_\\d+$", names(tWD), value = TRUE)),
                  variable.name = "VISIT")
      levels(tLD$VISIT) <- 0:3
      setnames(tLD, "TOT_0", "BSLN")
      tLD <- tLD[order(tLD$ID), ]
      tLD <- tLD[!is.na(tLD$TOT), ]
      LD[, c("TRT", "BSLN", "VISIT", "TIME", "TOT", "CFB") := NULL]
      LD <- merge(tLD, LD[!duplicated(LD$ID), ], by = "ID")
    }
  }
  
  
  # 3. OTHER DATASETS __________________________________________________________
  
  # MAX.VISIT <- 3 # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! REMOVE
  
  # LOCF2 (long dataset)
  LOCF2 <- LD
  LOCF2 <- merge(expand.grid(ID = unique(LD$ID), VISIT = 0:MAX.VISIT), LD,           
                 by = c("ID", "VISIT"), all.x = TRUE)
  LOCF2 <- LOCF2[order(LOCF2$ID, LOCF2$VISIT), ]
  
  # Fill in LOCF for TOT
  LOCF2$TOT <- ave(LOCF2$TOT, LOCF2$ID, FUN = function(x) {
    na.locf(x, na.rm = FALSE)
  })
  
  # Fill in LOCF for CFB
  LOCF2$CFB <- ave(LOCF2$CFB, LOCF2$ID, FUN = function(x) {
    na.locf(x, na.rm = FALSE)
  })
  
  # Fill in LOCF for TIME
  LOCF2$TIME <- ifelse(is.na(LOCF2$TIME),
                       tapply(LOCF2$TIME, LOCF2$VISIT, mean, na.rm = TRUE)[as.character(LOCF2$VISIT)], LOCF2$TIME)
  
  # Fill in open cells for fixed variables
  for (var in setdiff(names(LD), c("ID", "TOT", "CFB", "TIME", "VISIT"))) { 
    LOCF2[[var]] <- ave(LOCF2[[var]], LOCF2$ID, FUN = function(x) {
      fv <- x[!is.na(x)][1]
      ifelse(is.na(x), fv, x)
    })
  }
  
  LOCF2[c("ID", "VISIT")] <- lapply(LOCF2[c("ID", "VISIT")], as.factor)
  
  # Wide dataset
  LD$VISIT <- as.factor(LD$VISIT)
  
  if (MAX.VISIT == 6){
    WD <- LD %>%
      pivot_wider(
        id_cols = c(ID, TRT, SEX, AGE, RACE, RILUSE, ONSET, FAMILIAL, DISDUR, DIAGDELAY,
                    DELFRS, ESCORIAL, FVCL, FVCP, CENSOR, STIME),
        names_from = VISIT,
        values_from = c(TOT, CFB, TIME),
        values_fill = NA,
        names_glue = "V{VISIT}.{.value}"
      )
    
    WD <- WD[, !names(WD) %in% c("VNA.TOT", "V0.CFB", "VNA.CFB", "VNA.TIME")]
    
    # reorder columns based on prefixes and VISIT
    ord <- paste0("V", 1:6, ".CFB")
    ord <- c(ord, unlist(lapply(c("TOT", "TIME"), function(x) paste0("V", 0:6, ".", x))))
    ord <- c(unlist(lapply(c("TIME", "TOT"), function(x) paste0("V", 0:6, ".", x))),
             paste0("V", 1:6, ".CFB"))
    ord <- c(setdiff(colnames(WD), ord), ord)
    WD <- WD[, ord]
  } else {
    
    WD <- LD %>%
      pivot_wider(
        id_cols = c(ID, TRT, CNTRY, ONSET, AGE, RACE, SEX, RILUSE, DISDUR, DIAGDELAY, 
                    DELFRS, SVCL, SVCP, MITOS, BSLN, SCRN, SCRN.TIME, AGE_SYMP, TRICALS),
        names_from = VISIT,
        values_from = c(TIME, TOT, CFB),
        values_fill = NA,
        names_glue = "V{VISIT}.{.value}"
      )
  }
  
  # LOCF
  LOCF.cols <- WD[, which(names(WD) == "V0.TOT"):which(names(WD) == paste0("V", MAX.VISIT, ".CFB"))]
  LOCF <- as.data.frame(t(apply(LOCF.cols, 1, function(x) na.locf(x, na.rm = FALSE))))
  LOCF <- cbind(WD[, which(names(WD) == "ID"):which(names(WD) == paste0("V", MAX.VISIT, ".TIME"))], LOCF)
  
  # ZERO
  ZERO <- WD
  ZERO[, which(names(ZERO) == "V0.TOT"):which(names(ZERO) == paste0("V", MAX.VISIT, ".TOT"))] <-
    lapply(ZERO[, which(names(ZERO) == "V0.TOT"):which(names(ZERO) == paste0("V", MAX.VISIT, ".TOT"))],      
           function(x) {x[is.na(x)] <- 0; x})
  ZERO[, which(names(ZERO) == "V1.CFB"):which(names(ZERO) == paste0("V", MAX.VISIT, ".CFB"))] <-
    ZERO[, c(paste0("V", 1:MAX.VISIT, ".TOT"))] - ZERO$V0.TOT                          
  
  # COMP
  COMP <- WD[complete.cases(WD[, paste0("V", 1:MAX.VISIT, ".TOT")]), ]  
  
  
  rm(list = setdiff(ls(), c(lsf.str(), c("WD", "LD", "LOCF", "LOCF2", "ZERO", "COMP"), "TE", "TRT.cenrate")))
  list(WD = WD, LD = LD, LOCF = LOCF, LOCF2 = LOCF2, ZERO = ZERO, COMP = COMP)
}

