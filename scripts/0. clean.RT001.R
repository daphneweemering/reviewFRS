########################### ALSFRS-R ANALYSIS REVIEW ###########################

#### 0. libraries ####
library(readxl)
library(lme4)
library(tidyverse)


######## RT001 DATA FOR EXAMPLE SCENARIOS ########
setwd('/Users/dweemeri/surfdrive - Weemering, D.N. (Daphne)@surfdrive.surf.nl/data/RT001')
D <- read_excel('study_results.xlsx')

# remove ineligible patients
D <- D[!is.na(D$SCRN_VSFHGT_EU), ]

# age
D$BRT_MNTH <- match(D$SCRN_DEMMON, month.abb) # month in numbers
D$BRT_DATE <- as.Date(paste(D$SCRN_DEMDAY, '-' ,D$BRT_MNTH,'-', D$SCRN_DEMYEA, sep = ''), format = '%d-%m-%Y') 
D$AGE <- as.numeric(as.Date(D$SCRN_DOVDAT, format = '%d-%m-%Y') - D$BRT_DATE) / 365.25

# disease duration
D1 <- read_excel('added_data.xlsx')
colnames(D)[1] <- 'ID'
D <- merge(D, D1, by = 'ID', all.x = T)
# D[15, ]$SCRN_PREDAT <- '27-06-2020'
D$DISDUR <- as.numeric(as.Date(D$SCRN_DOVDAT, format = '%d-%m-%Y') - as.Date(D$date_onset, format = '%d-%m-%Y')) / (365.25/12)

# age at symptom onset
D$AGE_SYMP <- as.numeric(as.Date(D$date_onset, format = '%d-%m-%Y') - as.Date(D$BRT_DATE)) / 365.25

# diagnostic delay
D$DIAGDELAY <- as.numeric(as.Date(D$SCRN_PREDAT, format = "%d-%m-%Y") - as.Date(D$date_onset, format = '%d-%m-%Y')) / (365.25/12)
D$DIAGDELAY <- ifelse(D$DIAGDELAY <= 0, 1, D$DIAGDELAY) # one patient with negative diagnostic delay, one patient with zero DIAGDELAY. Both set to 1 month

# riluzole use
D2 <- read_excel('concomitant_medication.xlsx')
D2 <- subset(D2, D2$CMRNAM == 'riluzole' | D2$CMRNAM == 'Riluzole' |
               D2$CMRNAM == 'Riluzol' | D2$CMRNAM == 'riluzol')

colnames(D2)[1] <- 'ID'
D2 <- merge(D[, c('ID', 'SCRN_SVCOUN')], D2[, c('ID', 'CMRNAM')], by = 'ID', all.x = T)
D$RILUSE <- ifelse(is.na(D2$CMRNAM) == T, 0, 1)

colnames(D)[colnames(D) == 'SCRN_SVCOUN'] <- 'CNTRY'
colnames(D)[colnames(D) == 'bulbar'] <- 'ONSET'
colnames(D)[colnames(D) == 'SCRN_DEMRAC'] <- 'RACE'
colnames(D)[colnames(D) == 'SCRN_DEMSEX'] <- 'SEX'
colnames(D)[colnames(D) == 'SCRN_SVCSVC'] <- 'SVCL'
colnames(D)[colnames(D) == 'SCRN_SVCSVP'] <- 'SVCP'

# treatment assignment
D3 <- read_excel('randomization.xlsx')
D3 <- D3[!is.na(D3$TRX), ]
D$TRT <- ifelse(D3$TRX == 'Active', 1, 0)

## ALSFRS-R ##
alsfrs <- list(unique(grep(paste(c('ALRSPE', 'ALRSAL', 'ALRSWA', 'ALRHAN', 'ALRCFN', 'ALRCFY',
                                   'ALRDRE', 'ALRTUR', 'ALRWAL', 'SCRN_ALRCLI', 'ALRDYS', 
                                   'ALRORT', 'ALRRES'), collapse = "|"), 
                           names(D), value = T)))

# remove text
pos <- 1
for(i in alsfrs){
  alsfrs[[pos]] <- lapply(D[, i], function(x) substr(x, 1, 1))
  pos <- pos + 1
}

alsfrs <- rapply(alsfrs, function(x) ifelse(grepl('#', x), NA, x), how = 'replace')
alsfrs <- as.data.frame(rapply(alsfrs, as.numeric, how = 'replace'))

# move alsfrs-r data to correct measurement for 118-004 and 118-008
# 118-004
alsfrs[31, c(which(colnames(alsfrs) == 'V3_ALRSPE'):which(colnames(alsfrs) == 'V3_ALRRES'))] <-
  alsfrs[31, c(which(colnames(alsfrs) == 'V5_ALRSPE'):which(colnames(alsfrs) == 'V5_ALRRES'))]
alsfrs[31, c(which(colnames(alsfrs) == 'V5_ALRSPE'):which(colnames(alsfrs) == 'V5_ALRRES'))] <- NA

# 118-008
alsfrs[34, c(which(colnames(alsfrs) == 'V4_ALRSPE'):which(colnames(alsfrs) == 'V4_ALRRES'))] <-
  alsfrs[34, c(which(colnames(alsfrs) == 'V5_ALRSPE'):which(colnames(alsfrs) == 'V5_ALRRES'))]
alsfrs[34, c(which(colnames(alsfrs) == 'V5_ALRSPE'):which(colnames(alsfrs) == 'V5_ALRRES'))] <- NA

# get total alsfrs-r scores
lst <- list('SCRN', 'V2', 'V3', 'V4', 'V5')
totals <- as.data.frame(do.call(cbind, 
                                lapply(lst, function(x) rowSums(alsfrs[, grepl(x, names(alsfrs))], na.rm = T))))
totals[totals == 0] <- NA
names(totals) <- c('SCRN.TOT', 'V0.TOT', 'V1.TOT', 'V2.TOT', 'V3.TOT')

# change from baseline scores
cfb <- as.data.frame(totals[, 3:5] - totals[, 2])

# cfb <- as.data.frame(apply(totals, 2, function(x) x - totals[,2]))
# names(cfb) <- paste('c', names(totals), sep = '')
names(cfb) <- c("V1.CFB", "V2.CFB", "V3.CFB")

# get time 
time <- D[, grep('ALRDAT', names(D), value = T)]
time <- as.data.frame(apply(time, 2, function(x) (as.Date(x, '%d-%m-%Y') - as.Date(D$V2_ALRDAT, '%d-%m-%Y')) / (365.25/12)))
time[time > 100] <- NA


# put times in right place
time[31, 'V3_ALRDAT'] <- time[31, 'V5_ALRDAT']; time[31, 'V5_ALRDAT'] <- NA
time[34, 'V4_ALRDAT'] <- time[34, 'V5_ALRDAT']; time[34, 'V5_ALRDAT'] <- NA

# drop OLE times 
time <- time[-c(6:8)]

# names(time) <- paste('t', names(totals), sep = '')
names(time) <- c('SCRN.TIME', 'V0.TIME', 'V1.TIME', 'V2.TIME', 'V3.TIME')

# FRS slope
#D$DISDUR <- ifelse(D$DISDUR == 0, 0.0001, D$DISDUR)
D$DELFRS <- (totals$V0.TOT - 48) / D$DISDUR

# severity (MiToS)
MITOS <- within(alsfrs, {
  movement      <- as.integer(alsfrs$V2_ALRWAL <= 1 | alsfrs$V2_ALRDRE <= 1)
  swallowing    <- as.integer(alsfrs$V2_ALRSWA <= 1)
  communicating <- as.integer(alsfrs$V2_ALRSPE <= 1 & alsfrs$V2_ALRHAN <= 1)
  breathing     <- as.integer(alsfrs$V2_ALRDYS <= 1 | alsfrs$V2_ALRRES <= 2)
  staging       <- movement + swallowing + communicating + breathing 
})

D$MITOS <- MITOS[, 'staging']


# data frame (wide, long)
WD <- data.frame(D[, c('ID', 'TRT', 'CNTRY', 'ONSET', 'AGE', 'RACE', 'SEX', 'RILUSE', 
                       'DISDUR', 'DIAGDELAY', 'DELFRS', 'SVCL', 'SVCP', 'MITOS', 'AGE_SYMP')], totals, cfb, time)

WD[, c('SVCL', 'SVCP')] <- lapply(WD[, c('SVCL', 'SVCP')], as.numeric)

# TRICALS risk score
WD$TRICALS <- (0.474 * ((WD$SVCP/100)^-1) + ((WD$SVCP/100)^-0.5)) -
  (2.376 * ((WD$DIAGDELAY/10)^-0.5) + log(WD$DIAGDELAY/10)) - 
  (1.839 * ((-WD$DELFRS + 0.1)^-0.5)) - (0.264 * (WD$AGE_SYMP/100)^-2) + 
  (0.271 * WD$ONSET)

#### OLD ####
# LD <- gather(WD, V, value, V0.TOT:V4.TIME) %>%
#   separate(V, c('VISIT', 'col')) %>%
#   arrange(ID) %>%
#   spread(col, value)
#############

LD <- gather(WD, V, value, SCRN.TOT:V3.TIME) %>%
  separate(V, c('VISIT', 'col')) %>%
  arrange(ID) %>%
  spread(col, value)

# add baseline FRS column
extr <- lapply(split(LD, LD$ID), function(df) c(BSLN = df$TOT[2], SCRN = df$TOT[1], SCRN.TIME = df$TIME[1]))
extr <- do.call("rbind", extr)
extr <- data.frame(ID = rownames(extr), extr)
rownames(extr) <- NULL

LD <- merge(LD, extr, by = "ID", all.x = T)
LD <- LD[!LD$VISIT %in% "SCRN", ]
# LD <- LD[!LD$VISIT %in% c("SCRN", "V0"), ]

LD$VISIT <- as.factor(gsub("V", "", LD$VISIT))

LD$CFB <- ifelse(LD$VISIT == 0, 0, LD$CFB)

#### OLD ####
# BSLN <- as.data.frame(cbind(D$ID, LD[!duplicated(LD$ID),]$TOT))
# names(BSLN) <- c('ID', 'BSLN')
# LD <- merge(LD, BSLN, by = 'ID', all.x = T)
# LD$BSLN <- as.numeric(LD$BSLN)
# LD <- LD[!LD$VISIT == 'V0',]
#############


#### Missing data scenarios 
# LOCF wide
LOCF <- WD[, which(names(WD) == "SCRN.TOT"):which(names(WD) == "V3.CFB")]
LOCF[-1] <- t(apply(LOCF[-1], 1, zoo::na.locf))
LOCF <- cbind(WD[, which(names(WD) == "ID"):which(names(WD) == "MITOS")], LOCF, 
               WD[, which(names(WD) == "SCRN.TIME"):which(names(WD) == "V3.TIME")])
LOCF$V2.TIME <- ifelse(is.na(LOCF$V2.TIME) == T, mean(LOCF$V2.TIME, na.rm = T), LOCF$V2.TIME)
LOCF$V3.TIME <- ifelse(is.na(LOCF$V3.TIME) == T, mean(LOCF$V3.TIME, na.rm = T), LOCF$V3.TIME)

# LOCF long
LOCF2 <- LD
LOCF2$TOT <- c(NA, LD$TOT[!is.na(LD$TOT)])[cumsum(!is.na(LD$TOT)) + 1]
LOCF2$CFB <- c(NA, LOCF2$CFB[!is.na(LOCF2$CFB)])[cumsum(!is.na(LOCF2$CFB)) + 1]

# zero imputation wide
ZERO <- WD
ZERO[, which(names(ZERO) == "SCRN.TOT"):which(names(ZERO) == "V3.TOT")] <- 
  lapply(ZERO[, which(names(ZERO) == "SCRN.TOT"):which(names(ZERO) == "V3.TOT")], 
         function(x) {x[is.na(x)] <- 0; x})

ZERO[, which(names(ZERO) == "V1.CFB"):which(names(ZERO) == "V3.CFB")] <-
  ZERO[, c("V1.TOT", "V2.TOT", "V3.TOT")] - ZERO$V0.TOT

ZERO$V2.TIME <- ifelse(is.na(ZERO$V2.TIME) == T, mean(ZERO$V2.TIME, na.rm = T), ZERO$V2.TIME)
ZERO$V3.TIME <- ifelse(is.na(ZERO$V3.TIME) == T, mean(ZERO$V3.TIME, na.rm = T), ZERO$V3.TIME)

# complete case analysis (wide)
COMP <- WD[!(is.na(WD$V0.TOT) | is.na(WD$V1.TOT) | is.na(WD$V2.TOT) | is.na(WD$V3.TOT)),]

# remove missing rows from long dataset and keep one with missing rows for potential imputation
LD <- LD[!is.na(LD$TOT), ]

# rm(list = setdiff(ls(), c("D" ,"WD", "LD", "LD2", "LOCF", "LOCF2", "ZERO", "ZERO2", "COMP")))
rm(list = setdiff(ls(), c(lsf.str(), c("WD", "LD", "LD2", "LOCF", "LOCF2", "ZERO", "COMP"), "TE", "TRT.cenrate")))

