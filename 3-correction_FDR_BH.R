# PART 4 - Calculate Benjamini-Hochberg (B-H) false discovery rates (FDRs)
# of the 37 disease IVW MR results

# Read in IVW MR results table (Supplementary File 1e) created in PART 2 script
mydata <- read.table("Supp_File_1e.txt", header = TRUE, stringsAsFactors = FALSE)

### Calculate B-H estimate (FDR) for all 37 disease outcomes (including 5 in Martin et al. (2021))
for (FDR in c(0.1, 0.05)) {
  # For BMI (all disease outcomes)
  mydata$rank <- rank(mydata$BMI.p)
  mydata$BH <- (mydata$rank/nrow(mydata))*FDR
  mydata$diff <- mydata$BH - mydata$BMI.p
  maxP <- max(mydata[mydata$diff > 0,]$BMI.p)
  mydata$BMI.BH_sig <- 0
  mydata[mydata$BMI.p <= maxP,]$BMI.BH_sig <- 1
  
  # For BFP/FA/UFA (only BMI-significant disease outcomes)
  bmidata <- mydata[mydata$BMI.p < 0.05,]
  
  bmidata$rank <- rank(bmidata$BFP.p)
  bmidata$BH <- (bmidata$rank/nrow(bmidata))*FDR
  bmidata$diff <- bmidata$BH - bmidata$BFP.p
  maxP <- max(bmidata[bmidata$diff > 0,]$BFP.p)
  bmidata$BFP.BH_sig <- 0
  bmidata[bmidata$BFP.p <= maxP,]$BFP.BH_sig <- 1
  
  bmidata$rank <- rank(bmidata$FA.p)
  bmidata$BH <- (bmidata$rank/nrow(bmidata))*FDR
  bmidata$diff <- bmidata$BH - bmidata$FA.p
  maxP <- max(bmidata[bmidata$diff > 0,]$FA.p)
  bmidata$FA.BH_sig <- 0
  bmidata[bmidata$FA.p <= maxP,]$FA.BH_sig <- 1
  
  bmidata$rank <- rank(bmidata$UFA.p)
  bmidata$BH <- (bmidata$rank/nrow(bmidata))*FDR
  bmidata$diff <- bmidata$BH - bmidata$UFA.p
  maxP <- max(bmidata[bmidata$diff > 0,]$UFA.p)
  bmidata$UFA.BH_sig <- 0
  bmidata[bmidata$UFA.p <= maxP,]$UFA.BH_sig <- 1
  
  bmidata <- bmidata[c("Outcome", "BFP.BH_sig", "FA.BH_sig", "UFA.BH_sig")]
  mydata <- merge(mydata, bmidata, all.x = TRUE)
  
  colnames(mydata)[colnames(mydata) == "BMI.BH_sig"] <- paste0("BMI.BH_sig.", FDR)
  colnames(mydata)[colnames(mydata) == "BFP.BH_sig"] <- paste0("BFP.BH_sig.", FDR)
  colnames(mydata)[colnames(mydata) == "FA.BH_sig"] <- paste0("FA.BH_sig.", FDR)
  colnames(mydata)[colnames(mydata) == "UFA.BH_sig"] <- paste0("UFA.BH_sig.", FDR)
}

# Print results
print("All 37 diseases - FDR = 0.1")
print(paste0("BMI: ", nrow(mydata[!is.na(mydata$BMI.BH_sig.0.1) & mydata$BMI.BH_sig.0.1 == 1,]), " out of ", nrow(mydata[!is.na(mydata$BMI.BH_sig.0.1),])))
print(paste0("BFP: ", nrow(mydata[!is.na(mydata$BFP.BH_sig.0.1) & mydata$BFP.BH_sig.0.1 == 1,]), " out of ", nrow(mydata[!is.na(mydata$BFP.BH_sig.0.1),])))
print(paste0("FA: ", nrow(mydata[!is.na(mydata$FA.BH_sig.0.1) & mydata$FA.BH_sig.0.1 == 1,]), " out of ", nrow(mydata[!is.na(mydata$FA.BH_sig.0.1),])))
print(paste0("UFA: ", nrow(mydata[!is.na(mydata$UFA.BH_sig.0.1) & mydata$UFA.BH_sig.0.1 == 1,]), " out of ", nrow(mydata[!is.na(mydata$UFA.BH_sig.0.1),])))

print("All 37 diseases - FDR = 0.05")
print(paste0("BMI: ", nrow(mydata[!is.na(mydata$BMI.BH_sig.0.05) & mydata$BMI.BH_sig.0.05 == 1,]), " out of ", nrow(mydata[!is.na(mydata$BMI.BH_sig.0.05),])))
print(paste0("BFP: ", nrow(mydata[!is.na(mydata$BFP.BH_sig.0.05) & mydata$BFP.BH_sig.0.05 == 1,]), " out of ", nrow(mydata[!is.na(mydata$BFP.BH_sig.0.05),])))
print(paste0("FA: ", nrow(mydata[!is.na(mydata$FA.BH_sig.0.05) & mydata$FA.BH_sig.0.05 == 1,]), " out of ", nrow(mydata[!is.na(mydata$FA.BH_sig.0.05),])))
print(paste0("UFA: ", nrow(mydata[!is.na(mydata$UFA.BH_sig.0.05) & mydata$UFA.BH_sig.0.05 == 1,]), " out of ", nrow(mydata[!is.na(mydata$UFA.BH_sig.0.05),])))

### Calculate B-H estimate (FDR) for 32 disease outcomes (not including 5 in Martin et al. (2021))
mydata <- mydata[!grepl("BH.sig", colnames(mydata))]
mydata <- mydata[!mydata$Outcome %in% c("Type 2 diabetes", "Coronary artery disease", "Hypertension", "Stroke", "Polycystic ovary syndrome"),]

for (FDR in c(0.1, 0.05)) {
  # For BMI (all traits)
  mydata$rank <- rank(mydata$BMI.p)
  mydata$BH <- (mydata$rank/nrow(mydata))*FDR
  mydata$diff <- mydata$BH - mydata$BMI.p
  maxP <- max(mydata[mydata$diff > 0,]$BMI.p)
  mydata$BMI.BH_sig <- 0
  mydata[mydata$BMI.p <= maxP,]$BMI.BH_sig <- 1
  
  # For BFP/FA/UFA (only BMI-significant traits)
  bmidata <- mydata[mydata$BMI.p < 0.05,]
  
  bmidata$rank <- rank(bmidata$BFP.p)
  bmidata$BH <- (bmidata$rank/nrow(bmidata))*FDR
  bmidata$diff <- bmidata$BH - bmidata$BFP.p
  maxP <- max(bmidata[bmidata$diff > 0,]$BFP.p)
  bmidata$BFP.BH_sig <- 0
  bmidata[bmidata$BFP.p <= maxP,]$BFP.BH_sig <- 1
  
  bmidata$rank <- rank(bmidata$FA.p)
  bmidata$BH <- (bmidata$rank/nrow(bmidata))*FDR
  bmidata$diff <- bmidata$BH - bmidata$FA.p
  maxP <- max(bmidata[bmidata$diff > 0,]$FA.p)
  bmidata$FA.BH_sig <- 0
  bmidata[bmidata$FA.p <= maxP,]$FA.BH_sig <- 1
  
  bmidata$rank <- rank(bmidata$UFA.p)
  bmidata$BH <- (bmidata$rank/nrow(bmidata))*FDR
  bmidata$diff <- bmidata$BH - bmidata$UFA.p
  maxP <- max(bmidata[bmidata$diff > 0,]$UFA.p)
  bmidata$UFA.BH_sig <- 0
  bmidata[bmidata$UFA.p <= maxP,]$UFA.BH_sig <- 1
  
  bmidata <- bmidata[c("Outcome", "BFP.BH_sig", "FA.BH_sig", "UFA.BH_sig")]
  mydata <- merge(mydata, bmidata, all.x = TRUE)
  
  colnames(mydata)[colnames(mydata) == "BMI.BH_sig"] <- paste0("BMI.BH_sig.", FDR)
  colnames(mydata)[colnames(mydata) == "BFP.BH_sig"] <- paste0("BFP.BH_sig.", FDR)
  colnames(mydata)[colnames(mydata) == "FA.BH_sig"] <- paste0("FA.BH_sig.", FDR)
  colnames(mydata)[colnames(mydata) == "UFA.BH_sig"] <- paste0("UFA.BH_sig.", FDR)
}

# Print results
print("32 diseases not in Martin et al. (2021) - FDR = 0.1")
print(paste0("BMI: ", nrow(mydata[!is.na(mydata$BMI.BH_sig.0.1) & mydata$BMI.BH_sig.0.1 == 1,]), " out of ", nrow(mydata[!is.na(mydata$BMI.BH_sig.0.1),])))
print(paste0("BFP: ", nrow(mydata[!is.na(mydata$BFP.BH_sig.0.1) & mydata$BFP.BH_sig.0.1 == 1,]), " out of ", nrow(mydata[!is.na(mydata$BFP.BH_sig.0.1),])))
print(paste0("FA: ", nrow(mydata[!is.na(mydata$FA.BH_sig.0.1) & mydata$FA.BH_sig.0.1 == 1,]), " out of ", nrow(mydata[!is.na(mydata$FA.BH_sig.0.1),])))
print(paste0("UFA: ", nrow(mydata[!is.na(mydata$UFA.BH_sig.0.1) & mydata$UFA.BH_sig.0.1 == 1,]), " out of ", nrow(mydata[!is.na(mydata$UFA.BH_sig.0.1),])))

print("32 diseases not in Martin et al. (2021) - FDR = 0.05")
print(paste0("BMI: ", nrow(mydata[!is.na(mydata$BMI.BH_sig.0.05) & mydata$BMI.BH_sig.0.05 == 1,]), " out of ", nrow(mydata[!is.na(mydata$BMI.BH_sig.0.05),])))
print(paste0("BFP: ", nrow(mydata[!is.na(mydata$BFP.BH_sig.0.05) & mydata$BFP.BH_sig.0.05 == 1,]), " out of ", nrow(mydata[!is.na(mydata$BFP.BH_sig.0.05),])))
print(paste0("FA: ", nrow(mydata[!is.na(mydata$FA.BH_sig.0.05) & mydata$FA.BH_sig.0.05 == 1,]), " out of ", nrow(mydata[!is.na(mydata$FA.BH_sig.0.05),])))
print(paste0("UFA: ", nrow(mydata[!is.na(mydata$UFA.BH_sig.0.05) & mydata$UFA.BH_sig.0.05 == 1,]), " out of ", nrow(mydata[!is.na(mydata$UFA.BH_sig.0.05),])))

