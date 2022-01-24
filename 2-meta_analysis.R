# PART 3 - Meta-analyse published GWAS and FinnGen IVW MR results together
# for appropriate matching outcomes (DO NOT include UK Biobank traits)
## Generates Supplementary File 1e, 1f and 1h

# Load required packages
library(metafor)

# Read in file of diseases to be meta-analysed together (with meta-analysed trait name)
# - do not include UK Biobank traits
# e.g. Trait     Meta_name
#      CAD_GWAS  CAD_meta
#      CAD_FG    CAD_meta
#      PCOS_GWAS PCOS_meta
#      ...       ...
metanames <- read.table("disease_meta_analysis.txt", header = TRUE, stringsAsFactors = FALSE)

# Read in MR results output produced in PART 1 script
bmidata <- read.table("MR_results_BMI.txt", header = TRUE, stringsAsFactors = FALSE)
bmidata$cluster <- "BMI"

bfpdata <- read.table("MR_results_BFP.txt", header = TRUE, stringsAsFactors = FALSE)
bfpdata$cluster <- "BFP"

fadata <- read.table("MR_results_FA.txt", header = TRUE, stringsAsFactors = FALSE)
fadata$cluster <- "FA"

ufadata <- read.table("MR_results_UFA.txt", header = TRUE, stringsAsFactors = FALSE)
ufadata$cluster <- "UFA"

alldata <- rbind(bmidata, bfpdata, fadata, ufadata)

# Merge trait results with meta-analysis name
metadata <- merge(alldata, metanames, by = "Trait")

### Run non-UK Biobank IVW MR meta-analyses
# Make blank dataframe to fill with meta-analysis results
metacols <- c("OR", "se", "ORlci", "ORuci", "p", "Q", "I2", "hetP", "df", "tau2")
myresults <- data.frame(matrix(nrow = length(unique(metadata$Meta_name)), ncol = (length(metacols)*4)+1))
colnames(myresults) <- c("Outcome", paste0("BMI.", metacols), paste0("BFP.", metacols), paste0("FA.", metacols), paste0("UFA.", metacols))
myresults$Outcome <- unique(metadata$Meta_name)

for (i in unique(metadata$Meta_name)) {
  for (j in c("BMI", "BFP", "FA", "UFA")) {
    subdata <- metadata[metadata$Meta_name == i & metadata$cluster == j & metadata$Analysis == "IVW",]
    meta <- rma(yi = Beta, sei = SE, data = subdata, method = "REML")
    myresults[myresults$Outcome == i, colnames(myresults) == paste0(j, ".OR")] <- exp(meta$beta[,1])
    myresults[myresults$Outcome == i, colnames(myresults) == paste0(j, ".se")] <- meta$se
    myresults[myresults$Outcome == i, colnames(myresults) == paste0(j, ".p")] <- meta$pval
    myresults[myresults$Outcome == i, colnames(myresults) == paste0(j, ".Q")] <- meta$QE
    myresults[myresults$Outcome == i, colnames(myresults) == paste0(j, ".I2")] <- meta$I2
    myresults[myresults$Outcome == i, colnames(myresults) == paste0(j, ".hetP")] <- meta$QEp
    myresults[myresults$Outcome == i, colnames(myresults) == paste0(j, ".df")] <- meta$k - 1
    myresults[myresults$Outcome == i, colnames(myresults) == paste0(j, ".tau2")] <- meta$tau2
    
    myresults[myresults$Outcome == i, colnames(myresults) == paste0(j, ".ORlci")] <- exp(meta$beta[,1] - (qnorm(0.975)*meta$se))
    myresults[myresults$Outcome == i, colnames(myresults) == paste0(j, ".ORuci")] <- exp(meta$beta[,1] + (qnorm(0.975)*meta$se))
  }
}

myresults$Study <- "Meta-analysis"
metaresults <- myresults

# Add system column
metaresults$System <- NA
metaresults[metaresults$Outcome %in% c("Type 2 diabetes", "Hypertension", "PCOS", "CAD", "Stroke",
                                     "PAD", "HeartFailure", "AtrialFibrillation", "CKD", "VTE",
                                     "DVT", "PulmonaryEmbolism", "AbdominalAneurysm"),]$System <- "CardiovascularAndMetabolic"
metaresults[metaresults$Outcome %in% c("Gout", "Osteoarthritis", "Osteoporosis", "RheumatoidArthritis"),]$System <- "Musculoskeletal"
metaresults[metaresults$Outcome %in% c("Gallstones", "GORD"),]$System <- "Gastrointestinal"
metaresults[metaresults$Outcome %in% c("Alzheimers", "Depression", "MultipleSclerosis", "Parkinsons"),]$System <- "Nervous"
metaresults[metaresults$Outcome %in% c("Psoriasis"),]$System <- "Integumentary"
metaresults[metaresults$Outcome %in% c("Adult-onsetAsthma"),]$System <- "Respiratory"
metaresults[metaresults$Outcome %in% c("Barretts", "BreastCancer", "Myeloma", "ColorectalCancer",
                                     "EndometrialCancer", "LungCancer", "Meningioma", "OvarianCancer",
                                     "PancreaticCancer", "ProstateCancer", "RenalCancer", "ThyroidCancer"),]$System <- "Cancer"

metaresults <- metaresults[c("System", "Outcome", "BMI.Q", "BMI.hetP", "BMI.I2", "BFP.Q", "BFP.hetP", "BFP.I2",
                             "FA.Q", "FA.hetP", "FA.I2", "UFA.Q", "UFA.hetP", "UFA.I2")]
metaresults <- metaresults[order(metaresults$System, metaresults$Outcome),]
write.table(metaresults, file = "Supp_File_1f.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# Combine meta-analysis results with individual disease outcome results with no meta-analysis
# (i.e. those only available in one published GWAS or FinnGen)
singdata <- alldata[!(alldata$Trait %in% metanames$Trait),]
singdata <- singdata[singdata$Analysis == "IVW",]

# Calculate odds ratios and associated 95% confidence intervals
singdata$OR <- exp(singdata$Beta)
singdata$ORlci <- exp(singdata$lowerCI)
singdata$ORuci <- exp(singdata$upperCI)
singdata <- singdata[c("cluster", "Trait", "OR", "ORlci", "ORuci", "P")]

allsingdata <- data.frame(unique(singdata$Trait))
colnames(allsingdata) <- "Outcome"
for (i in c("BMI", "BFP", "FA", "UFA")) {
  subsingdata <- singdata[singdata$cluster == i,]
  subsingdata <- subsingdata[,-1]
  colnames(subsingdata) <- c("Outcome", paste0(i, ".", c("OR", "ORlci", "ORuci", "p")))
  allsingdata <- merge(allsingdata, subsingdata, by = "Outcome")
}

allsingdata$Study <- "PublishedGWAS"
if (nrow(allsingdata[grepl("_FG", allsingdata$Outcome),]) > 0) {
  allsingdata[grepl("_FG", allsingdata$Outcome),]$Study <- "FinnGen"
}
if (nrow(allsingdata[grepl("_UKB", allsingdata$Outcome),]) > 0) {
  allsingdata[grepl("_UKB", allsingdata$Outcome),]$Study <- "UKBiobank"
}

myresults <- myresults[colnames(allsingdata)]
allresults <- rbind(myresults, allsingdata)

# Add system column
allresults$System <- NA
allresults[allresults$Outcome %in% c("Type 2 diabetes", "Hypertension", "PCOS", "CAD", "Stroke",
                                   "PAD", "HeartFailure", "AtrialFibrillation", "CKD", "VTE",
                                   "DVT", "PulmonaryEmbolism", "AbdominalAneurysm"),]$System <- "CardiovascularAndMetabolic"
allresults[allresults$Outcome %in% c("Gout", "Osteoarthritis", "Osteoporosis", "RheumatoidArthritis"),]$System <- "Musculoskeletal"
allresults[allresults$Outcome %in% c("Gallstones", "GORD"),]$System <- "Gastrointestinal"
allresults[allresults$Outcome %in% c("Alzheimers", "Depression", "MultipleSclerosis", "Parkinsons"),]$System <- "Nervous"
allresults[allresults$Outcome %in% c("Psoriasis"),]$System <- "Integumentary"
allresults[allresults$Outcome %in% c("Adult-onsetAsthma"),]$System <- "Respiratory"
allresults[allresults$Outcome %in% c("Barretts", "BreastCancer", "Myeloma", "ColorectalCancer",
                                   "EndometrialCancer", "LungCancer", "Meningioma", "OvarianCancer",
                                   "PancreaticCancer", "ProstateCancer", "RenalCancer", "ThyroidCancer"),]$System <- "Cancer"

allresults <- allresults[,c(ncol(allresults),1,(ncol(allresults)-1),2:(ncol(allresults)-2))]

nonukbresults <- allresults[allresults$Study != "UKBiobank",]
nonukbresults <- nonukbresults[order(nonukbresults$System, nonukbresults$Outcome, nonukbresults$Study),]
write.table(nonukbresults, file = "Supp_File_1e.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

ukbresults <- allresults[allresults$Study == "UKBiobank",]
ukbresults <- ukbresults[,-3]
ukbresults <- ukbresults[order(ukbresults$System, ukbresults$Outcome),]
write.table(ukbresults, file = "Supp_File_1h.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

## Find the total number of diseases with BMI p<0.05
print(paste0("Total number of diseases with BMI IVW p<0.05: ", nrow(nonukbresults[nonukbresults$BMI.p < 0.05,])))
## Find the total number of diseases with BFP p<0.05 (out of those with BMI p<0.05)
bmisig <- nonukbresults[nonukbresults$BMI.p < 0.05,]
print(paste0("Total number of diseases with BFP IVW also p<0.05: ", nrow(bmisig[bmisig$BFP.p < 0.05,])))

## Find the total number of IVW tests out of 37*4 with p<0.05
ivwp <- c(nonukbresults$BMI.p, nonukbresults$BFP.p, nonukbresults$FA.p, nonukbresults$UFA.p)
print(paste0("Total number of IVW tests with p<0.05: ", length(ivwp[ivwp < 0.05])))
print(paste0("Number expected by chance to have p<0.05 out of 37*4: ", 37*4*0.05))

## Sensitivity analyses (including disease sub-types)
alldata$Study <- "PublishedGWAS"
if (nrow(alldata[grepl("_FG", alldata$Trait),]) > 0) {
  alldata[grepl("_FG", alldata$Trait),]$Study <- "FinnGen"
}
if (nrow(alldata[grepl("_UKB", alldata$Trait),]) > 0) {
  alldata[grepl("_UKB", alldata$Trait),]$Study <- "UKBiobank"
}
# Find the number of weighted median MR results directionally consistent with IVW MR results
for (i in c("BMI", "BFP", "FA", "UFA")) {
  subdata <- alldata[alldata$cluster == i,]
  subdata <- subdata[!grepl("_UKB", subdata$Trait),]
  subdata$Label <- paste0(subdata$Trait, ".", subdata$Study)
  
  ivwdata <- subdata[subdata$Analysis == "IVW",]
  ivwdata <- ivwdata[c("Label", "Beta", "P")]
  colnames(ivwdata)[2:3] <- paste0("IVW.", colnames(ivwdata)[2:3])
  
  meddata <- subdata[subdata$Analysis == "Weighted_median",]
  meddata <- meddata[c("Label", "Beta", "P")]
  colnames(meddata)[2:3] <- paste0("Med.", colnames(meddata)[2:3])
  
  subdata <- merge(ivwdata, meddata, by = "Label")
  ivwmed <- nrow(subdata[(subdata$IVW.Beta > 0 & subdata$Med.Beta > 0) | (subdata$IVW.Beta < 0 & subdata$Med.Beta < 0),])
  print(paste0(i, ": ", ivwmed, " of out ", nrow(subdata), " diseases directionally consistent between IVW and weighted median"))
  
  # Find the number of weighted median MR results (with p<0.05) directionally consistent with IVW MR results
  subdata <- subdata[subdata$Med.P < 0.05,]
  ivwmed <- nrow(subdata[(subdata$IVW.Beta > 0 & subdata$Med.Beta > 0) | (subdata$IVW.Beta < 0 & subdata$Med.Beta < 0),])
  print(paste0(i, ": ", ivwmed, " of out ", nrow(subdata), " diseases directionally consistent between IVW and weighted median (p<0.05)"))
}

# Find the number of MR-Egger results directionally consistent with IVW MR results
for (i in c("BMI", "BFP", "FA", "UFA")) {
  subdata <- alldata[alldata$cluster == i,]
  subdata <- subdata[!grepl("_UKB", subdata$Trait),]
  subdata$Label <- paste0(subdata$Trait, ".", subdata$Study)
  
  ivwdata <- subdata[subdata$Analysis == "IVW",]
  ivwdata <- ivwdata[c("Label", "Beta", "P")]
  colnames(ivwdata)[2:3] <- paste0("IVW.", colnames(ivwdata)[2:3])
  
  eggdata <- subdata[subdata$Analysis == "MR_Egger",]
  eggdata <- eggdata[c("Label", "Beta", "P")]
  colnames(eggdata)[2:3] <- paste0("Egg.", colnames(eggdata)[2:3])
  
  subdata <- merge(ivwdata, eggdata, by = "Label")
  ivwegg <- nrow(subdata[(subdata$IVW.Beta > 0 & subdata$Egg.Beta > 0) | (subdata$IVW.Beta < 0 & subdata$Egg.Beta < 0),])
  print(paste0(i, ": ", ivwegg, " of out ", nrow(subdata), " diseases directionally consistent between IVW and MR-Egger"))
  
  # Find the number of MR-Egger results (with p<0.05) directionally consistent with IVW MR results
  subdata <- subdata[subdata$Egg.P < 0.05,]
  ivwegg <- nrow(subdata[(subdata$IVW.Beta > 0 & subdata$Egg.Beta > 0) | (subdata$IVW.Beta < 0 & subdata$Egg.Beta < 0),])
  print(paste0(i, ": ", ivwegg, " of out ", nrow(subdata), " diseases directionally consistent between IVW and MR-Egger (p<0.05)"))
}

# Find the number of UK Biobank IVW results directionally consistent with published GWAS and/or FinnGen IVW results
for (i in c("BMI", "BFP", "FA", "UFA")) {
  subnonukb <- nonukbresults[c("Outcome", paste0(i, c(".OR", ".p")))]
  colnames(subnonukb)[2:3] <- paste0("Non.", c("OR", "p"))
  
  subukb <- ukbresults[c("Outcome", paste0(i, c(".OR", ".p")))]
  colnames(subukb)[2:3] <- paste0("UKB.", c("OR", "p"))
  
  subdata <- merge(subnonukb, subukb, by = "Outcome")
  ukbnonukb <- nrow(subdata[(subdata$Non.OR > 1 & subdata$UKB.OR > 1) | (subdata$Non.OR < 1 & subdata$UKB.OR < 1),])
  print(paste0(i, ": ", ukbnonukb, " of out ", nrow(subdata), " diseases directionally consistent between published GWAS/FinnGen and UKBiobank"))
  
  # Find the number of UK Biobank IVW results (with p<0.05) directionally consistent with published GWAS and/or FinnGen IVW results
  subdata <- subdata[subdata$UKB.p < 0.05,]
  ukbnonukb <- nrow(subdata[(subdata$Non.OR > 1 & subdata$UKB.OR > 1) | (subdata$Non.OR < 1 & subdata$UKB.OR < 1),])
  print(paste0(i, ": ", ukbnonukb, " of out ", nrow(subdata), " diseases directionally consistent between between published GWAS/FinnGen and UK Biobank (p<0.05)"))
}

