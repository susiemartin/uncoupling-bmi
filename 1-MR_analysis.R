# PART 2 - Run Mendelian Randomisation for each disease and sub-disease
# using IVW, weighted median, MR-Egger and penalised weighted median methods
## Generates Supplementary File 1g(i-iv)

# Define input arguments
args <- commandArgs(TRUE)

test_trait <- args[1] #name of exposure (BMI, BFP, FA or UFA)
snp_info <- args[2] #location of SNP INFO file
summ_stats_file <- args[3] #location of file containing individual BetaYG and seBetaYG for the disease outcomes

# Load required packages
library(MendelianRandomization)

# Load SNP INFO and disease outcome data
snp_info <- read.table(snp_info, header = TRUE, sep = "\t", strip.white = TRUE)
snp_info_specific <- snp_info[which(snp_info$Test_trait == test_trait),]
results <- read.table(summ_stats_file, header = TRUE, sep = "\t", strip.white = TRUE)

# Merge the snp_info and results file and flip to trait-raising allele
stats <- merge(results, snp_info_specific, by = "SNP")
stats$BetaYG <- ifelse(as.character(stats$Trait_raising) == as.character(stats$A1), stats$BetaYG, -stats$BetaYG)

# Conversion function
S4_to_dataframe <- function(s4obj) {
  nms <- slotNames(s4obj)
  
  lst <- lapply(nms, function(nm) slot(s4obj, nm))
  as.data.frame(setNames(lst, nms))
}

# Run MR analysis per outcome trait
allivwdata <- NULL
allmeddata <- NULL
alleggdata <- NULL
allpendata <- NULL

for (trait in unique(stats$Trait)) {
  substats <- stats[stats$Trait == trait,]
  # Create MR input object
  MRInputObject <- mr_input(bx = substats$BetaXG, bxse = substats$seBetaXG, by = substats$BetaYG, byse = substats$seBetaYG,
                            exposure = test_trait, outcome = trait)
  
  # IVW method and heterogeneity test
  ivwdata <- mr_ivw(MRInputObject)
  ivwdata <- S4_to_dataframe(ivwdata)
  allivwdata <- rbind(allivwdata, ivwdata)
  
  # Weighted median method
  meddata <- mr_median(MRInputObject, weighting = "weighted")
  meddata <- S4_to_dataframe(meddata)
  allmeddata <- rbind(allmeddata, meddata)
  
  # MR-Egger method and heterogeneity test
  eggdata <- mr_egger(MRInputObject)
  eggdata <- S4_to_dataframe(eggdata)
  alleggdata <- rbind(alleggdata, eggdata)
  
  # Penalised weighted median method
  pendata <- mr_median(MRInputObject, weighting = "penalized")
  pendata <- S4_to_dataframe(pendata)
  allpendata <- rbind(allpendata, pendata)
}

# Create table of heterogeneity statistics
savedata <- as.data.frame(matrix(nrow = length(unique(stats$Trait))*4, ncol = 13, NA))
colnames(savedata) <- c("Exposure", "Trait", "Analysis", "Beta", "SE", "lowerCI", "upperCI", "P", "Egger_intercept", "Intercept_P", "het_Q", "het_P", "I2_Egger")
savedata$Exposure <- test_trait
savedata$Trait <- rep(unique(stats$Trait), each = 4)
savedata$Analysis <- c("IVW", "Weighted_median", "MR_Egger", "Penalised_weighted_median")

for (i in unique(savedata$Trait)) {
  savedata[savedata$Trait == i & savedata$Analysis == "IVW",]$Beta <- allivwdata[allivwdata$Outcome == i,]$Estimate[1]
  savedata[savedata$Trait == i & savedata$Analysis == "IVW",]$SE <- allivwdata[allivwdata$Outcome == i,]$StdError[1]
  savedata[savedata$Trait == i & savedata$Analysis == "IVW",]$lowerCI <- allivwdata[allivwdata$Outcome == i,]$CILower[1]
  savedata[savedata$Trait == i & savedata$Analysis == "IVW",]$upperCI <- allivwdata[allivwdata$Outcome == i,]$CIUpper[1]
  savedata[savedata$Trait == i & savedata$Analysis == "IVW",]$P <- allivwdata[allivwdata$Outcome == i,]$Pvalue[1]
  savedata[savedata$Trait == i & savedata$Analysis == "IVW",]$het_Q <- allivwdata[allivwdata$Outcome == i,]$Heter.Stat[1]
  savedata[savedata$Trait == i & savedata$Analysis == "IVW",]$het_P <- allivwdata[allivwdata$Outcome == i,]$Heter.Stat[2]
  
  savedata[savedata$Trait == i & savedata$Analysis == "Weighted_median",]$Beta <- allmeddata[allmeddata$Outcome == i,]$Estimate
  savedata[savedata$Trait == i & savedata$Analysis == "Weighted_median",]$SE <- allmeddata[allmeddata$Outcome == i,]$StdError
  savedata[savedata$Trait == i & savedata$Analysis == "Weighted_median",]$lowerCI <- allmeddata[allmeddata$Outcome == i,]$CILower
  savedata[savedata$Trait == i & savedata$Analysis == "Weighted_median",]$upperCI <- allmeddata[allmeddata$Outcome == i,]$CIUpper
  savedata[savedata$Trait == i & savedata$Analysis == "Weighted_median",]$P <- allmeddata[allmeddata$Outcome == i,]$Pvalue
  
  savedata[savedata$Trait == i & savedata$Analysis == "MR_Egger",]$Beta <- alleggdata[alleggdata$Outcome == i,]$Estimate[1]
  savedata[savedata$Trait == i & savedata$Analysis == "MR_Egger",]$SE <- alleggdata[alleggdata$Outcome == i,]$StdError[1]
  savedata[savedata$Trait == i & savedata$Analysis == "MR_Egger",]$lowerCI <- alleggdata[alleggdata$Outcome == i,]$CILower[1]
  savedata[savedata$Trait == i & savedata$Analysis == "MR_Egger",]$upperCI <- alleggdata[alleggdata$Outcome == i,]$CIUpper[1]
  savedata[savedata$Trait == i & savedata$Analysis == "MR_Egger",]$P <- alleggdata[alleggdata$Outcome == i,]$Pvalue[1]
  savedata[savedata$Trait == i & savedata$Analysis == "MR_Egger",]$Egger_intercept <- alleggdata[alleggdata$Outcome == i,]$Intercept[1]
  savedata[savedata$Trait == i & savedata$Analysis == "MR_Egger",]$Intercept_P <- alleggdata[alleggdata$Outcome == i,]$Pvalue.Int[1]
  savedata[savedata$Trait == i & savedata$Analysis == "MR_Egger",]$het_Q <- alleggdata[alleggdata$Outcome == i,]$Heter.Stat[1]
  savedata[savedata$Trait == i & savedata$Analysis == "MR_Egger",]$het_P <- alleggdata[alleggdata$Outcome == i,]$Heter.Stat[2]
  savedata[savedata$Trait == i & savedata$Analysis == "MR_Egger",]$I2_Egger <- alleggdata[alleggdata$Outcome == i,]$I.sq[1]
  
  savedata[savedata$Trait == i & savedata$Analysis == "Penalised_weighted_median",]$Beta <- allpendata[allpendata$Outcome == i,]$Estimate
  savedata[savedata$Trait == i & savedata$Analysis == "Penalised_weighted_median",]$SE <- allpendata[allpendata$Outcome == i,]$StdError
  savedata[savedata$Trait == i & savedata$Analysis == "Penalised_weighted_median",]$lowerCI <- allpendata[allpendata$Outcome == i,]$CILower
  savedata[savedata$Trait == i & savedata$Analysis == "Penalised_weighted_median",]$upperCI <- allpendata[allpendata$Outcome == i,]$CIUpper
  savedata[savedata$Trait == i & savedata$Analysis == "Penalised_weighted_median",]$P <- allpendata[allpendata$Outcome == i,]$Pvalue
}

write.table(savedata, file = paste0("MR_results_", test_trait, ".txt"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# Add study origin label
savedata$Study <- "PublishedGWAS"
if (nrow(savedata[grepl("_FG", savedata$Trait),]) > 0) {
  savedata[grepl("_FG", savedata$Trait),]$Study <- "FinnGen"
}
if (nrow(savedata[grepl("_UKB", savedata$Trait),]) > 0) {
  savedata[grepl("_UKB", savedata$Trait),]$Study <- "UKBiobank"
}
# Remove UK Biobank study results from this table
savedata <- savedata[savedata$Study != "UKBiobank",]

# Add system column (including to sub-diseases here)
savedata$System <- NA
savedata[savedata$Trait %in% c("Type 2 diabetes", "Hypertension", "PCOS", "CAD", "Stroke",
                                     "PAD", "HeartFailure", "AtrialFibrillation", "CKD", "VTE",
                                     "DVT", "PulmonaryEmbolism", "AbdominalAneurysm"),]$System <- "CardiovascularAndMetabolic"
savedata[grepl("Stroke", savedata$Trait),]$System <- "CardiovascularAndMetabolic"
savedata[savedata$Trait %in% c("Gout", "Osteoarthritis", "Osteoporosis", "RheumatoidArthritis"),]$System <- "Musculoskeletal"
savedata[grepl("Osteoarthritis", savedata$Trait),]$System <- "Musculoskeletal"
savedata[savedata$Trait %in% c("Gallstones", "GORD"),]$System <- "Gastrointestinal"
savedata[savedata$Trait %in% c("Alzheimers", "Depression", "MultipleSclerosis", "Parkinsons"),]$System <- "Nervous"
savedata[savedata$Trait %in% c("Psoriasis"),]$System <- "Integumentary"
savedata[savedata$Trait %in% c("Adult-onsetAsthma", "Child-onsetAsthma"),]$System <- "Respiratory"
savedata[savedata$Trait %in% c("Barretts", "BreastCancer", "Myeloma", "ColorectalCancer",
                                     "EndometrialCancer", "LungCancer", "Meningioma", "OvarianCancer",
                                     "PancreaticCancer", "ProstateCancer", "RenalCancer", "ThyroidCancer"),]$System <- "Cancer"
savedata[grepl("Cancer", savedata$Trait),]$System <- "Cancer"

# Calculate odds ratios and associated 95% confidence intervals
savedata$OR <- exp(savedata$Beta)
savedata$ORlci <- exp(savedata$lowerCI)
savedata$ORuci <- exp(savedata$upperCI)
savedata <- savedata[c("System", "Trait", "Study", "Analysis", "OR", "ORlci", "ORuci", "P", "Egger_intercept", "Intercept_P", "het_Q", "het_P", "I2_Egger")]

savedata <- savedata[order(savedata$System, savedata$Trait, savedata$Study),]

write.table(savedata, file = paste0("Supp_file_1_g_", test_trait, ".txt"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

