library(openxlsx)
library(gplots)
library(dplyr)

setwd('P:/AG_Hertel/Fluxomics/Microbiome_Modelling/MetaPhlAn4')

# read net secretions, direct contributions of the predictMicrobeContributions function and normalized abundances
dnetsecr <- read.csv('net_secretion_MPH4.csv')
dsecr <- read.csv('predictContributions_MPH4.csv')
dnormcov <- read.csv('normalizedCoverage_MPH4.csv')
# Data cleaning
dsecr <- select(dsecr, -contains("_"), -contains("."))
dnormcov <- select(dnormcov, -contains("_"), -contains("."))
interscolnames <- intersect(colnames(dsecr), colnames(dnormcov))

dsecr <- dsecr[, c("X", interscolnames)]
dnormcov <- dnormcov[, c("ID", interscolnames)]

smplsrem <- colnames(dsecr)[colSums(dsecr != 0) == 0]


dsecr <- dsecr[, !names(dsecr) %in% smplsrem]
dnormcov <- dnormcov[, !names(dnormcov) %in% smplsrem]
dnetsecr <- dnetsecr[, c("Net.secretion", colnames(dsecr)[2:ncol(dsecr)])]

# readout metabolites and species

dnetsecr <- dnetsecr[rowSums(dnetsecr[, -1] != 0) > 0, ]
mets <- gsub("^EX_|\\[fe\\]$", "", dnetsecr$Net.secretion)


panspecies <- dnormcov$ID
species <- gsub("^pan", "", panspecies)

effects <- data.frame(matrix(ncol = length(species) * 9,
                             nrow = length(mets)))

colnames_new <- unlist(lapply(species, function(species_name) {
  c(paste("e", species_name, sep = "_"),
    paste("e_", species_name, ".lb", sep = ""),
    paste("e_", species_name, ".ub", sep = ""),
    paste("d", species_name, sep = "_"),
    paste("d_", species_name, ".lb", sep = ""),
    paste("d_", species_name, ".ub", sep = ""),
    paste("t", species_name, sep = "_"),
    paste("t_", species_name, ".lb", sep = ""),
    paste("t_", species_name, ".ub", sep = ""))
}))

colnames(effects) <- colnames_new

for (l in seq_along(mets)) {
  for (j in seq_along(panspecies)) {
    
    rows <- grepl(paste0("_", mets[l], "$"), dsecr$X)
    numeric_cols <- sapply(dsecr, is.numeric)
    sums <- colSums(dsecr[rows, numeric_cols], na.rm = TRUE)
    
    df <- data.frame(
      ID = names(sums),
      val = sums
    )
    colnames(df)[2] <- mets[l]
    
    subspec <- dnormcov[dnormcov$ID == panspecies[j], -1]
    names_Mj <- colnames(subspec)[which(subspec > 0)]
    subspec_bin <- as.numeric(subspec[1, ] > 0)
    
    subspec_long <- data.frame(
      ID = colnames(subspec),
      pres = subspec_bin
    )
    colnames(subspec_long)[2] <- panspecies[j]
    
    df <- merge(df, subspec_long, by = "ID")
    
    ########## Ecological effect
    rows_lj <- grepl(paste0("_", mets[l], "$"), dsecr$X) &
      !grepl(paste0("^", species[j]), dsecr$X)
    
    sums_eco <- colSums(dsecr[rows_lj, numeric_cols], na.rm = TRUE)
    
    df_eco <- data.frame(
      ID = names(sums_eco),
      val = sums_eco
    )
    colnames(df_eco)[2] <- mets[l]
    df_eco <- merge(df_eco, subspec_long, by = "ID")
    
    mod_eco <- lm(df_eco[, mets[l]] ~ df_eco[, panspecies[j]])
    e_coef_val <- coef(mod_eco)[2]
    e_ci <- confint(mod_eco)[2, ]
    
    ########## Direct effect (FIXED)
    idx <- which(dsecr$X == paste0(species[j], "_", mets[l]))
    
    if (length(idx) > 0) {
      
      selected_values <- as.numeric(
        unlist(dsecr[idx, names_Mj, drop = FALSE])
      )
      
      selected_values <- selected_values[!is.na(selected_values)]
      n <- length(selected_values)
      
      if (n == 1) {
        mean_val <- selected_values
        ci_lower <- mean_val
        ci_upper <- mean_val
      } else if (n > 1) {
        mean_val <- mean(selected_values)
        sd_val <- sd(selected_values)
        err <- qt(0.975, df = n - 1) * sd_val / sqrt(n)
        ci_lower <- mean_val - err
        ci_upper <- mean_val + err
      } else {
        mean_val <- NA
        ci_lower <- NA
        ci_upper <- NA
      }
      
    } else {
      mean_val <- NA
      ci_lower <- NA
      ci_upper <- NA
    }
    
    ########## Total effect
    mod_total <- lm(df[, mets[l]] ~ df[, panspecies[j]])
    t_coef_val <- coef(mod_total)[2]
    t_ci <- confint(mod_total)[2, ]
    
    ########## Store
    base <- (j - 1) * 9
    effects[l, base + 1] <- e_coef_val
    effects[l, base + 2] <- e_ci[1]
    effects[l, base + 3] <- e_ci[2]
    effects[l, base + 4] <- mean_val
    effects[l, base + 5] <- ci_lower
    effects[l, base + 6] <- ci_upper
    effects[l, base + 7] <- t_coef_val
    effects[l, base + 8] <- t_ci[1]
    effects[l, base + 9] <- t_ci[2]
  }
  
  print(l)
}

effects$VMHID <- mets
effects <- effects[, c("VMHID", setdiff(names(effects), "VMHID"))]

write.xlsx(
  effects,
  'C:/Users/faesslerd/ecological_direct_total_effects_M4_T0_1837.xlsx',
  rowNames = FALSE
)
