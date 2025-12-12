I am inlibrary(openxlsx)
library(gplots)
library(dplyr)

setwd('P:/AG_Hertel/Fluxomics/Microbiome_Modelling/MetaPhlAn4')

mph4_t0 <- read.csv("E:/BW/MPH4_t0.csv")


# read net secretions, direct contributions of the predictMicrobeContributions function and normalized abundances
dnetsecr <- read.csv('net_secretion_MPH4.csv')
dsecr <- read.csv('predictContributions_MPH4.csv')
dnormcov <- read.csv('normalizedCoverage_MPH4.csv')


# remove duplicates
mets <- gsub("^EX_|\\[fe\\]$", "", dnetsecr$Net.secretion)
dnetsecr <- dnetsecr[ , gsub("^", "S", mph4_t0$proband)]



interscolnames <- intersect(colnames(dsecr), colnames(dnetsecr))
interscolnames <- intersect(interscolnames, colnames(dnormcov))

dsecr <- dsecr[, c("X", interscolnames)]
dnormcov <- dnormcov[, c("ID", interscolnames)]
#dsecr[,2:ncol(dsecr)] <- (-1)*dsecr[,2:ncol(dsecr)]


all_zero_rows <- apply(dnormcov[, 2:ncol(dnormcov)], 1, function(x) all(x == 0))

# See which rows are all zeros
which(all_zero_rows)


# Print result
panspecies <- dnormcov$ID
species <- gsub("^pan", "", panspecies)
####
effects <- data.frame(matrix(ncol = length(species)*9, nrow = length(mets)))

# Create a vector for the column names
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

# Loop to calculate the average direct species net production effect for each metabolite (l) and microbe (j)
for(l in 1:length(mets)) {
  for(j in 1:length(panspecies)) {
    
    ### metabolite l
    rows <- grepl(paste0("_", mets[l], "$"), dsecr$X)
    numeric_cols <- sapply(dsecr, is.numeric)
    sums <- colSums(dsecr[rows, numeric_cols], na.rm = TRUE)
    
    df <- data.frame(
      ID = names(sums),
      j = sums
    )
    colnames(df)[2] <- mets[l]
    
    ### species j
    
    subspec <- dnormcov[dnormcov$ID == panspecies[j], -1]
    names_Mj <- colnames(subspec)[which(subspec > 0)]
    species_colnames <- names(subspec)
    subspec <- as.numeric(subspec[1, ])
    subspec_bin <- as.numeric(subspec > 0)
    
    subspec_long <- data.frame(
      ID = species_colnames,
      value = subspec_bin
    )
    colnames(subspec_long)[2] <- panspecies[j]
    
    df <- merge(df, subspec_long, by ="ID")
    
    ########## Ecological effect (e_species)
    rows_lj <- grepl(paste0("_", mets[l], "$"), dsecr$X) & 
      !grepl(paste0("^", species[j]), dsecr$X)
    
    sums_eco <- colSums(dsecr[rows_lj, numeric_cols], na.rm = TRUE)
    
    df_eco <- data.frame(
      ID = names(sums_eco),
      j = sums_eco
    )
    colnames(df_eco)[2] <- mets[l]
    df_eco <- merge(df_eco, subspec_long, by ="ID")
    
    mod_eco <- lm(df_eco[, mets[l]] ~ df_eco[, panspecies[j]])
    e_coef_val <- summary(mod_eco)$coef[2, 1]
    e_confint_vals <- confint(mod_eco)[2, ]  # 95% CI by default
    e_ci_lower <- e_confint_vals[1]
    e_ci_upper <- e_confint_vals[2]
    
    ########## Direct effect (d_species)
    if(any(which(dsecr$X == paste0(species[j], "_", mets[l])))) {
      abs_Mj <- sum(subspec > 0)
      
      # Extract the relevant data
      selected_values <- dsecr[which(dsecr$X == paste0(species[j], "_", mets[l])), names_Mj]
      flattened_values <- as.numeric(unlist(selected_values))  # flatten if it's a matrix or data frame
      
      n <- length(flattened_values)
      mean_val <- mean(flattened_values)
      sd_val <- sd(flattened_values)
      
      # 95% CI using t-distribution
      error_margin <- qt(0.975, df = n - 1) * sd_val / sqrt(n)
      ci_lower <- mean_val - error_margin
      ci_upper <- mean_val + error_margin
    } else {
      # If the condition is not met, set the direct effect values to NA
      mean_val <- NA
      ci_lower <- NA
      ci_upper <- NA
    }
    
    ########## Total effect (t_species)
    mod_total <- lm(df[, mets[l]] ~ df[, panspecies[j]])
    t_coef_val <- round(summary(mod_total)$coef[2, 1], 9)
    t_confint_vals <- confint(mod_total)[2, ]  # 95% CI by default
    t_ci_lower <- t_confint_vals[1]
    t_ci_upper <- t_confint_vals[2]
    
    ########## Save results in the effects dataframe
    effects[l, (j - 1) * 9 + 1] <- e_coef_val  # Ecological effect
    effects[l, (j - 1) * 9 + 2] <- e_ci_lower  # Ecological effect lower bound
    effects[l, (j - 1) * 9 + 3] <- e_ci_upper  # Ecological effect upper bound
    effects[l, (j - 1) * 9 + 4] <- mean_val  # Direct effect (mean)
    effects[l, (j - 1) * 9 + 5] <- ci_lower  # Direct effect lower bound
    effects[l, (j - 1) * 9 + 6] <- ci_upper  # Direct effect upper bound
    effects[l, (j - 1) * 9 + 7] <- t_coef_val  # Total effect
    effects[l, (j - 1) * 9 + 8] <- t_ci_lower  # Total effect lower bound
    effects[l, (j - 1) * 9 + 9] <- t_ci_upper  # Total effect upper bound
    
  }
  print(l)
}

effects$VMHID <- mets
# Assuming your dataframe is named df
effects <- effects[, c("VMHID", setdiff(names(effects), "VMHID"))]

write.xlsx(effects, 'E:/BW/ecological_direct_total_effects_M4_T0_1838.xlsx', rowNames=FALSE)