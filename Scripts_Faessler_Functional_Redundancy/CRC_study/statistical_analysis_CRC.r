library(ggplot2)
library(data.table)
library(gridExtra)
library(reshape2)
library(tidyr)
library(dplyr)
library(cowplot)
library(openxlsx)

setwd(r'(F:\Functional_Redundancy_CSBJ\R\ReSubmission\CRC)')
dsecr <- read.csv('Microbe_Secretion.csv')
d <- read.csv('RedundancyMeasures.csv')
VMH_metabolites <- fread('VMH_Metabolites.tsv')

mets <- names(d)[grep("^KL_Sample", names(d))]
mets <- gsub(".*KL_Sample_", "", mets)

d$healthstat <- as.factor(d$healthstat)
measures <- c("KL_Sample", "KL_Reference", "KL_Abundance")

##############################################
#### Interrelations redundancy measures   ####
##############################################

dc <- d
Rsqtable <- data.frame(VMHId=character(),
                       Producing.Species=double(),
                       Sample_Reference=double(),
                       Sample_Abundance=double(),
                       Reference_Abundance=double(),
                       stringsAsFactors=FALSE)

plot_list <- vector("list", length(mets))
cnt <- 0

coefs <- c()
for(i in mets){
  cnt <- cnt+1
  frm12 <- as.formula(paste(paste0("KL_Sample_",i), paste(paste0("KL_Reference_",i), collapse=" + "), sep=" ~ "))
  frm13 <- as.formula(paste(paste0("KL_Sample_",i), paste(paste0("KL_Abundance_",i), collapse=" + "), sep=" ~ "))
  frm23 <- as.formula(paste(paste0("KL_Reference_",i), paste(paste0("KL_Abundance_",i), collapse=" + "), sep=" ~ "))
  
  fit1 <- lm(frm12, data = d)
  coefs <- c(coefs, as.numeric(fit1$coefficients[1]))
  fit2 <- lm(frm13, data = d)
  fit3 <- lm(frm23, data = d)

  Rsqtable[nrow(Rsqtable)+1, 1] <- i
  Rsqtable[nrow(Rsqtable), 2:ncol(Rsqtable)] <- c(mean(d[, paste0("N_",i)]), summary(fit1)$r.squared, summary(fit2)$r.squared, summary(fit3)$r.squared)
  print(i)
}

RSQT <- Rsqtable
RSQT$Metabolite <- VMH_metabolites$fullName[match(RSQT$VMHId, VMH_metabolites$abbreviation)]
RSQT$KeggId <- VMH_metabolites$keggId[match(RSQT$VMHId, VMH_metabolites$abbreviation)]
  

Mapping <- RSQT[, c("Metabolite", "VMHId", "KeggId")]
ref <- c()

for(i in mets){
  metexch <- dsecr$X[grepl(paste0("_", i, "$"), dsecr$X)]
  print(length(metexch))
  ref <- c(ref, length(metexch))
}

RSQT$NoReference <- ref
colnames(RSQT)
new_order <- c("Metabolite", "VMHId", "KeggId", "NoReference", "Producing.Species", "Sample_Reference", "Sample_Abundance", "Reference_Abundance") # Define the new order
RSQT <- RSQT[, new_order]

#write.xlsx(RSQT, "Tables/S4.xlsx", rowNames = FALSE)

combined_data <- bind_rows(
  transform(d, Category = "Glycine", X = KL_Sample_gly, Y = KL_Reference_gly),
  transform(d, Category = "Tryptamine", X = KL_Sample_trypta, Y = KL_Reference_trypta),
  transform(d, Category = "Glucosamine", X = KL_Sample_gam, Y = KL_Reference_gam),
)


##########################################################
#### ### Species diversity <-> functional redundancy  ####
##########################################################
div_red <- data.frame(Metabolite=character(),
                      VMHId=character(),
                      KeggId=character(),
                      No.Reference=double(),
                      Producing.Species=double(),
                      Total.Sum.Prod=double(),
                      median.Interdependencies=double(),
                      Estimate.Sample.CI=double(),
                      pvalue.KL_Sample=double(),
                      Rsq.Sample=double(),
                      Estimate.Reference.CI=double(),
                      pvalue.KL_Reference=double(),
                      Rsq.Reference=double(),
                      Estimate.Abundance.CI=double(),
                      pvalue.KL_Abundance=double(),
                      Rsq.Abundance=double(),
                      stringsAsFactors=FALSE)

for(i in mets){
  modelSamp <- lm(d[, paste0("KL_Sample_", i)] ~ d[, "shannon"])
  rsquared_samp <- summary(modelSamp)$r.squared
  
  modelRef <- lm(d[, paste0("KL_Reference_", i)] ~ d[, "shannon"])
  rsquared_ref <- summary(modelRef)$r.squared
  
  modelAbu <- lm(d[, paste0("KL_Abundance_", i)] ~ d[, "shannon"])
  rsquared_abu <- summary(modelAbu)$r.squared
  
  # Round coefficients to 3 decimal places and format confidence intervals to 5 decimal places
  modelSampv <- paste0(round(summary(modelSamp)$coef[2, 1], 5), " (",
                       formatC(confint(modelSamp)[2,1], format = "f", digits = 5), ", ", 
                       formatC(confint(modelSamp)[2,2], format = "f", digits = 5), ")")
  
  modelRefv <- paste0(round(summary(modelRef)$coef[2, 1], 5), " (",
                      formatC(confint(modelRef)[2,1], format = "f", digits = 5), ", ", 
                      formatC(confint(modelRef)[2,2], format = "f", digits = 5), ")")
  
  modelAbuv <- paste0(round(summary(modelAbu)$coef[2, 1], 5), " (",
                      formatC(confint(modelAbu)[2,1], format = "f", digits = 5), ", ", 
                      formatC(confint(modelAbu)[2,2], format = "f", digits = 5), ")")
  
  NoReference <- RSQT$NoReference[match(i, RSQT$VMHId)]
  ProSpecies <- RSQT$Producing.Species[match(i, RSQT$VMHId)]
  snB <- mean(d[, paste0("sNB_", i)])
  med <- median(d[, paste0("KLI_", i)], na.rm=TRUE)
  
  div_red[nrow(div_red)+1, ] <- c(Mapping$Metabolite[match(i, Mapping$VMHId)], i, Mapping$KeggId[match(i, Mapping$VMHId)], NoReference, ProSpecies, snB, med,
                                  modelSampv, summary(modelSamp)$coef[2, 4], rsquared_samp,
                                  modelRefv, summary(modelRef)$coef[2, 4], rsquared_ref,
                                  modelAbuv, summary(modelAbu)$coef[2, 4], rsquared_abu)
}


pvaluesdiv_red <- c(div_red$pvalue.KL_Sample, div_red$pvalue.KL_Reference, div_red$pvalue.KL_Abundance)
pvaluesdiv_red <- as.numeric(pvaluesdiv_red)
Bonfpvaldiv_red <- p.adjust(pvaluesdiv_red, method="bonferroni")

nonsigSample <- div_red$Metabolite[which(div_red$BonfKL_Sample>=0.05)]
nonsigRef <- div_red$Metabolite[which(div_red$BonfKL_Reference>=0.05)]
div_red$Metabolite[div_red$BonfKL_Sample<0.05 & sign(as.numeric(div_red$Estimate.KL_Sample))==1]

div_red$Metabolite[div_red$BonfKL_Abundance<0.05 & sign(as.numeric(div_red$Estimate.KL_Abundance))==1]
div_red$Metabolite[div_red$BonfKL_Abundance<0.05]

div_red$BonfKL_Sample <- Bonfpvaldiv_red[1:45]
div_red$BonfKL_Reference <- Bonfpvaldiv_red[46:90]
div_red$BonfKL_Abundance <- Bonfpvaldiv_red[91:135]

columns_to_numeric <- c(4:7, 10, 12, 13, 15, 16)
div_red[, columns_to_numeric] <- lapply(div_red[, columns_to_numeric], as.numeric)

div_red <- div_red[, c(1:9, 17, 10:16, 18, 19:ncol(div_red))]

# Step 2: Move column 18 after column 14
div_red <- div_red[, c(1:14, 18, 15:16, 17, 19:ncol(div_red))]


#########################################################
#### faecal concentrations <-> functional redundancy  ####
#########################################################
set.seed(123)
# Initialize the list to store valid metabolites
metsana <- c()

# Check each metabolite and filter based on the log transformation condition
for (i in mets) {
  j <- ifelse(grepl("^[[:digit:]]+", i), paste0("X", i), i)
  if (!j %in% names(d)) next  # Skip if 'j' is not a column in 'd'
  
  logvi <- log(d[[j]])
  logvi[is.infinite(logvi)] <- NA
  if (sum(is.na(logvi)) < 0.5 * length(logvi)) {
    metsana <- c(metsana, j)
  }
}

# Remove the prefixed 'X' from valid metabolites
metsana <- gsub("^X", "", metsana)

# Number of permutations
n_perm <- 10000

# Function to compute p-value for permutation test
permute_test <- function(model_formula, data, var_name, n_perm) {
  # Fit the original model
  model <- lm(model_formula, data = data)
  
  # Get the observed p-value for the original variable
  observed_p_value <- coef(summary(model))[2, 4]  # 2,4 is the p-value of the second coefficient
  
  # Initialize vector to store permuted p-values
  perm_p_values <- numeric(n_perm)
  
  # Perform permutation test
  for (i in 1:n_perm) {
    # Permute the variable of interest (e.g., KL_Sample, KL_Reference, etc.)
    data$perm_var <- sample(data[[var_name]])  # Create permuted variable
    
    # Update the model formula to use the permuted variable
    perm_formula <- update(model_formula, paste(". ~ . -", var_name, "+ perm_var"))
    
    # Fit the model with permuted data
    perm_model <- lm(perm_formula, data = data)
    
    # Extract the p-value for the permuted variable (perm_var is usually the last row)
    perm_p_values[i] <- coef(summary(perm_model))["perm_var", "Pr(>|t|)"]
  }
  
  # Calculate p-value based on the permutation distribution
  # Count the number of permuted p-values that are as extreme as or more extreme than the observed p-value
  p_value <- mean(perm_p_values <= observed_p_value)  # For one-tailed test (lower)
  return(p_value)
}

# Initialize the parameter table
paramtab_faec <- data.frame(
  Metabolite = character(),
  VMHId = character(),
  KeggId = character(),
  Estimate.Sample.CI = character(),
  pvalue.KL_Sample = double(),
  pval.perm.sample = double(),  # Store permutation p-value for Sample
  Estimate.Reference.CI = character(),
  pvalue.KL_Reference = double(),
  pval.perm.reference = double(),  # Store permutation p-value for Reference
  Estimate.Abundance.CI = character(),
  pvalue.KL_Abundance = double(),
  pval.perm.abundance = double(),  # Store permutation p-value for Abundance
  stringsAsFactors = FALSE
)

# Populate the parameter table
for (i in metsana) {
  j <- ifelse(grepl("^[[:digit:]]+", i), paste0("X", i), i)
  
  logvi <- log(d[[j]])
  logvi[is.infinite(logvi)] <- NA
  
  # Fit the models
  modelKL_Sample <- lm(logvi ~ d[[paste0("KL_Sample_", i)]] + d$Age + d$BMI + d$Gender + d$healthstat)
  modelKL_Reference <- lm(logvi ~ d[[paste0("KL_Reference_", i)]] + d$Age + d$BMI + d$Gender + d$healthstat)
  modelKL_Abundance <- lm(logvi ~ d[[paste0("KL_Abundance_", i)]] + d$Age + d$BMI + d$Gender + d$healthstat)
  
  # Run permutation tests
  pval.perm.sample <- permute_test(logvi ~ d[[paste0("KL_Sample_", i)]] + d$Age + d$BMI + d$Gender + d$healthstat, d, paste0("KL_Sample_", i), n_perm)
  pval.perm.reference <- permute_test(logvi ~ d[[paste0("KL_Reference_", i)]] + d$Age + d$BMI + d$Gender + d$healthstat, d, paste0("KL_Reference_", i), n_perm)
  pval.perm.abundance <- permute_test(logvi ~ d[[paste0("KL_Abundance_", i)]] + d$Age + d$BMI + d$Gender + d$healthstat, d, paste0("KL_Abundance_", i), n_perm)
  
  # Format coefficients with confidence intervals
  modelSampf <- paste0(round(coef(summary(modelKL_Sample))[2, 1], 5), " (",
                       formatC(confint(modelKL_Sample)[2, 1], format = "f", digits = 5), ", ", 
                       formatC(confint(modelKL_Sample)[2, 2], format = "f", digits = 5), ")")
  
  modelReff <- paste0(round(coef(summary(modelKL_Reference))[2, 1], 5), " (",
                      formatC(confint(modelKL_Reference)[2, 1], format = "f", digits = 5), ", ", 
                      formatC(confint(modelKL_Reference)[2, 2], format = "f", digits = 5), ")")
  
  modelAbuf <- paste0(round(coef(summary(modelKL_Abundance))[2, 1], 5), " (",
                      formatC(confint(modelKL_Abundance)[2, 1], format = "f", digits = 5), ", ", 
                      formatC(confint(modelKL_Abundance)[2, 2], format = "f", digits = 5), ")")
  
  # Add row to paramtab_faec with permutation p-values
  paramtab_faec <- rbind(paramtab_faec, data.frame(
    Metabolite = Mapping$Metabolite[match(i, Mapping$VMHId)],
    VMHId = i,
    KeggId = Mapping$KeggId[match(i, Mapping$VMHId)],  # Added KeggId
    Estimate.Sample.CI = modelSampf,
    pvalue.KL_Sample = summary(modelKL_Sample)$coef[2, 4],
    pval.perm.sample = pval.perm.sample,
    Estimate.Reference.CI = modelReff,
    pvalue.KL_Reference = summary(modelKL_Reference)$coef[2, 4],
    pval.perm.reference = pval.perm.reference,
    Estimate.Abundance.CI = modelAbuf,
    pvalue.KL_Abundance = summary(modelKL_Abundance)$coef[2, 4],
    pval.perm.abundance = pval.perm.abundance,
    stringsAsFactors = FALSE
  ))
  print(i)
}


# Adjust p-values for multiple testing
pvalues <- c(paramtab_faec$pvalue.KL_Sample, paramtab_faec$pvalue.KL_Reference, paramtab_faec$pvalue.KL_Abundance)
Bonfpval <- p.adjust(pvalues, method = "bonferroni")

paramtab_faec$BonfKL_Sample <- Bonfpval[1:nrow(paramtab_faec)]
paramtab_faec$BonfKL_Reference <- Bonfpval[(nrow(paramtab_faec) + 1):(2 * nrow(paramtab_faec))]
paramtab_faec$BonfKL_Abundance <- Bonfpval[(2 * nrow(paramtab_faec) + 1):(3 * nrow(paramtab_faec))]

# Filter significant metabolites
KL_Samplesig <- paramtab_faec$VMHId[paramtab_faec$BonfKL_Sample < 0.05]
KL_Referencesig <- paramtab_faec$VMHId[paramtab_faec$BonfKL_Reference < 0.05]
KL_Abundancesig <- paramtab_faec$VMHId[paramtab_faec$BonfKL_Abundance < 0.05]

pvalues_permutation <- c(paramtab_faec$pval.perm.sample, paramtab_faec$pval.perm.reference, paramtab_faec$pval.perm.abundance)
Bonfpval_permutation <- p.adjust(pvalues_permutation, method = "bonferroni")

paramtab_faec$BonfKL_permuation_Sample <- Bonfpval_permutation[1:nrow(paramtab_faec)]
paramtab_faec$BonfKL_permuation_Reference <- Bonfpval_permutation[(nrow(paramtab_faec) + 1):(2 * nrow(paramtab_faec))]
paramtab_faec$BonfKL_permuation_Abundance <- Bonfpval_permutation[(2 * nrow(paramtab_faec) + 1):(3 * nrow(paramtab_faec))]


# Filter significant metabolites (permuation)
KL_Samplesig <- paramtab_faec$VMHId[paramtab_faec$BonfKL_Sample < 0.05]
KL_Referencesig <- paramtab_faec$VMHId[paramtab_faec$BonfKL_Reference < 0.05]
KL_Abundancesig <- paramtab_faec$VMHId[paramtab_faec$BonfKL_Abundance < 0.05]


# Check the structure of 'dt'
## Plot
dt <- d

# Prepare columns for plotting
cols <- unique(c(KL_Samplesig, KL_Referencesig, KL_Abundancesig))
cols <- gsub("^(\\d+)", "X\\1", cols)  # Ensure variables starting with digits are prefixed with "X"

# Log-transform the relevant columns
dt[cols] <- lapply(dt[cols], function(x) {
  logx <- log(x)
  logx[is.infinite(logx)] <- NA
  return(logx)
})

paramtab_faec <- paramtab_faec %>%
  select(
    Metabolite,
    VMHId,
    KeggId,
    Estimate.Sample.CI,
    pvalue.KL_Sample,
    pval.perm.sample,
    BonfKL_Sample,
    BonfKL_permuation_Sample,  # Fixed typo here
    Estimate.Reference.CI,
    pvalue.KL_Reference,
    pval.perm.reference,
    BonfKL_Reference,
    BonfKL_permuation_Reference,  # Fixed typo here
    Estimate.Abundance.CI,
    pvalue.KL_Abundance,
    pval.perm.abundance,
    BonfKL_Abundance,
    BonfKL_permuation_Abundance  # Fixed typo here
  )


#################################################
#### Functional redundancy <-> Healthstatus  ####
#################################################
# Initialize the parameter table with KEGGId column
# Initialize the parameter table with the desired column order
paramtabHealthstat <- data.frame(
  Metabolite = character(),
  VMHId = character(),
  KeggId = character(),  # Added KeggId
  Estimate.Sample.CI = character(),
  pvalue.KL_Sample = double(),
  BonfKL_Sample = double(),  # Bonferroni correction for Sample
  Estimate.Reference.CI = character(),
  pvalue.KL_Reference = double(),
  BonfKL_Reference = double(),  # Bonferroni correction for Reference
  Estimate.Abundance.CI = character(),
  pvalue.KL_Abundance = double(),
  BonfKL_Abundance = double(),  # Bonferroni correction for Abundance
  Estimate.Interdependency.CI = character(),
  pvalue.Interdependency = double(),
  BonfKL_I = double(),  # Bonferroni correction for Interdependency
  stringsAsFactors = FALSE
)

# Initialize an empty vector to collect all p-values
all_pvalues <- numeric()

# Loop through each metabolite
for (i in mets) {
  # Fit the models
  modelSamp <- lm(d[, paste0("KL_Sample_", i)] ~ d[, "healthstat"] + d[, "Age"] + d[, "BMI"] + d[, "Gender"])
  modelRef <- lm(d[, paste0("KL_Reference_", i)] ~ d[, "healthstat"] + d[, "Age"] + d[, "BMI"] + d[, "Gender"])
  modelAbu <- lm(d[, paste0("KL_Abundance_", i)] ~ d[, "healthstat"] + d[, "Age"] + d[, "BMI"] + d[, "Gender"])
  logvi <- log(d[, paste0("KLI_", i)])
  logvi[which(logvi == -Inf)] <- NaN
  modelI <- lm(logvi ~ d[, "healthstat"] + d[, "Age"] + d[, "BMI"] + d[, "Gender"])
  
  # Format coefficients with confidence intervals
  modelSampv <- paste0(format(summary(modelSamp)$coef[2, 1], scientific = FALSE, digits = 3), " (",
                       format(confint(modelSamp)[2, 1], scientific = FALSE, digits = 3), ",", 
                       format(confint(modelSamp)[2, 2], scientific = FALSE, digits = 3), ")")
  
  modelRefv <- paste0(format(summary(modelRef)$coef[2, 1], scientific = FALSE, digits = 3), " (",
                      format(confint(modelRef)[2, 1], scientific = FALSE, digits = 3), ",", 
                      format(confint(modelRef)[2, 2], scientific = FALSE, digits = 3), ")")
  
  modelAbuv <- paste0(format(summary(modelAbu)$coef[2, 1], scientific = FALSE, digits = 3), " (",
                      format(confint(modelAbu)[2, 1], scientific = FALSE, digits = 3), ",", 
                      format(confint(modelAbu)[2, 2], scientific = FALSE, digits = 3), ")")
  
  modelIv <- paste0(format(summary(modelI)$coef[2, 1], scientific = FALSE, digits = 3), " (",
                    format(confint(modelI)[2, 1], scientific = FALSE, digits = 3), ",", 
                    format(confint(modelI)[2, 2], scientific = FALSE, digits = 3), ")")
  
  # Collect p-values for each model
  pvalues <- c(
    summary(modelSamp)$coef[2, 4], 
    summary(modelRef)$coef[2, 4], 
    summary(modelAbu)$coef[2, 4], 
    summary(modelI)$coef[2, 4]
  )
  
  # Add p-values to the list
  all_pvalues <- c(all_pvalues, pvalues)
  
  # Add row to the parameter table without Bonferroni correction (we will correct later)
  paramtabHealthstat[nrow(paramtabHealthstat) + 1, ] <- c(
    Mapping$Metabolite[match(i, Mapping$VMHId)], 
    i, 
    Mapping$KeggId[match(i, Mapping$VMHId)],  # Added KEGGId
    modelSampv, 
    summary(modelSamp)$coef[2, 4], 
    NA,  # Placeholder for Bonferroni corrected p-value for Sample
    modelRefv, 
    summary(modelRef)$coef[2, 4], 
    NA,  # Placeholder for Bonferroni corrected p-value for Reference
    modelAbuv, 
    summary(modelAbu)$coef[2, 4], 
    NA,  # Placeholder for Bonferroni corrected p-value for Abundance
    modelIv, 
    summary(modelI)$coef[2, 4], 
    NA   # Placeholder for Bonferroni corrected p-value for Interdependency
  )
}

# Apply Bonferroni correction to the collected p-values (for 4*45 = 180 tests)
Bonfpval <- p.adjust(all_pvalues, method = "bonferroni")

# Now, assign the corrected p-values to the paramtabHealthstat dataframe
start_index <- 1  # To start assigning corrected p-values
for (i in 1:nrow(paramtabHealthstat)) {
  # Corrected p-values for each model (Sample, Reference, Abundance, Interdependency)
  paramtabHealthstat$pvalue.KL_Sample[i] <- paramtabHealthstat$pvalue.KL_Sample[i]
  paramtabHealthstat$BonfKL_Sample[i] <- Bonfpval[start_index]
  start_index <- start_index + 1
  
  paramtabHealthstat$pvalue.KL_Reference[i] <- paramtabHealthstat$pvalue.KL_Reference[i]
  paramtabHealthstat$BonfKL_Reference[i] <- Bonfpval[start_index]
  start_index <- start_index + 1
  
  paramtabHealthstat$pvalue.KL_Abundance[i] <- paramtabHealthstat$pvalue.KL_Abundance[i]
  paramtabHealthstat$BonfKL_Abundance[i] <- Bonfpval[start_index]
  start_index <- start_index + 1
  
  paramtabHealthstat$pvalue.Interdependency[i] <- paramtabHealthstat$pvalue.Interdependency[i]
  paramtabHealthstat$BonfKL_I[i] <- Bonfpval[start_index]
  start_index <- start_index + 1
}


modelI_g_spec <- lm(d[, "spec"] ~ d[, "healthstat"] + d[, "Age"] + d[, "BMI"] + d[, "Gender"])
summary(modelI_g_spec)
confint(modelI_g_spec)


modelI_g_shannon <- lm(d[, "shannon"] ~ d[, "healthstat"] + d[, "Age"] + d[, "BMI"] + d[, "Gender"])
summary(modelI_g_shannon)
confint(modelI_g_shannon)

columns_to_numeric <- c(5,6,8,9,11,12,14,15)

paramtabHealthstat[, columns_to_numeric] <- lapply(paramtabHealthstat[, columns_to_numeric], as.numeric)


## Summary stats of pyruvate
modelRef <- lm(d[, paste0("KL_Reference_", "pyr")] ~ d[, "healthstat"]  + d[, "Age"] + d[, "BMI"] + d[, "Gender"])
coef(summary(modelRef))[2,4]
confint(modelRef)

## Global redundancy measure
mets_I <- gsub("^", "KLI_", mets)

d$median_mets_I <- log(apply(d[, mets_I], 1, median))


modelI_g <- lm(d[, "median_mets_I"] ~ d[, "healthstat"] + d[, "Age"] + d[, "BMI"] + d[, "Gender"])
summary(modelI_g)
confint(modelI_g)
