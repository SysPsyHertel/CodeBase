library(ggplot2)
library(data.table)
library(gridExtra)
library(reshape2)
library(tidyr)
library(dplyr)
library(cowplot)
library(vegan)
library(sandwich)
library(emmeans)

## Individual maximum secretion fluxes
dsecr <- read.csv('MicrobeContributions_Fluxes.csv')

## Read generated redundancy measures
d <- read.csv('RedundancyMeasures_IBD.csv')

## read file of the cobratoolbox to map VMH to KEGGIds
VMH_metabolites <- fread('VMH_Metabolites.tsv')

VMH_metabolites$keggId[VMH_metabolites$abbreviation == "h2s"] <- "C00283"
VMH_metabolites$keggId[VMH_metabolites$abbreviation == "mnl"] <- "C00392"
VMH_metabolites$keggId[VMH_metabolites$abbreviation == "rbflvrd"] <- "C01007"

mets <- names(d)[grep("^KL_Sample", names(d))]
mets <- gsub(".*KL_Sample_", "", mets)

measures <- c("KL_Sample", "KL_Reference", "KL_Abundance")

d$Stratification <- factor(d$Stratification, levels = c("Healthy", "IBD_nondysbiotic", "IBD_dysbiotic"))

##############################################

## ANOVA
# Perform ANOVA for each diversity metric
modelrichness <- aov(d$spec ~ factor(d$Stratification))
# Extract the sum of squares
SS_model <- sum(modelrichness$coefficients^2)  # Model sum of squares
SS_residuals <- sum(residuals(modelrichness)^2)  # Residuals sum of squares
SS_total <- SS_model + SS_residuals

# Calculate R-squared
R_squared <- SS_model / SS_total

modelShannon <- aov(d$shannon ~ factor(d$Stratification))
SS_model <- sum(modelShannon$coefficients^2)  # Model sum of squares
SS_residuals <- sum(residuals(modelShannon)^2)  # Residuals sum of squares
SS_total <- SS_model + SS_residuals

# Calculate R-squared
R_squared <- SS_model / SS_total

modelSimpson <- aov(d$simpson ~ factor(d$Stratification))
modelEvenness <- aov(d$evenness ~ factor(d$Stratification))
modelrichness <- aov(d$spec ~ factor(d$Stratification))

# Extract the values for species richness
richness_df <- summary(modelrichness)[[1]]$Df[1]
richness_F <- summary(modelrichness)[[1]]$`F value`[1]
summary(modelrichness)[[1]]$`Pr(>F)`[1]


# Extract the values for Shannon
shannon_df <- summary(modelShannon)[[1]]$Df[1]
shannon_F <- summary(modelShannon)[[1]]$`F value`[1]
shannon_p <- summary(modelShannon)[[1]]$`Pr(>F)`[1]
summary(modelShannon)[[1]]$`Pr(>F)`[1]

# Extract the values for Simpson
simpson_df <- summary(modelSimpson)[[1]]$Df[1]
simpson_F <- summary(modelSimpson)[[1]]$`F value`[1]
simpson_p <- summary(modelSimpson)[[1]]$`Pr(>F)`[1]

# Extract the values for Evenness
evenness_df <- summary(modelEvenness)[[1]]$Df[1]
evenness_F <- summary(modelEvenness)[[1]]$`F value`[1]
evenness_p <- summary(modelEvenness)[[1]]$`Pr(>F)`[1]

# Create a table with all the information
result_table <- data.frame(
  Metric = c("Shannon", "Simpson", "Evenness"),
  Degrees_of_Freedom = c(shannon_df, simpson_df, evenness_df),
  F_value = c(shannon_F, simpson_F, evenness_F),
  p_value = c(shannon_p, simpson_p, evenness_p)
)

# Print the table
print(result_table)

anova_table <- data.frame(   
  Metabolite = character(),
  VMHId = character(),
  KeggId = character(),
  NoReference = double(),
  Producing.Species = double(),
  totalSum = double(),
  Interdependency = double(),
  Anova.Fval.Sample = double(),
  Anova.pval.Sample = double(),
  Anova.pvalbonf.Sample = double(),
  Anova.Fval.Reference = double(),
  Anova.pval.Reference = double(),
  Anova.pvalbonf.Reference = double(),
  Anova.Fval.Abundance = double(),
  Anova.pval.Abundance = double(),
  Anova.pvalbonf.Abundance = double(),
  stringsAsFactors = FALSE
)

for (i in mets) {
  
  Metabolite <- VMH_metabolites$fullName[match(i, VMH_metabolites$abbreviation)]
  VMHId <- i
  KeggId <- VMH_metabolites$keggId[match(i, VMH_metabolites$abbreviation)]
  metexch <- dsecr$X[grepl(paste0("_", i, "$"), dsecr$X)]
  NoReference <- length(metexch)
  Producing.Species <- mean(d[, paste0("N_", i)], na.rm = TRUE)
  snB <- mean(d[, paste0("N_", i)])
  Interdependency <- median(d[, paste0("KLI_", i)], na.rm = TRUE)
  
  # ANOVA models
  modelSamp <- aov(d[, paste0("KL_Sample_", i)] ~ factor(d$Stratification))
  modelRef <- aov(d[, paste0("KL_Reference_", i)] ~ factor(d$Stratification))
  modelAbu <- aov(d[, paste0("KL_Abundance_", i)] ~ factor(d$Stratification))
  
  # Extracting F-values and p-values
  Anova.Fval.Sample <- summary(modelSamp)[[1]]$`F value`[1]
  Anova.pval.Sample <- summary(modelSamp)[[1]]$`Pr(>F)`[1]
  
  Anova.Fval.Reference <- summary(modelRef)[[1]]$`F value`[1]
  Anova.pval.Reference <- summary(modelRef)[[1]]$`Pr(>F)`[1]
  
  Anova.Fval.Abundance <- summary(modelAbu)[[1]]$`F value`[1]
  Anova.pval.Abundance <- summary(modelAbu)[[1]]$`Pr(>F)`[1]
  
  # Bonferroni placeholders (to be updated after all calculations)
  Anova.pvalbonf.Sample <- NA
  Anova.pvalbonf.Reference <- NA
  Anova.pvalbonf.Abundance <- NA
  
  # Add the row to the anova_table
  anova_table <- rbind(anova_table, data.frame(
    Metabolite = Metabolite,
    VMHId = VMHId,
    KeggId = KeggId,
    NoReference = NoReference,
    Producing.Species = Producing.Species,
    totalSum = snB,
    Interdependency = Interdependency,
    Anova.Fval.Sample = Anova.Fval.Sample,
    Anova.pval.Sample = Anova.pval.Sample,
    Anova.pvalbonf.Sample = Anova.pvalbonf.Sample,
    Anova.Fval.Reference = Anova.Fval.Reference,
    Anova.pval.Reference = Anova.pval.Reference,
    Anova.pvalbonf.Reference = Anova.pvalbonf.Reference,
    Anova.Fval.Abundance = Anova.Fval.Abundance,
    Anova.pval.Abundance = Anova.pval.Abundance,
    Anova.pvalbonf.Abundance = Anova.pvalbonf.Abundance
  ))
}

pvalsbonf <- p.adjust(c(anova_table$Anova.pval.Sample, anova_table$Anova.pval.Reference, anova_table$Anova.pval.Abundance), method="bonferroni")

anova_table$Anova.pvalbonf.Sample <- pvalsbonf[1:134]
anova_table$Anova.pvalbonf.Reference <- pvalsbonf[135:268]
anova_table$Anova.pvalbonf.Abundance <- pvalsbonf[269:402]

sig_sample <- anova_table$VMHId[which(anova_table$Anova.pvalbonf.Sample < 0.05)]
sig_reference <- anova_table$VMHId[which(anova_table$Anova.pvalbonf.Reference < 0.05)]
sig_abundance <- anova_table$VMHId[which(anova_table$Anova.pvalbonf.Abundance < 0.05)]

sig_all <- unique(c(sig_sample, sig_reference, sig_abundance))
nonsig <- setdiff(anova_table$VMHId, sig_all)
nonsig

#write.xlsx(anova_table, "Tables/S9_anova.xlsx", rowNames = FALSE)

# Perform ANOVA
library(emmeans)
anova_table_contrasts <- data.frame(
  Metabolite = character(),
  VMHId = character(),
  KeggId = character(),
  Estimate.Sample.CI = character(),
  pvalue.KL_Sample = double(),
  Bonferroni.pvalue.KL_Sample = double(),
  Estimate.Reference.CI = character(),
  pvalue.KL_Reference = double(),
  Bonferroni.pvalue.KL_Reference = double(),
  Estimate.Abundance.CI = character(),
  pvalue.KL_Abundance = double(),
  Bonferroni.pvalue.KL_Abundance = double(),
  stringsAsFactors = FALSE
)

# Collect all p-values before applying Bonferroni correction
pvals_sample <- c()
pvals_reference <- c()
pvals_abundance <- c()

for (i in mets) {
  
  Metabolite <- VMH_metabolites$fullName[match(i, VMH_metabolites$abbreviation)]
  VMHId <- i
  KeggId <- VMH_metabolites$keggId[match(i, VMH_metabolites$abbreviation)]
  
  # ANOVA models
  modelSamp <- aov(d[, paste0("KL_Sample_", i)] ~ factor(d$Stratification))
  modelRef <- aov(d[, paste0("KL_Reference_", i)] ~ factor(d$Stratification))
  modelAbu <- aov(d[, paste0("KL_Abundance_", i)] ~ factor(d$Stratification))
  
  # Compute estimated means
  emmeansSamp <- emmeans(modelSamp, ~ Stratification)
  emmeansRef <- emmeans(modelRef, ~ Stratification)
  emmeansAbu <- emmeans(modelAbu, ~ Stratification)
  
  # Define and compute contrasts
  contrast_testSamp <- contrast(emmeansSamp, list("IBD vs Healthy" = c(-1, 0.5, 0.5)))
  contrast_testRef <- contrast(emmeansRef, list("IBD vs Healthy" = c(-1, 0.5, 0.5)))
  contrast_testAbu <- contrast(emmeansAbu, list("IBD vs Healthy" = c(-1, 0.5, 0.5)))
  
  # Extracting estimates, confidence intervals, and p-values
  estimate_sample <- summary(contrast_testSamp)$estimate
  lower_sample <- confint(contrast_testSamp)$lower.CL
  upper_sample <- confint(contrast_testSamp)$upper.CL
  pvalue_sample <- summary(contrast_testSamp)$p.value
  
  estimate_reference <- summary(contrast_testRef)$estimate
  lower_reference <- confint(contrast_testRef)$lower.CL
  upper_reference <- confint(contrast_testRef)$upper.CL
  pvalue_reference <- summary(contrast_testRef)$p.value
  
  estimate_abundance <- summary(contrast_testAbu)$estimate
  lower_abundance <- confint(contrast_testAbu)$lower.CL
  upper_abundance <- confint(contrast_testAbu)$upper.CL
  pvalue_abundance <- summary(contrast_testAbu)$p.value
  
  # Format as "estimate (CI_lower, CI_upper)"
  Estimate.Sample.CI <- paste0(round(estimate_sample, 5), " (",
                               formatC(lower_sample, format = "f", digits = 5), ", ", 
                               formatC(upper_sample, format = "f", digits = 5), ")")
  
  Estimate.Reference.CI <- paste0(round(estimate_reference, 5), " (",
                                  formatC(lower_reference, format = "f", digits = 5), ", ", 
                                  formatC(upper_reference, format = "f", digits = 5), ")")
  
  Estimate.Abundance.CI <- paste0(round(estimate_abundance, 5), " (",
                                  formatC(lower_abundance, format = "f", digits = 5), ", ", 
                                  formatC(upper_abundance, format = "f", digits = 5), ")")
  
  
  # Add the row without Bonferroni correction
  anova_table_contrasts <- rbind(anova_table_contrasts, data.frame(
    Metabolite = Metabolite,
    VMHId = VMHId,
    KeggId = KeggId,
    Estimate.Sample.CI = Estimate.Sample.CI,
    pvalue.KL_Sample = pvalue_sample,
    Estimate.Reference.CI = Estimate.Reference.CI,
    pvalue.KL_Reference = pvalue_reference,
    Estimate.Abundance.CI = Estimate.Abundance.CI,
    pvalue.KL_Abundance = pvalue_abundance,
    stringsAsFactors = FALSE
  ))
  
  # Collect all p-values for correction
  pvals_sample <- c(pvals_sample, pvalue_sample)
  pvals_reference <- c(pvals_reference, pvalue_reference)
  pvals_abundance <- c(pvals_abundance, pvalue_abundance)
  
  print(i)
}


#write.xlsx(anova_table_contrasts, "Tables/S10_anova_contrasts.xlsx", rowNames = FALSE)

contrasts_sig_sample <- anova_table_contrasts[anova_table_contrasts$VMHId %in% sig_sample, ]
contrasts_sig_reference <- anova_table_contrasts[anova_table_contrasts$VMHId %in% sig_reference, ]
contrasts_sig_abundance <- anova_table_contrasts[anova_table_contrasts$VMHId %in% sig_abundance, ]

# ANOVA models
model_species_richness <- aov(d[, "spec"] ~ factor(d$Stratification))
model_N_h2s <- aov(d[, "N_h2s"] ~ factor(d$Stratification))
model_N_glcn<- aov(d[, "N_glcn"] ~ factor(d$Stratification))

# Compute estimated means
emmeans_richness <- emmeans(model_species_richness, ~ Stratification)
emmeans_N_h2s <- emmeans(model_N_h2s, ~ Stratification)
emmeans_N_glcn <- emmeans(model_N_glcn, ~ Stratification)

# Define and compute contrasts
contrast_richness <- contrast(emmeans_richness, list("IBD vs Healthy" = c(-1, 0.5, 0.5)))
summary(contrast_richness)$estimate
confint(contrast_richness)$lower.CL
confint(contrast_richness)$upper.CL
summary(contrast_richness)$p.value

contrast_Nh2s <- contrast(emmeans_N_h2s, list("IBD vs Healthy" = c(-1, 0.5, 0.5)))
summary(contrast_Nh2s)$estimate
confint(contrast_Nh2s)$lower.CL
confint(contrast_Nh2s)$upper.CL
summary(contrast_Nh2s)$p.value

contrast_Nglcn <- contrast(emmeans_N_glcn, list("IBD vs Healthy" = c(-1, 0.5, 0.5)))
summary(contrast_Nglcn)$estimate
confint(contrast_Nglcn)$lower.CL
confint(contrast_Nglcn)$upper.CL
summary(contrast_Nglcn)$p.value

# Initialize the results table with the new structure
regr_table <- data.frame(   
  Metabolite = character(),
  VMHId = character(),
  KeggId = character(),
  
  # Sample estimates, p-values, and tests
  Estimate.Sample.nb.CI = character(),  
  pvalue.Sample.nb = double(),          
  Estimate.Sample.db.CI = character(),  
  pvalue.Sample.db = double(),          
  pval.Sample.WaldTest = double(),  # Moved here
  pval.Sample.Bonferroni = double(),  # Moved here
  
  # Reference estimates, p-values, and tests
  Estimate.Reference.nb.CI = character(),  
  pvalue.Reference.nb = double(),          
  Estimate.Reference.db.CI = character(),  
  pvalue.Reference.db = double(),          
  pval.Reference.WaldTest = double(),  # Moved here
  pval.Reference.Bonferroni = double(),  # Moved here
  
  # Abundance estimates, p-values, and tests
  Estimate.Abundance.nb.CI = character(),  
  pvalue.Abundance.nb = double(),          
  Estimate.Abundance.db.CI = character(),  
  pvalue.Abundance.db = double(),          
  pval.Abundance.WaldTest = double(),  # Moved here
  pval.Abundance.Bonferroni = double(),  # Moved here
  
  stringsAsFactors = FALSE
)

# Iterate over metabolites
for (i in mets) {
  
  # Retrieve metabolite info
  Metabolite <- VMH_metabolites$fullName[match(i, VMH_metabolites$abbreviation)]
  VMHId <- i
  KeggId <- VMH_metabolites$keggId[match(i, VMH_metabolites$abbreviation)]
  
  # Fit models for sample, reference, and abundance
  modelSamp <- lm(d[, paste0("KL_Sample_", i)] ~ d[, "shannon"] + factor(d$Stratification))
  modelRef  <- lm(d[, paste0("KL_Reference_", i)] ~ d[, "shannon"] + factor(d$Stratification))
  modelAbu  <- lm(d[, paste0("KL_Abundance_", i)] ~ d[, "shannon"] + factor(d$Stratification))
  
  # Function to extract estimates, CI, and p-values
  extract_CI_pvalue <- function(model, var) {
    estimate <- summary(model)$coef[var, 1]
    lower_CI <- confint(model)[var, 1]
    upper_CI <- confint(model)[var, 2]
    p_value  <- summary(model)$coef[var, 4]
    CI_string <- paste0(round(estimate, 5), " (", 
                        formatC(lower_CI, format = "f", digits = 5), ", ", 
                        formatC(upper_CI, format = "f", digits = 5), ")")
    return(list(CI = CI_string, p_value = p_value))
  }
  
  # Extract values for each condition
  sample_nb <- extract_CI_pvalue(modelSamp, "factor(d$Stratification)IBD_nondysbiotic")
  sample_db <- extract_CI_pvalue(modelSamp, "factor(d$Stratification)IBD_dysbiotic")
  ref_nb    <- extract_CI_pvalue(modelRef, "factor(d$Stratification)IBD_nondysbiotic")
  ref_db    <- extract_CI_pvalue(modelRef, "factor(d$Stratification)IBD_dysbiotic")
  abu_nb    <- extract_CI_pvalue(modelAbu, "factor(d$Stratification)IBD_nondysbiotic")
  abu_db    <- extract_CI_pvalue(modelAbu, "factor(d$Stratification)IBD_dysbiotic")
  
  # Perform Wald test
  vartest <- c("factor(d$Stratification)IBD_dysbiotic", "factor(d$Stratification)IBD_nondysbiotic")
  waldtestS <- car::linearHypothesis(modelSamp, vartest, vcov = vcovHC(modelSamp, type = "HC"))
  waldtestR <- car::linearHypothesis(modelRef, vartest, vcov = vcovHC(modelRef, type = "HC"))
  waldtestA <- car::linearHypothesis(modelAbu, vartest, vcov = vcovHC(modelAbu, type = "HC"))
  
  pval.Sample.WaldTest    <- waldtestS$`Pr(>F)`[2]
  pval.Reference.WaldTest <- waldtestR$`Pr(>F)`[2]
  pval.Abundance.WaldTest <- waldtestA$`Pr(>F)`[2]
  
  # Bonferroni correction
  n_tests <- length(mets)
  pval.Sample.Bonferroni    <- NA
  pval.Reference.Bonferroni <- NA
  pval.Abundance.Bonferroni <- NA
  
  # Append results to the table
  regr_table <- rbind(regr_table, data.frame(
    Metabolite = Metabolite,
    VMHId = VMHId,
    KeggId = KeggId,
    
    # Sample-based results
    Estimate.Sample.nb.CI = sample_nb$CI,
    pvalue.Sample.nb = sample_nb$p_value,
    Estimate.Sample.db.CI = sample_db$CI,
    pvalue.Sample.db = sample_db$p_value,
    pval.Sample.WaldTest = pval.Sample.WaldTest,
    pval.Sample.Bonferroni = pval.Sample.Bonferroni,
    
    # Reference-based results
    Estimate.Reference.nb.CI = ref_nb$CI,
    pvalue.Reference.nb = ref_nb$p_value,
    Estimate.Reference.db.CI = ref_db$CI,
    pvalue.Reference.db = ref_db$p_value,
    pval.Reference.WaldTest = pval.Reference.WaldTest,
    pval.Reference.Bonferroni = pval.Reference.Bonferroni,
    
    # Abundance-based results
    Estimate.Abundance.nb.CI = abu_nb$CI,
    pvalue.Abundance.nb = abu_nb$p_value,
    Estimate.Abundance.db.CI = abu_db$CI,
    pvalue.Abundance.db = abu_db$p_value,
    pval.Abundance.WaldTest = pval.Abundance.WaldTest,
    pval.Abundance.Bonferroni = pval.Abundance.Bonferroni
  ))
}


pvalsbonf <- p.adjust(c(regr_table$pval.Sample.WaldTest, regr_table$pval.Reference.WaldTest, regr_table$pval.Abundance.WaldTest), method="bonferroni")

regr_table$pval.Sample.Bonferroni <- pvalsbonf[1:134]
regr_table$pval.Reference.Bonferroni <- pvalsbonf[135:268]
regr_table$pval.Abundance.Bonferroni <- pvalsbonf[269:402]



sig_sample <- regr_table$VMHId[which(regr_table$pval.Sample.Bonferroni < 0.05)]
sig_reference <- regr_table$VMHId[which(regr_table$pval.Reference.Bonferroni < 0.05)]
sig_abundance <- regr_table$VMHId[which(regr_table$pval.Abundance.Bonferroni  < 0.05)]

sig_all <- unique(c(sig_sample, sig_reference, sig_abundance))
nonsig <- setdiff(anova_table$VMHId, sig_all)