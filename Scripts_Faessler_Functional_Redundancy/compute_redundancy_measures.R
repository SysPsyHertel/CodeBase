library(ggplot2)
library(dplyr)
library(vegan)
library(philentropy)
library(haven)


setwd(r'(C:\Users\faesslerd\Documents\Projects\CRC_Redundancy\Final_Vers)')

dsecr <- read.csv('Microbe_Secretion.csv')
dmet <- read.csv('metadata.csv')
dnormcov <- read.csv('normCoverage_CRC.csv')
faec_mets  <- read.csv('metabolite_data.csv')
agr <- read.csv('AGORA1_metabolites_abbr.csv')

sampnames <- colnames(dnormcov)[2:ncol(dnormcov)]
sampnames <- setdiff(sampnames, sampnames[grep("[0-9]_1$", sampnames)])

smplsrem <- colnames(dnormcov)[colSums(dnormcov != 0) < 25]
sampnames <- setdiff(sampnames,smplsrem)

div <- data.frame(ID=sampnames)
# Function to extract the substring after the last underscore
extract_last <- function(string) {
  parts <- strsplit(string, "_")[[1]]
  substring <- tail(parts, n = 1)
  
  # Check if the substring ends with "_L" or "_R"
  if (substring == "_L" || substring == "_R" || substring == "L" || substring == "R") {
    substring <- paste(parts[length(parts) - 1], substring, sep = "_")
  } else if (endsWith(substring, "_L") || endsWith(substring, "_R")) {
    substring <- substring
  }
  return(substring)
}

# Apply the function to the vector of strings
result <- sapply(dsecr$X, extract_last)
mets <- unique(result)

# Calculate redundancy measures
for (i in mets){
  
  metexch <- dsecr$X[grepl(paste0("_", i, "$"), dsecr$X)]
  # Sample-based, Reference-based, and Abundance-based 
  kl_sample <- c()
  kl_reference <- c()
  kl_abundance <- c()
  kl_I <- c()
  
  #Number of species
  nc <- c()
  
  #Sum of species that can produce the metabolite of interest
  sNB <- c()
  
  #Sum of exchange fluxes 
  sIEX <- c()
  
  for(j in sampnames){
    
    #readout exchange reactions
    IEX <- (-1)*dsecr[which(dsecr$X %in% metexch), j]
    
    if(all(IEX == 0) | sum(IEX < 0) > 0){
      kln_sample <- NaN
      kln_reference <- NaN
      kln_abundance <- NaN
      Ki_In <- NaN
      
    } else {
      
      nb <- c()
      for(k in 1:length(IEX)){
        nbz <- which(dnormcov$X == paste0("pan",gsub(paste0("_", i, "$"), "" ,metexch[k])))
        nb <- c(nb, dnormcov[nbz, j])
      }
      
      leftspec <- rep(0, sum(dnormcov[, j]>0)-length(IEX[IEX>0]))
      samplelevel <- c(IEX[IEX>0], leftspec)
      
      P <- samplelevel/sum(samplelevel)
      Q <- rep(1/length(samplelevel), length(samplelevel))
      x <- rbind(P,Q)
      kln_sample <- (-1)*suppressMessages(KL(x, test.na = TRUE, unit = "log", est.prob = NULL, epsilon = 1e-05))
      gcn_sample <- (-1)*dineq::gini.wtd(samplelevel)
      
      P2 <- IEX/sum(IEX)
      Q2 <- rep(1/length(IEX), length(IEX))
      x2 <- rbind(P2,Q2)
      kln_reference <- (-1)*suppressMessages(KL(x2, test.na = TRUE, unit = "log", est.prob = NULL, epsilon = 1e-05))
      gcn_reference <- (-1)*dineq::gini.wtd(IEX)
      
      P3 <- IEX/sum(IEX)
      Q3 <- nb
      x3 <- rbind(P3,Q3)
      
      kln_abundance <-  (-1)*suppressMessages(KL(x3, test.na = TRUE, unit = "log", est.prob = NULL, epsilon = 1e-05))
      
      P4 <- c(IEX[IEX>0]/sum(IEX))
      Q4 <- nb[IEX>0]/sum(nb)
      x3 <- rbind(P4,Q4)
      kl_In <-  suppressMessages(KL(x3, test.na = TRUE, unit = "log", est.prob = NULL, epsilon = 1e-05))
      
    }
    
    ncn <- sum(IEX>0)
    sIEXn <- sum(IEX)
    
    #Sum of species that can produce the metabolite of interest
    sNBn <- sum(nb)
    
    kl_sample <- c(kl_sample, kln_sample)
    kl_reference <- c(kl_reference, kln_reference)
    kl_abundance <- c(kl_abundance, kln_abundance)
    kl_I <- c(kl_I, kl_In)
    nc <- c(nc, ncn)
    sNB <- c(sNB, sNBn)
    sIEX <- c(sIEX, sIEXn)
    
  }
  brkp <- 612*0.8
  if(sum(!is.na(kl_sample))>brkp){
    
    div[, paste0("KL_Sample_", i)] <- kl_sample
    div[, paste0("KL_Reference_", i)] <- kl_reference
    div[, paste0("KL_Abundance_", i)] <- kl_abundance
    div[, paste0("KLI_", i)] <- kl_I

    div[, paste0("N_", i)] <- nc
    div[, paste0("sIEX_", i)] <- sIEX
    div[, paste0("sNB_", i)] <- sNB
    div[, paste0("sIEX_", i)] <- sIEX
  }
  print(i)
}

dtn <- as.data.frame(t(dnormcov))
colnames(dtn) <- dtn[1, ]

dtn <- dtn[-c(1), ]
dtn$ID <- rownames(dtn)

rownames(dtn) <- NULL
dtn <- dtn[, c(ncol(dtn), 1:ncol(dtn)-1)]
dtn[2:ncol(dtn)] <- sapply(dtn[2:ncol(dtn)], as.numeric)

# ecological diversity measures
divsh <- c()
nospec <- c()

for(i in 1:nrow(dtn)){
  
  specpres <- sum(dtn[i,2:ncol(dtn)] > 0)
  sh <- vegan::diversity(dtn[i,2:ncol(dtn)], index="shannon")
  nospec <- c(nospec, specpres)
  divsh <- c(divsh, sh)
  
}

divdf <- data.frame(ID = dtn$ID, spec = nospec, shannon = divsh)
dmet$healthstat <- 0
dmet$healthstat[dmet$Stratification == "CRC"] = 1

colnames(dmet)[1] <- "ID"

dn <- merge(dmet, divdf, by ="ID")
d <- merge(dn, div, by="ID")

d[13:ncol(d)] <- sapply(d[13:ncol(d)], as.numeric)
d$BMI <- as.numeric(d$BMI)
d$Alcohol <- as.numeric(d$Alcohol)


# faec_mets <-> VMH IDs
agr <- agr[-which(agr$keggId %in% agr$keggId[duplicated(agr$keggId)]), ]
b <- gsub("_.*","", faec_mets$X)

for (i in 1:length(b)){
  if(b[i] %in% gsub("_", "", agr$keggId)){
    b[i] = agr$abbreviation[match(b[i], gsub("_", "", agr$keggId))]
  }
}

faec_mets$X = b
faec_mets <- faec_mets[which(faec_mets$X %in% mets), ]
sampn <- gsub("X", "", colnames(faec_mets)[2:ncol(faec_mets)])

colnames(faec_mets)[2:ncol(faec_mets)] <- gsub("^", "sample_", sampn)
faec_mets <- data.frame(t(faec_mets))
colnames(faec_mets) <- faec_mets["X", ]

faec_mets <- faec_mets[-c(1), ]
faec_mets$ID <- rownames(faec_mets)
faec_mets <- faec_mets[, c(ncol(faec_mets), 1:ncol(faec_mets)-1)]
rownames(faec_mets) <- NULL

d <- merge(d, faec_mets, by="ID", all.x = TRUE)
d[12:ncol(d)] <- sapply(d[12:ncol(d)], as.numeric)

# Get column names that start with "KL_" or "KLI_"
selected_columns <- grep("^KL_|^KLI_", names(d), value = TRUE)

# Subset the dataframe to keep only these columns
d <- d[, c("ID", selected_columns)]
colnames(d)[1] <- "Sample"

#as.data.frame(mets)

pattern <- paste0("(", paste(mets, collapse = "|"), ")$")
filtered_dsecr <- dsecr %>% filter(grepl(pattern, X))
write.csv(filtered_dsecr, 'Tables/MicrobeSecretions.csv', row.names = FALSE)
write.csv(d, 'Tables/RedundancyMeasures.csv', row.names = FALSE)
