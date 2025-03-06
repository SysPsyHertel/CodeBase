library(vegan)
library(haven)
library(stringr)
library(openxlsx)
library(ggplot2)

## This file computes functional redundancy measures

# Invidiual maximum secretion fluxes computed by predictmicrobecontributions function
# of the cobratoolbox, generated with compute_individual_secretion_fluxes.m 
dsecr <- read.csv('Microbe_Secretion.csv')

# Metadata of the colorectal cancer (CRC) study
dmet <- read.csv('metadata.csv')

# Normalized abundances for the MicrobiomeModellingToolbox pipeline
dnormcov <- read.csv('normCoverage_CRC.csv')

# Generated maximum secretion fluxes by the MicrobiomeModellingToolbox pipeline
dnetsecr <- read.csv('inputDiet_net_secretion_fluxes.csv')

# Faecal concentrations of the CRC study
faec_mets  <- read.csv('metabolite_data.csv')

# AGORA <-> KEGGIds file
agr <- read.csv('AGORA1_metabolites_abbr.csv')

sampnames <- colnames(dnormcov)[2:ncol(dnormcov)]
sampnames <- setdiff(sampnames, sampnames[grep("[0-9]_1$", sampnames)])

smplsrem <- colnames(dnormcov)[colSums(dnormcov != 0) < 25]
sampnames <- setdiff(sampnames,smplsrem)

div <- data.frame(ID=sampnames)

metsec <- gsub("EX_", "", dnetsecr$Net.secretion)
metsec <- gsub("\\[fe\\]", "", metsec)

# Collapse the vector into a single regex pattern
metsec_pattern <- paste(c(metsec, "Lkynr"), collapse = "|")
last_matching_values <- str_extract(dsecr$X, paste0("(", paste(metsec_pattern, collapse = "|"), ")$"))
mets <- unique(last_matching_values)
cd <- c()
for (i in mets){

  metexch <- dsecr$X[grepl(paste0("_", i, "$"), dsecr$X)]
  #Sample-based, Reference-based, and Abundance-based
  kl_sample <- c()
  kl_reference <- c()
  kl_abundance <- c()
  KL_I <- c()

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
      KL_In <- NaN

    } else {

      nb <- c()
      for(k in 1:length(IEX)){
        re <- gsub(paste0("_", i, "$"), "" ,metexch[k])
        re <- gsub("_IEX$", "", re)
        nbz <- which(dnormcov$X == paste0("pan", re))
        nb <- c(nb, dnormcov[nbz, j])
      }

      abundance <- dnormcov[, j]
      f <- rep(0, length(abundance))
      for(s in 1:length(abundance)){

        secr_spec <- paste0(gsub("^pan", "", dnormcov$X[s]), "_", i)
        if(secr_spec %in% metexch){
          ind <- which(metexch == secr_spec)
          f[s] <- IEX[ind]
        } else{
          f[s] <- 0
        }
      }
      abundance <- abundance/sum(abundance)
      result <- FunRed::fredundancy(f, abundance, n_reference=length(IEX))
      kln_sample <- result$sample_based
      kln_reference <- result$reference_based
      kln_abundance <-  result$abundance_based
      KL_In <- result$interdependency
    }

    ncn <- sum(IEX>0)
    sIEXn <- sum(IEX)

    #Sum of species that can produce the metabolite of interest
    sNBn <- sum(nb)
    kl_sample <- c(kl_sample, kln_sample)
    kl_reference <- c(kl_reference, kln_reference)
    kl_abundance <- c(kl_abundance, kln_abundance)
    KL_I <- c(KL_I, KL_In)
    nc <- c(nc, ncn)
    sNB <- c(sNB, sNBn)
    sIEX <- c(sIEX, sIEXn)
  }
  brkp <- 612*0.8
  if(sum(!is.na(kl_sample))>brkp){

  div[, paste0("KL_Sample_", i)] <- kl_sample

  div[, paste0("KL_Reference_", i)] <- kl_reference

  div[, paste0("KL_Abundance_", i)] <- kl_abundance

  div[, paste0("KLI_", i)] <- KL_I

  #div[, paste0("KL_Abundance2_", i)] <- kl_abundance2
  #div[, paste0("GC_Abundance2_", i)] <- gc_abundance2
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
divsimp <- c()
nospec <- c()
evenness <- c()


for(i in 1:nrow(dtn)){

  specpres <- sum(dtn[i,2:ncol(dtn)] > 0)
  sh <- vegan::diversity(dtn[i,2:ncol(dtn)], index="shannon")
  simp <- vegan::diversity(dtn[i,2:ncol(dtn)], index="simpson")
  even <- sh/log(specpres)
  nospec <- c(nospec, specpres)
  divsh <- c(divsh, sh)
  divsimp <- c(divsimp, simp)
  evenness <- c(evenness, even)

}

divdf <- data.frame(ID = dtn$ID, spec = nospec, shannon = divsh, simpson = divsimp, evenness = evenness)
dmet$healthstat <- 0
dmet$healthstat[dmet$Stratification == "CRC"] = 1

colnames(dmet)[1] <- "ID"

dn <- merge(dmet, divdf, by ="ID")
d <- merge(dn, div, by="ID")
d <- merge(d, dtn, by="ID")

d[13:ncol(d)] <- sapply(d[13:ncol(d)], as.numeric)
d$BMI <- as.numeric(d$BMI)
d$Alcohol <- as.numeric(d$Alcohol)

# net secretion
#dnetsecr[nrow(dnetsecr)+1, ] <- c("Metabolic.Richness", colSums(dnetsecr[,2:ncol(dnetsecr)]>0))
dnetsecr <- data.frame(t(dnetsecr))

colnames(dnetsecr) <- dnetsecr[1,]
dnetsecr <- dnetsecr[-1,]

dnetsecr$ID <- rownames(dnetsecr)
row.names(dnetsecr) <- NULL

dt <- merge(d, dnetsecr, by ="ID")

dt <- dt[, colnames(dt)[c(1:16, (ncol(dt)-1):ncol(dt), 17:(ncol(dt)-2))]]

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

d <- merge(dt, faec_mets, by="ID", all.x = TRUE)
d[13:ncol(d)] <- sapply(d[13:ncol(d)], as.numeric)

write.csv(d, 'RedundancyMeasures.csv', row.names = FALSE)

#S1
mets <- names(d)[grep("^KL_Sample", names(d))]
mets <- gsub(".*KL_Sample_", "", mets)
pattern <- paste0("(", paste(mets, collapse = "|"), ")$")
filtered_dsecr <- dsecr %>% filter(grepl(pattern, X))
write.xlsx(filtered_dsecr, 'Tables/S1_MicrobeSecretionsCRC.xlsx', rowNames = FALSE)

#S2
# Get column names that start with "KL_" or "KLI_"
selected_columns <- grep("^KL_|^KLI_", names(d), value = TRUE)

# Subset the dataframe to keep only these columns
df <- d[, c("ID", selected_columns)]
colnames(df)[1] <- "Sample"
df[] <- lapply(df, function(x) {
  if (is.numeric(x)) {
    x[is.nan(x)] <- NA
  }
  return(x)
})

write.xlsx(df, 'Tables/S2_RedundancyMeasuresCRC.xlsx', rowNames = FALSE)
