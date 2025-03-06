library(ggplot2)
library(vegan)
library(philentropy)
library(nnet)
library(haven)
library(ggplot2)

setwd(r'(F:\Functional_Redundancy_CSBJ\R\ReSubmission\IBD)')

# Invidiual maximum secretion fluxes computed by predictmicrobecontributions function
# of the cobratoolbox
dsecr <- read.csv('MicrobeContributions_Fluxes.csv')
dsecr[,2:ncol(dsecr)] <- (-1)*dsecr[,2:ncol(dsecr)]

# Metadata of the IBD study
dmet <- read.csv('metadata_IBD.csv')

# Normalized abundances for computing microbiome community models
dnormcov <- read.csv('normalizedCoverage_IBD.csv')

# File tax contributions
dtax <- read.csv('Taxonomy_MicrobeContributions_Fluxes.csv')
mets <- dtax$Metabolite

string_counts <- table(mets)
sorted_counts <- sort(string_counts, decreasing = TRUE)
mets <- names(sorted_counts[sorted_counts>5])

sampnames <- colnames(dnormcov)[2:ncol(dnormcov)]

smplsrem <- colnames(dnormcov)[colSums(dnormcov != 0) < 25]
sampnames <- setdiff(sampnames,smplsrem)

div <- data.frame(ID=sampnames)

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
        nbz <- which(dnormcov$X == re)
        nb <- c(nb, dnormcov[nbz, j])
      }

      abundance <- dnormcov[, j]
      abundance <- abundance/sum(abundance)
      f <- rep(0, length(abundance))
      for(s in 1:length(abundance)){

        secr_spec <- paste0(gsub("^pan", "", dnormcov$X[s]), "_IEX_", i)
        if(secr_spec %in% metexch){
          ind <- which(metexch == secr_spec)
          f[s] <- IEX[ind]
        } else{
          f[s] <- 0
        }
      }

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
  brkp <- 102*0.8
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
colnames(dmet)[1] <- "ID"

dn <- merge(dmet, divdf, by ="ID")
d <- merge(dn, div, by="ID")
d <- merge(d, dtn, by="ID")

d[3:ncol(d)] <- sapply(d[3:ncol(d)], as.numeric)

#as.data.frame(mets)
write.csv(d, 'RedundancyMeasures_IBD.csv', row.names = FALSE)

#S3
# Get column names that start with "KL_" or "KLI_"
selected_columns <- grep("^KL_|^KLI_", names(d), value = TRUE)

# Subset the dataframe to keep only these columns
df <- d[, c("ID", selected_columns)]
colnames(df)[1] <- "Sample"
df[is.na(df)] <- "NA"


write.xlsx(df, 'Tables/S3_RedundancyMeasuresIBD.xlsx', rowNames = FALSE)
