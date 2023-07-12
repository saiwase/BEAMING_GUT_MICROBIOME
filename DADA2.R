# Script for DADA2
#-----------------------------------------------------------------------------------------------------
# Load packages ----
#-----------------------------------------------------------------------------------------------------
library(dplyr)
library(writexl)
library(rlist)
library(readxl)
library(stringr)
library(dada2)
library(DECIPHER)
library(phangorn)

#-----------------------------------------------------------------------------------------------------
# DADA2 for BEAMING samples  ----
#-----------------------------------------------------------------------------------------------------
# followed the DADA2 tutorial: https://benjjneb.github.io/dada2/tutorial.html

path <- "/path/to/fastq.files/"   # define the path
list.files(path)  #list the file names

# perform some string manipulation
fnFs <- sort(list.files(path,pattern="_R1_001.fastq",full.names=TRUE))
fnRs <- sort(list.files(path,pattern="_R2_001.fastq",full.names=TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names 
length(sample.names)
sample.names[duplicated(sample.names)]

#-----------------------------------------------------------------------------------------------------
# Quality control  ----
#-----------------------------------------------------------------------------------------------------
plotQualityProfile(fnFs[40:50])  
plotQualityProfile(fnRs[40:50])

#-----------------------------------------------------------------------------------------------------
# Filter and trim  ----
#-----------------------------------------------------------------------------------------------------
filtFs <- file.path(path,"filtered",paste0(sample.names,"_F_filt.fastq.gz"))
filtRs <- file.path(path,"filtered",paste0(sample.names,"_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# standard filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and maxEE=2
# truncLen must be large enough to maintain 20 + biological.length.variation nucleotides between them
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(265,210),
                     maxN=0, maxEE=c(2,2), truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)
# save the result
out_data <- out 
df <- data.frame(sample.names, out_data, check.names = T)
print(df)
write_xlsx(df,"./DADA2_output.xlsx")

#-----------------------------------------------------------------------------------------------------
# Learn the Error Rates  ----
#-----------------------------------------------------------------------------------------------------
# "learnErrors method" learns this error model from the data, by alternating estimation of the error rates and inference of sample composition until they converge on a jointly consistent solution.
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs,multithread = TRUE)

# check the error rates
plotErrors(errF, nominalQ = TRUE)

#-----------------------------------------------------------------------------------------------------
# Sample Inference  ----
#-----------------------------------------------------------------------------------------------------
# apply the core sample inference algorithm to the filtered and trimmed sequence data.
dadaFs <-dada(filtFs, err = errF, multithread = TRUE)
dadaRs <-dada(filtRs, err = errR, multithread = TRUE)

# inspecting the result
dadaFs[[1]] 

#-----------------------------------------------------------------------------------------------------
# Merge paired reads  ----
#-----------------------------------------------------------------------------------------------------
# merge the forward and reverse reads together to obtain the full denoised sequences
# by default, only output if the forward and reverse reads overlap by at least 12 bases
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# inspeect the result
head(mergers[[1]])

#-----------------------------------------------------------------------------------------------------
# Construct sequence table  ----
#-----------------------------------------------------------------------------------------------------
# construct an amplicon sequence variant table (ASV) table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)  

table(nchar(getSequences(seqtab)))

#-----------------------------------------------------------------------------------------------------
# Remove chimeras  ----
#-----------------------------------------------------------------------------------------------------
# identify how many chimeric reads in the total sequence
seqtab.nochim <- removeBimeraDenovo(seqtab,method="consensus", multithread=TRUE, verbose = TRUE)

# dimention of the seqtable after removing chimera
dim(seqtab.nochim)  

# ratio of reads that are not chimeric
sum(seqtab.nochim)/sum(seqtab) 

#-----------------------------------------------------------------------------------------------------
# Track reads through the pipeline  ----
#-----------------------------------------------------------------------------------------------------
getN <- function(x)sum(getUniques(x))
track <-cbind(out, sapply(dadaFs, getN),sapply(dadaRs, getN),sapply(mergers, getN),rowSums(seqtab.nochim))

colnames(track) <- c("input","filtered","denoisedF","denisedR","mergers","nochim")
rownames(track) <-sample.names
head(track) 

track_data <- track
df2 <- data.frame(sample.names,track_data, check.names = T)
print(df2)
write_xlsx(df2,"./DADA2_output.xlsx")

#-----------------------------------------------------------------------------------------------------
# Assign Taxonomy  ----
#-----------------------------------------------------------------------------------------------------
# use the updated SILVA training set (version 132) https://github.com/itsmisterbrown/updated_16S_dbs
taxa <- assignTaxonomy(seqtab.nochim,
                       "/path/to/silva_nr99_v138_wSpecies_train_set_BPB_091820_updated.fa.gz", multithread = TRUE)

# inspect the taxonomic assignments
taxa.print <- taxa
View(taxa.print) 
rownames(taxa.print) <-NULL  
head(taxa.print)

# save the taxonomy list
write.csv(taxa, file="ASVs_taxonomy.csv") 
saveRDS(taxa, "ASVs_taxonomy.rds") 

# ASVs count table
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
count.asv.tab <- t(seqtab.nochim)
row.names(count.asv.tab) <- sub(">", "", asv_headers)

# save the ASV count
write.csv(count.asv.tab, file="ASVs_counts.csv")
saveRDS(count.asv.tab, file="ASVs_counts.rds")

