# Scripts for pre-processing phyloseq object
#-----------------------------------------------------------------------------------------------------
# Load packages  ----
#-----------------------------------------------------------------------------------------------------
library(dplyr)
library(writexl)
library(writexl)
library(decontam)
library(ggplot2)
library(ggrepel)
library(cluster)
library(RColorBrewer)
library(microbiome)

#-----------------------------------------------------------------------------------------------------
# Remove the taxa apart from bacteria  ----
#-----------------------------------------------------------------------------------------------------
BEAMING_phy <- readRDS("BEAMING_phy.RDS") # this phyloseq includes 16S rRNA samples for BEAMING study (540 samples from SA and Nigerian samples & pc/nc)
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 5464 taxa and 540 samples ]
#sample_data() Sample Data:       [ 540 samples by 75 sample variables ]
#tax_table()   Taxonomy Table:    [ 5464 taxa by 7 taxonomic ranks ]
#refseq()      DNAStringSet:      [ 5464 reference sequences ]

## any samples with no OTU
any(taxa_sums(BEAMING_phy) == 0)  # TRUE
taxa_sums(BEAMING_phy) 

# remove OTUs with 0 or 1 (singletone)
BEAMING_phy2 <- prune_taxa(taxa_sums(BEAMING_phy) > 1, BEAMING_phy)
any(taxa_sums(BEAMING_phy2) == 0) # FALSE

# filter out non-bacteria, chloroplasts and mitochondria 
a = which(tax_table(BEAMING_phy2)[, "Kingdom"]! = "Bacteria")
b = which(tax_table(BEAMING_phy2)[, "Family"] == "mitochondria")
c = which(tax_table(BEAMING_phy2)[, "Class"] == "Chloroplast")
d = c(a, b, c)  
d  # [1]  300  835 1121 1568 1677 1798 3225 3862 3986
e = taxa_names(BEAMING_phy2)[-d]
BEAMING_phy3 = prune_taxa(e, BEAMING_phy2)

# check
ntaxa(BEAMING_phy)  # 5464 - before filtering
ntaxa(BEAMING_phy2)  # 4350 - after removing singletone
ntaxa(BEAMING_phy3) # 4341 - after removing non-bacteria

#-----------------------------------------------------------------------------------------------------
# Inspect NC and PC for possible contamination ----
#-----------------------------------------------------------------------------------------------------
sample_data(BEAMING_phy3)$Sample_or_Control <-  case_when((sample_data(BEAMING_phy3)$PID %in% c("M33-NC1", "M33-NC2", "M33-NC3", "M35-NC1", "M35-NC2", "M35-NC3", 
                                                                                                "M37-NC0", "M37-NC1", "M37-NC2", "M38-NC0", "M38-NC1", "M38-NC2")) ~"NegControl", 
                                                          (sample_data(BEAMING_phy3)$PID %in% c("M33-PC", "M35-PC", "M37-PC", "M38-PC")) ~ "PosControl", 
                                                          TRUE ~ "Sample")

### subset NC samples  ------
NC_phy <- subset_samples(BEAMING_phy3, Sample_or_Control == "NegControl")
otu_table(NC_phy)

# prune OTUs
NC_phy2 <- prune_taxa(taxa_sums(NC_phy) > 1, NC_phy)
NC_phy2

sample_sums(NC_phy2)
## M33-NC1 M33-NC2 M33-NC3 M35-NC1 M35-NC2 M35-NC3 M37-NC0 M37-NC1 M37-NC2 M38-NC0 M38-NC1 M38-NC2 
## 228     157      61     159      47      99      17      28       0       0       0       1 

# standarize the reads count
total = median(sample_sums(NC_phy2))
standf = function(x, t = total) round(t * (x / sum(x)))
NC_phy2.std = transform_sample_counts(NC_phy2, standf)

# plot
plot_bar(NC_phy2.std, x = "PID", fill = "Phylum")
plot_bar(NC_phy2.std, x = "PID", fill = "Class")
plot_bar(NC_phy2.std, x = "PID", fill = "Order")
plot_bar(NC_phy2.std, x = "PID", fill = "Family")
plot_bar(NC_phy2.std, x = "PID", fill = "Genus")
plot_bar(NC_phy2.std, x = "PID", fill = "Species")


### subset PC samples ------
PC_phy <- subset_samples(BEAMING_phy3, Sample_or_Control == "PosControl")
otu_table(PC_phy) 

# prune OTUs
PC_phy2 <- prune_taxa(taxa_sums(PC_phy) > 1, PC_phy)
PC_phy2

sample_sums(PC_phy2)

# standarize the reads count
total = median(sample_sums(PC_phy2))
standf = function(x, t = total) round(t * (x / sum(x)))
PC_phy2.std = transform_sample_counts(PC_phy2, standf) # standarized phyloseq object

sample_sums(PC_phy2) # before standarize
sample_sums(PC_phy2.std) # after standarize

# plot
plot_bar(PC_phy2.std, x = "PID", fill = "Phylum")
plot_bar(PC_phy2.std, x = "PID", fill = "Class")
plot_bar(PC_phy2.std, x = "PID", fill = "Order")
plot_bar(PC_phy2.std, x = "PID", fill = "Family")
plot_bar(PC_phy2.std, x = "PID", fill = "Genus")
plot_bar(PC_phy2.std, x = "PID", fill = "Species")

#-----------------------------------------------------------------------------------------------------
# Remove potential contamination by "decontam"  ----
#-----------------------------------------------------------------------------------------------------
# followed decontam tutorial: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
BEAMING_phy3 

## ----see-meta-table--------------------------------------------------------
head(sample_data(BEAMING_phy3))
## ----see-depths------------------------------------------------------------
# annotate PC/NC samples
sample_data(BEAMING_phy3)$Sample_or_Control <-  case_when((sample_data(BEAMING_phy3)$PID %in% c("M33-NC1", "M33-NC2", "M33-NC3", "M35-NC1", "M35-NC2", "M35-NC3", 
                                                                                                "M37-NC0", "M37-NC1", "M37-NC2", "M38-NC0", "M38-NC1", "M38-NC2")) ~"NegControl", 
                                                          (sample_data(BEAMING_phy3)$PID %in% c("M33-PC", "M35-PC", "M37-PC", "M38-PC")) ~ "PosControl", 
                                                          TRUE ~ "Sample")

df <- as.data.frame(sample_data(BEAMING_phy3))
df$LibrarySize <- sample_sums(BEAMING_phy3)
df <- df[order(df$LibrarySize), ] # change the order
df$Index <- seq(nrow(df))

sample_sums(BEAMING_phy3) # check indivisual read count
ggplot(data = df, aes(x = Index, y = LibrarySize, color = Sample_or_Control)) + geom_point()  

## ----Identify Contaminants - Prevalence------------------------------------------------------------
sample_data(BEAMING_phy3)$is.neg <- sample_data(BEAMING_phy3)$Sample_or_Control == "NegControl"
contamdf.prev <- isContaminant(BEAMING_phy3, method = "prevalence", neg = "is.neg")
table(contamdf.prev$contaminant) # whether they are contaminant or not (true/false)
## FALSE  TRUE 
## 4339     2 

rownames(contamdf.prev)[contamdf.prev$contaminant == "TRUE"] # list of contaminant ASVs
tax_table(BEAMING_phy3) [17] # Staphylococcus sciuri
tax_table(BEAMING_phy3) [2633] # Cyanobacteriia  (class) Chloroplast (order)

## ----prevalence-05---------------------------------------------------------
# stricter filtering with threshold 0.5
contamdf.prev05 <- isContaminant(BEAMING_phy3, method = "prevalence", neg = "is.neg", threshold = 0.5)
table(contamdf.prev05$contaminant)
# FALSE  TRUE 
# 4335     6 
head(which(contamdf.prev05$contaminant)) # [1]   14   17   50  963 1776 2633 

rownames(contamdf.prev05)[contamdf.prev05$contaminant == "TRUE"] # list of contaminant ASVs

tax_table(BEAMING_phy3) [14]  # Kocuria marina"
tax_table(BEAMING_phy3) [17]  # Staphylococcus sciuri
tax_table(BEAMING_phy3) [50]  # Bifidobacterium breve
tax_table(BEAMING_phy3) [963] # Staphylococcus sciuri
tax_table(BEAMING_phy3) [1776]  # Halomonas nitritophilus
tax_table(BEAMING_phy3) [2633]  # Cyanobacteriia  (class) Chloroplast (order)

## ----see-prev-05-----------------------------------------------------------
# make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(BEAMING_phy3, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "NegControl", ps.pa) 
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Sample", ps.pa)

# make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos = taxa_sums(ps.pa.pos), pa.neg = taxa_sums(ps.pa.neg), 
                    contaminant = contamdf.prev$contaminant)
df.pa[1:10, ]

# with annotation
ggplot(data = df.pa, aes(x = pa.neg, y = pa.pos, color = contaminant)) + 
  geom_point() + xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") + 
  geom_text_repel(aes(label = rownames(df.pa)))

## ----remove----------------------------------------------------------------
BEAMING_phy3  #4341 taxa (before removing the contaminant)

BEAMING_phy3.noncontam <- prune_taxa(!contamdf.prev05$contaminant, BEAMING_phy3)
BEAMING_phy3.noncontam # 4335 taxa (after removing the contaminant)

#-----------------------------------------------------------------------------------------------------
# Pre Processing the phyloseq object ----
#----------------------------------------------------------------------------------------------------
# prune samples with <2000 reads, PC/NC, and samples that did not fall into the inclusion criteria
reads <- sample_sums(BEAMING_phy3.noncontam)
length(which(reads<2000)) # 91 mostly W1 samples
Samples_toRemove <- dput(names(which(reads<2000)))
Samples_toRemove <- append(Samples_toRemove, c("M33-PC", "M35-PC", "M37-PC", 
                                              "N0240BBAS2", "N0267BV5", "N0319BV5", "N0314BV5")) 

BEAMING_phy4 <- subset_samples(BEAMING_phy3.noncontam, !(Sample_ID2 %in% Samples_toRemove))  

#-----------------------------------------------------------------------------------------------------
#  standardise and filtering, prepare for the analysis ----
#-----------------------------------------------------------------------------------------------------
## standardise
total = median(sample_sums(BEAMING_phy4))
standf = function(x, t = total) round(t * (x / sum(x)))
Phy.std = transform_sample_counts(BEAMING_phy4, standf) # standarized phyloseq object

sample_sums(BEAMING_phy4)
sample_sums(Phy.std) 

## filter the standardised phyloseq ( = Phy.f_std)
# the filter below retains only OTUs that are present at at least 10 counts at least 2% of samples OR that have a total relative abundance of at least 0.1% of the total number of reads/sample
Phy.f_std = filter_taxa(Phy.std, function(x) sum(x > 10) > (0.2*length(x)) | sum(x) > 0.001*total, TRUE) 
ntaxa(Phy.std) #  4335 
ntaxa(Phy.f_std) #  2084 

#saveRDS(Phy.f_std, "Phy.f_std.RDS")
#Phy.f_std <- readRDS("Phy.f_std.RDS") # use this for downstream analysis


#-----------------------------------------------------------------------------------------------------
# PAM clustering (determine the k value ) ----
#-----------------------------------------------------------------------------------------------------
Phy.f_std <- readRDS("Phy.f_std.RDS")

# cluster into PAMs
otu_table(Phy.f_std) <- t(otu_table(Phy.f_std))

phy.obj.rel <- Phy.f_std
phy.obj.rel <-transform_sample_counts(physeq = phy.obj.rel, fun = function(x) x/sum(x))
otu_table(phy.obj.rel)

jsd_dist <- phyloseq::distance(phy.obj.rel, method = "jsd")
ord = ordinate(phy.obj.rel, method = "PCoA", distance = jsd_dist)
plot_scree(ord) + xlim(as.character(seq(1, 20))) + ggtitle("PCoA-jsd ordination eigenvalues")

evs <- ord$value$Eigenvalues
print(evs[1:20])
print(tail(evs))

# remove those below magnitude of largest negative eignevalue
h_sub10 <- hist(evs[10:length(evs)], 100)
plot(h_sub10$mids, h_sub10$count, log = "y", type = 'h', lwd = 10, lend = 2)

# gap statistic for cluter number
NDIM <- 7 #11 falls below mag of largest neg eig
x <- ord$vectors[, 1:NDIM]  # rows = sample, cols = MDS axes, entries = value
pamPCoA = function(x, k) {
  list(cluster = pam(x[, 1:2], k, cluster.only = TRUE))
}
gs = clusGap(x, FUN = pamPCoA, d.power = 2, K.max = 12, B = 500) #Tibshirani def of power

# ElbowPlot ---
# k = 3 is the optimal
p1 = plot_clusgap(gs) + scale_x_continuous(breaks = c(seq(0, 12, 2))) + 
  theme_bw() + theme() + ggtitle("") + theme(axis.text = element_text(size = 13, face = "bold"), 
                                            axis.title = element_text(size = 16, face = "bold"), 
                                            legend.title = element_text(size = 14), 
                                            legend.text = element_text(size = 12))
p1

# k = 3
K <- 3
x <- ord$vectors[, 1:NDIM]
clust <- as.factor(pam(x, k = K, cluster.only = T))
# add to sample data
sample_data(phy.obj.rel)$PAM_3 <- clust
PAMs <- as.character(seq(K))

# k = 4
K <- 4
x <- ord$vectors[, 1:NDIM]
clust <- as.factor(pam(x, k = K, cluster.only = T))
# add to sample data
sample_data(phy.obj.rel)$PAM_4 <- clust
PAMs <- as.character(seq(K))

# k = 5
K <- 5
x <- ord$vectors[, 1:NDIM]
clust <- as.factor(pam(x, k = K, cluster.only = T))
# add to sample data
sample_data(phy.obj.rel)$PAM_5 <- clust
PAMs <- as.character(seq(K))

# k = 6
K <- 6
x <- ord$vectors[, 1:NDIM]
clust <- as.factor(pam(x, k = K, cluster.only = T))
# add to sample data
sample_data(phy.obj.rel)$PAM_6 <- clust
PAMs <- as.character(seq(K))

# k = 7
K <- 7
x <- ord$vectors[, 1:NDIM]
clust <- as.factor(pam(x, k = K, cluster.only = T))
# add to sample data
sample_data(phy.obj.rel)$PAM_7 <- clust
PAMs <- as.character(seq(K))

### visually evaluate clustering -----
### set up the color
PAMColors <- brewer.pal(7, "Paired")[c(1, 3, 2, 5, 4, 6, 7)] # Length 6 for consistency with pre-revision CST + coloration
names(PAMColors) <- PAMs
PAMColorScale <- scale_colour_manual(name = "PAM_3", values = PAMColors[1:7])
PAMFillScale <- scale_fill_manual(name = "PAM_3", values = PAMColors[1:7])
PAMColorScale <- scale_colour_manual(name = "PAM_4", values = PAMColors[1:7])
PAMFillScale <- scale_fill_manual(name = "PAM_4", values = PAMColors[1:7])
PAMColorScale <- scale_colour_manual(name = "PAM_5", values = PAMColors[1:7])
PAMFillScale <- scale_fill_manual(name = "PAM_5", values = PAMColors[1:7])
PAMColorScale <- scale_colour_manual(name = "PAM_6", values = PAMColors[1:7])
PAMFillScale <- scale_fill_manual(name = "PAM_6", values = PAMColors[1:7])
PAMColorScale <- scale_colour_manual(name = "PAM_7", values = PAMColors[1:7])
PAMFillScale <- scale_fill_manual(name = "PAM_7", values = PAMColors[1:7])

### plot
# k = 3
plot_ordination(phy.obj.rel, ord, color = "PAM_3") + PAMColorScale
plot_ordination(phy.obj.rel, ord, axes = c(3, 4), color = "PAM_3") + PAMColorScale
# k = 4
plot_ordination(phy.obj.rel, ord, color = "PAM_4") + PAMColorScale
# k = 5
plot_ordination(phy.obj.rel, ord, color = "PAM_5") + PAMColorScale
# k = 6
plot_ordination(phy.obj.rel, ord, color = "PAM_6") + PAMColorScale
# k = 7
plot_ordination(phy.obj.rel, ord, color = "PAM_7") + PAMColorScale

### add the cluster annotation (k = 3)
K <- 3
x <- ord$vectors[, 1:NDIM]
clust <- as.factor(pam(x, k = K, cluster.only = T))

sample_data(Phy.f_std)$Pam3_bray <- clust   # add the annotation to phyloseq object
levels(sample_data(Phy.f_std)$Pam3_bray) <- c("1" = "Cluster1", "2" = "Cluster2", "3" = "Cluster3")

# save the PAM cluter info
# saveRDS(Phy.f_std, "Phy.f_std.RDS")

#-----------------------------------------------------------------------------------------------------
# Shannon diversity info ----
#-----------------------------------------------------------------------------------------------------
Phy.f_std <- readRDS("Phy.f_std.RDS")
tab <- microbiome::alpha(Phy.f_std, index = "shannon")
sample_data(Phy.f_std)$Shannon <- tab$diversity_shannon 

# save the shannon diversity info info
# saveRDS(Phy.f_std, "Phy.f_std.RDS")