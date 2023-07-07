### Script for heatmaps - Fig 1A & FigS 4C
library(ape)
library(phyloseq)
library(RColorBrewer)
library(pheatmap)

# Function for annotation  
# This customized function tax.lab() was created by Katie Lennard available  at https://gist.github.com/kviljoen/97d36c689c5c9b9c39939c7a100720b9
#ARGUMENTS:
#otus = otu table of interest (often subset)
#physeq = phyloseq object of interst (for tax annotation)
#labrow: TRUE = annotate rows with taxonomy info, FALSE = annotate rows with OTU info.
#if OTUs have been merged, TRUE, else FALSE
tax.lab <- function(otus, physeq, labrow = TRUE, merged = FALSE){
  tax <- data.frame(tax_table(physeq))
  tax <- tax[c(rownames(otus)), ]
  tax[is.na(tax)] <- ""
  if(merged == FALSE){
    tax$unique <- rep("", dim(tax)[1])#make unique names
    for(i in 1:dim(tax)[1]){
      x = 7#column 7 = Species
      tax$unique[i] <- ifelse(tax[i, x] == ""|is.na(tax[i, x]) == TRUE, 
                              ifelse(tax[i, x-1] == ""|is.na(tax[i, x-1]) == TRUE, 
                                     ifelse(tax[i, x-2] == ""|is.na(tax[i, x-2]) == TRUE, 
                                            ifelse(tax[i, x-3] == ""|is.na(tax[i, x-3]) == TRUE, #phylum
                                                   ifelse(tax[i, x-4] == ""|is.na(tax[i, x-4]) == TRUE, as.character(tax[i, x-5]), 
                                                          as.character(tax[i, x-4])), #class
                                                   as.character(tax[i, x-3])), #order
                                            as.character(tax[i, x-2])), #family
                                     as.character(tax[i, x-1])), #genus
                              as.character(tax[i, x]))#species
    }
    #Custom cleanup - species:
    #For OTUs with species-level annotation: get the first character of the Genus and append to species.
    #first remove all "[]"brackets
    tax$Genus = sub("\\[", "", tax$Genus)
    tax$Genus = sub("\\]", "", tax$Genus)
    for(i in 1:dim(tax)[1]){
      if(is.na(tax$Species[i]) == FALSE & tax$Species[i]!= ""){#get OTUs with species-level annotation
        get = substr(tax$Genus[i], 1, 1) # Extract substrings in a character vector.
        tax$unique[i] = sub(tax$unique[i], paste0(get, ".", tax$unique[i]), tax$unique[i]) # sub(old, new, string)
      }	
    }
    tax$unique <- make.unique(tax$unique) #Make Character Strings Unique
    
    if(labrow == TRUE){   #labrow: TRUE = annotate rows with taxonomy info
      labs = tax$unique
    }else if(labrow == FALSE){   #labrow: FALSE = annotate rows with OTU info.
      labs = rownames(otus)
    }
  }else if(length(merged)>0){
    if(labrow == TRUE){
      labs = tax[, merged]
    }else if(labrow == FALSE){
      labs = rownames(otus)
    }
    
  }
  return(labs)
}


### Fig.1A  -------
# heatmap function
heatmap_stand_otu <- function(physeq) {
  set.seed(2)
  glom_random_tree = rtree(ntaxa(physeq), rooted = TRUE, tip.label = taxa_names(physeq))
  glom_tree = merge_phyloseq(physeq, glom_random_tree)
  diss <- phyloseq::distance(glom_tree, method = "bray", type = "samples")
  #top 30 most abundant taxa for plotting
  prune_30 = prune_taxa(names(sort(taxa_sums(physeq), TRUE))[1:30], physeq) 
  otus <- otu_table(prune_30)
  #log2 transform otus for better heatmap visualisation colour scale
  zs <- which(otus == 0)#find zero entries 
  if(length(zs)>0){#if there are zero values, transform
    otus[zs] <- 0.1 #set zs to 0.1
    otus <- log2(otus)
  } else{otus <- log2(otus)}
  plot <- otus #data to be plotted
  
  #get metadata for annotation 
  meta_hm <- data.frame(row.names = sample_names(physeq), sample_data(physeq)$Pam3_bray, sample_data(physeq)$Status2, sample_data(physeq)$Study_site)
  names(meta_hm) <- c("Community cluster", "HIV exposure status", "Study site")
  #taxa annotation
  labs <- tax.lab(physeq, otus = otus)
  #color
  col.pal <- brewer.pal(9, "YlOrRd")
  annot.color.col <- list("Community cluster" = c("Cluster 1" = "#FDCA89", "Cluster 2" = "#F1948A", "Cluster 3" = "#E06666"), "HIV exposure status" = c("iHEU" = "maroon1", "iHUU" = "#F0E161"), "Study site" = c("South Africa" = "red", "Nigeria" = "blue"))
  #heatmap
  hm <- pheatmap(plot, annotation = meta_hm, color = col.pal, clustering_method = "complete", clustering_distance_cols = diss, show_colnames = FALSE, labels_row = labs, annotation_colors = annot.color.col, fontsize = 11)
  return(hm)
} 

# load phyloseq object
Phy.f_std <- readRDS("Phy.f_std.RDS")
Phy.f_std_W1 <- subset_samples(Phy.f_std, Visit == "Week 1")

# merge the same species with different ASVs
Phy.f_std_W1_merge1 = merge_taxa(Phy.f_std_W1, c("ASV12", "ASV22")) #Staphylococcus saprophyticus
Phy.f_std_W1_merge2 = merge_taxa(Phy.f_std_W1_merge1, c("ASV18", "ASV68")) #Staphylococcus haemolyticus

# plot heatmap
Fig.1A <- heatmap_stand_otu(Phy.f_std_W1_merge2)
Fig.1A  



### Fig.S4C  -------
heatmap_stand_otu <- function(physeq) {
  set.seed(2)
  glom_random_tree = rtree(ntaxa(physeq), rooted = TRUE, tip.label = taxa_names(physeq))
  glom_tree = merge_phyloseq(physeq, glom_random_tree)
  diss <- phyloseq::distance(glom_tree, method = "bray", type = "samples")
  #top 25 most abundant taxa for plotting
  prune_25 = prune_taxa(names(sort(taxa_sums(physeq), TRUE))[1:24], physeq)  # 24 instead of 25 as two ASVs will be merged together
  otus <- otu_table(prune_25)
  #log2 transform otus for better heatmap visualisation colour scale
  zs <- which(otus == 0)#find zero entries 
  if(length(zs)>0){#if there are zero values, transform
    otus[zs] <- 0.1 #set zs to 0.1
    otus <- log2(otus)
  } else{otus <- log2(otus)}
  plot <- otus #data to be plotted
  #get metadata for annotation 
  meta_hm <- data.frame(row.names = sample_names(physeq), sample_data(physeq)$Pam3_bray, sample_data(physeq)$Status2, sample_data(physeq)$Visit, sample_data(physeq)$Study_site)
  names(meta_hm) <- c("Community cluster", "HIV exposure status", "Visit", "Study site")
  #taxa annotation 
  labs <- tax.lab(physeq, otus = otus)
  #add subspecies info for B.longum 
  labs[labs == "B.longum"] <- "B.longum (subsp. infantis)" # ASV1
  labs[labs == "B.longum.1"] <- "B.longum (subsp. longum)" # ASV4
  #color
  col.pal <- brewer.pal(9, "YlOrRd")
  annot.color.col <- list("Community cluster" = c("Cluster 1" = "#8E999F", "Cluster 2" = "#F1948A", "Cluster 3" = "#C5E1A5"), "HIV exposure status" = c("iHEU" = "#B05A7A", "iHUU" = "#F99417"), "Visit" = c("Week 1" = "#4FC3F7", "Week 15" = "#48648F"), "Study site" = c("South Africa" = "#E87D72", "Nigeria" = "#54BCC2"))
  #heatmap
  hm <- pheatmap(plot, annotation = meta_hm, color = col.pal, clustering_method = "complete", clustering_distance_cols = diss, show_colnames = FALSE, labels_row = labs, annotation_colors = annot.color.col, fontsize = 13)
  return(hm)
} 

# load phyloseq object
Phy.f_std <- readRDS("Phy.f_std.RDS")

# merge the same species with different ASVs
Phy.f_std_merge = merge_taxa(Phy.f_std, c("ASV12", "ASV22")) # Staphylococcus saprophyticus

# plot heatmap
Fig.S4C <- heatmap_stand_otu(Phy.f_std_merge) 
Fig.S4C
