### script for ANCOMBC analysis
library(phyloseq)
library(ANCOMBC)
library(dplyr)
library(microViz)
library(tibble)


### function for annotation --------
# This customized function tax.lab() was created by Katie Lennard available at https://gist.github.com/kviljoen/97d36c689c5c9b9c39939c7a100720b9
#ARGUMENTS:
#otus = otu table of interest (often subset)
#physeq = phyloseq object of interst (for tax annotation)
#labrow: TRUE = annotate rows with taxonomy info, FALSE = annotate rows with OTU info.
#if OTUs have been merged, TRUE, else FALSE
tax.lab <- function(otus, physeq, labrow=TRUE, merged=FALSE){
  tax <- data.frame(tax_table(physeq))
  tax <- tax[c(rownames(otus)),]
  tax[is.na(tax)] <- ""
  if(merged==FALSE){
    tax$unique <- rep("",dim(tax)[1])#make unique names
    for(i in 1:dim(tax)[1]){
      x = 7#column 7 = Species
      tax$unique[i] <- ifelse(tax[i,x]==""|is.na(tax[i,x])==TRUE,
                              ifelse(tax[i,x-1]==""|is.na(tax[i,x-1])==TRUE,
                                     ifelse(tax[i,x-2]==""|is.na(tax[i,x-2])==TRUE,
                                            ifelse(tax[i,x-3]==""|is.na(tax[i,x-3])==TRUE,#phylum
                                                   ifelse(tax[i,x-4]==""|is.na(tax[i,x-4])==TRUE, as.character(tax[i,x-5]),
                                                          as.character(tax[i,x-4])),#class
                                                   as.character(tax[i,x-3])),#order
                                            as.character(tax[i,x-2])),#family
                                     as.character(tax[i,x-1])),#genus
                              as.character(tax[i,x]))#species
    }
    #Custom cleanup - species:
    #For OTUs with species-level annotation: get the first character of the Genus and append to species.
    #first remove all "[]"brackets
    tax$Genus = sub("\\[","",tax$Genus)
    tax$Genus = sub("\\]","",tax$Genus)
    for(i in 1:dim(tax)[1]){
      if(is.na(tax$Species[i])==FALSE & tax$Species[i]!=""){#get OTUs with species-level annotation
        get = substr(tax$Genus[i],1,1) # Extract substrings in a character vector.
        tax$unique[i] = sub(tax$unique[i],paste0(get,".",tax$unique[i]),tax$unique[i]) # sub(old, new, string)
      }	
    }
    tax$unique <- make.unique(tax$unique) #Make Character Strings Unique
    
    if(labrow==TRUE){   #labrow: TRUE = annotate rows with taxonomy info
      labs = tax$unique
    }else if(labrow==FALSE){   #labrow: FALSE = annotate rows with OTU info.
      labs = rownames(otus)
    }
  }else if(length(merged)>0){
    if(labrow==TRUE){
      labs=tax[,merged]
    }else if(labrow==FALSE){
      labs = rownames(otus)
    }
    
  }
  return(labs)
}


### load phyloseq --------
Phy.f_std <- readRDS("Phy.f_std.RDS")
sample_data(Phy.f_std)$Visit
sample_data(Phy.f_std)$Status2 <-factor(sample_data(Phy.f_std)$Status2, levels = c("iHUU", "iHEU"))
levels(sample_data(Phy.f_std)$Feeding_W15_strictver) 

# subset by Site & Visit
p1 <- subset_samples(Phy.f_std, Study_site == "South Africa" & Visit == "Week 1")
p1 <- subset_samples(Phy.f_std, Study_site == "South Africa" & Visit == "Week 15")
p1 <- subset_samples(Phy.f_std, Study_site == "Nigeria" & Visit == "Week 1")
p1 <- subset_samples(Phy.f_std, Study_site == "Nigeria" & Visit == "Week 15")

otu = otu_table(p1)
labs <- tax.lab(p1, otu = otu)
head(labs)

p2 <- tax_mutate(p1, Species = labs) 

# for SA W1 / Nigerian samples
out = ancombc(data = p2, 
              tax_level = "Species", phyloseq = NULL, 
              formula = "Status2", 
              p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
              group = "Status2", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
              max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE, 
              n_cl = 1, verbose = TRUE)

# adjust with mode of feeding (i.e., for SA W15)
out = ancombc(data = p2, 
              tax_level = "Species", phyloseq = NULL, 
              formula = "Status2 + Feeding_W15_strictver", 
              p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
              group = "Status2", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
              max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE, 
              n_cl = 1, verbose = TRUE)

res = out$res

### Visualization --------
res = out$res
colnames(res$lfc) 

df_lfc = data.frame(res$lfc * res$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")
df_se = data.frame(res$se * res$diff_abn, check.names = FALSE) %>% 
  rownames_to_column("taxon_id")
colnames(df_se)[-1] = paste0(colnames(df_se)[-1], "SE")

df_fig_SS = df_lfc %>% 
  dplyr::left_join(df_se, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, Status2iHEU, Status2iHEUSE) %>%
  dplyr::filter(Status2iHEU > 0.5 | Status2iHEU < -0.5) %>% 
  dplyr::arrange(desc(Status2iHEU)) %>%
  dplyr::mutate(direct = ifelse(Status2iHEU > 0, "Positive LFC", "Negative LFC"))
df_fig_SS$taxon_id = factor(df_fig_SS$taxon_id, levels = df_fig_SS$taxon_id)
df_fig_SS$direct = factor(df_fig_SS$direct, 
                          levels = c("Positive LFC", "Negative LFC"))

tax_list <- read.csv("tax_list_BEAMING.csv") # get taxa list
names(tax_list)
tax_list <- tibble(tax_list) %>% mutate(taxon_id = X)
tax_list  <- tax_list[-1]

data <- left_join(df_fig_SS, tax_list, by = "taxon_id")
data$Species[is.na(data$Species)] = "(unclassified)"
data2 <- data %>%
  mutate(genspec = paste0(Genus, " ", Species)) %>% 
  relocate(Kingdom, Phylum, Class, Order, Family, Genus, Species, genspec) 

# add annotation
CT_W1_HIVexposure <- data2
CT_W15_HIVexposure <- data2
CT_W1_HIVexposure$tag <- "CT.W1"
CT_W15_HIVexposure$tag <- "CT.W15"

data3 <- rbind(CT_W1_HIVexposure, CT_W15_HIVexposure)
View(data3)

df <- data.frame(data3, check.names = T)
print(df)
#write_xlsx(df, "ANCOMBC.output.xlsx")
