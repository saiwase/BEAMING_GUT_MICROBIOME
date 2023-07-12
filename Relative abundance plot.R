### scripts for relative abundance plot
library(phyloseq)
library(ggplot2)
library(maditr)
library(RColorBrewer)

# Function for annotation  
# This customized function tax.lab() was created by Katie Lennard available at https://gist.github.com/kviljoen/97d36c689c5c9b9c39939c7a100720b9
#ARGUMENTS:
#otus = otu table of interest (often subset)
#physeq = phyloseq object of interest (for tax annotation)
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

### Fig3 (Beta diversity) -------
Phy.f_std <- readRDS("Phy.f_std.RDS")

barplot_stag_BEAMING_sort <- function(physeq) {
  p1 <- tax_glom(physeq, taxrank = "Species") #agglumerate at species level
  p1_30 = prune_taxa(names(sort(taxa_sums(p1), TRUE))[1:30], p1) #use top 30 most abundant taxa for plotting
  p2 <- transform_sample_counts(p1_30, function(x) x/sum(x)) #get abundance in %
  p3 <- psmelt(p2)  #create dataframe from phyloseq object
  p3$Species <- as.character(p3$Species) #convert to character
  p4 <- p3[, c("Sample", "Species","CST_pam_3", "Abundance")]
  p4 <- dcast(p3, Sample + CST_pam_3 ~ Species, value.var = "Abundance", fun.aggregate = sum)
  p4[p4$CST_pam_3 == "Cluster 1",] <- p4[p4$CST_pam_3 == "Cluster 1",][order(p4[p4$CST_pam_3 == "Cluster 1", "coli"]),]
  p4[p4$CST_pam_3 == "Cluster 2",] <- p4[p4$CST_pam_3 == "Cluster 2",][order(p4[p4$CST_pam_3 == "Cluster 2", "longum"]),]
  p4[p4$CST_pam_3 == "Cluster 3",] <- p4[p4$CST_pam_3 == "Cluster 3",][order(p4[p4$CST_pam_3 == "Cluster 3", "faecalis"]),]
  #reorder p3
  p3$Sample <- factor(p3$Sample, levels = c(p4$Sample))
  otus <- otu_table(p1_30)
  tax_table(p1_30)
  #annotation 
  labs <- tax.lab(physeq, otus = otus)
  p3$Species <- as.factor(p3$Species)
  levels(p3$Species) <- c("aerofaciens" = "C.aerofaciens", "atypica" = "V.atypica", "aureus" = "S.aureus", "bifidum" = "B.bifidum", "breve" = "B.breve", "caprae" = "S.caprae",
                          "carniphila" = "K.carniphila", "catenulatum" = "B.catenulatum", "coli" = "E.coli", "copri" = "P.copri", "dispar" = "V.dispar", "equorum" = "S.equorum",
                          "erythropolis" =  "R.erythropolis", "faecalis" = "E.faecalis", "faecium" = "E.faecium", "gallolyticus" = "S.gallolyticus", "gasseri" = "L.gasseri",
                          "haemolyticus" = "S.haemolyticus", "longum" = "B.longum", "lutetiensis" = "S.lutetiensis", "luteus" = "M.luteus", "palustris" = "K.palustris",
                          "pentosaceus" = "P.pentosaceus", "pneumoniae" = "K.pneumoniae", "putida" = "P.putida", "quasipneumoniae" = "K.quasipneumoniae", "salivarius" = "S.salivarius",
                          "saprophyticus" = "S.saprophyticus", "variicola"= "K.variicola", "vulgatus" = "B.vulgatus")
  
  
  #set color
  colourCount = length(unique(p3$Genus))
  getPalette = colorRampPalette(brewer.pal(8, "Accent"))  
  pal <- c("#7FC97F", "#CBB1C3", "#FDDA8D", "#83A3A7", "#D01487", "#8EC293", 
          "#BB5B19","#DAB6B1", "#FEE992","#FDCA89",  "#EC0877","#C9482C",
          "#9DBBA8", "#E9BA9E","#E3EA9C", "#FEF897", "#4B61AA", "#E01D5E", 
          "#ACB5BC", "#90603F","#F8BE8B", "#77479F", "#D43345", "#7B6352",
          "#BBAED1", "#A65E2C", "#B3C7A1", "#A32D93","#5380AC", "#666666")
  
  pal <- c("steelblue3", "#8EC293", "#FEE992", "azure4", "brown3", 
          "#CBB1C3", "tomato2","lightgreen", "cyan","#FDCA89",  
          "#EC0877","sienna1","lightseagreen", "lightcoral","#E3EA9C", 
          "yellow", "blue3", "#E01D5E", "chartreuse3", "chocolate4",
          "darkgoldenrod1", "#77479F", "honeydew3", "yellow4", "violet", 
          "tan3", "seagreen", "#A32D93","lightskyblue", "#666666")
  
  #plot
  barplot_species <- ggplot(data = p3, aes(x = Sample, y = Abundance, fill = Species)) 
  barplot_species <- barplot_species + geom_bar(aes(), stat = "identity", position = "stack") + 
    facet_wrap(Visit ~ Study_site, scales = "free") + guides(fill=guide_legend(nrow=3)) + 
    scale_fill_manual(values = pal) + ylab("Relative abundance") + xlab("Participants") + 
    theme(legend.position = "bottom", strip.background = element_rect(fill = "lightgrey", color = "black"), 
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.title.y = element_text(size = 17, colour = "black",face = "bold"),
          axis.title.x = element_blank(),
          strip.text = element_text(size = 17, face = "bold"), 
          axis.text.y = element_text(size = 17, face = "bold"), 
          text = element_text(size = 15, face = "bold"))
  barplot_species
  return(barplot_species)
}


# run the function
Fig2_A <- barplot_stag_BEAMING_sort(Phy.f_std)
Fig2_A

