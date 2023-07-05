
# this is for a script to publish BEAMING paper.

# I will be using Phyloseq object which are standadaized and then filtered.
# previous plots were generated using a plhyloseq pbejct with just filteing, so I need to update that part as well.
# let's call it "Phy.f_std"

# start with BEAMING_phy4 <- readRDS("BEAMING_phy4.RDS")
# BEAMING_phy4: BEAMING samples. Removed the potential contamination. PC/NCs are also removed. reads are also filtered

# fix:
# the standadization process
# tree
# figure out put has to be high enough
# Cluster 1, 2, 3
# Visit: Week 1 and Week 15

#-----------------------------------------------------------------------------------------------------
# LOAD PACKAGE ----
#-----------------------------------------------------------------------------------------------------
library(dplyr)
library(writexl)
library(rlist)
library(readxl)
library(stringr)
library(ggplot2)
library(ape)
library(vegan)
library(MMUPHin)
library(curatedMetagenomicData)
library(phyloseq)
library(ggpubr)
library(DESeq2)
library(microbiome)
library(DESeq2)
library(MMUPHin)
library(cowplot)
library(ShortRead)
library(decontam); packageVersion("decontam")

source("microbiome_custom_functions_tutorial.R") 


#-----------------------------------------------------------------------------------------------------
#  updating the meta data, prepare for the analysis ----
#-----------------------------------------------------------------------------------------------------
BEAMING_phy4 <- readRDS("BEAMING_phy4.RDS")


## Standardise (=Phy.std)
total = median(sample_sums(BEAMING_phy4))
standf = function(x, t=total) round(t * (x / sum(x)))
Phy.std = transform_sample_counts(BEAMING_phy4, standf) # standarized phyloseq object
#saveRDS(Phy.std, "Phy.std.RDS")
#Phy.std <- readRDS("Phy.std.RDS")
ntaxa(Phy.std) # 4335 - quite different number compared to the merged phyloseq (it used to be 1873 after filter)
sample_sums(Phy.std) 


## Filter the Standardised phyloseq (= Phy.f_std)
#The filter below retains only OTUs that are present at at least 10 counts at least 2% of samples OR that have a total relative abundance of at least 0.1% of the total number of reads/sample
Phy.f_std = filter_taxa(Phy.std, function(x) sum(x > 10) > (0.2*length(x)) | sum(x) > 0.001*total, TRUE) 
ntaxa(Phy.f_std) #  2084 
sample_sums(Phy.f_std) 
#saveRDS(Phy.f_std, "Phy.f_std.RDS")
#Phy.f_std <- readRDS("Phy.f_std.RDS")

# > Phy.f_std_std
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2084 taxa and 442 samples ]
# sample_data() Sample Data:       [ 442 samples by 82 sample variables ]
# tax_table()   Taxonomy Table:    [ 2084 taxa by 7 taxonomic ranks ]
# refseq()      DNAStringSet:      [ 2084 reference sequences ]



### checking the metadata
Phy.f_std <- readRDS("Phy.f_std.RDS")
new_meta <- tibble(data.frame(sample_data(Phy.f_std))) %>% View()


Phy.f <- readRDS("Phy.f.RDS")
old_meta <- tibble(data.frame(sample_data(Phy.f)))
names(old_meta)
old_meta_extract <- old_meta %>% select(PID, Visit, Sample_ID)


setdiff(new_meta$Sample_ID, old_meta$Sample_ID) # character(0)
length(new_meta$PID)
intersect(new_meta$Sample_ID, old_meta$Sample_ID)# 442


setdiff(new_meta$wflz, old_meta$wflz)
setdiff(new_meta$ChosenFeedingMeth, old_meta$ChosenFeedingMeth)
setdiff(new_meta$TAG3, old_meta$TAG3)


# update the variables all together
sample_data(Phy.f_std) <- sample_data(Phy.f)




levels(sample_data(Phy.f_std)$CST_pam_3) <- c("Cluster1"="Cluster 1","Cluster2"="Cluster 2", "Cluster3"="Cluster 3")
levels(sample_data(Phy.f_std)$CST_pam_3_bray) <- c("Cluster1"="Cluster 1","Cluster2"="Cluster 2", "Cluster3"="Cluster 3")


#-----------------------------------------------------------------------------------------------------
#  BEAMING paper figures ----
#-----------------------------------------------------------------------------------------------------
### Fig 1-A (HEATMAP) ----
### Fig 1-A 
## this is the HEATMAP function (w/o using tax_glom.kv) - for Fig 1 (A) - refer code below instead!

# remove the unnecessary annotation (i.e. feeding status)
library(RColorBrewer)
library(pheatmap)
heatmap_stand_otu <- function(physeq) {
  #calculate weighted unifrac for clustering
  set.seed(2)
  glom_random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
  glom_tree = merge_phyloseq(physeq,glom_random_tree)
  diss <- phyloseq::distance(glom_tree,method = "bray", type = "samples")
  #use top 30 most abundant taxa for plotting
  glom_30 = prune_taxa(names(sort(taxa_sums(physeq), TRUE))[1:30], physeq) 
  otus <- otu_table(glom_30)
  #log2 transform otus for better heatmap visualisation colour scale
  zs <- which(otus==0)#find zero entries 
  if(length(zs)>0){#if there are zero values, transform
    otus[zs] <- 0.1 #set zs to 0.1
    otus <- log2(otus)
  } else{otus <- log2(otus)}
  plot <- otus #data to be plotted
  #get metadata for annotation 
  meta_hm <- data.frame(row.names=sample_names(physeq),sample_data(physeq)$CST_pam_3_bray,sample_data(physeq)$Status2,sample_data(physeq)$Study_site)
  names(meta_hm) <- c("Community cluster", "HIV exposure status","Study site")
  #annotation 
  labs <- tax.lab(physeq,otus=otus)
  #color
  col.pal <- brewer.pal(9,"YlOrRd")
  annot.color.col <- list("Community cluster"=c("Cluster 1"="#FDCA89","Cluster 2"="#F1948A","Cluster 3"="#E06666"),"HIV exposure status" =c("iHEU" ="maroon1", "iHUU"="#F0E161"),"Study site" = c("South Africa" = "red", "Nigeria" = "blue"))
  #heatmap
  hm <- pheatmap(plot,annotation = meta_hm, color = col.pal, clustering_method = "complete", clustering_distance_cols = diss, show_colnames = FALSE, labels_row = labs, annotation_colors = annot.color.col, fontsize = 11)
  return(hm)
} 


# HEATMAP (2') - using "Phy.f_std"(all sites / Only birth)
Phy.f_std <- readRDS("Phy.f_std.RDS")
Phy.f_std_birth <- subset_samples(Phy.f_std, Visit =="Week 1") # 204 samples
levels(sample_data(Phy.f_std_birth)$"Study_site") <- c("CT"="South Africa","Nigeria"="Nigeria") # need to change
levels(sample_data(Phy.f_std_birth)$Status2) <- c("HEU"="iHEU","HU"="iHUU") 


Fig1_A <- heatmap_stand_otu(Phy.f_std_birth)
Fig1_A  # using PAM clustering.

#ggsave(file = "./figure_output/fig1_heatmap_birth.pdf", plot = Fig1_A, dpi = 800, width = 9, height = 10) # in this way, I can control the size of the figure
dev.off()



### Fig 1-A (HEATMAP) ---- top 20 instead (*)
### Fig 1-A 
## this is the HEATMAP function (w/o using tax_glom.kv) - for Fig 1 (A)

# remove the unnecessary annotation (i.e. feeding status)
library(RColorBrewer)
library(pheatmap)
heatmap_stand_otu <- function(physeq) {
  #calculate weighted unifrac for clustering
  set.seed(2)
  glom_random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
  glom_tree = merge_phyloseq(physeq,glom_random_tree)
  diss <- phyloseq::distance(glom_tree,method = "bray", type = "samples")
  #use top 20 most abundant taxa for plotting
  glom_20 = prune_taxa(names(sort(taxa_sums(physeq), TRUE))[1:18], physeq) # 18 instead of 20 as I merge two bacteria
  otus <- otu_table(glom_20)
  #log2 transform otus for better heatmap visualisation colour scale
  zs <- which(otus==0)#find zero entries 
  if(length(zs)>0){#if there are zero values, transform
    otus[zs] <- 0.1 #set zs to 0.1
    otus <- log2(otus)
  } else{otus <- log2(otus)}
  plot <- otus #data to be plotted
  #get metadata for annotation 
  meta_hm <- data.frame(row.names=sample_names(physeq),sample_data(physeq)$CST_pam_3_bray,sample_data(physeq)$Status2,sample_data(physeq)$Study_site)
  names(meta_hm) <- c("Community cluster", "HIV exposure status","Study site")
  #annotation 
  labs <- tax.lab(physeq,otus=otus)
  labs[labs == "B.longum"] <- "B.longum (subsp. infantis)" # ASV1
  labs[labs == "B.longum.1"] <- "B.longum (subsp. longum)" # ASV4
  #color
  col.pal <- brewer.pal(9,"YlOrRd")
  annot.color.col <- list("Community cluster"=c("Cluster 1"="#8E999F","Cluster 2"="#F1948A","Cluster 3"="#C5E1A5"),"HIV exposure status" =c("iHEU" ="#B05A7A", "iHUU"="#F99417"),"Study site" = c("South Africa" = "#E87D72", "Nigeria" = "#54BCC2"))
  #heatmap
  hm <- pheatmap(plot,annotation = meta_hm, color = col.pal, clustering_method = "complete", clustering_distance_cols = diss, show_colnames = FALSE, labels_row = labs, annotation_colors = annot.color.col, fontsize = 13)
  return(hm)
} 


# HEATMAP (2') - using "Phy.f_std"(all sites / Only birth)
Phy.f_std <- readRDS("Phy.f_std.RDS")
Phy.f_std_birth <- subset_samples(Phy.f_std, Visit =="Week 1") # 204 samples
levels(sample_data(Phy.f_std_birth)$"Study_site") <- c("CT"="South Africa","Nigeria"="Nigeria") # need to change
levels(sample_data(Phy.f_std_birth)$Status2) <- c("HEU"="iHEU","HU"="iHUU") 


# fine tuning - merge "Staphylococcus saprophyticus" and  Staphylococcus haemolyticus so that the annotation would be better
Phy.f_std_birth_merged.phy1 = merge_taxa(Phy.f_std_birth, c("ASV12","ASV22")) #Staphylococcus saprophyticus
Phy.f_std_birth_merged.phy2 = merge_taxa(Phy.f_std_birth_merged.phy1, c("ASV18","ASV68")) #Staphylococcus haemolyticus

Fig1_A <- heatmap_stand_otu(Phy.f_std_birth_merged.phy2)
Fig1_A  # using PAM clustering.

ggsave(file = "./figure_output/fig1_heatmap_birth.top20.bray.v3.pdf", plot = Fig1_A, dpi = 800, width = 9, height = 10) # in this way, I can control the size of the figure
dev.off()


### for better annotation
glom_20 = prune_taxa(names(sort(taxa_sums(Phy.f_std_birth), TRUE))[1:20], Phy.f_std_birth) 
otus <- otu_table(glom_20)
rownames(tax_table(glom_20)) # these are the top 20 taxa (for birth data)
# [1] "ASV1"  "ASV2"  "ASV3"  "ASV4"  "ASV5"  "ASV6"  "ASV7"  "ASV8"  "ASV9"  "ASV10" "ASV11" "ASV12" "ASV13" "ASV18" "ASV22"
# [16] "ASV24" "ASV30" "ASV31" "ASV37" "ASV68"

#
dput(rownames(tax_table(glom_20))) # output as vector
my.ASVlist <- c("ASV1", "ASV2", "ASV3", "ASV4", "ASV5", "ASV6", "ASV7", "ASV8", 
                "ASV9", "ASV10", "ASV11", "ASV12", "ASV13", "ASV18", "ASV22", 
                "ASV24", "ASV30", "ASV31", "ASV37", "ASV68")

## Get ASV245 info
top.W1.ASV <- subset_taxa(Phy.f_std, rownames(tax_table(Phy.f_std)) %in% my.ASVlist)
top.W1.ASV_seq <- refseq(top.W1.ASV)
ASV <- as.data.frame(top.W1.ASV_seq)
ASV$ASV <- rownames(ASV)
ASV <-  tibble(ASV)

# get annotation
tax_list <- read.csv("/Users/saoriiwase/Desktop/tax_list_BEAMING.csv") # obtain the tax_list (for full tax info)
names(tax_list)
tax_list <- tibble(tax_list) %>% mutate(ASV = X)

ASV.v2 <- left_join(ASV, tax_list, by = "ASV") # merge
ASV.v3 <- ASV.v2 %>% dplyr::rename("Sequence" =  "x") %>% select(-X)
# View(ASV.v3)


df <- data.frame(ASV.v3, check.names = T)
print(df)
write_xlsx(df,"./figure_output/fig1_W1.top20.ASV.Seq.list.xlsx")






### Fig 1-B (alpha diversity)----
## alpha diversity - birth only and compare by site - for Fig 1 (B). Change the plot theme.
Phy.f_std <- readRDS("Phy.f_std.RDS")
Phy.f_std_birth <- subset_samples(Phy.f_std, Visit =="Week 1")
levels(sample_data(Phy.f_std_birth)$"Study_site") <- c("CT"="South Africa","Nigeria"="Nigeria") # need to change the name
p1 <- plot_richness(Phy.f_std_birth, x="Study_site", measures="Shannon",color = "Study_site",
                   # title = "Alpha Diversity at birth"
                   )+xlab("Study site")+
  geom_boxplot(alpha=0.6)+ theme_bw()+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
        strip.text = element_blank(), 
        axis.text.y=element_text(size=16, face="bold"),
        text = element_text(size=17, face="bold"))+
  labs(y= "Shannon Index")

# and add statistics
# https://rpubs.com/lconteville/713954
a_my_comparisons <- list( c("South Africa", "Nigeria"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

Fig1_B <- p1 + stat_compare_means(method = "wilcox.test", paired = F,
                                comparisons = a_my_comparisons, 
                                label = "p.signif", symnum.args = symnum.args) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
Fig1_B  #significantly different - CT increases the alpha diversity over the time whereas decased in Nigerian infants

Fig1_B

ggsave(file = "./figure_output/fig1_alpha_diversity_birth_SAvsNigeria.color.pdf", plot = Fig1_B, dpi = 500, width = 7, height = 7) # in this way, I can control the size of the figure
dev.off()


### alpha diversity, adjusting for mother's maternal weight, age and gravity 
# https://bookdown.org/rwnahhas/RMPH/mlr-distinctions.html # introduction to regression models
str(sample_data(Phy.f_std_birth)$Study_site)
str(sample_data(Phy.f_std_birth)$Mum_weight)
str(sample_data(Phy.f_std_birth)$Mum_Age)
str(sample_data(Phy.f_std_birth)$Num_Preg)
str(sample_data(Phy.f_std_birth)$Shannon)

Df <- data.frame(sample_data(Phy.f_std_birth))

# SequenceLot/ Education/ Occupation.v2/ Houses/ Refrigerator/ Running_water/ Marital_status/ Num_Preg/ Mum_weight/ Ges_age/ DeliverVag
Df2 <- Df %>% select(Shannon, Study_site,SequenceLot, Education, Occupation.v2, Houses, Refrigerator,
                     Running_water, Marital_status, Num_Preg, Mum_weight, Mum_Age,Ges_age, DeliverVag) %>% tibble()%>% na.omit()

Df2 <- Df %>% select(Shannon,Mum_Age,Mum_weight,Num_Preg,Study_site,SequenceLot) %>% tibble()%>% na.omit()

# run Levene’s test.
# if the test came back insignificantly, implies that the variance is homogenous, and we can proceed with our ANOVA. 
leveneTest(Shannon~Study_site,Df2) # not significant -  homogeneity of variance (i.e., the difference between the variances is zero)
library(psych)
describeBy(Df2$Shannon, Df2$Study_site) # we can summarize the skew too!!!

ancova_model1 <- aov(Shannon ~ Study_site+Mum_weight+Mum_Age+Num_Preg, data = Df2)
summary(ancova_model1)
car::Anova(ancova_model1, type=3)
summary.lm(ancova_model1) #if you were to look at the fit2 object through the summary.lm command, which produces the output in the style of a linear model (i.e., OLS) and also uses type III errors, you would get the same correct information in the output as via the Anova command from the car package.


ancova_model2 <- aov(Shannon ~ Study_site+SequenceLot, data = Df2)
ancova_model2 <- aov(Shannon ~ Study_site+Mum_weight+Mum_Age+Num_Preg+SequenceLot, data = Df2)
car::Anova(ancova_model2, type=3)
summary(ancova_model2)
summary.lm(ancova_model2)   # Study_site was still significantly associated with Shannon diversity, after controllling for mum's weight, age and gravidity.
# Study site: P< 2.2e-16 *** Sum sq 52.845/ F-val 128.6807



ancova_model3 <- aov(Shannon ~ Study_site+SequenceLot+Education+ Occupation.v2+ Houses+ Refrigerator+ Running_water+ Marital_status+ Num_Preg+ Mum_weight+ Ges_age+ DeliverVag, data = Df2)
summary(ancova_model3)
SequenceLot/ Education/ Occupation.v2/ Houses/ Refrigerator/ Running_water/ Marital_status/ Num_Preg/ Mum_weight/ Ges_age/ DeliverVag


### Fig 1-C (Beta-diversity)----
### Beta-diversity (1') PCoA plot comparing by site for each time point (Birth samples only) - jsd Fig 1 (c)
# for PcOA plot using Birth data
Phy.f_std <- readRDS("Phy.f_std.RDS")
Phy.f_std_birth <- subset_samples(Phy.f_std, Visit=="Week 1") 
levels(sample_data(Phy.f_std_birth)$"Study_site") <- c("CT"="South Africa","Nigeria"="Nigeria")
ord.BC <- ordinate(Phy.f_std_birth, method = "PCoA", distance = "bray", k=3, trymax=1000) # try PcOA as NMDS does not work
color = c("Study_site")
shape = c("CST_pam_3_bray")
#title= c("PCoA of 16S microbiome at birth (Jensen-Shannon distance; k=2)")  
MDS = plot_ordination(Phy.f_std_birth, ord.BC, color = color, shape = shape#, title = title
                      )
Fig1_C  = MDS+theme_bw() +theme(axis.text=element_text(size=17, face="bold"),
                               axis.title=element_text(size=17,face="bold"), 
                               legend.title=element_text(size=17 ,face="bold"),
                               legend.text = element_text(size =16))+
  labs(color=color)+geom_point(size=3)+
  guides(color = guide_legend(title = "Study site"))+  # change the ledgend labee title
  guides(shape = guide_legend(title = "Community cluster")) # change the shape title

Fig1_C 

ggsave(file = "./figure_output/fig1_birth_PCoA_k2.resize.bray.pdf", plot = Fig1_C, dpi = 800, width = 9, height = 6) # in this way, I can control the size of the figure
dev.off()



#### PERMANOVA for beta-diversity W1 samples comapred by study site
names(sample_data(Phy.f_std_birth))

Phy.f_std_birth <- subset_samples(Phy.f_std, Visit =="Week 1")

diss <- phyloseq::distance(Phy.f_std_birth, "bray", parallel=TRUE)  # jsd = Jensen–Shannon divergence / bray = Bray Curtis 
Phy.f_std_birth_df <- data.frame(row.names=sample_names(Phy.f_std_birth),sample_data(Phy.f_std_birth)) #create dataframe with metadata

# Whether "Study_site" differ significantly from each other
set.seed(2)
adonis2(diss ~ Study_site, data=Phy.f_std_birth_df, permutations = 999, by = "terms", na.action=na.exclude)  #0.001

# Now, adjust Beta-diversity for maternal weight, age and gravity 
# https://forum.qiime2.org/t/how-to-properly-design-formulas-in-adonis/17996
#https://statisticsbyjim.com/regression/interaction-effects/
str(sample_data(Phy.f_std_birth)$Study_site)
str(sample_data(Phy.f_std_birth)$Mum_weight)

str(sample_data(Phy.f_std_birth)$Num_Preg)


str(sample_data(Phy.f_std_birth)$Mum_Age)

sample_data(Phy.f_std_birth)$Education <- as.factor(sample_data(Phy.f_std_birth)$Education)
sample_data(Phy.f_std_birth)$Occupation <-  as.factor(sample_data(Phy.f_std_birth)$Occupation) # 6 is unemployment

sample_data(Phy.f_std_birth)$Occupation.v2 <- forcats::fct_recode(sample_data(Phy.f_std_birth)$Occupation,
                                                                  "0" = "0",
                                                                  "0" = "1",
                                                                  "0" = "2",
                                                                  "0" = "3", 
                                                                  "0" = "4",
                                                                  "0" = "5", 
                                                                  "1" = "6", 
                                                                  "0" = "7", 
                                                                  "0" = "8") # since the latest demographic table shows Unemployed or not, I made it binary

sample_data(Phy.f_std_birth)$Houses #0 or 1

sample_data(Phy.f_std_birth)$Refrigerator <- as.factor(sample_data(Phy.f_std_birth)$Refrigerator)
sample_data(Phy.f_std_birth)$Running_water <- as.factor(sample_data(Phy.f_std_birth)$Running_water)
sample_data(Phy.f_std_birth)$Marital_status <- as.factor(sample_data(Phy.f_std_birth)$Marital_status)

# 1=Married; 2=Living together; 3=Widowed; 4=Divorced; 5=Separated; 6=Single; 7=Other
sample_data(Phy.f_std_birth)$Marital_status.v2 <-  forcats::fct_recode(sample_data(Phy.f_std_birth)$Marital_status,
                                                          "1" = "1",
                                                          "1" = "2", 
                                                          "0" = "3",
                                                          "0" ="5", 
                                                          "0" = "6")

sample_data(Phy.f_std_birth)$DeliverVag
sample_data(Phy.f_std_birth)$Ges_age_cat
sample_data(Phy.f_std_birth)$Ges_age

names(sample_data(Phy.f_std_birth))

set.seed(2)
# adjust for SequenceLot/ Education/ Occupation.v2/ Houses/ Refrigerator/ Running_water/ Marital_status/ Num_Preg/ Mum_weight/ Mum_Age/ Ges_age/ DeliverVag
adonis2(diss ~ Study_site +SequenceLot+ Education+ Occupation.v2+ Houses+ Refrigerator+ Running_water+ Marital_status+ Num_Preg+ Mum_weight+ Mum_Age +Ges_age+ DeliverVag , data=Phy.f_std_birth_df, permutations = 999, by = "terms", na.action=na.exclude)


adonis2(diss ~ Study_site + Mum_weight + Mum_Age + Num_Preg+ SequenceLot, data=Phy.f_std_birth_df, permutations = 999, by = "terms", na.action=na.exclude)
adonis2(diss ~ Study_site / Mum_weight / Mum_Age / Num_Preg, data=Phy.f_std_birth_df, permutations = 999, by = "terms", na.action=na.exclude)

adonis2(diss ~ Study_site * Mum_weight * Mum_Age * Num_Preg, data=Phy.f_std_birth_df, permutations = 999, by = "terms", na.action=na.exclude) #R2: 0.08557; p: 0.001
adonis2(diss ~ Study_site / Mum_weight / Mum_Age / Num_Preg, strata = Phy.f_std_birth_df$Study_site, data=Phy.f_std_birth_df, permutations = 999, by = "terms", na.action=na.exclude)



# https://stats.stackexchange.com/questions/188519/adonis-in-vegan-order-of-variables-or-use-of-strata/238962
# https://stats.stackexchange.com/questions/312302/adonis-in-vegan-order-of-variables-non-nested-with-one-degree-of-freedom-for


### Fig1-D (DESeq2) -----

## DESeq2 (4) - Site comparison Birth samples (CT vs Nigeria) - sgnificance 0.01 - fig 1 (D)
Phy.f_std_birth <- subset_samples(Phy.f_std, Visit=="Birth")
sample_data(Phy.f_std_birth)$Study_site  # Levels: CT Nigeria


sample_data(Phy.f_std_birth)$Study_site <- relevel(sample_data(Phy.f_std_birth)$Study_site, ref="South Africa")

#Convert the phyloseq object to a DESeqDataSet and run DESeq2:
ds = phyloseq_to_deseq2(Phy.f_std_birth, ~ Study_site)
#ds$Study_site <- relevel(ds$Study_site, ref = "South.Africa")
ds = DESeq(ds, fitType = "local")  # get an error...rror in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : every gene contains at least one zero, cannot compute log geometric means
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(ds), 1, gm_mean)
ds = estimateSizeFactors(ds, geoMeans = geoMeans)
ds = DESeq(ds, fitType="local")
res = results(ds)
res = res[order(res$padj, na.last = NA), ]
alpha = 0.01 #this is your significance level, adjust if desired
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Phy.f_std_birth)[rownames(sigtab), ], "matrix"))
head(sigtab)
sigtab %>% View()

resultsNames(ds) 
rownames(sigtab) 

#only show bugs annotated at the Species level
#sigtabspec = subset(sigtab, !is.na(Species))
sigtabgen = subset(sigtab, !is.na(Genus))
#sigtabord = subset(sigtab, !is.na(Order))
head(sigtabgen)

#write.csv(sigtabgen, file = "./figure_output/Sfig-sigtabgen_birth_p0.01_SAvsNigeria.csv")


#sort for pretty plotting
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

#Alternative: Phylum level, for plotting purposes - don't do both
#x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
#x = sort(x, TRUE)
#sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

#fix taxrank and subset to highly significant asvs
sigtabgen.hs <- sigtabgen[which(sigtabgen$padj < 0.005),]

#subset to greater magnitude of FC
sigtabgen$Species[is.na(sigtabgen$Species)] = "(unclassified)"
sigtabgen$Species[sigtabgen$Species == "NA"] <- "unknown"
sigtabgen$genspec <- paste(sigtabgen$Genus, sigtabgen$Species, sep = " ")
sigtab.merged.hs <- sigtabgen[which(sigtabgen$padj < 0.005),]

sigtabgen %>%
  subset(log2FoldChange > 0.5) -> sigtabgen.hi
sigtabgen %>%
  subset(log2FoldChange < -0.5) -> sigtabgen.low

sigtab.merged.hs %>%
  subset(log2FoldChange > 0.5) -> sigtabgen.hi

sigtab.merged.hs %>%
  subset(log2FoldChange < -0.5) -> sigtabgen.low

sigtabgen.fc <- rbind.data.frame(sigtabgen.hi, sigtabgen.low)                          
View(sigtabgen.fc)

#plot the results (Genus currently set, and colored by Order, set as desired)
Fig1_D <- ggplot(sigtabgen.fc, aes(y=Genus, x=log2FoldChange, fill=genspec)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_vline(xintercept = 0.5, color = "indianred3", size = 0.5, lty=2) +
  geom_vline(xintercept = -0.5, color = "indianred3", size = 0.5, lty=2) +
  geom_point(pch = 21, size=3)+ theme_bw()+
  theme(axis.text.y = element_text(size = 14, colour = "black",face = "bold"),
        axis.title.y = element_text(size = 17, colour = "black",face = "bold"),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 17, colour = "black",face = "bold"),
        strip.background =element_rect(fill="white"),
        strip.text.x = element_text(colour = "mediumpurple4", size = 15, hjust = 0.5, vjust = 0.5, face = "bold"),
        legend.text = element_text(size = 14, colour = "black"),
        legend.title = element_text(size=15)) +  theme(legend.position="right", legend.title=element_blank())+
  guides(fill = guide_legend(nrow = 12)) + guides(fill=guide_legend(ncol = 2))+
  labs(x=expression(log[2]~fold~change))  #scale_fill_manual(values=custom_colors_36) + #+ ggtitle("DESeq2 SA vs Nigeria at birth, p=0.01, w/o melting the lowest taxonomy")

Fig1_D # something wrong with DEseq - it looks like opposite

ggsave(file = "./figure_output/fig1_DESeq2_birth_SAvsNigeria.pdf", plot = Fig1_D, dpi = 900, width = 15, height = 7) # in this way, I can control the size of the figure
dev.off()









# Combine for Figure1
Fig1_A
Fig1_B
Fig1_C
Fig1_D

library(cowplot)
cowplot::plot_grid(Fig1_B, Fig1_C, Fig1_D, nrow=2)

### Figure 2 ###
### Fig 2-A (relative abundance)-----
### Simiar to the figure above but with re-ordering the barplot - fig 2 (A)
##BARPLOT - ORDERED AND CLUSTERED
library(maditr)
Phy.f_std <- readRDS("Phy.f_std.RDS")
class(Phy.f_std@sam_data$CST_pam_3)
View(Phy.f_std@tax_table)

barplot_stag_BEAMING_sort <- function(physeq) {
  p1 <- tax_glom(physeq, taxrank = 'Species') #agglumerate at species level
  p1_30 = prune_taxa(names(sort(taxa_sums(p1), TRUE))[1:30], p1) #use top 30 most abundant taxa for plotting
  p2 <- transform_sample_counts(p1_30, function(x) x/sum(x)) #get abundance in %
  p3 <- psmelt(p2)  #create dataframe from phyloseq object
  p3$Species <- as.character(p3$Species) #convert to character
  #sort based on abundance of dominant bacteria in each cluster
  p4 <- p3[,c("Sample", "Species","CST_pam_3", "Abundance")]
  p4 <- dcast(p3, Sample + CST_pam_3 ~ Species, value.var = "Abundance", fun.aggregate = sum)
  p4[p4$CST_pam_3 =="Cluster 1",] <- p4[p4$CST_pam_3=="Cluster 1",][order(p4[p4$CST_pam_3=="Cluster 1","coli"]),]
  p4[p4$CST_pam_3 =="Cluster 2",] <- p4[p4$CST_pam_3=="Cluster 2",][order(p4[p4$CST_pam_3=="Cluster 2","longum"]),]
  p4[p4$CST_pam_3 =="Cluster 3",] <- p4[p4$CST_pam_3=="Cluster 3",][order(p4[p4$CST_pam_3=="Cluster 3","faecalis"]),]
  #reorder p3
  p3$Sample <- factor(p3$Sample, levels=c(p4$Sample))
  otus <- otu_table(p1_30)
  tax_table(p1_30)
  #annotation 
  labs <- tax.lab(physeq,otus=otus)
  
  p3$Species <- as.factor(p3$Species)# this is to manually editt the annotation - very figure specific and order needs to be exactly the same.
  levels(p3$Species) <- c("aerofaciens"="C.aerofaciens","atypica"= "V.atypica", "aureus" = "S.aureus","bifidum" = "B.bifidum","breve"= "B.breve","caprae" = "S.caprae",
                          "carniphila" = "K.carniphila","catenulatum"= "B.catenulatum", "coli"= "E.coli", "copri"= "P.copri","dispar" = "V.dispar","equorum" = "S.equorum",
                          "erythropolis" =  "R.erythropolis","faecalis" = "E.faecalis","faecium" = "E.faecium","gallolyticus" = "S.gallolyticus", "gasseri" = "L.gasseri",
                          "haemolyticus" ="S.haemolyticus", "longum"="B.longum","lutetiensis"= "S.lutetiensis", "luteus"= "M.luteus","palustris" ="K.palustris",
                          "pentosaceus" ="P.pentosaceus", "pneumoniae" = "K.pneumoniae","putida"= "P.putida", "quasipneumoniae" ="K.quasipneumoniae", "salivarius" = "S.salivarius",
                          "saprophyticus" = "S.saprophyticus",  "variicola"=  "K.variicola","vulgatus" = "B.vulgatus")
  
  
  #set color palette to accommodate the number of species
  colourCount = length(unique(p3$Genus))
  getPalette = colorRampPalette(brewer.pal(8, "Accent"))  
  pal<- c("#7FC97F", "#CBB1C3", "#FDDA8D", "#83A3A7", "#D01487", "#8EC293", 
          "#BB5B19","#DAB6B1", "#FEE992","#FDCA89",  "#EC0877","#C9482C",
          "#9DBBA8", "#E9BA9E","#E3EA9C", "#FEF897", "#4B61AA", "#E01D5E", 
          "#ACB5BC", "#90603F","#F8BE8B", "#77479F", "#D43345", "#7B6352",
          "#BBAED1", "#A65E2C", "#B3C7A1", "#A32D93","#5380AC", "#666666")
  
  pal<- c("steelblue3", "#8EC293", "#FEE992", "azure4", "brown3", 
          "#CBB1C3", "tomato2","lightgreen", "cyan","#FDCA89",  
          "#EC0877","sienna1","lightseagreen", "lightcoral","#E3EA9C", 
          "yellow", "blue3", "#E01D5E", "chartreuse3", "chocolate4",
          "darkgoldenrod1", "#77479F", "honeydew3", "yellow4", "violet", 
          "tan3", "seagreen", "#A32D93","lightskyblue", "#666666")
  #annotation 
  #plot
  barplot_species <- ggplot(data=p3, aes(x=Sample, y=Abundance, fill=Species)) 
  barplot_species <- barplot_species + geom_bar(aes(), stat="identity", position="stack") + 
    facet_wrap(Visit ~ Study_site, scales = "free") + guides(fill=guide_legend(nrow=3)) + 
    scale_fill_manual(values=pal) + ylab("Relative abundance") +xlab("Participants") + 
    #ggtitle("Bacterial composition at each time pooint by study site") + 
    theme(legend.position="bottom", strip.background = element_rect(fill="lightgrey", color = "black"), 
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          panel.background=element_blank(), panel.border=element_blank(),
          axis.title.y = element_text(size = 17, colour = "black",face = "bold"),
          axis.title.x = element_blank(),
          strip.text = element_text(size=17, face="bold"), 
          axis.text.y=element_text(size=17, face="bold"), 
          text = element_text(size=15, face="bold"))
  barplot_species
  return(barplot_species)
}


levels(sample_data(Phy.f_std)$Visit) <- c("Birth"="Week 1","W15"="Week 15")
Fig2_A <- barplot_stag_BEAMING_sort(Phy.f_std)
Fig2_A

ggsave(file = "./figure_output/Fig2_barplot_birth_W15.pdf", plot = Fig2_A, dpi = 900, width = 20, height = 12) # in this way, I can control the size of the figure
dev.off()



### Fig 2-B (alpha diversity) ------
## alpha diversity by Visit - separate by TP and facet by Study site  - fig 2 (B) version 2
# join the PID by geom_line

p1 <- plot_richness(Phy.f_std, x="Visit", measures="Shannon",color = "Study_site")+
  geom_boxplot(alpha=0.6)+ geom_line(aes(group=PID2),size = 0.1) +
  facet_wrap(~Study_site)+theme_bw()+
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
        strip.text = element_text(size=17, face="bold"), axis.text.y=element_text(size=16, face="bold"), text = element_text(size=17, face="bold"))+ labs(y= "Shannon Index")

# and add statistics
# https://rpubs.com/lconteville/713954
a_my_comparisons <- list( c("Week 1", "Week 15"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
Fig2_B <- p1 + stat_compare_means(method = "wilcox.test", 
                                comparisons = a_my_comparisons, 
                                label = "p.signif", symnum.args = symnum.args) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

Fig2_B  

ggsave(file = "./figure_output/fig2_alpha_diversity_transition_v2.withLine.color.pdf", plot = Fig2_B, dpi = 500, width = 7, height = 7) # in this way, I can control the size of the figure
dev.off()




### Fig 2-C (Beta-diversity) -----
### Beta-diversity (3') PCoA plot comparing by site (include both birth and W15) - no annotation of hiv exposure buy annotate with TP  - fig 2 (C)

set.seed(20) #select any start seed, in this case 2 (ensures reproducibility)
ord.BC <- ordinate(Phy.f_std, method = "PCoA", distance = "bray", k=2, trymax=1000) # try PcOA as NMDS does not work
ord.BC

color = c("Study_site")
shape = c("Visit")
#title= c("PCoA of 16S microbiome (all visit), Bray, k=2")  
MDS = plot_ordination(Phy.f_std, ord.BC, color = color, shape=shape)
Fig2_C  = MDS +theme_bw() +theme(axis.text=element_text(size=17, face="bold"),
                                 axis.title=element_text(size=17,face="bold"), 
                                 legend.title=element_text(size=17 ,face="bold"),
                                 legend.text = element_text(size =16))+ 
  labs(color=color, shape=shape)+geom_point(size=2)+ 
  guides(color = guide_legend(title = "Study site")) 

Fig2_C

ggsave(file = "./figure_output/fig2_birth.W15_PCoA_k2.resize.bray.pdf", plot = Fig2_C, dpi = 500, width = 8.5, height = 5) # in this way, I can control the size of the figure
dev.off()


# check with the adonis2
diss <- phyloseq::distance(Phy.f_std, "bray", parallel=TRUE)   # bray = bray
Phy.f_std_df <- data.frame(row.names=sample_names(Phy.f_std),sample_data(Phy.f_std)) 
set.seed(2)

# Whether "Study_site" differ significantly from each other
adonis2(diss ~ Study_site, data=Phy.f_std_df, permutations = 999, by = "terms", na.action=na.exclude)  
adonis2(diss ~ Visit, data=Phy.f_std_df, permutations = 999, by = "terms", na.action=na.exclude)  



### Fig 2-D (Alluvial Plot) -----
library(msm)
Phy.f_std <- readRDS("Phy.f_std.RDS")
# sample_data(Phy.f_std)$DeliverVag = as.factor(sample_data(Phy.f_std)$DeliverVag)


levels(sample_data(Phy.f_std)$Visit) 
levels(sample_data(Phy.f_std)$CST_pam_3_bray) 
# sample_data(Phy.f_std)$Visit <- relevel(sample_data(Phy.f_std)$Visit, ref = "Birth")

sample_data(Phy.f_std)

df<- sample_data(Phy.f_std)
PID.with.2TP <- df$PID2[duplicated(df$PID2)]
physeq <- subset_samples(Phy.f_std, PID2 %in% PID.with.2TP) # only sample with 2TP

df = sample_data(physeq)[,c("Visit","CST_pam_3_bray","PID", "Status2", "Study_site")]
df$CST_pam_3_bray = as.character(unlist(df$CST_pam_3_bray))

df$CST_pam_3_bray[df$CST_pam_3_bray=="Cluster 1"] <- "1" #Change like this because the states apparently are required to be 1, 2, 3 as oppoased to C1, C2, C3 (not sure if they're assumed to be rank ordered, i.e. 1 < 2 <3)
df$CST_pam_3_bray[df$CST_pam_3_bray=="Cluster 2"] <- "2"
df$CST_pam_3_bray[df$CST_pam_3_bray=="Cluster 3"] <- "3"
df$CST_pam_3_bray = as.factor(df$CST_pam_3_bray)



statetable.msm(CST_pam_3_bray, PID, data=df) # transition of the cluster from one to another
## to
## from  1  2  3
## 1 34 23  1
## 2  2 17  0
## 3  9 72  5


library(ggalluvial)
## Alluvial Plots (1) - transition of PAM clustering -- fig 2 (D)
head(df)
df$CST_pam_3_bray <- forcats::fct_recode(df$CST_pam_3_bray, "Cluster 1 " = "1", "Cluster 2 " = "2","Cluster 3 " = "3")

Fig2_D <- ggplot(as.data.frame(df),
             aes(x=Visit, stratum=CST_pam_3_bray, alluvium=PID, fill=CST_pam_3_bray)) +
  geom_alluvium(aes(fill = CST_pam_3_bray), width = 1/16) +
  geom_stratum(width = 1/3, alpha = .8) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Week 1", "Week 15"), expand = c(.05, .05,.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  #ggtitle("Transition of PAM clustering (#1-#3) over time") +
  theme_bw() + 
  geom_text(stat = "stratum", aes(label = after_stat(stratum)))+
  facet_wrap(~Study_site) +
  theme(axis.text.y=element_text(size=17, face="bold"), 
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"),
        axis.title=element_text(size=18,face="bold"), 
        legend.title=element_text(size=17, face="bold"),
        legend.text = element_text(size =17),
        strip.text = element_text(size=18, face="bold"))+
  guides(fill=guide_legend(title="Community cluster"))+
  labs(y= "Number of infants")+scale_fill_manual(values=c("#8E999F","#F1948A","#C5E1A5"))

ggsave(file = "./figure_output/fig2_AlluvialPlots_PAMclustering.bray.v3.pdf", plot = Fig2_D, dpi = 800, width = 10, height = 10) 
dev.off()






###### Figure 3 (A) (alpha diversity) -----
## alpha diversity by Status - separate by TP and facet by Study site & TP - Figure 3 (A)
levels(sample_data(Phy.f_std)$"Study_site") <- c("CT"="South Africa","Nigeria"="Nigeria") # need to change the name
levels(sample_data(Phy.f_std)$Status2) <- c("HEU"="iHEU","HU"="iHUU") # need to change the name

p1 <- plot_richness(Phy.f_std, x="Status2", measures="Shannon",color = "Status2")+
  geom_boxplot(alpha=0.6,aes(fill=sample_data(Phy.f_std)$Status2))+ facet_wrap(Study_site~Visit)+theme_bw()+ xlab("HIV exposure status")+
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
        strip.text = element_text(size=17, face="bold"), axis.text.y=element_text(size=16, face="bold"),
        text = element_text(size=17, face="bold"))+ 
  scale_color_manual(values=c("#B05A7A","#F99417"))+
  scale_fill_manual(values=c("#B05A7A","#F99417"))+
  labs(y= "Shannon Index")

#A84448
# F99417 orange
# FF78F0

# and add statistics
# https://rpubs.com/lconteville/713954
a_my_comparisons <- list( c("iHEU", "iHUU"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
Fig3_A <- p1 + stat_compare_means(method = "wilcox.test", 
                                comparisons = a_my_comparisons, 
                                label = "p.signif", symnum.args = symnum.args, size  = 5)+ scale_y_continuous(expand = expansion(mult = c(0, 0.1))) 

Fig3_A
ggsave(file = "./figure_output/fig3_alpha_exposure_comparison.v2.pdf", plot = Fig3_A, dpi = 500, width = 7, height = 7) # in this way, I can control the size of the figure
dev.off()


#### PERMANOVA (HIV exposure)
set.seed(20)
Phy.f_std_SA <- subset_samples(Phy.f_std, Study_site=="South Africa")
diss <- phyloseq::distance(Phy.f_std_SA, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_SA_df <- data.frame(row.names=sample_names(Phy.f_std_SA),sample_data(Phy.f_std_SA)) #create dataframe with metadata
# Whether "SequenceLot" differ significantly from each other
adonis(diss ~ Status2, data=Phy.f_std_SA_df, permutations = 999, by = "terms", na.omit)  # p=0.719



set.seed(20)
Phy.f_std_Nigeria <- subset_samples(Phy.f_std, Study_site=="Nigeria")
diss <- phyloseq::distance(Phy.f_std_Nigeria, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_Nigeria_df <- data.frame(row.names=sample_names(Phy.f_std_Nigeria),sample_data(Phy.f_std_Nigeria)) #create dataframe with metadata
# Whether "SequenceLot" differ significantly from each other
adonis(diss ~ Status2, data=Phy.f_std_Nigeria_df, permutations = 999, by = "terms", na.omit)  # p=0.719





####### Figure 3 (B) Heatmap ------
library(RColorBrewer)
library(pheatmap)
heatmap_stand_otu <- function(physeq) {
  #calculate weighted unifrac for clustering
  set.seed(2)
  glom_random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
  glom_tree = merge_phyloseq(physeq,glom_random_tree)
  diss <- phyloseq::distance(glom_tree,method = "bray", type = "samples")
  #use top 25 most abundant taxa for plotting
  glom_25 = prune_taxa(names(sort(taxa_sums(physeq), TRUE))[1:24], physeq)  # 24 instead of 25 as I mergerd two ASV together
  otus <- otu_table(glom_25)
  #log2 transform otus for better heatmap visualisation colour scale
  zs <- which(otus==0)#find zero entries 
  if(length(zs)>0){#if there are zero values, transform
    otus[zs] <- 0.1 #set zs to 0.1
    otus <- log2(otus)
  } else{otus <- log2(otus)}
  plot <- otus #data to be plotted
  #get metadata for annotation 
  meta_hm <- data.frame(row.names=sample_names(physeq),sample_data(physeq)$CST_pam_3_bray,sample_data(physeq)$Status2,sample_data(physeq)$Visit, sample_data(physeq)$Study_site)
  names(meta_hm) <- c("Community cluster", "HIV exposure status","Visit","Study site")
  #annotation 
  labs <- tax.lab(physeq,otus=otus)
  labs[labs == "B.longum"] <- "B.longum (subsp. infantis)" # ASV1
  labs[labs == "B.longum.1"] <- "B.longum (subsp. longum)" # ASV4
  #color
  col.pal <- brewer.pal(9,"YlOrRd")
  annot.color.col <- list("Community cluster"=c("Cluster 1"="#8E999F","Cluster 2"="#F1948A","Cluster 3"="#C5E1A5"),"HIV exposure status" =c("iHEU" ="#B05A7A", "iHUU"="#F99417"),"Visit" = c("Week 1"="#4FC3F7","Week 15" = "#48648F"),"Study site" = c("South Africa" = "#E87D72", "Nigeria" = "#54BCC2"))
  #heatmap
  hm <- pheatmap(plot, annotation = meta_hm, color = col.pal, clustering_method = "complete", clustering_distance_cols = diss, show_colnames = FALSE, labels_row = labs, annotation_colors = annot.color.col, fontsize = 13)
  return(hm)
} 


# here is my "Phy.f_std" phyloseq  with PAM clustering annotation
head(sample_data(Phy.f_std)$CST_pam_3_bray) 

# HEATMAP (1) - using "Phy.f_std"(all sites / all visits) - Figure 3
# levels(sample_data(Phy.f_std)$"Study_site") <- c("CT"="South Africa","Nigeria"="Nigeria") # need to change the annotation
# levels(sample_data(Phy.f_std)$Status2) <- c("HEU"="iHEU","HU"="iHUU") # need to change the annotation
levels(sample_data(Phy.f_std)$"Mode_of_feeding") <- c("Exclusive_breastfeeding"="Exclusive breastfeeding","Mixed_feeding"="Mixed feeding") 

# for better annotation
merged_phy1 = merge_taxa(Phy.f_std, c("ASV12","ASV22")) # Staphylococcus saprophyticus

Fig3_B <- heatmap_stand_otu(merged_phy1) 
Fig3_B

ggsave(file = "./figure_output/fig3_heatmap_AllData.top25.bray.v3.pdf", plot = Fig3_B, dpi = 800, width = 9, height = 7) 
dev.off()


### for better annotation
glom_25 = prune_taxa(names(sort(taxa_sums(Phy.f_std), TRUE))[1:25], Phy.f_std) 
otus <- otu_table(glom_25)
rownames(tax_table(glom_25)) # these are the top 20 taxa (for birth data)
# [1] "ASV1"  "ASV2"  "ASV3"  "ASV4"  "ASV5"  "ASV6"  "ASV7"  "ASV8"  "ASV9"  "ASV10" "ASV11" "ASV12"
# [13] "ASV13" "ASV16" "ASV18" "ASV19" "ASV20" "ASV21" "ASV22" "ASV24" "ASV26" "ASV30" "ASV31" "ASV37"
# [25] "ASV45"

dput(rownames(tax_table(glom_25))) # output as vector
my.ASVlist <- c("ASV1", "ASV2", "ASV3", "ASV4", "ASV5", "ASV6", "ASV7", "ASV8", 
                "ASV9", "ASV10", "ASV11", "ASV12", "ASV13", "ASV16", "ASV18", 
                "ASV19", "ASV20", "ASV21", "ASV22", "ASV24", "ASV26", "ASV30", 
                "ASV31", "ASV37", "ASV45")

## Get ASV info
top.W1W15.ASV <- subset_taxa(Phy.f_std, rownames(tax_table(Phy.f_std)) %in% my.ASVlist)
top.W1W15.ASV_seq <- refseq(top.W1W15.ASV)
ASV <- as.data.frame(top.W1W15.ASV_seq)
ASV$ASV <- rownames(ASV)
ASV <-  tibble(ASV)

# get annotation
tax_list <- read.csv("/Users/saoriiwase/Desktop/tax_list_BEAMING.csv") # obtain the tax_list (for full tax info)
names(tax_list)
tax_list <- tibble(tax_list) %>% mutate(ASV = X)

ASV.v2 <- left_join(ASV, tax_list, by = "ASV") # merge
ASV.v3 <- ASV.v2 %>% dplyr::rename("Sequence" =  "x") %>% select(-X)
# View(ASV.v3)

df <- data.frame(ASV.v3, check.names = T)
print(df)
#write_xlsx(df,"./figure_output/fig3_W1W15.top25.ASV.Seq.list.xlsx")




###### Figure 3 (C) DEseq2 --------
## DESeq2 (8') - Exposure status (HEU vs HU) comparison by Study site (CT  BIRTH) - p= 0.01 --- Fig 3 (B)
#levels(sample_data(Phy.f_std)$"Study_site") <- c("CT"="South Africa","Nigeria"="Nigeria") 
levels(sample_data(Phy.f_std)$"Status2")
Phy.f_std_CT_birth <- subset_samples(Phy.f_std, Study_site=="South Africa" &  Visit == "Week 1")

sample_data(Phy.f_std_CT_birth)
Phy.f_std_CT_birth # 63 samples

#Convert the phyloseq object to a DESeqDataSet and run DESeq2 (use "Phy.f_std_CT_birth"): 
ds = phyloseq_to_deseq2(Phy.f_std_CT_birth, ~ Status2)
ds = DESeq(ds, fitType = "local")  # get an error...rror in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : every gene contains at least one zero, cannot compute log geometric means
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(ds), 1, gm_mean)
ds = estimateSizeFactors(ds, geoMeans = geoMeans)
ds = DESeq(ds, fitType="local")
res = results(ds)
res = res[order(res$padj, na.last = NA), ]
alpha = 0.01 #this is your significance level, adjust if desired
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Phy.f_std_CT_birth)[rownames(sigtab), ], "matrix"))
head(sigtab)
sigtab %>% View()

#only show bugs annotated at the Species level
#sigtabspec = subset(sigtab, !is.na(Species))
sigtabgen = subset(sigtab, !is.na(Genus))
#sigtabord = subset(sigtab, !is.na(Order))
head(sigtabgen)

#write.csv(sigtabgen, file = "./figure_output/sigtabgen_CT_birth_p0.01_HIVexposure.csv")


#sort for pretty plotting
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

#Alternative: Phylum level, for plotting purposes - don't do both
#x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
#x = sort(x, TRUE)
#sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

#fix taxrank and subset to highly significant asvs
sigtabgen.hs <- sigtabgen[which(sigtabgen$padj < 0.005),]

#subset to greater magnitude of FC
sigtabgen$Species[is.na(sigtabgen$Species)] = "(unclassified)"
sigtabgen$genspec <- paste(sigtabgen$Genus, sigtabgen$Species, sep = " ")
sigtab.merged.hs <- sigtabgen[which(sigtabgen$padj < 0.005),]

sigtabgen %>%
  subset(log2FoldChange > 0.5) -> sigtabgen.hi
sigtabgen %>%
  subset(log2FoldChange < -0.5) -> sigtabgen.low

sigtab.merged.hs %>%
  subset(log2FoldChange > 0.5) -> sigtabgen.hi

sigtab.merged.hs %>%
  subset(log2FoldChange < -0.5) -> sigtabgen.low

sigtabgen.fc <- rbind.data.frame(sigtabgen.hi, sigtabgen.low)  
DEseq2_CT_Birth_HICexposure <- sigtabgen.fc
View(sigtabgen.fc)

#plot the results (Genus currently set, and colored by Order, set as desired)
p1.1 <- ggplot(sigtabgen.fc, aes(y=Genus, x=log2FoldChange, fill=genspec)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_vline(xintercept = 0.5, color = "indianred3", size = 0.5, lty=2) +
  geom_vline(xintercept = -0.5, color = "indianred3", size = 0.5, lty=2) +
  geom_point(pch = 21, size=3) + theme_bw()+
  theme(axis.text.y = element_text(size = 11, colour = "black",face = "bold"),
        axis.title.y = element_text(size = 15, colour = "black",face = "bold"),
        axis.text.x = element_text(size = 9, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black",face = "bold"),
        strip.background =element_rect(fill="white"),
        strip.text.x = element_text(colour = "mediumpurple4", size = 15, hjust = 0.5, vjust = 0.5, face = "bold"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size=0)) + 
  #scale_fill_manual(values=custom_colors_36) +
  guides(fill = guide_legend(nrow = 12)) + 
  labs(x=expression(log[2]~fold~change)) +  ggtitle("South Africa (Birth visit)")


ggsave(file = "./figure_output/fig3_DESeq2_HEUvsHU_CTbirth.pdf", plot = p1.1, dpi = 100, width = 10, height = 7) # in this way, I can control the size of the figure
dev.off()





## DESeq2 (8'') - Exposure status (HEU vs HU) comparison by Study site (CT  W15) - p= 0.01  - figure 3
levels(sample_data(Phy.f_std)$"Study_site") <- c("CT"="South Africa","Nigeria"="Nigeria") # need to change the name
Phy.f_std_CT_W15 <- subset_samples(Phy.f_std, Study_site=="South Africa" &  Visit == "Week 15")

sample_data(Phy.f_std_CT_W15)
Phy.f_std_CT_W15 # 66 samples

#Convert the phyloseq object to a DESeqDataSet and run DESeq2 (use "Phy.f_std_CT_W15"):
ds = phyloseq_to_deseq2(Phy.f_std_CT_W15, ~ Status2)
ds = DESeq(ds, fitType = "local")  # get an error...rror in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : every gene contains at least one zero, cannot compute log geometric means
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(ds), 1, gm_mean)
ds = estimateSizeFactors(ds, geoMeans = geoMeans)
ds = DESeq(ds, fitType="local")
res = results(ds)
res = res[order(res$padj, na.last = NA), ]
alpha = 0.01 #this is your significance level, adjust if desired
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Phy.f_std_CT_W15)[rownames(sigtab), ], "matrix"))
head(sigtab)
sigtab %>% View()

#only show bugs annotated at the Species level
#sigtabspec = subset(sigtab, !is.na(Species))
sigtabgen = subset(sigtab, !is.na(Genus))
#sigtabord = subset(sigtab, !is.na(Order))
head(sigtabgen)

write.csv(sigtabgen, file = "./figure_output/sigtabgen_CT_W15_p0.01_HIVexposure.csv")


#sort for pretty plotting
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

#Alternative: Phylum level, for plotting purposes - don't do both
#x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
#x = sort(x, TRUE)
#sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

#fix taxrank and subset to highly significant asvs
sigtabgen.hs <- sigtabgen[which(sigtabgen$padj < 0.005),]

#subset to greater magnitude of FC
sigtabgen$Species[is.na(sigtabgen$Species)] = "(unclassified)"
sigtabgen$genspec <- paste(sigtabgen$Genus, sigtabgen$Species, sep = " ")
sigtab.merged.hs <- sigtabgen[which(sigtabgen$padj < 0.005),]

sigtabgen %>%
  subset(log2FoldChange > 0.5) -> sigtabgen.hi
sigtabgen %>%
  subset(log2FoldChange < -0.5) -> sigtabgen.low

sigtab.merged.hs %>%
  subset(log2FoldChange > 0.5) -> sigtabgen.hi

sigtab.merged.hs %>%
  subset(log2FoldChange < -0.5) -> sigtabgen.low

sigtabgen.fc <- rbind.data.frame(sigtabgen.hi, sigtabgen.low)                          
View(sigtabgen.fc)
DEseq2_CT_W15_HICexposure <- sigtabgen.fc

#plot the results (Genus currently set, and colored by Order, set as desired)
p1.1 <- ggplot(sigtabgen.fc, aes(y=Genus, x=log2FoldChange, fill=genspec)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_vline(xintercept = 0.5, color = "indianred3", size = 0.5, lty=2) +
  geom_vline(xintercept = -0.5, color = "indianred3", size = 0.5, lty=2) +
  geom_point(pch = 21, size=3) + theme_bw()+
  theme(axis.text.y = element_text(size = 14, colour = "black",face = "bold"),
        axis.title.y = element_text(size = 17, colour = "black",face = "bold"),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 17, colour = "black",face = "bold"),
        strip.background =element_rect(fill="white"),
        strip.text.x = element_text(colour = "mediumpurple4", size = 15, hjust = 0.5, vjust = 0.5, face = "bold"),
        legend.text = element_text(size = 14, colour = "black"),
        legend.title = element_text(size=0)) + 
  #scale_fill_manual(values=custom_colors_36) +
  guides(fill = guide_legend(nrow = 12)) + 
  labs(x=expression(log[2]~fold~change)) #+  ggtitle("South Africa (W15 visit)")

ggsave(file = "./figure_output/fig3_DESeq2_HEUvsHU_CTW15.pdf", plot = p1.1, dpi = 100, width = 10, height = 7) # in this way, I can control the size of the figure
dev.off()




## DESeq2 (9') - Exposure status (HEU vs HU) comparison by Study site (Nigeria Birth) - p= 0.01  - Figure 3
Phy.f_std_Nigeria_birth <- subset_samples(Phy.f_std, Study_site=="Nigeria" & Visit =="Week 1")
sample_data(Phy.f_std_Nigeria_birth)
Phy.f_std_Nigeria_birth # 141 samples

#Convert the phyloseq object to a DESeqDataSet and run DESeq2 (use "Phy.f_std_Nigeria_birth"):
ds = phyloseq_to_deseq2(Phy.f_std_Nigeria_birth, ~ Status2)
ds = DESeq(ds, fitType = "local")  # get an error...rror in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : every gene contains at least one zero, cannot compute log geometric means
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(ds), 1, gm_mean)
ds = estimateSizeFactors(ds, geoMeans = geoMeans)
ds = DESeq(ds, fitType="local")
res = results(ds)
res = res[order(res$padj, na.last = NA), ]
alpha = 0.01 #this is your significance level, adjust if desired
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Phy.f_std_Nigeria_birth)[rownames(sigtab), ], "matrix"))
head(sigtab)
sigtab %>% View()

#only show bugs annotated at the Species level
#sigtabspec = subset(sigtab, !is.na(Species))
sigtabgen = subset(sigtab, !is.na(Genus))
#sigtabord = subset(sigtab, !is.na(Order))
head(sigtabgen)

write.csv(sigtabgen, file = "./figure_output/sigtabgen_Nigeria_birth_p0.01_HIVexposure.csv")


#sort for pretty plotting
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

#Alternative: Phylum level, for plotting purposes - don't do both
#x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
#x = sort(x, TRUE)
#sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

#fix taxrank and subset to highly significant asvs
sigtabgen.hs <- sigtabgen[which(sigtabgen$padj < 0.005),]

#subset to greater magnitude of FC
sigtabgen$Species[is.na(sigtabgen$Species)] = "(unclassified)"
sigtabgen$genspec <- paste(sigtabgen$Genus, sigtabgen$Species, sep = " ")
sigtab.merged.hs <- sigtabgen[which(sigtabgen$padj < 0.005),]

sigtabgen %>%
  subset(log2FoldChange > 0.5) -> sigtabgen.hi
sigtabgen %>%
  subset(log2FoldChange < -0.5) -> sigtabgen.low

sigtab.merged.hs %>%
  subset(log2FoldChange > 0.5) -> sigtabgen.hi

sigtab.merged.hs %>%
  subset(log2FoldChange < -0.5) -> sigtabgen.low

sigtabgen.fc <- rbind.data.frame(sigtabgen.hi, sigtabgen.low)  
DEseq2_Nigeria_birth_HICexposure <- sigtabgen.fc
View(sigtabgen.fc)

#plot the results (Genus currently set, and colored by Order, set as desired)
p1.1 <- ggplot(sigtabgen.fc, aes(y=Genus, x=log2FoldChange, fill=genspec)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_vline(xintercept = 0.5, color = "indianred3", size = 0.5, lty=2) +
  geom_vline(xintercept = -0.5, color = "indianred3", size = 0.5, lty=2) +
  geom_point(pch = 21, size=3) + theme_bw()+
  theme(axis.text.y = element_text(size = 11, colour = "black",face = "bold"),
        axis.title.y = element_text(size = 15, colour = "black",face = "bold"),
        axis.text.x = element_text(size = 9, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black",face = "bold"),
        strip.background =element_rect(fill="white"),
        strip.text.x = element_text(colour = "mediumpurple4", size = 15, hjust = 0.5, vjust = 0.5, face = "bold"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size=0)) + 
  #scale_fill_manual(values=custom_colors_36) +
  guides(fill = guide_legend(nrow = 12)) + 
  labs(x=expression(log[2]~fold~change)) +  ggtitle("Nigeria (Birth visit)")


ggsave(file = "./figure_output/fig3_DESeq2_HEUvsHU_NigeriaBirth.pdf", plot = p1.1, dpi = 100, width = 10, height = 7) # in this way, I can control the size of the figure
dev.off()




## DESeq2 (9'') - Exposure status (HEU vs HU) comparison by Study site (Nigeria W15) - p= 0.01  --- Figure 3
Phy.f_std_Nigeria_W15 <- subset_samples(Phy.f_std, Study_site=="Nigeria" & Visit =="Week 15")
sample_data(Phy.f_std_Nigeria_W15)
Phy.f_std_Nigeria_W15 # 172 samples

#Convert the phyloseq object to a DESeqDataSet and run DESeq2 (use "Phy.f_std_Nigeria_W15"):
ds = phyloseq_to_deseq2(Phy.f_std_Nigeria_W15, ~ Status2)

ds = DESeq(ds, fitType = "local")  # get an error...rror in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : every gene contains at least one zero, cannot compute log geometric means
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(ds), 1, gm_mean)
ds = estimateSizeFactors(ds, geoMeans = geoMeans)

ds = DESeq(ds, fitType="local")
res = results(ds)
res = res[order(res$padj, na.last = NA), ]
alpha = 0.01 #this is your significance level, adjust if desired
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Phy.f_std_Nigeria_W15)[rownames(sigtab), ], "matrix"))
head(sigtab)
sigtab %>% View()

#only show bugs annotated at the Species level
#sigtabspec = subset(sigtab, !is.na(Species))
sigtabgen = subset(sigtab, !is.na(Genus))
#sigtabord = subset(sigtab, !is.na(Order))
head(sigtabgen)

write.csv(sigtabgen, file = "./figure_output/sigtabgen_Nigeria_W15_p0.01_HIVexposure.csv")


#sort for pretty plotting
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

#Alternative: Phylum level, for plotting purposes - don't do both
#x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
#x = sort(x, TRUE)
#sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

#fix taxrank and subset to highly significant asvs
sigtabgen.hs <- sigtabgen[which(sigtabgen$padj < 0.005),]

#subset to greater magnitude of FC
sigtabgen$Species[is.na(sigtabgen$Species)] = "(unclassified)"
sigtabgen$genspec <- paste(sigtabgen$Genus, sigtabgen$Species, sep = " ")
sigtab.merged.hs <- sigtabgen[which(sigtabgen$padj < 0.005),]

sigtabgen %>%
  subset(log2FoldChange > 0.5) -> sigtabgen.hi
sigtabgen %>%
  subset(log2FoldChange < -0.5) -> sigtabgen.low

sigtab.merged.hs %>%
  subset(log2FoldChange > 0.5) -> sigtabgen.hi

sigtab.merged.hs %>%
  subset(log2FoldChange < -0.5) -> sigtabgen.low

sigtabgen.fc <- rbind.data.frame(sigtabgen.hi, sigtabgen.low)      
DEseq2_Nigeria_W15_HICexposure <- sigtabgen.fc
View(sigtabgen.fc)

#plot the results (Genus currently set, and colored by Order, set as desired)
p1.1 <- ggplot(sigtabgen.fc, aes(y=Genus, x=log2FoldChange, fill=genspec)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_vline(xintercept = 0.5, color = "indianred3", size = 0.5, lty=2) +
  geom_vline(xintercept = -0.5, color = "indianred3", size = 0.5, lty=2) +
  geom_point(pch = 21, size=3) + theme_bw()+
  theme(axis.text.y = element_text(size = 11, colour = "black",face = "bold"),
        axis.title.y = element_text(size = 15, colour = "black",face = "bold"),
        axis.text.x = element_text(size = 9, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black",face = "bold"),
        strip.background =element_rect(fill="white"),
        strip.text.x = element_text(colour = "mediumpurple4", size = 15, hjust = 0.5, vjust = 0.5, face = "bold"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size=0)) + 
  #scale_fill_manual(values=custom_colors_36) +
  guides(fill = guide_legend(nrow = 12)) + 
  labs(x=expression(log[2]~fold~change)) +  ggtitle("Nigeria (W15 visit)")


ggsave(file = "./figure_output/fig3_DESeq2_HEUvsHU_NigeriaW15.pdf", plot = p1.1, dpi = 100, width = 10, height = 7) # in this way, I can control the size of the figure
dev.off()


### DEseq2 figure Fig3 (B) - combined ver
### For Fig 3B, I want to comobine the figures so that I would be easier to see and compare (instead of put 4 figures separately).
# DEseq2 results for each site and visit were saved (above). Need to add "tag" so that it would be easier to hadle
DEseq2_CT_Birth_HICexposure   # DEseq2 result for CT birth HEU vs HUU
DEseq2_CT_Birth_HICexposure$tag <- "CT_W1_HEUvsHUU"

DEseq2_CT_W15_HICexposure  #  DEseq2 result for CT W15 HEU vs HUU
DEseq2_CT_W15_HICexposure$tag <- "CT_W15_HEUvsHUU"

DEseq2_Nigeria_birth_HICexposure   # DEseq2 result for Nigeria birth HEU vs HUU
DEseq2_Nigeria_birth_HICexposure$tag <- "Nigeria_W1_HEUvsHUU"

DEseq2_Nigeria_W15_HICexposure
DEseq2_Nigeria_W15_HICexposure$tag <- "Nigeria_W15_HEUvsHUU"

head(DEseq2_CT_Birth_HICexposure)
# now combine these 4 datasets together
data <- rbind(DEseq2_CT_Birth_HICexposure, DEseq2_CT_W15_HICexposure, DEseq2_Nigeria_birth_HICexposure, DEseq2_Nigeria_W15_HICexposure)
View(data)

length(data$genspec) # 33 ASVs
length(unique(data$genspec)) # 28 ASVs
data$genspec[duplicated(data$genspec)] 
# [1] "Prevotella copri"          "Prevotella copri"          "Prevotella copri"         
# [4] "Collinsella aerofaciens"   "Escherichia-Shigella coli"

str(data)
sample_data(Phy.f_std)$Status2  # Levels:iHEU iHUU
levels(sample_data(Phy.f_std_birth)$"Study_site") <- c("CT"="South Africa","Nigeria"="Nigeria")
data$tag <- as.factor(data$tag)
levels(data$tag) <- c("CT_W1_HEUvsHUU" = "South Africa at Week 1", 
                      "CT_W15_HEUvsHUU" = "South Africa at Week  15",
                      "Nigeria_W1_HEUvsHUU" = "Nigeria at Week 1",
                      "Nigeria_W15_HEUvsHUU" ="Nigeria at Week 15")

# save the dataframe for supp info
df <- data.frame(data)
write.csv(df, file = "./figure_output/sigtabgen_DESeq2_HEUvsHU_CombinedVer_p0.01_HIVexposure.csv")

# now, let's plot the DEseq result (HEU vs HUU)
df$tag[df$tag == "CT_W1_HEUvsHUU"] <- "South Africa at Week 1"
df$tag[df$tag == "CT_W15_HEUvsHUU"] <- "South Africa at Week 15"
df$tag[df$tag == "Nigeria_W1_HEUvsHUU"] <- "Nigeria at Week 1"
df$tag[df$tag == "Nigeria_W15_HEUvsHUU"] <- "Nigeria at Week 15"

Fig3_C <- ggplot(df, aes(y=Genus, x=log2FoldChange, fill=genspec)) + 
  facet_wrap(~tag, ncol = 4)+
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_vline(xintercept = 0.5, color = "indianred3", size = 0.5, lty=2) +
  geom_vline(xintercept = -0.5, color = "indianred3", size = 0.5, lty=2) +
  geom_point(pch = 21, size=3) + theme_bw()+
  theme(axis.text.y = element_text(size = 12.5, colour = "black",face = "bold"),
        axis.title.y = element_text(size = 17, colour = "black",face = "bold"),
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.title.x = element_text(size = 17, colour = "black",face = "bold"),
        strip.background =element_rect(fill="white"),
        strip.text.x = element_text(colour = "black", size = 17, hjust = 0.5, vjust = 0.5, face = "bold"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size=0)) + theme(legend.position="bottom", legend.title=element_blank())+
  #scale_fill_manual(values=custom_colors_36) +
  #guides(fill = guide_legend(ncol = 6)) + 
  labs(x=expression(log[2]~fold~change)) + 
  #ggtitle("Bacteria in iHUU (compared to iHEU as base)") + 
  theme(plot.title = element_text(size = 17, face = "bold"), 
        panel.spacing = unit(1.5, "lines")
        )

Fig3_C


ggsave(file = "./figure_output/fig3_DESeq2_HEUvsHU_CombinedVer.pdf", plot = Fig3_C, dpi = 800, width = 18, height = 8) # in this way, I can control the size of the figure
dev.off()


### Beta-diversity, comparing HIV exposure--------------
set.seed(20) #select any start seed, in this case 2 (ensures reproducibility)
ord.BC <- ordinate(Phy.f_std, method = "PCoA", distance = "bray", k=2, trymax=1000) # try PcOA as NMDS does not work


color = c("Status2")
shape = c("Study_site")
#title= c("PCoA of 16S microbiome (all visit), Jensen-Shannon distance, k=2")  
MDS = plot_ordination(Phy.f_std, ord.BC, color = color, shape=shape)
Fig2S  = MDS +theme_bw() +theme(axis.text=element_text(size=17, face="bold"),
                                 axis.title=element_text(size=17,face="bold"), 
                                 legend.title=element_text(size=17 ,face="bold"),
                                 legend.text = element_text(size =16))+ 
  labs(color=color, shape=shape)+geom_point(size=2)+ 
  guides(color = guide_legend(title = "Exposure status"), shape = guide_legend(title = "Study site")) +
  scale_color_manual(values=c("#B05A7A","#F99417"))

Fig2S

ggsave(file = "./figure_output/fig3_birth.W15_PCoA_k2_exposure.resize.bray.color.v2.pdf", plot = Fig2S, dpi = 600, width = 8.5, height = 5) # in this way, I can control the size of the figure
dev.off()

set.seed(20)
diss <- phyloseq::distance(Phy.f_std, "bray", parallel=TRUE) 
Phy.f_std_df <- data.frame(row.names=sample_names(Phy.f_std),sample_data(Phy.f_std)) #create dataframe with metadata
# Whether "SequenceLot" differ significantly from each other
adonis2(diss ~ Status2, data=Phy.f_std_df, permutations = 999, by = "terms", na.omit)  # p=



########## Figure 4


#### Figure 4 A (scatter plot) -------
Tetanus <- read_excel("/Users/saoriiwase/Desktop/BEAMING study/BEAMING_microbiome_data_analysis/for paper/Tetanus data/*Tetanus IgG ELISA data_5_3_27 CLEAN_27Jan.xlsx", 
                      sheet = "Tetanus analysis_long")

Tetanus # A tibble: 602 × 7
names(Tetanus)
## [1] "PID"             "Status2"         "Study_site"      "Visit"           "Sample_ID"      
## [6] "Baby_titer"      "Mum_Titer_birth"

length(unique(Tetanus$PID))  # 301 unique PID
Tetanus$PID[duplicated(Tetanus$Sample_ID)]  # character(0) - if it list something, then fix it

# dealing with data types
str(Tetanus) # Baby_titer  and Mum's titer needs to be number. No < or  > sign in the colum (substitute)
Tetanus$Baby_titer[which(Tetanus$Baby_titer == "<0.1")] <- "0.099"  # replace the <0.1
Tetanus$Baby_titer <- as.numeric(Tetanus$Baby_titer)  # change it to numeric
Tetanus$Baby_titer <- round(Tetanus$Baby_titer,3)  # change the number of digit

Tetanus$Status2 <- as.factor(Tetanus$Status2)

Tetanus %>% View() # check if it worked.

stool_PID_info <- data.frame(sample_data(Phy.f_std))$PID # list of PIDs that we used for analysis

### Let's filter only PIDs that have both stool library and tetanus titer 
Tetanus # original. # A tibble: 602 × 7
Tetanus2 <- Tetanus %>% filter(PID %in% stool_PID_info) # select only PIDs that are involved in the stool samples. # A tibble: 578 × 7

# extract only birth data
Tetanus_birth <- Tetanus2 %>% filter(Visit == "Birth") # A tibble: 289 × 7
#Tetanus_birth %>% View()

### Scatter plots ###
# Plot showing Placental transfer (with spearman analysis) - fig4 (A)
levels(Tetanus_birth$Status2) <- c("HEU"="iHEU","HU"="iHUU") # need to change the name

#with Pearson correlation plot (Susan was not happy about using this method as some of the data was satulated)
Fig4_A <- ggscatter(data=Tetanus_birth, x = "Baby_titer", y = "Mum_titer_birth",   
                add = "reg.line", conf.int = TRUE, size=2) +
  stat_cor(method = "pearson", label.x = 0, label.y = 6.2, size = 5.5)+ 
  facet_wrap(~Status2)+ theme(axis.text=element_text(size=17, face="bold"),
                              axis.title=element_text(size=17,face="bold"),
                              strip.text = element_text(size=17,face="bold"))+
  labs(x= "Anti-TT IgG at W1 (IU/ml) [infant]") +
  labs(y= "Anti-TT IgG (IU/ml) [mother]") 
  
  
#ggsave(file = "./figure_output/Fig4_correlation_plot_BabyvsMum.pdf", plot = Fig4_A, dpi = 500, width =7 , height = 6) # in this way, I can control the size of the figure
dev.off()

# with "Spearman correlation" plot  - better?
ggscatter(data=Tetanus_birth, x = "Baby_titer", y = "Mum_titer_birth",   
          add = "reg.line", conf.int = TRUE, size=2) +
  stat_cor(method = "spearman", label.x = 0, label.y = 6.2, size = 5.5)+ 
  facet_wrap(~Status2)+ theme(axis.text=element_text(size=17, face="bold"),
                              axis.title=element_text(size=17,face="bold"),
                              strip.text = element_text(size=17,face="bold"))+
  labs(x= "Anti-TT IgG at W1 (IU/ml) [infant]") +
  labs(y= "Anti-TT IgG (IU/ml) [mother]")  

Tetanus_birth_iHEU <- Tetanus_birth %>% filter(Status2 == "iHEU")
Tetanus_birth_iHUU <- Tetanus_birth %>% filter(Status2 == "iHUU")

# check the correlation coefficients with different method
# http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r
res1 <-cor.test(Tetanus_birth_iHEU$Baby_titer, Tetanus_birth_iHEU$Mum_titer_birth,  method = "spearman")
res1  # rho = 0.7156186, p-value < 2.2e-16
# Extract the p-value
p_value1 <- res1$p.value  # 6.062285e-23
p.adjust(p_value1, method = "BH")
p.adjust(p_value1, method = "bonferroni")

res2 <-cor.test(Tetanus_birth_iHUU$Baby_titer, Tetanus_birth_iHUU$Mum_titer_birth,  method = "spearman")
res2 # rho = 0.9457199, p-value < 2.2e-16
# Extract the p-value
p_value2 <- res2$p.value  # 1.493413e-26
p.adjust(p_value2, method = "BH")
p.adjust(p_value2, method = "bonferroni")





#### Figure 4 (A) - version 2 - color by status
Fig4_A.v2 <- ggscatter(data=Tetanus_birth, x = "Baby_titer", y = "Mum_titer_birth",   
                    add = "reg.line", conf.int = TRUE, size=2, color = "Status2") +
  #stat_cor(method = "pearson", label.x = 0, label.y = 6.2, size = 5.5)+ 
  # facet_wrap(~Status2)+ 
  theme(axis.text=element_text(size=17, face="bold"),
        axis.title=element_text(size=17,face="bold"),
        strip.text = element_text(size=17,face="bold"),
        legend.title=element_text(size=17),
        legend.text = element_text(size =17))+
  guides(color = guide_legend(title = "HIV exposure status")) +
  scale_color_manual(values=c("#B05A7A","#F99417"))+
  scale_fill_manual(values=c("#B05A7A","#F99417"))+
  labs(x= "Anti-TT IgG at W1 (IU/ml) [infant]") +
  labs(y= "Anti-TT IgG (IU/ml) [mother]") 


Fig4_A.v2
ggsave(file = "/Users/saoriiwase/Desktop/IHMC 9th CONGRESS 2022 KOBE//Fig4_correlation_plot_BabyvsMum.v3.pdf", plot = Fig4_A.v2, dpi = 700, width =6 , height = 6) # in this way, I can control the size of the figure
dev.off()

ggsave(file = "./figure_output/Fig4_correlation_plot_BabyvsMum.v3.pdf", plot = Fig4_A.v2, dpi = 700, width =6 , height = 6) # in this way, I can control the size of the figure
dev.off()


### Figure 4 (B) (box plot HEU vs HU) -------
# get metadata from the phyloseq object
meta <- tibble(data.frame(sample_data(Phy.f_std)))



#levels(meta$Study_site) <- c("CT"="South Africa","Nigeria"="Nigeria") 
#levels(meta$Status2) <- c("HEU" = "iHEU", "HU" = "iHUU")

names(meta)$Status2 <- meta$Status2 
p1 <- ggplot(data=meta, aes(x=Status2, y=Baby_titer,color = Status2))+ theme_bw()+
  geom_boxplot(alpha=0.6, aes(fill=Status2)) +geom_jitter(width=0.3, alpha=0.2)+ 
  facet_wrap(Visit~Study_site,scales = "free_y")+
  labs(x= "HIV exposure status") +labs(y= "Anti-TT IgG (IU/ml)") + 
  #ggtitle("Baby Tetanus titer iHEU vs iHUU") + 
  #guides(color = guide_legend(title = "Exposure status"))  + #ylim(0, 7.0) + 
  theme(axis.title=element_text(size=17,face="bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
        strip.text = element_text(size=17, face="bold"), 
        axis.text.y=element_text(size=17, face="bold"),
        text = element_text(size=17, face="bold"),
        legend.title=element_blank())+ 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  scale_color_manual(values=c("#B05A7A","#F99417"))+
  scale_fill_manual(values=c("#B05A7A","#F99417"))
  

# add stats
a_my_comparisons <- list( c("iHEU", "iHUU"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
Fig4_B <- p1 + stat_compare_means(method = "wilcox.test", 
                                comparisons = a_my_comparisons, 
                                label = "p.signif", symnum.args = symnum.args
) 


# adjust for multiple comparison...
anno_df <- compare_means(Baby_titer ~ Status2, meta, method = "wilcox.test", paired = FALSE,
                         group.by = c("Visit","Study_site"), ref.group = NULL,
                         p.adjust.method = "BH") # this shows the p-value was nnot significantly significant after adjusted for multiple comaprison.

anno_df$y.position <- c(6, 6,6,6) 
anno_df$xmin <- c(1,1,1,1)
anno_df$xmax <- c(2,2,2,2)

# I want to use adjdusted p-value and indicate so on the plot
# Define the significance levels and asterisks
alpha_levels <- c(0.05, 0.01, 0.001)
asterisks <- c("*", "**", "***")
asterisks <- c("***", "**", "*")

# Add the significance level and asterisks to the data frame
anno_df$p.adj.signif <- ifelse(anno_df$p.adj <= alpha_levels[1], asterisks[1],
                          ifelse(anno_df$p.adj <= alpha_levels[2], asterisks[2],
                                 ifelse(anno_df$p.adj <= alpha_levels[3], asterisks[3], "ns")))



Fig4_B <- p1 + stat_pvalue_manual(
  anno_df,  label = "p.adj.signif", tip.length = 0.01,
  bracket.nudge.y = 0.5, label.size = 4)


Fig4_B
ggsave(file = "./figure_output/fig4_Tetanus_barplot_HEUvsHU.color.BH.pdf", plot = Fig4_B, dpi = 500, width = 8, height = 7) # in this way, I can control the size of the figure
dev.off()


# find out the p-value to indicate in the manuscript
meta_W15_Nig_HEU <- meta %>% 
  select(Baby_titer,Visit, Study_site, Status2) %>% 
  filter(Study_site == "Nigeria") %>% 
  filter(Visit == "Week 15") %>% 
  filter(Status2 == "iHEU")

meta_W15_Nig_HUU <- meta %>% 
  select(Baby_titer,Visit, Study_site, Status2) %>% 
  filter(Study_site == "Nigeria") %>% 
  filter(Visit == "Week 15") %>% 
  filter(Status2 == "iHUU")

wilcox.test(meta_W15_Nig_HEU$Baby_titer, meta_W15_Nig_HUU$Baby_titer) #p-value = 0.04408


meta_W15_SA_HEU <- meta %>% 
  select(Baby_titer,Visit, Study_site, Status2) %>% 
  filter(Study_site == "South Africa") %>% 
  filter(Visit == "Week 15") %>% 
  filter(Status2 == "iHEU")

meta_W15_SA_HUU <- meta %>% 
  select(Baby_titer,Visit, Study_site, Status2) %>% 
  filter(Study_site == "South Africa") %>% 
  filter(Visit == "Week 15") %>% 
  filter(Status2 == "iHUU")

wilcox.test(meta_W15_SA_HEU$Baby_titer, meta_W15_SA_HUU$Baby_titer) #p-value = 0.1462




###### Figure 4B - version 2 - color by status
p1 <- ggplot(data=meta, aes(x=Status2, y=Baby_titer,color = Status2))+ theme_bw()+
  geom_boxplot() +geom_jitter(width=0.3, alpha=0.2)+ 
  facet_wrap(Study_site ~ Visit,scales = "free_y")+
  labs(x= "HIV exposure status") +labs(y= "Anti-TT IgG (IU/ml)") + 
  #ggtitle("Baby Tetanus titer iHEU vs iHUU") + 
  guides(color = guide_legend(title = ""))  + #ylim(0, 7.0) + 
  theme(axis.title=element_text(size=17,face="bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
        strip.text = element_text(size=17, face="bold"), 
        axis.text.y=element_text(size=17, face="bold"),
        text = element_text(size=17, face="bold"))+ 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))


# add stats
a_my_comparisons <- list( c("iHEU", "iHUU"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
Fig4_B.v2 <- p1 + stat_compare_means(method = "wilcox.test", 
                                     comparisons = a_my_comparisons, 
                                     label = "p.signif", symnum.args = symnum.args
) 
Fig4_B.v2
ggsave(file = "./figure_output/fig4_Tetanus_barplot_HEUvsHU.v2.pdf", plot = Fig4_B.v2, dpi = 500, width = 8, height = 7) # in this way, I can control the size of the figure
dev.off()



#### Figure 4C (corr. plot - shannon) --------
# tab <- microbiome::alpha(Phy.f_std, index = "shannon")
# sample_data(Phy.f_std)$Shannon <- tab$diversity_shannon 
#levels(df$Study_site) <- c("CT" = "South Africa", "Nigeria" = "Nigeria")

df <- data.frame(sample_data(Phy.f_std))
library(tibble)
sample_data(Phy.f_std)$Shannon

Fig4_C <- ggscatter(data=df, x = "Shannon", y = "Baby_titer", 
               add = "reg.line", conf.int = TRUE, size=2) +
  stat_cor(method = "pearson", size = 5.5)+ 
  facet_wrap(Visit~Study_site, nrow=2, scales = "free_y")  + theme_bw()+ 
  theme(axis.text=element_text(size=17,face = "bold"),
        axis.title=element_text(size=17,face="bold"), 
        strip.text = element_text(size=17, face="bold")) +ylim(0,5)+
  labs(x= "Shannon diversity") +labs(y= "Anti-TT IgG (IU/ml)") #+ ggtitle("Scatter plot of Tetanus titer vs Shannon diversity")
  
  
ggsave(file = "./figure_output/fig4_correlation_plot_Shannon&Tetanus.pdf", plot = Fig4_C, dpi = 500, width =7 , height = 7) # in this way, I can control the size of the figure
dev.off()


### Heather told me that I should plot W15 titer vs  W1 shannon /  W15 titer vs  W15 shannon
names(df)
W15_titer <- df %>% filter(Visit == "Week 15") %>%
  select(PID, Baby_titer) %>% 
  dplyr::rename("W15_titer" =  "Baby_titer")

df.2 <- left_join(df, W15_titer, by = "PID")
df.2 

Fig4_C.v2 <- ggscatter(data = df.2, x = "Shannon", y = "W15_titer", color = "Study_site",
          add = "reg.line", conf.int = TRUE, size=2) +
  stat_cor(method = "spearman", size = 5,  label.x = 1.4, label.y = 4.7)+ ylim(0,5)+xlim(0,4.1)+
  facet_wrap(Visit~Study_site, nrow=2, scales = "free_y")  + theme_bw()+ 
  theme(axis.text=element_text(size=17,face = "bold"),
        legend.title=element_text(size=17, face = "bold"),
        legend.text = element_text(size =17),
        axis.title=element_text(size=17,face="bold"), 
        strip.text = element_text(size=17, face="bold"))+
  labs(fill = "Study site")+
  labs(color = "Study site")+
  labs(x= "Shannon diversity") +labs(y= "Anti-TT IgG at W15 (IU/ml)") #+ ggtitle("Scatter plot of Tetanus titer vs Shannon diversity")

ggsave(file = "./figure_output/fig4_correlation_plot_Shannon&Tetanus.V2.spearman.color.pdf", plot = Fig4_C.v2, dpi = 800, width =9 , height = 8) # in this way, I can control the size of the figure
dev.off()

### Figure 4 (Deseq/ only in text) -----
Phy.f_std <- readRDS("Phy.f_std.RDS")
sample_data(Phy.f_std) # Baby_titer (either at birth/W15) and Mum_titer (sampled at birth visit)
levels(sample_data(Phy.f_std)$"Study_site") <- c("CT"="South Africa","Nigeria"="Nigeria") # need to change the name


# select only the PIDs that have vaccine data.
library(tidyr)
library(dplyr)
sdf <- tibble(data.frame(sample_data(Phy.f_std)))

df <- sdf %>% filter(Visit == "Week 1")
df <- sdf %>% filter(Visit == "Week 15")
median(sdf$Baby_titer, na.rm = T)  # 1.4

names(sdf)
class(sdf$Baby_titer) #numeric
class(sdf$Mum_titer) #numeric

# select only PIDs that have tetanus titer data
Tetanus <- sdf %>% filter(Baby_titer > 0 | Mum_titer > 0)  # A tibble: 430 × 81
View(Tetanus) # looks extracted accordingly.

Tetanus_birth <-  sdf %>% filter(Visit == "Birth")    # placental transfer
Tetanus_W15 <-  sdf %>% filter(Visit == "W15")    


# TT response3 - categorize into two *using median
sdf_all <- sdf %>% filter(Baby_titer > 0 | Mum_titer > 0) 
median(sdf_all$Baby_titer) #1.4

sample_data(Phy.f_std)$TTresponse3 <- case_when((sample_data(Phy.f_std)$Baby_titer >= 1.4) ~"High_responder", 
                                            (sample_data(Phy.f_std)$Baby_titer <1.4) ~ "Low_responder")

sample_data(Phy.f_std)$TTresponse3 <- as.factor(sample_data(Phy.f_std)$TTresponse3) # Levels: High_responder Low_responder




#### DESeq2 for Tetanus - p= 0.05 --- 
#levels(sample_data(Phy.f_std)$"Study_site") <- c("CT"="South Africa","Nigeria"="Nigeria") 
Phy.f_std_Tetanus <- subset_samples(Phy.f_std, TTresponse3 == "High_responder"  |  TTresponse3 == "Low_responder" )
sample_data(Phy.f_std_Tetanus)$TTresponse3 <- as.factor(sample_data(Phy.f_std_Tetanus)$TTresponse3)   # Levels: High_responder Low_responder
data.frame(sample_data(Phy.f_std_Tetanus)) %>% tibble() %>% View() 

sample_data(Phy.f_std_Tetanus)
Phy.f_std_Tetanus 

Phy.f_std_Tetanus <- subset_samples(Phy.f_std_Tetanus, Visit == "Week 15")

#Convert the phyloseq object to a DESeqDataSet and run DESeq2 (use "Phy.f_std_Tetanus"): 
ds = phyloseq_to_deseq2(Phy.f_std_Tetanus, ~ TTresponse3)  #  got warning
ds = DESeq(ds, fitType = "local")  
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(ds), 1, gm_mean)
ds = estimateSizeFactors(ds, geoMeans = geoMeans)
ds = DESeq(ds, fitType="local")
res = results(ds)
res = res[order(res$padj, na.last = NA), ]
alpha = 0.01 #this is your significance level, adjust if desired
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Phy.f_std_Tetanus)[rownames(sigtab), ], "matrix"))
head(sigtab)
sigtab %>% View()

#only show bugs annotated at the Species level
#sigtabspec = subset(sigtab, !is.na(Species))
sigtabgen = subset(sigtab, !is.na(Genus))
#sigtabord = subset(sigtab, !is.na(Order))
head(sigtabgen)

write.csv(sigtabgen, file = "./figure_output/sigtabgen_TTrespoonse3_p0.05.csv")
# Enterococcus faecalis
# Rhodococcus erythropolis (somehow this disappeared in the figure)

#sort for pretty plotting
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

#Alternative: Phylum level, for plotting purposes - don't do both
#x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
#x = sort(x, TRUE)
#sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

#fix taxrank and subset to highly significant asvs
sigtabgen.hs <- sigtabgen[which(sigtabgen$padj < 0.005),]

#subset to greater magnitude of FC
sigtabgen$Species[is.na(sigtabgen$Species)] = "(unclassified)"
sigtabgen$genspec <- paste(sigtabgen$Genus, sigtabgen$Species, sep = " ")
sigtab.merged.hs <- sigtabgen[which(sigtabgen$padj < 0.005),]

sigtabgen %>%
  subset(log2FoldChange > 0.5) -> sigtabgen.hi
sigtabgen %>%
  subset(log2FoldChange < -0.5) -> sigtabgen.low

sigtab.merged.hs %>%
  subset(log2FoldChange > 0.5) -> sigtabgen.hi

sigtab.merged.hs %>%
  subset(log2FoldChange < -0.5) -> sigtabgen.low

sigtabgen.fc <- rbind.data.frame(sigtabgen.hi, sigtabgen.low)                          
View(sigtabgen.fc)

#plot the results (Genus currently set, and colored by Order, set as desired)
Fig4_x <- ggplot(sigtabgen.fc, aes(y=Genus, x=log2FoldChange, fill=genspec)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_vline(xintercept = 0.5, color = "indianred3", size = 0.5, lty=2) +
  geom_vline(xintercept = -0.5, color = "indianred3", size = 0.5, lty=2) +
  geom_point(pch = 21, size=3) + 
  theme(axis.text.y = element_text(size = 11, colour = "black",face = "bold"),
        axis.title.y = element_text(size = 15, colour = "black",face = "bold"),
        axis.text.x = element_text(size = 9, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black",face = "bold"),
        strip.background =element_rect(fill="white"),
        strip.text.x = element_text(colour = "mediumpurple4", size = 15, hjust = 0.5, vjust = 0.5, face = "bold"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size=0)) + 
  #scale_fill_manual(values=custom_colors_36) +
  guides(fill = guide_legend(nrow = 12)) +
  labs(x=expression(log[2]~fold~change)) +   ggtitle("Tetanus High responder (base) vs Low responder")


ggsave(file = "./figure_output/fig4_TTrespnse3.pdf", plot = Fig4_x, dpi = 500, width = 10, height = 7) # in this way, I can control 




#-----------------------------------------------------------------------------------------------------
# Ven diagram  ----
#-----------------------------------------------------------------------------------------------------
# https://r-charts.com/part-whole/ggvenndiagram/
#install.packages("ggVennDiagram")
library(ggVennDiagram)

# List of items
x <- list(A = 1:5, B = 2:7)

# 2D Venn diagram
ggVennDiagram::ggVennDiagram(x)


set.seed(20190708)
genes <- paste("gene",1:1000,sep="")
x <- list(
  A = sample(genes,300), 
  B = sample(genes,525), 
  C = sample(genes,440),
  D = sample(genes,350)
)

ggVennDiagram(x)

df <- tibble(data.frame(sample_data(Phy.f_std)))
select_SA_W1 <- df %>% filter(Study_site == "South Africa" & Visit =="Week 1") %>% pull(PID)
select_SA_W15 <- df %>% filter(Study_site == "South Africa" & Visit =="Week 15") %>% pull(PID)

SA <- list(
  'Week 1' = select_SA_W1,
  'Week 15' = select_SA_W15
)

Ven1 <- ggVennDiagram(SA, set_color ="white", set_size =6, label_size = 6, label_alpha = 0) + 
  scale_color_manual(values = c("Week 1" = "yellow","Week 15" ="steelblue"))+
  theme_void()+
  theme(legend.position = "none") + scale_fill_gradient(low = "#FBC5C5", high = "#FF8B8B")



df <- tibble(data.frame(sample_data(Phy.f_std)))
select_Nig_W1 <- df %>% filter(Study_site == "Nigeria" & Visit =="Week 1") %>% pull(PID)
select_Nig_W15 <- df %>% filter(Study_site == "Nigeria" & Visit =="Week 15") %>% pull(PID)

Nig <- list(
  'Week 1' = select_Nig_W1,
  'Week 15' = select_Nig_W15
)

Ven2 <- ggVennDiagram(Nig, set_color ="white", set_size =6, label_size = 6, label_alpha = 0) + 
  scale_color_manual(values = c("Week 1" = "yellow","Week 15" ="steelblue"))+
  theme_void()+
  theme(legend.position = "none") + scale_fill_gradient(low = "#B1D7B4", high = "#319DA0")


#let'scombine ven diagrams
library(ggpubr)
Ven1
Ven2


Ven_all <- ggarrange(Ven1,NULL, Ven2 + rremove("x.text"), 
                     labels = c(), nrow = 1, widths = c(1, 0.05, 1))

ggsave(file = "./figure_output/SuppX_Ven.diagram_combine.pdf", plot = Ven_all, dpi = 800, width = 10 , height = 5) # in this way, I can control the size of the figure
dev.off()


#-----------------------------------------------------------------------------------------------------
#  PAM clustering (determine the k value ) ----
#-----------------------------------------------------------------------------------------------------
### PAM clustering
library(cluster)
Phy.f_std <- readRDS("Phy.f_std.RDS")
Phy.f_std_test <- Phy.f_std

# Creating CSTs the same way Bryan did, for consistency
# He used a modified version of Susan Holmes’ script: http://statweb.stanford.edu/~susan/papers/Pregnancy/PNAS_Vaginal_Analysis.html

#Cluster into CSTs
otu_table(Phy.f_std_test) <- t(otu_table(Phy.f_std_test))

phy.obj.rel <- Phy.f_std_test
phy.obj.rel <-transform_sample_counts(physeq = phy.obj.rel, fun = function(x) x/sum(x))
otu_table(phy.obj.rel)


#jsd_dist <- phyloseq::distance(phy.obj.rel, method="jsd")
jsd_dist <- phyloseq::distance(phy.obj.rel, method="bray")

ord = ordinate(phy.obj.rel, method = "PCoA", distance = jsd_dist)
plot_scree(ord) + xlim(as.character(seq(1,20))) + ggtitle("PCoA-jsd ordination eigenvalues")

evs <- ord$value$Eigenvalues
print(evs[1:20])

print(tail(evs))

#remove those below magnitude of largest negative eignevalue
h_sub10 <- hist(evs[10:length(evs)], 100)
plot(h_sub10$mids, h_sub10$count, log="y", type='h', lwd=10, lend=2)

#gap statistic for cluter number
NDIM <- 7 #11 falls below mag of largest neg eig
x <- ord$vectors[,1:NDIM]  # rows=sample, cols=MDS axes, entries = value
pamPCoA = function(x, k) {
  list(cluster = pam(x[,1:2], k, cluster.only = TRUE))
}
gs = clusGap(x, FUN = pamPCoA, d.power = 2, K.max = 12, B = 500) #Tibshirani def of power

# Either k=3 or k=6 is the best value.
Elbow_plot <- plot_clusgap(gs) + scale_x_continuous(breaks=c(seq(0, 12, 2)))+ 
  theme_bw()+ theme() + ggtitle("") + theme(axis.text=element_text(size=17, face="bold"),
                                            axis.title=element_text(size=17,face="bold"), 
                                            legend.title=element_text(size=17),
                                            legend.text = element_text(size =17))

#ggsave(file = "./figure_output/SuppX.ElbowPlot.pdf", plot = Elbow_plot, dpi = 500, width = 5, height = 5) # in this way, I can control the size of the figure
dev.off()

#let's do k=3 (k=6 is too much?)
K <- 3
x <- ord$vectors[,1:NDIM]
clust <- as.factor(pam(x, k=K, cluster.only=T))
#add to sample data
sample_data(phy.obj.rel)$CST_pam_3 <- clust
CSTs <- as.character(seq(K))


#let's do k=4 (k=6 is too much?)
K <- 4
x <- ord$vectors[,1:NDIM]
clust <- as.factor(pam(x, k=K, cluster.only=T))
#add to sample data
sample_data(phy.obj.rel)$CST_pam_4 <- clust
CSTs <- as.character(seq(K))

#let's do k=5 (k=6 is too much?)
K <- 5
x <- ord$vectors[,1:NDIM]
clust <- as.factor(pam(x, k=K, cluster.only=T))
#add to sample data
sample_data(phy.obj.rel)$CST_pam_5 <- clust
CSTs <- as.character(seq(K))

#let's do k=6 (k=6 is too much?)
K <- 6
x <- ord$vectors[,1:NDIM]
clust <- as.factor(pam(x, k=K, cluster.only=T))
#add to sample data
sample_data(phy.obj.rel)$CST_pam_6 <- clust
CSTs <- as.character(seq(K))

#let's do k=7 (k=7 is too much?)
K <- 7
x <- ord$vectors[,1:NDIM]
clust <- as.factor(pam(x, k=K, cluster.only=T))
#add to sample data
sample_data(phy.obj.rel)$CST_pam_7 <- clust
CSTs <- as.character(seq(K))



#visually evaluate clustering
library(RColorBrewer)
CSTColors <- brewer.pal(8,"Paired")[c(1,3,2,5,4,6,7,8)] # Length 6 for consistency with pre-revision CST+ coloration
names(CSTColors) <- CSTs

CSTColorScale <- scale_colour_manual(name = "CST_pam_3", values = CSTColors[1:8])
CSTFillScale <- scale_fill_manual(name = "CST_pam_3", values = CSTColors[1:8])

CSTColorScale <- scale_colour_manual(name = "CST_pam_4", values =CSTColors[1:8])
CSTFillScale <- scale_fill_manual(name = "CST_pam_4", values = CSTColors[1:8])

CSTColorScale <- scale_colour_manual(name = "CST_pam_5", values = CSTColors[1:8])
CSTFillScale <- scale_fill_manual(name = "CST_pam_5", values = CSTColors[1:8])

CSTColorScale <- scale_colour_manual(name = "CST_pam_6", values = CSTColors[1:8])
CSTFillScale <- scale_fill_manual(name = "CST_pam_6", values = CSTColors[1:8])

CSTColorScale <- scale_colour_manual(name = "CST_pam_7", values = CSTColors[1:8])
CSTFillScale <- scale_fill_manual(name = "CST_pam_7", values = CSTColors[1:8])


PAM_3 <- plot_ordination(phy.obj.rel, ord, color="CST_pam_3") + 
  CSTColorScale + theme_bw()+
  theme(axis.text=element_text(size=17, face="bold"),
        axis.title=element_text(size=17,face="bold"), 
        legend.title=element_text(size=17,face="bold"),
        legend.text = element_text(size =16), legend.position="bottom") +geom_point(size=2)+ 
  guides(color = guide_legend(title = "PAM cluster (k=3)"))  

#ggsave(file = "./figure_output/SuppX.PAM3.pdf", plot = PAM_3, dpi = 500, width = 8, height = 7) # in this way, I can control the size of the figure
dev.off()


plot_ordination(phy.obj.rel, ord, axes=c(3,4), color="CST_pam_3") + CSTColorScale

plot_ordination(phy.obj.rel, ord, color="CST_pam_4") + CSTColorScale

plot_ordination(phy.obj.rel, ord, color="CST_pam_5") + CSTColorScale

PAM_6 <- plot_ordination(phy.obj.rel, ord, color="CST_pam_6") + CSTColorScale+ 
  theme_bw()+ theme(axis.text=element_text(size=17, face="bold"),
                    axis.title=element_text(size=17,face="bold"), 
                    legend.title=element_text(size=17,face="bold"),
                    legend.text = element_text(size =16), legend.position="bottom") +geom_point(size=2)+ 
  guides(color = guide_legend(title = "PAM cluster (k=6)"))  

#ggsave(file = "./figure_output/SuppX.PAM6.pdf", plot = PAM_6, dpi = 500, width = 8, height = 7) # in this way, I can control the size of the figure
dev.off()


PAM_7 <- plot_ordination(phy.obj.rel, ord, color="CST_pam_7") + CSTColorScale+
  theme_bw()+ theme(axis.text=element_text(size=17, face="bold"),
                    axis.title=element_text(size=17,face="bold"), 
                    legend.title=element_text(size=17,face="bold"),
                    legend.text = element_text(size =16), legend.position="bottom") +geom_point(size=2)+ 
  guides(color = guide_legend(title = "PAM cluster (k=7)")) 



ggsave(file = "./figure_output/SuppX.PAM7.pdf", plot = PAM_7, dpi = 500, width = 8, height = 7) # in this way, I can control the size of the figure
dev.off()


# let's add the cluster annotation (k=3)
sample_data(phy.obj.rel)$CST_pam_3 <- clust

sample_data(Phy.f_std)$CST_pam_3 <- clust   # add the annotation to Phy.f_std object too.
levels(sample_data(Phy.f_std)$CST_pam_3) <- c("1"="Cluster1","2"="Cluster2","3"="Cluster3")


## add cluster annotation (k=3) using brey instead 21MAr2023
sample_data(phy.obj.rel)$CST_pam_3_bray <- clust

sample_data(Phy.f_std)$CST_pam_3_bray <- clust   # add the annotation to Phy.f_std object too.
levels(sample_data(Phy.f_std)$CST_pam_3_bray) <- c("1"="Cluster1","2"="Cluster2","3"="Cluster3")

data <- tibble(data.frame(sample_data(Phy.f_std)))
data_2 <- data %>% select(PID, Study_site, Status2,CST_pam_3,CST_pam_3_bray)
View(data_2) # can see some annotation are different whether we use "bray" or "Jensen-shannon" distance




# let's add the cluster annotation (k=6)
sample_data(phy.obj.rel)$CST_pam_6 <- clust

sample_data(Phy.f_std)$CST_pam_6 <- clust   # add the annotation to Phy.f_std object too.
levels(sample_data(Phy.f_std)$CST_pam_6) <- c("1"="Cluster1","2"="Cluster2","3"="Cluster3","4"="Cluster4","5"="Cluster5","6"="Cluster6")


### let's coombine PCoA plots together
PAM_3
PAM_6
PAM_7
library(ggpubr)
PAM_all <- ggarrange(PAM_3, PAM_6, PAM_7 + rremove("x.text"), 
          labels = c(), nrow = 1)

ggsave(file = "./figure_output/SuppX.PAM_combine.pdf", plot = PAM_all, dpi = 900, width = 17, height = 6) # in this way, I can control the size of the figure
dev.off()

# save the PAM cluter info...
saveRDS(Phy.f_std, "Phy.f_std.RDS")


#-----------------------------------------------------------------------------------------------------
#  Sequence Lot  ----
#-----------------------------------------------------------------------------------------------------
### 
### Beta-diversity (11) PCoA plot of SA comparing by sequenceLot - use "jsd"


Phy.f_std_CT <- subset_samples(Phy.f_std, Study_site=="South Africa")
levels(sample_data(Phy.f_std_CT)$SequenceLot) <- c("1"="Batch 1","5"="Batch 5")

set.seed(20) #select any start seed

ord.BC <- ordinate(Phy.f_std_CT, method = "PCoA", distance = "jsd", k=2, trymax=1000) 
color = c("SequenceLot")
shape = c("Visit")
# title=c("PCoA of 16S microbiome (CT samples), Jensen–Shannon distance, k=2")  
MDS = plot_ordination(Phy.f_std_CT, ord.BC, color = color, shape=shape
                      #,title = title
                      )

SequenceLot_SA  = MDS+
  theme_bw()+
  theme(axis.text=element_text(size=17, face="bold"),
                    axis.title=element_text(size=17,face="bold"), 
                    legend.title=element_text(size=17 ,face="bold"),
                    legend.text = element_text(size =16))+ 
  guides(color = guide_legend(title = "Miseq run")) +
  labs(color=color, shape=shape)+geom_point(size=3)

SequenceLot_SA

ggsave(file = "./figure_output/SuppX_PCoA_jsd_k2_SA.batch.pdf", plot = SequenceLot_SA, dpi = 500, width = 8, height = 7) # in this way, I can control the size of the figure
dev.off()




### PERMANOVA 
set.seed(20)
Phy.f_std_CT_W1 <- subset_samples(Phy.f_std, Study_site=="South Africa" & Visit =="Week 1")
diss <- phyloseq::distance(Phy.f_std_CT_W1, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_CT_W1_df <- data.frame(row.names=sample_names(Phy.f_std_CT_W1),sample_data(Phy.f_std_CT_W1)) #create dataframe with metadata
# Whether "SequenceLot" differ significantly from each other
adonis(diss ~ SequenceLot, data=Phy.f_std_CT_W1_df, permutations = 999, by = "terms", na.omit)  # p=0.719


Phy.f_std_CT_W15 <- subset_samples(Phy.f_std, Study_site=="South Africa" & Visit =="Week 15")
diss <- phyloseq::distance(Phy.f_std_CT_W15, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_CT_W15_df <- data.frame(row.names=sample_names(Phy.f_std_CT_W15),sample_data(Phy.f_std_CT_W15)) #create dataframe with metadata
# Whether "SequenceLot" differ significantly from each other
adonis(diss ~ SequenceLot, data=Phy.f_std_CT_W15_df, permutations = 999, by = "terms", na.omit)  # p=0.524



### Beta-diversity (12) PCoA plot of Nigeria comparing by sequenceLot - use "jsd"
Phy.f_std_Nigeria <- subset_samples(Phy.f_std, Study_site=="Nigeria")
levels(sample_data(Phy.f_std_Nigeria)$SequenceLot) <- c("6"="Batch 6","7"="Batch 7")
set.seed(20) #select any start seed

ord.BC <- ordinate(Phy.f_std_Nigeria, method = "PCoA", distance = "jsd", k=2, trymax=1000) 
color = c("SequenceLot")
shape = c("Visit")
# title=c("PCoA of 16S microbiome (Nigeria samples), Jensen–Shannon distance, k=2")  
MDS = plot_ordination(Phy.f_std_Nigeria, ord.BC, color = color, shape=shape #, title = title
                      )
SequenceLot_Nig  = MDS+
  theme_bw()+
  theme(axis.text=element_text(size=17, face="bold"),
        axis.title=element_text(size=17,face="bold"), 
        legend.title=element_text(size=17 ,face="bold"),
        legend.text = element_text(size =16))+ 
  guides(color = guide_legend(title = "Miseq run")) +
  labs(color=color, shape=shape)+geom_point(size=3)

SequenceLot_Nig

ggsave(file = "./figure_output/SuppX_PCoA_jsd_k2_Nigeria.batch.pdf", plot = SequenceLot_Nig, dpi = 500, width = 8, height = 7) # in this way, I can control the size of the figure
dev.off()


### PERMANOVA 
set.seed(20)
Phy.f_std_Nigeria_W1 <- subset_samples(Phy.f_std, Study_site=="Nigeria" & Visit == "Week 1")
diss <- phyloseq::distance(Phy.f_std_Nigeria_W1, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_Nigeria_W1_df <- data.frame(row.names=sample_names(Phy.f_std_Nigeria_W1),sample_data(Phy.f_std_Nigeria_W1)) #create dataframe with metadata
# Whether "SequenceLot" differ significantly from each other
adonis(diss ~ SequenceLot, data=Phy.f_std_Nigeria_W1_df, permutations = 999, by = "terms", na.omit)  # p=0.085

set.seed(20)
Phy.f_std_Nigeria_W15 <- subset_samples(Phy.f_std, Study_site=="Nigeria" & Visit == "Week 15")
diss <- phyloseq::distance(Phy.f_std_Nigeria_W15, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_Nigeria_W15_df <- data.frame(row.names=sample_names(Phy.f_std_Nigeria_W15),sample_data(Phy.f_std_Nigeria_W15)) #create dataframe with metadata
# Whether "SequenceLot" differ significantly from each other
adonis(diss ~ SequenceLot, data=Phy.f_std_Nigeria_W15_df, permutations = 999, by = "terms", na.omit)  # 0.067


# let's plot PcOA plot (sequence Lot) together
library(ggpubr)
SequenceLot_SA
SequenceLot_Nig

PAM_all <- ggarrange(SequenceLot_SA, SequenceLot_Nig + rremove("x.text"), 
                     labels = c(), nrow = 1)

ggsave(file = "./figure_output/SuppX_PCoA_jsd_k2_batch.combine.pdf", plot = PAM_all, dpi = 900, width = 15, height = 7) # in this way, I can control the size of the figure
dev.off()


#-----------------------------------------------------------------------------------------------------
#  Relative abundance by PAM cluster groups (birth)   ----
#-----------------------------------------------------------------------------------------------------
### FigS  (relative abundance)----- related to Fig1 (phylum level compared by PAM cluster)
### Simiar to the figure above but with re-ordering the barplot - 
##BARPLOT - ORDERED AND CLUSTERED
library(maditr)
Phy.f_std <- readRDS("Phy.f_std.RDS")
class(Phy.f_std@sam_data$CST_pam_3)
View(Phy.f_std@tax_table)
Species

barplot_stag_BEAMING_sort <- function(physeq) {
  p1 <- tax_glom(physeq, taxrank = 'Phylum') #agglumerate at species level
  p1_30 = prune_taxa(names(sort(taxa_sums(p1), TRUE))[1:30], p1) #use top 30 most abundant taxa for plotting
  p2 <- transform_sample_counts(p1_30, function(x) x/sum(x)) #get abundance in %
  p3 <- psmelt(p2)  #create dataframe from phyloseq object
  p3$Phylum <- as.character(p3$Phylum) #convert to character
  #sort based on abundance of dominant bacteria in each cluster
  # p4 <- p3[,c("Sample", "Species","CST_pam_3", "Abundance")]
  # p4 <- dcast(p3, Sample + CST_pam_3 ~ Species, value.var = "Abundance", fun.aggregate = sum)
  #p4[p4$CST_pam_3 =="Cluster 1",] <- p4[p4$CST_pam_3=="Cluster 1",][order(p4[p4$CST_pam_3=="Cluster 1","coli"]),]
  # p4[p4$CST_pam_3 =="Cluster 2",] <- p4[p4$CST_pam_3=="Cluster 2",][order(p4[p4$CST_pam_3=="Cluster 2","longum"]),]
  # p4[p4$CST_pam_3 =="Cluster 3",] <- p4[p4$CST_pam_3=="Cluster 3",][order(p4[p4$CST_pam_3=="Cluster 3","faecalis"]),]
  #reorder p3
  p3$Sample <- factor(p3$Sample, levels=c(p4$Sample))
  otus <- otu_table(p1_30)
  tax_table(p1_30)
  #annotation 
  labs <- tax.lab(physeq,otus=otus)
  
  
  #set color palette to accommodate the number of species
  colourCount = length(unique(p3$Phylum))
  getPalette = colorRampPalette(brewer.pal(8, "Accent"))  
  pal<- c("#7FC97F", "#CBB1C3", "#FDDA8D", "#83A3A7", "#D01487", "#8EC293", 
          "#BB5B19","#DAB6B1", "#FEE992","#FDCA89",  "#EC0877","#C9482C",
          "#9DBBA8", "#E9BA9E","#E3EA9C", "#FEF897", "#4B61AA", "#E01D5E", 
          "#ACB5BC", "#90603F","#F8BE8B", "#77479F", "#D43345", "#7B6352",
          "#BBAED1", "#A65E2C", "#B3C7A1", "#A32D93","#5380AC", "#666666")
  
  pal<- c("steelblue3", "#8EC293", "#FEE992", "azure4", "brown3", 
          "#CBB1C3", "tomato2","lightgreen", "cyan","#FDCA89",  
          "#EC0877","sienna1","lightseagreen", "lightcoral","#E3EA9C", 
          "yellow", "blue3", "#E01D5E", "chartreuse3", "chocolate4",
          "darkgoldenrod1", "#77479F", "honeydew3", "yellow4", "violet", 
          "tan3", "seagreen", "#A32D93","lightskyblue", "#666666")
  #annotation 
  #plot
  barplot_species <- ggplot(data=p3, aes(x=Sample, y=Abundance, fill=Phylum)) 
  barplot_species <- barplot_species + geom_bar(aes(), stat="identity", position="stack") + 
    facet_wrap(Study_site ~ CST_pam_3, scales = "free") + guides(fill=guide_legend(nrow=2)) + 
    scale_fill_manual(values=pal) + ylab("Relative abundance") +xlab("Participants") + 
    #ggtitle("Bacterial composition at each time pooint by study site") + 
    theme(legend.position="bottom", strip.background = element_rect(fill="lightgrey", color = "black"), 
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          panel.background=element_blank(), panel.border=element_blank(),
          axis.title.y = element_text(size = 17, colour = "black",face = "bold"),
          axis.title.x = element_blank(),
          strip.text = element_text(size=17, face="bold"), 
          axis.text.y=element_text(size=17, face="bold"), 
          text = element_text(size=15, face="bold"))
  barplot_species
  return(barplot_species)
}


levels(sample_data(Phy.f_std)$Visit) <- c("Birth"="Week 1","W15"="Week 15")
Phy.f_std_birth <- subset_samples(Phy.f_std, Visit == "Week 1")
FigS <- barplot_stag_BEAMING_sort(Phy.f_std_birth)

ggsave(file = "./figure_output/FigS_barplot_birth_PhylumLevel.PAMcluster.pdf", plot = FigS, dpi = 900, width = 20, height = 12) # in this way, I can control the size of the figure
dev.off()



#-----------------------------------------------------------------------------------------------------
#  Relative abundance by PAM cluster groups (W15)   ----
#-----------------------------------------------------------------------------------------------------
### FigS  (relative abundance)----- related to Fig2 (phylum level compared by PAM cluster)
### Simiar to the figure above but with re-ordering the barplot - 
##BARPLOT - ORDERED AND CLUSTERED
library(maditr)
Phy.f_std <- readRDS("Phy.f_std.RDS")
class(Phy.f_std@sam_data$CST_pam_3)
View(Phy.f_std@tax_table)
Species

barplot_stag_BEAMING_sort <- function(physeq) {
  p1 <- tax_glom(physeq, taxrank = 'Phylum') #agglumerate at species level
  p1_30 = prune_taxa(names(sort(taxa_sums(p1), TRUE))[1:30], p1) #use top 30 most abundant taxa for plotting
  p2 <- transform_sample_counts(p1_30, function(x) x/sum(x)) #get abundance in %
  p3 <- psmelt(p2)  #create dataframe from phyloseq object
  p3$Phylum <- as.character(p3$Phylum) #convert to character
  #sort based on abundance of dominant bacteria in each cluster
  # p4 <- p3[,c("Sample", "Species","CST_pam_3", "Abundance")]
  # p4 <- dcast(p3, Sample + CST_pam_3 ~ Species, value.var = "Abundance", fun.aggregate = sum)
  #p4[p4$CST_pam_3 =="Cluster 1",] <- p4[p4$CST_pam_3=="Cluster 1",][order(p4[p4$CST_pam_3=="Cluster 1","coli"]),]
  # p4[p4$CST_pam_3 =="Cluster 2",] <- p4[p4$CST_pam_3=="Cluster 2",][order(p4[p4$CST_pam_3=="Cluster 2","longum"]),]
  # p4[p4$CST_pam_3 =="Cluster 3",] <- p4[p4$CST_pam_3=="Cluster 3",][order(p4[p4$CST_pam_3=="Cluster 3","faecalis"]),]
  #reorder p3
  p3$Sample <- factor(p3$Sample, levels=c(p4$Sample))
  otus <- otu_table(p1_30)
  tax_table(p1_30)
  #annotation 
  labs <- tax.lab(physeq,otus=otus)
  
  
  #set color palette to accommodate the number of species
  colourCount = length(unique(p3$Phylum))
  getPalette = colorRampPalette(brewer.pal(8, "Accent"))  
  pal<- c("#7FC97F", "#CBB1C3", "#FDDA8D", "#83A3A7", "#D01487", "#8EC293", 
          "#BB5B19","#DAB6B1", "#FEE992","#FDCA89",  "#EC0877","#C9482C",
          "#9DBBA8", "#E9BA9E","#E3EA9C", "#FEF897", "#4B61AA", "#E01D5E", 
          "#ACB5BC", "#90603F","#F8BE8B", "#77479F", "#D43345", "#7B6352",
          "#BBAED1", "#A65E2C", "#B3C7A1", "#A32D93","#5380AC", "#666666")
  
  pal<- c("steelblue3", "#8EC293", "#FEE992", "azure4", "brown3", 
          "#CBB1C3", "tomato2","lightgreen", "cyan","#FDCA89",  
          "#EC0877","sienna1","lightseagreen", "lightcoral","#E3EA9C", 
          "yellow", "blue3", "#E01D5E", "chartreuse3", "chocolate4",
          "darkgoldenrod1", "#77479F", "honeydew3", "yellow4", "violet", 
          "tan3", "seagreen", "#A32D93","lightskyblue", "#666666")
  #annotation 
  #plot
  barplot_species <- ggplot(data=p3, aes(x=Sample, y=Abundance, fill=Phylum)) 
  barplot_species <- barplot_species + geom_bar(aes(), stat="identity", position="stack") + 
    facet_wrap(Study_site ~ CST_pam_3, scales = "free") + guides(fill=guide_legend(nrow=2)) + 
    scale_fill_manual(values=pal) + ylab("Relative abundance") +xlab("Participants") + 
    #ggtitle("Bacterial composition at each time pooint by study site") + 
    theme(legend.position="bottom", strip.background = element_rect(fill="lightgrey", color = "black"), 
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          panel.background=element_blank(), panel.border=element_blank(),
          axis.title.y = element_text(size = 17, colour = "black",face = "bold"),
          axis.title.x = element_blank(),
          strip.text = element_text(size=17, face="bold"), 
          axis.text.y=element_text(size=17, face="bold"), 
          text = element_text(size=15, face="bold"))
  barplot_species
  return(barplot_species)
}


levels(sample_data(Phy.f_std)$Visit) <- c("Birth"="Week 1","W15"="Week 15")
Phy.f_std_W15 <- subset_samples(Phy.f_std, Visit == "Week 15")
FigS <- barplot_stag_BEAMING_sort(Phy.f_std_W15)
FigS

ggsave(file = "./figure_output/FigS_barplot_W15_PhylumLevel.PAMcluster.pdf", plot = FigS, dpi = 900, width = 20, height = 12) # in this way, I can control the size of the figure
dev.off()

#-----------------------------------------------------------------------------------------------------
#  Relative abundance with matching PID   ----
#-----------------------------------------------------------------------------------------------------
### Simiar to the figure above but with re-ordering the barplot - Supp.figX; correlate with fig 2 (A)
##BARPLOT - ORDERED AND CLUSTERED
## only samples with both time point, and each colum with the same PID

barplot_stag_BEAMING_sort.modi <- function(physeq) {
  p1 <- tax_glom(physeq, taxrank = 'Species') #agglumerate at species level
  p1_30 = prune_taxa(names(sort(taxa_sums(p1), TRUE))[1:30], p1) #use top 30 most abundant taxa for plotting
  p2 <- transform_sample_counts(p1_30, function(x) x/sum(x)) #get abundance in %
  p3 <- psmelt(p2)  #create dataframe from phyloseq object
  p3$Species <- as.character(p3$Species) #convert to character
  #sort based on abundance of dominant bacteria in each cluster
  p4 <- p3[,c("Sample", "Species","CST_pam_3", "Abundance")]
  p4 <- dcast(p3, Sample + CST_pam_3 ~ Species, value.var = "Abundance", fun.aggregate = sum)
  #p4[p4$CST_pam_3 =="Cluster1",] <- p4[p4$CST_pam_3=="Cluster1",][order(p4[p4$CST_pam_3=="Cluster1","coli"]),]
  #p4[p4$CST_pam_3 =="Cluster2",] <- p4[p4$CST_pam_3=="Cluster2",][order(p4[p4$CST_pam_3=="Cluster2","longum"]),]
  #p4[p4$CST_pam_3 =="Cluster3",] <- p4[p4$CST_pam_3=="Cluster3",][order(p4[p4$CST_pam_3=="Cluster3","faecalis"]),]
  #reorder p3
  p3$Sample <- factor(p3$Sample, levels=c(p4$Sample))
  otus <- otu_table(p1_30)
  #annotation 
  labs <- tax.lab(physeq,otus=otus)
  #set color palette to accommodate the number of species
  colourCount = length(unique(p3$Genus))
  getPalette = colorRampPalette(brewer.pal(8, "Accent")) 
  
  pal<- c("steelblue3", "#8EC293", "#FEE992", "azure4", "brown3", 
          "#CBB1C3", "tomato2","lightgreen", "cyan","#FDCA89",  
          "#EC0877","sienna1","lightseagreen", "lightcoral","#E3EA9C", 
          "yellow", "blue3", "#E01D5E", "chartreuse3", "chocolate4",
          "darkgoldenrod1", "#77479F", "honeydew3", "yellow4", "violet", 
          "tan3", "seagreen", "#A32D93","lightskyblue", "#666666")
  #annotation 
  #plot
  barplot_species <- ggplot(data=p3, aes(x=Sample, y=Abundance, fill=Species)) 
  barplot_species <- barplot_species + geom_bar(aes(), stat="identity", position="stack") + 
    facet_wrap(Visit ~ Study_site, scales = "free")  + guides(fill=guide_legend(nrow=3)) + 
    scale_fill_manual(values=pal) + ylab("Relative abundance") +xlab("Participants") + 
    #ggtitle("Bacterial composition at each time pooint by study site (onyl PIDs with 2TP)")+ 
    theme(legend.position="bottom") + theme(axis.text.x=element_text(size=7, angle=90,hjust=1,vjust=0.2))+
    theme(legend.position="bottom",strip.background = element_rect(fill="lightgrey", color = "black"),
          axis.ticks.x = element_blank(),panel.border=element_blank(),
          strip.text = element_text(size=16, face="bold"), 
          axis.text.y=element_text(size=16, face="bold"), 
          text = element_text(size=16, face="bold"))
  return(barplot_species)
}



library(maditr)
#levels(sample_data(Phy.f_std)$"Study_site") <- c("CT"="South Africa","Nigeria"="Nigeria") # change the site name

df<- sample_data(Phy.f_std)
PID.with.2TP <- df$PID2[duplicated(df$PID2)] # filter only the samples that both 2 TP are available
Phy.f_std_2TP <- subset_samples(Phy.f_std,  PID2 %in% PID.with.2TP)

Fig2_A_supp <- barplot_stag_BEAMING_sort.modi(Phy.f_std_2TP)
Fig2_A_supp


ggsave(file = "./figure_output/SuppX_barplot_matchingPID.pdf", plot = Fig2_A_supp, dpi = 900, width = 20, height = 12) # in this way, I can control the size of the figure
dev.off()



#-----------------------------------------------------------------------------------------------------
# Transition of bacteria (Deseq) over time SA / Nigeria  ----
#-----------------------------------------------------------------------------------------------------
### SA
### Deseq2 result 
Phy.f_std <- readRDS("Phy.f_std.RDS")
Phy.f_std_SA <- subset_samples(Phy.f_std, Study_site == "South Africa")

Phy.f_std_SA_glom <- tax_glom.kv(Phy.f_std_SA)

sample_data(Phy.f_std_SA)$Visit  # Levels: Week 1 Week 15

levels(sample_data(Phy.f_std_SA)$Visit) <- c("Week 1"="W1", "Week 15" = "W15") # [1] "W1"  "W15"



#Convert the phyloseq object to a DESeqDataSet and run DESeq2:
ds = phyloseq_to_deseq2(Phy.f_std_SA, ~ Visit)
ds = DESeq(ds, fitType = "local")  # get an error...rror in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : every gene contains at least one zero, cannot compute log geometric means
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(ds), 1, gm_mean)
ds = estimateSizeFactors(ds, geoMeans = geoMeans)
ds = DESeq(ds, fitType="local")
res = results(ds)
res = res[order(res$padj, na.last = NA), ]
alpha = 0.01 #this is your significance level, adjust if desired
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Phy.f_std_SA)[rownames(sigtab), ], "matrix"))
head(sigtab)
sigtab %>% View()

resultsNames(ds) 
rownames(sigtab) 

#only show bugs annotated at the Species level
#sigtabspec = subset(sigtab, !is.na(Species))
sigtabgen = subset(sigtab, !is.na(Genus))
#sigtabord = subset(sigtab, !is.na(Order))
head(sigtabgen)

#write.csv(sigtabgen, file = "./figure_output/Sfig-sigtabgen_p0.01_Supp_SA.W1vsW15.deseq.csv")

#sort for pretty plotting
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

#Alternative: Phylum level, for plotting purposes - don't do both
#x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
#x = sort(x, TRUE)
#sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

#fix taxrank and subset to highly significant asvs
sigtabgen.hs <- sigtabgen[which(sigtabgen$padj < 0.005),]

#subset to greater magnitude of FC
sigtabgen$Species[is.na(sigtabgen$Species)] = "(unclassified)"
sigtabgen$Species[sigtabgen$Species == "NA"] <- "unknown"
sigtabgen$genspec <- paste(sigtabgen$Genus, sigtabgen$Species, sep = " ")
sigtab.merged.hs <- sigtabgen[which(sigtabgen$padj < 0.005),]

sigtabgen %>%
  subset(log2FoldChange > 0.5) -> sigtabgen.hi
sigtabgen %>%
  subset(log2FoldChange < -0.5) -> sigtabgen.low

sigtab.merged.hs %>%
  subset(log2FoldChange > 0.5) -> sigtabgen.hi

sigtab.merged.hs %>%
  subset(log2FoldChange < -0.5) -> sigtabgen.low

sigtabgen.fc <- rbind.data.frame(sigtabgen.hi, sigtabgen.low)                          
View(sigtabgen.fc)

#plot the results (Genus currently set, and colored by Order, set as desired)
SuppFig.SA.overtime.deseq <- ggplot(sigtabgen.fc, aes(y=Genus, x=log2FoldChange, fill=genspec)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_vline(xintercept = 0.5, color = "indianred3", size = 0.5, lty=2) +
  geom_vline(xintercept = -0.5, color = "indianred3", size = 0.5, lty=2) +
  geom_point(pch = 21, size=3)+ theme_bw()+
  theme(axis.text.y = element_text(size = 14, colour = "black",face = "bold"),
        axis.title.y = element_text(size = 17, colour = "black",face = "bold"),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 17, colour = "black",face = "bold"),
        strip.background =element_rect(fill="white"),
        strip.text.x = element_text(colour = "mediumpurple4", size = 15, hjust = 0.5, vjust = 0.5, face = "bold"),
        legend.text = element_text(size = 14, colour = "black"),
        legend.title = element_text(size=15),
        plot.title = element_text(size = 17, face = "bold")) +  theme(legend.position="right", legend.title=element_blank())+
  guides(fill = guide_legend(nrow = 12)) + guides(fill=guide_legend(ncol = 2))+
  labs(x=expression(log[2]~fold~change))+ ggtitle("South Africa (Week 1 vs Week 15)")  #scale_fill_manual(values=custom_colors_36) + #+ ggtitle("DESeq2 SA vs Nigeria at birth, p=0.01, w/o melting the lowest taxonomy")

SuppFig.SA.overtime.deseq # something wrong with DEseq - it looks like opposite

ggsave(file = "./figure_output/SuppFig.SA.overtime.deseq.pdf", plot = SuppFig.SA.overtime.deseq, dpi = 900, width = 15, height = 7) # in this way, I can control the size of the figure
dev.off()






### Nigeria

Phy.f_std <- readRDS("Phy.f_std.RDS")
Phy.f_std_Nig <- subset_samples(Phy.f_std, Study_site == "Nigeria")

Phy.f_std_Nig_glom <- tax_glom.kv(Phy.f_std_Nig)

sample_data(Phy.f_std_Nig)$Visit  # Levels: Week 1 Week 15

levels(sample_data(Phy.f_std_Nig)$Visit) <- c("Week 1"="W1", "Week 15" = "W15") # [1] "W1"  "W15"



#Convert the phyloseq object to a DESeqDataSet and run DESeq2:
ds = phyloseq_to_deseq2(Phy.f_std_Nig, ~ Visit)
ds = DESeq(ds, fitType = "local")  # get an error...rror in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : every gene contains at least one zero, cannot compute log geometric means
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(ds), 1, gm_mean)
ds = estimateSizeFactors(ds, geoMeans = geoMeans)
ds = DESeq(ds, fitType="local")
res = results(ds)
res = res[order(res$padj, na.last = NA), ]
alpha = 0.01 #this is your significance level, adjust if desired
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Phy.f_std_Nig)[rownames(sigtab), ], "matrix"))
head(sigtab)
sigtab %>% View()

resultsNames(ds) 
rownames(sigtab) 

#only show bugs annotated at the Species level
#sigtabspec = subset(sigtab, !is.na(Species))
sigtabgen = subset(sigtab, !is.na(Genus))
#sigtabord = subset(sigtab, !is.na(Order))
head(sigtabgen)

#write.csv(sigtabgen, file = "./figure_output/Sfig-sigtabge_p0.01_Supp_Nig.W1vsW15.deseq.csv")

#sort for pretty plotting
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

#Alternative: Phylum level, for plotting purposes - don't do both
#x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
#x = sort(x, TRUE)
#sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

#fix taxrank and subset to highly significant asvs
sigtabgen.hs <- sigtabgen[which(sigtabgen$padj < 0.005),]

#subset to greater magnitude of FC
sigtabgen$Species[is.na(sigtabgen$Species)] = "(unclassified)"
sigtabgen$Species[sigtabgen$Species == "NA"] <- "unknown"
sigtabgen$genspec <- paste(sigtabgen$Genus, sigtabgen$Species, sep = " ")
sigtab.merged.hs <- sigtabgen[which(sigtabgen$padj < 0.005),]

sigtabgen %>%
  subset(log2FoldChange > 0.5) -> sigtabgen.hi
sigtabgen %>%
  subset(log2FoldChange < -0.5) -> sigtabgen.low

sigtab.merged.hs %>%
  subset(log2FoldChange > 0.5) -> sigtabgen.hi

sigtab.merged.hs %>%
  subset(log2FoldChange < -0.5) -> sigtabgen.low

sigtabgen.fc <- rbind.data.frame(sigtabgen.hi, sigtabgen.low)                          
View(sigtabgen.fc)

#plot the results (Genus currently set, and colored by Order, set as desired)
SuppFig.Nig.overtime.deseq <- ggplot(sigtabgen.fc, aes(y=Genus, x=log2FoldChange, fill=genspec)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_vline(xintercept = 0.5, color = "indianred3", size = 0.5, lty=2) +
  geom_vline(xintercept = -0.5, color = "indianred3", size = 0.5, lty=2) +
  geom_point(pch = 21, size=3)+ theme_bw()+
  theme(axis.text.y = element_text(size = 14, colour = "black",face = "bold"),
        axis.title.y = element_text(size = 17, colour = "black",face = "bold"),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 17, colour = "black",face = "bold"),
        strip.background =element_rect(fill="white"),
        strip.text.x = element_text(colour = "mediumpurple4", size = 15, hjust = 0.5, vjust = 0.5, face = "bold"),
        legend.text = element_text(size = 14, colour = "black"),
        legend.title = element_text(size=15),
        plot.title = element_text(size = 17, face = "bold")) +  theme(legend.position="right", legend.title=element_blank())+
  guides(fill = guide_legend(nrow = 12)) + guides(fill=guide_legend(ncol = 2))+
  labs(x=expression(log[2]~fold~change)) + ggtitle("Nigeria (Week 1 vs Week 15)") #+ #scale_fill_manual(values=custom_colors_36) + 

SuppFig.Nig.overtime.deseq # something wrong with DEseq - it looks like opposite

ggsave(file = "./figure_output/SuppFig.Nig.overtime.deseq.pdf", plot = SuppFig.Nig.overtime.deseq, dpi = 900, width = 15, height = 7) # in this way, I can control the size of the figure
dev.off()



#-----------------------------------------------------------------------------------------------------
#  Phylum comparison   ----
#-----------------------------------------------------------------------------------------------------
## Phylom comparison - Birth (bar plot)  - Supp figure?
library(ggpubr)
library(rstatix)
Phy.f_std <- readRDS("Phy.f_std.RDS") 
Phy.f_std <- microbiome::transform(Phy.f_std, 'clr')

#levels(sample_data(Phy.f_std)$"Study_site") <- c("CT"="South Africa","Nigeria"="Nigeria") # need to change the name
Phy.f_std_birth <- subset_samples(Phy.f_std, Visit =="Week 1")

#Agglomerate to phylum-level and rename
ps_phylum <- phyloseq::tax_glom(Phy.f_std_birth, "Phylum")
phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::otu_table(ps_phylum)[1:5, 1:5]

#Melt and plot
bxp <- phyloseq::psmelt(ps_phylum) %>%
  ggplot(data = ., aes(x = Study_site, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .1) +
  labs(x = "", y = "Abundance (clr transformed) \n") +
  facet_wrap(~ OTU, scales = "free_y", nrow = 3) + theme_bw() + guides(color = guide_legend(title = "Phylum")) +
  theme(axis.text = element_text(size=15, face="bold"),
       axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=15, face="bold"),
       axis.title = element_text(size=17,face="bold"), 
       legend.position = "none",
       strip.text = element_text(size=17, face="bold"))


# stats test- https://www.datanovia.com/en/blog/how-to-add-p-values-onto-a-grouped-ggplot-using-the-ggpubr-r-package/
stat.test <- phyloseq::psmelt(ps_phylum) %>%
  group_by(OTU) %>% 
  wilcox_test(Abundance ~ Study_site) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

stat.result <- data.frame(stat.test)

dat <- data.frame(stat.result, check.names = T)
print(dat)
write_xlsx(dat,"./figure_output/SuppX_phylum_comparison_Birth.clr_stats.xlsx")


# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_y_position() %>%
  add_x_position()

# Show adjusted p-values and significance levels
# Hide ns (non-significant)
Phylum_W1 <- bxp + stat_pvalue_manual(stat.test,  label = "{p.adj.signif}", 
                               tip.length = 0, hide.ns = TRUE) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(file = "./figure_output/SuppX_phylum_comparison_Birth.clr.pdf", plot = Phylum_W1, dpi = 900, width =14 , height = 9) # in this way, I can control the size of the figure
dev.off()


phyloseq::psmelt(ps_phylum) %>%
  group_by(OTU) %>% summarise(mean = mean(Abundance), n = n())　 %>% View()


### Phylom comparison - W15 (bar plot)  - Supp figure?
library(ggpubr)
library(rstatix)
Phy.f_std <- readRDS("Phy.f_std.RDS")
Phy.f_std <- microbiome::transform(Phy.f_std, 'clr')
#levels(sample_data(Phy.f_std)$"Study_site") <- c("CT"="South Africa","Nigeria"="Nigeria") # need to change the name
Phy.f_std_W15 <- subset_samples(Phy.f_std, Visit =="Week 15")

#Agglomerate to phylum-level and rename
ps_phylum <- phyloseq::tax_glom(Phy.f_std_W15, "Phylum")
phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::otu_table(ps_phylum)[1:5, 1:5]

df <- phyloseq::psmelt(ps_phylum)
df$OTU <- as.factor(df$OTU)
df$Study_site <- as.factor(df$Study_site)
df$Visit
df$Abundance

#Melt and plot - [W15 samples] 
bxp <- phyloseq::psmelt(ps_phylum) %>%
  ggplot(data = ., aes(x = Study_site, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .1) +
  labs(x = "", y = "Abundance (clr transformed)\n") +
  facet_wrap(~ OTU, scales = "free_y", nrow = 3) + theme_bw() + guides(color = guide_legend(title = "Phylum")) +
  theme(axis.text = element_text(size=15, face="bold"),
        axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=15, face="bold"),
        axis.title = element_text(size=17,face="bold"), 
        legend.position = "none",
        strip.text = element_text(size=17, face="bold"))



# stats test- https://www.datanovia.com/en/blog/how-to-add-p-values-onto-a-grouped-ggplot-using-the-ggpubr-r-package/
stat.test <- df %>% #filter(!Phylum %in% c("Chloroflexi","Spirochaetota","Synergistota")) %>% 
  group_by(OTU) %>%   wilcox_test(Abundance ~ Study_site) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

stat.result <- data.frame(stat.test)
dat <- data.frame(stat.result, check.names = T)
print(dat)
write_xlsx(dat,"./figure_output/SuppX_phylum_comparison_W15_clt.stats.xlsx")



stat.test <- stat.test %>%
  add_y_position() %>%
  add_x_position()

# Show adjusted p-values and significance levels
# Hide ns (non-significant)
Phylum_W15 <- bxp + stat_pvalue_manual(stat.test,  label = "{p.adj.signif}", 
                               tip.length = 0, hide.ns = TRUE) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(file = "./figure_output/SuppX_phylum_comparison_W15.clt.pdf", plot = Phylum_W15, dpi = 900, width =14 , height = 9) # in this way, I can control the size of the figure
dev.off()

phyloseq::psmelt(ps_phylum) %>%
  group_by(OTU) %>% summarise(mean = mean(Abundance), n = n())　 %>% View()


### let's combine the phylom comparison
library(ggpubr)
Phylum_W1
Phylum_W15

PAM_all <- ggarrange(Phylum_W1, Phylum_W15 ,
                     labels = c(), ncol = 1
                     )

ggsave(file = "./figure_output/SuppX_phylum_comparison.clt_combine.pdf", plot = PAM_all, dpi = 900, width =14 , height = 22) # in this way, I can control the size of the figure
dev.off()



### Phylum comparison- lower taxomony classification - W15 (bar plot)  - Supp figure?
## feeding duration
library(ggpubr)
library(rstatix)
Phy.f_std <- readRDS("Phy.f_std.RDS") 
Phy.f_std <- microbiome::transform(Phy.f_std, 'clr')

#levels(sample_data(Phy.f_std)$"Study_site") <- c("CT"="South Africa","Nigeria"="Nigeria") # need to change the name
Phy.f_std_W15 <- subset_samples(Phy.f_std, Visit =="Week 15" & Study_site == "South Africa")

Phy.f_std_W15_subset = subset_taxa(Phy.f_std_W15, Order %in% c("Bifidobacteriales", "Lactobacillales", "Bacteroidales",
                                                               "Enterobacterales","Staphylococcales", "Coriobacteriales","Veillonellales-Selenomonadales"))


#Agglomerate to phylum-level and rename
ps_phylum <- phyloseq::tax_glom(Phy.f_std_W15_subset, "Family")
phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Family"]
phyloseq::otu_table(ps_phylum)[1:5, 1:5]

df <- phyloseq::psmelt(ps_phylum)
df$OTU <- as.factor(df$OTU)
df$Study_site <- as.factor(df$Study_site)
df$Visit
df$Abundance
df$ FeedingDuration_cat

#Melt and plot - [W15 samples] 
bxp <- phyloseq::psmelt(ps_phylum) %>% 
  ggplot(data = ., aes(x =  FeedingDuration_cat, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .1) +
  labs(x = "", y = "Abundance (clt)\n") +
  facet_wrap(~ OTU, scales = "free_y", nrow = 3) + theme_bw() + guides(color = guide_legend(title = "Phylum")) +
  theme(axis.text = element_text(size=15, face="bold"),
        axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=15, face="bold"),
        axis.title = element_text(size=17,face="bold"), 
        legend.position = "none",
        strip.text = element_text(size=8, face="bold")) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

bxp

a_my_comparisons <- list(c("upto2M", "morethan2M"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
bxp_stat <- bxp + stat_compare_means(method = "wilcox.test", 
                                              comparisons = a_my_comparisons, 
                                              label = "p.signif", symnum.args = symnum.args) 
bxp_stat 





#-----------------------------------------------------------------------------------------------------
#  PCoA - actual age ()   ----
#-----------------------------------------------------------------------------------------------------
names(tibble(data.frame(sample_data(Phy.f_std))))

# N0348B V5 date - need to check.
# N0348B DOB was 2014/01/13
# N0348B Birth visit 2014/01/25 (gernder 1)
# N0348B V5 visit 2014/01/23  (gernder 0)


# consider to create a new annotation for visits (Visit3)...
data <- tibble(data.frame(sample_data(Phy.f_std))) %>%  
  select(PID, Visit, Visit2, Infant_DOB, DateofVisit) %>%
  mutate(age_days = as.Date(DateofVisit) - as.Date(Infant_DOB)) %>% 
  mutate(Visit3 = case_when(age_days <= 3 ~ "Baseline", 
                            age_days >= 4 & age_days <= 20 ~ "Day4-7", 
                            Visit == "W15" ~ "Week 15" 
  )) %>% 
  mutate(Visit4 = case_when(age_days <= 2 ~ "Day0_2", 
                            age_days >= 3 & age_days <= 8 ~ "Day3_8", 
                            age_days >= 9 & age_days <= 20 ~ "Day9_above", 
                            Visit == "Week 15" & age_days <= 111 ~ "Week15",  
                            Visit == "Week 15" & age_days >= 112 ~ "Week16_above"
  ))

sample_data(Phy.f_std)$Visit3 <- data$Visit3
sample_data(Phy.f_std)$Visit4 <- data$Visit4
sample_data(Phy.f_std)$age_days <- as.numeric(data$age_days, units="days") # otherwise it won't allow me to plot on the ggplot2.

data


df <- data.frame(data, check.names = T)
print(df)
write_xlsx(df,"/Users/saoriiwase/Desktop/BEAMING_actual.days.since.DOB.xlsx")



### beta-diversity with actual days (gradient) in the plot - all
to_remove <- c("N0348B")  # this PID's date of collection is wrong
Phy.f_std_pruned <- prune_samples(!(sample_data(Phy.f_std)$PID %in% to_remove), Phy.f_std)

# levels(sample_data(Phy.f_std_pruned)$"Study_site") <- c("CT"="South Africa","Nigeria"="Nigeria")
ord.BC <- ordinate(Phy.f_std_pruned, method = "PCoA", distance = "jsd", k=3, trymax=1000) # try PcOA as NMDS does not work
color = c("age_days")
shape = c("Study_site")
MDS = plot_ordination(Phy.f_std_pruned, ord.BC, color = color, shape = shape)
Fig2_C_supp  = MDS+theme_bw() +
  scale_color_viridis_c(option = "plasma")+
  theme(axis.text=element_text(size=17, face="bold"),
        axis.title=element_text(size=17,face="bold"), 
        legend.title=element_text(size=17 ,face="bold"),
        legend.text = element_text(size =16))+
  labs(color= "Age (days)")+geom_point(size=3)+
  guides(shape = guide_legend(title = "Study site")) # change the shape title



ggsave(file = "./figure_output/SuppX_PCoA_jsd_k2_age.pdf", plot = Fig2_C_supp, dpi = 500, width = 9, height = 7) # in this way, I can control the size of the figure
dev.off()



### beta-diversity with actual days (gradient) in the plot - all & with joinning line
to_remove <- c("N0348B")  # this PID's date of collection is wrong
Phy.f_std_pruned <- prune_samples(!(sample_data(Phy.f_std)$PID %in% to_remove), Phy.f_std)

#levels(sample_data(Phy.f_std_pruned)$"Study_site") <- c("CT"="South Africa","Nigeria"="Nigeria")
ord.BC <- ordinate(Phy.f_std_pruned, method = "PCoA", distance = "jsd", k=3, trymax=1000) # try PcOA as NMDS does not work
color = c("age_days")
shape = c("CST_pam_3")
#title= c("PCoA of 16S microbiome at birth (Jensen-Shannon distance; k=3)")  
MDS = plot_ordination(Phy.f_std_pruned, ord.BC, color = color, shape = shape#, title = title
)
Fig2_C_supp_p2 = MDS+theme_bw() +
  scale_color_viridis_c(option = "plasma")+
  theme(axis.text=element_text(size=17, face="bold"),
        axis.title=element_text(size=17,face="bold"), 
        legend.title=element_text(size=17, face="bold"),
        legend.text = element_text(size =16),
        strip.text = element_text(size=17, face="bold"))+
  labs(color= "Age (days)")+geom_point(size=3)+ facet_wrap(~Study_site)+
  guides(shape = guide_legend(title = "Community cluster"))+ # change the shape title
  geom_line(aes(group=PID2),size = 0.1, color = "grey",arrow=arrow(length=unit(0.3,"cm"))) 

ggsave(file = "./figure_output/SuppX_PCoA_jsd_k2_age_p2.pdf", plot = Fig2_C_supp_p2, dpi = 800, width = 14, height = 7) # in this way, I can control the size of the figure
dev.off()

#AEDBCE


### let's combine the "actual age on PCoA plot"
library(ggpubr)
Fig2_C_supp
Fig2_C_supp_p2

PAM_all <- ggarrange(Fig2_C_supp, Fig2_C_supp_p2 + rremove("x.text"), 
                     labels = c(), nrow = 1, widths = c(0.4, 0.6), heights = c(1, 1))


ggsave(file = "./figure_output/SuppX_PCoA_jsd_k2_age_combine.pdf", plot = PAM_all, dpi = 900, width =17 , height = 5) # in this way, I can control the size of the figure
dev.off()




#-----------------------------------------------------------------------------------------------------
#  alpha diversity - actual age ()   ----
#-----------------------------------------------------------------------------------------------------
data <- tibble(data.frame(sample_data(Phy.f_std))) %>%  
  select(PID, Visit, Visit2, Infant_DOB, DateofVisit) %>%
  mutate(age_days = as.Date(DateofVisit) - as.Date(Infant_DOB)) %>% 
  mutate(Visit5 = case_when(age_days <= 3 ~ "D0-3", 
                            age_days >= 4 & age_days <= 20~ "> D3", 
                            Visit == "Week 15" & age_days <= 105 ~ "D105",  
                            Visit == "Week 15" & age_days >= 106 ~ "> D105"
  )) %>% select(PID, age_days,Visit, Visit5) 

sample_data(Phy.f_std)$Visit5 <- data$Visit5


###### alpha diversity - 
levels(sample_data(Phy.f_std)$"Study_site") <- c("CT"="South Africa","Nigeria"="Nigeria") # need to change the name

sample_data(Phy.f_std)$Visit5 <- as.factor(sample_data(Phy.f_std)$Visit5)
sample_data(Phy.f_std)$Visit5 <-factor(sample_data(Phy.f_std)$Visit5, levels = c("D0-3","> D3", "D105","> D105"))
# Levels: Day0_3 Above_Day3 Week15 Above_W15

p1 <- plot_richness(Phy.f_std, x="Visit5", measures="Shannon",
                    #title = "Alpha Diversity transition over time"
)+
  geom_boxplot(alpha=0.6)+ 
  facet_wrap(~Study_site)+theme_bw()+ xlab("Collection time (days since birth)")+
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
        strip.text = element_text(size=17, face="bold"), axis.text.y=element_text(size=16, face="bold"),text = element_text(size=17, face="bold"))

# and add statistics
# https://rpubs.com/lconteville/713954
a_my_comparisons <- list( c("D0-3", "> D3"), c("> D3", "D105"),c("D105", "> D105"),c("D0-3", "D105"),c("> D3", "> D105"),c("D0-3", "> D105"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
FigS <- p1 + stat_compare_means(method = "wilcox.test", hide.ns = FALSE,
                        comparisons = a_my_comparisons, 
                        label = "p.signif", symnum.args = symnum.args) 



FigS
ggsave(file = "./figure_output/SuppX_alpha_actual.Days.pdf", plot = FigS, dpi = 700, width = 7, height = 7) # in this way, I can control the size of the figure
dev.off()


### Compare the meconium samples vs stool at W1
data <- tibble(data.frame(sample_data(Phy.f_std))) %>%  
  select(PID, Visit, Visit2, Infant_DOB, DateofVisit) %>%
  mutate(age_days = as.Date(DateofVisit) - as.Date(Infant_DOB)) %>% 
  mutate(Visit5 = case_when(age_days <= 3 ~ "Meconium", 
                            age_days >= 4 & age_days <= 20~ "W1 stool", 
                            Visit == "Week 15" & age_days <= 105 ~ "D105",  
                            Visit == "Week 15" & age_days >= 106 ~ "> D105"
  )) %>% select(PID, age_days,Visit, Visit5) 

sample_data(Phy.f_std)$Visit5 <- data$Visit5
sample_data(Phy.f_std)$Visit5 <- as.factor(sample_data(Phy.f_std)$Visit5)

Phy.f_std_W1 <- subset_samples(Phy.f_std, Visit5 == c("Meconium"))

p1 <- plot_richness(Phy.f_std_W1, x="Study_site", measures="Shannon", color = "Study_site")+ 
  geom_boxplot(alpha=0.6)+theme_bw()+ xlab("Study site")+
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
        strip.text = element_blank(), axis.text.y=element_text(size=16, face="bold"),text = element_text(size=17, face="bold"))+  labs(y= "Shannon Index")


a_my_comparisons <- list( c("South Africa", "Nigeria"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
FigS <- p1 + stat_compare_means(method = "wilcox.test", hide.ns = FALSE,
                                comparisons = a_my_comparisons, 
                                label = "p.signif", symnum.args = symnum.args)+  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) 

FigS
ggsave(file = "./figure_output/SuppX_alpha_meconium.color.pdf", plot = FigS, dpi = 700, width = 6, height = 7) # in this way, I can control the size of the figure
dev.off()





#-----------------------------------------------------------------------------------------------------
#  Alluvial plots by status2 (supp fig)  ----
#-----------------------------------------------------------------------------------------------------
#### Alluvial plots
#library(msm)
#install.packages("msm")
Phy.f_std <- readRDS("Phy.f_std.RDS")
# sample_data(Phy.f_std)$DeliverVag = as.factor(sample_data(Phy.f_std)$DeliverVag)


levels(sample_data(Phy.f_std)$Visit) 
levels(sample_data(Phy.f_std)$CST_pam_3_bray) 
# sample_data(Phy.f_std)$Visit <- relevel(sample_data(Phy.f_std)$Visit, ref = "Birth")

sample_data(Phy.f_std)

df<- sample_data(Phy.f_std)
PID.with.2TP <- df$PID2[duplicated(df$PID2)]
physeq <- subset_samples(Phy.f_std, PID2 %in% PID.with.2TP) # only sample with 2TP

df = sample_data(physeq)[,c("Visit","CST_pam_3_bray","PID", "Status2", "Study_site")]
df$CST_pam_3_bray = as.character(unlist(df$CST_pam_3_bray))

df$CST_pam_3_bray[df$CST_pam_3_bray=="Cluster 1"] <- "1" #Change like this because the states apparently are required to be 1, 2, 3 as oppoased to C1, C2, C3 (not sure if they're assumed to be rank ordered, i.e. 1 < 2 <3)
df$CST_pam_3_bray[df$CST_pam_3_bray=="Cluster 2"] <- "2"
df$CST_pam_3_bray[df$CST_pam_3_bray=="Cluster 3"] <- "3"
df$CST_pam_3_bray = as.factor(df$CST_pam_3_bray)



statetable.msm(CST_pam_3_bray, PID, data=df) # transition of the cluster from one to another
## to
## from  1  2  3
## 1 34 23  1
## 2  2 17  0
## 3  9 72  5

#install.packages("ggalluvial")
library(ggalluvial)
## Alluvial Plots (2) - transition of PAM clustering color by Status - maybe add as supp fig?
head(df)
df$CST_pam_3_bray <- forcats::fct_recode(df$CST_pam_3_bray, "Cluster 1 " = "1", "Cluster 2 " = "2","Cluster 3 " = "3")


Fig3_supp_Alluvial <- ggplot(data = as.data.frame(df),
       aes(x=Visit, stratum=CST_pam_3_bray, alluvium=PID)) +
  facet_wrap(~Study_site) + 
  geom_alluvium(aes(fill = Status2),width = 1/16) +
  scale_fill_manual(values=c("#B05A7A","#F99417"))+
  geom_stratum(width = 1/2.5) +
  geom_text(stat = "stratum", size = 5,
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Week 1", "Week 15"),
                   expand = c(0.15, 0.05)) +
  #scale_fill_brewer(type = "qual", palette = "Set1") +
  theme_bw()+
  guides(fill=guide_legend(title="Exposure status"))+
    theme(axis.text.y=element_text(size=17, face="bold"), 
          axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"),
          axis.title=element_text(size=18,face="bold"), 
          legend.title=element_text(size=17, face="bold"),
          legend.text = element_text(size =17),
          strip.text = element_text(size=18, face="bold"))+
  labs(y= "Number of infants")

Fig3_supp_Alluvial

ggsave(file = "./figure_output/SuppX_AlluvialPlots_PAMclustering_Status.bray.color.v2.pdf", plot = Fig3_supp_Alluvial, dpi = 800, width = 11, height = 11) 

dev.off()


#-----------------------------------------------------------------------------------------------------
#  Tetanus titer supp (supp fig)  ----
#-----------------------------------------------------------------------------------------------------
Tetanus <- read_excel("/Users/saoriiwase/Desktop/BEAMING study/BEAMING_microbiome_data_analysis/for paper/Tetanus data/*Tetanus IgG ELISA data_5_3_27 CLEAN_27Jan.xlsx", 
                      sheet = "Tetanus analysis_long")

Tetanus # A tibble: 602 × 7

# dealing with data types
str(Tetanus) # Baby_titer  and Mum's titer needs to be number. No < or  > sign in the colum (substitute)
Tetanus$Baby_titer[which(Tetanus$Baby_titer == "<0.1")] <- "0.099"  # replace the <0.1
Tetanus$Baby_titer <- as.numeric(Tetanus$Baby_titer)  # change it to numeric
Tetanus$Baby_titer <- round(Tetanus$Baby_titer,3)  # change the number of digit

Tetanus$Status2 <- as.factor(Tetanus$Status2)

Tetanus2 <- Tetanus %>% filter(PID %in% stool_PID_info) # select only PIDs that are involved in the stool samples. # A tibble: 578 × 7

# extract only birth data
Tetanus_birth <- Tetanus2 %>% filter(Visit == "Birth") # A tibble: 289 × 7
# extract only W15 data
Tetanus_W15 <- Tetanus2 %>% filter(Visit == "W15") # A tibble: 289 × 7

# extract mum's data - there is no mum's titer data from CT.
Tetanus_mum <-  Tetanus2 %>% filter(Visit == "Birth" & Mum_titer_birth >0) # A tibble: 199 × 7



### Fig S?-A (mum's titer comparison w/ w/o HIV) -------
### barplot for Mum (Nigeria site)
levels(Tetanus_mum$Status2) <- c("HEU"="w/ HIV","HU"="w/o HIV") 

p1 <- ggplot(data=Tetanus_mum, aes(x=Status2, y=Mum_titer_birth, color = Status2, fill = Status2))+
  geom_boxplot(alpha = 0.5) + theme_bw()+ geom_jitter(width=0.3, alpha=0.2)+ #facet_wrap(~Study_site,scales = "free_y")+
  theme(axis.text=element_text(size=17, face="bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"),
        axis.title=element_text(size=17,face="bold"),
        strip.text = element_text(size=17,face="bold"),
        legend.position = "none") +
  scale_color_manual(values=c("#B05A7A","#F99417"))+
  scale_fill_manual(values=c("#B05A7A","#F99417"))+

  labs(x= "Mother's HIV status") +
  labs(y= "Maternal anti-TT IgG (IU/ml)")  #+ ggtitle("Infant mothers' Tetanus titer comparison (only Nigeria)")


# add stats
a_my_comparisons <- list( c("w/ HIV", "w/o HIV"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
Fig4_A_supp.mum <- p1 + stat_compare_means(method = "wilcox.test", 
                                comparisons = a_my_comparisons, 
                                label = "p.signif", symnum.args = symnum.args) 
Fig4_A_supp.mum


compare_means(Mum_titer_birth ~ Status2, Tetanus_mum, method = "wilcox.test", 
              paired = FALSE, ref.group = NULL)

compare_means(Mum_titer_birth ~ Status2, Tetanus_mum, method = "wilcox.test", 
              paired = FALSE,
              group.by = NULL, ref.group = NULL,
              p.adjust.method = "BH") # still not significant even after the multiple comparison adjustment by BH


#ggsave(file = "./figure_output/SuppX_Tetatus_Mum.color.pdf", plot = Fig4_A_supp.mum, dpi = 600, width =5 , height = 6) # in this way, I can control the size of the figure
dev.off()



#### Fig S?-B plot infant titer by overtime & Site (longituidinal comparison) -------
# Supp Figure X

# https://www.datanovia.com/en/blog/how-to-add-p-values-onto-a-grouped-ggplot-using-the-ggpubr-r-package/
meta <- meta %>% filter(Baby_titer >0)
levels(meta$Visit) <- c("Birth" = "Week 1", "W15" = "Week 15")

bxp <- ggboxplot(
  meta, x = "Study_site", y = "Baby_titer", 
  color = "Visit",fill = "Visit",alpha = 0.7, palette = c("#4FC3F7", "#48648F"))  # palette = c("#B3E5BE", "#7DB9B6")

bxp # boxplot


TestFrag1 <- compare_means(Baby_titer ~ Visit, meta, method = "wilcox.test", paired = FALSE,
                           group.by = "Study_site", ref.group = NULL,
                           p.adjust.method = "BH")
# need to add info about position into the stat result
TestFrag1$y.position <- c(6.34, 6.44) 
TestFrag1$x <- c(1, 2) 
TestFrag1$xmin <- c(0.8,1.8)
TestFrag1$xmax <- c(1.2,2.2)

# add the TestFrag1 to the boxplot - we get "bxp.complex"
bxp.complex <- bxp + stat_pvalue_manual(
  TestFrag1,  label = "p.adj", tip.length = 0.01,label.size = 5)


TestFrag2 <- compare_means(Baby_titer ~ Study_site, meta, method = "wilcox.test", paired = FALSE,
                           group.by = "Visit", ref.group = NULL,
                           p.adjust.method = "BH")

# need to add info about position into the stat result
TestFrag2$y.position <- c(7.2, 6.5) 
TestFrag2$xmin <- c(1.2,0.8)
TestFrag2$xmax <- c(2.2,1.8)

# add stat2 on the boxplot - we get "bxp.complex2"
bxp.complex2 <- bxp.complex + 
  stat_pvalue_manual(
    TestFrag2,  label = "p.adj", tip.length = 0.01,
    bracket.nudge.y = 0.5, label.size = 5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Display the plot
Tetanus_supp.mum <- bxp.complex2 + theme_bw()+
  theme(axis.title=element_text(size=17,face="bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
        strip.text = element_text(size=17, face="bold"), 
        axis.text.y=element_text(size=17, face="bold"),
        legend.title=element_text(size=17, face="bold"),
        legend.text = element_text(size =17),
        text = element_text(size=16))+
  labs(x= "Study site") +labs(y= "Anti-TT IgG (IU/ml)") 

Tetanus_supp.mum

ggsave(file = "./figure_output/SuppX_Tetatus_boxplot.W1vsW15.color.BHadj.v2.pdf", plot = Tetanus_supp.mum, dpi = 600, width =6 , height = 6) # in this way, I can control the size of the figure
dev.off()



#### Fig S?-C scatter plot W15 titer vs mum titer (Nigeria)-----
# Supp Figure X

# plot comparing Baby W15 titer & Mum's titer
# use "Tetanus_W15"
View(Tetanus_W15)
Tetanus_W15v2 <- Tetanus_W15 %>% filter(Mum_titer_birth > 0 & Study_site == "Nigeria")
levels(Tetanus_W15v2$Status2) <- c("HEU" = "iHEU", "HU" = "iHUU")

Fig4_A_supp <-ggscatter(data=Tetanus_W15v2, x = "Baby_titer", y = "Mum_titer_birth",   
                add = "reg.line", conf.int = TRUE, size=2, color = "Status2") +
  stat_cor(method = "spearman", label.x = 2.5, label.y = 5.3, size = 5)+ 
  facet_wrap(~Status2)+ 
  theme(axis.text=element_text(size=17,face = "bold"),
        axis.title=element_text(size=17,face="bold"), 
        strip.text = element_text(size=17, face="bold"),
        legend.title=element_blank(),
        legend.text=element_blank()) +
  labs(x= "Anti-TT IgG (IU/ml) [infant]") +
  labs(y= "Anti-TT IgG (IU/ml) [mother]") +
  scale_color_manual(values=c("#B05A7A","#F99417"))+
    scale_fill_manual(values=c("#B05A7A","#F99417"))+
    labs(color= "HIV exposure status")+
    labs(fill= "HIV exposure status")
  # ggtitle("Baby W15 vs mum TT titer comparison (Nigeria)")

ggsave(file = "./figure_output/SuppX_correlation_plot_BabyW15vsMum.spearman.color.pdf", plot = Fig4_A_supp, dpi =700, width =8 , height = 6) # in this way, I can control the size of the figure
dev.off()


Tetanus_W15v2_iHEU <- Tetanus_W15v2 %>% filter(Status2 == "iHEU")
Tetanus_W15v2_iHUU <- Tetanus_W15v2 %>% filter(Status2 == "iHUU")


res1 <-cor.test(Tetanus_W15v2_iHEU$Mum_titer_birth, Tetanus_W15v2_iHEU$Baby_titer,  method = "spearman")
res1  # rho =0.04513498, p-value 0.5991
# Extract the p-value
p_value1 <- res1$p.value  # 0.5991237
p.adjust(p_value1, method = "BH")
p.adjust(p_value1, method = "bonferroni")


res2 <-cor.test(Tetanus_W15v2_iHUU$Mum_titer_birth, Tetanus_W15v2_iHUU$Baby_titer,  method = "spearman")
res2 # rho = 0.3144678, p-value 0.02183
# Extract the p-value
p_value2 <- res2$p.value  # 0.02183144
p.adjust(p_value2, method = "BH")
p.adjust(p_value2, method = "bonferroni")




#### Fig S?-D scatter plot W1vs W15 titer -----
# Supp Figure X

# Plot comparing birth vs W15 titer - Supp figure
Tetanus_birth_data <- Tetanus_birth %>% select(PID, Baby_titer) %>% mutate(Baby_Birthtiter = Baby_titer) %>% select(-Baby_titer)

Tetanus_W15_data <- Tetanus_W15 %>% select(PID, Status2, Study_site,Baby_titer) %>% mutate(Baby_W15titer= Baby_titer) %>% select(-Baby_titer)

data <- left_join(Tetanus_birth_data, Tetanus_W15_data,Study_site, by ="PID")  # A tibble: 289 × 3
View(data)


levels(data$Status2) <- c("HEU" = "iHEU", "HU" = "iHUU")
data$Study_site <- as.factor(data$Study_site)
data$Study_site<-factor(data$Study_site, levels = c("South Africa", "Nigeria"))


Tetanus_supp.W1W15.Status <- ggscatter(data=data, x = "Baby_Birthtiter", y = "Baby_W15titer",   
                add = "reg.line", conf.int = TRUE, size=3) +
  stat_cor(method = "pearson", label.x = 3.5, label.y = 6.2, size = 5)+ 
  facet_wrap(Study_site~Status2)+ 
  theme(axis.text.x=element_text(angle=0,hjust=1,vjust=1,size=17,face = "bold"),
        axis.title=element_text(size=17,face="bold"), strip.text = element_text(size=17,face="bold")) +
  labs(x= "Anti-TT IgG at Week 1 (IU/ml)") +labs(y= "Anti-TT IgG at Week 15 (IU/ml)") #+
  # ggtitle("Baby birth titer vs Baby W15 titer (Nigeria)")

ggsave(file = "./figure_output/SuppX_correlation_plot_W1vsW15.pdf", plot = Tetanus_supp.W1W15.Status, dpi = 800, width =8 , height = 8) # in this way, I can control the size of the figure
dev.off()




#### Fig S?-E scatter plot W1vs W15 titer (both sites combined)-----
# Supp Figure X
# box plot HEU vs HU - SuppX
levels(meta$Study_site) <- c("CT"="South Africa","Nigeria"="Nigeria") 

p1 <- ggplot(data=meta, aes(x=Status2, y=Baby_titer))+ theme_bw()+
  geom_boxplot(aes(color = Status2, fill = Status2) ,alpha = 0.5) +geom_jitter(aes(color = Status2, fill = Status2,alpha = 0.5),width=0.3, alpha=0.2)+ facet_wrap(~Visit,scales = "free_y")+
  theme(axis.title=element_text(size=17,face="bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17,face ="bold"),
        axis.text.y=element_text(size=17,face ="bold"),
        strip.text = element_text(size=17,face="bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=16))  +
  labs(x= "HIV exposure status") +labs(y= "Anti-TT IgG (IU/ml)") + 
  labs(color= "HIV exposure status") +  labs(fill= "HIV exposure status") + 
  # ggtitle("Baby Tetanus titer HEU vs HUU") + 
  ylim(0, 6.7) + 
  scale_color_manual(values=c("#B05A7A","#F99417"))+
  scale_fill_manual(values=c("#B05A7A","#F99417"))

# adjust for multiple comparison...
anno_df <- compare_means(Baby_titer ~ Status2, meta, method = "wilcox.test", paired = FALSE,
                         group.by = "Visit", ref.group = NULL,
                         p.adjust.method = "BH")
anno_df$y.position <- c(6, 6) 
anno_df$xmin <- c(1,1)
anno_df$xmax <- c(2,2)

Tetanus_supp.HEUvsHU.allcombine <- p1 + stat_pvalue_manual(
  anno_df,  label = "p.signif", tip.length = 0.01,
  bracket.nudge.y = 0.5, label.size = 5)



Tetanus_supp.HEUvsHU.allcombine
ggsave(file = "./figure_output/SuppX_Tetanus_barplot_HEUvsHU_BothSites.combined.color.BH.pdf", plot = Tetanus_supp.HEUvsHU.allcombine, dpi = 100, width = 7, height = 7) # in this way, I can control the size of the figure
dev.off()



# find out the median value
median(meta$Baby_titer[meta$Visit == "Week 1" & meta$Status2 == "iHEU"])  # 1.21 # iHEU W1
median(meta$Baby_titer[meta$Visit == "Week 1" & meta$Status2 == "iHUU"])  # 2.2 # iHUU W1


median(meta$Baby_titer[meta$Visit == "Week 15" & meta$Status2 == "iHEU"])  # 1.3 # iHEU W15
median(meta$Baby_titer[meta$Visit == "Week 15" & meta$Status2 == "iHUU"])  # 2.5 # iHUU W15

median(meta$Baby_titer[meta$Visit == "Week 1" & meta$Study_site == "South Africa"])  # 1 
median(meta$Baby_titer[meta$Visit == "Week 1" & meta$Study_site == "Nigeria"]) #1.5


median(meta$Baby_titer[meta$Visit == "Week 15" & meta$Study_site == "South Africa"])  # 1.9 
median(meta$Baby_titer[meta$Visit == "Week 15" & meta$Study_site == "Nigeria"]) #1.6




### Tetanus titer vs mode of feeding
sample_data(Phy.f_std)$Feeding_W15_strictver
sample_data(Phy.f_std)$Feeding_W15_lenientver
sample_data(Phy.f_std)$FeedingDuration_cat
meta_SA <- tibble(data.frame(sample_data(Phy.f_std))) %>% filter(Study_site == "South Africa" & Baby_titer >= 0)
meta_SA$Baby_titer <- as.numeric(meta_SA$Baby_titer)

p1 <- ggplot(data=meta_SA, aes(x=Feeding_W15_lenientver, y=Baby_titer))+ theme_bw()+
  geom_boxplot() +geom_jitter(width=0.3, alpha=0.2)+ facet_wrap(Study_site ~Visit, scales = "free_y")+
  theme(axis.title=element_text(size=17,face="bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17,face ="bold"),
        strip.text = element_text(size=17,face="bold"))  +
  labs(x= "Mode of feeding during the study") +labs(y= "Anti-TT IgG (IU/ml)") + 
  # ggtitle("Baby Tetanus titer HEU vs HUU") + 
  guides(color = guide_legend(title = "Study site"))  + 
  ylim(0, 7.0) + theme(legend.title = element_text(face = "bold"))

# add stats
a_my_comparisons <- list( c("Exclusive_breastfeeding", "Mixed_feeding"))
a_my_comparisons <- list(c("Exclusive_breastfeeding", "Switched_to_formula"))
a_my_comparisons <- list(c("upto2M", "morethan2M"))


symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
Tetanus.vs.feeding <- p1 + stat_compare_means(method = "wilcox.test", 
                                                           comparisons = a_my_comparisons, 
                                                           label = "p.signif", symnum.args = symnum.args) 
Tetanus.vs.feeding # mixed-feeding group had a higher antibody titer?


# try with scatter plot (Tetanus vs duration of feeding)
meta_SA <- tibble(data.frame(sample_data(Phy.f_std))) %>% filter(Study_site == "South Africa" & Baby_titer >= 0)
meta_SA$Baby_titer <- as.numeric(meta_SA$Baby_titer)

ggscatter(data=meta_SA, x = "Feeding_Duration_days", y = "Baby_titer",   
                    add = "reg.line", conf.int = TRUE, size=2) +
  stat_cor(method = "pearson", label.x = 3, label.y = 6.2, size = 5.5)+ 
  facet_wrap(Study_site~Visit)+ theme(axis.text=element_text(size=17, face="bold"),
                              axis.title=element_text(size=17,face="bold"),
                              strip.text = element_text(size=17,face="bold"))+
  labs(x= "Duration of feeding") +
  labs(y= "Anti-TT IgG (IU/ml) [infant at W15]") 


# scatter plot (Tetanus vs gestational age)
meta <- tibble(data.frame(sample_data(Phy.f_std))) %>% filter(Baby_titer >= 0)
names(meta)
meta$Ges_age <- as.numeric(meta$Ges_age)
meta$Baby_titer <- as.numeric(meta$Baby_titer)

ggscatter(data=meta, x = "Ges_age", y = "Baby_titer",   
          add = "reg.line", conf.int = TRUE, size=2) +
  stat_cor(method = "pearson"#, label.x = 3, label.y = 6.2, size = 5.5
           )+ 
  facet_wrap(Study_site~Visit)+ theme(axis.text=element_text(size=17, face="bold"),
                                      axis.title=element_text(size=17,face="bold"),
                                      strip.text = element_text(size=17,face="bold"))+
  labs(x= "Gestational age") +
  labs(y= "Anti-TT IgG (IU/ml) [infant at W15]") 


# scatter plot (Tetanus vs growth score)
meta <- tibble(data.frame(sample_data(Phy.f_std))) %>% filter(Baby_titer >= 0)
names(meta)
meta$Ges_age <- as.numeric(meta$Ges_age)
meta$Baby_titer <- as.numeric(meta$Baby_titer)
meta$wfaz <- as.numeric(meta$wfaz)

ggscatter(data=meta, x = "wfaz", y = "Baby_titer",   
          add = "reg.line", conf.int = TRUE, size=2) +
  stat_cor(method = "pearson"#, label.x = 3, label.y = 6.2, size = 5.5
  )+ 
  facet_wrap(Study_site~Visit)+ theme(axis.text=element_text(size=17, face="bold"),
                                      axis.title=element_text(size=17,face="bold"),
                                      strip.text = element_text(size=17,face="bold"))+
  labs(x= " wfaz") +
  labs(y= "Anti-TT IgG (IU/ml)") 


# box plot (Tetanus vs mode of delivery)
meta <- tibble(data.frame(sample_data(Phy.f_std))) %>% filter(Baby_titer >= 0)
names(meta)
meta$Ges_age <- as.numeric(meta$Ges_age)
meta$Baby_titer <- as.numeric(meta$Baby_titer)
meta$wfaz <- as.numeric(meta$wfaz)
meta$DeliverVag

p1 <- ggplot(data=meta, aes(x=DeliverVag, y=Baby_titer))+ theme_bw()+
  geom_boxplot() +geom_jitter(width=0.3, alpha=0.2)+ facet_wrap(Study_site ~Visit, scales = "free_y")+
  theme(axis.title=element_text(size=17,face="bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17,face ="bold"),
        strip.text = element_text(size=17,face="bold"))  +
  labs(x= "Mode of delivery") +labs(y= "Anti-TT IgG (IU/ml)") + 
  # ggtitle("Baby Tetanus titer HEU vs HUU") + 
  guides(color = guide_legend(title = "Study site"))  + 
  ylim(0, 7.0) + theme(legend.title = element_text(face = "bold"))

# add stats
a_my_comparisons <- list( c("0", "1"))

symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
Tetanus.vs.delivery <- p1 + stat_compare_means(method = "wilcox.test", 
                                              comparisons = a_my_comparisons, 
                                              label = "p.signif", symnum.args = symnum.args) 
Tetanus.vs.delivery # mixed-feeding group had a higher antibody titer?





#-----------------------------------------------------------------------------------------------------
#  Mode of feeding (supp fig)  ----
#-----------------------------------------------------------------------------------------------------
feeding_data <- tibble(data.frame(sample_data(Phy.f_std))) %>% select(PID, Visit, Study_site, Status2, 
                                                                  Mode_of_feeding,  
                                                                  Mode_of_feeding_birth, Mode_of_feeding_W15, 
                                                                  Feeding_W15_strictver, Feeding_W15_lenientver, 
                                                                  Feeding_Duration_days, Baby_titer) 

View(feeding_data)
str(feeding_data)

feeding_data2 <- feeding_data %>% mutate(FeedingDuration_cat = 
                                           case_when(
                                             (Feeding_Duration_days < 60 ) ~ "upto2M",
                                             (Feeding_Duration_days >= 69) ~ "morethan2M"
                                           )) 

sample_data(Phy.f_std)$FeedingDuration_cat <- feeding_data2$FeedingDuration_cat
sample_data(Phy.f_std)$FeedingDuration_cat <- as.factor(sample_data(Phy.f_std)$FeedingDuration_cat)


### How the alpha-diversity transitioned over time depending on the mode of feeding
Phy.f_std_SA <- subset_samples(Phy.f_std, Study_site == "South Africa")
p1 <- plot_richness(Phy.f_std_SA, x= "Visit", measures="Shannon",
                    # title = "Alpha Diversity at birth"
)+xlab("Mode of feeding")+ facet_wrap(~Feeding_W15_strictver)+
  geom_boxplot(alpha=0.6)+ geom_line(aes(group=PID2),size = 0.1) +
  theme_bw()+ theme(legend.position="none", axis.text.x=element_text(size=17, face="bold"), 
                    strip.text = element_text(size=17, face="bold"), axis.text.y=element_text(size=16, face="bold"),
                    text = element_text(size=17, face="bold"))

a_my_comparisons <- list(c("Week 1", "Week 15"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
Supp_feeding.strict <- p1 + stat_compare_means(method = "wilcox.test",  comparisons = a_my_comparisons,label = "p.signif", symnum.args = symnum.args) 
Supp_feeding.strict




#Question - will the duration of breastfeeding matters in terms of microbiota profile at W15?
#### alpha - Feeding (1) ----- Fig S_A (or B...)
Phy.f_std_SA_W15 <- subset_samples(Phy.f_std, Study_site == "South Africa" & Visit =="Week 15")
levels(sample_data(Phy.f_std_SA_W15)$Feeding_W15_strictver) # [1] "Exclusive_breastfeeding" "Mixed_feeding"  
levels(sample_data(Phy.f_std_SA_W15)$Feeding_W15_strictver) <- c("Exclusive_breastfeeding"="EBF","Mixed_feeding"="MF")



# alpha diversity - mode of feeding - Feeding_W15_strictver (1)
p1 <- plot_richness(Phy.f_std_SA_W15, x= "Feeding_W15_strictver", measures="Shannon",
                    # title = "Alpha Diversity at birth"
)+xlab("Mode of feeding")+
  geom_boxplot(alpha=0.6)+
  theme_bw()+ theme(legend.position="none", axis.text.x=element_text(size=17, face="bold"), 
                    strip.text = element_text(size=17, face="bold"), axis.text.y=element_text(size=16, face="bold"),
                    text = element_text(size=17, face="bold"))

a_my_comparisons <- list(c("EBF", "MF"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
Supp_feeding.strict <- p1 + stat_compare_means(method = "wilcox.test",  comparisons = a_my_comparisons,label = "p.signif", symnum.args = symnum.args) 
Supp_feeding.strict



ggsave(file = "./figure_output/Supp_feeding.strict.pdf", plot = Supp_feeding.strict, dpi = 600, width = 6, height = 7) # in this way, I can control the size of the figure
dev.off()



# alpha diversity - mode of feeding (2) - Feeding_W15_lenientver
levels(sample_data(Phy.f_std_SA_W15)$Feeding_W15_lenientver)  
# [1] "Breastfeeding"           "Exclusive_breastfeeding" "Switched_to_formula"   

sample_data(Phy.f_std_SA_W15)$Feeding_W15_lenientver <-factor(sample_data(Phy.f_std_SA_W15)$Feeding_W15_lenientver, levels = c("Exclusive_breastfeeding", "Breastfeeding", "Switched_to_formula"))
levels(sample_data(Phy.f_std_SA_W15)$Feeding_W15_lenientver) <- c("Exclusive_breastfeeding"="EBF","Breastfeeding"="BF", "Switched_to_formula" = "FF" )


p1 <- plot_richness(Phy.f_std_SA_W15, x= "Feeding_W15_lenientver", measures="Shannon",
                    # title = "Alpha Diversity at birth"
)+xlab("Mode of feeding")+
  geom_boxplot(alpha=0.6)+
  theme_bw()+ theme(legend.position="none", axis.text.x=element_text(size=17, face="bold"), 
                    strip.text = element_text(size=17, face="bold"), axis.text.y=element_text(size=16, face="bold"),
                    text = element_text(size=17, face="bold"))

a_my_comparisons <- list(c("EBF", "BF"), c("BF", "FF"), c("EBF", "FF"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
Supp_feeding.lenient <- p1 + stat_compare_means(method = "wilcox.test",  comparisons = a_my_comparisons,label = "p.signif", symnum.args = symnum.args) 
Supp_feeding.lenient

ggsave(file = "./figure_output/Supp_feeding.lenient.pdf", plot = Supp_feeding.lenient, dpi = 600, width = 8, height = 7) # in this way, I can control the size of the figure
dev.off()


# alpha diversity - mode of feeding (3) - FeedingDuration_cat
feeding_data2 <- feeding_data %>% mutate(FeedingDuration_cat = 
                                           case_when(
                                             (Feeding_Duration_days < 60 ) ~ "<2 months",
                                             (Feeding_Duration_days >= 60) ~ "3 months"
                                           )) 

sample_data(Phy.f_std)$FeedingDuration_cat <- feeding_data2$FeedingDuration_cat
sample_data(Phy.f_std)$FeedingDuration_cat <- as.factor(sample_data(Phy.f_std)$FeedingDuration_cat)
sample_data(Phy.f_std)$FeedingDuration_cat  <-factor(sample_data(Phy.f_std)$FeedingDuration_cat, levels = c("3 months", "<2 months"))


Phy.f_std_SA_W15 <- subset_samples(Phy.f_std, Study_site == "South Africa" & Visit =="Week 15")

levels(sample_data(Phy.f_std_SA_W15)$FeedingDuration_cat)  


p1 <- plot_richness(Phy.f_std_SA_W15, x= "FeedingDuration_cat", measures="Shannon",
                    # title = "Alpha Diversity at birth"
)+xlab("Duration of breast feeding")+
  geom_boxplot(alpha=0.6)+
  theme_bw()+ theme(legend.position="none", axis.text.x=element_text(size=17, face="bold"), 
                    strip.text = element_text(size=17, face="bold"), axis.text.y=element_text(size=16, face="bold"),
                    text = element_text(size=17, face="bold"))

a_my_comparisons <- list(c("<2 months", "3 months"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
Supp_feeding.duration <- p1 + stat_compare_means(method = "wilcox.test",  comparisons = a_my_comparisons,label = "p.signif", symnum.args = symnum.args) 


Supp_feeding.duration

ggsave(file = "./figure_output/Supp_feeding.duration.pdf", plot = Supp_feeding.duration, dpi = 600, width = 6, height = 7) # in this way, I can control the size of the figure
dev.off()




# alpha diversity - mode of feeding (4) - comparing "exclusively breastfed" SA and Nigeria
Phy.f_std_EBF_W15 <- subset_samples(Phy.f_std, Feeding_W15_strictver == "Exclusive_breastfeeding" & Visit =="Week 15")

p1 <- plot_richness(Phy.f_std_EBF_W15, x= "Study_site", measures="Shannon",color = "Study_site",
                    # title = "Alpha Diversity at birth"
)+xlab("Study site")+
  geom_boxplot(alpha=0.6)+
  theme_bw()+ theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
                    strip.text = element_blank(),  axis.text.y=element_text(size=16, face="bold"),
                    text = element_text(size=17, face="bold"))+ labs(y= "Shannon Index")

a_my_comparisons <- list(c("South Africa", "Nigeria"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
Supp_feeding.EBF <- p1 + stat_compare_means(method = "wilcox.test",  comparisons = a_my_comparisons,label = "p.signif", symnum.args = symnum.args) +  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
Supp_feeding.EBF


ggsave(file = "./figure_output/Supp_feeding.EBF.color.pdf", plot = Supp_feeding.EBF, dpi = 600, width = 6, height = 7) # in this way, I can control the size of the figure
dev.off()


### Beta-diversity ------
#### Beta-diversity and feeding - stricter

Phy.f_std_SA_W15 <- subset_samples(Phy.f_std, Study_site == "South Africa" & Visit =="Week 15")
levels(sample_data(Phy.f_std_SA_W15)$Feeding_W15_strictver) <- c("Exclusive_breastfeeding"="EBF","Mixed_feeding"="MF")

ord.BC <- ordinate(Phy.f_std_SA_W15, method = "PCoA", distance = "jsd", k=3, trymax=1000) # try PcOA as NMDS does not work
color = c("Feeding_W15_strictver")
#title= c("PCoA of 16S microbiome at birth (Jensen-Shannon distance; k=3)")  
MDS = plot_ordination(Phy.f_std_SA_W15, ord.BC, color = color #, shape = shape#, title = title
)
Supp_feeding.Beta.strict <- MDS+theme_bw() +theme(axis.text=element_text(size=17, face="bold"),
                      axis.title=element_text(size=17,face="bold"), 
                      strip.text = element_text(size=17,face="bold"),
                      legend.title=element_text(size=17,face="bold"),
                      legend.text = element_text(size =17),
                      legend.position="bottom")+ 
  labs(color=color)+geom_point(size=3)+ 
  guides(color = guide_legend(title = "Mode of feeding"))+  # change the ledgend labee title
 # guides(shape = guide_legend(title = "Study site"))+ # change the shape title
  annotate("text", x=- 0.28, y= -0.3, label= "p = 0.03", size = 5)

Supp_feeding.Beta.strict

ggsave(file = "./figure_output/Supp_feeding.Beta.strict.pdf", plot = Supp_feeding.Beta.strict, dpi = 700, width =6, height = 6) # in this way, I can control the size of the figure
dev.off()


###look at the PERMANOVA (Feeding_W15_strictver)
sample_data(Phy.f_std_SA_W15)$Feeding_W15_strictver  # Levels: EBF MF

diss <- phyloseq::distance(Phy.f_std_SA_W15, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_SA_W15_df <- data.frame(row.names=sample_names(Phy.f_std_SA_W15),sample_data(Phy.f_std_SA_W15)) #create dataframe with metadata

# Whether "ARV intake" differ significantly from each other
set.seed(2)
adonis2(diss ~ Feeding_W15_strictver, data=Phy.f_std_SA_W15_df, permutations = 999, by = "terms", na.omit)  #p = 0.03


###  Beta-diversity (Feeding_W15_lenientver)

levels(sample_data(Phy.f_std_SA_W15)$Feeding_W15_lenientver)  
# [1] "Breastfeeding"           "Exclusive_breastfeeding" "Switched_to_formula"   

sample_data(Phy.f_std_SA_W15)$Feeding_W15_lenientver <-factor(sample_data(Phy.f_std_SA_W15)$Feeding_W15_lenientver, levels = c("Exclusive_breastfeeding", "Breastfeeding", "Switched_to_formula"))
levels(sample_data(Phy.f_std_SA_W15)$Feeding_W15_lenientver) <- c("Exclusive_breastfeeding"="EBF","Breastfeeding"="BF", "Switched_to_formula" = "FF" )

ord.BC <- ordinate(Phy.f_std_SA_W15, method = "PCoA", distance = "jsd", k=3, trymax=1000) # try PcOA as NMDS does not work
color = c("Feeding_W15_lenientver")
#title= c("PCoA of 16S microbiome at birth (Jensen-Shannon distance; k=3)")  
MDS = plot_ordination(Phy.f_std_SA_W15, ord.BC, color = color #, shape = shape, title = title
)
Supp_feeding.Beta.lenient <- MDS+theme_bw() +theme(axis.text=element_text(size=17, face="bold"),
                                                  axis.title=element_text(size=17,face="bold"), 
                                                  strip.text = element_text(size=17,face="bold"),
                                                  legend.title=element_text(size=17,face="bold"),
                                                  legend.text = element_text(size =17),
                                                  legend.position="bottom")+ 
  labs(color=color)+geom_point(size=3)+ 
  guides(color = guide_legend(title = "Mode of feeding"))+  # change the ledgend labee title
  guides(shape = guide_legend(title = "Study site"))+ # change the shape title
  annotate("text", x=- 0.28, y= -0.3, label= "p = 0.006", size = 5)

Supp_feeding.Beta.lenient

ggsave(file = "./figure_output/Supp_feeding.Beta.lenient.pdf", plot = Supp_feeding.Beta.lenient, dpi = 700, width =7.5, height = 6) # in this way, I can control the size of the figure
dev.off()


###look at the PERMANOVA (Feeding_W15_lenientver)
sample_data(Phy.f_std_SA_W15)$Feeding_W15_lenientver  # Levels: EBF BF FF

diss <- phyloseq::distance(Phy.f_std_SA_W15, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_SA_W15_df <- data.frame(row.names=sample_names(Phy.f_std_SA_W15),sample_data(Phy.f_std_SA_W15)) #create dataframe with metadata

# Whether "ARV intake" differ significantly from each other
set.seed(2)
test <- adonis2(diss ~ Feeding_W15_lenientver, data=Phy.f_std_SA_W15_df, permutations = 999, by = "terms", na.omit)  #p = 0.006
str(test)


###do the pair-wise comparison!!!
## since the adonis oveall for "Feeding_W15_lenientver" was significant, let's find out which pairs are significantly different
# https://www.youtube.com/watch?v=1ETBgbXl-BM
sample_data(Phy.f_std_SA_W15)$Feeding_W15_lenientver # EBF BF FF
pairwise_p <- numeric()  # create a vector
set.seed(2)

Phy.f_std_pair1 <- subset_samples(Phy.f_std_SA_W15, Feeding_W15_lenientver == "EBF" |Feeding_W15_lenientver == "BF" ) 
diss <- phyloseq::distance(Phy.f_std_pair1, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_pair1_df <- data.frame(row.names=sample_names(Phy.f_std_pair1),sample_data(Phy.f_std_pair1)) #create dataframe with metadata
pair1_test <- adonis2(diss ~ Feeding_W15_lenientver, data=Phy.f_std_pair1_df, permutations = 999, by = "terms", na.omit)  #p= 0.387
pairwise_p["EBF.vs.BF"] <- pair1_test$`Pr(>F)`[1]


Phy.f_std_pair2 <- subset_samples(Phy.f_std_SA_W15, Feeding_W15_lenientver == "BF" |Feeding_W15_lenientver == "FF" ) 
diss <- phyloseq::distance(Phy.f_std_pair2, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_pair2_df <- data.frame(row.names=sample_names(Phy.f_std_pair2),sample_data(Phy.f_std_pair2)) #create dataframe with metadata
pair2_test <- adonis2(diss ~ Feeding_W15_lenientver, data=Phy.f_std_pair2_df, permutations = 999, by = "terms", na.omit)  #p= 0.056
pairwise_p["BF.vs.FF"] <-pair2_test$`Pr(>F)`[1]


Phy.f_std_pair3 <- subset_samples(Phy.f_std_SA_W15, Feeding_W15_lenientver == "EBF" |Feeding_W15_lenientver == "FF" ) 
diss <- phyloseq::distance(Phy.f_std_pair3, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_pair3_df <- data.frame(row.names=sample_names(Phy.f_std_pair3),sample_data(Phy.f_std_pair3)) #create dataframe with metadata
pair3_test <- adonis2(diss ~ Feeding_W15_lenientver, data=Phy.f_std_pair3_df, permutations = 999, by = "terms", na.omit)  #p= 0.002
pairwise_p["EBF.vs.FF"] <- pair3_test$`Pr(>F)`[1]

pairwise_p #
# EBF.vs.BF  BF.vs.FF EBF.vs.FF 
# 0.387     0.046     0.002 

p.adjust(pairwise_p, method = "BH") # Normal.vs.Wasting and  Overweight.vs.Wasting became not significnat after adjusting for p-value.
# EBF.vs.BF  BF.vs.FF EBF.vs.FF 
# 0.387     0.069     0.006 




#### Beta-diversity and feeding (2) - plot FeedingDuration_cat
Phy.f_std_SA_W15 <- subset_samples(Phy.f_std, Study_site == "South Africa" & Visit =="Week 15")
ord.BC <- ordinate(Phy.f_std_SA_W15, method = "PCoA", distance = "jsd", k=3, trymax=1000) # try PcOA as NMDS does not work
color = c("FeedingDuration_cat")
#title= c("PCoA of 16S microbiome at birth (Jensen-Shannon distance; k=3)")  
MDS = plot_ordination(Phy.f_std_SA_W15, ord.BC, color = color #, shape = shape, title = title
)
Supp_feeding.Beta.Duration_cat <- MDS+theme_bw() +theme(axis.text=element_text(size=13, face="bold"),
                      axis.title=element_text(size=17, face="bold"), 
                      legend.title=element_text(size=17, face="bold"),
                      legend.text = element_text(size =17),
                      legend.position="bottom")+ 
  labs(color=color)+geom_point(size=3)+
  guides(color = guide_legend(title = "Mode of feeding"))+  # change the ledgend labee title
  guides(shape = guide_legend(title = "Study site"))+ # change the shape title
  annotate("text", x=- 0.28, y= -0.3, label= "p = 0.001", size = 5)


Supp_feeding.Beta.Duration_cat

ggsave(file = "./figure_output/Supp_feeding.Beta.Duration_cat.pdf", plot = Supp_feeding.Beta.Duration_cat, dpi = 700, width =6, height = 6) # in this way, I can control the size of the figure
dev.off()


###look at the PERMANOVA 
Phy.f_std_SA_W15 <- subset_samples(Phy.f_std, Study_site == "South Africa" & Visit =="Week 15")
sample_data(Phy.f_std_SA_W15)$FeedingDuration_cat  # Levels:  3 months <2 months

diss <- phyloseq::distance(Phy.f_std_SA_W15, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_SA_W15_df <- data.frame(row.names=sample_names(Phy.f_std_SA_W15),sample_data(Phy.f_std_SA_W15)) #create dataframe with metadata

# Whether "ARV intake" differ significantly from each other
set.seed(2)
adonis2(diss ~ FeedingDuration_cat, data=Phy.f_std_SA_W15_df, permutations = 999, by = "terms", na.omit)  #0.001
adonis2(diss ~ FeedingDuration_cat*Status2, data=Phy.f_std_SA_W15_df, permutations = 999, by = "terms", na.omit) 




#### Beta-diversity and feeding (3) -gradient with "duration of breastfeeding"
Phy.f_std_SA_W15 <- subset_samples(Phy.f_std, Study_site == "South Africa" & Visit =="W15")
ord.BC <- ordinate(Phy.f_std_SA_W15, method = "PCoA", distance = "jsd", k=3, trymax=1000) # try PcOA as NMDS does not work
color = c("Feeding_Duration_days")
shape = c("Study_site")
#title= c("PCoA of 16S microbiome at birth (Jensen-Shannon distance; k=3)")  
MDS = plot_ordination(Phy.f_std_SA_W15, ord.BC, color = color, shape = shape#, title = title
)
MDS+theme_bw() +scale_color_viridis_c(option = "plasma")+
  theme(axis.text=element_text(size=13, face="bold"),
        axis.title=element_text(size=16,face="bold"), 
        legend.title=element_text(size=14),
        legend.text = element_text(size =12))+
  labs(color= "Duration of EBF (days)")+geom_point(size=3)+
  guides(shape = guide_legend(title = "Study site"))



### Heather told me to look at beta diversity of EBF infants and compare by site
Phy.f_std_EBF_W15 <- subset_samples(Phy.f_std, Feeding_W15_strictver == "Exclusive_breastfeeding" & Visit =="Week 15")
ord.BC <- ordinate(Phy.f_std_EBF_W15, method = "PCoA", distance = "bray", k=3, trymax=1000) # try PcOA as NMDS does not work
color = c("Study_site")
#title= c("PCoA of 16S microbiome at birth (Jensen-Shannon distance; k=3)")  
MDS = plot_ordination(Phy.f_std_EBF_W15, ord.BC, color = color)

plot <- MDS+theme_bw() +theme(axis.text=element_text(size=15, face="bold"),
                      axis.title=element_text(size=15,face="bold"), 
                      strip.text = element_text(size=15,face="bold"),
                      legend.title=element_text(size=15,face="bold"),
                      legend.text = element_text(size =15),
                      legend.position="bottom")+ 
  labs(color=color)+geom_point(size=2)+ 
  guides(color = guide_legend(title = "Study site")) # change the color title

plot
ggsave(file = "./figure_output/Supp_feeding.EBFcomparison.resize.bray.pdf", plot = plot, dpi = 700, width =6, height = 4) # in this way, I can control the size of the figure
dev.off()


###look at the PERMANOVA 
diss <- phyloseq::distance(Phy.f_std_EBF_W15, "bray", parallel=TRUE)  
Phy.f_std_EBF_W15_df <- data.frame(row.names=sample_names(Phy.f_std_EBF_W15),sample_data(Phy.f_std_EBF_W15)) #create dataframe with metadata

# Whether "Study site" still differ significantly from each other even only comparing EBF infants.
set.seed(2)
adonis2(diss ~ Study_site, data=Phy.f_std_EBF_W15_df, permutations = 999, by = "terms", na.omit)  #0.001 - still significant



# now combine the three plots...
library(ggpubr)
Supp_feeding.strict
Supp_feeding.lenient
Supp_feeding.duration
Supp_feeding.EBF
Supp_feeding.Beta.strict
Supp_feeding.Beta.Duration_cat


top.4.panel <- ggarrange(Supp_feeding.strict, NULL,Supp_feeding.lenient, NULL, Supp_feeding.duration,
                         NULL,NULL,NULL,NULL,NULL,
                         Supp_feeding.EBF, labels = c(), nrow = 3, ncol =5 , widths = c(1, 0.07, 1, 0.07, 1),heights = c(1, 0.07, 1))
ggsave(file = "./figure_output/SuppX_Feeding.top.4.panel.pdf", plot = top.4.panel, dpi = 900, width = 12 , height = 9) # in this way, I can control the size of the figure
dev.off()


bottom.2.panel <- ggarrange(Supp_feeding.Beta.strict, NULL,Supp_feeding.Beta.Duration_cat,NULL, labels = c(), nrow = 1, ncol =4, widths = c(1, 0.1, 1,0.2))
ggsave(file = "./figure_output/SuppX_Feeding.bottom.2.panel.pdf", plot = bottom.2.panel, dpi = 900, width = 10 , height = 5) # in this way, I can control the size of the figure
dev.off()


###### Relative abundace -"Mode of feeding"

barplot_stag_BEAMING_sort <- function(physeq) {
  p1 <- tax_glom(physeq, taxrank = 'Species') #agglumerate at species level
  p1_30 = prune_taxa(names(sort(taxa_sums(p1), TRUE))[1:30], p1) #use top 30 most abundant taxa for plotting
  p2 <- transform_sample_counts(p1_30, function(x) x/sum(x)) #get abundance in %
  p3 <- psmelt(p2)  #create dataframe from phyloseq object
  p3$Species <- as.character(p3$Species) #convert to character
  #sort based on abundance of dominant bacteria in each cluster
  p4 <- p3[,c("Sample", "Species","CST_pam_3", "Abundance")]
  p4 <- dcast(p3, Sample + CST_pam_3 ~ Species, value.var = "Abundance", fun.aggregate = sum)
  p4[p4$CST_pam_3 =="Cluster 1",] <- p4[p4$CST_pam_3=="Cluster 1",][order(p4[p4$CST_pam_3=="Cluster 1","coli"]),]
  p4[p4$CST_pam_3 =="Cluster 2",] <- p4[p4$CST_pam_3=="Cluster 2",][order(p4[p4$CST_pam_3=="Cluster 2","longum"]),]
  #p4[p4$CST_pam_3 =="Cluster 3",] <- p4[p4$CST_pam_3=="Cluster 3",][order(p4[p4$CST_pam_3=="Cluster 3","faecalis"]),]
  #reorder p3
  p3$Sample <- factor(p3$Sample, levels=c(p4$Sample))
  otus <- otu_table(p1_30)
  tax_table(p1_30)
  #annotation 
  labs <- tax.lab(physeq,otus=otus)
  #set color palette to accommodate the number of species
  colourCount = length(unique(p3$Genus))
  getPalette = colorRampPalette(brewer.pal(8, "Accent"))  
  pal<- c("#7FC97F", "#CBB1C3", "#FDDA8D", "#83A3A7", "#D01487", "#8EC293", 
          "#BB5B19","#DAB6B1", "#FEE992","#FDCA89",  "#EC0877","#C9482C",
          "#9DBBA8", "#E9BA9E","#E3EA9C", "#FEF897", "#4B61AA", "#E01D5E", 
          "#ACB5BC", "#90603F","#F8BE8B", "#77479F", "#D43345", "#7B6352",
          "#BBAED1", "#A65E2C", "#B3C7A1", "#A32D93","#5380AC", "#666666")
  
  pal<- c("steelblue3", "#8EC293", "#FEE992", "azure4", "brown3", 
          "#CBB1C3", "tomato2","lightgreen", "cyan","#FDCA89",  
          "#EC0877","sienna1","lightseagreen", "lightcoral","#E3EA9C", 
          "yellow", "blue3", "#E01D5E", "chartreuse3", "chocolate4",
          "darkgoldenrod1", "#77479F", "honeydew3", "yellow4", "violet", 
          "tan3", "seagreen", "#A32D93","lightskyblue", "#666666")
  #annotation 
  #plot
  barplot_species <- ggplot(data=p3, aes(x=Sample, y=Abundance, fill=Species)) 
  barplot_species <- barplot_species + geom_bar(aes(), stat="identity", position="stack") + 
    facet_wrap(~ FeedingDuration_cat, scales = "free") + guides(fill=guide_legend(nrow=3)) + 
    scale_fill_manual(values=pal) + ylab("Relative abundance") +xlab("Participants") + 
    #ggtitle("Bacterial composition at each time pooint by study site") + 
    theme(legend.position="bottom", strip.background = element_rect(fill="lightgrey", color = "black"), 
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          panel.background=element_blank(), panel.border=element_blank(),
          axis.title.y = element_text(size = 17, colour = "black",face = "bold"),
          axis.title.x = element_blank(),
          strip.text = element_text(size=17, face="bold"),
          axis.text.y=element_text(size=17, face="bold"), 
          text = element_text(size=15, face="bold"))
  barplot_species
  return(barplot_species)
}



Phy.f_std_SA_W15 <- subset_samples(Phy.f_std, Study_site == "South Africa" & Visit =="Week 15")
sample_data(Phy.f_std_SA_W15)$FeedingDuration_cat # Levels: morethan2M upto2M
sample_data(Phy.f_std_SA_W15)$FeedingDuration_cat  <-factor(sample_data(Phy.f_std_SA_W15)$FeedingDuration_cat, levels = c("morethan2M","upto2M")) # cahnge the order


Supp_feeding.relative.plot <- barplot_stag_BEAMING_sort(Phy.f_std_SA_W15)


Supp_feeding.relative.plot

ggsave(file = "./figure_output/Supp_feeding.relative.plot.pdf", plot = Supp_feeding.relative.plot, dpi = 900, width = 18, height = 9) # in this way, I can control the size of the figure
dev.off()


### Deseq2 result - check whats abundant between EBF vs FF
Phy.f_std <- readRDS("Phy.f_std.RDS")
Phy.f_std_birth <- subset_samples(Phy.f_std, Visit=="Week 15" & Study_site == "South Africa")
Phy.f_std_birth_feeding <- subset_samples(Phy.f_std_birth, Feeding_W15_lenientver =="Exclusive_breastfeeding" | Feeding_W15_lenientver =="Switched_to_formula")

Phy.f_std_birth_feeding_glom <- tax_glom.kv(Phy.f_std_birth_feeding)

sample_data(Phy.f_std_birth_feeding)$FeedingDuration_cat  #  3 months <2 months
sample_data(Phy.f_std_birth_feeding)$FeedingDuration_cat  <-factor(sample_data(Phy.f_std_birth_feeding)$FeedingDuration_cat, levels = c("upto2M","morethan2M")) # cahnge the order



#Convert the phyloseq object to a DESeqDataSet and run DESeq2:
ds = phyloseq_to_deseq2(Phy.f_std_birth_feeding, ~ FeedingDuration_cat)
#ds$Study_site <- relevel(ds$Study_site, ref = "South.Africa")
ds = DESeq(ds, fitType = "local")  # get an error...rror in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : every gene contains at least one zero, cannot compute log geometric means
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(ds), 1, gm_mean)
ds = estimateSizeFactors(ds, geoMeans = geoMeans)
ds = DESeq(ds, fitType="local")
res = results(ds)
res = res[order(res$padj, na.last = NA), ]
alpha = 0.01 #this is your significance level, adjust if desired
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Phy.f_std_birth_feeding)[rownames(sigtab), ], "matrix"))
head(sigtab)
sigtab %>% View()

resultsNames(ds) 
rownames(sigtab) 

#only show bugs annotated at the Species level
#sigtabspec = subset(sigtab, !is.na(Species))
sigtabgen = subset(sigtab, !is.na(Genus))
#sigtabord = subset(sigtab, !is.na(Order))
head(sigtabgen)

#write.csv(sigtabgen, file = "./figure_output/Sfig-sigtabgen_Supp_feeding.deseq.csv")

#sort for pretty plotting
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

#Alternative: Phylum level, for plotting purposes - don't do both
#x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
#x = sort(x, TRUE)
#sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

#fix taxrank and subset to highly significant asvs
sigtabgen.hs <- sigtabgen[which(sigtabgen$padj < 0.005),]

#subset to greater magnitude of FC
sigtabgen$Species[is.na(sigtabgen$Species)] = "(unclassified)"
sigtabgen$Species[sigtabgen$Species == "NA"] <- "unknown"
sigtabgen$genspec <- paste(sigtabgen$Genus, sigtabgen$Species, sep = " ")
sigtab.merged.hs <- sigtabgen[which(sigtabgen$padj < 0.005),]

sigtabgen %>%
  subset(log2FoldChange > 0.5) -> sigtabgen.hi
sigtabgen %>%
  subset(log2FoldChange < -0.5) -> sigtabgen.low

sigtab.merged.hs %>%
  subset(log2FoldChange > 0.5) -> sigtabgen.hi

sigtab.merged.hs %>%
  subset(log2FoldChange < -0.5) -> sigtabgen.low

sigtabgen.fc <- rbind.data.frame(sigtabgen.hi, sigtabgen.low)                          
View(sigtabgen.fc)

#plot the results (Genus currently set, and colored by Order, set as desired)
SuppFig.feeding.deseq <- ggplot(sigtabgen.fc, aes(y=Genus, x=log2FoldChange, fill=genspec)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_vline(xintercept = 0.5, color = "indianred3", size = 0.5, lty=2) +
  geom_vline(xintercept = -0.5, color = "indianred3", size = 0.5, lty=2) +
  geom_point(pch = 21, size=3)+ theme_bw()+
  theme(axis.text.y = element_text(size = 14, colour = "black",face = "bold"),
        axis.title.y = element_text(size = 17, colour = "black",face = "bold"),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 17, colour = "black",face = "bold"),
        strip.background =element_rect(fill="white"),
        strip.text.x = element_text(colour = "mediumpurple4", size = 15, hjust = 0.5, vjust = 0.5, face = "bold"),
        legend.text = element_text(size = 14, colour = "black"),
        legend.title = element_text(size=15)) +  theme(legend.position="right", legend.title=element_blank())+
  guides(fill = guide_legend(nrow = 12)) + guides(fill=guide_legend(ncol = 2))+
  labs(x=expression(log[2]~fold~change))  #scale_fill_manual(values=custom_colors_36) + #+ ggtitle("DESeq2 SA vs Nigeria at birth, p=0.01, w/o melting the lowest taxonomy")

SuppFig.feeding.deseq # something wrong with DEseq - it looks like opposite

ggsave(file = "./figure_output/Supp_feeding.deseq.pdf", plot = SuppFig.feeding.deseq, dpi = 900, width = 15, height = 7) # in this way, I can control the size of the figure
dev.off()


#-----------------------------------------------------------------------------------------------------
#  CD4 count & Viral load ----
#-----------------------------------------------------------------------------------------------------
# CD4 count

# take a look at the related variables.
dat <- tibble(data.frame(sample_data(Phy.f_std))) 
dat$Last_ViralLoad[dat$Last_ViralLoad == "<20"] <- "19"
dat$Last_ViralLoad[dat$Last_ViralLoad == "<35"] <- "19"
dat$Last_ViralLoad <- as.numeric(dat$Last_ViralLoad)

dat2 <- dat %>% 
  select(PID, Study_site, Visit, Status2, Last_CD4count, Last_ViralLoad) %>% 
  mutate(Last_CD4count_cat = 
           case_when(
             (Last_CD4count < 250) ~ "CD4_below250",
             (Last_CD4count >= 250) ~ "CD4_above250")) %>% 
  mutate(Last_ViralLoad_cat = 
           case_when((Last_ViralLoad <= 19) ~ "VL_undetectable",
                     (Last_ViralLoad > 250) ~ "VL_detectable"))

sample_data(Phy.f_std)$Last_CD4count_cat <- dat2$Last_CD4count_cat  # adding "Last_CD4count_cat"
sample_data(Phy.f_std)$Last_ViralLoad_cat <- dat2$Last_ViralLoad_cat   # adding "Last_ViralLoad_cat"



### alpha diversity ---- CD4 count 
sample_data(Phy.f_std)$Last_CD4count_cat <- as.factor(sample_data(Phy.f_std)$Last_CD4count_cat) # Levels: CD4_above250 CD4_below250
sample_data(Phy.f_std)$Last_ViralLoad_cat <- as.factor(sample_data(Phy.f_std)$Last_ViralLoad_cat) # Levels: VL_detectable VL_undetectable

Phy.f_std_HEU <- subset_samples(Phy.f_std, Status2 == "iHEU") 
Phy.f_std_HEU_v2 <- subset_samples(Phy.f_std_HEU, Last_CD4count_cat == "CD4_above250" | Last_CD4count_cat == "CD4_below250") 

levels(sample_data(Phy.f_std_HEU_v2)$Last_CD4count_cat) <- c("CD4_above250"=">250", "CD4_below250" = "<250")


p1 <- plot_richness(Phy.f_std_HEU_v2, x= "Last_CD4count_cat", measures="Shannon",
                    # title = "Alpha Diversity at birth"
)+xlab("Maternal CD4 count")+
  geom_boxplot(alpha=0.6)+ facet_wrap(Study_site ~ Visit)+
  theme_bw()+ theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
                    strip.text = element_text(size=17, face="bold"), axis.text.y=element_text(size=16, face="bold"),
                    text = element_text(size=17, face="bold"))

# and add statistics
# https://rpubs.com/lconteville/713954
a_my_comparisons <- list(c(">250", "<250"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
SuppFig.CD4.alpha <- p1 + stat_compare_means(method = "wilcox.test", 
                                comparisons = a_my_comparisons, 
                                label = "p.signif", symnum.args = symnum.args)+ 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) # to give a bit of space above the p-value
  #significantly different - CT increases the alpha diversity over the time whereas decased in Nigerian infants


SuppFig.CD4.alpha
ggsave(file = "./figure_output/SuppFig.CD4.alpha.pdf", plot = SuppFig.CD4.alpha, dpi = 700, width = 7, height = 8) # in this way, I can control the size of the figure
dev.off()


# Beta-diversity - CD4
levels(sample_data(Phy.f_std_HEU_v2)$Last_CD4count_cat) <- c("CD4_above250"=">250", "CD4_below250" = "<250")
ord.BC <- ordinate(Phy.f_std_HEU_v2, method = "PCoA", distance = "jsd", k=3, trymax=1000) # try PcOA as NMDS does not work
color = c("Last_CD4count_cat")
shape = c("Study_site")
#title= c("PCoA of 16S microbiome at birth (Jensen-Shannon distance; k=3)")  
MDS = plot_ordination(Phy.f_std_HEU_v2, ord.BC, color = color, shape = shape#, title = title
)
SuppFig.CD4.beta <- MDS+theme_bw() +theme(axis.text=element_text(size=13, face="bold"),
                      axis.title=element_text(size=16,face="bold"), 
                      legend.title=element_text(size=14),
                      legend.text = element_text(size =12))+ 
  labs(color=color)+geom_point(size=3)+
  guides(color = guide_legend(title = "Maternal CD4 count"))+  # change the ledgend labee title
  guides(shape = guide_legend(title = "Study site")) # change the shape title


SuppFig.CD4.beta
ggsave(file = "./figure_output/SuppFig.CD4.beta.pdf", plot = SuppFig.CD4.beta, dpi = 700, width = 7, height = 6) # in this way, I can control the size of the figure
dev.off()


### PERMANOVA - CD4 count
set.seed(20)
diss <- phyloseq::distance(Phy.f_std_HEU_v2, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_HEU_v2_df <- data.frame(row.names=sample_names(Phy.f_std_HEU_v2),sample_data(Phy.f_std_HEU_v2)) #create dataframe with metadata
# Whether "Last_CD4count_cat" differ significantly from each other
adonis2(diss ~ Last_CD4count_cat, data=Phy.f_std_HEU_v2_df, permutations = 999, by = "terms", na.omit)  # 

# only in SA
Phy.f_std_HEU_v2_SA <- subset_samples(Phy.f_std_HEU_v2, Study_site == "South Africa")
set.seed(20)
diss <- phyloseq::distance(Phy.f_std_HEU_v2_SA, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_HEU_v2_SA_df <- data.frame(row.names=sample_names(Phy.f_std_HEU_v2_SA),sample_data(Phy.f_std_HEU_v2_SA)) #create dataframe with metadata
# Whether "Last_CD4count_cat" differ significantly from each other
adonis2(diss ~ Last_CD4count_cat, data=Phy.f_std_HEU_v2_SA_df, permutations = 999, by = "terms", na.omit)  # p=0.271


# only in Nigeria
Phy.f_std_HEU_v2_Nig <- subset_samples(Phy.f_std_HEU_v2, Study_site == "Nigeria")
set.seed(20)
diss <- phyloseq::distance(Phy.f_std_HEU_v2_Nig, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_HEU_v2_Nig_df <- data.frame(row.names=sample_names(Phy.f_std_HEU_v2_Nig),sample_data(Phy.f_std_HEU_v2_Nig)) #create dataframe with metadata
# Whether "Last_CD4count_cat" differ significantly from each other
adonis2(diss ~ Last_CD4count_cat, data=Phy.f_std_HEU_v2_Nig_df, permutations = 999, by = "terms", na.omit)  # p= 0.498



### Viral load info
### alpha diversity ---- Viral Load
sample_data(Phy.f_std)$Last_ViralLoad_cat <- as.factor(sample_data(Phy.f_std)$Last_ViralLoad_cat) # Levels: VL_detectable VL_undetectable

Phy.f_std_HEU <- subset_samples(Phy.f_std, Status2 == "HEU") 
Phy.f_std_HEU_v2 <- subset_samples(Phy.f_std_HEU, Last_ViralLoad_cat == "VL_detectable" | Last_ViralLoad_cat == "VL_undetectable") 

levels(sample_data(Phy.f_std_HEU_v2)$Last_ViralLoad_cat) <- c("VL_detectable"="Detactable", "VL_undetectable" = "Undetectable")


p1 <- plot_richness(Phy.f_std_HEU_v2, x= "Last_ViralLoad_cat", measures="Shannon",
                    # title = "Alpha Diversity at birth"
)+xlab("Maternal viral load")+
  geom_boxplot(alpha=0.6)+ facet_wrap(Study_site ~ Visit)+
  theme_bw()+ theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
                    strip.text = element_text(size=17, face="bold"), axis.text.y=element_text(size=16, face="bold"),
                    text = element_text(size=17, face="bold"))+ scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# and add statistics
# https://rpubs.com/lconteville/713954
a_my_comparisons <- list(c("Detactable", "Undetectable"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
SuppFig.VL.alpha <- p1 + stat_compare_means(method = "wilcox.test", 
                                comparisons = a_my_comparisons, 
                                label = "p.signif", symnum.args = symnum.args) 
SuppFig.VL.alpha  #significantly different - CT increases the alpha diversity over the time whereas decased in Nigerian infants


SuppFig.VL.alpha
ggsave(file = "./figure_output/SuppFig.VL.alpha.pdf", plot = SuppFig.VL.alpha, dpi = 700, width = 7, height = 8) # in this way, I can control the size of the figure
dev.off()


# for PcOA plot  -  Viral Load
Phy.f_std_HEU <- subset_samples(Phy.f_std, Status2 == "iHEU") 
Phy.f_std_HEU_v2 <- subset_samples(Phy.f_std_HEU, Last_ViralLoad_cat == "VL_detectable" | Last_ViralLoad_cat == "VL_undetectable") 
levels(sample_data(Phy.f_std_HEU_v2)$Last_ViralLoad_cat) <- c("VL_detectable"="Detactable", "VL_undetectable" = "Undetectable")

ord.BC <- ordinate(Phy.f_std_HEU_v2, method = "PCoA", distance = "jsd", k=3, trymax=1000) # try PcOA as NMDS does not work
color = c("Last_ViralLoad_cat")
shape = c("Study_site")
#title= c("PCoA of 16S microbiome at birth (Jensen-Shannon distance; k=3)")  
MDS = plot_ordination(Phy.f_std_HEU_v2, ord.BC, color = color, shape = shape#, title = title
)
SuppFig.VL.beta <- MDS+theme_bw() +theme(axis.text=element_text(size=13, face="bold"),
                                          axis.title=element_text(size=16,face="bold"), 
                                          legend.title=element_text(size=14),
                                          legend.text = element_text(size =12))+ 
  labs(color=color)+geom_point(size=3)+
  guides(color = guide_legend(title = "Maternal viral load"))+  # change the ledgend labee title
  guides(shape = guide_legend(title = "Study site")) # change the shape title


SuppFig.VL.beta
ggsave(file = "./figure_output/SuppFig.VL4.beta.pdf", plot = SuppFig.VL.beta, dpi = 700, width = 7, height = 6) # in this way, I can control the size of the figure
dev.off()


### PERMANOVA - VL count
set.seed(20)
diss <- phyloseq::distance(Phy.f_std_HEU_v2, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_HEU_v2_df <- data.frame(row.names=sample_names(Phy.f_std_HEU_v2),sample_data(Phy.f_std_HEU_v2)) #create dataframe with metadata
# Whether "Last_ViralLoad_cat" differ significantly from each other
adonis2(diss ~ Last_ViralLoad_cat, data=Phy.f_std_HEU_v2_df, permutations = 999, by = "terms", na.omit)  # p=0.706

# only in SA
Phy.f_std_HEU_v2_SA <- subset_samples(Phy.f_std_HEU_v2, Study_site == "South Africa")
set.seed(20)
diss <- phyloseq::distance(Phy.f_std_HEU_v2_SA, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_HEU_v2_SA_df <- data.frame(row.names=sample_names(Phy.f_std_HEU_v2_SA),sample_data(Phy.f_std_HEU_v2_SA)) #create dataframe with metadata
# Whether "Last_ViralLoad_cat" differ significantly from each other
adonis2(diss ~ Last_ViralLoad_cat, data=Phy.f_std_HEU_v2_SA_df, permutations = 999, by = "terms", na.omit)  # p= 0.479


# only in Nigeria
Phy.f_std_HEU_v2_Nig <- subset_samples(Phy.f_std_HEU_v2, Study_site == "Nigeria")
set.seed(20)
diss <- phyloseq::distance(Phy.f_std_HEU_v2_Nig, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_HEU_v2_Nig_df <- data.frame(row.names=sample_names(Phy.f_std_HEU_v2_Nig),sample_data(Phy.f_std_HEU_v2_Nig)) #create dataframe with metadata
# Whether "Last_ViralLoad_cat" differ significantly from each other
adonis2(diss ~ Last_ViralLoad_cat, data=Phy.f_std_HEU_v2_Nig_df, permutations = 999, by = "terms", na.omit)  # p= 0.789


#-----------------------------------------------------------------------------------------------------
#  Mode of delivery ----
#-----------------------------------------------------------------------------------------------------
### alpha diversity

# Alpha Diversity by Mode of delivery (9) 
sample_data(Phy.f_std)
sample_data(Phy.f_std)$DeliverVag <- as.factor(sample_data(Phy.f_std)$DeliverVag)
levels(sample_data(Phy.f_std)$"DeliverVag") <- c("0"="C-section","1"="Vaginal birth")

Phy.f_std_Nigeria <- subset_samples(Phy.f_std, Study_site == "Nigeria")

p1 <- plot_richness(Phy.f_std_Nigeria, x="DeliverVag", measures="Shannon",
                     #title = "Alpha diversirty of delivery mode and Visit (Nigeria)"
                    )+
                    xlab("Mode of delivery")+
    geom_boxplot(alpha=0.6)+ theme_bw()+ 
    theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
          strip.text = element_text(size=17, face="bold"), axis.text.y=element_text(size=16, face="bold"),
          text = element_text(size=17, face="bold")) + facet_wrap(~Visit)

symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

SuppFig.delivery <- p1 + stat_compare_means(method = "wilcox.test", 
                        comparisons = list(c("C-section", "Vaginal birth")),
                        label = "p.signif", symnum.args = symnum.args)


SuppFig.delivery
ggsave(file = "./figure_output/SuppFig.delivery.alpha.pdf", plot = SuppFig.delivery, dpi = 700, width = 7, height = 6) # in this way, I can control the size of the figure
dev.off()


### compare only vaginal birth infants by country (Anna's suggestion)
sample_data(Phy.f_std)$DeliverVag <- as.factor(sample_data(Phy.f_std)$DeliverVag)
levels(sample_data(Phy.f_std)$"DeliverVag") <- c("0"="C-section","1"="Vaginal birth")

Phy.f_std_VagiBirth <- subset_samples(Phy.f_std, DeliverVag == "Vaginal birth")

p1 <- plot_richness(Phy.f_std_VagiBirth, x= "Study_site", measures="Shannon",color = "Study_site",
                    #title = "Alpha diversirty of delivery mode and Visit (Nigeria)"
)+
  xlab("Study site")+
  geom_boxplot(alpha=0.6)+ theme_bw()+ facet_wrap(~ Visit)+
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
        axis.text.y=element_text(size=16, face="bold"),
        text = element_text(size=17, face="bold"))+   labs(y= "Shannon Index")

symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

p1.1 <- p1 + stat_compare_means(method = "wilcox.test", 
                                            comparisons = list(c("South Africa", "Nigeria")),
                                            label = "p.signif", symnum.args = symnum.args)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))



ggsave(file = "./figure_output/SuppFig.VagiBirth.alpha.color.pdf", plot = p1.1, dpi = 700, width = 5, height = 6) # in this way, I can control the size of the figure
dev.off()

### Beta diversity
# for PcOA plot  -  delivery
sample_data(Phy.f_std)$DeliverVag <- as.factor(sample_data(Phy.f_std)$DeliverVag)
levels(sample_data(Phy.f_std)$"DeliverVag") <- c("0"="C-section","1"="Vaginal birth")

Phy.f_std_v2 <- subset_samples(Phy.f_std, DeliverVag == "C-section" | DeliverVag == "Vaginal birth")
Phy.f_std_v3 <-  subset_samples(Phy.f_std_v2, Study_site == "Nigeria")

ord.BC <- ordinate(Phy.f_std_v3, method = "PCoA", distance = "jsd", k=3, trymax=1000) # try PcOA as NMDS does not work
color = c("DeliverVag")
shape = c("CST_pam_3")
#title= c("PCoA of 16S microbiome at birth (Jensen-Shannon distance; k=3)")  
MDS = plot_ordination(Phy.f_std_v2, ord.BC, color = color, shape = shape, #title = title
)
SuppFig.delivery.beta <- MDS+theme_bw() +theme(axis.text=element_text(size=13, face="bold"),
                                         axis.title=element_text(size=16,face="bold"), 
                                         legend.title=element_text(size=14),
                                         legend.text = element_text(size =12))+ 
  labs(color=color)+geom_point(size=3)+
  guides(color = guide_legend(title = "Mode of delivery"))+  # change the ledgend labee title
  guides(shape = guide_legend(title = "Community cluster"))+ # change the shape title
  annotate("text", x= 0.15, y= -0.35, label= "p = 0.485", size = 5)

SuppFig.delivery.beta
ggsave(file = "./figure_output/SuppFig.delivery.beta.pdf", plot = SuppFig.delivery.beta, dpi = 700, width = 7, height = 6) # in this way, I can control the size of the figure
dev.off()


#PAMANOVA  only in Nigeria
set.seed(20)
diss <- phyloseq::distance(Phy.f_std_v3, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_v3_df <- data.frame(row.names=sample_names(Phy.f_std_v3),sample_data(Phy.f_std_v3)) #create dataframe with metadata
# Whether "Last_ViralLoad_cat" differ significantly from each other
adonis2(diss ~ DeliverVag, data=Phy.f_std_v3_df, permutations = 999, by = "terms", na.omit)  # p=0.485


#PAMANOVA  only in Nigeria birth
set.seed(20)
Phy.f_std_v3_W1 <- subset_samples(Phy.f_std_v3, Visit == "Week 1")
diss <- phyloseq::distance(Phy.f_std_v3_W1, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_v3_W1_df <- data.frame(row.names=sample_names(Phy.f_std_v3_W1),sample_data(Phy.f_std_v3_W1)) #create dataframe with metadata
# Whether "Last_ViralLoad_cat" differ significantly from each other
adonis2(diss ~ DeliverVag, data=Phy.f_std_v3_W1_df, permutations = 999, by = "terms", na.omit)  # p=0.566



#PAMANOVA  only in Nigeria W15
set.seed(20)
Phy.f_std_v3_W15 <- subset_samples(Phy.f_std_v3, Visit == "Week 15")
diss <- phyloseq::distance(Phy.f_std_v3_W15, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_v3_W15_df <- data.frame(row.names=sample_names(Phy.f_std_v3_W15),sample_data(Phy.f_std_v3_W15)) #create dataframe with metadata
# Whether "Last_ViralLoad_cat" differ significantly from each other
adonis2(diss ~ DeliverVag, data=Phy.f_std_v3_W15_df, permutations = 999, by = "terms", na.omit)  # p= 0.24




### only Vaginally delivered infants 
levels(sample_data(Phy.f_std)$"DeliverVag") <- c("0"="C-section","1"="Vaginal birth")
Phy.f_std_VagiBirth <- subset_samples(Phy.f_std, DeliverVag == "Vaginal birth")

ord.BC <- ordinate(Phy.f_std_VagiBirth, method = "PCoA", distance = "bray", k=3, trymax=1000) # try PcOA as NMDS does not work
color = c("Study_site")
#title= c("PCoA of 16S microbiome at birth (Jensen-Shannon distance; k=3)")  
MDS = plot_ordination(Phy.f_std_VagiBirth, ord.BC, color = color)
plot <- MDS+theme_bw() +theme(axis.text=element_text(size=13, face="bold"),
                                               axis.title=element_text(size=16,face="bold"), 
                                                strip.text = element_text(size=15, face="bold"), 
                                               legend.title=element_text(size=14),
                                               legend.text = element_text(size =12))+ facet_wrap(~Visit)+
  labs(color=color)+geom_point(size=2)+
  guides(color = guide_legend(title = "Study site"))  # change the ledgend labee title
 
plot
ggsave(file = "./figure_output/SuppFig.VagiBirth.beta.resize.bray.pdf", plot = plot, dpi = 700, width = 11, height = 4) # in this way, I can control the size of the figure
dev.off()

#PAMANOVA  only vaginally delivered infants
set.seed(20)
Phy.f_std_VagiBirth <- subset_samples(Phy.f_std, DeliverVag == "Vaginal birth" & Visit == "Week 1") # cahnge Week 1 or Week 15

diss <- phyloseq::distance(Phy.f_std_VagiBirth, "bray", parallel=TRUE)  
Phy.f_std_VagiBirth_df <- data.frame(row.names=sample_names(Phy.f_std_VagiBirth),sample_data(Phy.f_std_VagiBirth)) #create dataframe with metadata
# Whether "Study_site" differ significantly from each other
adonis2(diss ~ Study_site, data=Phy.f_std_VagiBirth_df, permutations = 999, by = "terms", na.omit)  # p= 0.24



####### Relative abundance
### Simiar to the figure above but with re-ordering the barplot - C-section vs Vaginal birth
##BARPLOT - ORDERED AND CLUSTERED
library(maditr)
Phy.f_std <- readRDS("Phy.f_std.RDS")
sample_data(Phy.f_std)$DeliverVag <- as.factor(sample_data(Phy.f_std)$DeliverVag)
levels(sample_data(Phy.f_std)$"DeliverVag") <- c("0"="C-section","1"="Vaginal birth")

Phy.f_std_v2 <- subset_samples(Phy.f_std, DeliverVag == "C-section" | DeliverVag == "Vaginal birth")
Phy.f_std_v3 <-  subset_samples(Phy.f_std_v2, Study_site == "Nigeria" & Visit == "Week 1")



barplot_stag_BEAMING_sort <- function(physeq) {
  p1 <- tax_glom(physeq, taxrank = 'Species') #agglumerate at species level
  p1_30 = prune_taxa(names(sort(taxa_sums(p1), TRUE))[1:30], p1) #use top 30 most abundant taxa for plotting
  p2 <- transform_sample_counts(p1_30, function(x) x/sum(x)) #get abundance in %
  p3 <- psmelt(p2)  #create dataframe from phyloseq object
  p3$Species <- as.character(p3$Species) #convert to character
  #sort based on abundance of dominant bacteria in each cluster
  p4 <- p3[,c("Sample", "Species","CST_pam_3", "Abundance")]
  p4 <- dcast(p3, Sample + CST_pam_3 ~ Species, value.var = "Abundance", fun.aggregate = sum)
  p4[p4$CST_pam_3 =="Cluster 1",] <- p4[p4$CST_pam_3=="Cluster 1",][order(p4[p4$CST_pam_3=="Cluster 1","coli"]),]
  p4[p4$CST_pam_3 =="Cluster 2",] <- p4[p4$CST_pam_3=="Cluster 2",][order(p4[p4$CST_pam_3=="Cluster 2","longum"]),]
  p4[p4$CST_pam_3 =="Cluster 3",] <- p4[p4$CST_pam_3=="Cluster 3",][order(p4[p4$CST_pam_3=="Cluster 3","faecalis"]),]
  #reorder p3
  p3$Sample <- factor(p3$Sample, levels=c(p4$Sample))
  otus <- otu_table(p1_30)
  tax_table(p1_30)
  #annotation 
  labs <- tax.lab(physeq,otus=otus)
  #set color palette to accommodate the number of species
  colourCount = length(unique(p3$Genus))
  getPalette = colorRampPalette(brewer.pal(8, "Accent"))  
  pal<- c("#7FC97F", "#CBB1C3", "#FDDA8D", "#83A3A7", "#D01487", "#8EC293", 
          "#BB5B19","#DAB6B1", "#FEE992","#FDCA89",  "#EC0877","#C9482C",
          "#9DBBA8", "#E9BA9E","#E3EA9C", "#FEF897", "#4B61AA", "#E01D5E", 
          "#ACB5BC", "#90603F","#F8BE8B", "#77479F", "#D43345", "#7B6352",
          "#BBAED1", "#A65E2C", "#B3C7A1", "#A32D93","#5380AC", "#666666")
  
  pal<- c("steelblue3", "#8EC293", "#FEE992", "azure4", "brown3", 
          "#CBB1C3", "tomato2","lightgreen", "cyan","#FDCA89",  
          "#EC0877","sienna1","lightseagreen", "lightcoral","#E3EA9C", 
          "yellow", "blue3", "#E01D5E", "chartreuse3", "chocolate4",
          "darkgoldenrod1", "#77479F", "honeydew3", "yellow4", "violet", 
          "tan3", "seagreen", "#A32D93","lightskyblue", "#666666")
  #annotation 
  #plot
  barplot_species <- ggplot(data=p3, aes(x=Sample, y=Abundance, fill=Species)) 
  barplot_species <- barplot_species + geom_bar(aes(), stat="identity", position="stack") + 
    facet_wrap(~ DeliverVag, scales = "free") + guides(fill=guide_legend(nrow=3)) + 
    scale_fill_manual(values=pal) + ylab("Relative abundance") +xlab("Participants") + 
    #ggtitle("Bacterial composition at each time pooint by study site") + 
    theme(legend.position="bottom", strip.background = element_rect(fill="lightgrey", color = "black"), 
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          panel.background=element_blank(), panel.border=element_blank(),
          axis.title.y = element_text(size = 17, colour = "black",face = "bold"),
          axis.title.x = element_blank(),
          strip.text = element_text(size=17, face="bold"), 
          axis.text.y=element_text(size=17, face="bold"), 
          text = element_text(size=15, face="bold"))
  barplot_species
  return(barplot_species)
}



SuppFig.delivery.relative.plot <- barplot_stag_BEAMING_sort(Phy.f_std_v3)
SuppFig.delivery.relative.plot

ggsave(file = "./figure_output/SuppFig.delivery.relative.NigW1.plot.pdf", plot = SuppFig.delivery.relative.plot, dpi = 900, width = 20, height = 8) # in this way, I can control the size of the figure
dev.off()


##### Deseq2
Phy.f_std <- readRDS("Phy.f_std.RDS")
sample_data(Phy.f_std)$DeliverVag <- as.factor(sample_data(Phy.f_std)$DeliverVag)
#levels(sample_data(Phy.f_std)$"DeliverVag") <- c("0"="C-section","1"="Vaginal birth")

Phy.f_std_v2 <- subset_samples(Phy.f_std, DeliverVag == "0" | DeliverVag == "1")
Phy.f_std_v3 <-  subset_samples(Phy.f_std_v2, Study_site == "Nigeria" & Visit == "Week 1")
sample_data(Phy.f_std_v3)$DeliverVag  # Levels: 0 1


#Convert the phyloseq object to a DESeqDataSet and run DESeq2:
ds = phyloseq_to_deseq2(Phy.f_std_v3, ~ DeliverVag)
ds = DESeq(ds, fitType = "local")  # get an error...rror in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : every gene contains at least one zero, cannot compute log geometric means
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(ds), 1, gm_mean)
ds = estimateSizeFactors(ds, geoMeans = geoMeans)
ds = DESeq(ds, fitType="local")
res = results(ds)
res = res[order(res$padj, na.last = NA), ]
alpha = 0.01 #this is your significance level, adjust if desired
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Phy.f_std_v3)[rownames(sigtab), ], "matrix"))
head(sigtab)
sigtab %>% View()

resultsNames(ds) 
rownames(sigtab) 

#only show bugs annotated at the Species level
#sigtabspec = subset(sigtab, !is.na(Species))
sigtabgen = subset(sigtab, !is.na(Genus))
#sigtabord = subset(sigtab, !is.na(Order))
head(sigtabgen)

#write.csv(sigtabgen, file = "./figure_output/Sfig-sigtabgen_p0.01_delivery.csv")


#sort for pretty plotting
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

#Alternative: Phylum level, for plotting purposes - don't do both
#x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
#x = sort(x, TRUE)
#sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

#fix taxrank and subset to highly significant asvs
sigtabgen.hs <- sigtabgen[which(sigtabgen$padj < 0.005),]

#subset to greater magnitude of FC
sigtabgen$Species[is.na(sigtabgen$Species)] = "(unclassified)"
sigtabgen$Species[sigtabgen$Species == "NA"] <- "unknown"
sigtabgen$genspec <- paste(sigtabgen$Genus, sigtabgen$Species, sep = " ")
sigtab.merged.hs <- sigtabgen[which(sigtabgen$padj < 0.005),]

sigtabgen %>%
  subset(log2FoldChange > 0.5) -> sigtabgen.hi
sigtabgen %>%
  subset(log2FoldChange < -0.5) -> sigtabgen.low

sigtab.merged.hs %>%
  subset(log2FoldChange > 0.5) -> sigtabgen.hi

sigtab.merged.hs %>%
  subset(log2FoldChange < -0.5) -> sigtabgen.low

sigtabgen.fc <- rbind.data.frame(sigtabgen.hi, sigtabgen.low)                          
View(sigtabgen.fc)

#plot the results (Genus currently set, and colored by Order, set as desired)
SuppFig.delivery.deseq <- ggplot(sigtabgen.fc, aes(y=Genus, x=log2FoldChange, fill=genspec)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_vline(xintercept = 0.5, color = "indianred3", size = 0.5, lty=2) +
  geom_vline(xintercept = -0.5, color = "indianred3", size = 0.5, lty=2) +
  geom_point(pch = 21, size=3)+ theme_bw()+
  theme(axis.text.y = element_text(size = 14, colour = "black",face = "bold"),
        axis.title.y = element_text(size = 17, colour = "black",face = "bold"),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 17, colour = "black",face = "bold"),
        strip.background =element_rect(fill="white"),
        strip.text.x = element_text(colour = "mediumpurple4", size = 15, hjust = 0.5, vjust = 0.5, face = "bold"),
        legend.text = element_text(size = 14, colour = "black"),
        legend.title = element_text(size=15), plot.title = element_text(size = 16)) + 
          theme(legend.position="right", legend.title=element_blank())+
  guides(fill = guide_legend(nrow = 12),
         legend.title = element_text(size=15)) + guides(fill=guide_legend(ncol = 2))+
  labs(x=expression(log[2]~fold~change))+ ggtitle("C-section vs Vaginal birth (Nigeria Week 1)")  #scale_fill_manual(values=custom_colors_36) +

SuppFig.delivery.deseq # something wrong with DEseq - it looks like opposite

ggsave(file = "./figure_output/SuppFig.delivery.deseq.pdf", plot = SuppFig.delivery.deseq, dpi = 900, width = 15, height = 6) # in this way, I can control the size of the figure
dev.off()



#-----------------------------------------------------------------------------------------------------
# Growth related plots ----
#-----------------------------------------------------------------------------------------------------
### Growth transition over time
Phy.f_std <- readRDS("Phy.f_std.RDS")
meta <- tibble(data.frame(sample_data(Phy.f_std)))
names(meta)

meta2 <-meta %>% select(PID,PID2, Visit, Study_site, Status2, Mode_of_feeding,  
                        Mode_of_feeding_birth, Mode_of_feeding_W15,  
                        Feeding_W15_strictver, Feeding_W15_lenientver,  
                        Feeding_Duration_days, FeedingDuration_cat,
                        InfWeight, InfHeight, InfHeadCircum,
                        wflz, wfaz, lfaz, bmiAgeZ, hca) 




### transition of growth score over time (wfaz)
p1 <- ggplot(data=meta2, aes(x=Visit, y=wfaz, color = Study_site))+ theme_bw()+
  geom_boxplot() +geom_jitter(width=0.3, alpha=0.2)+ 
  facet_wrap(~Study_site)+
  labs(x= "Visit") +labs(y= "Weight-for-age Z score") + 
  #ggtitle("Baby Tetanus titer iHEU vs iHUU") + 
  guides(color = guide_legend(title = "Study site"))  + ylim(0, 5) + 
  theme(axis.title=element_text(size=17,face="bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
        strip.text = element_text(size=17, face="bold"), 
        axis.text.y=element_text(size=17, face="bold"),
        legend.text = element_text(size = 17, colour = "black"),
        legend.title = element_text(size = 17, colour = "black", face="bold"))+ 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))


# add stats
a_my_comparisons <- list( c("Week 1", "Week 15"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
Z.score.transition.wfaz <- p1 + stat_compare_means(method = "wilcox.test", 
                                                   comparisons = a_my_comparisons, 
                                                   label = "p.signif", symnum.args = symnum.args
) 

Z.score.transition.wfaz
ggsave(file = "./figure_output/SuppFig.Z.score.transition.wfaz.pdf", plot = Z.score.transition.wfaz, dpi = 900, width = 8, height = 6) # in this way, I can control the size of the figure
dev.off()


### transition of growth score over time (lfaz)
p1 <- ggplot(data=meta2, aes(x=Visit, y=lfaz, color = Study_site))+ theme_bw()+
  geom_boxplot() +geom_jitter(width=0.3, alpha=0.2)+ 
  facet_wrap(~Study_site)+
  labs(x= "Visit") +labs(y= "Length-for-age Z score") + 
  #ggtitle("Baby Tetanus titer iHEU vs iHUU") + 
  guides(color = guide_legend(title = "Study site"))  + ylim(0, 5) + 
  theme(axis.title=element_text(size=17,face="bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
        strip.text = element_text(size=17, face="bold"), 
        axis.text.y=element_text(size=17, face="bold"),
        legend.text = element_text(size = 17, colour = "black"),
        legend.title = element_text(size = 17, colour = "black", face="bold"))+ 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))


# add stats
a_my_comparisons <- list( c("Week 1", "Week 15"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
Z.score.transition.lfaz <- p1 + stat_compare_means(method = "wilcox.test", 
                                                   comparisons = a_my_comparisons, 
                                                   label = "p.signif", symnum.args = symnum.args
) 

Z.score.transition.lfaz
ggsave(file = "./figure_output/SuppFig.Z.score.transition.lfaz.pdf", plot = Z.score.transition.lfaz, dpi = 900, width = 8, height = 6) # in this way, I can control the size of the figure
dev.off()

### transition of growth score over time (wflz)
p1 <- ggplot(data=meta2, aes(x=Visit, y=wflz, color = Study_site))+ theme_bw()+
  geom_boxplot() +geom_jitter(width=0.3, alpha=0.2)+ 
  facet_wrap(~Study_site)+
  labs(x= "Visit") +labs(y= "Weight-for-length Z score") + 
  #ggtitle("Baby Tetanus titer iHEU vs iHUU") + 
  guides(color = guide_legend(title = "Study site"))  + #ylim(0, 7.0) + 
  theme(axis.title=element_text(size=17,face="bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
        strip.text = element_text(size=17, face="bold"), 
        axis.text.y=element_text(size=17, face="bold"),
        legend.text = element_text(size = 17, colour = "black"),
        legend.title = element_text(size = 17, colour = "black", face="bold"))+ 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))


# add stats
a_my_comparisons <- list( c("Week 1", "Week 15"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
Z.score.transition.wflz <- p1 + stat_compare_means(method = "wilcox.test", 
                                                 comparisons = a_my_comparisons, 
                                                 label = "p.signif", symnum.args = symnum.args
) 

Z.score.transition.wflz
ggsave(file = "./figure_output/SuppFig.Z.score.transition.wflz.pdf", plot = Z.score.transition.wflz, dpi = 900, width = 8, height = 6) # in this way, I can control the size of the figure
dev.off()








### Mode of breastfeeding vs growth
Phy.f_std <- readRDS("Phy.f_std.RDS")
meta <- tibble(data.frame(sample_data(Phy.f_std)))
names(meta)

meta2 <-meta %>% select(PID, Visit, Study_site, Status2, Mode_of_feeding,  
                               Mode_of_feeding_birth, Mode_of_feeding_W15,  
                               Feeding_W15_strictver, Feeding_W15_lenientver,  
                               Feeding_Duration_days, FeedingDuration_cat,
                              InfWeight, InfHeight, InfHeadCircum,
                               wflz, wfaz, lfaz, bmiAgeZ, hca,DeliverVag) 

meta2_SA_W15<- meta2  %>% filter(Visit == "Week 15" & Study_site == "South Africa")
length(meta2_SA_W15$Feeding_W15_lenientver[meta2_SA_W15$Feeding_W15_lenientver == "Breastfeeding"]) #13
length(meta2_SA_W15$Feeding_W15_lenientver[meta2_SA_W15$Feeding_W15_lenientver == "Exclusive_breastfeeding"]) #40
length(meta2_SA_W15$Feeding_W15_lenientver[meta2_SA_W15$Feeding_W15_lenientver == "Switched_to_formula"]) #13

# Nigeria - all exclusively breastfed

meta2_W15 <- meta2  %>% filter(Visit == "Week 15")
meta2_W15$Feeding_W15_lenientver
meta2_W15$FeedingDuration_cat
meta2_W15$InfWeight <- as.numeric(meta2_W15$InfWeight)
levels(meta2_W15$Feeding_W15_strictver) <- c("Exclusive_breastfeeding"="Exclusive breastfeeding","Mixed_feeding"="Mixed feeding")

### weight atW15
p1 <- ggplot(data=meta2_W15, aes(x=Study_site, y=InfWeight, color = Feeding_W15_strictver))+ theme_bw()+
  geom_boxplot() +geom_jitter(width=0.3, alpha=0.2)+ 
  #facet_wrap(~Feeding_W15_strictver,scales = "free_y")+
  labs(x= "Study site") +labs(y= "Weight at 15 weeks of age") + 
  #ggtitle("Baby Tetanus titer iHEU vs iHUU") + 
  guides(color = guide_legend(title = "Mode of feeding"))  + #ylim(0, 7.0) + 
  theme(axis.title=element_text(size=17,face="bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
        strip.text = element_text(size=17, face="bold"), 
        axis.text.y=element_text(size=17, face="bold"),
        text = element_text(size=17, face="bold"))+ 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# add stats
a_my_comparisons <- list( c("South Africa", "Nigeria"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
growth.weight.boxplot <- p1 + stat_compare_means(method = "wilcox.test", 
                                               comparisons = a_my_comparisons, 
                                               label = "p.signif", symnum.args = symnum.args
) 

growth.weight.boxplot

ggsave(file = "./figure_output/SuppFig.weight.feeding.boxplot.pdf", plot = growth.weight.boxplot, dpi = 500, width = 8, height = 7) # in this way, I can control the size of the figure
dev.off()

### Height atW15
p1 <- ggplot(data=meta2_W15, aes(x=Study_site, y=InfHeight, color = Feeding_W15_strictver))+ theme_bw()+
  geom_boxplot() +geom_jitter(width=0.3, alpha=0.2)+ 
  #facet_wrap(~Feeding_W15_strictver,scales = "free_y")+
  labs(x= "Study site") +labs(y= "Height at 15 weeks of age") + 
  #ggtitle("Baby Tetanus titer iHEU vs iHUU") + 
  guides(color = guide_legend(title = "Mode of feeding"))  + #ylim(0, 7.0) + 
  theme(axis.title=element_text(size=17,face="bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
        strip.text = element_text(size=17, face="bold"), 
        axis.text.y=element_text(size=17, face="bold"),
        text = element_text(size=17, face="bold"))+ 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# add stats
a_my_comparisons <- list( c("South Africa", "Nigeria"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
growth.height.boxplot <- p1 + stat_compare_means(method = "wilcox.test", 
                                                 comparisons = a_my_comparisons, 
                                                 label = "p.signif", symnum.args = symnum.args
) 

growth.height.boxplot

ggsave(file = "./figure_output/SuppFig.height.feeding.boxplot.pdf", plot = growth.height.boxplot, dpi = 500, width = 8, height = 7) # in this way, I can control the size of the figure
dev.off()






### length-for-age (lfa)
p1 <- ggplot(data=meta2_W15, aes(x=Study_site, y=lfaz, color = Feeding_W15_strictver))+ theme_bw()+
  geom_boxplot() +geom_jitter(width=0.3, alpha=0.2)+ 
  #facet_wrap(~Feeding_W15_strictver,scales = "free_y")+
  labs(x= "Study site") +labs(y= "Length-for-age z score") + 
  #ggtitle("Baby Tetanus titer iHEU vs iHUU") + 
  guides(color = guide_legend(title = "Mode of feeding"))  + #ylim(0, 7.0) + 
  theme(axis.title=element_text(size=17,face="bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
        strip.text = element_text(size=17, face="bold"), 
        axis.text.y=element_text(size=17, face="bold"),
        text = element_text(size=17, face="bold"))+ 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_segment(aes(x=0.8,xend=1.2,y=2.5,yend=2.5), color = "black", size = 0.15)+
  annotate("text", x= 1, y= 2.8, label= "ns", size = 4)# p=xxx

# add stats
a_my_comparisons <- list( c("South Africa", "Nigeria"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
growth.lfaz.boxplot <- p1 + stat_compare_means(method = "wilcox.test", 
                                               comparisons = a_my_comparisons, tip.length = 0,
                                               label = "p.signif", symnum.args = symnum.args
) 

growth.lfaz.boxplot

ggsave(file = "./figure_output/SuppFig.lfaz.feeding.boxplot.pdf", plot = growth.lfaz.boxplot, dpi = 500, width = 8, height = 7) # in this way, I can control the size of the figure
dev.off()




### Weight-for-age (wfa)

p1 <- ggplot(data=meta2_W15, aes(x=Study_site, y=wfaz, color = Feeding_W15_strictver))+ theme_bw()+
  geom_boxplot() +geom_jitter(width=0.3, alpha=0.2)+ 
  #facet_wrap(~Feeding_W15_strictver,scales = "free_y")+
  labs(x= "Study site") +labs(y= "Weight-for-age z score") + 
  #ggtitle("Baby Tetanus titer iHEU vs iHUU") + 
  guides(color = guide_legend(title = "Mode of feeding"))  + #ylim(0, 7.0) + 
  theme(axis.title=element_text(size=17,face="bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
        strip.text = element_text(size=17, face="bold"), 
        axis.text.y=element_text(size=17, face="bold"),
        text = element_text(size=17, face="bold"))+ 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  geom_segment(aes(x=0.8,xend=1.2,y=3.1,yend=3.1), color = "black", size = 0.15)+
  annotate("text", x= 1, y= 3.3, label= "ns", size = 4)# p=xxx

# add stats
a_my_comparisons <- list( c("South Africa", "Nigeria"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
growth.wfaz.boxplot <- p1 + stat_compare_means(method = "wilcox.test", 
                                  comparisons = a_my_comparisons, tip.length = 0,
                                  label = "p.signif", symnum.args = symnum.args
) 

growth.wfaz.boxplot

ggsave(file = "./figure_output/SuppFig.wfaz.feeding.boxplot.pdf", plot = growth.wfaz.boxplot, dpi = 500, width = 8, height = 7) # in this way, I can control the size of the figure
dev.off()


# wflz
### Weight-for-length (wfl)

p1 <- ggplot(data=meta2_W15, aes(x=Study_site, y=wflz, color = Feeding_W15_strictver))+ theme_bw()+
  geom_boxplot() +geom_jitter(width=0.3, alpha=0.2)+ 
  #facet_wrap(~Feeding_W15_strictver,scales = "free_y")+
  labs(x= "Study site") +labs(y= "Weight-for-length z score") + 
  #ggtitle("Baby Tetanus titer iHEU vs iHUU") + 
  guides(color = guide_legend(title = "Mode of feeding"))  + #ylim(0, 7.0) + 
  theme(axis.title=element_text(size=17,face="bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
        strip.text = element_text(size=17, face="bold"), 
        axis.text.y=element_text(size=17, face="bold"),
        text = element_text(size=17, face="bold"))+ 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  geom_segment(aes(x=0.8,xend=1.2,y=5,yend=5), color = "black", size = 0.15)+
  annotate("text", x= 1, y= 5.3, label= "*", size = 4)# p=0.045


# add stats
a_my_comparisons <- list( c("South Africa", "Nigeria"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
growth.wflz.boxplot <- p1 + stat_compare_means(method = "wilcox.test", tip.length = 0,
                                               comparisons = a_my_comparisons, 
                                               label = "p.signif", symnum.args = symnum.args
) 

growth.wflz.boxplot

ggsave(file = "./figure_output/SuppFig.wflz.feeding.boxplot.pdf", plot = growth.wflz.boxplot, dpi = 500, width = 8, height = 7) # in this way, I can control the size of the figure
dev.off()



### let's combine these three plots
library(ggpubr)

growth.wfaz.boxplot
growth.wflz.boxplot
growth.lfaz.boxplot

# not work nicely as the y-axis are different...
feedinv.vs.growth.plots <- ggarrange(growth.wfaz.boxplot,growth.wflz.boxplot, growth.lfaz.boxplot +rremove("x.text"),
                                     common.legend = TRUE,labels = c(), nrow = 1)


### growth score vs mode of delivery
meta2$DeliverVag
levels(meta2$"DeliverVag") <- c("0"="C-section","1"="Vaginal birth")

### Weight-for-length (wfl)
p1 <- ggplot(data=meta2, aes(x=Study_site, y=wflz, color = DeliverVag))+ theme_bw()+
  geom_boxplot() +geom_jitter(width=0.3, alpha=0.2)+ 
  #facet_wrap(~Feeding_W15_strictver,scales = "free_y")+
  labs(x= "Study site") +labs(y= "Weight-for-length z score") + 
  #ggtitle("Baby Tetanus titer iHEU vs iHUU") + 
  guides(color = guide_legend(title = "Mode of delivery"))  + #ylim(0, 7.0) + 
  theme(axis.title=element_text(size=17,face="bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
        strip.text = element_text(size=17, face="bold"), 
        axis.text.y=element_text(size=17, face="bold"),
        text = element_text(size=17, face="bold"))+ 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  geom_segment(aes(x=1.8,xend=2.2,y=9,yend=9), color = "black", size = 0.15)+
  annotate("text", x= 2, y= 9.2, label= "**", size = 4)# p=xxx

# add stats
a_my_comparisons <- list( c("South Africa", "Nigeria"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
delivery.wflz.boxplot <- p1 + stat_compare_means(method = "wilcox.test", tip.length = 0,
                                               comparisons = a_my_comparisons, 
                                               label = "p.signif", symnum.args = symnum.args
) 

delivery.wflz.boxplot

ggsave(file = "./figure_output/SuppFig.wflz.delivery.boxplot.pdf", plot = delivery.wflz.boxplot, dpi = 500, width = 8, height = 7) # in this way, I can control the size of the figure
dev.off()



### Weight-for-length (wfa)
p1 <- ggplot(data=meta2, aes(x=Study_site, y=wfaz, color = DeliverVag))+ theme_bw()+
  geom_boxplot() +geom_jitter(width=0.3, alpha=0.2)+ 
  #facet_wrap(~Feeding_W15_strictver,scales = "free_y")+
  labs(x= "Study site") +labs(y= "Weight-for-zge z score") + 
  #ggtitle("Baby Tetanus titer iHEU vs iHUU") + 
  guides(color = guide_legend(title = "Mode of delivery"))  + #ylim(0, 7.0) + 
  theme(axis.title=element_text(size=17,face="bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
        strip.text = element_text(size=17, face="bold"), 
        axis.text.y=element_text(size=17, face="bold"),
        text = element_text(size=17, face="bold"))+ 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  geom_segment(aes(x=1.8,xend=2.2,y=9,yend=9), color = "black", size = 0.15)+
  annotate("text", x= 2, y= 9.3, label= "ns", size = 4)# p=xxx

# add stats
a_my_comparisons <- list( c("South Africa", "Nigeria"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
delivery.wfaz.boxplot <- p1 + stat_compare_means(method = "wilcox.test", tip.length = 0,
                                                 comparisons = a_my_comparisons, 
                                                 label = "p.signif", symnum.args = symnum.args
) 

delivery.wfaz.boxplot

ggsave(file = "./figure_output/SuppFig.wfaz.delivery.boxplot.pdf", plot = delivery.wfaz.boxplot, dpi = 500, width = 8, height = 7) # in this way, I can control the size of the figure
dev.off()




### Length-for-length (lfa)
p1 <- ggplot(data=meta2, aes(x=Study_site, y=lfaz, color = DeliverVag))+ theme_bw()+
  geom_boxplot() +geom_jitter(width=0.3, alpha=0.2)+ 
  #facet_wrap(~Feeding_W15_strictver,scales = "free_y")+
  labs(x= "Study site") +labs(y= "Length-for-zge z score") + 
  #ggtitle("Baby Tetanus titer iHEU vs iHUU") + 
  guides(color = guide_legend(title = "Mode of delivery"))  + #ylim(0, 7.0) + 
  theme(axis.title=element_text(size=17,face="bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
        strip.text = element_text(size=17, face="bold"), 
        axis.text.y=element_text(size=17, face="bold"),
        text = element_text(size=17, face="bold"))+ 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  geom_segment(aes(x=1.8,xend=2.2,y=9,yend=9), color = "black", size = 0.15)+
  annotate("text", x= 2, y= 9.3, label= "*", size = 4)# p=xxx

# add stats
a_my_comparisons <- list( c("South Africa", "Nigeria"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
delivery.lfaz.boxplot <- p1 + stat_compare_means(method = "wilcox.test", tip.length = 0,
                                                 comparisons = a_my_comparisons, 
                                                 label = "p.signif", symnum.args = symnum.args
) 

delivery.lfaz.boxplot

ggsave(file = "./figure_output/SuppFig.lfaz.delivery.boxplot.pdf", plot = delivery.lfaz.boxplot, dpi = 500, width = 8, height = 7) # in this way, I can control the size of the figure
dev.off()


#### growth vs microbiota
my_data <- tibble(data.frame(sample_data(Phy.f_std))) %>% select(PID, Study_site, Visit, wfaz, wflz) 
my_data

SA_birth <- my_data %>% filter(Study_site =="South Africa" & Visit == "Birth")
SA_W15 <- my_data %>% filter(Study_site =="South Africa" & Visit == "W15")
median(SA_birth$wfaz) # -0.39
barplot(SA_birth$wfaz)
barplot(SA_W15$wfaz)

Nigeria_birth <- my_data %>% filter(Study_site =="Nigeria" & Visit == "Birth")
Nigeria_W15 <- my_data %>% filter(Study_site =="Nigeria" & Visit == "W15")

median(Nigeria_birth$wfaz) # -0.52
barplot(Nigeria_birth$wfaz)
barplot(Nigeria_W15$wfaz)


my_data.v2 <- my_data %>% mutate(wfaz_groups = 
                                   case_when(
                                     (wfaz > 2 ) ~ "Overweight",
                                     (wfaz <= 2 & wfaz >= -2) ~ "Normal",
                                     (wfaz < -2) ~ "Wasting")) 


my_data.v2 <- my_data %>% mutate(wfaz_groups = 
                                   case_when(
                                     (wfaz  >= -2) ~ "Not wasting",
                                     (wfaz < -2) ~ "Wasting")) %>% 
  mutate(wflz_groups = 
           case_when(
             (wflz  >= -2) ~ "Not wasting",
             (wflz < -2) ~ "Wasting"))
  
  



my_data.v2$wfaz_groups <- as.factor(my_data.v2$wfaz_groups)
my_data.v2$wflz_groups <- as.factor(my_data.v2$wflz_groups)

sample_data(Phy.f_std)$wfaz_groups <- my_data.v2$wfaz_groups  # add the annotation into the phyloseqobject
sample_data(Phy.f_std)$wflz_groups <- my_data.v2$wflz_groups 

tibble(data.frame(sample_data(Phy.f_std))) %>% 
  select(PID, Study_site, Visit, lfaz, wfaz, wfaz_groups, wflz,wflz_groups) %>% View()


### alpha diversity with "wfaz_groups" ---- not sure why there are some missing values...
p1 <- plot_richness(Phy.f_std, x= "wfaz_groups", measures="Shannon",
                    # title = "Alpha Diversity at birth"
)+xlab("Nutritional Status (wfaz)")+
  geom_boxplot(alpha=0.6)+ facet_wrap(Study_site ~ Visit)+
  theme_bw()+ theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
                    strip.text = element_text(size=17, face="bold"), axis.text.y=element_text(size=16, face="bold"),
                    text = element_text(size=17, face="bold"))+ scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# and add statistics
# https://rpubs.com/lconteville/713954
a_my_comparisons <- list( c("Not wasting", "Wasting"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
Nutritional.wfaz <- p1 + stat_compare_means(method = "wilcox.test", 
                                comparisons = a_my_comparisons, 
                                label = "p.signif", symnum.args = symnum.args) 
Nutritional.wfaz  #significantly different - CT increases the alpha diversity over the time whereas decased in Nigerian infants


ggsave(file = "./figure_output/Supp_Nutritional.wfaz.pdf", plot = Nutritional.wfaz, dpi = 800, width = 7, height = 7) # in this way, I can control the size of the figure
dev.off()





### alpha diversity with "wflz_groups" ----
Phy.f_std_v2 <- subset_samples(Phy.f_std, wflz_groups == "Not wasting"|  wflz_groups == "Wasting")
p1 <- plot_richness(Phy.f_std_v2, x= "wflz_groups", measures="Shannon",
                    # title = "Alpha Diversity at birth"
)+xlab("Nutritional Status (wflz)")+
  geom_boxplot(alpha=0.6)+ facet_wrap(Study_site ~ Visit)+
  theme_bw()+ theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
                    strip.text = element_text(size=17, face="bold"), axis.text.y=element_text(size=16, face="bold"),
                    text = element_text(size=17, face="bold"))+ scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# and add statistics
# https://rpubs.com/lconteville/713954
a_my_comparisons <- list( c("Not wasting", "Wasting"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
Nutritional.wflz <- p1 + stat_compare_means(method = "wilcox.test", 
                                            comparisons = a_my_comparisons, 
                                            label = "p.signif", symnum.args = symnum.args) 
Nutritional.wflz  #significantly different - CT increases the alpha diversity over the time whereas decased in Nigerian infants


ggsave(file = "./figure_output/Supp_Nutritional.wflz.pdf", plot = Nutritional.wflz, dpi = 800, width = 7, height = 7) # in this way, I can control the size of the figure
dev.off()

### Scattter Plot with continuous variale with loop function - automatically generate & save plots
my_variable <- c("wfaz", "lfaz", "wflz")

for(variable in my_variable){
  
  df <- tibble(data.frame(sample_data(Phy.f_std)))
  
  final.plot <-   ggscatter(data=df, x = variable, y = "Shannon", 
                            add = "reg.line", conf.int = TRUE, size=2) +
    stat_cor(method = "pearson", size = 5.5)+ 
    facet_wrap(Study_site~ Visit, nrow=2, scales = "free_y")  + theme_bw()+ 
    theme(axis.text=element_text(size=17,face = "bold"),
          axis.title=element_text(size=17,face="bold"), 
          strip.text = element_text(size=17, face="bold")) +ylim(0,5.5)+ 
    geom_vline(xintercept = -2, color = "red", size=0.3)+
    labs(x= variable) +labs(y= "Shannon diversity")
  
  pdf(paste0("./figure_output/Supp_Nutritional.",variable, ".scatter.pdf"))
  print(final.plot)
  dev.off()
}





### Beta diversirty ----- wfaz_groups
levels(sample_data(Phy.f_std)$"Study_site") <- c("CT"="South Africa","Nigeria"="Nigeria")

set.seed(20) #select any start seed, in this case 2 (ensures reproducibility)
ord.BC <- ordinate(Phy.f_std, method = "PCoA", distance = "jsd", k=2, trymax=1000) # try PcOA as NMDS does not work
ord.BC

color = c("wfaz_groups")
shape = c("Study_site")
# title= c("PCoA of 16S microbiome (all visit), Jensen-Shannon distance, k=2")  
MDS = plot_ordination(Phy.f_std, ord.BC, color = color, shape=shape #, title = title
                      )
Nutritional.wfaz.beta = MDS +theme_bw() +theme(axis.text=element_text(size=13, face="bold"),
                                            axis.title=element_text(size=16,face="bold"), 
                                            legend.title=element_text(size=14),
                                            legend.text = element_text(size =12))+ 
  labs(color=color, shape=shape)+ geom_point(size=3)+ 
  guides(color = guide_legend(title = "Nutritional status (wfaz)"))+  # change the ledgend labee title
  guides(shape = guide_legend(title = "Study site")) # change the shape title
  #annotate("text", x= 0.15, y= -0.35, label= "p = XXX", size = 5)


Nutritional.wfaz.beta
ggsave(file = "./figure_output/SuppFig.Nutritional.wfaz.beta.pdf", plot = Nutritional.wfaz.beta, dpi = 800, width = 8, height = 6) # in this way, I can control the size of the figure
dev.off()



### beta-diversity with gradient
Phy.f_std_v2 <- subset_samples(Phy.f_std, Visit == "Week1")

ord.BC <- ordinate(Phy.f_std_v2, method = "PCoA", distance = "jsd", k=3, trymax=1000) # try PcOA as NMDS does not work
color = c("wfaz")
shape = c("Study_site")
MDS = plot_ordination(Phy.f_std_v2, ord.BC, color = color, shape = shape)
Nutritional.wfaz.beta.gradient  = MDS+theme_bw() +
  scale_color_viridis_c(option = "plasma")+
  theme(axis.text=element_text(size=17, face="bold"),
        axis.title=element_text(size=17,face="bold"), 
        legend.title=element_text(size=17 ,face="bold"),
        legend.text = element_text(size =16))+
  labs(color= "Weight-for-age")+geom_point(size=3)+
  guides(shape = guide_legend(title = "Study site")) # change the shape title

Nutritional.wfaz.beta.gradient

Nutritional.wfaz.beta.gradient
ggsave(file = "./figure_output/SuppFig.Nutritional.wfaz.beta.gradient.pdf", plot = Nutritional.wfaz.beta.gradient, dpi = 800, width = 8, height = 6) # in this way, I can control the size of the figure
dev.off()



### adonis (1) overall
set.seed(2)
sample_data(Phy.f_std)$wfaz_groups

diss <- phyloseq::distance(Phy.f_std, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_df <- data.frame(row.names=sample_names(Phy.f_std),sample_data(Phy.f_std)) #create dataframe with metadata

# Whether "wfaz_groups" differ significantly from each other
adonis2(diss ~ wfaz_groups, data=Phy.f_std_df, permutations = 999, by = "terms", na.omit)  # 0.022. significantly different...
adonis2(diss ~ wfaz_groups + Visit, data=Phy.f_std_df, permutations = 999, by = "terms", na.omit)  # 0.022. significantly different...
adonis2(diss ~ wfaz_groups * Visit, data=Phy.f_std_df, permutations = 999, by = "terms", na.omit)  # * means covariate
adonis2(diss ~ wfaz_groups * Visit * Study_site, data=Phy.f_std_df, permutations = 999, by = "terms", na.omit)  # * means covariate

adonis2(diss ~ wfaz, data=Phy.f_std_df, permutations = 999, by = "terms", na.omit)  # * means covariate



### adonis (2) only birth
sample_data(Phy.f_std)$wfaz_groups
Phy.f_std_birth <- subset_samples(Phy.f_std, Visit == "Week 1") 

diss <- phyloseq::distance(Phy.f_std_birth, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_birth_df <- data.frame(row.names=sample_names(Phy.f_std_birth),sample_data(Phy.f_std_birth)) #create dataframe with metadata
adonis2(diss ~ wfaz_groups, data=Phy.f_std_birth_df, permutations = 999, by = "terms", na.omit) #p=0.477

adonis2(diss ~ wfaz, data=Phy.f_std_birth_df, permutations = 999, by = "terms", na.omit) #p=0.477



### adonis (2) only W15
sample_data(Phy.f_std)$wfaz_groups
Phy.f_std_W15 <- subset_samples(Phy.f_std, Visit == "Week 15") 

diss <- phyloseq::distance(Phy.f_std_W15, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_W15_df <- data.frame(row.names=sample_names(Phy.f_std_W15),sample_data(Phy.f_std_W15)) #create dataframe with metadata
adonis2(diss ~ wfaz_groups, data=Phy.f_std_W15_df, permutations = 999, by = "terms", na.omit) #p=0.477





### Beta diversirty ----- wflz_groups
### Beta-diversity (3') PCoA plot comparing by site (include both birth and W15) - no annotation of hiv exposure buy annotate with TP - jsd - fig 2 (C)
sample_data(Phy.f_std)$wflz_groups  # Levels: Not wasting Wasting

sample_data(Phy.f_std)$"Study_site"
Phy.f_std_v2 <- subset_samples(Phy.f_std, wflz_groups == "Not wasting" | wflz_groups == "Wasting")


set.seed(20) #select any start seed, in this case 2 (ensures reproducibility)
ord.BC <- ordinate(Phy.f_std_v2, method = "PCoA", distance = "jsd", k=2, trymax=1000) # try PcOA as NMDS does not work
ord.BC

color = c("wflz_groups")
shape = c("Study_site")
# title= c("PCoA of 16S microbiome (all visit), Jensen-Shannon distance, k=2")  
MDS = plot_ordination(Phy.f_std_v2, ord.BC, color = color, shape=shape #, title = title
)
Nutritional.wflz.beta = MDS +theme_bw() +theme(axis.text=element_text(size=13, face="bold"),
                                               axis.title=element_text(size=16,face="bold"), 
                                               legend.title=element_text(size=14),
                                               legend.text = element_text(size =12))+ 
  labs(color=color, shape=shape)+ geom_point(size=3)+ 
  guides(color = guide_legend(title = "Nutritional status (wflz)"))+  # change the ledgend labee title
  guides(shape = guide_legend(title = "Study site")) # change the shape title
#annotate("text", x= 0.15, y= -0.35, label= "p = XXX", size = 5)


Nutritional.wflz.beta
ggsave(file = "./figure_output/SuppFig.Nutritional.wflz.beta.pdf", plot = Nutritional.wflz.beta, dpi = 800, width = 8, height = 6) # in this way, I can control the size of the figure
dev.off()


### adonis (1) overall
set.seed(2)
sample_data(Phy.f_std_v2)$wflz_groups

diss <- phyloseq::distance(Phy.f_std_v2, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_v2_df <- data.frame(row.names=sample_names(Phy.f_std_v2),sample_data(Phy.f_std_v2)) #create dataframe with metadata

# Whether "wflz_groups" differ significantly from each other
adonis2(diss ~ wflz_groups, data=Phy.f_std_v2_df, permutations = 999, by = "terms", na.omit)  # 0.022. significantly different...
adonis2(diss ~ wflz_groups + Visit, data=Phy.f_std_v2_df, permutations = 999, by = "terms", na.omit)  # 0.022. significantly different...
adonis2(diss ~ wflz_groups * Visit, data=Phy.f_std_v2_df, permutations = 999, by = "terms", na.omit)  # * means covariate
adonis2(diss ~ wflz_groups * Visit * Study_site, data=Phy.f_std_v2_df, permutations = 999, by = "terms", na.omit)  # * means covariate




### Beta-diversity with loop function - automatically save
# for W1 dataset
Phy.f_std_v2 <- subset_samples(Phy.f_std, Visit == "Week 1")
my_variable <- c("wfaz", "lfaz", "wflz")

for(variable in my_variable){
  
  ord.BC <- ordinate(Phy.f_std_v2, method = "PCoA", distance = "jsd", k=3, trymax=1000) # try PcOA as NMDS does not work
  final.plot = plot_ordination(Phy.f_std_v2, ord.BC, color = variable, shape = "Study_site")+ theme_bw() +
    scale_color_viridis_c(option = "plasma")+
    theme(axis.text=element_text(size=17, face="bold"),
          axis.title=element_text(size=17,face="bold"), 
          legend.title=element_text(size=17 ,face="bold"),
          legend.text = element_text(size =16))+
    labs(color= variable)+geom_point(size=3)+
    guides(shape = FALSE) # turn off the shape ledend
  
  
  ggsave(file = paste0("./figure_output/Supp_Nutritional.W1.",variable, ".beta.pdf"), plot = final.plot, dpi = 800, width = 7, height = 6) # in this way, I can control the size of the figure
}


# for W15 dataset
Phy.f_std_v3 <- subset_samples(Phy.f_std, Visit == "Week 15")
my_variable <- c("wfaz", "lfaz", "wflz")

for(variable in my_variable){
  
  ord.BC <- ordinate(Phy.f_std_v3, method = "PCoA", distance = "jsd", k=3, trymax=1000) # try PcOA as NMDS does not work
  final.plot = plot_ordination(Phy.f_std_v3, ord.BC, color = variable, shape = "Study_site")+ theme_bw() +
    scale_color_viridis_c(option = "plasma")+
    theme(axis.text=element_text(size=17, face="bold"),
          axis.title=element_text(size=17,face="bold"), 
          legend.title=element_text(size=17 ,face="bold"),
          legend.text = element_text(size =16))+
    labs(color= variable)+geom_point(size=3)+
    guides(shape = FALSE) # turn off the shape ledend
  
  
  ggsave(file = paste0("./figure_output/Supp_Nutritional.W15.",variable, ".beta.pdf"), plot = final.plot, dpi = 800, width = 7, height = 6) # in this way, I can control the size of the figure
}


### adonis (2) only birth
set.seed(2)
sample_data(Phy.f_std)$wfaz_groups
Phy.f_std_birth <- subset_samples(Phy.f_std, Visit == "Week 1") 
Phy.f_std_birth <- subset_samples(Phy.f_std_birth, wflz >0 |wflz <= 0 ) 

diss <- phyloseq::distance(Phy.f_std_birth, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_birth_df <- data.frame(row.names=sample_names(Phy.f_std_birth),sample_data(Phy.f_std_birth)) #create dataframe with metadata

adonis2(diss ~ wfaz, data=Phy.f_std_birth_df, permutations = 999, by = "terms", na.omit) #p=0.377/ R2: 0.00511
adonis2(diss ~ lfaz, data=Phy.f_std_birth_df, permutations = 999, by = "terms", na.omit) #p=0.282/ R2: 0.0056
adonis2(diss ~ wflz, data=Phy.f_std_birth_df, permutations = 999, by = "terms", na.omit) #p = 0.646/ R2: 0.0048



### adonis (2) only W15
sample_data(Phy.f_std)$wfaz_groups
Phy.f_std_W15 <- subset_samples(Phy.f_std, Visit == "Week 15") 
Phy.f_std_W15 <- subset_samples(Phy.f_std_W15, wflz >0 |wflz <= 0 ) 

set.seed(2)
diss <- phyloseq::distance(Phy.f_std_W15, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_W15_df <- data.frame(row.names=sample_names(Phy.f_std_W15),sample_data(Phy.f_std_W15)) #create dataframe with metadata

adonis2(diss ~ wfaz, data=Phy.f_std_W15_df, permutations = 999, by = "terms", na.omit) #p= 0.001/ R2: 0.02503
adonis2(diss ~ lfaz, data=Phy.f_std_W15_df, permutations = 999, by = "terms", na.omit) #p= 0.065 / R2: 0.00944
adonis2(diss ~ wflz, data=Phy.f_std_W15_df, permutations = 999, by = "terms", na.omit) #p= 0.23 / R2: 0.00539



### adonis (2) only W15 & Nigeria
sample_data(Phy.f_std)$wfaz_groups
Phy.f_std_W15 <- subset_samples(Phy.f_std, Visit == "Week 15" & Study_site == "Nigeria") 
#Phy.f_std_W15 <- subset_samples(Phy.f_std_W15, wflz >0 |wflz <= 0 ) 

set.seed(2)
diss <- phyloseq::distance(Phy.f_std_W15, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_W15_df <- data.frame(row.names=sample_names(Phy.f_std_W15),sample_data(Phy.f_std_W15)) #create dataframe with metadata

adonis2(diss ~ wfaz, data=Phy.f_std_W15_df, permutations = 999, by = "terms", na.omit) #p= 0.477/ R2: 0.00532
adonis2(diss ~ lfaz, data=Phy.f_std_W15_df, permutations = 999, by = "terms", na.omit) #p= 0.705/ R2: 0.00344
adonis2(diss ~ wflz, data=Phy.f_std_W15_df, permutations = 999, by = "terms", na.omit) #p= 0.505/ R2: 0.00472


### adonis (2) only W15 & SA
sample_data(Phy.f_std)$wfaz_groups
Phy.f_std_W15 <- subset_samples(Phy.f_std, Visit == "Week 15" & Study_site == "South Africa") 
#Phy.f_std_W15 <- subset_samples(Phy.f_std_W15, wflz >0 |wflz <= 0 ) 

set.seed(2)
diss <- phyloseq::distance(Phy.f_std_W15, "jsd", parallel=TRUE)  # jsd = Jensen–Shannon divergence 
Phy.f_std_W15_df <- data.frame(row.names=sample_names(Phy.f_std_W15),sample_data(Phy.f_std_W15)) #create dataframe with metadata

adonis2(diss ~ wfaz, data=Phy.f_std_W15_df, permutations = 999, by = "terms", na.omit) #p= 0.405/ R2: 0.01582
adonis2(diss ~ lfaz, data=Phy.f_std_W15_df, permutations = 999, by = "terms", na.omit) #p= 0.121/ R2: 0.02242
adonis2(diss ~ wflz, data=Phy.f_std_W15_df, permutations = 999, by = "terms", na.omit) #p= 0.7/ R2: 0.01185




df <- tibble(data.frame(sample_data(Phy.f_std))) %>% select(PID, Study_site, Status2, Visit,
                                                            Diarrhea_byW15,Diarrhea_days, Cough_byW15, 
                                                            Vitamin_byW15, Antibiotics_ever) %>% View()
names(df)


########## Alpha diversity & different measures
### FigS (alpha diversity)---- corresponds to Fig 1(B)
Phy.f_std <- readRDS("Phy.f_std.RDS")
Phy.f_std_birth <- subset_samples(Phy.f_std, Visit =="Week 1")
levels(sample_data(Phy.f_std_birth)$"Study_site") <- c("CT"="South Africa","Nigeria"="Nigeria") # need to change the name
p1 <- plot_richness(Phy.f_std_birth, x="Study_site",measures=c("Shannon","Simpson", "Chao1", "Observed")
                    # title = "Alpha Diversity at birth"
)+xlab("Study site")+
  geom_boxplot(alpha=0.6)+ theme_bw()+ theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
                                             strip.text = element_text(size=17, face="bold"), axis.text.y=element_text(size=16, face="bold"),
                                             text = element_text(size=17, face="bold"))

# and add statistics
# https://rpubs.com/lconteville/713954
a_my_comparisons <- list( c("South Africa", "Nigeria"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
p1.1<- p1 + stat_compare_means(method = "wilcox.test", 
                                  comparisons = a_my_comparisons, 
                                  label = "p.signif", symnum.args = symnum.args) 
p1.1  #significantly different - CT increases the alpha diversity over the time whereas decased in Nigerian infants

p1.1

#ggsave(file = "./figure_output/figS_alpha_diversity_birth_SAvsNigeria.differnt.matrics.pdf", plot = Fig1_B, dpi = 500, width = 7, height = 7) # in this way, I can control the size of the figure
dev.off()






### FigS  (Beta-diversity)---- Use differnt matrics [Bray curtis]
### 
Phy.f_std <- readRDS("Phy.f_std.RDS")
#Phy.f_std_birth <- subset_samples(Phy.f_std, Visit=="Week 1") 
levels(sample_data(Phy.f_std)$"Study_site") <- c("CT"="South Africa","Nigeria"="Nigeria")
ord.BC <- ordinate(Phy.f_std, method = "PCoA", distance = "bray", k=3, trymax=1000) # try PcOA as NMDS does not work
color = c("Study_site")
shape = c("Visit")
#title= c("PCoA of 16S microbiome at birth (Jensen-Shannon distance; k=2)")  
MDS = plot_ordination(Phy.f_std, ord.BC, color = color, shape = shape#, title = title
)
p1  = MDS+theme_bw() +theme(axis.text=element_text(size=17, face="bold"),
                                axis.title=element_text(size=17,face="bold"), 
                                legend.title=element_text(size=17 ,face="bold"),
                                legend.text = element_text(size =16))+
  labs(color=color)+geom_point(size=3)+
  guides(color = guide_legend(title = "Study site"))+  # change the ledgend labee title
  guides(shape = guide_legend(title = "Community cluster")) # change the shape title

p1 

#ggsave(file = "./figure_output/fig1_birth_PCoA_jsd_k2.pdf", plot = Fig1_C, dpi = 500, width = 9, height = 7) # in this way, I can control the size of the figure
#dev.off()


######################################################
###### checking sequencing data----------
######################################################
Phy.f  <- readRDS("Phy.f.RDS")
Phy.f  <- readRDS("Phy.f_std.RDS")

sample_sums(Phy.f)
sample_sums(Phy.f_std)
df <- data.frame(otu_table(Phy.f_std))
View(df)

sdt1 = data.table(as(sample_data(Phy.f), "data.frame"),
                  TotalReads = sample_sums(Phy.f), keep.rownames = TRUE)
median(sdt1$TotalReads) # [1] 53719
sum(sdt1$TotalReads) # [1] 26248372
length(sdt1$TotalReads) #[1] 442
26248372/442 #[1] 59385.46


sdt2 = data.table(as(sample_data(Phy.f_std), "data.frame"),
                  TotalReads = sample_sums(Phy.f_std), keep.rownames = TRUE)
median(sdt2$TotalReads) # [1] 53871
sum(sdt2$TotalReads) #[1] 23784821
length(sdt2$TotalReads)  # 442
23784821/442   #[1] 53811.81


setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth

ggplot(data=sdt1, aes(x=Visit3, y=TotalReads))+ theme_bw()+
  geom_boxplot() +geom_jitter(width=0.3, alpha=0.2)+ 
  facet_wrap(~Study_site,scales = "free_y")+
  labs(x= "Visit") +labs(y= "Read depth (before standardization)") + 
  #ggtitle("Baby Tetanus titer iHEU vs iHUU") + 
  guides(color = guide_legend(title = "Study site"))  + #ylim(0, 7.0) + 
  theme(axis.title=element_text(size=17,face="bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
        strip.text = element_text(size=17, face="bold"), 
        axis.text.y=element_text(size=17, face="bold"),
        text = element_text(size=17, face="bold"))+ 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))


ggplot(data=sdt2, aes(x=Visit3, y=TotalReads))+ theme_bw()+
  geom_boxplot() +geom_jitter(width=0.3, alpha=0.2)+ 
  facet_wrap(~Study_site,scales = "free_y")+
  labs(x= "Visit") +labs(y= "Read depth (after standardization)") + 
  #ggtitle("Baby Tetanus titer iHEU vs iHUU") + 
  guides(color = guide_legend(title = "Study site"))  + #ylim(0, 7.0) + 
  theme(axis.title=element_text(size=17,face="bold"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=17, face="bold"), 
        strip.text = element_text(size=17, face="bold"), 
        axis.text.y=element_text(size=17, face="bold"),
        text = element_text(size=17, face="bold"))+ 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))





