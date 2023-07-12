### scripts for alpha- and beta-diversity plots
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(vegan)

### Fig1B (Alpha diversity) -------
# load phyloseq object
Phy.f_std <- readRDS("Phy.f_std.RDS")
# filter only W1 samples
Phy.f_std_W1 <- subset_samples(Phy.f_std, Visit == "Week 1")

p1 = plot_richness(Phy.f_std_W1, x = "Study_site", measures = "Shannon", color = "Study_site") + 
  geom_boxplot(alpha = 0.6) + theme_bw() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 17, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"), 
        strip.text = element_blank(), 
        text = element_text(size = 17, face = "bold")) + xlab("Study site") + labs(y = "Shannon Index")

a_my_comparisons <- list( c("South Africa", "Nigeria"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
Fig1_B = p1 + stat_compare_means(method = "wilcox.test", paired = F, comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args) + 　
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

Fig1_B  


### alpha diversity - adjusted for sequencing batch or demographic factors that significantly differed between countries
Df <- data.frame(sample_data(Phy.f_std_W1))
names(Df)

Df2 <- Df %>% select(Shannon, Study_site, SequenceLot, Education, Occupation, Houses, Refrigerator, 
                     Running_water, Marital_status, Num_Preg, Mum_weight, Mum_Age, Ges_age, DeliverVag) %>% tibble()%>% na.omit()

ancova_model1 <- aov(Shannon ~ Study_site + SequenceLot, data = Df2)
summary(ancova_model1)
summary.lm(ancova_model1)

ancova_model2 <- aov(Shannon ~ Study_site + Mum_weight + Mum_Age + Num_Preg, data = Df2)
summary(ancova_model2)
summary.lm(ancova_model2) 

ancova_model3 <- aov(Shannon ~ Study_site + Mum_weight + Mum_Age + Num_Preg + SequenceLot, data = Df2)
summary(ancova_model3)
summary.lm(ancova_model3)   

ancova_model4 <- aov(Shannon ~ Study_site + SequenceLot + Education + Occupation + Houses + Refrigerator + Running_water + Marital_status + Num_Preg + Mum_weight + Ges_age + DeliverVag, data = Df2)
summary(ancova_model4)
summary.lm(ancova_model4)  


### Fig1C (Beta diversity) -------
Phy.f_std <- readRDS("Phy.f_std.RDS")
Phy.f_std_W1 <- subset_samples(Phy.f_std, Visit == "Week 1") 
ord.BC <- ordinate(Phy.f_std_W1, method = "PCoA", distance = "bray", k = 3, trymax = 1000)
color = c("Study_site")
shape = c("Pam3_bray")

p1 = plot_ordination(Phy.f_std_W1, ord.BC, color = color, shape = shape)

Fig1C = p1 + geom_point(size = 3) + theme_bw() + 
  theme(axis.text = element_text(size = 17, face = "bold"), 
        axis.title = element_text(size = 17, face = "bold"), 
        legend.title = element_text(size = 17 , face = "bold"), 
        legend.text = element_text(size = 16)) + 
  labs(color = color) + 
  guides(color = guide_legend(title = "Study site")) + 
  guides(shape = guide_legend(title = "Community cluster")) 

Fig1C 


### PERMANOVA
Phy.f_std <- readRDS("Phy.f_std.RDS")
Phy.f_std_W1 <- subset_samples(Phy.f_std, Visit == "Week 1") 
diss <- phyloseq::distance(Phy.f_std_W1, "bray", parallel = TRUE) 
Phy.f_std_W1_df <- data.frame(row.names = sample_names(Phy.f_std_W1), sample_data(Phy.f_std_W1))

# adjusted for sequencing batch or demographic factors that significantly differed between countries
sample_data(Phy.f_std_W1)$Education <- as.factor(sample_data(Phy.f_std_W1)$Education)
sample_data(Phy.f_std_W1)$Occupation <-  as.factor(sample_data(Phy.f_std_W1)$Occupation)
sample_data(Phy.f_std_W1)$Houses <- as.factor(sample_data(Phy.f_std_W1)$Houses)
sample_data(Phy.f_std_W1)$Refrigerator <- as.factor(sample_data(Phy.f_std_W1)$Refrigerator)
sample_data(Phy.f_std_W1)$Running_water <- as.factor(sample_data(Phy.f_std_W1)$Running_water)
sample_data(Phy.f_std_W1)$Marital_status <- as.factor(sample_data(Phy.f_std_W1)$Marital_status)

set.seed(2)
adonis2(diss ~ Study_site, data = Phy.f_std_W1_df, permutations = 999, by = "terms", na.action = na.exclude)
adonis2(diss ~ Study_site + SequenceLot, data = Phy.f_std_W1_df, permutations = 999, by = "terms", na.action = na.exclude)
adonis2(diss ~ Study_site + SequenceLot + Education + Occupation + Houses + Refrigerator + Running_water + Marital_status + Num_Preg + Mum_weight + Mum_Age + Ges_age + DeliverVag, data = Phy.f_std_W1_df, permutations = 999, by = "terms", na.action = na.exclude)


### Fig2A (Alpha diversity) -------
Phy.f_std <- readRDS("Phy.f_std.RDS")
p1 <- plot_richness(Phy.f_std, x = "Visit", measures = "Shannon", color = "Study_site") + 
  geom_boxplot(alpha = 0.6) + 
  geom_line(aes(group = PID2), size = 0.1) + # join by PID
  facet_wrap(~Study_site) + theme_bw() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 17, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"), 
        strip.text = element_text(size = 17, face = "bold"), 
        text = element_text(size = 17, face = "bold")) + labs(y = "Shannon Index")

a_my_comparisons <- list( c("Week 1", "Week 15"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
Fig2A = p1 + stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, 
                                label = "p.signif", symnum.args = symnum.args) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

Fig2A  


### Fig2B (Beta diversity) -------
Phy.f_std <- readRDS("Phy.f_std.RDS")
set.seed(20) 
ord.BC <- ordinate(Phy.f_std, method = "PCoA", distance = "bray", k = 2, trymax = 1000)
color = c("Study_site")
shape = c("Visit")

p1 = plot_ordination(Phy.f_std, ord.BC, color = color, shape = shape)
Fig2B = p1 + theme_bw() + 
  theme(axis.text = element_text(size = 17, face = "bold"), 
        axis.title = element_text(size = 17, face = "bold"), 
        legend.title = element_text(size = 17 , face = "bold"), 
        legend.text = element_text(size = 16)) + 
  labs(color = color, shape = shape) + 
  geom_point(size = 2) + 
  guides(color = guide_legend(title = "Study site")) 

Fig2B


### PERMANOVA
diss <- phyloseq::distance(Phy.f_std, "bray", parallel = TRUE)
Phy.f_std_df <- data.frame(row.names = sample_names(Phy.f_std), sample_data(Phy.f_std)) 
set.seed(2)
adonis2(diss ~ Study_site, data = Phy.f_std_df, permutations = 999, by = "terms", na.action = na.exclude)  
adonis2(diss ~ Visit, data = Phy.f_std_df, permutations = 999, by = "terms", na.action = na.exclude)  


### FigS1 (Alpha diversity) -------
Phy.f_std <- readRDS("Phy.f_std.RDS")
data <- tibble(data.frame(sample_data(Phy.f_std))) %>%  
  mutate(Visit5 = case_when(age_days <= 3 ~ "Meconium", 
                            age_days >= 4 & age_days <= 20 ~ "W1 stool", 
                            Visit == "Week 15" ~ "W15 stool")) %>% select(PID, age_days, Visit, Visit5) 

sample_data(Phy.f_std)$Visit5 <- data$Visit5
sample_data(Phy.f_std)$Visit5 <- as.factor(sample_data(Phy.f_std)$Visit5)

# filter only Meconium samples
Phy.f_std_W1 <- subset_samples(Phy.f_std, Visit5 == "Meconium")
p1 = plot_richness(Phy.f_std_W1, x = "Study_site", measures = "Shannon", color = "Study_site") + 
  geom_boxplot(alpha = 0.6) + theme_bw() + xlab("Study site") + 
  theme(legend.position = "none", 
        strip.text = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 17, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"), 
        text = element_text(size = 17, face = "bold")) + 
  labs(y = "Shannon Index")

a_my_comparisons <- list( c("South Africa", "Nigeria"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
FigS1 = p1 + stat_compare_means(method = "wilcox.test", hide.ns = FALSE, 
                                comparisons = a_my_comparisons, 
                                label = "p.signif", symnum.args = symnum.args) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) 

FigS1


### FigS2A (Alpha diversity) -------
Phy.f_std <- readRDS("Phy.f_std.RDS")
Phy.f_std_EBF_W15 <- subset_samples(Phy.f_std, Feeding_W15_strictver == "Exclusive_breastfeeding" & Visit == "Week 15")

p1 = plot_richness(Phy.f_std_EBF_W15, x = "Study_site", measures = "Shannon", color = "Study_site") + 
  xlab("Study site") + 
  geom_boxplot(alpha = 0.6) + 
  theme_bw() + theme(legend.position = "none", 
                    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 17, face = "bold"), 
                    axis.text.y = element_text(size = 16, face = "bold"), 
                    strip.text = element_blank(),  
                    text = element_text(size = 17, face = "bold")) + 
  labs(y = "Shannon Index")

a_my_comparisons <- list(c("South Africa", "Nigeria"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
FigS2A = p1 + stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

FigS2A


### FigS2B (Beta diversity) -------
Phy.f_std <- readRDS("Phy.f_std.RDS")
Phy.f_std_EBF_W15 <- subset_samples(Phy.f_std, Feeding_W15_strictver == "Exclusive_breastfeeding" & Visit == "Week 15")

ord.BC = ordinate(Phy.f_std_EBF_W15, method = "PCoA", distance = "bray", k = 3, trymax = 1000)
color = c("Study_site")

p1 = plot_ordination(Phy.f_std_EBF_W15, ord.BC, color = color)

FigS2B = p1 + theme_bw() + theme(axis.text = element_text(size = 15, face = "bold"), 
                              axis.title = element_text(size = 15, face = "bold"), 
                              strip.text = element_text(size = 15, face = "bold"), 
                              legend.title = element_text(size = 15, face = "bold"), 
                              legend.text = element_text(size = 15), 
                              legend.position = "bottom") + 
  labs(color = color) + geom_point(size = 2) + 
  guides(color = guide_legend(title = "Study site")) 

FigS2B

### PERMANOVA 
diss <- phyloseq::distance(Phy.f_std_EBF_W15, "bray", parallel = TRUE)  
Phy.f_std_EBF_W15_df <- data.frame(row.names = sample_names(Phy.f_std_EBF_W15), sample_data(Phy.f_std_EBF_W15))
set.seed(2)
adonis2(diss ~ Study_site, data = Phy.f_std_EBF_W15_df, permutations = 999, by = "terms", na.omit)


### FigS3A (Alpha diversity) -------
Phy.f_std <- readRDS("Phy.f_std.RDS")
sample_data(Phy.f_std)$DeliverVag <- as.factor(sample_data(Phy.f_std)$DeliverVag)
levels(sample_data(Phy.f_std)$DeliverVag) <- c("0" = "C-section", "1" = "Vaginal birth")

# filter only vaginal birth 
Phy.f_std_VagiBirth <- subset_samples(Phy.f_std, DeliverVag == "Vaginal birth")

p1 = plot_richness(Phy.f_std_VagiBirth, x = "Study_site", measures = "Shannon", color = "Study_site") + 
  xlab("Study site") + 
  geom_boxplot(alpha = 0.6) + 
  theme_bw() + facet_wrap(~ Visit) + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 17, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"), 
        text = element_text(size = 17, face = "bold")) + 
  labs(y = "Shannon Index")

symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
FigS3A <- p1 + stat_compare_means(method = "wilcox.test", comparisons = list(c("South Africa", "Nigeria")), 
                                label = "p.signif", symnum.args = symnum.args) + scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

FigS3A


### FigS3B (Beta diversity) -------
Phy.f_std <- readRDS("Phy.f_std.RDS")
levels(sample_data(Phy.f_std)$"DeliverVag") <- c("0" = "C-section", "1" = "Vaginal birth")

# filter only vaginal birth 
Phy.f_std_VagiBirth <- subset_samples(Phy.f_std, DeliverVag == "Vaginal birth")
ord.BC <- ordinate(Phy.f_std_VagiBirth, method = "PCoA", distance = "bray", k = 3, trymax = 1000) 
color = c("Study_site")

p1 = plot_ordination(Phy.f_std_VagiBirth, ord.BC, color = color)
FigS3B = p1 + theme_bw() + theme(axis.text = element_text(size = 13, face = "bold"), 
                              axis.title = element_text(size = 16, face = "bold"), 
                              strip.text = element_text(size = 15, face = "bold"), 
                              legend.title = element_text(size = 14), 
                              legend.text = element_text(size = 12)) + facet_wrap(~Visit) + 
  labs(color = color) + geom_point(size = 2) + 
  guides(color = guide_legend(title = "Study site"))

FigS3B


### PERMANOVA 
Phy.f_std_VagiBirth <- subset_samples(Phy.f_std, DeliverVag == "Vaginal birth" & Visit = "Week 1")
diss <- phyloseq::distance(Phy.f_std_VagiBirth, "bray", parallel = TRUE)  
Phy.f_std_VagiBirth_df <- data.frame(row.names = sample_names(Phy.f_std_VagiBirth), sample_data(Phy.f_std_VagiBirth)) 
set.seed(2)
adonis2(diss ~ Study_site, data = Phy.f_std_VagiBirth_df, permutations = 999, by = "terms", na.omit)


### FigS4A (Alpha diversity) -------
Phy.f_std <- readRDS("Phy.f_std.RDS")
p1 = plot_richness(Phy.f_std, x = "Status2", measures = "Shannon", color = "Status2") + 
  geom_boxplot(alpha = 0.6, aes(fill = sample_data(Phy.f_std)$Status2)) + 
  facet_wrap(Study_site~Visit) + theme_bw() + 
  xlab("HIV exposure status") + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 17, face = "bold"), 
        axis.text.y = element_text(size = 16, face = "bold"), 
        strip.text = element_text(size = 17, face = "bold"), 
        text = element_text(size = 17, face = "bold")) + 
  scale_color_manual(values = c("#B05A7A", "#F99417")) + 
  scale_fill_manual(values = c("#B05A7A", "#F99417")) + 
  labs(y = "Shannon Index")

a_my_comparisons <- list( c("iHEU", "iHUU"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
FigS4A = p1 + stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, 
                                  label = "p.signif", symnum.args = symnum.args, size = 5) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) 

FigS4A


### FigS4B (Beta diversity) -------
Phy.f_std <- readRDS("Phy.f_std.RDS")
set.seed(20) 
ord.BC <- ordinate(Phy.f_std, method = "PCoA", distance = "bray", k=2, trymax=1000)

color = c("Status2")
shape = c("Study_site")
p1 = plot_ordination(Phy.f_std, ord.BC, color = color, shape=shape)
FigS4B = p1 + theme_bw() + theme(axis.text=element_text(size=17, face="bold"),
                                axis.title=element_text(size=17,face="bold"), 
                                legend.title=element_text(size=17 ,face="bold"),
                                legend.text = element_text(size =16)) +  
  labs(color=color, shape=shape) + geom_point(size=2) +  
  guides(color = guide_legend(title = "Exposure status"), shape = guide_legend(title = "Study site")) + 
  scale_color_manual(values=c("#B05A7A","#F99417"))

FigS4B


### PERMANOVA 
diss <- phyloseq::distance(Phy.f_std, "bray", parallel=TRUE) 
Phy.f_std_df <- data.frame(row.names=sample_names(Phy.f_std),sample_data(Phy.f_std))
set.seed(20)
adonis2(diss ~ Status2, data=Phy.f_std_df, permutations = 999, by = "terms", na.omit)  
