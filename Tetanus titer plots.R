library(readxl)
library(phyloseq)
library(dplyr)
library(ggpubr)

### Prepare the dataset
# load Tetanus titer data (including mum's titer)
Tetanus <- read_excel("Tetanus titer_long.xlsx")

Tetanus$Baby_titer[which(Tetanus$Baby_titer == "<0.1")] <- "0.099"  
Tetanus$Baby_titer <- as.numeric(Tetanus$Baby_titer)
Tetanus$Baby_titer <- round(Tetanus$Baby_titer, 3)
Tetanus$Status2 <- as.factor(Tetanus$Status2)

# extract PIDs that have microbiome data available
Phy.f_std <- readRDS("Phy.f_std.RDS")
stool_PID_info <- data.frame(sample_data(Phy.f_std))$PID

# subset the Tetanus data based on the PID list 
Tetanus2 <- Tetanus %>% filter(PID %in% stool_PID_info) 

# subset only birth data
Tetanus_birth <- Tetanus2 %>% filter(Visit == "Birth")


### Fig4A (Scatter plots) -------
Fig4A = ggscatter(data = Tetanus_birth, x = "Baby_titer", y = "Mum_titer_birth", 
                       add = "reg.line", conf.int = TRUE, size = 2, color = "Status2") + 
  theme(axis.text = element_text(size = 17, face = "bold"), 
        axis.title = element_text(size = 17, face = "bold"), 
        strip.text = element_text(size = 17, face = "bold"), 
        legend.position = "none") + 
  scale_color_manual(values = c("#B05A7A", "#F99417")) + 
  scale_fill_manual(values = c("#B05A7A", "#F99417")) + 
  labs(x = "Anti-TT IgG at W1 (IU/ml) [infant]") + 
  labs(y = "Anti-TT IgG (IU/ml) [mother]")
Fig4A


### Fig4B (box plots) -------
Phy.f_std <- readRDS("Phy.f_std.RDS")
meta <- data.frame(sample_data(Phy.f_std)) %>% tibble()  %>% filter(Baby_titer >0)
meta$Baby_titer # Tetanus titer of infants

p1 = ggplot(data = meta, aes(x = Status2, y = Baby_titer)) + 
  theme_bw() + 
  geom_boxplot(aes(color = Status2, fill = Status2), alpha = 0.5) + 
  geom_jitter(aes(color = Status2, fill = Status2, alpha = 0.5), width = 0.3, alpha = 0.2) + 
  facet_wrap(~Visit, scales = "free_y") + 
  theme(axis.title = element_text(size = 17, face = "bold"), 
        strip.text = element_text(size = 17, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 17, face = "bold"), 
        axis.text.y = element_text(size = 17, face = "bold"), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 16)) + 
  labs(x = "HIV exposure status") + 
  labs(y = "Anti-TT IgG (IU/ml)") + 
  labs(color = "HIV exposure status") + 
  labs(fill = "HIV exposure status") + 
  ylim(0, 6.7) + 
  scale_color_manual(values = c("#B05A7A", "#F99417")) + 
  scale_fill_manual(values = c("#B05A7A", "#F99417"))

# adjust for multiple comparison
anno_df <- compare_means(Baby_titer ~ Status2, meta, method = "wilcox.test", paired = FALSE, 
                         group.by = "Visit", ref.group = NULL, 
                         p.adjust.method = "BH")
anno_df$y.position <- c(6, 6) 
anno_df$xmin <- c(1, 1)
anno_df$xmax <- c(2, 2)

Fig4B = p1 + stat_pvalue_manual(anno_df, label = "p.signif", 
                                 tip.length = 0.01, bracket.nudge.y = 0.5, label.size = 5)
Fig4B


# median value
median(meta$Baby_titer[meta$Visit == "Week 1" & meta$Status2 == "iHEU"]) # 1.21 
median(meta$Baby_titer[meta$Visit == "Week 1" & meta$Status2 == "iHUU"]) # 2.2 

median(meta$Baby_titer[meta$Visit == "Week 15" & meta$Status2 == "iHEU"]) # 1.3 
median(meta$Baby_titer[meta$Visit == "Week 15" & meta$Status2 == "iHUU"]) # 2.5

median(meta$Baby_titer[meta$Visit == "Week 1" & meta$Study_site == "South Africa"]) # 1 
median(meta$Baby_titer[meta$Visit == "Week 1" & meta$Study_site == "Nigeria"]) #1.5

median(meta$Baby_titer[meta$Visit == "Week 15" & meta$Study_site == "South Africa"]) # 1.9 
median(meta$Baby_titer[meta$Visit == "Week 15" & meta$Study_site == "Nigeria"]) # 1.6


### Fig5A (correlation plot) -------
Phy.f_std <- readRDS("Phy.f_std.RDS")
meta <- data.frame(sample_data(Phy.f_std))

W15_titer <- meta %>% filter(Visit == "Week 15") %>% # select Tetatnus titer at W15
  select(PID, Baby_titer) %>% dplyr::rename("W15_titer" = "Baby_titer") %>% filter(W15_titer > 0)

meta2 <- left_join(meta, W15_titer, by = "PID") # update the metadata
 
Fig5A = ggscatter(data = meta2, x = "Shannon", y = "W15_titer", color = "Study_site", 
                       add = "reg.line", conf.int = TRUE, size = 2) + 
  stat_cor(method = "spearman", size = 5, label.x = 1.4, label.y = 4.7) + 
  facet_wrap(Visit~Study_site, nrow = 2, scales = "free_y") + 
  theme_bw() + 
  theme(axis.text = element_text(size = 17, face = "bold"), 
        axis.title = element_text(size = 17, face = "bold"), 
        legend.title = element_text(size = 17, face = "bold"), 
        legend.text = element_text(size = 17), 
        strip.text = element_text(size = 17, face = "bold")) + 
  ylim(0, 5) + xlim(0, 4.1) + 
  labs(fill = "Study site") + 
  labs(color = "Study site") + 
  labs(x = "Shannon diversity") + 
  labs(y = "Anti-TT IgG at W15 (IU/ml)")
Fig5A


### FigS5A (box plot) -------
meta <- data.frame(sample_data(Phy.f_std)) %>% tibble()  %>% filter(Baby_titer >0) # filter only samples that Tetanus titer is available

# boxplot
bxp = ggboxplot(meta, x = "Study_site", y = "Baby_titer", color = "Visit", fill = "Visit", alpha = 0.7, palette = c("#4FC3F7", "#48648F"))

# stat1
TestFrag1 <- compare_means(Baby_titer ~ Visit, meta, method = "wilcox.test", paired = FALSE, 
                           group.by = "Study_site", ref.group = NULL, p.adjust.method = "BH")

# add position info
TestFrag1$y.position <- c(6.34, 6.44) 
TestFrag1$x <- c(1, 2) 
TestFrag1$xmin <- c(0.8, 1.8)
TestFrag1$xmax <- c(1.2, 2.2)

# add stat1
bxp.complex = bxp + stat_pvalue_manual(TestFrag1, label = "p.adj", 
                                        tip.length = 0.01, label.size = 5)

# stat2
TestFrag2 <- compare_means(Baby_titer ~ Study_site, meta, method = "wilcox.test", paired = FALSE, 
                           group.by = "Visit", ref.group = NULL, p.adjust.method = "BH")

# add position info
TestFrag2$y.position <- c(7.2, 6.5) 
TestFrag2$xmin <- c(1.2, 0.8)
TestFrag2$xmax <- c(2.2, 1.8)

# add stat2
bxp.complex2 = bxp.complex + stat_pvalue_manual(TestFrag2, label = "p.adj", 
                                                 tip.length = 0.01, bracket.nudge.y = 0.5, label.size = 5) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Display the plot
FigS5A = bxp.complex2 + theme_bw() + 
  theme(axis.title = element_text(size = 17, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 17, face = "bold"), 
        axis.text.y = element_text(size = 17, face = "bold"), 
        strip.text = element_text(size = 17, face = "bold"), 
        legend.title = element_text(size = 17, face = "bold"), 
        legend.text = element_text(size = 17)) + 
  labs(x = "Study site") + labs(y = "Anti-TT IgG (IU/ml)") 
FigS5A


### FigS5B (scatter plot) -------
Tetanus <- read_excel("Tetanus titer_long.xlsx")

Tetanus$Baby_titer[which(Tetanus$Baby_titer == "<0.1")] <- "0.099"  
Tetanus$Baby_titer <- as.numeric(Tetanus$Baby_titer)  # change it to numeric
Tetanus$Baby_titer <- round(Tetanus$Baby_titer, 3)  # change the number of digit
Tetanus$Status2 <- as.factor(Tetanus$Status2)

# extract PIDs that have microbiome data available
Phy.f_std <- readRDS("Phy.f_std.RDS")
stool_PID_info <- data.frame(sample_data(Phy.f_std))$PID

# subset the Tetanus data based on the PID list 
Tetanus2 <- Tetanus %>% filter(PID %in% stool_PID_info) 

# subset only W15 data in Nigerian cohort
Tetanus_W15 <- Tetanus2 %>%
  filter(Visit == "W15" & Study_site == "Nigeria" & Mum_titer_birth > 0 & Baby_titer >0)
levels(Tetanus_W15$Status2) <- c("HEU" = "iHEU", "HU" = "iHUU")

FigS5B = ggscatter(data = Tetanus_W15, x = "Baby_titer", y = "Mum_titer_birth", 
                        add = "reg.line", conf.int = TRUE, size = 2, color = "Status2") + 
  stat_cor(method = "spearman", label.x = 2.5, label.y = 5.3, size = 5) + 
  facet_wrap(~Status2) + 
  theme(axis.text = element_text(size = 17, face = "bold"), 
        axis.title = element_text(size = 17, face = "bold"), 
        strip.text = element_text(size = 17, face = "bold"), 
        legend.position = "none") + 
  labs(color = "HIV exposure status") + 
  labs(fill = "HIV exposure status") + 
scale_color_manual(values = c("#B05A7A", "#F99417")) + 
  scale_fill_manual(values = c("#B05A7A", "#F99417")) + 
  labs(x = "Anti-TT IgG (IU/ml) [infant]") + 
  labs(y = "Anti-TT IgG (IU/ml) [mother]")
FigS5B


# stat
Tetanus_W15_iHEU <- Tetanus_W15 %>% filter(Status2 == "iHEU")
res1 <- cor.test(Tetanus_W15_iHEU$Mum_titer_birth, Tetanus_W15_iHEU$Baby_titer, method = "spearman")
p_value1 <- res1$p.value  
p.adjust(p_value1, method = "BH")

Tetanus_W15_iHUU <- Tetanus_W15 %>% filter(Status2 == "iHUU")
res2 <- cor.test(Tetanus_W15_iHUU$Mum_titer_birth, Tetanus_W15v2_iHUU$Baby_titer, method = "spearman")
p_value2 <- res2$p.value
p.adjust(p_value2, method = "BH")


### FigS5C (box plot) -------
Tetanus <- read_excel("Tetanus titer_long.xlsx")

Tetanus$Baby_titer[which(Tetanus$Baby_titer == "<0.1")] <- "0.099"  
Tetanus$Baby_titer <- as.numeric(Tetanus$Baby_titer)  # change it to numeric
Tetanus$Baby_titer <- round(Tetanus$Baby_titer, 3)  # change the number of digit
Tetanus$Status2 <- as.factor(Tetanus$Status2)

# extract PIDs that have microbiome data available
Phy.f_std <- readRDS("Phy.f_std.RDS")
stool_PID_info <- data.frame(sample_data(Phy.f_std))$PID

# subset the Tetanus data based on the PID list 
Tetanus2 <- Tetanus %>% filter(PID %in% stool_PID_info) 

# extract only mum's tetanus titer (measured at birth)
Tetanus_mum <-  Tetanus2 %>% filter(Visit == "Birth" & Mum_titer_birth >0)
levels(Tetanus_mum$Status2) <- c("HEU" = "w/ HIV", "HU" = "w/o HIV") # change the annotation

p1 = ggplot(data = Tetanus_mum, aes(x = Status2, y = Mum_titer_birth, color = Status2, fill = Status2)) + 
  geom_boxplot(alpha = 0.5) + theme_bw() + geom_jitter(width = 0.3, alpha = 0.2) + 
  theme(axis.text.y = element_text(size = 17, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 17, face = "bold"), 
        axis.title = element_text(size = 17, face = "bold"), 
        strip.text = element_text(size = 17, face = "bold"), 
        legend.position = "none") + 
  scale_color_manual(values = c("#B05A7A", "#F99417")) + 
  scale_fill_manual(values = c("#B05A7A", "#F99417")) + 
  labs(x = "Mother's HIV status") + 
  labs(y = "Maternal anti-TT IgG (IU/ml)")

compare_means(Mum_titer_birth ~ Status2, Tetanus_mum, method = "wilcox.test", paired = FALSE, ref.group = NULL)
compare_means(Mum_titer_birth ~ Status2, Tetanus_mum, method = "wilcox.test", paired = FALSE, ref.group = NULL, group.by = NULL, p.adjust.method = "BH") 

# stats
a_my_comparisons <- list( c("w/ HIV", "w/o HIV"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
FigS5C = p1 + stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args) 
FigS5C


### FigS5D (box plot) -------
Phy.f_std <- readRDS("Phy.f_std.RDS")
meta <- tibble(data.frame(sample_data(Phy.f_std)))

p1 = ggplot(data = meta, aes(x = Status2, y = Baby_titer, color = Status2)) + 
  theme_bw() + 
  geom_boxplot(alpha = 0.6, aes(fill = Status2)) + geom_jitter(width = 0.3, alpha = 0.2) + 
  facet_wrap(Visit~Study_site, scales = "free_y") + 
  labs(x = "HIV exposure status") + labs(y = "Anti-TT IgG (IU/ml)") + 
  theme(strip.text = element_text(size = 17, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 17, face = "bold"), 
        axis.text.y = element_text(size = 17, face = "bold"), 
        axis.title = element_text(size = 17, face = "bold"), 
        legend.title = element_blank()) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  scale_color_manual(values = c("#B05A7A", "#F99417")) + 
  scale_fill_manual(values = c("#B05A7A", "#F99417"))

# stats - adjust for multiple comparison
anno_df <- compare_means(Baby_titer ~ Status2, meta, method = "wilcox.test", paired = FALSE, 
                         group.by = c("Visit", "Study_site"), ref.group = NULL, p.adjust.method = "BH")

anno_df$y.position <- c(6, 6, 6, 6) 
anno_df$xmin <- c(1, 1, 1, 1)
anno_df$xmax <- c(2, 2, 2, 2)

# Use adjdusted p-value and add on the plot
alpha_levels <- c(0.05, 0.01, 0.001)
asterisks <- c("*", "**", "***")
asterisks <- c("***", "**", "*")
anno_df$p.adj.signif <- ifelse(anno_df$p.adj <= alpha_levels[1], asterisks[1], 
                               ifelse(anno_df$p.adj <= alpha_levels[2], asterisks[2], 
                                      ifelse(anno_df$p.adj <= alpha_levels[3], asterisks[3], "ns")))

FigS5D = p1 + stat_pvalue_manual(anno_df, label = "p.adj.signif", tip.length = 0.01, bracket.nudge.y = 0.5, label.size = 4)
FigS5D
