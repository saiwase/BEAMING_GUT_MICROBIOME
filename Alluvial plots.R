### scripts for alluvial plots
library(phyloseq)
library(ggplot2)
library(msm)
library(ggalluvial)

### Fig2C -------
# load phyloseq object
Phy.f_std <- readRDS("Phy.f_std.RDS")
df<- sample_data(Phy.f_std)

# PIDs that have both time points
PID.with.2TP <- df$PID2[duplicated(df$PID2)]
physeq <- subset_samples(Phy.f_std, PID2 %in% PID.with.2TP) # subset the phyloseq object

df = sample_data(physeq)[, c("Visit", "Pam3_bray", "PID", "Status2", "Study_site")]
df$Pam3_bray = as.character(unlist(df$Pam3_bray))

df$Pam3_bray[df$Pam3_bray == "Cluster 1"] <- "1"
df$Pam3_bray[df$Pam3_bray == "Cluster 2"] <- "2"
df$Pam3_bray[df$Pam3_bray == "Cluster 3"] <- "3"
df$Pam3_bray = as.factor(df$Pam3_bray)

statetable.msm(Pam3_bray, PID, data = df) # transition of the cluster from one to another

df$Pam3_bray <- forcats::fct_recode(df$Pam3_bray, "Cluster 1 " = "1", "Cluster 2 " = "2", "Cluster 3 " = "3")

Fig2C <- ggplot(as.data.frame(df), 
                 aes(x = Visit, stratum = Pam3_bray, alluvium = PID, fill = Pam3_bray)) + 
  geom_alluvium(aes(fill = Pam3_bray), width = 1/16) + 
  geom_stratum(width = 1/3, alpha = .8) + 
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) + 
  scale_x_discrete(limits = c("Week 1", "Week 15"), expand = c(.05, .05, .05, .05)) + 
  scale_fill_brewer(type = "qual", palette = "Set1") + 
  theme_bw() + facet_wrap(~Study_site) + 
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) + 
  theme(axis.text.y = element_text(size = 17, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 17, face = "bold"), 
        axis.title = element_text(size = 18, face = "bold"), 
        legend.title = element_text(size = 17, face = "bold"), 
        legend.text = element_text(size = 17), 
        strip.text = element_text(size = 18, face = "bold")) + 
  guides(fill = guide_legend(title = "Community cluster")) + 
  labs(y = "Number of infants") + scale_fill_manual(values = c("#8E999F", "#F1948A", "#C5E1A5"))

Fig2C


### Fig4SD -------
Phy.f_std <- readRDS("Phy.f_std.RDS")
df<- sample_data(Phy.f_std)
PID.with.2TP <- df$PID2[duplicated(df$PID2)]
physeq <- subset_samples(Phy.f_std, PID2 %in% PID.with.2TP) # only sample with 2TP

df = sample_data(physeq)[, c("Visit", "Pam3_bray", "PID", "Status2", "Study_site")]
df$Pam3_bray = as.character(unlist(df$Pam3_bray))

df$Pam3_bray[df$Pam3_bray == "Cluster 1"] <- "1"
df$Pam3_bray[df$Pam3_bray == "Cluster 2"] <- "2"
df$Pam3_bray[df$Pam3_bray == "Cluster 3"] <- "3"
df$Pam3_bray = as.factor(df$Pam3_bray)

statetable.msm(Pam3_bray, PID, data = df)

df$Pam3_bray <- forcats::fct_recode(df$Pam3_bray, "Cluster 1 " = "1", "Cluster 2 " = "2", "Cluster 3 " = "3")

Fig4SD = ggplot(data = as.data.frame(df), 
                             aes(x = Visit, stratum = Pam3_bray, alluvium = PID)) + 
  facet_wrap(~Study_site) + 
  geom_alluvium(aes(fill = Status2), width = 1/16) + 
  scale_fill_manual(values = c("#B05A7A", "#F99417")) + 
  geom_stratum(width = 1/2.5) + theme_bw() + 
  geom_text(stat = "stratum", size = 5, aes(label = after_stat(stratum))) + 
  scale_x_discrete(limits = c("Week 1", "Week 15"), expand = c(0.15, 0.05)) + 
  guides(fill = guide_legend(title = "Exposure status")) + 
  theme(axis.text.y = element_text(size = 17, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 17, face = "bold"), 
        axis.title = element_text(size = 18, face = "bold"), 
        legend.title = element_text(size = 17, face = "bold"), 
        legend.text = element_text(size = 17), 
        strip.text = element_text(size = 18, face = "bold")) + 
  labs(y = "Number of infants")

Fig4SD

