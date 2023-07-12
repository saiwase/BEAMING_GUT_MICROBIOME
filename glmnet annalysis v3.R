### script for glmnet analysis
library(microbiome)
library(phyloseq)
library(readxl)
library(xlsx)
library(dplyr)
library(stringr)
library(glmnet)

### create a function for glmnet
glmnet.fun <- function(){
  # <Tetanus data>
  Tetanus <- read_excel("Tetanus titer data.xlsx")
  # setup the data in the right format
  Tetanus$Birth_Titer[which(Tetanus$Birth_Titer == "<0.1")] <- "0.099"
  Tetanus$W15_Titer[which(Tetanus$W15_Titer == "<0.1")] <- "0.099"
  Tetanus$Birth_Titer <- as.numeric(Tetanus$Birth_Titer)
  Tetanus$W15_Titer <- as.numeric(Tetanus$W15_Titer)
  Tetanus$Birth_Titer <- round(Tetanus$Birth_Titer, 3) 
  Tetanus$W15_Titer <- round(Tetanus$W15_Titer, 3) 
  Tetanus$Status2 <- as.factor(Tetanus$Status2)
  Tetanus$Study_site <- as.factor(Tetanus$Study_site) 
  
  Tetanus <- Tetanus %>% 
    dplyr::select(PID, Birth_Titer, W15_Titer) # extract titer info
  
  # <Meta general>
  Phy.f_std <- readRDS("Phy.f_std.RDS")
  current_meta <- tibble(data.frame(sample_data(Phy.f_std)))
  meta_general <- current_meta  %>% dplyr::select(PID, Status2, Study_site) %>% unique()
  
  # <Rank-transformed ASV >
  Phy.f_std <- readRDS("Phy.f_std.RDS")
  meta <- data.frame(sample_data(Phy.f_std))
  PID.of.interest <- meta$PID[duplicated(meta$PID)] # PIDs with two time points
  phy_country <- subset_samples(Phy.f_std, Study_site == Country  & DeliverVag == "1")
  phy_country <- do.call("subset_samples", list(quote(phy_country), substitute(PID %in% PID.of.interest)))
  
  length(sample_data(phy_country)$PID[duplicated(sample_data(phy_country)$PID)])
  
  # select the top taxa
  N <- 50
  top_taxa <- names(sort(taxa_sums(phy_country), decreasing = TRUE))[1:N]
  phy_country_only.top <- prune_taxa(top_taxa, phy_country)
  
  ## get the top taxa
  phy_country_only.top_Birth <- subset_samples(phy_country_only.top, Visit == "Week 1")
  rownames(otu_table(phy_country_only.top_Birth))
  
  topASV <- t(abundances(phy_country_only.top_Birth)[1:N, ])  
  topASV[1:10, 1:10]
  
  # apply::MARGIN = 2: the manipulation is performed on columns/ apply::MARGIN = 1: the manipulation is performed on rows
  topASV_birth <<- t(apply(topASV, 1, rank))
  topASV_birth[1:10, 1:10] 
  row.names(topASV_birth) = meta(phy_country_only.top_Birth)$PID
  
  ### do the same thing for W15
  phy_country_only.top_W15 <- subset_samples(phy_country_only.top, Visit == "Week 15")
  rownames(otu_table(phy_country_only.top_W15))
  
  topASV <- t(abundances(phy_country_only.top_W15)[1:N, ])  
  topASV[1:10, 1:10]
  
  # apply::MARGIN = 2: the manipulation is performed on columns/ apply::MARGIN = 1: the manipulation is performed on rows
  topASV_W15 <<- t(apply(topASV, 1, rank))
  topASV_W15[1:10, 1:10] 
  row.names(topASV_W15) = meta(phy_country_only.top_W15)$PID
  
  # check
  identical(colnames(topASV_birth), colnames(topASV_W15))
  identical(rownames(topASV_birth), rownames(topASV_W15))
  
  ### get taxa names
  tax_list <- read.csv("tax_list_BEAMING.csv")
  tax_list <- tibble(tax_list)
  colnames(tax_list)[1] <- "ASV.ID"
  tax.name <- tax_list %>% 
    mutate(Genus_species = paste0(Genus, "_", Species)) %>% dplyr::select(ASV.ID, Genus_species)
  
  # change the W1 dataset names
  col.name <- colnames(topASV_birth)
  col.name2 <- tibble(col.name)
  col.name3 <- left_join(col.name2, tax.name, by = c("col.name" = "ASV.ID"))
  col.name4 <- col.name3 %>% mutate(new.name = paste0(col.name, "_", Genus_species))
  col.name5 <- col.name4 %>% pull(new.name)
  col.name6 <- str_replace_all(col.name5, "ASV", "W1_ASV")
  colnames(topASV_birth) <- col.name6
  topASV_birth[1:10, 1:10] 
  
  # change the W15 dataset names
  col.name <- colnames(topASV_W15)
  col.name2 <- tibble(col.name)
  col.name3 <- left_join(col.name2, tax.name, by = c("col.name" = "ASV.ID"))
  col.name4 <- col.name3 %>% mutate(new.name = paste0(col.name, "_", Genus_species))
  col.name5 <- col.name4 %>% pull(new.name)
  col.name6 <- str_replace_all(col.name5, "ASV", "W15_ASV")
  colnames(topASV_W15) <- col.name6
  topASV_W15[1:10, 1:10]
  
  identical(sample_data(phy_country_only.top_Birth)$PID, sample_data(phy_country_only.top_W15)$PID) # TRUE
  df1 <- data.frame(row.names = meta(phy_country_only.top_Birth)$PID, 
                    PID = meta(phy_country_only.top_Birth)$PID)
  
  ### select(1): either dataset for W1 (general metadata + ASVs at W1) OR dataset of W15 (general metadata + ASVs at W15)
  #abundance_data_v2 <- cbind(df1, topASV_birth)
  abundance_data_v2 <- cbind(df1, topASV_W15)
  
  df_abundance <- tibble(abundance_data_v2)
  
  # <Combine>
  mydata1 <- left_join(meta_general, Tetanus, by = "PID") 
  mydata2 <- left_join(mydata1, df_abundance, by = "PID") 
  
  mydata2$Status2 <- as.factor(mydata2$Status2)
  mydata2$Status2 <- factor(mydata2$Status2, levels = c("iHUU", "iHEU"))
  
  mydata2$Study_site <- as.factor(mydata2$Study_site)
  
  # glmnet
  ASVs <- names(df_abundance) # lsit of ASVs
  
  # subset
  subset <- mydata2 %>% filter (Study_site == Country) %>% 
    dplyr::select(-PID, -Birth_Titer, -Study_site) %>% na.omit()  # select(2): whether take the "Birth_Titer" out or not
  
  names(subset)
  subset$W15_Titer <- as.numeric(subset$W15_Titer)
  
  # define the response and predictors
  y <- subset$W15_Titer
  x <- model.matrix(W15_Titer ~ ., data = subset)[, -1]
  
  set.seed(29)
  cvfit <- cv.glmnet(x, y, alpha = 1, family = "gaussian", type.measure = "mae", nfolds = 10) 
  cvfit$lambda.min
  cvfit$lambda.1se
  coef(cvfit)
  
  #plot(cvfit)
  best_lasso_coef.subset <<- coef(cvfit, s = c(cvfit$lambda.1se, cvfit$lambda.min))# this only works for numerical
  beta.glmnet <<- coef(cvfit, s = cvfit$lambda.min)  # for plot 
  
  print(best_lasso_coef.subset)
  
}

# run the function ------
# in the function, choose (1) whether to use top ASVs dataset from W1 or W15 (2) whether take the "Birth_Titer" out of the equation
Country <- "Nigeria"
glmnet.fun()

Country <- "South Africa"
glmnet.fun()


# save the glmnet output ------
output <- best_lasso_coef.subset 
df <- as.data.frame(as.matrix(output), check.names = T)
colnames(df) <- c("lambda.1se", "lambda.min")
df <- round(df, 5)
print(df)
#write.xlsx(df, "Glmnet.output.xlsx", col.names = TRUE, row.names = TRUE) # save the result


### Bar plot based on the glmnet results ------
# FigS5B/FigS6A/FigS6B
df <- data.frame(ASV.ID = rownames(beta.glmnet)[-1], 
                 Beta.coef = as.numeric(beta.glmnet)[-1]) %>% tibble()

df$ASV.ID <- str_replace_all(df$ASV.ID, "W1_", "")
df$ASV.ID <- str_replace_all(df$ASV.ID, "W15_", "")
df$ASV.ID <- str_replace_all(df$ASV.ID, "`", "")

df$ASV <- sapply(strsplit(df$ASV.ID, "\\_"), `[`, 1)

# add taxa names
tax_list <- read.csv("tax_list_BEAMING.csv") # tax table 
tax_list <- tibble(tax_list)
colnames(tax_list)[1] <- "ASV"
tax_list <- tax_list %>% 
  mutate(Genus_species = paste0(Genus, "_", Species)) %>%
  relocate(ASV, Genus_species)

df2 <- left_join(df, tax_list, by = "ASV") # add taxa info

# change the description for plotting
df3.2 <- df2 %>% 
  mutate(Covariate = ifelse(is.na(df2$Kingdom), ASV.ID, paste0(Genus, " ", Species, " [", ASV, "]"))) %>% 
  mutate(Covariate2 = ifelse(is.na(df2$Kingdom), ASV.ID, paste0(Genus, " ", Species)))

df3.2$Covariate2 <- str_replace(df3.2$Covariate2, "NA", "") 
df3.2$Covariate2 <- str_replace(df3.2$Covariate2, "Status2iHEU", "Exposure status (iHEU)")

p1 = df3.2%>% filter(Beta.coef > 0 |Beta.coef < 0) %>% 
  mutate(color = ifelse(Beta.coef < 0, "#blue", "#red")) %>% 
  arrange(desc(abs(Beta.coef)))%>%
  ggplot(., aes(x = reorder(Covariate2,  + Beta.coef), y = Beta.coef, fill = Family)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw() + coord_flip() + 
  ylab("Coefficients") + 
  xlab("Selected variables") + 
  theme(axis.text.x = element_text(size = 17, face = "bold"), 
        axis.text.y = element_text(size = 17, face = "bold"), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 17, face = "bold"), 
        axis.title = element_text(size = 17, colour = "black", face = "bold")) + 
  geom_hline(yintercept = 0)

p1
