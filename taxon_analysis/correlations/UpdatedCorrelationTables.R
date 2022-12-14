# Additional analyses for Wyatt Paper
# By Cliff Bueno de Mesquita
# 1. Make Table S8 (and input table for Figure 3a, Figure S9) for Wyatt Paper
# Correlations for each taxonomic level
# Important: set relative = FALSE in mctoolsr summarize code
# Counts NOT log2 transformed
# 2. Remake Figure S5e-f with new guilds

#### Setup ####
suppressMessages(library(dplyr))
suppressMessages(library(rlang))
suppressMessages(library(tibble))
library(mctoolsr)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
`%notin%` <- Negate(`%in%`)
save_pheatmap_pdf <- function(x, filename, width = 2.5, height = 2.5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

setwd("~/Documents/GitHub/SF_microbe_methane/taxon_analysis/")

# Data Imports From Wyatt
# Import OTU Table
SilvaOTUs <- read.table("../silvaOTUs/Silva_OTU_VSTcpm.txt", sep='\t')
row.names(SilvaOTUs) <- SilvaOTUs[,"OTU"]                                
SilvaOTUs <- SilvaOTUs[,-1]
otu_V <- SilvaOTUs

# Sort OTU table                                                                      
otu_V <- otu_V[order(otu_V$Consensus.lineage),]                                                        

# Make new top level plotting var (should be in PRE-PROCESS ? )
otu_V$Taxonomy <- ifelse(otu_V$Phylum == "Proteobacteria", paste(otu_V$Class), paste(otu_V$Phylum))
otu_V <- data.frame(otu_V)

# METADATA
metaDB <-read.table("../data/meta/SF_sal_meta_FIX3.txt", sep="\t", header=TRUE)

# Prune metadata to only iTag samples
# Get Sample names in OTU table               
OTU_samps <- data.frame('Sample'= colnames(otu_V))

# Merge site order and Samples
Meta_iTag <- merge(metaDB, OTU_samps, by='Sample')
rownames(Meta_iTag) <- Meta_iTag$Sample

# add any extra variables as needed
CH4 <- Meta_iTag$CH4_ug_m2_h
Meta_iTag$CH4_logn1 <- log10(CH4 - (1.05*min(CH4)))

# New Code by Cliff
dim(otu_V)
dim(Meta_iTag)
sum(colnames(otu_V) %in% Meta_iTag$Sample)

# make mctoolsr object
data_loaded <- otu_V %>%
  select(1:168) %>%
  t() %>%
  as.data.frame() %>%
  arrange(rownames(.)) %>%
  t() %>%
  as.data.frame()
map_loaded <- Meta_iTag
taxonomy_loaded <- otu_V %>%
  select(170:175) %>%
  set_names(c("taxonomy1", "taxonomy2", "taxonomy3", "taxonomy4", "taxonomy5", "taxonomy6"))
sum(colnames(data_loaded) != rownames(map_loaded))
sum(rownames(data_loaded) != rownames(taxonomy_loaded))
input <- list()
input$data_loaded <- data_loaded
input$map_loaded <- map_loaded
input$taxonomy_loaded <- taxonomy_loaded
input <- filter_taxa_from_input(input,
                                taxa_to_remove = "PH",
                                at_spec_level = 3)



#### Correlations ####
# Run Pearson correlation at each taxonomic level

# Domain
dom <- summarize_taxonomy(input, level = 1, relative = F, report_higher_tax = FALSE)
dom_t <- as.data.frame(t(dom))
cor.df.dom <- as.data.frame(matrix(NA, nrow = ncol(dom_t), ncol = 2)) %>%
  set_names(c("Domain", "Domain_r"))
for (i in 1:ncol(dom_t)) {
  m <- cor.test(dom_t[,i], input$map_loaded$CH4_logn1, method = "pearson")
  cor.df.dom$Domain[i] <- names(dom_t)[i]
  cor.df.dom$Domain_r[i] <- round(m$estimate, digits = 2)
}
cor.test(dom_t$Archaea, input$map_loaded$CH4_logn1, method = "pearson")
cor.test(dom_t$Bacteria, input$map_loaded$CH4_logn1, method = "pearson")

# Phylum
phy <- summarize_taxonomy(input, level = 2, relative = F, report_higher_tax = FALSE)
phy_t <- as.data.frame(t(phy))
cor.df.phy <- as.data.frame(matrix(NA, nrow = ncol(phy_t), ncol = 2)) %>%
  set_names(c("Phylum", "Phylum_r"))
for (i in 1:ncol(phy_t)) {
  m <- cor.test(phy_t[,i], input$map_loaded$CH4_logn1, method = "pearson")
  cor.df.phy$Phylum[i] <- names(phy_t)[i]
  cor.df.phy$Phylum_r[i] <- round(m$estimate, digits = 2)
}

# Class
cla <- summarize_taxonomy(input, level = 3, relative = F, report_higher_tax = FALSE)
cla_t <- as.data.frame(t(cla))
cor.df.cla <- as.data.frame(matrix(NA, nrow = ncol(cla_t), ncol = 2)) %>%
  set_names(c("Class", "Class_r"))
for (i in 1:ncol(cla_t)) {
  m <- cor.test(cla_t[,i], input$map_loaded$CH4_logn1, method = "pearson")
  cor.df.cla$Class[i] <- names(cla_t)[i]
  cor.df.cla$Class_r[i] <- round(m$estimate, digits = 2)
}

# Order
ord <- summarize_taxonomy(input, level = 4, relative = F, report_higher_tax = FALSE)
ord_t <- as.data.frame(t(ord))
cor.df.ord <- as.data.frame(matrix(NA, nrow = ncol(ord_t), ncol = 2)) %>%
  set_names(c("Order", "Order_r"))
for (i in 1:ncol(ord_t)) {
  m <- cor.test(ord_t[,i], input$map_loaded$CH4_logn1, method = "pearson")
  cor.df.ord$Order[i] <- names(ord_t)[i]
  cor.df.ord$Order_r[i] <- round(m$estimate, digits = 2)
}

# Family
fam <- summarize_taxonomy(input, level = 5, relative = F, report_higher_tax = FALSE)
fam_t <- as.data.frame(t(fam))
cor.df.fam <- as.data.frame(matrix(NA, nrow = ncol(fam_t), ncol = 2)) %>%
  set_names(c("Family", "Family_r"))
for (i in 1:ncol(fam_t)) {
  m <- cor.test(fam_t[,i], input$map_loaded$CH4_logn1, method = "pearson")
  cor.df.fam$Family[i] <- names(fam_t)[i]
  cor.df.fam$Family_r[i] <- round(m$estimate, digits = 2)
}

# Genus
gen <- summarize_taxonomy(input, level = 6, relative = F, report_higher_tax = FALSE)
gen_t <- as.data.frame(t(gen))
cor.df.gen <- as.data.frame(matrix(NA, nrow = ncol(gen_t), ncol = 2)) %>%
  set_names(c("Genus", "Genus_r"))
for (i in 1:ncol(gen_t)) {
  m <- cor.test(gen_t[,i], input$map_loaded$CH4_logn1, method = "pearson")
  cor.df.gen$Genus[i] <- names(gen_t)[i]
  cor.df.gen$Genus_r[i] <- round(m$estimate, digits = 2)
}



#### Tables ####
# Make Table S8 and Wyatt table for Figure 3
tax <- input$taxonomy_loaded %>%
  set_names(c("Domain", "Phylum", "Class", "Order", "Family", "Genus")) %>%
  filter(Class != "PH") %>%
  group_by(Genus) %>%
  slice(n = 1)
TableS8 <- cor.df.gen %>%
  left_join(., tax, by = "Genus") %>%
  left_join(., cor.df.fam, by = "Family") %>%
  left_join(., cor.df.ord, by = "Order") %>%
  left_join(., cor.df.cla, by = "Class") %>%
  left_join(., cor.df.phy, by = "Phylum") %>%
  select(Domain, Phylum, Class, Order, Family, Genus, 
         Phylum_r, Class_r, Order_r, Family_r, Genus_r) %>%
  arrange(Domain, Phylum, Class, Order, Family, Genus)
# write.table(TableS8, "TableS8.txt", sep = "\t", row.names = F)
# Careful - a few names got changed to dates
# Row 56, 282, 572

# Now need to mimic Wyatt's output to feed into Figure 3
tax2 <- otu_V %>%
  select(169:176) %>%
  filter(Class != "PH") %>%
  group_by(Genus) %>%
  slice(n = 1) %>%
  filter(Genus %in% tax$Genus) %>%
  select(Consensus.lineage, Taxonomy, Genus) 
wyatt_table <- TableS8 %>%
  left_join(., tax2, by = "Genus") %>%
  rename(Kingdom = Domain) %>%
  mutate(Kingdom_r = ifelse(Kingdom == "Archaea",
                            0.53,
                            -0.53)) %>%
  mutate(Taxonomy_r = ifelse(Phylum == "Proteobacteria", 
                             Class_r, 
                             Phylum_r)) %>%
  select(Consensus.lineage, Kingdom, Phylum, Class, Order, Family, Genus, Taxonomy,
         Kingdom_r, Phylum_r, Class_r, Order_r, Family_r, Genus_r, Taxonomy_r)
#write.table(wyatt_table, "correlations/CH4corrRanks_Silva_genus_CB.txt", sep = "\t")



#### Delta ####
# Now need to do for the Delta for Figure S9
# Filter to Delta and rerun same code as above
input_d <- filter_data(input,
                       filter_cat = "SALTgroup",
                       keep_vals = c("FW", "Oligo"))

#### _Correlations ####
# Domain
dom <- summarize_taxonomy(input_d, level = 1, relative = FALSE, report_higher_tax = FALSE)
dom_t <- as.data.frame(t(dom))
cor.df.dom <- as.data.frame(matrix(NA, nrow = ncol(dom_t), ncol = 2)) %>%
  set_names(c("Domain", "Domain_r"))
for (i in 1:ncol(dom_t)) {
  m <- cor.test(dom_t[,i], input_d$map_loaded$CH4_logn1, method = "pearson")
  cor.df.dom$Domain[i] <- names(dom_t)[i]
  cor.df.dom$Domain_r[i] <- round(m$estimate, digits = 2)
}
cor.test(dom_t$Archaea, input_d$map_loaded$CH4_logn1, method = "pearson")
cor.test(dom_t$Bacteria, input_d$map_loaded$CH4_logn1, method = "pearson")

# Phylum
phy <- summarize_taxonomy(input_d, level = 2, relative = F, report_higher_tax = FALSE)
phy_t <- as.data.frame(t(phy))
cor.df.phy <- as.data.frame(matrix(NA, nrow = ncol(phy_t), ncol = 2)) %>%
  set_names(c("Phylum", "Phylum_r"))
for (i in 1:ncol(phy_t)) {
  m <- cor.test(phy_t[,i], input_d$map_loaded$CH4_logn1, method = "pearson")
  cor.df.phy$Phylum[i] <- names(phy_t)[i]
  cor.df.phy$Phylum_r[i] <- round(m$estimate, digits = 2)
}

# Class
cla <- summarize_taxonomy(input_d, level = 3, relative = F, report_higher_tax = FALSE)
cla_t <- as.data.frame(t(cla))
cor.df.cla <- as.data.frame(matrix(NA, nrow = ncol(cla_t), ncol = 2)) %>%
  set_names(c("Class", "Class_r"))
for (i in 1:ncol(cla_t)) {
  m <- cor.test(cla_t[,i], input_d$map_loaded$CH4_logn1, method = "pearson")
  cor.df.cla$Class[i] <- names(cla_t)[i]
  cor.df.cla$Class_r[i] <- round(m$estimate, digits = 2)
}

# Order
ord <- summarize_taxonomy(input_d, level = 4, relative = F, report_higher_tax = FALSE)
ord_t <- as.data.frame(t(ord))
cor.df.ord <- as.data.frame(matrix(NA, nrow = ncol(ord_t), ncol = 2)) %>%
  set_names(c("Order", "Order_r"))
for (i in 1:ncol(ord_t)) {
  m <- cor.test(ord_t[,i], input_d$map_loaded$CH4_logn1, method = "pearson")
  cor.df.ord$Order[i] <- names(ord_t)[i]
  cor.df.ord$Order_r[i] <- round(m$estimate, digits = 2)
}

# Family
fam <- summarize_taxonomy(input_d, level = 5, relative = F, report_higher_tax = FALSE)
fam_t <- as.data.frame(t(fam))
cor.df.fam <- as.data.frame(matrix(NA, nrow = ncol(fam_t), ncol = 2)) %>%
  set_names(c("Family", "Family_r"))
for (i in 1:ncol(fam_t)) {
  m <- cor.test(fam_t[,i], input_d$map_loaded$CH4_logn1, method = "pearson")
  cor.df.fam$Family[i] <- names(fam_t)[i]
  cor.df.fam$Family_r[i] <- round(m$estimate, digits = 2)
}

# Genus
gen <- summarize_taxonomy(input_d, level = 6, relative = F, report_higher_tax = FALSE)
gen_t <- as.data.frame(t(gen))
cor.df.gen <- as.data.frame(matrix(NA, nrow = ncol(gen_t), ncol = 2)) %>%
  set_names(c("Genus", "Genus_r"))
for (i in 1:ncol(gen_t)) {
  m <- cor.test(gen_t[,i], input_d$map_loaded$CH4_logn1, method = "pearson")
  cor.df.gen$Genus[i] <- names(gen_t)[i]
  cor.df.gen$Genus_r[i] <- round(m$estimate, digits = 2)
}

#### _Tables ####
# Don't need to make Table S8, just need to mimic Wyatt's output to feed into Figure S9
tax <- input_d$taxonomy_loaded %>%
  set_names(c("Domain", "Phylum", "Class", "Order", "Family", "Genus")) %>%
  filter(Class != "PH") %>%
  group_by(Genus) %>%
  slice(n = 1)
tax2 <- otu_V %>%
  select(169:176) %>%
  filter(Class != "PH") %>%
  group_by(Genus) %>%
  slice(n = 1) %>%
  filter(Genus %in% tax$Genus) %>%
  select(Consensus.lineage, Taxonomy, Genus) 
wyatt_table <- cor.df.gen %>%
  left_join(., tax, by = "Genus") %>%
  left_join(., cor.df.fam, by = "Family") %>%
  left_join(., cor.df.ord, by = "Order") %>%
  left_join(., cor.df.cla, by = "Class") %>%
  left_join(., cor.df.phy, by = "Phylum") %>%
  select(Domain, Phylum, Class, Order, Family, Genus, 
         Phylum_r, Class_r, Order_r, Family_r, Genus_r) %>%
  arrange(Domain, Phylum, Class, Order, Family, Genus) %>%
  left_join(., tax2, by = "Genus") %>%
  rename(Kingdom = Domain) %>%
  mutate(Kingdom_r = ifelse(Kingdom == "Archaea",
                            0.46,
                            -0.46)) %>%
  mutate(Taxonomy_r = ifelse(Phylum == "Proteobacteria", 
                             Class_r, 
                             Phylum_r)) %>%
  select(Consensus.lineage, Kingdom, Phylum, Class, Order, Family, Genus, Taxonomy,
         Kingdom_r, Phylum_r, Class_r, Order_r, Family_r, Genus_r, Taxonomy_r)
#write.table(wyatt_table, "correlations/CH4corrRanks_Silva_Delta_genus_CB.txt", sep = "\t")



#### S10 Genera ####
# Need to get log2 guild genera abundances for LASSO for Figure S10b
# For Delta only
# Need to get guild table to merge those guilds with the genera table from above
# Filter to just significant LASSO guilds
# Then further filter to top 20 genera of those
setwd("~/Documents/GitHub/SF_microbe_methane/guild_analysis/")
sig_guilds <- c("SRB", "MOB_I", "MOB_IIa", "AOB", "FeRB", "Anamx", "MOB_II",
                "SOxB", "CH4_ac", "FeOB")
guilds <- read.table("Silva_OTU_Guild_taxa_counts.txt", sep = "\t") %>%
  select(OTU, Guild, input_d$map_loaded$Sample) %>%
  filter(Guild %in% sig_guilds)
rs <- data.frame("rs" = rowSums(guilds[,3:75]),
                 "OTU" = guilds$OTU) %>%
  filter(rs > 0)
input_d_g <- filter_taxa_from_input(input = input_d,
                                    taxa_IDs_to_keep = rs$OTU)
gen_g <- summarize_taxonomy(input_d_g, level = 6, relative = F, report_higher_tax = F)  
top20 <- plot_taxa_bars(gen_g,
                        input_d_g$map_loaded,
                        "Sample",
                        num_taxa = 20,
                        data_only = TRUE)
levels(top20$taxon)[1:20]
gen_g_top <- gen_g %>%
  filter(row.names(.) %in% levels(top20$taxon)[1:20])
gen_g_top[, 1:73] <- log(gen_g_top[1:73], 2)
gen_g_top_t <- gen_g_top %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample") %>%
  select(Sample, everything()) %>%
  mutate_all(function(x) ifelse(is.infinite(x), 0, x))
#write.table(gen_g_top_t, 
#            "~/Documents/GitHub/SF_microbe_methane/data/SEM_data/SEM_base_log2_10_delta_genera_soil_data.txt", 
#            sep = "\t")

# List of Genera and their guilds (for coloring)
gen_otu <- input_d_g$taxonomy_loaded %>%
  rownames_to_column(var = "OTU") %>%
  select(OTU, taxonomy6)
guild_otu <- guilds %>%
  select(OTU, Guild) %>%
  filter(OTU %in% gen_otu$OTU)
gen_guild <- gen_otu %>%
  left_join(., guild_otu, by = "OTU") %>%
  select(taxonomy6, Guild) %>%
  group_by(taxonomy6) %>%
  slice(n = 1) %>%
  set_names(c("Genus", "Guild")) %>%
  filter(Genus %in% levels(top20$taxon)[1:20])
#write.table(gen_guild, "Top20_SigGuild_Genera.txt", sep = "\t")

# Even though wanted top 20 genera in Delta, need to make table of those for all data too?
# In original LASSO module, load all data and Delta data

#### Figure S5e-f ####
# Need to import guilds, Salinity, SO4, CH4, CO2 and log transform
# Make correlation matrix
# Plot heatmap of Pearson r values with guild colorbars
# Do for All sites and Delta sites
# Need RdBu matplotlib color palette

#### _All ####
setwd("~/Documents/GitHub/SF_microbe_methane/taxon_analysis/")

# Data Imports From Wyatt
# Import OTU Table
SilvaOTUs <- read.table("../silvaOTUs/Silva_OTU_VSTcpm.txt", sep='\t')
row.names(SilvaOTUs) <- SilvaOTUs[,"OTU"]                                
SilvaOTUs <- SilvaOTUs[,-1]
otu_V <- SilvaOTUs

# Sort OTU table                                                                      
otu_V <- otu_V[order(otu_V$Consensus.lineage),]                                                        

# Make new top level plotting var (should be in PRE-PROCESS ? )
otu_V$Taxonomy <- ifelse(otu_V$Phylum == "Proteobacteria", paste(otu_V$Class), paste(otu_V$Phylum))
otu_V <- data.frame(otu_V)

# METADATA
metaDB <-read.table("../data/meta/SF_sal_meta_FIX3.txt", sep="\t", header=TRUE)

# Prune metadata to only iTag samples
# Get Sample names in OTU table               
OTU_samps <- data.frame('Sample'= colnames(otu_V))

# Get aggregated guild counts by sample
guilds <- read.table("../guild_analysis/Silva_OTU_Guild_abundT_counts.txt", sep = "\t") %>%
  select(1:17) %>%
  rownames_to_column(var = "Sample")

# Merge site order and Samples,  log transform, merge guilds
Meta_iTag <- merge(metaDB, OTU_samps, by='Sample') %>%
  select(Sample, Location, SALTgroup, CH4_ug_m2_h, SO4, Salinity.x, CO2_mg_m2_h) %>%
  column_to_rownames(var = "Sample") %>%
  mutate(Sample = rownames(.)) %>%
  mutate(Salinity = log10(Salinity.x),
         SO4 = log10(SO4),
         CH4_ug_m2_h = log10(CH4_ug_m2_h - (1.05*min(CH4_ug_m2_h))),
         CO2_mg_m2_h = log10(CO2_mg_m2_h)) %>%
  left_join(., guilds, by = "Sample") %>%
  select(-Salinity.x) %>%
  select(Sample, Location, SALTgroup, Salinity, SO4, CH4_ug_m2_h, CO2_mg_m2_h,
         CH4_mix, CH4_me, everything())

# Get colors for colorbar. Note will need to make the env. vars white
Guild_cols <- read.table("../data/colors/Guild_color_palette.txt", sep='\t') %>%
  select(Guild, G_index, color) %>%
  set_names(c("Guild", "Index", "color")) %>%
  mutate(Index = rev(Index)) %>%
  add_row(Guild = "CH4_me", Index = 16, color = "#FDC086") %>%
  add_row(Guild = "CH4_mix", Index = 17, color = "#FFFF99") %>%
  add_row(Guild = "CO2_mg_m2_h", Index = 18, color = "white") %>%
  add_row(Guild = "CH4_ug_m2_h", Index = 19, color = "white") %>%
  add_row(Guild = "SO4", Index = 20, color = "white") %>%
  add_row(Guild = "Salinity", Index = 21, color = "white") %>%
  arrange(-Index) %>%
  mutate(Guild = gsub("MeOB", "ANME", Guild))

# Cor matrix
vars <- Meta_iTag %>%
  select(4:24)

sum(names(vars) != Guild_cols$Guild)

cm <- cor(vars)
cm[cm == 1.00] <- 1

mycolors <- colorRampPalette(brewer.pal(11, "RdBu"))(100)

ann_rows <- data.frame(row.names = rownames(cm), 
                       "colorbar.r" = Guild_cols$Guild)
ann_cols <- data.frame(row.names = rownames(cm), 
                       "colorbar.c" = Guild_cols$Guild)
ann_colors <- list(colorbar.r = c(Salinity = Guild_cols$color[1],
                                  SO4 = Guild_cols$color[2],
                                  CH4_ug_m2_h = Guild_cols$color[3],
                                  CO2_mg_m2_h = Guild_cols$color[4],
                                  CH4_mix = Guild_cols$color[5],
                                  CH4_me = Guild_cols$color[6],
                                  CH4_H2 = Guild_cols$color[7],
                                  CH4_ac = Guild_cols$color[8],
                                  MOB_I = Guild_cols$color[9],
                                  MOB_II = Guild_cols$color[10],
                                  MOB_IIa = Guild_cols$color[11],
                                  ANME = Guild_cols$color[12],
                                  AOA = Guild_cols$color[13],
                                  AOB = Guild_cols$color[14],
                                  NOB = Guild_cols$color[15],
                                  Anamx = Guild_cols$color[16],
                                  SOxB = Guild_cols$color[17],
                                  SRB_syn = Guild_cols$color[18],
                                  SRB = Guild_cols$color[19],
                                  FeOB = Guild_cols$color[20],
                                  FeRB = Guild_cols$color[21]),
                   colorbar.c = c(Salinity = Guild_cols$color[1],
                                  SO4 = Guild_cols$color[2],
                                  CH4_ug_m2_h = Guild_cols$color[3],
                                  CO2_mg_m2_h = Guild_cols$color[4],
                                  CH4_mix = Guild_cols$color[5],
                                  CH4_me = Guild_cols$color[6],
                                  CH4_H2 = Guild_cols$color[7],
                                  CH4_ac = Guild_cols$color[8],
                                  MOB_I = Guild_cols$color[9],
                                  MOB_II = Guild_cols$color[10],
                                  MOB_IIa = Guild_cols$color[11],
                                  ANME = Guild_cols$color[12],
                                  AOA = Guild_cols$color[13],
                                  AOB = Guild_cols$color[14],
                                  NOB = Guild_cols$color[15],
                                  Anamx = Guild_cols$color[16],
                                  SOxB = Guild_cols$color[17],
                                  SRB_syn = Guild_cols$color[18],
                                  SRB = Guild_cols$color[19],
                                  FeOB = Guild_cols$color[20],
                                  FeRB = Guild_cols$color[21]))

pheatmap(cm,
         legend = F,
         scale = "none",
         color = rev(mycolors),
         cluster_cols = F,
         cluster_rows = F,
         angle_col = 90,
         border_color = NA,
         annotation_row = ann_rows,
         annotation_col = ann_cols,
         annotation_colors = ann_colors,
         annotation_legend = F,
         annotation_names_row = F,
         annotation_names_col = F,
         display_numbers = T,
         number_format = ifelse(cm == 1, "%.0f", "%.2f"),
         fontsize = 3,
         fontsize_number = 2,
         number_color = ifelse(abs(cm) > 0.5, "white", "black"),
         filename = "../figs/FigureS5e.png", width = 2.25, height = 2.25)

e <- pheatmap(cm,
              legend = F,
              scale = "none",
              color = rev(mycolors),
              cluster_cols = F,
              cluster_rows = F,
              angle_col = 90,
              border_color = NA,
              annotation_row = ann_rows,
              annotation_col = ann_cols,
              annotation_colors = ann_colors,
              annotation_legend = F,
              annotation_names_row = F,
              annotation_names_col = F,
              display_numbers = T,
              number_format = ifelse(cm == 1, "%.0f", "%.2f"),
              fontsize = 3,
              fontsize_number = 2,
              number_color = ifelse(abs(cm) > 0.5, "white", "black"))
save_pheatmap_pdf(e, "../figs/FigureS5e.pdf")


#### _Delta ####
Meta_iTag_d <- Meta_iTag %>%
  filter(SALTgroup == "FW" | SALTgroup == "Oligo")

# Cor matrix
vars_d <- Meta_iTag_d %>%
  select(4:24)

sum(names(vars_d) != Guild_cols$Guild)

cm_d <- cor(vars_d)

pheatmap(cm_d,
         legend = F,
         scale = "none",
         color = rev(mycolors),
         cluster_cols = F,
         cluster_rows = F,
         angle_col = 90,
         border_color = NA,
         annotation_row = ann_rows,
         annotation_col = ann_cols,
         annotation_colors = ann_colors,
         annotation_legend = F,
         annotation_names_row = F,
         annotation_names_col = F,
         display_numbers = T,
         number_format = ifelse(cm_d == 1, "%.0f", "%.2f"),
         fontsize = 3,
         fontsize_number = 2,
         number_color = ifelse(abs(cm_d) > 0.5, "white", "black"),
         filename = "../figs/FigureS5f.png", width = 2.25, height = 2.25)

f <- pheatmap(cm_d,
              legend = F,
              scale = "none",
              color = rev(mycolors),
              cluster_cols = F,
              cluster_rows = F,
              angle_col = 90,
              border_color = NA,
              annotation_row = ann_rows,
              annotation_col = ann_cols,
              annotation_colors = ann_colors,
              annotation_legend = F,
              annotation_names_row = F,
              annotation_names_col = F,
              display_numbers = T,
              number_format = ifelse(cm_d == 1, "%.0f", "%.2f"),
              fontsize = 3,
              fontsize_number = 2,
              number_color = ifelse(abs(cm_d) > 0.5, "white", "black"))
save_pheatmap_pdf(f, "../figs/FigureS5f.pdf")
