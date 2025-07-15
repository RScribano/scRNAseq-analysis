# packages ####
#tidyverse
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(harmony)

#single cell analysis
library(scran)
library(scater) 
library(tidyomics)
library(tidySingleCellExperiment)
library(bluster)
library(clusterProfiler)
library(org.Hs.eg.db)
library(EnhancedVolcano)

# pseudobulk DE
library(edgeR)
library(Matrix)
library(SummarizedExperiment)
library(clusterProfiler)


# These steps were carried out on the CCDL server. SCE object obtained by
# 1. converting the Seurat object to SCE
# 2. Obtaining 2000 highly variable genes
# 3. Running PCA
# 4. Running UMAP (seed == 42) with dimereduct == PCA
# 5. Adding metadata to colData (individual, timepoint, MRD_risk, relapse)

# import data ####

# These steps were carried out on the CCDL server. SCE object obtained by
# 1. converting the Seurat object to SCE
# 2. Obtaining 2000 highly variable genes
# 3. Running PCA
# 4. Running UMAP (seed == 42) with dimereduct == PCA
# 5. Adding metadata to colData (individual, timepoint, MRD_risk, relapse)

rds_path <- file.path("scRNAseq data","20250415_SCE.RDS")
sce<-readRDS(rds_path)


# Data wrangling####
coldata_df <- as.data.frame(colData(sce))

# new column to distinguish MRD and healthy

coldata_df <- coldata_df |>
  mutate(
    sample_type = case_when(
      individual == "Healthy" ~ "Healthy",
      timepoint %in% c("D8", "D15", "D33", "D78") ~ "MRD timepoints",
      timepoint == "Dx" ~ "Diagnosis",
      timepoint == "Rx" ~ "Relapse"
      
    )
  ) |>
  mutate(
    sample_type = factor(sample_type, levels = c("Healthy","Diagnosis", "MRD timepoints", "Relapse"),
                         ordered= TRUE)
  )

colData(sce) <- DataFrame(coldata_df)

# Fix labels MRD_43 ####
# in cartridge F tags 6 and 7 have been swapped around so that Dx is swapped with D33. we need to change that 

# 1. Swap sample tags
tags <- sce$ident[sce$orig.ident == "cart_f"]

tags_new <- dplyr::case_when(tags == "SampleTag07_hs" ~ "SampleTag06_hs",
                             tags == "SampleTag06_hs" ~ "SampleTag07_hs")

sce$ident[sce$orig.ident == "cart_f"] <- tags_new

# Data wrangling####
coldata_df <- as.data.frame(colData(sce))

# new column to distinguish MRD and healthy

coldata_df <- coldata_df |>
  mutate(
    sample_type = case_when(
      individual == "Healthy" ~ "Healthy",
      timepoint %in% c("D8", "D15", "D33", "D78") ~ "MRD timepoints",
      timepoint == "Dx" ~ "Diagnosis",
      timepoint == "Rx" ~ "Relapse"
      
    )
  ) |>
  mutate(
    sample_type = factor(sample_type, levels = c("Healthy","Diagnosis", "MRD timepoints", "Relapse"),
                         ordered= TRUE)
  )

colData(sce) <- DataFrame(coldata_df)

# distribution of total counts per cartridge ####
ggplot(as.data.frame(colData(sce)), aes(x = orig.ident,y = log10(nCount_RNA) ))+
  geom_violin()

# import gates ####
gated_path <- file.path("scRNAseq data", "gated", "20250423_concat.csv")
gated_df <- read.csv(gated_path) 

# restore sample names 
gated_df <- gated_df |>
  dplyr::select(-Sample_Name) |>
  dplyr::rename(Sample_Name = OmiqFileIndex) |>
  dplyr::mutate(Sample_Name = gsub(x = Sample_Name ,pattern = "\\.csv", replacement = ""))


# order sce by sample name (the order of cells within samples is the same as in the gated file)

ord <- order(factor(sce$Sample_Name, levels = unique(gated_df$Sample_Name)))
sce <- sce[,ord]


# add the gate information to the SCE
sce$gate <- gated_df$OmiqFilter

# find marker genes between blasts (annotated on proteins) and B-cells ####

sce_b <- sce[,sce$Cell_Type_Experimental == "B"]

# for diagnosis
markers_dx <- findMarkers(sce_b[, sce_b$sample_type == "Diagnosis"], test.type = "t", 
                       direction = "up", 
                       groups = sce_b[, sce_b$sample_type == "Diagnosis"]$gate)

# MRD timepoints
markers_MRD <- findMarkers(sce_b[, sce_b$sample_type == "MRD timepoints"], test.type = "t", 
                          direction = "up", 
                          groups = sce_b[, sce_b$sample_type == "MRD timepoints"]$gate)

# Relapse
markers_Rx <- findMarkers(sce_b[, sce_b$sample_type == "Relapse"], test.type = "t", 
                           direction = "up", 
                           groups = sce_b[, sce_b$sample_type == "Relapse"]$gate)



# store and save as adataframes for later use
markers_dx_df <- as.data.frame(markers_dx[[1]])
markers_mrd_df <- as.data.frame(markers_MRD[[1]])
markers_rx_df <- as.data.frame(markers_Rx[[1]])


# transform gene names to column

markers_dx_df |>
  rownames_to_column(var = "gene") -> markers_dx_df

markers_mrd_df |>
  rownames_to_column(var = "gene") -> markers_mrd_df

markers_rx_df |>
  rownames_to_column(var = "gene") -> markers_rx_df

# save csv files
write.csv(markers_dx_df, "20250710_markers_dx.csv")
write.csv(markers_mrd_df, "20250710_markers_mrd.csv")
write.csv(markers_rx_df, "20250710_markers_rx.csv")

hist(markers_Rx[[1]]$p.value)

# plot differentially expressed at ANY timepoint ####
markers_dx_df |>
  filter(FDR<0.05, summary.logFC >0.5) |>
  pull(gene) -> genes_dx

markers_mrd_df |>
  filter(FDR<0.05, summary.logFC >0.5) |>
  pull(gene) -> genes_mrd

markers_rx_df |>
  filter(FDR<0.05, summary.logFC >0.5) |>
  pull(gene) -> genes_rx

# create common list

deg_any <- 
  c(genes_dx, genes_mrd, genes_rx) |>
  unique()


# boxplots markers at any timepoint
sce_b |>
  join_features(features = deg_any) |>
  filter(sample_type == "Diagnosis") |>
  ggplot(aes(x = .feature , y=.abundance_logcounts)) +
  geom_boxplot(position = position_dodge(), aes(fill = gate))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_cartesian(ylim = c(0,4))


sce_b |>
  join_features(features = deg_any) |>
  filter(sample_type == "MRD timepoints") |>
  ggplot(aes(x = .feature , y=.abundance_logcounts)) +
  geom_boxplot(position = position_dodge(), aes(fill = gate))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_cartesian(ylim = c(0,4))

sce_b |>
  join_features(features = deg_any) |>
  filter(sample_type == "Relapse") |>
  ggplot(aes(x = .feature , y=.abundance_logcounts)) +
  geom_boxplot(position = position_dodge(), aes(fill = gate))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_cartesian(ylim = c(0,4))
  

# pathway analysis markers ####


# diagnosis
sig_dx <- markers_dx_df |>
  filter(FDR < 0.05, logFC.Unfiltered>0.3)

genes_dx <-  bitr(sig_dx$gene, fromType = "SYMBOL", 
                          toType = "ENTREZID", 
                          OrgDb = org.Hs.eg.db)

sig_dx <- sig_dx |>
  rename(SYMBOL = gene) |>
  left_join(genes_dx, by = "SYMBOL")

  
ego <- enrichGO(gene         = sig_dx$ENTREZID,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",          # "BP"=Biological Process, "MF", "CC"
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2,
                readable     = TRUE)


# pseudobulk differential expression analysis - RELAPSED ####

# in this section we are investigating in a more statistical sound way the 
# differences between different timepoints. Note that we are doing this on
# all cells, not only those that relapsed



# create column for timepoint_individual 
sce_b$timepoint_patient <- paste(sce_b$sample_type, 
                                        sce_b$relapse, sce_b$individual)

# aggregate pseudobulk counts
pb <- scuttle::aggregateAcrossCells(sce_b[, sce_b$sample_type != "Healthy"],
                                    ids = sce_b[, sce_b$sample_type != "Healthy"]$timepoint_patient )
counts <- counts(pb)


# set up differential expression analysis

#1 = diagnosis
#2= MRD timepoints
#3 = Relapse timepoint


sce_b <- sce_b[,sce_b$sample_type != "Healthy"]
sce_b_relapsed <- sce_b[, sce_b$relapse == "yes"&
                      sce_b$gate == "Leukemia"]


# Add QC for this subset

sce_b_relapsed$sum <- NULL
sce_b_relapsed$detected <- NULL
sce_b_relapsed$total <- NULL
sce_b_relapsed$subsets_mito_sum <- NULL
sce_b_relapsed$subsets_mito_detected <- NULL
sce_b_relapsed$subsets_mito_sum <- NULL

# QC
sce_b_relapsed <- scater::addPerCellQC(sce_b_relapsed)


# create column for timepoint_individual 
sce_b_relapsed$timepoint_patient <- paste(sce_b_relapsed$sample_type, 
                                          sce_b_relapsed$individual)




# aggregate pseudobulk counts
pb <- scuttle::aggregateAcrossCells(sce_b_relapsed, ids = sce_b_relapsed$timepoint_patient )
counts <- counts(pb)


# set up differential expression analysis ####

#1 = diagnosis
#2= MRD timepoints
#3 = Relapse timepoint

group <- factor(c(1,1,1,1,1,2,2,2,2,3,3,3,3))
patient <- factor(c('MRD_01', 'MRD_10', 'MRD_11', 'MRD_27',
                    'MRD_45', 'MRD_01', 'MRD_10', 'MRD_11', 
                    'MRD_27', 'MRD_01', 'MRD_10', 'MRD_11',
                    'MRD_27'))

# create dgelist 
y <- DGEList(counts = counts,
             group = group)

# filter lowly expressed genes

keep <- filterByExpr(y)
y<- y[keep, , keep.lib.sizes = FALSE]


# create design
design <- model.matrix(~patient + group)

# normalize
y <- normLibSizes(y)

# estimate dispersion
y <- estimateDisp(y, design)

# fit model
fit <- glmFit(y, design)

# comparisons

# To test MRD vs Diagnosis
lrt_mrd <- glmLRT(fit, coef = "group2")

# To test Relapse vs Diagnosis
lrt_relapse <- glmLRT(fit, coef = "group3")


# MRD vs relapse
contrast <- makeContrasts(group3 - group2, levels = design)
lrt_mrd_relapse <- glmLRT(fit, contrast = contrast)


# contrasts
topde_mrd <- as.data.frame(topTags(lrt_mrd,adjust.method = 'hochberg',n = Inf))
topde_rx <- as.data.frame(topTags(lrt_relapse,adjust.method = 'hochberg',n = Inf))
topde_mrd_rx <- as.data.frame(topTags(lrt_mrd_relapse,adjust.method = 'hochberg',n = Inf))

# volcano plots

EnhancedVolcano(toptable = topde_mrd,
                lab = rownames(topde_mrd),
                x= "logFC",
                y = "PValue",
                pCutoff = 0.05,
                pCutoffCol = "FWER",
                FCcutoff = 1,
                title = 'Relapsed - Dx vs. MRD')


EnhancedVolcano(toptable = topde_rx,
                lab = rownames(topde_rx),
                x= "logFC",
                y = "PValue",
                pCutoff = 0.05,
                pCutoffCol = "FWER",
                FCcutoff = 1,
                title = 'Relapsed - Dx vs. Rx')


EnhancedVolcano(toptable = topde_mrd_rx,
                lab = rownames(topde_mrd_rx),
                x= "logFC",
                y = "PValue",
                pCutoff = 0.05,
                pCutoffCol = "FWER",
                FCcutoff = 1,
                title = 'Relapsed - MRD vs. Rx')


# now i want to plot these genes at single cell level in relapsed patients
deg_relapsed_any <- c("C16orf54", "ENSG00000287784", "ENSG00000287784", "MS4A1")

sce_b |>
  join_features(features = deg_relapsed_any) |>
  filter(gate == "Leukemia",
         relapse == "yes") |>
  ggplot(aes(x = .feature , y=.abundance_logcounts)) +
  geom_jitter(position = position_jitter(),alpha = 0.4)+
  geom_violin(aes(fill = sample_type)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  

sce_b |>
  join_features(features = deg_relapsed_any) |>
  filter(gate == "Leukemia",
         relapse == "yes") |>
  group_by(.feature, sample_type) |>
  summarise(mean_expr = mean(.abundance_logcounts, na.rm = TRUE),
            sem = sd(.abundance_logcounts)/n()) |>
  ggplot(aes(x = .feature , y=mean_expr)) +
  geom_bar(stat = "identity", aes( fill = sample_type), 
           position = position_dodge()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text("mean logcounts"),
        axis.title.x = element_text("Gene ID"))


# Dimension reduction on blasts ####

sce_blasts <- sce_b[, sce_b$gate == "Leukemia"]

highvar <- getTopHVGs(sce_blasts)

reducedDim(sce_blasts, "PCA") <- NULL
reducedDim(sce_blasts, "UMAP") <- NULL

# PCA
sce_blasts <- runPCA(sce_blasts, 
                     ncomponents = 20, 
                     ntop = 500,
                     exprs_values = "logcounts")



# UMAP

sce_blasts <- runUMAP(sce_blasts, 
                    dimred = "PCA")


sce_blasts |>
  ggplot(aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = relapse)) +
  facet_wrap(~sample_type) +
  theme_bw(base_size = 22)



# tSNE
sce_blasts<-runTSNE(sce_blasts,
                    dimred = "PCA")

sce_blasts |>
  ggplot(aes(x = TSNE1, y = TSNE2)) +
  geom_point(aes(color = individual)) +
  facet_wrap(~sample_type) +
  theme_bw(base_size = 22)


# clustering ON BLASTS only ####

blasts_pca <- reducedDim(sce_blasts, "PCA")

cluster.louvein <- clusterCells(sce_blasts, use.dimred = "PCA",
                                BLUSPARAM = SNNGraphParam(cluster.fun = "louvain"))

cluster.leiden <- clusterCells(sce_blasts, use.dimred = "PCA",
                                BLUSPARAM = SNNGraphParam(cluster.fun = "leiden"))

#piechart total count of cells per cluster 

label_counts <- as.data.frame(table(sce_blasts$label)) %>%
  rename(label = Var1, count = Freq)

ggplot(label_counts, aes(x = "", y = count, fill = label)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "Cluster Label Distribution (Pie Chart)") +
  theme_void() +
  theme(legend.title = element_blank())


# umap colored by labels
colLabels(sce_blasts) <- cluster.louvein

sce_blasts |>
  ggplot(aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = label)) +
  facet_wrap(~sample_type) +
  theme_bw(base_size = 22)

# contingency table invidual vs. cluster label

tbl <- table(sce_blasts$individual, sce_blasts$label) |>
  as.matrix() |>
  pheatmap::pheatmap()


# which plot is the cool one? abundance overtime

clusters_count <-sce_blasts |>
  group_by(sample_type, label, relapse) |>
  summarise(abundance =  n(), .groups = "drop")

# plot proportion
clusters_count |>
  group_by(sample_type, relapse) |>
  mutate(prop = abundance/sum(abundance)) |>
  ggplot(aes(x = sample_type, y = prop)) +
  geom_bar(stat = "identity",
    position = position_dodge(), aes(fill = relapse)) +
  facet_wrap(~label, scales = "free")+
  theme_bw()


# find markers

# cluster 1 

topgene_clusters <- findMarkers(sce_blasts,
                                test.type = "wilcox", 
                                direction = "up", lfc = 1,
                                pval.type = "any")

top_markers <- topgene_clusters[[3]] |>
  as.data.frame() |>
  filter(FDR< 0.05) |>
  rownames()


plotExpression(sce_blasts, features = top_markers[16:31], x = "label", color_by = "label")

top_markers <- topgene_clusters[[3]] |>
  as.data.frame() |>
  filter(FDR< 0.05) |>
  write.csv("20250714_Cluster3.csv")




 # explorative plots
sce_blasts |>
  join_features(features = head(highvar)) |>
  ggplot(aes(x = .feature , y=.abundance_logcounts)) +
  geom_violin(aes(fill = relapse), position = position_dodge()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text("mean logcounts"),
        axis.title.x = element_text("Gene ID")) +
  facet_wrap(~sample_type)

sce_blasts |>
  join_features(features = head(highvar)) |>
  ggplot(aes(x = .abundance_logcounts, y = nCount_RNA)) +
  geom_point(aes(color = .feature))




