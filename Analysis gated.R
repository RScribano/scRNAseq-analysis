# packages ####
#tidyverse
library(dplyr)
library(tidyr)
library(purrr)

#single cell analysis
library(scran)
library(scater) 

#plotting
library(ggplot2)
library(ggridges)
library(scales)
library(patchwork)
library(UCell)


# import gated data ####
gated_path <- file.path("gated", "20250423_concat.csv")
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

# re-make UMAP df. Now it will retain the same order as the one matching the gated df

umap_df <- as.data.frame(reducedDim(sce, "UMAP"))
umap_df$signature_mean <- colData(sce)$signature_mean
umap_df$signature_sum <- colData(sce)$signature_sum
umap_df$individual <- colData(sce)$individual
umap_df$timepoint <- colData(sce)$timepoint
umap_df$sample_type <- colData(sce)$sample_type
umap_df$relapse <- colData(sce)$relapse
umap_df$MRD_risk <- colData(sce)$MRD_risk
umap_df$Cell_Type_Experimental <- colData(sce)$Cell_Type_Experimental
umap_df$gate <- colData(sce)$gate

# write as csv
write.csv(umap_df, "20250424_gated_umapd_df.csv")

# Basic plot gates for sanity check
umap_df |> 
  dplyr::filter(individual == "MRD_43") |>
  ggplot(aes(x=UMAP1, y=UMAP2, color = gate))+
    geom_point(alpha = 0.6, size = 0.7)+
  facet_wrap(~timepoint)

# Contingency table cell type experimental gate ####
# Comparison between the automated cell type annotation from Seurat (reference-based) and the manual gating 

table <- as.data.frame( table(umap_df$Cell_Type_Experimental, umap_df$gate))
colnames(table) <- c("Celltype.Azimuth", "Manual.gate", "N.cells")

table <- table |>
  dplyr:: group_by(Celltype.Azimuth) |>
  dplyr:: mutate( total_cells = sum(N.cells),
                  percent = 100*N.cells/total_cells)

table$percent <- round(table$percent, digits = 2)


# cell counts
ggplot(table, aes(x = Manual.gate , y = Celltype.Azimuth, fill = N.cells)) +
  geom_tile() +
  geom_text(aes(label = N.cells), color = "black") +
  scale_fill_gradient(low = "grey80", high = "steelblue") +
  theme_minimal()

# normalized

ggplot(table, aes(x = Manual.gate , y = Celltype.Azimuth, fill = percent)) +
  geom_tile() +
  geom_text(aes(label = percent), color = "black") +
  scale_fill_gradient(low = "grey80", high = "steelblue") +
  theme_minimal()


# Signature scores (Ucell) ####

# check if signature genes are in "detected" (genes detected above thresholds)

# instantiate signatures
signatures <- list(
  ComplexI= c("CYCS","CYC1", "CS","NDUFA4","NDUFC1","NDUFS2","NDUFS5","NDUFV1","NDUFB3", 
                     "NDUFA13","NDUFB8","NDUFB10","NDUFA3","NDUFB9","NDUFB11","NDUFS1","NDUFA11","NDUFS8", 
                     "NDUFA2","NDUFA9","MRPL3","MRPS21","MRPL17","MRPS30","MRPL16"),
  
  Mito_ETC = c("CYCS", "CYC1", "CS", "COX5A", "ATP5A", "ATP5B"),
  
  Pyrimidine_synthesis = c("DHODH", "CAD", "UMPS")
)

#check if genes are in the sce (=passed QC steps)
complexi_detected <- signatures$ComplexI[signatures$ComplexI%in% rownames(sce)]
mito_detected <- signatures$Mito_ETC[signatures$Mito_ETC%in% rownames(sce)]
pyr_detected <- signatures$Pyrimidine_synthesis[signatures$Pyrimidine_synthesis%in% rownames(sce)]



# some genes in the mito ( "ATP5A", "ATP5B") were not detected. edit the signautres accordingly
signatures$Mito_ETC <- mito_detected

# run UCell
sce <- ScoreSignatures_UCell(sce, features = signatures,
                             assay = 'logcounts', name = NULL)

# export as dataframe column
umap_df$UCell_signature_score <- t(as.data.frame(assay(altExp(sce, "UCell"))))


