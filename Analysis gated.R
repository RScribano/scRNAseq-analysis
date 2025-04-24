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
library(ggpubr)
library(pheatmap)


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

UCell_scores <- as.data.frame(t(as.data.frame(assay(altExp(sce, "UCell")))))


# export as datafr# export as datafr# export as dataframe column

signatures_df <- data.frame(
  
  sample_name = colData(sce)$Sample_Name,
  timepoint = colData(sce)$timepoint,
  individual = colData(sce)$individual,
  MRD_risk = colData(sce)$MRD_risk,
  sample_type = colData(sce)$sample_type,
  relapse = colData(sce)$relapse,
  gate = colData(sce)$gate,
  ComplexI = UCell_scores$ComplexI,
  Mito_ETC = UCell_scores$Mito_ETC,
  Pyrimidine_synthesis = UCell_scores$Pyrimidine_synthesis,
  stringsAsFactors = FALSE
)

write.csv(signatures_df, "20250424_signatures.csv")


# plot signature #### 

# convert to long format for better handling during plotting 

signatures_longer <- signatures_df |>
  filter(gate == "Leukemia") |>
  tidyr::pivot_longer(cols = c( "ComplexI","Mito_ETC","Pyrimidine_synthesis"),
               names_to = "signature")


# Violin plot 1 - individual, by sample type
p0<- signatures_longer |>
  dplyr::filter(sample_type != "Healthy") |>
  ggplot(aes(x=individual, y=value))+
  geom_violin(aes(fill = signature))+
  geom_jitter(size = 0.5, aes(fill = signature, color = signature), position = position_jitterdodge(), alpha = 0.2)+
  theme_gray(base_size = 22)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  facet_grid(sample_type~signature)

ggsave("20250424_gray_overall violin plot signatures.png",
       plot = p0, width = 1600,height =1100, unit = "px",
       scale = 5)


# divided by relapse 

p1<-signatures_longer |>
  dplyr::filter(
    sample_type %in% c("Diagnosis", 'MRD timepoints'),
    relapse == "no") |>
  ggplot(aes(x=individual, y=value))+
  geom_violin(aes(fill = signature))+
  geom_jitter(aes(fill = signature, color = signature), size = 0.5, position = position_jitterdodge(), alpha = 0.2)+
  theme_gray(base_size = 22)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  facet_grid(sample_type~signature)+
  ggtitle("Complete Remission")

p2 <- signatures_longer |>
  dplyr::filter(
    sample_type %in% c("Diagnosis", 'MRD timepoints'),
    relapse == "yes") |>
  ggplot(aes(x=individual, y=value))+
  geom_violin(aes(fill = signature))+
  geom_jitter(aes(fill = signature, color = signature),size = 0.5, position = position_jitterdodge(), alpha = 0.2)+
  theme_gray(base_size = 22)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  facet_grid(sample_type~signature)+
  ggtitle("Relapsed")

grid<-wrap_plots(p1,p2, ncol = 2)

ggsave("20250424_no relapse_signature score by outcome.png",
       plot = grid, width = 2100,height =1226, unit = "px",
       scale = 5)

# heatmaps single genes  ###

# the part below could DEFINETLY be automated with a function + using purrr
# but for now let's do mancave style

# one heatmpa per gene of each signature 

# Complex I signatures #####

signatures_df_leukemia <- signatures_df |>
  filter(gate == "Leukemia",
         sample_type != "Healthy")

complexi_df <-  data.frame(
  individual = signatures_df_leukemia$individual,
  timepoint = signatures_df_leukemia$timepoint,
  sample_type = signatures_df_leukemia$sample_type,
  sample_name = signatures_df_leukemia$sample_name,
  relapse = signatures_df_leukemia$relapse
)

# individual expression of genes in the signatures
for (gene in signatures$ComplexI){
  print(gene)
  complexi_df[[gene]] <- as.vector(assay(sce[, (sce$gate == "Leukemia")&(sce$sample_type != "Healthy")], "counts")[gene,]) # !! RAW COUNTS!!
}

# convert to longer
complexi_longer <- complexi_df |>
  pivot_longer(cols = 6:30,
                names_to = "genes",
               values_to = 'value')

# annotation matrix
annot <- complexi_longer |>
  filter(sample_type != "Healthy") |>
  group_by(sample_name) |>
  summarise(
    sample_name = first(sample_name),
    sample_type = first(sample_type),  # or use `mode` or any other appropriate function
    relapse = first(relapse),  # Assuming you want the first value, or aggregate them
    invididual = first(individual),
    .groups = 'drop'  # This drops the grouping after summarizing
  ) |>
  tibble::column_to_rownames("sample_name")


# mean expression per patient
complexi_longer <- complexi_longer |>
      group_by(genes, sample_name, individual) |>
      summarise(mean = mean(value), .groups = "drop") 

# convert to matrix
complexi_mat <- complexi_longer |>
  select(c("sample_name", "genes", "mean")) |>
  pivot_wider(names_from = sample_name, values_from = mean) |>
  tibble::column_to_rownames("genes") |>
  as.matrix()

# annotation colors
annot_colors <- list(
  relapse = c("yes"="red","no"="blue"),
  sample_type = c("Diagnosis" = "lightblue", "MRD timepoints" = "purple", "Relapse" = "rosybrown")
  )

#heatmap
pheatmap(complexi_mat,  # Only numeric columns for the heatmap
         scale = "none",  # Optional: scales data by columns (variables)
         clustering_distance_rows = "euclidean",  # Distance metric for rows (patients)
         clustering_distance_cols = "euclidean",  # Distance metric for columns (variables)
         clustering_method = "complete",  # Clustering method (complete linkage, etc.)
         show_rownames = TRUE,  # Show row names (patient IDs)
         show_colnames = TRUE,  # Show column names (variable names)
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Color scale
        annotation_col = annot,#annotation_row = annot,
        annotation_colors = annot_colors,
          # Define the colors for annotations
         main = "Complex I genes",  # Title of the heatmap
         fontsize_row = 8,  # Decrease font size of row names (patient IDs) to make more space
         angle_row = 45,  # Rotate row names to avoid squishing
         cellheight = 8,  # Adjust the height of each cell for better spacing of row names
         
)




