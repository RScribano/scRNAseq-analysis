# packages ####
#tidyverse
library(dplyr)
library(tidyr)
library(purrr)

# stats
library(vegan)
library(PMCMRplus)
library(FSA)

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
library(ggrepel)

# for general purpose

summary_table <- coldata_df |>
  filter(sample_type != "Healthy") |>
  select(c("sample_type", "relapse", "individual")) |>
  group_by(sample_type, relapse) |>
  summarise(
    patients = n_distinct(individual),
    .groups = "drop") |>
  as.data.frame()


summary_cells <- coldata_df |>
  filter(sample_type != "Healthy") |>
  group_by(sample_type, relapse) |>
  summarise(
    cells = n(),
    .groups = "drop"
  )
colnames(summary_table) <- c("sample_type", "relapse")

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
p0<- umap_df |>
  filter(sample_type != "Healthy") |>
  ggplot(aes(x = UMAP1, y= UMAP2))+
  geom_point(aes(fill = Cell_Type_Experimental, color = Cell_Type_Experimental), alpha = 0.6, size = 0.7)+
  theme_void(base_size = 28)+ 
  guides(color = guide_legend(override.aes = list(size = 6)))+ #scale legend
  ggtitle("Seurat-Azimuth reference-based cell type annotation")

p1 <- umap_df |>
  filter(sample_type != "Healthy") |>
  ggplot(aes(x = UMAP1, y= UMAP2))+
  geom_point(aes(fill = gate, color = gate), alpha = 0.6, size = 0.7)+
  guides(color = guide_legend(override.aes = list(size = 6)))+
  theme_void(base_size = 28)+
  ggtitle("Manual gate on CITE-seq markers")

grid <- wrap_plots(p0,p1, ncol = 1)

ggsave("20250425_celltype vs gate umap.png",
       plot = grid, width = 1800,height =1400, unit = "px",
       scale = 3)


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
p0<- ggplot(table, aes(x = Manual.gate , y = Celltype.Azimuth, fill = N.cells)) +
  geom_tile() +
  geom_text(aes(label = N.cells), color = "black") +
  scale_fill_gradient(low = "grey80", high = "steelblue") +
  ggtitle("Absolute counts (N cells)")+
  theme_minimal()

# normalized

p1<-ggplot(table, aes(x = Manual.gate , y = Celltype.Azimuth, fill = percent)) +
  geom_tile() +
  geom_text(aes(label = percent), color = "black") +
  scale_fill_gradient(low = "grey80", high = "steelblue") +
  ggtitle("Normalized by cell-type (azimuth)")+
  theme_minimal()

grid <- wrap_plots(p0,p1, ncol = 2, axis_titles = "collect")

grid

# Signature scores (Ucell) ####

# check if signature genes are in "detected" (genes detected above thresholds)

# instantiate signatures
signatures <- list(
  ComplexI= c("NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4", 
              "NDUFA5", "NDUFA6", "NDUFA7", "NDUFA8", 
              "NDUFA9", "NDUFB1", "NDUFB2", "NDUFB3", 
              "NDUFB4", "NDUFB5", "NDUFB6", "NDUFB7", 
              "NDUFB8", "NDUFB9", "NDUFC1", "NDUFC2", 
              "NDUFS1", "NDUFS2", "NDUFS3", "NDUFV1", 
              "NDUFS4", "NDUFS5", "NDUFS6", "NDUFS8", 
              "NDUFV2", "NDUFV3", "NDUFS7"),
  
  Mito_ETC = c("CYCS", "CYC1", "CS", "COX5A", "ATP5A", "ATP5B"),
  
  Pyrimidine_synthesis = c("DHODH", "CAD", "UMPS")
)

#check if genes are in the sce (=passed QC steps)
complexi_detected <- signatures$ComplexI[signatures$ComplexI%in% rownames(sce)]
print(length(complexi_detected) == length(signatures$ComplexI))
mito_detected <- signatures$Mito_ETC[signatures$Mito_ETC%in% rownames(sce)]
print(length(mito_detected) == length(signatures$Mito_ETC))
pyr_detected <- signatures$Pyrimidine_synthesis[signatures$Pyrimidine_synthesis%in% rownames(sce)]
print(length(pyr_detected) == length(signatures$Pyrimidine_synthesis))


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

write.csv(signatures_df, "20250425_signatures.csv")





# calculate the mean signature score for each sample

signatures_mean <- signatures_df |>
  filter(gate == "Leukemia") |>
  select("sample_name", "sample_type", "relapse","ComplexI", "Mito_ETC", "Pyrimidine_synthesis") |>
  group_by(sample_name,) |>
  summarise(sample_type = first(sample_type),
            relapse = first(relapse),
            ComplexI = mean(ComplexI, na.rm = TRUE),
            Mito_ETC = mean(Mito_ETC, na.rm = TRUE),
            Pyrimidine_synthesis = mean(Pyrimidine_synthesis, na.rm = TRUE),
            .groups = "drop"
  )

write.csv(signatures_mean, "20250425_mean_signatures.csv")

# statistical tests signature scores ####
# wilcox test with post-hoc comparisons

features <- c("ComplexI", "Mito_ETC", "Pyrimidine_synthesis")

# Empty list to store results
all_results <- list()

for (feat in features) {
  
  sub_data <- signatures_mean %>%
    select(sample_type, all_of(feat)) %>%
    drop_na()
  
  # Get the unique groups
  groups <- unique(sub_data$sample_type)
  
  # Only run if there are at least 2 groups
  if (length(groups) >= 2) {
    
    # Generate all pairwise comparisons
    combs <- combn(groups, 2, simplify = FALSE)
    
    temp_results <- lapply(combs, function(pair) {
      group1 <- pair[1]
      group2 <- pair[2]
      
      # Subset data
      data1 <- sub_data[sub_data$sample_type == group1, feat, drop = TRUE]
      data2 <- sub_data[sub_data$sample_type == group2, feat, drop = TRUE]
      
      # Only run test if both groups have at least 2 samples
      if (length(data1) >= 2 && length(data2) >= 2) {
        test <- wilcox.test(data1, data2)
        
        data.frame(
          feature = feat,
          sample_type_1 = group1,
          sample_type_2 = group2,
          p_value = test$p.value
        )
      } else {
        NULL  # Skip if not enough samples
      }
    })
    
    # Remove NULLs
    temp_results <- bind_rows(temp_results)
    
    # Adjust p-values if there are results
    if (nrow(temp_results) > 0) {
      temp_results$adj_p_value <- p.adjust(temp_results$p_value, method = "holm")
      all_results[[feat]] <- temp_results
    }
    
  } else {
    message(paste("Skipping", feat, ": not enough valid groups"))
  }
}

# Combine all results into a single dataframe
pairwise_results_df <- bind_rows(all_results)
print(pairwise_results_df)

# export results
write.csv(pairwise_results_df, "20250428_pairwise results.csv")


# wilcox test with post-hoc comparisons 
# COMPLETE REMISSION ONLY

features <- c("ComplexI", "Mito_ETC", "Pyrimidine_synthesis")

# Empty list to store results
all_results <- list()

for (feat in features) {
  
  sub_data <- signatures_mean %>%
    filter(relapse == "no") |>
    select(sample_type, all_of(feat)) %>%
    drop_na()
  
  # Get the unique groups
  groups <- unique(sub_data$sample_type)
  
  # Only run if there are at least 2 groups
  if (length(groups) >= 2) {
    
    # Generate all pairwise comparisons
    combs <- combn(groups, 2, simplify = FALSE)
    
    temp_results <- lapply(combs, function(pair) {
      group1 <- pair[1]
      group2 <- pair[2]
      
      # Subset data
      data1 <- sub_data[sub_data$sample_type == group1, feat, drop = TRUE]
      data2 <- sub_data[sub_data$sample_type == group2, feat, drop = TRUE]
      
      # Only run test if both groups have at least 2 samples
      if (length(data1) >= 2 && length(data2) >= 2) {
        test <- wilcox.test(data1, data2)
        
        data.frame(
          feature = feat,
          sample_type_1 = group1,
          sample_type_2 = group2,
          p_value = test$p.value
        )
      } else {
        NULL  # Skip if not enough samples
      }
    })
    
    # Remove NULLs
    temp_results <- bind_rows(temp_results)
    
    # Adjust p-values if there are results
    if (nrow(temp_results) > 0) {
      temp_results$adj_p_value <- p.adjust(temp_results$p_value, method = "holm")
      all_results[[feat]] <- temp_results
    }
    
  } else {
    message(paste("Skipping", feat, ": not enough valid groups"))
  }
}

# Combine all results into a single dataframe
pairwise_results_cr_df <- bind_rows(all_results)
print(pairwise_results_cr_df)

# export results
write.csv(pairwise_results_cr_df, "20250428_pairwise results_CR.csv")

# wilcox test with post-hoc comparisons 
# RELAPSED ONLY

features <- c("ComplexI", "Mito_ETC", "Pyrimidine_synthesis")

# Empty list to store results
all_results <- list()

for (feat in features) {
  
  sub_data <- signatures_mean %>%
    filter(relapse == "yes") |>
    select(sample_type, all_of(feat)) %>%
    drop_na()
  
  # Get the unique groups
  groups <- unique(sub_data$sample_type)
  
  # Only run if there are at least 2 groups
  if (length(groups) >= 2) {
    
    # Generate all pairwise comparisons
    combs <- combn(groups, 2, simplify = FALSE)
    
    temp_results <- lapply(combs, function(pair) {
      group1 <- pair[1]
      group2 <- pair[2]
      
      # Subset data
      data1 <- sub_data[sub_data$sample_type == group1, feat, drop = TRUE]
      data2 <- sub_data[sub_data$sample_type == group2, feat, drop = TRUE]
      
      # Only run test if both groups have at least 2 samples
      if (length(data1) >= 2 && length(data2) >= 2) {
        test <- wilcox.test(data1, data2)
        
        data.frame(
          feature = feat,
          sample_type_1 = group1,
          sample_type_2 = group2,
          p_value = test$p.value
        )
      } else {
        NULL  # Skip if not enough samples
      }
    })
    
    # Remove NULLs
    temp_results <- bind_rows(temp_results)
    
    # Adjust p-values if there are results
    if (nrow(temp_results) > 0) {
      temp_results$adj_p_value <- p.adjust(temp_results$p_value, method = "holm")
      all_results[[feat]] <- temp_results
    }
    
  } else {
    message(paste("Skipping", feat, ": not enough valid groups"))
  }
}

# Combine all results into a single dataframe
pairwise_results_rx_df <- bind_rows(all_results)
print(pairwise_results_rx_df)

# export results
write.csv(pairwise_results_rx_df, "20250428_pairwise results_rx.csv")


# plot signature #### 

# add to the umap df
for (gene in pyr_detected){
  print(gene)
  umap_df[[gene]] <- as.vector(assay(sce, "logcounts")[gene,])
}

# plot for each gene divded by sample type ####

#labels for faceting

labels <- labeller(
  relapse = c(yes = "relapsed", no = "CR")
)

#plots
plots <- purrr::map(pyr_detected, function(gene_name){
  umap_df |>
    filter(sample_type != "Healthy") |>
    filter(gate == "Leukemia") |>
    ggplot(aes(x = UMAP1, y=  UMAP2, color = .data[[gene_name]], alpha = .data[[gene_name]]))+
    geom_point(size = 0.6)+
    scale_color_gradient(low= "grey70", high = "blue")+
    guides(alpha = "none")+
    ggtitle(gene_name)+
    facet_grid(cols = vars(sample_type), rows = vars(relapse), labeller = labels)+
    theme_minimal(base_size = 20)
})

# save plots 
purrr::walk2(
  .x = plots,
  .y = pyr_detected,
  .f = ~ggsave(
    filename = paste0("20250507_", .y, ".png"),
    plot = .x,
    width = 1100,
    height = 680,
    unit = "px",
    dpi = 72,
    scale = 1,
    limitsize = FALSE,
    bg = "white"
  )
)

# UMAP signautres

# complexI

umap_df |>
  
  filter(gate == "Leukemia",
         sample_type != "Healthy") |>
  ggplot( aes(x= UMAP1, y= UMAP2, color =ComplexI))+
  geom_point(size = 0.6)+
  scale_color_gradient(low= "grey70", high = "green")+
  guides(alpha = "none")+
  ggtitle("Complex I signature score")+
  facet_grid(cols = vars(sample_type), rows = vars(relapse), labeller = labels)+
  theme_gray(base_size = 20)

# Mito ETC
umap_df |>
  
  filter(gate == "Leukemia",
         sample_type != "Healthy") |>
  ggplot( aes(x= UMAP1, y= UMAP2, color =Mito_ETC))+
  geom_point(size = 0.6)+
  scale_color_gradient(low= "grey70", high = "orange")+
  guides(alpha = "none")+
  ggtitle("Mitochondria ETC signature score")+
  facet_grid(cols = vars(sample_type), rows = vars(relapse), labeller = labels)+
  theme_gray(base_size = 20)

#pyrimidine synthesis

umap_df |>
  
  filter(gate == "Leukemia",
         sample_type != "Healthy") |>
  ggplot( aes(x= UMAP1, y= UMAP2, color =Pyrimidine_synthesis))+
  geom_point(size = 0.6)+
  scale_color_gradient(low= "grey70", high = "purple")+
  guides(alpha = "none")+
  ggtitle("Pyrimidine synthesis signature score")+
  facet_grid(cols = vars(sample_type), rows = vars(relapse), labeller = labels)+
  theme_gray(base_size = 20)


  
  

# convert to long format for better handling during plotting 

signatures_longer <- signatures_mean |>
  tidyr::pivot_longer(cols = c( "ComplexI","Mito_ETC","Pyrimidine_synthesis"),
               names_to = "signature")


# Violin plot 1 - individual, by sample type

signatures_longer |>
  dplyr::filter(sample_type != "Healthy") |>
  ggplot(aes(x=sample_type, y=value))+
  geom_violin(aes(fill = signature), alpha = 0.7)+
  geom_jitter(aes(fill = signature), size = 1.1, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75))+
  theme_gray(base_size = 24)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  guides(fill = "none", color = "none")+
  ggtitle("Per-sample mean UCell signature score")+
  labs(y = "score")+
  facet_wrap(~signature, ncol = 1, axis.labels = "all_y")


# divide by relapse
p1 <- signatures_longer |>
  dplyr::filter(
    sample_type %in% c("Diagnosis", "MRD timepoints"),
    relapse == "yes",
    ) |>
  
  ggplot(aes(x=sample_type, y=value))+
  geom_violin(aes(fill = signature), alpha = 0.7)+
  geom_jitter(aes(fill = signature),, size = 1.1, position = position_jitterdodge())+
  theme_gray(base_size = 24)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "none")+
  ggtitle("Relapsed")+
  facet_wrap(~signature, axis.labels = "all")

p2 <- signatures_longer |>
  dplyr::filter(
    sample_type != "Healthy",
    relapse == "no") |>
  
  ggplot(aes(x=sample_type, y=value))+
  geom_violin(aes(fill = signature), alpha = 0.7)+
  geom_jitter(aes(fill = signature),, size = 1.1, position = position_jitterdodge())+
  theme_gray(base_size = 26)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  ggtitle("Complete remission")+
  facet_wrap(~signature, axis.labels = "all")
  

grid <- wrap_plots(p1,p2, 
                   axes = "collect", 
                   axis_titles = "collect",
                   ncol = 2,
                   tag_level = "new", 
                   guides = "collect")
grid


ggsave("20250425_no relapse_signature score by outcome.png",
       plot = grid, width = 1100,height = 800, unit = "px",
       scale = 4)


# another way


signatures_longer |>
  dplyr::filter(
    sample_type %in% c("Diagnosis", "MRD timepoints"),
    relapse %in% c("yes", "no")
  ) |>

  ggplot(aes(x = sample_type, y = value)) +
    geom_violin(aes(fill = signature), alpha = 0.7) +
    geom_jitter(aes(fill = signature), size = 1.1, position = position_jitterdodge()) +
    facet_grid(relapse~signature) +
    theme_gray(base_size = 24) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    labs(title = "Signature score by relapse outcome", x = "Sample Type", y = "Signature score")
# heatmaps single genes  ###

# the part below could DEFINETLY be automated with a function + using purrr
# but for now let's do mancave style

# one heatmpa per gene of each signature 

# Complex I signature heatmap #####

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
  complexi_df[[gene]] <- as.vector(assay(sce[, (sce$gate == "Leukemia")&(sce$sample_type != "Healthy")], "logcounts")[gene,]) # !! RAW COUNTS!!
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
      summarise(
        mean = mean(value),
        .groups = "drop") |>
  # calculate z score
    group_by(genes) |>
  mutate(global_mean = mean(mean),
         global_sd = sd(mean))|>
  mutate(
    z_score = (mean -global_mean)/global_sd
  )

# convert to matrix
complexi_mat <- complexi_longer |>
  select(c("sample_name", "genes", "z_score")) |>
  pivot_wider(names_from = sample_name, values_from = z_score) |>
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

# Mito_ETC heatmap ####



mito_df <-  data.frame(
  individual = signatures_df_leukemia$individual,
  timepoint = signatures_df_leukemia$timepoint,
  sample_type = signatures_df_leukemia$sample_type,
  sample_name = signatures_df_leukemia$sample_name,
  relapse = signatures_df_leukemia$relapse
)

# individual expression of genes in the signatures
for (gene in signatures$Mito_ETC){
  print(gene)
  mito_df[[gene]] <- as.vector(assay(sce[, (sce$gate == "Leukemia")&(sce$sample_type != "Healthy")], "logcounts")[gene,]) # !! RAW COUNTS!!
}

# convert to longer
mito_longer <- mito_df |>
  pivot_longer(cols = 6:9,
               names_to = "genes",
               values_to = 'value')

# annotation matrix
annot <- mito_longer |>
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
mito_longer <- mito_longer |>
  group_by(genes, sample_name, individual) |>
  summarise(
    mean = mean(value),
    .groups = "drop") |>
  # calculate z score
  group_by(genes) |>
  mutate(global_mean = mean(mean),
         global_sd = sd(mean))|>
  mutate(
    z_score = (mean -global_mean)/global_sd
  )

# convert to matrix
mito_mat <- mito_longer |>
  select(c("sample_name", "genes", "z_score")) |>
  pivot_wider(names_from = sample_name, values_from = z_score) |>
  tibble::column_to_rownames("genes") |>
  as.matrix()

# annotation colors
annot_colors <- list(
  relapse = c("yes"="red","no"="blue"),
  sample_type = c("Diagnosis" = "lightblue", "MRD timepoints" = "purple", "Relapse" = "rosybrown")
)

#heatmap
pheatmap(mito_mat,  # Only numeric columns for the heatmap
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
         main = "Mitochondrial genes - ETC",  # Title of the heatmap
         fontsize_row = 8,  # Decrease font size of row names (patient IDs) to make more space
         angle_row = 45,  # Rotate row names to avoid squishing
         cellheight = 8,  # Adjust the height of each cell for better spacing of row names
         
)

# Pyrimidine synthesis genes


pyr_df <-  data.frame(
  individual = signatures_df_leukemia$individual,
  timepoint = signatures_df_leukemia$timepoint,
  sample_type = signatures_df_leukemia$sample_type,
  sample_name = signatures_df_leukemia$sample_name,
  relapse = signatures_df_leukemia$relapse
)

# individual expression of genes in the signatures
for (gene in signatures$Pyrimidine_synthesis){
  print(gene)
  pyr_df[[gene]] <- as.vector(assay(sce[, (sce$gate == "Leukemia")&(sce$sample_type != "Healthy")], "logcounts")[gene,]) # !! RAW COUNTS!!
}

# convert to longer
pyr_longer <- pyr_df |>
  pivot_longer(cols = 6:8,
               names_to = "genes",
               values_to = 'value')

# annotation matrix
annot <- pyr_longer |>
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
pyr_longer <- pyr_longer |>
  group_by(genes, sample_name, individual) |>
  summarise(
    mean = mean(value),
    .groups = "drop") |>
  # calculate z score
  group_by(genes) |>
  mutate(global_mean = mean(mean),
         global_sd = sd(mean))|>
  mutate(
    z_score = (mean -global_mean)/global_sd
  )

# convert to matrix
pyr_mat <- pyr_longer |>
  select(c("sample_name", "genes", "z_score")) |>
  pivot_wider(names_from = sample_name, values_from = z_score) |>
  tibble::column_to_rownames("genes") |>
  as.matrix()

# annotation colors
annot_colors <- list(
  relapse = c("yes"="red","no"="blue"),
  sample_type = c("Diagnosis" = "lightblue", "MRD timepoints" = "purple", "Relapse" = "rosybrown")
)

#heatmap
pheatmap(pyr_mat,  # Only numeric columns for the heatmap
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
         main = "De-novo pyrimidine synthesis genes",  # Title of the heatmap
         fontsize_row = 8,  # Decrease font size of row names (patient IDs) to make more space
         angle_row = 45,  # Rotate row names to avoid squishing
         cellheight = 8,  # Adjust the height of each cell for better spacing of row names
)


# Violin pltos single-cell expression of genes COMPLEX I ####


# one grid per signature
plots <- purrr::map(colnames(complexi_df[6:36]), function(gene) {
  complexi_df |>
    filter(sample_type != "Healthy") |>
    ggplot(aes(x= sample_type, y = .data[[gene]], fill = sample_type))+
    geom_violin(aes(fill= sample_type))+
    ggtitle(gene)+
    #geom_jitter(position = position_jitter(), aes(fill = sample_type,color = sample_type), alpha = 0.3)+
    theme_gray(base_size = 18)
})

grid<-wrap_plots(plots, ncol = 7, axes = "collect_x", axis_titles = "collect", guides = "collect")
ggsave("20250425_COmplexI violin grid.png",
       plot = grid, width = 2100,height = 1300, unit = "px",
       scale = 4)

grid
# Violin pltos single-cell expression of genes MITO-ETC ####


# one grid per signature
plots <- purrr::map(colnames(mito_df[6:9]), function(gene) {
  mito_df |>
    filter(sample_type != "Healthy") |>
    ggplot(aes(x= sample_type, y = .data[[gene]], fill = sample_type))+
    geom_violin(aes(fill= sample_type))+
    ggtitle(gene)+
    labs(y = "gene expression")+
    #geom_jitter(position = position_jitter(), aes(fill = sample_type,color = sample_type), alpha = 0.3)+
    theme_gray(base_size = 24)
})

grid<-wrap_plots(plots, ncol = 2, axes = "collect_x", axis_titles = "collect", guides = "collect")
ggsave("20250425_mito etc violin grid.png",
       plot = grid, width = 1600,height =1226, unit = "px",
       scale = 5)


# Violin pltos single-cell expression of genesPYRIMIDINE SYNTHESIS ####


# one grid per signature
plots <- purrr::map(colnames(pyr_df[6:8]), function(gene) {
  pyr_df |>
    filter(sample_type != "Healthy") |>
    ggplot(aes(x= sample_type, y = .data[[gene]], fill = sample_type))+
    geom_violin(aes(fill= sample_type))+
    ggtitle(gene)+
    labs(y = "gene expression")+
    #geom_jitter(position = position_jitter(), aes(fill = sample_type,color = sample_type), alpha = 0.3)+
    theme_gray(base_size = 26)
})

grid<-wrap_plots(plots, ncol = 3, axes = "collect_x", axis_titles = "collect", guides = "collect")
ggsave("20250425_pyr_synthesis violin grid.png",
       plot = grid, width = 1850,height =500, unit = "px",
       scale = 5)


