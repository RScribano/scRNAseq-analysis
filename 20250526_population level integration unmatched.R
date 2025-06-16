# packages ####
#tidyverse
library(dplyr)
library(tidyr)
library(purrr)
library(tidyverse)

# stats
library(vegan)
library(PMCMRplus)
library(FSA)

#single cell analysis
library(scran)
library(scater) 

#plotting
library(ggplot2)
library(pheatmap)

# wrangling ####
original_tags <- rownames(altExp(sce, "ADT"))

new_tags <- original_tags |>
  sub("-.*", "", x = _) |>
  sub(":.*", "", x= _)

rowData(altExp(sce, "ADT"))$antigen <- new_tags

rownames(altExp(sce, "ADT")) <- rowData(altExp(sce, "ADT"))$antigen

# extract adt expression as dataframe
adt_df <- as.data.frame(t(logcounts(altExp(sce, "ADT"))))

# extract genes for dtd, CD34
adt_df['tdt'] <- assay(sce, "logcounts")["DNTT",]
adt_df['cd34'] <- assay(sce, "logcounts")["CD34",]

# add key for sce metadata 
adt_df$Sample_Name <- coldata_df$Sample_Name[]

#export .csv
write.csv(adt_df, "20250526_adt_exp_supplemented.csv")

# import CITEseq data mapped to heatlhy BM cytof ####
adt_mapped_path <- file.path('..', 'MRD_project', 'Population level integration', 'CITEseq mapped', "20250526_citeseq_mapped.csv")
adt_mapped_df <- read.csv(adt_mapped_path)

# check if they are aligned
all.equal(adt_mapped_df$cd34, assay(sce, "logcounts")["CD34", ])

# add to umap_df
umap_df$celltype_cytof <- adt_mapped_df$Pop_names

# Comparison with Azimuth ####
table <- as.data.frame( table(umap_df$Cell_Type_Experimental, umap_df$celltype_cytof))
colnames(table) <- c("Celltype.Azimuth", "Predicted.cytof", "N.cells")

# create contingency table
table <- table |>
  dplyr:: group_by(Celltype.Azimuth) |>
  dplyr:: mutate( total_cells = sum(N.cells),
                  percent = 100*N.cells/total_cells)

table$percent <- round(table$percent, digits = 2)

# plot contingency table
ggplot(table, aes(x = Predicted.cytof , y = Celltype.Azimuth, fill = percent)) +
  geom_tile() +
  geom_text(aes(label = percent), color = "black") +
  scale_fill_gradient(low = "grey80", high = "steelblue") +
  ggtitle("Normalized by cell-type (azimuth)")+
  theme_minimal(base_size = 22)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Comparison with manual gate ####
table <- as.data.frame( table(umap_df$gate, umap_df$celltype_cytof))
colnames(table) <- c("Manual_gate.blasts", "Predicted.cytof", "N.cells")

# create contingency table
table <- table |>
  dplyr:: group_by(Manual_gate.blasts) |>
  dplyr:: mutate( total_cells = sum(N.cells),
                  percent = 100*N.cells/total_cells)

table$percent <- round(table$percent, digits = 2)

# plot contingency table
ggplot(table, aes(x = Predicted.cytof , y = Manual_gate.blasts, fill = percent)) +
  geom_tile() +
  geom_text(aes(label = percent), color = "black") +
  scale_fill_gradient(low = "grey80", high = "steelblue") +
  ggtitle("Normalized by manual gate")+
  theme_minimal(base_size = 22)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  

# plot UMAP with predicted cell type ####

celltype_palette <- c(
  "Mature_non-B"     = "#E41A1C",  # red
  "Pre-BI"           = "#377EB8",  # blue
  "Pre-ProB"         = "#4DAF4A",  # green
  "Mature-BII"       = "#984EA3",  # purple
  "Immature-BI"      = "#FF7F00",  # orange
  "Mature-BI"        = "#A65628",  # brown
  "Immature-BII"     = "#F781BF",  # pink
  "Pre-BII"          = "#999999",  # grey
  "late_progenitors" = "#66C2A5",  # teal
  "early_non-B"      = "#FFD92F"   # yellow
)

# all cells
ggplot(umap_df, aes(x = UMAP1, y = UMAP2))+
  geom_point(aes(color = celltype_cytof), size = 0.7, alpha = 0.3)+
  theme_minimal(base_size = 22) +
  scale_color_manual(values = celltype_palette)+
  guides(color = guide_legend(override.aes = list(size=6)))+
  facet_wrap(~sample_type)
  

# healthy only
umap_df |>
  filter(sample_type == "Healthy") |>
  ggplot(aes(x = UMAP1, y = UMAP2))+
  geom_point(aes(color = celltype_cytof), size = 0.7, alpha = 0.3)+
  theme_minimal(base_size = 22) +
  scale_color_manual(values = celltype_palette)+
  guides(color = guide_legend(override.aes = list(size=6)))


# barplot abundance of predicted populations in different sampletypes ####


# add sample type information
adt_mapped_df$sample_type <- umap_df$sample_type

# crete abundance df
abundance_df <- adt_mapped_df |>
  group_by(sample_type, Pop_names) |>
  summarise(count = n(), .groups = "drop") |>
  group_by(sample_type) |>
  mutate(freq = count/sum(count))

ggplot(abundance_df, aes(x = Pop_names, y = freq)) +
  geom_bar(stat = "identity", position = position_dodge(), aes(fill = Pop_names)) +
  scale_fill_manual(values = celltype_palette) +
  facet_wrap(~sample_type) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Heatmaps expression of CITE-seq markers//genes used in developmental classifier in predicted populations ####

#columns to exclude 
exclude <- c(
  "mah_Pre.BII",
  "mah_Mature.BII",
  "mah_Immature.BI",
  "mah_Mature_non.B",
  "mah_early_non.B",
  "mah_Immature.BII",
  "mah_Pro.BII",
  "mah_late_progenitors",
  "mah_CLP",
  "mah_Mature.BI",
  "mah_Pre.ProB",
  "mah_HSC",
  "mah_Pre.BI",
  "mah_Pro.BI",
  "Sample_Name",
  "plate",
  "sample_type"
)

# cretae df with scaled expression of markers
markers_df <- adt_mapped_df |>
  select(-any_of(exclude)) |>
  pivot_longer(cols = -c("Pop_names"),
               values_to = "expression",
               names_to = "marker") |>
  # get populations means
  group_by(Pop_names, marker) |>
  summarise(pop_mean = mean(expression),
            .groups = "drop") |>
  # get global mean for each marker
  left_join(
    adt_mapped_df|>
      select(-any_of(exclude)) |>
      pivot_longer(cols = -c("Pop_names"),
                   values_to = "expression",
                   names_to = "marker") |>
      group_by(marker) |>
      summarise(global_mean = mean(expression),
                global_sd = sd(expression)),
    by = "marker"
  ) |>
  
  # calculate zscore
  mutate(z_score = (pop_mean - global_mean)/global_sd)


# turn to matrix
markers_mtx <- markers_df |>
  select(-c('global_mean', "global_sd", "pop_mean")) |>
  pivot_wider(names_from = "marker",
              id_cols = "Pop_names",
              values_from = "z_score"
              ) |>
  column_to_rownames("Pop_names") |>
  as.matrix()


# create heatmap
pheatmap(mat = markers_mtx,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete")


## violinplot markers by predicted population #### 

adt_mapped_df |>
  select(-any_of(exclude)) |>
  pivot_longer(cols = -c("Pop_names"),
               names_to = "marker",
               values_to = "expression") |>
  ggplot(aes(x = Pop_names, y = expression)) +
  geom_violin(aes(fill = Pop_names), position = position_dodge()) +
  scale_fill_manual(values = celltype_palette) +
  facet_wrap(~marker, ncol = 2, scales = "free_y") +
  theme_grey(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
