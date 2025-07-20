# packages ####
#tidyverse
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(zellkonverter)
library(BoneMarrowMap)


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
library(Seurat)

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


#

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




# Assign populations from reference (bonemarrowMap) ####

# load seurat object (BMM is deisgned to work with seurat objects..)

seurat_path <- file.path("merged_carts.RDS")
sce_seurat <- readRDS(seurat_path)

# obtain patient
meta <- as.data.frame(sce_seurat@meta.data)

meta <- meta |>
  separate(Sample_Name, into = c("individual","timepoint"), sep = "_")

sce_seurat@meta.data <- sce_seurat@meta.data |>
  separate(Sample_Name, into = c("individual","timepoint"), sep = "_")
  

# Load Symphony reference
ref_path <- file.path('BoneMarrowMap_SymphonyReference.rds')
ref <- readRDS(ref_path)

ref$save_uwot_path <- file.path('BoneMarrowMap_uwot_model.uwot')

# map to reference 

sce_seurat <- map_Query(
  exp_query = sce_seurat[['RNA']]@counts, 
  metadata_query = sce_seurat@meta.data,
  ref_obj = ref,
  vars = 'individual',
)

# caclualte mapping error
sce_seurat <- sce_seurat |> calculate_MappingError(query = _, reference = ref, 
                                                   MAD_threshold = 2.5,
                                                   threshold_by_donor = TRUE,
                                                   donor_key = "individual")

# plot mapping error distirbutions
sce_seurat@meta.data %>% 
  ggplot(aes(x = mapping_error_score, fill = mapping_error_QC)) + 
  geom_histogram(bins = 200) + facet_wrap(.~get("individual"))


# Cell type assignemnt

sce_seurat <- predict_CellTypes(
  query_obj = sce_seurat,
  ref_obj = ref,
  initial_label = 'initial_CellType_BoneMarrowMap', #before QC filtering
  final_label = 'predicted_CellType_BoneMarrowMap',
  
)



# transfer to SCE object
DimPlot(subset(sce_seurat, mapping_error_QC == 'Pass'), reduction = 'umap_projected', group.by = c('predicted_CellType_BoneMarrowMap'), 
        raster=FALSE, label=TRUE, label.size = 4)

# Extract UMAP coordinates from Seurat
umap_coords <- Embeddings(sce_seurat, reduction = "umap_projected")

# Confirm it's a matrix (should be, but good to check)
umap_coords <- as.matrix(umap_coords)

# Assign it to the SCE object's reducedDims
reducedDim(sce, "UMAP_projected") <- umap_coords

# attach cell type labels 
sce$predicted_BMM<-sce_seurat$predicted_CellType_BoneMarrowMap
sce$predicted_BMM_proba <- sce_seurat$predicted_CellType_BoneMarrowMap_prob

sce$predicted_BMM[is.na(sce$predicted_BMM)] <- "Unassigned"
sce$predicted_BMM_proba[is.na(sce$predicted_BMM_proba)] <- 0


saveRDS(sce, "20250716_sce_projected.rds")

sce <- readRDS("20250716_sce_projected.rds")

# plot UMAP healthy with assinged celltype

sce |>
   ggplot(aes(x = umap_1, y =  umap_2)) +
  geom_point(aes(color = predicted_BMM), size = 0.6) +
  facet_wrap(~individual) +
  theme_bw(base_size = 22) +
  theme(
    legend.text  = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

sce |>
  join_features("CD4") |>
  ggplot(aes(x = umap_1, y =  umap_2)) +
  geom_point(aes(color = .abundance_logcounts ), size = 0.6, alpha = 0.6) +
  facet_wrap(~individual) +
  theme_bw(base_size = 22) +
  scale_color_viridis_c(trans = "log10")+
  theme(
    legend.text  = element_text(size = 14),
    legend.title = element_text(size = 16)
  )



# downstream checks populations blasts #####

# used to see the proportion of populations within each gate
table_1 <- as.data.frame(table(sce[,sce$sample_type != "Healthy"]$gate, 
                               sce[,sce$sample_type != "Healthy"]$predicted_BMM,
                               sce[,sce$sample_type != "Healthy"]$individual))

proportions_pops <- table_1 |>
  group_by(Var1, Var3) |>
  mutate(tot_pop = sum(Freq)) |>
  ungroup() |>
  group_by(Var1,Var2, Var3) |>
  mutate(rel_freq = Freq/tot_pop) |>
  ungroup() |>
  rename(
    Population_BMM = Var2,
    Gate = Var1,
    individual = Var3
  )

# used to see the proportion of different populations within the blasts gate 
# across timepoints and relapse status
table_2 <- as.data.frame(table(sce[,sce$sample_type != "Healthy"]$gate, 
                               sce[,sce$sample_type != "Healthy"]$predicted_BMM,
                               sce[,sce$sample_type != "Healthy"]$sample_type,
                               sce[,sce$sample_type != "Healthy"]$relapse))

proportions_pops_2 <- table_2 |>
  filter(Var1 == "Leukemia",
         Var3 != "Healthy") |>
  group_by(Var3, Var4) |>
  mutate(tot = sum(Freq)) |>
  ungroup() |>
  group_by(Var2, Var3, Var4) |>
  mutate(rel_freq = Freq/tot) |>
  ungroup() |>
  rename(
    gate = Var1,
    Population_BMM = Var2,
    sample_type = Var3,
    relapse = Var4
  )
  
proportions_pops_2 |>
  ggplot(aes(x = Population_BMM, y = rel_freq))+
  geom_bar(stat = "identity", aes(fill = relapse)) +
  facet_wrap(~ sample_type) +
  theme(
    axis.text.x  = element_text(size = 9,
                                angle = 90,
                                hjust = 1)
  )
  
  
  
# plot proportions populations in gates
proportions_pops |>
  ggplot(aes(x = Gate, y = Freq )) +
  geom_bar(position = position_stack(), aes(fill = Population_BMM),
           stat = "identity") +
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12)
  )+
  facet_wrap(~individual, scales = 'free') 


# propotion of populations by relapse and sample type within blasts
proportions_pops_2 |>
  ggplot(aes(x = Population_BMM, y = rel_freq))+
  geom_bar(stat = "identity", aes(fill = relapse)) +
  facet_wrap(~ sample_type) +
  theme(
    axis.text.x  = element_text(size = 9,
                                angle = 90,
                                hjust = 1)
  )

# subset to developmental stages B cells only

b_lymphopoiesis_cells <- c(
  "HSC",
  "LMPP",
  "CLP",
  "MLP",
  "MLP-II",
  "Pre-ProB",
  "Pro-B VDJ",
  "Pro-B Cycling",
  "Large Pre-B",
  "Small Pre-B",
  "Immature B",
  "Mature B",
  "Plasma Cell"
)

# same as above but only B cell development
proportions_pops_b <- table_2 |>
  filter(Var1 == "Leukemia" &
         Var3 != "Healthy" &
         Var2 %in% b_lymphopoiesis_cells) |>
  group_by(Var3, Var4) |>
  mutate(tot = sum(Freq)) |>
  ungroup() |>
  group_by(Var2, Var3, Var4) |>
  mutate(rel_freq = Freq/tot) |>
  ungroup() |>
  rename(
    gate = Var1,
    Population_BMM = Var2,
    sample_type = Var3,
    relapse = Var4
  ) |>
  mutate(Population_BMM = factor(Population_BMM, levels = b_lymphopoiesis_cells, ordered = TRUE))

# visualization: alternative 1 
proportions_pops_b |>
  ggplot(aes(x = Population_BMM, y = rel_freq))+
  geom_bar(stat = "identity", aes(fill = relapse), position = position_dodge()) +
  facet_wrap(~ sample_type) +
  theme(
    axis.text.x  = element_text(size = 9,
                                angle = 90,
                                hjust = 1)
  )

# visualziation alternative 2
proportions_pops_b |>
  ggplot(aes(x = Population_BMM, y = rel_freq))+
  geom_bar(stat = "identity", aes(fill = sample_type), position = position_dodge()) +
  facet_wrap(~ relapse) +
  theme(
    axis.text.x  = element_text(size = 9,
                                angle = 90,
                                hjust = 1)
  )

# proportion of B developmental populations in blasts of relapsed vs. nonrelapsed 
# at different timepoints

proportions_pops |>
  ggplot(aes(x = Gate, y = Freq )) +
  geom_bar(position = position_stack(), aes(fill = Population_BMM),
           stat = "identity") +
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12)
  )+
  facet_wrap(~individual, scales = 'free') 


  

  
# Plots expression of ADTs ####

# these plots were made t ocheck the expression of ADTs in the populations
# predicted using BMM and to check if the gates (leukemia) were re-imported correctly



# extract UMAP projected coordinates
umap_projected <- reducedDim(sce, "UMAP_projected") |>
  as.data.frame()

# make tags readable
original_tags <- rownames(altExp(sce, "ADT"))

new_tags <- original_tags |>
  sub("-.*", "", x = _) |>
  sub(":.*", "", x= _)

tag_map <- setNames(new_tags, original_tags)



## dataframe with ADT expression
adt_df <- assay(altExp(sce, "ADT"), "logcounts") |>
  t() |>
  as.data.frame() |>
  dplyr::mutate(CellType_BMM  = sce$predicted_BMM,
                umap_1 = umap_projected$umap_1,
                umap_2 = umap_projected$umap_2,
                UMAP_1 = as.data.frame(reducedDim(sce, "UMAP"))$"UMAP1",
                UMAP_2 = as.data.frame(reducedDim(sce, "UMAP"))$"UMAP2",
                gate = sce$gate,
                sample_type = sce$sample_type,
                individual = sce$individual,
                proba = sce$predicted_BMM_proba) |>
  pivot_longer(cols = -c(CellType_BMM, umap_1, umap_2, gate, 
                         sample_type, individual, UMAP_1, UMAP_2, proba), 
               values_to = "expression", names_to = "ADT_marker") |>
  mutate(ADT_marker = recode(ADT_marker, !!!tag_map))

# Violinplot predicted popultions
adt_df |>
  ggplot(aes(x = ADT_marker, y = expression)) +
  geom_violin(aes(fill = ADT_marker), position = position_dodge()) +
  scale_y_log10()+
  facet_wrap(~CellType_BMM) +
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


adt_df |>
  filter(sample_type != "Healthy",
         gate == "Leukemia",
         ADT_marker %in% c("IgM")) |>
  mutate(expression = as.numeric(expression)) |>
  ggplot(aes(x = umap_1, y = umap_2)) +
  geom_point(aes(color = expression), size = 0.8, alpha  = 0.6) +
  scale_color_viridis_c(trans = "log10")+
  theme_bw() +
  ggtitle("IgM")+
  theme(
    plot.title = element_text(size = 28, face = "bold")
  )
  

# color BLASTS by predoicted probability BMM mapping
adt_df |>
  filter(sample_type != "Healthy",
         gate == "Leukemia") |>
  ggplot(aes(x = umap_1, y = umap_2)) +
  geom_point(aes(color = proba), size = 0.8, alpha  = 0.6) +
  scale_color_viridis_c()+
  theme_bw()+
  facet_wrap(~individual)

# only cells with >= 75% proba
adt_df |>
  filter(sample_type != "Healthy",
         gate == "Leukemia",
         proba >= 0.75) |>
  ggplot(aes(x = umap_1, y = umap_2)) +
  geom_point(aes(color = CellType_BMM), size = 0.8, alpha  = 0.6) +
  theme_bw()+
  facet_wrap(~individual)

# distribution proba BMM for each population
pops <- unique(sce$predicted_BMM)
pops_1 <- pops[1:11]
pops_2 <-pops[12:23]
pops_3 <- pops[24:35]
pops_4 <- pops[36:47]
pops_5 <- pops[48:55]

pops_list <- list(pops_1, pops_2, pops_3, pops_4, pops_5)

for (i in seq_along(pops_list)) {
  current_pops <- pops_list[[i]]
  
  # Filter and plot
  p <- adt_df |>
    filter(sample_type != "Healthy",
           ADT_marker %in% c("IgM"),
           CellType_BMM %in% current_pops) |>
    ggplot(aes(x = CellType_BMM, y = proba)) +
    geom_violin(position = position_dodge(), aes(fill = gate)) +
    theme(axis.text.x = element_text(size = 7, angle = 90, face = "italic")) +
    ggtitle(paste("Violin plot for pops group", i))
  
  # Save the plot
  ggsave(filename = paste0("violin_plot_group_", i, ".png"), plot = p, width = 6, height = 4, dpi = 300)
}


# violin plots ADT expression ####
adt_df |>
  filter(sample_type != "Healthy",
         
         ADT_marker %in% c("CD10", "CD45", "CD19")) |>
  ggplot(aes(x = ADT_marker, y = expression)) +
  geom_violin(position = position_dodge(), aes(fill = gate)) +
  scale_y_log10()+
  theme_bw(base_size = 22)+
  facet_wrap(~sample_type)

# by individual
adt_df |>
  filter(sample_type != "Healthy",
         ADT_marker %in% c("CD10", "CD45", "CD19", "IgM")) |>
  ggplot(aes(x = ADT_marker, y = expression)) +
  geom_violin(position = position_dodge(), aes(fill = gate)) +
  scale_y_log10()+
  theme_bw(base_size = 22)+
  facet_wrap(~individual)


# UMAPs (on our UMAP coordinates, not projected)

sub_adt <- adt_df |>
  filter(ADT_marker %in% c("CD10", "CD19", "CD45", "IgM"),
                  sample_type != "Healthy")

sub_adt |>
  ggplot(aes( x= UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = expression), size = 0.8, alpha = 0.6) +
  facet_wrap(~ADT_marker) +
  theme_minimal() +
  scale_color_viridis_c(trans = "log10")

ggsave("plot.png")


# UMAPs ADT expression predicted populations BMM ####
# Loop over each marker name
for (marker in unique(adt_df$ADT_marker)) {
  
  # Filter data for the current marker
  plot_data <- adt_df %>% filter(ADT_marker == marker)
  
  # Build the plot
  p <- ggplot(plot_data, aes(x = umap_1, y = umap_2)) +
    geom_point(aes(color = expression), alpha = 0.5, size = 0.7) +
    scale_color_viridis_c() +
    ggtitle(paste("UMAP for", marker)) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid = element_line(color = "grey90")
    ) +
    scale_color_viridis_c(trans = "log10")
  
  
  # Save to file
  ggsave(
    filename = paste0("ADT_UMAP_", marker, ".png"),
    plot = p,
    width = 5,
    height = 4,
    dpi = 300,
  )
}


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

sce_blasts <- sce[,sce$gate == "Leukemia"]
sce_blasts <- sce_blasts[,sce_blasts$sample_type != "Healthy"]

# create column for timepoint_individual 

sce_blasts$timepoint_patient <- paste(sce_blasts$sample_type, 
                                      sce_blasts$relapse, sce_blasts$individual)

# aggregate pseudobulk counts - relapsed patients
pb <- scuttle::aggregateAcrossCells(sce_blasts[, sce_blasts$relapse == "yes"],
                                    ids = sce_blasts[, sce_blasts$relapse == "yes"]$timepoint_patient )

counts_relapsed <- counts(pb)


# set up differential expression analysis

#1 = diagnosis
#2= MRD timepoints
#3 = Relapse timepoint

group <- factor(c(1,1,1,1,1,2,2,2,2,3,3,3,3))
patient <- factor(c('MRD_01', 'MRD_10', 'MRD_11', 'MRD_27',
                    'MRD_45', 'MRD_01', 'MRD_10', 'MRD_11', 
                    'MRD_27', 'MRD_01', 'MRD_10', 'MRD_11',
                    'MRD_27'))


# create dgelist 
y <- DGEList(counts = counts_relapsed,
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
topde_mrd <- as.data.frame(topTags(lrt_mrd,adjust.method = 'holm',n = Inf))
topde_rx <- as.data.frame(topTags(lrt_relapse,adjust.method = 'holm',n = Inf))
topde_mrd_rx <- as.data.frame(topTags(lrt_mrd_relapse,adjust.method = 'holm',n = Inf))

# volcano plots

EnhancedVolcano(toptable = topde_mrd,
                lab = rownames(topde_mrd),
                x= "logFC",
                y = "PValue",
                pCutoff = 0.25,
                pCutoffCol = "FWER",
                FCcutoff = 1,
                title = 'Relapsed - Dx vs. MRD')


EnhancedVolcano(toptable = topde_rx,
                lab = rownames(topde_rx),
                x= "logFC",
                y = "PValue",
                pCutoff = 0.25,
                pCutoffCol = "FWER",
                FCcutoff = 1,
                title = 'Relapsed - Dx vs. Rx')


EnhancedVolcano(toptable = topde_mrd_rx,
                lab = rownames(topde_mrd_rx),
                x= "logFC",
                y = "PValue",
                pCutoff = 0.25,
                pCutoffCol = "FWER",
                FCcutoff = 1,
                title = 'Relapsed - MRD vs. Rx')


# now i want to plot these genes at single cell level in relapsed patients ####
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



# DEG timepoints (Dx vs MRD) for NOT relapsed ####

# DEG leukemia vs. unfiltered ####

sce$patient_gate <- paste(sce$individual, sce$gate)

pb <- scuttle::aggregateAcrossCells(sce[, sce$sample_type != 'Healthy'],
                                    ids = sce[, sce$sample_type != 'Healthy']$patient_gate )

counts <- counts(pb)




# set up differential expression analysis

#1 = leukemia
#2= unfiltered

group <- c(
  1, 2,  # MRD_01
  1, 2,  # MRD_09
  1, 2,  # MRD_10
  1, 2,  # MRD_11
  1, 2,  # MRD_20
  1, 2,  # MRD_25
  1, 2,  # MRD_27
  1, 2,  # MRD_30
  1, 2,  # MRD_43
  1, 2,  # MRD_45
  1, 2   # MRD_53
)

patient <- c(
  "MRD_01", "MRD_01",
  "MRD_09", "MRD_09",
  "MRD_10", "MRD_10",
  "MRD_11", "MRD_11",
  "MRD_20", "MRD_20",
  "MRD_25", "MRD_25",
  "MRD_27", "MRD_27",
  "MRD_30", "MRD_30",
  "MRD_43", "MRD_43",
  "MRD_45", "MRD_45",
  "MRD_53", "MRD_53"
)



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
lrt <- glmLRT(fit, coef = "group")

# contrasts
topde <- as.data.frame(topTags(lrt,adjust.method = 'holm',n = Inf))

# volcano plots

EnhancedVolcano(toptable = topde,
                lab = rownames(topde),
                x= "logFC",
                y = "PValue",
                pCutoff = 0.25,
                pCutoffCol = "FWER",
                FCcutoff = 1,
                title = 'CR - Dx vs. MRD')

# Log-CPM values for PCA
logCPM <- cpm(y, log=TRUE)

# Run PCA
pca <- prcomp(t(logCPM), scale. = TRUE)

# Create metadata
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  PC3= pca$x[,3] ,
  PC4 = pca$x[,4],
  group = group
)

# Plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA of Pseudobulk Samples", x = "PC1", y = "PC2") +
  theme(text = element_text(size = 14))

ggplot(pca_df, aes(x = PC3, y = PC4, color = group)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA of Pseudobulk Samples", x = "PC3", y = "PC4") +
  theme(text = element_text(size = 14))


# pathway analysis

topde_filtered <- topde %>%
  filter(FWER < 0.05, abs(logFC) > 1) %>%
  rownames_to_column("SYMBOL")

# Step 2: Map SYMBOL to ENTREZID
topde_entrez <- bitr(topde_filtered$SYMBOL,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)

# Step 3: Join back safely
topde_merged <- left_join(topde_filtered, topde_entrez, by = "SYMBOL")


ego <- enrichGO(gene         = topde_merged$ENTREZID,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",          # "BP"=Biological Process, "MF", "CC"
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2,
                readable     = TRUE)

clusterProfiler::dotplot(ego)


# create new identifier mixing relapse status and timepoint


pb <- scuttle::aggregateAcrossCells(sce_blasts,
                                    ids = sce_blasts[, sce_blasts$relapse == "no"]$timepoint_patient )

counts_cr <- counts(pb)


# set up differential expression analysis

#1 = diagnosis
#2= MRD timepoints

group <- factor(c(1,1,1,1,1,1,2,2,2,2,2))
patient <- factor(c("MRD_09", "MRD_20", "MRD_25", "MRD30", 
                    "MRD_43", "MRD_53","MRD_09", "MRD_20", "MRD_25", 
                    "MRD_43", "MRD_53" ))


# create dgelist 
y <- DGEList(counts = counts_cr,
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

# contrasts
topde_mrd <- as.data.frame(topTags(lrt_mrd,adjust.method = 'holm',n = Inf))


# volcano plots

EnhancedVolcano(toptable = topde_mrd,
                lab = rownames(topde_mrd),
                x= "logFC",
                y = "PValue",
                pCutoff = 0.25,
                pCutoffCol = "FWER",
                FCcutoff = 1,
                title = 'CR - Dx vs. MRD')




# quick patwhay analysis

topde_mrd <- topde_mrd |>
  rownames_to_column("gene") |>
  filter(FWER <= 0.25,
         abs(logFC) > 1)


genes_dx <-  bitr(topde_mrd$gene, fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db)

topde_mrd <- topde_mrd |>
  rename(SYMBOL = gene) |>
  left_join(genes_dx, by = "SYMBOL")


ego <- enrichGO(gene         = topde_mrd$ENTREZID,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",          # "BP"=Biological Process, "MF", "CC"
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2,
                readable     = TRUE)

clusterProfiler::dotplot(ego)




# DEG Dx relapse vs Not relapsed ####
sce_blasts$relapse_timepoint_patient <- paste(sce_blasts$sample_type, 
                                              sce_blasts$relapse,
                                              sce_blasts$individual)

pb <- scuttle::aggregateAcrossCells(sce_blasts[,sce_blasts$sample_type=="Diagnosis"],
                                    ids = sce_blasts[,sce_blasts$sample_type=="Diagnosis"]$relapse_timepoint_patient )

counts_dx <- counts(pb)

# set up differential expression analysis

#1 = CR
#2= Rx

group <- factor(c(1,1,1,1,1,1,2,2,2,2,2))

# create dgelist 
y <- DGEList(counts = counts_dx,
             group = group)

# filter lowly expressed genes

keep <- filterByExpr(y)
y<- y[keep, , keep.lib.sizes = FALSE]


# create design
design <- model.matrix(~group)

# normalize
y <- normLibSizes(y)

# estimate dispersion
y <- estimateDisp(y, design)

# fit model
fit <- glmFit(y, design)


# To test MRD vs Diagnosis
lrt_dx <- glmLRT(fit, coef = "group2")

# contrasts
topde_dx <- as.data.frame(topTags(lrt_dx,adjust.method = 'BH',n = Inf))

topde_dx |>
  filter(FWER <= 0.25,
         abs(logFC) >0.5) |>
  write.csv("20250715_DX_CR_Rx.csv")


EnhancedVolcano(toptable = topde_dx,
                lab = rownames(topde_dx),
                x= "logFC",
                y = "PValue",
                pCutoff = 0.25,
                pCutoffCol = "FDR",
                FCcutoff = 0.5,
                title = 'Dx: CR vs Relapsed')

# Log-CPM values for PCA
logCPM <- cpm(y, log=TRUE)

# Run PCA
pca <- prcomp(t(logCPM), scale. = TRUE)

# Create metadata
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  group = group
)

# Plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA of Pseudobulk Samples", x = "PC1", y = "PC2") +
  theme(text = element_text(size = 14))

# DEG MRD relapse vs. Not relapsed ####
pb <- scuttle::aggregateAcrossCells(sce_blasts[,sce_blasts$sample_type=="MRD timepoints"],
                                    ids = sce_blasts[,sce_blasts$sample_type=="MRD timepoints"]$relapse_timepoint_patient )

counts_mrd <- counts(pb)

# set up differential expression analysis

#1 = CR
#2= Rx

group <- factor(c(1,1,1,1,1,2,2,2,2))

# create dgelist 
y <- DGEList(counts = counts_mrd,
             group = group)

# filter lowly expressed genes

keep <- filterByExpr(y)
y<- y[keep, , keep.lib.sizes = FALSE]


# create design
design <- model.matrix(~group)

# normalize
y <- normLibSizes(y)

# estimate dispersion
y <- estimateDisp(y, design)

# fit model
fit <- glmFit(y, design)


# To test MRD vs Diagnosis
lrt_mrd <- glmLRT(fit, coef = "group2")

# contrasts
topde_mrd <- as.data.frame(topTags(lrt_mrd,adjust.method = 'BH',n = Inf))

topde_mrd |>
  filter(FDR <= 0.25,
         abs(logFC) >0.5) |>
  write.csv("20250718_MRD_CR_Rx.csv")


EnhancedVolcano(toptable = topde_mrd,
                lab = rownames(topde_mrd),
                x= "logFC",
                y = "PValue",
                pCutoff = 0.25,
                pCutoffCol = "FDR",
                FCcutoff = 1,
                title = 'MRD: CR vs Relapsed')


# pathway analysis 

topde_mrd <- topde_mrd |>
  rownames_to_column("gene") |>
  filter(FDR <= 0.25,
         abs(logFC) > 1)


genes_mrd <-  bitr(topde_mrd$gene, fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db)

#symbols <- topde_mrd$gene

#unmapped <- setdiff(symbols, genes_mrd$SYMBOL)

topde_mrd <- topde_mrd |>
  rename(SYMBOL = gene) |>
  left_join(genes_mrd, by = "SYMBOL")


ego <- enrichGO(gene         = topde_mrd$ENTREZID,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",          # "BP"=Biological Process, "MF", "CC"
                pAdjustMethod = "holm",
                pvalueCutoff = 0.25,
                qvalueCutoff = 0.1,
                readable     = TRUE)

ereactome <- enrichPathway(gene = topde_mrd$ENTREZID, 
                                 organism = "human", 
                                 pvalueCutoff = 0.05, 
                                 readable = TRUE)  # readable = gene symbols in output

clusterProfiler::dotplot(ego)
clusterProfiler::dotplot(ereactome)
clusterProfiler::(ereactome)

# PCA plot
# Log-CPM values for PCA
logCPM <- cpm(y, log=TRUE)

# Run PCA
pca <- prcomp(t(logCPM), scale. = TRUE)

# Create metadata
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  group = group
)

# Plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA of Pseudobulk Samples", x = "PC1", y = "PC2") +
  theme(text = element_text(size = 14))

# pathway analysis




# Dimension reduction on blasts ####



highvar <- getTopHVGs(sce_blasts)

reducedDim(sce_blasts, "PCA") <- NULL
reducedDim(sce_blasts, "UMAP") <- NULL
reducedDim(sce_blasts, "TSNE") <- NULL


# PCA
sce_blasts <- runPCA(sce_blasts, 
                     ncomponents = 20, 
                     ntop = 500,
                     exprs_values = "logcounts")



# UMAP

sce_blasts <- runUMAP(sce_blasts, 
                    dimred = "PCA")

# plot UMAp colored by relapse
sce_blasts |>
  ggplot(aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = relapse), size = 0.8, alpha = 0.6) +
  facet_wrap(~sample_type) +
  theme_bw(base_size = 22)

# plot colored by individual
sce_blasts |>
  ggplot(aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = individual), size = 0.8, alpha = 0.6) +
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

colLabels(sce_blasts) <- cluster.louvein


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

sce_blasts |>
  ggplot(aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = label)) +
  facet_wrap(~sample_type) +
  theme_bw(base_size = 22)

# contingency table invidual vs. cluster label

tbl <- table(sce_blasts$individual, sce_blasts$label) |>
  as.matrix() |>
  pheatmap::pheatmap()

# normalized
tbl <- table(sce_blasts$individual, sce_blasts$label)

# Normalize by column (i.e., by label)
normalized_tbl <- prop.table(tbl, margin = 2)

# Plot the heatmap
pheatmap::pheatmap(normalized_tbl)


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


# characterization clusters ####

# find markers for all clusters
topgene_clusters <- findMarkers(sce_blasts,
                                test.type = "wilcox", 
                                direction = "up", lfc = 1,
                                pval.type = "any")


# Create an empty list to store results
entrez_dfs <- list()

# Loop over all 14 clusters
for (i in 1:14) {
  entrez_dfs[[i]] <- topgene_clusters[[i]] |>
    as.data.frame() |>
    filter(FDR < 0.05) |>
    rownames_to_column(var = "gene_symbol") |>
    pull(gene_symbol) |>
    bitr(fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
}

for (i in seq_along(entrez_dfs)) {
  write.csv(entrez_dfs[[i]],
            file = paste0("entrez_df_", i, ".csv"),
            row.names = FALSE)
}


# pahway enrichment analysis 
enrich_results <- list()

# run for each cluster
for (i in seq_along(entrez_dfs)) {
  enrich_results[[i]] <- enrichPathway(
    gene = entrez_dfs[[i]]$ENTREZID,
    organism = "human",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE
  )
}

# enrichment of BMM predicted populations in clusters ####

clusters_pops <- as.data.frame(table(sce_blasts$label,
                               sce_blasts$predicted_BMM))

# all populations
clusters_pops_df <- clusters_pops |>
  group_by(Var1) |>
  mutate(tot_pop = sum(Freq)) |>
  ungroup() |>
  group_by(Var1,Var2) |>
  mutate(rel_freq = Freq/tot_pop) |>
  ungroup() |>
  rename(
    Population_BMM = Var2,
    cluster = Var1)

# only B cell development trajectory

clusters_b_df <- clusters_pops |>
  filter(Var2 %in% b_lymphopoiesis_cells) |>
  group_by(Var1) |>
  mutate(tot_pop = sum(Freq)) |>
  ungroup() |>
  group_by(Var1,Var2) |>
  mutate(rel_freq = Freq/tot_pop) |>
  ungroup() |>
  rename(
    Population_BMM = Var2,
    cluster = Var1)


# barplots all pops
clusters_pops_df |>
  ggplot(aes(x = cluster, y = rel_freq))+
  geom_bar(stat = 'identity',
           aes(fill = Population_BMM))+
  theme(axis.text.x = element_text(size = 17,face = 'bold'))


# barplots B cell development 
clusters_b_df |>
  ggplot(aes(x = cluster, y = rel_freq))+
  geom_bar(stat = 'identity',
           aes(fill = Population_BMM))+
  theme(axis.text.x = element_text(size = 17,face = 'bold'))

# umap
sce_blasts |>
  filter(predicted_BMM %in% b_lymphopoiesis_cells) |>
  ggplot(aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = predicted_BMM), size = 0.8, alpha = 0.6) +
  facet_wrap(~sample_type) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_bw(base_size = 23)

# Cell cycle assignment ####


