# packages ####
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scran)
library(scater) 
library(scales)
library(patchwork)

# import data ####

# These steps were carreid out on the CCDL server. SCE object obtained by
# 1. converting the Seurat object to SCE
# 2. Obtaining 2000 highly variable genes
# 3. Running PCA
# 4. Running UMAP (seed == 42) with dimereduct == PCA
# 5. Adding metadata to colData (individual, timepoint, MRD_risk, relapse)

rds_path <- file.path("20250415_SCE.RDS")
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


# Gene signature ####

# calculate mean gene expression per cell of gene signature
gene_signature_path <- file.path("gene list checking in MRD.csv")
gene_signature <- read.csv(gene_signature_path)

# check that genes are in the sce
gene_signature_detected <- gene_signature$gene.list[gene_signature$gene.list%in% rownames(sce)]

# extract expression of target genes
expr_mat_signature <- logcounts(sce)[gene_signature_detected,]

# calculate mean across rows (genes)
mean_exp <- colMeans(expr_mat_signature)
sum_exp <- colSums(expr_mat_signature)

# append to coldata 
colData(sce)$signature_mean <- mean_exp
colData(sce)$signature_sum <- sum_exp

# Plots ####

# everything

# colored by cell type
p<- plotReducedDim(sce, 
                   "UMAP",
                   color_by = "Cell_Type_Experimental",
                   other_fields = c("timepoint", 
                                    "sample_type"),
                   point_size = 1,
                   theme_size = 22)

# colored by BCR paired chain (binary)
p<- plotReducedDim(sce, 
                   "UMAP",
                   color_by = "BCR_Paired_Chains",
                   other_fields = c("timepoint", 
                                    "sample_type"),
                   point_size = 1,
                   theme_size = 22)

# by expression of the genes in the gene signature
p<- plotReducedDim(sce, 
                   "UMAP",
                   color_by = "BCR_Paired_Chains",
                   other_fields = c("timepoint", 
                                    "sample_type"),
                   point_size = 1,
                   theme_size = 22)

bcr_palette <- c("True"="green",
                 "False"= "red")

p+
  facet_wrap(~sample_type)+
  scale_color_manual(values = bcr_palette)+
  guides(colour = guide_legend(override.aes = list(size = 5)))



# excluding healthy samples
p_1 <-  plotReducedDim(
  sce_patients,
  "UMAP",
  color_by = "Cell_Type_Experimental",
  other_fields = c("timepoint", "relapse", "individual"),
  point_size = 1,
  theme_size = 22
)

p_1 <-  plotReducedDim(
  sce_patients,
  "UMAP",
  color_by = "signature_mean",
  other_fields = c("timepoint", "relapse", "individual"),
  point_size = 1,
  theme_size = 22
)

# by timepoint



p_1+
  facet_wrap(~timepoint)+
  guides(colour = guide_legend(override.aes = list(size = 5)))+
  scale_colour_gradient(low = "yellow", high = "blue")

# by patient

p_1+
  facet_wrap(~individual)+
  guides(colour = guide_legend(override.aes = list(size = 5)))


# manual plotting with ggplot ####
umap_df <- as.data.frame(reducedDim(sce, "UMAP"))
umap_df$signature_mean <- colData(sce)$signature_mean
umap_df$signature_sum <- colData(sce)$signature_sum
umap_df$individual <- colData(sce)$individual
umap_df$timepoint <- colData(sce)$timepoint
umap_df$sample_type <- colData(sce)$sample_type
umap_df$relapse <- colData(sce)$relapse
umap_df$MRD_risk <- colData(sce)$MRD_risk


# plot all samples colored by sample (for sanity check)
umap_df |>
  ggplot(aes(x=UMAP1, y= UMAP2, color = individual))+
  geom_point(size = 0.8)+
  scale_color_discrete()+
  labs(color = "Patient/Donor")+
  theme_light(
    base_size = 26
  )+
  guides(colour = guide_legend(override.aes = list(size = 5)))+ #rescale legend
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot = element_rect(fill = "white", color = NA)# Remove minor gridlines
  )



# plot all ut healthy by timepoint colored by mean signature expression

# all samples by sample type
umap_df |>
  ggplot(aes(x =UMAP1, y = UMAP2, color = signature_mean ))+
  geom_point(size= 0.8)+
  scale_color_distiller(palette = "Spectral")+
  facet_wrap(~sample_type)+
  labs(color = "Mean expression signature")+
  theme_light(
    base_size = 26
  )+
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot = element_rect(fill = "white", color = NA)# Remove minor gridlines
  )

# only patients by timepoint 
umap_df |>
  filter(sample_type!="Healthy") |>
  ggplot(aes(x =UMAP1, y = UMAP2, color = signature_mean ))+
    geom_point(size= 1, shape = 16)+
  scale_color_distiller(palette = "Spectral")+
    facet_wrap(~timepoint)+
    labs(color = "Mean expression signature")+
  theme_light(
    base_size = 26
  )+
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
  )




# color by sum expression signature

umap_df |>
  ggplot(aes(x =UMAP1, y = UMAP2, color = signature_sum ))+
  geom_point(size= 0.8)+
  scale_color_viridis_c()+
  facet_wrap(~sample_type)+
  labs(color = "Sum expression signature")+
  theme_light(
    base_size = 26
  )+
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot = element_rect(fill = "white", color = NA)# Remove minor gridlines
  )
  

umap_df |>
  filter(
    sample_type != "Healthy"
  ) |>
  ggplot(aes(x =UMAP1, y = UMAP2, color = signature_sum ))+
  geom_point(size= 0.8)+
  scale_color_viridis_c()+
  facet_wrap(~timepoint)+
  labs(color = "Sum expression signature")+
  theme_light(
    base_size = 26
  )+
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot = element_rect(fill = "white", color = NA)# Remove minor gridlines
  )

# Color by each gene of the signature individually ####

# add to the umap df
for (gene in gene_signature_detected){
  print(gene)
  umap_df[[gene]] <- as.vector(assay(sce, "logcounts")[gene,])
}

# use patchwork to plot the UMAP based on the genes from the gene signature

# create list of plots

plots <- purrr::map(gene_signature_detected, function(gene_name){
  umap_df |>
    filter(sample_type != "Healthy") |>
    ggplot(aes(x = UMAP1, y=  UMAP2, color = .data[[gene_name]]))+
    geom_point(size = 0.6, alpha = 0.8)+
    scale_color_distiller(palette = "Spectral")+
    ggtitle(gene_name)+
    theme_minimal(base_size = 18)
})
grid<-wrap_plots(plots, ncol = 6, axes = "collect")
ggsave("20250415_gene signature patchwork.png",
       plot = grid, width = 2100,height =1226, unit = "px",
       scale = 5)

