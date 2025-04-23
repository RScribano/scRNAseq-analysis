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


# Fix labels MRD_43 ####
# in cartridge F tags 6 and 7 have been swapped around so that Dx is swapped with D33. we need to change that 
tags <- sce$ident[sce$orig.ident == "cart_f"]

tags_new <- dplyr::case_when(tags == "SampleTag07_hs" ~ "SampleTag06_hs",
                             tags == "SampleTag06_hs" ~ "SampleTag07_hs")

sce$ident[sce$orig.ident == "cart_f"] <- tags_new

# now we need t oswitch around timepoint and sample type

timepoints <- as.character(sce$timepoint[sce$individual == "MRD_43"])
sample_type <- as.character(sce$sample_type[sce$individual == "MRD_43"])
sample_names <- as.character(sce$Sample_Name[sce$individual == "MRD_43"])


timepoints_new <- dplyr::case_when(
  timepoints == "Dx" ~ "D33",
  timepoints == "D33" ~ "Dx",
  TRUE ~ timepoints
)

sampletype_new <- dplyr::case_when(
  sample_type == "Diagnosis" ~ "MRD timepoints",
  sample_type == "MRD timepoints" ~ "Diagnosis",
  TRUE ~ sample_type
)

samplenames_new <- dplyr::case_when(
  sample_names == "MRD43_Dx" ~ "MRD43_D33",
  sample_names == "MRD43_D33" ~ "MRD43_Dx",
  TRUE ~ sample_names
)
  
sce$timepoint[sce$individual == "MRD_43"] <- timepoints_new
sce$sample_type[sce$individual == "MRD_43"] <- sampletype_new
sce$Sample_Name[sce$individual == "MRD_43"] <- samplenames_new


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



# make UMAP df ####

umap_df <- as.data.frame(reducedDim(sce, "UMAP"))
umap_df$signature_mean <- colData(sce)$signature_mean
umap_df$signature_sum <- colData(sce)$signature_sum
umap_df$individual <- colData(sce)$individual
umap_df$timepoint <- colData(sce)$timepoint
umap_df$sample_type <- colData(sce)$sample_type
umap_df$relapse <- colData(sce)$relapse
umap_df$MRD_risk <- colData(sce)$MRD_risk
umap_df$Cell_Type_Experimental <- colData(sce)$Cell_Type_Experimental


# add to the umap df
for (gene in gene_signature_detected){
  print(gene)
  umap_df[[gene]] <- as.vector(assay(sce, "logcounts")[gene,])
}

#save as .csv
write.csv(umap_df, "20250423_umapd_df.csv")

# manual plotting with ggplot ####
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



# use patchwork to plot the UMAP based on the genes from the gene signature

# create list of plots

plots <- purrr::map(gene_signature_detected, function(gene_name){
  umap_df |>
    filter(sample_type != "Healthy") |>
    ggplot(aes(x = UMAP1, y=  UMAP2, color = .data[[gene_name]], alpha = .data[[gene_name]]))+
    geom_point(size = 0.6)+
    scale_color_distiller(palette = "Spectral")+
    ggtitle(gene_name)+
    theme_minimal(base_size = 18)
})
grid<-wrap_plots(plots, ncol = 6, axes = "collect")
ggsave("20250415_gene signature patchwork_1.png",
       plot = grid, width = 2100,height =1226, unit = "px",
       scale = 5)


# two layers
plots <- purrr::map(gene_signature_detected, function(gene_name) {
  umap_df |>
    filter(sample_type != "Healthy") |>
    ggplot(aes(x = UMAP1, y = UMAP2)) +
    
    # Layer 1: all points, grey
    geom_point(color = "grey90", size = 0.6) +
    
    # Layer 2: only gene > 0, with color + alpha
    geom_point(
      
      aes(color = .data[[gene_name]]),
      size = 0.6
    ) +
    
    scale_color_distiller(palette = "Blues", direction = 1) +  # or use scale_color_steps()
    guides(alpha = "none") +  # Hide alpha legend if too noisy
    ggtitle(gene_name) +
    theme_minimal(base_size = 18)
})

grid<-wrap_plots(plots, ncol = 6, axes = "collect")
ggsave("20250415_gene signature patchwork_3.png",
       plot = grid, width = 2100,height =1226, unit = "px",
       scale = 5)

# Cell signature scoring using UCell ####
signature <- list(
  metabolic_genes= c("CYCS","CYC1", "CS","NDUFA4","NDUFC1","NDUFS2","NDUFS5","NDUFV1","NDUFB3", 
                      "NDUFA13","NDUFB8","NDUFB10","NDUFA3","NDUFB9","NDUFB11","NDUFS1","NDUFA11","NDUFS8", 
                      "NDUFA2","NDUFA9","MRPL3","MRPS21","MRPL17","MRPS30","MRPL16")
)
sce <- ScoreSignatures_UCell(sce, features = signature,
                             assay = 'logcounts', name = NULL)

# export as dataframe column

umap_df$UCell_signature_score <- t(as.data.frame(assay(altExp(sce, "UCell"))))

# plottini..

ggplot(data = umap_df, aes(x= UMAP1, y=UMAP2))+
  geom_point(size = 0.7, aes(color = UCell_signature_score))+
  scale_color_gradient(low = "grey70", high = "purple")+
  theme_minimal(base_size = 22)+
  facet_grid(~sample_type, cols = 2, rows = NULL)

#barplot cell signature score in B-cells only



umap_df |> 
  group_by(sample_type, Cell_Type_Experimental)|>
  summarise(mean_UCell_score = 
              mean(UCell_signature_score),
            se = sd(UCell_signature_score)/sqrt(n()),
            .groups = "drop") |>
  ggplot(aes(x = Cell_Type_Experimental, y =mean_UCell_score ))+
  geom_bar(aes(fill = Cell_Type_Experimental), stat = "identity", position =position_dodge())+
  geom_errorbar(aes(ymin = mean_UCell_score - se, 
                    ymax = mean_UCell_score + se),
                position = position_dodge(),
                width = 0.2)+
  scale_fill_brewer(palette = "Set3")+
  facet_wrap(~sample_type)+
  theme_minimal(base_size = 22)+
  theme( axis.text.x = element_blank())
 

# working on ADT expression ####

#fixing names 
original_tags <- rownames(altExp(sce, "ADT"))

new_tags <- original_tags |>
  sub("-.*", "", x = _) |>
  sub(":.*", "", x= _)
  
rowData(altExp(sce, "ADT"))$antigen <- new_tags

rownames(altExp(sce, "ADT")) <- rowData(altExp(sce, "ADT"))$antigen

# extract adt expression as dataframe
adt_df <- as.data.frame(t(logcounts(altExp(sce, "ADT"))))

# add key for sce metadata 
adt_df$Sample_Name <- coldata_df$Sample_Name

#export just in case
write.csv(adt_df, "20250415_adt_exp.csv")

# plot distributions ADT #### 

# merge with metadata (we didn't save it already merged to spare memory)

adt_df$UMAP1 <- umap_df$UMAP1
adt_df$UMAP2 <- umap_df$UMAP2
adt_df$individual <- colData(sce)$individual
adt_df$timepoint <- colData(sce)$timepoint
adt_df$sample_type <- colData(sce)$sample_type
adt_df$relapse <- colData(sce)$relapse
adt_df$MRD_risk <- colData(sce)$MRD_risk



# plot distirbutions by sample type


#pivot to long format
adt_long <- adt_df |>
  tidyr::pivot_longer(
    cols = -c("individual", "timepoint","sample_type", "relapse", "MRD_risk", "Sample_Name"),
    names_to = "marker",
    values_to = "logcount"
    )

# joyplot/ridges
  ggplot(adt_long, aes(x = logcount, y= marker, fill = marker))+
    geom_density_ridges(panel_scaling = TRUE)+
    facet_grid(~sample_type)
  
  
  # joyplot/ridges
  
  adt_long |>
    filter(sample_type != "Healthy")|>
  ggplot(aes(x = logcount, y= marker, fill = marker))+
    geom_density_ridges(panel_scaling = TRUE)+
    theme_minimal(base_size = 24)+
    labs(rota)
    facet_grid(~timepoint)


# UMAP colored by ADT ####


# plot umaps of all patients (not healthy)
plots <- purrr::map(colnames(adt_df[1:15]), function(adt) {
  adt_df |>
    filter(sample_type != "Healthy") |>
    ggplot(aes(x = UMAP1, y = UMAP2)) +
             
             # first layer grey
    geom_point(size = 0.6, color = "grey80", alpha = 0.5) +
      
    geom_point(size = 0.6,
               aes(color = .data[[adt]], alpha = .data[[adt]]))+
    scale_color_distiller(palette = "Reds", direction = 1) +
    guides(alpha = "none") + 
    ggtitle(adt) +
    theme_minimal(base_size = 18)
    })

grid<-wrap_plots(plots, ncol = 6, axes = "collect")
ggsave("20250415_umap_adt_1.png",
       plot = grid, width = 2100,height =1226, unit = "px",
       scale = 5)

# plot for each gene divded by sample type ####
plots <- purrr::map(gene_signature_detected, function(gene_name){
  umap_df |>
    ggplot(aes(x = UMAP1, y=  UMAP2, color = .data[[gene_name]], alpha = .data[[gene_name]]))+
    geom_point(size = 0.6)+
    scale_color_gradient(low= "grey70", high = "blue")+
    guides(alpha = "none")+
    ggtitle(gene_name)+
    facet_wrap(~sample_type, ncol = 4)+
    theme_minimal(base_size = 20)
})

purrr::walk2(
  .x = plots,
  .y = gene_signature_detected,
  .f = ~ggsave(
    filename = paste0("20250416_", .y, ".png"),
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


ggplot(data = umap_df,aes(x = UMAP1, y=  UMAP2))+
  geom_point(size = 0.6, aes(color = CYCS, alpha = CYCS))+
  scale_color_gradient(low = "grey70", high = "blue")+
  guides(alpha = "none")+
  ggtitle("CYCS")+
  facet_wrap(~ sample_type, ncol = 4)+
  theme_minimal(base_size = 18)


