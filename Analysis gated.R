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
sce <- sce[,order]

sce$gate <- gated_df$OmiqFilter

umap_df$gate <- gated_df$OmiqFilter

umap_df |> 
  dplyr::filter(sample_type != "Healthy") |>
  ggplot(aes(x=UMAP1, y=UMAP2, color = gate))+
    geom_point(alpha = 0.6)+
  facet_wrap(~sample_type)
