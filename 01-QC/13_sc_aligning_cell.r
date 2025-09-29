# Single-cell RNA-seq analysis - integrated analysis

library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(tidyverse)
library(ggplot2)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(dplyr)
library(patchwork)
library(clustree)
library(Seurat)
library(harmony)
resultdir<-"/hsfscqjf1/ST_CQ/P23Z32300N0001/hemingmin/10.sc_merged1/results/"

setwd(resultdir)



load(file = "seurat_filtered.RData")


filtered_seurat<- SCTransform(filtered_seurat) # 


seurat_integrated<-RunPCA(object = filtered_seurat)

seurat_integrated <-seurat_integrated %>% RunHarmony(group.by.vars="sample", plot_convergence = TRUE)


save(seurat_integrated,file = "seurat_integrated_seurat.RData")




