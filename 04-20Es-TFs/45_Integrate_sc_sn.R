# scRNA-seq and snRNA-seq Integrated and distance analysis

.libPaths(c("/hsfscqjf2/ST_CQ/Reference/software/envs/monocle3/lib/R/library",
            "/hsfscqjf2/ST_CQ/Reference/software/envs/RTest/lib/R/library"))
suppressPackageStartupMessages({
  library(Matrix)
  library(cowplot)
  library(RCurl)
  library(Seurat)
  library(patchwork)
  library(cluster)
  library(future)
  library(ggsci)
  library(scales)
  library(RColorBrewer)
  library(tidyverse)
  library(ggplot2)
  library(ComplexHeatmap)
  library(harmony)})
pal=c("Cell_morphogenesis"='#1f77b4',
      "Epithelial_1"='#279e68',
      "Epithelial_2"='#98df8a',
      "Epithelial_3"='#b5bd61',
      "Cuticle_1"='#d62728',
      "Cuticle_2"='#e377c2',
      "Cuticle_3"='#ff9896',
      "Immune"='#f7b6d2',
      "Apoptosis"='#aa40fc',
      "Metabolic_process"='#c49c94',
      "Axon_development"='#ffbb78',
      "Ciliated_cell"='#ff7f0e')
pal2=c("scRNA"="#D45651",
       "snRNA"="#A1CFFA")
order =c( "Ciliated_cell",
          "Axon_development",
          "Metabolic_process",
          "Apoptosis",
          "Immune",
          "Cuticle_3",
          "Cuticle_2",
          "Cuticle_1",
          "Epithelial_3",
          "Epithelial_2",
          "Epithelial_1",
          "Cell_morphogenesis"
)
order_sn <- c("Cell_morphogenesis","Epithelial_1","Cuticle_1","Cuticle_2","Immune",
              "Apoptosis")


snRNA<-readRDS("./sc_snRNA/anno.rds")
snRNA@meta.data$celltype <-factor(snRNA@meta.data$celltype,levels = order_sn)
Idents(snRNA) <- "celltype"

scRNA<-readRDS("./sc_snRNA/anno_final.rds")
scRNA@meta.data$celltype <-scRNA@meta.data$celltype_v1
scRNA@meta.data$celltype <- factor(scRNA@meta.data$celltype,levels = rev(order))
Idents(scRNA) <- "celltype"

combined <- merge(scRNA, y = snRNA, add.cell.ids = c("scRNA", "snRNA"))
combined@meta.data$orig.ident[is.na(combined@meta.data$orig.ident)] <- "scRNA"
combined@meta.data$orig.ident[combined@meta.data$orig.ident=="SeuratProject"] <- "snRNA"
combined@meta.data <- combined@meta.data[, colSums(is.na(combined@meta.data)) == 0]

combined <- NormalizeData(object = combined, normalization.method = "LogNormalize", scale.factor = 10000)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(combined)
combined <- ScaleData(combined, features = all.genes)

combined <- RunPCA(combined, features = VariableFeatures(object = combined),verbose = FALSE,reduction.name = "pca")

combined <- RunHarmony(combined,reduction="pca",group.by.vars = "orig.ident",reduction.save = "harmony")
combined <- RunUMAP(combined, dims = 1:30,reduction = "harmony",reduction.name = "umap")
saveRDS(combined,file="./sc_snRNA/combined.rds")

data <- as.data.frame(combined[["umap"]]@cell.embeddings)
data$celltype <- combined$celltype
data$celltype <- factor(data$celltype,levels = rev(order))
data$group <- combined$orig.ident
colnames(data)[1:2] <- c("UMAP_1","UMAP_2")
pdf("./Umap_celltype.pdf",width = 10,height = 8)
ggplot(data, aes(x = UMAP_1, 
                 y = UMAP_2, 
                 fill = celltype,
                 color = celltype)) +
  geom_point(size = 1) +
  theme_classic()+
  theme(axis.text = element_text(colour = 'black', size = 12))+
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal)
dev.off()

pdf("./Umap_group.pdf",width = 10,height = 8)
ggplot(data, aes(x = UMAP_1, 
                 y = UMAP_2, 
                 fill = group,
                 color = group)) +
  geom_point(size = 1) +
  theme_classic()+
  theme(axis.text = element_text(colour = 'black', size = 12))+
  scale_fill_manual(values = pal2) +
  scale_color_manual(values = pal2)
dev.off()

# Extract the PCA coordinates
pca_coords <- Embeddings(combined, "pca")

meta_data <- combined@meta.data
meta_data$pca1 <- pca_coords[, 1]
meta_data$pca2 <- pca_coords[, 2]

# timepoint groups
sn_time <- c("0h","10min","20min", "30min","40min","50min",
             "60min","4h","6h")
time_points <- c("0h","10min","20min", "30min","40min","50min",
                 "60min","4h","6h","L5D1", "L5D2", "L5D3", "L5D4", 
                 "L5D5", "L5D6","L5D7","WD1","WD2","P6")
sc_time <- c("L5D1", "L5D2", "L5D3", "L5D4", 
             "L5D5", "L5D6","L5D7","WD1","WD2","P6")
# Compute the centroid coordinates at each time point
centroids <- sapply(time_points, function(tp) {
  subset <- meta_data[meta_data$timepoint == tp, ]
  c(mean(subset$pca1), mean(subset$pca2))  
})
centroids <- t(centroids)
rownames(centroids) <- time_points

#lorentzian distances
final_dis <-c()
for (i in sn_time){
  sub_time <- c(i,sc_time)
  centroids2 <- centroids[sub_time,]
  distance_matrix <- distance(centroids2, method = "lorentzian")
  rownames(distance_matrix) <- rownames(centroids2)
  colnames(distance_matrix) <- rownames(centroids2)
  distances <- distance_matrix[i, -1,drop=F]
  distances <- t(distances)
  final_dis <- cbind(final_dis,distances)
}
write.csv(final_dis,file = "./final_lorentzian_pca.csv")