# Single-cell RNA-seq analysis - doublet removal

library(tidyverse)
library(ggplot2)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(dplyr)
library(Seurat)
library(patchwork)
library(clustree)
library(ggsci)
library(scales)
library(RColorBrewer)
library(DoubletFinder)

rm_doublet <- function(sample_name,input=NULL,dim.usage=30,auto="true") {

if (auto=="true") {
inpath <- paste0(input, sample_name,"/04.Matrix/FilterMatrix")#
}

sc_obj <-CreateSeuratObjsc_objt(counts = Read10X(data.dir = inpath, gene.column = 1),
                    projsc_objt = sample_name)
sc_obj <- NormalizeData(sc_obj)
sc_obj <- FindVariableFeatures(sc_obj, selsc_objtion.method = "vst", nfeatures = 3000)
sc_obj <- ScaleData(sc_obj)
sc_obj <- RunPCA(sc_obj)
sc_obj <- RunUMAP(sc_obj, dims = 1:dim.usage)

##DoubletFinder
Find_doublet <- function(data){
sweep.res.list <- paramSweep_v3(data, PCs = 1:dim.usage, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
nExp_poi <- round(0.05*ncol(data))
p<-as.numeric(as.vsc_objtor(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
return(data)
}

##
sc_obj<-Find_doublet(sc_obj)
# sc_obj<-subset(sc_obj,subset=doublet_info=="Singlet")
sc_obj@meta.data$library = sample_name 


c <- grep("pANN_",colnames(sc_obj@meta.data))
sc_obj@meta.data <- sc_obj@meta.data[,-c]

print(paste0(sample_name," cells: ", length(sc_obj@meta.data$orig.ident)," is read!"," Time:",format(Sys.time(), "%Y%m%d %X"),sep = " "))
return(sc_obj)
}



merge_seob <- function(infile,input,batch=NULL){
    seob_list <- lapply(infile, rm_doublet, input=input)
    n <- length(seob_list)

    all <- merge(seob_list[[1]], seob_list[c(2:n)])
    
    if (!is.null(batch)) {
        all$batch <- batch
    }
    return(all)
}


input<-"/hsfscqjf2/ST_CQ/Reference/Projsc_objt_data/SWU_Bombyx/scRNA/"
sample_name <- list.files(path=input)
ob1 <- merge_seob(infile=sample_name,input=input)

resultdir<-"/hsfscqjf1/ST_CQ/P23Z32300N0001/hemingmin/10.sc_merged1/results/"
figuredir<-"/hsfscqjf1/ST_CQ/P23Z32300N0001/hemingmin/10.sc_merged1/results/figures/"
save(ob1,file=paste0(resultdir,"merged_seurat.RData"))
