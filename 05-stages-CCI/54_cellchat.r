rm(list= ls())
ls()


library(dplyr)
library(Seurat)
library(ggplot2)
library(CellChat)
library(patchwork)

resultdir<-"/hsfscqjf1/ST_CQ/P23Z32300N0005/hemingmin/bom/10.sc_merged1/results/12_CCI/stage_filter/"
figuredir<-resultdir
setwd(resultdir)

timepoint<-c("Stage_1","Stage_2","Stage_3","Stage_4","Stage_5")
for (i in timepoint){
  load(file=paste0(i,"_exp_data_cci.Rdata"))
    cellchat <- createCellChat(object = exp_t, 
                               meta = meta_t, 
                               group.by = "celltype")
    cellchat <- addMeta(cellchat, meta = meta_t)


cellchat <- setIdent(cellchat, ident.use = "celltype")

levels(cellchat@idents)

groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.human 

CellChatDB.use <-CellChatDB  #use all CellChatDB for cell-cell communication analysis

cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
# subset the expression data of signaling genes for saving computation cost


future::plan("multisession", workers = 10)


cellchat <- identifyOverExpressedGenes(cellchat)


cellchat <- identifyOverExpressedInteractions(cellchat)


cellchat <- projectData(cellchat,PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE)

cellchat <- filterCommunication(cellchat, min.cells = 10)

df_cellchat_net <- subsetCommunication(cellchat) 

write.csv(df_cellchat_net,file=paste0(resultdir,i,"_df_cellchat_net.csv"))

cellchat <- computeCommunProbPathway(cellchat)

cellchat_netP <- subsetCommunication(cellchat, slot.name = "netP") 

write.csv(cellchat_netP,file=paste0(resultdir,i,"_cellchat_netP.csv"))
cellchat <- aggregateNet(cellchat)  
saveRDS(cellchat,file=paste0(resultdir,i,"_cellchat.rds"))
                         
                         }


