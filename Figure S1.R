############Figure S1A############
library(CellChat)
library(ggalluvial)
library(svglite)
library(multtest)
library(dplyr)
library(mindr)
library(tidyverse)
.libPaths(c("~/SeuratObjectV5.0.0", .libPaths()))
unloadNamespace("Seurat")
library(Seurat)
library(SeuratObject)
library(SeuratData)

pbmc_1 <- subset(x = pbmc, idents = c( "Mono"  ,  "NCMono" , "Bcell" ,  
                                       "cDC" ,"Plasma" ,"pDC"))
TP <- merge(expanded_TCR, y = pbmc_1, add.cell.ids = c("T", "P"), project = "TP")
splitlist = SplitObject(TP, split.by = "group_day")
data.input_0d = splitlist$MTT1_0d@assays$RNA@data
identity_0d = data.frame(group = splitlist$MTT1_0d@active.ident, row.names = names(splitlist$MTT1_0d@active.ident))
unique(identity_0d$group)
cellchat_0d = createCellChat(data.input_0d)

cellchat_0d <- addMeta(cellchat_0d, meta = identity_0d, meta.name = "labels")
cellchat_0d <- setIdent(cellchat_0d, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_0d@idents) # show factor levels of the cell labels

cellchat_0d <- addMeta(cellchat_0d, meta = identity_0d, meta.name = "labels")
cellchat_0d <- setIdent(cellchat_0d, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_0d@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_0d@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat_0d@DB <- CellChatDB.use 

cellchat_0d <- subsetData(cellchat_0d) 
cellchat_0d <- identifyOverExpressedGenes(cellchat_0d)
cellchat_0d <- identifyOverExpressedInteractions(cellchat_0d)
cellchat_0d <- projectData(cellchat_0d, PPI.human)  

cellchat_0d <- computeCommunProb(cellchat_0d, raw.use = TRUE)
cellchat_0d <- filterCommunication(cellchat_0d, min.cells = 10)
df_0d_net = subsetCommunication(cellchat_0d)
write.csv(df_0d_net, "net_0d.csv")
cellchat_0d <- computeCommunProbPathway(cellchat_0d)
df_0d_netp = subsetCommunication(cellchat_0d, slot.name = "netp")
write.csv(df_0d_netp, "netp_0d.csv")
cellchat_0d <- aggregateNet(cellchat_0d)

data.input_7d = splitlist$MTT1_7d@assays$RNA@data
identity_7d = data.frame(group = splitlist$MTT1_7d@active.ident, row.names = names(splitlist$MTT1_7d@active.ident))
unique(identity_7d$group)
cellchat_7d = createCellChat(data.input_7d)

cellchat_7d <- addMeta(cellchat_7d, meta = identity_7d, meta.name = "labels")
cellchat_7d <- setIdent(cellchat_7d, ident.use = "labels") 
levels(cellchat_7d@idents) 

cellchat_7d <- addMeta(cellchat_7d, meta = identity_7d, meta.name = "labels")
cellchat_7d <- setIdent(cellchat_7d, ident.use = "labels") 
levels(cellchat_7d@idents) 
groupSize <- as.numeric(table(cellchat_7d@idents)) 

CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat_7d@DB <- CellChatDB.use 


cellchat_7d <- subsetData(cellchat_7d)
cellchat_7d <- identifyOverExpressedGenes(cellchat_7d)
cellchat_7d <- identifyOverExpressedInteractions(cellchat_7d)
cellchat_7d <- projectData(cellchat_7d, PPI.human)  

cellchat_7d <- computeCommunProb(cellchat_7d, raw.use = TRUE)
cellchat_7d <- filterCommunication(cellchat_7d, min.cells = 10)
df_7d_net = subsetCommunication(cellchat_7d)
write.csv(df_7d_net, "net_7d.csv")
cellchat_7d <- computeCommunProbPathway(cellchat_7d)
df_7d_netp = subsetCommunication(cellchat_7d, slot.name = "netp")
write.csv(df_7d_netp, "netp_7d.csv")
cellchat_7d <- aggregateNet(cellchat_7d)
cellchat_0d <- filterCommunication(cellchat_0d, min.cells = 10)
cellchat_7d <- filterCommunication(cellchat_7d, min.cells = 10)
cellchat_0d <- netAnalysis_computeCentrality(cellchat_0d, slot.name = "netP")
cellchat_7d <- netAnalysis_computeCentrality(cellchat_7d, slot.name = "netP")
object.list <- list(MTT_0d = cellchat_0d, MTT_7d = cellchat_7d)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

####netAnalysis_signalingRole_scatter
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) 
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)+ 
    scale_y_continuous(limits = c(0,2.5)) +  scale_x_continuous(limits = c(0,2))
} 
patchwork::wrap_plots(plots = gg)

####number and strength
gg1 <- netVisual_heatmap(cellchat, comparison = c(1,2))
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat,comparison = c(1,2) ,measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
