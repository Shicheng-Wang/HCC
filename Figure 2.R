library(Seurat)
############Figure 2A############
library(ggplot2)
library(ggVolcano)

Bcell.markers <- FindMarkers(Bcell, ident.1 = "MTT1_7d", ident.2 = "MTT1_0d", group.by = "group_day",
                             only.pos = F, min.pct = 0.25, logfc.threshold = 0)

data <- add_regulate(Bcell.markers, log2FC_name = "avg_log2FC",
                     fdr_name = "p_val_adj",log2FC = 0.5, fdr = 0.05)
gene <- rownames(data)
data <- cbind(gene,data)

ggvolcano(data, x = "log2FoldChange", y = "padj", label = "gene", add_line = F,
          custom_label = c("IGHG1","IGHG4", "IGHG3", "IGHA1", 
                          "JCHAIN","IGHA2","XBP1","MZB1", "SEC11C", "HSP90B1","IGKC", "CD69",
                         "AFF3",
                        "AL158823.1",
                       "CXCR4",
                      "DUSP1"
    ),  output = FALSE)

############Figure 2B############
library(msigdbr) 
library(fgsea)
library(Seurat)

DefaultAssay(Bcell) <- "RNA"
markers <- FindMarkers(Bcell, ident.1 = "MTT1_7d", ident.2 = "MTT1_0d",group.by = "group_day",
                       min.pct = 0.25, logfc.threshold = 0)  
markers$genes = rownames(markers)
cluster.genes<- markers %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes, avg_log2FC)
ranks<- deframe(cluster.genes)
####GSEA####
mdb_c5 <- msigdbr(species = "Homo sapiens", category = "C5")
fgsea_sets = mdb_c5 [grep("GOBP",mdb_c5 $gs_name),] %>% split(x = .$gene_symbol, f = .$gs_name)
length(fgsea_sets)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
write.csv(fgseaRes[,1:7], "Bcell_GOBP_7dvs0d.csv")


############Figure 2C############

Bcell <- NormalizeData(Bcell, normalization.method = "LogNormalize", scale.factor = 10000)
Bcell <- FindVariableFeatures(Bcell, selection.method = "vst", nfeatures = 2000)
Bcell <- ScaleData(Bcell, features = rownames(Bcell))
Bcell <- RunPCA(Bcell, features = VariableFeatures(object = Bcell))
Bcell <- JackStraw(Bcell, num.replicate = 100)
Bcell <- ScoreJackStraw(Bcell, dims = 1:20)
ElbowPlot(Bcell)
Bcell <- FindNeighbors(Bcell, dims = 1:15) 
Bcell <- FindClusters(Bcell, resolution = 0.4) 
Bcell <- RunUMAP(Bcell, dims = 1:20)
DimPlot(Bcell, reduction = "umap", label = TRUE,  label.size = 3) 


############Figure 2E############
library(ggsignif)
VlnPlot(NK, features = c("FCGR3A"),  group.by = "group_day",raster=FALSE, y.max = 6)  + NoLegend()+
  geom_signif(comparisons =  list(c("MTT1_0d", "MTT1_7d")),y_position= 5.1,  tip_length = 0.04, vjust=0.2)+
  geom_signif(comparisons =  list(c("MTT1_0d", "MTT2_30d")),y_position= 5.5,  tip_length = 0.04, vjust=0.2) |
  VlnPlot(NK, features = c("GZMB"), group.by = "group_day",raster=FALSE, y.max = 6)+NoLegend()+
  geom_signif(comparisons =  list(c("MTT1_0d", "MTT1_7d")),y_position= 5.1,  tip_length = 0.04, vjust=0.2)+
  geom_signif(comparisons =  list(c("MTT1_0d", "MTT2_30d")),y_position= 5.5,  tip_length = 0.04, vjust=0.2) |
  VlnPlot(NK, features = c("PRF1"), group.by = "group_day",raster=FALSE, y.max = 6)+NoLegend()+
  geom_signif(comparisons =  list(c("MTT1_0d", "MTT1_7d")),y_position= 5.1,  tip_length = 0.04, vjust=0.2)+
  geom_signif(comparisons =  list(c("MTT1_0d", "MTT2_30d")),y_position= 5.5,  tip_length = 0.04, vjust=0.2)


############Figure 2F############
gene_cell_exp <- AverageExpression(cDC,
                                   features = c("IL12A","IL15","IL18","CD86","CD1D"                                   ),
                                   group.by = 'group_day',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
write.csv(gene_cell_exp,"cDC mature.csv")

############Figure 2F############
library(msigdbr) 
library(fgsea)
library(Seurat)

DefaultAssay(DCs) <- "RNA"
markers <- FindMarkers(DCs, ident.1 = "MTT1_7d", ident.2 = "MTT1_0d",group.by = "group_day",
                       min.pct = 0.25, logfc.threshold = 0)  
markers$genes = rownames(markers)
cluster.genes<- markers %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes, avg_log2FC)
ranks<- deframe(cluster.genes)
####GSEA####
mdb_c2 <- msigdbr(species = "Homo sapiens", category = "C2")
fgsea_sets = mdb_c2 [grep("KEGG",mdb_c2 $gs_name),] %>% split(x = .$gene_symbol, f = .$gs_name)
length(fgsea_sets)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
write.csv(fgseaRes[,1:7], "DCs_KEGG_7dvs0d.csv")
mdb_c2 <- msigdbr(species = "Homo sapiens", category = "C2")
fgsea_sets = mdb_c2 [grep("REACTOME",mdb_c2 $gs_name),] %>% split(x = .$gene_symbol, f = .$gs_name)
length(fgsea_sets)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
write.csv(fgseaRes[,1:7], "DCs_REACTOME_7dvs0d.csv")

############Figure 2G############
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

splitlist = SplitObject(pbmc, split.by = "group_day")
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
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") # use Secreted Signaling for cell-cell communication analysis
cellchat_0d@DB <- CellChatDB.use # set the used database in the object

cellchat_0d <- subsetData(cellchat_0d) # subset the expression data of signaling genes for saving computation cost
cellchat_0d <- identifyOverExpressedGenes(cellchat_0d)
cellchat_0d <- identifyOverExpressedInteractions(cellchat_0d)
cellchat_0d <- projectData(cellchat_0d, PPI.human)  

cellchat_0d <- computeCommunProb(cellchat_0d, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
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
cellchat_7d <- setIdent(cellchat_7d, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_7d@idents) # show factor levels of the cell labels

cellchat_7d <- addMeta(cellchat_7d, meta = identity_7d, meta.name = "labels")
cellchat_7d <- setIdent(cellchat_7d, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_7d@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_7d@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") # use Secreted Signaling for cell-cell communication analysis
cellchat_7d@DB <- CellChatDB.use # set the used database in the object


cellchat_7d <- subsetData(cellchat_7d) # subset the expression data of signaling genes for saving computation cost
cellchat_7d <- identifyOverExpressedGenes(cellchat_7d)
cellchat_7d <- identifyOverExpressedInteractions(cellchat_7d)
cellchat_7d <- projectData(cellchat_7d, PPI.human)  

cellchat_7d <- computeCommunProb(cellchat_7d, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_7d <- filterCommunication(cellchat_7d, min.cells = 10)
df_7d_net = subsetCommunication(cellchat_7d)
write.csv(df_7d_net, "net_7d.csv")
cellchat_7d <- computeCommunProbPathway(cellchat_7d)
df_7d_netp = subsetCommunication(cellchat_7d, slot.name = "netp")
write.csv(df_7d_netp, "netp_7d.csv")
cellchat_7d <- aggregateNet(cellchat_7d)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_0d <- filterCommunication(cellchat_0d, min.cells = 10)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_7d <- filterCommunication(cellchat_7d, min.cells = 10)
object.list <- list(MTT_0d = cellchat_0d, MTT_7d = cellchat_7d)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
netVisual_bubble(cellchat, sources.use = c(10), targets.use = c(2,4), 
                       max.dataset = 2, title.name = "Increased signaling in MTT_7D", comparison = c(1, 2), angle.x = 45)











