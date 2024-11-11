library(Seurat)
############Figure 3A############
library(msigdbr) 
library(fgsea)
library(Seurat)

DefaultAssay(DCs) <- "RNA"
markers <- FindMarkers(DCs, ident.1 = "MTT1_7d", ident.2 = "MTT1_0d",group.by = "group_day",
                       min.pct = 0.1, logfc.threshold = 0)  
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


############Figure 3B############
gene_cell_exp <- AverageExpression(cDC,
                                   features = c("IL12A","IL15","IL18","CD86","CD1D"                                   ),
                                   group.by = 'group_day',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
write.csv(gene_cell_exp,"cDC mature.csv")


############Figure 3C############
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


############Figure 3D############
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


############Figure 3E############
library(msigdbr) 
library(fgsea)
library(Seurat)

DefaultAssay(Bcell) <- "RNA"
markers <- FindMarkers(Bcell, ident.1 = "MTT1_7d", ident.2 = "MTT1_0d",group.by = "group_day",
                       min.pct = 0.1, logfc.threshold = 0)  
markers$genes = rownames(markers)
cluster.genes<- markers %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes, avg_log2FC)
ranks<- deframe(cluster.genes)
####GSEA####
mdb_c5 <- msigdbr(species = "Homo sapiens", category = "C5")
fgsea_sets = mdb_c5 [grep("GOBP",mdb_c5 $gs_name),] %>% split(x = .$gene_symbol, f = .$gs_name)
length(fgsea_sets)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
write.csv(fgseaRes[,1:7], "Bcell_GOBP_7dvs0d.csv")

############Figure 3F############
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

############Figure 3H############

library(msigdbr) 
library(fgsea)
library(Seurat)

DefaultAssay(Bcell) <- "RNA"
markers <- FindMarkers(Tcells, ident.1 = "MTT1_7d", ident.2 = "MTT1_0d",group.by = "group_day",
                       min.pct = 0.1, logfc.threshold = 0)  
markers$genes = rownames(markers)
cluster.genes<- markers %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes, avg_log2FC)
ranks<- deframe(cluster.genes)
####GSEA####
mdb_c5 <- msigdbr(species = "Homo sapiens", category = "C5")
fgsea_sets = mdb_c5 [grep("GOBP",mdb_c5 $gs_name),] %>% split(x = .$gene_symbol, f = .$gs_name)
length(fgsea_sets)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
write.csv(fgseaRes[,1:7], "Bcell_GOBP_7dvs0d.csv")
