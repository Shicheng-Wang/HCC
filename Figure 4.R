############Figure 4A############

library(multtest)
library(Seurat)
library(dplyr)
library(mindr)
library(tidyverse)
library(ggsci)
mycolor <- c(pal_npg()(10),pal_jama()(7), pal_uchicago()(7), pal_jama()(7))

Tcell <- NormalizeData(Tcell, normalization.method = "LogNormalize", scale.factor = 10000)
Tcell <- FindVariableFeatures(Tcell, selection.method = "vst", nfeatures = 2000)
DefaultAssay(Tcell) <- "integrated"
Tcell <- ScaleData(Tcell, features = rownames(Tcell))
Tcell <- RunPCA(Tcell, features = VariableFeatures(object = Tcell))
Tcell <- JackStraw(Tcell, num.replicate = 100)
Tcell <- ScoreJackStraw(Tcell, dims = 1:20)
ElbowPlot(Tcell)
Tcell <- FindNeighbors(Tcell, dims = 1:15) 
Tcell <- FindClusters(Tcell, resolution = 0.7) 
Tcell <- RunUMAP(Tcell, dims = 1:30) 
DimPlot(Tcell, reduction = "umap", label = TRUE,cols = mycolor, raster=FALSE) + NoLegend()

############Figure 4B############
gene_cell_exp <- AverageExpression(Tcell,
                                   features = c("CD4","CD8A",
                                                "LTB","TNFSF10","TNFRSF4","GATA3","IL7R",
                                                "TCF7","CCR7","SELL","LEF1",
                                                "GZMB","PRF1","FGFBP2","CX3CR1","ZNF683","GNLY",
                                                "TRDC","TRGC1","TRGC2",
                                                "GZMK","GZMA","CMC1","EOMES",
                                                "SLC4A10","KLRB1","CXCR6","NCR3","CEBPD",
                                                "FOXP3","IL2RA","CTLA4","TNFRSF18","LIMS1",
                                                "XCL1","XCL2",
                                                "HSPA1A","HSPA1B","DUSP1","HSP90AA1","DNAJB1",
                                                "PDCD1","LAG3","TOX","TIGIT","HAVCR2"),
                                   group.by = 'ident',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)

library(ComplexHeatmap)
col_cluster <- setNames(c(rep("#9ECABE",2), rep("#F6F5B4",5),
                          rep("gray",4), rep("#E3AD68",6),
                          rep("#9ECABE",3), rep("#F6F5B4",4),
                          rep("gray",5), rep("#E3AD68",5),
                          rep("#9ECABE",2), rep("#F6F5B4",5),
                          rep("gray",5)),rownames(marker_exp))

row_info = rowAnnotation(foo = anno_text(rownames(marker_exp), 
                                         location = 0, 
                                         just = "left",
                                         gp = gpar(fill = col_cluster, 
                                                   col = "black"),
                                         width = max_text_width(rownames(marker_exp))*1.2))
Heatmap(marker_exp, cluster_rows = F, cluster_columns = F, show_column_names = F, show_row_names = T,
        column_title = NULL, heatmap_legend_param = list( title=' '), 
        col = colorRampPalette(c("#87CEFA","white","red"))(100), border = 'black',rect_gp = gpar(col = "black", lwd = 1),
        row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),)+row_info

############Figure 4C-D############
library(ggplot2)
library(UCell)
library(ggsignif)
markers <- list()
markers$naive <- c("TCF7", "CCR7", "SELL", "LEF1")
markers$dysfunction <- c("PDCD1", "CTLA4", "TIGIT", "HAVCR2", "LAG3", "LAYN","TOX")
markers$cytotoxic <- c("PRF1","IFNG", "GZMA", "GZMB", "GZMH", "GNLY",  "NKG7",
                       "KLRK1","KLRB1","KLRD1","CTSW","CST7")
Tcell_Ucell <- AddModuleScore_UCell(Tcell, features = markers)
signature.names <- paste0(names(markers), "_UCell")
library(ggsignif)
library(ggsci)
mycolor <- c(pal_npg()(10),pal_jama()(7), pal_uchicago()(7), pal_jama()(7))
library(ggrastr)
library(dplyr)
data<- FetchData(Tcell_Ucell,vars = c("cluster","naive_UCell"))
data$cellid <- case_when(data$cluster ==unique(data$cluster)[1] ~ "CD4_LIMS1",
                         data$cluster ==unique(data$cluster)[2] ~ 'CD4_GZMB',
                         data$cluster ==unique(data$cluster)[3] ~ 'T_LTB',
                         data$cluster ==unique(data$cluster)[4] ~ 'CD4_CCR7',
                         data$cluster ==unique(data$cluster)[5] ~ 'MAIT',
                         data$cluster ==unique(data$cluster)[6] ~ 'GDT',
                         data$cluster ==unique(data$cluster)[7] ~ 'CD8_GZMK',
                         data$cluster ==unique(data$cluster)[8] ~ 'CD8_GZMB',
                         data$cluster ==unique(data$cluster)[9] ~ 'CD4_FOXP3',
                         data$cluster ==unique(data$cluster)[10] ~ 'CD8_CCR7',
                         data$cluster ==unique(data$cluster)[11] ~ 'CD8_XCL1',
                         data$cluster ==unique(data$cluster)[12] ~ 'CD8_HSPA1A'
)
colors <- mycolor

ggplot(data, aes(x=cellid,y=naive_UCell,fill=cellid,color=cellid) ) +
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(size=1),
        axis.text.y = element_text(size=10),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')+
  labs(x=NULL,y=NULL,title = "naive_UCell")

############Figure 4G############
library(multtest)
library(Seurat)
library(dplyr)
library(mindr)
library(tidyverse)
library(ggsci)
mycolor <- c(pal_npg()(10),pal_jama()(7), pal_uchicago()(7), pal_jama()(7))
DimPlot(Tcell, reduction = "umap", label = TRUE,cols = mycolor, raster=FALSE) + NoLegend()

FeaturePlot(Tcell, features = c("CX3CR1","ADGRG1"), label = T, ncol = 2,raster=FALSE,
            pt.size = (1)) &scale_colour_gradient(low="#87CEFA", high="red")

