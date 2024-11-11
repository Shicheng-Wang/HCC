############Figure 2C############

library(multtest)
library(Seurat)
library(dplyr)
library(mindr)
library(tidyverse)
library(ggsci)
mycolor <- c(pal_npg()(10),pal_jama()(7), pal_uchicago()(7), pal_jama()(7))
pbmc <- subset(pbmc, subset = nFeature_RNA >300 & nFeature_RNA < 4000 & percent.mt < 35& nCount_RNA > 1000)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt",features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:15) 
pbmc <- FindClusters(pbmc, resolution = 0.4) 
pbmc <- RunUMAP(pbmc, dims = 1:20) 
DimPlot(pbmc, reduction = "umap", label = TRUE,cols = mycolor, raster=FALSE) 


############Figure 2D############
markers = c("PTPRC","NKG7","CD3E","CD4","CD14","CD8B","GPC3","TTR","AFP"
                            ,"FCGR3A","MS4A1", "CD37","MKI67",'MGP', 'VWF',"CD1C",
                            "CD79A","MZB1","ITGA2B","ITGB3","CLEC4C","IL3RA")

DotPlot(pbmc, features = markers,
             assay='RNA' ) + 
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ 
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) 
