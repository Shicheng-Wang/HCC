############Figure 5A############
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
