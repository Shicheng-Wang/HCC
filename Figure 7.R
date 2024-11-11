############Figure 7A############
library(msigdbr) 
library(fgsea)
library(Seurat)

DefaultAssay(Tcell) <- "RNA"
markers <- FindMarkers(Tcell, ident.1 = "CD8_Cytotoxic", min.pct = 0.25, logfc.threshold = 0)  
markers$genes = rownames(markers)
cluster.genes<- markers %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes, avg_log2FC)
ranks<- deframe(cluster.genes)
####GSEA####
mdb_c2 <- msigdbr(species = "Homo sapiens", category = "C2")
fgsea_sets = mdb_c2 [grep("KEGG",mdb_c2 $gs_name),] %>% split(x = .$gene_symbol, f = .$gs_name)
length(fgsea_sets)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
write.csv(fgseaRes[,1:7], "CD8_GZMB_KEGG.csv")

markers <- FindMarkers(Tcell, ident.1 = "CD4_Cytotoxic", min.pct = 0.25, logfc.threshold = 0)  
markers$genes = rownames(markers)
cluster.genes<- markers %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes, avg_log2FC)
ranks<- deframe(cluster.genes)
####GSEA####
mdb_c2 <- msigdbr(species = "Homo sapiens", category = "C2")
fgsea_sets = mdb_c2 [grep("KEGG",mdb_c2 $gs_name),] %>% split(x = .$gene_symbol, f = .$gs_name)
length(fgsea_sets)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
write.csv(fgseaRes[,1:7], "CD4_GZMB_KEGG.csv")

############Figure 7B############
markers = c("CX3CR1","CXCR1","CXCR2","CXCR3","CXCR4","CXCR5","CXCR6","CCR1","CCR2"
            ,"CCR3","CCR4", "CCR5","CCR6",'CCR7', 'CCR8',"CCR9", "CCR10","XCR1")
DotPlot(Tcell, features = markers,
        assay='RNA' ) + 
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ 
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) 


############Figure 7C############
VlnPlot(pbmc, features = "CX3CL1")
