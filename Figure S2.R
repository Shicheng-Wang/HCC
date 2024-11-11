############Figure S2############
markers = c("CD274","ARG1","MRC1","CCL8","CD163","CLEC10A","CD9")

DotPlot(pbmc, features = markers,group.by = "sample", assay='RNA' ) + theme(panel.grid = element_blank(), 
                                                                            axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ 
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) 
