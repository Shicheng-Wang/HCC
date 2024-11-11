############Figure 9A############
library(immunarch)
data <-repLoad("TCR", .coding = TRUE ,)
tc1 <- trackClonotypes(data$data, list(1, 50), .col = "aa")
vis(tc1)

############Figure 9D-G############
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
cellchat_0d <- setIdent(cellchat_0d, ident.use = "labels") 
levels(cellchat_0d@idents) 

cellchat_0d <- addMeta(cellchat_0d, meta = identity_0d, meta.name = "labels")
cellchat_0d <- setIdent(cellchat_0d, ident.use = "labels")
levels(cellchat_0d@idents) 
groupSize <- as.numeric(table(cellchat_0d@idents)) 

CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") 
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
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact")
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
    scale_y_continuous(limits = c(0,12)) +  scale_x_continuous(limits = c(0,6))
} 
patchwork::wrap_plots(plots = gg)

####number and strength
gg1 <- netVisual_heatmap(cellchat, comparison = c(1,2))
gg2 <- netVisual_heatmap(cellchat,comparison = c(1,2) ,measure = "weight")
gg1 + gg2

####information flow
gg3 <- rankNet(cellchat, mode = "comparison",comparison = c(2,1), stacked = F, do.stat = TRUE)
gg3

####MHC-I pathway network
pathways.show <- c("MHC-I") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
} 


############Figure 9H############
mdb_c2 <- msigdbr(species = "Homo sapiens", category = "C5")## 定义基因集，选取C2
fgsea_sets = mdb_c2 [grep("GO",mdb_c2 $gs_name),] %>% split(x = .$gene_symbol, f = .$gs_name)
length(fgsea_sets)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000) #运行fgsea
plotEnrichment(fgsea_sets[["GOMF_MHC_CLASS_I_PROTEIN_BINDING"]],ranks) + labs(title="GOMF_MHC_CLASS_I_PROTEIN_BINDING")
plotEnrichment(fgsea_sets[["GOMF_ANTIGEN_BINDING"]],ranks) + labs(title="GOMF_ANTIGEN_BINDING")
plotEnrichment(fgsea_sets[["GOMF_MHC_PROTEIN_BINDING"]],ranks) + labs(title="GOMF_MHC_PROTEIN_BINDING")
plotEnrichment(fgsea_sets[["GOMF_PEPTIDE_ANTIGEN_BINDING"]],ranks) + labs(title="GOMF_PEPTIDE_ANTIGEN_BINDING")
plotEnrichment(fgsea_sets[["GOBP_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY"]],ranks) + 
                                labs(title="GOBP_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY")
plotEnrichment(fgsea_sets[["GOCC_MHC_CLASS_II_PROTEIN_COMPLEX"]],ranks) + labs(title="GOCC_MHC_CLASS_II_PROTEIN_COMPLEX")

############Figure 9I############

library(ggplot2)
library(UCell)
library(ggsignif)
markers <- list()
markers$cytotoxic <- c("PRF1","IFNG", "GZMA", "GZMB", "GZMH", "GNLY",  "NKG7",
                       "KLRK1","KLRB1","KLRD1","CTSW","CST7")
expanded_TCR_UCell <- AddModuleScore_UCell(Tcell, features = markers)
signature.names <- paste0(names(markers), "_UCell")
library(ggsignif)
library(ggsci)
mycolor <- c(pal_npg()(10),pal_jama()(7), pal_uchicago()(7), pal_jama()(7))
library(ggrastr)
library(dplyr)
data<-  FetchData(expanded_TCR_UCell,vars = c("group_day","cytotoxic_UCell"))
data$cellid <- case_when(data$group_day ==unique(data$group_day)[1] ~ "MTT1_0d",
                         data$group_day ==unique(data$group_day)[2] ~ 'MTT1_7d',
                         data$group_day ==unique(data$group_day)[3] ~ 'MTT2_30d'
)
colors <- mycolor

ggplot(data, aes(x=cellid,y=cytotoxic_UCell,fill=cellid,color=cellid) ) +
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(size=1),
        axis.text.y = element_text(size=10),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')+
  labs(x=NULL,y=NULL,title = "cytotoxic_UCell")+ 
  geom_jitter_rast(col="#00000033", pch=19,cex=2, position = position_jitter(0.2))+
  geom_boxplot(position=position_dodge(0),outliers = F,)+
  scale_fill_manual(values = colors)+
  geom_boxplot(position=position_dodge(0),color='black',
               outlier.colour  = NA,outlier.fill=NA,outlier.shape=NA)+
  geom_signif(comparisons =  list(c("MTT1_0d", "MTT1_7d")),y_position= 0.78,  tip_length = 0.04, vjust=0.2)+
  geom_signif(comparisons =  list(c("MTT1_0d", "MTT2_30d")),y_position= 0.83,  tip_length = 0.04, vjust=0.2) 

