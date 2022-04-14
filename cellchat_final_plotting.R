### Code used for Cellchat analysis and plotting split by regions
#Starting with object 'clean_demux_full.rds' which contains all QC filtered, demuxed, immune cells
library(CellChat)
test_list = SplitObject(ln,split.by = 'region')
#DimPlot(test_list[[2]])
#DimPlot(test_list[[3]])
cellchat_outer = createCellChat(object = test_list[[3]],group.by = 'ident',assay = 'RNA')
cellchat_inner = createCellChat(object = test_list[[2]],group.by = 'ident',assay = 'RNA')
cellchat_mid = createCellChat(object = test_list[[1]],group.by = 'ident',assay = 'RNA')
DB = CellChatDB.mouse

cellchat_inner@DB <- DB
cellchat_inner <- subsetData(cellchat_inner)
future::plan("multiprocess",workers=6)
cellchat_inner<- identifyOverExpressedGenes(cellchat_inner)
cellchat_inner <- identifyOverExpressedInteractions(cellchat_inner)
beepr::beep(2)

cellchat_outer@DB <- DB
cellchat_outer <- subsetData(cellchat_outer)
future::plan("multiprocess",workers=6)
cellchat_outer<- identifyOverExpressedGenes(cellchat_outer)
cellchat_outer <- identifyOverExpressedInteractions(cellchat_outer)
beepr::beep(2)

cellchat_mid@DB <- DB
cellchat_mid <- subsetData(cellchat_mid)
future::plan("multiprocess",workers=6)
cellchat_mid<- identifyOverExpressedGenes(cellchat_mid)
cellchat_mid <- identifyOverExpressedInteractions(cellchat_mid)
beepr::beep(2)


#cellchat_inner<-computeCommunProb(cellchat_inner,type = "triMean",raw.use = TRUE)
cellchat_inner<-computeCommunProb(cellchat_inner,type = "triMean",raw.use = TRUE,population.size = TRUE)
cellchat_inner <- computeCommunProbPathway(cellchat_inner)
cellchat_inner <- aggregateNet(cellchat_inner)
cellchat_inner <- netAnalysis_computeCentrality(cellchat_inner,slot.name = "netP")

#cellchat_outer<-computeCommunProb(cellchat_outer,type = "triMean",raw.use = TRUE)
cellchat_outer<-computeCommunProb(cellchat_outer,type = "triMean",raw.use = TRUE,population.size = TRUE)
cellchat_outer <- computeCommunProbPathway(cellchat_outer)
cellchat_outer <- aggregateNet(cellchat_outer)
cellchat_outer <- netAnalysis_computeCentrality(cellchat_outer,slot.name = "netP")

cellchat_mid<-computeCommunProb(cellchat_mid,type = "triMean",raw.use = TRUE,population.size = TRUE)
cellchat_mid <- computeCommunProbPathway(cellchat_mid)
cellchat_mid <- aggregateNet(cellchat_mid)
cellchat_mid <- netAnalysis_computeCentrality(cellchat_mid,slot.name = "netP")

object.list <- list(OUTER = cellchat_outer,INNER = cellchat_inner)
cellchat_beeg <- mergeCellChat(object.list, add.names = names(object.list))

#####Circle plot for CSF signaling pathway
pathways.show = 'CSF'
weight.max <- getMaxWeight(object.list,slot.name = c("netP"),attribute = pathways.show)
vertex.receiver = seq(1,2)
par(mfrow = (2,1),xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]],signaling = pathways.show, vertex.receiver = vertex.receiver, 
                      edge.weight.max = weight.max[[1]], edge.width.max = 10,signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("CSF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]),targets.use = c(1,2,3,4,5,6,8,9,13),show.legend = TRUE)
}

