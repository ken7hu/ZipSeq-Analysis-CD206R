ln<-subset(ln,idents = c(13),invert = TRUE)
ln<- SCTransform(ln, vars.to.regress = "percent.mt", verbose = TRUE)
ln <- RunPCA(ln,verbose = FALSE)
ln <- FindNeighbors(ln, dims = 1:17, verbose = FALSE)
ln <- FindClusters(ln, verbose = FALSE,resolution = 0.6)
ElbowPlot(ln,reduction = "pca",ndims = 40)
ln <- ln %>% RunUMAP( dims = 1:17, verbose = FALSE,n.neighbors = 30,min.dist = 0.2)
DimPlot(pbmc_sub,group.by = 'ADT_maxID')
DefaultAssay(pbmc)<-'RNA'
DimPlot(ln)

pbmc.markers <- FindAllMarkers(ln, only.pos = TRUE, min.pct = 0.15, logfc.threshold = 0.2)
#pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
#top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
b_plasma_harmony_0.2_Markers_padj0.1_Top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DotPlot(ln,
        features = unique(b_plasma_harmony_0.2_Markers_padj0.1_Top5$gene),
        cols = "RdBu", assay = "RNA") + theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()

ln<-subset(ln,idents = c(9),invert = TRUE)
# 0: apoe
# 1: prolif macs
# 2: ms4a7 
# 3: inflamm mono
# 4: ly6c lo mono
#5: nos2 macc
#6: cd4 t cells
#7:ly6c lo mono
#8: retnla mac
#9: dying
#10: naive cd8
#11: exhausted cd8
#12: NK
#13:DC
new.cluster.ids <- c("Apoe Mac.", "Prolif. Mac.", "Ms4a7 Mac.", "Inflamm. Mono.", "Interm. Mac.", "Nos2 Mac.", 
                     "CD4 T", "MoDC", "Retnla Mac.","Tcf7+ CD8 T","Exhaust CD8 T",'NK','DC')
names(new.cluster.ids) <- levels(ln)
ln <- RenameIdents(ln, new.cluster.ids)
DimPlot(ln, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(magrittr)
#Starting with object 'clean_demux_full.rds' which contains all QC filtered, demuxed, immune cells
# Cells 
ln = readRDS("D:/Krummel Lab/SC_Experiments/20210112_CD206R/for_kelly/final_figs/20210527_full_clean_mito_removed.rds")
#############Plot whole object (Figure 6B):####################
b <- DimPlot(ln,pt.size = 1.2)
b+ylab('UMAP 2')+xlab('UMAP 1')+theme(legend.text = element_text(size=20),axis.title = element_text(size=24))
#### For Bubble plot of DE for whole object: (Supp Fig )#########
DefaultAssay(ln)<-'RNA'
markers <- FindAllMarkers(ln, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.20)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DotPlot(ln,
        features = unique(top10$gene),
        cols = "RdBu", assay = "RNA") + theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()
#####Nebulosa plotting (Supp Fig 7a)##########
library(Nebulosa)
p <- plot_density(ln, c("Cd4","Cd8b1","Tcf7","Pdcd1","Ly6c2","Apoe","Nos2","Ms4a7","Top2a"))
p + plot_layout(ncol = 3)
###get the original colors
library(scales)
my_color_palette <- hue_pal()(length(levels(Idents(ln))))
my_color_palette <- hue_pal()(13)

########################
#for plotting umap of monomacs
#########################
mono.macs = subset(ln,idents = c('Apoe Mac.','Prolif. Mac.','Ms4a7 Mac.','Inflamm. Mono.','Interm. Mac.','Nos2 Mac.','Retnla Mac.'))
DefaultAssay(mono.macs) = 'SCT'
mono.macs <- RunPCA(mono.macs,verbose = FALSE)
ElbowPlot(mono.macs,reduction = "pca",ndims = 40)
mono.macs <- mono.macs %>% RunUMAP( dims = 1:17, verbose = FALSE,n.neighbors = 50,min.dist = 0.2,seed.use = 12345)
###retrieve this monomac object:
mac_sub = mono.macs
m1 <- DimPlot(mac_sub,label = TRUE,label.size = 6,pt.size = 1.5,cols = my_color_palette[c(1:6,9)])
m1 = m1+theme(legend.text = element_text(size=20),axis.title = element_blank())
m2 <- DimPlot(mac_sub,pt.size = 1.5,group.by = 'region')
m2= m2+theme(legend.text = element_text(size=20),axis.title = element_blank())
####For the CD8 T subset
library(patchwork)
cd8.t = subset(ln,idents = c('Tcf7+ CD8 T','Exhaust CD8 T'))
cd8.t$region = factor(cd8.t$region,levels = c('Outer','Mid','Inner'))
DefaultAssay(cd8.t) = 'SCT'
cd8.t <- RunPCA(cd8.t,verbose = FALSE)
ElbowPlot(cd8.t,reduction = "pca",ndims = 40)
cd8.t <- cd8.t %>% RunUMAP( dims = 1:8, verbose = FALSE,n.neighbors = 20,min.dist = 0.2,seed.use = 12345)
################################
### For plotting Figure 6D
m1 <- DimPlot(cd8.t,label = TRUE,label.size = 6,pt.size = 2,cols = c(my_color_palette[10:11]))
m1 = m1+theme(legend.text = element_text(size=20),axis.title = element_blank(),legend.position = 'None')+scale_x_continuous(limits = c(-8,8),expand = c(0, 0)) +
  scale_y_continuous(limits = c(-3,3),expand = c(0, 0))
#m2 <- DimPlot(cd8.t,pt.size = 3,group.by = 'region')
#m2= m2+theme(legend.text = element_text(size=20),axis.title = element_blank())
#m1 + m2 + plot_layout(ncol = 1)
##############playing with contour plot version:
#need to extract umap coordinates
coords = cd8.t@reductions$umap@cell.embeddings
coords = as.data.frame(coords)
coords$region = cd8.t$region
ghost = rgb(0,0,0,0.4)
#For generating figure 6D
m2 = ggplot(data = coords,aes(x=UMAP_1,y=UMAP_2))+geom_density_2d_filled(data = subset(coords,subset = region == 'Outer'),alpha=0.4,bins=9)+scale_fill_manual(values = brewer.pal(n=9,name = "Reds"))+
  geom_point(data = coords,aes(color = region),size=2)+
  scale_color_manual(values = c('red',ghost,ghost))
m2 = m2+theme_classic()+theme(legend.position = "None",legend.text = element_text(size=20),axis.title = element_blank())+scale_x_continuous(limits = c(-8,8),expand = c(0, 0)) +
  scale_y_continuous(limits = c(-3,3),expand = c(0, 0))
m3 = ggplot(data = coords,aes(x=UMAP_1,y=UMAP_2))+geom_density_2d_filled(data = subset(coords,subset = region == 'Mid'),alpha=0.5,bins=9)+scale_fill_manual(values = brewer.pal(n=9,name = "Greens"))+
  geom_point(data = coords,aes(color = region),size=2)+
  scale_color_manual(values = c(ghost,'green',ghost))
m3 = m3+theme_classic()+theme(legend.position = "None",legend.text = element_text(size=20),axis.title = element_blank())+scale_x_continuous(limits = c(-8,8),expand = c(0, 0)) +
  scale_y_continuous(limits = c(-3,3),expand = c(0, 0))
m4 = ggplot(data = coords,aes(x=UMAP_1,y=UMAP_2))+geom_density_2d_filled(data = subset(coords,subset = region == 'Inner'),alpha=0.5,bins=9)+scale_fill_manual(values = brewer.pal(n=9,name = "Blues"))+
  geom_point(data = coords,aes(color = region),size=2)+
  scale_color_manual(values = c(ghost,ghost,'blue'))
m4 = m4+theme_classic()+theme(legend.position = "None",legend.text = element_text(size=20),axis.title = element_blank())+scale_x_continuous(limits = c(-8,8),expand = c(0, 0)) +
  scale_y_continuous(limits = c(-3,3),expand = c(0, 0))
m1 + m2 + m3+m4+plot_layout(ncol = 4)






pdf(file = "~/Krummel lab/20210112_CD206R/for_kelly/final_figs/cd8t_umaps_densities.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 4,useDingbats = FALSE) # The height of the plot in inches

m1 + m2 + m3+m4+plot_layout(ncol = 4)

dev.off()
###retrieve this cd8 t object:
cd8.t = readRDS("~/Krummel lab/20210112_CD206R/for_kelly/final_figs/20210628_cd8t_final.rds")
