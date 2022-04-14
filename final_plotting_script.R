library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(magrittr)
#Starting with object 'clean_demux_full.rds' which contains all QC filtered, demuxed, immune cells

#############Plot whole object (Figure 6B):####################
b <- DimPlot(ln,pt.size = 1.2)
b+ylab('UMAP 2')+xlab('UMAP 1')+theme(legend.text = element_text(size=20),axis.title = element_text(size=24))

######## Plot stacked barchart based on regional frequency (figure 6C) ######
freqs = interaction(Idents(ln),ln$region) %>% table() %>% as.data.frame()
freqs$region = c(rep('Inner',13),rep('Mid',13),rep('Outer',13))
freqs$cluster = rep(new.cluster.ids,3)
freqs$cluster = factor(freqs$cluster,levels = new.cluster.ids)
freqs$region = factor(freqs$region,levels = c('Outer','Mid','Inner'))
ggplot(freqs, aes(fill=region, y=Freq, x=cluster)) + 
  geom_bar(position="fill", stat="identity")+theme_classic()+
  theme(axis.text.x = element_text(size = 14,angle=45,hjust =1),axis.text.y = element_text(size=14),legend.text = element_text(size=14),legend.title = element_blank(),
        axis.title = element_blank())+scale_fill_manual(values = c('red','green','blue'))

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

####For the CD8 T subset #############
library(patchwork)
cd8.t = subset(ln,idents = c('Tcf7+ CD8 T','Exhaust CD8 T'))
cd8.t$region = factor(cd8.t$region,levels = c('Outer','Mid','Inner'))
DefaultAssay(cd8.t) = 'SCT'
cd8.t <- RunPCA(cd8.t,verbose = FALSE)
ElbowPlot(cd8.t,reduction = "pca",ndims = 40)
cd8.t <- cd8.t %>% RunUMAP( dims = 1:8, verbose = FALSE,n.neighbors = 20,min.dist = 0.2,seed.use = 12345)

###for the Mono/Mac subset (excluding the MoDC and DC) ##############
mac_sub<-subset(ln,idents = c('Inflamm. Mono.','Interm. Mac.','Retnla Mac.','Apoe Mac.','Nos2 Mac.','Ms4a7 Mac.','Prolif. Mac.'))
DefaultAssay(mac_sub)<-'SCT'
mac_sub <- RunPCA(mac_sub,verbose = FALSE)
ElbowPlot(mac_sub,reduction = "pca",ndims = 50)
mac_sub <- RunUMAP(mac_sub, dims = 1:18, min.dist = 0.15,  verbose = TRUE,assay = "SCT",seed.use = 12345L)
DimPlot(mac_sub,label=TRUE)
################################


###############For generating figure 6D
m1 <- DimPlot(cd8.t,label = TRUE,label.size = 6,pt.size = 2,cols = c(my_color_palette[10:11]))
m1 = m1+theme(legend.text = element_text(size=20),axis.title = element_blank(),legend.position = 'None')+scale_x_continuous(limits = c(-8,8),expand = c(0, 0)) +
  scale_y_continuous(limits = c(-3,3),expand = c(0, 0))
#m2 <- DimPlot(cd8.t,pt.size = 3,group.by = 'region')
#m2= m2+theme(legend.text = element_text(size=20),axis.title = element_blank())
#m1 + m2 + plot_layout(ncol = 1)
##### for splitting out UMAPs by region: 
#need to extract umap coordinates
coords = cd8.t@reductions$umap@cell.embeddings
coords = as.data.frame(coords)
coords$region = cd8.t$region
ghost = rgb(0,0,0,0.4)
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
m1 + m2 + m3 + m4 + plot_layout(ncol = 4)

############ For generating Figure 6E
m1 <- DimPlot(mac_sub,label = TRUE,label.size = 6,pt.size = 1.5,cols = my_color_palette[c(1:6,9)])
m1 = m1+theme(legend.text = element_text(size=20),axis.title = element_blank())
##### for splitting out UMAPs by region: 
#need to extract umap coordinates
coords = mac_sub@reductions$umap@cell.embeddings
coords = as.data.frame(coords)
coords$region = factor(mac_sub$region,levels = c('Outer','Mid','Inner'))
ghost = rgb(0,0,0,0.4)
m2 = ggplot(data = coords,aes(x=UMAP_1,y=UMAP_2))+geom_density_2d_filled(data = subset(coords,subset = region == 'Outer'),alpha=0.4,bins=9)+scale_fill_manual(values = brewer.pal(n=9,name = "Reds"))+
  geom_point(data = coords,aes(color = region),size=0.5)+
  scale_color_manual(values = c('red',ghost,ghost))
m2 = m2+theme_classic()+theme(legend.position = "None",legend.text = element_text(size=20),axis.title = element_blank())+scale_x_continuous(limits = c(-9,8),expand = c(0, 0)) +
  scale_y_continuous(limits = c(-4,4),expand = c(0, 0))
m3 = ggplot(data = coords,aes(x=UMAP_1,y=UMAP_2))+geom_density_2d_filled(data = subset(coords,subset = region == 'Mid'),alpha=0.5,bins=9)+scale_fill_manual(values = brewer.pal(n=9,name = "Greens"))+
  geom_point(data = coords,aes(color = region),size=0.5)+
  scale_color_manual(values = c(ghost,'green',ghost))
m3 = m3+theme_classic()+theme(legend.position = "None",legend.text = element_text(size=20),axis.title = element_blank())+scale_x_continuous(limits = c(-9,8),expand = c(0, 0)) +
  scale_y_continuous(limits = c(-4,4),expand = c(0, 0))
m4 = ggplot(data = coords,aes(x=UMAP_1,y=UMAP_2))+geom_density_2d_filled(data = subset(coords,subset = region == 'Inner'),alpha=0.5,bins=9)+scale_fill_manual(values = brewer.pal(n=9,name = "Blues"))+
  geom_point(data = coords,aes(color = region),size=0.5)+
  scale_color_manual(values = c(ghost,ghost,'blue'))
m4 = m4+theme_classic()+theme(legend.position = "None",legend.text = element_text(size=20),axis.title = element_blank())+scale_x_continuous(limits = c(-9,8),expand = c(0, 0)) +
  scale_y_continuous(limits = c(-4,4),expand = c(0, 0))

m1 + m2 + m3 + m4 + plot_layout(ncol = 4)

######## Signature score generation or the CD8 T object
sig = list(Glycolysis_Zenith_mouseB16$X1)
cd8.t <- AddModuleScore(object  = cd8.t,features = sig,ctrl = 50,name = c('glycolysis'),assay = 'RNA')
VlnPlot(cd8.t,features = 'glycolysis1',group.by = 'region')
sig = list(t_ex$X1)
cd8.t <- AddModuleScore(object  = cd8.t,features = sig,ctrl = 50,name = c('exhaustion'),assay = 'RNA')
VlnPlot(cd8.t,features = 'exhaustion1',group.by = 'region')
###add antigen presentation and glycolysis score to monomacs
GO_term_summary_antigen_presentation <- read_excel("C:/Users/KHu/Downloads/GO_term_summary_antigen_presentation.xlsx")
ag = GO_term_summary_antigen_presentation$Symbol %>% unique()
sig = list(ag)
mac_sub <- AddModuleScore(object  = mac_sub,features = sig,ctrl = 50,name = c('antigen_presentation'),assay = 'RNA')
sig = list(Glycolysis_Zenith_mouseB16$X1)
mac_sub <- AddModuleScore(object  = mac_sub,features = sig,ctrl = 50,name = c('glycolysis'),assay = 'RNA')
######### Monocle3 analysis on Mono/Mac subset ########
library(monocle3)
library(SeuratWrappers)
cds <- as.cell_data_set(mac_sub)
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
#integrated.sub <- subset(as.Seurat(cds), monocle3_partitions == 1)
#cds <- as.cell_data_set(integrated.sub)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
cds <- order_cells(cds)
######## generate Supp Fig 7B
pp = plot_cells(
  cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE, cell_size = 2,trajectory_graph_segment_size = 1.5
)
pp+theme(legend.title = element_text(size=20),legend.text = element_text(size=20),axis.title = element_blank(),axis.text = element_blank())

cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(mac_sub[["SCT"]])

##### extract pseudotime to put as metadata
sudo = cds@principal_graph_aux@listData$UMAP$pseudotime
mac_sub<-AddMetaData(mac_sub,metadata = sudo,col.name = 'pseudotime')
VlnPlot(mac_sub,features = c('pseudotime'),group.by = 'region')

#### FINAL PLOTTING:
###### for generating Supplementary Figure 7D
plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="region",
                         min_expr=0.5,label_by_short_name = FALSE)+scale_color_manual(values= c('blue','green','red'))
###### for generating Supplementary Figure 7C
plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="ident",
                         min_expr=0.5,label_by_short_name = FALSE)+scale_color_manual(values= my_color_palette[c(1:6,9)])


######## VlnPlots w lines
####For a signature/module score####################
####################################################
#########################   for generating figure 6F (the cd8 t violin plots)  ############
gene = 'glycolysis1'
region = c('Outer','Mid','Inner')
ln_sub = cd8.t
Idents(ln_sub)<-factor(ln_sub$region,levels = c('Outer','Mid','Inner'))
summary = data.frame(region)
summary$mean = c(mean(ln_sub[[gene]][ln_sub$region == 'Outer',1]),mean(ln_sub[[gene]][ln_sub$region == 'Mid',1]),mean(ln_sub[[gene]][ln_sub$region == 'Inner',1]))
summary$ident = summary$region
v1 = VlnPlot(ln_sub,features = c(gene))+  
  geom_path(data = summary,aes(x=region,y=mean,group=1),size=1.2,linetype = 2)+ggtitle('CD8 T Glycolysis')+
  geom_point(data = summary,aes(x=region,y=mean,group=1),shape=23,fill='white',color = 'black',size=2,stroke=2)+
  labs(y='Signature Score')+theme(axis.title.x = element_blank(),axis.text.x = element_text(size=20),legend.position = 'None',axis.text.y = element_text(size = 20),
                                  axis.title.y = element_text(size=20))+scale_fill_manual(values = c('red','green','blue'))
##########################  for generating figure 6F  ###############
gene = 'exhaustion1'
summary = data.frame(region)
summary$mean = c(mean(ln_sub[[gene]][ln_sub$region == 'Outer',1]),mean(ln_sub[[gene]][ln_sub$region == 'Mid',1]),mean(ln_sub[[gene]][ln_sub$region == 'Inner',1]))
summary$ident = summary$region
v2 = VlnPlot(ln_sub,features = c(gene))+  
  geom_path(data = summary,aes(x=region,y=mean,group=1),size=1.2,linetype = 2)+ggtitle('CD8 T Exhaustion')+
  geom_point(data = summary,aes(x=region,y=mean,group=1),shape=23,fill='white',color = 'black',size=2,stroke=2)+
  labs(y='Signature Score')+theme(axis.title.x = element_blank(),axis.text.x = element_text(size=20),legend.position = 'None',axis.text.y = element_text(size = 20),
                                  axis.title.y = element_text(size=20))+scale_fill_manual(values = c('red','green','blue'))

v1+v2+plot_layout(ncol=2)

###################### for generating bottom half of figure 6F (the monomac violin plots)
gene = 'pseudotime'
region = c('Outer','Mid','Inner')
ln_sub = mac_sub
Idents(ln_sub)<-factor(ln_sub$region,levels = c('Outer','Mid','Inner'))
summary = data.frame(region)
summary$mean = c(mean(ln_sub[[gene]][ln_sub$region == 'Outer',1]),mean(ln_sub[[gene]][ln_sub$region == 'Mid',1]),mean(ln_sub[[gene]][ln_sub$region == 'Inner',1]))
summary$ident = summary$region
v1 = VlnPlot(ln_sub,features = c(gene))+  
  geom_path(data = summary,aes(x=region,y=mean,group=1),size=1.2,linetype = 2)+ggtitle('CD8 T Glycolysis')+
  geom_point(data = summary,aes(x=region,y=mean,group=1),shape=23,fill='white',color = 'black',size=2,stroke=2)+
  labs(y='Signature Score')+theme(axis.title.x = element_blank(),axis.text.x = element_text(size=20),legend.position = 'None',axis.text.y = element_text(size = 20),
                                  axis.title.y = element_text(size=20))+scale_fill_manual(values = c('red','green','blue'))
##########################  for generating figure 6F  ###############
gene = 'glycolysis1'
region = c('Outer','Mid','Inner')
ln_sub = mac_sub
Idents(ln_sub)<-factor(ln_sub$region,levels = c('Outer','Mid','Inner'))
summary = data.frame(region)
summary$mean = c(mean(ln_sub[[gene]][ln_sub$region == 'Outer',1]),mean(ln_sub[[gene]][ln_sub$region == 'Mid',1]),mean(ln_sub[[gene]][ln_sub$region == 'Inner',1]))
summary$ident = summary$region
v2 = VlnPlot(ln_sub,features = c(gene))+  
  geom_path(data = summary,aes(x=region,y=mean,group=1),size=1.2,linetype = 2)+ggtitle('CD8 T Exhaustion')+
  geom_point(data = summary,aes(x=region,y=mean,group=1),shape=23,fill='white',color = 'black',size=2,stroke=2)+
  labs(y='Signature Score')+theme(axis.title.x = element_blank(),axis.text.x = element_text(size=20),legend.position = 'None',axis.text.y = element_text(size = 20),
                                  axis.title.y = element_text(size=20))+scale_fill_manual(values = c('red','green','blue'))

v1+v2+plot_layout(ncol=2)
######## scatterplots for figure 6I
####################SIGNATURE 1##########################\
lymph_sub$exhaust_signature = lymph_sub$naive_ex2
lymph_sub$naive_signature = lymph_sub$naive_ex2
lymph_sub$stemlike_signature = lymph_sub$stem_term1
lymph_sub$terminal_signature = lymph_sub$stem_term2
#################
samp_mean <- function(x, i) {
  mean(x[i])
}
#########################For gene signature score of interest in cd8 T cells ##############################################
var1 = 'exhaustion1'
ln_sub = lymph_sub
region = c('Outer','Mid','Inner')
summary = data.frame(region)
summary$region = factor(summary$region,levels = c('Outer','Mid','Inner'))
summary$mean = c(mean(ln_sub[[var1]][ln_sub$region == 'Outer',1]),mean(ln_sub[[var1]][ln_sub$region == 'Mid',1]),mean(ln_sub[[var1]][ln_sub$region == 'Inner',1]))
#######calculate bootstrap 95% CI#############
results <- boot(data=ln_sub[[var1]][lymph_sub$region=='Outer',1], statistic=samp_mean,R=10000)
q = quantile(results$t,c(0.025,0.975))
results <- boot(data=ln_sub[[var1]][lymph_sub$region=='Mid',1], statistic=samp_mean,R=10000)
q = rbind(q,quantile(results$t,c(0.025,0.975))) 
results <- boot(data=ln_sub[[var1]][lymph_sub$region=='Inner',1], statistic=samp_mean,R=10000)
q = rbind(q,quantile(results$t,c(0.025,0.975))) 
q = as.data.frame(q)
summary = cbind(summary,q)
#########################For gene of interest in cd8 T cells ##############################################
var1 = 'Csf1'
ln_sub = lymph_sub
region = c('Outer','Mid','Inner')
summary = data.frame(region)
summary$region = factor(summary$region,levels = c('Outer','Mid','Inner'))
gene_expression = GetAssayData(ln_sub,assay = 'RNA',slot = 'scale.data')
gene_expression = gene_expression[row.names(gene_expression)==var1,]
summary$mean = c(mean(gene_expression[ln_sub$region == 'Outer']),mean(gene_expression[ln_sub$region == 'Mid']),mean(gene_expression[ln_sub$region == 'Inner']))

results <- boot(data=gene_expression[ln_sub$region=='Outer'], statistic=samp_mean,R=10000)
q = quantile(results$t,c(0.025,0.975))
results <- boot(data=gene_expression[ln_sub$region=='Mid'], statistic=samp_mean,R=10000)
q = rbind(q,quantile(results$t,c(0.025,0.975))) 
results <- boot(data=gene_expression[ln_sub$region=='Inner'], statistic=samp_mean,R=10000)
q = rbind(q,quantile(results$t,c(0.025,0.975))) 
q = as.data.frame(q)
summary = cbind(summary,q)
#########################For gene signature score of interest in mono/macs ##############################################
var2 = 'antigen_presentaion1'
ln_sub = mac_sub
region = c('Outer','Mid','Inner')
summary2 = data.frame(region)
summary2$region = factor(summary2$region,levels = c('Outer','Mid','Inner'))
summary2$mean = c(mean(ln_sub[[var2]][ln_sub$region == 'Outer',1]),mean(ln_sub[[var2]][ln_sub$region == 'Mid',1]),mean(ln_sub[[var2]][ln_sub$region == 'Inner',1]))

results <- boot(data=ln_sub[[var2]][ln_sub$region=='Outer',1], statistic=samp_mean,R=10000)
q = quantile(results$t,c(0.025,0.975))
results <- boot(data=ln_sub[[var2]][ln_sub$region=='Mid',1], statistic=samp_mean,R=10000)
q = rbind(q,quantile(results$t,c(0.025,0.975))) 
results <- boot(data=ln_sub[[var2]][ln_sub$region=='Inner',1], statistic=samp_mean,R=10000)
q = rbind(q,quantile(results$t,c(0.025,0.975))) 
q = as.data.frame(q)
colnames(q)<-c('low1','high2')
summary2 = cbind(summary2,q)
##########################For gene  of interest in mono/macs###############################
var2 = 'Ms4a7'
ln_sub = mac_sub
region = c('Outer','Mid','Inner')
summary2 = data.frame(region)
summary2$region = factor(summary2$region,levels = c('Outer','Mid','Inner'))
gene_expression = GetAssayData(ln_sub,assay = 'RNA',slot = 'scale.data')
gene_expression = gene_expression[row.names(gene_expression)==var2,]
summary2$mean = c(mean(gene_expression[ln_sub$region == 'Outer']),mean(gene_expression[ln_sub$region == 'Mid']),mean(gene_expression[ln_sub$region == 'Inner']))

results <- boot(data=gene_expression[ln_sub$region=='Outer'], statistic=samp_mean,R=10000)
q = quantile(results$t,c(0.025,0.975))
results <- boot(data=gene_expression[ln_sub$region=='Mid'], statistic=samp_mean,R=10000)
q = rbind(q,quantile(results$t,c(0.025,0.975))) 
results <- boot(data=gene_expression[ln_sub$region=='Inner'], statistic=samp_mean,R=10000)
q = rbind(q,quantile(results$t,c(0.025,0.975))) 
q = as.data.frame(q)
summary2 = cbind(summary2,q)

#Combine summary 1(from cd8 t) and summary 2 (from mono/mac) 
##########
summary_full = cbind(summary,summary2)
colnames(summary_full)<-c('region1','mean1','low1','high1','region2','mean2','low2','high2')
#######
#calculate limits of the plot:
xrange = max(summary_full$high1)-min(summary_full$low1)
yrange = max(summary_full$high2)-min(summary_full$low2)

g1<-ggplot(data = summary_full,aes(x=mean1,y=mean2,color = region1))+geom_point(size=3)+
  geom_errorbar(aes(ymin = low2,ymax = high2),size=1,width = xrange*0.05)+geom_errorbarh(aes(xmin =low1,xmax = high1),size=1,height = yrange*0.05 )+
  theme_bw()+theme(axis.text = element_text(size=18),axis.title = element_text(size=20),
                   legend.title = element_blank(),legend.text = element_text(size=14))+xlab(var1)+ylab(var2)

g4<-ggplot(data = summary_full,aes(x=mean1,y=mean2,color = region1))+geom_point(size=3)+
  geom_errorbar(aes(ymin = low2,ymax = high2),size=1,width = xrange*0.05)+geom_errorbarh(aes(xmin =low1,xmax = high1),size=1,height = yrange*0.05 )+
  theme_bw()+theme(axis.text = element_text(size=10),axis.title = element_text(size=16),
                   legend.title = element_blank(),legend.text = element_text(size=14))+xlab(var1)+ylab(var2)+ylim(0,12)

ggplot(data = summary_full,aes(x=mean1,y=mean2,color = region1))+geom_point(size=3)+
  geom_errorbar(aes(ymin = low2,ymax = high2),size=1,width = xrange*0.05)+geom_errorbarh(aes(xmin =low1,xmax = high1),size=1,height = yrange*0.05 )+
  theme_bw()+theme(axis.text = element_text(size=16),axis.title = element_text(size=16),
                   legend.title = element_blank(),legend.text = element_text(size=14))+xlab(var1)+ylab(var2)
g15<-g15+theme(axis.text = element_text(size=10))
gbig = arrangeGrob(g1,g2,g3,g4,g5,g6,ncol = 3)
gbig = arrangeGrob(g10,g11,g12,g13,g14,g15,ncol = 3)
gbig = arrangeGrob(g16,g17,g18,ncol = 3)

g1+g2+g3+g4+g5+g6+plot_layout(ncol = 3)
ggsave(gbig, file = "~/Krummel lab/20210112_CD206R/for_kelly/tam_vs_exhaust.pdf",height= 2.6,width=12)
ggsave(g18, file = "~/Krummel lab/20210112_CD206R/for_kelly/m1_v_exhaust.pdf",height= 2.6,width=4)

ggsave(gbig, file = "~/Krummel lab/20210112_CD206R/for_kelly/2x3grid.pdf",height= 6,width=12)
gbig+theme(axis.title.x = element_text(size=10))



