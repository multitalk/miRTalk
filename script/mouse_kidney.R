library(miRTalk)
load("rawdata.rda")
samplename <- unique(sc_meta$Time.point)
sc_data <- rev_gene(data = sc_data, data_type = "count",species = "Mouse",geneinfo = geneinfo)

a2 <- sc_meta[sc_meta$Time.point == samplename[1],]
a1 <- sc_data[,a2$id]
obj <- create_miRTalk(sc_data = a1,sc_celltype = as.character(a2$Annotation.Level.2), species = "Mouse",if_normalize = F)
obj <- find_miRNA(obj, mir_info = mir_info)
obj <- find_hvtg(obj)
obj <- find_miRTalk(obj, mir2tar = mir2tar,use_n_cores = 4)
obj1 <- obj

a2 <- sc_meta[sc_meta$Time.point != samplename[1],]
a1 <- sc_data[,a2$id]
obj <- create_miRTalk(sc_data = a1,sc_celltype = as.character(a2$Annotation.Level.2), species = "Mouse",if_normalize = F)
obj <- find_miRNA(obj, mir_info = mir_info)
obj <- find_hvtg(obj)
obj <- find_miRTalk(obj, mir2tar = mir2tar,use_n_cores = 4)
obj2 <- obj

# plot

# tsne plot
library(Seurat)
obj_seurat <- CreateSeuratObject(counts = sc_data)
obj_seurat <- ScaleData(obj_seurat)
obj_seurat <- FindVariableFeatures(obj_seurat)
obj_seurat <- RunPCA(obj_seurat)
obj_seurat <- RunTSNE(object = obj_seurat,dims = 1:30)
obj_seurat <- RunUMAP(object = obj_seurat,dims = 1:30)

sc_meta$condition <- sc_meta$Time.point
sc_meta[sc_meta$condition != "Uninjured",]$condition <- "UUO"
obj_seurat$condition <- sc_meta$condition
p1 <- DimPlot(object = obj_seurat,reduction = "tsne",group.by = "condition",pt.size = 2)+NoLegend()

# chord plot
plot_miRTalk_chord(obj1)
plot_miRTalk_chord(obj2)
plot_miR_heatmap(obj2,grid_color = "black")
plot_miR2tar_circle(object = obj2,celltype_sender = "Injured_Vascular_Smooth_Muscle_Cells",celltype_receiver = "Myofibroblast",celltype_color = "NO")

emt_score <- gsva_ch["hallmark_epithelial_mesenchymal_transition", ]
obj_seurat$emt <- as.numeric(emt_score)
library(viridis)
FeaturePlot(obj_seurat,features = "emt",reduction = "tsne",pt.size = 3)+scale_color_viridis(option = "C")

y1 <- gsva_c2["reactome_collagen_formation",]
y2 <- gsva_c2["reactome_collagen_biosynthesis_and_modifying_enzymes",]
y3 <- gsva_c5["gobp_negative_regulation_of_epithelial_to_mesenchymal_transition",]
x1 <- sc_data_rev["Mef2c",]
x2 <- sc_data_rev["Rcan2",]
x3 <- sc_data_rev["Npy1r",]
x4 <- sc_data_rev["Rgs7bp",]
x4 <- sc_data_rev["Gprc5c",]

# Mef2c
d1 <- data.frame(x = x1, y = y1, celltype = sc_meta$Annotation.Level.2, stringsAsFactors = F)
d1[d1$celltype != "Myofibroblast",]$celltype <- "Other"
ggplot(d1)+geom_point(aes(x,y),size = 1)+geom_smooth(aes(x,y),method = "lm")+scale_color_manual(values = c("#1ab33b", "black"))
correlation::cor_test(d1,x = "x",y = "y")
summary(lm(d1$y~d1$x))
wilcox.test(d1[d1$celltype != "Other",]$y,d1[d1$celltype == "Other",]$y,alternative = "greater")
# Mef2c
d1 <- data.frame(x = x1, y = y2, celltype = sc_meta$Annotation.Level.2, stringsAsFactors = F)
d1[d1$celltype != "Myofibroblast",]$celltype <- "Other"
ggplot(d1)+geom_point(aes(x,y),size = 1)+geom_smooth(aes(x,y),method = "lm")+scale_color_manual(values = c("#1ab33b", "black"))
correlation::cor_test(d1,x = "x",y = "y")
summary(lm(d1$y~d1$x))
wilcox.test(d1[d1$celltype != "Other",]$y,d1[d1$celltype == "Other",]$y,alternative = "greater")
# Rcan2
d1 <- data.frame(x = x2, y = y1, celltype = sc_meta$Annotation.Level.2, stringsAsFactors = F)
d1[d1$celltype != "Myofibroblast",]$celltype <- "Other"
ggplot(d1)+geom_point(aes(x,y),size = 1)+geom_smooth(aes(x,y),method = "lm")+scale_color_manual(values = c("#1ab33b", "black"))
correlation::cor_test(d1,x = "x",y = "y")
summary(lm(d1$y~d1$x))
# Npy1r
d1 <- data.frame(x = x3, y = y1, celltype = sc_meta$Annotation.Level.2, stringsAsFactors = F)
d1[d1$celltype != "Myofibroblast",]$celltype <- "Other"
ggplot(d1)+geom_point(aes(x,y),size = 1)+geom_smooth(aes(x,y),method = "lm")+scale_color_manual(values = c("#1ab33b", "black"))
correlation::cor_test(d1,x = "x",y = "y")
summary(lm(d1$y~d1$x))
# Rgs7bp
d1 <- data.frame(x = x4, y = y1, celltype = sc_meta$Annotation.Level.2, stringsAsFactors = F)
d1[d1$celltype != "Myofibroblast",]$celltype <- "Other"
ggplot(d1)+geom_point(aes(x,y),size = 1)+geom_smooth(aes(x,y),method = "lm")+scale_color_manual(values = c("#1ab33b", "black"))
correlation::cor_test(d1,x = "x",y = "y")
summary(lm(d1$y~d1$x))
# Gprc5c
d1 <- data.frame(x = x5, y = y1, celltype = sc_meta$Annotation.Level.2, stringsAsFactors = F)
d1[d1$celltype != "Myofibroblast",]$celltype <- "Other"
ggplot(d1)+geom_point(aes(x,y),size = 1)+geom_smooth(aes(x,y),method = "lm")+scale_color_manual(values = c("#1ab33b", "black"))
correlation::cor_test(d1,x = "x",y = "y")
summary(lm(d1$y~d1$x))


obj_seurat$score <- as.numeric(y1)
p1 <- FeaturePlot(obj_seurat,features = "score",reduction = "tsne",pt.size = 3)+scale_color_viridis(option = "D")+NoLegend()
obj_seurat$score <- as.numeric(y2)
p1 <- FeaturePlot(obj_seurat,features = "score",reduction = "tsne",pt.size = 3)+scale_color_viridis(option = "D")+NoLegend()
plot_miR2tar_heatmap(object = obj2,celltype_sender = "Injured_Vascular_Smooth_Muscle_Cells",celltype_receiver = "Parietal_Epithelial_Cells")
plot_miRTalk_circle(obj2,show_type = "score")
plot_miRTalk_circle_simple(obj2,celltype = "Injured_Vascular_Smooth_Muscle_Cells",miRNA = "mmu-miR-27a-3p",if_show_autocrine = T,show_type = "score")
obj_seurat$score <- as.numeric(y3)
p1 <- FeaturePlot(obj_seurat,features = "score",reduction = "tsne",pt.size = 4,order = T)+scale_color_viridis(option = "D")+NoLegend()
############### bubble
cci_sim <- get_miRTalk_cci(obj2)
cci1 <- cci_sim[cci_sim$celltype_sender == "Injured_Vascular_Smooth_Muscle_Cells",]
cci1 <- cci1[cci1$miRNA == "mmu-miR-27a-3p",]

y_name <-  unique(cci1$celltype_receiver)
x_name <- unique(cci1$target_gene)

y_name <- data.frame(names = y_name, y = 6:1,stringsAsFactors = F)
x_name <- data.frame(names = x_name, x = 1:length(x_name),stringsAsFactors = F)

cci1$x <- 0
cci1$y <- 0

for (i in 1:nrow(cci1)) {
  d1 <- cci1$celltype_receiver[i]
  d1 <- y_name[y_name$names == d1,]$y
  cci1$y[i] <- d1

  d1 <- cci1$target_gene[i]
  d1 <- x_name[x_name$names == d1,]$x
  cci1$x[i] <- d1
}

cci12 <- cci1[1,]
cci12$x <- 1
cci12$y <- 1
cci12$score <- 0
cci12$prob <- min(cci1$prob)
cci1 <- rbind(cci1,cci12)
p1 <- ggplot2::ggplot(data = cci1) + ggplot2::geom_point(ggplot2::aes(x = x, y = y, color = prob,
                                                                      size = score)) + ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
                                                                                                                            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))+scale_color_viridis(option = "C")+
  ggplot2::scale_x_continuous(name = "targets", breaks = x_name$x, labels = x_name$names, limits = c(1, max(x_name$x)))+
  scale_size_continuous(range = c(2,10))
p1
ggsave(filename = "miR-27a-3p_bubble.pdf",device = "pdf",plot = p1,width = 7,height = 7)

library(Nebulosa)
plot_density(object = obj_seurat,features = 'Mir27a',reduction = 'tsne')
plot_density(object = obj_seurat,features = 'Gata3',reduction = 'tsne')
plot_density(object = obj_seurat,features = 'Epha4',reduction = 'tsne')

Idents(obj_seurat) <- factor(sc_meta$Annotation.Level.2, levels = c("Injured Vascular Smooth Muscle Cells",
    "Mesangial Cells", "Myofibroblast", "Parietal Epithelial Cells", "Pericytes", "Vascular Smooth Muscle Cells"))
VlnPlot(obj_seurat,features = "Mef2c",pt.size = 0)
VlnPlot(obj_seurat,features = "Rcan2",pt.size = 0)
VlnPlot(obj_seurat,features = "Npy1r",pt.size = 0)
# Gata3 Rrad Mef2c
VlnPlot(obj_seurat_subset,features = "Gata3",pt.size = 0)
VlnPlot(obj_seurat_subset,features = "Rrad",pt.size = 0)
VlnPlot(obj_seurat_subset,features = "Mef2c",pt.size = 0)
# monocle3
obj_cds <- as.cell_data_set(obj_seurat)
obj_cds <- cluster_cells(cds = obj_cds, reduction_method = "UMAP")
sc_meta$cluster <- 1
sc_meta[sc_meta$Annotation.Level.2 == "Injured Vascular Smooth Muscle Cells",]$cluster <- 2
sc_meta[sc_meta$Annotation.Level.2 == "Mesangial Cells",]$cluster <- 3
sc_meta[sc_meta$Annotation.Level.2 == "Parietal Epithelial Cells",]$cluster <- 4
sc_meta[sc_meta$Annotation.Level.2 == "Pericytes",]$cluster <- 5
sc_meta[sc_meta$Annotation.Level.2 == "Myofibroblast",]$cluster <- 6
sc_meta_clu <- as.integer(sc_meta$cluster)
names(sc_meta_clu) <- sc_meta$id
obj_cds@clusters@listData[["UMAP"]][["clusters"]] <- factor(sc_meta_clu)
# obj_cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- obj_cds@int_colData@listData[["reducedDims"]]@listData[["TSNE"]]
obj_cds <- learn_graph(obj_cds, use_partition = FALSE)
obj_cds <- order_cells(obj_cds)
plot_cells(
  cds = obj_cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,cell_size = 1
)
