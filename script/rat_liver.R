library(miRTalk)
obj_con <- readRDS("obj_con.rds")
obj_fat <- readRDS("obj_fat.rds")
obj_condition <- list()

a2 <- obj_con@meta
a1 <- obj_con@data$data
obj <- create_miRTalk(sc_data = a1,sc_celltype = a2$celltype,
                      species = "Rat",condition = rep("control", ncol(a1)), evbiog = evbiog, risc = risc, if_normalize = F)
obj <- find_hvtg(object = obj)
obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
obj <- find_miRTalk(obj, if_filter_miRNA = F, use_n_cores = 16)
obj_condition[[1]] <- obj

a2 <- obj_fat@meta
a1 <- obj_fat@data$data
obj <- create_miRTalk(sc_data = a1,sc_celltype = a2$celltype,
                      species = "Rat",condition = rep("fat", ncol(a1)), evbiog = evbiog, risc = risc, if_normalize = F)
obj <- find_hvtg(object = obj)
obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
obj <- find_miRTalk(obj, if_filter_miRNA = F, use_n_cores = 16)
obj_condition[[2]] <- obj
obj_con <- obj_condition[[1]]
obj_fat <- obj_condition[[2]]
obj_con_data <- obj_con@data$data
obj_fat_data <- obj_fat@data$data

obj_con_meta <- obj_con@meta
obj_fat_meta <- obj_fat@meta

plot_miRTalk_circle(object = obj_con)
plot_miRTalk_circle(object = obj_fat)

plot_target_heatmap(object = obj_fat, celltype = "hepatocyte")

obj_con_meta$mir223 <- obj_con_data["Mir223",]
obj_fat_meta$mir223 <- obj_fat_data["Mir223",]
#obj_data1 <- obj_seurat[["RNA"]]@data
#obj_meta$condition <- as.character(obj_meta$Group)
#obj_meta$celltype <- as.character(obj_meta$celltypenew)
#obj_meta[obj_meta$celltype %in% c("pDC", "cDC"),]$celltype <- "DC"
#obj_seurat$celltype <- obj_meta$celltype
#DimPlot(obj_seurat,reduction = "tsne", group.by = "celltype",label = T,label.box = T)+NoLegend()

cci_cdl <- obj_con@cci
cci_fdl <- obj_fat@cci
mirnames <- unique(c(cci_cdl$miR_gene,cci_fdl$miR_gene))
mir_change <- data.frame(mirgene = mirnames, stringsAsFactors = F)
mir_change$exp_control <- 0
mir_change$exp_fat <- 0
mir_change$wilcol_p <- 1
celltypes1 <- unique(c(obj_con_meta$celltype, obj_fat_meta$celltype))
celltypes1 <- celltypes1[order(celltypes1)]
celltypes2 <- unique(c(cci_cdl$celltype_sender,cci_cdl$celltype_receiver,cci_fdl$celltype_sender,cci_fdl$celltype_receiver))
celltypes2 <- celltypes2[order(celltypes2)]
celltypes <- data.frame(celltype_seurat = celltypes1, celltype_mirtalk = celltypes2, stringsAsFactors = F)
res_all <- data.frame()
for (k in 1:nrow(celltypes)) {
  mir_change1 <- mir_change
  celltype_seurat <- celltypes$celltype_seurat[k]
  celltype_mirtalk <- celltypes$celltype_mirtalk[k]
  for (i in 1:nrow(mir_change1)) {
    obj_con_meta$gene <- as.numeric(as.numeric(obj_con_data[mir_change1$mirgene[i],]))
    obj_fat_meta$gene <- as.numeric(as.numeric(obj_fat_data[mir_change1$mirgene[i],]))
    obj_meta1 <- obj_con_meta[obj_con_meta$celltype == celltype_seurat,]
    obj_meta2 <- obj_fat_meta[obj_fat_meta$celltype == celltype_seurat,]
    mir_change1$exp_control[i] <- mean(obj_meta1$gene)
    mir_change1$exp_fat[i] <- mean(obj_meta2$gene)
    if (mean(obj_meta1$gene) < mean(obj_meta2$gene)) {
      mir_change1$wilcol_p[i] <- wilcox.test(obj_meta1$gene, obj_meta2$gene, alternative = "less")$p.value
    }
    if (mean(obj_meta1$gene) > mean(obj_meta2$gene)) {
      mir_change1$wilcol_p[i] <- wilcox.test(obj_meta1$gene, obj_meta2$gene, alternative = "greater")$p.value
    }
  }
  mir_change1 <- mir_change1[mir_change1$wilcol_p < 0.05,]
  if (nrow(mir_change1) > 0) {
    mir_change1$celltype <- celltype_mirtalk
    res_all <- rbind(res_all, mir_change1)
  }
}

miRnames <- unique(cci_fdl$miRNA)
miRnames <- miRnames[order(miRnames)]
celltypes <- unique(cci_fdl$celltype_sender)
celltypes <- celltypes[order(celltypes)]

res_plot <- as.data.frame(matrix(0, nrow = length(celltypes), ncol = length(miRnames)))
colnames(res_plot) <- miRnames
rownames(res_plot) <- celltypes

res_plot1 <- res_plot
for (i in 1:nrow(res_plot1)) {
  celltypename <- rownames(res_plot1)[i]
  cci_cdl1 <- cci_cdl[cci_cdl$celltype_sender == celltypename, ]
  for (j in 1:ncol(res_plot1)) {
    miRname <- colnames(res_plot1)[j]
    cci_cdl2 <- cci_cdl1[cci_cdl1$miRNA == miRname,]
    if (nrow(cci_cdl2) > 0) {
      cci_cdl2 <- unique(cci_cdl2$EVmiR_score)
      res_plot1[i,j] <- cci_cdl2
    }
  }
}
res_plot_n <- res_plot1

res_plot1 <- res_plot
for (i in 1:nrow(res_plot1)) {
  celltypename <- rownames(res_plot1)[i]
  cci_fdl1 <- cci_fdl[cci_fdl$celltype_sender == celltypename, ]
  for (j in 1:ncol(res_plot1)) {
    miRname <- colnames(res_plot1)[j]
    cci_fdl2 <- cci_fdl1[cci_fdl1$miRNA == miRname,]
    if (nrow(cci_fdl2) > 0) {
      cci_fdl2 <- unique(cci_fdl2$EVmiR_score)
      res_plot1[i,j] <- cci_fdl2
    }
  }
}
res_plot_f <- res_plot1

rownames(res_plot_n) <- paste0("normal", rownames(res_plot_n))
rownames(res_plot_f) <- paste0("fat", rownames(res_plot_f))
res_plot <- rbind(res_plot_n, res_plot_f)

heat_col <- viridis::viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "D")

heatmaply::heatmaply(x = as.matrix(res_plot), colors = heat_col, dendrogram = "none",margins = c(60,100,40,20),scale = "row",
                     titleX = FALSE, branches_lwd = 0.1,heatmap_layers = theme(axis.line=element_blank()), label_names = c("Sender","miRNA","EVmiR_score"))

res_all1 <- res_all[res_all$wilcol_p < 0.05,]
cci_fdl1 <- cci_fdl[cci_fdl$miR_gene %in% res_all1$mirgene,]
miRnames <- unique(cci_fdl1$miRNA)

exp_bulk <- data.frame()

load("GSE56071_mirna.rda")
exp_data <- as.data.frame(exp_data)
exp_data$gene <- exp_gene$miRNA_ID
exp_data <- exp_data[which(!is.na(exp_data$GSM1354615)),]
exp_data <- exp_data[which(!is.na(exp_data$GSM1354616)),]
exp_data <- exp_data[which(!is.na(exp_data$GSM1354617)),]
exp_data1 <- exp_data[grep("rno-miR-24", exp_data$gene),]
exp_data1 <- exp_data1["3637",]
exp_bulk <- rbind(exp_bulk, exp_data1)

exp_data1 <- exp_data[grep("rno-miR-21", exp_data$gene),]
exp_data1 <- exp_data1["4893",]
exp_bulk <- rbind(exp_bulk, exp_data1)

exp_data1 <- exp_data[grep("rno-miR-223", exp_data$gene),]
exp_data1 <- exp_data1["9236",]
exp_bulk <- rbind(exp_bulk, exp_data1)

exp_data1 <- exp_data[grep("rno-miR-290", exp_data$gene),]
exp_bulk <- rbind(exp_bulk, exp_data1)

res_all1 <- res_all1[res_all1$celltype != "B_cell",]
res_all1 <- res_all1[res_all1$celltype != "natural_killer_cell",]
res_all1 <- res_all1[res_all1$celltype != "hepatocyte",]
res_all1 <- res_all1[order(res_all1$mirgene),]
res_all1 <- res_all1[-c(4,5),]


for (i in 1:nrow(res_all1)) {
  obj_con_meta$gene <- as.numeric(obj_con_data[res_all1$mirgene[i],])
  obj_fat_meta$gene <- as.numeric(obj_fat_data[res_all1$mirgene[i],])
  x_normal <- obj_con_meta[obj_con_meta$celltype == res_all1$celltype[i],]$gene
  x_tumor <- obj_fat_meta[obj_fat_meta$celltype == res_all1$celltype[i],]$gene
  data_plot <- data.frame(x = c("normal", "tumor"), y = c(mean(x_normal), mean(x_tumor)),
                          ymin = c(mean(x_normal) - sd(x_normal)/sqrt(length(x_normal)), mean(x_tumor) - sd(x_tumor)/sqrt(length(x_tumor))),
                          ymax = c(mean(x_normal) + sd(x_normal)/sqrt(length(x_normal)), mean(x_tumor) + sd(x_tumor)/sqrt(length(x_tumor))),stringsAsFactors = F)
  p1 <- ggplot(data_plot)+geom_bar(aes(x = x, y=y),stat = "identity")+geom_errorbar(aes(x = x, ymin = ymin, ymax = ymax))
  ggsave(filename = paste0("difference_sc/",res_all1$celltype[i],"_",res_all1$mirgene[i], "_change.png"),plot = p1,device = "png",width = 6,height = 6)
}

res_path <- get_miRTalk_pathway(object = obj_fat,gene2path = gene2path,mir2path = mir2path,miRNA = "rno-miR-21-5p",targetgenes = "Smad7")
res_path1 <- get_miRTalk_pathway(object = obj_fat,gene2path = gene2path,mir2path = mir2path,miRNA = "rno-miR-223-3p",targetgenes = "Vim")

gene2path_rat <- gene2path[gene2path$species == "Rat",]
gene2path_rat1 <- gene2path_rat[grep("APOPTOSIS", gene2path_rat$term),]
gene2path_rat1$pair <- paste0(gene2path_rat1$term, "_", gene2path_rat1$type)
termnames <- unique(gene2path_rat1$term)
gene2path_rat2 <- gene2path_rat1[gene2path_rat1$term == "APOPTOSIS",]

genesets_apop <- list(set1 = gene2path_rat2[gene2path_rat2$type == "KEGG",]$gene,
                      set2 = gene2path_rat2[gene2path_rat2$type == "REACTOME",]$gene,
                      set3 = gene2path_rat2[gene2path_rat2$type == "WIKIPATHWAYS",]$gene)

obj_con_seurat <- CreateSeuratObject(obj_con_data)
obj_fat_seurat <- CreateSeuratObject(obj_fat_data)
obj_seurat_com <- merge(obj_con_seurat, obj_fat_seurat)
obj_seurat_com <- AddModuleScore(obj_seurat_com, genesets_apop)
obj_seurat_meta <- obj_seurat_com@meta.data
obj_meta_com <- rbind(obj_con_meta, obj_fat_meta)
obj_seurat_meta$celltype <- obj_meta_com$celltype
obj_seurat_meta$condition <- obj_meta_com$condition
obj_seurat_meta1 <- obj_seurat_meta[obj_seurat_meta$celltype == "hepatocyte",]

obj_seurat_meta1$Cluster3 <- obj_seurat_meta1$Cluster3-(min(obj_seurat_meta1$Cluster3))
x_normal <- obj_seurat_meta1[obj_seurat_meta1$condition == "control",]$Cluster3
x_tumor <- obj_seurat_meta1[obj_seurat_meta1$condition == "fat",]$Cluster3
data_plot <- data.frame(x = c("normal", "tumor"), y = c(mean(x_normal), mean(x_tumor)),
                        ymin = c(mean(x_normal) - sd(x_normal)/sqrt(length(x_normal)), mean(x_tumor) - sd(x_tumor)/sqrt(length(x_tumor))),
                        ymax = c(mean(x_normal) + sd(x_normal)/sqrt(length(x_normal)), mean(x_tumor) + sd(x_tumor)/sqrt(length(x_tumor))),stringsAsFactors = F)
p1 <- ggplot(data_plot)+geom_bar(aes(x = x, y=y),stat = "identity")+geom_errorbar(aes(x = x, ymin = ymin, ymax = ymax))

cci_fdl1 <- cci_fdl[cci_fdl$miR_gene %in% res_all1$mirgene,]
cci_fdl1 <- cci_fdl1[,c("miRNA","target_gene")]
cci_fdl1 <- unique(cci_fdl1)
res_MiTI_xlsx <- data.frame(miRNA = unique(cci_injured_new1$miRNA), target_gene = "NO",stringsAsFactors = F)
res_MiTI_xlsx <- res_MiTI_xlsx[order(res_MiTI_xlsx$miRNA),]
for (i in 1:nrow(res_MiTI_xlsx)) {
  d1 <- cci_injured_new1[cci_injured_new1$miRNA == res_MiTI_xlsx$miRNA[i],]
  d1 <- unique(d1$target_gene)
  d2 <- d1[1]
  if (length(d1) > 1) {
    for (j in 2:length(d1)) {
      d2 <- paste0(d2, ", ", d1[j])
    }
  }
  res_MiTI_xlsx$target_gene[i] <- d2
}
