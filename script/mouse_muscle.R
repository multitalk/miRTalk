library(miRTalk)
# https://doi.org/10.1038/s41587-022-01517-6
sc_meta <- readRDS("st_meta.rds")
sc_data <- readRDS("st_data.rds")
obj_condition <- list()
samplename <- unique(sc_meta$condition)
samplename <- samplename[-1]
for (i in 1:length(samplename)) {
  print(i)
  a2 <- sc_meta[sc_meta$condition %in% samplename[i], ]
  a1 <- sc_data[,a2$cell]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = a2$cluster,
                        species = "Mouse", condition = a2$condition, evbiog = evbiog, risc = risc)
  obj <- find_hvtg(object = obj)
  obj <- find_miRNA(object = obj,mir_info = mir_info,mir2tar = mir2tar)
  obj <- find_miRTalk(obj, use_n_cores = 16)
  res_condition[[i]] <- obj
}
names(obj_condition) <- samplename

# uninjured, dpi_2, dpi_5, dpi_7
# Myofiber, Injury locus, Injury border

# data plot

st_meta1 <- st_meta[st_meta$condition == "uninjured",]
p1 <- ggplot()+geom_point(data = st_meta1, aes(x = x, y = y, color = cluster), size = 4)+
  scale_color_manual(values = "#2e6198")+NoLegend()
ggsave(filename = "plot_uninjured_data.pdf",plot = p1,device = "pdf",width = 8,height = 6)

st_meta1 <- st_meta[st_meta$condition == "dpi_2",]
p1 <- ggplot()+geom_point(data = st_meta1, aes(x = x, y = y, color = cluster), size = 4)+
  scale_color_manual(values = c("#d8a837","#c42a2e","#2e6198"))+NoLegend()
ggsave(filename = "plot_dpi_2_data.pdf",plot = p1,device = "pdf",width = 8,height = 6)

st_meta1 <- st_meta[st_meta$condition == "dpi_5",]
p1 <- ggplot()+geom_point(data = st_meta1, aes(x = x, y = y, color = cluster), size = 4)+
  scale_color_manual(values = c("#d8a837","#c42a2e","#2e6198"))+NoLegend()
ggsave(filename = "plot_dpi_5_data.pdf",plot = p1,device = "pdf",width = 8,height = 6)


st_meta1 <- st_meta[st_meta$condition == "dpi_7",]
p1 <- ggplot()+geom_point(data = st_meta1, aes(x = x, y = y, color = cluster), size = 4)+
  scale_color_manual(values = c("#d8a837","#c42a2e","#2e6198"))+NoLegend()
ggsave(filename = "plot_dpi_7_data.pdf",plot = p1,device = "pdf",width = 8,height = 6)

#miR_score
obj2 <- obj_condition[["dpi_2"]]
obj5 <- obj_condition[["dpi_5"]]
obj7 <- obj_condition[["dpi_7"]]

obj2_cci <- obj2@cci
obj5_cci <- obj5@cci
obj7_cci <- obj7@cci

obj_cci <- rbind(obj2_cci, obj5_cci, obj7_cci)
obj_cci_raw <- obj@cci
obj@cci <- obj_cci

cci <- unique(obj_cci[, c("celltype_sender", "miRNA", "EVmiR_score","condition")])
cci$celltype_sender <- paste0(cci$condition, "_", cci$celltype_sender)
cci <- cci[,-4]
cci <- reshape2::dcast(data = cci, formula = celltype_sender ~ miRNA, fun.aggregate = mean, value.var = "EVmiR_score", fill = 0)
rownames(cci) <- cci$celltype_sender
cci <- cci[, -1]
cci <- cci[c(2, 5, 1, 4, 3, 6),]
heat_col <- viridis::viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "D")
heatmaply::heatmaply(x = as.matrix(cci), colors = heat_col, dendrogram = "none", xlab = "miRNA", ylab = "Senders", main = "EVmiR_score in senders", margins = c(60,100,40,20),
                     titleX = FALSE, branches_lwd = 0.1, fontsize_row = 10, fontsize_col = 10, labCol = colnames(cci), labRow = rownames(cci),
                     heatmap_layers = theme(axis.line=element_blank()), label_names = c("Sender","miRNA","EVmiR_score"))

st_obj <- CreateSeuratObject(st_data)
st_obj <- NormalizeData(st_obj)
st_obj <- ScaleData(st_obj)
st_ndata <- st_obj[["RNA"]]@data
st_sdata <- st_obj[["RNA"]]@scale.data
st_meta$mir145a <- st_ndata["Mir145a",]

a1 <- st_meta[st_meta$condition == "uninjured",]
a2 <- st_meta[st_meta$condition == "dpi_2",]
a3 <- st_meta[st_meta$condition == "dpi_5",]
a4 <- st_meta[st_meta$condition == "dpi_7",]
boxplot(a1$mir145a,a2$mir145a,a3$mir145a,a4$mir145a)

# Myofiber, Injury locus, Injury border
a2_myofiber <- a2[a2$cluster == "Myofiber",]
a2_locus <- a2[a2$cluster == "Injury locus",]
a2_border <- a2[a2$cluster == "Injury border",]

a3_myofiber <- a3[a3$cluster == "Myofiber",]
a3_locus <- a3[a3$cluster == "Injury locus",]
a3_border <- a3[a3$cluster == "Injury border",]

a4_myofiber <- a4[a4$cluster == "Myofiber",]
a4_border <- a4[a4$cluster == "Injury border",]

boxplot(a1$mir145a, a2_locus$mir145a, a2_border$mir145a, a2_myofiber$mir145a,a3_locus$mir145a, a3_border$mir145a, a3_myofiber$mir145a)

cci <- unique(obj_cci[, c("celltype_sender","celltype_receiver", "miRNA", "target_gene", "condition")])
cci$celltype_sender <- paste0(cci$condition, "_", cci$celltype_sender)
cci$celltype_receiver <- paste0(cci$condition, "_", cci$celltype_receiver)
cci$pair <- paste0(cci$miRNA, "_", cci$target_gene)
cci <- cci[,-c(3,4,5)]
cci <- reshape2::dcast(data = cci, formula = celltype_sender ~ celltype_receiver, value.var = "pair", fill = 0)
rownames(cci) <- cci$celltype_sender
cci <- cci[, -1]
cci <- cci[c(2, 5, 1, 4, 3, 6),]
cci <- cci[,c(2,1,3, 5,4,6,7,8)]

heat_col <- viridis::viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "D")
heatmaply::heatmaply(x = as.matrix(cci), colors = heat_col, dendrogram = "none", xlab = "miRNA", ylab = "Senders", main = "EVmiR_score in senders", margins = c(60,100,40,20),
                     titleX = FALSE, branches_lwd = 0.1, fontsize_row = 10, fontsize_col = 10, labCol = colnames(cci), labRow = rownames(cci),
                     heatmap_layers = theme(axis.line=element_blank()), label_names = c("Sender","miRNA","EVmiR_score"))


# Bone marrow; Liver; Skin; Spleen

mirnames <- unique(obj_cci$miR_gene)
mir_change <- data.frame(mirgene = mirnames, stringsAsFactors = F)
mir_change$exp_day0 <- 0
mir_change$exp_day2 <- 0
mir_change$wilcol_p <- 1
celltypes1 <- unique(st_meta$cluster)
celltypes1 <- celltypes1[order(celltypes1)]
celltypes2 <- unique(c(obj_cci$celltype_sender, obj_cci$celltype_receiver))
celltypes2 <- celltypes2[order(celltypes2)]
celltypes <- data.frame(celltype_seurat = celltypes1, celltype_mirtalk = celltypes2, stringsAsFactors = F)
res_all <- data.frame()
for (k in 1:nrow(celltypes)) {
  mir_change1 <- mir_change
  celltype_seurat <- celltypes$celltype_seurat[k]
  celltype_mirtalk <- celltypes$celltype_mirtalk[k]
  for (i in 1:nrow(mir_change1)) {
    st_meta$gene <- as.numeric(st_ndata[mir_change1$mirgene[i],])
    obj_meta1 <- st_meta[st_meta$cluster == "Myofiber" & st_meta$condition == "uninjured",]
    obj_meta2 <- st_meta[st_meta$cluster == celltype_seurat & st_meta$condition == "dpi_2",]
    mir_change1$exp_day0[i] <- mean(obj_meta1$gene)
    mir_change1$exp_day2[i] <- mean(obj_meta2$gene)
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
obj_cci1 <- obj_cci[obj_cci$miR_gene %in% res_all$mirgene,]
obj_cci1 <- unique(obj_cci1[,c("miR_gene","miRNA")])
colnames(obj_cci1)[1] <- "mirgene"
res_all1 <- merge(res_all,obj_cci1)
mirnames <- unique(res_all1$miRNA)

bulk_data <- readRDS("rawdata.rds")
bulk_data1 <- bulk_data
bulk_data1[,c(2:9)] <- log2(bulk_data1[,c(2:9)]+1)
bulk_data2 <- as.data.frame(bulk_data1)
bulk_data2$day0 <- rowMeans(bulk_data2[,c(2,3)])
bulk_data2$day2 <- rowMeans(bulk_data2[,c(4,5)])
bulk_data2$day5 <- rowMeans(bulk_data2[,c(6,7)])
bulk_data2$day7 <- rowMeans(bulk_data2[,c(8,9)])

bulk_data2 <- bulk_data2[bulk_data2$miRNA %in% mirnames,]
res_all1 <- res_all1[res_all1$miRNA %in% bulk_data2$miRNA,]
write.csv(res_all1, file = "res_miR_sc.csv")

# difference_sc
for (i in 1:nrow(res_all1)) {
  st_meta$gene <- as.numeric(st_ndata[res_all1$mirgene[i],])
  obj_meta1 <- st_meta[st_meta$cluster == "Myofiber" & st_meta$condition == "uninjured",]
  obj_meta2 <- st_meta[st_meta$cluster == res_all1$celltype[i] & st_meta$condition == "dpi_2",]
  x_normal <- obj_meta1$gene
  x_tumor <- obj_meta2$gene
  data_plot <- data.frame(x = c("normal", "tumor"), y = c(mean(x_normal), mean(x_tumor)),
                          ymin = c(mean(x_normal) - sd(x_normal)/sqrt(length(x_normal)), mean(x_tumor) - sd(x_tumor)/sqrt(length(x_tumor))),
                          ymax = c(mean(x_normal) + sd(x_normal)/sqrt(length(x_normal)), mean(x_tumor) + sd(x_tumor)/sqrt(length(x_tumor))),stringsAsFactors = F)
  p1 <- ggplot(data_plot)+geom_bar(aes(x = x, y=y),stat = "identity")+geom_errorbar(aes(x = x, ymin = ymin, ymax = ymax))
  ggsave(filename = paste0("",res_all1$celltype[i],"_",res_all1$mirgene[i], "_change.png"),plot = p1,device = "png",width = 6,height = 6)
}

# difference_bulk
for (i in 1:nrow(bulk_data2)) {
  a1 <- bulk_data2[i,]
  x_normal <- a1$day0
  x_tumor <- a1$day2
  data_plot <- data.frame(x = c("normal", "tumor"), y = c(mean(x_normal), mean(x_tumor)),stringsAsFactors = F)
  p1 <- ggplot(data_plot)+geom_bar(aes(x = x, y=y),stat = "identity")
  ggsave(filename = paste0("exp_",bulk_data2$miRNA[i], "_change.png"),plot = p1,device = "png",width = 6,height = 6)
}

st_meta$mir142a <- st_ndata["Mir142",]
st_meta$mir145a <- st_ndata["Mir145a",]
st_meta$Nmrk2 <- st_ndata["Nmrk2",]
st_meta$Ms4a6c <- st_ndata["Ms4a6c",]

library(viridis)

# Mir142
value_min <- min(st_meta$mir142a)
value_max <- max(st_meta$mir142a)
conditionnames <- unique(st_meta$condition)
for(i in 1:4){
  st_meta1 <- st_meta[st_meta$condition == conditionnames[i],]
  p1 <- ggplot()+geom_point(data = st_meta1, aes(x = x, y = y, color = mir142a), size = 4)+scale_color_viridis(limits = c(value_min,value_max),option = "C")+NoLegend()
  ggsave(filename = paste0("plot_", conditionnames[i], "_mir142a.pdf"),plot = p1,device = "pdf",width = 8,height = 6)
}

# Mir145
value_min <- min(st_meta$mir145a)
value_max <- max(st_meta$mir145a)
conditionnames <- unique(st_meta$condition)
for(i in 1:4){
  st_meta1 <- st_meta[st_meta$condition == conditionnames[i],]
  p1 <- ggplot()+geom_point(data = st_meta1, aes(x = x, y = y, color = mir145a), size = 4)+scale_color_viridis(limits = c(value_min,value_max))+NoLegend()
  ggsave(filename = paste0("plot_", conditionnames[i], "_mir145a.pdf"),plot = p1,device = "pdf",width = 8,height = 6)
}

# Nmrk2
value_min <- min(st_meta$Nmrk2)
value_max <- max(st_meta$Nmrk2)
conditionnames <- unique(st_meta$condition)
for(i in 1:4){
  st_meta1 <- st_meta[st_meta$condition == conditionnames[i],]
  p1 <- ggplot()+geom_point(data = st_meta1, aes(x = x, y = y, color = Nmrk2), size = 4)+scale_color_viridis(limits = c(value_min,value_max))+NoLegend()
  ggsave(filename = paste0("plot_", conditionnames[i], "_Nmrk2.pdf"),plot = p1,device = "pdf",width = 8,height = 6)
}

# Ms4a6c
value_min <- min(st_meta$Ms4a6c)
value_max <- max(st_meta$Ms4a6c)
conditionnames <- unique(st_meta$condition)
for(i in 1:4){
  st_meta1 <- st_meta[st_meta$condition == conditionnames[i],]
  p1 <- ggplot()+geom_point(data = st_meta1, aes(x = x, y = y, color = Ms4a6c), size = 4)+scale_color_viridis(limits = c(value_min,value_max))+NoLegend()
  ggsave(filename = paste0("plot_", conditionnames[i], "_Ms4a6c.pdf"),plot = p1,device = "pdf",width = 8,height = 6)
}

gene2path_mouse <- gene2path[gene2path$species == "Mouse",]
res_path_gene <- gene2path_mouse[gene2path_mouse$gene == "Nmrk2",]
res_path_gene$term <- tolower(res_path_gene$term)
mir2path_mouse <- mir2path
res_path_mir <- mir2path_mouse[mir2path_mouse$mir == "hsa-miR-142-5p",]
res_path_mir$term <- tolower(res_path_mir$term)

d1 <- res_path_mir[res_path_mir$term %in% res_path_gene$term,]

res_path1 <- get_miRTalk_pathway(object = obj,gene2path = gene2path,mir2path = mir2path,miRNA = "mmu-miR-142a-5p",targetgenes = "Nmrk2")

a1 <- geneinfo1[geneinfo1$term == "NAD METABOLISM, SIRTUINS AND AGING",]
a1 <- a1[a1$species == "Mouse",]$gene
a2 <- geneinfo1[geneinfo1$term == "NEGATIVE REGULATION OF MYOTUBE DIFFERENTIATION",]
a2 <- a2[a2$species == "Mouse",]$gene
geneset1 <- list(gene1 = a1, gene2 = a2)
st_obj <- AddModuleScore(st_obj,features = geneset1)
st_meta$genemodule1 <- st_obj$Cluster1
st_meta$genemodule2 <- st_obj$Cluster2



# miR-142a-5p:Nmrk2

res_plot <- data.frame()
# uninjured
d1_day <- st_meta[st_meta$condition == "dpi_2",]
a1 <- d1_day[d1_day$cluster == "Myofiber",]$mir142a
a2 <- d1_day[d1_day$cluster == "Myofiber",]$Nmrk2
data_plot <- data.frame(x = c("a_aMyofiber", "a_aMyofiber"), y = c(mean(a1), mean(a2)), celltype = c("mirna", "target"),
                        ymin = c(mean(a1) - sd(a1)/sqrt(length(a1)), mean(a2) - sd(a2)/sqrt(length(a2))),
                        ymax = c(mean(a1) + sd(a1)/sqrt(length(a1)), mean(a2) + sd(a2)/sqrt(length(a2))),stringsAsFactors = F)
res_plot <- rbind(res_plot, data_plot)

a1 <- d1_day[d1_day$cluster == "Injury border",]$mir142a
a2 <- d1_day[d1_day$cluster == "Injury border",]$Nmrk2
data_plot <- data.frame(x = c("a_border", "a_border"), y = c(mean(a1), mean(a2)), celltype = c("mirna", "target"),
                        ymin = c(mean(a1) - sd(a1)/sqrt(length(a1)), mean(a2) - sd(a2)/sqrt(length(a2))),
                        ymax = c(mean(a1) + sd(a1)/sqrt(length(a1)), mean(a2) + sd(a2)/sqrt(length(a2))),stringsAsFactors = F)
res_plot <- rbind(res_plot, data_plot)

a1 <- d1_day[d1_day$cluster == "Injury locus",]$mir142a
a2 <- d1_day[d1_day$cluster == "Injury locus",]$Nmrk2
data_plot <- data.frame(x = c("a_locus", "a_locus"), y = c(mean(a1), mean(a2)), celltype = c("mirna", "target"),
                        ymin = c(mean(a1) - sd(a1)/sqrt(length(a1)), mean(a2) - sd(a2)/sqrt(length(a2))),
                        ymax = c(mean(a1) + sd(a1)/sqrt(length(a1)), mean(a2) + sd(a2)/sqrt(length(a2))),stringsAsFactors = F)
res_plot <- rbind(res_plot, data_plot)

p1 <- ggplot(res_plot[res_plot$celltype == "mirna",])+geom_bar(aes(x = x, y=y),stat = "identity")+geom_errorbar(aes(x = x, ymin = ymin, ymax = ymax))

p1 <- ggplot(res_plot[res_plot$celltype == "target",])+geom_bar(aes(x = x, y=y),stat = "identity")+geom_errorbar(aes(x = x, ymin = ymin, ymax = ymax))


# miR-145a-5p:Ms4a6c

res_plot <- data.frame()
# uninjured
d1_day <- st_meta[st_meta$condition == "dpi_2",]
a1 <- d1_day[d1_day$cluster == "Myofiber",]$mir145a
a2 <- d1_day[d1_day$cluster == "Myofiber",]$Ms4a6c
data_plot <- data.frame(x = c("a_aMyofiber", "a_aMyofiber"), y = c(mean(a1), mean(a2)), celltype = c("mirna", "target"),
                        ymin = c(mean(a1) - sd(a1)/sqrt(length(a1)), mean(a2) - sd(a2)/sqrt(length(a2))),
                        ymax = c(mean(a1) + sd(a1)/sqrt(length(a1)), mean(a2) + sd(a2)/sqrt(length(a2))),stringsAsFactors = F)
res_plot <- rbind(res_plot, data_plot)

a1 <- d1_day[d1_day$cluster == "Injury border",]$mir145a
a2 <- d1_day[d1_day$cluster == "Injury border",]$Ms4a6c
data_plot <- data.frame(x = c("a_border", "a_border"), y = c(mean(a1), mean(a2)), celltype = c("mirna", "target"),
                        ymin = c(mean(a1) - sd(a1)/sqrt(length(a1)), mean(a2) - sd(a2)/sqrt(length(a2))),
                        ymax = c(mean(a1) + sd(a1)/sqrt(length(a1)), mean(a2) + sd(a2)/sqrt(length(a2))),stringsAsFactors = F)
res_plot <- rbind(res_plot, data_plot)

a1 <- d1_day[d1_day$cluster == "Injury locus",]$mir145a
a2 <- d1_day[d1_day$cluster == "Injury locus",]$Ms4a6c
data_plot <- data.frame(x = c("a_locus", "a_locus"), y = c(mean(a1), mean(a2)), celltype = c("mirna", "target"),
                        ymin = c(mean(a1) - sd(a1)/sqrt(length(a1)), mean(a2) - sd(a2)/sqrt(length(a2))),
                        ymax = c(mean(a1) + sd(a1)/sqrt(length(a1)), mean(a2) + sd(a2)/sqrt(length(a2))),stringsAsFactors = F)
res_plot <- rbind(res_plot, data_plot)

p1 <- ggplot(res_plot[res_plot$celltype == "mirna",])+geom_bar(aes(x = x, y=y),stat = "identity")+geom_errorbar(aes(x = x, ymin = ymin, ymax = ymax))

p1 <- ggplot(res_plot[res_plot$celltype == "target",])+geom_bar(aes(x = x, y=y),stat = "identity")+geom_errorbar(aes(x = x, ymin = ymin, ymax = ymax))


res_path1 <- get_miRTalk_pathway(object = obj,gene2path = gene2path,mir2path = mir2path,miRNA = "mmu-miR-145a-5p",targetgenes = "Ms4a6c")

mir2path1 <- mir2path[mir2path$mir == "hsa-miR-145-5p", ]
mir2path1$term <- tolower(mir2path1$term)
mir2path1 <- mir2path1[grep(pattern = "muscle", x = mir2path1$term),]

mir_terms <- unique(mir2path1$term)
mir_terms <- mir_terms[c(5,6,21,22,26,27,28)]
gene2path1 <- gene2path[gene2path$species == "Mouse",]
gene2path1$term <- tolower(gene2path1$term)
gene2path1 <- gene2path1[gene2path1$term %in% mir_terms,]

mir_terms <- unique(gene2path1$term)
gene2sets1 <- list()
for(i in 1:length(mir_terms)){
  gene2sets1[[i]] <- gene2path1[gene2path1$term == mir_terms[i],]$gene
}
st_obj <- AddModuleScore(st_obj, features = gene2sets1)
st_meta$cluster1 <- st_obj$Cluster1
st_meta$cluster2 <- st_obj$Cluster2
st_meta$cluster3 <- st_obj$Cluster3
st_meta$cluster4 <- st_obj$Cluster4
st_meta$cluster5 <- st_obj$Cluster5
st_meta$cluster6 <- st_obj$Cluster6
st_meta$cluster7 <- st_obj$Cluster7

# cluster 1 and 2 selected
st_meta$x <- st_meta$Ms4a6c/max(st_meta$Ms4a6c)
st_meta$y1 <- st_meta$cluster1 - min(st_meta$cluster1)
st_meta$y1 <- st_meta$y1/max(st_meta$y1)

st_meta$y2 <- st_meta$cluster2 - min(st_meta$cluster2)
st_meta$y2 <- st_meta$y2/max(st_meta$y2)

a1 <- correlation(st_meta[,c("x","y1")],method = "spearman")
a1$p
a1$r

a1 <- correlation(st_meta[,c("x","y2")],method = "spearman")
a1$p
a1$r

p1 <- ggplot()+geom_point(data = st_meta, aes(x = x, y = y1))+geom_smooth(data = st_meta,aes(x = x, y = y1),method = "lm")
ggsave(filename = paste0("plot_cor_Ms4a6c_1.pdf"),plot = p1,device = "pdf",width = 6,height = 6)

p1 <- ggplot()+geom_point(data = st_meta, aes(x = x, y = y2))+geom_smooth(data = st_meta,aes(x = x, y = y2),method = "lm")
ggsave(filename = paste0("plot_cor_Ms4a6c_2.pdf"),plot = p1,device = "pdf",width = 6,height = 6)


#res_circulating--------------------------------------------------------------------------
res_circulating <- get_miRTalk_circulating_score(obj)
res_circulating$celltype_receiver <- paste0(res_circulating$condition,"_",res_circulating$celltype_receiver)
res_circulating <- unique(res_circulating[,c("miRNA","score","celltype_receiver")])

res_cir_plot <- reshape2::dcast(data = res_circulating, formula = miRNA ~ celltype_receiver, value.var = "score", fun.aggregate = mean, fill = 0)
rownames(res_cir_plot) <- res_cir_plot$miRNA
res_cir_plot <- res_cir_plot[,-1]
heat_col <- viridis::viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "D")
heatmaply::heatmaply(x = as.matrix(res_cir_plot), colors = heat_col, limits = c(0,1),dendrogram = "none",  margins = c(60,100,40,20),
                     titleX = FALSE, branches_lwd = 0.1, fontsize_row = 10, fontsize_col = 10, labCol = colnames(res_cir_plot), labRow = rownames(res_cir_plot),
                     heatmap_layers = theme(axis.line=element_blank()))

res_tissues <- unique(res_circulating[ ,c("miRNA","tissue_TarBase")])
tissues_names <- res_tissues$tissue_TarBase
tissues_names <- strsplit(x = tissues_names, split = "; ")
tissues_names_all <- character()
for (i in 1:length(tissues_names)) {
  tissues_names_all <- c(tissues_names_all, tissues_names[[i]])
}
tissues_names_all <- unique(tissues_names_all)
tissues_names_all <- tissues_names_all[order(tissues_names_all)]
tissues_names_all <- tissues_names_all[-5]


#---all
tissues_names <- mir_info$tissue_TarBase
tissues_names <- strsplit(x = tissues_names, split = "; ")
tissues_names_all <- character()
for (i in 1:length(tissues_names)) {
  tissues_names_all <- c(tissues_names_all, tissues_names[[i]])
}
tissues_names_all <- unique(tissues_names_all)
tissues_names_all <- tissues_names_all[order(tissues_names_all)]
tissues_names_all <- tissues_names_all[-5]

#-----MiTIs
res_all1 <- res_all[res_all$mirgene %in% c("Mir142", "Mir145a", "Mir15b", "Mir221", "Mir674"),]

cci_fdl1 <- obj2_cci[obj2_cci$miR_gene %in% res_all1$mirgene,]
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