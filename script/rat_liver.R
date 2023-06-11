library(miRTalk)
library(Nebulosa)
load('sc_data_raw.rda')
load('sc_meta.rda')
load('sc_meta_specimen.rda')
sc_data_raw <- rev_gene(sc_data_raw,data_type = 'count',species = 'Rat',geneinfo = geneinfo)
res_condition <- list()
samplename <- unique(sc_meta_specimen$Status)
for (i in 1:length(samplename)) {
  print(i)
  a3 <- sc_meta_specimen[sc_meta_specimen$Status == samplename[i],]
  a2 <- sc_meta[sc_meta$Specimen %in% a3$Specimen, ]
  a1 <- sc_data_raw[,a2$Id]
  obj <- create_miRTalk(sc_data = a1,sc_celltype = a2$Name_Celltype,species = "Rat",if_normalize = F)
  obj <- find_miRNA(obj, mir_info = mir_info)
  obj <- find_hvtg(obj)
  obj <- find_miRTalk(obj, mir2tar = mir2tar,use_n_cores = 8)
  res_condition[[i]] <- obj
}
names(res_condition) <- samplename

obj_con <- res_condition[["Normal"]]
obj_fat <- res_condition[["NAFLD (high-fat diet, 8 weeks)"]]

plot_miRTalk_circle(object = obj_con)
plot_miRTalk_circle(object = obj_fat)

# analysis
celltype = NULL
miRNA = NULL
text_size = 10
if_horizontal = TRUE
object <- obj_con
# check input
if (!is(object, "miRTalk")) {
  stop("Invalid class for object: must be 'miRTalk'!")
}
cci <- object@cci
if (nrow(cci) == 0) {
  stop("No cci found in object!")
}
celltype_sender <- unique(cci$celltype_sender)
miR_name <- unique(cci$miRNA)
if (!is.null(celltype[1])) {
  if (all(celltype %in% celltype_sender)) {
    stop("Please input the right celltype name!")
  }
  cci <- cci[cci$celltype_sender %in% celltype, ]
}
if (!is.null(miRNA[1])) {
  if (!all(miRNA %in% miR_name)) {
    stop("Please input the right miRNA name!")
  }
  cci <- cci[cci$miRNA %in% miRNA, ]
}
cci <- unique(cci[, c("celltype_sender", "miRNA", "miRNA_activity")])
if (length(unique(cci$celltype_sender)) < 2 | length(unique(cci$miRNA)) < 2) {
  stop("Limited celltype_sender and miRNA to plot!")
}
cci <- reshape2::dcast(data = cci, formula = celltype_sender ~ miRNA, value.var = "miRNA_activity", fill = 0)
rownames(cci) <- cci$celltype_sender
cci <- cci[, -1]
if (!if_horizontal) {
  cci <- as.data.frame(t(cci))
}

cci_new <- as.matrix(cci)

object <- obj_fat
# check input
if (!is(object, "miRTalk")) {
  stop("Invalid class for object: must be 'miRTalk'!")
}
cci <- object@cci
if (nrow(cci) == 0) {
  stop("No cci found in object!")
}
celltype_sender <- unique(cci$celltype_sender)
miR_name <- unique(cci$miRNA)
if (!is.null(celltype[1])) {
  if (all(celltype %in% celltype_sender)) {
    stop("Please input the right celltype name!")
  }
  cci <- cci[cci$celltype_sender %in% celltype, ]
}
if (!is.null(miRNA[1])) {
  if (!all(miRNA %in% miR_name)) {
    stop("Please input the right miRNA name!")
  }
  cci <- cci[cci$miRNA %in% miRNA, ]
}
cci <- unique(cci[, c("celltype_sender", "miRNA", "miRNA_activity")])
if (length(unique(cci$celltype_sender)) < 2 | length(unique(cci$miRNA)) < 2) {
  stop("Limited celltype_sender and miRNA to plot!")
}
cci <- reshape2::dcast(data = cci, formula = celltype_sender ~ miRNA, value.var = "miRNA_activity", fill = 0)
rownames(cci) <- cci$celltype_sender
cci <- cci[, -1]
if (!if_horizontal) {
  cci <- as.data.frame(t(cci))
}
cci <- as.matrix(cci)
cci <- rbind(cci_new, cci)
heatmaply::heatmaply(x = as.matrix(cci),dendrogram = "none", xlab = "", ylab = "", main = "miRNA activity in senders", margins = c(60,100,40,20),
                     grid_color = "white", grid_width = 0.00001, titleX = FALSE, branches_lwd = 0.1, fontsize_row = text_size,
                     fontsize_col = text_size, labCol = colnames(cci), labRow = rownames(cci), heatmap_layers = theme(axis.line=element_blank()))

data_con <- obj_con@data$data
data_fat <- obj_fat@data$data
meta_con <- obj_con@meta
meta_fat <- obj_fat@meta

summary(data_con["Mir21",meta_con[meta_con$celltype %in% c("dendritic cell", "granulocyte",
                                                           "Kupffer cell", "macrophage", "monocyte", "natural killer cell"),]$cell])
a1 <- data_con["Mir21",meta_con[meta_con$celltype %in% c("dendritic cell", "granulocyte",
                                                         "Kupffer cell", "macrophage", "monocyte", "natural killer cell"),]$cell]
summary(data_fat["Mir21",meta_fat[meta_fat$celltype %in% c("dendritic cell", "granulocyte",
                                                           "Kupffer cell", "macrophage", "monocyte", "natural killer cell"),]$cell])
a2 <- data_fat["Mir21",meta_fat[meta_fat$celltype %in% c("dendritic cell", "granulocyte",
                                                         "Kupffer cell", "macrophage", "monocyte", "natural killer cell"),]$cell]
res1 <- data.frame(gene = c(a1,a2), group = c(rep("normal",length(a1)), rep("fatty",length(a2))),stringsAsFactors = F)

summary(data_con["Smad7",meta_con[meta_con$celltype == "hepatocyte",]$cell])

a1 <- data_con["Smad7",meta_con[meta_con$celltype == "hepatocyte",]$cell]

summary(data_fat["Smad7",meta_fat[meta_fat$celltype == "hepatocyte",]$cell])

a2 <- data_fat["Smad7",meta_fat[meta_fat$celltype == "hepatocyte",]$cell]

res1 <- data.frame(gene = c(a1,a2), group = c(rep("normal",length(a1)), rep("fatty",length(a2))),stringsAsFactors = F)

a1 <- data_con["Col1a1",meta_con[meta_con$celltype == "hepatocyte",]$cell]
a2 <- data_fat["Col1a1",meta_fat[meta_fat$celltype == "hepatocyte",]$cell]
res1 <- data.frame(gene = c(a1,a2), group = c(rep("normal",length(a1)), rep("fatty",length(a2))),stringsAsFactors = F)

a1 <- data_con["Mir223",meta_con[meta_con$celltype == "granulocyte",]$cell]
a2 <- data_fat["Mir223",meta_fat[meta_fat$celltype == "granulocyte",]$cell]
res1 <- data.frame(gene = c(a1,a2), group = c(rep("normal",length(a1)), rep("fatty",length(a2))),stringsAsFactors = F)

a1 <- data_con["Vim",meta_con[meta_con$celltype == "hepatocyte",]$cell]
a2 <- data_fat["Vim",meta_fat[meta_fat$celltype == "hepatocyte",]$cell]
res1 <- data.frame(gene = c(a1,a2), group = c(rep("normal",length(a1)), rep("fatty",length(a2))),stringsAsFactors = F)

plot_miR2tar_chord(object = obj_con, celltype_sender = c("DC", "granulocyte", "Kupffer", "macrophage", "monocyte"),
                   celltype_receiver = "hepatocyte", miRNA = "rno-miR-21-5p", celltype_color = "NO")
plot_miR2tar_chord(object = obj_con, celltype_sender = c("DC", "granulocyte", "Kupffer", "macrophage", "monocyte"),
                   celltype_receiver = "hepatocyte", miRNA = "rno-miR-21-5p",
                   celltype_color = c("#16A28D", "#00468B", "#E2746E", "#D2070F", "#688C31", "#5178A8"))
plot_miR2tar_chord(object = obj_con, celltype_sender = c("DC", "granulocyte", "Kupffer", "macrophage", "monocyte"),
                   celltype_receiver = "hepatocyte", miRNA = "rno-miR-21-5p",
                   celltype_color = "NO")
plot_miR2tar_chord(object = obj_fat, celltype_sender = c("dendritic_cell", "granulocyte", "Kupffer_cell", "macrophage", "monocyte"),
                   celltype_receiver = "hepatocyte", miRNA = "rno-miR-21-5p",
                   celltype_color = c("#16A28D", "#00468B", "#E2746E", "#D2070F", "#688C31", "#5178A8"))

