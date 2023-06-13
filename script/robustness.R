library(miRTalk)
# function
cell_sub <- function(x,y){
  celltype <- as.data.frame(table(x[,2]),stringsAsFactors = F)
  celltype$num <- round(celltype$Freq*y)
  res_sub <- NULL
  for (i in 1:nrow(celltype)) {
    celltype1 <- x[x[,2] == celltype$Var1[i],]
    set.seed(i)
    cell_id <- sample(x = 1:celltype$Freq[i],size = celltype$num[i])
    celltype1 <- celltype1[cell_id,]
    res_sub <- rbind(res_sub, celltype1)
  }
  return(res_sub)
}
rename_chr_name <- function(x){
  x <- strsplit(x, split = " ")
  x_new <- NULL
  for (i in 1:length(x)) {
    x1 <- x[[i]]
    if (length(x1) > 1) {
      x2 <- x1[1]
      for (j in 2:length(x1)) {
        x2 <- paste(x2, x1[j], sep = "_")
      }
      x1 <- x2
    }
    x_new <- c(x_new, x1)
  }
  x <- strsplit(x_new, split = "-")
  x_new <- NULL
  for (i in 1:length(x)) {
    x1 <- x[[i]]
    if (length(x1) > 1) {
      x2 <- x1[1]
      for (j in 2:length(x1)) {
        x2 <- paste(x2, x1[j], sep = "_")
      }
      x1 <- x2
    }
    x_new <- c(x_new, x1)
  }
  return(x_new)
}
res_sub_all <- list()
# bladder cancer
d1 <- Read10X_h5("BLCA_GSE130001_expression.h5")
d2 <- read.table('BLCA_GSE130001_CellMetainfo_table.tsv',header = T,sep = '\t',stringsAsFactors = F)
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
sampleid <- unique(d2$Sample)
colnames(d1) <- rename_chr_name(colnames(d1))
d2$Cell <- rename_chr_name(d2$Cell)
d2$Celltype..major.lineage. <- rename_chr_name(d2$Celltype..major.lineage.)
res_tmp <- NULL
res_mir2tar <- data.frame()
for (i in 1:length(sampleid)) {
  a2 <- d2[d2$Sample == sampleid[i],]
  # 1.0
  print(paste0("bladder_",i,"_1.0"))
  a3 <- a2[,c("Cell","Celltype..major.lineage.")]
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]), species = "Human",if_normalize = F)
  obj@data$data <- a1
  obj <- find_miRNA(object = obj,mir_info = mir_info)
  obj <- find_hvtg(object = obj)
  obj <- find_miRTalk(obj, mir2tar = mir2tar)
  obj <- obj@cci
  obj$sample <- i
  obj$sample_size <- 1.0
  res_tmp <- rbind(res_tmp, obj)
  # 0.9
  print(paste0("bladder_",i,"_0.9"))
  a3 <- cell_sub(x = a2[,c("Cell","Celltype..major.lineage.")],y = 0.9)
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]), species = "Human",if_normalize = F)
  obj@data$data <- a1
  obj <- find_miRNA(object = obj,mir_info = mir_info)
  obj <- find_hvtg(object = obj)
  obj <- find_miRTalk(obj, mir2tar = mir2tar)
  obj <- obj@cci
  obj$sample <- i
  obj$sample_size <- 0.9
  res_tmp <- rbind(res_tmp, obj)
  # 0.8
  print(paste0("bladder_",i,"_0.8"))
  a3 <- cell_sub(x = a2[,c("Cell","Celltype..major.lineage.")],y = 0.8)
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]), species = "Human",if_normalize = F)
  obj@data$data <- a1
  obj <- find_miRNA(object = obj,mir_info = mir_info)
  obj <- find_hvtg(object = obj)
  obj <- find_miRTalk(obj, mir2tar = mir2tar)
  obj <- obj@cci
  obj$sample <- i
  obj$sample_size <- 0.8
  res_tmp <- rbind(res_tmp, obj)
  # 0.7
  print(paste0("bladder_",i,"_0.7"))
  a3 <- cell_sub(x = a2[,c("Cell","Celltype..major.lineage.")],y = 0.7)
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]), species = "Human",if_normalize = F)
  obj@data$data <- a1
  obj <- find_miRNA(object = obj,mir_info = mir_info)
  obj <- find_hvtg(object = obj)
  obj <- find_miRTalk(obj, mir2tar = mir2tar)
  obj <- obj@cci
  obj$sample <- i
  obj$sample_size <- 0.7
  res_tmp <- rbind(res_tmp, obj)
  # 0.6
  print(paste0("bladder_",i,"_0.6"))
  a3 <- cell_sub(x = a2[,c("Cell","Celltype..major.lineage.")],y = 0.6)
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]), species = "Human",if_normalize = F)
  obj@data$data <- a1
  obj <- find_miRNA(object = obj,mir_info = mir_info)
  obj <- find_hvtg(object = obj)
  obj <- find_miRTalk(obj, mir2tar = mir2tar)
  obj <- obj@cci
  obj$sample <- i
  obj$sample_size <- 0.6
  res_tmp <- rbind(res_tmp, obj)
  # 0.5
  print(paste0("bladder_",i,"_0.5"))
  a3 <- cell_sub(x = a2[,c("Cell","Celltype..major.lineage.")],y = 0.5)
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]), species = "Human",if_normalize = F)
  obj@data$data <- a1
  obj <- find_miRNA(object = obj,mir_info = mir_info)
  obj <- find_hvtg(object = obj)
  obj <- find_miRTalk(obj, mir2tar = mir2tar)
  obj <- obj@cci
  obj$sample <- i
  obj$sample_size <- 0.5
  res_tmp <- rbind(res_tmp, obj)
}

# chol
load("sc_data_raw.rda")
load("sc_meta.rda")
d1 <- sc_data_raw
d2 <- sc_meta
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
d2[d2$Specimen %in% c("ICC24_1", "ICC24_2"),]$Specimen <- "ICC24"
sampleid <- unique(d2$Specimen)

res_tmp <- NULL
for (i in 1:length(sampleid)) {
  a2 <- d2[d2$Specimen == sampleid[i],]
  # 1.0
  print(paste0("chol_",i,"_1.0"))
  a3 <- a2[,c("Id","Celltype")]
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]), species = "Human")
  obj <- find_miRNA(object = obj,mir_info = mir_info)
  obj <- find_hvtg(object = obj)
  obj <- find_miRTalk(obj, mir2tar = mir2tar)
  obj <- obj@cci
  obj$sample <- i
  obj$sample_size <- 1.0
  res_tmp <- rbind(res_tmp, obj)
  # 0.9
  print(paste0("chol_",i,"_0.9"))
  a3 <- cell_sub(x = a2[,c("Id","Celltype")],y = 0.9)
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]), species = "Human")
  obj <- find_miRNA(object = obj,mir_info = mir_info)
  obj <- find_hvtg(object = obj)
  obj <- find_miRTalk(obj, mir2tar = mir2tar)
  obj <- obj@cci
  obj$sample <- i
  obj$sample_size <- 0.9
  res_tmp <- rbind(res_tmp, obj)
  # 0.8
  print(paste0("chol_",i,"_0.8"))
  a3 <- cell_sub(x = a2[,c("Id","Celltype")],y = 0.8)
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]), species = "Human")
  obj <- find_miRNA(object = obj,mir_info = mir_info)
  obj <- find_hvtg(object = obj)
  obj <- find_miRTalk(obj, mir2tar = mir2tar)
  obj <- obj@cci
  obj$sample <- i
  obj$sample_size <- 0.8
  res_tmp <- rbind(res_tmp, obj)
  # 0.7
  print(paste0("chol_",i,"_0.7"))
  a3 <- cell_sub(x = a2[,c("Id","Celltype")],y = 0.7)
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]), species = "Human")
  obj <- find_miRNA(object = obj,mir_info = mir_info)
  obj <- find_hvtg(object = obj)
  obj <- find_miRTalk(obj, mir2tar = mir2tar)
  obj <- obj@cci
  obj$sample <- i
  obj$sample_size <- 0.7
  res_tmp <- rbind(res_tmp, obj)
  # 0.6
  print(paste0("chol_",i,"_0.6"))
  a3 <- cell_sub(x = a2[,c("Id","Celltype")],y = 0.6)
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]), species = "Human")
  obj <- find_miRNA(object = obj,mir_info = mir_info)
  obj <- find_hvtg(object = obj)
  obj <- find_miRTalk(obj, mir2tar = mir2tar)
  obj <- obj@cci
  obj$sample <- i
  obj$sample_size <- 0.6
  res_tmp <- rbind(res_tmp, obj)
  # 0.5
  print(paste0("chol_",i,"_0.5"))
  a3 <- cell_sub(x = a2[,c("Id","Celltype")],y = 0.5)
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]), species = "Human")
  obj <- find_miRNA(object = obj,mir_info = mir_info)
  obj <- find_hvtg(object = obj)
  obj <- find_miRTalk(obj, mir2tar = mir2tar)
  obj <- obj@cci
  obj$sample <- i
  obj$sample_size <- 0.5
  res_tmp <- rbind(res_tmp, obj)
}

# ovarian
d1 <- readRDS('GSE165897_UMIcounts_HGSOC.rds')
d2 <- readRDS('GSE165897_meta.rds')
d2 <- d2[d2$treatment_phase == "treatment-naive",]
d1 <- d1[,d2$cell]
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
sampleid <- unique(d2$patient_id)
sampleid <- sampleid[c(1,2,3,5,8,9)]
res_tmp <- NULL
for (i in 1:length(sampleid)) {
  a2 <- d2[d2$patient == sampleid[i],]
  # 1.0
  print(paste0("ovarian_",i,"_1.0"))
  a3 <- a2[,c("cell","cell_subtype")]
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]), species = "Human")
  obj <- find_miRNA(object = obj,mir_info = mir_info)
  obj <- find_hvtg(object = obj)
  obj <- find_miRTalk(obj, mir2tar = mir2tar)
  obj <- obj@cci
  obj$sample <- i
  obj$sample_size <- 1.0
  res_tmp <- rbind(res_tmp, obj)
  # 0.9
  print(paste0("ovarian_",i,"_0.9"))
  a3 <- cell_sub(x = a2[,c("cell","cell_subtype")],y = 0.9)
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]), species = "Human")
  obj <- find_miRNA(object = obj,mir_info = mir_info)
  obj <- find_hvtg(object = obj)
  obj <- find_miRTalk(obj, mir2tar = mir2tar)
  obj <- obj@cci
  obj$sample <- i
  obj$sample_size <- 0.9
  res_tmp <- rbind(res_tmp, obj)
  # 0.8
  print(paste0("ovarian_",i,"_0.8"))
  a3 <- cell_sub(x = a2[,c("cell","cell_subtype")],y = 0.8)
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]), species = "Human")
  obj <- find_miRNA(object = obj,mir_info = mir_info)
  obj <- find_hvtg(object = obj)
  obj <- find_miRTalk(obj, mir2tar = mir2tar)
  obj <- obj@cci
  obj$sample <- i
  obj$sample_size <- 0.8
  res_tmp <- rbind(res_tmp, obj)
  # 0.7
  print(paste0("ovarian_",i,"_0.7"))
  a3 <- cell_sub(x = a2[,c("cell","cell_subtype")],y = 0.7)
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]), species = "Human")
  obj <- find_miRNA(object = obj,mir_info = mir_info)
  obj <- find_hvtg(object = obj)
  obj <- find_miRTalk(obj, mir2tar = mir2tar)
  obj <- obj@cci
  obj$sample <- i
  obj$sample_size <- 0.7
  res_tmp <- rbind(res_tmp, obj)
  # 0.6
  print(paste0("ovarian_",i,"_0.6"))
  a3 <- cell_sub(x = a2[,c("cell","cell_subtype")],y = 0.6)
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]), species = "Human")
  obj <- find_miRNA(object = obj,mir_info = mir_info)
  obj <- find_hvtg(object = obj)
  obj <- find_miRTalk(obj, mir2tar = mir2tar)
  obj <- obj@cci
  obj$sample <- i
  obj$sample_size <- 0.6
  res_tmp <- rbind(res_tmp, obj)
  # 0.5
  print(paste0("ovarian_",i,"_0.5"))
  a3 <- cell_sub(x = a2[,c("cell","cell_subtype")],y = 0.5)
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]), species = "Human")
  obj <- find_miRNA(object = obj,mir_info = mir_info)
  obj <- find_hvtg(object = obj)
  obj <- find_miRTalk(obj, mir2tar = mir2tar)
  obj <- obj@cci
  obj$sample <- i
  obj$sample_size <- 0.5
  res_tmp <- rbind(res_tmp, obj)
}
