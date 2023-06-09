library(miRTalk)
library(Seurat)

res_sub_all <- list()
# bladder cancer
d1 <- Read10X_h5("BLCA_GSE130001_expression.h5")
d2 <- read.table('BLCA_GSE130001_CellMetainfo_table.tsv',header = T,sep = '\t',stringsAsFactors = F)
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
sampleid <- unique(d2$Sample)

res_tmp <- NULL
res_mir2tar <- data.frame()
for (i in 1:length(sampleid)) {
  a2 <- d2[d2$Sample == sampleid[i],]
  # 1.0
  print(paste0("bladder_",i,"_1.0"))
  a3 <- a2[,c("Cell","Celltype..major.lineage.")]
  a1 <- d1[,a3[,1]]
  obj <- create_miRTalk(sc_data = a1, sc_celltype = as.character(a3[,2]), species = "Human",if_normalize = F)
  colnames(a1) <- colnames(obj@data$data)
  obj <- find_miRNA(object = obj,mir_info = mir_info)
  obj <- find_hvtg(object = obj)
  obj <- find_miRTalk(obj, mir2tar = mir2tar, use_n_cores = 16)
  obj <- obj@cci
  obj$sample <- sampleid[i]
  obj$sample_size <- 1.0
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
  obj <- find_miRTalk(obj, mir2tar = mir2tar, use_n_cores = 16)
  obj <- obj@cci
  obj$sample <- sampleid[i]
  obj$sample_size <- 1.0
  res_tmp <- rbind(res_tmp, obj)
}
# ovarian
d1 <- readRDS('GSE165897_UMIcounts_HGSOC.rds')
d2 <- readRDS('GSE165897_meta.rds')
d2 <- d2[d2$treatment_phase == "treatment-naive",]
d1 <- d1[,d2$cell]
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
sampleid <- unique(d2$patient_id)

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
  obj <- find_miRTalk(obj, mir2tar = mir2tar, use_n_cores = 16)
  obj <- obj@cci
  obj$sample <- sampleid[i]
  obj$sample_size <- 1.0
  res_tmp <- rbind(res_tmp, obj)
}

d1 <- readRDS('GSE165897_UMIcounts_HGSOC.rds')
d2 <- readRDS('GSE165897_meta.rds')
d2 <- d2[d2$treatment_phase == "post-NACT",]
d1 <- d1[,d2$cell]
d1 <- rev_gene(d1,data_type = 'count',species = 'Human',geneinfo = geneinfo)
sampleid <- unique(d2$patient_id)

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
  obj <- find_miRTalk(obj, mir2tar = mir2tar, use_n_cores = 16)
  obj <- obj@cci
  obj$sample <- sampleid[i]
  obj$sample_size <- 1.0
  res_tmp <- rbind(res_tmp, obj)
}
